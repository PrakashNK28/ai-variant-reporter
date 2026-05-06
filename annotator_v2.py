# app.py
# SpectralG — AI-Powered Clinical Variant Interpretation
# Version 2.0 — Universal input + Full ACMG engine + Clinical PDF output
# Upgrades: PP5 removed | gnomAD SAS | Confidence levels | 1-page summary | All input types

import streamlit as st
import json
import pandas as pd
import os
from datetime import date

# ── IMPORTS (with graceful fallback for missing optional deps) ────────────────
from parser import parse_vcf
from annotator import annotate_all, filter_rare_variants, rank_variants, apply_acmg_classification
from report_generator import generate_report, generate_word_report
from input_handler import (
    parse_plain_text, parse_pdf_report, parse_csv_excel,
    lookup_rsid, parse_clinical_info
)

try:
    from pdf_generator import generate_pdf_download
    PDF_AVAILABLE = True
except ImportError:
    PDF_AVAILABLE = False

# ── CONFIG ───────────────────────────────────────────────────────────────────
MAX_VARIANTS = 10

st.set_page_config(
    page_title="SpectralG",
    page_icon="🧬",
    layout="wide"
)

# ── TITLE ────────────────────────────────────────────────────────────────────
st.title("🧬 SpectralG")
st.markdown(
    "**AI-powered clinical variant interpretation** — from any input to clinical report in minutes.\n\n"
    "*Ensembl VEP · ACMG/AMP 2015 + 2023 · Claude AI · gnomAD SAS*"
)
st.caption("Research tool — requires review by a qualified clinical geneticist before clinical use")
st.divider()

# ── SIDEBAR ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("🔑 Claude API Key")
    user_api_key = st.text_input(
        "Enter your Claude API key",
        type="password",
        placeholder="sk-ant-api03-...",
        help="Get a key at console.anthropic.com — required for AI report generation"
    )
    if user_api_key:
        os.environ["ANTHROPIC_API_KEY"] = user_api_key
        st.success("✅ API key loaded")
    else:
        st.info("ℹ️ Enter API key to generate AI interpretation")

    st.divider()

    # ── REPORT TYPE ──────────────────────────────────────
    st.header("📋 Report Settings")
    report_type = st.selectbox(
        "Test / Report Type",
        ["Clinical Exome (WES)", "Whole Genome Sequencing (WGS)",
         "Disease Panel", "Carrier Screening",
         "Pharmacogenomics (PGx)", "CNV / Chromosomal",
         "Oncology Somatic Panel", "Prenatal",
         "Mitochondrial Panel"],
        index=0,
        help="Determines which sections appear in the output report"
    )

    report_language = st.selectbox(
        "Report Language",
        ["English", "Tamil (தமிழ்)", "Malayalam (മലയാളം)",
         "Kannada (ಕನ್ನಡ)", "Telugu (తెలుగు)", "Hindi (हिंदी)"],
        index=0
    )

    st.divider()

    # ── PATIENT INFO ─────────────────────────────────────
    st.header("👤 Patient Information")
    patient_id = st.text_input("Sample / Report ID", value="VC-2026-001")
    patient_name = st.text_input("Patient Name (optional)", value="",
                                 help="Leave blank for de-identified report")

    # Age
    age = st.text_input("Age", value="",
                        placeholder="e.g. 34 years / 6 months")

    # Sex — NEVER inferred from name
    sex = st.selectbox(
        "Biological Sex",
        ["Not provided", "Female", "Male", "Indeterminate / Not specified"],
        index=0,
        help="Enter only if explicitly documented in the clinical referral. "
             "Never assumed from patient name."
    )

    indication = st.text_area(
        "Clinical Indication",
        value="",
        placeholder="e.g. Bilateral renal cysts, hypertension, family history of ESRD",
        height=80
    )

    clinical_features = st.text_area(
        "Clinical Features / History",
        value="",
        placeholder="Presenting symptoms, relevant investigations, family history...",
        height=80
    )

    referring_clinician = st.text_input("Referring Clinician", value="")

    # Genotype-phenotype correlation — explicit
    gp_correlation = st.selectbox(
        "Genotype-Phenotype Correlation",
        ["Not assessed", "Present", "Absent", "Partial"],
        index=0,
        help="State explicitly whether the clinical features match the identified variant"
    )

    st.divider()

    st.header("⚙️ Advanced Settings")
    show_raw = st.checkbox("Show raw variant data", False)
    fetch_sas = st.checkbox("Fetch gnomAD SAS (slower)", True,
                            help="Fetches South Asian subpopulation frequency separately. "
                                 "Recommended for Indian patients.")

    st.divider()
    st.warning(
        "⚠️ Do NOT upload real patient data to public apps.\n\n"
        "For clinical use, deploy locally or on secure servers."
    )
    st.caption("SpectralG v2.0 | ACMG/AMP 2015 + 2023 | PP5 not applied")


# ── BUILD CLINICAL INFO DICT ─────────────────────────────────────────────────
clinical_info = parse_clinical_info({
    "patient_id":     patient_id,
    "patient_name":   patient_name or "[Not provided]",
    "age":            age or "Not provided",
    "sex":            sex,
    "indication":     indication,
    "clinical_features": clinical_features,
    "family_history": "",
    "report_type":    report_type,
    "referring_clinician": referring_clinician or "Not provided",
    "genotype_phenotype_correlation": gp_correlation,
    "gp_narrative": f"Genotype-phenotype correlation: {gp_correlation}. "
                    f"Clinical features: {clinical_features[:200] if clinical_features else 'not provided.'}",
})

# ── INPUT METHOD SELECTOR ────────────────────────────────────────────────────
st.subheader("📥 Input Method")

input_method = st.radio(
    "Choose how to provide variant data:",
    ["VCF File", "Plain Text / Gene Name", "PDF Lab Report", "CSV / Excel", "SNP / rsID"],
    horizontal=True,
    help="All methods feed into the same annotation and classification pipeline."
)

source_variants = None

# ── INPUT METHOD: VCF ────────────────────────────────────────────────────────
if input_method == "VCF File":
    col1, col2 = st.columns([3,1])
    with col1:
        uploaded_file = st.file_uploader("Upload VCF file (.vcf)", type=["vcf"])
    with col2:
        st.write("")
        st.write("")
        if st.button("📋 Use Demo Data"):
            st.session_state["demo"] = True

    if uploaded_file:
        source_variants = parse_vcf(uploaded_file)
    elif st.session_state.get("demo"):
        demo_vcf = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t43071077\t.\tA\tT\t50\tPASS\t.
17\t7674220\t.\tC\tT\t45\tPASS\t.
7\t117548628\t.\tA\tG\t60\tPASS\t.
"""
        source_variants = parse_vcf(demo_vcf.splitlines(keepends=True))

# ── INPUT METHOD: PLAIN TEXT ─────────────────────────────────────────────────
elif input_method == "Plain Text / Gene Name":
    st.info("💡 Enter any variant description: gene name, HGVS notation, SNP description, or mixed text")
    text_input = st.text_area(
        "Variant description",
        placeholder=(
            "Examples:\n"
            "• BRCA1 c.5266dup p.Gln1756fs Heterozygous Pathogenic female 38yo hereditary cancer\n"
            "• Gene: PKD1 | Variant: c.10255C>T (p.Arg3419Ter) | Zygosity: Heterozygous | LP\n"
            "• MAP1A c.5983C>T nonsense heterozygous ADHD 12yo"
        ),
        height=120
    )
    if st.button("🔍 Parse Variant", type="primary") and text_input.strip():
        with st.spinner("Parsing variant description..."):
            source_variants = parse_plain_text(text_input, user_api_key)
        if source_variants and not source_variants[0].get("error"):
            st.success(f"✅ Extracted {len(source_variants)} variant(s)")
        else:
            err = source_variants[0].get("error","Unknown error") if source_variants else "Parse failed"
            st.error(f"❌ {err}")
            source_variants = None

# ── INPUT METHOD: PDF ─────────────────────────────────────────────────────────
elif input_method == "PDF Lab Report":
    st.info("📄 Upload a lab report PDF — variant data will be extracted automatically using AI")
    pdf_file = st.file_uploader("Upload PDF report", type=["pdf"])
    if pdf_file:
        with st.spinner("Extracting variant data from PDF..."):
            source_variants = parse_pdf_report(pdf_file.read(), user_api_key)
        if source_variants and not source_variants[0].get("error"):
            st.success(f"✅ Extracted {len(source_variants)} variant(s) from PDF")
            st.info("📋 Review extracted variants below — edit clinical info in sidebar if needed")
        else:
            err = source_variants[0].get("error","Extraction failed") if source_variants else "Extraction failed"
            st.error(f"❌ {err}")
            source_variants = None

# ── INPUT METHOD: CSV / EXCEL ─────────────────────────────────────────────────
elif input_method == "CSV / Excel":
    st.info("📊 Upload a spreadsheet with columns: Gene, HGVS_c, HGVS_p, Zygosity, Classification")
    csv_file = st.file_uploader("Upload CSV or Excel file", type=["csv","xlsx","xls"])
    if csv_file:
        with st.spinner("Parsing spreadsheet..."):
            source_variants = parse_csv_excel(csv_file.read(), csv_file.name)
        if source_variants and not source_variants[0].get("error"):
            st.success(f"✅ Loaded {len(source_variants)} variant(s) from file")
        else:
            err = source_variants[0].get("error","Parse failed") if source_variants else "Parse failed"
            st.error(f"❌ {err}")
            source_variants = None

# ── INPUT METHOD: SNP / rsID ──────────────────────────────────────────────────
elif input_method == "SNP / rsID":
    rsid_input = st.text_input(
        "Enter rsID",
        placeholder="rs80357914",
        help="Enter a dbSNP rsID to fetch gene, position, and ClinVar data automatically"
    )
    if st.button("🔍 Look Up SNP", type="primary") and rsid_input.strip():
        with st.spinner(f"Fetching {rsid_input} from dbSNP..."):
            source_variants = lookup_rsid(rsid_input.strip())
        if source_variants and not source_variants[0].get("error"):
            st.success(f"✅ Found: {source_variants[0].get('gene','Unknown')} | rsID: {rsid_input}")
        else:
            err = source_variants[0].get("error","Not found") if source_variants else "Lookup failed"
            st.error(f"❌ {err}")
            source_variants = None


# ── PIPELINE ─────────────────────────────────────────────────────────────────
if source_variants:

    # Validate
    valid = [v for v in source_variants if not v.get("error") and isinstance(v, dict)]
    if not valid:
        st.error("❌ No valid variants to process")
        st.stop()

    if len(valid) > MAX_VARIANTS:
        st.warning(f"⚠️ Processing first {MAX_VARIANTS} variants only.")
        valid = valid[:MAX_VARIANTS]

    st.success(f"✅ {len(valid)} variant(s) ready for annotation")

    if show_raw:
        st.subheader("Raw Variant Data")
        st.json(valid)

    # ── ANNOTATION PIPELINE ───────────────────────────────
    try:
        progress = st.progress(0, text="Starting annotation pipeline...")

        # Step 1: VEP annotation (skip if manual entry)
        progress.progress(20, "Step 1/4: Annotating via Ensembl VEP...")
        needs_vep = [v for v in valid if v.get("source") not in ("plain_text","csv_excel","rsid_lookup")]
        manual = [v for v in valid if v.get("source") in ("plain_text","csv_excel","rsid_lookup")]

        if needs_vep:
            annotated_vep = annotate_all(needs_vep)
        else:
            annotated_vep = []

        annotated = annotated_vep + manual

        # Step 2: gnomAD SAS enrichment
        if fetch_sas:
            progress.progress(40, "Step 2/4: Fetching gnomAD South Asian frequencies...")
            from annotator import enrich_gnomad_sas
            annotated = enrich_gnomad_sas(annotated)

        # Step 3: Filter + Rank
        progress.progress(60, "Step 3/4: Filtering and ranking variants...")
        before = len(annotated)
        filtered = filter_rare_variants(annotated)
        ranked = rank_variants(filtered)
        removed = before - len(filtered)

        # Step 4: ACMG classification
        progress.progress(80, "Step 4/4: Applying ACMG criteria (PP5 not applied)...")
        final = apply_acmg_classification(ranked)

        progress.progress(100, "✅ Annotation complete")
        progress.empty()

    except Exception as e:
        st.error(f"❌ Annotation pipeline failed: {str(e)}")
        st.stop()

    if removed > 0:
        st.info(f"Filtered out {removed} common variant(s) (gnomAD AF > 1%)")
    st.success("✅ Annotation complete")

    # ── SUMMARY DASHBOARD ────────────────────────────────
    st.subheader("📊 Summary Dashboard")

    m1, m2, m3, m4, m5 = st.columns(5)
    total = len(annotated)
    high  = len([v for v in annotated if v.get("priority") == "HIGH"])
    path  = len([v for v in annotated if v.get("acmg") in ["Pathogenic","Likely Pathogenic"]])
    vus   = len([v for v in annotated if v.get("acmg") == "VUS"])
    benign= len([v for v in annotated if v.get("acmg") in ["Likely Benign","Benign"]])

    m1.metric("Total Variants", total)
    m2.metric("🔴 HIGH Priority", high)
    m3.metric("⚠️ Pathogenic / LP", path)
    m4.metric("🟡 VUS", vus)
    m5.metric("🟢 Benign / LB", benign)
    st.divider()

    # ACMG pie chart
    import plotly.express as px
    acmg_counts = {}
    for v in annotated:
        label = v.get("acmg","VUS")
        acmg_counts[label] = acmg_counts.get(label,0) + 1

    if acmg_counts:
        col_chart, col_info = st.columns([2,1])
        with col_chart:
            fig = px.pie(
                names=list(acmg_counts.keys()),
                values=list(acmg_counts.values()),
                title="ACMG Classification Distribution",
                color=list(acmg_counts.keys()),
                color_discrete_map={
                    "Pathogenic":"#c62828","Likely Pathogenic":"#e53935",
                    "VUS":"#f9a825","Likely Benign":"#43a047","Benign":"#1b5e20"
                }
            )
            fig.update_layout(height=280, margin=dict(t=40,b=0,l=0,r=0))
            st.plotly_chart(fig, use_container_width=True)
        with col_info:
            st.markdown("**Note on PP5:**")
            st.info(
                "PP5 (external assertion) has been **removed** from this version "
                "per ACMG 2023 guidance (Biesecker & Harrison). "
                "Classifications are based on primary evidence only."
            )

    st.divider()

    # ── ANNOTATED VARIANT TABLE ───────────────────────────
    st.subheader("Annotated Variants")

    def _fmt_af(v):
        af = v.get("gnomad_af",{})
        if isinstance(af, dict):
            sas = af.get("south_asian")
            glb = af.get("global")
            if sas is not None:
                return f"SAS:{sas:.5f}"
            if glb is not None:
                return f"Global:{glb:.5f}"
        try:
            return f"{float(af):.5f}"
        except:
            return "N/A"

    table_rows = []
    for v in final:
        ann = v.get("annotation",{})
        table_rows.append({
            "Gene":       v.get("gene","Unknown"),
            "HGVS c.":    v.get("hgvsc","") or ann.get("hgvsc",""),
            "HGVS p.":    v.get("hgvsp","") or ann.get("hgvsp",""),
            "Consequence":v.get("consequence","?"),
            "Impact":     ann.get("impact","?"),
            "gnomAD SAS": _fmt_af(v),
            "ACMG":       v.get("acmg","VUS"),
            "Confidence": v.get("confidence_level","N/A"),
            "Priority":   v.get("priority","NA"),
            "PP5 Applied":"No (ACMG 2023)",
            "Evidence":   ", ".join(v.get("acmg_evidence",[])),
        })

    df = pd.DataFrame(table_rows)

    def _colour_row(row):
        priority = row.get("Priority","")
        if priority == "HIGH":
            return ["background-color:#ffebee;color:#c62828"] * len(row)
        elif priority == "MEDIUM":
            return ["background-color:#fff3e0;color:#e65100"] * len(row)
        elif priority == "LOW":
            return ["background-color:#e8f5e9;color:#2e7d32"] * len(row)
        return [""] * len(row)

    styled_df = df.style.apply(_colour_row, axis=1)
    st.dataframe(styled_df, use_container_width=True)
    st.caption("🟢 Low priority | 🟡 Medium priority | 🔴 High priority | gnomAD SAS = South Asian subpopulation")
    st.divider()

    # ── ACMG DETAIL EXPANDER ──────────────────────────────
    if final:
        with st.expander("🔬 View Full ACMG Criteria Evidence (Top Variant)", expanded=False):
            top = final[0]
            criteria = top.get("acmg_criteria_table",[])
            if criteria:
                crit_rows = []
                for c in criteria:
                    crit_rows.append({
                        "Criterion": c["code"],
                        "Weight":    c.get("weight",""),
                        "Applied":   "✓ YES" if c.get("applied") else "No",
                        "Evidence":  c.get("evidence","")[:200],
                    })
                st.dataframe(pd.DataFrame(crit_rows), use_container_width=True)
                st.caption("PP5 not applied per ACMG 2023 guidance.")
            else:
                st.info("ACMG criteria table not available for this variant.")

        with st.expander("📋 Evidence Panel (3billion-style)", expanded=False):
            panel = final[0].get("evidence_panel",{})
            if panel:
                for label, content in panel.items():
                    st.markdown(f"**{label}**")
                    st.markdown(content)
                    st.divider()
            else:
                st.info("Evidence panel not available.")

    st.divider()

    # ── REPORT GENERATION ─────────────────────────────────
    st.subheader("📋 Generate Clinical Report")

    col_gen1, col_gen2 = st.columns(2)
    with col_gen1:
        generate_btn = st.button("🤖 Generate AI Report", type="primary",
                                  disabled=not (user_api_key or os.environ.get("ANTHROPIC_API_KEY")))
    with col_gen2:
        if not (user_api_key or os.environ.get("ANTHROPIC_API_KEY")):
            st.warning("API key required for AI report generation")

    if generate_btn:
        try:
            with st.spinner("Generating clinical report..."):
                report = generate_report(final, patient_id, clinical_info, report_language)
        except Exception as e:
            st.error(f"❌ Report generation failed: {str(e)}")
            st.stop()

        st.subheader("Clinical Interpretation")
        st.markdown(report)
        st.divider()

        # ── DOWNLOADS ─────────────────────────────────────
        st.subheader("⬇️ Download Report")
        dl1, dl2, dl3, dl4 = st.columns(4)

        with dl1:
            st.download_button(
                "📄 Download TXT",
                data=report,
                file_name=f"{patient_id}_report.txt",
                mime="text/plain"
            )

        with dl2:
            st.download_button(
                "📊 Download JSON",
                data=json.dumps(final, indent=2, default=str),
                file_name=f"{patient_id}_variants.json",
                mime="application/json"
            )

        with dl3:
            try:
                word_file = generate_word_report(final, report, patient_id, report_language)
                with open(word_file, "rb") as f:
                    st.download_button(
                        "📝 Download Word",
                        data=f.read(),
                        file_name=word_file,
                        mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
                    )
            except Exception as e:
                st.error(f"Word generation failed: {e}")

        with dl4:
            if PDF_AVAILABLE:
                try:
                    pdf_bytes, pdf_filename = generate_pdf_download(
                        final, clinical_info, report, patient_id
                    )
                    st.download_button(
                        "📕 Download Clinical PDF",
                        data=pdf_bytes,
                        file_name=pdf_filename,
                        mime="application/pdf",
                        type="primary"
                    )
                except Exception as e:
                    st.error(f"PDF generation failed: {e}")
            else:
                st.warning("Install reportlab for PDF output: `pip install reportlab`")

    # Quick PDF without AI (for fast testing)
    elif PDF_AVAILABLE:
        if st.button("📕 Generate PDF Report (no AI text)", type="secondary"):
            try:
                with st.spinner("Generating clinical PDF..."):
                    pdf_bytes, pdf_filename = generate_pdf_download(
                        final, clinical_info, "", patient_id
                    )
                st.download_button(
                    "⬇️ Download Clinical PDF",
                    data=pdf_bytes,
                    file_name=pdf_filename,
                    mime="application/pdf",
                )
                st.success("✅ PDF generated — includes ACMG table, evidence panel, and summary box")
            except Exception as e:
                st.error(f"PDF generation failed: {str(e)}")

# ── FOOTER ────────────────────────────────────────────────────────────────────
st.divider()
st.caption(
    "SpectralG v2.0 | Built using Streamlit · Ensembl VEP · Anthropic Claude | India\n\n"
    "ACMG/AMP 2015 + 2023 updates + ACMG/ACGS-2024v1.2 | PP5 not applied (ACMG 2023) | "
    "gnomAD SAS subpopulation used for Indian patients"
)