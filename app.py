# app.py
# Streamlit UI for AI Variant Report Generator

import streamlit as st
import json
import pandas as pd
import os
from parser import parse_vcf
from annotator import annotate_all, filter_rare_variants, rank_variants, apply_acmg_classification
from report_generator import generate_report

# ── CONFIG ──────────────────────────────────────────────
MAX_VARIANTS = 10

st.set_page_config(
    page_title="SpectralG",
    #page_icon="🧬",
    layout="centered"
)

#st.title("🧬 SpectralG")
st.markdown("AI-powered clinical variant interpretation — from VCF to clinical report in minutes.")
st.caption("Ensembl VEP + ACMG + Claude AI | Research tool — requires clinical geneticist review")
st.divider()

# ── SIDEBAR ─────────────────────────────────────────────
with st.sidebar:
    st.header("Patient Info")
    patient_id = st.text_input("Sample ID", value="SAMPLE_001")

    st.divider()

    st.header("🔑 Claude API Key")
    user_api_key = st.text_input(
        "Enter your Claude API key",
        type="password",
        placeholder="sk-ant-api03-...",
        help="Get a free key at console.anthropic.com"
    )
    if user_api_key:
        os.environ["ANTHROPIC_API_KEY"] = user_api_key
        st.success("✅ API key loaded")
    else:
        st.info("ℹ️ Enter your Claude API key to generate AI reports")

    st.divider()

    st.header("Settings")
    show_raw = st.checkbox("Show raw variant data", False)

    st.divider()
    st.header("🌐 Report Language")
    report_language = st.selectbox(
        "Generate report in:",
        ["English", "Tamil (தமிழ்)","Malayalam (മലയാളം)","Kannada (ಕನ್ನಡ)", "Telugu (తెలుగు)", "Hindi (हिंदी)"],
        index=0
    )

    st.divider()

    st.warning(
        "⚠️ Do NOT upload real patient data to public apps.\n\n"
        "For clinical use, deploy locally or on secure servers."
    )

# ── FILE INPUT ──────────────────────────────────────────
uploaded_file = st.file_uploader(
    "Upload VCF file (.vcf)",
    type=["vcf"]
)

if st.button("📋 Use Demo Data"):
    st.session_state["demo"] = True

# ── DATA SOURCE ─────────────────────────────────────────
source = None

if uploaded_file:
    source = uploaded_file

elif st.session_state.get("demo"):
    demo_vcf = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t43071077\t.\tA\tT\t50\tPASS\t.
17\t7674220\t.\tC\tT\t45\tPASS\t.
7\t117548628\t.\tA\tG\t60\tPASS\t.
"""
    source = demo_vcf.splitlines(keepends=True)

# ── PIPELINE ────────────────────────────────────────────
if source:

    # Step 1: Parse
    try:
        with st.spinner("Step 1/3: Parsing VCF..."):
            variants = parse_vcf(source)

        if not variants:
            st.error("❌ No variants found. Invalid VCF.")
            st.stop()

    except Exception as e:
        st.error(f"❌ Parsing failed: {str(e)}")
        st.stop()

    # Limit variants
    if len(variants) > MAX_VARIANTS:
        st.warning(f"Only first {MAX_VARIANTS} variants will be processed.")
        variants = variants[:MAX_VARIANTS]

    st.success(f"✅ {len(variants)} variants loaded")

    if show_raw:
        st.subheader("Raw Variants")
        st.json(variants)

    # Step 2: Annotation Pipeline
    try:
        with st.spinner("Step 2/3: Annotating + Filtering + Ranking..."):

            # ✅ STEP 2.1: Annotate
            annotated = annotate_all(variants)

            # ✅ STEP 2.2: Filter (gnomAD)
            before = len(annotated)
            filtered = filter_rare_variants(annotated)
            st.info(f"Filtered {before - len(filtered)} common variants (AF > 1%)")

            # ✅ STEP 2.3: Rank
            ranked = rank_variants(filtered)

            # ✅ STEP 2.4: ACMG
            final = apply_acmg_classification(ranked)

    except Exception as e:
        st.error(f"❌ Annotation failed: {str(e)}")
        st.stop()

    st.success("✅ Annotation complete")

    # ── SUMMARY DASHBOARD ───────────────────────────────
    st.subheader("📊 Summary")

    col1, col2, col3, col4 = st.columns(4)

    total = len(annotated)
    high = len([v for v in annotated if v.get("priority") == "HIGH"])
    medium = len([v for v in annotated if v.get("priority") == "MEDIUM"])
    pathogenic = len([v for v in annotated if v.get("acmg") in [
        "Pathogenic", "Likely pathogenic", "Likely Pathogenic"
    ]])

    col1.metric("Total Variants", total)
    col2.metric("🔴 HIGH Priority", high)
    col3.metric("🟡 MEDIUM Priority", medium)
    col4.metric("⚠️ Pathogenic/LP", pathogenic)

    st.divider()

    # ── ACMG PIE CHART ──────────────────────────────────
    acmg_counts = {}
    for v in annotated:
        label = v.get("acmg", "VUS")
        acmg_counts[label] = acmg_counts.get(label, 0) + 1

    if acmg_counts:
        import plotly.express as px
        fig = px.pie(
            names=list(acmg_counts.keys()),
            values=list(acmg_counts.values()),
            title="ACMG Classification Distribution",
            color=list(acmg_counts.keys()),
            color_discrete_map={
                "Pathogenic": "#c62828",
                "Likely pathogenic": "#e53935",
                "Likely Pathogenic": "#e53935",
                "VUS": "#f9a825",
                "Likely Benign": "#43a047",
                "Benign": "#1b5e20"
            }
        )
        fig.update_layout(height=300, margin=dict(t=40, b=0, l=0, r=0))
        st.plotly_chart(fig, use_container_width=True)

    st.divider()

    # 🔥 TOP PRIORITY SECTION
    # ── GENE FUNCTION CARDS ─────────────────────────────
    st.subheader("🧬 Identified Genes")

    gene_info = {
        "TP53": "Tumour suppressor gene. Mutations found in over 50% of human cancers. Associated with Li-Fraumeni syndrome.",
        "BRCA1": "DNA repair gene. Pathogenic variants significantly increase risk of breast and ovarian cancer.",
        "BRCA2": "DNA repair gene. Associated with breast, ovarian, and pancreatic cancer risk.",
        "CFTR": "Chloride channel gene. Mutations cause cystic fibrosis — a serious lung and digestive disease.",
        "HBB": "Haemoglobin beta gene. Mutations cause sickle cell disease and beta-thalassaemia.",
        "EGFR": "Epidermal growth factor receptor. Frequently mutated in lung cancer. Targetable by specific drugs.",
        "KRAS": "Cell signalling gene. One of the most commonly mutated genes in pancreatic and colorectal cancer.",
        "MLH1": "DNA mismatch repair gene. Mutations cause Lynch syndrome — hereditary colorectal cancer.",
        "PTEN": "Tumour suppressor. Mutations cause Cowden syndrome — increased risk of breast and thyroid cancer.",
        "APC":  "Colorectal cancer suppressor. Mutations cause familial adenomatous polyposis.",
    }

    unique_genes = list(set([
        v.get("gene", "Unknown") for v in annotated
        if v.get("gene") not in ["Unknown", "Intergenic", None]
    ]))

    if unique_genes:
        gene_cols = st.columns(len(unique_genes))
        for i, gene in enumerate(unique_genes):
            with gene_cols[i]:
                info = gene_info.get(gene, f"{gene} — clinical significance requires expert review.")
                priority_variants = [v for v in annotated if v.get("gene") == gene]
                top_priority = priority_variants[0].get("priority", "LOW") if priority_variants else "LOW"
                colour = "🔴" if top_priority == "HIGH" else "🟡" if top_priority == "MEDIUM" else "🟢"
                st.info(f"{colour} **{gene}**\n\n{info}")
                
                
    st.subheader("🔥 Top Priority Variants")
    top = [v for v in final if v.get("priority") == "HIGH"]

    if top:
        st.success(f"{len(top)} high-priority variants found")
    else:
        st.warning("No high-priority variants detected")

    # ── TABLE ───────────────────────────────────────────
    st.subheader("Annotated Variants")

    table = []
    for v in final:
        ann = v.get("annotation", {})

        table.append({
            "Gene": v.get("gene", "Unknown"),
            "Position": f"{v['chrom']}:{v['pos']}",
            "Change": f"{v['ref']}>{v['alt']}",
            "Impact": ann.get("impact", "?"),
            "Consequence": ann.get("consequence", "?"),
            "ClinVar": v.get("clinvar", "Unknown"),
            "SIFT": str(ann.get("sift", "N/A")),
            "PolyPhen": str(ann.get("polyphen", "N/A")),
            "gnomAD_AF": v.get("gnomad_af", "N/A"),
            "Score": v.get("score", 0),
            "Priority": v.get("priority", "NA"),
            "ACMG": v.get("acmg", "NA"),
            "Evidence": ", ".join(v.get("acmg_evidence", [])),
        })

    import pandas as pd

    df = pd.DataFrame(table)

    def colour_priority(row):
        if row["Priority"] == "HIGH":
            return ["background-color: #ffebee; color: #c62828"] * len(row)
        elif row["Priority"] == "MEDIUM":
            return ["background-color: #fff3e0; color: #e65100"] * len(row)
        elif row["Priority"] == "LOW":
            return ["background-color: #e8f5e9; color: #2e7d32"] * len(row)
        else:
            return [""] * len(row)

    styled_df = df.style.apply(colour_priority, axis=1)
    st.dataframe(styled_df, use_container_width=True)
    st.divider()

    # Step 3: Report
    if st.button("🤖 Generate AI Report", type="primary"):

        try:
            with st.spinner("Step 3/3: Generating report..."):
                report = generate_report(final, patient_id)

        except Exception as e:
            st.error(f"❌ Report generation failed: {str(e)}")
            st.stop()

        st.subheader("📋 Clinical Report")
        st.markdown(report)
        st.divider()

        # Downloads
        col1, col2, col3 = st.columns(3)

        with col1:
            st.download_button(
                "⬇️ Download TXT",
                data=report,
                file_name=f"{patient_id}_report.txt",
                mime="text/plain"
            )

        with col2:
            st.download_button(
                "⬇️ Download JSON",
                data=json.dumps(annotated, indent=2),
                file_name=f"{patient_id}_variants.json",
                mime="application/json"
            )

        with col3:
            from report_generator import generate_word_report
            word_file = generate_word_report (annotated,report, patient_id,report_language)
            with open(word_file, "rb") as f:
                st.download_button(
                    "⬇️ Download Word",
                    data=f.read(),
                    file_name=word_file,
                    mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
                )

# ── FOOTER ─────────────────────────────────────────────
st.divider()
st.caption("SpectralG | Built using Streamlit · Ensembl VEP · anthropic | India")