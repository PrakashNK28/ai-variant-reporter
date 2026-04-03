# app.py
# Streamlit UI for AI Variant Report Generator

import streamlit as st
import json
from parser import parse_vcf
from annotator import annotate_all, filter_rare_variants, rank_variants, apply_acmg_classification
from report_generator import generate_report

# ── CONFIG ──────────────────────────────────────────────
MAX_VARIANTS = 10

st.set_page_config(
    page_title="AI Variant Report Generator",
    page_icon="🧬",
    layout="centered"
)

st.title("🧬 AI Variant Report Generator")
st.markdown("Upload a VCF file to generate a clinical genetics report.")
st.caption("VEP + ClinVar + AI | Research tool only")
st.divider()

# ── SIDEBAR ─────────────────────────────────────────────
with st.sidebar:
    st.header("Patient Info")
    patient_id = st.text_input("Sample ID", value="SAMPLE_001")

    st.divider()

    st.header("Settings")
    show_raw = st.checkbox("Show raw variant data", False)

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

    # 🔥 TOP PRIORITY
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

    st.dataframe(table, use_container_width=True)
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
        col1, col2 = st.columns(2)

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
                data=json.dumps(final, indent=2),
                file_name=f"{patient_id}_variants.json",
                mime="application/json"
            )

# ── FOOTER ─────────────────────────────────────────────
st.divider()
st.caption("Built with Streamlit | VEP | AI")