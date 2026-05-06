# report_generator.py
# AI-powered clinical report generator with robust fallback

import os
import json
import anthropic
from dotenv import load_dotenv
from pathlib import Path

# ── LOAD ENV ────────────────────────────────────────────
load_dotenv(dotenv_path=Path.home() / ".env", override=True)
load_dotenv()


def generate_report(variants, patient_id="SAMPLE_001", clinical_info=None, language="English"):
    """
    Generate a clinical variant interpretation report using AI (Claude).
    Falls back to rule-based report if API fails.
    """

    # ── Prepare structured variant summary ───────────────
    variant_summary = []

    for v in variants:
        ann = v.get("annotation", {})

        variant_summary.append({
            "gene": v.get("gene", "Unknown"),
            "position": f"chr{v['chrom']}:{v['pos']}",
            "change": f"{v['ref']} > {v['alt']}",
            "consequence": ann.get("consequence", "unknown"),
            "impact": ann.get("impact", "UNKNOWN"),
            "clinvar": v.get("clinvar", "Unknown"),
            "gnomad_af": v.get("gnomad_af", "N/A"),
            "sift": ann.get("sift", "N/A"),
            "polyphen": ann.get("polyphen", "N/A"),
            "protein_change": ann.get("hgvsp", ""),
            "acmg": v.get("acmg", "VUS"),
            "acmg_evidence": v.get("acmg_evidence", []),
            "priority": v.get("priority", "LOW"),
            "score": v.get("score", 0)
        })

    # ── PROMPT ───────────────────────────────────────────
    # ── PROMPT ───────────────────────────────────────────
    # Extract clinical context if provided
    if clinical_info is None:
        clinical_info = {}

    sex = clinical_info.get("sex", "Sex: not provided in referral")
    indication = clinical_info.get("indication", "Not provided")
    report_type = clinical_info.get("report_type", "Clinical WES")
    gp = clinical_info.get("genotype_phenotype_correlation", "Not assessed")

    prompt = f"""
You are an expert clinical geneticist using SpectralG, an AI-powered variant interpretation tool.

Generate a clinical variant interpretation report.

Patient ID: {patient_id}
Report Type: {report_type}
Sex: {sex}
Clinical Indication: {indication}
Genotype-Phenotype Correlation: {gp}
Variants Analysed: {len(variants)}

Variant Data:
{json.dumps(variant_summary, indent=2)}

Instructions:

1. EXECUTIVE SUMMARY
- State key finding clearly
- Highlight high-risk variants
- State genotype-phenotype correlation explicitly as present/absent/partial/not assessed

2. VARIANT INTERPRETATION
For each variant:
- Gene name and function
- ACMG classification with evidence
- Clinical relevance in plain language

3. CLINICAL RECOMMENDATIONS
- Next steps for the ordering clinician
- Do not say treatment is mandatory — say "standard-of-care per guideline, to be implemented by treating team"

4. DISCLAIMER
- AI-assisted, requires clinical geneticist review
- PP5 not applied per ACMG 2023

Keep report under 600 words.
Use professional medical language.
Generate the entire report in {language}.
If language is Tamil, Malayalam, Kannada, Telugu or Hindi — use that language throughout
but keep gene names, variant positions, and ACMG terms in English.
"""

    # ── AI GENERATION ────────────────────────────────────
    ai_text = None

    try:
        api_key = os.getenv("ANTHROPIC_API_KEY")

        if api_key:
            client = anthropic.Anthropic(api_key=api_key)

            response = client.messages.create(
                model="claude-haiku-4-5-20251001",
                max_tokens=800,
                messages=[{"role": "user", "content": prompt}]
            )

            ai_text = response.content[0].text

    except Exception as e:
        print(f"AI error: {e}")

    # ── FALLBACK LOGIC ───────────────────────────────────
    if not ai_text:
        ai_text = generate_rule_based_report(variants, patient_id)

    return ai_text


# ── RULE-BASED FALLBACK (SMART) ─────────────────────────
def generate_rule_based_report(variants, patient_id):
    """
    Smarter fallback using ACMG + impact
    """

    report = f"# Clinical Variant Report — {patient_id}\n\n"

    report += "## Executive Summary\n"
    high = [v for v in variants if v.get("priority") == "HIGH"]

    if high:
        report += f"{len(high)} high-priority variants detected.\n\n"
    else:
        report += "No high-risk variants detected.\n\n"

    report += "## Variant Interpretation\n"

    for v in variants:
        ann = v.get("annotation", {})

        gene = v.get("gene", "Unknown")
        acmg = v.get("acmg", "VUS")
        impact = ann.get("impact", "UNKNOWN")

        # 🧠 Interpretation logic
        if acmg == "Pathogenic":
            interpretation = "Likely disease-causing variant."
        elif acmg == "Likely pathogenic":
            interpretation = "Possibly disease-associated variant."
        elif acmg == "VUS":
            interpretation = "Uncertain clinical significance."
        else:
            interpretation = "Likely benign variant."

        report += f"""
Gene: {gene}
Position: chr{v['chrom']}:{v['pos']}
Change: {v['ref']} > {v['alt']}
Impact: {impact}
ACMG: {acmg}

Interpretation:
{interpretation}

-----------------------------------
"""

    report += "\n## Clinical Recommendations\n"
    report += "Consider clinical correlation, family history analysis, and functional validation if needed.\n"

    report += "\n## Disclaimer\n"
    report += "This is an AI-assisted research report and must be reviewed by a clinical geneticist.\n"

    return report


def generate_word_report(variants, report_text, patient_id, language="English"):
    from docx import Document
    from docx.shared import Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from datetime import date

    doc = Document()

    # Title
    title = doc.add_heading("SpectralG — Clinical Variant Report", 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # Patient info
    doc.add_paragraph(f"Sample ID: {patient_id}")
    doc.add_paragraph(f"Date: {date.today().strftime('%d %B %Y')}")
    doc.add_paragraph(f"Tool: SpectralG | Ensembl VEP + ACMG + Claude AI")
    doc.add_paragraph(f"Variants Analysed: {len(variants)}")
    doc.add_paragraph("")

    # Variant table
    doc.add_heading("Variant Summary", level=1)
    table = doc.add_table(rows=1, cols=6)
    table.style = "Table Grid"

    headers = ["Gene", "Position", "Change", "Impact", "ACMG", "Priority"]
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        for para in cell.paragraphs:
            for run in para.runs:
                run.bold = True

    for v in variants:
        ann = v.get("annotation", {})
        row = table.add_row().cells
        row[0].text = str(v.get("gene", "Unknown"))
        row[1].text = f"chr{v['chrom']}:{v['pos']}"
        row[2].text = f"{v['ref']}>{v['alt']}"
        row[3].text = str(ann.get("impact", "?"))
        row[4].text = str(v.get("acmg", "VUS"))
        row[5].text = str(v.get("priority", "LOW"))

    doc.add_paragraph("")

    # AI Report text
    doc.add_heading("Clinical Interpretation", level=1)
    for line in report_text.split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("## "):
            doc.add_heading(line.replace("## ", ""), level=2)
        elif line.startswith("# "):
            doc.add_heading(line.replace("# ", ""), level=1)
        else:
            doc.add_paragraph(line)

    # Disclaimer
    doc.add_heading("Disclaimer", level=1)
    p = doc.add_paragraph(
        "SpectralG is an AI-assisted research tool. "
        "All reports must be reviewed by a qualified clinical geneticist "
        "before clinical use. Do not use for diagnosis without professional review."
    )
    p.runs[0].italic = True

    filename = f"{patient_id}_SpectralG_report.docx"
    doc.save(filename)
    return filename