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


def generate_report(variants, patient_id="SAMPLE_001"):
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
    prompt = f"""
You are an expert clinical geneticist.

Analyze the following variant data and generate a clinical report.

Patient ID: {patient_id}
Variants analyzed: {len(variants)}

Variant Data:
{json.dumps(variant_summary, indent=2)}

Instructions:

1. EXECUTIVE SUMMARY
- Highlight high-risk variants

2. VARIANT INTERPRETATION
- Gene
- ACMG classification
- Clinical relevance

3. CLINICAL RECOMMENDATIONS

4. DISCLAIMER

Keep report under 400 words.
"""

    # ── AI GENERATION ────────────────────────────────────
    ai_text = None

    try:
        api_key = os.getenv("ANTHROPIC_API_KEY")

        if api_key:
            client = anthropic.Anthropic(api_key=api_key)

            response = client.messages.create(
                model="claude-3-haiku-20240307",
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