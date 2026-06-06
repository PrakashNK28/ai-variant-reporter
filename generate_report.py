# generate_report.py
# SpectralG Commercial Report Generator
# Usage: python3 generate_report.py --json variants.json --patient "Patient Name" --id "VC-2026-001"
#
# This is the paid service workflow — uses claude-sonnet-4-6 for higher quality output.
# temperature=0.1 enforced for all clinical reports (anti-hallucination).

import json
import os
import argparse
from pathlib import Path
from dotenv import load_dotenv
import anthropic

load_dotenv(dotenv_path=Path.home() / ".env", override=True)


def load_prompt():
    """Load SpectralG system prompt from spectralg_prompt.txt if it exists."""
    prompt_file = Path(__file__).parent / "spectralg_prompt.txt"
    if prompt_file.exists():
        return prompt_file.read_text()
    # Default system prompt if file not found
    return (
        "You are an expert clinical geneticist generating professional variant "
        "interpretation reports for SpectralG. Apply ACMG/AMP 2015 + 2023 criteria. "
        "PP5 is NEVER applied. State confidence levels for every classification. "
        "Never say 'diagnosis confirmed' — use 'findings suggest' or 'consistent with'. "
        "Never infer sex from patient name. Never say treatment is mandatory. "
        "Every claim must be traceable to the variant data provided. "
        "Never invent statistics, patient numbers, or literature citations."
    )


def format_gnomad(af):
    """Format gnomAD frequency for report text."""
    if isinstance(af, dict):
        sas = af.get("south_asian")
        glb = af.get("global")
        if sas is not None:
            return f"SAS: {sas:.6f}"
        if glb is not None:
            return f"Global: {glb:.6f}"
        return "Not in gnomAD"
    try:
        return f"{float(af):.6f}"
    except Exception:
        return "Not in gnomAD"


def generate_commercial_report(variants_json_path, patient_id,
                                patient_name, clinical_info=""):
    """
    Generate a commercial-grade SpectralG report from a JSON variants file.

    Args:
        variants_json_path: path to annotated variants JSON file
        patient_id:         report/sample ID string
        patient_name:       patient name or [De-identified]
        clinical_info:      clinical indication and context string

    Returns:
        report_text: complete report as string
    """
    # Load variants
    with open(variants_json_path) as f:
        variants = json.load(f)

    if not variants:
        raise ValueError("No variants found in JSON file.")

    system_prompt = load_prompt()

    # Build variant summary blocks
    variant_blocks = []
    for i, v in enumerate(variants):
        ann     = v.get("annotation", {})
        ct      = v.get("acmg_criteria_table", [])
        applied = [c["code"] for c in ct
                   if c.get("applied") and c["code"] not in {"PP5", "BP6"}]

        block = f"""
VARIANT {i+1}:
Gene:                  {v.get('gene', 'Unknown')}
Position:              chr{v.get('chrom', '?')}:{v.get('pos', '?')}
Change:                {v.get('ref', '?')} > {v.get('alt', '?')}
HGVS c.:               {v.get('hgvsc', 'Not available')}
HGVS p.:               {v.get('hgvsp', 'Not available')}
Zygosity:              {v.get('zygosity', 'Unknown')}
Consequence:           {ann.get('consequence', 'unknown')}
VEP Impact:            {ann.get('impact', 'UNKNOWN')}
SIFT:                  {ann.get('sift', 'N/A')}
PolyPhen-2:            {ann.get('polyphen', 'N/A')}
gnomAD AF:             {format_gnomad(v.get('gnomad_af'))}
ClinVar:               {v.get('clinvar', 'Unknown')}
ACMG Classification:   {v.get('acmg', 'VUS')}
Confidence Level:      {v.get('confidence_level', 'Limited')}
Applied Criteria:      {', '.join(applied) if applied else 'None'}
Priority:              {v.get('priority', 'LOW')}
Score:                 {v.get('score', 0)}/10
"""
        variant_blocks.append(block)

    # Summary counts
    high_priority = [v for v in variants if v.get("priority") == "HIGH"]
    pathogenic    = [v for v in variants if v.get("acmg") in
                     ("Pathogenic", "Likely Pathogenic")]

    user_message = f"""Generate a complete SpectralG clinical variant interpretation report.

Patient ID:          {patient_id}
Patient Name:        {patient_name}
Clinical Information:{clinical_info if clinical_info else 'Not provided'}
Total Variants:      {len(variants)}
High Priority:       {len(high_priority)}
Pathogenic / LP:     {len(pathogenic)}

VARIANT DATA (use ONLY this — do not invent any additional information):
{''.join(variant_blocks)}

Generate the complete report with these sections:

1. EXECUTIVE SUMMARY (2-3 sentences, key finding, genotype-phenotype correlation)
2. VARIANT INTERPRETATIONS (one block per variant with all fields above)
3. CLINICAL RECOMMENDATIONS (standard-of-care per guideline — never say mandatory)
4. GENOTYPE-PHENOTYPE CORRELATION (Present / Absent / Partial / Not assessed)
5. METHODS NOTE (Ensembl VEP, ACMG/AMP 2015 + 2023, gnomAD SAS, PP5 not applied)
6. DISCLAIMER (research tool, requires clinical geneticist review)

Rules:
- PP5 is NOT applied — do not mention as evidence
- State confidence level for every classification
- Never say "diagnosis confirmed" — use "findings suggest"
- Never infer sex from patient name
- Never invent statistics or literature citations
- If a field says Not available — write Not available, never guess
"""

    # Call Claude Sonnet at temperature=0.1 (factual mode)
    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        raise EnvironmentError(
            "ANTHROPIC_API_KEY not found. "
            "Set it in ~/.env or as environment variable."
        )

    client   = anthropic.Anthropic(api_key=api_key)
    response = client.messages.create(
        model       = "claude-sonnet-4-6",  # Sonnet for commercial quality
        max_tokens  = 2000,
        temperature = 0.1,                   # Factual mode — prevents hallucination
        system      = system_prompt,
        messages    = [{"role": "user", "content": user_message}]
    )

    report_text = response.content[0].text

    # Save report to file
    output_file = f"{patient_id}_SpectralG_commercial_report.txt"
    with open(output_file, "w") as f:
        f.write(report_text)

    print(f"✅ Commercial report generated: {output_file}")
    print(f"   Word count:  {len(report_text.split())}")
    print(f"   Variants:    {len(variants)}")
    print(f"   High priority: {len(high_priority)}")
    print(f"   Pathogenic/LP: {len(pathogenic)}")

    return report_text


# ── CLI ENTRYPOINT ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SpectralG Commercial Report Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 generate_report.py "
            "--json variants.json --patient 'Patient A' --id VC-2026-001\n"
            "  python3 generate_report.py "
            "--json output.json --patient '[De-identified]' --id LAB-001 "
            "--clinical 'Hereditary breast cancer risk assessment'"
        )
    )
    parser.add_argument("--json",     required=True,
                        help="Path to annotated variants JSON file")
    parser.add_argument("--patient",  required=True,
                        help="Patient name or [De-identified]")
    parser.add_argument("--id",       required=True,
                        help="Report / Sample ID")
    parser.add_argument("--clinical", default="",
                        help="Clinical indication and context (optional)")
    args = parser.parse_args()

    generate_commercial_report(
        args.json, args.id, args.patient, args.clinical
    )