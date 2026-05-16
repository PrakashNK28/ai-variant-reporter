# generate_report.py
# SpectralG Commercial Report Generator
# Usage: python3 generate_report.py --json variants.json --patient "Patient Name" --id "VC-2026-001"

import json
import sys
import os
import argparse
from pathlib import Path
from dotenv import load_dotenv
import anthropic

load_dotenv(dotenv_path=Path.home() / ".env", override=True)


def load_prompt():
    prompt_file = Path(__file__).parent / "spectralg_prompt.txt"
    if prompt_file.exists():
        return prompt_file.read_text()
    return ""


def generate_commercial_report(variants_json_path, patient_id,
                                patient_name, clinical_info=""):
    """
    Generate a commercial-grade SpectralG report from a JSON variants file.
    This is the paid service workflow.
    """
    # Load variants
    with open(variants_json_path) as f:
        variants = json.load(f)

    system_prompt = load_prompt()

    # Build variant summary for Claude
    variant_blocks = []
    for i, v in enumerate(variants):
        ann = v.get("annotation", {})
        ct = v.get("acmg_criteria_table", [])
        applied = [c["code"] for c in ct
                   if c.get("applied") and c["code"] not in {"PP5", "BP6"}]

        af = v.get("gnomad_af", {})
        if isinstance(af, dict):
            sas = af.get("south_asian")
            glb = af.get("global")
            af_str = f"SAS: {sas:.6f}" if sas else f"Global: {glb:.6f}" if glb else "Not in gnomAD"
        else:
            af_str = "Not in gnomAD"

        block = f"""
VARIANT {i+1}:
Gene: {v.get('gene', 'Unknown')}
Position: chr{v.get('chrom')}:{v.get('pos')}
Change: {v.get('ref')} > {v.get('alt')}
HGVS c.: {v.get('hgvsc', 'Not available')}
HGVS p.: {v.get('hgvsp', 'Not available')}
Zygosity: {v.get('zygosity', 'Unknown')}
Consequence: {ann.get('consequence', 'unknown')}
Impact: {ann.get('impact', 'UNKNOWN')}
SIFT: {ann.get('sift', 'N/A')}
PolyPhen: {ann.get('polyphen', 'N/A')}
gnomAD AF: {af_str}
ClinVar: {v.get('clinvar', 'Unknown')}
ACMG Classification: {v.get('acmg', 'VUS')}
Confidence Level: {v.get('confidence_level', 'Limited')}
Applied Criteria: {', '.join(applied) if applied else 'None'}
Priority: {v.get('priority', 'LOW')}
Score: {v.get('score', 0)}/10
"""
        variant_blocks.append(block)

    user_message = f"""
Generate a complete SpectralG clinical variant interpretation report.

Patient ID: {patient_id}
Patient Name: {patient_name}
Clinical Information: {clinical_info if clinical_info else 'Not provided'}
Total Variants: {len(variants)}
High Priority: {len([v for v in variants if v.get('priority') == 'HIGH'])}
Pathogenic/LP: {len([v for v in variants if v.get('acmg') in ('Pathogenic', 'Likely Pathogenic')])}

VARIANT DATA:
{''.join(variant_blocks)}

Generate the complete report following the mandatory structure in your instructions.
Be specific, evidence-based, and clinically appropriate.
State confidence levels for every classification.
PP5 is NOT applied — do not mention it as evidence.
"""

    client = anthropic.Anthropic(api_key=os.getenv("ANTHROPIC_API_KEY"))

    response = client.messages.create(
        model="claude-sonnet-4-6",  # Use Sonnet for commercial reports — higher quality
        max_tokens=2000,
        system=system_prompt,
        messages=[{"role": "user", "content": user_message}]
    )

    report_text = response.content[0].text

    # Save report
    output_file = f"{patient_id}_SpectralG_commercial_report.txt"
    with open(output_file, "w") as f:
        f.write(report_text)

    print(f"✅ Commercial report generated: {output_file}")
    print(f"Word count: {len(report_text.split())}")
    return report_text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SpectralG Commercial Report Generator")
    parser.add_argument("--json", required=True, help="Path to variants JSON file")
    parser.add_argument("--patient", required=True, help="Patient name")
    parser.add_argument("--id", required=True, help="Report/Sample ID")
    parser.add_argument("--clinical", default="", help="Clinical information")
    args = parser.parse_args()

    generate_commercial_report(args.json, args.id, args.patient, args.clinical)