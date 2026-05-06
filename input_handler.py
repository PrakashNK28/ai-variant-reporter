# input_handler.py
# SpectralG v2.0 — Universal input handler
# Handles: plain text, PDF, CSV/Excel, rsID lookup, clinical info parsing

import re
import os
import json
import requests
import io
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(dotenv_path=Path.home() / ".env", override=True)


# ── PLAIN TEXT PARSER ────────────────────────────────────────────────────────
def parse_plain_text(text, api_key=None):
    """
    Parse a plain text variant description using Claude AI.
    Returns a list of variant dicts compatible with the annotation pipeline.
    """
    try:
        import anthropic

        key = api_key or os.getenv("ANTHROPIC_API_KEY", "")
        if not key:
            return [{"error": "No API key — cannot parse plain text"}]

        client = anthropic.Anthropic(api_key=key)

        prompt = f"""Extract genetic variant information from this text and return ONLY a JSON array.

Text: {text}

Return a JSON array where each element has these fields:
- gene: gene symbol (string, e.g. "BRCA1")
- hgvsc: HGVS cDNA notation (string, e.g. "c.5266dup") or ""
- hgvsp: HGVS protein notation (string, e.g. "p.Gln1756fs") or ""
- zygosity: "Heterozygous" or "Homozygous" or "Hemizygous" or "Unknown"
- consequence: variant type (string, e.g. "frameshift_variant") or "unknown"
- acmg: classification if stated ("Pathogenic","Likely Pathogenic","VUS","Likely Benign","Benign") or "VUS"
- chrom: chromosome number as string (e.g. "17") or ""
- pos: chromosomal position as integer or 0
- ref: reference allele or ""
- alt: alternate allele or ""
- source: "plain_text"

Return ONLY the JSON array. No explanation. No markdown."""

        response = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=500,
            messages=[{"role": "user", "content": prompt}]
        )

        raw = response.content[0].text.strip()
        # Clean JSON
        raw = re.sub(r"```json|```", "", raw).strip()
        variants = json.loads(raw)

        # Ensure source field
        for v in variants:
            v["source"] = "plain_text"
            if not v.get("annotation"):
                v["annotation"] = {
                    "gene": v.get("gene", "Unknown"),
                    "consequence": v.get("consequence", "unknown"),
                    "impact": "UNKNOWN",
                    "sift": None,
                    "polyphen": None,
                    "hgvsp": v.get("hgvsp", "")
                }

        return variants

    except json.JSONDecodeError:
        return [{"error": "Could not parse AI response as JSON"}]
    except Exception as e:
        return [{"error": f"Plain text parsing failed: {str(e)}"}]


# ── PDF LAB REPORT PARSER ────────────────────────────────────────────────────
def parse_pdf_report(pdf_bytes, api_key=None):
    """
    Extract variant data from a PDF lab report using Claude AI vision.
    """
    try:
        import anthropic
        import base64

        key = api_key or os.getenv("ANTHROPIC_API_KEY", "")
        if not key:
            return [{"error": "No API key — cannot parse PDF"}]

        client = anthropic.Anthropic(api_key=key)

        # Encode PDF as base64
        pdf_b64 = base64.standard_b64encode(pdf_bytes).decode("utf-8")

        prompt = """This is a genetic laboratory report PDF.
Extract all genetic variants mentioned and return ONLY a JSON array.

Each element should have:
- gene: gene symbol
- hgvsc: HGVS cDNA notation or ""
- hgvsp: HGVS protein notation or ""
- zygosity: "Heterozygous" or "Homozygous" or "Hemizygous" or "Unknown"
- consequence: variant type or "unknown"
- acmg: ACMG classification stated in report or "VUS"
- chrom: chromosome number as string or ""
- pos: position as integer or 0
- ref: reference allele or ""
- alt: alternate allele or ""
- source: "pdf_report"

Return ONLY the JSON array. No explanation."""

        response = client.messages.create(
            model="claude-haiku-4-5-20251001",
            max_tokens=1000,
            messages=[{
                "role": "user",
                "content": [
                    {
                        "type": "document",
                        "source": {
                            "type": "base64",
                            "media_type": "application/pdf",
                            "data": pdf_b64
                        }
                    },
                    {"type": "text", "text": prompt}
                ]
            }]
        )

        raw = response.content[0].text.strip()
        raw = re.sub(r"```json|```", "", raw).strip()
        variants = json.loads(raw)

        for v in variants:
            v["source"] = "pdf_report"
            if not v.get("annotation"):
                v["annotation"] = {
                    "gene": v.get("gene", "Unknown"),
                    "consequence": v.get("consequence", "unknown"),
                    "impact": "UNKNOWN",
                    "sift": None,
                    "polyphen": None,
                    "hgvsp": v.get("hgvsp", "")
                }

        return variants

    except Exception as e:
        return [{"error": f"PDF parsing failed: {str(e)}"}]


# ── CSV / EXCEL PARSER ───────────────────────────────────────────────────────
def parse_csv_excel(file_bytes, filename):
    """
    Parse a CSV or Excel spreadsheet of variants.
    Expected columns: Gene, HGVS_c, HGVS_p, Zygosity, Classification
    """
    try:
        import pandas as pd

        if filename.endswith(".csv"):
            df = pd.read_csv(io.BytesIO(file_bytes))
        else:
            df = pd.read_excel(io.BytesIO(file_bytes))

        # Normalise column names
        df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

        variants = []
        for _, row in df.iterrows():
            v = {
                "gene":        str(row.get("gene", "Unknown")),
                "hgvsc":       str(row.get("hgvs_c", row.get("hgvsc", ""))),
                "hgvsp":       str(row.get("hgvs_p", row.get("hgvsp", ""))),
                "zygosity":    str(row.get("zygosity", "Unknown")),
                "consequence": str(row.get("consequence", "unknown")),
                "acmg":        str(row.get("classification",
                                row.get("acmg", "VUS"))),
                "chrom":       str(row.get("chromosome", row.get("chrom", ""))),
                "pos":         int(row.get("position", row.get("pos", 0)) or 0),
                "ref":         str(row.get("ref", "")),
                "alt":         str(row.get("alt", "")),
                "source":      "csv_excel"
            }
            v["annotation"] = {
                "gene":        v["gene"],
                "consequence": v["consequence"],
                "impact":      "UNKNOWN",
                "sift":        None,
                "polyphen":    None,
                "hgvsp":       v["hgvsp"]
            }
            variants.append(v)

        return variants if variants else [{"error": "No data rows found in file"}]

    except Exception as e:
        return [{"error": f"CSV/Excel parsing failed: {str(e)}"}]


# ── rsID LOOKUP ──────────────────────────────────────────────────────────────
def lookup_rsid(rsid):
    """
    Look up a dbSNP rsID and return variant data.
    """
    try:
        api_key = os.getenv("NCBI_API_KEY", "")

        # Query dbSNP via NCBI E-utilities
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            "db": "snp",
            "id": rsid.replace("rs", ""),
            "retmode": "json",
            "api_key": api_key
        }

        r = requests.get(url, params=params, timeout=15)
        if not r.ok:
            return [{"error": f"NCBI lookup failed: HTTP {r.status_code}"}]

        data = r.json()
        result = data.get("result", {})
        snp_id = rsid.replace("rs", "")
        snp_data = result.get(snp_id, {})

        if not snp_data:
            return [{"error": f"rsID {rsid} not found in dbSNP"}]

        # Extract basic info
        gene = "Unknown"
        genes = snp_data.get("genes", [])
        if genes:
            gene = genes[0].get("name", "Unknown")

        chrom = snp_data.get("chr", "")
        pos = int(snp_data.get("chrpos", 0) or 0)

        variant = {
            "gene":        gene,
            "chrom":       str(chrom),
            "pos":         pos,
            "ref":         "",
            "alt":         "",
            "hgvsc":       "",
            "hgvsp":       "",
            "zygosity":    "Unknown",
            "consequence": "unknown",
            "acmg":        "VUS",
            "source":      "rsid_lookup",
            "rsid":        rsid,
            "annotation": {
                "gene":        gene,
                "consequence": "unknown",
                "impact":      "UNKNOWN",
                "sift":        None,
                "polyphen":    None,
                "hgvsp":       ""
            }
        }

        return [variant]

    except Exception as e:
        return [{"error": f"rsID lookup failed: {str(e)}"}]


# ── CLINICAL INFO PARSER ─────────────────────────────────────────────────────
def parse_clinical_info(raw_dict):
    """
    Validate and structure clinical info from sidebar inputs.
    Applies SpectralG tone rules:
    - Never infer sex from name
    - Never say 'not provided' by inference
    - Flag missing data explicitly
    """
    sex_raw = raw_dict.get("sex", "Not provided")

    # Never infer sex from patient name
    if sex_raw in ["Not provided", "Not specified", ""]:
        sex_final = "Sex: not provided in referral"
    else:
        sex_final = sex_raw

    return {
        "patient_id":    raw_dict.get("patient_id", "SAMPLE_001"),
        "patient_name":  raw_dict.get("patient_name", "[Not provided]"),
        "age":           raw_dict.get("age", "Not provided"),
        "sex":           sex_final,
        "indication":    raw_dict.get("indication", "Not provided"),
        "clinical_features": raw_dict.get("clinical_features", "Not provided"),
        "family_history":    raw_dict.get("family_history", "Not provided"),
        "report_type":       raw_dict.get("report_type", "Clinical WES"),
        "referring_clinician": raw_dict.get("referring_clinician", "Not provided"),
        "genotype_phenotype_correlation": raw_dict.get(
            "genotype_phenotype_correlation", "Not assessed"),
        "gp_narrative":  raw_dict.get("gp_narrative", ""),
        "action_points": raw_dict.get("action_points", []),
        "mdt_plan":      raw_dict.get("mdt_plan", []),
    }