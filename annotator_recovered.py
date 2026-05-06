# annotator.py
# Handles VEP annotation + ClinVar + gnomAD + scoring + ACMG

import requests
import time
import os 

# ── CACHE ───────────────────────────────────────────────
vep_cache = {}
# ── GENE LOOKUP BY POSITION (Ensembl overlap API) ───────
def get_gene_from_ncbi(chrom, pos):
    """
    Find gene name at a chromosomal position using
    Ensembl overlap API - more reliable than NCBI for this.
    """
    try:
        # Ensembl overlap endpoint - finds genes at a position
        url = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{pos}-{pos}"

        headers = {"Accept": "application/json"}
        params = {"feature": "gene", "content-type": "application/json"}

        response = requests.get(
            url,
            headers=headers,
            params=params,
            timeout=15
        )

        if response.status_code != 200:
            print(f"⚠️ Ensembl overlap failed: {response.status_code}")
            return "Unknown"

        data = response.json()

        if not data:
            # No gene exactly at this position - search wider window
            url2 = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{max(1,pos-50000)}-{pos+50000}"
            response2 = requests.get(
                url2,
                headers=headers,
                params=params,
                timeout=15
            )
            if response2.status_code == 200:
                data = response2.json()

        if not data:
            print(f"⚠️ No gene found at chr{chrom}:{pos}")
            return "Intergenic"

        # Get the gene with external name (symbol)
        for item in data:
            gene_name = item.get("external_name")
            if gene_name:
                print(f"✅ Ensembl overlap found: {gene_name} at chr{chrom}:{pos}")
                return gene_name

        return "Unknown"

    except Exception as e:
        print(f"Gene lookup error: {e}")
        return "Unknown"
    
# ── HGVS API (PRIMARY) ──────────────────────────────────
def call_vep_hgvs(chrom, pos, ref, alt):
    key = f"{chrom}-{pos}-{ref}-{alt}"

    if key in vep_cache:
        return vep_cache[key]

    try:
        # ✅ FIX: DEFINE HGVS
        hgvs = f"{chrom}:g.{pos}{ref}>{alt}"

        url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs}"

        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        params = {"canonical": 1, "sift": 1, "polyphen": 1, "numbers": 1}

        response = requests.get(url, headers=headers, params=params, timeout=20)

        if response.status_code != 200:
            return None

        data = response.json()
        vep_cache[key] = data

        return data

    except Exception as e:
        print(f"HGVS error: {e}")
        return None


# ── REGION API (FALLBACK) ───────────────────────────────
def call_vep_region(chrom, pos, ref, alt):
    key = f"{chrom}-{pos}-{ref}-{alt}-region"

    if key in vep_cache:
        return vep_cache[key]

    try:
        region = f"{chrom}:{pos}:{ref}/{alt}"
        url = "https://rest.ensembl.org/vep/human/region"

        headers = {"Content-Type": "application/json"}

        response = requests.post(
            url,
            headers=headers,
            json={"variants": [region]},
            timeout=15
        )

        if not response.ok:
            return None

        data = response.json()
        vep_cache[key] = data

        return data

    except Exception as e:
        print(f"Region error: {e}")
        return None


# ── MAIN ANNOTATION ─────────────────────────────────────
def annotate_variant(variant):
    chrom = variant["chrom"]
    pos = variant["pos"]
    ref = variant["ref"]
    alt = variant["alt"]

    # Try HGVS first
    data = call_vep_hgvs(chrom, pos, ref, alt)

    # Try region fallback
    if not data:
        print(f"⚠️ HGVS failed → using region fallback")
        data = call_vep_region(chrom, pos, ref, alt)

    # If BOTH VEP calls failed — go straight to Ensembl overlap
    if not data:
        print(f"⚠️ Both VEP calls failed → using Ensembl overlap for chr{chrom}:{pos}")
        gene_name = get_gene_from_ncbi(chrom, pos)
        variant["gene"] = gene_name
        variant["annotation"] = {
            "gene": gene_name,
            "consequence": "unknown",
            "impact": "UNKNOWN",
            "sift": None,
            "polyphen": None,
            "hgvsp": ""
        }
        variant["gnomad_af"] = None
        variant["clinvar"] = "Unknown"
        return variant

    vep_data = data[0]
    transcripts = vep_data.get("transcript_consequences", [])

    # Get gene from canonical transcript first, then any transcript
    gene = "Unknown"
    chosen_transcript = {}

    for t in transcripts:
        if t.get("canonical") == 1:
            gene = t.get("gene_symbol", "Unknown")
            chosen_transcript = t
            break

    if gene == "Unknown" and transcripts:
        gene = transcripts[0].get("gene_symbol", "Unknown")
        chosen_transcript = transcripts[0]

    # If still no gene from VEP — use Ensembl overlap
    if gene in ["Unknown", "Intergenic/Unknown", "", None]:
        print(f"⚠️ VEP returned no gene → Ensembl overlap for chr{chrom}:{pos}")
        gene = get_gene_from_ncbi(chrom, pos)

    annotation = {
        "gene": gene,
        "consequence": chosen_transcript.get("consequence_terms", ["unknown"])[0] if chosen_transcript else "unknown",
        "impact": chosen_transcript.get("impact", "UNKNOWN") if chosen_transcript else "UNKNOWN",
        "sift": chosen_transcript.get("sift_score"),
        "polyphen": chosen_transcript.get("polyphen_score"),
        "hgvsp": chosen_transcript.get("hgvsp", "")
    }

    variant["annotation"] = annotation
    variant["gene"] = gene
    variant["gnomad_af"] = extract_gnomad_af(vep_data)
    variant["clinvar"] = extract_clinvar(vep_data)

    return variant


# ── GNOMAD ──────────────────────────────────────────────
def extract_gnomad_af(vep_data):
    try:
        for var in vep_data.get("colocated_variants", []):
            if "frequencies" in var:
                for allele, data in var["frequencies"].items():
                    if "gnomad" in data:
                        return data["gnomad"]
    except:
        pass
    return None


# ── CLINVAR ─────────────────────────────────────────────
def extract_clinvar(vep_data):
    try:
        for var in vep_data.get("colocated_variants", []):
            if "clinical_significance" in var:
                return ",".join(var["clinical_significance"])
    except:
        pass
    return "Unknown"


# ── BULK ────────────────────────────────────────────────
def annotate_all(variants):
    annotated = []
    for v in variants:
        annotated.append(annotate_variant(v))
        time.sleep(1)
    return annotated


# ── FILTER ──────────────────────────────────────────────
def filter_rare_variants(variants, threshold=0.01):
    return [
        v for v in variants
        if v.get("gnomad_af") is None or float(v.get("gnomad_af", 1)) <= threshold
    ]


# ── SCORING (UPGRADED) ──────────────────────────────────
def score_variant(v):
    score = 0
    ann = v.get("annotation", {})

    # 🔥 GENE PRIORITY
    if v.get("gene") in ["TP53", "BRCA1", "BRCA2"]:
        score += 3

    # 🔥 IMPACT
    if ann.get("impact") == "HIGH":
        score += 3
    elif ann.get("impact") == "MODERATE":
        score += 2

    # 🔥 SIFT
    try:
        if ann.get("sift") is not None and float(ann.get("sift")) < 0.05:
            score += 2
    except:
        pass

    # 🔥 CLINVAR
    if v.get("clinvar") in ["Pathogenic", "Likely pathogenic"]:
        score += 4

    # 🔥 GNOMAD
    try:
        if v.get("gnomad_af") and float(v["gnomad_af"]) < 0.01:
            score += 2
    except:
        pass

    return score


def rank_variants(variants):
    for v in variants:
        score = score_variant(v)

        if score >= 6:
            priority = "HIGH"
        elif score >= 3:
            priority = "MEDIUM"
        else:
            priority = "LOW"

        v["score"] = score
        v["priority"] = priority

    return sorted(variants, key=lambda x: x["score"], reverse=True)


# ── ACMG (UPGRADED) ─────────────────────────────────────
def evaluate_acmg_evidence(v):
    evidence = []
    ann = v.get("annotation", {})

    sift = ann.get("sift")
    gnomad_af = v.get("gnomad_af")
    clinvar = v.get("clinvar")

    # PP3
    try:
        if sift is not None and float(sift) < 0.05:
            evidence.append("PP3")
    except:
        pass

    # PM2
    try:
        if gnomad_af and float(gnomad_af) < 0.01:
            evidence.append("PM2")
    except:
        pass

    # PS1
    if clinvar == "Pathogenic":
        evidence.append("PS1")

    return list(set(evidence))


def combine_acmg(evidence):
    ev = set(evidence)

    if "PS1" in ev and "PM2" in ev:
        return "Likely pathogenic"
    elif len(ev) >= 3:
        return "Likely pathogenic"
    else:
        return "VUS"


def apply_acmg_classification(variants):
    for v in variants:
        ev = evaluate_acmg_evidence(v)
        v["acmg_evidence"] = ev
        v["acmg"] = combine_acmg(ev)

    return variants

# ── NCBI FALLBACK GENE LOOKUP ───────────────────────────
# ── NCBI FALLBACK GENE LOOKUP ───────────────────────────

