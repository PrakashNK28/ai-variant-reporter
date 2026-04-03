# annotator.py
# Handles VEP annotation + ClinVar + gnomAD + scoring + ACMG

import requests
import time

# ── CACHE ───────────────────────────────────────────────
vep_cache = {}

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

    # 🔥 TRY HGVS
    data = call_vep_hgvs(chrom, pos, ref, alt)

    # 🔥 FALLBACK
    if not data:
        print("⚠️ HGVS failed → using region fallback")
        data = call_vep_region(chrom, pos, ref, alt)

    if not data:
        return variant

    vep_data = data[0]

    transcript = vep_data.get("transcript_consequences", [{}])[0]

    # ✅ FIX: SAFE GENE HANDLING
    gene = transcript.get("gene_symbol")
    if not gene:
        gene = "Intergenic/Unknown"

    annotation = {
        "gene": gene,
        "consequence": transcript.get("consequence_terms", ["unknown"])[0],
        "impact": transcript.get("impact", "UNKNOWN"),
        "sift": transcript.get("sift_score"),
        "polyphen": transcript.get("polyphen_score"),
        "hgvsp": transcript.get("hgvsp", "")
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