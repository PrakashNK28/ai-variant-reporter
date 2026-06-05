# annotator.py
# SpectralG — Clinical Variant Annotator v2.1
# Research-grade: full ACMG table, HGVS notation, evidence panel, gnomAD SAS
#
# Framework: ACMG/AMP 2015 (Richards et al., Genet Med 2015)
# PP5/BP6 NOT applied per ACMG 2023 (Biesecker & Harrison)

import re
import requests
import time
import os
from pathlib import Path
from dotenv import load_dotenv

# Indian Variant Registry — import with fallback if file missing
try:
    from indian_variant_registry import lookup_indian_evidence
    REGISTRY_AVAILABLE = True
except ImportError:
    REGISTRY_AVAILABLE = False
    print("⚠️ Indian variant registry not found — registry lookup disabled")

load_dotenv(dotenv_path=Path.home() / ".env", override=True)

# ── CACHE ─────────────────────────────────────────────────────────────────────
vep_cache = {}

# ── KNOWN DISEASE GENES (PP2 criterion) ──────────────────────────────────────
KNOWN_DISEASE_GENES = {
    "TP53","BRCA1","BRCA2","MLH1","MSH2","MSH6","PMS2",
    "APC","PTEN","STK11","CDH1","PALB2","CHEK2","ATM",
    "CFTR","HBB","HBA1","HBA2","LDLR","APOB","PCSK9",
    "PKD1","PKD2","TSC1","TSC2","NF1","NF2","VHL",
    "RB1","WT1","RUNX1","FLT3","KRAS","NRAS","BRAF",
    "EGFR","ALK","RET","MEN1","GJB2","GJB6","SLC26A4",
    "MYO7A","OTOF","HEXA","HEXB","GBA","ASPA","ARSA",
    "DMD","MECP2","FMR1","SMN1","DMPK","MAP1A","OCRL",
    "MYBPC3","MYH7","KCNQ1","SCN5A","SCN1A","KCNQ2",
    "DEPDC5","CYP2D6","CYP2C19","DPYD","TPMT"
}

HIGH_IMPACT_CONSEQUENCES = {
    "stop_gained","frameshift_variant","splice_acceptor_variant",
    "splice_donor_variant","start_lost","stop_lost",
    "transcript_ablation","transcript_amplification"
}

MODERATE_IMPACT_CONSEQUENCES = {
    "missense_variant","inframe_insertion","inframe_deletion",
    "protein_altering_variant","regulatory_region_variant"
}


# ── VEP HGVS (PRIMARY) ────────────────────────────────────────────────────────
def call_vep_hgvs(chrom, pos, ref, alt):
    key = f"{chrom}-{pos}-{ref}-{alt}"
    if key in vep_cache:
        return vep_cache[key]
    try:
        hgvs    = f"{chrom}:g.{pos}{ref}>{alt}"
        url     = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs}"
        headers = {"Content-Type":"application/json","Accept":"application/json"}
        params  = {"canonical":1,"sift":1,"polyphen":1,"numbers":1,"hgvs":1,"domains":1}
        r = requests.get(url, headers=headers, params=params, timeout=20)
        if r.status_code != 200:
            return None
        data = r.json()
        vep_cache[key] = data
        return data
    except Exception as e:
        print(f"VEP HGVS error: {e}")
        return None


# ── VEP REGION (FALLBACK) ─────────────────────────────────────────────────────
def call_vep_region(chrom, pos, ref, alt):
    key = f"{chrom}-{pos}-{ref}-{alt}-region"
    if key in vep_cache:
        return vep_cache[key]
    try:
        region  = f"{chrom} {pos} . {ref} {alt} . . ."
        url     = "https://rest.ensembl.org/vep/human/region"
        headers = {"Content-Type":"application/json","Accept":"application/json"}
        r = requests.post(url, headers=headers,
                          json={"variants":[region],"canonical":1,
                                "sift":1,"polyphen":1,"hgvs":1},
                          timeout=20)
        if not r.ok:
            return None
        data = r.json()
        vep_cache[key] = data
        return data
    except Exception as e:
        print(f"VEP region error: {e}")
        return None


# ── ENSEMBL OVERLAP (GENE NAME FALLBACK) ──────────────────────────────────────
def get_gene_from_ensembl_overlap(chrom, pos):
    """Find gene name at position when VEP returns no gene symbol."""
    try:
        url     = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{pos}-{pos}"
        headers = {"Accept":"application/json"}
        params  = {"feature":"gene","content-type":"application/json"}
        r = requests.get(url, headers=headers, params=params, timeout=15)
        if r.status_code == 200:
            for item in r.json():
                name = item.get("external_name")
                if name:
                    print(f"✅ Overlap found: {name} at chr{chrom}:{pos}")
                    return name
        # Widen window
        url2 = (f"https://rest.ensembl.org/overlap/region/human/"
                f"{chrom}:{max(1,pos-50000)}-{pos+50000}")
        r2 = requests.get(url2, headers=headers, params=params, timeout=15)
        if r2.status_code == 200:
            for item in r2.json():
                name = item.get("external_name")
                if name:
                    print(f"✅ Overlap (window) found: {name} at chr{chrom}:{pos}")
                    return name
    except Exception as e:
        print(f"Overlap error: {e}")
    return "Unknown"


# ── DIRECT CLINVAR LOOKUP ─────────────────────────────────────────────────────
# Defined BEFORE annotate_variant which calls it
def lookup_clinvar_direct(chrom, pos, ref, alt):
    """
    Direct ClinVar lookup via NCBI E-utilities.
    Called when VEP does not return ClinVar colocated data.
    """
    try:
        api_key = os.getenv("NCBI_API_KEY", "")
        base    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        term    = f"{chrom}[chr] AND {pos}[chrpos37] AND human[orgn]"

        r = requests.get(f"{base}esearch.fcgi", params={
            "db": "clinvar", "term": term,
            "retmax": 3, "retmode": "json", "api_key": api_key
        }, timeout=10)
        if not r.ok:
            return "Unknown"

        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return "Unknown"

        r2 = requests.get(f"{base}esummary.fcgi", params={
            "db": "clinvar", "id": ids[0],
            "retmode": "json", "api_key": api_key
        }, timeout=10)
        if not r2.ok:
            return "Unknown"

        result = r2.json().get("result", {}).get(ids[0], {})
        sig    = result.get("clinical_significance", {})
        if isinstance(sig, dict):
            return sig.get("description", "Unknown")
        return str(sig) if sig else "Unknown"

    except Exception as e:
        print(f"ClinVar direct lookup error: {e}")
        return "Unknown"


# ── EXTRACT HGVS NOTATION ─────────────────────────────────────────────────────
def extract_hgvs(transcript):
    """
    Extract HGVS c. and p. from VEP transcript consequence.
    Strips Ensembl transcript prefix if present.
    Returns (hgvsc, hgvsp) — never empty strings.
    """
    hgvsc = transcript.get("hgvsc") or transcript.get("hgvs_c") or ""
    hgvsp = transcript.get("hgvsp") or transcript.get("hgvs_p") or ""
    if hgvsc and ":" in hgvsc:
        hgvsc = hgvsc.split(":")[-1]
    if hgvsp and ":" in hgvsp:
        hgvsp = hgvsp.split(":")[-1]
    return (hgvsc or "Not available (VEP missing)",
            hgvsp or "Not available (VEP missing)")


# ── EXTRACT GNOMAD FREQUENCIES ────────────────────────────────────────────────
def extract_gnomad_af(vep_data):
    """
    Extract gnomAD frequencies from VEP colocated_variants.
    Prioritises South Asian (SAS) subpopulation for Indian patients.
    """
    try:
        for var in vep_data.get("colocated_variants", []):
            freqs = var.get("frequencies", {})
            if not freqs:
                continue
            for allele, fd in freqs.items():
                if not isinstance(fd, dict):
                    continue
                sas = fd.get("gnomad_sas") or fd.get("gnomadg_sas") or fd.get("sas")
                glb = fd.get("gnomad")     or fd.get("gnomadg")    or fd.get("af")
                if sas is not None or glb is not None:
                    return {
                        "global":      glb,
                        "south_asian": sas,
                        "display":     sas if sas is not None else glb,
                        "source":      "VEP colocated_variants"
                    }
    except Exception as e:
        print(f"gnomAD extraction error: {e}")
    return None


# ── EXTRACT CLINVAR FROM VEP ──────────────────────────────────────────────────
def extract_clinvar(vep_data):
    try:
        for var in vep_data.get("colocated_variants", []):
            if "clinical_significance" in var:
                sigs = var["clinical_significance"]
                return ", ".join(sigs) if isinstance(sigs, list) else str(sigs)
    except Exception as e:
        print(f"ClinVar extraction error: {e}")
    return "Unknown"


# ── BUILD FULL ACMG CRITERIA TABLE ────────────────────────────────────────────
def build_acmg_criteria_table(variant, annotation, gnomad_af, clinvar):
    """
    Build complete ACMG/AMP 2015 criteria table.
    ALL 28 criteria included — applied True or False.
    PP5 and BP6 always False per ACMG 2023.
    Each entry: {code, weight, applied, evidence}
    """
    gene        = annotation.get("gene", "Unknown")
    consequence = annotation.get("consequence", "unknown")
    sift        = annotation.get("sift")
    polyphen    = annotation.get("polyphen")

    # Resolve AF — prefer South Asian
    af = None
    af_str = "Not available in gnomAD"
    if isinstance(gnomad_af, dict):
        af      = gnomad_af.get("south_asian") or gnomad_af.get("global")
        sas_val = gnomad_af.get("south_asian")
        glb_val = gnomad_af.get("global")
        af_str  = (f"gnomAD SAS: {sas_val if sas_val is not None else 'N/A'} | "
                   f"gnomAD Global: {glb_val if glb_val is not None else 'N/A'}")
    elif gnomad_af is not None:
        af     = float(gnomad_af)
        af_str = f"gnomAD AF: {af:.6f}"

    # ── PVS1 ─────────────────────────────────────────────────────────────────
    pvs1 = consequence in HIGH_IMPACT_CONSEQUENCES
    pvs1_ev = (
        f"Loss-of-function consequence '{consequence}' identified. "
        f"Gene: {gene}. PVS1 applied when LOF is established disease mechanism."
        if pvs1 else
        f"Consequence '{consequence}' is not loss-of-function. PVS1 not applied."
    )

    # ── PS1 ──────────────────────────────────────────────────────────────────
    ps1    = False
    ps1_ev = ("Requires comparison to established pathogenic amino acid changes. "
              "Cannot determine computationally — manual curation required. Not applied.")

    # ── PS2 ──────────────────────────────────────────────────────────────────
    ps2    = False
    ps2_ev = "Parental testing data not available. PS2 not applied."

    # ── PS3 ──────────────────────────────────────────────────────────────────
    ps3    = False
    ps3_ev = "Functional study data not available computationally. PS3 not applied."

    # ── PS4 — uses Indian Variant Registry ───────────────────────────────────
    indian_ev = variant.get("indian_registry")
    if indian_ev and indian_ev.get("observation_count", 0) >= 5:
        ps4    = True
        ps4_ev = (
            f"Variant observed in {indian_ev['observation_count']} "
            f"Indian patients via SpectralG Indian Variant Registry. "
            f"Consensus: {indian_ev['consensus']}. "
            f"PS4 applied — requires expert review before clinical reporting."
        )
    elif indian_ev and indian_ev.get("observation_count", 0) >= 3:
        ps4    = False
        ps4_ev = (
            f"Variant observed in {indian_ev['observation_count']} Indian patients "
            f"(minimum 5 required for PS4). Consensus: {indian_ev['consensus']}. "
            f"Accumulating evidence — PS4 not yet applied."
        )
    else:
        ps4    = False
        ps4_ev = (
            "Case-control data not available computationally. "
            "Not yet observed in SpectralG Indian Variant Registry. "
            "PS4 not applied."
        )

    # ── PM1 ──────────────────────────────────────────────────────────────────
    pm1    = False
    pm1_ev = ("Domain-level annotation requires ClinGen/UniProt data. "
              "Not determined from VEP alone. PM1 not applied.")
    try:
        from gene_specific_rules import get_gene_rule
        gene_rule      = get_gene_rule(gene)
        hotspot_codons = gene_rule.get("hotspot_codons", [])
        hgvsp_str      = variant.get("hgvsp", "") or annotation.get("hgvsp", "")
        codon_num      = None
        match = re.search(r'[A-Za-z]+(\d+)[A-Za-z]', hgvsp_str)
        if match:
            codon_num = int(match.group(1))
        if hotspot_codons and codon_num and codon_num in hotspot_codons:
            pm1    = True
            pm1_ev = (
                f"Codon {codon_num} in {gene} is a defined mutational hotspot "
                f"per ClinGen VCEP. Hotspot codons: {hotspot_codons}. "
                f"PM1 applied per {gene_rule.get('source','ClinGen VCEP')}."
            )
        elif hotspot_codons:
            pm1    = False
            pm1_ev = (
                f"Codon {codon_num} not in hotspot list for {gene} {hotspot_codons}. PM1 not applied."
                if codon_num else
                f"Codon not extractable from HGVS notation. PM1 not applied."
            )
    except Exception as e:
        pm1    = False
        pm1_ev = f"PM1 gene-specific check unavailable: {e}. Not applied."

    # ── PM2 ──────────────────────────────────────────────────────────────────
    if af is None:
        pm2    = True
        pm2_ev = f"Absent from gnomAD ({af_str}). PM2 applied — absence supports pathogenicity."
    elif af < 0.0001:
        pm2    = True
        pm2_ev = f"Extremely rare ({af_str}). AF {af:.6f} < 0.0001. PM2 applied."
    elif af < 0.01:
        pm2    = True
        pm2_ev = f"Rare in population ({af_str}). AF {af:.5f} < 0.01. PM2 applied."
    else:
        pm2    = False
        pm2_ev = f"AF {af:.5f} ({af_str}) exceeds rare threshold. PM2 not applied."

    # ── PM3 ──────────────────────────────────────────────────────────────────
    pm3    = False
    pm3_ev = "Phase and family data not available. PM3 not applied."

    # ── PM4 ──────────────────────────────────────────────────────────────────
    pm4    = consequence in {"inframe_insertion", "inframe_deletion"}
    pm4_ev = (
        f"In-frame indel '{consequence}' causes protein length change. PM4 applied."
        if pm4 else
        f"Consequence '{consequence}' is not an in-frame indel. PM4 not applied."
    )

    # ── PM5 ──────────────────────────────────────────────────────────────────
    pm5    = False
    pm5_ev = ("Requires variant database comparison at amino acid position. "
              "Cannot determine computationally. PM5 not applied.")

    # ── PM6 ──────────────────────────────────────────────────────────────────
    pm6    = False
    pm6_ev = "Parental data not available. PM6 not applied."

    # ── PP1 ──────────────────────────────────────────────────────────────────
    pp1    = False
    pp1_ev = "Family segregation data not available. PP1 not applied."

    # ── PP2 ──────────────────────────────────────────────────────────────────
    pp2    = (consequence == "missense_variant") and (gene in KNOWN_DISEASE_GENES)
    pp2_ev = (
        f"Missense variant in {gene} — gene with established missense disease mechanism. PP2 applied."
        if pp2 else
        f"Gene '{gene}' not in curated disease gene list or not missense. PP2 not applied."
    )

    # ── PP3 ──────────────────────────────────────────────────────────────────
    sift_dam = sift is not None and float(sift) <= 0.05
    poly_dam = polyphen is not None and float(polyphen) >= 0.85
    pp3      = sift_dam or poly_dam
    sift_str = (f"SIFT {sift:.3f} ({'deleterious' if sift_dam else 'tolerated'})"
                if sift is not None else "SIFT N/A")
    poly_str = (f"PolyPhen {polyphen:.3f} ({'damaging' if poly_dam else 'benign'})"
                if polyphen is not None else "PolyPhen N/A")
    pp3_ev   = (
        f"{sift_str} | {poly_str}. "
        f"PP3 {'applied — at least one tool predicts damaging' if pp3 else 'not applied'}."
    )

    # ── PP4 ──────────────────────────────────────────────────────────────────
    pp4    = False
    pp4_ev = "Clinical phenotype not provided. PP4 not applied."

    # ── PP5 — INTENTIONALLY NOT APPLIED ──────────────────────────────────────
    pp5    = False
    pp5_ev = ("PP5 INTENTIONALLY NOT APPLIED per ACMG 2023 (Biesecker & Harrison). "
              "External laboratory assertions excluded as independent evidence. "
              f"ClinVar data documented separately: {clinvar}.")

    # ── BA1 ──────────────────────────────────────────────────────────────────
    ba1    = af is not None and float(af) > 0.05
    ba1_ev = (
        f"AF {af:.4f} ({af_str}) > 5% threshold. BA1 applied — common variant, likely benign."
        if ba1 else
        f"AF ({af_str}) does not exceed 5% threshold. BA1 not applied."
    )

    # ── BS1 ──────────────────────────────────────────────────────────────────
    bs1    = af is not None and 0.01 < float(af) <= 0.05
    bs1_ev = (
        f"AF {af:.4f} ({af_str}) > 1% threshold. BS1 applied."
        if bs1 else
        f"AF ({af_str}) does not trigger BS1 threshold. Not applied."
    )

    # ── BS2 ──────────────────────────────────────────────────────────────────
    bs2    = False
    bs2_ev = "Healthy adult observation data not available. BS2 not applied."

    # ── BS3 ──────────────────────────────────────────────────────────────────
    bs3    = False
    bs3_ev = "Functional study data not available. BS3 not applied."

    # ── BS4 ──────────────────────────────────────────────────────────────────
    bs4    = False
    bs4_ev = "Family data not available. BS4 not applied."

    # ── BP1 ──────────────────────────────────────────────────────────────────
    bp1    = False
    bp1_ev = ("Requires ClinGen gene-disease validity data — cannot determine computationally. "
              "BP1 not applied.")

    # ── BP2 ──────────────────────────────────────────────────────────────────
    bp2    = False
    bp2_ev = "Phase and family data not available. BP2 not applied."

    # ── BP3 ──────────────────────────────────────────────────────────────────
    bp3    = False
    bp3_ev = "Repeat region annotation not available from VEP alone. BP3 not applied."

    # ── BP4 ──────────────────────────────────────────────────────────────────
    sift_tol = sift is not None and float(sift) > 0.05
    poly_ben = polyphen is not None and float(polyphen) < 0.45
    bp4      = sift_tol and poly_ben
    bp4_ev   = (
        f"{sift_str} | {poly_str}. "
        f"BP4 {'applied — both tools predict benign/tolerated' if bp4 else 'not applied'}."
    )

    # ── BP5 ──────────────────────────────────────────────────────────────────
    bp5    = False
    bp5_ev = "Alternate molecular diagnosis data not available. BP5 not applied."

    # ── BP6 — NOT APPLIED ────────────────────────────────────────────────────
    bp6    = False
    bp6_ev = "BP6 NOT APPLIED per ACMG 2023 — analogous to PP5. Excluded."

    # ── BP7 ──────────────────────────────────────────────────────────────────
    bp7    = consequence == "synonymous_variant"
    bp7_ev = (
        f"Synonymous variant ('{consequence}') — no amino acid change. BP7 applied."
        if bp7 else
        f"Consequence '{consequence}' is not synonymous. BP7 not applied."
    )

    return [
        {"code":"PVS1","weight":"Very Strong Pathogenic",  "applied":pvs1,"evidence":pvs1_ev},
        {"code":"PS1", "weight":"Strong Pathogenic",       "applied":ps1, "evidence":ps1_ev},
        {"code":"PS2", "weight":"Strong Pathogenic",       "applied":ps2, "evidence":ps2_ev},
        {"code":"PS3", "weight":"Strong Pathogenic",       "applied":ps3, "evidence":ps3_ev},
        {"code":"PS4", "weight":"Strong Pathogenic",       "applied":ps4, "evidence":ps4_ev},
        {"code":"PM1", "weight":"Moderate Pathogenic",     "applied":pm1, "evidence":pm1_ev},
        {"code":"PM2", "weight":"Moderate Pathogenic",     "applied":pm2, "evidence":pm2_ev},
        {"code":"PM3", "weight":"Moderate Pathogenic",     "applied":pm3, "evidence":pm3_ev},
        {"code":"PM4", "weight":"Moderate Pathogenic",     "applied":pm4, "evidence":pm4_ev},
        {"code":"PM5", "weight":"Moderate Pathogenic",     "applied":pm5, "evidence":pm5_ev},
        {"code":"PM6", "weight":"Moderate Pathogenic",     "applied":pm6, "evidence":pm6_ev},
        {"code":"PP1", "weight":"Supporting Pathogenic",   "applied":pp1, "evidence":pp1_ev},
        {"code":"PP2", "weight":"Supporting Pathogenic",   "applied":pp2, "evidence":pp2_ev},
        {"code":"PP3", "weight":"Supporting Pathogenic",   "applied":pp3, "evidence":pp3_ev},
        {"code":"PP4", "weight":"Supporting Pathogenic",   "applied":pp4, "evidence":pp4_ev},
        {"code":"PP5", "weight":"NOT APPLIED (ACMG 2023)", "applied":pp5, "evidence":pp5_ev},
        {"code":"BA1", "weight":"Stand-alone Benign",      "applied":ba1, "evidence":ba1_ev},
        {"code":"BS1", "weight":"Strong Benign",           "applied":bs1, "evidence":bs1_ev},
        {"code":"BS2", "weight":"Strong Benign",           "applied":bs2, "evidence":bs2_ev},
        {"code":"BS3", "weight":"Strong Benign",           "applied":bs3, "evidence":bs3_ev},
        {"code":"BS4", "weight":"Strong Benign",           "applied":bs4, "evidence":bs4_ev},
        {"code":"BP1", "weight":"Supporting Benign",       "applied":bp1, "evidence":bp1_ev},
        {"code":"BP2", "weight":"Supporting Benign",       "applied":bp2, "evidence":bp2_ev},
        {"code":"BP3", "weight":"Supporting Benign",       "applied":bp3, "evidence":bp3_ev},
        {"code":"BP4", "weight":"Supporting Benign",       "applied":bp4, "evidence":bp4_ev},
        {"code":"BP5", "weight":"Supporting Benign",       "applied":bp5, "evidence":bp5_ev},
        {"code":"BP6", "weight":"NOT APPLIED (ACMG 2023)", "applied":bp6, "evidence":bp6_ev},
        {"code":"BP7", "weight":"Supporting Benign",       "applied":bp7, "evidence":bp7_ev},
    ]


# ── COMBINE ACMG TABLE INTO CLASSIFICATION ────────────────────────────────────
def combine_acmg_from_table(criteria_table):
    """
    Combine criteria into classification per Richards et al. 2015, Table 5.
    Returns (classification, confidence_level, evidence_list).
    """
    applied = {c["code"] for c in criteria_table
               if c.get("applied") and c["code"] not in {"PP5","BP6"}}

    pvs = "PVS1" in applied
    ps  = len(applied & {"PS1","PS2","PS3","PS4"})
    pm  = len(applied & {"PM1","PM2","PM3","PM4","PM5","PM6"})
    pp  = len(applied & {"PP1","PP2","PP3","PP4"})
    ba  = "BA1" in applied
    bs  = len(applied & {"BS1","BS2","BS3","BS4"})
    bp  = len(applied & {"BP1","BP2","BP3","BP4","BP5","BP7"})

    if ba:
        cls, conf = "Benign", "High"
    elif ((pvs and ps >= 1) or (pvs and pm >= 2) or
          (pvs and pp >= 2) or (ps >= 2)):
        cls, conf = "Pathogenic", "High"
    elif ((pvs and pm == 1) or (ps == 1 and pm >= 3) or
          (ps == 1 and pp >= 4) or (pm >= 3 and pp >= 2)):
        cls, conf = "Likely Pathogenic", "Moderate"
    elif bs >= 2 or (bs == 1 and bp >= 1):
        cls, conf = ("Benign" if bs >= 2 else "Likely Benign"), "Moderate"
    elif bp >= 2:
        cls, conf = "Likely Benign", "Moderate"
    else:
        cls  = "VUS"
        n    = len(applied)
        conf = "Moderate" if n >= 3 else "Limited" if n >= 1 else "Insufficient"

    return cls, conf, sorted(applied)


# ── BUILD EVIDENCE PANEL ──────────────────────────────────────────────────────
def build_evidence_panel(variant, annotation, gnomad_af, clinvar):
    """
    Build structured evidence panel (3billion/Sbimon-style).
    Always fully populated — uses explicit N/A statements when data missing.
    """
    gene        = annotation.get("gene", "Unknown")
    consequence = annotation.get("consequence", "unknown")
    impact      = annotation.get("impact", "UNKNOWN")
    sift        = annotation.get("sift")
    polyphen    = annotation.get("polyphen")
    hgvsc       = variant.get("hgvsc", "Not available")
    hgvsp       = variant.get("hgvsp", "Not available")

    if isinstance(gnomad_af, dict):
        sas     = gnomad_af.get("south_asian")
        glb     = gnomad_af.get("global")
        sas_str = f"{sas:.6f}" if sas is not None else "Not available in VEP response"
        glb_str = f"{glb:.6f}" if glb is not None else "Not available in VEP response"
    else:
        sas_str = glb_str = "Not available in VEP response"

    sift_str = (f"{sift:.3f} ({'Deleterious ≤0.05' if sift <= 0.05 else 'Tolerated >0.05'})"
                if sift is not None else "Not available")
    poly_str = (f"{polyphen:.3f} ({'Probably Damaging ≥0.85' if polyphen >= 0.85 else 'Not Damaging <0.85'})"
                if polyphen is not None else "Not available")

    indian_reg = variant.get("indian_registry")
    if indian_reg:
        indian_str = (
            f"SpectralG Indian Variant Registry: "
            f"{indian_reg.get('observation_count', 0)} observations in Indian patients. "
            f"Consensus: {indian_reg.get('consensus', 'No data')}. "
            f"Registry provides PS4 evidence when ≥5 observations present."
        )
    else:
        indian_str = (
            "Not yet observed in SpectralG Indian Variant Registry. "
            "Registry grows with each interpreted Indian patient case. "
            "Submit deidentified cases to contribute to Indian population data."
        )

    return {
        "Population Data": (
            f"gnomAD South Asian (SAS): {sas_str} | "
            f"gnomAD Global: {glb_str} | "
            "South Asian frequency prioritised for Indian patients per SpectralG design."
        ),
        "Predicted Consequence": (
            f"Gene: {gene} | Consequence: {consequence} | "
            f"VEP Impact: {impact} | HGVS c.: {hgvsc} | HGVS p.: {hgvsp}"
        ),
        "Computational Evidence": (
            f"SIFT: {sift_str} | PolyPhen-2: {poly_str} | "
            "PP3 applied when SIFT ≤0.05 OR PolyPhen ≥0.85 (ACMG 2015). "
            "BP4 applied when both tools predict tolerated/benign."
        ),
        "ClinVar Significance": (
            f"{clinvar} | "
            "Documented for reference only. PP5 NOT applied per ACMG 2023 guidance. "
            "ClinVar assertion does not count as independent pathogenic evidence."
        ),
        "Segregation Data": (
            "Not provided. Parental testing recommended to assess de novo vs inherited status. "
            "De novo confirmation would apply PS2 (Strong Pathogenic evidence)."
        ),
        "Previously Reported": (
            "ClinVar colocated variant data used where available from VEP response. "
            "Manual PubMed literature review recommended for comprehensive assessment."
        ),
        "Functional Studies": (
            "No functional study data available computationally. "
            "PS3/BS3 criteria require wet-lab functional evidence — manual curation required."
        ),
        "Indian Population Evidence": indian_str,
        "Sanger Validation": (
            "Sanger sequencing confirmation recommended before clinical reporting "
            "per standard laboratory protocols and ACMG/AMP reporting guidelines."
        ),
    }


# ── MAIN ANNOTATION FUNCTION ──────────────────────────────────────────────────
def annotate_variant(variant):
    """
    Complete annotation pipeline for a single variant.
    Steps:
    1. VEP HGVS (primary)
    2. VEP Region (fallback)
    3. Ensembl Overlap (gene name fallback)
    4. Extract HGVS c./p., gnomAD SAS, ClinVar
    5. Direct ClinVar NCBI lookup if VEP returns nothing
    6. Check Indian Variant Registry
    7. Build full ACMG criteria table (28 criteria, always populated)
    8. Build evidence panel (9 categories, always populated)
    """
    chrom = str(variant.get("chrom", ""))
    pos   = int(variant.get("pos", 0))
    ref   = str(variant.get("ref", ""))
    alt   = str(variant.get("alt", ""))

    # ── Steps 1 + 2: VEP calls ───────────────────────────────────────────────
    data = call_vep_hgvs(chrom, pos, ref, alt)
    if not data:
        print(f"⚠️ HGVS failed → region fallback chr{chrom}:{pos}")
        data = call_vep_region(chrom, pos, ref, alt)

    if not data:
        # ── Step 3: Overlap fallback (VEP completely failed) ─────────────────
        print(f"⚠️ Both VEP failed → Ensembl Overlap chr{chrom}:{pos}")
        gene_name  = get_gene_from_ensembl_overlap(chrom, pos)
        hgvsc      = "Not available (VEP failed)"
        hgvsp      = "Not available (VEP failed)"
        annotation = {
            "gene": gene_name, "consequence": "unknown",
            "impact": "UNKNOWN", "sift": None, "polyphen": None,
            "hgvsc": hgvsc, "hgvsp": hgvsp
        }
        gnomad_af = None
        clinvar   = "Unknown"
        variant.update({
            "annotation": annotation, "gene": gene_name,
            "hgvsc": hgvsc, "hgvsp": hgvsp,
            "gnomad_af": None, "clinvar": "Unknown",
            "vep_status": "failed"
        })

    else:
        # ── VEP succeeded — extract all data ─────────────────────────────────
        vep_data    = data[0]
        transcripts = vep_data.get("transcript_consequences", [])
        most_severe = vep_data.get("most_severe_consequence", "unknown")

        gene_name          = "Unknown"
        hgvsc = hgvsp      = "Not available (VEP missing)"
        sift_score         = None
        polyphen_score     = None
        chosen_consequence = most_severe
        chosen_impact      = "UNKNOWN"

        # Priority 1: canonical transcript
        for t in transcripts:
            if t.get("canonical") == 1:
                gene_name          = t.get("gene_symbol", "Unknown")
                hgvsc, hgvsp       = extract_hgvs(t)
                sift_score         = t.get("sift_score")
                polyphen_score     = t.get("polyphen_score")
                chosen_consequence = t.get("consequence_terms", [most_severe])[0]
                chosen_impact      = t.get("impact", "UNKNOWN")
                break

        # Priority 2: first transcript if no canonical found
        if gene_name == "Unknown" and transcripts:
            t                  = transcripts[0]
            gene_name          = t.get("gene_symbol", "Unknown")
            hgvsc, hgvsp       = extract_hgvs(t)
            sift_score         = t.get("sift_score")
            polyphen_score     = t.get("polyphen_score")
            chosen_consequence = t.get("consequence_terms", [most_severe])[0]
            chosen_impact      = t.get("impact", "UNKNOWN")

        # Priority 3: Overlap fallback if gene still unknown
        if gene_name in ("Unknown", "", None):
            print(f"⚠️ VEP no gene → Overlap chr{chrom}:{pos}")
            gene_name = get_gene_from_ensembl_overlap(chrom, pos)

        # Infer impact from consequence if VEP did not return it
        if chosen_impact == "UNKNOWN":
            if chosen_consequence in HIGH_IMPACT_CONSEQUENCES:
                chosen_impact = "HIGH"
            elif chosen_consequence in MODERATE_IMPACT_CONSEQUENCES:
                chosen_impact = "MODERATE"
            else:
                chosen_impact = "LOW"

        annotation = {
            "gene": gene_name, "consequence": chosen_consequence,
            "impact": chosen_impact, "sift": sift_score,
            "polyphen": polyphen_score,
            "hgvsc": hgvsc, "hgvsp": hgvsp
        }

        gnomad_af = extract_gnomad_af(vep_data)
        clinvar   = extract_clinvar(vep_data)

        # Step 5: Direct ClinVar NCBI lookup if VEP returned nothing
        if clinvar == "Unknown":
            print(f"ClinVar not in VEP → direct NCBI lookup chr{chrom}:{pos}")
            clinvar = lookup_clinvar_direct(chrom, pos, ref, alt)

        variant.update({
            "annotation": annotation, "gene": gene_name,
            "hgvsc": hgvsc, "hgvsp": hgvsp,
            "gnomad_af": gnomad_af, "clinvar": clinvar,
            "vep_status": "success"
        })

    # ── Step 6: Indian Variant Registry lookup ───────────────────────────────
    # This runs AFTER both VEP success and VEP failure paths
    current_gene  = variant.get("gene", "Unknown")
    current_hgvsc = variant.get("hgvsc", "")
    if REGISTRY_AVAILABLE and current_gene not in ("Unknown", ""):
        try:
            indian_evidence = lookup_indian_evidence(current_gene, current_hgvsc)
            if indian_evidence and indian_evidence.get("observation_count", 0) >= 3:
                variant["indian_registry"] = indian_evidence
                print(f"✅ Indian registry hit: {current_gene} {current_hgvsc} — "
                      f"{indian_evidence['observation_count']} observations — "
                      f"consensus: {indian_evidence['consensus']}")
            else:
                variant["indian_registry"] = None
        except Exception as e:
            print(f"Registry lookup error: {e}")
            variant["indian_registry"] = None
    else:
        variant["indian_registry"] = None

    # ── Steps 7 + 8: ACMG table + evidence panel ─────────────────────────────
    ann = variant["annotation"]
    gaf = variant.get("gnomad_af")
    cv  = variant.get("clinvar", "Unknown")

    criteria_table          = build_acmg_criteria_table(variant, ann, gaf, cv)
    acmg_cls, confidence, evidence_list = combine_acmg_from_table(criteria_table)

    variant.update({
        "acmg":                acmg_cls,
        "confidence_level":    confidence,
        "acmg_evidence":       evidence_list,
        "acmg_criteria_table": criteria_table,
        "evidence_panel":      build_evidence_panel(variant, ann, gaf, cv)
    })

    return variant


# ── BULK ANNOTATION ───────────────────────────────────────────────────────────
def annotate_all(variants):
    """Annotate list of variants with 1.1s rate limit between Ensembl API calls."""
    annotated = []
    for i, v in enumerate(variants):
        print(f"Annotating {i+1}/{len(variants)}: chr{v.get('chrom')}:{v.get('pos')}")
        annotated.append(annotate_variant(v))
        time.sleep(1.1)
    return annotated


def enrich_gnomad_sas(variants):
    """Optional SAS enrichment pass — no-op if already populated."""
    return variants


# ── FILTER RARE VARIANTS ──────────────────────────────────────────────────────
def filter_rare_variants(variants, threshold=0.01):
    filtered = []
    for v in variants:
        af_data = v.get("gnomad_af")
        if af_data is None:
            filtered.append(v)
            continue
        af = (af_data.get("south_asian") or af_data.get("global")
              if isinstance(af_data, dict) else af_data)
        try:
            if float(af) <= threshold:
                filtered.append(v)
        except (TypeError, ValueError):
            filtered.append(v)
    return filtered


# ── SCORE AND RANK VARIANTS ───────────────────────────────────────────────────
def score_variant(v):
    score  = 0
    ann    = v.get("annotation", {})
    gene   = v.get("gene", "")
    impact = ann.get("impact", "").upper()
    score += {"HIGH":3, "MODERATE":2, "LOW":1}.get(impact, 0)
    try:
        if ann.get("sift") is not None and float(ann["sift"]) <= 0.05:
            score += 2
    except (TypeError, ValueError):
        pass
    try:
        if ann.get("polyphen") is not None and float(ann["polyphen"]) >= 0.85:
            score += 2
    except (TypeError, ValueError):
        pass
    cv = v.get("clinvar", "").lower()
    if "pathogenic" in cv and "likely" not in cv:
        score += 3
    elif "likely pathogenic" in cv:
        score += 2
    if gene in KNOWN_DISEASE_GENES:
        score += 1
    try:
        af_data = v.get("gnomad_af")
        af      = (af_data.get("south_asian") or af_data.get("global")
                   if isinstance(af_data, dict) else af_data)
        if af is not None and float(af) < 0.001:
            score += 1
    except (TypeError, ValueError):
        pass
    return min(score, 10)


def rank_variants(variants):
    for v in variants:
        s          = score_variant(v)
        v["score"] = s
        v["priority"] = "HIGH" if s >= 7 else "MEDIUM" if s >= 4 else "LOW"
    return sorted(variants, key=lambda x: x.get("score", 0), reverse=True)


# ── APPLY ACMG CLASSIFICATION ─────────────────────────────────────────────────
def apply_acmg_classification(variants):
    """Ensure every variant has ACMG table and evidence panel populated."""
    for v in variants:
        if not v.get("acmg_criteria_table"):
            ann = v.get("annotation", {})
            ct  = build_acmg_criteria_table(
                v, ann, v.get("gnomad_af"), v.get("clinvar", "Unknown")
            )
            cls, conf, ev = combine_acmg_from_table(ct)
            v.update({
                "acmg":                cls,
                "confidence_level":    conf,
                "acmg_evidence":       ev,
                "acmg_criteria_table": ct
            })
        if not v.get("evidence_panel"):
            ann = v.get("annotation", {})
            v["evidence_panel"] = build_evidence_panel(
                v, ann, v.get("gnomad_af"), v.get("clinvar", "Unknown")
            )
    return variants