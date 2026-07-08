# annotator.py
# SpectralG — Clinical Variant Annotator v2.2
# Research-grade: full ACMG table, HGVS notation, evidence panel, gnomAD SAS
#
# Framework: ACMG/AMP 2015 (Richards et al., Genet Med 2015)
# PP5/BP6 NOT applied per ACMG 2023 (Biesecker & Harrison)
#
# v2.2 changes vs v2.1:
#   - GAP 1: filter_rare_variants() now returns (rare, common) — common variants
#             are classified as Benign/LB and included in report, never silently dropped
#   - GAP 2: BA1/BS1 applied correctly before any filtering decision
#   - GAP 3: SAS AF used for BA1 check independently of global AF
#   - GAP 4: REVEL score fetched from MyVariant.info — PP3 now REVEL-primary
#   - GAP 5: CADD phred score fetched from MyVariant.info
#   - GAP 6: MANE Select transcript preferred in VEP params and selection logic
#   - GAP 7: BS2 applied when gnomAD homozygote count >= 50 in any population
#   - GAP 8: BS4 wired to clinical_evidence sidebar input
#   - GAP 9: BP7 extended to highly conservative missense (REVEL < 0.3 + SIFT tol + PP ben)
#   - GAP 10: NM_ RefSeq transcript ID extracted and stored
#   - GAP 11: rsID from VCF column 3 surfaced in report
#   - GAP 12: gene_specific_rules PM1 hotspot logic retained + HL-VCEP BA1 note
#   - GAP 13: ClinVar lookup now uses chrpos38 + rsID strategy
#   - GAP 14: gnomAD homozygote count extracted and stored
#   - GAP 15: version field added to variant dict for audit trail

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

SPECTRALG_VERSION = "v2.2"


# ── SAFE SCALAR EXTRACTION ────────────────────────────────────────────────────
def _safe_float(value):
    """
    Safely convert a value to float.
    Handles: float, int, str, list (takes max), None.
    VEP and MyVariant.info can return frequency/score values as lists
    when multiple alleles or transcripts are present — this helper
    ensures we always get a single comparable scalar or None.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except (ValueError, TypeError):
            return None
    if isinstance(value, list):
        # Filter out None values, take max (most conservative for pathogenicity)
        scalars = [_safe_float(x) for x in value if x is not None]
        scalars = [x for x in scalars if x is not None]
        return max(scalars) if scalars else None
    return None

# ── CACHE ─────────────────────────────────────────────────────────────────────
vep_cache = {}
myvariant_cache = {}

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


# ── ZYGOSITY FROM GENOTYPE ────────────────────────────────────────────────────
def zygosity_from_genotype(genotype, chrom=""):
    """
    Convert raw VCF GT field into a clinically meaningful zygosity label.
    Homozygosity matters for ACMG/AMP: BS2 fires on many SAS homozygotes
    in gnomAD; PM3 note surfaces for recessive disease genes.
    """
    if not genotype or genotype in (".", "./.", ".|."):
        return "Not provided"
    sep = "/" if "/" in genotype else "|" if "|" in genotype else None
    if not sep:
        return "Not provided"
    alleles = genotype.split(sep)
    if len(alleles) == 1:
        return "Hemizygous"
    if len(alleles) != 2:
        return "Not provided"
    a1, a2 = alleles[0], alleles[1]
    if a1 == "0" and a2 == "0":
        return "Homozygous reference (not a variant)"
    if (a1 == "0") != (a2 == "0"):
        chrom_clean = str(chrom).replace("chr", "").upper()
        if chrom_clean in ("X", "Y"):
            return "Heterozygous (verify sex — chrX/Y het call)"
        return "Heterozygous"
    if a1 == a2 and a1 != "0":
        return "Homozygous"
    if a1 != a2 and a1 != "0" and a2 != "0":
        return "Compound heterozygous (two different ALT alleles, this position)"
    return "Not provided"


# ── MYVARIANT.INFO — REVEL + CADD (GAP 4 + 5) ────────────────────────────────
def fetch_myvariant_scores(chrom, pos, ref, alt):
    """
    Fetch REVEL and CADD phred scores from MyVariant.info.
    These are not reliably returned by Ensembl VEP REST API.

    REVEL threshold for PP3 (supporting): ≥ 0.5  (Pejaver et al. 2022)
    REVEL threshold for PP3 (moderate):   ≥ 0.7  (ClinGen SVI)
    CADD threshold widely cited:          ≥ 15   (Kircher et al. 2014)

    Returns dict: {revel, cadd_phred} — values may be None if unavailable.
    """
    cache_key = f"mv-{chrom}-{pos}-{ref}-{alt}"
    if cache_key in myvariant_cache:
        return myvariant_cache[cache_key]

    result = {"revel": None, "cadd_phred": None}
    try:
        # MyVariant.info hg38 HGVS format
        vid = f"chr{chrom}:g.{pos}{ref}>{alt}"
        url = f"https://myvariant.info/v1/variant/{vid}"
        r   = requests.get(url, params={
            "fields": "dbnsfp.revel.score,cadd.phred",
            "assembly": "hg38"
        }, timeout=10)
        if r.ok:
            data = r.json()
            # REVEL — use _safe_float throughout to handle list/dict/scalar
            dbnsfp     = data.get("dbnsfp", {})
            if isinstance(dbnsfp, list) and dbnsfp:
                dbnsfp = dbnsfp[0]  # take first if list of dicts
            revel_data = dbnsfp.get("revel", {}) if isinstance(dbnsfp, dict) else {}
            if isinstance(revel_data, dict):
                result["revel"] = _safe_float(revel_data.get("score"))
            elif isinstance(revel_data, list):
                scores = []
                for x in revel_data:
                    if isinstance(x, dict):
                        scores.append(_safe_float(x.get("score")))
                    else:
                        scores.append(_safe_float(x))
                scores = [s for s in scores if s is not None]
                result["revel"] = max(scores) if scores else None
            else:
                result["revel"] = _safe_float(revel_data)
            # CADD — same treatment
            cadd_data = data.get("cadd", {})
            if isinstance(cadd_data, list):
                cadd_data = cadd_data[0] if cadd_data else {}
            if isinstance(cadd_data, dict):
                result["cadd_phred"] = _safe_float(cadd_data.get("phred"))
            else:
                result["cadd_phred"] = _safe_float(cadd_data)
        print(f"MyVariant: REVEL={result['revel']}, CADD={result['cadd_phred']} "
              f"for chr{chrom}:{pos}")
    except Exception as e:
        print(f"MyVariant.info error for chr{chrom}:{pos}: {e}")

    myvariant_cache[cache_key] = result
    return result


# ── VEP HGVS (PRIMARY) — with MANE Select (GAP 6) ────────────────────────────
def call_vep_hgvs(chrom, pos, ref, alt):
    key = f"{chrom}-{pos}-{ref}-{alt}"
    if key in vep_cache:
        return vep_cache[key]
    try:
        hgvs    = f"{chrom}:g.{pos}{ref}>{alt}"
        url     = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs}"
        headers = {"Content-Type":"application/json","Accept":"application/json"}
        params  = {
            "canonical":      1,
            "sift":           1,
            "polyphen":       1,
            "numbers":        1,
            "hgvs":           1,
            "domains":        1,
            "check_existing": 1,
            "total_length":   1,
            "mane":           1,   # GAP 6: request MANE Select flag
            "refseq":         1,   # GAP 10: request RefSeq NM_ IDs
        }
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
def call_vep_region(chrom, pos, ref, alt, info_fields=None):
    info_fields = info_fields or {}
    svtype = info_fields.get("SVTYPE")
    end    = info_fields.get("END")
    key    = f"{chrom}-{pos}-{ref}-{alt}-{svtype}-{end}-region"
    if key in vep_cache:
        return vep_cache[key]
    try:
        if svtype and end:
            region = (f"{chrom} {pos} . {ref} {alt} . . "
                      f"SVTYPE={svtype};END={end}")
            print(f"[ VEP REGION ] SV detected (SVTYPE={svtype}, END={end})")
        else:
            region = f"{chrom} {pos} . {ref} {alt} . . ."
        url     = "https://rest.ensembl.org/vep/human/region"
        headers = {"Content-Type":"application/json","Accept":"application/json"}
        r = requests.post(url, headers=headers,
                          json={"variants":[region],"canonical":1,
                                "sift":1,"polyphen":1,"hgvs":1,
                                "check_existing":1,"mane":1,"refseq":1},
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


# ── EXTRACT HGVS NOTATION ─────────────────────────────────────────────────────
def extract_hgvs(transcript):
    hgvsc = transcript.get("hgvsc") or transcript.get("hgvs_c") or ""
    hgvsp = transcript.get("hgvsp") or transcript.get("hgvs_p") or ""
    if hgvsc and ":" in hgvsc:
        hgvsc = hgvsc.split(":")[-1]
    if hgvsp and ":" in hgvsp:
        hgvsp = hgvsp.split(":")[-1]
    return (hgvsc or "Not available (VEP missing)",
            hgvsp or "Not available (VEP missing)")


# ── EXTRACT REFSEQ TRANSCRIPT ID (GAP 10) ─────────────────────────────────────
def extract_refseq_id(transcript):
    """
    Extract the best RefSeq (NM_) transcript ID from a VEP transcript.
    Priority: MANE Select > canonical RefSeq > first RefSeq available.
    """
    # Direct refseq_transcript_ids field (returned when refseq=1 in params)
    refseq_ids = transcript.get("refseq_transcript_ids", [])
    if refseq_ids:
        # Prefer NM_ over NR_/XM_
        nm_ids = [x for x in refseq_ids if x.startswith("NM_")]
        if nm_ids:
            return nm_ids[0]
        return refseq_ids[0]
    # Fallback: try to extract from hgvsc prefix
    hgvsc_raw = transcript.get("hgvsc", "")
    if hgvsc_raw and ":" in hgvsc_raw:
        prefix = hgvsc_raw.split(":")[0]
        if prefix.startswith("NM_") or prefix.startswith("NR_"):
            return prefix
    return ""


# ── EXTRACT GNOMAD FREQUENCIES + HOMOZYGOTE COUNT (GAP 14) ───────────────────
def extract_gnomad_af(vep_data):
    """
    Extract gnomAD AF and SAS homozygote count from VEP colocated_variants.
    Prioritises South Asian (SAS) subpopulation.
    Returns dict or None.
    """
    try:
        for var in vep_data.get("colocated_variants", []):
            freqs = var.get("frequencies", {})
            if not freqs:
                continue
            for allele, fd in freqs.items():
                if not isinstance(fd, dict):
                    continue
                sas_raw = fd.get("gnomad_sas") or fd.get("gnomadg_sas") or fd.get("sas")
                glb_raw = fd.get("gnomad")     or fd.get("gnomadg")    or fd.get("af")
                # Use _safe_float to handle cases where VEP returns a list
                # (multi-allelic positions can produce list-valued frequencies)
                sas = _safe_float(sas_raw)
                glb = _safe_float(glb_raw)
                if sas is not None or glb is not None:
                    # Homozygote counts — field names vary by VEP version
                    hom_sas_raw = (fd.get("gnomad_sas_nhomalt")
                                   or fd.get("nhomalt_sas")
                                   or fd.get("sas_nhomalt")
                                   or 0)
                    hom_glb_raw = (fd.get("gnomad_nhomalt")
                                   or fd.get("nhomalt")
                                   or 0)
                    hom_sas = int(hom_sas_raw) if isinstance(hom_sas_raw, (int,float)) else 0
                    hom_glb = int(hom_glb_raw) if isinstance(hom_glb_raw, (int,float)) else 0
                    return {
                        "global":             glb,
                        "south_asian":        sas,
                        "display":            sas if sas is not None else glb,
                        "source":             "VEP colocated_variants",
                        "sas_homozygotes":    hom_sas,
                        "global_homozygotes": hom_glb,
                    }
    except Exception as e:
        print(f"gnomAD extraction error: {e}")
    return None


# ── EXTRACT CLINVAR FROM VEP ──────────────────────────────────────────────────
def extract_clinvar(vep_data):
    try:
        for var in vep_data.get("colocated_variants", []):
            if "clin_sig" in var:
                sigs = var["clin_sig"]
                return ", ".join(sigs) if isinstance(sigs, list) else str(sigs)
    except Exception as e:
        print(f"ClinVar extraction error: {e}")
    return "Unknown"


# ── DIRECT CLINVAR LOOKUP — chrpos38 + rsID strategy (GAP 13) ────────────────
def _fetch_clinvar_summary(clinvar_id, base, api_key):
    """Fetch and parse a single ClinVar esummary record."""
    try:
        r2 = requests.get(f"{base}esummary.fcgi", params={
            "db": "clinvar", "id": clinvar_id,
            "retmode": "json", "api_key": api_key
        }, timeout=10)
        if not r2.ok:
            return "Unknown"
        result = r2.json().get("result", {}).get(clinvar_id, {})
        if not result:
            return "Unknown"
        sig = result.get("clinical_significance", {})
        if isinstance(sig, dict):
            return sig.get("description", "Unknown")
        return str(sig) if sig else "Unknown"
    except Exception:
        return "Unknown"


def lookup_clinvar_direct(chrom, pos, ref, alt, rsid=None):
    """
    Direct ClinVar lookup via NCBI E-utilities.
    Strategy 1: rsID lookup (most reliable — resolves across assemblies).
    Strategy 2: GRCh38 coordinate lookup (chrpos38).
    Strategy 3: GRCh37 coordinate lookup (legacy fallback).
    """
    base    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    api_key = os.getenv("NCBI_API_KEY", "")

    # Strategy 1 — rsID
    if rsid and rsid not in (".", "", None):
        try:
            r = requests.get(f"{base}esearch.fcgi", params={
                "db": "clinvar", "term": f"{rsid}[RS]",
                "retmax": 3, "retmode": "json", "api_key": api_key
            }, timeout=10)
            if r.ok:
                ids = r.json().get("esearchresult", {}).get("idlist", [])
                if ids:
                    classification = _fetch_clinvar_summary(ids[0], base, api_key)
                    if classification != "Unknown":
                        print(f"✅ ClinVar via rsID {rsid}: {classification}")
                        return classification
        except Exception as e:
            print(f"ClinVar rsID lookup error: {e}")

    # Strategy 2 — GRCh38 coordinate (chrpos38)
    try:
        term = f"{chrom}[chr] AND {pos}[chrpos38] AND human[orgn]"
        r = requests.get(f"{base}esearch.fcgi", params={
            "db": "clinvar", "term": term,
            "retmax": 3, "retmode": "json", "api_key": api_key
        }, timeout=10)
        if r.ok:
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                classification = _fetch_clinvar_summary(ids[0], base, api_key)
                if classification != "Unknown":
                    print(f"✅ ClinVar via GRCh38 chr{chrom}:{pos}: {classification}")
                    return classification
    except Exception as e:
        print(f"ClinVar GRCh38 lookup error: {e}")

    # Strategy 3 — GRCh37 coordinate fallback
    try:
        term = f"{chrom}[chr] AND {pos}[chrpos37] AND human[orgn]"
        r = requests.get(f"{base}esearch.fcgi", params={
            "db": "clinvar", "term": term,
            "retmax": 3, "retmode": "json", "api_key": api_key
        }, timeout=10)
        if r.ok:
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                classification = _fetch_clinvar_summary(ids[0], base, api_key)
                if classification != "Unknown":
                    print(f"✅ ClinVar via GRCh37 chr{chrom}:{pos}: {classification}")
                    return classification
    except Exception as e:
        print(f"ClinVar GRCh37 lookup error: {e}")

    return "Unknown"


# ── BUILD FULL ACMG CRITERIA TABLE ────────────────────────────────────────────
def build_acmg_criteria_table(variant, annotation, gnomad_af, clinvar):
    """
    Build complete ACMG/AMP 2015 criteria table — 28 criteria always populated.
    PP5/BP6 always False per ACMG 2023.
    v2.2: REVEL-primary PP3, BS2 via homozygote count, BP7 extended,
          BA1 checks SAS independently, BS4 wired to clinical evidence.
    """
    gene        = annotation.get("gene", "Unknown")
    consequence = annotation.get("consequence", "unknown")
    sift        = annotation.get("sift")
    polyphen    = annotation.get("polyphen")
    revel       = variant.get("revel")
    cadd        = variant.get("cadd_phred")
    zygosity    = variant.get("zygosity", "Not provided")

    # ── Resolve AF — prefer South Asian ──────────────────────────────────────
    af     = None
    af_sas = None
    af_glb = None
    af_str = "Not available in gnomAD"
    hom_sas = 0

    if isinstance(gnomad_af, dict):
        af_sas  = gnomad_af.get("south_asian")
        af_glb  = gnomad_af.get("global")
        hom_sas = gnomad_af.get("sas_homozygotes", 0) or 0
        af      = af_sas if af_sas is not None else af_glb
        af_str  = (f"gnomAD SAS: {af_sas if af_sas is not None else 'N/A'} | "
                   f"gnomAD Global: {af_glb if af_glb is not None else 'N/A'}")
    elif gnomad_af is not None:
        af     = _safe_float(gnomad_af)
        af_str = f"gnomAD AF: {af:.6f}" if af is not None else "gnomAD AF: not parseable"

    # ── Computational tool strings ────────────────────────────────────────────
    sift_dam = sift is not None and _safe_float(sift) <= 0.05
    poly_dam = polyphen is not None and _safe_float(polyphen) >= 0.85
    sift_tol = sift is not None and _safe_float(sift) > 0.05
    poly_ben = polyphen is not None and _safe_float(polyphen) < 0.45

    sift_str  = (f"SIFT {sift:.3f} ({'deleterious' if sift_dam else 'tolerated'})"
                 if sift is not None else "SIFT N/A")
    poly_str  = (f"PolyPhen {polyphen:.3f} ({'damaging' if poly_dam else 'benign'})"
                 if polyphen is not None else "PolyPhen N/A")
    revel_dam = revel is not None and _safe__safe_float(revel) >= 0.5
    revel_str = (f"REVEL {revel:.3f} ({'≥0.5 damaging' if revel_dam else '<0.5 not damaging'})"
                 if revel is not None else "REVEL not available (manual lookup required)")
    cadd_str  = (f"CADD phred {cadd:.1f} ({'≥15 potentially deleterious' if cadd and _safe_float(cadd) >= 15 else '<15 likely tolerated'})"
                 if cadd is not None else "CADD not available")

    # ── PVS1 ─────────────────────────────────────────────────────────────────
    pvs1 = consequence in HIGH_IMPACT_CONSEQUENCES
    pvs1_ev = (
        f"Loss-of-function consequence '{consequence}' identified in {gene}. "
        f"PVS1 applied when LOF is established disease mechanism."
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

    # ── PS4 — Indian Variant Registry ────────────────────────────────────────
    indian_ev = variant.get("indian_registry")
    if indian_ev and indian_ev.get("observation_count", 0) >= 5:
        ps4    = True
        ps4_ev = (
            f"Variant observed in {indian_ev['observation_count']} Indian patients "
            f"via SpectralG Indian Variant Registry. Consensus: {indian_ev['consensus']}. "
            f"PS4 applied — requires expert review before clinical reporting."
        )
    else:
        ps4    = False
        ps4_ev = ("Case-control data not available computationally. "
                  "Not yet observed in SpectralG Indian Variant Registry. PS4 not applied.")

    # ── PM1 — Hotspot/ClinGen VCEP (GAP 12 retained) ─────────────────────────
    pm1    = False
    pm1_ev = ("Domain-level annotation requires ClinGen/UniProt data. "
              "Not determined from VEP alone. PM1 not applied.")
    try:
        from gene_specific_rules import get_gene_rule
        gene_rule      = get_gene_rule(gene)
        hotspot_codons = gene_rule.get("hotspot_codons", [])
        hgvsp_str      = variant.get("hgvsp","") or annotation.get("hgvsp","")
        codon_num      = None
        m = re.search(r'[A-Za-z]+(\d+)[A-Za-z*]', hgvsp_str)
        if m:
            codon_num = int(m.group(1))
        if hotspot_codons and codon_num and codon_num in hotspot_codons:
            pm1    = True
            pm1_ev = (
                f"Codon {codon_num} in {gene} is a defined mutational hotspot "
                f"per ClinGen VCEP. Hotspot codons: {hotspot_codons}. "
                f"PM1 applied per {gene_rule.get('source','ClinGen VCEP')}."
            )
        elif hotspot_codons and codon_num:
            pm1_ev = (f"Codon {codon_num} not in hotspot list {hotspot_codons} for {gene}. "
                      f"PM1 not applied.")
    except Exception:
        pass  # gene_specific_rules not available — pm1 stays False

    # ── PM2 — Absent/rare in population ──────────────────────────────────────
    # Uses SAS AF preferentially per SpectralG design
    if af is None:
        pm2    = True
        pm2_ev = f"Absent from gnomAD ({af_str}). PM2 applied."
    elif _safe_float(af) < 0.0001:
        pm2    = True
        pm2_ev = f"Extremely rare ({af_str}). AF {_safe_float(af):.6f} < 0.0001. PM2 applied."
    elif _safe_float(af) < 0.01:
        pm2    = True
        pm2_ev = f"Rare ({af_str}). AF {_safe_float(af):.5f} < 0.01. PM2 applied."
    else:
        pm2    = False
        pm2_ev = f"AF {_safe_float(af):.5f} ({af_str}) exceeds rare threshold. PM2 not applied."

    # ── PM3 — In trans / zygosity note ───────────────────────────────────────
    pm3 = False
    if zygosity == "Homozygous":
        pm3_ev = (
            "Variant is HOMOZYGOUS. PM3 requires a second, different pathogenic "
            "variant in trans — homozygosity alone does not satisfy PM3. However, "
            "homozygosity is clinically significant for autosomal recessive genes. "
            "Review family history for consanguinity. PM3 not applied."
        )
    else:
        pm3_ev = "Phase and family data not available. PM3 not applied."

    # ── PM4 — Protein length change ───────────────────────────────────────────
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

    # ── PP3 — REVEL-primary (GAP 4) ──────────────────────────────────────────
    # Per Pejaver et al. 2022 (AJHG): REVEL ≥ 0.5 = PP3 supporting
    # Fallback: if REVEL unavailable, require BOTH SIFT deleterious AND PolyPhen damaging
    if revel is not None:
        pp3 = _safe__safe_float(revel) >= 0.5
        pp3_ev = (
            f"{revel_str} | {sift_str} | {poly_str} | {cadd_str}. "
            f"PP3 {'applied — REVEL ≥ 0.5 (Pejaver et al. 2022)' if pp3 else 'not applied — REVEL < 0.5'}. "
            f"REVEL is the primary computational evidence tool per ClinGen SVI guidance."
        )
    else:
        # Fallback without REVEL: both SIFT AND PolyPhen must agree
        pp3 = sift_dam and poly_dam
        pp3_ev = (
            f"REVEL not available — using SIFT+PolyPhen as fallback. "
            f"{sift_str} | {poly_str} | {cadd_str}. "
            f"PP3 {'applied — both tools predict damaging' if pp3 else 'not applied — requires both SIFT deleterious AND PolyPhen damaging when REVEL unavailable'}. "
            f"Manual REVEL lookup recommended at dbnsfp.genetics.utah.edu."
        )

    # ── PP4 ──────────────────────────────────────────────────────────────────
    pp4    = False
    pp4_ev = "Clinical phenotype not provided. PP4 not applied."

    # ── PP5 — INTENTIONALLY NOT APPLIED ──────────────────────────────────────
    pp5    = False
    pp5_ev = ("PP5 INTENTIONALLY NOT APPLIED per ACMG 2023 (Biesecker & Harrison). "
              "External laboratory assertions excluded as independent evidence. "
              f"ClinVar data documented separately: {clinvar}.")

    # ── BA1 — Stand-alone Benign: SAS AF > 5% OR Global AF > 5% (GAP 2+3) ──
    # ClinGen HL-VCEP and most VCEPs: BA1 threshold = >5% in ANY population
    ba1_sas = af_sas is not None and _safe_float(af_sas) > 0.05
    ba1_glb = af_glb is not None and _safe_float(af_glb) > 0.05
    ba1     = ba1_sas or ba1_glb
    if ba1_sas:
        ba1_ev = (f"gnomAD SAS AF = {af_sas:.4f} ({af_sas*100:.2f}%) exceeds 5% BA1 threshold "
                  f"in South Asian population. BA1 applied — variant is common in Indian ancestry. "
                  f"Note: Global AF may be lower due to population stratification.")
    elif ba1_glb:
        ba1_ev = (f"gnomAD Global AF = {af_glb:.4f} ({af_glb*100:.2f}%) exceeds 5% BA1 threshold. "
                  f"BA1 applied — common variant, likely benign.")
    else:
        ba1_ev = f"AF ({af_str}) does not exceed 5% threshold in any population. BA1 not applied."

    # ── BS1 — Strong Benign: AF 1-5% ─────────────────────────────────────────
    bs1_sas = af_sas is not None and 0.01 < _safe_float(af_sas) <= 0.05
    bs1_glb = af_glb is not None and 0.01 < _safe_float(af_glb) <= 0.05
    bs1     = bs1_sas or bs1_glb
    if bs1_sas:
        bs1_ev = (f"gnomAD SAS AF = {af_sas:.4f} ({af_sas*100:.2f}%) in 1-5% range. "
                  f"BS1 applied — relatively common in South Asian population.")
    elif bs1_glb:
        bs1_ev = f"gnomAD Global AF = {af_glb:.4f} in 1-5% range. BS1 applied."
    else:
        bs1_ev = f"AF ({af_str}) does not trigger BS1 threshold (1-5%). Not applied."

    # ── BS2 — Observed in healthy adults / gnomAD homozygotes (GAP 7) ────────
    # gnomAD homozygotes in SAS population = variant present in unaffected individuals
    # Threshold: ≥ 50 SAS homozygotes = strong evidence variant is not fully penetrant
    if hom_sas >= 50:
        bs2    = True
        bs2_ev = (
            f"{hom_sas} South Asian homozygotes observed in gnomAD v4.1. "
            f"Variant present in healthy unaffected individuals — strong evidence "
            f"against fully penetrant autosomal recessive disease causing variants. "
            f"BS2 applied per ACMG/AMP 2015 criterion."
        )
    elif hom_sas > 0:
        bs2    = False
        bs2_ev = (
            f"{hom_sas} SAS homozygotes in gnomAD (below threshold of 50 for BS2). "
            f"Insufficient evidence for BS2 application. Not applied."
        )
    else:
        bs2    = False
        bs2_ev = ("No gnomAD homozygote count data available. "
                  "Healthy adult observation data not available. BS2 not applied.")

    # ── BS3 ──────────────────────────────────────────────────────────────────
    bs3    = False
    bs3_ev = "Functional study data not available. BS3 not applied."

    # ── BS4 — Lack of family segregation (GAP 8) ─────────────────────────────
    # Wired to clinical_evidence sidebar input
    clinical_ev = variant.get("clinical_evidence", {})
    family_input = clinical_ev.get("family_affected", "Unknown")
    if "Unaffected family members carry variant" in family_input:
        bs4    = True
        bs4_ev = ("Unaffected family members carry this variant — lack of segregation "
                  "with disease. BS4 applied per clinical evidence provided.")
    else:
        bs4    = False
        bs4_ev = "Family segregation data not provided. BS4 not applied."

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

    # ── BP4 — Both SIFT tolerated AND PolyPhen benign AND REVEL < 0.5 ────────
    if revel is not None:
        bp4 = sift_tol and poly_ben and _safe__safe_float(revel) < 0.5
        bp4_ev = (
            f"{sift_str} | {poly_str} | {revel_str}. "
            f"BP4 {'applied — all three tools predict benign/tolerated' if bp4 else 'not applied — one or more tools do not support benign'}."
        )
    else:
        bp4 = sift_tol and poly_ben
        bp4_ev = (
            f"{sift_str} | {poly_str}. REVEL not available. "
            f"BP4 {'applied — SIFT tolerated + PolyPhen benign (REVEL unavailable)' if bp4 else 'not applied'}."
        )

    # ── BP5 ──────────────────────────────────────────────────────────────────
    bp5    = False
    bp5_ev = "Alternate molecular diagnosis data not available. BP5 not applied."

    # ── BP6 — NOT APPLIED (ACMG 2023) ────────────────────────────────────────
    bp6    = False
    bp6_ev = "BP6 NOT APPLIED per ACMG 2023 — analogous to PP5. Excluded."

    # ── BP7 — Synonymous OR highly conservative missense (GAP 9) ─────────────
    is_synonymous = consequence == "synonymous_variant"
    # Conservative missense: SIFT tolerated + PolyPhen benign + REVEL < 0.3
    # Val→Ile (Grantham 29), Ala→Gly (60), Ser→Thr (58) are classic conservative pairs
    is_conservative_missense = (
        consequence == "missense_variant"
        and sift_tol
        and poly_ben
        and (revel is None or _safe__safe_float(revel) < 0.3)
        and not pp3  # only apply BP7 when PP3 is also not firing
    )
    bp7    = is_synonymous or is_conservative_missense
    if is_synonymous:
        bp7_ev = f"Synonymous variant — no amino acid change. BP7 applied."
    elif is_conservative_missense:
        bp7_ev = (
            f"Highly conservative missense: {sift_str}, {poly_str}, {revel_str}. "
            f"Multiple benign computational signals. BP7 applied per ClinGen SVI guidance."
        )
    else:
        bp7_ev = f"Consequence '{consequence}' does not qualify for BP7. Not applied."

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
          (pvs and pm >= 1 and pp >= 1) or
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
    Build structured evidence panel.
    v2.2: includes REVEL, CADD, homozygote count, transcript ID, rsID.
    """
    gene        = annotation.get("gene", "Unknown")
    consequence = annotation.get("consequence", "unknown")
    impact      = annotation.get("impact", "UNKNOWN")
    sift        = annotation.get("sift")
    polyphen    = annotation.get("polyphen")
    hgvsc       = variant.get("hgvsc", "Not available")
    hgvsp       = variant.get("hgvsp", "Not available")
    revel       = variant.get("revel")
    cadd        = variant.get("cadd_phred")
    rsid        = variant.get("id", "Not in VCF")
    refseq_id   = variant.get("refseq_id", "")
    zygosity    = variant.get("zygosity", "Not provided")

    if isinstance(gnomad_af, dict):
        sas     = gnomad_af.get("south_asian")
        glb     = gnomad_af.get("global")
        hom_sas = gnomad_af.get("sas_homozygotes", 0)
        sas_str = f"{sas:.6f} ({sas*100:.3f}%)" if sas is not None else "Not available"
        glb_str = f"{glb:.6f} ({glb*100:.3f}%)" if glb is not None else "Not available"
        hom_str = f"{hom_sas} SAS homozygotes in gnomAD" if hom_sas else "Homozygote count not available"
    else:
        sas_str = glb_str = "Not available in VEP response"
        hom_str = "Not available"

    sift_str  = (f"{sift:.3f} ({'Deleterious ≤0.05' if sift <= 0.05 else 'Tolerated >0.05'})"
                 if sift is not None else "Not available")
    poly_str  = (f"{polyphen:.3f} ({'Probably Damaging ≥0.85' if polyphen >= 0.85 else 'Not Damaging <0.85'})"
                 if polyphen is not None else "Not available")
    revel_str = (f"{revel:.3f} ({'≥0.5 damaging' if revel >= 0.5 else '<0.5 not damaging'})"
                 if revel is not None else "Not available — manual lookup recommended at dbnsfp.genetics.utah.edu")
    cadd_str  = (f"{cadd:.1f} ({'≥15 potentially deleterious' if cadd >= 15 else '<15 likely tolerated'})"
                 if cadd is not None else "Not available")

    transcript_display = ""
    if refseq_id:
        transcript_display = f"{refseq_id}:{hgvsc}"
    else:
        transcript_display = hgvsc

    indian_reg = variant.get("indian_registry")
    indian_str = (
        f"SpectralG Indian Variant Registry: "
        f"{indian_reg.get('observation_count',0)} observations. "
        f"Consensus: {indian_reg.get('consensus','No data')}. "
        f"PS4 evidence when ≥5 observations."
        if indian_reg else
        "Not yet observed in SpectralG Indian Variant Registry."
    )

    return {
        "Zygosity": (
            f"{zygosity}. "
            + (
                "CLINICALLY SIGNIFICANT: Homozygous variants in recessive genes "
                "affect both alleles. Review consanguinity history. BS2 may apply "
                "if gnomAD homozygote count is high."
                if zygosity == "Homozygous" else
                "Hemizygous — confirm patient sex and chromosome "
                "(typical for male chrX/Y or mitochondrial variants)."
                if zygosity == "Hemizygous" else
                "Standard heterozygous call."
                if "Heterozygous" in zygosity else
                "Zygosity not determined from available genotype data."
            )
        ),
        "Population Data": (
            f"gnomAD SAS: {sas_str} | gnomAD Global: {glb_str} | "
            f"{hom_str} | "
            f"South Asian frequency prioritised for Indian patients per SpectralG design."
        ),
        "Predicted Consequence": (
            f"Gene: {gene} | rsID: {rsid} | Consequence: {consequence} | "
            f"VEP Impact: {impact} | "
            f"HGVS c.: {transcript_display} | HGVS p.: {hgvsp}"
        ),
        "Computational Evidence": (
            f"SIFT: {sift_str} | PolyPhen-2: {poly_str} | "
            f"REVEL: {revel_str} | CADD phred: {cadd_str} | "
            "PP3 applied when REVEL ≥ 0.5 (Pejaver et al. 2022 / ClinGen SVI). "
            "BP4 applied when SIFT tolerated + PolyPhen benign + REVEL < 0.5."
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
    v2.2: adds REVEL/CADD fetch, MANE Select preference, RefSeq ID extraction,
          ClinVar rsID+chrpos38 lookup, gnomAD homozygote extraction.
    """
    chrom = str(variant.get("chrom",""))
    pos   = int(variant.get("pos", 0))
    ref   = str(variant.get("ref",""))
    alt   = str(variant.get("alt",""))
    rsid  = variant.get("id", "")  # GAP 11: rsID from VCF column 3

    # Tag with version for audit trail (GAP 15)
    variant["spectralg_version"] = SPECTRALG_VERSION

    # ── Fetch REVEL + CADD from MyVariant.info (GAP 4 + 5) ──────────────────
    mv = fetch_myvariant_scores(chrom, pos, ref, alt)
    variant["revel"]      = mv.get("revel")
    variant["cadd_phred"] = mv.get("cadd_phred")
    time.sleep(0.3)  # rate limit MyVariant.info

    # ── VEP calls ────────────────────────────────────────────────────────────
    data = call_vep_hgvs(chrom, pos, ref, alt)
    if not data:
        print(f"⚠️ HGVS failed → region fallback chr{chrom}:{pos}")
        data = call_vep_region(chrom, pos, ref, alt,
                               info_fields=variant.get("info_fields"))

    if not data:
        # VEP completely failed
        parser_gene = variant.get("gene")
        has_real_parser_gene = parser_gene and parser_gene != "Unknown"
        if has_real_parser_gene:
            print(f"⚠️ Both VEP calls failed — keeping parser gene '{parser_gene}'")
            gene_name = parser_gene
        else:
            print(f"⚠️ Both VEP failed → Ensembl Overlap chr{chrom}:{pos}")
            gene_name = get_gene_from_ensembl_overlap(chrom, pos)

        annotation = {
            "gene": gene_name, "consequence": "unknown",
            "impact": "UNKNOWN", "sift": None, "polyphen": None,
            "hgvsc": "Not available (VEP failed)",
            "hgvsp": "Not available (VEP failed)"
        }
        zygosity = zygosity_from_genotype(variant.get("genotype"), chrom)
        gnomad_af = None
        clinvar   = lookup_clinvar_direct(chrom, pos, ref, alt, rsid)
        variant.update({
            "annotation": annotation, "gene": gene_name,
            "hgvsc": "Not available (VEP failed)",
            "hgvsp": "Not available (VEP failed)",
            "refseq_id": "", "rsid": rsid,
            "gnomad_af": None, "clinvar": clinvar,
            "vep_status": "failed", "zygosity": zygosity
        })

    else:
        vep_data    = data[0]
        transcripts = vep_data.get("transcript_consequences", [])
        most_severe = vep_data.get("most_severe_consequence", "unknown")

        # ── SV handling (unchanged) ───────────────────────────────────────────
        info_fields = variant.get("info_fields", {})
        is_sv = bool(info_fields.get("SVTYPE"))

        if is_sv and transcripts:
            affected_genes = sorted({
                t.get("gene_symbol") for t in transcripts if t.get("gene_symbol")
            })
            variant["sv_affected_genes"] = affected_genes
            parser_gene = variant.get("gene")
            if parser_gene and parser_gene != "Unknown":
                gene_name = parser_gene
                matching = next((t for t in transcripts
                                 if t.get("gene_symbol") == parser_gene), None)
                if matching:
                    hgvsc, hgvsp = extract_hgvs(matching)
                    refseq_id    = extract_refseq_id(matching)
                    sift_score   = matching.get("sift_score")
                    polyphen_score = matching.get("polyphen_score")
                    chosen_consequence = matching.get("consequence_terms",[most_severe])[0]
                    chosen_impact      = matching.get("impact","UNKNOWN")
                else:
                    hgvsc = hgvsp = "Not available (gene not in VEP transcript list)"
                    refseq_id = ""
                    sift_score = polyphen_score = None
                    chosen_consequence = most_severe
                    chosen_impact = "UNKNOWN"

                annotation = {
                    "gene": gene_name, "consequence": chosen_consequence,
                    "impact": chosen_impact, "sift": sift_score,
                    "polyphen": polyphen_score, "hgvsc": hgvsc, "hgvsp": hgvsp,
                    "sv_affected_genes": affected_genes
                }
                gnomad_af = extract_gnomad_af(vep_data)
                clinvar   = extract_clinvar(vep_data)
                if clinvar == "Unknown":
                    clinvar = lookup_clinvar_direct(chrom, pos, ref, alt, rsid)
                zygosity = zygosity_from_genotype(variant.get("genotype"), chrom)
                variant.update({
                    "annotation": annotation, "gene": gene_name,
                    "hgvsc": hgvsc, "hgvsp": hgvsp,
                    "refseq_id": refseq_id, "rsid": rsid,
                    "gnomad_af": gnomad_af, "clinvar": clinvar,
                    "vep_status": "success", "zygosity": zygosity
                })
                ann = variant["annotation"]
                gaf = variant.get("gnomad_af")
                cv  = variant.get("clinvar","Unknown")
                ct  = build_acmg_criteria_table(variant, ann, gaf, cv)
                cls, conf, ev = combine_acmg_from_table(ct)
                variant.update({
                    "acmg": cls, "confidence_level": conf,
                    "acmg_evidence": ev, "acmg_criteria_table": ct,
                    "evidence_panel": build_evidence_panel(variant, ann, gaf, cv)
                })
                return variant

        # ── Standard point-mutation transcript picking ─────────────────────────
        gene_name = "Unknown"
        hgvsc = hgvsp = "Not available (VEP missing)"
        refseq_id = ""
        sift_score = polyphen_score = None
        chosen_consequence = most_severe
        chosen_impact = "UNKNOWN"

        DISTAL = {"upstream_gene_variant","downstream_gene_variant","intergenic_variant"}

        # Priority 0: MANE Select transcript (GAP 6)
        for t in transcripts:
            if t.get("mane_select") or t.get("mane") == "MANE Select":
                terms = set(t.get("consequence_terms",[]))
                if terms and terms.issubset(DISTAL):
                    continue
                gene_name          = t.get("gene_symbol","Unknown")
                hgvsc, hgvsp       = extract_hgvs(t)
                refseq_id          = extract_refseq_id(t)
                sift_score         = t.get("sift_score")
                polyphen_score     = t.get("polyphen_score")
                chosen_consequence = t.get("consequence_terms",[most_severe])[0]
                chosen_impact      = t.get("impact","UNKNOWN")
                print(f"✅ MANE Select transcript used for {gene_name}")
                break

        # Priority 1: canonical transcript (if MANE not found)
        if gene_name == "Unknown":
            for t in transcripts:
                if t.get("canonical") == 1:
                    terms = set(t.get("consequence_terms",[]))
                    if terms and terms.issubset(DISTAL):
                        continue
                    gene_name          = t.get("gene_symbol","Unknown")
                    hgvsc, hgvsp       = extract_hgvs(t)
                    refseq_id          = extract_refseq_id(t)
                    sift_score         = t.get("sift_score")
                    polyphen_score     = t.get("polyphen_score")
                    chosen_consequence = t.get("consequence_terms",[most_severe])[0]
                    chosen_impact      = t.get("impact","UNKNOWN")
                    break

        # Priority 2: first within-gene transcript
        if gene_name == "Unknown" and transcripts:
            within = next((t for t in transcripts
                           if not set(t.get("consequence_terms",[])).issubset(DISTAL)), None)
            t = within or transcripts[0]
            gene_name          = t.get("gene_symbol","Unknown")
            hgvsc, hgvsp       = extract_hgvs(t)
            refseq_id          = extract_refseq_id(t)
            sift_score         = t.get("sift_score")
            polyphen_score     = t.get("polyphen_score")
            chosen_consequence = t.get("consequence_terms",[most_severe])[0]
            chosen_impact      = t.get("impact","UNKNOWN")

        # Priority 3: Overlap or parser gene fallback
        if gene_name in ("Unknown","",None):
            parser_gene = variant.get("gene")
            if parser_gene and parser_gene not in ("Unknown","",None):
                print(f"⚠️ VEP no gene — keeping parser gene '{parser_gene}'")
                gene_name = parser_gene
            else:
                print(f"⚠️ VEP no gene → Overlap chr{chrom}:{pos}")
                gene_name = get_gene_from_ensembl_overlap(chrom, pos)

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
            "hgvsc": hgvsc, "hgvsp": hgvsp,
        }
        gnomad_af = extract_gnomad_af(vep_data)
        clinvar   = extract_clinvar(vep_data)

        # GAP 13: direct ClinVar fallback if VEP returned nothing
        if clinvar == "Unknown":
            print(f"ClinVar not in VEP → direct NCBI lookup chr{chrom}:{pos} rsID={rsid}")
            clinvar = lookup_clinvar_direct(chrom, pos, ref, alt, rsid)

        zygosity = zygosity_from_genotype(variant.get("genotype"), chrom)
        variant.update({
            "annotation": annotation, "gene": gene_name,
            "hgvsc": hgvsc, "hgvsp": hgvsp,
            "refseq_id": refseq_id, "rsid": rsid,
            "gnomad_af": gnomad_af, "clinvar": clinvar,
            "vep_status": "success", "zygosity": zygosity
        })

    # ── Indian Variant Registry (step 6) ──────────────────────────────────────
    current_gene  = variant.get("gene","Unknown")
    current_hgvsc = variant.get("hgvsc","")
    if REGISTRY_AVAILABLE and current_gene not in ("Unknown",""):
        try:
            indian_ev = lookup_indian_evidence(current_gene, current_hgvsc)
            variant["indian_registry"] = (
                indian_ev if indian_ev and indian_ev.get("observation_count",0) >= 3
                else None
            )
        except Exception as e:
            print(f"Registry lookup error: {e}")
            variant["indian_registry"] = None
    else:
        variant["indian_registry"] = None

    # ── ACMG table + evidence panel ───────────────────────────────────────────
    ann = variant["annotation"]
    gaf = variant.get("gnomad_af")
    cv  = variant.get("clinvar","Unknown")

    ct = build_acmg_criteria_table(variant, ann, gaf, cv)
    cls, conf, ev = combine_acmg_from_table(ct)
    variant.update({
        "acmg":                cls,
        "confidence_level":    conf,
        "acmg_evidence":       ev,
        "acmg_criteria_table": ct,
        "evidence_panel":      build_evidence_panel(variant, ann, gaf, cv)
    })
    return variant


# ── BULK ANNOTATION ───────────────────────────────────────────────────────────
def annotate_all(variants):
    """Annotate list with rate limiting — 1.1s VEP + 0.3s MyVariant."""
    annotated = []
    for i, v in enumerate(variants):
        print(f"Annotating {i+1}/{len(variants)}: chr{v.get('chrom')}:{v.get('pos')}")
        annotated.append(annotate_variant(v))
        time.sleep(1.1)
    return annotated


def enrich_gnomad_sas(variants):
    """No-op — SAS now extracted in annotate_variant."""
    return variants


# ── FILTER — returns (rare, common) NEVER silently drops variants (GAP 1) ────
def filter_rare_variants(variants, threshold=0.01):
    """
    Split variants into rare (clinically significant) and common (Benign/LB).

    CRITICAL CHANGE v2.2:
    Common variants (AF > threshold in SAS or global) are NO LONGER dropped.
    They are tagged with filter_reason and returned in the 'common' list.
    Both lists are returned so the caller can include common variants in the
    report as Benign/Likely Benign rather than silently omitting them.

    A blank report for a case with only common variants is a clinical
    patient safety error — it implies no variants were found, when in fact
    variants were found and correctly classified as Benign.

    Returns:
        (rare_variants, common_variants) — two lists
    """
    rare   = []
    common = []
    for v in variants:
        af_data = v.get("gnomad_af")
        if af_data is None:
            rare.append(v)
            continue

        # Use SAS AF preferentially (SpectralG design principle)
        af_sas = af_data.get("south_asian") if isinstance(af_data, dict) else None
        af_glb = af_data.get("global")      if isinstance(af_data, dict) else None
        af     = af_sas if af_sas is not None else af_glb

        if af is None:
            try:
                af = _safe_float(af_data)
            except (TypeError, ValueError):
                rare.append(v)
                continue

        try:
            if _safe_float(af) > threshold:
                v["filtered_common"] = True
                v["filter_reason"]   = (
                    f"Common variant — AF {_safe_float(af):.4f} "
                    f"({'SAS' if af_sas is not None else 'Global'}) "
                    f"> {threshold} threshold. "
                    f"BA1 applies if SAS AF > 5%. "
                    f"Classified and reported as Benign/Likely Benign."
                )
                common.append(v)
            else:
                rare.append(v)
        except (TypeError, ValueError):
            rare.append(v)

    return rare, common


def score_variant(v):
    """Score and prioritise variant for clinical ranking."""
    score  = 0
    ann    = v.get("annotation", {})
    gene   = v.get("gene","")
    impact = ann.get("impact","").upper()
    score += {"HIGH":3,"MODERATE":2,"LOW":1}.get(impact, 0)
    try:
        if ann.get("sift") is not None and _safe_float(ann["sift"]) <= 0.05:
            score += 2
    except (TypeError,ValueError): pass
    try:
        if ann.get("polyphen") is not None and _safe_float(ann["polyphen"]) >= 0.85:
            score += 2
    except (TypeError,ValueError): pass
    # REVEL bonus
    revel = v.get("revel")
    try:
        if revel is not None and _safe__safe_float(revel) >= 0.7:
            score += 2
        elif revel is not None and _safe__safe_float(revel) >= 0.5:
            score += 1
    except (TypeError,ValueError): pass
    cv = v.get("clinvar","").lower()
    if "pathogenic" in cv and "likely" not in cv: score += 3
    elif "likely pathogenic" in cv:               score += 2
    if gene in KNOWN_DISEASE_GENES:               score += 1
    try:
        af_data = v.get("gnomad_af")
        af = (af_data.get("south_asian") or af_data.get("global")
              if isinstance(af_data,dict) else af_data)
        if af is not None and _safe_float(af) < 0.001: score += 1
    except (TypeError,ValueError): pass
    # Criteria count boost
    criteria_count = len(v.get("acmg_evidence", []))
    if criteria_count >= 4:   score += 2
    elif criteria_count >= 2: score += 1
    return min(score, 10)


def rank_variants(variants):
    for v in variants:
        s          = score_variant(v)
        v["score"] = s
        v["priority"] = "HIGH" if s >= 7 else "MEDIUM" if s >= 4 else "LOW"
    return sorted(variants, key=lambda x: x.get("score",0), reverse=True)


def apply_acmg_classification(variants):
    """Ensure every variant has ACMG table and evidence panel."""
    for v in variants:
        if not v.get("acmg_criteria_table"):
            ann = v.get("annotation",{})
            ct  = build_acmg_criteria_table(
                v, ann, v.get("gnomad_af"), v.get("clinvar","Unknown")
            )
            cls, conf, ev = combine_acmg_from_table(ct)
            v.update({
                "acmg":                cls,
                "confidence_level":    conf,
                "acmg_evidence":       ev,
                "acmg_criteria_table": ct
            })
        if not v.get("evidence_panel"):
            ann = v.get("annotation",{})
            v["evidence_panel"] = build_evidence_panel(
                v, ann, v.get("gnomad_af"), v.get("clinvar","Unknown")
            )
    return variants
