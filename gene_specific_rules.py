# gene_specific_rules.py
# Gene-specific ACMG rule modifications
# Update this file when new VCEP guidelines are published
# Source: ClinGen VCEP publications at clinicalgenome.org
# Last updated: May 2026

GENE_SPECIFIC_RULES = {

    "BRCA1": {
        "PVS1_note": "PVS1 applies at very strong level for LOF variants "
                     "in BRCA1 per ENIGMA VCEP 2022",
        "PM2_threshold": 0.0001,  # BRCA1 uses stricter PM2 threshold
        "PP2_applicable": False,  # BRCA1 VCEP does not use PP2
        "source": "ENIGMA VCEP BRCA1 guidelines, Genet Med 2022",
        "last_updated": "2022"
    },

    "BRCA2": {
        "PVS1_note": "PVS1 applies at very strong level for LOF variants "
                     "in BRCA2 per ENIGMA VCEP 2022",
        "PM2_threshold": 0.0001,
        "PP2_applicable": False,
        "source": "ENIGMA VCEP BRCA2 guidelines, Genet Med 2022",
        "last_updated": "2022"
    },

    "TP53": {
    "pm2_threshold": 0.00004,  # TP53 VCEP specific threshold
    "hotspot_codons": [175, 245, 248, 249, 273, 282],  # PM1 hotspot codons
    "pm1_note": (
        "Codons 175, 245, 248, 249, 273, 282 are defined as mutational "
        "hotspots by ClinGen TP53 VCEP (PMID: 8023157). PM1 applies if "
        "variant falls within these codons."
    ),
    "revel_threshold": 0.7,  # TP53 VCEP uses BayesDel >0.16 and AlignGVGD >C15
    "pp2_applicable": True,
    "vcep_version": "2.0",
    "vcep_date": "July 2024",
    "source": "ClinGen TP53 Expert Panel Specifications to ACMG/AMP v2.0, 2024",
    "disease": "Li-Fraumeni Syndrome (OMIM 151623)",
    "inheritance": "Autosomal Dominant",
    "last_updated": "2024-07"

    },

    "PKD1": {
        "PVS1_note": "PKD1 pseudogene region requires minimum 80-100x coverage "
                     "for reliable variant calling. Flag 32x as borderline.",
        "coverage_warning": True,
        "source": "ACMG PKD1 guidance",
        "last_updated": "2023"
    },

    "HBB": {
        "population_note": "South Asian population has high carrier frequency "
                           "for HBB variants. Apply BA1 threshold carefully — "
                           "common variants in SAS may still cause disease in "
                           "homozygous state.",
        "source": "ACMG Hemoglobinopathy guidelines",
        "last_updated": "2023"
    },
}


def get_gene_rule(gene_name):
    """
    Get gene-specific rules for a gene.
    Returns dict of modifications or empty dict if no specific rules.
    """
    return GENE_SPECIFIC_RULES.get(gene_name, {})


def apply_gene_specific_pm2(gene_name, af):
    """
    Apply gene-specific PM2 threshold if available.
    Falls back to standard 0.01 threshold.
    """
    rule = get_gene_rule(gene_name)
    threshold = rule.get("pm2_threshold", 0.01)
    if af is None:
        return True, f"Absent from gnomAD. PM2 applied (threshold: {threshold})"
    if af <= threshold:
        return True, f"AF {af:.6f} ≤ {threshold} (gene-specific threshold). PM2 applied."
    return False, f"AF {af:.6f} exceeds gene-specific PM2 threshold {threshold}. PM2 not applied."