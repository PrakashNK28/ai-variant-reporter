# indian_variant_registry.py
# SpectralG — Indian Variant Evidence Registry
# Crowdsources Indian population evidence for VUS reclassification
#
# Every interpreted Indian patient case adds one observation.
# Consensus classification after ≥3 observations.
# PS4 applied when ≥5 observations with consistent classification.

import json
from datetime import date
from pathlib import Path

REGISTRY_FILE = Path(__file__).parent / "indian_variant_registry.json"


# ── LOAD / SAVE ───────────────────────────────────────────────────────────────

def load_registry():
    """Load existing registry or create empty one."""
    if REGISTRY_FILE.exists():
        try:
            return json.loads(REGISTRY_FILE.read_text())
        except Exception:
            pass
    return {
        "variants": [],
        "metadata": {
            "created":     str(date.today()),
            "description": "SpectralG Indian Variant Evidence Registry",
            "total_entries": 0
        }
    }


def save_registry(registry):
    """Save registry to disk."""
    try:
        REGISTRY_FILE.write_text(json.dumps(registry, indent=2))
    except Exception as e:
        print(f"Registry save error: {e}")


# ── ADD OBSERVATION ───────────────────────────────────────────────────────────

def add_variant_observation(
    gene,
    hgvsc,
    hgvsp,
    population,
    phenotype,
    expert_classification,
    spectralg_classification,
    sanger_confirmed=False,
    lab_id="anonymous",
    notes=""
):
    """
    Add a new variant observation to the Indian registry.
    Call this after delivering and verifying a report.

    Args:
        gene:                    Gene symbol e.g. "BRCA1"
        hgvsc:                   HGVS cDNA notation e.g. "c.5266dupC"
        hgvsp:                   HGVS protein notation e.g. "p.Gln1756fs"
        population:              e.g. "Indian-South Asian", "Indian-Telugu"
        phenotype:               Clinical indication e.g. "Hereditary breast cancer"
        expert_classification:   Final classification by clinician
        spectralg_classification: What SpectralG classified it as
        sanger_confirmed:        Whether Sanger confirmation was done
        lab_id:                  Anonymous lab identifier
        notes:                   Any additional notes
    """
    registry    = load_registry()
    variant_key = f"{gene}:{hgvsc}"
    existing    = None

    for v in registry["variants"]:
        if v["key"] == variant_key:
            existing = v
            break

    observation = {
        "population":              population,
        "phenotype":               phenotype,
        "expert_classification":   expert_classification,
        "spectralg_classification": spectralg_classification,
        "sanger_confirmed":        sanger_confirmed,
        "lab_id":                  lab_id,
        "date":                    str(date.today()),
        "notes":                   notes
    }

    if existing:
        existing["observations"].append(observation)
        existing["observation_count"] += 1
        existing["consensus"]          = calculate_consensus(existing["observations"])
        existing["last_updated"]        = str(date.today())
    else:
        registry["variants"].append({
            "key":               variant_key,
            "gene":              gene,
            "hgvsc":             hgvsc,
            "hgvsp":             hgvsp,
            "observation_count": 1,
            "consensus":         expert_classification,
            "last_updated":      str(date.today()),
            "observations":      [observation]
        })

    registry["metadata"]["total_entries"] = len(registry["variants"])
    save_registry(registry)
    total = next(
        (v["observation_count"] for v in registry["variants"]
         if v["key"] == variant_key), 1
    )
    print(f"✅ Registry updated: {gene} {hgvsc} — total observations: {total}")
    return total


# ── CONSENSUS CALCULATION ─────────────────────────────────────────────────────

def calculate_consensus(observations):
    """
    Calculate consensus classification from multiple observations.
    Requires minimum 3 observations for a meaningful consensus.
    Expert classification takes priority over SpectralG classification.
    """
    if len(observations) < 3:
        return f"Preliminary ({len(observations)} observation(s) — minimum 3 required)"

    classifications = [
        o["expert_classification"] for o in observations
        if o["expert_classification"] not in ("Unknown", "VUS", "Not classified", "")
    ]

    if not classifications:
        return "VUS — insufficient Indian population evidence"

    from collections import Counter
    counts      = Counter(classifications)
    most_common = counts.most_common(1)[0]
    total       = len(observations)
    agreement   = most_common[1] / total

    if agreement >= 0.8:
        return (f"{most_common[0]} "
                f"(Indian evidence: {most_common[1]}/{total} observations, "
                f"{agreement*100:.0f}% agreement)")
    elif agreement >= 0.6:
        return (f"Likely {most_common[0]} "
                f"(Indian evidence: moderate confidence, "
                f"{most_common[1]}/{total} observations)")
    else:
        return "VUS — conflicting Indian population evidence"


# ── LOOKUP ────────────────────────────────────────────────────────────────────

def lookup_indian_evidence(gene, hgvsc):
    """
    Check if a variant has Indian population evidence in the registry.
    Called during annotation to supplement gnomAD SAS data.

    Returns evidence dict or None.
    """
    if not gene or not hgvsc:
        return None

    registry    = load_registry()
    variant_key = f"{gene}:{hgvsc}"

    for v in registry["variants"]:
        if v["key"] == variant_key:
            return {
                "found":             True,
                "gene":              v["gene"],
                "hgvsc":             v["hgvsc"],
                "hgvsp":             v.get("hgvsp", ""),
                "observation_count": v["observation_count"],
                "consensus":         v["consensus"],
                "last_updated":      v.get("last_updated", ""),
                "observations":      v["observations"]
            }

    return None


# ── STATISTICS ────────────────────────────────────────────────────────────────

def get_registry_stats():
    """Print summary statistics of the Indian Variant Registry."""
    registry  = load_registry()
    variants  = registry["variants"]
    metadata  = registry["metadata"]

    print(f"\n{'='*55}")
    print(f"SpectralG Indian Variant Registry")
    print(f"{'='*55}")
    print(f"Registry created:        {metadata.get('created', 'Unknown')}")
    print(f"Unique variants:         {len(variants)}")
    print(f"Total observations:      "
          f"{sum(v['observation_count'] for v in variants)}")

    if variants:
        from collections import Counter
        genes = Counter(v["gene"] for v in variants)
        print(f"\nTop genes:")
        for gene, count in genes.most_common(5):
            print(f"  {gene}: {count} variant(s)")

        reclassified = [
            v for v in variants
            if "Indian evidence" in v.get("consensus", "")
        ]
        print(f"\nVariants with consensus: {len(reclassified)}")

        sanger_confirmed = sum(
            1 for v in variants
            for o in v["observations"]
            if o.get("sanger_confirmed")
        )
        print(f"Sanger-confirmed obs:    {sanger_confirmed}")
    else:
        print("\nRegistry is empty.")
        print("It grows automatically with each SpectralG report delivered.")
        print("Use add_variant_observation() to add your first entry.")

    print(f"{'='*55}\n")


# ── EXPORT ────────────────────────────────────────────────────────────────────

def export_registry_csv():
    """Export registry to CSV for sharing or analysis."""
    import csv
    registry    = load_registry()
    output_file = Path(__file__).parent / "indian_variant_registry_export.csv"

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "hgvsc", "hgvsp", "observation_count",
            "consensus", "last_updated"
        ])
        for v in registry["variants"]:
            writer.writerow([
                v["gene"], v["hgvsc"], v.get("hgvsp", ""),
                v["observation_count"], v["consensus"],
                v.get("last_updated", "")
            ])

    print(f"✅ Registry exported to: {output_file}")
    return str(output_file)