#!/usr/bin/env python3
# guideline_monitor.py
# SpectralG — Automated Guideline Update Monitor
# Run monthly: python3 guideline_monitor.py
# Or set up as a cron job: 0 9 1 * * python3 /path/to/guideline_monitor.py
#
# Checks:
# 1. PubMed for new ACMG/ClinGen VCEP publications
# 2. Ensembl VEP current version
# 3. gnomAD current version
# 4. ClinVar last update date
# Sends summary report to terminal and optionally email

import requests
import json
import os
import sys
from datetime import date, datetime, timedelta
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(dotenv_path=Path.home() / ".env", override=True)

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
NCBI_API_KEY  = os.getenv("NCBI_API_KEY", "")
STATE_FILE    = Path(__file__).parent / "guideline_monitor_state.json"
REPORT_FILE   = Path(__file__).parent / "guideline_update_report.txt"

# How many days back to search for new publications
LOOKBACK_DAYS = 35  # slightly more than a month to avoid missing anything


# ── LOAD / SAVE STATE ─────────────────────────────────────────────────────────
def load_state():
    """Load previous check state to detect new items."""
    if STATE_FILE.exists():
        return json.loads(STATE_FILE.read_text())
    return {
        "last_run": None,
        "known_pubmed_ids": [],
        "last_ensembl_version": None,
        "last_gnomad_version": None,
    }


def save_state(state):
    STATE_FILE.write_text(json.dumps(state, indent=2))


# ── CHECK 1: NEW VCEP / ACMG PUBLICATIONS ON PUBMED ──────────────────────────
def check_pubmed_guidelines():
    """
    Search PubMed for new ACMG, ClinGen VCEP, and variant classification
    guideline publications in the last LOOKBACK_DAYS days.
    Returns list of new article dicts.
    """
    print("Checking PubMed for new guideline publications...")

    cutoff = (date.today() - timedelta(days=LOOKBACK_DAYS)).strftime("%Y/%m/%d")

    # Search terms cover ACMG guidelines, ClinGen VCEP, and variant classification
    searches = [
        (
            "ACMG variant classification guidelines",
            '("ACMG"[Title/Abstract] OR "ClinGen"[Title/Abstract]) AND '
            '("variant classification"[Title/Abstract] OR "VCEP"[Title/Abstract] OR '
            '"variant interpretation"[Title/Abstract]) AND '
            f'("{cutoff}"[Date - Publication] : "3000"[Date - Publication])'
        ),
        (
            "Genetics in Medicine ACMG guidelines",
            '("Genetics in Medicine"[Journal]) AND '
            '("variant"[Title/Abstract] OR "ACMG"[Title/Abstract]) AND '
            f'("{cutoff}"[Date - Publication] : "3000"[Date - Publication])'
        ),
        (
            "gnomAD population frequencies update",
            '"gnomAD"[Title/Abstract] AND '
            f'("{cutoff}"[Date - Publication] : "3000"[Date - Publication])'
        ),
    ]

    all_articles = []
    seen_ids = set()

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    for search_name, query in searches:
        try:
            # Search
            r = requests.get(f"{base}esearch.fcgi", params={
                "db":      "pubmed",
                "term":    query,
                "retmax":  10,
                "retmode": "json",
                "sort":    "pub+date",
                "api_key": NCBI_API_KEY
            }, timeout=15)

            if not r.ok:
                print(f"  PubMed search failed for '{search_name}': {r.status_code}")
                continue

            ids = r.json().get("esearchresult", {}).get("idlist", [])
            new_ids = [i for i in ids if i not in seen_ids]
            seen_ids.update(new_ids)

            if not new_ids:
                continue

            # Fetch summaries
            r2 = requests.get(f"{base}esummary.fcgi", params={
                "db":      "pubmed",
                "id":      ",".join(new_ids),
                "retmode": "json",
                "api_key": NCBI_API_KEY
            }, timeout=15)

            if not r2.ok:
                continue

            result = r2.json().get("result", {})
            for pmid in new_ids:
                art = result.get(pmid, {})
                title = art.get("title", "")
                if title:
                    all_articles.append({
                        "pmid":    pmid,
                        "title":   title[:200],
                        "journal": art.get("source", ""),
                        "date":    art.get("pubdate", ""),
                        "url":     f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        "search":  search_name,
                    })

        except Exception as e:
            print(f"  PubMed error for '{search_name}': {e}")

    return all_articles


# ── CHECK 2: ENSEMBL VEP CURRENT VERSION ──────────────────────────────────────
def check_ensembl_version():
    """
    Check current Ensembl REST API version.
    Returns version string or None.
    """
    print("Checking Ensembl VEP version...")
    try:
        r = requests.get(
            "https://rest.ensembl.org/info/software",
            headers={"Accept": "application/json"},
            timeout=10
        )
        if r.ok:
            data = r.json()
            version = data.get("release", data.get("version", "Unknown"))
            print(f"  Ensembl current version: {version}")
            return str(version)
    except Exception as e:
        print(f"  Ensembl version check failed: {e}")
    return None


# ── CHECK 3: GNOMAD CURRENT VERSION ───────────────────────────────────────────
def check_gnomad_version():
    """
    Check gnomAD current version via their API.
    Returns version string or None.
    """
    print("Checking gnomAD version...")
    try:
        # gnomAD GraphQL API
        query = '{ meta { version } }'
        r = requests.post(
            "https://gnomad.broadinstitute.org/api",
            json={"query": query},
            headers={"Content-Type": "application/json"},
            timeout=10
        )
        if r.ok:
            version = (r.json().get("data", {})
                               .get("meta", {})
                               .get("version", "Unknown"))
            print(f"  gnomAD current version: {version}")
            return str(version)
    except Exception as e:
        print(f"  gnomAD version check failed: {e}")
    return None


# ── CHECK 4: CLINVAR LAST UPDATED ─────────────────────────────────────────────
def check_clinvar_update():
    """
    Check ClinVar database info for last update date.
    Returns date string or None.
    """
    print("Checking ClinVar update status...")
    try:
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi",
            params={"db": "clinvar", "retmode": "json",
                    "api_key": NCBI_API_KEY},
            timeout=10
        )
        if r.ok:
            data = r.json()
            db_info = data.get("einforesult", {}).get("dbinfo", {})
            last_update = db_info.get("lastupdatedate", "Unknown")
            print(f"  ClinVar last updated: {last_update}")
            return last_update
    except Exception as e:
        print(f"  ClinVar check failed: {e}")
    return None


# ── CHECK 5: CLINGEN VCEP AFFILIATIONS ───────────────────────────────────────
def check_clingen_vcep():
    """
    Search PubMed specifically for ClinGen VCEP gene-specific guideline papers.
    Returns list of new VCEP publications.
    """
    print("Checking for new ClinGen VCEP gene-specific guidelines...")

    cutoff = (date.today() - timedelta(days=LOOKBACK_DAYS)).strftime("%Y/%m/%d")

    try:
        query = (
            '"ClinGen"[Title/Abstract] AND '
            '"expert panel"[Title/Abstract] AND '
            '"variant curation"[Title/Abstract] AND '
            f'("{cutoff}"[Date - Publication] : "3000"[Date - Publication])'
        )

        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        r = requests.get(f"{base}esearch.fcgi", params={
            "db":      "pubmed",
            "term":    query,
            "retmax":  5,
            "retmode": "json",
            "api_key": NCBI_API_KEY
        }, timeout=15)

        if not r.ok:
            return []

        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        r2 = requests.get(f"{base}esummary.fcgi", params={
            "db": "pubmed", "id": ",".join(ids),
            "retmode": "json", "api_key": NCBI_API_KEY
        }, timeout=15)

        if not r2.ok:
            return []

        result = r2.json().get("result", {})
        vcep_papers = []
        for pmid in ids:
            art = result.get(pmid, {})
            title = art.get("title", "")
            if title:
                vcep_papers.append({
                    "pmid":    pmid,
                    "title":   title[:200],
                    "journal": art.get("source", ""),
                    "date":    art.get("pubdate", ""),
                    "url":     f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                })

        return vcep_papers

    except Exception as e:
        print(f"  ClinGen VCEP check failed: {e}")
        return []


# ── GENERATE REPORT ───────────────────────────────────────────────────────────
def generate_report(state, articles, vcep_papers,
                    ensembl_ver, gnomad_ver, clinvar_date):
    """
    Generate a structured update report.
    Highlights new items vs previous check.
    """
    today = date.today().strftime("%d %B %Y")
    report_lines = [
        "=" * 70,
        f"SPECTRALG GUIDELINE MONITOR REPORT",
        f"Generated: {today}",
        f"Lookback period: {LOOKBACK_DAYS} days",
        "=" * 70,
        "",
    ]

    # ── Version changes ───────────────────────────────────────────────────────
    report_lines.append("DATABASE AND TOOL VERSIONS")
    report_lines.append("-" * 40)

    prev_ensembl = state.get("last_ensembl_version")
    if ensembl_ver:
        if prev_ensembl and ensembl_ver != prev_ensembl:
            report_lines.append(
                f"🔴 ENSEMBL VEP UPDATED: {prev_ensembl} → {ensembl_ver}"
            )
            report_lines.append(
                "   ACTION REQUIRED: Check Ensembl release notes for VEP changes"
            )
            report_lines.append(
                "   URL: https://ensembl.info/category/02-ensembl-vep"
            )
        else:
            report_lines.append(f"✅ Ensembl VEP: version {ensembl_ver} (no change)")

    prev_gnomad = state.get("last_gnomad_version")
    if gnomad_ver:
        if prev_gnomad and gnomad_ver != prev_gnomad:
            report_lines.append(
                f"🔴 GNOMAD UPDATED: {prev_gnomad} → {gnomad_ver}"
            )
            report_lines.append(
                "   ACTION REQUIRED: Check new population frequencies impact on AF thresholds"
            )
        else:
            report_lines.append(f"✅ gnomAD: version {gnomad_ver} (no change)")

    if clinvar_date:
        report_lines.append(f"ℹ️  ClinVar last updated: {clinvar_date}")

    report_lines.append("")

    # ── New VCEP papers ───────────────────────────────────────────────────────
    report_lines.append("NEW CLINGEN VCEP GENE-SPECIFIC GUIDELINES")
    report_lines.append("-" * 40)

    prev_ids = set(state.get("known_pubmed_ids", []))
    new_vcep = [p for p in vcep_papers if p["pmid"] not in prev_ids]

    if new_vcep:
        report_lines.append(
            f"🔴 {len(new_vcep)} NEW VCEP PAPER(S) FOUND — ACTION REQUIRED"
        )
        for p in new_vcep:
            report_lines.append(f"\n  Title: {p['title']}")
            report_lines.append(f"  Journal: {p['journal']} ({p['date']})")
            report_lines.append(f"  URL: {p['url']}")
            report_lines.append(
                "  ACTION: Read paper, identify affected gene, "
                "update gene_specific_rules.py"
            )
    else:
        report_lines.append("✅ No new ClinGen VCEP papers found this period")

    report_lines.append("")

    # ── New guideline publications ────────────────────────────────────────────
    report_lines.append("NEW ACMG / VARIANT CLASSIFICATION PUBLICATIONS")
    report_lines.append("-" * 40)

    new_articles = [a for a in articles if a["pmid"] not in prev_ids]

    if new_articles:
        report_lines.append(f"ℹ️  {len(new_articles)} new publication(s) found:")
        for a in new_articles[:10]:
            report_lines.append(f"\n  [{a['search']}]")
            report_lines.append(f"  Title: {a['title']}")
            report_lines.append(f"  Journal: {a['journal']} ({a['date']})")
            report_lines.append(f"  URL: {a['url']}")
    else:
        report_lines.append("✅ No new guideline publications found this period")

    report_lines.append("")

    # ── Action summary ────────────────────────────────────────────────────────
    report_lines.append("ACTION SUMMARY")
    report_lines.append("-" * 40)

    actions_needed = []
    if ensembl_ver and prev_ensembl and ensembl_ver != prev_ensembl:
        actions_needed.append(
            "Update Ensembl VEP version string in report_generator.py"
        )
    if gnomad_ver and prev_gnomad and gnomad_ver != prev_gnomad:
        actions_needed.append(
            "Update gnomAD version string in report_generator.py and pdf_generator.py"
        )
    if new_vcep:
        for p in new_vcep:
            actions_needed.append(
                f"Read VCEP paper and update gene_specific_rules.py: {p['title'][:80]}"
            )

    if actions_needed:
        report_lines.append("🔴 ACTIONS REQUIRED:")
        for i, action in enumerate(actions_needed, 1):
            report_lines.append(f"  {i}. {action}")
    else:
        report_lines.append(
            "✅ No immediate actions required — SpectralG is current"
        )

    report_lines.append("")
    report_lines.append(
        "SpectralG Guideline Monitor | Prakash NK | "
        f"Next run recommended: {(date.today() + timedelta(days=30)).strftime('%d %B %Y')}"
    )
    report_lines.append("=" * 70)

    return "\n".join(report_lines)


# ── MAIN ──────────────────────────────────────────────────────────────────────
def run_monitor():
    print("=" * 70)
    print("SpectralG Guideline Monitor")
    print(f"Date: {date.today().strftime('%d %B %Y')}")
    print("=" * 70)
    print()

    # Load previous state
    state = load_state()
    if state.get("last_run"):
        print(f"Last run: {state['last_run']}")
    else:
        print("First run — establishing baseline")
    print()

    # Run all checks
    articles    = check_pubmed_guidelines()
    vcep_papers = check_clingen_vcep()
    ensembl_ver = check_ensembl_version()
    gnomad_ver  = check_gnomad_version()
    clinvar_date = check_clinvar_update()

    print()

    # Generate report
    report = generate_report(
        state, articles, vcep_papers,
        ensembl_ver, gnomad_ver, clinvar_date
    )

    # Print to terminal
    print(report)

    # Save to file
    REPORT_FILE.write_text(report)
    print(f"\n✅ Report saved: {REPORT_FILE}")

    # Update state
    all_new_ids = (
        [a["pmid"] for a in articles] +
        [p["pmid"] for p in vcep_papers]
    )
    known_ids = list(set(state.get("known_pubmed_ids", []) + all_new_ids))

    new_state = {
        "last_run":             date.today().isoformat(),
        "known_pubmed_ids":     known_ids,
        "last_ensembl_version": ensembl_ver or state.get("last_ensembl_version"),
        "last_gnomad_version":  gnomad_ver  or state.get("last_gnomad_version"),
    }
    save_state(new_state)
    print("✅ State saved — next run will detect new items only")


if __name__ == "__main__":
    run_monitor()