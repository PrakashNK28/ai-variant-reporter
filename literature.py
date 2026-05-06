# literature.py
# Fetches recent PubMed literature for identified genes

import requests
import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(dotenv_path=Path.home() / ".env", override=True)


def fetch_pubmed_articles(gene_name, max_results=3):
    """
    Search PubMed for recent articles about a gene and pathogenicity.
    Returns list of article summaries.
    """
    try:
        api_key = os.getenv("NCBI_API_KEY", "")
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

        # Search for recent clinical articles about this gene
        search_term = (
            f"{gene_name}[gene] AND "
            f"(pathogenic[title/abstract] OR "
            f"variant[title/abstract] OR "
            f"mutation[title/abstract]) AND "
            f"humans[mesh]"
        )

        r = requests.get(f"{base}esearch.fcgi", params={
            "db": "pubmed",
            "term": search_term,
            "retmax": max_results,
            "retmode": "json",
            "sort": "relevance",
            "api_key": api_key
        }, timeout=10)

        if not r.ok:
            return []

        ids = r.json().get("esearchresult", {}).get("idlist", [])

        if not ids:
            return []

        # Fetch article summaries
        r2 = requests.get(f"{base}esummary.fcgi", params={
            "db": "pubmed",
            "id": ",".join(ids),
            "retmode": "json",
            "api_key": api_key
        }, timeout=10)

        if not r2.ok:
            return []

        result = r2.json().get("result", {})
        articles = []

        for pmid in ids:
            article = result.get(pmid, {})
            title = article.get("title", "")
            authors = article.get("authors", [])
            first_author = authors[0].get("name", "") if authors else ""
            pub_date = article.get("pubdate", "")
            source = article.get("source", "")

            if title:
                articles.append({
                    "pmid": pmid,
                    "title": title[:150],
                    "first_author": first_author,
                    "year": pub_date[:4],
                    "journal": source,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })

        return articles

    except Exception as e:
        print(f"PubMed error for {gene_name}: {e}")
        return []


def fetch_literature_for_variants(variants):
    """
    Fetch PubMed articles for all unique genes in variant list.
    Returns dict: gene_name → list of articles
    """
    gene_literature = {}

    unique_genes = list(set([
        v.get("gene", "Unknown") for v in variants
        if v.get("gene") not in ["Unknown", "Intergenic", None]
    ]))

    for gene in unique_genes:
        print(f"📚 Fetching PubMed articles for {gene}...")
        articles = fetch_pubmed_articles(gene, max_results=3)
        gene_literature[gene] = articles

    return gene_literature