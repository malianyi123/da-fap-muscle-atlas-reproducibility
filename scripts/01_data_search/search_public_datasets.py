#!/usr/bin/env python3
"""Query NCBI GEO DataSets for configured Phase 1 search terms.

This helper is intentionally lightweight. It creates a reproducible search log
for broad discovery, while the curated inclusion table remains
data/metadata/manual_dataset_seed.csv.
"""

from __future__ import annotations

import csv
import logging
import sys
import time
from pathlib import Path
from urllib.parse import urlencode

import requests
import yaml


NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    config_path = root / "config/config.yaml"
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Missing config file: {config_path}") from exc
    except yaml.YAMLError as exc:
        raise RuntimeError(f"Could not parse {config_path}: {exc}") from exc


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase1_ncbi_geo_search.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def get_json(url: str, params: dict[str, str | int]) -> dict:
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    return response.json()


def search_term(term: str, retmax: int = 10) -> list[dict[str, str]]:
    search = get_json(
        NCBI_ESEARCH,
        {
            "db": "gds",
            "term": term,
            "retmax": retmax,
            "retmode": "json",
            "sort": "relevance",
        },
    )
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    summary = get_json(
        NCBI_ESUMMARY,
        {
            "db": "gds",
            "id": ",".join(ids),
            "retmode": "json",
        },
    )
    records = []
    for uid in ids:
        item = summary.get("result", {}).get(uid, {})
        accession = item.get("accession", "")
        title = item.get("title", "")
        gds_type = item.get("gdstype", "")
        taxon = item.get("taxon", "")
        records.append(
            {
                "search_term": term,
                "uid": uid,
                "accession": accession,
                "title": title,
                "gds_type": gds_type,
                "taxon": taxon,
                "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}" if accession else "",
            }
        )
    return records


def write_rows(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["search_term", "uid", "accession", "title", "gds_type", "taxon", "url"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote NCBI GEO search log: %s", output_path)


def main() -> int:
    root = project_root()
    setup_logging(root)

    try:
        config = read_config(root)
        terms = config["phase1"]["search_terms"]
        output_path = root / config["phase1"]["search_output"]
        rows: list[dict[str, str]] = []
        for term in terms:
            logging.info("Searching GEO DataSets: %s", term)
            rows.extend(search_term(term))
            time.sleep(0.35)
        write_rows(rows, output_path)
        logging.info("Completed %d search terms and %d result rows.", len(terms), len(rows))
        return 0
    except requests.HTTPError as exc:
        logging.exception("NCBI request failed: %s", exc)
        return 1
    except Exception as exc:
        logging.exception("Search failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

