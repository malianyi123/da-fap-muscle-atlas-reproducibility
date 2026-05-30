#!/usr/bin/env python3
"""Validate curated Phase 1 dataset records and build dataset_metadata.csv."""

from __future__ import annotations

import csv
import logging
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import yaml


REQUIRED_COLUMNS = [
    "dataset_id",
    "repository",
    "organism",
    "species",
    "tissue",
    "muscle_type",
    "disease_or_condition",
    "age_group",
    "sample_size",
    "number_of_cells_if_single_cell",
    "platform",
    "data_type",
    "raw_data_available",
    "processed_matrix_available",
    "clinical_metadata_available",
    "key_phenotypes",
    "publication_title",
    "publication_year",
    "DOI",
    "download_url",
    "inclusion_decision",
    "reason_for_decision",
]


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
    log_file = root / "results/logs/phase1_dataset_inventory.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def read_records(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            missing = [col for col in REQUIRED_COLUMNS if col not in (reader.fieldnames or [])]
            extra = [col for col in (reader.fieldnames or []) if col not in REQUIRED_COLUMNS]
            if missing:
                raise RuntimeError(f"Missing required columns in {path}: {missing}")
            if extra:
                logging.warning("Ignoring non-standard columns in %s: %s", path, extra)
            return [{col: row.get(col, "").strip() for col in REQUIRED_COLUMNS} for row in reader]
    except FileNotFoundError as exc:
        raise RuntimeError(f"Missing curated seed file: {path}") from exc


def validate_records(records: list[dict[str, str]]) -> None:
    if not records:
        raise RuntimeError("Dataset inventory is empty.")

    seen: set[str] = set()
    for i, row in enumerate(records, start=2):
        dataset_id = row["dataset_id"]
        if not dataset_id:
            raise RuntimeError(f"Row {i} has an empty dataset_id.")
        if dataset_id in seen:
            raise RuntimeError(f"Duplicate dataset_id: {dataset_id}")
        seen.add(dataset_id)
        if row["inclusion_decision"] not in {"include", "exclude", "maybe"}:
            raise RuntimeError(
                f"Row {i} ({dataset_id}) has invalid inclusion_decision: {row['inclusion_decision']}"
            )
        if not row["download_url"].startswith(("http://", "https://")):
            raise RuntimeError(f"Row {i} ({dataset_id}) has invalid download_url.")


def write_records(records: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
        writer.writeheader()
        writer.writerows(records)
    logging.info("Wrote validated inventory: %s", path)


def count_matching(records: list[dict[str, str]], include_only: bool, terms: list[str]) -> int:
    total = 0
    for row in records:
        if include_only and row["inclusion_decision"] != "include":
            continue
        data_type = row["data_type"].lower()
        if any(term in data_type for term in terms):
            total += 1
    return total


def write_summary(records: list[dict[str, str]], path: Path) -> None:
    include_records = [row for row in records if row["inclusion_decision"] == "include"]
    decision_counts = Counter(row["inclusion_decision"] for row in records)
    organism_counts = Counter(row["organism"] for row in records)
    data_type_counts = Counter(row["data_type"] for row in records)

    sc_or_sn = count_matching(records, True, ["scrna-seq", "snrna-seq"])
    bulk = count_matching(records, True, ["bulk"])
    spatial = count_matching(records, True, ["spatial"])

    top_includes = "\n".join(
        f"- {row['dataset_id']}: {row['data_type']} ({row['species']}); {row['reason_for_decision']}"
        for row in include_records
    )

    content = f"""# Phase 1 Dataset Inventory Summary

Generated: {datetime.now().isoformat(timespec="seconds")}

## Counts

- Total candidate records: {len(records)}
- Inclusion decisions: {dict(decision_counts)}
- Organisms: {dict(organism_counts)}
- Data types: {dict(data_type_counts)}

## Phase Gate

- Included scRNA-seq/snRNA-seq candidates: {sc_or_sn}
- Included bulk validation candidates: {bulk}
- Included spatial transcriptomics candidates: {spatial}

Minimum discovery gate status: {"PASS" if sc_or_sn >= 1 and bulk >= 1 else "FAIL"}

Spatial status: {"available for supportive analysis" if spatial >= 1 else "fallback required"}

## Included Datasets

{top_includes}

## Caution

Human spatial transcriptomics specifically for skeletal muscle aging/degeneration was not confirmed in this Phase 1 pass. Current spatial support is mouse dystrophic/regeneration centered. Treat spatial conclusions as cross-species mechanistic support unless a human spatial dataset is added later.
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    logging.info("Wrote summary: %s", path)


def main() -> int:
    root = project_root()
    setup_logging(root)

    try:
        config = read_config(root)
        seed_path = root / config["phase1"]["metadata_seed"]
        output_path = root / config["phase1"]["metadata_output"]
        summary_path = root / config["phase1"]["summary_output"]

        records = read_records(seed_path)
        validate_records(records)
        write_records(records, output_path)
        write_summary(records, summary_path)
        logging.info("Phase 1 dataset inventory build completed successfully.")
        return 0
    except Exception as exc:
        logging.exception("Phase 1 dataset inventory build failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

