#!/usr/bin/env python3
"""Create compact main-text candidate tables from source result tables."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
METADATA_PATH = PROJECT_ROOT / "data" / "metadata" / "dataset_metadata.csv"


CORE_REGULATORS = ["PPARG", "SMAD3", "STAT3", "CEBPA", "JUN", "FOS", "KLF4", "KLF5"]


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input not found: {path.relative_to(PROJECT_ROOT)}")


def role_from_row(row: pd.Series) -> str:
    data_type = str(row.get("data_type", "")).lower()
    dataset_id = str(row.get("dataset_id", ""))
    decision = str(row.get("inclusion_decision", "")).lower()
    if "spatial" in data_type:
        return "spatial validation"
    if "bulk" in data_type:
        return "bulk validation/meta-validation"
    if "scrna" in data_type or "snrna" in data_type or "single" in data_type:
        if decision == "include":
            return "single-cell discovery/integration"
        return "single-cell support or sensitivity"
    if dataset_id.startswith("GSE254300"):
        return "paraspinal extension context"
    return "supporting evidence"


def build_dataset_overview() -> Path:
    require_file(METADATA_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    selected = metadata[metadata["inclusion_decision"].isin(["include", "maybe"])].copy()
    selected["role_in_project"] = selected.apply(role_from_row, axis=1)
    columns = [
        "dataset_id",
        "repository",
        "species",
        "tissue",
        "muscle_type",
        "disease_or_condition",
        "age_group",
        "sample_size",
        "number_of_cells_if_single_cell",
        "platform",
        "data_type",
        "publication_year",
        "DOI",
        "inclusion_decision",
        "role_in_project",
        "reason_for_decision",
    ]
    table = selected[[col for col in columns if col in selected.columns]].sort_values(
        ["inclusion_decision", "data_type", "dataset_id"]
    )
    out = TABLES_DIR / "main_table_1_dataset_overview.csv"
    table.to_csv(out, index=False)
    return out


def build_virtual_knockdown_table() -> Path:
    path = TABLES_DIR / "phase8b_grn_virtual_knockdown_ranking.csv"
    require_file(path)
    ranking = pd.read_csv(path)
    table = ranking[ranking["regulator"].isin(CORE_REGULATORS)].copy()
    table["rank_among_core_candidates"] = table["perturbation_sensitivity_score"].rank(
        ascending=False, method="first"
    ).astype(int)
    table = table.sort_values("rank_among_core_candidates")
    columns = [
        "rank_among_core_candidates",
        "regulator",
        "predicted_delta_da_fap_program",
        "predicted_delta_adipogenic",
        "predicted_delta_fibrotic_ecm",
        "predicted_delta_inflammatory_remodeling",
        "perturbation_sensitivity_score",
        "interpretation",
    ]
    out = TABLES_DIR / "main_table_2_core_virtual_knockdown_regulators.csv"
    table[[col for col in columns if col in table.columns]].to_csv(out, index=False)
    return out


def build_bulk_meta_table() -> Path:
    path = TABLES_DIR / "phase10b_bulk_random_effects_meta_analysis.csv"
    require_file(path)
    meta = pd.read_csv(path)
    priority = {
        "DA_FAP": 0,
        "adipogenic": 1,
        "fibrotic_ecm": 2,
        "inflammatory_remodeling": 3,
        "SMAD3_TGFB_like": 4,
        "STAT3_IL6_like": 5,
        "PPARG_adipogenic_regulon_like": 6,
        "OXPHOS_mito": 7,
        "myogenesis": 8,
    }
    meta["_order"] = meta["module"].map(priority).fillna(99)
    table = meta.sort_values(["_order", "module"]).drop(columns="_order")
    out = TABLES_DIR / "main_table_3_bulk_meta_validation_summary.csv"
    table.to_csv(out, index=False)
    return out


def main() -> int:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    try:
        outputs = [
            build_dataset_overview(),
            build_virtual_knockdown_table(),
            build_bulk_meta_table(),
        ]
    except Exception as exc:
        print(f"Failed to build main tables: {exc}", file=sys.stderr)
        return 1

    print("Main tables written:")
    for path in outputs:
        print(f"- {path.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
