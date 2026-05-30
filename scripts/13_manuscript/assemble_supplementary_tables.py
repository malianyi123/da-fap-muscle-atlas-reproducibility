#!/usr/bin/env python3
"""Create stable manuscript-facing supplementary table files."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import yaml


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def copy_table(root: Path, source: str, target: str) -> None:
    src = root / source
    dst = root / "results/tables" / target
    if not src.exists():
        raise RuntimeError(f"Missing required table: {src}")
    shutil.copyfile(src, dst)


def write_module_table(root: Path) -> None:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    rows = []
    for module, genes in config["gene_modules"].items():
        rows.append({"module": module, "genes": ";".join(genes), "n_genes": len(genes)})
    pd.DataFrame(rows).to_csv(root / "results/tables/supplementary_table_5_da_fap_gene_modules.csv", index=False)


def main() -> int:
    root = project_root()
    (root / "results/tables").mkdir(parents=True, exist_ok=True)
    copy_table(root, "data/metadata/dataset_metadata.csv", "supplementary_table_1_dataset_metadata.csv")
    copy_table(root, "results/tables/phase3_qc_summary.csv", "supplementary_table_2_qc_thresholds_retained_cells.csv")
    copy_table(root, "results/tables/phase3_marker_score_summary.csv", "supplementary_table_3_major_cell_type_marker_scores.csv")
    copy_table(root, "results/tables/phase4_fap_subtype_markers.csv", "supplementary_table_4_fap_subtype_markers.csv")
    write_module_table(root)
    copy_table(root, "results/tables/phase6_ligand_receptor_scores.csv", "supplementary_table_6_cell_cell_communication_results.csv")
    copy_table(root, "results/tables/phase7_regulator_correlations.csv", "supplementary_table_7_regulator_activity_proxy_results.csv")
    copy_table(root, "results/tables/phase8b_grn_virtual_knockdown_ranking.csv", "supplementary_table_8_grn_virtual_knockdown_ranking.csv")
    copy_table(root, "results/tables/phase9_spatial_colocalization_stats.csv", "supplementary_table_9_spatial_colocalization_statistics.csv")
    copy_table(root, "results/tables/phase10b_bulk_multi_cohort_effect_sizes.csv", "supplementary_table_10_bulk_validation_effect_sizes.csv")
    copy_table(root, "results/tables/phase4_fap_subtype_module_means.csv", "supplementary_table_11_subtype_module_programs.csv")
    copy_table(root, "results/tables/phase15_leave_one_dataset_out_sensitivity.csv", "supplementary_table_12_leave_one_dataset_out_sensitivity.csv")
    copy_table(root, "results/tables/phase15_module_gene_sensitivity.csv", "supplementary_table_13_module_gene_sensitivity.csv")
    copy_table(root, "results/tables/phase15_artifact_correlation_checks.csv", "supplementary_table_14_artifact_correlation_checks.csv")
    copy_table(root, "results/tables/phase15_integration_method_sensitivity.csv", "supplementary_table_15_integration_method_sensitivity.csv")
    copy_table(root, "results/tables/phase10b_bulk_random_effects_meta_analysis.csv", "supplementary_table_16_bulk_random_effects_meta_analysis.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
