#!/usr/bin/env python3
"""Phase 3: QC, human integration, marker-based annotation, and initial figures."""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy import sparse


MARKER_SETS = {
    "myonuclei_muscle_fiber": ["MYH1", "MYH2", "MYH7", "ACTA1", "DES", "TNNT3", "CKM"],
    "satellite_cell": ["PAX7", "MYF5", "MYOD1", "VCAM1"],
    "fap_stromal": ["PDGFRA", "DCN", "LUM", "COL1A1", "COL3A1", "PI16", "DPP4", "CXCL12"],
    "endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "CLDN5"],
    "smooth_muscle_pericyte": ["RGS5", "PDGFRB", "ACTA2", "TAGLN", "MYH11"],
    "macrophage_monocyte": ["LYZ", "CD68", "C1QA", "C1QB", "APOE", "SPP1", "FCGR3A"],
    "t_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "CD2"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "MZB1"],
    "schwann_cell": ["SOX10", "MPZ", "PLP1", "S100B"],
    "adipocyte": ["ADIPOQ", "PLIN1", "FABP4", "LPL"],
}


REFERENCE_MAP = [
    ("fap_stromal", ["fibro", "fap", "stromal", "mesenchymal", "tenocyte"]),
    ("myonuclei_muscle_fiber", ["myofiber", "myonuc", "muscle fiber", "skeletal muscle cell", "myocyte"]),
    ("satellite_cell", ["satellite", "muscle stem"]),
    ("endothelial", ["endothelial", "vascular endothelial"]),
    ("smooth_muscle_pericyte", ["smooth muscle", "pericyte", "smc"]),
    ("macrophage_monocyte", ["macrophage", "monocyte", "myeloid"]),
    ("t_cell", ["t cell", "nk cell", "lymphocyte"]),
    ("b_cell", ["b cell", "plasma"]),
    ("schwann_cell", ["schwann"]),
    ("adipocyte", ["adipocyte"]),
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase3_qc_integration.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def read_manifest(root: Path, config: dict) -> pd.DataFrame:
    path = root / config["phase2"]["standardized_manifest"]
    if not path.exists():
        raise RuntimeError(f"Missing standardized manifest: {path}")
    manifest = pd.read_csv(path)
    return manifest[manifest["status"] == "written"].copy()


def add_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    return adata


def filter_qc(adata: ad.AnnData, thresholds: dict) -> tuple[ad.AnnData, dict[str, int | float | str]]:
    before = adata.n_obs
    add_qc_metrics(adata)
    keep = (
        (adata.obs["n_genes_by_counts"] >= thresholds["min_genes"])
        & (adata.obs["total_counts"] >= thresholds["min_counts"])
        & (adata.obs["pct_counts_mt"] <= thresholds["max_pct_mito"])
    )
    filtered = adata[keep].copy()
    summary = {
        "dataset_id": str(adata.obs["dataset_id"].iloc[0]),
        "species": str(adata.obs["species"].iloc[0]),
        "cells_before_qc": int(before),
        "cells_after_qc": int(filtered.n_obs),
        "genes": int(filtered.n_vars),
        "median_genes": float(np.median(filtered.obs["n_genes_by_counts"])) if filtered.n_obs else 0,
        "median_counts": float(np.median(filtered.obs["total_counts"])) if filtered.n_obs else 0,
        "median_pct_mito": float(np.median(filtered.obs["pct_counts_mt"])) if filtered.n_obs else 0,
        "threshold_min_genes": thresholds["min_genes"],
        "threshold_min_counts": thresholds["min_counts"],
        "threshold_max_pct_mito": thresholds["max_pct_mito"],
    }
    return filtered, summary


def downsample_by_dataset(adata: ad.AnnData, max_cells: int, random_seed: int) -> ad.AnnData:
    if max_cells <= 0:
        return adata
    sampled_indices: list[str] = []
    rng = np.random.default_rng(random_seed)
    for dataset_id, idx in adata.obs.groupby("dataset_id", observed=True).indices.items():
        names = np.asarray(adata.obs_names)[idx]
        if len(names) > max_cells:
            names = rng.choice(names, size=max_cells, replace=False)
        sampled_indices.extend(names.tolist())
        logging.info("Dataset %s: kept %d cells for integration.", dataset_id, len(names))
    return adata[sampled_indices].copy()


def normalize_for_analysis(adata: ad.AnnData) -> ad.AnnData:
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def score_markers(adata: ad.AnnData) -> ad.AnnData:
    genes = set(adata.var_names)
    for cell_type, markers in MARKER_SETS.items():
        present = [gene for gene in markers if gene in genes]
        if len(present) >= 2:
            sc.tl.score_genes(adata, present, score_name=f"score_{cell_type}", use_raw=False)
        else:
            adata.obs[f"score_{cell_type}"] = 0.0
            logging.warning("Only %d markers available for %s.", len(present), cell_type)
    return adata


def map_reference_label(label: str) -> str | None:
    text = str(label).lower()
    if text in {"", "nan", "none", "unknown"}:
        return None
    for target, tokens in REFERENCE_MAP:
        if any(token in text for token in tokens):
            return target
    return None


def annotate_major_cell_types(adata: ad.AnnData, use_reference: bool) -> ad.AnnData:
    score_cols = [f"score_{cell_type}" for cell_type in MARKER_SETS]
    score_matrix = adata.obs[score_cols].to_numpy()
    best_idx = np.argmax(score_matrix, axis=1)
    best_scores = score_matrix[np.arange(score_matrix.shape[0]), best_idx]
    marker_labels = np.array(list(MARKER_SETS.keys()), dtype=object)[best_idx]
    marker_labels[best_scores <= 0] = "unassigned"

    labels = marker_labels.astype(object)
    if use_reference and "reference_cell_type" in adata.obs.columns:
        reference = adata.obs["reference_cell_type"].astype(str).copy()
        unknown_mask = reference.str.lower().isin(["", "nan", "none", "unknown"])
        for fallback_col in ["cell_annotation", "cell_annotation.1", "annotation", "CellType", "cell_type"]:
            if fallback_col in adata.obs.columns:
                fallback = adata.obs[fallback_col].astype(str)
                reference.loc[unknown_mask] = fallback.loc[unknown_mask]
                unknown_mask = reference.str.lower().isin(["", "nan", "none", "unknown"])
        ref_labels = reference.map(map_reference_label)
        ref_mask = ref_labels.notna().to_numpy()
        labels[ref_mask] = ref_labels[ref_mask].to_numpy()
    adata.obs["major_cell_type"] = pd.Categorical(labels)
    return adata


def run_integration(adata: ad.AnnData, config: dict) -> tuple[ad.AnnData, str]:
    n_top = int(config["phase3"]["n_top_hvg"])
    n_pcs = int(config["phase3"]["n_pcs"])
    batch_key = "dataset_id"

    logging.info("Selecting highly variable genes.")
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, batch_key=batch_key, flavor="seurat_v3")
    except Exception as exc:
        logging.warning("seurat_v3 HVG failed (%s); falling back to seurat flavor.", exc)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, batch_key=batch_key, flavor="seurat")

    n_hvg = int(adata.var["highly_variable"].sum())
    if n_hvg > 0:
        logging.info("Subsetting to %d highly variable genes for integration.", n_hvg)
        adata = adata[:, adata.var["highly_variable"]].copy()

    logging.info("Scaling and computing PCA.")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1), svd_solver="arpack")

    integration_rep = "X_pca"
    try:
        import harmonypy as hm

        if adata.obs[batch_key].nunique() > 1:
            logging.info("Running harmonypy integration by %s.", batch_key)
            harmony = hm.run_harmony(adata.obsm["X_pca"], adata.obs, [batch_key], verbose=True)
            corrected = harmony.Z_corr
            expected_shape = adata.obsm["X_pca"].shape
            if corrected.shape == expected_shape[::-1]:
                corrected = corrected.T
            if corrected.shape != expected_shape:
                raise RuntimeError(f"Harmony returned shape {corrected.shape}; expected {adata.obsm['X_pca'].shape}")
            adata.obsm["X_pca_harmony"] = corrected
            integration_rep = "X_pca_harmony"
    except Exception as exc:
        logging.warning("Harmony integration unavailable or failed: %s. Using PCA.", exc)

    sc.pp.neighbors(adata, n_pcs=min(n_pcs, adata.obsm[integration_rep].shape[1]), use_rep=integration_rep)
    sc.tl.umap(adata, random_state=int(config["phase3"]["random_seed"]))
    sc.tl.leiden(adata, resolution=0.8, key_added="leiden_phase3")
    return adata, integration_rep


def save_umap(adata: ad.AnnData, color: str, output_prefix: Path, title: str) -> None:
    sc.pl.umap(
        adata,
        color=color,
        frameon=False,
        size=5,
        title=title,
        show=False,
        legend_loc="right margin",
    )
    for suffix in [".png", ".pdf"]:
        plt.savefig(output_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_composition_plot(composition: pd.DataFrame, output_prefix: Path) -> None:
    plot_df = composition.copy()
    plot_df["fraction"] = plot_df["fraction"].astype(float)
    pivot = plot_df.pivot(index="dataset_id", columns="major_cell_type", values="fraction").fillna(0)
    ax = pivot.plot(kind="bar", stacked=True, figsize=(9, 4.8), width=0.8, colormap="tab20")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("")
    ax.legend(title="Major cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(output_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def summarize_marker_scores(adata: ad.AnnData) -> pd.DataFrame:
    rows = []
    for cell_type in adata.obs["major_cell_type"].cat.categories:
        subset = adata.obs[adata.obs["major_cell_type"] == cell_type]
        row = {"major_cell_type": cell_type, "n_cells": subset.shape[0]}
        for marker_type in MARKER_SETS:
            row[f"mean_score_{marker_type}"] = float(subset[f"score_{marker_type}"].mean()) if subset.shape[0] else 0
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    root = project_root()
    setup_logging(root)

    try:
        config = read_config(root)
        sc.settings.figdir = str(root / config["paths"]["figures"])
        manifest = read_manifest(root, config)
        species = config["phase3"]["selected_species_for_integration"]
        human_rows = manifest[manifest["species"] == species]
        if human_rows.empty:
            raise RuntimeError(f"No standardized objects found for species {species}.")

        objects: list[ad.AnnData] = []
        qc_rows: list[dict[str, int | float | str]] = []
        for row in human_rows.itertuples(index=False):
            path = root / row.output_path
            logging.info("Loading standardized object: %s", path)
            sample = sc.read_h5ad(path)
            sample.var_names_make_unique()
            filtered, summary = filter_qc(sample, config["phase3"]["qc_thresholds"])
            qc_rows.append(summary)
            if filtered.n_obs > 0:
                objects.append(filtered)

        if not objects:
            raise RuntimeError("All human objects were empty after QC.")

        logging.info("Concatenating %d human objects.", len(objects))
        adata = ad.concat(objects, join="inner", label="phase2_object", keys=[str(x.obs["dataset_id"].iloc[0]) for x in objects])
        adata.obs_names_make_unique()
        adata = downsample_by_dataset(
            adata,
            int(config["phase3"]["max_cells_per_dataset_for_integration"]),
            int(config["phase3"]["random_seed"]),
        )
        adata = normalize_for_analysis(adata)
        adata = score_markers(adata)
        adata = annotate_major_cell_types(
            adata,
            bool(config["phase3"]["annotation"]["use_reference_labels_when_available"]),
        )
        adata, integration_rep = run_integration(adata, config)
        adata.uns["phase3_integration_representation"] = integration_rep

        output_path = root / config["phase3"]["integrated_object"]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        logging.info("Writing integrated object: %s", output_path)
        adata.write_h5ad(output_path, compression="gzip")

        qc_table = pd.DataFrame(qc_rows)
        qc_table.to_csv(root / config["phase3"]["qc_summary_table"], index=False)

        comp = (
            adata.obs.groupby(["dataset_id", "major_cell_type"], observed=True)
            .size()
            .reset_index(name="n_cells")
        )
        totals = comp.groupby("dataset_id", observed=True)["n_cells"].transform("sum")
        comp["fraction"] = comp["n_cells"] / totals
        comp.to_csv(root / config["phase3"]["composition_table"], index=False)

        marker_summary = summarize_marker_scores(adata)
        marker_summary.to_csv(root / config["phase3"]["marker_score_table"], index=False)

        figure_dir = root / config["paths"]["figures"]
        figure_dir.mkdir(parents=True, exist_ok=True)
        save_umap(adata, "dataset_id", figure_dir / "phase3_umap_by_dataset", "Initial human atlas integration")
        save_umap(adata, "major_cell_type", figure_dir / "phase3_umap_by_major_cell_type", "Marker-supported major cell types")
        save_umap(adata, "leiden_phase3", figure_dir / "phase3_umap_by_leiden", "Phase 3 Leiden clusters")
        save_composition_plot(comp, figure_dir / "phase3_cell_type_composition_by_dataset")

        logging.info("Phase 3 QC/integration/annotation completed.")
        return 0
    except Exception as exc:
        logging.exception("Phase 3 failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
