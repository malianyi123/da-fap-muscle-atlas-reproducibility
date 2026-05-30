#!/usr/bin/env python3
"""Phase 9 spatial signature support using GSE225766 mouse Visium data.

This is a cross-species spatial-support workflow. It scores ortholog-like mouse
versions of the human DA-FAP modules in Visium spots, summarizes condition-level
patterns, and estimates spot-level co-localization between DA-FAP and ECM,
adipogenic, inflammatory, macrophage, and muscle atrophy signatures.
"""

from __future__ import annotations

import logging
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy.stats import kruskal, mannwhitneyu, spearmanr


GEO_SUPPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225766/suppl/GSE225766_RAW.tar"

MOUSE_MODULES = {
    "DA_FAP": [
        "Pparg", "Cebpa", "Fabp4", "Lpl", "Cd36",
        "Col1a1", "Col3a1", "Col6a1", "Fn1", "Postn", "Tnc", "Thbs2", "Lox", "Acta2",
        "Il6", "Ccl2", "Cxcl12", "Cxcl14", "Icam1", "Vcam1", "Spp1", "Mmp2", "Mmp14", "Timp1",
    ],
    "adipogenic": ["Pparg", "Cebpa", "Fabp4", "Lpl", "Adipoq", "Plin1", "Cd36"],
    "fibrotic_ecm": ["Col1a1", "Col3a1", "Col6a1", "Fn1", "Postn", "Tnc", "Thbs2", "Lox", "Acta2"],
    "inflammatory_remodeling": ["Il6", "Ccl2", "Cxcl12", "Cxcl14", "Icam1", "Vcam1", "Spp1", "Mmp2", "Mmp14", "Timp1"],
    "fap_stromal": ["Pdgfra", "Dcn", "Lum", "Col1a1", "Col3a1", "Pi16", "Dpp4", "Cxcl12"],
    "macrophage": ["Lyz2", "Cd68", "Adgre1", "C1qa", "C1qb", "Apoe", "Spp1"],
    "atrophic_muscle": ["Fbxo32", "Trim63", "Foxo3", "Mstn", "Casp3"],
    "myofiber": ["Acta1", "Des", "Myh1", "Myh2", "Myh4", "Myh7", "Tnnt3"],
}


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase9_spatial_validation.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def ensure_data(raw_dir: Path) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    h5_files = list(raw_dir.glob("*_filtered_feature_bc_matrix.h5"))
    if h5_files:
        logging.info("Found %d extracted Visium matrices.", len(h5_files))
        return
    tar_path = raw_dir / "GSE225766_RAW.tar"
    if not tar_path.exists():
        logging.info("Downloading %s", GEO_SUPPL_URL)
        try:
            urllib.request.urlretrieve(GEO_SUPPL_URL, tar_path)
        except Exception as exc:
            logging.warning("urllib download failed (%s); retrying with curl.", exc)
            subprocess.run(["curl", "-L", "--fail", "-o", str(tar_path), GEO_SUPPL_URL], check=True)
    logging.info("Extracting %s", tar_path)
    with tarfile.open(tar_path, "r") as tar:
        tar.extractall(raw_dir)


def parse_sample(stem: str) -> dict[str, str | int]:
    # Example: GSM7055903_MDX1_filtered_feature_bc_matrix
    sample_token = stem.replace("_filtered_feature_bc_matrix", "")
    gsm, label = sample_token.split("_", 1)
    if label.startswith("WT"):
        condition = "WT"
        condition_group = "WT"
        day = np.nan
    elif label.startswith("MDX"):
        condition = "MDX"
        condition_group = "MDX"
        day = np.nan
    elif "Day" in label and "PostInjury" in label:
        day_text = label.replace("Day", "").replace("_PostInjury", "")
        day = int(day_text)
        condition = f"injury_day{day}"
        condition_group = "injury"
    else:
        condition = label
        condition_group = "other"
        day = np.nan
    return {
        "sample_id": sample_token,
        "gsm": gsm,
        "sample_label": label,
        "condition": condition,
        "condition_group": condition_group,
        "post_injury_day": day,
        "species": "Mus musculus",
        "data_type": "spatial transcriptomics",
        "dataset_id": "GSE225766",
    }


def load_positions(raw_dir: Path, sample_token: str) -> pd.DataFrame:
    path = raw_dir / f"{sample_token}_tissue_positions_list.csv.gz"
    if not path.exists():
        raise RuntimeError(f"Missing tissue positions file for {sample_token}: {path}")
    pos = pd.read_csv(
        path,
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    return pos.set_index("barcode")


def score_modules(adata: ad.AnnData) -> pd.DataFrame:
    coverage_rows = []
    for module, genes in MOUSE_MODULES.items():
        present = [gene for gene in genes if gene in adata.var_names]
        coverage_rows.append({
            "dataset_id": "GSE225766",
            "module": module,
            "n_genes": len(genes),
            "n_present": len(present),
            "present_genes": ";".join(present),
        })
        if len(present) >= 2:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"{module}_score", use_raw=False)
        else:
            adata.obs[f"{module}_score"] = 0.0
            logging.warning("Too few genes detected for %s: %s", module, present)
    return pd.DataFrame(coverage_rows)


def load_spatial_dataset(raw_dir: Path) -> ad.AnnData:
    objects: list[ad.AnnData] = []
    for h5_path in sorted(raw_dir.glob("*_filtered_feature_bc_matrix.h5")):
        sample_info = parse_sample(h5_path.stem)
        sample_token = sample_info["sample_id"]
        logging.info("Reading %s", h5_path.name)
        sample = sc.read_10x_h5(h5_path)
        sample.var_names_make_unique()
        pos = load_positions(raw_dir, sample_token)
        shared = sample.obs_names.intersection(pos.index)
        if len(shared) == 0:
            raise RuntimeError(f"No matching barcodes between matrix and tissue positions for {sample_token}")
        sample = sample[shared].copy()
        pos = pos.loc[sample.obs_names]
        for col in pos.columns:
            sample.obs[col] = pos[col].values
        sample = sample[sample.obs["in_tissue"].astype(int) == 1].copy()
        sample.obsm["spatial"] = sample.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        for key, value in sample_info.items():
            sample.obs[key] = value
        sample.obs_names = [f"{sample_token}:{barcode}" for barcode in sample.obs_names]
        objects.append(sample)
    if not objects:
        raise RuntimeError(f"No Visium h5 files found in {raw_dir}")
    return ad.concat(objects, join="outer", label="sample_batch", fill_value=0)


def summarize_scores(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_cols = [f"{module}_score" for module in MOUSE_MODULES]
    spot_summary = adata.obs[["sample_id", "condition", "condition_group", *score_cols]].copy()
    sample_summary = (
        spot_summary
        .groupby(["sample_id", "condition", "condition_group"], observed=True)[score_cols]
        .mean()
        .reset_index()
    )
    stats_rows = []
    for module in MOUSE_MODULES:
        col = f"{module}_score"
        groups = [values[col].dropna().to_numpy() for _, values in sample_summary.groupby("condition_group", observed=True) if len(values) > 0]
        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
            h_stat, h_p = kruskal(*groups)
        else:
            h_stat, h_p = np.nan, np.nan
        wt = sample_summary[sample_summary["condition_group"] == "WT"][col]
        mdx = sample_summary[sample_summary["condition_group"] == "MDX"][col]
        injury = sample_summary[sample_summary["condition_group"] == "injury"][col]
        for contrast, other in [("MDX_vs_WT", mdx), ("injury_vs_WT", injury)]:
            if len(wt) and len(other):
                u_stat, u_p = mannwhitneyu(other, wt, alternative="two-sided")
                delta = other.mean() - wt.mean()
            else:
                u_stat, u_p, delta = np.nan, np.nan, np.nan
            stats_rows.append({
                "dataset_id": "GSE225766",
                "module": module,
                "contrast": contrast,
                "test_level": "sample_mean",
                "n_wt_samples": len(wt),
                "n_comparison_samples": len(other),
                "comparison_minus_wt": delta,
                "mannwhitney_u": u_stat,
                "mannwhitney_pvalue": u_p,
                "kruskal_h_across_groups": h_stat,
                "kruskal_p_across_groups": h_p,
            })

    coloc_rows = []
    partner_cols = [
        "adipogenic_score",
        "fibrotic_ecm_score",
        "inflammatory_remodeling_score",
        "fap_stromal_score",
        "macrophage_score",
        "atrophic_muscle_score",
        "myofiber_score",
    ]
    for sample_id, sub in spot_summary.groupby("sample_id", observed=True):
        for partner in partner_cols:
            rho, pval = spearmanr(sub["DA_FAP_score"], sub[partner])
            coloc_rows.append({
                "dataset_id": "GSE225766",
                "sample_id": sample_id,
                "partner_signature": partner.replace("_score", ""),
                "spearman_rho": rho,
                "spearman_pvalue": pval,
                "n_spots": sub.shape[0],
                "test_level": "spot_level_within_sample",
            })
    return sample_summary, pd.DataFrame(stats_rows), pd.DataFrame(coloc_rows)


def save_figures(root: Path, adata: ad.AnnData, sample_summary: pd.DataFrame, coloc: pd.DataFrame) -> None:
    fig_dir = root / "results/figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font_scale=0.95)

    order = ["WT", "MDX", "injury"]
    plt.figure(figsize=(7.2, 4.4))
    sns.boxplot(data=sample_summary, x="condition_group", y="DA_FAP_score", order=order, color="#dce7ef", width=0.55)
    sns.stripplot(data=sample_summary, x="condition_group", y="DA_FAP_score", order=order, hue="condition", size=7, linewidth=0.4, edgecolor="black")
    plt.xlabel("")
    plt.ylabel("Spatial DA-FAP signature score\n(sample mean)")
    plt.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase9_spatial_gse225766_da_score_by_condition").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    corr_summary = (
        coloc.groupby("partner_signature", observed=True)["spearman_rho"]
        .median()
        .reindex(["fap_stromal", "fibrotic_ecm", "adipogenic", "inflammatory_remodeling", "macrophage", "atrophic_muscle", "myofiber"])
        .dropna()
    )
    plt.figure(figsize=(6.6, 3.6))
    sns.barplot(x=corr_summary.values, y=corr_summary.index, color="#8097b1")
    plt.axvline(0, color="#333333", linewidth=0.8)
    plt.xlabel("Median within-sample Spearman rho with DA-FAP score")
    plt.ylabel("")
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase9_spatial_gse225766_da_colocalization").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    top_sample = sample_summary.sort_values("DA_FAP_score", ascending=False).iloc[0]["sample_id"]
    obs = adata.obs[adata.obs["sample_id"] == top_sample].copy()
    plt.figure(figsize=(5.8, 5.4))
    scatter = plt.scatter(
        obs["pxl_col_in_fullres"],
        -obs["pxl_row_in_fullres"],
        c=obs["DA_FAP_score"],
        s=8,
        cmap="magma",
        linewidths=0,
    )
    plt.title(str(top_sample))
    plt.axis("equal")
    plt.axis("off")
    cbar = plt.colorbar(scatter, fraction=0.046, pad=0.02)
    cbar.set_label("DA-FAP score")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase9_spatial_gse225766_top_sample_da_map").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        raw_dir = root / config["phase9"]["spatial_raw_dir"]
        ensure_data(raw_dir)
        adata = load_spatial_dataset(raw_dir)
        logging.info("Loaded spatial object: %d spots x %d genes", adata.n_obs, adata.n_vars)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        coverage = score_modules(adata)

        sample_summary, stats, coloc = summarize_scores(adata)
        table_dir = root / "results/tables"
        table_dir.mkdir(parents=True, exist_ok=True)
        sample_summary.to_csv(root / config["phase9"]["spatial_summary_table"], index=False)
        coverage.to_csv(table_dir / "phase9_spatial_module_gene_coverage.csv", index=False)
        stats.to_csv(table_dir / "phase9_spatial_condition_stats.csv", index=False)
        coloc.to_csv(table_dir / "phase9_spatial_colocalization_stats.csv", index=False)
        adata.write_h5ad(root / config["phase9"]["spatial_object"], compression="gzip")
        save_figures(root, adata, sample_summary, coloc)
        logging.info("Spatial validation completed.")
        return 0
    except Exception as exc:
        logging.exception("Spatial validation failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
