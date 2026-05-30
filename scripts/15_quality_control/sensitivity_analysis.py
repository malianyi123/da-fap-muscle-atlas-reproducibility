#!/usr/bin/env python3
"""Phase 15 DA-FAP robustness and artifact sensitivity analyses."""

from __future__ import annotations

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
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import adjusted_rand_score, roc_auc_score


MODULES = {
    "adipogenic": ["PPARG", "CEBPA", "FABP4", "LPL", "ADIPOQ", "PLIN1", "CD36"],
    "fibrotic_ecm": ["COL1A1", "COL3A1", "COL6A1", "FN1", "POSTN", "TNC", "THBS2", "LOX", "ACTA2"],
    "inflammatory_remodeling": ["IL6", "CCL2", "CXCL12", "CXCL14", "ICAM1", "VCAM1", "SPP1", "MMP2", "MMP14", "TIMP1"],
}

ARTIFACT_MODULES = {
    "stress_response": ["HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "DUSP1", "IER2", "ATF3", "JUN", "FOS"],
    "cell_cycle": ["MKI67", "TOP2A", "PCNA", "STMN1", "HMGB2", "TYMS", "UBE2C", "MCM5"],
    "apoptosis": ["BAX", "BCL2", "CASP3", "CASP8", "JUN", "FOS", "ATF3"],
}


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase15_sensitivity_analysis.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = np.nanstd(values)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def as_dense(x) -> np.ndarray:
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def prepare_lognorm(fap: ad.AnnData) -> ad.AnnData:
    x = fap.layers["counts"].copy() if "counts" in fap.layers else fap.X.copy()
    tmp = ad.AnnData(X=x, obs=fap.obs.copy(), var=fap.var.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    return tmp


def module_score(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs), []
    x = as_dense(adata[:, present].X)
    zx = (x - np.nanmean(x, axis=0)) / np.nanstd(x, axis=0)
    zx = np.nan_to_num(zx)
    return np.nanmean(zx, axis=1), present


def leave_one_dataset_out(fap: ad.AnnData) -> pd.DataFrame:
    features = ["homeostatic_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z"]
    obs = fap.obs.copy()
    rows = []
    for heldout in sorted(obs["dataset_id"].astype(str).unique()):
        train = obs["dataset_id"].astype(str).ne(heldout)
        test = ~train
        y_train = obs.loc[train, "da_fap_flag"].astype(int)
        y_test = obs.loc[test, "da_fap_flag"].astype(int)
        if y_train.nunique() < 2 or y_test.nunique() < 2:
            rows.append({"heldout_dataset": heldout, "status": "skipped_one_class", "n_test": int(test.sum())})
            continue
        model = LogisticRegression(max_iter=1000, class_weight="balanced")
        model.fit(obs.loc[train, features], y_train)
        prob = model.predict_proba(obs.loc[test, features])[:, 1]
        auc = roc_auc_score(y_test, prob)
        pred = prob >= 0.5
        rows.append(
            {
                "heldout_dataset": heldout,
                "status": "tested",
                "n_test": int(test.sum()),
                "da_fap_fraction_observed": float(y_test.mean()),
                "predicted_da_fap_fraction": float(pred.mean()),
                "auc": float(auc),
                "spearman_predicted_probability_vs_da_score": float(spearmanr(prob, obs.loc[test, "da_fap_score_z"]).correlation),
            }
        )
    return pd.DataFrame(rows)


def module_gene_sensitivity(lognorm: ad.AnnData, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    original = np.asarray(lognorm.obs["da_fap_score_z"], dtype=float)
    rows = []
    for dropped_module in [None, *MODULES.keys()]:
        used_modules = [m for m in MODULES if m != dropped_module]
        scores = []
        present_map = {}
        for module in used_modules:
            score, present = module_score(lognorm, MODULES[module])
            scores.append(zscore(score))
            present_map[module] = present
        alt = zscore(np.nanmean(np.vstack(scores), axis=0))
        rows.append(
            {
                "analysis": "drop_module" if dropped_module else "all_modules_recomputed",
                "removed": dropped_module or "none",
                "n_iterations": 1,
                "spearman_with_original_da_score": float(spearmanr(alt, original).correlation),
                "top_quartile_overlap": float(np.mean((alt >= np.quantile(alt, 0.75)) & (original >= np.quantile(original, 0.75))) / 0.25),
                "genes_used": ";".join(gene for genes in present_map.values() for gene in genes),
            }
        )

    present_by_module = {module: [gene for gene in genes if gene in lognorm.var_names] for module, genes in MODULES.items()}
    for i in range(100):
        scores = []
        genes_used = []
        for module, genes in present_by_module.items():
            if len(genes) <= 2:
                chosen = genes
            else:
                chosen = sorted(rng.choice(genes, size=max(2, int(np.ceil(0.8 * len(genes)))), replace=False).tolist())
            score, present = module_score(lognorm, chosen)
            scores.append(zscore(score))
            genes_used.extend(present)
        alt = zscore(np.nanmean(np.vstack(scores), axis=0))
        rows.append(
            {
                "analysis": "bootstrap_80pct_genes",
                "removed": f"iteration_{i + 1}",
                "n_iterations": 100,
                "spearman_with_original_da_score": float(spearmanr(alt, original).correlation),
                "top_quartile_overlap": float(np.mean((alt >= np.quantile(alt, 0.75)) & (original >= np.quantile(original, 0.75))) / 0.25),
                "genes_used": ";".join(genes_used),
            }
        )
    return pd.DataFrame(rows)


def artifact_checks(root: Path, fap: ad.AnnData, lognorm: ad.AnnData) -> pd.DataFrame:
    phase3 = sc.read_h5ad(root / "results/objects/phase3_human_integrated_initial.h5ad", backed="r")
    qc_cols = [col for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if col in phase3.obs.columns]
    obs = fap.obs.copy()
    obs = obs.join(phase3.obs.loc[obs.index.intersection(phase3.obs_names), qc_cols], how="left")
    rows = []
    for col in qc_cols:
        rho, pval = spearmanr(obs["da_fap_score_z"], obs[col], nan_policy="omit")
        rows.append({"feature": col, "feature_type": "qc_covariate", "n_present_genes": np.nan, "spearman_rho_with_da_fap_score": rho, "pvalue": pval})

    for module, genes in ARTIFACT_MODULES.items():
        score, present = module_score(lognorm, genes)
        rho, pval = spearmanr(lognorm.obs["da_fap_score_z"], score, nan_policy="omit")
        rows.append({"feature": module, "feature_type": "artifact_module", "n_present_genes": len(present), "spearman_rho_with_da_fap_score": rho, "pvalue": pval})

    mt_genes = [gene for gene in lognorm.var_names if str(gene).startswith("MT-")]
    if mt_genes:
        score, present = module_score(lognorm, mt_genes)
        rho, pval = spearmanr(lognorm.obs["da_fap_score_z"], score, nan_policy="omit")
        rows.append({"feature": "mitochondrial_gene_module", "feature_type": "artifact_module", "n_present_genes": len(present), "spearman_rho_with_da_fap_score": rho, "pvalue": pval})
    return pd.DataFrame(rows)


def artifact_adjusted_da_score(root: Path, fap: ad.AnnData) -> pd.DataFrame:
    phase3 = sc.read_h5ad(root / "results/objects/phase3_human_integrated_initial.h5ad", backed="r")
    qc_cols = [col for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if col in phase3.obs.columns]
    obs = fap.obs.copy()
    obs = obs.join(phase3.obs.loc[obs.index.intersection(phase3.obs_names), qc_cols], how="left")
    obs = obs.dropna(subset=qc_cols + ["da_fap_score_z", "da_fap_flag", "dataset_id"]).copy()
    if obs.empty:
        return pd.DataFrame()
    design = obs[qc_cols].copy()
    for col in ["n_genes_by_counts", "total_counts"]:
        design[col] = np.log1p(design[col])
    design = pd.concat([design, pd.get_dummies(obs["dataset_id"].astype(str), prefix="dataset", drop_first=True)], axis=1)
    model = LinearRegression()
    model.fit(design, obs["da_fap_score_z"])
    fitted = model.predict(design)
    residual = obs["da_fap_score_z"] - fitted
    obs["artifact_adjusted_da_fap_residual"] = residual
    da = obs[obs["da_fap_flag"].astype(int).eq(1)]["artifact_adjusted_da_fap_residual"]
    other = obs[obs["da_fap_flag"].astype(int).eq(0)]["artifact_adjusted_da_fap_residual"]
    try:
        auc = roc_auc_score(obs["da_fap_flag"].astype(int), residual)
    except Exception:
        auc = np.nan
    rows = [
        {
            "analysis": "linear_residualization",
            "covariates": ";".join(qc_cols + ["dataset_id"]),
            "n_cells": obs.shape[0],
            "spearman_residual_vs_original_da_score": float(spearmanr(residual, obs["da_fap_score_z"]).correlation),
            "residual_auc_for_phase4_da_fap_label": auc,
            "mean_residual_da_fap": da.mean(),
            "mean_residual_other_fap": other.mean(),
            "da_minus_other_residual": da.mean() - other.mean(),
            "top_quartile_overlap_residual_vs_original": float(
                np.mean((residual >= np.quantile(residual, 0.75)) & (obs["da_fap_score_z"] >= np.quantile(obs["da_fap_score_z"], 0.75))) / 0.25
            ),
        }
    ]
    return pd.DataFrame(rows)


def integration_method_sensitivity(fap: ad.AnnData) -> pd.DataFrame:
    work = fap.copy()
    rows = []
    original_subtype = work.obs["fap_subtype"].astype(str)
    for rep in ["X_pca", "X_pca_harmony"]:
        if rep not in work.obsm:
            continue
        key = f"neighbors_{rep}"
        cluster_key = f"leiden_{rep}"
        sc.pp.neighbors(work, use_rep=rep, key_added=key)
        sc.tl.leiden(work, resolution=0.9, key_added=cluster_key, adjacency=work.obsp[f"{key}_connectivities"])
        cluster_means = work.obs.groupby(cluster_key, observed=True)[["homeostatic_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z", "da_fap_score_z"]].mean()
        da_cut = cluster_means["da_fap_score_z"].quantile(0.75)
        subtype_by_cluster = {}
        for cluster, row in cluster_means.iterrows():
            if row["da_fap_score_z"] >= da_cut and max(row["adipogenic_z"], row["fibrotic_ecm_z"], row["inflammatory_remodeling_z"]) > 0:
                subtype = "DA-FAP"
            elif row["homeostatic_z"] >= max(row["adipogenic_z"], row["fibrotic_ecm_z"], row["inflammatory_remodeling_z"]):
                subtype = "homeostatic FAP"
            elif row["inflammatory_remodeling_z"] >= max(row["adipogenic_z"], row["fibrotic_ecm_z"]):
                subtype = "inflammatory FAP"
            elif row["adipogenic_z"] >= row["fibrotic_ecm_z"]:
                subtype = "adipogenic FAP"
            else:
                subtype = "ECM-remodeling FAP"
            subtype_by_cluster[str(cluster)] = subtype
        subtype = work.obs[cluster_key].astype(str).map(subtype_by_cluster).astype(str)
        rows.append(
            {
                "integration_representation": rep,
                "n_clusters": int(work.obs[cluster_key].nunique()),
                "adjusted_rand_vs_phase4_subtypes": float(adjusted_rand_score(original_subtype, subtype)),
                "da_fap_fraction": float((subtype == "DA-FAP").mean()),
                "da_fap_jaccard_vs_phase4": float(
                    np.sum((subtype == "DA-FAP") & (original_subtype == "DA-FAP")) /
                    np.sum((subtype == "DA-FAP") | (original_subtype == "DA-FAP"))
                ),
            }
        )
    return pd.DataFrame(rows)


def save_figures(root: Path, loo: pd.DataFrame, module_sens: pd.DataFrame, artifacts: pd.DataFrame, integration: pd.DataFrame, adjusted: pd.DataFrame) -> None:
    fig_dir = root / "results/figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font_scale=0.95)

    tested = loo[loo["status"].eq("tested")].copy()
    if not tested.empty:
        plt.figure(figsize=(7.2, 4.2))
        sns.barplot(data=tested, x="auc", y="heldout_dataset", color="#6f93a8")
        plt.axvline(0.5, color="#333333", linewidth=0.8)
        plt.xlabel("Held-out dataset AUC")
        plt.ylabel("")
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase15_leave_one_dataset_out_auc").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()

    plt.figure(figsize=(7.4, 4.4))
    plot_df = module_sens.copy()
    plot_df["group"] = np.where(plot_df["analysis"].eq("bootstrap_80pct_genes"), "bootstrap 80% genes", plot_df["removed"])
    sns.boxplot(data=plot_df, x="group", y="spearman_with_original_da_score", color="#dce7ef")
    sns.stripplot(data=plot_df[~plot_df["analysis"].eq("bootstrap_80pct_genes")], x="group", y="spearman_with_original_da_score", color="#b22222", size=6)
    plt.xticks(rotation=35, ha="right")
    plt.xlabel("")
    plt.ylabel("Spearman rho with original DA-FAP score")
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase15_module_gene_sensitivity").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    plt.figure(figsize=(7.4, 4.6))
    art = artifacts.sort_values("spearman_rho_with_da_fap_score")
    sns.barplot(data=art, x="spearman_rho_with_da_fap_score", y="feature", hue="feature_type", dodge=False)
    plt.axvline(0, color="#333333", linewidth=0.8)
    plt.xlabel("Spearman rho with DA-FAP score")
    plt.ylabel("")
    plt.legend(title="", frameon=False)
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase15_artifact_correlation_checks").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    if not integration.empty:
        plt.figure(figsize=(6.5, 4.2))
        sns.barplot(data=integration, x="integration_representation", y="da_fap_jaccard_vs_phase4", color="#8ba6b8")
        plt.ylim(0, 1)
        plt.xlabel("")
        plt.ylabel("DA-FAP Jaccard vs Phase 4 labels")
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase15_integration_method_sensitivity").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()

    if not adjusted.empty:
        labels = ["Residual vs original rho", "Residual DA-FAP AUC", "Top-quartile overlap"]
        values = [
            adjusted.iloc[0]["spearman_residual_vs_original_da_score"],
            adjusted.iloc[0]["residual_auc_for_phase4_da_fap_label"],
            adjusted.iloc[0]["top_quartile_overlap_residual_vs_original"],
        ]
        plt.figure(figsize=(6.8, 4.0))
        sns.barplot(x=values, y=labels, color="#8aa6a3")
        plt.xlim(0, 1)
        plt.xlabel("Robustness statistic")
        plt.ylabel("")
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase15_artifact_adjusted_da_score").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        seed = int(config["phase4_8"]["random_seed"])
        fap = sc.read_h5ad(root / config["phase4_8"]["fap_object"])
        lognorm = prepare_lognorm(fap)
        for col in ["da_fap_score_z", "fap_subtype", "dataset_id"]:
            lognorm.obs[col] = fap.obs[col]

        loo = leave_one_dataset_out(fap)
        module_sens = module_gene_sensitivity(lognorm, seed)
        artifacts = artifact_checks(root, fap, lognorm)
        adjusted = artifact_adjusted_da_score(root, fap)
        integration = integration_method_sensitivity(fap)

        table_dir = root / "results/tables"
        table_dir.mkdir(parents=True, exist_ok=True)
        loo.to_csv(table_dir / "phase15_leave_one_dataset_out_sensitivity.csv", index=False)
        module_sens.to_csv(table_dir / "phase15_module_gene_sensitivity.csv", index=False)
        artifacts.to_csv(table_dir / "phase15_artifact_correlation_checks.csv", index=False)
        adjusted.to_csv(table_dir / "phase15_artifact_adjusted_da_score.csv", index=False)
        integration.to_csv(table_dir / "phase15_integration_method_sensitivity.csv", index=False)
        save_figures(root, loo, module_sens, artifacts, integration, adjusted)
        logging.info("Sensitivity analysis completed.")
        return 0
    except Exception as exc:
        logging.exception("Sensitivity analysis failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
