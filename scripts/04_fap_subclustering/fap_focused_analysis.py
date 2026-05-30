#!/usr/bin/env python3
"""FAP-focused Phase 4-8 analysis.

This script performs an executable, cautious computational pass:
- subset FAP-like stromal cells from the human integrated atlas;
- rebuild a full-gene FAP AnnData object from standardized datasets;
- score homeostatic, adipogenic, fibrotic/ECM, and inflammatory modules;
- subcluster and define DA-FAP-like states;
- infer pseudotime with DPT;
- score DA-FAP-centered ligand-receptor pairs;
- rank candidate regulators and virtual-KO proxy perturbations.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr


MODULES = {
    "homeostatic": ["PDGFRA", "PI16", "DPP4", "CXCL12", "DCN", "LUM"],
    "adipogenic": ["PPARG", "CEBPA", "FABP4", "LPL", "ADIPOQ", "PLIN1", "CD36"],
    "fibrotic_ecm": ["COL1A1", "COL3A1", "COL6A1", "FN1", "POSTN", "TNC", "THBS2", "LOX", "ACTA2"],
    "inflammatory_remodeling": ["IL6", "CCL2", "CXCL12", "CXCL14", "ICAM1", "VCAM1", "SPP1", "MMP2", "MMP14", "TIMP1"],
}

REGULATORS = ["PPARG", "SMAD3", "STAT3", "CEBPA", "JUN", "FOS", "KLF4", "KLF5", "RUNX1", "RUNX2", "FOXO1", "FOXO3", "HIF1A", "TEAD1", "YAP1", "WWTR1"]

LR_PAIRS = [
    ("TGFB1", "TGFBR1", "TGFB"),
    ("TGFB1", "TGFBR2", "TGFB"),
    ("PDGFA", "PDGFRA", "PDGF"),
    ("PDGFB", "PDGFRB", "PDGF"),
    ("IL6", "IL6R", "IL6"),
    ("IL6", "IL6ST", "IL6"),
    ("CCL2", "CCR2", "CCL"),
    ("CXCL12", "CXCR4", "CXCL12"),
    ("SPP1", "CD44", "SPP1"),
    ("TNF", "TNFRSF1A", "TNF"),
    ("JAG1", "NOTCH1", "NOTCH"),
    ("WNT5A", "FZD4", "WNT"),
    ("COL1A1", "ITGA1", "ECM_integrin"),
    ("COL1A1", "ITGB1", "ECM_integrin"),
    ("FN1", "ITGA5", "ECM_integrin"),
    ("FN1", "ITGB1", "ECM_integrin"),
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase4_8_fap_analysis.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(arr)
    return (arr - np.nanmean(arr)) / sd


def expression_vector(adata: ad.AnnData, gene: str) -> np.ndarray | None:
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].X
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def load_full_fap_object(root: Path, config: dict) -> ad.AnnData:
    phase3 = sc.read_h5ad(root / config["phase3"]["integrated_object"])
    fap_obs = phase3.obs[phase3.obs["major_cell_type"].astype(str) == "fap_stromal"].copy()
    if fap_obs.empty:
        raise RuntimeError("No FAP/stromal cells found in Phase 3 object.")
    fap_obs["global_umap_1"] = phase3.obsm["X_umap"][phase3.obs_names.isin(fap_obs.index), 0]
    fap_obs["global_umap_2"] = phase3.obsm["X_umap"][phase3.obs_names.isin(fap_obs.index), 1]

    objects: list[ad.AnnData] = []
    for dataset_id in fap_obs["dataset_id"].astype(str).unique():
        path = root / config["phase2"]["standard_objects_dir"] / f"{dataset_id}.h5ad"
        if not path.exists():
            logging.warning("Missing standardized source object for %s", dataset_id)
            continue
        full = sc.read_h5ad(path)
        names = fap_obs.index[fap_obs["dataset_id"].astype(str) == dataset_id]
        names = names.intersection(full.obs_names)
        if len(names) == 0:
            continue
        sample = full[names].copy()
        sample.obs = sample.obs.join(fap_obs.loc[names, ["major_cell_type", "leiden_phase3", "global_umap_1", "global_umap_2"]], how="left", rsuffix="_phase3")
        objects.append(sample)
        logging.info("Loaded %d FAP/stromal cells from %s", sample.n_obs, dataset_id)

    if not objects:
        raise RuntimeError("Could not rebuild FAP object from standardized datasets.")
    fap = ad.concat(objects, join="inner")
    fap.obs_names_make_unique()
    fap.var_names_make_unique()

    max_cells = int(config["phase4_8"]["max_fap_cells_per_dataset"])
    rng = np.random.default_rng(int(config["phase4_8"]["random_seed"]))
    selected: list[str] = []
    for dataset_id, idx in fap.obs.groupby("dataset_id", observed=True).indices.items():
        names = np.asarray(fap.obs_names)[idx]
        if len(names) > max_cells:
            names = rng.choice(names, max_cells, replace=False)
        selected.extend(names.tolist())
    fap = fap[selected].copy()
    return fap


def normalize_score_integrate(fap: ad.AnnData, config: dict) -> ad.AnnData:
    if "counts" not in fap.layers:
        fap.layers["counts"] = fap.X.copy()
    sc.pp.normalize_total(fap, target_sum=1e4)
    sc.pp.log1p(fap)

    present_rows = []
    for module, genes in MODULES.items():
        present = [gene for gene in genes if gene in fap.var_names]
        present_rows.append({"module": module, "n_genes": len(genes), "n_present": len(present), "present_genes": ";".join(present)})
        if len(present) >= 2:
            sc.tl.score_genes(fap, present, score_name=f"{module}_score", use_raw=False)
        else:
            fap.obs[f"{module}_score"] = 0.0
            logging.warning("Too few genes for module %s: %s", module, present)
        fap.obs[f"{module}_z"] = zscore(fap.obs[f"{module}_score"])

    fap.obs["da_fap_score"] = (
        fap.obs["adipogenic_z"] + fap.obs["fibrotic_ecm_z"] + fap.obs["inflammatory_remodeling_z"]
    ) / 3.0
    fap.obs["da_fap_score_z"] = zscore(fap.obs["da_fap_score"])
    fap.uns["module_gene_coverage"] = pd.DataFrame(present_rows).to_dict("list")

    try:
        sc.pp.highly_variable_genes(fap, n_top_genes=int(config["phase4_8"]["n_top_hvg"]), batch_key="dataset_id", flavor="seurat_v3")
    except Exception as exc:
        logging.warning("seurat_v3 HVG failed (%s); using seurat flavor.", exc)
        sc.pp.highly_variable_genes(fap, n_top_genes=int(config["phase4_8"]["n_top_hvg"]), batch_key="dataset_id", flavor="seurat")
    forced_genes = sorted({gene for genes in MODULES.values() for gene in genes} | set(REGULATORS))
    forced_present = [gene for gene in forced_genes if gene in fap.var_names]
    fap.var.loc[forced_present, "highly_variable"] = True
    logging.info("Forced retention of %d module/regulator genes for downstream interpretation.", len(forced_present))
    fap = fap[:, fap.var["highly_variable"]].copy()

    sc.pp.scale(fap, max_value=10)
    sc.tl.pca(fap, n_comps=min(int(config["phase4_8"]["n_pcs"]), fap.n_vars - 1), svd_solver="arpack")
    try:
        import harmonypy as hm

        harmony = hm.run_harmony(fap.obsm["X_pca"], fap.obs, ["dataset_id"], verbose=True)
        corrected = harmony.Z_corr
        if corrected.shape == fap.obsm["X_pca"].shape[::-1]:
            corrected = corrected.T
        if corrected.shape == fap.obsm["X_pca"].shape:
            fap.obsm["X_pca_harmony"] = corrected
            use_rep = "X_pca_harmony"
        else:
            logging.warning("Unexpected Harmony shape %s; using PCA.", corrected.shape)
            use_rep = "X_pca"
    except Exception as exc:
        logging.warning("Harmony failed: %s; using PCA.", exc)
        use_rep = "X_pca"

    sc.pp.neighbors(fap, use_rep=use_rep)
    sc.tl.umap(fap, random_state=int(config["phase4_8"]["random_seed"]))
    sc.tl.leiden(fap, resolution=0.9, key_added="fap_leiden")
    fap.uns["fap_integration_representation"] = use_rep
    return fap


def assign_fap_states(fap: ad.AnnData, config: dict) -> ad.AnnData:
    cluster_means = fap.obs.groupby("fap_leiden", observed=True)[
        ["homeostatic_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z", "da_fap_score_z"]
    ].mean()
    da_cut = cluster_means["da_fap_score_z"].quantile(float(config["phase4_8"]["da_fap_quantile"]))
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
        elif row["fibrotic_ecm_z"] > 0.5:
            subtype = "fibrotic/ECM-remodeling FAP"
        else:
            subtype = "ECM-remodeling FAP"
        subtype_by_cluster[str(cluster)] = subtype
    fap.obs["fap_subtype"] = fap.obs["fap_leiden"].astype(str).map(subtype_by_cluster).astype("category")
    fap.obs["da_fap_flag"] = (fap.obs["fap_subtype"].astype(str) == "DA-FAP").astype(int)
    fap.uns["fap_cluster_state_means"] = cluster_means.reset_index().to_dict("list")
    return fap


def run_pseudotime(fap: ad.AnnData) -> ad.AnnData:
    try:
        root_score = zscore(fap.obs["homeostatic_score"]) - zscore(fap.obs["da_fap_score"])
        fap.uns["iroot"] = int(np.argmax(root_score))
        sc.tl.diffmap(fap)
        sc.tl.dpt(fap)
        fap.obs["fap_pseudotime"] = fap.obs["dpt_pseudotime"].astype(float)
    except Exception as exc:
        logging.warning("DPT failed: %s; using PC1-oriented fallback.", exc)
        pc1 = fap.obsm["X_pca"][:, 0]
        corr = spearmanr(pc1, fap.obs["da_fap_score_z"]).correlation
        if np.isfinite(corr) and corr < 0:
            pc1 = -pc1
        fap.obs["fap_pseudotime"] = (pc1 - pc1.min()) / (pc1.max() - pc1.min())
    return fap


def differential_markers(fap: ad.AnnData, out_path: Path) -> None:
    try:
        sc.tl.rank_genes_groups(fap, groupby="fap_subtype", method="wilcoxon", pts=True)
        rows = []
        result = fap.uns["rank_genes_groups"]
        groups = result["names"].dtype.names
        for group in groups:
            for rank, gene in enumerate(result["names"][group][:100], start=1):
                rows.append(
                    {
                        "fap_subtype": group,
                        "rank": rank,
                        "gene": gene,
                        "score": result["scores"][group][rank - 1],
                        "logfoldchange": result["logfoldchanges"][group][rank - 1],
                        "pvals_adj": result["pvals_adj"][group][rank - 1],
                    }
                )
        pd.DataFrame(rows).to_csv(out_path, index=False)
    except Exception as exc:
        logging.warning("Marker test failed: %s", exc)


def build_lr_expression(root: Path, config: dict, fap: ad.AnnData) -> pd.DataFrame:
    phase3 = sc.read_h5ad(root / config["phase3"]["integrated_object"], backed="r")
    obs = phase3.obs[["dataset_id", "major_cell_type"]].copy()
    phase3.file.close()
    obs["communication_group"] = obs["major_cell_type"].astype(str)
    fap_groups = fap.obs[["fap_subtype"]].copy()
    obs.loc[obs.index.intersection(fap_groups.index), "communication_group"] = fap_groups.loc[obs.index.intersection(fap_groups.index), "fap_subtype"].astype(str)

    genes = sorted({g for pair in LR_PAIRS for g in pair[:2] if g})
    pieces = []
    for dataset_id in obs["dataset_id"].astype(str).unique():
        path = root / config["phase2"]["standard_objects_dir"] / f"{dataset_id}.h5ad"
        if not path.exists():
            continue
        full = sc.read_h5ad(path)
        names = obs.index[(obs["dataset_id"].astype(str) == dataset_id) & obs.index.isin(full.obs_names)]
        present_genes = [g for g in genes if g in full.var_names]
        if len(names) == 0 or len(present_genes) == 0:
            continue
        mini = full[names, present_genes].copy()
        mini.obs = obs.loc[names].copy()
        pieces.append(mini)
    if not pieces:
        return pd.DataFrame()
    lr = ad.concat(pieces, join="outer")
    sc.pp.normalize_total(lr, target_sum=1e4)
    sc.pp.log1p(lr)
    rows = []
    for group, idx in lr.obs.groupby("communication_group", observed=True).indices.items():
        sub = lr[list(idx), :]
        mean_expr = np.asarray(sub.X.mean(axis=0)).ravel()
        for gene, val in zip(sub.var_names, mean_expr):
            rows.append({"communication_group": group, "gene": gene, "mean_log_expr": float(val), "n_cells": sub.n_obs})
    return pd.DataFrame(rows)


def score_lr_pairs(expr: pd.DataFrame) -> pd.DataFrame:
    if expr.empty:
        return expr
    lookup = expr.set_index(["communication_group", "gene"])["mean_log_expr"].to_dict()
    groups = sorted(expr["communication_group"].unique())
    rows = []
    for sender in groups:
        for receiver in groups:
            for ligand, receptor, pathway in LR_PAIRS:
                ligand_expr = lookup.get((sender, ligand), 0.0)
                receptor_expr = lookup.get((receiver, receptor), 0.0)
                if ligand_expr <= 0 or receptor_expr <= 0:
                    continue
                rows.append(
                    {
                        "sender": sender,
                        "receiver": receiver,
                        "ligand": ligand,
                        "receptor": receptor,
                        "pathway": pathway,
                        "ligand_mean_log_expr": ligand_expr,
                        "receptor_mean_log_expr": receptor_expr,
                        "interaction_score": ligand_expr * receptor_expr,
                    }
                )
    return pd.DataFrame(rows).sort_values("interaction_score", ascending=False)


def regulator_and_vko(fap: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    module_cols = ["da_fap_score_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z"]
    reg_rows = []
    vko_rows = []
    edge_rows = []
    target_genes = sorted({gene for genes in MODULES.values() for gene in genes if gene in fap.var_names})
    for tf in REGULATORS:
        tf_expr = expression_vector(fap, tf)
        if tf_expr is None:
            reg_rows.append({"regulator": tf, "status": "missing"})
            continue
        row = {"regulator": tf, "status": "tested", "mean_expression": float(np.mean(tf_expr))}
        vko = {"regulator": tf, "interpretation": "proxy virtual inhibition; negative predicted_delta means module attenuation"}
        attenuation = []
        for col in module_cols:
            corr, pval = spearmanr(tf_expr, fap.obs[col])
            corr = 0.0 if not np.isfinite(corr) else float(corr)
            row[f"spearman_{col}"] = corr
            row[f"p_{col}"] = float(pval) if np.isfinite(pval) else np.nan
            delta = -corr
            vko[f"predicted_delta_{col}"] = delta
            if col == "da_fap_score_z":
                vko["perturbation_sensitivity_score"] = corr
            attenuation.append(corr)
        vko["mean_positive_module_correlation"] = float(np.mean([max(x, 0) for x in attenuation]))
        reg_rows.append(row)
        vko_rows.append(vko)
        for target in target_genes:
            target_expr = expression_vector(fap, target)
            if target_expr is None:
                continue
            corr, pval = spearmanr(tf_expr, target_expr)
            if np.isfinite(corr) and abs(corr) >= 0.15:
                edge_rows.append({"regulator": tf, "target": target, "spearman": float(corr), "pvalue": float(pval) if np.isfinite(pval) else np.nan})

    reg = pd.DataFrame(reg_rows)
    vko = pd.DataFrame(vko_rows).sort_values("perturbation_sensitivity_score", ascending=False)
    edges = pd.DataFrame(edge_rows).sort_values("spearman", key=lambda s: s.abs(), ascending=False) if edge_rows else pd.DataFrame()
    return reg, vko, edges


def save_umap(fap: ad.AnnData, color: str, out_prefix: Path, title: str) -> None:
    sc.pl.umap(fap, color=color, frameon=False, size=10, title=title, show=False, legend_loc="right margin")
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_score_heatmap(fap: ad.AnnData, out_prefix: Path) -> None:
    cols = ["homeostatic_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z", "da_fap_score_z"]
    data = fap.obs.groupby("fap_subtype", observed=True)[cols].mean()
    plt.figure(figsize=(7, 3.6))
    sns.heatmap(data, cmap="vlag", center=0, linewidths=0.5, cbar_kws={"label": "Mean z-score"})
    plt.ylabel("")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_abundance_plot(comp: pd.DataFrame, out_prefix: Path) -> None:
    pivot = comp.pivot(index="dataset_id", columns="fap_subtype", values="fraction").fillna(0)
    ax = pivot.plot(kind="bar", stacked=True, figsize=(8, 4.2), colormap="tab20")
    ax.set_ylabel("Fraction of FAP-like cells")
    ax.set_xlabel("")
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_lr_heatmap(lr: pd.DataFrame, out_prefix: Path) -> None:
    if lr.empty:
        return
    centered = lr[(lr["sender"].eq("DA-FAP")) | (lr["receiver"].eq("DA-FAP"))].copy()
    if centered.empty:
        centered = lr.head(40).copy()
    centered["pair"] = centered["ligand"] + "-" + centered["receptor"]
    top = centered.sort_values("interaction_score", ascending=False).head(40)
    mat = top.pivot_table(index="pair", columns="sender", values="interaction_score", aggfunc="max").fillna(0)
    plt.figure(figsize=(8, max(4, 0.18 * mat.shape[0])))
    sns.heatmap(mat, cmap="mako", linewidths=0.2, cbar_kws={"label": "LR score"})
    plt.ylabel("")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_vko_heatmap(vko: pd.DataFrame, out_prefix: Path) -> None:
    if vko.empty:
        return
    cols = [c for c in vko.columns if c.startswith("predicted_delta_")]
    data = vko.set_index("regulator")[cols].head(12)
    data.columns = [c.replace("predicted_delta_", "").replace("_z", "") for c in data.columns]
    plt.figure(figsize=(7, 4.4))
    sns.heatmap(data, cmap="vlag", center=0, linewidths=0.4, cbar_kws={"label": "Predicted delta after inhibition"})
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def save_network(edges: pd.DataFrame, out_prefix: Path) -> None:
    if edges.empty:
        return
    top = edges.head(40)
    graph = nx.DiGraph()
    for row in top.itertuples(index=False):
        graph.add_edge(row.regulator, row.target, weight=abs(row.spearman), sign=np.sign(row.spearman))
    plt.figure(figsize=(7, 5.5))
    pos = nx.spring_layout(graph, seed=13, k=0.65)
    regulators = set(top["regulator"])
    colors = ["#d62728" if node in regulators else "#7f7f7f" for node in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=colors, node_size=[650 if n in regulators else 280 for n in graph.nodes], alpha=0.9)
    nx.draw_networkx_edges(graph, pos, arrows=True, width=[1 + 3 * graph[u][v]["weight"] for u, v in graph.edges], edge_color="#555555", alpha=0.55)
    nx.draw_networkx_labels(graph, pos, font_size=8)
    plt.axis("off")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig(out_prefix.with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        fig_dir = root / config["paths"]["figures"]
        table_dir = root / config["paths"]["tables"]
        fig_dir.mkdir(parents=True, exist_ok=True)
        table_dir.mkdir(parents=True, exist_ok=True)

        fap = load_full_fap_object(root, config)
        fap = normalize_score_integrate(fap, config)
        fap = assign_fap_states(fap, config)
        fap = run_pseudotime(fap)

        object_path = root / config["phase4_8"]["fap_object"]
        object_path.parent.mkdir(parents=True, exist_ok=True)
        logging.info("Writing FAP object: %s", object_path)
        fap.write_h5ad(object_path, compression="gzip")
        fap.obs.to_csv(root / config["phase4_8"]["fap_metadata"], compression="gzip")

        module_coverage = pd.DataFrame(fap.uns["module_gene_coverage"])
        module_coverage.to_csv(table_dir / "phase4_da_fap_module_gene_coverage.csv", index=False)
        comp = fap.obs.groupby(["dataset_id", "fap_subtype"], observed=True).size().reset_index(name="n_cells")
        comp["fraction"] = comp["n_cells"] / comp.groupby("dataset_id", observed=True)["n_cells"].transform("sum")
        comp.to_csv(table_dir / "phase4_fap_subtype_composition.csv", index=False)
        fap.obs.groupby("fap_subtype", observed=True)[["homeostatic_z", "adipogenic_z", "fibrotic_ecm_z", "inflammatory_remodeling_z", "da_fap_score_z", "fap_pseudotime"]].mean().to_csv(table_dir / "phase4_fap_subtype_module_means.csv")
        differential_markers(fap, table_dir / "phase4_fap_subtype_markers.csv")

        da = fap.obs[fap.obs["da_fap_flag"] == 1]["fap_pseudotime"]
        other = fap.obs[fap.obs["da_fap_flag"] == 0]["fap_pseudotime"]
        stat, pval = mannwhitneyu(da, other, alternative="two-sided") if len(da) and len(other) else (np.nan, np.nan)
        pd.DataFrame(
            [
                {
                    "comparison": "DA-FAP vs non-DA-FAP pseudotime",
                    "n_da_fap": len(da),
                    "n_other_fap": len(other),
                    "median_da_fap_pseudotime": float(np.median(da)) if len(da) else np.nan,
                    "median_other_pseudotime": float(np.median(other)) if len(other) else np.nan,
                    "mannwhitney_u": stat,
                    "pvalue": pval,
                }
            ]
        ).to_csv(table_dir / "phase5_pseudotime_stats.csv", index=False)

        expr = build_lr_expression(root, config, fap)
        expr.to_csv(table_dir / "phase6_lr_group_mean_expression.csv", index=False)
        lr = score_lr_pairs(expr)
        lr.to_csv(table_dir / "phase6_ligand_receptor_scores.csv", index=False)

        reg, vko, edges = regulator_and_vko(fap)
        reg.to_csv(table_dir / "phase7_regulator_correlations.csv", index=False)
        vko.to_csv(table_dir / "phase8_virtual_ko_proxy_ranking.csv", index=False)
        edges.to_csv(table_dir / "phase7_tf_target_proxy_edges.csv", index=False)

        save_umap(fap, "fap_subtype", fig_dir / "phase4_fap_umap_by_subtype", "FAP-like cell states")
        save_umap(fap, "da_fap_score_z", fig_dir / "phase4_fap_umap_da_score", "DA-FAP composite score")
        save_umap(fap, "fap_pseudotime", fig_dir / "phase5_fap_umap_pseudotime", "FAP pseudotime")
        save_score_heatmap(fap, fig_dir / "phase4_fap_module_heatmap")
        save_abundance_plot(comp, fig_dir / "phase4_fap_subtype_composition_by_dataset")
        save_lr_heatmap(lr, fig_dir / "phase6_da_fap_ligand_receptor_heatmap")
        save_vko_heatmap(vko, fig_dir / "phase8_virtual_ko_proxy_heatmap")
        save_network(edges, fig_dir / "phase7_regulator_target_proxy_network")

        logging.info("FAP-focused Phase 4-8 analysis completed.")
        return 0
    except Exception as exc:
        logging.exception("FAP-focused analysis failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
