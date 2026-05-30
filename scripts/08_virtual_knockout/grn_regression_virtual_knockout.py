#!/usr/bin/env python3
"""GRN-regression in silico knockdown for DA-FAP candidate regulators.

This workflow is stronger than the earlier correlation-only proxy. It trains
directed TF-to-target regression models from FAP-like expression data and then
simulates TF knockdown in DA-FAP cells by perturbing candidate regulator
expression and predicting downstream target/module changes.

The analysis remains computational and should be described as GRN-regression
virtual perturbation, not as experimental causality.
"""

from __future__ import annotations

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
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor


MODULES = {
    "adipogenic": ["PPARG", "CEBPA", "FABP4", "LPL", "ADIPOQ", "PLIN1", "CD36"],
    "fibrotic_ecm": ["COL1A1", "COL3A1", "COL6A1", "FN1", "POSTN", "TNC", "THBS2", "LOX", "ACTA2"],
    "inflammatory_remodeling": ["IL6", "CCL2", "CXCL12", "CXCL14", "ICAM1", "VCAM1", "SPP1", "MMP2", "MMP14", "TIMP1"],
}

REGULATORS = [
    "PPARG", "SMAD3", "STAT3", "CEBPA", "JUN", "FOS", "KLF4", "KLF5",
    "RUNX1", "RUNX2", "FOXO1", "FOXO3", "HIF1A", "TEAD1", "YAP1", "WWTR1",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase8b_grn_virtual_knockout.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def as_dense(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def expression_matrix(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [gene for gene in genes if gene in adata.var_names]
    return as_dense(adata[:, present].X).astype(np.float32)


def stratified_sample(obs: pd.DataFrame, group_col: str, max_cells: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected: list[str] = []
    groups = obs[group_col].astype(str).fillna("unknown")
    per_group = max(50, int(np.ceil(max_cells / groups.nunique())))
    for _, names in obs.groupby(groups, observed=True).groups.items():
        arr = np.asarray(list(names))
        if len(arr) > per_group:
            arr = rng.choice(arr, size=per_group, replace=False)
        selected.extend(arr.tolist())
    if len(selected) > max_cells:
        selected = rng.choice(np.asarray(selected), size=max_cells, replace=False).tolist()
    return np.asarray(selected)


def prepare_lognorm_fap(root: Path, config: dict) -> ad.AnnData:
    path = root / config["phase4_8"]["fap_object"]
    if not path.exists():
        raise RuntimeError(f"Missing FAP object: {path}")
    fap = sc.read_h5ad(path)
    if "counts" in fap.layers:
        fap.X = fap.layers["counts"].copy()
    sc.pp.normalize_total(fap, target_sum=1e4)
    sc.pp.log1p(fap)
    return fap


def target_genes(root: Path, fap: ad.AnnData) -> list[str]:
    module_genes = {gene for genes in MODULES.values() for gene in genes}
    marker_path = root / "results/tables/phase4_fap_subtype_markers.csv"
    marker_genes: set[str] = set()
    if marker_path.exists():
        markers = pd.read_csv(marker_path)
        marker_genes = set(
            markers.loc[markers["fap_subtype"].astype(str).eq("DA-FAP")]
            .sort_values("rank")
            .head(80)["gene"].astype(str)
        )
    genes = sorted((module_genes | marker_genes | set(REGULATORS)).intersection(set(fap.var_names)))
    return genes


def fit_grn_models(
    fap: ad.AnnData,
    train_names: np.ndarray,
    regulators: list[str],
    targets: list[str],
    seed: int,
) -> tuple[dict[str, tuple[ExtraTreesRegressor, list[str]]], pd.DataFrame, np.ndarray]:
    train = fap[train_names].copy()
    x_tf = expression_matrix(train, regulators)
    models: dict[str, tuple[ExtraTreesRegressor, list[str]]] = {}
    edge_rows = []
    logging.info("Training GRN regressors for %d targets using %d cells and %d regulators.", len(targets), train.n_obs, len(regulators))
    for i, target in enumerate(targets, start=1):
        y = expression_matrix(train, [target]).ravel()
        if np.nanstd(y) < 1e-4:
            continue
        predictors = [regulator for regulator in regulators if regulator != target]
        if not predictors:
            continue
        predictor_indices = [regulators.index(regulator) for regulator in predictors]
        model = ExtraTreesRegressor(
            n_estimators=220,
            max_depth=8,
            min_samples_leaf=20,
            random_state=seed + i,
            n_jobs=-1,
        )
        model.fit(x_tf[:, predictor_indices], y)
        models[target] = (model, predictors)
        for regulator, importance in zip(predictors, model.feature_importances_):
            if importance <= 0:
                continue
            rho = spearmanr(x_tf[:, regulators.index(regulator)], y).correlation
            if not np.isfinite(rho):
                rho = 0.0
            edge_rows.append(
                {
                    "regulator": regulator,
                    "target": target,
                    "feature_importance": float(importance),
                    "spearman_sign": float(np.sign(rho)),
                    "signed_grn_weight": float(importance * np.sign(rho)),
                    "spearman_rho": float(rho),
                    "target_is_module_gene": target in {gene for genes in MODULES.values() for gene in genes},
                }
            )
    edges = pd.DataFrame(edge_rows)
    return models, edges, x_tf


def simulate_knockdown(
    fap: ad.AnnData,
    da_names: np.ndarray,
    regulators: list[str],
    targets: list[str],
    models: dict[str, tuple[ExtraTreesRegressor, list[str]]],
    train_tf: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    da = fap[da_names].copy()
    x_da = expression_matrix(da, regulators)
    knockdown_floor = np.quantile(train_tf, 0.05, axis=0)
    target_rows = []
    module_rows = []

    for reg_idx, regulator in enumerate(regulators):
        x_pert = x_da.copy()
        x_pert[:, reg_idx] = np.minimum(x_pert[:, reg_idx], knockdown_floor[reg_idx])
        target_delta: dict[str, float] = {}
        for target in targets:
            if target == regulator and target in fap.var_names:
                original = expression_matrix(da, [target]).ravel()
                delta = float(np.mean(np.minimum(original, knockdown_floor[reg_idx]) - original))
                target_delta[target] = delta
                target_rows.append(
                    {
                        "regulator": regulator,
                        "target": target,
                        "mean_predicted_delta_log_expression": delta,
                        "mean_predicted_baseline_log_expression": float(np.mean(original)),
                        "n_da_fap_cells": da.n_obs,
                        "effect_type": "direct_knockdown",
                    }
                )
                continue
            model_record = models.get(target)
            if model_record is None:
                continue
            model, predictors = model_record
            predictor_indices = [regulators.index(gene) for gene in predictors]
            base_pred = model.predict(x_da[:, predictor_indices])
            pert_pred = model.predict(x_pert[:, predictor_indices])
            delta = float(np.mean(pert_pred - base_pred))
            target_delta[target] = delta
            target_rows.append(
                    {
                        "regulator": regulator,
                        "target": target,
                        "mean_predicted_delta_log_expression": delta,
                        "mean_predicted_baseline_log_expression": float(np.mean(base_pred)),
                        "n_da_fap_cells": da.n_obs,
                        "effect_type": "grn_regression_prediction",
                    }
                )

        per_module = {}
        for module, genes in MODULES.items():
            present = [gene for gene in genes if gene in target_delta]
            value = float(np.mean([target_delta[gene] for gene in present])) if present else np.nan
            per_module[module] = value
        da_delta = float(np.nanmean([per_module["adipogenic"], per_module["fibrotic_ecm"], per_module["inflammatory_remodeling"]]))
        module_rows.append(
            {
                "regulator": regulator,
                "perturbation": "GRN-regression simulated knockdown to fifth-percentile expression in DA-FAP cells",
                "predicted_delta_da_fap_program": da_delta,
                "predicted_delta_adipogenic": per_module["adipogenic"],
                "predicted_delta_fibrotic_ecm": per_module["fibrotic_ecm"],
                "predicted_delta_inflammatory_remodeling": per_module["inflammatory_remodeling"],
                "perturbation_sensitivity_score": -da_delta,
                "interpretation": "negative delta indicates predicted attenuation; computational prediction only",
            }
        )
    ranking = pd.DataFrame(module_rows).sort_values("perturbation_sensitivity_score", ascending=False)
    return ranking, pd.DataFrame(target_rows)


def save_figures(root: Path, ranking: pd.DataFrame, target_effects: pd.DataFrame, edges: pd.DataFrame) -> None:
    fig_dir = root / "results/figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font_scale=0.95)

    heat = ranking.set_index("regulator")[
        ["predicted_delta_adipogenic", "predicted_delta_fibrotic_ecm", "predicted_delta_inflammatory_remodeling", "predicted_delta_da_fap_program"]
    ]
    plt.figure(figsize=(8.6, 5.8))
    sns.heatmap(heat, cmap="vlag", center=0, linewidths=0.25, cbar_kws={"label": "Predicted delta"})
    plt.xlabel("")
    plt.ylabel("")
    plt.title("GRN-regression virtual knockdown")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase8b_grn_virtual_knockdown_heatmap").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    top = ranking.head(10).copy()
    plt.figure(figsize=(7.2, 4.6))
    sns.barplot(data=top, y="regulator", x="perturbation_sensitivity_score", color="#7393a7")
    plt.axvline(0, color="#333333", linewidth=0.8)
    plt.xlabel("Predicted DA-FAP attenuation score")
    plt.ylabel("")
    sns.despine()
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase8b_grn_virtual_knockdown_ranking").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    top_regs = ranking.head(6)["regulator"].tolist()
    top_edges = (
        edges[edges["regulator"].isin(top_regs)]
        .assign(abs_weight=lambda x: x["signed_grn_weight"].abs())
        .sort_values("abs_weight", ascending=False)
        .head(55)
    )
    graph = nx.DiGraph()
    for row in top_edges.itertuples(index=False):
        graph.add_edge(row.regulator, row.target, weight=row.signed_grn_weight)
    if graph.number_of_edges():
        plt.figure(figsize=(9, 7))
        pos = nx.spring_layout(graph, seed=11, k=0.55)
        node_colors = ["#b85c5c" if node in top_regs else "#d9e6ef" for node in graph.nodes]
        nx.draw_networkx_nodes(graph, pos, node_color=node_colors, edgecolors="#333333", linewidths=0.5, node_size=650)
        widths = [0.8 + 4 * abs(graph[u][v]["weight"]) for u, v in graph.edges]
        colors = ["#b22222" if graph[u][v]["weight"] > 0 else "#315f99" for u, v in graph.edges]
        nx.draw_networkx_edges(graph, pos, width=widths, edge_color=colors, alpha=0.65, arrows=True, arrowsize=10)
        nx.draw_networkx_labels(graph, pos, font_size=8)
        plt.axis("off")
        plt.title("Top GRN edges for prioritized perturbation nodes")
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase8b_grn_virtual_knockdown_network").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        seed = int(config["phase4_8"]["random_seed"])
        fap = prepare_lognorm_fap(root, config)
        regulators = [gene for gene in REGULATORS if gene in fap.var_names]
        targets = target_genes(root, fap)
        if len(regulators) < 3 or len(targets) < 10:
            raise RuntimeError("Insufficient regulators or target genes for GRN perturbation.")

        train_names = stratified_sample(fap.obs, "fap_subtype", max_cells=10000, seed=seed)
        da_obs = fap.obs[fap.obs["fap_subtype"].astype(str).eq("DA-FAP")]
        if da_obs.empty:
            raise RuntimeError("No DA-FAP cells found for perturbation simulation.")
        da_names = stratified_sample(da_obs, "dataset_id", max_cells=6000, seed=seed + 1)

        models, edges, train_tf = fit_grn_models(fap, train_names, regulators, targets, seed)
        ranking, target_effects = simulate_knockdown(fap, da_names, regulators, targets, models, train_tf)

        table_dir = root / "results/tables"
        table_dir.mkdir(parents=True, exist_ok=True)
        ranking.to_csv(table_dir / "phase8b_grn_virtual_knockdown_ranking.csv", index=False)
        target_effects.to_csv(table_dir / "phase8b_grn_virtual_knockdown_target_effects.csv", index=False)
        edges.to_csv(table_dir / "phase8b_grn_edges.csv", index=False)
        save_figures(root, ranking, target_effects, edges)
        logging.info("GRN virtual knockdown completed with %d models and %d edges.", len(models), len(edges))
        return 0
    except Exception as exc:
        logging.exception("GRN virtual knockdown failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
