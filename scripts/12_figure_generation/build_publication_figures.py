#!/usr/bin/env python3
"""Assemble manuscript-facing figures from completed analysis outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def save_all(fig: plt.Figure, path_base: Path, dpi: int = 300) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(path_base.with_suffix(suffix), bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def add_image(ax: plt.Axes, path: Path, label: str, title: str | None = None) -> None:
    if not path.exists():
        ax.text(0.5, 0.5, f"Missing\n{path.name}", ha="center", va="center", fontsize=11)
        ax.set_axis_off()
        return
    img = mpimg.imread(path)
    ax.imshow(img)
    ax.set_axis_off()
    ax.text(
        0.01,
        0.98,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    if title:
        ax.set_title(title, fontsize=11, pad=6)


def figure1(root: Path) -> None:
    meta = pd.read_csv(root / "data/metadata/dataset_metadata.csv")
    included = meta[meta["inclusion_decision"].eq("include")].copy()
    counts = included.groupby("data_type").size().sort_values(ascending=False)

    fig = plt.figure(figsize=(12, 6.2), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.45], wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    colors = ["#4f6d7a", "#c0d6df", "#dd6e42", "#e8dab2"]
    ax1.barh(counts.index, counts.values, color=colors[: len(counts)])
    ax1.invert_yaxis()
    ax1.set_xlabel("Included public resources")
    ax1.set_ylabel("")
    ax1.set_title("A. Dataset inventory", loc="left", fontweight="bold")
    for i, value in enumerate(counts.values):
        ax1.text(value + 0.05, i, str(value), va="center", fontsize=10)
    sns.despine(ax=ax1)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title("B. Computational workflow", loc="left", fontweight="bold")
    ax2.set_axis_off()
    steps = [
        ("Public skeletal muscle atlases", "single-cell / snRNA-seq"),
        ("FAP-like stromal subset", "module scoring and subclustering"),
        ("DA-FAP state definition", "adipogenic + ECM + inflammatory programs"),
        ("Trajectory and communication", "pseudotime, LR and niche signals"),
        ("GRN virtual perturbation", "candidate TF nodes and simulated knockdown"),
        ("Spatial and bulk support", "Visium signatures and aging RNA-seq"),
    ]
    y_positions = list(reversed([0.12 + i * 0.15 for i in range(len(steps))]))
    for i, ((header, sub), y) in enumerate(zip(steps, y_positions)):
        ax2.add_patch(plt.Rectangle((0.08, y - 0.045), 0.76, 0.085, facecolor="#f5f7f9", edgecolor="#8fa3ad", linewidth=1))
        ax2.text(0.11, y + 0.012, header, fontsize=10.5, fontweight="bold", va="center")
        ax2.text(0.11, y - 0.018, sub, fontsize=9, color="#46545c", va="center")
        if i < len(steps) - 1:
            ax2.annotate("", xy=(0.46, y_positions[i + 1] + 0.048), xytext=(0.46, y - 0.052), arrowprops=dict(arrowstyle="->", color="#59656f", linewidth=1.2))
    save_all(fig, root / "results/figures/figure1_study_design_dataset_map")


def contact_sheet(root: Path, output: str, panels: list[tuple[str, str, str]], ncols: int = 2, figsize: tuple[float, float] = (12, 9)) -> None:
    nrows = (len(panels) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, facecolor="white")
    axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
    for ax, (label, title, relpath) in zip(axes, panels):
        add_image(ax, root / relpath, label, title)
    for ax in axes[len(panels):]:
        ax.set_axis_off()
    fig.subplots_adjust(wspace=0.04, hspace=0.16)
    save_all(fig, root / f"results/figures/{output}")


def figure8(root: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 6.2), facecolor="white")
    ax.set_axis_off()
    boxes = [
        ("Homeostatic\nFAP-like stroma", 0.17, 0.70, "#eaf3ea"),
        ("Inflammatory\nactivation", 0.50, 0.70, "#eef0f7"),
        ("Adipogenic and\nECM remodeling", 0.83, 0.70, "#f8efe5"),
        ("DA-FAP-like\ncell state", 0.83, 0.34, "#f6e0df"),
        ("Fibro-adipogenic\ndegenerative niche", 0.50, 0.34, "#eadff0"),
        ("MRI-visible fatty/fibrotic\ndegeneration framework", 0.17, 0.34, "#e8edf2"),
    ]
    width, height = 0.25, 0.16
    for text, x, y, color in boxes:
        ax.add_patch(plt.Rectangle((x - width / 2, y - height / 2), width, height, facecolor=color, edgecolor="#59656f", linewidth=1.1))
        ax.text(x, y, text, ha="center", va="center", fontsize=10.5, fontweight="bold")

    arrow = dict(arrowstyle="->", linewidth=1.5, color="#3f4b53", shrinkA=2, shrinkB=2)
    ax.annotate("", xy=(0.50 - width / 2, 0.70), xytext=(0.17 + width / 2, 0.70), arrowprops=arrow)
    ax.annotate("", xy=(0.83 - width / 2, 0.70), xytext=(0.50 + width / 2, 0.70), arrowprops=arrow)
    ax.annotate("", xy=(0.83, 0.34 + height / 2), xytext=(0.83, 0.70 - height / 2), arrowprops=arrow)
    ax.annotate("", xy=(0.50 + width / 2, 0.34), xytext=(0.83 - width / 2, 0.34), arrowprops=arrow)
    ax.annotate("", xy=(0.17 + width / 2, 0.34), xytext=(0.50 - width / 2, 0.34), arrowprops=arrow)

    ax.text(0.18, 0.53, "Loss of regenerative balance", ha="center", fontsize=9.5, color="#46545c")
    ax.text(0.52, 0.53, "PPARG / CEBPA, ECM-integrin,\nCXCL12-CXCR4, macrophage signals", ha="center", fontsize=9.5, color="#46545c")
    ax.text(0.50, 0.16, "Computational perturbation prioritizes candidate sensitivity nodes; experimental validation remains required.", ha="center", fontsize=10, color="#222222")
    ax.text(0.50, 0.08, "Interpretation: DA-FAPs are modeled as a degeneration-associated FAP-like state, not a newly defined cell type.", ha="center", fontsize=10, color="#222222")
    save_all(fig, root / "results/figures/figure8_mechanistic_model")


def main() -> int:
    root = project_root()
    sns.set_theme(style="white", font="DejaVu Sans", font_scale=1.0)
    figure1(root)
    contact_sheet(
        root,
        "figure2_integrated_single_cell_atlas",
        [
            ("A", "Datasets", "results/figures/phase3_umap_by_dataset.png"),
            ("B", "Major cell types", "results/figures/phase3_umap_by_major_cell_type.png"),
            ("C", "Cell type composition", "results/figures/phase3_cell_type_composition_by_dataset.png"),
        ],
        ncols=2,
        figsize=(12, 9),
    )
    contact_sheet(
        root,
        "figure3_fap_subclustering_da_fap_definition",
        [
            ("A", "FAP subtypes", "results/figures/phase4_fap_umap_by_subtype.png"),
            ("B", "DA-FAP module score", "results/figures/phase4_fap_umap_da_score.png"),
            ("C", "Subtype module programs", "results/figures/phase4_fap_module_heatmap.png"),
            ("D", "Subtype composition", "results/figures/phase4_fap_subtype_composition_by_dataset.png"),
        ],
        ncols=2,
        figsize=(12, 10),
    )
    contact_sheet(
        root,
        "figure4_trajectory_cell_communication",
        [
            ("A", "FAP pseudotime", "results/figures/phase5_fap_umap_pseudotime.png"),
            ("B", "DA-FAP-centered ligand-receptor scores", "results/figures/phase6_da_fap_ligand_receptor_heatmap.png"),
        ],
        ncols=2,
        figsize=(12, 5.4),
    )
    contact_sheet(
        root,
        "figure5_regulatory_network_virtual_knockout",
        [
            ("A", "GRN-regression knockdown network", "results/figures/phase8b_grn_virtual_knockdown_network.png"),
            ("B", "Predicted module changes", "results/figures/phase8b_grn_virtual_knockdown_heatmap.png"),
            ("C", "Perturbation ranking", "results/figures/phase8b_grn_virtual_knockdown_ranking.png"),
        ],
        ncols=2,
        figsize=(12, 9),
    )
    contact_sheet(
        root,
        "figure6_spatial_da_fap_niche_support",
        [
            ("A", "Spatial DA-FAP score by condition", "results/figures/phase9_spatial_gse225766_da_score_by_condition.png"),
            ("B", "Spatial co-localization", "results/figures/phase9_spatial_gse225766_da_colocalization.png"),
            ("C", "High DA-FAP sample map", "results/figures/phase9_spatial_gse225766_top_sample_da_map.png"),
        ],
        ncols=2,
        figsize=(12, 9),
    )
    contact_sheet(
        root,
        "figure7_bulk_validation",
        [
            ("A", "DA-FAP meta-analysis", "results/figures/phase10b_bulk_meta_da_fap_forest.png"),
            ("B", "Module effect sizes", "results/figures/phase10b_bulk_meta_module_effect_heatmap.png"),
            ("C", "GSE111017 sarcopenia subseries", "results/figures/phase10b_gse111017_sarcopenia_subseries_scores.png"),
        ],
        ncols=2,
        figsize=(12, 9),
    )
    figure8(root)
    contact_sheet(
        root,
        "figure9_sensitivity_analysis",
        [
            ("A", "Leave-one-dataset-out", "results/figures/phase15_leave_one_dataset_out_auc.png"),
            ("B", "Module gene sensitivity", "results/figures/phase15_module_gene_sensitivity.png"),
            ("C", "Artifact correlations", "results/figures/phase15_artifact_correlation_checks.png"),
            ("D", "QC/dataset-adjusted score", "results/figures/phase15_artifact_adjusted_da_score.png"),
            ("E", "Integration method sensitivity", "results/figures/phase15_integration_method_sensitivity.png"),
        ],
        ncols=2,
        figsize=(12, 12),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
