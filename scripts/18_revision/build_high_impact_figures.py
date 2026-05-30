#!/usr/bin/env python3
"""Redesign Figure 1 and Figure 5 for the upgraded submission."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = ROOT / "results" / "figures"
TABLE_DIR = ROOT / "results" / "tables"


BLUE = "#1f4e79"
TEAL = "#2a7f8f"
RED = "#b92b2b"
GOLD = "#b7791f"
GRAY = "#5f6b73"
LIGHT_BLUE = "#eaf2f8"
LIGHT_TEAL = "#e8f4f2"
LIGHT_GOLD = "#fbf1df"
LIGHT_RED = "#f8e7e7"
LIGHT_GRAY = "#f4f6f7"


def save_all(fig: plt.Figure, path_base: Path) -> None:
    path_base.parent.mkdir(parents=True, exist_ok=True)
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(path_base.with_suffix(suffix), bbox_inches="tight", dpi=350)
    plt.close(fig)


def clean_axis(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def card(ax: plt.Axes, xy, wh, title, subtitle="", fill=LIGHT_GRAY, edge=BLUE, title_size=12, sub_size=9.5):
    x, y = xy
    w, h = wh
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        linewidth=1.5,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(box)
    ax.text(x + w * 0.05, y + h * 0.66, title, ha="left", va="center", fontsize=title_size, fontweight="bold", color="#0d2233")
    if subtitle:
        ax.text(x + w * 0.05, y + h * 0.34, subtitle, ha="left", va="center", fontsize=sub_size, color="#38464f")


def arrow(ax: plt.Axes, start, end, color=BLUE, lw=1.8):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=16, color=color, linewidth=lw))


def figure1() -> None:
    fig, ax = plt.subplots(figsize=(14, 7.5), facecolor="white")
    clean_axis(ax)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.5, 0.97, "Figure 1. Study design and evidence architecture", ha="center", va="top", fontsize=19, fontweight="bold", color=BLUE)

    # Data layer
    ax.text(0.05, 0.88, "Public-data inputs", fontsize=13, fontweight="bold", color=BLUE)
    card(ax, (0.05, 0.72), (0.25, 0.12), "Human single-cell/snRNA-seq", "121,520 cells/nuclei\n4 public resources", LIGHT_BLUE)
    card(ax, (0.375, 0.72), (0.25, 0.12), "Focused FAP compartment", "30,065 FAP-like cells\n5 stromal states", LIGHT_TEAL, TEAL)
    card(ax, (0.70, 0.72), (0.25, 0.12), "External validation", "Mouse Visium + 6 human bulk\naging/sarcopenia contrasts", LIGHT_GOLD, GOLD)
    arrow(ax, (0.305, 0.78), (0.365, 0.78))
    arrow(ax, (0.63, 0.78), (0.69, 0.78))

    # Discovery and interrogation layer
    ax.text(0.05, 0.62, "Analysis chain", fontsize=13, fontweight="bold", color=BLUE)
    steps = [
        ("Atlas integration", "QC, Harmony,\ncell annotation", LIGHT_BLUE, BLUE),
        ("DA-FAP definition", "adipogenic + ECM +\ninflammatory modules", LIGHT_RED, RED),
        ("State dynamics", "pseudotime +\nligand-receptor signals", LIGHT_TEAL, TEAL),
        ("GRN perturbation", "virtual knockdown +\nnegative controls", LIGHT_GOLD, GOLD),
        ("Validation", "spatial support +\nbulk meta-analysis", LIGHT_BLUE, BLUE),
    ]
    x0s = [0.05, 0.235, 0.42, 0.605, 0.79]
    for x0, (title, sub, fill, edge) in zip(x0s, steps):
        card(ax, (x0, 0.43), (0.155, 0.135), title, sub, fill, edge, title_size=10.5, sub_size=8.6)
    for i in range(len(x0s) - 1):
        arrow(ax, (x0s[i] + 0.158, 0.50), (x0s[i + 1] - 0.008, 0.50), lw=1.5)

    # Robustness layer
    ax.text(0.05, 0.33, "Robustness layer added for higher-impact submission", fontsize=13, fontweight="bold", color=BLUE)
    robustness = [
        ("Null-label permutation", "tests DA-FAP-specificity\nof regulator effects"),
        ("Matched-expression controls", "tests expression-driven\nfalse positives"),
        ("Bootstrap ranking stability", "top-5/top-10 overlap\nand rank correlation"),
        ("Donor-level aggregation", "checks pseudoreplication\nwhere metadata permit"),
    ]
    xs = [0.05, 0.285, 0.52, 0.755]
    for x0, (title, sub) in zip(xs, robustness):
        card(ax, (x0, 0.15), (0.19, 0.12), title, sub, "#f7f7f7", "#7a8790", title_size=9.9, sub_size=8.2)

    ax.text(
        0.5,
        0.055,
        "Interpretation guardrail: DA-FAPs are modeled as a degeneration-associated FAP-like state; GRN perturbation nominates candidate regulators, not experimental causality.",
        ha="center",
        va="center",
        fontsize=10.2,
        color="#26333b",
    )
    save_all(fig, FIG_DIR / "figure1_study_design_dataset_map")


def figure5() -> None:
    sns.set_theme(style="white", font="DejaVu Sans")
    ranking = pd.read_csv(TABLE_DIR / "phase8b_grn_virtual_knockdown_ranking.csv")
    ranking = ranking.copy()
    ranking["sensitivity"] = ranking["perturbation_sensitivity_score"]
    perm = pd.read_csv(TABLE_DIR / "supplementary_table_17_grn_label_permutation.csv")
    controls = pd.read_csv(TABLE_DIR / "supplementary_table_18_matched_expression_controls.csv")
    boot = pd.read_csv(TABLE_DIR / "supplementary_table_19_regulator_bootstrap_stability.csv")
    edges = pd.read_csv(TABLE_DIR / "phase8b_grn_edges.csv")

    perm_fdr = dict(zip(perm["regulator"].astype(str), perm["perm_fdr"].astype(float)))
    top = ranking.head(10).copy()
    top["perm_fdr"] = top["regulator"].map(perm_fdr)
    top["passes_null"] = top["perm_fdr"] < 0.05

    fig = plt.figure(figsize=(15.8, 11.0), facecolor="white")
    gs = fig.add_gridspec(2, 3, height_ratios=[0.95, 1.05], width_ratios=[1.08, 1.02, 1.0], hspace=0.56, wspace=0.34)
    fig.suptitle("Figure 5. GRN-regression perturbation screen with negative controls", fontsize=18, fontweight="bold", color=BLUE, y=0.985)

    # A: workflow
    ax = fig.add_subplot(gs[0, 0])
    clean_axis(ax)
    ax.set_title("A  Perturbation framework", loc="left", fontsize=13, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    workflow = [
        ("FAP expression", "30,065 cells"),
        ("GRN regressors", "110 target models"),
        ("Simulated TF knockdown", "5th percentile in DA-FAPs"),
        ("Control analyses", "null labels | matched controls | bootstrap"),
    ]
    ys = [0.78, 0.56, 0.34, 0.12]
    for (title, sub), y in zip(workflow, ys):
        card(ax, (0.08, y), (0.80, 0.13), title, sub, LIGHT_BLUE if y > 0.2 else LIGHT_GOLD, BLUE if y > 0.2 else GOLD, 10.5, 8.0)
    for y1, y2 in zip(ys[:-1], ys[1:]):
        arrow(ax, (0.48, y1 - 0.005), (0.48, y2 + 0.14), lw=1.4)

    # B: ranking
    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(len(top))[::-1]
    colors = [RED if flag else "#9aa6ad" for flag in top["passes_null"]]
    ax.barh(y, top["sensitivity"], color=colors, edgecolor="#45525a", height=0.72)
    ax.axvline(0, color="#333333", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(top["regulator"], fontsize=9.5)
    ax.set_xlabel("Predicted DA-FAP attenuation")
    ax.set_title("B  Candidate ranking", loc="left", fontsize=13, fontweight="bold")
    for yi, (_, row) in zip(y, top.iterrows()):
        if row["passes_null"]:
            ax.text(row["sensitivity"] + 0.002, yi, "FDR<0.05", va="center", fontsize=7.8, color=RED)
    sns.despine(ax=ax)

    # C: module heatmap
    ax = fig.add_subplot(gs[0, 2])
    heat = top.set_index("regulator")[
        [
            "predicted_delta_adipogenic",
            "predicted_delta_fibrotic_ecm",
            "predicted_delta_inflammatory_remodeling",
        ]
    ]
    heat.columns = ["Adipogenic", "Fibrotic/ECM", "Inflamm."]
    sns.heatmap(
        heat,
        cmap="vlag",
        center=0,
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": "Predicted delta", "shrink": 0.72},
        ax=ax,
    )
    ax.set_title("C  Program-level response", loc="left", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=35, labelsize=9)

    # D: permutation null
    ax = fig.add_subplot(gs[1, 0])
    plot_perm = perm.sort_values("observed_sensitivity", ascending=True).tail(10)
    yp = np.arange(len(plot_perm))
    ax.barh(yp, plot_perm["observed_sensitivity"], color=[RED if f < 0.05 else "#9aa6ad" for f in plot_perm["perm_fdr"]], edgecolor="#45525a")
    ax.scatter(plot_perm["null_p95"], yp, marker="|", s=160, color="#222222", label="null 95th pct")
    ax.scatter(plot_perm["null_mean"], yp, s=22, color="#222222", label="null mean")
    ax.set_yticks(yp)
    ax.set_yticklabels(plot_perm["regulator"], fontsize=9.5)
    ax.set_xlabel("Sensitivity")
    ax.set_title("D  DA-FAP null-label permutation", loc="left", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    sns.despine(ax=ax)

    # E: matched-expression controls
    ax = fig.add_subplot(gs[1, 1])
    subset = controls[controls["type"].isin(["regulator", "matched_control"])].copy()
    order = ["regulator", "matched_control"]
    sns.boxplot(data=subset, x="type", y="sensitivity", order=order, ax=ax, color="#d9e1e6", width=0.5, fliersize=0)
    sns.stripplot(data=subset, x="type", y="sensitivity", order=order, ax=ax, hue="type", palette=[BLUE, "#9aa6ad"], dodge=False, size=5, alpha=0.85, legend=False)
    control_95 = subset.loc[subset["type"].eq("matched_control"), "sensitivity"].quantile(0.95)
    ax.axhline(control_95, ls="--", color=RED, lw=1)
    ax.text(1.04, control_95 + 0.001, "control 95th", color=RED, fontsize=8)
    ax.set_xticklabels(["Candidate\nregulators", "Matched-expression\ncontrols"])
    ax.set_xlabel("")
    ax.set_ylabel("Predicted DA-FAP attenuation")
    ax.set_title("E  Expression-matched controls", loc="left", fontsize=13, fontweight="bold")
    sns.despine(ax=ax)

    # F: bootstrap and network summary
    ax = fig.add_subplot(gs[1, 2])
    clean_axis(ax)
    ax.set_title("F  Stable prioritized nodes", loc="left", fontsize=13, fontweight="bold")
    freq = boot[boot["record_type"].eq("regulator_selection_frequency")].copy()
    freq = freq.sort_values("top5_selection_freq", ascending=False).head(7)

    G = nx.DiGraph()
    top_regs = list(freq["regulator"].astype(str).head(5))
    edge_subset = (
        edges[edges["regulator"].isin(top_regs)]
        .sort_values("feature_importance", ascending=False)
        .groupby("regulator")
        .head(3)
    )
    for reg in top_regs:
        G.add_node(reg, kind="reg")
    for _, row in edge_subset.iterrows():
        G.add_node(row["target"], kind="target")
        G.add_edge(row["regulator"], row["target"], weight=float(row["feature_importance"]))
    pos = {}
    for i, reg in enumerate(top_regs):
        pos[reg] = (0.18, 0.86 - i * 0.17)
    targets = [n for n, d in G.nodes(data=True) if d["kind"] == "target"]
    for i, target in enumerate(targets):
        pos[target] = (0.72, 0.90 - i * (0.78 / max(1, len(targets) - 1)))
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=True, arrowstyle="-|>", arrowsize=10, width=0.8, alpha=0.45, edge_color="#5d6d75")
    nx.draw_networkx_nodes(G, pos, nodelist=top_regs, node_color=BLUE, node_size=480, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=targets, node_color="#d7dee2", edgecolors="#7f8b93", node_size=220, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7.6, font_color="white", labels={reg: reg for reg in top_regs}, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=6.5, font_color="#1f2933", labels={target: target for target in targets}, ax=ax)
    ax.text(
        0.02,
        0.02,
        "Bootstrap B=300: top-5 Jaccard 0.98; top-10 Jaccard 0.90; Spearman 0.98.\n"
        "Ridge surrogate vs ExtraTrees: Spearman 0.91.",
        transform=ax.transAxes,
        fontsize=8.4,
        color="#26333b",
        va="bottom",
    )

    save_all(fig, FIG_DIR / "figure5_regulatory_network_virtual_knockout")


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    figure1()
    figure5()
    print("Redesigned Figure 1 and Figure 5.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
