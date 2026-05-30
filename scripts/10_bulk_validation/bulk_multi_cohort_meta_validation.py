#!/usr/bin/env python3
"""Multi-cohort bulk validation and random-effects meta-analysis."""

from __future__ import annotations

import csv
import gzip
import logging
import re
import sys
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy.stats import mannwhitneyu, norm, spearmanr, ttest_ind


MODULES = {
    "DA_FAP": ["PPARG", "CEBPA", "FABP4", "LPL", "CD36", "COL1A1", "COL3A1", "COL6A1", "FN1", "POSTN", "TNC", "THBS2", "LOX", "ACTA2", "IL6", "CCL2", "CXCL12", "CXCL14", "ICAM1", "VCAM1", "SPP1", "MMP2", "MMP14", "TIMP1"],
    "adipogenic": ["PPARG", "CEBPA", "FABP4", "LPL", "ADIPOQ", "PLIN1", "CD36"],
    "fibrotic_ecm": ["COL1A1", "COL3A1", "COL6A1", "FN1", "POSTN", "TNC", "THBS2", "LOX", "ACTA2"],
    "inflammatory_remodeling": ["IL6", "CCL2", "CXCL12", "CXCL14", "ICAM1", "VCAM1", "SPP1", "MMP2", "MMP14", "TIMP1"],
    "PPARG_like": ["PPARG", "CEBPA", "FABP4", "LPL", "CD36"],
    "SMAD3_TGFB": ["SMAD3", "TGFB1", "TGFBR1", "TGFBR2", "COL1A1", "COL3A1", "FN1", "SERPINE1"],
    "STAT3_IL6": ["STAT3", "IL6", "IL6R", "IL6ST", "SOCS3", "JUNB"],
    "myogenesis": ["PAX7", "MYOD1", "MYOG", "MYH1", "MYH2", "MYH7", "ACTA1", "DES"],
    "OXPHOS_mito": ["NDUFA1", "NDUFB8", "COX4I1", "COX5A", "ATP5F1A", "UQCRC1", "SDHB"],
}


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase10b_bulk_meta_validation.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0).fillna(0)


def score_modules(expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.DataFrame(index=expr.columns)
    coverage_rows = []
    for module, genes in MODULES.items():
        present = [gene for gene in genes if gene in expr.index]
        coverage_rows.append({"module": module, "n_genes": len(genes), "n_present": len(present), "present_genes": ";".join(present)})
        if present:
            scores[f"{module}_score"] = zscore_df(expr.loc[present]).mean(axis=0)
        else:
            scores[f"{module}_score"] = 0.0
    return scores, pd.DataFrame(coverage_rows)


def deduplicate_genes(expr: pd.DataFrame) -> pd.DataFrame:
    expr = expr[expr.index.notna()].copy()
    expr.index = expr.index.astype(str).str.split("///").str[0].str.strip()
    expr = expr[(expr.index != "") & (expr.index != "nan")]
    means = expr.mean(axis=1, skipna=True)
    return expr.assign(__mean=means).sort_values("__mean", ascending=False).drop(columns="__mean").groupby(level=0).first()


def read_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.readlines()
    meta_rows: list[list[str]] = []
    table_lines: list[str] = []
    in_table = False
    for line in lines:
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            in_table = False
            continue
        if in_table:
            table_lines.append(line)
        elif line.startswith("!Sample_"):
            meta_rows.append(next(csv.reader([line], delimiter="\t")))
    accessions = []
    titles = []
    for row in meta_rows:
        if row[0] == "!Sample_geo_accession":
            accessions = row[1:]
        if row[0] == "!Sample_title":
            titles = row[1:]
    metadata = pd.DataFrame(index=accessions)
    metadata["sample_id"] = accessions
    if titles:
        metadata["title"] = titles
    char_count = 0
    for row in meta_rows:
        key = row[0]
        values = row[1:]
        if key == "!Sample_source_name_ch1":
            metadata["source_name"] = values
        elif key == "!Sample_characteristics_ch1":
            parsed_keys = []
            parsed_vals = []
            for value in values:
                if ":" in value:
                    k, v = value.split(":", 1)
                    parsed_keys.append(k.strip().lower().replace(" ", "_").replace("/", "_").replace("(", "").replace(")", ""))
                    parsed_vals.append(v.strip())
                else:
                    parsed_keys.append(f"characteristic_{char_count}")
                    parsed_vals.append(value)
            if len(set(parsed_keys)) == 1:
                metadata[parsed_keys[0]] = parsed_vals
            else:
                metadata[f"characteristic_{char_count}"] = values
            char_count += 1
    expr = None
    if table_lines and len(table_lines) > 1:
        expr = pd.read_csv(StringIO("".join(table_lines)), sep="\t")
        if expr.shape[0] == 0 or expr.shape[1] <= 1:
            expr = None
    return metadata, expr


def load_gpl570_symbols(path: Path) -> pd.Series:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("ID\t"))
    annot = pd.read_csv(StringIO("".join(lines[start:])), sep="\t")
    annot["symbol"] = annot["Gene symbol"].astype(str).str.split("///").str[0].str.strip()
    annot = annot[annot["symbol"].ne("") & annot["symbol"].ne("nan")]
    return annot.set_index("ID")["symbol"]


def load_ensg_symbols(root: Path) -> pd.Series:
    path = root / "data/raw/bulk_validation/GSE164471/GSE164471_GESTALT_Muscle_ENSG_counts_annotated.csv.gz"
    annot = pd.read_csv(path, usecols=["Tracking_ID", "HG19en82 Gene Name"])
    annot["Tracking_ID"] = annot["Tracking_ID"].astype(str).str.split(".").str[0]
    annot["symbol"] = annot["HG19en82 Gene Name"].astype(str).str.split("///").str[0].str.strip()
    annot = annot[annot["symbol"].ne("") & annot["symbol"].ne("nan")]
    return annot.drop_duplicates("Tracking_ID").set_index("Tracking_ID")["symbol"]


def load_microarray_dataset(root: Path, dataset_id: str, symbol_map: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    path = root / f"data/raw/bulk_validation/{dataset_id}/{dataset_id}_series_matrix.txt.gz"
    metadata, expr_table = read_geo_series_matrix(path)
    if expr_table is None:
        raise RuntimeError(f"No expression table found in {dataset_id}")
    expr_table = expr_table.rename(columns={expr_table.columns[0]: "ID_REF"}).set_index("ID_REF")
    expr_table = expr_table.apply(pd.to_numeric, errors="coerce")
    expr_table.index = expr_table.index.map(symbol_map)
    expr = deduplicate_genes(expr_table)
    if np.nanmax(expr.to_numpy()) > 50:
        expr = np.log2(expr + 1)
    metadata = metadata.loc[expr.columns].copy()
    metadata["dataset_id"] = dataset_id
    if dataset_id == "GSE25941":
        metadata["contrast_group"] = metadata["age"].astype(str).str.lower().map({"young": "control", "old": "case", "older": "case"})
        metadata["contrast_label"] = "old_or_older_vs_young"
        metadata["phenotype"] = metadata["age"]
    else:
        metadata["contrast_group"] = metadata["age_group"].astype(str).str.lower().map({"young": "control", "old": "case", "older": "case"})
        metadata["contrast_label"] = "old_or_older_vs_young"
        metadata["phenotype"] = metadata["age_group"]
    metadata["platform"] = "GPL570 microarray"
    return expr, metadata, "microarray"


def load_gse111_subseries(root: Path, dataset_id: str, ensg_map: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    raw_dir = root / "data/raw/bulk_validation" / dataset_id
    matrix_path = raw_dir / f"{dataset_id}_series_matrix.txt.gz"
    metadata, _ = read_geo_series_matrix(matrix_path)
    expr_file = next(raw_dir.glob("*lcpmy1*.csv.gz"))
    expr = pd.read_csv(expr_file, index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    title_to_gsm = {}
    for gsm, row in metadata.iterrows():
        clean = re.sub(r"\s+\[[^\]]+\]$", "", str(row["title"]))
        title_to_gsm[clean] = gsm
    expr = expr.rename(columns=title_to_gsm)
    expr = expr.loc[:, [c for c in expr.columns if c in metadata.index]]
    expr.index = expr.index.astype(str).str.split(".").str[0].map(ensg_map)
    expr = deduplicate_genes(expr)
    metadata = metadata.loc[expr.columns].copy()
    metadata["dataset_id"] = dataset_id
    metadata["contrast_group"] = metadata["sarcopenia_status"].astype(str).str.lower().map({"yes": "case", "no": "control"})
    metadata["contrast_label"] = "sarcopenia_vs_control"
    metadata["phenotype"] = metadata["sarcopenia_status"]
    metadata["platform"] = "RNA-seq logCPM"
    return expr, metadata, "rnaseq_logcpm"


def load_gse254300(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    path = root / "data/raw/bulk_validation/GSE254300/GSE254300_fpkm.xlsx"
    expr = pd.read_excel(path)
    expr = expr.rename(columns={expr.columns[0]: "GeneID"}).set_index("GeneID")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = np.log2(expr + 1)
    expr = deduplicate_genes(expr)
    metadata = pd.DataFrame(index=expr.columns)
    metadata["sample_id"] = expr.columns
    metadata["dataset_id"] = "GSE254300"
    metadata["scoliosis_type"] = metadata.index.to_series().str.extract(r"^(CS|IS)", expand=False).map({"CS": "congenital_scoliosis", "IS": "idiopathic_scoliosis"})
    metadata["side"] = metadata.index.to_series().str[-1].map({"A": "A_side", "T": "T_side"})
    metadata["contrast_group"] = metadata["side"].map({"A_side": "control", "T_side": "case"})
    metadata["contrast_label"] = "paraspinal_T_side_vs_A_side"
    metadata["phenotype"] = metadata["scoliosis_type"] + "_" + metadata["side"]
    metadata["platform"] = "RNA-seq FPKM"
    return expr, metadata, "rnaseq_fpkm"


def add_scores(expr: pd.DataFrame, metadata: pd.DataFrame, dataset_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores, coverage = score_modules(expr)
    scores = metadata.join(scores, how="inner")
    scores["dataset_id"] = dataset_id
    coverage.insert(0, "dataset_id", dataset_id)
    return scores.reset_index(drop=True), coverage


def hedges_g(case: pd.Series, control: pd.Series) -> tuple[float, float]:
    case = pd.to_numeric(case, errors="coerce").dropna()
    control = pd.to_numeric(control, errors="coerce").dropna()
    n1, n0 = len(case), len(control)
    if n1 < 2 or n0 < 2:
        return np.nan, np.nan
    s1, s0 = case.std(ddof=1), control.std(ddof=1)
    pooled = np.sqrt(((n1 - 1) * s1**2 + (n0 - 1) * s0**2) / (n1 + n0 - 2))
    if pooled == 0 or not np.isfinite(pooled):
        return np.nan, np.nan
    d = (case.mean() - control.mean()) / pooled
    j = 1 - 3 / (4 * (n1 + n0) - 9)
    g = d * j
    se = np.sqrt((n1 + n0) / (n1 * n0) + (g**2) / (2 * (n1 + n0 - 2)))
    return float(g), float(se)


def summarize_effects(score_rows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset_id, contrast_label), sub in score_rows.groupby(["dataset_id", "contrast_label"], observed=True):
        case = sub[sub["contrast_group"].eq("case")]
        control = sub[sub["contrast_group"].eq("control")]
        if case.empty or control.empty:
            continue
        for module in MODULES:
            col = f"{module}_score"
            g, se = hedges_g(case[col], control[col])
            try:
                _, p_t = ttest_ind(case[col], control[col], equal_var=False, nan_policy="omit")
            except Exception:
                p_t = np.nan
            try:
                _, p_u = mannwhitneyu(case[col].dropna(), control[col].dropna(), alternative="two-sided")
            except Exception:
                p_u = np.nan
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "contrast_label": contrast_label,
                    "module": module,
                    "case_n": len(case),
                    "control_n": len(control),
                    "case_mean": case[col].mean(),
                    "control_mean": control[col].mean(),
                    "case_minus_control": case[col].mean() - control[col].mean(),
                    "hedges_g": g,
                    "se": se,
                    "ci_low": g - 1.96 * se if np.isfinite(se) else np.nan,
                    "ci_high": g + 1.96 * se if np.isfinite(se) else np.nan,
                    "welch_pvalue": p_t,
                    "mannwhitney_pvalue": p_u,
                }
            )
    return pd.DataFrame(rows)


def random_effects_meta(effect_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module, sub in effect_df.dropna(subset=["hedges_g", "se"]).groupby("module", observed=True):
        sub = sub[sub["se"] > 0].copy()
        if sub.empty:
            continue
        yi = sub["hedges_g"].to_numpy(float)
        vi = sub["se"].to_numpy(float) ** 2
        wi = 1 / vi
        fixed = np.sum(wi * yi) / np.sum(wi)
        q = np.sum(wi * (yi - fixed) ** 2)
        df = len(yi) - 1
        c = np.sum(wi) - (np.sum(wi**2) / np.sum(wi))
        tau2 = max(0.0, (q - df) / c) if c > 0 and df > 0 else 0.0
        wre = 1 / (vi + tau2)
        pooled = np.sum(wre * yi) / np.sum(wre)
        se = np.sqrt(1 / np.sum(wre))
        z = pooled / se if se > 0 else np.nan
        p = 2 * norm.sf(abs(z)) if np.isfinite(z) else np.nan
        i2 = max(0.0, (q - df) / q) * 100 if q > 0 and df > 0 else 0.0
        rows.append(
            {
                "module": module,
                "n_contrasts": len(yi),
                "random_effect_hedges_g": pooled,
                "se": se,
                "ci_low": pooled - 1.96 * se,
                "ci_high": pooled + 1.96 * se,
                "pvalue": p,
                "tau2": tau2,
                "i2_percent": i2,
            }
        )
    return pd.DataFrame(rows).sort_values("random_effect_hedges_g", ascending=False)


def load_existing_gse164471(root: Path) -> pd.DataFrame:
    path = root / "data/processed/bulk_validation/GSE164471/GSE164471_module_scores.csv"
    scores = pd.read_csv(path)
    scores["contrast_group"] = np.where(scores["age"] >= 65, "case", np.where(scores["age"] < 50, "control", "middle"))
    scores["contrast_label"] = "old_65_plus_vs_young_under_50"
    scores["phenotype"] = scores["age_group"]
    scores["platform"] = "RNA-seq counts logCPM"
    scores = scores[scores["contrast_group"].isin(["case", "control"])].copy()
    return scores


def save_figures(root: Path, effects: pd.DataFrame, meta: pd.DataFrame, score_rows: pd.DataFrame) -> None:
    fig_dir = root / "results/figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", font_scale=0.95)

    da = effects[effects["module"].eq("DA_FAP")].copy()
    da["label"] = da["dataset_id"] + " (" + da["contrast_label"] + ")"
    da = da.sort_values("hedges_g")
    fig, ax = plt.subplots(figsize=(8.8, max(4, 0.45 * len(da) + 1.8)))
    ax.errorbar(da["hedges_g"], da["label"], xerr=1.96 * da["se"], fmt="o", color="#2d5f7f", ecolor="#7c93a6", capsize=3)
    pooled = meta.loc[meta["module"].eq("DA_FAP")]
    if not pooled.empty:
        row = pooled.iloc[0]
        ax.axvline(row["random_effect_hedges_g"], color="#b22222", linestyle="--", linewidth=1.2, label="random-effects pooled")
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.set_xlabel("Hedges g: degeneration/old/sarcopenia - control")
    ax.set_ylabel("")
    ax.legend(frameon=False, loc="lower right")
    sns.despine(ax=ax)
    fig.tight_layout()
    for suffix in [".png", ".pdf"]:
        fig.savefig((fig_dir / "phase10b_bulk_meta_da_fap_forest").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close(fig)

    heat = effects.pivot_table(index="dataset_id", columns="module", values="hedges_g", aggfunc="first")
    ordered_cols = ["DA_FAP", "adipogenic", "fibrotic_ecm", "inflammatory_remodeling", "PPARG_like", "SMAD3_TGFB", "STAT3_IL6", "myogenesis", "OXPHOS_mito"]
    heat = heat[[c for c in ordered_cols if c in heat.columns]]
    plt.figure(figsize=(10.5, 4.8))
    sns.heatmap(heat, cmap="vlag", center=0, linewidths=0.3, cbar_kws={"label": "Hedges g"})
    plt.xlabel("")
    plt.ylabel("")
    plt.title("Bulk validation effect sizes")
    plt.tight_layout()
    for suffix in [".png", ".pdf"]:
        plt.savefig((fig_dir / "phase10b_bulk_meta_module_effect_heatmap").with_suffix(suffix), bbox_inches="tight", dpi=300)
    plt.close()

    keep = score_rows[score_rows["dataset_id"].isin(["GSE111006", "GSE111010", "GSE111016"])].copy()
    if not keep.empty:
        plt.figure(figsize=(8.5, 4.6))
        sns.boxplot(data=keep, x="dataset_id", y="DA_FAP_score", hue="contrast_group", palette={"control": "#dce7ef", "case": "#d98b73"})
        sns.stripplot(data=keep, x="dataset_id", y="DA_FAP_score", hue="contrast_group", dodge=True, color="#333333", alpha=0.55, size=3, legend=False)
        plt.xlabel("")
        plt.ylabel("DA-FAP score")
        plt.legend(title="", frameon=False)
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase10b_gse111017_sarcopenia_subseries_scores").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        raw_dir = root / config["phase10"]["bulk_raw_dir"]
        gpl570 = load_gpl570_symbols(raw_dir / "GPL570.annot.gz")
        ensg_map = load_ensg_symbols(root)

        all_scores = [load_existing_gse164471(root)]
        all_coverage = []
        for dataset_id in ["GSE25941", "GSE38718"]:
            expr, metadata, _ = load_microarray_dataset(root, dataset_id, gpl570)
            scores, coverage = add_scores(expr, metadata, dataset_id)
            all_scores.append(scores)
            all_coverage.append(coverage)
            out_dir = root / "data/processed/bulk_validation" / dataset_id
            out_dir.mkdir(parents=True, exist_ok=True)
            expr.to_csv(out_dir / f"{dataset_id}_gene_symbol_expression.csv.gz", compression="gzip")
            scores.to_csv(out_dir / f"{dataset_id}_module_scores.csv", index=False)
            logging.info("Processed %s: %d genes x %d samples", dataset_id, expr.shape[0], expr.shape[1])

        for dataset_id in ["GSE111006", "GSE111010", "GSE111016"]:
            expr, metadata, _ = load_gse111_subseries(root, dataset_id, ensg_map)
            scores, coverage = add_scores(expr, metadata, dataset_id)
            all_scores.append(scores)
            all_coverage.append(coverage)
            out_dir = root / "data/processed/bulk_validation" / dataset_id
            out_dir.mkdir(parents=True, exist_ok=True)
            expr.to_csv(out_dir / f"{dataset_id}_gene_symbol_logCPM.csv.gz", compression="gzip")
            scores.to_csv(out_dir / f"{dataset_id}_module_scores.csv", index=False)
            logging.info("Processed %s: %d genes x %d samples", dataset_id, expr.shape[0], expr.shape[1])

        expr, metadata, _ = load_gse254300(root)
        scores, coverage = add_scores(expr, metadata, "GSE254300")
        all_scores.append(scores)
        all_coverage.append(coverage)
        out_dir = root / "data/processed/bulk_validation/GSE254300"
        out_dir.mkdir(parents=True, exist_ok=True)
        expr.to_csv(out_dir / "GSE254300_gene_symbol_log2FPKM.csv.gz", compression="gzip")
        scores.to_csv(out_dir / "GSE254300_module_scores.csv", index=False)

        score_rows = pd.concat(all_scores, ignore_index=True, sort=False)
        effects = summarize_effects(score_rows)
        # Primary meta-analysis excludes the paraspinal scoliosis extension because it is not aging/sarcopenia.
        primary = effects[~effects["dataset_id"].eq("GSE254300")].copy()
        meta = random_effects_meta(primary)

        table_dir = root / "results/tables"
        table_dir.mkdir(parents=True, exist_ok=True)
        score_rows.to_csv(table_dir / "phase10b_bulk_multi_cohort_module_scores.csv", index=False)
        effects.to_csv(table_dir / "phase10b_bulk_multi_cohort_effect_sizes.csv", index=False)
        meta.to_csv(table_dir / "phase10b_bulk_random_effects_meta_analysis.csv", index=False)
        pd.concat(all_coverage, ignore_index=True).to_csv(table_dir / "phase10b_bulk_multi_cohort_gene_coverage.csv", index=False)
        save_figures(root, primary, meta, score_rows)
        logging.info("Bulk multi-cohort meta-validation completed.")
        return 0
    except Exception as exc:
        logging.exception("Bulk meta-validation failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
