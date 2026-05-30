#!/usr/bin/env python3
"""Phase 10 bulk validation on GSE164471 healthy-aging skeletal muscle RNA-seq."""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy.stats import mannwhitneyu, spearmanr


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


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase10_bulk_validation.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def parse_sample_name(name: str) -> dict[str, str | int]:
    base = name.replace("_COUNT", "")
    age = re.search(r"AGE(\d+)", base)
    sex = re.search(r"_([MF])_GROUP", base)
    group = re.search(r"GROUP([A-Z0-9_]+?)_", base)
    age_value = int(age.group(1)) if age else np.nan
    return {
        "sample_id": base,
        "age": age_value,
        "sex": sex.group(1) if sex else "unknown",
        "age_group": "young_20_49" if age_value < 50 else ("middle_50_64" if age_value < 65 else "old_65_plus"),
        "dataset_id": "GSE164471",
    }


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0).fillna(0)


def module_score(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(0.0, index=expr.columns), []
    z = zscore_df(expr.loc[present])
    return z.mean(axis=0), present


def main() -> int:
    root = project_root()
    setup_logging(root)
    try:
        config = read_config(root)
        raw_dir = root / config["phase10"]["bulk_raw_dir"] / "GSE164471"
        count_path = raw_dir / "GSE164471_GESTALT_Muscle_ENSG_counts_annotated.csv.gz"
        if not count_path.exists():
            raise RuntimeError(f"Missing GSE164471 counts file: {count_path}")

        logging.info("Reading %s", count_path)
        counts = pd.read_csv(count_path)
        sample_cols = [c for c in counts.columns if c.endswith("_COUNT")]
        if "gene_name" in counts.columns:
            gene_col = "gene_name"
        elif "Gene" in counts.columns:
            gene_col = "Gene"
        else:
            # The provided file uses a trailing HGNC symbol column named gene_name in current GEO export.
            candidates = [c for c in counts.columns if counts[c].astype(str).str.fullmatch(r"[A-Z0-9.-]+").mean() > 0.5]
            gene_col = candidates[-1]
        counts = counts[counts[gene_col].notna()].copy()
        counts[gene_col] = counts[gene_col].astype(str)
        expr_counts = counts.set_index(gene_col)[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        expr_counts = expr_counts.groupby(expr_counts.index).sum()

        lib_sizes = expr_counts.sum(axis=0)
        log_cpm = np.log2(expr_counts.div(lib_sizes, axis=1) * 1e6 + 1)
        log_cpm.columns = [c.replace("_COUNT", "") for c in log_cpm.columns]

        metadata = pd.DataFrame([parse_sample_name(c) for c in sample_cols]).set_index("sample_id")
        log_cpm = log_cpm.loc[:, metadata.index]

        score_rows = metadata.copy()
        coverage_rows = []
        for module, genes in MODULES.items():
            score, present = module_score(log_cpm, genes)
            score_rows[f"{module}_score"] = score
            coverage_rows.append({"dataset_id": "GSE164471", "module": module, "n_genes": len(genes), "n_present": len(present), "present_genes": ";".join(present)})

        long_rows = []
        for module in MODULES:
            score_col = f"{module}_score"
            corr, pval = spearmanr(score_rows["age"], score_rows[score_col])
            young = score_rows[score_rows["age"] < 50][score_col]
            old = score_rows[score_rows["age"] >= 65][score_col]
            u, p_u = mannwhitneyu(old, young, alternative="two-sided") if len(young) and len(old) else (np.nan, np.nan)
            long_rows.append(
                {
                    "dataset_id": "GSE164471",
                    "module": module,
                    "n_samples": score_rows.shape[0],
                    "spearman_age_rho": corr,
                    "spearman_age_pvalue": pval,
                    "young_n": len(young),
                    "old_n": len(old),
                    "young_mean": young.mean(),
                    "old_mean": old.mean(),
                    "old_minus_young": old.mean() - young.mean(),
                    "mannwhitney_u": u,
                    "mannwhitney_pvalue": p_u,
                }
            )

        out_dir = root / config["phase10"]["bulk_processed_dir"] / "GSE164471"
        out_dir.mkdir(parents=True, exist_ok=True)
        log_cpm.to_csv(out_dir / "GSE164471_logCPM_gene_symbol.csv.gz", compression="gzip")
        score_rows.to_csv(out_dir / "GSE164471_module_scores.csv", index=True)
        pd.DataFrame(coverage_rows).to_csv(root / "results/tables/phase10_bulk_module_gene_coverage.csv", index=False)
        results = pd.DataFrame(long_rows)
        results.to_csv(root / config["phase10"]["bulk_results_table"], index=False)

        fig_dir = root / config["paths"]["figures"]
        fig_dir.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(8, 4.4))
        sns.boxplot(data=score_rows.reset_index(), x="age_group", y="DA_FAP_score", order=["young_20_49", "middle_50_64", "old_65_plus"], color="#d9e6f2")
        sns.stripplot(data=score_rows.reset_index(), x="age_group", y="DA_FAP_score", order=["young_20_49", "middle_50_64", "old_65_plus"], color="#333333", size=4, alpha=0.7)
        plt.xlabel("")
        plt.ylabel("Bulk DA-FAP score")
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase10_bulk_gse164471_da_score_by_age").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()

        plt.figure(figsize=(5.5, 4.2))
        sns.regplot(data=score_rows.reset_index(), x="age", y="DA_FAP_score", scatter_kws={"s": 28, "alpha": 0.75}, line_kws={"color": "#b22222"})
        plt.xlabel("Age")
        plt.ylabel("Bulk DA-FAP score")
        sns.despine()
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase10_bulk_gse164471_da_score_age_correlation").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()

        heat = score_rows[[f"{m}_score" for m in MODULES]].copy()
        heat.columns = list(MODULES)
        heat = heat.sort_index(key=lambda idx: score_rows.loc[idx, "age"])
        plt.figure(figsize=(11, 4.8))
        sns.heatmap(heat.T, cmap="vlag", center=0, xticklabels=False, cbar_kws={"label": "Module score"})
        plt.xlabel("Samples ordered by age")
        plt.ylabel("")
        plt.tight_layout()
        for suffix in [".png", ".pdf"]:
            plt.savefig((fig_dir / "phase10_bulk_gse164471_module_heatmap").with_suffix(suffix), bbox_inches="tight", dpi=300)
        plt.close()

        logging.info("Bulk validation completed.")
        return 0
    except Exception as exc:
        logging.exception("Bulk validation failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

