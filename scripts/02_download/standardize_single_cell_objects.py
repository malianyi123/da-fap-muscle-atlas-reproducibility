#!/usr/bin/env python3
"""Load downloaded single-cell objects and save standardized AnnData files."""

from __future__ import annotations

import csv
import gzip
import logging
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse


STANDARD_COLUMNS = [
    "dataset_id",
    "sample_id",
    "donor_id",
    "species",
    "tissue",
    "muscle_type",
    "condition",
    "age_group",
    "sex",
    "disease_status",
    "platform",
    "data_type",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase2_standardize_single_cell.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def read_dataset_metadata(root: Path) -> dict[str, dict[str, str]]:
    with (root / "data/metadata/dataset_metadata.csv").open("r", encoding="utf-8", newline="") as handle:
        return {row["dataset_id"]: row for row in csv.DictReader(handle)}


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {col.lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def infer_age_group(value: object) -> str:
    text = str(value).lower()
    if any(token in text for token in ["old", "aged", "elder", "sarcopen"]):
        return "old_or_degenerative"
    if any(token in text for token in ["young", "adult"]):
        return "young_or_adult"
    return "unknown"


def add_standard_metadata(
    adata: ad.AnnData,
    dataset_id: str,
    inventory: dict[str, str],
    sample_id: str | None = None,
) -> ad.AnnData:
    obs = adata.obs

    sample_col = first_existing_column(obs, ["sample_id", "sample", "Sample", "orig.ident", "library_id"])
    donor_col = first_existing_column(obs, ["donor_id", "donor", "Donor", "individual", "subject"])
    sex_col = first_existing_column(obs, ["sex", "Sex", "gender"])
    disease_col = first_existing_column(obs, ["disease", "disease_status", "condition", "Condition"])
    age_col = first_existing_column(obs, ["age_group", "age", "development_stage", "Age"])
    tissue_col = first_existing_column(obs, ["tissue", "Tissue"])

    obs["dataset_id"] = dataset_id
    obs["sample_id"] = obs[sample_col].astype(str) if sample_col else (sample_id or dataset_id)
    obs["donor_id"] = obs[donor_col].astype(str) if donor_col else obs["sample_id"].astype(str)
    obs["species"] = inventory.get("species", "unknown")
    obs["tissue"] = obs[tissue_col].astype(str) if tissue_col else inventory.get("tissue", "unknown")
    obs["muscle_type"] = inventory.get("muscle_type", "unknown")
    obs["condition"] = obs[disease_col].astype(str) if disease_col else inventory.get("disease_or_condition", "unknown")
    obs["age_group"] = obs[age_col].map(infer_age_group) if age_col else inventory.get("age_group", "unknown")
    obs["sex"] = obs[sex_col].astype(str) if sex_col else "unknown"
    obs["disease_status"] = obs[disease_col].astype(str) if disease_col else inventory.get("disease_or_condition", "unknown")
    obs["platform"] = inventory.get("platform", "unknown")
    obs["data_type"] = inventory.get("data_type", "unknown")

    for col in STANDARD_COLUMNS:
        obs[col] = obs[col].astype("category") if col not in {"muscle_type", "platform"} else obs[col].astype(str)
    return adata


def ensure_gene_symbols(adata: ad.AnnData) -> ad.AnnData:
    for col in ["feature_name", "gene_symbols", "gene_symbol", "name", "features"]:
        if col in adata.var.columns:
            values = adata.var[col].astype(str)
            if values.notna().sum() > 0:
                adata.var_names = values
                break
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()
    return adata


def preserve_counts_if_reasonable(adata: ad.AnnData) -> ad.AnnData:
    if sparse.issparse(adata.X):
        data = adata.X.data
    else:
        data = np.asarray(adata.X).ravel()
        data = data[np.isfinite(data)]
    if data.size and np.nanmax(data) > 20:
        adata.layers["counts"] = adata.X.copy()
    return adata


def load_cellxgene_human(raw_root: Path, inventory: dict[str, str]) -> ad.AnnData | None:
    dataset_id = "CELLxGENE_E-MTAB-13874_human_ageing"
    h5ad_paths = sorted((raw_root / dataset_id).glob("*.h5ad"))
    if not h5ad_paths:
        logging.warning("No CELLxGENE H5AD files found for %s", dataset_id)
        return None
    objects: list[ad.AnnData] = []
    for h5ad_path in h5ad_paths:
        logging.info("Reading %s", h5ad_path)
        sample = sc.read_h5ad(h5ad_path)
        sample = ensure_gene_symbols(sample)
        sample = preserve_counts_if_reasonable(sample)
        sample = add_standard_metadata(sample, dataset_id, inventory, sample_id=h5ad_path.stem)
        sample.obs["source_dataset_title"] = h5ad_path.stem
        sample.obs["reference_cell_type"] = "unknown"
        for col in ["cell_type", "author_cell_type", "annotation", "CellType"]:
            if col in sample.obs.columns:
                sample.obs["reference_cell_type"] = sample.obs[col].astype(str)
                break
        sample.obs_names = [f"{h5ad_path.stem}_{barcode}" for barcode in sample.obs_names]
        objects.append(sample)
    return ad.concat(objects, label="source_dataset_title", keys=[p.stem for p in h5ad_paths], join="outer")


def read_gse167_sample_metadata(raw_root: Path) -> pd.DataFrame:
    xlsx = raw_root / "GSE167186_snRNA/GSE167186_SimplifiedMetadataSheet.xlsx"
    if not xlsx.exists():
        return pd.DataFrame()
    meta = pd.read_excel(xlsx)
    meta.columns = [str(col).strip() for col in meta.columns]
    return meta


def match_gse167_metadata(sample_id: str, sample_meta: pd.DataFrame) -> dict[str, str]:
    if sample_meta.empty:
        return {}
    for col in sample_meta.columns:
        matches = sample_meta[sample_meta[col].astype(str).str.contains(sample_id, case=False, regex=False, na=False)]
        if not matches.empty:
            return matches.iloc[0].astype(str).to_dict()
    return {}


def load_gse167(raw_root: Path, inventory: dict[str, str]) -> ad.AnnData | None:
    dataset_id = "GSE167186_snRNA"
    extracted = raw_root / dataset_id / "extracted"
    h5_files = sorted(extracted.glob("*_filtered_feature_bc_matrix.h5"))
    if not h5_files:
        logging.warning("No GSE167186 10x H5 files found in %s", extracted)
        return None
    sample_meta = read_gse167_sample_metadata(raw_root)
    objects: list[ad.AnnData] = []
    for h5_path in h5_files:
        match = re.search(r"_(HM\\d+)_", h5_path.name)
        sample_id = match.group(1) if match else h5_path.stem
        logging.info("Reading %s", h5_path.name)
        sample = sc.read_10x_h5(h5_path)
        sample.var_names_make_unique()
        sample.obs_names = [f"{sample_id}_{barcode}" for barcode in sample.obs_names]
        sample = add_standard_metadata(sample, dataset_id, inventory, sample_id=sample_id)
        matched = match_gse167_metadata(sample_id, sample_meta)
        for col, value in matched.items():
            safe_col = f"gse167_{col}"
            sample.obs[safe_col] = str(value)
        sample.obs["reference_cell_type"] = "unknown"
        sample.layers["counts"] = sample.X.copy()
        objects.append(sample)
    return ad.concat(objects, label="sample_id_from_file", keys=[x.obs["sample_id"].iloc[0] for x in objects], join="outer")


def read_expression_csv(path: Path) -> pd.DataFrame:
    opener = gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r", encoding="utf-8")
    with opener as handle:
        return pd.read_csv(handle, sep=None, engine="python", index_col=0)


def load_gse130(raw_root: Path, inventory: dict[str, str]) -> ad.AnnData | None:
    dataset_id = "GSE130646"
    files = sorted((raw_root / dataset_id / "extracted").glob("*Counts.csv.gz"))
    if not files:
        logging.warning("No GSE130646 count CSV files found.")
        return None
    objects: list[ad.AnnData] = []
    for path in files:
        sample_id = path.name.replace("_Counts.csv.gz", "").replace("GSM", "GSM")
        logging.info("Reading %s", path.name)
        matrix = read_expression_csv(path)
        matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        sample = ad.AnnData(X=sparse.csr_matrix(matrix.T.values))
        sample.var_names = matrix.index.astype(str)
        sample.obs_names = [f"{sample_id}_{cell}" for cell in matrix.columns.astype(str)]
        sample = ensure_gene_symbols(sample)
        sample = add_standard_metadata(sample, dataset_id, inventory, sample_id=sample_id)
        sample.obs["reference_cell_type"] = "unknown"
        sample.layers["counts"] = sample.X.copy()
        objects.append(sample)
    return ad.concat(objects, join="outer")


def load_gse143(raw_root: Path, inventory: dict[str, str]) -> ad.AnnData | None:
    dataset_id = "GSE143704"
    matrix_path = raw_root / dataset_id / "GSE143704_DeMicheli_HumanMuscleAtlas_rawdata.txt.gz"
    meta_path = raw_root / dataset_id / "GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt.gz"
    if not matrix_path.exists() or not meta_path.exists():
        logging.warning("GSE143704 raw matrix or metadata missing.")
        return None
    logging.info("Reading %s", matrix_path.name)
    matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    logging.info("Reading %s", meta_path.name)
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    sample = ad.AnnData(X=sparse.csr_matrix(matrix.T.values))
    sample.var_names = matrix.index.astype(str)
    sample.obs_names = matrix.columns.astype(str)
    shared = sample.obs_names.intersection(meta.index.astype(str))
    if len(shared) > 0:
        sample = sample[shared].copy()
        sample.obs = sample.obs.join(meta.loc[shared].astype(str), how="left")
    sample = ensure_gene_symbols(sample)
    sample = add_standard_metadata(sample, dataset_id, inventory)
    sample.obs["reference_cell_type"] = "unknown"
    for col in ["cell_type", "CellType", "cluster", "annotation"]:
        if col in sample.obs.columns:
            sample.obs["reference_cell_type"] = sample.obs[col].astype(str)
            break
    sample.layers["counts"] = sample.X.copy()
    return sample


def load_gse138(raw_root: Path, inventory: dict[str, str]) -> ad.AnnData | None:
    dataset_id = "GSE138826"
    matrix_path = raw_root / dataset_id / "GSE138826_expression_matrix.txt.gz"
    if not matrix_path.exists():
        logging.warning("GSE138826 expression matrix missing.")
        return None
    logging.warning(
        "Deferring GSE138826: downloaded matrix is a large whitespace-delimited mouse expression table "
        "without companion cell metadata; it should be converted separately for mouse-only support analysis."
    )
    return None


LOADERS = {
    "CELLxGENE_E-MTAB-13874_human_ageing": load_cellxgene_human,
    "GSE167186_snRNA": load_gse167,
    "GSE143704": load_gse143,
    "GSE130646": load_gse130,
    "GSE138826": load_gse138,
}


def write_manifest(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["dataset_id", "status", "n_obs", "n_vars", "species", "output_path", "message"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    root = project_root()
    setup_logging(root)

    try:
        config = read_config(root)
        raw_root = root / config["paths"]["data_raw"]
        out_dir = root / config["phase2"]["standard_objects_dir"]
        out_dir.mkdir(parents=True, exist_ok=True)
        metadata = read_dataset_metadata(root)
        rows: list[dict[str, str]] = []

        for dataset_id in config["phase2"]["selected_single_cell_datasets"]:
            loader = LOADERS.get(dataset_id)
            if loader is None:
                continue
            output_path = out_dir / f"{dataset_id}.h5ad"
            if output_path.exists() and dataset_id != "GSE138826":
                try:
                    existing = sc.read_h5ad(output_path, backed="r")
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "status": "written",
                            "n_obs": existing.n_obs,
                            "n_vars": existing.n_vars,
                            "species": str(existing.obs["species"].iloc[0]) if existing.n_obs else metadata[dataset_id]["species"],
                            "output_path": str(output_path.relative_to(root)),
                            "message": "existing object reused",
                        }
                    )
                    existing.file.close()
                    logging.info("Reusing existing standardized object: %s", output_path)
                    continue
                except Exception:
                    logging.warning("Existing object could not be read; rebuilding: %s", output_path)
            try:
                adata = loader(raw_root, metadata[dataset_id])
                if adata is None:
                    message = "required files not found"
                    if dataset_id == "GSE138826":
                        message = "deferred: mouse text matrix requires separate sparse conversion and is excluded from human Phase 3 integration"
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "status": "missing",
                            "n_obs": 0,
                            "n_vars": 0,
                            "species": metadata[dataset_id]["species"],
                            "output_path": "",
                            "message": message,
                        }
                    )
                    continue
                for col in STANDARD_COLUMNS:
                    if col not in adata.obs.columns:
                        adata.obs[col] = "unknown"
                logging.info("Writing standardized object: %s (%d cells, %d genes)", output_path, adata.n_obs, adata.n_vars)
                adata.write_h5ad(output_path, compression="gzip")
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "status": "written",
                        "n_obs": adata.n_obs,
                        "n_vars": adata.n_vars,
                        "species": str(adata.obs["species"].iloc[0]),
                        "output_path": str(output_path.relative_to(root)),
                        "message": "ok",
                    }
                )
            except Exception as exc:
                logging.exception("Failed to standardize %s: %s", dataset_id, exc)
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "status": "failed",
                        "n_obs": 0,
                        "n_vars": 0,
                        "species": metadata.get(dataset_id, {}).get("species", "unknown"),
                        "output_path": "",
                        "message": str(exc),
                    }
                )

        write_manifest(rows, root / config["phase2"]["standardized_manifest"])
        failures = [row for row in rows if row["status"] == "failed"]
        return 1 if failures else 0
    except Exception as exc:
        logging.exception("Phase 2 standardization failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
