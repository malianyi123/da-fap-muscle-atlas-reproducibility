#!/usr/bin/env python3
"""Download selected Phase 2 public single-cell/snRNA-seq processed matrices."""

from __future__ import annotations

import csv
import logging
import sys
import tarfile
from pathlib import Path

import requests
import yaml
from tqdm import tqdm


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_config(root: Path) -> dict:
    with (root / "config/config.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase2_download.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def download_file(url: str, output_path: Path) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.stat().st_size > 0:
        logging.info("Already present: %s", output_path)
        return "present"

    tmp_path = output_path.with_suffix(output_path.suffix + ".part")
    logging.info("Downloading %s -> %s", url, output_path)
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        total = int(response.headers.get("content-length", 0))
        with tmp_path.open("wb") as handle, tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            desc=output_path.name,
            leave=False,
        ) as progress:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    progress.update(len(chunk))
    tmp_path.replace(output_path)
    return "downloaded"


def get_cellxgene_assets(config: dict) -> list[dict[str, str]]:
    collection_id = config["phase2"]["cellxgene"]["collection_id"]
    max_bytes = int(config["phase2"]["cellxgene"].get("max_h5ad_asset_bytes", 0))
    preferred = {config["phase2"]["cellxgene"]["preferred_dataset_title"]}
    preferred.update(config["phase2"]["cellxgene"].get("supportive_dataset_titles", []))
    url = f"https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    collection = response.json()

    assets: list[dict[str, str]] = []
    for dataset in collection.get("datasets", []):
        title = dataset.get("title", "")
        if title not in preferred:
            continue
        for asset in dataset.get("assets", []):
            if asset.get("filetype") == "H5AD":
                size = int(asset.get("filesize", 0))
                if max_bytes and size > max_bytes:
                    logging.info("Skipping large CELLxGENE asset %s (%d bytes > limit %d).", title, size, max_bytes)
                    continue
                safe_title = title.lower().replace(" ", "_").replace("(", "").replace(")", "")
                assets.append(
                    {
                        "dataset_id": "CELLxGENE_E-MTAB-13874_human_ageing",
                        "source_dataset_title": title,
                        "file_name": f"{safe_title}.h5ad",
                        "url": asset["url"],
                        "source": "CELLxGENE",
                    }
                )
    if not assets:
        raise RuntimeError("No CELLxGENE H5AD assets matched configured dataset titles.")
    return assets


def extract_tar(tar_path: Path, output_dir: Path) -> None:
    marker = output_dir / ".extract_complete"
    if marker.exists():
        logging.info("Archive already extracted: %s", tar_path)
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.info("Extracting %s -> %s", tar_path, output_dir)
    with tarfile.open(tar_path, "r") as archive:
        archive.extractall(output_dir)
    marker.write_text("ok\n", encoding="utf-8")


def write_manifest(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["dataset_id", "source", "source_dataset_title", "file_name", "url", "local_path", "status"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote download manifest: %s", path)


def main() -> int:
    root = project_root()
    setup_logging(root)

    try:
        config = read_config(root)
        raw_root = root / config["paths"]["data_raw"]
        manifest_rows: list[dict[str, str]] = []

        for asset in get_cellxgene_assets(config):
            dataset_dir = raw_root / asset["dataset_id"]
            local_path = dataset_dir / asset["file_name"]
            status = download_file(asset["url"], local_path)
            manifest_rows.append({**asset, "local_path": str(local_path.relative_to(root)), "status": status})

        for dataset_id, entry in config["phase2"]["geo_downloads"].items():
            dataset_dir = raw_root / dataset_id
            for file_name in entry["files"]:
                url = f"{entry['base_url'].rstrip('/')}/{file_name}"
                local_path = dataset_dir / file_name
                status = download_file(url, local_path)
                manifest_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "source": "GEO",
                        "source_dataset_title": "",
                        "file_name": file_name,
                        "url": url,
                        "local_path": str(local_path.relative_to(root)),
                        "status": status,
                    }
                )
                if file_name.endswith(".tar"):
                    extract_tar(local_path, dataset_dir / "extracted")

        write_manifest(manifest_rows, root / config["phase2"]["download_manifest"])
        logging.info("Phase 2 download completed.")
        return 0
    except Exception as exc:
        logging.exception("Phase 2 download failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
