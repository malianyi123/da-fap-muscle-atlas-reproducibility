#!/usr/bin/env python3
"""Validate the Phase 0 project structure."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import yaml


REQUIRED_DIRS = [
    "environment",
    "data/raw",
    "data/processed",
    "data/metadata",
    "results/figures",
    "results/tables",
    "results/logs",
    "results/objects",
    "scripts/00_setup",
    "scripts/01_data_search",
    "scripts/02_download",
    "scripts/03_sc_qc_integration",
    "scripts/04_fap_subclustering",
    "scripts/05_trajectory_velocity",
    "scripts/06_cell_communication",
    "scripts/07_regulatory_network",
    "scripts/08_virtual_knockout",
    "scripts/09_spatial_validation",
    "scripts/10_bulk_validation",
    "scripts/11_subtyping",
    "scripts/12_figure_generation",
    "scripts/13_manuscript",
    "notebooks",
    "manuscript",
    "config",
]

REQUIRED_FILES = [
    "README.md",
    "environment/environment.yml",
    "environment/sessionInfo.txt",
    "config/config.yaml",
    "manuscript/main_text.md",
    "manuscript/abstract.md",
    "manuscript/cover_letter.md",
    "manuscript/supplementary_methods.md",
]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def setup_logging(root: Path) -> None:
    log_file = root / "results/logs/phase0_structure_check.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )


def main() -> int:
    root = project_root()
    setup_logging(root)
    config_path = root / "config/config.yaml"

    try:
        with config_path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except FileNotFoundError:
        logging.error("Missing config file: %s", config_path)
        return 1
    except yaml.YAMLError as exc:
        logging.error("Could not parse %s: %s", config_path, exc)
        return 1

    missing_dirs = [item for item in REQUIRED_DIRS if not (root / item).is_dir()]
    missing_files = [item for item in REQUIRED_FILES if not (root / item).is_file()]

    if missing_dirs or missing_files:
        for item in missing_dirs:
            logging.error("Missing directory: %s", item)
        for item in missing_files:
            logging.error("Missing file: %s", item)
        return 1

    logging.info("Project '%s' structure is complete for Phase 0.", config["project"]["name"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

