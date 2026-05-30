#!/usr/bin/env python3
"""Write a compact Python package version manifest."""

from __future__ import annotations

import importlib.metadata
import platform
from pathlib import Path


PACKAGES = [
    "anndata",
    "scanpy",
    "scvi-tools",
    "scvelo",
    "squidpy",
    "celloracle",
    "decoupler",
    "gseapy",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "seaborn",
    "scikit-learn",
]


def version_for(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    out = root / "environment/python_package_versions.tsv"
    rows = ["package\tversion", f"python\t{platform.python_version()}"]
    rows.extend(f"{package}\t{version_for(package)}" for package in PACKAGES)
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
