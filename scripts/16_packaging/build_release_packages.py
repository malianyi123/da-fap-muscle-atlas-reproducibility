#!/usr/bin/env python3
"""Build data, script, and figure/table release packages for the DA-FAP project."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import tarfile
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "release_packages"
PACKAGE_DATE = "2026-05-29"

EXCLUDED_PARTS = {
    ".git",
    ".venv",
    "__pycache__",
    ".ipynb_checkpoints",
    "release_packages",
}
EXCLUDED_NAMES = {".DS_Store", ".gitkeep"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}


@dataclass(frozen=True)
class PackageSpec:
    name: str
    title: str
    description: str
    sources: tuple[str, ...]
    readme: str


def data_readme() -> str:
    return f"""# DA-FAP 数据资料包

生成日期：{PACKAGE_DATE}

本资料包用于复现 DA-FAP 骨骼肌退变项目的分析和图表生成。它保留了项目相对路径，解压到项目根目录后即可与脚本包配合使用。

## 内容

- `data/metadata/`：数据集清单、下载记录和标准化对象清单。
- `data/raw/`：本项目已下载的公开数据输入，包括单细胞/单核 RNA-seq、空间转录组和 bulk validation 数据。
- `data/processed/`：标准化后的元数据、bulk 表达矩阵和中间数据。
- `results/objects/`：关键 AnnData 对象，包括标准化单细胞对象、整合 atlas、FAP 子集对象和空间打分对象。
- `config/config.yaml`：所有分析脚本读取的路径与参数配置。
- `environment/`：Python 环境和版本记录。

## 使用方式

1. 解压本资料包。
2. 将 `da_fap_python_scripts_package` 解压到同一项目根目录，允许覆盖相同的 `config/` 与 `environment/` 文件。
3. 使用脚本包中的 `scripts/16_packaging/run_release_workflow.py` 重新生成分析结果。

注意：虚拟敲除为基于公开单细胞表达与 GRN 回归的计算预测，不能替代实验敲除因果验证。
"""


def scripts_readme() -> str:
    return f"""# DA-FAP Python 脚本包

生成日期：{PACKAGE_DATE}

本脚本包用于分析 DA-FAP 数据资料包，并重新生成文章所需主图、补图、主表、补表和核心中间结果。

## 主要入口

推荐从项目根目录运行：

```bash
python scripts/16_packaging/run_release_workflow.py
```

查看可运行步骤：

```bash
python scripts/16_packaging/run_release_workflow.py --list-steps
```

只重跑后半部分分析和图表：

```bash
python scripts/16_packaging/run_release_workflow.py \\
  --steps phase4_to_8_fap_core phase8_grn_virtual_knockdown phase10_bulk_meta \\
          phase15_sensitivity phase12_publication_figures phase13_supplementary_tables
```

## 环境

首选使用 `environment/environment.yml` 创建环境；若使用现有虚拟环境，请确认 `scanpy`、`anndata`、`numpy`、`pandas`、`scikit-learn`、`scipy`、`matplotlib`、`seaborn`、`networkx`、`statsmodels` 等依赖可用。

所有分析脚本均从 `config/config.yaml` 读取路径，避免硬编码绝对路径。
"""


def figures_readme() -> str:
    return f"""# DA-FAP 图片及表格资料包

生成日期：{PACKAGE_DATE}

本资料包用于论文撰写和投稿整理，包含文章主图、候选补图、主结果表、补充表和图注文档。

## 内容

- `results/figures/figure1_*` 至 `figure9_*`：主文候选图，提供 PNG/PDF/SVG。
- `results/figures/phase*`：分析阶段图，可作为补充图或审稿回复材料。
- `results/tables/supplementary_table_*`：补充表 1-16。
- `results/tables/phase*`：主分析和敏感性分析的源结果表。
- `manuscript/figure_legends.md`：图注草稿。
- `manuscript/main_text.md`、`abstract.md`、`supplementary_methods.md`：配套文稿材料。

图表使用白底、统一字体和谨慎表述；虚拟扰动相关图表应表述为“计算预测的候选扰动敏感节点”。
"""


PACKAGE_SPECS = [
    PackageSpec(
        name="da_fap_data_package",
        title="数据资料包",
        description="Public input data, processed matrices, standardized metadata, and AnnData objects.",
        sources=(
            "config",
            "environment",
            "data",
            "results/objects",
        ),
        readme=data_readme(),
    ),
    PackageSpec(
        name="da_fap_python_scripts_package",
        title="Python 脚本包",
        description="Modular analysis scripts and workflow driver for regenerating figures and tables.",
        sources=(
            "README.md",
            "config",
            "environment",
            "scripts",
        ),
        readme=scripts_readme(),
    ),
    PackageSpec(
        name="da_fap_figures_tables_package",
        title="图片及表格资料包",
        description="Publication-ready main figures, supplementary figures, source tables, and manuscript files.",
        sources=(
            "README.md",
            "results/figures",
            "results/tables",
            "manuscript",
            "data/metadata/dataset_metadata.csv",
        ),
        readme=figures_readme(),
    ),
]


def should_include(path: Path) -> bool:
    rel = path.relative_to(PROJECT_ROOT)
    if any(part in EXCLUDED_PARTS for part in rel.parts):
        return False
    if path.name in EXCLUDED_NAMES:
        return False
    if path.suffix in EXCLUDED_SUFFIXES:
        return False
    return path.is_file()


def iter_source_files(source: str) -> list[Path]:
    path = PROJECT_ROOT / source
    if not path.exists():
        return []
    if path.is_file():
        return [path] if should_include(path) else []
    return sorted(p for p in path.rglob("*") if should_include(p))


def copy_or_link(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(package_dir: Path) -> tuple[int, int]:
    rows = []
    for path in sorted(package_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(package_dir).as_posix()
        if rel in {"MANIFEST.tsv", "SHA256SUMS.txt"}:
            continue
        rows.append(
            {
                "relative_path": rel,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )

    manifest_path = package_dir / "MANIFEST.tsv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    sha_path = package_dir / "SHA256SUMS.txt"
    with sha_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['sha256']}  {row['relative_path']}\n")

    return len(rows), sum(int(row["size_bytes"]) for row in rows)


def make_archive(package_dir: Path) -> Path:
    archive_path = package_dir.with_suffix(".tar.gz")
    if archive_path.exists():
        archive_path.unlink()
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(package_dir, arcname=package_dir.name)
    return archive_path


def write_top_level_readme(index_rows: list[dict[str, str]]) -> None:
    readme_path = RELEASE_ROOT / "README.md"
    lines = [
        "# DA-FAP release packages",
        "",
        f"Generated on {PACKAGE_DATE}.",
        "",
        "This directory contains three separated deliverable packages requested for manuscript preparation:",
        "",
        "1. `da_fap_data_package`: data inputs and processed AnnData objects.",
        "2. `da_fap_python_scripts_package`: analysis scripts and the workflow driver.",
        "3. `da_fap_figures_tables_package`: article-ready figures, tables, and manuscript support files.",
        "",
        "Each package includes `README.md`, `MANIFEST.tsv`, and `SHA256SUMS.txt`.",
        "",
        "| package | files | uncompressed size | archive |",
        "|---|---:|---:|---|",
    ]
    for row in index_rows:
        lines.append(
            f"| `{row['package']}` | {row['file_count']} | {row['total_size_human']} | `{row['archive']}` |"
        )
    lines.append("")
    readme_path.write_text("\n".join(lines), encoding="utf-8")


def human_size(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024 or unit == "TB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{num_bytes} B"


def build_package(spec: PackageSpec) -> dict[str, str]:
    package_dir = RELEASE_ROOT / spec.name
    if package_dir.exists():
        shutil.rmtree(package_dir)
    archive_path = package_dir.with_suffix(".tar.gz")
    if archive_path.exists():
        archive_path.unlink()

    package_dir.mkdir(parents=True, exist_ok=True)
    (package_dir / "README.md").write_text(spec.readme, encoding="utf-8")

    copied: set[Path] = set()
    for source in spec.sources:
        for source_file in iter_source_files(source):
            rel = source_file.relative_to(PROJECT_ROOT)
            if rel in copied:
                continue
            copied.add(rel)
            copy_or_link(source_file, package_dir / rel)

    file_count, total_size = write_manifest(package_dir)
    archive = make_archive(package_dir)

    return {
        "package": spec.name,
        "title": spec.title,
        "description": spec.description,
        "directory": package_dir.relative_to(PROJECT_ROOT).as_posix(),
        "archive": archive.relative_to(PROJECT_ROOT).as_posix(),
        "file_count": str(file_count),
        "total_size_bytes": str(total_size),
        "total_size_human": human_size(total_size),
        "archive_size_bytes": str(archive.stat().st_size),
        "archive_size_human": human_size(archive.stat().st_size),
    }


def write_index(rows: list[dict[str, str]]) -> None:
    index_path = RELEASE_ROOT / "package_index.tsv"
    fieldnames = [
        "package",
        "title",
        "description",
        "directory",
        "archive",
        "file_count",
        "total_size_bytes",
        "total_size_human",
        "archive_size_bytes",
        "archive_size_human",
    ]
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    RELEASE_ROOT.mkdir(parents=True, exist_ok=True)
    index_rows = []
    for spec in PACKAGE_SPECS:
        print(f"Building {spec.name}...")
        index_rows.append(build_package(spec))
    write_index(index_rows)
    write_top_level_readme(index_rows)

    print("\nRelease packages built:")
    for row in index_rows:
        print(
            f"- {row['package']}: {row['file_count']} files, "
            f"{row['total_size_human']} uncompressed, {row['archive_size_human']} archive"
        )
    print(f"\nIndex: {RELEASE_ROOT / 'package_index.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
