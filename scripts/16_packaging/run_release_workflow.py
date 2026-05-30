#!/usr/bin/env python3
"""Run the DA-FAP release workflow from packaged data.

This driver assumes the data package and script package have been extracted
into the same project root, preserving the repository folder layout.
"""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
from pathlib import Path


WORKFLOW_STEPS = [
    ("structure_check", "scripts/00_setup/check_project_structure.py"),
    ("phase3_qc_integration", "scripts/03_sc_qc_integration/qc_integrate_annotate.py"),
    ("phase4_to_8_fap_core", "scripts/04_fap_subclustering/fap_focused_analysis.py"),
    ("phase8_grn_virtual_knockdown", "scripts/08_virtual_knockout/grn_regression_virtual_knockout.py"),
    ("phase9_spatial_validation", "scripts/09_spatial_validation/spatial_validation_gse225766.py"),
    ("phase10_bulk_primary", "scripts/10_bulk_validation/bulk_validation_gse164471.py"),
    ("phase10_bulk_meta", "scripts/10_bulk_validation/bulk_multi_cohort_meta_validation.py"),
    ("phase15_sensitivity", "scripts/15_quality_control/sensitivity_analysis.py"),
    ("phase12_publication_figures", "scripts/12_figure_generation/build_publication_figures.py"),
    ("main_tables", "scripts/13_manuscript/build_main_tables.py"),
    ("phase13_supplementary_tables", "scripts/13_manuscript/assemble_supplementary_tables.py"),
]


def find_project_root() -> Path:
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        if (parent / "config" / "config.yaml").exists():
            return parent
    raise FileNotFoundError("Could not locate project root containing config/config.yaml")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--steps",
        nargs="+",
        default=[name for name, _ in WORKFLOW_STEPS],
        help="Step names to run. Use --list-steps to inspect available names.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable to use for child scripts. Defaults to the current interpreter.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue running later steps if a step exits with a non-zero status.",
    )
    parser.add_argument(
        "--list-steps",
        action="store_true",
        help="List available workflow step names and exit.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.list_steps:
        for name, script in WORKFLOW_STEPS:
            print(f"{name}\t{script}")
        return 0

    root = find_project_root()
    logs_dir = root / "results" / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = logs_dir / f"release_workflow_{timestamp}.log"

    selected = set(args.steps)
    known = {name for name, _ in WORKFLOW_STEPS}
    unknown = sorted(selected - known)
    if unknown:
        raise ValueError(f"Unknown workflow step(s): {', '.join(unknown)}")

    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"Project root: {root}\n")
        log.write(f"Python: {args.python}\n")
        log.write(f"Dry run: {args.dry_run}\n\n")

        failures: list[str] = []
        for name, script_rel in WORKFLOW_STEPS:
            if name not in selected:
                continue
            script_path = root / script_rel
            if not script_path.exists():
                message = f"[{name}] missing script: {script_rel}"
                print(message)
                log.write(message + "\n")
                failures.append(name)
                if not args.continue_on_error:
                    break
                continue

            cmd = [args.python, str(script_path)]
            message = f"[{name}] {' '.join(cmd)}"
            print(message)
            log.write(message + "\n")
            if args.dry_run:
                continue

            result = subprocess.run(cmd, cwd=root, text=True, capture_output=True)
            log.write(result.stdout)
            if result.stderr:
                log.write("\n[stderr]\n")
                log.write(result.stderr)
            log.write(f"\n[{name}] exit_code={result.returncode}\n\n")
            log.flush()

            if result.returncode != 0:
                failures.append(name)
                print(f"[{name}] failed; see {log_path}")
                if not args.continue_on_error:
                    break

    if failures:
        print(f"Workflow finished with failures: {', '.join(failures)}")
        print(f"Log: {log_path}")
        return 1

    print(f"Workflow completed successfully. Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
