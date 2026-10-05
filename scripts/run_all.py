#!/usr/bin/env python
"""Rebuild the manuscript analysis from ingested section files.

Inputs are the public-data-derived ``*.raw.h5ad`` files under ``data/interim``
and their rows in ``data/ledger.csv``. The raw source data are not distributed
with this repository. This runner is deterministic and can be rerun safely.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(ROOT))

from sparta.io_ import Paths, admitted_slides, load_config  # noqa: E402


TEST_FILES = (
    "test_barrier.py",
    "test_counterfactual.py",
    "test_loaders.py",
    "test_statistics.py",
    "test_new_experiments.py",
)
SENSITIVITY_SLIDES = ("MEL01", "MEL03", "CSCC01", "CSCC03", "CSCC05")


def _pinned_versions() -> None:
    requirements = ROOT / "requirements.txt"
    for line in requirements.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, expected = line.split("==", maxsplit=1)
        actual = importlib.metadata.version(name)
        if actual != expected:
            raise SystemExit(
                f"Environment mismatch: {name}=={actual}; expected {name}=={expected}. "
                "Install the pinned requirements before reproducing results."
            )


def _command(script: str, *args: str) -> list[str]:
    return [sys.executable, str(SCRIPTS / script), *args]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild SPARTA manuscript results, figures, and the Word draft."
    )
    parser.add_argument("--dry-run", action="store_true", help="print the run plan only")
    parser.add_argument(
        "--include-review-analyses", action="store_true",
        help="rerun the optional manuscript checks in run_22 through run_26",
    )
    args = parser.parse_args()

    cfg = load_config()
    paths = Paths(cfg)
    slides = list(admitted_slides(paths))
    if not slides:
        raise SystemExit("No admitted sections found in data/ledger.csv.")

    missing = [str(paths.interim / f"{sid}.raw.h5ad") for sid in slides
               if not (paths.interim / f"{sid}.raw.h5ad").is_file()]
    if missing:
        raise SystemExit(
            "Ingest the public data before running the analysis; missing inputs:\n  "
            + "\n  ".join(missing)
        )

    if not args.dry_run:
        _pinned_versions()

    cscc = [sid for sid in slides if sid.startswith("CSCC")]
    melanoma = [sid for sid in slides if sid.startswith("MEL")]
    sensitivity = [sid for sid in SENSITIVITY_SLIDES if sid in slides]
    if len(sensitivity) < 3 or not any(s.startswith("CSCC") for s in sensitivity) \
            or not any(s.startswith("MEL") for s in sensitivity):
        raise SystemExit(
            "The published sensitivity grid requires the five-section selection "
            "MEL01 MEL03 CSCC01 CSCC03 CSCC05."
        )

    plan: list[list[str]] = []
    plan.extend([[sys.executable, str(ROOT / "tests" / test)] for test in TEST_FILES])
    plan.extend([
        _command("run_00c_admission_audit.py", "--slides", *slides),
        _command("run_batch.py", "--through", "05", "--slides", *slides,
                 "--timeout", "14400"),
        _command("run_06_validate.py", "--slides", *slides, "--dim", "decoupling"),
    ])
    if cscc:
        plan.append(_command("run_06_validate.py", "--slides", *cscc,
                             "--dim", "consistency", "--tag", "CSCC"))
    if melanoma:
        plan.append(_command("run_06_validate.py", "--slides", *melanoma,
                             "--dim", "consistency", "--tag", "MEL"))
    plan.extend([
        _command("run_07_screen.py", "--slides", *slides, "--no-plot"),
        _command("run_08_benchmark.py"),
        _command("run_09_benchmark_real.py", *melanoma),
        _command("run_10_benchmark_ext.py", *slides),
        _command("run_11_review_diagnostics.py", "--slides", *slides),
        _command("run_12_paper_stats.py"),
        _command("run_13_benchmark_tools.py", "--slides", *slides, "--no-plot"),
        _command("run_13b_ndomains_scan.py", "--slides", *slides, "--no-plot"),
        _command("run_14_shared_ecm_check.py", "--slides", *slides),
        _command("run_17_sensitivity.py", "--slides", *sensitivity, "--no-plot"),
        _command("run_18_table1.py"),
        _command("run_19_runtime.py", "--slides", *slides),
        _command("run_21_mesh_stats.py"),
    ])
    if args.include_review_analyses:
        plan.extend([
            _command("run_22_null_crosslink.py"),
            _command("run_23_stromal_intervention.py"),
            _command("run_24_naive_baseline.py"),
            _command("run_25_biological_validation.py"),
            _command("run_26_radius_sensitivity.py"),
        ])
    plan.extend([
        _command("run_15_figure1.py", "--dpi", "300"),
        _command("run_16_figures.py", "--which", "2", "--dpi", "300"),
        _command("run_20_pending_figures.py"),
    ])
    node = shutil.which("node")
    if node is None:
        raise SystemExit("Node.js is required to regenerate docs/manuscript_cbc_draft.docx.")
    plan.append([node, str(SCRIPTS / "generate_docx.js")])

    print(f"SPARTA full rerun: {len(slides)} admitted sections")
    print("CSCC13 remains excluded by the ledger and is never added by this runner.")
    if args.include_review_analyses:
        print("Optional checks enabled: run_22 through run_26.")
    else:
        print("Optional checks skipped; use --include-review-analyses to rerun run_22 through run_26.")
    for i, cmd in enumerate(plan, start=1):
        print(f"\n[{i:02d}/{len(plan):02d}] {' '.join(cmd)}", flush=True)
        if not args.dry_run:
            subprocess.run(cmd, cwd=ROOT, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
