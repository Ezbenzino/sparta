#!/usr/bin/env python
"""Package a local, DOI-ready CMPB code/results archive candidate.

The script does not publish to GitHub or Zenodo. It creates a ZIP from an
explicit allowlist, a SHA-256 manifest, and a human-readable note.
Raw expression data, intermediate arrays, local environments, and logs are
excluded; data sources remain discoverable through public accessions.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "release_candidate"
STAMP = "20261004"
ARCHIVE_NAME = f"SPARTA-CMPB-analysis-candidate-{STAMP}.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_info():
    def run(*args):
        try:
            return subprocess.check_output(["git", *args], cwd=ROOT,
                                           text=True, stderr=subprocess.DEVNULL).strip()
        except Exception:
            return None
    return {"commit": run("rev-parse", "HEAD"),
            "branch": run("branch", "--show-current"),
            "working_tree_status": run("status", "--short")}


def collect_files():
    files = set()
    top_files = ["README.md", "README.zh.md", "LICENSE", "CITATION.cff",
                 "pyproject.toml", "requirements.txt", "environment.yml"]
    for rel in top_files:
        p = ROOT / rel
        if p.is_file():
            files.add(p)
    for dirname in ("configs", "sparta", "tests"):
        base = ROOT / dirname
        files.update(p for p in base.rglob("*") if p.is_file()
                     and not any(part in {"__pycache__", ".pytest_cache"}
                                 for part in p.parts))
    scripts = ROOT / "scripts"
    files.update(p for p in scripts.glob("run_*.py") if p.is_file())
    files.update(p for p in (scripts / "generate_cmpb_figures.py",
                             scripts / "build_cmpb_docx.py",
                             scripts / "build_cmpb_release_candidate.py") if p.is_file())
    for rel in ("data/ledger.csv", "data/admission_audit.csv",
                "docs/manuscript_cmpb_draft.md",
                "docs/manuscript_cmpb_draft.docx",
                "docs/manuscript_cmpb_highlights.txt",
                "docs/figure_legends.md", "docs/release_checklist.md",
                "docs/review_for_journal.md", "docs/cmpb_readiness_review.md",
                "docs/cover_letter.md"):
        p = ROOT / rel
        if p.is_file():
            files.add(p)
    for dirname, pattern in (("data/external", "*.json"),
                             ("results/validation", "*.json"),
                             ("results/counterfactual", "*.json"),
                             ("results/figures/cmpb", "*")):
        base = ROOT / dirname
        if base.exists():
            files.update(p for p in base.rglob(pattern) if p.is_file())
    return sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    archive_path = OUT_DIR / ARCHIVE_NAME
    manifest_path = OUT_DIR / "manifest.json"
    info_path = OUT_DIR / "README_ARCHIVE.txt"
    files = collect_files()
    manifest = {
        "archive": ARCHIVE_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "local candidate; not deposited; no DOI",
        "git": git_info(),
        "scope": "Source, configs, analysis scripts, primary/external validation JSON, counterfactual JSON, CMPB figures and manuscript sources.",
        "runtime": {"python": sys.version.split()[0],
                    "platform": platform.platform(),
                    "packages": {}},
        "generation_commands": [
            "python scripts/run_27_spatial_null.py",
            "python scripts/run_28_s2_matched.py",
            "python scripts/run_29_ext_validation.py",
            "python scripts/run_30_s2_matched_ext.py",
            "python scripts/run_31_benchmark_ext_slides.py",
            "python scripts/run_32_prereadiness_audit.py",
            "python scripts/run_33_input_connectivity_sensitivity.py",
            "python scripts/run_34_benchmark_lambda_sensitivity.py",
            "python scripts/generate_cmpb_figures.py",
            "python scripts/build_cmpb_docx.py"],
        "random_seeds_and_configs": "See configs/default.yaml and the analysis scripts; per-analysis seeds and settings are recorded there or in output JSON metadata.",
        "excluded": ["raw GEO expression data", "intermediate arrays and AnnData files",
                     "local environments", "logs", "scratch scripts"],
        "reproduction_note": "Public source datasets must be downloaded from the accessions documented in the manuscript/README. Analysis outputs are included to support exact manuscript verification.",
        "files": [{"path": p.relative_to(ROOT).as_posix(),
                   "size_bytes": p.stat().st_size,
                   "sha256": sha256(p)} for p in files],
    }
    try:
        from importlib.metadata import PackageNotFoundError, version
        for package in ("numpy", "scipy", "pandas", "networkx", "scanpy",
                        "squidpy", "anndata", "matplotlib", "seaborn"):
            try:
                manifest["runtime"]["packages"][package] = version(package)
            except PackageNotFoundError:
                manifest["runtime"]["packages"][package] = None
    except Exception:
        manifest["runtime"]["packages"] = {"note": "Version inspection unavailable"}
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8")
    info_path.write_text(
        "SPARTA CMPB analysis archive candidate\n"
        "=======================================\n\n"
        "This is a local release candidate. It has not been deposited and has no DOI.\n"
        "Review the contents and exact manuscript/code state before publication.\n"
        "Raw public expression datasets and local intermediates are excluded.\n"
        "The SHA-256 manifest is included as manifest.json.\n",
        encoding="utf-8")
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=9) as zf:
        for p in files:
            zf.write(p, p.relative_to(ROOT).as_posix())
        zf.write(manifest_path, "manifest.json")
        zf.write(info_path, "README_ARCHIVE.txt")
    hash_path = OUT_DIR / f"{ARCHIVE_NAME}.sha256"
    hash_path.write_text(f"{sha256(archive_path)}  {ARCHIVE_NAME}\n",
                         encoding="ascii")
    print(f"Wrote {archive_path}")
    print(f"Files: {len(files)}; bytes: {archive_path.stat().st_size}")
    print(f"SHA-256: {sha256(archive_path)}")
    print(f"Git state: {manifest['git']}")


if __name__ == "__main__":
    main()
