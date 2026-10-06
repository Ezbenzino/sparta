"""Build the SPARTA submission package, manifest, checks and ZIP archive."""
from __future__ import annotations

import csv
import hashlib
import shutil
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / "submission_IS"
OUT = OUT_ROOT / "SPARTA_submission_package"
FIG = ROOT / "results/figures/is"
SFIG = FIG / "supplement"

MAIN_FIG_STEMS = [
    "Fig1_overview", "Fig2_synthetic", "Fig3_calibration", "Fig4_association",
    "Fig5_structure", "Fig6_parameters", "Fig7_codex", "Fig8_baselines",
]
SUPPLEMENT_FIG_STEMS = [
    "FigS1_parameter_sensitivity", "FigS2a_maps", "FigS2b_maps",
    "FigS3_calibration", "FigS4_domain_scans", "FigS5_radius",
    "FigS6_reproduction", "FigS7_runtime", "FigS8_intervention",
    "FigS9_batch_effect",
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_file(src: Path, dst: Path, required: bool = True) -> bool:
    if not src.exists():
        if required:
            raise FileNotFoundError(src)
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def copy_tree(src: Path, dst: Path, ignore_names: set[str] | None = None,
              ignore_suffixes: set[str] | None = None):
    ignore_names = ignore_names or set()
    ignore_suffixes = ignore_suffixes or set()
    for path in src.rglob("*"):
        if path.is_dir():
            continue
        rel = path.relative_to(src)
        if any(part in ignore_names for part in rel.parts):
            continue
        if path.name in ignore_names or path.suffix in ignore_suffixes:
            continue
        target = dst / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)


def collect_manuscript():
    manuscript = OUT / "manuscript"
    files = [
        ("manuscript_IS.md", ROOT / "docs_is/manuscript_IS.md"),
        ("assembled_manuscript.md", ROOT / "docs_is/build/assembled_after_extension.md"),
        ("online_resource_1.md", ROOT / "docs_is/online_resource_1.md"),
        ("references.md", ROOT / "docs_is/references.md"),
        ("cover_letter.md", ROOT / "docs_is/cover_letter.md"),
        ("title_page.md", ROOT / "docs_is/title_page.md"),
        ("suggested_reviewers.md", ROOT / "docs_is/suggested_reviewers.md"),
    ]
    for name, src in files:
        copy_file(src, manuscript / name, required=not name.startswith(("assembled", "suggested")))


def collect_tables():
    tables = OUT / "tables"
    for name in ["table1_sections.md", "table1_sections.csv", "table1_sections.xlsx"]:
        copy_file(ROOT / "docs" / name, tables / name)
    for name in ["codex_cores.csv", "graph_coverage_table.csv"]:
        copy_file(ROOT / "results/validation" / name, tables / name, required=False)


def collect_figures():
    figures = OUT / "figures"
    for stem in MAIN_FIG_STEMS:
        for suffix in [".png", ".pdf", ".tif", ".eps"]:
            copy_file(FIG / f"{stem}{suffix}", figures / "main" / f"{stem}{suffix}")
    for stem in SUPPLEMENT_FIG_STEMS:
        for suffix in [".png", ".pdf", ".tif"]:
            copy_file(SFIG / f"{stem}{suffix}", figures / "supplement" / f"{stem}{suffix}",
                      required=False)


def collect_code():
    code = OUT / "code"
    for name in ["README.md", "README.zh.md", "LICENSE", "requirements.txt",
                 "environment.yml", "pyproject.toml", ".gitignore"]:
        copy_file(ROOT / name, code / name, required=name not in {".gitignore"})
    copy_tree(ROOT / "sparta", code / "sparta",
              ignore_names={"__pycache__"},
              ignore_suffixes={".pyc", ".pyo"})
    copy_tree(ROOT / "scripts", code / "scripts",
              ignore_names={"__pycache__"},
              ignore_suffixes={".pyc", ".pyo"})
    copy_tree(ROOT / "configs", code / "configs")
    copy_tree(ROOT / "tests", code / "tests",
              ignore_names={"__pycache__", ".pytest_cache"},
              ignore_suffixes={".pyc", ".pyo"})


def collect_reports():
    reports = OUT / "reports"
    files = [
        (ROOT / "docs/melanoma_mapping_correction.md", "melanoma_mapping_correction.md"),
        (ROOT / "docs/response_to_reviewers_2026.md", "response_to_reviewers_2026.md"),
        (ROOT / "results/qc/extension_2026_qc_summary.csv", "extension_qc_summary.csv"),
        (ROOT / "results/qc/extension_2026_qc_summary.json", "extension_qc_summary.json"),
        (ROOT / "data/external/extension_2026/extension_registry.csv", "extension_registry.csv"),
        (ROOT / "data/external/extension_2026/README.md", "extension_README.md"),
        (ROOT / "results/logs/extension_2026_ingest.log", "extension_ingest.log"),
        (ROOT / "results/logs/clinical_term_scan.txt", "clinical_term_scan.txt"),
        (ROOT / "results/logs/reproducibility_log.txt", "reproducibility_log.txt"),
    ]
    for src, name in files:
        copy_file(src, reports / name, required=False)


def self_check() -> list[str]:
    checks = []
    for stem in MAIN_FIG_STEMS:
        for suffix in [".png", ".pdf", ".tif", ".eps"]:
            p = OUT / "figures/main" / f"{stem}{suffix}"
            checks.append(("PASS" if p.exists() else "FAIL", f"main figure {p.name}"))
    for stem in SUPPLEMENT_FIG_STEMS:
        for suffix in [".png", ".pdf"]:
            p = OUT / "figures/supplement" / f"{stem}{suffix}"
            checks.append(("PASS" if p.exists() else "FAIL", f"supplement figure {p.name}"))
    for rel in ["manuscript/manuscript_IS.md", "manuscript/assembled_manuscript.md",
                "code/sparta", "code/scripts", "code/README.md",
                "reports/extension_qc_summary.csv", "reports/reproducibility_log.txt"]:
        p = OUT / rel
        checks.append(("PASS" if p.exists() else "FAIL", rel))
    forbidden = list((OUT / "code").rglob("exploratory")) + list((OUT / "code").rglob("legacy"))
    checks.append(("PASS" if not forbidden else "FAIL", "no exploratory/legacy in clean code"))
    return [f"{status}: {name}" for status, name in checks]


def write_manifest_and_summary(checks: list[str]):
    rows = []
    for path in sorted(OUT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(OUT).as_posix()
        rows.append({"relative_path": rel, "size_bytes": path.stat().st_size, "md5": md5(path)})
    with (OUT / "MANIFEST.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["relative_path", "size_bytes", "md5"])
        writer.writeheader()
        writer.writerows(rows)

    failures = [line for line in checks if line.startswith("FAIL")]
    # MANIFEST cannot contain a final hash of itself; it lists all other content files.
    total_files = sum(1 for p in OUT.rglob("*") if p.is_file()) + 1  # adds PACKAGE_SUMMARY below
    summary = [
        "# SPARTA submission package summary",
        "",
        f"Built: {datetime.now().astimezone().isoformat()}",
        f"Files: {total_files}",
        f"Manifested content files: {len(rows)}",
        f"Self-check status: {'PASS' if not failures else 'FAIL'}",
        "",
        "## Self-checks",
        *checks,
        "",
        "## Known rendering limitation",
        "Pandoc, pdflatex and R were not available in PATH when this package was built; therefore rendered DOCX/PDF outputs are not included. The package contains source Markdown, assembled Markdown, vector/high-resolution figures, tables, clean code, QC reports and reproducibility logs.",
    ]
    (OUT / "PACKAGE_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    return failures


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    collect_manuscript()
    collect_tables()
    collect_figures()
    collect_code()
    collect_reports()
    checks = self_check()
    failures = write_manifest_and_summary(checks)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    zip_base = OUT_ROOT / "SPARTA_submission_IS"
    if zip_base.with_suffix(".zip").exists():
        zip_base.with_suffix(".zip").unlink()
    archive = shutil.make_archive(str(zip_base), "zip", root_dir=OUT_ROOT,
                                 base_dir=OUT.name)
    snapshot_dir = OUT_ROOT / "snapshots"
    snapshot_dir.mkdir(exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(archive, snapshot_dir / f"SPARTA_submission_IS_{timestamp}.zip")
    print(f"package: {OUT}")
    print(f"zip: {archive}")
    print(f"checks: {sum('PASS' in c for c in checks)} pass, {len(failures)} fail")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
