#!/usr/bin/env python
"""
run_54_correct_melanoma_mapping.py
==================================

Correct the primary-cohort melanoma patient mapping in ``data/ledger.csv`` and
write a traceable correction report.

Input:
    data/ledger.csv
Output:
    data/ledger.csv (corrected in place; Git provides before/after traceability)
    docs/melanoma_mapping_correction.md

The correction is idempotent. Section-level expression values, QC values, graph
objects and per-section estimates are not changed; only patient-level aggregation
and pooled patient counts depend on this mapping.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "ledger.csv"
REPORT = ROOT / "docs" / "melanoma_mapping_correction.md"

# GSE250636 GEO sample characteristics:
# GSM7983359 (sternum) and GSM7983366 (ribcage) are patient A;
# GSM7983364 (cecal nodule) and GSM7983365 (chest wall) are patient B.
CORRECT_MAPPING = {
    "MEL01": {"accession": "GSM7983359", "patient": "MEL_PtA", "geo_patient": "A"},
    "MEL02": {"accession": "GSM7983364", "patient": "MEL_PtB", "geo_patient": "B"},
    "MEL03": {"accession": "GSM7983365", "patient": "MEL_PtB", "geo_patient": "B"},
    "MEL04": {"accession": "GSM7983366", "patient": "MEL_PtA", "geo_patient": "A"},
}
NOTE_MARKER = "2026-10-05 patient mapping corrected from GSE250636 sample characteristics"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments; heavy imports are intentionally avoided."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ledger", default=str(LEDGER), help="path to data/ledger.csv")
    ap.add_argument("--report", default=str(REPORT), help="correction report path")
    return ap.parse_args()


def correct_ledger(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Apply the corrected patient mapping and return rows and change records."""
    if not path.exists():
        sys.exit(f"Ledger not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not fieldnames:
        sys.exit("Ledger has no header")

    rows_by_id = {row["slide_id"]: row for row in rows}
    changes: list[dict[str, str]] = []
    for slide, spec in CORRECT_MAPPING.items():
        row = rows_by_id.get(slide)
        if row is None:
            sys.exit(f"Expected slide {slide} in ledger")
        if row.get("accession") != spec["accession"]:
            sys.exit(
                f"{slide}: expected accession {spec['accession']}, "
                f"found {row.get('accession')}"
            )
        old_patient = row.get("patient", "")
        new_patient = spec["patient"]
        if old_patient != new_patient:
            row["patient"] = new_patient
        notes = row.get("notes", "")
        correction_note = (
            f"{NOTE_MARKER}: GEO patient {spec['geo_patient']} -> {new_patient}"
        )
        if NOTE_MARKER not in notes:
            row["notes"] = f"{notes} | {correction_note}" if notes else correction_note
        changes.append({
            "slide_id": slide,
            "accession": spec["accession"],
            "old_patient": old_patient,
            "new_patient": new_patient,
            "geo_patient": spec["geo_patient"],
        })

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows, changes


def write_report(path: Path, changes: list[dict[str, str]]) -> None:
    """Write a human-readable correction report for editors and reviewers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Melanoma patient-mapping correction",
        "",
        "## Summary",
        "",
        "The initial primary-cohort ledger assigned all four extracranial melanoma",
        "sections to one patient (`MEL_PtB`). The GSE250636 sample characteristics",
        "identify two patients among those sections: patient A contributed MEL01 and",
        "MEL04, and patient B contributed MEL02 and MEL03.",
        "",
        "The primary cohort therefore contains **19 sections from 8 patients**, not 7.",
        "Section-level measurements did not change; patient-level summaries, nested",
        "models, confidence intervals, sign-flip tests and pooled analyses with the",
        "Thrane replication cohort were recomputed.",
        "",
        "## Evidence and corrected mapping",
        "",
        "| Slide | GEO sample | GEO patient ID | Old ledger patient | Corrected patient |",
        "|---|---|---|---|---|",
    ]
    for change in changes:
        lines.append(
            f"| {change['slide_id']} | {change['accession']} | "
            f"{change['geo_patient']} | {change['old_patient']} | "
            f"{change['new_patient']} |"
        )
    lines += [
        "",
        "## Source",
        "",
        "- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250636",
        "- Sample characteristics field: `patient id`",
        "",
        "## Reproduction",
        "",
        "1. Run `python scripts/run_54_correct_melanoma_mapping.py`.",
        "2. Re-run patient-level and table-generation scripts:",
        "   `python scripts/run_40_graph_coverage_table.py`;",
        "   `python scripts/run_36_patient_level.py`;",
        "   `python scripts/run_48_replication_cohort.py --stats-only`;",
        "   `python scripts/run_18_table1.py`;",
        "   `python scripts/is_figures/facts.py`.",
        "",
        "Git supplies the exact before/after version of the ledger.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Correct the ledger and print a concise completion message."""
    args = parse_args()
    ledger = Path(args.ledger)
    report = Path(args.report)
    _rows, changes = correct_ledger(ledger)
    write_report(report, changes)
    for change in changes:
        print(
            f"{change['slide_id']}: {change['old_patient']} -> "
            f"{change['new_patient']} (GEO {change['geo_patient']})"
        )
    print(f"Wrote corrected ledger: {ledger}")
    print(f"Wrote correction report: {report}")


if __name__ == "__main__":
    main()
