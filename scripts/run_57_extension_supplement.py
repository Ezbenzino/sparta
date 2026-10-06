#!/usr/bin/env python
"""Build the 2026 extension cohort section of Online Resource 1."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "docs_is/online_resource_1.md"
REGISTRY = ROOT / "data/external/extension_2026/extension_registry.csv"
MARKER = "# S15 2026 extension cohort"


def fmt_int(value):
    return f"{int(float(value)):,}" if value not in (None, "") else "—"


def fmt_float(value, digits=1):
    return f"{float(value):.{digits}f}" if value not in (None, "") else "—"


def main():
    rows = list(csv.DictReader(REGISTRY.open(encoding="utf-8")))
    rows.sort(key=lambda r: (r["project_id"], r["slide_id"]))
    lines = [
        MARKER,
        "",
        "This section documents the 50 public sections added as the 2026 extension cohort. "
        "QC was performed at the native Visium spot or Slide-seqV2 bead level. Slide-seqV2 counts "
        "were subsequently summed into fixed 50-µm grid bins for structural analysis; the h5ad files "
        "retain the QC-passed bead-level counts.",
        "",
        "**Table S11** Extension cohort metadata. Patient relationship was not provided by GSE289745; "
        "those sections are marked as relationship unavailable and were treated as conservative "
        "section-level sampling units.",
        "",
        "| Project | Section | Patient identifier | Disease / stratum | Platform | Site |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        patient = row["patient_id"]
        if row["project_id"] == "GSE289745":
            patient = "relationship unavailable"
        disease = f"{row['cancer_type']} / {row['progression']}"
        lines.append(
            f"| {row['project_id']} | {row['slide_id']} | {patient} | {disease} | "
            f"{row['platform']} | {row.get('site') or '—'} |"
        )

    lines += [
        "",
        "**Table S12** Extension cohort QC. Raw and retained locations are Visium spots or Slide-seqV2 "
        "beads before the separate Slide-seqV2 50-µm grid aggregation. MT, mitochondrial UMI fraction; "
        "NA indicates that the supplied matrix did not contain mitochondrial genes.",
        "",
        "| Section | Platform | Raw locations | Retained locations | Retained genes | Median retained UMI | Median MT fraction | QC |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        mt = "NA" if row["mitochondria_assessable"] != "True" else fmt_float(
            row["median_mito_fraction"], 3)
        lines.append(
            f"| {row['slide_id']} | {row['platform']} | {fmt_int(row['raw_spots'])} | "
            f"{fmt_int(row['retained_spots'])} | {fmt_int(row['retained_genes'])} | "
            f"{fmt_float(row['retained_median_umi'], 0)} | {mt} | {row['qc_passed']} |"
        )
    lines.append("")

    text = SUPP.read_text(encoding="utf-8")
    if MARKER in text:
        text = text.split(MARKER, 1)[0].rstrip() + "\n\n"
    SUPP.write_text(text.rstrip() + "\n\n" + "\n".join(lines), encoding="utf-8")
    print(f"wrote extension supplement with {len(rows)} sections -> {SUPP}")


if __name__ == "__main__":
    main()
