#!/usr/bin/env python
"""
run_55_extension_ingest.py
==========================

Ingest and quality-control the 2026 SPARTA extension spatial-transcriptomics
datasets.

Inputs:
    data/external/extension_2026/GSE*/...
Outputs:
    data/interim/extension_2026/{slide}.h5ad
    data/external/extension_2026/extension_registry.csv
    data/external/extension_2026/README.md
    results/qc/extension_2026_qc_summary.csv
    results/qc/extension_2026_qc_summary.json
    results/logs/extension_2026_ingest.log
    data/ledger.csv (extension rows appended/updated)

The script is idempotent: existing extension ledger rows are replaced and all
derived files are regenerated.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform as platform_mod
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
EXT_ROOT = ROOT / "data" / "external" / "extension_2026"
INTERIM_DIR = ROOT / "data" / "interim" / "extension_2026"
QC_DIR = ROOT / "results" / "qc"
LOG_DIR = ROOT / "results" / "logs"
LEDGER_PATH = ROOT / "data" / "ledger.csv"
REGISTRY_PATH = EXT_ROOT / "extension_registry.csv"
README_PATH = EXT_ROOT / "README.md"

VISIUM_MIN_COUNTS = 500
SLIDESEQ_MIN_COUNTS = 100
MIN_CELLS_PER_GENE = 3
MAX_MITO_FRACTION = 0.20
VISIUM_MIN_SECTION_SPOTS = 1000
VISIUM_MIN_MEDIAN_UMI = 1500
# Extension metastatic deposits can be physically smaller than primary sections.
# Use a separate, explicitly reported section-size threshold and later run a
# sensitivity analysis excluding sections below the primary-cohort 1000-spot rule.
EXT_VISIUM_MIN_SECTION_SPOTS = 300
EXT_VISIUM_MIN_RETAINED_MEDIAN_UMI = 1500
SLIDESEQ_MIN_SECTION_BEADS = 1000
SLIDESEQ_MIN_MEDIAN_UMI = 100


@dataclass
class SampleSpec:
    """Static identity metadata for one extension section."""

    slide_id: str
    patient_id: str
    project_id: str
    accession: str
    platform: str
    cancer_type: str
    subtype: str
    progression: str
    site: str
    treatment: str
    replicate: str
    count_path: str
    spatial_path: str | None


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-h5ad", action="store_true",
                    help="run loading/QC and write tables without saving h5ad files")
    return ap.parse_args()


def label_from_count(project_id: str, path: Path) -> str:
    """Extract the submitter sample label from a downloaded count file name."""
    accession = path.name.split("_", 1)[0]
    rest = path.name[len(accession) + 1:]
    for suffix in (
        "_matrix.mtx.gz",
        "_filtered_feature_bc_matrix.h5",
    ):
        if rest.endswith(suffix):
            label = rest[: -len(suffix)]
            break
    else:
        raise ValueError(f"Cannot infer label from {path.name}")
    if project_id in {"GSE289745", "GSE300445"}:
        label = label.removesuffix("_SpT").removesuffix("_processed")
    return label


def make_standard_spec(project_id: str, count_path: Path) -> SampleSpec:
    """Build a sample spec for Visium datasets."""
    accession = count_path.name.split("_", 1)[0]
    label = label_from_count(project_id, count_path)
    spatial_path = None
    if project_id == "GSE289745":
        slide_id = f"CSCC289_{label}"
        spec = SampleSpec(
            slide_id=slide_id, patient_id=slide_id, project_id=project_id,
            accession=accession, platform="Visium", cancer_type="cscc",
            subtype="Cutaneous squamous cell carcinoma",
            progression="cutaneous cSCC (progression not otherwise specified)",
            site="skin", treatment="unknown", replicate=label,
            count_path=str(count_path), spatial_path=spatial_path,
        )
    elif project_id == "GSE300445":
        number = label.split("-", 1)[1]
        slide_id = f"MEL300_{number}"
        spec = SampleSpec(
            slide_id=slide_id, patient_id=slide_id, project_id=project_id,
            accession=accession, platform="Visium", cancer_type="melanoma",
            subtype="Primary cutaneous melanoma", progression="primary",
            site="skin", treatment="pre-adjuvant anti-PD-1",
            replicate=label, count_path=str(count_path), spatial_path=spatial_path,
        )
    elif project_id == "GSE316760":
        slide_id = f"MEL316_{label}"
        spec = SampleSpec(
            slide_id=slide_id, patient_id=slide_id, project_id=project_id,
            accession=accession, platform="Visium", cancer_type="melanoma",
            subtype="Primary cutaneous melanoma", progression="primary",
            site="skin", treatment="unknown", replicate=label,
            count_path=str(count_path), spatial_path=spatial_path,
        )
    elif project_id == "GSE320041":
        patient_label = re.sub(r"_(\d+)$", "", label)
        slide_id = f"MEL320_{label}"
        patient_id = f"MEL320_{patient_label}"
        spec = SampleSpec(
            slide_id=slide_id, patient_id=patient_id, project_id=project_id,
            accession=accession, platform="Visium", cancer_type="melanoma",
            subtype="Metastatic melanoma", progression="metastatic",
            site="metastatic melanoma (specific site not annotated)",
            treatment="unknown", replicate=label,
            count_path=str(count_path), spatial_path=spatial_path,
        )
    elif project_id == "GSE321832":
        slide_id = f"CSCC321_{label}"
        spec = SampleSpec(
            slide_id=slide_id, patient_id=slide_id, project_id=project_id,
            accession=accession, platform="Visium", cancer_type="cscc",
            subtype="Primary cutaneous squamous cell carcinoma",
            progression="primary", site="skin", treatment="unknown",
            replicate=label, count_path=str(count_path), spatial_path=spatial_path,
        )
    else:
        raise ValueError(f"Unsupported standard project {project_id}")
    return spec


def make_slideseq_spec(count_path: Path) -> SampleSpec:
    """Build a sample spec for a Slide-seqV2 section."""
    pattern = re.compile(
        r"^(?P<accession>GSM\d+)_(?P<label>.+)_slide_raw_counts\.csv\.gz$"
    )
    match = pattern.match(count_path.name)
    if not match:
        raise ValueError(f"Cannot parse Slide-seq file {count_path.name}")
    accession = match.group("accession")
    label = match.group("label")
    spatial_path = Path(
        str(count_path).replace("raw_counts.csv.gz", "spatial_info.csv.gz")
    )
    base = label.split("_rep", 1)[0]
    if base.startswith("MBM"):
        patient_id = base
        subtype = "Melanoma brain metastasis"
        site = "brain"
    elif base.startswith("ECM"):
        patient_id = base.replace("ECM", "MPM", 1)
        subtype = "Melanoma subcutaneous metastasis"
        site = "subcutaneous tissue"
    else:
        raise ValueError(f"Unknown Slide-seq label {label}")
    return SampleSpec(
        slide_id=label, patient_id=patient_id, project_id="GSE200278",
        accession=accession, platform="Slide-seqV2", cancer_type="melanoma",
        subtype=subtype, progression="metastatic", site=site,
        treatment="treatment-naive metastasis", replicate=label,
        count_path=str(count_path), spatial_path=str(spatial_path),
    )


def discover_samples() -> list[SampleSpec]:
    """Discover all 50 included extension sections."""
    specs: list[SampleSpec] = []
    for project_dir in sorted(path for path in EXT_ROOT.iterdir() if path.is_dir()):
        project_id = project_dir.name
        if project_id == "GSE200278":
            counts = sorted(project_dir.glob("*raw_counts.csv.gz"))
            specs.extend(make_slideseq_spec(path) for path in counts)
            continue
        counts = sorted([
            *project_dir.glob("*matrix.mtx.gz"),
            *project_dir.glob("*filtered_feature_bc_matrix.h5"),
        ])
        specs.extend(make_standard_spec(project_id, path) for path in counts)
    if len(specs) != 50:
        raise AssertionError(f"Expected 50 extension sections, found {len(specs)}")
    return sorted(specs, key=lambda spec: (spec.project_id, spec.slide_id))


def load_sample(spec: SampleSpec):
    """Load a sample using its platform-specific loader."""
    from sparta.extension_loaders import load_visium, read_slideseqv2

    if spec.platform == "Visium":
        return load_visium(spec.count_path, spec.slide_id, spec.project_id)
    if spec.platform == "Slide-seqV2":
        assert spec.spatial_path is not None
        return read_slideseqv2(
            spec.count_path, spec.spatial_path, spec.slide_id, spec.project_id
        )
    raise ValueError(f"Unsupported platform {spec.platform}")


def qc_sample(slide) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    """Run QC and return summary plus spot/gene keep arrays."""
    counts = slide.counts
    total = np.asarray(counts.sum(axis=0)).ravel()
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    mt_mask = np.char.startswith(slide.gene_names.astype(str), "MT-")
    mt_assessable = bool(mt_mask.any())
    if mt_assessable:
        mt_sum = np.asarray(counts[mt_mask].sum(axis=0)).ravel()
        mt_fraction = np.divide(
            mt_sum, total, out=np.zeros_like(total, dtype=float), where=total > 0
        )
    else:
        mt_fraction = np.full(counts.shape[1], np.nan, dtype=float)
    min_counts = (
        SLIDESEQ_MIN_COUNTS if slide.platform == "Slide-seqV2" else VISIUM_MIN_COUNTS
    )
    keep_umi = total >= min_counts
    keep_mito = (
        mt_fraction <= MAX_MITO_FRACTION if mt_assessable else np.ones(counts.shape[1], bool)
    )
    spot_keep = keep_umi & keep_mito
    counts_after_spots = counts[:, spot_keep]
    gene_present = np.asarray((counts_after_spots > 0).sum(axis=1)).ravel()
    gene_keep = gene_present >= MIN_CELLS_PER_GENE

    retained_total = total[spot_keep]
    retained_detected = detected[spot_keep]
    finite_coords = bool(np.isfinite(slide.coords_um).all())
    unique_barcodes = len(set(slide.barcodes)) == len(slide.barcodes)
    if slide.platform == "Slide-seqV2":
        enough_spots = spot_keep.sum() >= SLIDESEQ_MIN_SECTION_BEADS
        enough_umi = (
            np.median(retained_total) >= SLIDESEQ_MIN_MEDIAN_UMI
            if spot_keep.any() else False
        )
    else:
        enough_spots = spot_keep.sum() >= EXT_VISIUM_MIN_SECTION_SPOTS
        enough_umi = (
            np.median(retained_total) >= EXT_VISIUM_MIN_RETAINED_MEDIAN_UMI
            if spot_keep.any() else False
        )
    checks = {
        "finite_coordinates": finite_coords,
        "unique_barcodes": unique_barcodes,
        "enough_spots": bool(enough_spots),
        "enough_median_umi": bool(enough_umi),
        "matrix_dimensions_match": counts.shape[1] == len(slide.barcodes)
        and counts.shape[0] == len(slide.gene_names),
    }
    reasons = [name for name, passed in checks.items() if not passed]
    summary = {
        "platform": slide.platform,
        "raw_spots": int(counts.shape[1]),
        "raw_genes": int(counts.shape[0]),
        "retained_spots": int(spot_keep.sum()),
        "retained_genes": int(gene_keep.sum()),
        "raw_median_umi": float(np.median(total)),
        "retained_median_umi": float(np.median(retained_total)) if spot_keep.any() else np.nan,
        "raw_median_detected_genes": float(np.median(detected)),
        "retained_median_detected_genes": float(np.median(retained_detected))
        if spot_keep.sum() else np.nan,
        "spots_removed_umi": int((~keep_umi).sum()),
        "spots_removed_mito": int((~keep_mito).sum()) if mt_assessable else 0,
        "genes_removed_low_detection": int((~gene_keep).sum()),
        "mitochondria_assessable": mt_assessable,
        "median_mito_fraction": float(np.nanmedian(mt_fraction)) if mt_assessable else None,
        "max_mito_fraction": float(np.nanmax(mt_fraction)) if mt_assessable else None,
        "min_counts_cutoff": int(min_counts),
        "max_mito_fraction_cutoff": MAX_MITO_FRACTION if mt_assessable else None,
        "min_cells_per_gene": MIN_CELLS_PER_GENE,
        "qc_passed": not reasons,
        "qc_failure_reasons": "; ".join(reasons),
    }
    return summary, spot_keep, gene_keep


def subset_slide(slide, spot_keep: np.ndarray, gene_keep: np.ndarray):
    """Return a QC-filtered SpatialExpressionData object."""
    counts = slide.counts[gene_keep, :][:, spot_keep].tocsr()
    coords_pixel = slide.coords_pixel[spot_keep] if slide.coords_pixel is not None else None
    metadata = dict(slide.metadata)
    metadata.update({
        "qc_filtered": True,
        "qc_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    })
    from sparta.extension_loaders import SpatialExpressionData

    return SpatialExpressionData(
        slide_id=slide.slide_id,
        project_id=slide.project_id,
        platform=slide.platform,
        counts=counts,
        gene_ids=slide.gene_ids[gene_keep],
        gene_names=slide.gene_names[gene_keep],
        barcodes=slide.barcodes[spot_keep],
        coords_um=slide.coords_um[spot_keep],
        coords_pixel=coords_pixel,
        in_tissue=slide.in_tissue[spot_keep],
        scalefactors=slide.scalefactors,
        metadata=metadata,
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a list of dictionaries as CSV, preserving all keys."""
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def update_ledger(specs: list[SampleSpec], qc_rows: dict[str, dict[str, Any]]) -> None:
    """Replace extension rows in data/ledger.csv, preserving primary/replication rows."""
    with LEDGER_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        existing = list(reader)
    if not fieldnames:
        raise AssertionError("Ledger has no fieldnames")
    retained = [
        row for row in existing
        if row.get("status") not in {"extension", "extension_rejected"}
    ]
    new_rows: list[dict[str, str]] = []
    for spec in specs:
        qc = qc_rows[spec.slide_id]
        status = "extension" if qc["qc_passed"] else "extension_rejected"
        row = {key: "" for key in fieldnames}
        row.update({
            "slide_id": spec.slide_id,
            "patient": spec.patient_id,
            "replicate": spec.replicate,
            "cancer_type": spec.cancer_type,
            "platform": "slideseqv2" if spec.platform == "Slide-seqV2" else "visium",
            "source": f"GEO {spec.project_id}",
            "accession": spec.accession,
            "raw_path": f"data/external/extension_2026/{spec.project_id}",
            "format": "slideseqv2" if spec.platform == "Slide-seqV2" else "visium",
            "coord_source": "file",
            "has_counts": "True",
            "has_coords": "True",
            "has_image": "False",
            "n_spots": str(qc["retained_spots"]),
            "n_genes": str(qc["retained_genes"]),
            "median_umi": f"{qc['retained_median_umi']:.1f}",
            "treatment": spec.treatment,
            "site": spec.site,
            "status": status,
            "notes": (f"{spec.subtype}; {spec.progression}; QC min counts "
                      f"{qc['min_counts_cutoff']}, gene >= {MIN_CELLS_PER_GENE} spots"),
        })
        new_rows.append(row)
    with LEDGER_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(retained)
        writer.writerows(new_rows)


def cohort_counts(specs: list[SampleSpec]) -> dict[str, dict[str, int]]:
    """Count sections and patients by project."""
    out: dict[str, dict[str, int]] = {}
    for spec in specs:
        item = out.setdefault(spec.project_id, {"sections": 0, "patients": 0})
        item["sections"] += 1
    for project_id, item in out.items():
        item["patients"] = len({
            spec.patient_id for spec in specs if spec.project_id == project_id
        })
    return out


def write_readme(specs: list[SampleSpec], qc_rows: dict[str, dict[str, Any]]) -> None:
    """Write the extension dataset README."""
    counts = cohort_counts(specs)
    disease_summary: dict[str, int] = {}
    for spec in specs:
        disease_summary[spec.subtype] = disease_summary.get(spec.subtype, 0) + 1
    lines = [
        "# SPARTA 2026 extension datasets",
        "",
        "## Scope",
        "",
        "These public spatial-transcriptomics datasets extend the method-oriented",
        "SPARTA analysis across cutaneous SCC, primary cutaneous melanoma and",
        "metastatic melanoma. The ingestion uses expression matrices and spatial",
        "coordinates only. Treatment-response, prognosis and other clinical labels",
        "present in the source publications are not used.",
        "",
        "## Included datasets",
        "",
        "| Project | Platform | Sections | Patients | Main disease strata |",
        "|---|---|---:|---:|---|",
    ]
    for project_id in sorted(counts):
        strata = sorted({
            spec.subtype for spec in specs if spec.project_id == project_id
        })
        item = counts[project_id]
        lines.append(
            f"| {project_id} | "
            f"{'Slide-seqV2' if project_id == 'GSE200278' else 'Visium'} | "
            f"{item['sections']} | {item['patients']} | {'; '.join(strata)} |"
        )
    lines += [
        "",
        "## Disease-strand totals",
        "",
    ]
    for subtype, count in sorted(disease_summary.items()):
        lines.append(f"- {subtype}: {count} sections")
    lines += [
        "",
        "## QC rules",
        "",
        f"- Visium spots: total UMI >= {VISIUM_MIN_COUNTS}; extension sections require "
        f"at least {EXT_VISIUM_MIN_SECTION_SPOTS} retained spots and retained median UMI "
        f">= {EXT_VISIUM_MIN_RETAINED_MEDIAN_UMI}. The primary-cohort rule used "
        f"{VISIUM_MIN_SECTION_SPOTS} spots; sections below that threshold are retained "
        "here as small metastatic deposits and separately examined in a sensitivity analysis.",
        f"- Slide-seqV2 beads: total UMI >= {SLIDESEQ_MIN_COUNTS}; section-level minimum "
        f"{SLIDESEQ_MIN_SECTION_BEADS} retained beads and retained median UMI >= {SLIDESEQ_MIN_MEDIAN_UMI}.",
        f"- Genes: detected in at least {MIN_CELLS_PER_GENE} retained spots/beads.",
        "- Slide-seqV2 coordinates are treated as micrometres based on the 10-um bead scale; beads sharing one published coordinate are aggregated into one spatial point.",
        f"- Mitochondria: if mitochondrial genes are present, fraction <= {MAX_MITO_FRACTION:.2f}. "
        "Some processed Visium matrices exclude mitochondrial genes; this is recorded as not assessable.",
        "",
        "## Exclusions and unavailable samples",
        "",
        "- GSE320041: WU1457 and WU2109 are Visium HD 16-um samples and are not included in this 50-section batch.",
        "- GSE321832: only the two cutaneous SCC sections are included; two head-and-neck and one lung SCC sections are outside this skin-focused batch.",
        "- GSE300445: only the four Visium primary tumours are included; bulk RNA-seq and other study components are not used.",
        "- GSE200278: only the 16 sections with matched Slide-seqV2 data are included; samples without spatial sections are not used.",
        "",
        "## Reproduction",
        "",
        "```powershell",
        "& '.venv\\Scripts\\python.exe' scripts/run_55_extension_ingest.py",
        "```",
        "",
        "QC-filtered AnnData files are written to `data/interim/extension_2026/`.",
        "The detailed registry is `extension_registry.csv`; per-section QC summaries",
        "are in `results/qc/extension_2026_qc_summary.csv`.",
        "",
    ]
    README_PATH.write_text("\n".join(lines), encoding="utf-8")


def environment_versions() -> dict[str, str]:
    """Return key software versions for the ingest log."""
    import scipy

    try:
        from importlib.metadata import version
        anndata_version = version("anndata")
    except Exception as error:  # noqa: BLE001
        anndata_version = f"unavailable ({error})"
    return {
        "python": sys.version.replace("\n", " "),
        "executable": sys.executable,
        "platform": platform_mod.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "h5py": h5py.__version__,
        "anndata": anndata_version,
    }


def main() -> None:
    """Run extension ingestion end to end."""
    args = parse_args()
    t0 = time.time()
    for path in (INTERIM_DIR, QC_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)
    specs = discover_samples()
    qc_rows: dict[str, dict[str, Any]] = {}
    registry_rows: list[dict[str, Any]] = []
    log_lines = [
        f"Extension ingestion started {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        f"environment {json.dumps(environment_versions(), ensure_ascii=False)}",
    ]
    for i, spec in enumerate(specs, start=1):
        sample_t0 = time.time()
        print(f"[{i:02d}/{len(specs)}] {spec.slide_id} ({spec.platform})")
        slide = load_sample(spec)
        qc, spot_keep, gene_keep = qc_sample(slide)
        qc_rows[spec.slide_id] = qc
        log_lines.append(
            f"{spec.slide_id}: raw {slide.counts.shape}, retained "
            f"{qc['retained_spots']} spots/{qc['retained_genes']} genes, "
            f"passed={qc['qc_passed']}, seconds={time.time() - sample_t0:.1f}"
        )
        if qc["qc_passed"] and not args.no_h5ad:
            filtered = subset_slide(slide, spot_keep, gene_keep)
            if "group_barcodes" in filtered.metadata:
                filtered.metadata["group_barcodes_json"] = json.dumps(
                    filtered.metadata.pop("group_barcodes"), ensure_ascii=False
                )
            adata = filtered.to_anndata()
            h5ad_path = INTERIM_DIR / f"{spec.slide_id}.h5ad"
            adata.write_h5ad(h5ad_path, compression="gzip", compression_opts=6)
        registry = asdict(spec)
        registry.update(qc)
        registry_rows.append(registry)

    qc_csv = QC_DIR / "extension_2026_qc_summary.csv"
    write_csv(qc_csv, [qc_rows[spec.slide_id] | {"slide_id": spec.slide_id}
                       for spec in specs])
    write_csv(REGISTRY_PATH, registry_rows)
    (QC_DIR / "extension_2026_qc_summary.json").write_text(
        json.dumps(qc_rows, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    update_ledger(specs, qc_rows)
    write_readme(specs, qc_rows)
    passed = sum(row["qc_passed"] for row in qc_rows.values())
    elapsed = time.time() - t0
    log_lines += [
        f"Discovered sections: {len(specs)}",
        f"QC passed: {passed}/{len(specs)}",
        f"Elapsed seconds: {elapsed:.1f}",
    ]
    (LOG_DIR / "extension_2026_ingest.log").write_text(
        "\n".join(log_lines) + "\n", encoding="utf-8"
    )
    print(f"Ingested {len(specs)} sections; QC passed {passed}/{len(specs)}")
    print(f"QC table: {qc_csv}")
    print(f"Ledger updated: {LEDGER_PATH}")
    print(f"README: {README_PATH}")


if __name__ == "__main__":
    main()
