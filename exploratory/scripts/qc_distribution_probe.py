"""Probe QC distributions for representative extension samples."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sparta.extension_loaders import load_visium, read_slideseqv2  # noqa: E402

EXT = ROOT / "data" / "external" / "extension_2026"
SAMPLES = [
    ("GSE289745", EXT / "GSE289745" / "GSM8797975_S11_SpT_matrix.mtx.gz"),
    ("GSE316760", EXT / "GSE316760" / "GSM9459774_mel2_filtered_feature_bc_matrix.h5"),
    ("GSE200278", EXT / "GSE200278" / "GSM6025937_MBM05_rep3_slide_raw_counts.csv.gz"),
]


def summarize(slide) -> None:
    X = slide.counts  # genes x spots
    total = np.asarray(X.sum(axis=0)).ravel()
    detected = (X > 0).sum(axis=0)
    detected = np.asarray(detected).ravel()
    mt = np.char.startswith(slide.gene_names.astype(str), "MT-")
    mt_sum = np.asarray(X[mt].sum(axis=0)).ravel() if mt.any() else np.zeros(X.shape[1])
    mt_frac = np.divide(mt_sum, total, out=np.zeros_like(mt_sum, dtype=float), where=total > 0)
    print(f"\n===== {slide.slide_id} {slide.platform} =====")
    print("shape genes x spots:", X.shape, "nonzero:", X.nnz)
    for name, values in (("total UMI", total), ("detected genes", detected), ("mito fraction", mt_frac)):
        qs = np.percentile(values, [0, 1, 5, 25, 50, 75, 95, 99, 100])
        print(name, "quantiles 0/1/5/25/50/75/95/99/100:", np.round(qs, 4))
    for cutoff in (100, 200, 500, 1000):
        print(f"spots total >= {cutoff}: {int((total >= cutoff).sum())}/{X.shape[1]}")
    print("spots mt <= 0.20:", int((mt_frac <= 0.20).sum()) / X.shape[1])


def main() -> None:
    for project, count_path in SAMPLES:
        if project == "GSE200278":
            spatial = Path(str(count_path).replace("raw_counts.csv.gz", "spatial_info.csv.gz"))
            slide = read_slideseqv2(count_path, spatial, slide_id="probe_MBM05_rep3", project_id=project)
        else:
            slide = load_visium(count_path, slide_id="probe", project_id=project)
        summarize(slide)


if __name__ == "__main__":
    main()
