"""Inspect on-disk structure of the 2026 extension datasets.

The script reads HDF5 metadata, gzip headers and the first few CSV rows. It does
not load large count tables into memory.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "data" / "external" / "extension_2026"
GSES = ["GSE200278", "GSE289745", "GSE300445", "GSE316760", "GSE320041", "GSE321832"]


def first_gzip_lines(path: Path, n: int = 4) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [next(handle).rstrip("\n") for _ in range(n)]


def inspect_h5(path: Path) -> None:
    print(f"\nH5: {path.relative_to(EXT)}")
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        print("  attrs:", dict(matrix.attrs))
        for key, item in matrix.items():
            if isinstance(item, h5py.Dataset):
                print(f"  {key}: shape={item.shape}, first={item[:3]}")
        features = matrix["features"]
        for key, item in features.items():
            if isinstance(item, h5py.Dataset):
                print(f"  features/{key}: first={item[:5]}")
    scale = Path(str(path).replace("_filtered_feature_bc_matrix.h5", "_scalefactors_json.json.gz"))
    position = Path(str(path).replace("_filtered_feature_bc_matrix.h5", "_tissue_positions.csv.gz"))
    if scale.exists():
        with gzip.open(scale, "rt", encoding="utf-8") as handle:
            print("  scalefactors:", json.load(handle))
    if position.exists():
        print("  positions first lines:", first_gzip_lines(position, 3))


def inspect_mtx(path: Path) -> None:
    print(f"\nMTX group: {path.name}")
    print("  mtx header:", first_gzip_lines(path, 4))
    stem = path.name.replace("_matrix.mtx.gz", "")
    for suffix in ("_barcodes.tsv.gz", "_features.tsv.gz", "_scalefactors_json.json.gz"):
        candidate = path.with_name(stem + suffix)
        if candidate.exists():
            if suffix.endswith("json.gz"):
                with gzip.open(candidate, "rt", encoding="utf-8") as handle:
                    print("  scalefactors:", json.load(handle))
            else:
                print(f"  {suffix}:", first_gzip_lines(candidate, 3))
    position = path.with_name(stem + "_tissue_positions_list.csv.gz")
    if not position.exists():
        position = path.with_name(stem + "_tissue_positions.csv.gz")
    if position.exists():
        print("  positions:", first_gzip_lines(position, 3))


def inspect_standard() -> None:
    for gse in GSES[1:]:
        root = EXT / gse
        h5s = sorted(root.rglob("*filtered_feature_bc_matrix.h5"))
        mtxs = sorted(root.rglob("*matrix.mtx.gz"))
        print(f"\n===== {gse}: {len(h5s)} h5, {len(mtxs)} mtx groups =====")
        for path in h5s[:2]:
            inspect_h5(path)
        for path in mtxs[:2]:
            inspect_mtx(path)


def inspect_slideseq() -> None:
    root = EXT / "GSE200278"
    raws = sorted(root.rglob("*raw_counts.csv.gz"))
    print(f"\n===== GSE200278: {len(raws)} count files =====")
    for raw in raws[:4]:
        spatial = Path(str(raw).replace("raw_counts.csv.gz", "spatial_info.csv.gz"))
        print(f"\nSample: {raw.name}")
        raw_head = pd.read_csv(raw, nrows=3, compression="gzip")
        print("  raw columns:", list(raw_head.columns[:6]))
        print(raw_head.iloc[:2, :5].to_string(index=False))
        spatial_head = pd.read_csv(spatial, nrows=5, compression="gzip")
        print("  spatial columns:", list(spatial_head.columns))
        print(spatial_head.to_string(index=False))


def main() -> None:
    inspect_standard()
    inspect_slideseq()


if __name__ == "__main__":
    main()
