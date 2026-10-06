"""Loaders for the 2026 SPARTA extension spatial-transcriptomics datasets.

Supported inputs
-----------------
* 10x Visium / Visium CytAssist-style HDF5 feature matrices
  (``*_filtered_feature_bc_matrix.h5``).
* 10x Visium three-file Matrix Market groups
  (``*_matrix.mtx.gz`` plus barcodes and features).
* Slide-seqV2 bead tables
  (``*_raw_counts.csv.gz`` plus ``*_spatial_info.csv.gz``).

The loaders intentionally return raw counts and coordinates. Normalisation,
signature scoring, graph construction and filtering happen downstream so that
the same raw object can support QC and reproducible analysis.
"""
from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix, coo_matrix

# Visium capture area nominally contains 55-um-diameter spots. The per-image
# spot diameter in pixels is stored in scalefactors, giving pixels -> um.
VISIUM_SPOT_DIAMETER_UM = 55.0
SLIDESEQ_BEAD_DIAMETER_UM = 10.0


@dataclass
class SpatialExpressionData:
    """Raw spatial expression dataset with platform metadata."""

    slide_id: str
    project_id: str
    platform: str
    counts: csr_matrix
    gene_ids: np.ndarray
    gene_names: np.ndarray
    barcodes: np.ndarray
    coords_um: np.ndarray
    coords_pixel: np.ndarray | None
    in_tissue: np.ndarray
    scalefactors: dict[str, Any]
    metadata: dict[str, Any]

    @property
    def n_genes(self) -> int:
        return int(self.counts.shape[0])

    @property
    def n_spots(self) -> int:
        return int(self.counts.shape[1])

    def to_anndata(self):
        """Return an AnnData object with spots as observations and raw counts."""
        from anndata import AnnData

        obs = {
            "barcode": self.barcodes.astype(str),
            "in_tissue": self.in_tissue.astype(bool),
            "x_um": self.coords_um[:, 0],
            "y_um": self.coords_um[:, 1],
        }
        var = {
            "gene_id": self.gene_ids.astype(str),
            "gene_name": self.gene_names.astype(str),
        }
        adata = AnnData(
            X=self.counts.transpose().tocsr().astype(np.float32),
            obs=obs,
            var=var,
        )
        adata.obsm["spatial_um"] = self.coords_um.astype(np.float64)
        if self.coords_pixel is not None:
            adata.obsm["spatial_pixel"] = self.coords_pixel.astype(np.float64)
        adata.uns["spatial_scalefactors"] = self.scalefactors
        adata.uns["extension_metadata"] = self.metadata
        adata.layers["counts"] = adata.X.copy()
        var_index = self._unique_index(self.gene_ids, fallback=self.gene_names)
        adata.var_names = var_index
        adata.obs_names = self._unique_index(self.barcodes, fallback=self.barcodes)
        return adata

    @staticmethod
    def _unique_index(values: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
        """Return non-empty unique strings for AnnData index names."""
        source = values.astype(str)
        if fallback is not None and (not np.any(source) or len(set(source)) != len(source)):
            source = fallback.astype(str)
        out: list[str] = []
        seen: dict[str, int] = {}
        for value in source:
            base = value if value else "unknown"
            if base not in seen:
                seen[base] = 0
                out.append(base)
            else:
                seen[base] += 1
                out.append(f"{base}.{seen[base]}")
        return np.asarray(out)


def _decode(values: np.ndarray) -> np.ndarray:
    """Decode byte-string arrays from 10x HDF5 files."""
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def read_10x_h5(path: str | Path) -> tuple[csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Read a 10x HDF5 feature matrix as genes-by-barcodes CSR."""
    path = Path(path)
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        barcodes = _decode(matrix["barcodes"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        n_genes, n_barcodes = (int(x) for x in matrix["shape"][:])
        # HDF5 stores one CSR-like compressed axis per barcode (barcodes x genes);
        # the public shape attribute is reported as genes x barcodes.
        raw = csr_matrix(
            (data, indices, indptr), shape=(n_barcodes, n_genes)
        )
        counts = raw.transpose().tocsr()
        features = matrix["features"]
        gene_ids = _decode(features["id"][:])
        gene_names = _decode(features["name"][:])
    return counts, gene_ids, gene_names, barcodes


def _read_gz_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle]


def read_10x_mtx(
    path: str | Path,
) -> tuple[csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Read a three-file 10x Matrix Market group as genes-by-barcodes CSR."""
    mtx_path = Path(path)
    if not mtx_path.name.endswith("_matrix.mtx.gz"):
        raise ValueError(f"Expected a *_matrix.mtx.gz file, got {mtx_path}")
    prefix = mtx_path.name[: -len("_matrix.mtx.gz")]
    barcodes_path = mtx_path.with_name(prefix + "_barcodes.tsv.gz")
    features_path = mtx_path.with_name(prefix + "_features.tsv.gz")
    counts = csr_matrix(mmread(str(mtx_path))).tocsr()
    barcodes = np.asarray([line.strip() for line in _read_gz_lines(barcodes_path) if line.strip()])
    gene_ids: list[str] = []
    gene_names: list[str] = []
    for line in _read_gz_lines(features_path):
        fields = line.split("\t")
        gene_id = fields[0]
        gene_name = fields[1] if len(fields) > 1 else gene_id
        gene_ids.append(gene_id)
        gene_names.append(gene_name)
    return counts, np.asarray(gene_ids), np.asarray(gene_names), barcodes


def _spatial_paths(count_path: Path) -> tuple[Path, Path]:
    """Return scalefactor and tissue-position paths for a 10x count file."""
    if count_path.name.endswith("_filtered_feature_bc_matrix.h5"):
        stem = count_path.name[: -len("_filtered_feature_bc_matrix.h5")]
    elif count_path.name.endswith("_matrix.mtx.gz"):
        stem = count_path.name[: -len("_matrix.mtx.gz")]
    else:
        raise ValueError(f"Unsupported 10x count file: {count_path.name}")
    scale = count_path.with_name(stem + "_scalefactors_json.json.gz")
    for suffix in ("_tissue_positions.csv.gz", "_tissue_positions_list.csv.gz"):
        positions = count_path.with_name(stem + suffix)
        if positions.exists():
            return scale, positions
    raise FileNotFoundError(f"No tissue positions found for {count_path.name}")


def _read_positions(path: Path) -> dict[str, tuple[int, float, float]]:
    """Read old or new Visium tissue-position files."""
    positions: dict[str, tuple[int, float, float]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        first = next(reader)
        rows = reader if first and first[0] == "barcode" else [first, *reader]
        for row in rows:
            if len(row) < 6:
                continue
            barcode = row[0]
            in_tissue = int(row[1])
            y_pixel = float(row[4])  # pxl_row_in_fullres
            x_pixel = float(row[5])  # pxl_col_in_fullres
            positions[barcode] = (in_tissue, x_pixel, y_pixel)
    return positions


def _read_scalefactors(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return json.load(handle)


def _subset_barcodes(
    counts: csr_matrix,
    gene_ids: np.ndarray,
    gene_names: np.ndarray,
    barcodes: np.ndarray,
    ordered_barcodes: np.ndarray,
) -> tuple[csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Subset a genes-by-barcodes matrix to ordered, available barcodes."""
    index = {barcode: i for i, barcode in enumerate(barcodes)}
    columns = [index[barcode] for barcode in ordered_barcodes if barcode in index]
    if not columns:
        raise ValueError("No barcodes overlap with spatial positions")
    return counts[:, columns], gene_ids, gene_names, np.asarray(ordered_barcodes)


def load_visium(
    count_path: str | Path,
    slide_id: str,
    project_id: str,
) -> SpatialExpressionData:
    """Load a Visium HDF5 or Matrix Market sample."""
    count_path = Path(count_path)
    if count_path.name.endswith(".h5"):
        counts, gene_ids, gene_names, barcodes = read_10x_h5(count_path)
    else:
        counts, gene_ids, gene_names, barcodes = read_10x_mtx(count_path)
    scale_path, positions_path = _spatial_paths(count_path)
    scalefactors = _read_scalefactors(scale_path)
    positions = _read_positions(positions_path)
    ordered = np.asarray([barcode for barcode in barcodes if barcode in positions])
    counts, gene_ids, gene_names, barcodes = _subset_barcodes(
        counts, gene_ids, gene_names, barcodes, ordered
    )
    in_tissue = np.asarray([positions[barcode][0] for barcode in barcodes], dtype=int)
    x_pixel = np.asarray([positions[barcode][1] for barcode in barcodes], dtype=float)
    y_pixel = np.asarray([positions[barcode][2] for barcode in barcodes], dtype=float)
    pixel_per_um = float(scalefactors["spot_diameter_fullres"]) / VISIUM_SPOT_DIAMETER_UM
    coords_pixel = np.column_stack([x_pixel, y_pixel])
    coords_um = coords_pixel / pixel_per_um
    metadata = {
        "count_file": str(count_path),
        "scalefactors_file": str(scale_path),
        "positions_file": str(positions_path),
        "pixel_per_um": pixel_per_um,
        "coordinate_unit": "um",
    }
    return SpatialExpressionData(
        slide_id=slide_id,
        project_id=project_id,
        platform="Visium",
        counts=counts,
        gene_ids=gene_ids,
        gene_names=gene_names,
        barcodes=barcodes,
        coords_um=coords_um,
        coords_pixel=coords_pixel,
        in_tissue=in_tissue,
        scalefactors=scalefactors,
        metadata=metadata,
    )


def _read_slideseq_spatial(
    path: Path,
) -> tuple[dict[str, tuple[float, float]], dict[str, str]]:
    """Read Slide-seqV2 spatial coordinates and published RCTD labels."""
    coords: dict[str, tuple[float, float]] = {}
    labels: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            barcode = row.get("barcode") or next(iter(row.values()))
            coords[barcode] = (float(row["xcoord"]), float(row["ycoord"]))
            labels[barcode] = row.get("rctd_cell_type", "")
    return coords, labels


def read_slideseqv2(
    count_path: str | Path,
    spatial_path: str | Path,
    slide_id: str,
    project_id: str,
) -> SpatialExpressionData:
    """Load a Slide-seqV2 raw count CSV and matching spatial-info CSV."""
    count_path = Path(count_path)
    spatial_path = Path(spatial_path)
    spatial, labels = _read_slideseq_spatial(spatial_path)

    rows: list[int] = []
    cols: list[int] = []
    data: list[int] = []
    gene_names: list[str] = []
    with gzip.open(count_path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        all_barcodes = header[1:]
        # Multiple barcodes can share one localized coordinate in the published
        # table. Assign them to one spatial point; CSR construction sums counts.
        coord_to_column: dict[tuple[float, float], int] = {}
        ordered_barcodes: list[str] = []
        group_barcodes: list[list[str]] = []
        barcode_to_column: dict[str, int] = {}
        for barcode in all_barcodes:
            if barcode not in spatial:
                continue
            xy = spatial[barcode]
            if xy not in coord_to_column:
                coord_to_column[xy] = len(ordered_barcodes)
                ordered_barcodes.append(barcode)
                group_barcodes.append([])
            column = coord_to_column[xy]
            barcode_to_column[barcode] = column
            group_barcodes[column].append(barcode)
        for gene_i, row in enumerate(reader):
            gene_names.append(row[0])
            for barcode, value in zip(all_barcodes, row[1:]):
                if barcode not in barcode_to_column or value in {"", "0", "0.0"}:
                    continue
                number = float(value)
                if number != 0:
                    rows.append(gene_i)
                    cols.append(barcode_to_column[barcode])
                    data.append(int(number))
    barcodes = np.asarray(ordered_barcodes)
    gene_names_array = np.asarray(gene_names)
    counts = coo_matrix(
        (np.asarray(data), (np.asarray(rows), np.asarray(cols))),
        shape=(len(gene_names_array), len(barcodes)),
    ).tocsr()
    coords_um = np.asarray(
        [spatial[barcode] for barcode in barcodes], dtype=float
    )
    metadata = {
        "count_file": str(count_path),
        "spatial_file": str(spatial_path),
        "bead_diameter_um": SLIDESEQ_BEAD_DIAMETER_UM,
        "coordinate_unit": "um",
        "published_rctd_labels_available": True,
        "published_rctd_labels_not_used": True,
        "duplicate_coordinates_aggregated": True,
        "group_barcodes": group_barcodes,
        "rctd_cell_type": [labels[barcode] for barcode in barcodes],
    }
    return SpatialExpressionData(
        slide_id=slide_id,
        project_id=project_id,
        platform="Slide-seqV2",
        counts=counts,
        gene_ids=gene_names_array,
        gene_names=gene_names_array,
        barcodes=barcodes,
        coords_um=coords_um,
        coords_pixel=None,
        in_tissue=np.ones(len(barcodes), dtype=int),
        scalefactors={},
        metadata=metadata,
    )
