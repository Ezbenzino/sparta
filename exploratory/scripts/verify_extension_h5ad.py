"""Verify extension h5ad files independently of the ingestion path."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import anndata

ROOT = Path(__file__).resolve().parents[2]
H5AD_DIR = ROOT / "data" / "interim" / "extension_2026"
QC_JSON = ROOT / "results" / "qc" / "extension_2026_qc_summary.json"


def main() -> None:
    qc = json.loads(QC_JSON.read_text(encoding="utf-8"))
    h5ads = sorted(H5AD_DIR.glob("*.h5ad"))
    if len(h5ads) != 50:
        raise AssertionError(f"Expected 50 h5ad, found {len(h5ads)}")
    for path in h5ads:
        if path.stat().st_size <= 0:
            raise AssertionError(f"Empty h5ad: {path.name}")
        adata = anndata.read_h5ad(path)
        spec = qc[path.stem]
        if adata.n_obs != int(spec["retained_spots"]):
            raise AssertionError(f"{path.name}: obs {adata.n_obs} != QC {spec['retained_spots']}")
        if adata.n_vars != int(spec["retained_genes"]):
            raise AssertionError(f"{path.name}: vars {adata.n_vars} != QC {spec['retained_genes']}")
        if "spatial_um" not in adata.obsm:
            raise AssertionError(f"{path.name}: missing spatial_um")
        if adata.obsm["spatial_um"].shape != (adata.n_obs, 2):
            raise AssertionError(f"{path.name}: bad spatial shape")
    print(f"Verified {len(h5ads)} extension h5ad files")
    print("All retained dimensions and spatial coordinates match QC JSON")


if __name__ == "__main__":
    main()
