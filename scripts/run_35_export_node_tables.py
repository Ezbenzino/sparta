#!/usr/bin/env python
"""
run_35_export_node_tables.py -- export per-spot model inputs to lightweight .npz node tables
============================================================================================
Input : data/interim/{sid}.scored.h5ad  (M2 output; X = log1p(CP10k), layers["counts"] = UMI)
Output: data/interim/{sid}.nodes.npz     (one per ingested section)
        results/validation/node_tables_manifest.json

Why this exists (2026-10-04, IS submission)
-------------------------------------------
All downstream statistics (spatial null, patient-level model, shared-input ablation,
graph-disconnection table, spatial maps) only need: the rank-normalised signature
columns (*_n), the micrometre coordinates and a small marker-gene panel.  Exporting
them once makes every later step runnable with numpy/scipy/networkx only, i.e. on a
machine without scanpy/anndata/h5py, and lets reviewers re-run the statistics from
the archived node tables without the 2-3 GB of h5ad files.

Nothing here changes a model input: columns are copied verbatim from adata.obs.
Row order equals adata.obs order, which is the node order of {sid}.graph.npz.

Usage (Windows):
    .\.venv\Scripts\python.exe scripts\run_35_export_node_tables.py
    .\.venv\Scripts\python.exe scripts\run_35_export_node_tables.py --slides MEL01 CSCC05
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

GENES = [
    # immune
    "PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD8B", "CD4", "GZMB", "GZMK", "PRF1",
    "NKG7", "IFNG", "CXCL9", "CXCL10", "CXCL13", "CCL19", "CCL21", "CCL5",
    "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
    "CD68", "CD163", "LYZ", "MS4A1", "CD79A", "JCHAIN",
    # vessel
    "PECAM1", "VWF", "CDH5", "PLVAP", "ACKR1",
    # stroma / matrix
    "COL1A1", "COL1A2", "COL3A1", "FN1", "FAP", "ACTA2", "PDGFRB", "POSTN",
    "LOX", "LOXL1", "LOXL2", "PLOD2", "TGM2", "TGFB1", "TGFBI", "MMP2", "MMP9",
    # tumour lineage
    "KRT5", "KRT14", "KRT6A", "KRT17", "PMEL", "MLANA", "TYR", "SOX10",
    "EPCAM", "KRT8", "KRT18", "MUC1",
    # state
    "MKI67", "TOP2A", "HIF1A", "VEGFA", "CA9", "ABCB1", "ABCG2", "FCGRT",
]


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--slides", nargs="*", default=None,
                    help="default: every ledger row with status == ingested")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    import numpy as np
    import scipy.sparse as sp
    import scanpy as sc
    from sparta.io_ import Paths, load_config, save_json, stamp_run

    try:
        sys.stdout.reconfigure(errors="replace")
    except Exception:  # noqa: BLE001
        pass

    cfg = load_config(args.config)
    P = Paths(cfg)
    if args.slides:
        slides = args.slides
    else:
        with open(P.root / "data" / "ledger.csv", encoding="utf-8") as f:
            slides = [r["slide_id"] for r in csv.DictReader(f) if r.get("status") == "ingested"]

    manifest = {}
    for sid in slides:
        src = P.scored(sid)
        if not src.exists():
            print(f"{sid}: SKIP (missing {src.name})")
            continue
        ad = sc.read_h5ad(src)
        out = {}
        num_cols = []
        for c in ad.obs.columns:
            try:
                v = np.asarray(ad.obs[c].values, dtype=float)
            except (TypeError, ValueError):
                continue
            out[f"obs__{c}"] = v
            num_cols.append(c)
        out["obs_names"] = np.asarray(ad.obs_names.astype(str))
        for key in ("spatial_um", "spatial"):
            if key in ad.obsm:
                out[f"obsm__{key}"] = np.asarray(ad.obsm[key], dtype=float)
        # marker genes: log1p(CP10k) from X, raw UMI from layers["counts"]
        present = [g for g in GENES if g in ad.var_names]
        if present:
            idx = []
            for g in present:
                loc = ad.var_names.get_loc(g)
                if not isinstance(loc, (int, np.integer)):   # duplicated symbol -> first match
                    loc = int(np.flatnonzero(np.asarray(ad.var_names) == g)[0])
                idx.append(int(loc))
            X = ad.X[:, idx]
            X = X.toarray() if sp.issparse(X) else np.asarray(X)
            out["genes"] = np.asarray(present)
            out["gene_logcp10k"] = X.astype(np.float32)
            if "counts" in ad.layers:
                C = ad.layers["counts"][:, idx]
                C = C.toarray() if sp.issparse(C) else np.asarray(C)
                out["gene_counts"] = C.astype(np.float32)
        if "counts" in ad.layers:
            tot = ad.layers["counts"].sum(axis=1)
            out["total_umi"] = np.asarray(tot).ravel().astype(float)
        run_meta = ad.uns.get("sparta_run", {})
        try:
            run_meta = {k: (v if isinstance(v, (str, int, float, bool)) or v is None else str(v))
                        for k, v in dict(run_meta).items()}
        except Exception:  # noqa: BLE001
            run_meta = {"repr": str(run_meta)}
        meta = dict(slide=sid, n_obs=int(ad.n_obs), n_vars=int(ad.n_vars),
                    numeric_obs_columns=num_cols, genes_present=present,
                    genes_missing=[g for g in GENES if g not in present],
                    x_is_log1p_cp10k=bool("log1p" in ad.uns),
                    source_file=src.name, source_sha256=_sha256(src),
                    m2_run=run_meta)
        out["meta"] = np.array([json.dumps(meta, ensure_ascii=False)], dtype=object)
        dst = P.interim / f"{sid}.nodes.npz"
        np.savez_compressed(dst, **out)
        manifest[sid] = dict(file=dst.name, sha256=_sha256(dst), n_obs=int(ad.n_obs),
                             n_numeric_obs=len(num_cols), n_genes=len(present),
                             tumor_type=run_meta.get("tumor_type"))
        print(f"{sid}: {ad.n_obs} spots, {len(num_cols)} obs columns, {len(present)} genes -> {dst.name}")
        del ad

    save_json(P.validation("node_tables_manifest.json"),
              dict(slides=manifest, genes=GENES,
                   meta=stamp_run(cfg, {"module": "M35-export-node-tables"})))
    print(f"\nDone: {len(manifest)} node tables. Manifest: results/validation/node_tables_manifest.json")


if __name__ == "__main__":
    main()
