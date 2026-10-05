#!/usr/bin/env python
"""Neighborhood-weight sensitivity for the BANKSY-style descriptive benchmark.

This is not the official BANKSY implementation and is not a performance test.
It varies the neighborhood contribution while holding the eight-domain count,
PCA recipe, graph, seed, and cut-edge comparison fixed.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np
import scanpy as sc

from sparta.io_ import Paths, load_config, load_graph, load_json, save_json, stamp_run
from sparta.barrier import edge_pairs
from run_13_benchmark_tools import banksy_style_domains

LAM_GRID = (0.1, 0.3, 0.5, 0.7)
N_DOMAINS = 8
N_PCS = 20


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    s2 = load_json(P.validation("s2_matched_selection.json"))["per_slide"]
    rows = {}
    for sid in sorted(s2):
        adata = sc.read_h5ad(P.scored(sid))
        A, _, _, _, _, _ = load_graph(P.graph(sid))
        cut = load_json(P.mincut(sid))["cut_edges"]
        pairs = edge_pairs(A, upper_only=True)
        edge_id = {(int(min(u, v)), int(max(u, v))): i
                   for i, (u, v) in enumerate(pairs)}
        cut_idx = np.array([edge_id[(min(u, v), max(u, v))]
                            for u, v in cut
                            if (min(u, v), max(u, v)) in edge_id], int)
        rec = {"n_nodes": int(A.shape[0]), "runs": {}}
        for lam in LAM_GRID:
            labels = banksy_style_domains(adata, A, N_DOMAINS, lam,
                                          N_PCS, cfg["seed"])
            boundary = labels[pairs[:, 0]] != labels[pairs[:, 1]]
            if len(cut_idx):
                frac_boundary = float(boundary.mean())
                frac_cut_boundary = float(boundary[cut_idx].mean())
                enrichment = frac_cut_boundary / max(frac_boundary, 1e-12)
                precision = float(boundary[cut_idx].sum() /
                                  max(boundary.sum(), 1))
            else:
                frac_boundary = frac_cut_boundary = enrichment = precision = float("nan")
            rec["runs"][str(lam)] = {
                "n_domains": int(len(np.unique(labels))),
                "frac_boundary_edges": frac_boundary,
                "frac_boundary_within_cut": frac_cut_boundary,
                "enrichment": float(enrichment),
                "precision_of_boundary_for_cut": float(precision),
            }
        rows[sid] = rec
        vals = [rec["runs"][str(l)]["enrichment"] for l in LAM_GRID]
        print(f"{sid:<8} enrichment {min(vals):.2f}–{max(vals):.2f}")

    summary = {}
    for lam in LAM_GRID:
        key = str(lam)
        for cohort, ids in (
            ("all_primary", list(rows)),
            ("visium", [s for s in rows if s in {
                "MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"}]),
            ("legacy_st", [s for s in rows if s.startswith("CSCC") and
                            s not in {"CSCC01", "CSCC02", "CSCC03", "CSCC04"}]),
        ):
            enrich = [rows[s]["runs"][key]["enrichment"] for s in ids
                      if np.isfinite(rows[s]["runs"][key]["enrichment"])]
            precision = [rows[s]["runs"][key]["precision_of_boundary_for_cut"]
                         for s in ids
                         if np.isfinite(rows[s]["runs"][key]["precision_of_boundary_for_cut"])]
            summary[f"lambda_{key}_{cohort}"] = {
                "n": len(enrich),
                "median_enrichment": float(np.median(enrich)),
                "enrichment_range": [float(np.min(enrich)), float(np.max(enrich))],
                "median_precision": float(np.median(precision)),
            }
    out = {
        "params": {
            "implementation": "BANKSY-style augmented PCA/neighborhood-mean k-means; not official BANKSY",
            "n_domains": N_DOMAINS,
            "lambda_grid": list(LAM_GRID),
            "n_pcs": N_PCS,
            "seed": int(cfg["seed"]),
            "fixed_graph": True,
            "comparison": "minimum-cut edge overlap with domain-boundary edges",
        },
        "per_slide": rows,
        "summary": summary,
        "meta": stamp_run(cfg, {"module": "M34-benchmark-lambda-sensitivity"}),
    }
    save_json(P.validation("benchmark_lambda_sensitivity.json"), out)
    print(f"\nWrote {P.validation('benchmark_lambda_sensitivity.json')}")


if __name__ == "__main__":
    main()
