#!/usr/bin/env python
"""
run_26_radius_sensitivity.py —— 图构建 radius 稳健性
====================================================
审稿人会问："你改邻接半径，结论还稳吗？"

对每张切片用 3-4 个 radius 重建图、重新定义源汇、重算两个屏障。
检查两个主要结论的稳健性：
  (1) B_cell vs B_mab partial Spearman ρ（median 应在 +0.15~+0.25，大部分正）
  (2) 30% stromal reduction 后两个屏障是否同降

输出：results/validation/radius_sensitivity.json
"""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import csv
import numpy as np
import scanpy as sc

from sparta.io_ import Paths, load_config, patient_map, save_json, set_seed, stamp_run
from sparta.graph import build_graph_radius, define_source_sink
from sparta.barrier import (compute_b_cell, compute_b_mab, compute_b_meta,
                            scores_from_adata)
from sparta.validate import decoupling_stats


def _col(adata, name):
    return np.asarray(adata.obs[name + "_n"].values, float) if name + "_n" in adata.obs \
        else np.full(adata.n_obs, 0.5)


def main():
    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)
    q = cfg["validate"]["decouple_q"]

    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        slides = [(r["slide_id"], r["platform"]) for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    # radius grids
    radius_grid = {
        "visium": [100.0, 150.0, 200.0, 250.0],
        "legacy_st": [200.0, 300.0, 400.0],
    }

    print(f"Radius sensitivity ({len(slides)} slides)\n")

    all_rows = {}
    for sid, platform in slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
        except FileNotFoundError:
            continue
        coords = np.asarray(adata.obsm["spatial_um"], float)
        S = scores_from_adata(adata)
        radii = radius_grid.get(platform, [150.0])

        rows_sid = {}
        for r_um in radii:
            A, D = build_graph_radius(coords, radius_um=r_um)
            ss = define_source_sink(
                A, _col(adata, "Endothelial"), _col(adata, "T_NK"), _col(adata, "Malignant"),
                **{k: cfg["source_sink"][k] for k in ("q_vessel", "q_immune_nbr", "q_malig", "q_core")},
                fallback_border_coords=coords if cfg["source_sink"]["use_border_fallback"] else None,
            )
            source, sink, vessel = ss["source"], ss["sink"], ss["vessel"]

            # baseline barriers
            bc = compute_b_cell(A, S["ecm"], S["caf"], source, sink,
                                **cfg["barrier"]["b_cell"])["b_cell"]
            bm_res = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                                   **cfg["barrier"]["b_mab"])
            bm = bm_res["b_mab"]

            # vessel distance control
            dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                                D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]

            # coupling (use per-node fields for partial rho)
            from sparta.barrier import compute_b_cell_field
            bcf = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                       **cfg["barrier"]["b_cell"])["b_cell_field"]
            d_coup = decoupling_stats(bcf, bm, q=q, control=dv)

            # 30% stromal reduction
            bc_red = compute_b_cell(A, S["ecm"]*0.7, S["caf"]*0.7, source, sink,
                                    **cfg["barrier"]["b_cell"])["b_cell"]
            bm_red = compute_b_mab(A, S["ecm"]*0.7, S["crosslink"]*0.7, S["ag_target"], vessel,
                                   **cfg["barrier"]["b_mab"])
            _mask_red = np.asarray(bm_red["reachable"], dtype=bool)
            _mask_base = np.asarray(bm_res["reachable"], dtype=bool)
            bm_red_mean = float(np.nanmean(bm_red["b_mab"][_mask_red])) if _mask_red.any() else float("nan")
            bm_mean = float(np.nanmean(bm_res["b_mab"][_mask_base])) if _mask_base.any() else float("nan")
            bc_drop = (bc - bc_red) / bc * 100 if bc > 0 else 0
            bm_drop = (bm_mean - bm_red_mean) / bm_mean * 100 if bm_mean and bm_mean > 0 else 0

            rows_sid[r_um] = dict(
                n_nodes=int(A.shape[0]), n_edges=int(A.nnz // 2),
                rho_partial=d_coup["rho_partial"],
                bc_drop_30=bc_drop, bm_drop_30=bm_drop,
            )

        all_rows[sid] = dict(platform=platform, patient=pmap.get(sid, "?"), radii=rows_sid)
        # print compact
        med_rho = np.median([v["rho_partial"] for v in rows_sid.values()])
        med_bc = np.median([v["bc_drop_30"] for v in rows_sid.values()])
        med_bm = np.median([v["bm_drop_30"] for v in rows_sid.values()])
        print(f"{sid:<10} rho med={med_rho:+.3f}  bc_drop med={med_bc:5.1f}%  bm_drop med={med_bm:5.1f}%")

    # summary: across all slides and radii
    all_rhos = [v["rho_partial"] for sid in all_rows.values() for v in sid["radii"].values()]
    all_bc = [v["bc_drop_30"] for sid in all_rows.values() for v in sid["radii"].values()]
    all_bm = [v["bm_drop_30"] for sid in all_rows.values() for v in sid["radii"].values()]
    both_drop = sum(1 for sid in all_rows.values() for v in sid["radii"].values()
                    if v["bc_drop_30"] > 0 and v["bm_drop_30"] > 0)
    total_cells = len(all_rhos)

    summary = dict(
        n_slides=len(all_rows),
        median_rho_all=float(np.median(all_rhos)),
        frac_rho_positive=float(np.mean(np.array(all_rhos) > 0)),
        median_bc_drop=float(np.median(all_bc)),
        median_bm_drop=float(np.median(all_bm)),
        frac_both_drop=both_drop / total_cells,
        total_radius_configs=total_cells,
    )

    print("\n=== Summary ===")
    print(f"Across {len(all_rows)} slides x 3-4 radii = {total_cells} configurations:")
    print(f"  Median partial rho (B_cell field vs B_mab): {summary['median_rho_all']:+.3f}")
    print(f"  Fraction positive: {summary['frac_rho_positive']*100:.1f}%")
    print(f"  Median B_cell drop @30% stromal reduction: {summary['median_bc_drop']:.1f}%")
    print(f"  Median B_mAb drop @30% stromal reduction: {summary['median_bm_drop']:.1f}%")
    print(f"  Fraction where BOTH barriers drop: {summary['frac_both_drop']*100:.1f}%")

    out = dict(per_slide={k: {rr: vv for rr, vv in v["radii"].items()}
                          for k, v in all_rows.items()},
               summary=summary, meta=stamp_run(cfg, {"module": "M26-radius-sensitivity"}))
    p = P.validation("radius_sensitivity.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
