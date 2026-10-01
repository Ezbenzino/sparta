#!/usr/bin/env python
"""
run_25_biological_validation.py —— 屏障场与真实细胞分布的对齐验证
================================================================
审稿人会问：你的 B_cell 和 B_mAb 是数学构造，它们和切片上真实的
免疫细胞分布、肿瘤增殖对齐吗？

本脚本用已有 adata.obs 里的 per-spot signature 做 in silico ground truth：
  B_cell_field  vs  T_NK_n   (预期负相关：屏障高 = T 细胞少)
  B_cell_field  vs  CD8T_n   (同上，CD8 效应 T)
  B_mab         vs  Proliferation_n  (预期正相关：抗体屏障高 = 增殖不受控)

partial Spearman 控制距血管距离 (d_vessel)。
输出：results/validation/biological_validation.json
"""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import csv
import numpy as np
import scanpy as sc

from sparta.io_ import (Paths, load_config, load_graph, patient_map,
                        save_json, set_seed, stamp_run)
from sparta.barrier import (compute_b_cell_field, compute_b_mab, compute_b_meta,
                            scores_from_adata)
from sparta.validate import decoupling_stats


def main():
    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)
    q = cfg["validate"]["decouple_q"]

    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    print(f"Biological alignment validation ({len(slides)} slides)\n")
    print(f"{'slide':<10}{'Bc vs T_NK':>12}{'Bc vs CD8':>12}{'Bm vs Prolif':>14}")
    print("-" * 50)

    rows = {}
    for sid in slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        S = scores_from_adata(adata)

        # vessel distance as control
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]

        # barrier fields
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                  **cfg["barrier"]["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           **cfg["barrier"]["b_mab"])["b_mab"]

        # biological scores from adata.obs (some first-gen ST slides may lack CD8T_n)
        t_nk = np.asarray(adata.obs["T_NK_n"].values, float)
        if "CD8T_n" in adata.obs.columns:
            cd8 = np.asarray(adata.obs["CD8T_n"].values, float)
        else:
            cd8 = t_nk  # fallback: use T_NK as proxy
        prolif = np.asarray(adata.obs["Proliferation_n"].values, float)

        # partial Spearman (control = vessel distance)
        r_bc_tnk = decoupling_stats(bc, t_nk, q=q, control=dv)
        r_bc_cd8 = decoupling_stats(bc, cd8, q=q, control=dv)
        r_bm_prol = decoupling_stats(bm, prolif, q=q, control=dv)

        rows[sid] = dict(
            patient=pmap.get(sid, "?"),
            cohort="CSCC" if sid.startswith("CSCC") else "MEL",
            bc_vs_tnk_rho=r_bc_tnk["rho_partial"],
            bc_vs_tnk_p=r_bc_tnk["p_partial"],
            bc_vs_cd8_rho=r_bc_cd8["rho_partial"],
            bc_vs_cd8_p=r_bc_cd8["p_partial"],
            bm_vs_prolif_rho=r_bm_prol["rho_partial"],
            bm_vs_prolif_p=r_bm_prol["p_partial"],
        )
        print(f"{sid:<10}{r_bc_tnk['rho_partial']:>+12.3f}"
              f"{r_bc_cd8['rho_partial']:>+12.3f}{r_bm_prol['rho_partial']:>+14.3f}")

    # summary
    def med(key):
        return float(np.median([r[key] for r in rows.values()]))

    def n_neg(key, pkey):
        return int(sum(1 for r in rows.values() if r[key] < 0 and r[pkey] < 0.05))

    def n_pos(key, pkey):
        return int(sum(1 for r in rows.values() if r[key] > 0 and r[pkey] < 0.05))

    summary = dict(
        n_slides=len(rows),
        median_bc_vs_tnk=med("bc_vs_tnk_rho"),
        median_bc_vs_cd8=med("bc_vs_cd8_rho"),
        median_bm_vs_prolif=med("bm_vs_prolif_rho"),
        n_bc_tnk_neg_sig=n_neg("bc_vs_tnk_rho", "bc_vs_tnk_p"),
        n_bc_cd8_neg_sig=n_neg("bc_vs_cd8_rho", "bc_vs_cd8_p"),
        n_bm_prolif_pos_sig=n_pos("bm_vs_prolif_rho", "bm_vs_prolif_p"),
    )

    print("\n=== Summary ===")
    print(f"B_cell field vs T_NK signature: median rho = {summary['median_bc_vs_tnk']:+.3f}"
          f" (neg sig in {summary['n_bc_tnk_neg_sig']}/{summary['n_slides']})")
    print(f"B_cell field vs CD8T: median rho = {summary['median_bc_vs_cd8']:+.3f}"
          f" (neg sig in {summary['n_bc_cd8_neg_sig']}/{summary['n_slides']})")
    print(f"B_mAb vs Proliferation: median rho = {summary['median_bm_vs_prolif']:+.3f}"
          f" (pos sig in {summary['n_bm_prolif_pos_sig']}/{summary['n_slides']})")
    print("\nInterpretation:")
    print("  - B_cell vs T_NK/CD8 should be NEGATIVE: barrier-high spots have fewer T cells")
    print("  - B_mAb vs Prolif: positive would mean antibody-blocked spots retain proliferation")

    out = dict(per_slide=rows, summary=summary,
               meta=stamp_run(cfg, {"module": "M25-biological-validation"}))
    p = P.validation("biological_validation.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
