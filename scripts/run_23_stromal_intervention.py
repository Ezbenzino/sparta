#!/usr/bin/env python
"""
run_23_stromal_intervention.py —— In silico stromal co-targeting 预测
====================================================================
这是论文的 punchline：如果 stromal modification（TGF-β/LOX 抑制、抗纤维化）
把基质密度和交联度下调 X%，两个屏障会同时降多少？

做法：
  对每张切片，把 ecm 和 crosslink 各乘 (1 - reduction)，重算 B_cell 和 B_mAb，
  报告两个 barrier 各自的下降幅度。
  如果两个都显著下降 -> 支持"matrix-directed intervention 同时改善两种递送"的预测。

输出：results/validation/stromal_intervention.json
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.barrier import (compute_b_cell, compute_b_mab,  # noqa: E402
                            scores_from_adata)

REDUCTIONS = [0.0, 0.20, 0.30, 0.50]  # 干预强度


def main():
    import scanpy as sc

    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)

    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = cfg["barrier"]["b_mab"]

    # 从 ledger 拿 ingested 的
    import csv
    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    print(f"In silico stromal intervention  (reduction levels: {REDUCTIONS})")
    print(f"Slides ({len(slides)}): {slides}\n")

    rows = {}
    for sid in slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        S = scores_from_adata(adata)

        results = {}
        for red in REDUCTIONS:
            ecm = S["ecm"] * (1.0 - red)
            caf = S["caf"] * (1.0 - red)
            xl = S["crosslink"] * (1.0 - red)

            bc = compute_b_cell(A, ecm, caf, source, sink, **cfg_cell)
            bm = compute_b_mab(A, ecm, xl, S["ag_target"], vessel, **cfg_mab)
            bm_core = bm["b_mab"][bm["reachable"]]
            results[f"red_{int(red*100)}"] = dict(
                b_cell=float(bc["b_cell"]),
                b_mab_mean=float(np.nanmean(bm_core)) if len(bm_core) else float("nan"),
            )

        b0 = results["red_0"]
        rows[sid] = dict(
            patient=pmap.get(sid, "?"),
            cohort="CSCC" if sid.startswith("CSCC") else "MEL",
            b_cell_baseline=b0["b_cell"],
            b_mab_baseline=b0["b_mab_mean"],
            pct_drop_b_cell_at_30=100.0 * (1.0 - results["red_30"]["b_cell"] / b0["b_cell"]) if b0["b_cell"] else float("nan"),
            pct_drop_b_mab_at_30=100.0 * (1.0 - results["red_30"]["b_mab_mean"] / b0["b_mab_mean"]) if b0["b_mab_mean"] else float("nan"),
            pct_drop_b_cell_at_50=100.0 * (1.0 - results["red_50"]["b_cell"] / b0["b_cell"]) if b0["b_cell"] else float("nan"),
            pct_drop_b_mab_at_50=100.0 * (1.0 - results["red_50"]["b_mab_mean"] / b0["b_mab_mean"]) if b0["b_mab_mean"] else float("nan"),
            per_level=results,
        )
        print(f"{sid:<8} ({rows[sid]['cohort']:<3} P{pmap.get(sid,'?'):<4})  "
              f"B_cell drop@30%: {rows[sid]['pct_drop_b_cell_at_30']:5.1f}%   "
              f"B_mAb drop@30%: {rows[sid]['pct_drop_b_mab_at_30']:5.1f}%")

    # 汇总
    def _med(key, cohort=None):
        vals = [r[key] for r in rows.values()
                if (cohort is None or r["cohort"] == cohort)
                and np.isfinite(r[key])]
        return float(np.median(vals)) if vals else float("nan")

    summary = dict(
        n_slides=len(rows),
        n_patients=len({r["patient"] for r in rows.values()}),
        median_pct_drop_b_cell_30=_med("pct_drop_b_cell_at_30"),
        median_pct_drop_b_mab_30=_med("pct_drop_b_mab_at_30"),
        median_pct_drop_b_cell_50=_med("pct_drop_b_cell_at_50"),
        median_pct_drop_b_mab_50=_med("pct_drop_b_mab_at_50"),
        cscc=dict(
            median_drop_b_cell_30=_med("pct_drop_b_cell_at_30", "CSCC"),
            median_drop_b_mab_30=_med("pct_drop_b_mab_at_30", "CSCC"),
        ),
        mel=dict(
            median_drop_b_cell_30=_med("pct_drop_b_cell_at_30", "MEL"),
            median_drop_b_mab_30=_med("pct_drop_b_mab_at_30", "MEL"),
        ),
        both_drop_at_30_pct=float(np.mean([
            (r["pct_drop_b_cell_at_30"] > 0 and r["pct_drop_b_mab_at_30"] > 0)
            for r in rows.values()
            if np.isfinite(r["pct_drop_b_cell_at_30"]) and np.isfinite(r["pct_drop_b_mab_at_30"])
        ]) * 100),
    )

    print("\n=== Summary (median across slides) ===")
    print(f"At 30% stromal reduction:")
    print(f"  B_cell drops by {summary['median_pct_drop_b_cell_30']:.1f}%  (cSCC {summary['cscc']['median_drop_b_cell_30']:.1f}%, MEL {summary['mel']['median_drop_b_cell_30']:.1f}%)")
    print(f"  B_mAb drops by {summary['median_pct_drop_b_mab_30']:.1f}%  (cSCC {summary['cscc']['median_drop_b_mab_30']:.1f}%, MEL {summary['mel']['median_drop_b_mab_30']:.1f}%)")
    print(f"  Slides where BOTH barriers drop: {summary['both_drop_at_30_pct']:.0f}%")
    print(f"At 50% stromal reduction:")
    print(f"  B_cell drops by {summary['median_pct_drop_b_cell_50']:.1f}%,  B_mAb drops by {summary['median_pct_drop_b_mab_50']:.1f}%")

    out = dict(per_slide=rows, summary=summary,
               reductions=REDUCTIONS,
               meta=stamp_run(cfg, {"module": "M23-stromal-intervention"}))
    p = P.validation("stromal_intervention.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
