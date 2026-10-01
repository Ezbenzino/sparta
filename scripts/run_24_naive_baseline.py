#!/usr/bin/env python
"""
run_24_naive_baseline.py —— 为什么需要图？naive 基质评分 vs graph barrier
=========================================================================
审稿人一定会问："你为什么需要图？直接用 0.5*(ECM+CAF) 评分不就行了？"

本脚本定量回答：
  naive score = 0.5*(ecm + caf)  （节点级，无图结构、无 min-cut、无 screened Poisson）
  graph field = compute_b_cell_field(...)  （Dijkstra shortest-path cost）

比较两者与 B_mAb 的耦合强度（partial Spearman rho，控制距血管距离）。
如果 graph 的耦合不显著强于 naive，说明图结构没加值；
如果 graph 显著更强，说明 spatial arrangement 本身有信息。

输出：results/validation/naive_baseline.json
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.barrier import (compute_b_cell_field, compute_b_mab,  # noqa: E402
                            compute_b_meta, scores_from_adata)
from sparta.validate import decoupling_stats  # noqa: E402


def main():
    import scanpy as sc
    import csv

    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)

    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    q = cfg["validate"]["decouple_q"]

    print(f"Naive baseline vs graph barrier  ({len(slides)} slides)\n")
    print(f"{'slide':<10}{'graph rho_p':>12}{'naive rho_p':>13}{'delta':>9}")
    print("-" * 46)

    rows = {}
    for sid in slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        S = scores_from_adata(adata)

        # control: 距血管加权图距
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]

        # B_mAb 场
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           **cfg["barrier"]["b_mab"])["b_mab"]

        # graph B_cell field
        bf_graph = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                        **cfg["barrier"]["b_cell"])["b_cell_field"]

        # naive score: 0.5*(ecm+caf)，无图结构
        bf_naive = 0.5 * (S["ecm"] + S["caf"])

        d_graph = decoupling_stats(bf_graph, bm, q=q, control=dv)
        d_naive = decoupling_stats(bf_naive, bm, q=q, control=dv)

        rows[sid] = dict(
            patient=pmap.get(sid, "?"),
            cohort="CSCC" if sid.startswith("CSCC") else "MEL",
            graph_rho_partial=d_graph["rho_partial"],
            graph_p=d_graph["p_partial"],
            naive_rho_partial=d_naive["rho_partial"],
            naive_p=d_naive["p_partial"],
            delta_rho=d_graph["rho_partial"] - d_naive["rho_partial"],
        )
        print(f"{sid:<10}{d_graph['rho_partial']:>+12.3f}{d_naive['rho_partial']:>+13.3f}"
              f"{d_graph['rho_partial']-d_naive['rho_partial']:>+9.3f}")

    # 汇总
    graph_rhos = np.array([r["graph_rho_partial"] for r in rows.values()])
    naive_rhos = np.array([r["naive_rho_partial"] for r in rows.values()])
    deltas = graph_rhos - naive_rhos

    summary = dict(
        n_slides=len(rows),
        median_graph_rho=float(np.median(graph_rhos)),
        median_naive_rho=float(np.median(naive_rhos)),
        median_delta=float(np.median(deltas)),
        n_graph_higher=int((deltas > 0).sum()),
        n_graph_significant=int((np.array([r["graph_p"] for r in rows.values()]) < 0.05).sum()),
        n_naive_significant=int((np.array([r["naive_p"] for r in rows.values()]) < 0.05).sum()),
    )

    print("\n=== Summary ===")
    print(f"Median graph B_cell field vs B_mAb partial rho: {summary['median_graph_rho']:+.3f}")
    print(f"Median naive (0.5*(ECM+CAF)) vs B_mAb partial rho: {summary['median_naive_rho']:+.3f}")
    print(f"Median delta (graph - naive): {summary['median_delta']:+.3f}")
    print(f"Slides where graph coupling > naive: {summary['n_graph_higher']}/{summary['n_slides']}")
    print(f"Significant (p<0.05): graph {summary['n_graph_significant']}/{summary['n_slides']}, "
          f"naive {summary['n_naive_significant']}/{summary['n_slides']}")

    out = dict(per_slide=rows, summary=summary,
               meta=stamp_run(cfg, {"module": "M24-naive-baseline"}))
    p = P.validation("naive_baseline.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
