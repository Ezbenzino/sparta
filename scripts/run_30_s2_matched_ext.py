#!/usr/bin/env python
"""External S2 matched-selection counterfactual (2026-10-03 supplement).

Same logic as run_28_s2_matched.py (matched selection pressure on the
scattered control: 8 candidates per group, argmin kept, 100 groups) but run
only on the external validation slides and written to an independent summary
file so the main-cohort s2_matched_selection.json stays untouched.

Output: results/validation/s2_matched_ext.json
        + results/counterfactual/{slide}.s2matched.json (per-slide, same
          naming convention as run_28)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402
from sparta.counterfactual import s2_ring_breaking  # noqa: E402

GROUPS = 100          # 与 run_28 一致：每组 8 个分散候选、取剩余屏障最小者
K_FRAC = (0.20,)      # 主检验分数（稿件预注册口径）


def main():
    ap = argparse.ArgumentParser(description="外部验证 S2 匹配选择敏感性")
    ap.add_argument("--slides", nargs="+",
                    default=["BRCA01", "BRCA02", "LN01"])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cf = cfg["counterfactual"]
    pmap = patient_map(P)

    print(f"External S2 matched-selection  N_groups={GROUPS}, k_frac={K_FRAC}\n")
    summary = {}
    for sid in args.slides:
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_nodes(nodes)
        coords = np.asarray(coords_from_nodes(nodes), float)
        cfg_cell = cfg["barrier"]["b_cell"]
        s2 = s2_ring_breaking(A, S["ecm"], S["caf"], source, sink, cfg_cell,
                              coords=coords, k_list=tuple(cf["s2"]["k_list"]),
                              k_frac=K_FRAC, n_rand=cf["s2"]["n_rand"],
                              seed=cfg["seed"], low_q=cf["s2"]["low_q"],
                              match_selection=True, n_match_groups=GROUPS)
        k = int(next(iter(s2["per_k"])))
        r = s2["per_k"][k]
        summary[sid] = dict(
            patient=pmap.get(sid, "?"),
            cohort="EXT",
            n_cut_nodes=int(s2["n_cut_nodes"]),
            k=int(k),
            targeted=float(r["targeted"]),
            rand_in_cut_matched=float(r["rand_in_cut_mean"]),
            ratio_vs_in_cut_matched=float(r["ratio_vs_in_cut"]),
            p_vs_in_cut_matched=float(r["p_vs_in_cut"]),
            n_groups=int(r["n_rand_in_cut"]),
            match_selection=bool(r["match_selection"]),
        )
        save_json(P.counterfactual(sid).with_name(f"{sid}.s2matched.json"),
                  {"slide": sid, "s2_matched": s2,
                   "meta": stamp_run(cfg, {"module": "M30-ext-s2-matched",
                                           "slide": sid})})
        print(f"{sid:<7} k={k:>3}  targeted={r['targeted']:.4f}  "
              f"scatter(matched)={r['rand_in_cut_mean']:.4f}  "
              f"ratio={r['ratio_vs_in_cut']:.3f}x  p={r['p_vs_in_cut']:.3f}")

    ratios = np.array([v["ratio_vs_in_cut_matched"] for v in summary.values()])
    ps = np.array([v["p_vs_in_cut_matched"] for v in summary.values()])
    out = dict(
        per_slide=summary,
        summary=dict(
            n_slides=len(summary),
            median_ratio=float(np.median(ratios)),
            n_ratio_above_1=int((ratios > 1.0).sum()),
            n_significant_raw=int((ps < 0.05).sum()),
            p_values=[float(x) for x in ps],
        ),
        note="外部独立验证 S2 匹配选择：与 run_28 完全同口径（8 候选取最优、100 组）",
    )
    # BH 校正（与 run_28 一致的小样本步骤）
    ps_sorted = np.sort(ps)
    m = len(ps_sorted)
    bh_pass = [i for i, p in enumerate(ps_sorted) if p <= 0.05 * (i + 1) / m]
    n_bh = len(bh_pass)
    out["summary"]["n_significant_bh"] = n_bh

    p = P.validation("s2_matched_ext.json")
    save_json(p, out)
    print(f"\n=== 外部 S2 匹配选择汇总 ===")
    print(f"ratio>1: {out['summary']['n_ratio_above_1']}/{out['summary']['n_slides']}"
          f"，中位 {out['summary']['median_ratio']:.3f}x，"
          f"raw p<0.05: {out['summary']['n_significant_raw']}，"
          f"BH p<0.05: {n_bh}/{out['summary']['n_slides']}")
    print(f"已写出 {p}")


if __name__ == "__main__":
    raise SystemExit(main())
