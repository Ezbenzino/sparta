#!/usr/bin/env python
"""
run_28_s2_matched.py —— S2 匹配选择敏感性（2026-10-03 审查 P4）
=============================================================================

定向连续弧每片从 ≤8 个候选取"剩余屏障最小"的一段（argmin），而旧版对照 B
（屏障内分散）是均匀随机抽分散集，两臂候选搜索力度不对等，会把"连续优势"
系统性夸大。本脚本把相同的选择压力施加给对照 B：每组抽 8 个分散候选、
保留剩余屏障最小者，量化选择匹配后的效应量与 p。

- 只跑主检验 k_frac = 0.20（稿件预注册口径），组数 n_match_groups = 100，
  与主结果里对照 B 的 100 次同分辨率（p 下限 1/101）。
- 输出：results/validation/s2_matched_selection.json（汇总）
        + results/counterfactual/{slide}.s2matched.json（逐片，不覆盖主结果）
- 靶向弧的种子/候选顺序与 run_05 主跑一致（同一 seed），因此 targeted 值与
  主结果可逐片比对；只有对照 B 变了。

用法
----
    python scripts/run_28_s2_matched.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402
from sparta.counterfactual import s2_ring_breaking  # noqa: E402

GROUPS = 100          # 组数：每组 8 个分散候选、取剩余屏障最小者
K_FRAC = (0.20,)      # 主检验分数


def main():
    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cf = cfg["counterfactual"]
    pmap = patient_map(P)

    # 2026-10-05: primary sections only (the ledger now also lists the external arm as "ingested");
    # inputs come from the archived node tables, so Scanpy is not needed (run_35b verifies equality)
    sys.path.insert(0, str(Path(__file__).resolve().parent / "is_figures"))
    from isdata import PRIMARY_ORDER
    slides = list(PRIMARY_ORDER)

    print(f"S2 matched-selection sensitivity  N_groups={GROUPS}, k_frac={K_FRAC}\n")
    summary = {}
    for sid in slides:
        try:
            nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        filled: list[str] = []
        S = scores_from_nodes(nodes, missing_out=filled)
        coords = coords_from_nodes(nodes)
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
            cohort="CSCC" if sid.startswith("CSCC") else "MEL",
            n_cut_nodes=int(s2["n_cut_nodes"]),
            k=int(k),
            targeted=float(r["targeted"]),
            rand_in_cut_matched=float(r["rand_in_cut_mean"]),
            ratio_vs_in_cut_matched=float(r["ratio_vs_in_cut"]),
            p_vs_in_cut_matched=float(r["p_vs_in_cut"]),
            n_groups=int(r["n_rand_in_cut"]),
            match_selection=bool(r["match_selection"]),
            # S2 本身只用 ecm/caf/source/sink，不受 ag_target/efflux 填充影响；
            # 但记录下来便于排查"这张切片 M2 打分是否完整"。
            filled_keys=filled,
        )
        save_json(P.counterfactual(sid).with_name(f"{sid}.s2matched.json"),
                  {"slide": sid, "s2_matched": s2,
                   "meta": stamp_run(cfg, {"module": "M28-s2-matched", "slide": sid})})
        print(f"{sid:<8} k={k:>3}  targeted={r['targeted']:.4f}  "
              f"scatter(matched)={r['rand_in_cut_mean']:.4f}  "
              f"ratio={r['ratio_vs_in_cut']:.3f}x  p={r['p_vs_in_cut']:.3f}")

    ratios = np.array([v["ratio_vs_in_cut_matched"] for v in summary.values()])
    ps = np.array([v["p_vs_in_cut_matched"] for v in summary.values()])
    out = dict(
        per_slide=summary,
        summary=dict(
            n_slides=len(summary),
            n_ratio_above_1=int((ratios > 1.0).sum()),
            n_significant=int((ps < 0.05).sum()),
            median_ratio=float(np.median(ratios)),
            ratio_range=[float(np.min(ratios)), float(np.max(ratios))],
            groups=GROUPS,
            k_frac=list(K_FRAC),
            p_correction="(k+1)/(N+1)",
        ),
        meta=stamp_run(cfg, {"module": "M28-s2-matched"}),
    )
    save_json(P.validation("s2_matched_selection.json"), out)
    print(f"\n=== Summary (matched selection) ===")
    print(f"n={out['summary']['n_slides']}, ratio>1: "
          f"{out['summary']['n_ratio_above_1']}/{out['summary']['n_slides']}, "
          f"p<0.05: {out['summary']['n_significant']}/{out['summary']['n_slides']}, "
          f"median ratio {out['summary']['median_ratio']:.3f}")
    print(f"\nWrote {P.validation('s2_matched_selection.json')}")


if __name__ == "__main__":
    main()
