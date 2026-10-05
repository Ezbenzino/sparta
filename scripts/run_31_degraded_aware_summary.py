#!/usr/bin/env python
"""run_31_degraded_aware_summary.py —— 全队列 vs ag_target 充足子集的对照汇总
=====================================================================
CSCC10/14/15/16 这 4 张切片在第一代 ST 平台上 CD274/PDCD1LG2 命中 0–1 个，
B_mAb 的抗原吸收项被中性 0.5 填充（见 scores_from_adata）。这不影响 B_cell、
ECM、crosslink、hypoxia 等其他分量，但 B_mAb 本身少了抗原介导的结合消耗。

本脚本从已落盘的 5 个验证 JSON 中提取关键数字，同时报：
  - ALL        : 全部 19 张主队切片
  - AG_OK      : 剔除 ag_target-degraded 后的 15 张（MEL01 只缺 efflux，保留）
如果两套结论方向一致、量级接近，说明"两屏障耦合"不依赖抗原吸收项——
它由共享 ECM 通道驱动，这恰好是论文的核心主张。

输出：results/validation/degraded_aware_summary.json + 打印对照表
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

V = ROOT / "results" / "validation"

# 哪些切片因 ag_target 缺失而 B_mAb 退化（efflux-only 不影响 rho_partial）
AG_TARGET_DEGRADED = {"CSCC10", "CSCC14", "CSCC15", "CSCC16"}


def _stats(values, pvals=None):
    v = np.asarray(values, dtype=float)
    out = dict(
        n=int(len(v)),
        median=float(np.median(v)),
        q25=float(np.percentile(v, 25)),
        q75=float(np.percentile(v, 75)),
        n_pos=int((v > 0).sum()),
        n_neg=int((v < 0).sum()),
    )
    if pvals is not None:
        p = np.asarray(pvals, dtype=float)
        out["n_sig_p05"] = int((p < 0.05).sum())
        out["n_pos_sig"] = int(((v > 0) & (p < 0.05)).sum())
        out["n_neg_sig"] = int(((v < 0) & (p < 0.05)).sum())
    return out


def compare(name, per, value_key, pkey=None):
    all_slides = sorted(per.keys())
    ag_ok = [s for s in all_slides if s not in AG_TARGET_DEGRADED]
    all_vals = [per[s][value_key] for s in all_slides]
    ag_vals = [per[s][value_key] for s in ag_ok]
    all_p = [per[s][pkey] for s in all_slides] if pkey else None
    ag_p = [per[s][pkey] for s in ag_ok] if pkey else None
    s_all = _stats(all_vals, all_p)
    s_ag = _stats(ag_vals, ag_p)
    print(f"\n=== {name} ===")
    print(f"  ALL (n={s_all['n']}): median={s_all['median']:+.3f} "
          f"[{s_all['q25']:+.3f}, {s_all['q75']:+.3f}]  "
          f"pos/neg={s_all['n_pos']}/{s_all['n_neg']}", end="")
    if pkey:
        print(f"  sig(p<0.05) pos/neg={s_all['n_pos_sig']}/{s_all['n_neg_sig']}", end="")
    print()
    print(f"  AG_OK (n={s_ag['n']}): median={s_ag['median']:+.3f} "
          f"[{s_ag['q25']:+.3f}, {s_ag['q75']:+.3f}]  "
          f"pos/neg={s_ag['n_pos']}/{s_ag['n_neg']}", end="")
    if pkey:
        print(f"  sig(p<0.05) pos/neg={s_ag['n_pos_sig']}/{s_ag['n_neg_sig']}", end="")
    print()
    return {"all": s_all, "ag_ok": s_ag}


def main():
    out = {"ag_target_degraded_slides": sorted(AG_TARGET_DEGRADED),
           "note": "MEL01 只缺 efflux，不影响 rho_partial（control 是纯图距离），保留在两组中"}

    # --- R3 decoupling ---
    dec = json.load(open(V / "decoupling.json", encoding="utf-8"))
    out["decoupling_rho_partial"] = compare(
        "R3 decoupling: B_cell vs B_mAb (control=d_vessel)",
        dec, "rho_partial", "p_partial")
    out["decoupling_enrichment"] = compare(
        "R3 dissociation-zone enrichment vs chance",
        dec, "enrichment_vs_chance_r")

    # --- spatial null ---
    sn = json.load(open(V / "spatial_null_check.json", encoding="utf-8"))["per_slide"]
    out["spatial_null_real_rho"] = compare(
        "M3 spatial-null: real rho_partial",
        sn, "real_rho_partial")
    ps_all = [sn[s]["empirical_p_one_sided"] for s in sn]
    ps_ag = [sn[s]["empirical_p_one_sided"] for s in sn if s not in AG_TARGET_DEGRADED]
    print(f"\n=== M3 spatial-null: empirical p ===")
    print(f"  ALL  (n={len(ps_all)}): median p={np.median(ps_all):.3f}, "
          f"p<0.05 in {sum(p<0.05 for p in ps_all)}/{len(ps_all)}")
    print(f"  AG_OK(n={len(ps_ag)}): median p={np.median(ps_ag):.3f}, "
          f"p<0.05 in {sum(p<0.05 for p in ps_ag)}/{len(ps_ag)}")
    out["spatial_null_p"] = {
        "all": {"n": len(ps_all), "median_p": float(np.median(ps_all)),
                "n_sig": int(sum(p < 0.05 for p in ps_all))},
        "ag_ok": {"n": len(ps_ag), "median_p": float(np.median(ps_ag)),
                  "n_sig": int(sum(p < 0.05 for p in ps_ag))},
    }

    # --- biological validation ---
    bv = json.load(open(V / "biological_validation.json", encoding="utf-8"))["per_slide"]
    out["bc_vs_cd8"] = compare(
        "Bio: B_cell vs CD8T (expect NEGATIVE)",
        bv, "bc_vs_cd8_rho", "bc_vs_cd8_p")
    out["bm_vs_prolif"] = compare(
        "Bio: B_mAb vs Proliferation",
        bv, "bm_vs_prolif_rho", "bm_vs_prolif_p")

    # --- S2 matched selection ---
    s2 = json.load(open(V / "s2_matched_selection.json", encoding="utf-8"))["per_slide"]
    out["s2_ratio"] = compare(
        "S2 matched: targeted/scatter ratio (>1 = continuous arc works)",
        s2, "ratio_vs_in_cut_matched", "p_vs_in_cut_matched")

    path = V / "degraded_aware_summary.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
