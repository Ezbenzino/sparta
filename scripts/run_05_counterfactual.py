#!/usr/bin/env python
"""
run_05_counterfactual.py —— M5 反事实实验
==========================================

输入：{slide_id}.scored.h5ad + {slide_id}.graph.npz
输出：results/counterfactual/{slide_id}.json + 图件
上游模块：run_04_barrier.py
下游模块：run_06_validate.py

三个实验，分别封堵三个致命质疑：
  S1 空间重排  -> "分数是不是细胞比例的复杂重写？"
  S2 环带断裂  -> "为什么非要用空间数据？"
  S3 尺寸扫描  -> "两个屏障凭什么说是解耦的？"

成本极低、说服力极高。主线跑通后应立即完成，不要拖到最后。

用法
----
    python scripts/run_05_counterfactual.py --slide MEL01
    python scripts/run_05_counterfactual.py --slide MEL01 --n-perm 1000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed, stamp_run  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="M5 反事实实验",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slide", required=True)
    ap.add_argument("--n-perm", type=int, default=None, help="S1 置换次数")
    ap.add_argument("--n-rand", type=int, default=None, help="S2 随机对照次数")
    ap.add_argument("--k-absolute", action="store_true",
                    help="S2 强制用 config 的绝对 k_list（只用于复现 2026-08-26 之前的旧结果）")
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy 读取 .scored.h5ad。仅验证实验逻辑请跑 scripts/run_00_demo.py。")

    from sparta.barrier import scores_from_adata
    from sparta.counterfactual import s1_spatial_permutation, s2_ring_breaking, s3_size_scan

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cf = cfg["counterfactual"]

    adata = sc.read_h5ad(P.scored(args.slide))
    A, D, source, sink, vessel, _ = load_graph(P.graph(args.slide))
    S = scores_from_adata(adata)
    coords = np.asarray(adata.obsm["spatial_um"], float)
    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}

    out = {"slide": args.slide, "meta": stamp_run(cfg, {"module": "M5"})}

    # ---- S1 ----
    n_perm = args.n_perm or cf["s1"]["n_perm"]
    print(f"[S1] 空间重排对照（{n_perm} 次）...")
    out["s1"] = {}
    for mode in cf["s1"]["modes"]:
        kw = {}
        if mode == "follow":
            kw = dict(endothelial=S.get("ecm"), t_nk=S.get("caf"), malignant=S.get("ag_target"),
                      cfg_source_sink={k: cfg["source_sink"][k] for k in
                                       ("q_vessel", "q_immune_nbr", "q_malig", "q_core")})
            missing = []
            for name, col in (("endothelial", "Endothelial_n"), ("t_nk", "T_NK_n"),
                              ("malignant", "Malignant_n")):
                if col in adata.obs:
                    kw[name] = np.asarray(adata.obs[col].values, float)
                else:
                    missing.append(col)
            if missing:
                # 不要静默用 ecm/caf/ag_target 顶替——那会让 follow 模式测的东西
                # 与它声称测的东西完全不同，而且下游一点都看不出来。
                sys.exit(f"[S1] follow 模式需要 {missing}，但 {args.slide}.scored.h5ad 里没有。"
                         f"\n     请先重跑 run_02_score.py，或用 --modes fixed 只跑 fixed。")
        r = s1_spatial_permutation(A, S["ecm"], S["caf"], source, sink, cfg_cell,
                                   n_perm=n_perm, seed=cfg["seed"], mode=mode, **kw)
        print(f"     mode={mode}: 真实 {r['b_real']:.4f}｜零分布 {r['null'].mean():.4f}"
              f"±{r['null'].std():.4f}｜z={r['z']:.2f}｜p={r['p_emp']:.4f}")
        out["s1"][mode] = {k: v for k, v in r.items() if k != "null"}
        # 效应量用比值而不是 z：z 大可能只是零分布方差小。
        # 例：MEL01 fixed 的 z=14.9 看着很强，但真实 0.0308 / 零分布均值 0.0136
        # 只有 2.27x。论文里两个都要报，主报比值。
        nm = float(r["null"].mean())
        out["s1"][mode]["null_summary"] = dict(mean=nm, std=float(r["null"].std()),
                                               n=int(len(r["null"])))
        out["s1"][mode]["ratio_real_over_null"] = float(r["b_real"] / nm) if nm else None
        out["s1"][mode]["p_floor"] = 1.0 / (len(r["null"]) + 1)
        if mode == cf["s1"]["modes"][0]:
            s1_plot = r

    # ---- S2 ----
    n_rand = args.n_rand or cf["s2"]["n_rand"]
    # k 优先按割集比例取（真实切片必须这样，见 counterfactual.s2_ring_breaking 的说明）；
    # 只有 config 里没写 k_frac 时才回退到绝对 k_list。
    k_frac = cf["s2"].get("k_frac") if not args.k_absolute else None
    if k_frac:
        print(f"[S2] 环带断裂（k = 割集的 {k_frac}，随机对照 {n_rand} 次）...")
    else:
        print(f"[S2] 环带断裂（k={cf['s2']['k_list']} 绝对值，随机对照 {n_rand} 次）...")
    s2 = s2_ring_breaking(A, S["ecm"], S["caf"], source, sink, cfg_cell,
                          coords=coords, k_list=tuple(cf["s2"]["k_list"]),
                          k_frac=tuple(k_frac) if k_frac else None,
                          n_rand=n_rand, seed=cfg["seed"], low_q=cf["s2"]["low_q"])
    print(f"     割集 {s2['n_cut_nodes']} 节点（可动 {s2['n_cut_pool']} 个）")
    for k, r in sorted(s2["per_k"].items()):
        print(f"     k={k:>3}（占割集 {r['k_over_cut']*100:>4.1f}%）: 连续缺口 {r['targeted']:.4f}"
              f"｜屏障内分散 {r['rand_in_cut_mean']:.4f}"
              f"（{r['ratio_vs_in_cut']:.2f}x, p={r['p_vs_in_cut']:.3f}）"
              f"｜全局随机 {r['rand_global_mean']:.4f}")
    out["s2"] = s2

    # ---- S3 ----
    print(f"[S3] 分子尺寸扫描 {cf['s3']['radii_nm']} nm ...")
    s3 = s3_size_scan(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, cfg_mab,
                      radii_nm=tuple(cf["s3"]["radii_nm"]), core_idx=sink)
    for r_, m in zip(s3["radii_nm"], s3["mean_core"]):
        print(f"     r={r_:>5.1f} nm -> 瘤巢核心 B_mAb = {m:.3f}")
    out["s3"] = {k: v for k, v in s3.items() if k != "fields"}

    save_json(P.counterfactual(args.slide), out)
    print(f"[M5] 已写出 {P.counterfactual(args.slide)}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sparta.viz import (plot_permutation_null, plot_ring_breaking,
                                plot_size_scan, setup_cjk_font)
        setup_cjk_font()
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
        plot_permutation_null(s1_plot, ax=ax[0])
        plot_ring_breaking(s2, ax=ax[1])
        plot_size_scan(s3, ax=ax[2])
        fig.tight_layout()
        p = P.figure(f"counterfactual_{args.slide}.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        print(f"[M5] 已保存 {p}")


if __name__ == "__main__":
    main()
