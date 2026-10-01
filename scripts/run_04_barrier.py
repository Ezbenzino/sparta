#!/usr/bin/env python
"""
run_04_barrier.py —— M4 三分量屏障算子（核心）
============================================

输入：{slide_id}.scored.h5ad + {slide_id}.graph.npz
输出：{slide_id}.barrier.npz + {slide_id}.mincut.json + 屏障线叠加图
上游模块：run_03_graph.py
下游模块：run_05_counterfactual.py / run_06_validate.py

计算 B_cell（最小割）、B_mAb（有效阻抗+结合位点屏障）、B_meta（代谢庇护），
并输出最小割"封锁线"叠加图——本课题最具说服力的图件。

在跑这一步之前，请确认 `python tests/test_barrier.py` 全绿。
算子写错了在真实数据上察觉不到，只有合成图测试能发现。

用法
----
    python scripts/run_04_barrier.py --slide MEL01
    python scripts/run_04_barrier.py --slide MEL01 --r-nm 0.5   # 小分子对照
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed, stamp_run  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="M4 三分量屏障算子",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slide", required=True)
    ap.add_argument("--r-nm", type=float, default=None,
                    help="分子流体力学半径。不传则用 config（IgG=5.5）。0.5 为小分子对照")
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy 读取 .scored.h5ad。仅验证算子请跑 scripts/run_00_demo.py。")

    from sparta.barrier import (compute_b_cell, compute_b_mab, compute_b_meta,
                                save_barriers, scores_from_adata)

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    adata = sc.read_h5ad(P.scored(args.slide))
    A, D, source, sink, vessel, gmeta = load_graph(P.graph(args.slide))
    S = scores_from_adata(adata)
    bcfg = cfg["barrier"]

    print(f"[M4] {args.slide}：{A.shape[0]} 节点｜源 {len(source)}｜汇 {len(sink)}｜"
          f"血管 {len(vessel)}")

    # --- B_cell ---
    bc = compute_b_cell(A, S["ecm"], S["caf"], source, sink, **bcfg["b_cell"])
    print(f"[M4] B_cell = {bc['b_cell']:.4f}（最大流 {bc['max_flow']:.4f}，"
          f"最小割 {len(bc['cut_edges'])} 条边 / {len(bc['cut_nodes'])} 个节点）")

    # --- B_mAb ---
    mab_cfg = dict(bcfg["b_mab"])
    if args.r_nm is not None:
        mab_cfg["r_nm"] = args.r_nm
    bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **mab_cfg)
    n_unreach = int((~bm["reachable"]).sum())
    print(f"[M4] B_mAb (r={mab_cfg['r_nm']} nm)：均值 {np.nanmean(bm['b_mab']):.3f}，"
          f"范围 {np.nanmin(bm['b_mab']):.3f}–{np.nanmax(bm['b_mab']):.3f}")
    if n_unreach:
        print(f"[M4] 注意：{n_unreach} 个节点与血管不连通，其 B_mAb 已赋为最大值并标记，"
              f"统计时请按 reachable 掩码处理")

    # --- B_meta ---
    # D 必须传：不传就退回"跳数 × spacing_um"，其偏差随平台阵列几何反号
    # （Visium +11%、第一代 ST −11%），会在队列比较里冒充生物学差异。
    bt = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                        D=D, **bcfg["b_meta"])
    print(f"[M4] B_meta：均值 {np.nanmean(bt['b_meta']):.3f}；"
          f"距血管 {np.nanmin(bt['d_vessel_um']):.0f}–{np.nanmax(bt['d_vessel_um']):.0f} μm"
          f"（口径 {bt['d_vessel_mode']}）")
    if bt["d_vessel_mode"] != "weighted":
        print("[M4] ⚠ d_vessel_um 退回了跳数近似，请检查 graph.npz 里的距离矩阵")

    res = {**bc, **bm, **bt}
    save_barriers(P.barrier(args.slide), res, P.mincut(args.slide))
    save_json(P.interim / f"{args.slide}.barrier_meta.json",
              stamp_run(cfg, {"module": "M4", "slide": args.slide,
                              "b_cell": bc["b_cell"], "max_flow": bc["max_flow"],
                              "r_nm": mab_cfg["r_nm"], "n_unreachable": n_unreach,
                              "d_vessel_mode": bt["d_vessel_mode"]}))
    print(f"[M4] 已写出 {P.barrier(args.slide)}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sparta.viz import (plot_barrier_field, plot_barrier_line,
                                plot_decoupling, setup_cjk_font)
        setup_cjk_font()
        coords = np.asarray(adata.obsm["spatial_um"], float)
        fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
        plot_barrier_line(coords, bc["cut_edges"], node_color=S["caf"],
                          source=source, sink=sink, ax=ax[0],
                          title=f"{args.slide} 最小割封锁线 (B_cell={bc['b_cell']:.3f})")
        plot_barrier_field(coords, bm["b_mab"], ax=ax[1],
                           title=f"B_mAb (r={mab_cfg['r_nm']} nm)", cbar_label="B_mAb")
        plot_decoupling(0.5 * (S["ecm"] + S["caf"]), bm["b_mab"], ax=ax[2])
        fig.tight_layout()
        p = P.figure(f"barrier_{args.slide}.png")
        fig.savefig(p, dpi=140); plt.close(fig)
        print(f"[M4] 已保存 {p}")


if __name__ == "__main__":
    main()
