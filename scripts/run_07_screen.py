#!/usr/bin/env python
"""
run_07_screen.py —— 早期决策点：解离潜力筛查
================================================

输入：若干张已跑完 run_03_graph 的切片
输出：results/validation/screen_decision.json + screen_decision.png
上游模块：run_03_graph.py（不需要先跑 run_04）
下游模块：无 —— 这是一个**决策节点**，它的输出决定你论文的叙事走哪条分支

为什么要在建全库之前跑这个
--------------------------
本课题的中心主张是"T 细胞迁移屏障与抗体传质屏障可以解离"。但两者共享
同一个 ECM 项——ECM 既进了最小割的边容量，也进了扩散的边电导——这条耦合
通道无法用统计方法校正掉（详见 README 的「已知的核心风险」）。

只有**只影响抗体不影响细胞**的因素（靶抗原、交联度、分子尺寸）才能产生解离。
`dissociation_drivers()` 量化这些因素相对于共享 ECM 的贡献：

    解离潜力比 = (抗原贡献 + 交联贡献) / ECM 贡献

**> 1 才有望在全片尺度观察到解离。**

拿到 3–5 张切片就该跑这个，不要等到 15 张全建完。它只需几分钟，
却能在第三周告诉你叙事该怎么写——而不是第八个月。

用法
----
    python scripts/run_07_screen.py --slides MEL01 MEL02 SCC01
    python scripts/run_07_screen.py --slides MEL01 --r-nm 5.5 --no-plot
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed, stamp_run  # noqa: E402

THRESH_GOOD = 1.0      # 潜力比 > 1：主线叙事可行
THRESH_WEAK = 0.3      # 0.3–1：边缘，需谨慎并加强抗原/交联维度


def main():
    ap = argparse.ArgumentParser(
        description="早期决策点：解离潜力筛查",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=None)
    ap.add_argument("--synthetic", action="store_true",
                    help="用合成切片演示决策流程（不需要任何真实数据）。"
                         "建议在拿到真实数据之前先跑一次，看看输出长什么样")
    ap.add_argument("--r-nm", type=float, default=None, help="分子半径，默认取 config")
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if not (args.slides or args.synthetic):
        ap.error("需要 --slides 或 --synthetic 之一")

    if not args.synthetic:
        try:
            import scanpy as sc
        except ImportError:
            sys.exit("需要 scanpy 读取签名分数。pip install scanpy\n"
                     "想先看看输出长什么样，可以跑：python scripts/run_07_screen.py --synthetic")

    from sparta.barrier import (compute_b_cell_field, compute_b_mab,
                                compute_b_meta, scores_from_adata)
    from sparta.validate import decoupling_stats, dissociation_drivers

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    cfg_cell = cfg["barrier"]["b_cell"]
    cfg_mab = {k: v for k, v in cfg["barrier"]["b_mab"].items() if k != "r_nm"}
    r_nm = args.r_nm or cfg["barrier"]["b_mab"]["r_nm"]

    # 合成演示模式：现造三张性质不同的切片
    if args.synthetic:
        from sparta.synthetic import (make_dissociated_grid, make_ring_grid,
                                      make_scattered_grid)
        print("【合成演示模式】以下切片由 sparta.synthetic 现场生成，"
              "不含任何真实数据，也不构成任何研究结论。\n")
        synth = {
            "SYN_ring": make_ring_grid(n=25, width_spots=2.0, seed=0),
            "SYN_gap": make_ring_grid(n=25, width_spots=2.0, gap_spots=6, seed=0),
            "SYN_diss": make_dissociated_grid(n=25, seed=0),
        }
        args.slides = list(synth)

    rows = {}
    print(f"{'切片':<12}{'解离潜力':>10}{'ECM%':>8}{'抗原%':>8}{'交联%':>8}"
          f"{'ECM~交联':>10}{'校正ρ':>9}{'解离区%':>9}")
    print("-" * 78)

    for sid in args.slides:
        if args.synthetic:
            sl = synth[sid]
            A, source, sink, vessel = sl.A, sl.source, sl.sink, sl.vessel
            S = sl.scores
        else:
            try:
                adata = sc.read_h5ad(P.scored(sid))
                A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
            except FileNotFoundError as e:
                print(f"{sid:<12} 跳过：{e}")
                continue
            S = scores_from_adata(adata)

        dd = dissociation_drivers(A, S["ecm"], S["caf"], S["crosslink"],
                                  S["ag_target"], source, vessel,
                                  cfg_cell=cfg_cell, cfg_mab=cfg_mab, r_nm=r_nm)

        bcf = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg_cell)["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           r_nm=r_nm, **cfg_mab)["b_mab"]
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        dec = decoupling_stats(bcf, bm, q=cfg["validate"]["decouple_q"], control=dv)

        rows[sid] = {**{k: v for k, v in dd.items()},
                     "rho": dec["rho"], "rho_partial": dec["rho_partial"],
                     "frac_discordant_r": dec["frac_discordant_r"],
                     "n_spots": int(A.shape[0])}
        flag = "*" if dd["crosslink_collinear"] else " "
        print(f"{sid:<12}{dd['dissociation_potential']:>10.3f}"
              f"{dd['frac_shared_ecm']*100:>8.1f}{dd['frac_antigen']*100:>8.1f}"
              f"{dd['frac_crosslink']*100:>8.1f}"
              f"{dd['ecm_crosslink_corr']:>9.2f}{flag}"
              f"{dec['rho_partial']:>9.3f}{dec['frac_discordant_r']*100:>9.1f}")

    if not rows:
        sys.exit("没有任何切片可用。请先跑 run_01 → run_02 → run_03。")

    pots = np.array([r["dissociation_potential"] for r in rows.values()])
    med = float(np.median(pots))

    print("-" * 78)
    n_col = sum(1 for r in rows.values() if r["crosslink_collinear"])
    if n_col:
        print(f"* {n_col} 张切片的交联度与 ECM 共线（|r| >= 0.80），"
              f"其交联贡献已保守地并入共享桶——分不开就不声称分得开。")
    print(f"中位解离潜力比 = {med:.3f}   （n={len(pots)} 张切片）\n")
    print("=" * 76)
    if med > THRESH_GOOD:
        verdict = "GO"
        print("判定：GO —— 主线叙事可行")
        print("  只影响抗体的因素（抗原/交联）贡献超过共享 ECM，")
        print("  两道屏障有望在全片尺度上解离。按原计划建全库（15–20 张），")
        print("  重点报告 B_mAb 的独立贡献与增量价值。")
    elif med > THRESH_WEAK:
        verdict = "CAUTION"
        print("判定：CAUTION —— 边缘，需要加强")
        print("  解离信号存在但不占优势。建议：")
        print("  1) 优先纳入抗原表达异质性大的切片（先看各切片 Ag_target 的空间方差）；")
        print("  2) 在论文中同时报告解离与耦合两方面证据，不要只讲一边；")
        print("  3) S3 分子尺寸扫描的权重加大——它不依赖抗原，是更稳的解离证据。")
    else:
        verdict = "PIVOT"
        print("判定：PIVOT —— 建议调整叙事")
        print("  共享 ECM 主导了两个屏障，解离在全片尺度上大概率观察不到。")
        print("  这不是失败，是一个可发表的澄清性结论。建议把主张改为：")
        print()
        print('    "皮肤肿瘤中 T 细胞迁移屏障与抗体传质屏障高度耦合，')
        print('     且耦合来源是细胞外基质而非几何深度——因此针对基质的干预')
        print('     有望同时改善两类药物的递送。"')
        print()
        print("  支撑证据你已经有了：驱动分解（量化 ECM 的主导地位）、")
        print("  几何校正（证明耦合不是深度造成的）、S3 尺寸扫描（证明两者")
        print("  确实是不同的物理过程，只是被同一个因素驱动）。")
        print()
        print("  ⚠ 绝不要为了让解离出现而调 b_cell.b_ecm 或 b_mab.lam。")
        print("    这两个参数的相对大小只能依据文献确定。")
    print("=" * 76)

    out = dict(verdict=verdict, median_potential=med,
               per_slide=rows, meta=stamp_run(cfg, {"module": "screen", "r_nm": r_nm}))
    p = P.validation("screen_decision.json")
    save_json(p, out)
    print(f"\n已写出 {p}")

    if not args.no_plot and len(rows) >= 1:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sparta.viz import setup_cjk_font
        setup_cjk_font()

        sids = list(rows)
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        bottom = np.zeros(len(sids))
        for key, lab, c in (("frac_shared_ecm", "共享 ECM（同时挡两者）", "#8FA3B8"),
                            ("frac_antigen", "抗原 BSB（只挡抗体）", "#C0392B"),
                            ("frac_crosslink", "交联/尺寸（只挡抗体）", "#E67E22")):
            v = np.array([rows[s][key] for s in sids])
            ax[0].bar(sids, v, bottom=bottom, label=lab, color=c)
            bottom += v
        ax[0].set_ylabel("占 B_mAb 变异的比例")
        ax[0].set_title("驱动分解")
        ax[0].legend(fontsize=8)
        ax[0].tick_params(axis="x", rotation=45)

        ax[1].bar(sids, pots, color=["#2E8B57" if p > THRESH_GOOD else
                                     "#E67E22" if p > THRESH_WEAK else "#C0392B"
                                     for p in pots])
        ax[1].axhline(THRESH_GOOD, color="#2E8B57", ls="--", lw=1.2, label="GO 阈值 = 1.0")
        ax[1].axhline(THRESH_WEAK, color="#C0392B", ls=":", lw=1.2, label="PIVOT 阈值 = 0.3")
        ax[1].set_ylabel("解离潜力比")
        ax[1].set_title(f"决策：{verdict}（中位 {med:.3f}）")
        ax[1].legend(fontsize=8)
        ax[1].tick_params(axis="x", rotation=45)
        fig.tight_layout()
        pf = P.figure("screen_decision.png")
        fig.savefig(pf, dpi=140); plt.close(fig)
        print(f"已保存 {pf}")


if __name__ == "__main__":
    main()
