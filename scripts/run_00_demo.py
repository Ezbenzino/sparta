#!/usr/bin/env python
"""
run_00_demo.py —— 在合成数据上跑通全流程（不需要任何真实数据）
================================================================

输入：无（数据由 sparta.synthetic 现场生成）
输出：results/figures/demo_*.png
上游依赖：无
下游模块：这是 run_01–run_06 的预演，跑通它说明环境装对了

这个脚本存在的意义
------------------
在你还没拿到任何真实切片之前，先用它确认三件事：
  1. 环境装对了（跑得动）
  2. 三个屏障算子在"已知正确答案"的图上行为正确
  3. 你看得懂每一步在算什么、输出长什么样

强烈建议在开始处理真实数据之前，先把这个脚本从头到尾读一遍并跑一次。
它只依赖 numpy·scipy·networkx·matplotlib，不需要 scanpy。

用法
----
    python scripts/run_00_demo.py                 # 默认配置
    python scripts/run_00_demo.py --n 31 --fast   # 更大的网格 + 减少置换次数
    python scripts/run_00_demo.py --help
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sparta.barrier import (compute_b_cell, compute_b_cell_field,  # noqa: E402
                            compute_b_mab, compute_b_meta)
from sparta.counterfactual import s1_spatial_permutation, s2_ring_breaking, s3_size_scan  # noqa: E402
from sparta.graph import graph_report  # noqa: E402
from sparta.synthetic import (make_dissociated_grid, make_ring_grid,  # noqa: E402
                              make_scattered_grid)
from sparta.validate import decoupling_stats, dissociation_drivers  # noqa: E402

CFG_CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)
CFG_MAB = dict(g0=1.0, lam=3.0, xi0_nm=20.0, beta=3.0, kd_eff=0.5, kappa_w=1.0)
CFG_META = dict(spacing_um=100.0, d0_um=300.0, tau_um=100.0)


def _hdr(t):
    print("\n" + "=" * 66)
    print(t)
    print("=" * 66)


def main():
    ap = argparse.ArgumentParser(
        description="在合成数据上跑通 SPARTA 全流程，验证环境与算子行为",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--n", type=int, default=25, help="合成网格边长（节点数 = n²）")
    ap.add_argument("--width", type=float, default=2.0, help="屏障带厚度（以 spot 为单位）")
    ap.add_argument("--fast", action="store_true", help="减少置换/随机次数，快速过一遍")
    ap.add_argument("--outdir", default="results/figures", help="图件输出目录")
    ap.add_argument("--no-figures", action="store_true", help="只跑数值，不出图")
    args = ap.parse_args()

    n_perm = 60 if args.fast else 300
    n_rand = 30 if args.fast else 120
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- 步骤 1
    _hdr("步骤 1｜生成两张合成切片：环带 vs 散在（细胞组成完全相同，只有排布不同）")
    ring = make_ring_grid(n=args.n, width_spots=args.width, gap_spots=0, seed=0)
    scat = make_scattered_grid(n=args.n, width_spots=args.width, seed=0)
    rep = graph_report(ring.A)
    print(f"  节点 {rep['n_nodes']}，边 {rep['n_edges']}，平均度 {rep['mean_degree']:.1f}，"
          f"连通分量 {rep['n_components']}")
    print(f"  源（免疫入口）{len(ring.source)} 个，汇（瘤巢核心）{len(ring.sink)} 个")
    print(f"  高阻抗节点：环带 {len(ring.blocked_idx)} 个 vs 散在 {len(scat.blocked_idx)} 个")
    assert len(ring.blocked_idx) == len(scat.blocked_idx), "对照的组成不一致，无法比较"

    # ---------------------------------------------------------------- 步骤 2
    _hdr("步骤 2｜B_cell：源汇最小割 —— 屏障是拓扑性质吗？")
    rc = compute_b_cell(ring.A, ring.ecm, ring.caf, ring.source, ring.sink, **CFG_CELL)
    sc_ = compute_b_cell(scat.A, scat.ecm, scat.caf, scat.source, scat.sink, **CFG_CELL)
    print(f"  环带 B_cell = {rc['b_cell']:.4f}   (最大流 {rc['max_flow']:.3f}，"
          f"最小割边 {len(rc['cut_edges'])} 条)")
    print(f"  散在 B_cell = {sc_['b_cell']:.4f}   (最大流 {sc_['max_flow']:.3f})")
    print(f"  -> 比值 {rc['b_cell'] / sc_['b_cell']:.1f}x。同样多的阻力细胞，"
          f"连成一圈就挡得住，散开就挡不住。")

    # ---------------------------------------------------------------- 步骤 3
    _hdr("步骤 3｜B_mAb：尺寸依赖有效阻抗 —— 大分子和小分子面对同一道屏障吗？")
    igg = compute_b_mab(ring.A, ring.ecm, ring.crosslink, ring.ag_target,
                        ring.vessel, r_nm=5.5, **CFG_MAB)
    sml = compute_b_mab(ring.A, ring.ecm, ring.crosslink, ring.ag_target,
                        ring.vessel, r_nm=0.5, **CFG_MAB)
    core = ring.sink
    print(f"  瘤巢核心 B_mAb（IgG, r=5.5nm）  = {np.nanmean(igg['b_mab'][core]):8.3f}")
    print(f"  瘤巢核心 B_mAb（小分子, r=0.5nm）= {np.nanmean(sml['b_mab'][core]):8.3f}")
    print(f"  -> 同一张图、同一道基质，IgG 基本被挡在外面，小分子畅通无阻。")

    # ---------------------------------------------------------------- 步骤 4
    _hdr("步骤 4｜B_meta：血管图扩散距离 × 代谢状态")
    meta = compute_b_meta(ring.A, ring.hypoxia, ring.proliferation, ring.efflux,
                          ring.vessel, **CFG_META)
    print(f"  距血管距离范围 {np.nanmin(meta['d_vessel_um']):.0f}–"
          f"{np.nanmax(meta['d_vessel_um']):.0f} μm")
    print(f"  B_meta 范围 {np.nanmin(meta['b_meta']):.3f}–{np.nanmax(meta['b_meta']):.3f}")

    # ---------------------------------------------------------------- 步骤 5
    _hdr("步骤 5｜解耦性：这是本课题的中心主张")
    # B_cell 是切片级标量，逐 spot 版本用"局部阻力"作为代理量。
    # 真实分析中同样需要这样处理，并在方法里说明这是代理量而非严格的逐点最小割。
    b_cell_proxy = compute_b_cell_field(ring.A, ring.ecm, ring.caf, ring.source,
                                        **CFG_CELL)["b_cell_field"]
    dec_ring = decoupling_stats(b_cell_proxy, igg["b_mab"], q=0.75,
                                control=meta["d_vessel_um"])
    print(f"  [环带图] 原始 ρ = {dec_ring['rho']:.3f}｜"
          f"校正几何后 ρ = {dec_ring['rho_partial']:.3f}｜"
          f"解离区 {dec_ring['frac_discordant_r'] * 100:.1f}%")
    print("  注意两件事：")
    print("  (1) 原始 ρ 很高，是因为两个屏障都是'从血管出发的累积代价'，天然共享")
    print("      一个很强的几何成分（离血管越远两者都越大）。不校正就报告相关性，")
    print("      等于在测深度而不是测机制。论文里必须报告校正后的偏相关。")
    print("  (2) 环带图里 ECM/CAF/交联度是同一个数组，构造上不可能解离 —— 阴性对照。")

    # 换一张"已知存在解离"的合成切片：ECM 环带开大缺口（T 细胞进得去），
    # 但瘤巢外缘高抗原（抗体被结合位点屏障消耗掉）。即临床上的"热而无效"。
    diss = make_dissociated_grid(n=args.n, seed=0)
    bc_d = compute_b_cell(diss.A, diss.ecm, diss.caf, diss.source, diss.sink, **CFG_CELL)
    bm_d = compute_b_mab(diss.A, diss.ecm, diss.crosslink, diss.ag_target,
                         diss.vessel, r_nm=5.5, **CFG_MAB)
    bcf_d = compute_b_cell_field(diss.A, diss.ecm, diss.caf, diss.source,
                                 **CFG_CELL)["b_cell_field"]
    meta_d = compute_b_meta(diss.A, diss.hypoxia, diss.proliferation, diss.efflux,
                            diss.vessel, **CFG_META)
    dec_d = decoupling_stats(bcf_d, bm_d["b_mab"], q=0.75,
                             control=meta_d["d_vessel_um"])
    print()
    print(f"  [解离图] B_cell = {bc_d['b_cell']:.4f}（远低于环带图的 {rc['b_cell']:.2f}，"
          f"说明 T 细胞进得去）")
    print(f"  [解离图] 原始 ρ = {dec_d['rho']:.3f}｜"
          f"校正几何后 ρ = {dec_d['rho_partial']:.3f}｜"
          f"解离区 {dec_d['frac_discordant_r'] * 100:.1f}%")
    print("  -> 这个比例就是本课题最有临床意义的数字：这类区域 T 细胞进得去、")
    print("     抗体进不去，在现有框架下会被判为'热而无效'，本框架给出了")
    print("     可测量的机制解释。")
    print()
    dd = dissociation_drivers(diss.A, diss.ecm, diss.caf, diss.crosslink,
                              diss.ag_target, diss.source, diss.vessel,
                              cfg_cell=CFG_CELL, cfg_mab=CFG_MAB)
    print("  [驱动分解] B_mAb 的变异来自哪里：")
    print(f"    共享 ECM 通道（同时挡细胞和分子）  {dd['frac_shared_ecm']*100:5.1f}%")
    print(f"    抗原 BSB（只挡分子，不挡细胞）     {dd['frac_antigen']*100:5.1f}%")
    print(f"    交联/尺寸排阻（只挡分子）          {dd['frac_crosslink']*100:5.1f}%")
    print(f"    解离潜力比 = {dd['dissociation_potential']:.3f}")
    print()
    print("  ⚠ 这是本项目最重要的一个诊断。解离潜力比 > 1 才有望在全片尺度")
    print("    观察到解离。当前合成设定下只有 0.04——ECM 解释了 96% 的变异，")
    print("    两个屏障主要由同一个因素驱动。请务必读 README 的『已知的核心风险』")
    print("    一节，并在建库阶段就对每张真实切片跑一遍这个诊断。")

    # ---------------------------------------------------------------- 步骤 6
    _hdr("步骤 6｜S1 空间重排对照 —— 屏障是不是细胞比例的复杂重写？")
    s1 = s1_spatial_permutation(ring.A, ring.ecm, ring.caf, ring.source, ring.sink,
                                CFG_CELL, n_perm=n_perm, seed=1, mode="fixed")
    print(f"  真实 B_cell = {s1['b_real']:.4f}")
    print(f"  重排零分布  = {s1['null'].mean():.4f} ± {s1['null'].std():.4f} "
          f"（{s1['n_perm']} 次）")
    print(f"  z = {s1['z']:.1f}，经验 p = {s1['p_emp']:.4f}")
    print("  -> 组成完全没变，只是把细胞随机搬了家，屏障就没了。")

    # ---------------------------------------------------------------- 步骤 7
    _hdr("步骤 7｜S2 环带断裂 —— 为什么非要用空间数据？")
    s2 = s2_ring_breaking(ring.A, ring.ecm, ring.caf, ring.source, ring.sink,
                          CFG_CELL, coords=ring.coords, k_list=(5, 10),
                          n_rand=n_rand, seed=1)
    print(f"  基线 B_cell = {s2['b0']:.4f}，最小割节点 {s2['n_cut_nodes']} 个")
    for k, r in sorted(s2["per_k"].items()):
        print(f"  k={k:>2}｜连续缺口 {r['targeted']:.4f} ｜ 屏障内分散 "
              f"{r['rand_in_cut_mean']:.4f} ｜ 全局随机 {r['rand_global_mean']:.4f}")
        print(f"       剩余屏障比 vs 分散 {r['ratio_vs_in_cut']:.2f}x，p = {r['p_vs_in_cut']:.3f}")
    print("  -> 同样多的屏障物质，拿掉连续的一段 vs 到处拿掉一点，效果完全不同。")
    print("     这个结论只有空间数据能给出，bulk 数据完全无法企及。")

    # ---------------------------------------------------------------- 步骤 8
    _hdr("步骤 8｜S3 分子尺寸扫描 —— 两个屏障凭什么说是解耦的？")
    s3 = s3_size_scan(ring.A, ring.ecm, ring.crosslink, ring.ag_target,
                      ring.vessel, CFG_MAB, core_idx=core)
    for r, m in zip(s3["radii_nm"], s3["mean_core"]):
        bar = "█" * max(int((m - min(s3["mean_core"])) / 1.2), 0)
        print(f"  r = {r:>4.1f} nm  核心区 B_mAb = {m:8.3f}  {bar}")
    print(f"  -> IgG 相对小分子的核心区屏障增幅 = {s3['igg_vs_small_core_delta']:.2f}")

    # ---------------------------------------------------------------- 出图
    if not args.no_figures:
        _hdr("出图")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sparta.viz import (setup_cjk_font, plot_barrier_field, plot_barrier_line,
                                plot_permutation_null, plot_ring_breaking,
                                plot_size_scan, plot_decoupling)

        setup_cjk_font(verbose=True)

        fig, axes = plt.subplots(2, 3, figsize=(16.5, 10.5))
        plot_barrier_line(ring.coords, rc["cut_edges"], node_color=ring.ecm,
                          source=ring.source, sink=ring.sink, ax=axes[0, 0],
                          title=f"环带：最小割封锁线 (B_cell={rc['b_cell']:.3f})")
        plot_barrier_line(scat.coords, sc_["cut_edges"], node_color=scat.ecm,
                          source=scat.source, sink=scat.sink, ax=axes[0, 1],
                          title=f"散在对照 (B_cell={sc_['b_cell']:.3f})")
        plot_barrier_field(ring.coords, igg["b_mab"], ax=axes[0, 2],
                           title="B_mAb (IgG, r=5.5 nm)", cbar_label="B_mAb")
        plot_permutation_null(s1, ax=axes[1, 0])
        plot_ring_breaking(s2, ax=axes[1, 1])
        plot_size_scan(s3, ax=axes[1, 2])
        fig.tight_layout()
        p1 = out / "demo_overview.png"
        fig.savefig(p1, dpi=140); plt.close(fig)

        fig2, ax2 = plt.subplots(1, 3, figsize=(16, 5))
        plot_decoupling(b_cell_proxy, igg["b_mab"], ax=ax2[0])
        ax2[0].set_title(f"环带图：构造上不可解离 ({dec_ring['frac_discordant']*100:.0f}%)")
        plot_decoupling(bcf_d, bm_d["b_mab"], ax=ax2[1])
        ax2[1].set_title(f"解离图：T细胞进得去/抗体进不去 ({dec_d['frac_discordant']*100:.0f}%)")
        plot_barrier_field(ring.coords, meta["b_meta"], ax=ax2[2],
                           title="B_meta（代谢庇护）", cbar_label="B_meta")
        fig2.tight_layout()
        p2 = out / "demo_decoupling.png"
        fig2.savefig(p2, dpi=140); plt.close(fig2)
        print(f"  已保存 {p1}")
        print(f"  已保存 {p2}")

    _hdr("完成")
    print("  以上全部在合成数据上完成，没有使用任何真实数据，也没有产生任何研究结论。")
    print("  下一步：按 README 的『真实数据流程』准备切片，然后依次跑 run_01 – run_06。")


if __name__ == "__main__":
    main()
