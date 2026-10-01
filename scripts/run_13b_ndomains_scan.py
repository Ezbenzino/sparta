#!/usr/bin/env python
"""
run_13b_ndomains_scan.py —— 空间域数量 n_domains 的敏感性扫描（补充材料）
=========================================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz + {sid}.mincut.json
输出：results/validation/ndomains_sensitivity.json
      results/figures/supp_ndomains_scan.png
上游模块：run_13_benchmark_tools.py（复用其 banksy_style_domains，口径完全一致）
下游模块：无（成文用，对应 Results R6 的局限性段落与补充材料 S5）

为什么单独一个脚本
------------------
run_13 的产物是读回合并的：换 --n-domains 直接重跑会把主结果
（n_domains=8 的 19 张）覆盖掉。本脚本**只写自己的 JSON**，
且在 k=8 处与 benchmark_tools.json 逐张对数（同种子同配方，必须一致）。

扫描设计（预先指定，写进 Methods）
----------------------------------
1. 固定网格：n_domains ∈ {4, 6, 8, 10, 12}，全部切片同一组；
2. 按切片大小缩放的一臂：n_scaled = clip(round(n_nodes/200), 4, 16)。
   规则的动机：每域期望约 200 个 spot——Visium 默认的 8 域除以
   MEL01 的 ~1700 spot 恰好就是 200，等于"把 MEL01 的默认密度推广到全队列"。
   第一代 ST 的几百节点切片因此自动拿到更粗的分割，
   用于检验 R6 的判断：legacy 富集掉到 1.07x 是组合学的（8 域切小图，
   边界边占 42%），不是生物学的。

结论怎么读
----------
- precision（域边界预测割边的精度）若在整个网格上都低，
  R6 的核心论点（域分割给候选集、最小割做选择）就不依赖域数选择；
- 富集若随边界边占比上升而向 1 回落、缩放臂若让 legacy 富集回升，
  则组合学解释得到直接支持。

用法
----
    python scripts/run_13b_ndomains_scan.py --slides MEL01 ... CSCC16
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, load_json,  # noqa: E402
                        save_json, set_seed, stamp_run)

GRID = [4, 6, 8, 10, 12]
SCALE_DENOM = 200          # 每域期望 spot 数（缩放规则的斜率）
SCALE_MIN, SCALE_MAX = 4, 16

VISIUM = ["MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"]


def main():
    ap = argparse.ArgumentParser(
        description="n_domains 敏感性扫描（不覆盖 benchmark_tools.json）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", required=True)
    ap.add_argument("--lam", type=float, default=0.3, help="BANKSY 式邻域项权重 lambda")
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy。pip install scanpy scikit-learn")

    from run_13_benchmark_tools import banksy_style_domains
    from sparta.barrier import edge_pairs

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    out = {
        "params": dict(grid=GRID, lam=args.lam, n_pcs=args.n_pcs,
                       scale_rule=f"clip(round(n_nodes/{SCALE_DENOM}), {SCALE_MIN}, {SCALE_MAX})",
                       seed=cfg["seed"]),
        "per_slide": {},
        "meta": stamp_run(cfg, {"module": "M13b-ndomains-scan"}),
    }

    ref = load_json(P.validation("benchmark_tools.json")).get("per_slide", {})

    print(f"{'slide':<8}{'nodes':>6}{'scaled_k':>9}   " +
          "".join(f"k={k:<3}" for k in GRID))
    print("-" * 78)

    for sid in args.slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
            cut = load_json(P.mincut(sid))["cut_edges"]
        except FileNotFoundError as e:
            print(f"{sid:<8} 跳过：{e}")
            continue

        n_nodes = int(A.shape[0])
        k_scaled = int(np.clip(round(n_nodes / SCALE_DENOM), SCALE_MIN, SCALE_MAX))

        # PCA 只算一次（与 run_13 完全相同的配方），存回 obsm 后各 k 复用
        if "X_pca" not in adata.obsm:
            ad = adata.copy()
            sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
            ad = ad[:, ad.var.highly_variable].copy()
            sc.pp.scale(ad, max_value=10)
            sc.tl.pca(ad, n_comps=args.n_pcs, random_state=cfg["seed"])
            adata.obsm["X_pca"] = ad.obsm["X_pca"].copy()

        pairs = edge_pairs(A, upper_only=True)
        edge_id = {(int(min(u, v)), int(max(u, v))): i
                   for i, (u, v) in enumerate(pairs)}
        cut_idx = np.array([edge_id[(min(u, v), max(u, v))] for u, v in cut
                            if (min(u, v), max(u, v)) in edge_id], int)

        rec = {"n_nodes": n_nodes, "scaled_n_domains": k_scaled, "runs": {}}
        ks = sorted(set(GRID) | {k_scaled})
        for k in ks:
            lab = banksy_style_domains(adata, A, k, args.lam, args.n_pcs, cfg["seed"])
            is_boundary = lab[pairs[:, 0]] != lab[pairs[:, 1]]
            p_boundary = float(is_boundary.mean())
            if len(cut_idx):
                p_in_cut = float(is_boundary[cut_idx].mean())
                enrich = p_in_cut / max(p_boundary, 1e-12)
                precision = float(is_boundary[cut_idx].sum() / max(is_boundary.sum(), 1))
            else:
                p_in_cut = enrich = precision = np.nan
            rec["runs"][str(k)] = dict(
                n_domains=int(len(np.unique(lab))),
                frac_boundary_edges=p_boundary,
                frac_boundary_within_cut=p_in_cut,
                enrichment=float(enrich),
                precision_of_boundary_for_cut=precision)

        # 与 run_13 主结果对数（k=8 时两者必须逐张一致）
        if sid in ref and "8" in rec["runs"]:
            d = abs(rec["runs"]["8"]["enrichment"] - ref[sid]["enrichment"])
            flag = "" if d < 1e-9 else f"  ⚠ 与 benchmark_tools 不一致（Δ={d:.2e}）"
            print(f"{sid:<8}{n_nodes:>6}{k_scaled:>9}   " +
                  "".join(f"{rec['runs'][str(k)]['enrichment']:.2f} "
                          for k in GRID if str(k) in rec["runs"]) + flag)
        else:
            print(f"{sid:<8}{n_nodes:>6}{k_scaled:>9}   （无 k=8 参照或无割边）")

        out["per_slide"][sid] = rec

    if not out["per_slide"]:
        sys.exit("没有任何切片可扫描。")

    # ---- 汇总：固定网格各 k 的中位数（分队列）+ 缩放臂 ----
    legacy = [s for s in out["per_slide"] if s not in VISIUM]
    summ = {}
    for k in GRID:
        for name, grp in [("visium", VISIUM), ("legacy_st", legacy)]:
            es = [out["per_slide"][s]["runs"][str(k)]["enrichment"]
                  for s in grp if s in out["per_slide"]
                  and np.isfinite(out["per_slide"][s]["runs"][str(k)]["enrichment"])]
            ps = [out["per_slide"][s]["runs"][str(k)]["precision_of_boundary_for_cut"]
                  for s in grp if s in out["per_slide"]
                  and np.isfinite(out["per_slide"][s]["runs"][str(k)]["precision_of_boundary_for_cut"])]
            if es:
                summ[f"fixed_k{k}_{name}"] = dict(
                    median_enrichment=float(np.median(es)),
                    median_precision=float(np.median(ps)),
                    n=len(es))
    for name, grp in [("visium", VISIUM), ("legacy_st", legacy)]:
        es, ps = [], []
        for s in grp:
            if s not in out["per_slide"]:
                continue
            ks = out["per_slide"][s]["scaled_n_domains"]
            r = out["per_slide"][s]["runs"][str(ks)]
            if np.isfinite(r["enrichment"]):
                es.append(r["enrichment"]); ps.append(r["precision_of_boundary_for_cut"])
        if es:
            summ[f"scaled_{name}"] = dict(
                median_enrichment=float(np.median(es)),
                median_precision=float(np.median(ps)), n=len(es))
    out["summary"] = summ

    all_prec = [r["precision_of_boundary_for_cut"]
                for rec in out["per_slide"].values() for r in rec["runs"].values()
                if np.isfinite(r["precision_of_boundary_for_cut"])]
    print("-" * 78)
    for k, v in summ.items():
        print(f"{k:<22} 富集中位 {v['median_enrichment']:.2f}x   "
              f"precision 中位 {v['median_precision']*100:.1f}%   n={v['n']}")
    print(f"\n全网格 precision 范围：{min(all_prec)*100:.1f}%–{max(all_prec)*100:.1f}%"
          f"（{len(all_prec)} 个 切片×k 组合）")
    print("R6 的读法：若 precision 在整个网格上都远低于 100%，"
          "'域分割给候选集、最小割做选择'就不依赖域数选择。")

    p = P.validation("ndomains_sensitivity.json")
    save_json(p, out)
    print(f"\n已写出 {p}")

    if args.no_plot:
        return
    _plot(P, out, legacy)


def _plot(P, out, legacy):
    """补充材料图：富集热图（极性、中点 1.0）+ precision 热图（量级、单色相）。

    约定：队列用分组 + 直接标注（白色间隔行），不用颜色——色相留给量。
    文字颜色随背景深浅切换（深格用白字），色标带标签。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, to_rgb

    vis = [s for s in VISIUM if s in out["per_slide"]]
    leg = [s for s in legacy if s in out["per_slide"]]
    sids = vis + [None] + leg                       # None = 组间白色间隔行
    rows = [s for s in sids if s]
    cols = GRID + ["scaled"]

    def mat(key):
        M = np.full((len(sids), len(cols)), np.nan)
        for i, s in enumerate(sids):
            if s is None:
                continue
            rec = out["per_slide"][s]
            for j, k in enumerate(GRID):
                if str(k) in rec["runs"]:
                    M[i, j] = rec["runs"][str(k)][key]
            ks = rec["scaled_n_domains"]
            if str(ks) in rec["runs"]:
                M[i, len(GRID)] = rec["runs"][str(ks)][key]
        return M

    E, Q = mat("enrichment"), mat("precision_of_boundary_for_cut") * 100
    qmax = float(np.nanmax(Q))

    lo, hi = to_rgb("#2a78d6"), to_rgb("#eb6834")
    div = LinearSegmentedColormap.from_list(
        "div", [(lo[0]*.45, lo[1]*.45, lo[2]*.45), lo, (.90, .895, .88), hi,
                (hi[0]*.45, hi[1]*.45, hi[2]*.45)])
    base = to_rgb("#3E6B8A")
    seq = LinearSegmentedColormap.from_list(
        "seq", [tuple(1 - .96 * (1 - c) for c in base), base,
                tuple(c * .42 for c in base)])

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 8.0))
    gap = sids.index(None)
    for ax, M, title, cmap, norm, fmt, clabel in [
            (axes[0], E, "Cut edges on domain boundaries\n(fold enrichment, 1.0 = chance)",
             div, TwoSlopeNorm(vcenter=1.0, vmin=0.5, vmax=2.0), "%.2f",
             "fold enrichment"),
            (axes[1], Q, "Precision of boundary for cut\n(% of boundary edges on the cut)",
             seq, plt.Normalize(vmin=0, vmax=qmax), "%.0f", "precision (%)")]:
        im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([f"{c}" if isinstance(c, int) else "scaled" for c in cols])
        ax.set_xlabel("n_domains  (BANKSY-style, $\\lambda$=0.3)")
        ax.set_yticks(range(len(sids)))
        ax.set_yticklabels([s if s else "" for s in sids], fontsize=8)
        ax.tick_params(length=0)
        for i in range(len(sids)):
            for j in range(len(cols)):
                if np.isfinite(M[i, j]):
                    dark = float(norm(M[i, j])) > 0.62
                    ax.text(j, i, fmt % M[i, j], ha="center", va="center",
                            fontsize=7, color="#fcfcfb" if dark else "#0b0b0b")
        ax.axvline(len(GRID) - .5, color="#52514e", lw=1.2)          # scaled 分隔
        ax.add_patch(plt.Rectangle((GRID.index(8) - .5, -.5), 1, len(sids),
                                   fill=False, ec="#0b0b0b", lw=1.6))  # 主文默认 k=8
        for yy, lab in [(gap / 2 - .5, "Visium"),
                        (gap + (len(sids) - gap) / 2 - .5, "first-gen ST")]:
            ax.text(-.55, yy, lab, rotation=90, va="center", ha="right",
                    fontsize=9, color="#0b0b0b", clip_on=False)
        ax.set_title(title, fontsize=10)
        cb = fig.colorbar(im, ax=ax, shrink=.55)
        cb.set_label(clabel, fontsize=9)
    fig.suptitle("S5  Min-cut vs domain boundaries across the domain count\n"
                 "(boxed column = n_domains 8 used in the main text; 'scaled' = "
                 "clip(round(n_nodes/200), 4, 16))", fontsize=11)
    fig.subplots_adjust(left=.17, right=.97, top=.88, bottom=.07, wspace=.52)
    pf = P.figure("supp_ndomains_scan.png")
    fig.savefig(pf, dpi=200); plt.close(fig)
    print(f"已保存 {pf}")


if __name__ == "__main__":
    main()
