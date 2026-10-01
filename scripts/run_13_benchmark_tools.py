#!/usr/bin/env python
"""
run_13_benchmark_tools.py —— 与现有空间分析工具的正面对比
============================================================

输入：{sid}.scored.h5ad + {sid}.graph.npz + {sid}.mincut.json
输出：results/validation/benchmark_tools.json + results/figures/benchmark_tools.png
上游模块：run_03_graph.py / run_04_barrier.py
下游模块：无（成文用，对应 Results 的"与现有工具的关系"一节）

它回答审稿人一定会问的那句话
----------------------------
    "你的最小割相对 BANKSY 的空间域、Squidpy 的邻域富集，多了什么？"

回答分两半，都必须诚实地报出来：

**一半是"我们并不比它们更早看见这条带"**——如果最小割画出的封锁线
落在标准空间域分割给出的域边界上（富集倍数远大于 1），
那说明这条结构本身用现成工具也能看到。这一点要主动承认，
硬说"只有我们能发现"会被一眼识破。

**另一半才是本框架真正独有的**：域边界只告诉你"这里两侧不一样"，
它没有方向、没有容量、也无法被反事实操作。最小割额外给出：
  · **容量**——这条边让多少通量过得去（可与药物尺寸挂钩）；
  · **方向**——它是相对某一对源/汇（血管->瘤巢）定义的，换一对源汇就是另一条线；
  · **可反事实**——可以在割上开缺口重算（S2），域边界做不了这件事；
  · **可跨模态**——同一张图换边权语义即得到 B_mAb，域分割没有这个自由度。

方法说明（必须写进 Methods）
----------------------------
本脚本用的是 **BANKSY 式**（BANKSY-style）空间域分割，不是 BANKSY 官方实现：
把每个 spot 的主成分与其邻域均值主成分拼接成增广特征再聚类
（BANKSY 论文的核心思想），只依赖 scanpy + scikit-learn，不引入新依赖。
若你安装了官方 banksy/squidpy，可用 --domains-from 传入一列现成的域标签，
本脚本会直接用它来做同样的对比。Squidpy 的 nhood_enrichment 在装了 squidpy
时会一并计算并记录。

用法
----
    python scripts/run_13_benchmark_tools.py --slides MEL01 MEL02 MEL03 MEL04 \
                                             CSCC01 CSCC02 CSCC03 CSCC04
    python scripts/run_13_benchmark_tools.py --slides MEL01 --n-domains 10
    python scripts/run_13_benchmark_tools.py --slides MEL01 --domains-from banksy_labels
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, load_json, save_json, set_seed, stamp_run  # noqa: E402


def banksy_style_domains(adata, A, n_domains: int, lam: float, n_pcs: int, seed: int):
    """BANKSY 式空间域：把自身主成分与邻域均值主成分拼接后聚类。

    lam 是邻域项的权重（BANKSY 论文里的 lambda）：0 等于纯表达聚类，
    1 等于只看邻域。0.2–0.4 是常用区间，域会明显更连片。
    """
    import scanpy as sc
    from sklearn.cluster import KMeans

    ad = adata.copy()
    if "X_pca" not in ad.obsm:
        sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
        ad = ad[:, ad.var.highly_variable].copy()
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=n_pcs, random_state=seed)
    X = np.asarray(ad.obsm["X_pca"][:, :n_pcs], float)

    deg = np.maximum(np.asarray(A.sum(axis=1)).ravel(), 1.0)
    Xn = (A @ X) / deg[:, None]                       # 邻域均值
    Z = np.hstack([np.sqrt(1 - lam) * X, np.sqrt(lam) * Xn])
    lab = KMeans(n_clusters=n_domains, n_init=10, random_state=seed).fit_predict(Z)
    return lab


def main():
    ap = argparse.ArgumentParser(description="与现有空间工具的对比",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", required=True)
    ap.add_argument("--n-domains", type=int, default=8, help="空间域数量")
    ap.add_argument("--lam", type=float, default=0.3, help="BANKSY 式邻域项权重 lambda")
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--domains-from", default=None,
                    help="改用 adata.obs 里现成的一列域标签（例如官方 BANKSY 的输出）")
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy。pip install scanpy scikit-learn")

    from sparta.barrier import edge_pairs

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)

    try:
        import squidpy as sq
        has_sq = True
    except ImportError:
        has_sq = False
        print("[提示] 未安装 squidpy，跳过 nhood_enrichment；"
              "BANKSY 式域对比照常进行。")
        print("        补上它：pip install squidpy   然后重跑本脚本。")
        print("        产物里会记 params.squidpy=false，所以'没跑'这件事是留痕的，"
              "不会被误当成'跑了但没结果'。")
    # API 对照 squidpy 1.8.3 文档核实过（2026-08-27）：
    #   sq.gr.nhood_enrichment(adata, cluster_key=..., seed=...) 接受这两个参数，
    #   结果写在 adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore" / "count"]。
    # 若将来 squidpy 改了存放位置，下面那段 try/except 会把异常记进
    # rec["squidpy_error"] 而不是中断主流程——主对比不依赖 squidpy。

    out = {"params": dict(n_domains=args.n_domains, lam=args.lam, n_pcs=args.n_pcs,
                          domains_from=args.domains_from, squidpy=has_sq),
           "per_slide": {}, "meta": stamp_run(cfg, {"module": "M13-benchmark-tools"})}

    print(f"{'slide':<8}{'边数':>7}{'割边':>7}{'边界边%':>9}{'割边中边界%':>12}"
          f"{'富集':>7}{'precision':>10}{'recall':>8}")
    print("-" * 78)

    for sid in args.slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
            cut = load_json(P.mincut(sid))["cut_edges"]
        except FileNotFoundError as e:
            print(f"{sid:<8} 跳过：{e}")
            continue

        if args.domains_from and args.domains_from in adata.obs:
            lab = np.asarray(adata.obs[args.domains_from].values)
            lab = np.unique(lab, return_inverse=True)[1]
        else:
            lab = banksy_style_domains(adata, A, args.n_domains, args.lam,
                                       args.n_pcs, cfg["seed"])

        pairs = edge_pairs(A, upper_only=True)
        is_boundary = lab[pairs[:, 0]] != lab[pairs[:, 1]]
        edge_id = {(int(min(u, v)), int(max(u, v))): i for i, (u, v) in enumerate(pairs)}
        cut_idx = np.array([edge_id[(min(u, v), max(u, v))] for u, v in cut
                            if (min(u, v), max(u, v)) in edge_id], int)

        p_boundary = float(is_boundary.mean())
        if len(cut_idx):
            p_boundary_in_cut = float(is_boundary[cut_idx].mean())
            enrich = p_boundary_in_cut / max(p_boundary, 1e-12)
            # 把"域边界"当成对最小割的预测器：能召回多少割边、精度如何
            recall = float(is_boundary[cut_idx].sum() / max(len(cut_idx), 1))
            precision = float(is_boundary[cut_idx].sum() / max(is_boundary.sum(), 1))
        else:
            p_boundary_in_cut = enrich = recall = precision = np.nan

        rec = dict(n_edges=int(len(pairs)), n_cut_edges=int(len(cut_idx)),
                   n_domains=int(len(np.unique(lab))),
                   frac_boundary_edges=p_boundary,
                   frac_boundary_within_cut=p_boundary_in_cut,
                   enrichment=float(enrich), recall_of_cut_by_boundary=recall,
                   precision_of_boundary_for_cut=precision)

        if has_sq:
            try:
                adata.obs["_sparta_domain"] = [str(x) for x in lab]
                adata.obs["_sparta_domain"] = adata.obs["_sparta_domain"].astype("category")

                # 直接把 SPARTA 已经建好的图交给 squidpy，**不要让它按 radius 重建**。
                #
                # 2026-08-28 修的 bug：原来调
                #     sq.gr.spatial_neighbors(..., radius=cfg["graph"]["radius_um"])
                # 而跑第一代 ST 时若没带 --config configs/cscc_legacy_st.yaml，
                # cfg 就是 default.yaml，radius_um=150——可第一代 ST 的点间距是 200 μm，
                # 150 μm 半径**一个邻居都连不上**。空图下 nhood_enrichment 的置换零分布
                # 方差为 0，z 全是 NaN，而且 squidpy 不抛异常，
                # 于是 11 张切片安静地留下一串 NaN，看上去像"squidpy 在这个平台上不支持"。
                #
                # 用同一张图还有个好处：这本来就是更正确的对比口径——
                # 域边界与最小割是在**同一个邻接**上比的，不是两张不同的图。
                adata.obsp["spatial_connectivities"] = (A != 0).astype(float).tocsr()
                adata.obsp["spatial_distances"] = D.tocsr()
                adata.uns["spatial_neighbors"] = dict(
                    connectivities_key="spatial_connectivities",
                    distances_key="spatial_distances",
                    params=dict(n_neighbors=-1, coord_type="generic",
                                radius=None, transform=None,
                                source="sparta.graph (radius adjacency)"))
                n_edges_sq = int(adata.obsp["spatial_connectivities"].nnz // 2)
                if n_edges_sq == 0:
                    raise ValueError("交给 squidpy 的邻接图没有任何边")

                sq.gr.nhood_enrichment(adata, cluster_key="_sparta_domain", seed=cfg["seed"])
                z = np.asarray(adata.uns["_sparta_domain_nhood_enrichment"]["zscore"], float)
                off = z[~np.eye(len(z), dtype=bool)]
                if not np.isfinite(off).any():
                    # 全 NaN 必须留痕，不能当成"跑过了但没结果"
                    raise ValueError(
                        f"nhood_enrichment 的 z 全是 NaN（{n_edges_sq} 条边，"
                        f"{len(np.unique(lab))} 个域）——置换零分布方差为 0")
                rec["squidpy_n_edges"] = n_edges_sq
                rec["squidpy_nhood_z_offdiag_min"] = float(np.nanmin(off))
                rec["squidpy_nhood_z_offdiag_max"] = float(np.nanmax(off))
            except Exception as e:  # noqa: BLE001 —— squidpy 版本差异不应中断主流程
                rec["squidpy_error"] = f"{type(e).__name__}: {e}"
                print(f"         [squidpy] {sid} 失败：{type(e).__name__}: {e}")

        out["per_slide"][sid] = rec
        print(f"{sid:<8}{len(pairs):>7}{len(cut_idx):>7}{p_boundary*100:>8.1f}%"
              f"{p_boundary_in_cut*100:>11.1f}%{enrich:>7.2f}"
              f"{precision*100:>9.1f}%{recall*100:>7.1f}%")

    # 读回再合并（2026-08-27 修）：旧实现整文件覆盖，先跑 Visium 再跑 legacy 会把
    # 前一批的结果冲掉，写论文汇总时以为另一队列没跑过。分队列跑时各自合并进去。
    old = load_json(P.validation("benchmark_tools.json")) if P.validation("benchmark_tools.json").exists() else {}
    old_per = old.get("per_slide", {}) if isinstance(old, dict) else {}
    old_per.update(out["per_slide"])
    out["per_slide"] = old_per
    out["slides_in_this_run"] = list(args.slides)

    if not out["per_slide"]:
        sys.exit("没有任何切片可对比。")

    es = [r["enrichment"] for r in out["per_slide"].values() if np.isfinite(r["enrichment"])]
    ps = [r["precision_of_boundary_for_cut"] for r in out["per_slide"].values()
          if np.isfinite(r["precision_of_boundary_for_cut"])]
    print("-" * 78)
    print(f"割边落在域边界上的富集：中位 {np.median(es):.2f}x（范围 {min(es):.2f}–{max(es):.2f}x）")
    print(f"反过来，域边界对割边的精度：中位 {np.median(ps)*100:.1f}%")
    print()
    print("怎么写进论文（两句话都要写）：")
    print(f"  1) 最小割封锁线与标准空间域边界高度重合（富集 {np.median(es):.1f}x），"
          f"说明这条结构本身用现成工具也能看到——本框架不宣称它是新发现的结构；")
    print(f"  2) 但域边界只有'两侧不同'这一个属性：中位仅 {np.median(ps)*100:.0f}% 的域边界"
          f"落在割上，且它不带容量、不带源汇方向、无法做反事实操作，")
    print( "     也无法在同一张图上换边权语义得到抗体屏障。这四点是本框架的增量。")

    p = P.validation("benchmark_tools.json")
    save_json(p, out)
    print(f"\n已写出 {p}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sids = list(out["per_slide"])
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].bar(sids, [out["per_slide"][s]["enrichment"] for s in sids], color="#3E6B8A")
        ax[0].axhline(1.0, color="#C0392B", ls="--", lw=1.2, label="no enrichment")
        ax[0].set_ylabel("Min-cut edges on domain boundaries\n(fold enrichment)")
        ax[0].set_title("Does a standard spatial-domain method\nsee the same line?")
        ax[0].legend(fontsize=8); ax[0].tick_params(axis="x", rotation=45)

        w = 0.38
        xs = np.arange(len(sids))
        ax[1].bar(xs - w / 2, [out["per_slide"][s]["recall_of_cut_by_boundary"] * 100 for s in sids],
                  w, label="recall of cut by boundary", color="#6FA8A0")
        ax[1].bar(xs + w / 2, [out["per_slide"][s]["precision_of_boundary_for_cut"] * 100 for s in sids],
                  w, label="precision of boundary for cut", color="#D08C60")
        ax[1].set_xticks(xs); ax[1].set_xticklabels(sids, rotation=45)
        ax[1].set_ylabel("%")
        ax[1].set_title("Domain boundary as a predictor of the min-cut")
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        pf = P.figure("benchmark_tools.png")
        fig.savefig(pf, dpi=150); plt.close(fig)
        print(f"已保存 {pf}")


if __name__ == "__main__":
    main()
