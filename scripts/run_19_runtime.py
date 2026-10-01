#!/usr/bin/env python
"""
run_19_runtime.py —— 运行时 / 内存基准（JBHI 投稿补料）
=========================================================

输入：{sid}.scored.h5ad（图从坐标重建，与 run_03 同配方）
输出：results/validation/runtime_benchmark.json
上游模块：run_02/03/04/13
下游模块：docs/jbhi_repositioning.md 的对比总表

测什么
------
逐张切片对流水线的**每个阶段**单独计时（perf_counter，默认 3 次取最小，
squidpy 因置换检验昂贵只跑 1 次）：

    score   基因集打分（scores_from_adata）
    graph   空间图构建 + 源汇定义（build_graph_radius + define_source_sink）
    b_cell  源-汇最小割（networkx preflow-push）
    b_mab   屏蔽泊松扩散-吸收场（稀疏直接/迭代解）
    b_meta  代谢屏障（dijkstra 距离）
    banksy  BANKSY 式域分割（PCA+邻域均值+k-means，run_13 同款）
    squidpy 邻域富集（nhood_enrichment，同一张 SPARTA 图）

另用 tracemalloc 单独跑一遍 score+graph+三算子测 Python 层峰值内存
（计时与内存分开测，tracemalloc 的钩子会拖慢计时）。

公平性口径
----------
- Visium 用 default.yaml，第一代 ST 用 cscc_legacy_st.yaml（radius 150/300），
  与正式流水线完全一致，图重建半径与已有 graph.npz 的 meta 核对；
- banksy/squidpy 是**本仓库实现的对比臂**（k-means 式域分割、置换 z 值），
  不是官方 BANKSY/Squidpy 包的全部功能——表里必须这么标；
- STAGATE/GraphST 等 GPU 训练方法**不在此跑**，表中如实标注"未运行"，
  不编造数字。

用法
----
    python scripts/run_19_runtime.py [--repeats 3] [--no-squidpy]
"""
from __future__ import annotations

import argparse
import os
import platform
import sys
import time
import tracemalloc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, save_json, set_seed, stamp_run  # noqa: E402

VISIUM = ["MEL01", "MEL02", "MEL03", "MEL04",
          "CSCC01", "CSCC02", "CSCC03", "CSCC04"]
LEGACY = ["CSCC05", "CSCC06", "CSCC07", "CSCC08", "CSCC09", "CSCC10",
          "CSCC11", "CSCC12", "CSCC14", "CSCC15", "CSCC16"]
ALL_SLIDES = VISIUM + LEGACY


def config_path_for(sid: str) -> str | None:
    return None if sid in VISIUM else "configs/cscc_legacy_st.yaml"


def main():
    ap = argparse.ArgumentParser(description="运行时/内存基准",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slides", nargs="+", default=ALL_SLIDES)
    ap.add_argument("--repeats", type=int, default=3,
                    help="每阶段重复次数，取最小值（squidpy 恒为 1 次）")
    ap.add_argument("--no-squidpy", action="store_true")
    ap.add_argument("--config", default=None, help="仅当跑单一队列时手动指定")
    args = ap.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("需要 scanpy。")

    from run_13_benchmark_tools import banksy_style_domains
    from sparta.barrier import (compute_b_cell, compute_b_mab, compute_b_meta,
                                scores_from_adata)
    from sparta.graph import build_graph_radius, define_source_sink

    if not args.no_squidpy:
        try:
            import squidpy as sq
            has_sq = True
        except ImportError:
            has_sq = False
            print("[提示] 未装 squidpy，跳过其计时。")

    def timed(fn, n):
        best = float("inf")
        for _ in range(n):
            t0 = time.perf_counter()
            out = fn()
            best = min(best, time.perf_counter() - t0)
        return best, out

    env = dict(
        platform=platform.platform(), processor=platform.processor(),
        cpu_count=os.cpu_count(),
        python=platform.python_version(),
        versions={m.__name__: m.__version__ for m in
                  [np, sc, __import__("scipy"), __import__("networkx")]},
    )
    try:
        import psutil
        env["ram_gb"] = round(psutil.virtual_memory().total / 2**30, 1)
    except ImportError:
        env["ram_gb"] = None

    out = {"env": env, "params": dict(repeats=args.repeats, squidpy=not args.no_squidpy),
           "per_slide": {},
           "meta": stamp_run(load_config(args.config or config_path_for(args.slides[0])),
                             {"module": "M19-runtime"})}

    print(f"环境：{env['processor']} × {env['cpu_count']} 核，"
          f"RAM {env['ram_gb']} GB，numpy {env['versions']['numpy']}，"
          f"networkx {env['versions']['networkx']}")
    hdr = (f"{'slide':<8}{'nodes':>6}{'edges':>7}{'score':>7}{'graph':>7}"
           f"{'b_cell':>8}{'b_mab':>8}{'b_meta':>8}{'SPARTA':>8}{'banksy':>8}{'squidpy':>9}{'peakMB':>8}")
    print(hdr)
    print("-" * len(hdr))

    for sid in args.slides:
        cfg = load_config(args.config or config_path_for(sid))
        set_seed(cfg["seed"])
        P = Paths(cfg)
        adata = sc.read_h5ad(P.scored(sid))
        coords = np.asarray(adata.obsm["spatial_um"], float)

        # 与正式流水线同配方重建图，并与已存图核对半径
        radius = cfg["graph"]["radius_um"]
        try:
            _, _, _, _, _, gmeta = load_graph(P.graph(sid))
            if abs(float(gmeta.get("radius_um", radius)) - radius) > 1e-6:
                print(f"[警告] {sid} 配置 radius={radius} 与已存图 "
                      f"{gmeta.get('radius_um')} 不一致，计时以配置为准")
        except FileNotFoundError:
            pass

        def col(name):
            return (np.asarray(adata.obs[name + "_n"].values, float)
                    if name + "_n" in adata.obs
                    else np.full(adata.n_obs, 0.5))

        t_score, S = timed(lambda: scores_from_adata(adata), args.repeats)

        def build():
            A, D = build_graph_radius(coords, radius_um=radius)
            ss = define_source_sink(
                A, col("Endothelial"), col("T_NK"), col("Malignant"),
                **{k: cfg["source_sink"][k] for k in
                   ("q_vessel", "q_immune_nbr", "q_malig", "q_core")},
                fallback_border_coords=coords
                if cfg["source_sink"]["use_border_fallback"] else None)
            return A, D, ss

        t_graph, (A, D, ss) = timed(build, args.repeats)
        source, sink, vessel = ss["source"], ss["sink"], ss["vessel"]
        bcfg = cfg["barrier"]
        t_cell, _ = timed(lambda: compute_b_cell(
            A, S["ecm"], S["caf"], source, sink, **bcfg["b_cell"]), args.repeats)
        t_mab, _ = timed(lambda: compute_b_mab(
            A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **bcfg["b_mab"]),
            args.repeats)
        t_meta, _ = timed(lambda: compute_b_meta(
            A, S["hypoxia"], S["proliferation"], S["efflux"], vessel, D=D,
            **bcfg["b_meta"]), args.repeats)

        t_banksy, lab = timed(
            lambda: banksy_style_domains(adata, A, 8, 0.3, 20, cfg["seed"]),
            args.repeats)

        t_sq = None
        sq_err = None
        if not args.no_squidpy and has_sq:
            try:
                ad = adata.copy()
                ad.obs["_sparta_domain"] = pd_cat([str(x) for x in lab])
                ad.obsp["spatial_connectivities"] = (A != 0).astype(float).tocsr()
                ad.obsp["spatial_distances"] = D.tocsr()
                ad.uns["spatial_neighbors"] = dict(
                    connectivities_key="spatial_connectivities",
                    distances_key="spatial_distances",
                    params=dict(n_neighbors=-1, coord_type="generic",
                                radius=None, transform=None,
                                source="sparta.graph (radius adjacency)"))
                t0 = time.perf_counter()
                sq.gr.nhood_enrichment(ad, cluster_key="_sparta_domain",
                                       seed=cfg["seed"])
                t_sq = time.perf_counter() - t0
            except Exception as e:  # noqa: BLE001
                sq_err = f"{type(e).__name__}: {e}"

        # 内存峰值（单独一遍，不带计时）
        tracemalloc.start()
        scores_from_adata(adata)
        A2, D2 = build_graph_radius(coords, radius_um=radius)
        compute_b_cell(A2, S["ecm"], S["caf"], source, sink, **bcfg["b_cell"])
        compute_b_mab(A2, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                      **bcfg["b_mab"])
        compute_b_meta(A2, S["hypoxia"], S["proliferation"], S["efflux"],
                       vessel, D=D2, **bcfg["b_meta"])
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        rec = dict(
            platform=("visium" if sid in VISIUM else "legacy_st"),
            n_nodes=int(A.shape[0]), n_edges=int((A != 0).sum() // 2),
            radius_um=float(radius),
            t_score_s=t_score, t_graph_s=t_graph, t_b_cell_s=t_cell,
            t_b_mab_s=t_mab, t_b_meta_s=t_meta,
            t_sparta_core_s=t_score + t_graph + t_cell + t_mab + t_meta,
            t_banksy_style_s=t_banksy, t_squidpy_s=t_sq,
            peak_mem_mb=peak / 2**20,
        )
        if sq_err:
            rec["squidpy_error"] = sq_err
        out["per_slide"][sid] = rec
        print(f"{sid:<8}{rec['n_nodes']:>6}{rec['n_edges']:>7}"
              f"{t_score:>7.2f}{t_graph:>7.2f}{t_cell:>8.3f}{t_mab:>8.3f}"
              f"{t_meta:>8.2f}{rec['t_sparta_core_s']:>8.2f}{t_banksy:>8.2f}"
              f"{(t_sq if t_sq is not None else float('nan')):>9.2f}"
              f"{rec['peak_mem_mb']:>8.0f}")

    # ---- 汇总 ----
    recs = list(out["per_slide"].values())
    summ = {}
    for name, grp in [("visium", VISIUM), ("legacy_st", LEGACY),
                      ("all", list(out["per_slide"]))]:
        rs = [out["per_slide"][s] for s in grp if s in out["per_slide"]]
        if not rs:
            continue
        sq = [r["t_squidpy_s"] for r in rs if r["t_squidpy_s"] is not None]
        summ[name] = dict(
            n=len(rs),
            nodes=f"{min(r['n_nodes'] for r in rs)}–{max(r['n_nodes'] for r in rs)}",
            t_sparta_core_s_median=float(np.median([r["t_sparta_core_s"] for r in rs])),
            t_sparta_core_s_total=float(sum(r["t_sparta_core_s"] for r in rs)),
            t_banksy_style_s_median=float(np.median([r["t_banksy_style_s"] for r in rs])),
            t_squidpy_s_median=(float(np.median(sq)) if sq else None),
            peak_mem_mb_max=float(max(r["peak_mem_mb"] for r in rs)),
        )
    out["summary"] = summ

    print("-" * len(hdr))
    for name, s in summ.items():
        sq = (f"{s['t_squidpy_s_median']:.1f}s" if s["t_squidpy_s_median"] else "—")
        print(f"{name:<12} n={s['n']}  SPARTA 中位 {s['t_sparta_core_s_median']:.2f}s"
              f"（19 张合计口径 {s['t_sparta_core_s_total']:.1f}s）"
              f"  banksy 中位 {s['t_banksy_style_s_median']:.2f}s  squidpy 中位 {sq}"
              f"  峰值内存 ≤{s['peak_mem_mb_max']:.0f} MB")

    p = P.validation("runtime_benchmark.json")
    save_json(p, out)
    print(f"\n已写出 {p}")


def pd_cat(x):
    import pandas as pd
    return pd.Categorical(x)


if __name__ == "__main__":
    main()
