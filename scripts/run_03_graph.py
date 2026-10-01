#!/usr/bin/env python
"""
run_03_graph.py —— M3 空间图与源汇定义
========================================

输入：data/interim/{slide_id}.scored.h5ad
输出：data/interim/{slide_id}.graph.npz
上游模块：run_02_score.py
下游模块：run_04_barrier.py

建空间邻接图并定义血管、免疫入口（源）、瘤巢核心（汇）。

两个必看的输出：
  · 图健康报告中的**连通分量数**。图碎成很多块说明 radius_um 太小，
    有效阻抗与图距离都会失去意义。
  · 源/汇/血管的节点数。任一为 0 则下游算不出来。

用法
----
    python scripts/run_03_graph.py --slide MEL01
    python scripts/run_03_graph.py --slide SCC01 --radius-um 280   # 第一代 ST
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, save_json, set_seed, stamp_run  # noqa: E402


def _require_scanpy():
    """真实数据流程需要 scanpy；合成流程（run_00_demo）不需要。"""
    try:
        import scanpy  # noqa: F401
    except ImportError:
        sys.exit(
            "需要 scanpy 才能处理真实空间数据。\n"
            "  pip install scanpy squidpy\n"
            "如果只想先验证环境与算子，请跑 python scripts/run_00_demo.py（不需要 scanpy）。"
        )


def main():
    ap = argparse.ArgumentParser(description="M3 空间图与源汇定义",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--slide", required=True)
    ap.add_argument("--radius-um", type=float, default=None,
                    help="邻接半径。不传则用 config 的 graph.radius_um")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    _require_scanpy()
    import scanpy as sc
    from sparta.graph import build_graph_radius, define_source_sink, graph_report
    from sparta.io_ import save_graph

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    radius = args.radius_um or cfg["graph"]["radius_um"]

    adata = sc.read_h5ad(P.scored(args.slide))
    coords = np.asarray(adata.obsm["spatial_um"], float)

    A, D = build_graph_radius(coords, radius_um=radius)
    rep = graph_report(A)
    print(f"[M3] 图：{rep['n_nodes']} 节点，{rep['n_edges']} 边，平均度 {rep['mean_degree']:.1f}")
    print(f"[M3] 连通分量 {rep['n_components']} 个，最大分量占 {rep['largest_component_frac']*100:.1f}%")
    if rep["n_components"] > 1:
        print(f"[M3] 警告：图不连通。{rep['n_isolated']} 个孤立节点。"
              f"考虑增大 radius_um（当前 {radius}）。跨分量的有效阻抗无意义。")

    def col(name):
        return np.asarray(adata.obs[name + "_n"].values, float) if name + "_n" in adata.obs \
               else np.full(adata.n_obs, 0.5)

    ss = define_source_sink(
        A, col("Endothelial"), col("T_NK"), col("Malignant"),
        **{k: cfg["source_sink"][k] for k in
           ("q_vessel", "q_immune_nbr", "q_malig", "q_core")},
        fallback_border_coords=coords if cfg["source_sink"]["use_border_fallback"] else None,
    )
    print(f"[M3] 血管 {ss['n_vessel']} 个，源（免疫入口）{ss['n_source']} 个，"
          f"汇（瘤巢核心）{ss['n_sink']} 个")
    if ss["used_fallback"]:
        print("[M3] 注意：内皮信号不足，已启用『组织边缘作为血管代理源』的降级路径。"
              "这一替代必须在论文方法部分明确写出。")
    for k, n_ in (("源", ss["n_source"]), ("汇", ss["n_sink"]), ("血管", ss["n_vessel"])):
        if n_ == 0:
            sys.exit(f"[M3] {k}集为空，下游无法计算。请调整 source_sink 的分位数阈值。")

    save_graph(P.graph(args.slide), A, D, ss["source"], ss["sink"], ss["vessel"],
               meta={**rep, "radius_um": radius, "used_fallback": ss["used_fallback"],
                     **stamp_run(cfg, {"module": "M3", "slide": args.slide})})
    print(f"[M3] 已写出 {P.graph(args.slide)}")


if __name__ == "__main__":
    main()
