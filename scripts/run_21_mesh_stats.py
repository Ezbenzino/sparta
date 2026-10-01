#!/usr/bin/env python
"""
run_21_mesh_stats.py — 网孔尺寸排除统计（补全 60–65% 声明的可复现性缺口）
=========================================================================
从 interim/{slide}.scored.h5ad + {slide}.graph.npz 出发，重算逐边的有效网孔
尺寸 ξ 和排除标志（ξ < r_nm），汇总到 results/validation/mesh_stats.json。

不修改 barrier.py 的返回签名，不重跑管线，只读已有产物。
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import anndata as ad
import scipy.sparse as sp

from sparta.io_ import Paths, admitted_slides, load_config, save_json, stamp_run


def _pair_mean(arr, pairs):
    """逐边取两端点的均值。"""
    return 0.5 * (arr[pairs[:, 0]] + arr[pairs[:, 1]])


def _edge_pairs(A: sp.spmatrix, upper_only: bool = True):
    """从稀疏邻接矩阵取边对。"""
    A = sp.triu(A) if upper_only else A
    coo = A.tocoo()
    return np.column_stack([coo.row, coo.col])


def main():
    cfg = load_config()
    P = Paths(cfg)
    slides = admitted_slides(P)
    if not slides:
        sys.exit("台账里没有 status=ingested 的切片。")

    bmab_cfg = cfg["barrier"]["b_mab"]
    r_nm = bmab_cfg["r_nm"]
    xi0 = bmab_cfg["xi0_nm"]
    beta = bmab_cfg["beta"]
    print(f"[mesh] r_nm={r_nm}, xi0={xi0}, beta={beta}")
    print(f"[mesh] slides: {len(slides)}")

    per_slide = {}
    for sid in slides:
        scored_path = P.scored(sid)
        graph_path = P.graph(sid)
        if not scored_path.exists() or not graph_path.exists():
            print(f"  {sid}: SKIP (scored/graph not found)")
            continue

        adata = ad.read_h5ad(scored_path)
        z = np.load(graph_path, allow_pickle=True)
        shape = tuple(z["A_shape"])
        A = sp.coo_matrix((z["A_data"], (z["A_row"], z["A_col"])),
                          shape=shape).tocsr()

        crosslink = adata.obs["ECM_crosslink_n"].values.astype(float)
        pairs = _edge_pairs(A, upper_only=True)
        if len(pairs) == 0:
            print(f"  {sid}: SKIP (no edges)")
            continue

        xl_e = _pair_mean(crosslink, pairs)
        xi = xi0 * np.exp(-beta * xl_e)
        s = r_nm / np.maximum(xi, 1e-6)
        excluded = s >= 1.0
        frac_excluded = float(np.mean(excluded))
        xi_median = float(np.median(xi))
        xi_min = float(np.min(xi))
        xi_max = float(np.max(xi))

        per_slide[sid] = {
            "n_edges": int(len(pairs)),
            "r_nm": r_nm,
            "xi0_nm": xi0,
            "beta": beta,
            "xi_median_nm": round(xi_median, 3),
            "xi_min_nm": round(xi_min, 3),
            "xi_max_nm": round(xi_max, 3),
            "frac_size_excluded": round(frac_excluded, 4),
            "frac_size_excluded_pct": round(frac_excluded * 100, 1),
        }
        print(f"  {sid}: xi_med={xi_median:.2f} nm, "
              f"excluded={frac_excluded*100:.1f}%")

    fracs = [v["frac_size_excluded"] for v in per_slide.values()]
    summary = {
        "n_slides": len(per_slide),
        "r_nm": r_nm,
        "xi0_nm": xi0,
        "beta": beta,
        "frac_excluded_median": round(float(np.median(fracs)), 4),
        "frac_excluded_min": round(float(np.min(fracs)), 4),
        "frac_excluded_max": round(float(np.max(fracs)), 4),
        "frac_excluded_median_pct": round(float(np.median(fracs)) * 100, 1),
        "frac_excluded_range_pct": [
            round(float(np.min(fracs)) * 100, 1),
            round(float(np.max(fracs)) * 100, 1),
        ],
        "per_slide": per_slide,
        "meta": stamp_run(cfg, {"module": "M21-mesh-stats"}),
    }

    out = P.validation("mesh_stats.json")
    save_json(out, summary)
    print(f"\nDone: {out}")
    print(f"  median excluded: {summary['frac_excluded_median_pct']}%")
    print(f"  range: {summary['frac_excluded_range_pct']}")


if __name__ == "__main__":
    main()
