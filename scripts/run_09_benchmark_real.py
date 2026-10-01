# M9 真实数据 benchmark：最小割 vs 局部统计量对空间排布的敏感性
"""
在 4 张真实黑色素瘤 Visium 切片上，对比四类指标对【空间重排】的敏感性。

核心逻辑
--------
S1 空间重排 = 保持每个 spot 的分子状态（组成不变）、只打乱空间位置。
于是：
  · 若某指标在重排后显著下降 -> 它携带【排布/拓扑】信息（最小割应当如此）
  · 若某指标在重排后基本不变 -> 它只携带【组成】信息
  · 若某指标重排后仍上升或不变但无法给出屏障语义 -> 局部统计的局限

这里对 4 张真实切片分别计算：
  B_cell (最小割, 1/maxflow)
  CAF 密度 (top-CAF 比例, 纯组成)
  邻域 CAF 富集 (top-CAF spot 的 kNN 中 top-CAF 比例, 局部统计)
  Ripley's L 峰值 (对 top-CAF spot, 二阶统计)
并报告它们对空间重排（50 次）的 z 得分，同时给出 B_cell 与各局部量的
空间相关（验证 B_cell 捕获的是不同信息）。

运行：
  python scripts/run_09_benchmark_real.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import compute_b_cell, compute_b_cell_field, scores_from_adata

cfg = load_config(None)
P = Paths(cfg)
# 只跑黑色素瘤是**有意**的子集（这份对比只在黑色素瘤上做），不是漏改。
# 命令行可覆盖：python scripts/run_09_benchmark_real.py CSCC01 CSCC05
SLIDES = sys.argv[1:] or ["MEL01", "MEL02", "MEL03", "MEL04"]
N_PERM = 50
SEED = 0
cfg_cell = cfg["barrier"]["b_cell"]
OUT = Path(__file__).resolve().parents[1] / "results" / "validation"


def ripley_L_max(coords, idx, r_max=400.0):
    pts = coords[idx]
    n = len(pts)
    if n < 2:
        return 0.0
    xmin, ymin = coords.min(axis=0); xmax, ymax = coords.max(axis=0)
    area = max((xmax - xmin) * (ymax - ymin), 1.0)
    tree = cKDTree(pts)
    r_grid = np.linspace(50, r_max, 20)
    K = np.array([area * tree.count_neighbors(tree, r=r) / (n * (n - 1)) for r in r_grid])
    L = np.sqrt(np.maximum(K / np.pi, 0)) - r_grid
    return float(L.max())


def neighbor_enrichment(coords, is_high, k=10):
    hi = np.where(is_high)[0]
    if len(hi) == 0:
        return 0.0
    tree = cKDTree(coords)
    fracs = []
    for i in hi:
        d, nn = tree.query(coords[i], k=min(k + 1, len(coords)))
        nn = np.atleast_1d(nn)
        fracs.append(is_high[nn].mean())
    return float(np.mean(fracs))


print(f"{'切片':<7}{'指标':<24}{'真实':>10}{'重排均值':>10}{'重排sd':>9}{'z':>8}")
print("-" * 70)

summary = {}
for sid in SLIDES:
    adata = sc.read_h5ad(P.scored(sid))
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    S = scores_from_adata(adata)
    # QC 后坐标已标定为微米，存于 obsm['spatial_um']
    coords = np.asarray(adata.obsm["spatial_um"]).astype(float)

    caf = S["caf"]
    is_high = caf >= np.quantile(caf, 0.85)   # top 15% 视为高 CAF
    rng = np.random.default_rng(SEED)

    # ---- 真实值 ----
    b_real = compute_b_cell(A, S["ecm"], S["caf"], source, sink, **cfg_cell)["b_cell"]
    dens_real = is_high.mean()
    nb_real = neighbor_enrichment(coords, is_high)
    L_real = ripley_L_max(coords, np.where(is_high)[0])
    bfield = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg_cell)["b_cell_field"]

    # ---- 重排零分布（S1 fixed：打乱 spot 在网格上的位置）----
    b_perm, dens_perm, nb_perm, L_perm = [], [], [], []
    for _ in range(N_PERM):
        perm = rng.permutation(len(coords))
        ecm_p, caf_p = S["ecm"][perm], caf[perm]
        b_perm.append(compute_b_cell(A, ecm_p, caf_p, source, sink, **cfg_cell)["b_cell"])
        # 局部统计量基于"位置不变、标签重排"
        is_high_p = caf_p >= np.quantile(caf_p, 0.85)
        dens_perm.append(is_high_p.mean())
        nb_perm.append(neighbor_enrichment(coords, is_high_p))
        L_perm.append(ripley_L_max(coords, np.where(is_high_p)[0]))

    def zstat(real, perms):
        arr = np.array(perms)
        sd = arr.std()
        if sd == 0:
            return float("nan")
        return float((real - arr.mean()) / sd)

    rows = [
        ("B_cell(最小割)", b_real, b_perm),
        ("CAF密度(组成)", dens_real, dens_perm),
        ("邻域CAF富集", nb_real, nb_perm),
        ("Ripley's L峰值", L_real, L_perm),
    ]
    for name, real, perms in rows:
        z = zstat(real, perms)
        print(f"{sid:<7}{name:<24}{real:>10.4f}{np.mean(perms):>10.4f}"
              f"{np.std(perms):>9.4f}{z:>8.2f}")

    # B_cell 与局部量的空间相关
    rho_dens = spearmanr(bfield, caf).statistic
    summary[sid] = dict(
        b_cell=b_real, z_b_cell=zstat(b_real, b_perm),
        dens=dens_real, z_dens=zstat(dens_real, dens_perm),
        nbr=nb_real, z_nbr=zstat(nb_real, nb_perm),
        ripleyL=L_real, z_ripley=zstat(L_real, L_perm),
        rho_Bcell_vs_CAF=float(rho_dens))
    print("-" * 70)

print("\n关键判定：")
for sid, r in summary.items():
    topo = r["z_b_cell"]
    max_local_z = max(r["z_dens"], r["z_nbr"], r["z_ripley"])
    print(f"  {sid}: B_cell 对排布 z={topo:+.2f}  局部量最大 z={max_local_z:+.2f}  "
          f"B_cell~CAF 相关 ρ={r['rho_Bcell_vs_CAF']:+.2f}")

import json
with open(OUT / "benchmark_real.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)
print(f"\n已写出 {OUT / 'benchmark_real.json'}")
