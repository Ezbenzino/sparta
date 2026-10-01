# M10 扩展 benchmark：最小割 vs 局部统计量 vs 免疫排斥简单代理
"""
双队列 8 张切片，六类指标对【空间重排】的敏感性 + B_cell 空间场相关性。

新增两类免疫排斥代理：
  T-to-core distance : T_NK 高 spot 到肿瘤核心(汇)的最近欧氏距离均值 ——
                       经典"免疫排斥"的代理，越大表示 T 细胞被挡在外面
  T infiltration     : T_NK 高 spot 占比 —— 经典"免疫浸润"的代理

若 B_cell(最小割) 携带与这些代理不同的拓扑信息，且对重排敏感，
则证明最小割提供了 T 距离/浸润比例看不到的"屏障结构"信号。

运行：python scripts/run_10_benchmark_ext.py
"""
from __future__ import annotations
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sparta.io_ import Paths, admitted_slides, load_config, load_graph
from sparta.barrier import compute_b_cell, compute_b_cell_field, scores_from_adata

cfg = load_config(None)
P = Paths(cfg)

# 名单以台账为准（status == ingested），不再硬编码。
# 2026-08-27 之前这里写死了 8 张，队列扩到 19 张后它会**安静地**只算其中 8 张，
# 而 R5 引用的正是这个文件的数字。命令行可覆盖：
#     python scripts/run_10_benchmark_ext.py MEL01 CSCC05
SLIDES = sys.argv[1:] or admitted_slides(P)
N_PERM = 50
SEED = 0
cfg_cell = cfg["barrier"]["b_cell"]
OUT = Path(__file__).resolve().parents[1] / "results" / "validation"


def ripley_L_max(coords, idx, r_max=400.0):
    pts = coords[idx]; n = len(pts)
    if n < 2: return 0.0
    xmin, ymin = coords.min(axis=0); xmax, ymax = coords.max(axis=0)
    area = max((xmax - xmin) * (ymax - ymin), 1.0)
    tree = cKDTree(pts)
    r_grid = np.linspace(50, r_max, 20)
    K = np.array([area * tree.count_neighbors(tree, r=r) / (n * (n - 1)) for r in r_grid])
    L = np.sqrt(np.maximum(K / np.pi, 0)) - r_grid
    return float(L.max())


def neighbor_enrichment(coords, is_high, k=10):
    hi = np.where(is_high)[0]
    if len(hi) == 0: return 0.0
    tree = cKDTree(coords)
    fracs = [is_high[np.atleast_1d(tree.query(coords[i], k=min(k + 1, len(coords)))[1])].mean()
             for i in hi]
    return float(np.mean(fracs))


def t_to_core_dist(coords, t_high, sink):
    """T_NK-high spot 到最近肿瘤核心 spot 的距离均值（免疫排斥代理）。"""
    if len(t_high) == 0: return np.nan
    tree = cKDTree(coords[sink])
    d, _ = tree.query(coords[t_high])
    return float(np.mean(d))


print(f"{'切片':<7}{'指标':<26}{'真实':>9}{'重排均值':>9}{'z':>7}")
print("-" * 62)
summary = {}
for sid in SLIDES:
    adata = sc.read_h5ad(P.scored(sid))
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    S = scores_from_adata(adata)
    coords = np.asarray(adata.obsm["spatial_um"]).astype(float)
    caf = S["caf"]; tnk = np.asarray(adata.obs["T_NK_n"].values, dtype=float)
    is_high_caf = caf >= np.quantile(caf, 0.85)
    is_high_t = tnk >= np.quantile(tnk, 0.85)
    rng = np.random.default_rng(SEED)

    b_real = compute_b_cell(A, S["ecm"], S["caf"], source, sink, **cfg_cell)["b_cell"]
    dens_real = is_high_caf.mean()
    nb_real = neighbor_enrichment(coords, is_high_caf)
    L_real = ripley_L_max(coords, np.where(is_high_caf)[0])
    tdist_real = t_to_core_dist(coords, np.where(is_high_t)[0], sink)
    tinf_real = is_high_t.mean()
    bfield = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg_cell)["b_cell_field"]

    b_perm = []; dens_perm = []; nb_perm = []; L_perm = []; td_perm = []; ti_perm = []
    for _ in range(N_PERM):
        perm = rng.permutation(len(coords))
        caf_p, tnk_p = caf[perm], tnk[perm]
        b_perm.append(compute_b_cell(A, S["ecm"][perm], caf_p, source, sink, **cfg_cell)["b_cell"])
        ih_caf = caf_p >= np.quantile(caf_p, 0.85)
        ih_t = tnk_p >= np.quantile(tnk_p, 0.85)
        dens_perm.append(ih_caf.mean())
        nb_perm.append(neighbor_enrichment(coords, ih_caf))
        L_perm.append(ripley_L_max(coords, np.where(ih_caf)[0]))
        td_perm.append(t_to_core_dist(coords, np.where(ih_t)[0], sink))
        ti_perm.append(ih_t.mean())

    def zstat(real, perms):
        """置换 z 值。零方差或**近零方差**的零分布一律返回 NaN。

        为什么不能只判 sd == 0（2026-08-27 踩到）
        ----------------------------------------
        第一代 ST 是规则的 200 μm 交错阵列，Ripley's L 在这种格点上几乎不随位置置换
        变化：零分布的 sd 能小到 1e-14 量级，于是 z = (real-mean)/sd 爆到 1e15，
        而这个数会一路进到结果表和正文里，看上去还像个"极显著"。
        实际发生过：11 张第一代 ST 里 7 张的 z_ripley 是 1e14–1e15，另 4 张是 NaN。

        组成类指标（CAF 密度、T 浸润比例）也有同一个毛病的温和版：它们在置换下
        **本该完全不变**，但分位阈值处的并列会让计数抖动 1 个 spot，
        零分布只取两个值，z 就恰好是 ±1.0——那不是效应量，是退化。

        所以这里用相对判据：sd 小于均值尺度的 1e-9，或算出来的 |z| 超过 1e6，
        都判为退化并返回 NaN。宁可缺一个数，也不要让一个假的极端值进表。
        """
        arr = np.array(perms, dtype=float)
        if arr.size == 0 or not np.all(np.isfinite(arr)):
            return float("nan")
        mu, sd = float(arr.mean()), float(arr.std())
        if sd <= 1e-9 * max(abs(mu), 1.0):
            return float("nan")
        # 零分布只取一两个值时，z 不是效应量而是"落在哪一档"的记号
        # （组成类指标就是这样恰好给出 ±1.0 的）。同样判退化。
        if np.unique(arr).size <= 2:
            return float("nan")
        z = (real - mu) / sd
        return float(z) if np.isfinite(z) and abs(z) < 1e6 else float("nan")

    rows = [("B_cell (min-cut)", b_real, b_perm),
            ("CAF density (composition)", dens_real, dens_perm),
            ("CAF neighbor enrichment", nb_real, nb_perm),
            ("Ripley's L peak (CAF)", L_real, L_perm),
            ("T-to-core distance (exclusion)", tdist_real, td_perm),
            ("T infiltration fraction", tinf_real, ti_perm)]
    print(f"  --- {sid} ---")
    for name, real, perms in rows:
        z = zstat(real, perms)
        print(f"{sid:<7}{name:<26}{real:>9.3f}{np.mean(perms):>9.3f}{z:>7.2f}")

    rho_dens = spearmanr(bfield, caf).statistic
    rho_tdist = spearmanr(bfield, tnk).statistic
    summary[sid] = dict(
        b_cell=b_real, z_b_cell=zstat(b_real, b_perm),
        dens=dens_real, z_dens=zstat(dens_real, dens_perm),
        nbr=nb_real, z_nbr=zstat(nb_real, nb_perm),
        ripleyL=L_real, z_ripley=zstat(L_real, L_perm),
        tdist=tdist_real, z_tdist=zstat(tdist_real, td_perm),
        tinf=tinf_real, z_tinf=zstat(tinf_real, ti_perm),
        rho_Bfield_vs_CAF=float(rho_dens),
        rho_Bfield_vs_TNK=float(rho_tdist))

with open(OUT / "benchmark_ext.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

print("\n=== 汇总：B_cell 是否携带独立拓扑信息 ===")
for sid, r in summary.items():
    topo = r["z_b_cell"]
    others = [r["z_dens"], r["z_nbr"], r["z_ripley"], r["z_tdist"], r["z_tinf"]]
    print(f"  {sid}: B_cell z={topo:+.2f} | 其余 z={['%+.2f'%v for v in others]} | "
          f"B~CAF ρ={r['rho_Bfield_vs_CAF']:+.2f} B~TNK ρ={r['rho_Bfield_vs_TNK']:+.2f}")
print(f"\n已写出 {OUT / 'benchmark_ext.json'}")
