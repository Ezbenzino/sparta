# M31 外部验证 benchmark：最小割 vs 局部统计量（R7 增补）
"""
与 run_10_benchmark_ext.py 完全同口径（六类指标、50 次空间重排、同一退化
判据），但只跑外部独立验证切片（默认 BRCA01 BRCA02 LN01），输出到独立的
results/validation/benchmark_ext_external.json，**不覆盖**主队列的
benchmark_ext.json。

回答的问题：在 10x 公开乳腺癌/淋巴结切片上，B_cell（最小割）是否仍然
携带与 CAF 密度、邻域富集、Ripley's L、T 排斥距离/浸润比例不同的拓扑
信息？若外部三片的 z_b_cell 仍显著而组成/局部统计量不敏感，则 R7 的
"管线未改、结论复现"延伸到了 benchmark 层面。

运行：python scripts/run_31_benchmark_ext_slides.py [SLIDE ...]
"""
from __future__ import annotations
import sys, json
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

SLIDES = sys.argv[1:] or ["BRCA01", "BRCA02", "LN01"]
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
    if len(t_high) == 0: return np.nan
    tree = cKDTree(coords[sink])
    d, _ = tree.query(coords[t_high])
    return float(np.mean(d))


def zstat(real, perms):
    """与 run_10 完全一致：sd 退化、|z|>1e6、零分布 ≤2 个取值均判 NaN。"""
    arr = np.array(perms, dtype=float)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        return float("nan")
    mu, sd = float(arr.mean()), float(arr.std())
    if sd <= 1e-9 * max(abs(mu), 1.0):
        return float("nan")
    if np.unique(arr).size <= 2:
        return float("nan")
    z = (real - mu) / sd
    return float(z) if np.isfinite(z) and abs(z) < 1e6 else float("nan")


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

with open(OUT / "benchmark_ext_external.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

print("\n=== 外部验证：B_cell 是否携带独立拓扑信息 ===")
for sid, r in summary.items():
    topo = r["z_b_cell"]
    others = [r["z_dens"], r["z_nbr"], r["z_ripley"], r["z_tdist"], r["z_tinf"]]
    print(f"  {sid}: B_cell z={topo:+.2f} | 其余 z={['%+.2f'%v for v in others]} | "
          f"B~CAF ρ={r['rho_Bfield_vs_CAF']:+.2f} B~TNK ρ={r['rho_Bfield_vs_TNK']:+.2f}")
print(f"\n已写出 {OUT / 'benchmark_ext_external.json'}（主队列 benchmark_ext.json 未改动）")
