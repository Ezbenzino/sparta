#!/usr/bin/env python
"""
run_27_spatial_null.py — 保留空间自相关的偏相关空模型
=====================================================================
审稿人指出：R3 的 rho_partial 用标准 Spearman p 值，把每个 spot 当独立观察，
但 spot 在空间上自相关，有效样本量远小于 spot 数，p 值过于乐观。

做法（graph spectral phase randomization）：
  1. 对每张切片的图，算归一化拉普拉斯 L = I - D^{-1/2} A D^{-1/2}
  2. 特征分解 L = V Λ V^T
  3. 把 B_mAb 投影到谱域: f_hat = V^T B_mAb
  4. 保留幅度 |f_hat|（即空间自相关结构/variogram），随机化符号:
     f_hat_null = random_sign * |f_hat|
  5. 重构 B_mAb_null = V f_hat_null
     —— 这个 null 场和真实 B_mAb 有相同的空间自相关、相同的均值方差，
        但和 B_cell 的空间对齐被随机化了
  6. 在 null 场上重新算 rho_partial（控制 d_vessel），重复 N_PERM 次
  7. 经验 p = (null rho >= real rho 的次数 + 1) / (N_PERM + 1)
     —— 有限置换修正：纯 p=0 在 500 次下写成 1/501 ≈ 0.002，
        而不是 0（后者暗示"无限次置换后仍无一超越"）。

输出：results/validation/spatial_null_check.json
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.barrier import (compute_b_cell_field, compute_b_meta,  # noqa: E402
                            compute_b_mab, scores_from_adata)

N_PERM = 500          # 与稿件 M3 声明的 500 次置换口径一致
SEED = 20261003


def _rankz(v):
    from scipy.stats import rankdata
    return (rankdata(v) - np.mean(rankdata(v))) / np.std(rankdata(v))


def _residualize(v, control):
    """和 sparta.validate._residualize 完全一致：秩空间二次回归取残差。"""
    zv, zc = _rankz(v), _rankz(control)
    X = np.column_stack([np.ones_like(zc), zc, zc ** 2])
    beta, *_ = np.linalg.lstsq(X, zv, rcond=None)
    return zv - X @ beta


def partial_spearman(x, y, control):
    """和 decoupling_stats 的 rho_partial 同口径。"""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(control)
    x, y, c = x[ok], y[ok], control[ok]
    rx, ry = _residualize(x, c), _residualize(y, c)
    rho, _ = spearmanr(rx, ry)
    return float(rho)


def normalized_laplacian(A):
    """归一化图拉普拉斯 L = I - D^{-1/2} A D^{-1/2}，稀疏矩阵。"""
    import scipy.sparse as sp
    A = A.tocsr()
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).flatten()
    deg_safe = np.where(deg > 0, deg, 1.0)
    d_inv_sqrt = 1.0 / np.sqrt(deg_safe)
    D_inv_sqrt = sp.diags(d_inv_sqrt)
    L = sp.eye(n) - D_inv_sqrt @ A @ D_inv_sqrt
    return L


def main():
    cfg = load_config(None)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)
    rng = np.random.default_rng(SEED)

    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    with open(ledger, encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f)
                  if r.get("status") == "ingested"]

    print(f"Spatial null (spectral phase randomization)  "
          f"N={N_PERM} per slide, seed={SEED}\n")

    rows = {}
    for sid in slides:
        try:
            adata = sc.read_h5ad(P.scored(sid))
            A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        except FileNotFoundError as e:
            print(f"{sid}: SKIP ({e})")
            continue
        filled: list[str] = []
        S = scores_from_adata(adata, missing_out=filled)

        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"],
                            vessel, D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                  **cfg["barrier"]["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           **cfg["barrier"]["b_mab"])["b_mab"]

        n = A.shape[0]
        # 真实 rho_partial
        real_rho = partial_spearman(bc, bm, dv)

        # 图拉普拉斯特征分解（稠密，n 最大 ~2500，可接受）
        L = normalized_laplacian(A)
        L_dense = L.toarray()
        eigenvalues, V = np.linalg.eigh(L_dense)
        # eigenvalues 升序，V[:, i] 是第 i 个特征向量

        # 把 B_mAb 投影到谱域
        f_hat = V.T @ bm
        f_hat_abs = np.abs(f_hat)  # 保留功率谱/自相关结构

        null_rhos = []
        for _ in range(N_PERM):
            signs = rng.choice([-1.0, 1.0], size=n)
            f_hat_null = signs * f_hat_abs
            bm_null = V @ f_hat_null
            rho_null = partial_spearman(bc, bm_null, dv)
            null_rhos.append(rho_null)
        null_rhos = np.array(null_rhos)

        # 有限置换修正：空分布中 >= 真实值的个数 k，p = (k+1)/(N+1)。
        # 不加 +1 会把"500 次无一超越"写成 p=0，夸大证据强度。
        k_ge = int((null_rhos >= real_rho).sum())
        emp_p = (k_ge + 1) / (len(null_rhos) + 1)
        rows[sid] = dict(
            patient=pmap.get(sid, "?"),
            cohort="CSCC" if sid.startswith("CSCC") else "MEL",
            n_nodes=int(n),
            real_rho_partial=real_rho,
            null_mean=float(null_rhos.mean()),
            null_std=float(null_rhos.std()),
            null_ci95=[float(np.percentile(null_rhos, 2.5)),
                       float(np.percentile(null_rhos, 97.5))],
            empirical_p_one_sided=emp_p,
            n_perm=N_PERM,
            # 哪些签名被中性 0.5 填充。ag_target 被填会让 B_mAb 的结合吸收项退化成
            # 常数（k_abs 处处相同），efflux 被填会让 B_meta 的外排项退化——这两种
            # 切片上的 rho_partial 不是"生物学阴性"，而是"数据缺失"，必须在产物里标出来。
            filled_keys=filled,
            degraded=bool(filled),
        )
        tag = f"  [DEGRADED: {','.join(filled)}]" if filled else ""
        print(f"{sid:<8} real={real_rho:+.3f}  "
              f"null={null_rhos.mean():+.3f}±{null_rhos.std():.3f}  "
              f"p={emp_p:.3f}{tag}")

    # 汇总
    reals = np.array([r["real_rho_partial"] for r in rows.values()])
    nulls = np.array([r["null_mean"] for r in rows.values()])
    ps = np.array([r["empirical_p_one_sided"] for r in rows.values()])
    # BH-FDR 校正：19 张切片是一个检验族。用户 2026-10-03 指出 raw p<0.05=15/19
    # 未校正，BH 后应为 14/19（CSCC10 raw=0.048 → adj=0.061）。
    m = len(ps)
    order = np.argsort(ps)
    bh_adj = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = m - rank
        bh_adj[idx] = min(prev, ps[idx] * m / (i + 1))
        prev = bh_adj[idx]
    for sid, a in zip(rows.keys(), bh_adj):
        rows[sid]["bh_padj"] = float(a)

    summary = dict(
        median_real=float(np.median(reals)),
        median_null=float(np.median(nulls)),
        median_empirical_p=float(np.median(ps)),
        n_slides=len(rows),
        n_per_significant=int((ps < 0.05).sum()),
        n_sig_bh_p05=int((bh_adj < 0.05).sum()),
        n_patients=int(len(set(r["patient"] for r in rows.values()))),
        n_degraded=sum(1 for r in rows.values() if r["degraded"]),
        degraded_slides={sid: r["filled_keys"] for sid, r in rows.items() if r["degraded"]},
        p_correction="(k+1)/(N+1) finite-sample per slide; BH-FDR across the 19 slides",
    )
    print("\n=== Summary ===")
    print(f"Median real rho_partial: {summary['median_real']:+.3f}")
    print(f"Median null rho_partial (phase-randomized B_mAb): {summary['median_null']:+.3f}")
    print(f"Median one-sided empirical p: {summary['median_empirical_p']:.3f}")
    print(f"Slides raw p<0.05: {summary['n_per_significant']}/{summary['n_slides']}")
    print(f"Slides BH-FDR adj<0.05: {summary['n_sig_bh_p05']}/{summary['n_slides']}")
    if summary["n_degraded"]:
        print(f"⚠ Degraded slides (中性填充, 结果不可直接当生物学阴性): "
              f"{summary['n_degraded']}/{summary['n_slides']}")
        for sid, keys in summary["degraded_slides"].items():
            print(f"    {sid}: filled={keys}")

    out = dict(per_slide=rows, summary=summary,
               seed=SEED, n_perm=N_PERM,
               meta=stamp_run(cfg, {"module": "M27-spatial-null"}))
    p = P.validation("spatial_null_check.json")
    save_json(p, out)
    print(f"\nWrote {p}")


if __name__ == "__main__":
    main()
