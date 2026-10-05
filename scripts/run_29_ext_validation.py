#!/usr/bin/env python
"""External independent validation for the manuscript (2026-10-03 supplement).

Runs the two key confirmatory analyses — R3 spot-level decoupling and the
spectral phase-randomization null (M3) — on the external 10x Genomics public
Visium slides (breast cancer Block A S1/S2 + human lymph node), WITHOUT
touching the main result files (results/validation/decoupling.json and
spatial_null_check.json stay exactly as the 19-section main cohort).

Output: results/validation/ext_validation.json

The helper implementations below are copied verbatim from run_27_spatial_null.py
and run_06_validate.py so the external numbers are exactly comparable with the
main-cohort numbers.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from sparta.io_ import (Paths, load_config, load_graph, load_json, patient_map,  # noqa: E402
                        save_json, set_seed, stamp_run)
from sparta.barrier import (compute_b_cell_field, compute_b_meta,  # noqa: E402
                            compute_b_mab, scores_from_adata)
from sparta.validate import decoupling_stats  # noqa: E402

N_PERM = 500          # 与 run_27 完全一致（稿件 M3 声明 500 次置换）
SEED = 20261003


def _bh_adjust(p_values):
    """Benjamini–Hochberg adjusted p-values, returned in input order."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = np.clip(adjusted_sorted, 0.0, 1.0)
    return adjusted


# ---------------------------------------------------------------------------
# 以下三个 helper 与 run_27_spatial_null.py 逐字一致（2026-10-03 复制）。
# ---------------------------------------------------------------------------
def _rankz(v):
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
    ap = argparse.ArgumentParser(description="外部独立验证（decoupling + 空间空模型）")
    ap.add_argument("--slides", nargs="+",
                    default=["BRCA01", "BRCA02", "LN01"],
                    help="外部切片 ID（默认 3 片）")
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_seed(cfg["seed"])
    P = Paths(cfg)
    pmap = patient_map(P)
    rng = np.random.default_rng(SEED)

    rows = {}
    for sid in args.slides:
        adata = sc.read_h5ad(P.scored(sid))
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled: list[str] = []
        S = scores_from_adata(adata, missing_out=filled)
        with np.load(P.barrier(sid), allow_pickle=True) as z:
            b = {k: np.asarray(z[k]) for k in z.files}

        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"],
                            vessel, D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                                  **cfg["barrier"]["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           **cfg["barrier"]["b_mab"])["b_mab"]

        # R3 同款：decoupling_stats（control = 距血管图距离）
        d = decoupling_stats(bc, bm, q=cfg["validate"]["decouple_q"],
                             control=np.asarray(b["d_vessel_um"], float))
        dec = {k: v for k, v in d.items()
               if k not in ("idx_discordant", "idx_discordant_r")}

        # M3 同款：谱相位随机化空模型
        n = A.shape[0]
        real_rho = partial_spearman(bc, bm, dv)
        L = normalized_laplacian(A)
        eigenvalues, V = np.linalg.eigh(L.toarray())
        f_hat = V.T @ bm
        f_hat_abs = np.abs(f_hat)
        null_rhos = []
        for _ in range(N_PERM):
            signs = rng.choice([-1.0, 1.0], size=n)
            bm_null = V @ (signs * f_hat_abs)
            null_rhos.append(partial_spearman(bc, bm_null, dv))
        null_rhos = np.array(null_rhos)
        k_ge = int((null_rhos >= real_rho).sum())
        emp_p = (k_ge + 1) / (len(null_rhos) + 1)

        rows[sid] = dict(
            patient=pmap.get(sid, "?"),
            filled_keys=filled,
            degraded=bool(filled),
            decoupling=dec,
            spatial_null=dict(
                n_nodes=int(n),
                real_rho_partial=real_rho,
                null_mean=float(null_rhos.mean()),
                null_std=float(null_rhos.std()),
                null_ci95=[float(np.percentile(null_rhos, 2.5)),
                           float(np.percentile(null_rhos, 97.5))],
                empirical_p_one_sided=emp_p,
                n_perm=N_PERM,
            ),
        )
        print(f"[外部验证] {sid:<7} ρ={dec['rho']:+.3f} ρ_partial={dec['rho_partial']:+.3f} "
              f"(p_partial={dec['p_partial']:.2e}) 解离区 {dec['frac_discordant_r']*100:.1f}% "
              f"(随机 {dec['chance_discordant']*100:.2f}%) | 空模型 ρ={real_rho:+.3f} p={emp_p:.3f}")

    slide_ids = list(rows)
    ps = np.array([rows[sid]["spatial_null"]["empirical_p_one_sided"]
                   for sid in slide_ids], dtype=float)
    ps_bh = _bh_adjust(ps)
    for sid, p_adj in zip(slide_ids, ps_bh):
        rows[sid]["spatial_null"]["bh_padj"] = float(p_adj)
    rps = np.array([r["decoupling"]["rho_partial"] for r in rows.values()])
    pp = np.array([r["decoupling"].get("p_partial", 1.0) for r in rows.values()])
    n_deg = sum(1 for r in rows.values() if r["degraded"])
    summary = dict(
        n_slides=len(rows),
        n_decoupling_positive=int((rps > 0).sum()),
        n_decoupling_significant=int((pp < 0.05).sum()),
        n_null_significant=int((ps < 0.05).sum()),
        n_null_significant_bh=int((ps_bh < 0.05).sum()),
        median_rho_partial=float(np.median(rps)),
        median_null_p=float(np.median(ps)),
        n_degraded=n_deg,
        degraded_slides={sid: r["filled_keys"] for sid, r in rows.items() if r["degraded"]},
        p_correction="(k+1)/(N+1) finite-sample; BH-FDR across the three external sections",
        note=("External cross-dataset reproduction on 10x public Visium: two adjacent "
              "breast-cancer sections from one patient and one non-tumour lymph-node section. "
              "These sections are exploratory and are not pooled with the primary cohort."),
        caveat=("The breast sections are two adjacent sections from one patient; LN01 is a "
                "non-tumour control. All three were scored with the brca epithelial-marker "
                "profile for sink definition, which has not been histopathologically verified. "
                "The external arm supports cross-dataset reproduction of model-field coupling, "
                "not independent patient-level or functional validation."),
    )
    out = dict(summary=summary, per_slide=rows,
               run=stamp_run(cfg, {"module": "M29-ext-validation"}))
    p = P.validation("ext_validation.json")
    save_json(p, out)
    print(f"\n=== 外部验证汇总 ===")
    print(f"decoupling rho_partial 为正：{summary['n_decoupling_positive']}/{summary['n_slides']}"
          f"，中位 {summary['median_rho_partial']:+.3f}")
    print(f"空模型显著（p<0.05）：{summary['n_null_significant']}/{summary['n_slides']}"
          f"，中位 p={summary['median_null_p']:.3f}")
    print(f"BH-FDR 后显著（q<0.05）：{summary['n_null_significant_bh']}/{summary['n_slides']}")
    print(f"已写出 {p}")


if __name__ == "__main__":
    raise SystemExit(main())
