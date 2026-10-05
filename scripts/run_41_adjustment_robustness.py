#!/usr/bin/env python
"""
run_41_adjustment_robustness.py -- does the field association depend on how vessel distance is adjusted?
=======================================================================================================
Inputs : data/interim/{sid}.nodes.npz, {sid}.graph.npz
Output : results/validation/adjustment_robustness.json

The primary statistic adjusts both fields for a quadratic function of the rank-
transformed weighted graph distance to the nearest vessel (average ranks for ties).
Graph distances on regular spot lattices take few distinct values, so the covariate
is heavily tied.  This script recomputes the association under five alternative
adjustments and reports how often the sign and magnitude are preserved:

  rank_quadratic   the locked primary statistic
  rank_cubic       adds the cubic rank term
  spline_raw       cubic B-spline of the raw distance (um), 4 interior quantile knots
  shell_stratified Spearman within distance shells (deciles of distance), pooled by
                   size-weighted Fisher z
  unadjusted       plain Spearman of the two fields
It also reports the effect of breaking the distance ties at random (jitter of 1e-9 um),
which is how a platform with different floating-point summation could perturb the
covariate.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy.interpolate import BSpline  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from sparta.barrier import compute_b_cell_field, compute_b_mab, compute_b_meta  # noqa: E402
from sparta.io_ import Paths, admitted_slides, load_config, load_graph, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402
from sparta.spatial_stats import partial_spearman, rankz  # noqa: E402


def resid_on(Xd, v):
    zv = rankz(v)
    beta, *_ = np.linalg.lstsq(Xd, zv, rcond=None)
    return zv - Xd @ beta


def corr_ranks(a, b):
    return float(spearmanr(a, b)[0])


def bspline_design(x, n_interior=4, k=3):
    x = np.asarray(x, float)
    qs = np.quantile(x, np.linspace(0, 1, n_interior + 2)[1:-1])
    qs = np.unique(qs)
    lo, hi = x.min(), x.max() + 1e-9
    t = np.r_[[lo] * (k + 1), qs, [hi] * (k + 1)]
    B = BSpline.design_matrix(np.clip(x, lo, hi - 1e-12), t, k).toarray()
    return B


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    bcfg = cfg["barrier"]
    rng = np.random.default_rng(20261006)
    per = {}
    for sid in admitted_slides(P):
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_nodes(nodes)
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **bcfg["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **bcfg["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **bcfg["b_mab"])["b_mab"]
        r = {}
        r["rank_quadratic"] = partial_spearman(bc, bm, dv)
        zc = rankz(dv)
        X3 = np.column_stack([np.ones_like(zc), zc, zc ** 2, zc ** 3])
        r["rank_cubic"] = corr_ranks(resid_on(X3, bc), resid_on(X3, bm))
        B = bspline_design(dv)
        r["spline_raw"] = corr_ranks(resid_on(B, bc), resid_on(B, bm))
        edges = np.unique(np.quantile(dv, np.linspace(0, 1, 11)))
        shell = np.clip(np.searchsorted(edges, dv, side="right") - 1, 0, len(edges) - 2)
        zs, ws = [], []
        for s in np.unique(shell):
            m = shell == s
            if m.sum() >= 10 and np.ptp(bc[m]) > 0 and np.ptp(bm[m]) > 0:
                rr = spearmanr(bc[m], bm[m])[0]
                zs.append(np.arctanh(np.clip(rr, -0.999, 0.999)))
                ws.append(m.sum() - 3)
        r["shell_stratified"] = float(np.tanh(np.average(zs, weights=ws))) if zs else None
        r["unadjusted"] = corr_ranks(bc, bm)
        jit = [partial_spearman(bc, bm, dv + 1e-9 * rng.standard_normal(len(dv))) for _ in range(20)]
        r["tie_jitter_min"], r["tie_jitter_median"], r["tie_jitter_max"] = (float(np.min(jit)), float(np.median(jit)), float(np.max(jit)))
        r["n_unique_distance"] = int(len(np.unique(dv)))
        r["n_spots"] = int(len(dv))
        per[sid] = r
        print(f"{sid:<7} quad={r['rank_quadratic']:+.3f} cubic={r['rank_cubic']:+.3f} spline={r['spline_raw']:+.3f} "
              f"shell={r['shell_stratified']:+.3f} unadj={r['unadjusted']:+.3f} jitter=[{r['tie_jitter_min']:+.3f},{r['tie_jitter_max']:+.3f}] "
              f"uniqueD={r['n_unique_distance']}/{r['n_spots']}", flush=True)
    prim = [s for s in per if not s.startswith(("BRCA", "LN"))]
    summ = {}
    for k in ["rank_quadratic", "rank_cubic", "spline_raw", "shell_stratified", "unadjusted", "tie_jitter_min", "tie_jitter_median"]:
        v = np.array([per[s][k] for s in prim], float)
        summ[k] = dict(median=float(np.median(v)), n_positive=int((v > 0).sum()), n=len(v),
                       min=float(v.min()), max=float(v.max()))
    save_json(P.validation("adjustment_robustness.json"),
              dict(per_slide=per, summary_primary=summ, meta=stamp_run(cfg, {"module": "M41-adjustment-robustness"})))
    for k, v in summ.items():
        print(f"{k:<18} median={v['median']:+.3f} positive {v['n_positive']}/{v['n']} range [{v['min']:+.3f},{v['max']:+.3f}]")


if __name__ == "__main__":
    main()
