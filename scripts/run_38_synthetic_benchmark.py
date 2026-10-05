#!/usr/bin/env python
"""
run_38_synthetic_benchmark.py -- planted-barrier simulation study with a process-level ground truth
===================================================================================================
Output: results/validation/synthetic_benchmark.json  (+ per-tissue table synthetic_benchmark_tissues.csv)

What is benchmarked
-------------------
On Visium-like hexagonal lattices (100 um pitch, 150 um radius graph) we plant
a tumour nest, peripheral vessels and the SAME number K of matrix-rich spots in
seven geometries:

    closed      one continuous capsule around the nest (touching the tumour edge)
    closed_far  the same capsule 300 um further out, separated from the nest by
                ordinary stroma (blocks access without touching the tumour)
    gap05..gap40  the same capsule with a contiguous angular gap (5-40 % of the
                circumference); removed spots are re-placed at random in the
                outer stroma, so K is unchanged
    band        K spots scattered at random inside a wider peritumoural band
                (same peritumoural matrix density as a capsule, no designed
                continuity; near the site-percolation threshold)
    patches     K spots in four compact stromal patches
    scattered   K spots uniformly at random in the stroma

Observed scores are rank-normalised (truth + smooth background + Gaussian noise)
at three noise levels, as for real signature scores.

Ground truth is a stochastic process that is NOT the min-cut: W agents start at
vessel spots and perform a lazy random walk on the TRUE (noise-free) tissue,
accepting a move u->v with probability c_uv = sigmoid(a - b*ECM - c*CAF) (the
SPARTA edge capacity, i.e. the same local permeability law).  The fraction of
agents that reach the tumour core within T steps is the ground-truth access.

Each summary is computed from the OBSERVED (noisy) scores:
    SPARTA B_cell (1 / max-flow), SPARTA field (median shortest-path cost at
    the core), global ECM mean, peritumoural ECM mean, tumour-fibrosis
    neighbourhood-enrichment z (label permutation), a BANKSY-style domain
    boundary coverage, and Ripley's L of matrix-rich spots.
The benchmark asks which summary tracks ground-truth access across geometries,
and which separates a closed capsule from an equally dense random band.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from sparta.barrier import compute_b_cell, compute_b_cell_field  # noqa: E402
from sparta.counterfactual import s2_ring_breaking  # noqa: E402
from sparta.io_ import Paths, load_config, save_json, stamp_run  # noqa: E402

SEED = 20261004
CELL = dict(a=3.0, b_ecm=8.0, c_caf=4.0)       # identical to configs/default.yaml
GEOMS = ["closed", "gap05", "gap10", "gap20", "gap40", "closed_far", "band", "patches", "scattered"]


# --------------------------------------------------------------------------
# tissue generator
# --------------------------------------------------------------------------
def hex_lattice(R=1500.0, h=100.0):
    dy = h * np.sqrt(3) / 2
    pts = []
    for i, y in enumerate(np.arange(-R, R + 1e-9, dy)):
        off = 0.0 if i % 2 == 0 else h / 2
        xs = np.arange(-R - h, R + h, h) + off
        for x in xs:
            if x * x + y * y <= R * R:
                pts.append((x, y))
    xy = np.array(pts)
    pairs = np.array(sorted(cKDTree(xy).query_pairs(r=1.5 * h)))
    n = len(xy)
    A = sp.coo_matrix((np.ones(2 * len(pairs)),
                       (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                      shape=(n, n)).tocsr()
    return xy, A


def smooth_noise(xy, rng, length=250.0, n_feat=60):
    """Cheap stationary smooth random field via random Fourier features."""
    W = rng.normal(0, 1.0 / length, size=(n_feat, 2))
    b = rng.uniform(0, 2 * np.pi, n_feat)
    f = np.sqrt(2.0 / n_feat) * np.cos(xy @ W.T + b).sum(axis=1)
    return (f - f.mean()) / (f.std() + 1e-12)


class Tissue:
    pass


def make_tissue(xy, A, geom, rng, r_t=500.0, r_ring=620.0, ring_w=160.0, n_vessel=14):
    n = len(xy)
    r = np.linalg.norm(xy, axis=1)
    ang = np.mod(np.arctan2(xy[:, 1], xy[:, 0]), 2 * np.pi)
    tumour = r <= r_t
    core = np.flatnonzero(r <= 0.6 * r_t)
    ring = np.flatnonzero(np.abs(r - r_ring) <= ring_w / 2)
    K = len(ring)
    stroma_outer = np.flatnonzero(r >= r_ring + ring_w / 2 + 100)
    vessel_pool = np.flatnonzero((r >= r_ring + 250) & (r <= np.max(r) - 50))
    vessel = rng.choice(vessel_pool, size=n_vessel, replace=False)
    avoid = set(vessel.tolist())

    if geom == "closed_far":
        ring = np.flatnonzero(np.abs(r - (r_ring + 300.0)) <= ring_w / 2)
        stroma_outer = np.flatnonzero(r >= r_ring + 300.0 + ring_w / 2 + 100)
        vessel_pool = np.flatnonzero((r >= r_ring + 300.0 + 200) & (r <= np.max(r) - 50))
        vessel = rng.choice(vessel_pool, size=n_vessel, replace=False)
        avoid = set(vessel.tolist())
        fib = ring
        K = len(ring)
    elif geom == "closed":
        fib = ring
    elif geom.startswith("gap"):
        frac = int(geom[3:]) / 100.0
        start = rng.uniform(0, 2 * np.pi)
        da = np.mod(ang[ring] - start, 2 * np.pi)
        keep = ring[da >= 2 * np.pi * frac]
        pool = np.setdiff1d(stroma_outer, list(avoid))
        extra = rng.choice(pool, size=K - len(keep), replace=False)
        fib = np.r_[keep, extra]
    elif geom == "band":
        band = np.flatnonzero((r > r_t + 30) & (r <= r_ring + 230))
        band = np.setdiff1d(band, list(avoid))
        fib = rng.choice(band, size=min(K, len(band)), replace=False)
    elif geom == "patches":
        pool = np.setdiff1d(np.flatnonzero(r > r_t + 120), list(avoid))
        centres = rng.choice(pool, size=4, replace=False)
        per = [K // 4 + (1 if i < K % 4 else 0) for i in range(4)]
        fib = []
        for c_, m in zip(centres, per):
            d = np.linalg.norm(xy[pool] - xy[c_], axis=1)
            cand = pool[np.argsort(d)]
            cand = [x for x in cand if x not in fib]
            fib += cand[:m]
        fib = np.array(fib, int)
    elif geom == "scattered":
        pool = np.setdiff1d(np.flatnonzero(r > r_t + 120), list(avoid))
        fib = rng.choice(pool, size=K, replace=False)
    else:
        raise ValueError(geom)

    fib = np.setdiff1d(np.asarray(fib, int), vessel)
    ecm = 0.12 + 0.06 * np.clip(smooth_noise(xy, rng), -2, 2)
    ecm[fib] = 0.90
    caf = 0.12 + 0.06 * np.clip(smooth_noise(xy, rng), -2, 2)
    caf[fib] = 0.85
    mal = np.where(tumour, 0.9, 0.1)
    t = Tissue()
    t.xy, t.A, t.geom = xy, A, geom
    t.ecm_true, t.caf_true, t.mal = np.clip(ecm, 0, 1), np.clip(caf, 0, 1), mal
    t.fib, t.K, t.core, t.vessel, t.tumour = fib, K, core, vessel, tumour
    t.r, t.r_t = r, r_t
    return t


def observe(t, rng, noise):
    """Rank-normalised observed scores (what an expression signature would give)."""
    def obs(v):
        o = v + 0.5 * noise * smooth_noise(t.xy, rng, length=200) + noise * rng.standard_normal(len(v))
        return rankdata(o) / len(o)
    return obs(t.ecm_true), obs(t.caf_true)


# --------------------------------------------------------------------------
# ground truth: agent-based lazy random walk on the TRUE permeability
# --------------------------------------------------------------------------
def access_ground_truth(t, rng, W=3000, T=1500):
    A = t.A.tocsr()
    n = A.shape[0]
    deg = np.diff(A.indptr)
    maxd = deg.max()
    nb = np.full((n, maxd), -1, int)
    for i in range(n):
        nbrs = A.indices[A.indptr[i]:A.indptr[i + 1]]
        nb[i, :len(nbrs)] = nbrs
    def cap(u, v):
        z = CELL["a"] - CELL["b_ecm"] * 0.5 * (t.ecm_true[u] + t.ecm_true[v]) \
            - CELL["c_caf"] * 0.5 * (t.caf_true[u] + t.caf_true[v])
        return 1.0 / (1.0 + np.exp(-z))
    pos = rng.choice(t.vessel, size=W)
    arrived = np.zeros(W, bool)
    t_arr = np.full(W, np.inf)
    is_core = np.zeros(n, bool)
    is_core[t.core] = True
    for step in range(1, T + 1):
        act = ~arrived
        if not act.any():
            break
        p = pos[act]
        k = (rng.random(p.size) * deg[p]).astype(int)
        q = nb[p, k]
        acc = rng.random(p.size) < cap(p, q)
        p = np.where(acc, q, p)
        pos[act] = p
        hit = is_core[p]
        idx = np.flatnonzero(act)[hit]
        arrived[idx] = True
        t_arr[idx] = step
    return float(arrived.mean()), float(np.median(t_arr) if arrived.any() else np.inf)


# --------------------------------------------------------------------------
# summaries computed from observed scores
# --------------------------------------------------------------------------
def ripley_L(xy_pts, area, r=200.0):
    m = len(xy_pts)
    if m < 2:
        return 0.0
    tree = cKDTree(xy_pts)
    cnt = tree.count_neighbors(tree, r=r) - m
    K = area * cnt / (m * (m - 1))
    return float(np.sqrt(K / np.pi) - r)


def summaries(t, ecm_o, caf_o, rng, top_q=None):
    out = {}
    bc = compute_b_cell(t.A, ecm_o, caf_o, t.vessel, t.core, **CELL)
    out["sparta_bcell"] = bc["b_cell"]
    fld = compute_b_cell_field(t.A, ecm_o, caf_o, t.vessel, **CELL)["b_cell_field"]
    out["sparta_field_core"] = float(np.median(fld[t.core]))
    out["ecm_global_mean"] = float(ecm_o.mean())
    peri = (t.r > t.r_t) & (t.r <= t.r_t + 400)
    out["ecm_peritumoural_mean"] = float(ecm_o[peri].mean())
    # matrix-rich labels: top fraction matching the planted fraction (oracle K)
    q = 1.0 - t.K / len(ecm_o) if top_q is None else top_q
    fibl = ecm_o >= np.quantile(ecm_o, q)
    # neighbourhood enrichment z: tumour-fibrotic contacts vs label permutation
    Au = sp.triu(t.A, k=1).tocoo()
    tum = t.tumour
    def contacts(lbl):
        return float(np.sum((tum[Au.row] & lbl[Au.col]) | (lbl[Au.row] & tum[Au.col])))
    obs = contacts(fibl)
    nonT = np.flatnonzero(~tum)
    null = []
    for _ in range(200):
        lbl = np.zeros_like(fibl)
        lbl[rng.choice(nonT, size=int(fibl[~tum].sum()), replace=False)] = True
        lbl[tum] = fibl[tum]
        null.append(contacts(lbl))
    null = np.array(null)
    out["nhood_enrichment_z"] = float((obs - null.mean()) / (null.std() + 1e-12))
    # BANKSY-style domains: own + neighbour-mean features, k-means k=3
    from sklearn.cluster import KMeans
    Dg = sp.diags(1.0 / np.maximum(np.asarray(t.A.sum(1)).ravel(), 1))
    F = np.column_stack([ecm_o, caf_o, t.mal])
    Fn = Dg @ (t.A @ F)
    Z = np.column_stack([np.sqrt(1 - 0.3) * F, np.sqrt(0.3) * Fn])
    lab = KMeans(n_clusters=3, n_init=4, random_state=0).fit_predict(Z)
    fib_dom = int(np.argmax([ecm_o[lab == k].mean() for k in range(3)]))
    # coverage: fraction of tumour-boundary spots with a fibrotic-domain neighbour
    bnd = np.flatnonzero(tum & (np.asarray(t.A @ (~tum).astype(float)).ravel() > 0))
    isfib = (lab == fib_dom).astype(float)
    nbfib = np.asarray(t.A @ isfib).ravel()
    # extend one hop outwards so a capsule one spot away still counts
    nb2 = np.asarray(t.A @ (nbfib > 0).astype(float)).ravel()
    out["domain_boundary_coverage"] = float(np.mean((nbfib[bnd] > 0) | (nb2[bnd] > 0)))
    area = np.pi * (np.max(t.r) ** 2)
    out["ripley_L_200um"] = ripley_L(t.xy[fibl], area)
    # geometry only (added 2026-10-05): median distance from the vessel sources to the nest boundary
    out["vessel_boundary_distance"] = float(np.median(t.r[t.vessel] - t.r_t))
    return out


# direction: +1 if larger summary means LESS access (barrier-like)
DIRECTION = dict(sparta_bcell=+1, sparta_field_core=+1, ecm_global_mean=+1,
                 ecm_peritumoural_mean=+1, nhood_enrichment_z=+1,
                 domain_boundary_coverage=+1, ripley_L_200um=+1, vessel_boundary_distance=+1)


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    gt = (pos[:, None] > neg[None, :]).mean()
    eq = (pos[:, None] == neg[None, :]).mean()
    return float(gt + 0.5 * eq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--noise", nargs="*", type=float, default=[0.05, 0.15, 0.30])
    ap.add_argument("--s2-reps", type=int, default=10)
    args = ap.parse_args()

    cfg = load_config(None)
    P = Paths(cfg)
    rng = np.random.default_rng(SEED)
    xy, A = hex_lattice()
    print(f"lattice: {len(xy)} spots, {A.nnz // 2} edges")
    rows = []
    t0 = time.time()
    for noise in args.noise:
        for geom in GEOMS:
            for rep in range(args.reps):
                t = make_tissue(xy, A, geom, rng)
                acc, tmed = access_ground_truth(t, rng)
                ecm_o, caf_o = observe(t, rng, noise)
                s = summaries(t, ecm_o, caf_o, rng)
                rows.append(dict(noise=noise, geom=geom, rep=rep, K=int(t.K),
                                 access=acc, t_median=tmed, **s))
            print(f"noise={noise:.2f} {geom:<9} access={np.mean([r['access'] for r in rows if r['geom']==geom and r['noise']==noise]):.3f}  "
                  f"({time.time() - t0:.0f} s)", flush=True)

    keys = list(DIRECTION)
    res = dict(n_tissues=len(rows), geoms=GEOMS, noise_levels=args.noise, reps=args.reps,
               lattice_spots=int(len(xy)), per_noise={}, pooled={})

    def evaluate(sub):
        acc = np.array([r["access"] for r in sub])
        d = {}
        for k in keys:
            v = np.array([r[k] for r in sub], float) * DIRECTION[k]
            ok = np.isfinite(v)
            rho = spearmanr(v[ok], -acc[ok])[0] if ok.sum() > 3 else np.nan
            closed = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "closed"]
            band = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "band"]
            scat = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "scattered"]
            g05 = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "gap05"]
            far = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "closed_far"]
            d[k] = dict(spearman_vs_lost_access=float(rho),
                        auc_closed_vs_band=auc(closed, band),
                        auc_closed_vs_scattered=auc(closed, scat),
                        auc_closed_vs_gap05=auc(closed, g05),
                        auc_closedfar_vs_scattered=auc(far, scat) if far else None,
                        auc_closedfar_vs_band=auc(far, band) if far else None)
        return d

    for noise in args.noise:
        res["per_noise"][f"{noise:g}"] = evaluate([r for r in rows if r["noise"] == noise])
    res["pooled"] = evaluate(rows)
    res["access_by_geom"] = {g: dict(mean=float(np.mean([r["access"] for r in rows if r["geom"] == g])),
                                     sd=float(np.std([r["access"] for r in rows if r["geom"] == g])))
                             for g in GEOMS}

    # synthetic S2: contiguous vs selection-matched scattered removal on closed capsules
    s2 = []
    for rep in range(args.s2_reps):
        t = make_tissue(xy, A, "closed", rng)
        ecm_o, caf_o = observe(t, rng, 0.15)
        o = s2_ring_breaking(A, ecm_o, caf_o, t.vessel, t.core, CELL, coords=xy,
                             k_frac=(0.20,), n_rand=40, seed=int(rng.integers(1e9)),
                             low_q=0.05, match_selection=True, n_match_groups=30)
        kk = next(iter(o["per_k"]))
        s2.append(dict(ratio=o["per_k"][kk]["ratio_vs_in_cut"], p=o["per_k"][kk]["p_vs_in_cut"],
                       n_cut=o["n_cut_nodes"]))
    res["synthetic_s2_closed"] = dict(median_ratio=float(np.median([x["ratio"] for x in s2])),
                                      ratios=[x["ratio"] for x in s2], p=[x["p"] for x in s2],
                                      n_cut=[x["n_cut"] for x in s2])
    res["meta"] = stamp_run(cfg, {"module": "M38-synthetic-benchmark", "seed": SEED})
    save_json(P.validation("synthetic_benchmark.json"), res)
    import csv
    with open(P.validation("synthetic_benchmark_tissues.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\npooled:")
    for k, v in res["pooled"].items():
        print(f"  {k:<26} rho={v['spearman_vs_lost_access']:+.3f}  AUC closed/band={v['auc_closed_vs_band']:.3f}  "
              f"closed/scattered={v['auc_closed_vs_scattered']:.3f} closed/gap05={v['auc_closed_vs_gap05']:.3f}")
    print("access by geometry:", {g: round(v["mean"], 3) for g, v in res["access_by_geom"].items()})
    print("synthetic S2 median ratio:", res["synthetic_s2_closed"]["median_ratio"])


if __name__ == "__main__":
    main()
