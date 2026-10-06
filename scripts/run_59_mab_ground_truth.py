#!/usr/bin/env python
"""
run_59_mab_ground_truth.py -- process-level ground truth for the molecular operator
=================================================================================
Output: results/validation/mab_ground_truth.json  (+ per-tissue table mab_ground_truth_tissues.csv)

Why this experiment exists
--------------------------
The cellular operator (minimum cut) is checked against an agent-based access
process in run_38, and against withheld CD8+ T-cell positions in the CODEX
cohort. Until v2.3 the molecular operator (screened diffusion-absorption field)
had neither: the manuscript stated "the diffusion operator has no such test in
this study". This script closes that gap in silico, on the same simulated
tissues as run_38, with a ground truth that differs from the operator in four
specific ways:

  1. it is TRANSIENT -- particles are censused at finite time windows
     (200 / 800 / 2,500 steps), whereas the operator solves the steady state;
  2. absorption is SATURABLE in one arm (40 binding sites per spot, a
     nonlinearity the operator does not have);
  3. it runs on the TRUE noise-free tissue, while the operator reads the
     noisy rank-normalised observed scores;
  4. it is a delivery-direction process (vessels -> tissue), whereas the
     operator's phi is the escape probability of the killed walk.

What it does NOT test (stated in the manuscript): the local transport law
itself. The particle process shares the conductance and absorption functional
forms of Eqs. (2)-(3) by construction, so agreement validates the operator's
global summary, its steady-state and linear-absorption approximations and its
robustness to noisy inputs -- not the biological correctness of the law.

Design
------
Same hexagonal lattices, geometries and observation model as run_38, plus two
planted molecular fields: crosslinking co-planted with the matrix-rich spots
over a wide background (so that a spot's within-section rank approximates its
planted value, keeping the operator's rank-based pore size consistent with the
true physics; see the comment above plant_molecular) and patchy ligand
expression planted on the tumour nest.  Per tissue, W particles start at the
vessel spots and perform a conductance-weighted walk: at each step the particle
is absorbed at its current spot with probability kappa_i/(d_i + kappa_i)
(the backward-equation form of the screened Poisson equation; vessel spots are
Dirichlet sources and never absorb), otherwise moves to a neighbour chosen
with probability proportional to the edge conductance of Eq. (2) at the
specified molecular radius.  Core delivery = fraction of the initial particles
alive inside the tumour core at the census time.

Summaries computed from the OBSERVED scores are compared with delivery across
tissues exactly as in run_38: SPARTA B_mAb at the core, peritumoural ECM mean,
global ECM mean, the B_cell field at the core (does the cellular operator also
predict molecular delivery?) and the vessel-nest distance.  Core clipping at
the solver's phi floor (B_mAb = -log 1e-12) is tracked per run and reported by
geometry: at the IgG radius a closed capsule of top-rank crosslinking drives
the steady state below the floor and the operator can no longer rank those
tissues -- a boundary behaviour of the steady-state formulation that the
transient process does not share.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

import run_38_synthetic_benchmark as base  # noqa: E402
from sparta.barrier import compute_b_cell_field, compute_b_mab, edge_pairs  # noqa: E402
from sparta.io_ import Paths, load_config, save_json, stamp_run  # noqa: E402

SEED = 20261006
RADII = [0.5, 2.0, 5.5]            # small molecule / nanobody-scale / IgG (nm)
CHECKPOINTS = [200, 800, 2500]     # census times (steps)
NOISE = 0.15                       # mid noise level of run_38
SAT_CAP = 40                       # binding sites per spot, saturating arm
GEOMS = base.GEOMS

_MAB_CACHE: dict | None = None


def mab_cfg() -> dict:
    global _MAB_CACHE
    if _MAB_CACHE is None:
        _MAB_CACHE = load_config(None)["barrier"]["b_mab"]
    return _MAB_CACHE


def mab_kwargs() -> dict:
    m = mab_cfg()
    return {k: m[k] for k in ("g0", "lam", "xi0_nm", "beta", "kd_eff", "kappa_w", "g_floor")}


# --------------------------------------------------------------------------
# planted molecular fields
# --------------------------------------------------------------------------
# Crosslinking is planted with a WIDE background (roughly uniform on 0.02-0.78)
# so that the within-section rank of a spot approximately equals its planted
# value.  This makes the operator's rank-based reading of the observed scores
# (pore size 20*exp(-3*rank) nm) consistent with the true physics the particles
# experience (pore size 20*exp(-3*value) nm).  A narrow background would let the
# ring occupy the top ~10% of ranks while its absolute crosslinking stays low,
# and the operator would fully exclude edges that the true physics leaves open.
# With this planting the radius scan spans the regimes seen in real sections:
# at 0.5 nm no edge is fully excluded (partial steric attenuation only); at
# 2 nm ring-internal edges are excluded; at the IgG radius 5.5 nm every edge
# whose mean crosslinking rank exceeds ln(20/5.5)/3 = 0.43 is fully excluded --
# about half of all edges, as in real sections -- and a CLOSED capsule of
# top-rank crosslinking drives the steady state below the solver's 1e-12 phi
# floor, where B_mAb saturates at -log(1e-12) and can no longer rank tissues
# (tracked per run via frac_core_clipped; transient delivery still differs).
XL_RING = 0.88


def plant_molecular(t, rng):
    """Crosslinking co-planted with the matrix spots; ligand on the tumour."""
    sn = lambda: np.clip(base.smooth_noise(t.xy, rng), -2, 2)
    xl = np.clip(0.40 + 0.30 * sn(), 0.02, 0.78)
    xl[t.fib] = np.clip(XL_RING + 0.06 * sn()[t.fib], 0.80, 0.98)
    ag = np.clip(0.04 + 0.04 * sn(), 0, 1)
    tum = np.flatnonzero(t.tumour)
    ag[tum] = np.clip(0.15 + 0.10 * sn()[tum], 0, 1)
    hot = rng.choice(tum, size=max(1, int(0.4 * len(tum))), replace=False)
    ag[hot] = np.clip(0.55 + 0.20 * sn()[hot], 0, 1)
    t.xl_true, t.ag_true = xl, ag
    return t


def observe_field(t, v, rng, noise):
    o = v + 0.5 * noise * base.smooth_noise(t.xy, rng, length=200) + noise * rng.standard_normal(len(v))
    return rankdata(o) / len(o)


# --------------------------------------------------------------------------
# ground truth: transient, (optionally) saturable particle delivery
# --------------------------------------------------------------------------
class WalkMachinery:
    """Neighbour lists and per-radius log-conductances for one tissue."""

    def __init__(self, t, r_nm, mab_cfg):
        A = t.A.tocsr()
        n = A.shape[0]
        deg = np.diff(A.indptr)
        maxd = deg.max()
        nb = np.full((n, maxd), -1, int)
        for i in range(n):
            nbrs = A.indices[A.indptr[i]:A.indptr[i + 1]]
            nb[i, :len(nbrs)] = nbrs

        pairs = edge_pairs(A, upper_only=True)
        pu, pv = pairs[:, 0], pairs[:, 1]
        ecm_e = 0.5 * (t.ecm_true[pu] + t.ecm_true[pv])
        xl_e = 0.5 * (t.xl_true[pu] + t.xl_true[pv])
        xi = mab_cfg["xi0_nm"] * np.exp(-mab_cfg["beta"] * xl_e)
        s = r_nm / np.maximum(xi, 1e-6)
        phi = np.where(s < 1.0, (1.0 - s) ** 2, 0.0)
        g = mab_cfg["g0"] * np.exp(-mab_cfg["lam"] * ecm_e) * phi + mab_cfg["g_floor"]

        Wnb = np.zeros((n, maxd))
        for u, v, w in zip(pu, pv, g):
            k_u = int(np.where(nb[u] == v)[0][0])
            k_v = int(np.where(nb[v] == u)[0][0])
            Wnb[u, k_u] = w
            Wnb[v, k_v] = w
        self.nb, self.Wnb = nb, Wnb
        self.Lnb = np.log(Wnb + 1e-300)
        self.d = Wnb.sum(axis=1)
        self.maxd = maxd


def kappa_true(t):
    m = mab_cfg()
    kap = t.ag_true / (t.ag_true + m["kd_eff"] + 1e-9)
    kap[t.vessel] = 0.0                      # Dirichlet sources never absorb
    return kap


def particle_delivery(mach, t, rng, W=6000, checkpoints=CHECKPOINTS, saturate=False, cap=SAT_CAP):
    """Run the particle process; returns censuses, cumulative visits, delivery."""
    n, maxd = mach.nb.shape
    kap = kappa_true(t)
    kill_p = kap / (mach.d + kap + 1e-300)   # kappa/(d+kappa): backward-equation killing
    pos = rng.choice(t.vessel, size=W)
    alive = np.ones(W, bool)
    bound = np.zeros(n)
    visits = np.zeros(n, int)
    census = {T: np.zeros(n, int) for T in checkpoints}
    in_core = np.zeros(n, bool)
    in_core[t.core] = True
    T_max = checkpoints[-1]
    for step in range(1, T_max + 1):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        u = pos[idx]
        if saturate:
            frac = np.clip(bound[u] / cap, 0, 1)
            keff = kap[u] * (1.0 - frac)
            kp = keff / (mach.d[u] + keff + 1e-300)
        else:
            kp = kill_p[u]
        die = rng.random(u.size) < kp
        if saturate and die.any():
            np.add.at(bound, u[die], 1)
        idx = idx[~die]
        u = pos[idx]
        gum = -np.log(-np.log(rng.random((u.size, maxd)) + 1e-300) + 1e-300)
        k = np.argmax(mach.Lnb[u] + gum, axis=1)
        v = mach.nb[u, k]
        pos[idx] = v
        np.add.at(visits, v, 1)
        if step in census:
            cens = np.zeros(n, int)
            np.add.at(cens, v, 1)
            census[step] = cens
    delivery = {T: float(census[T][t.core].sum() / W) for T in checkpoints}
    return delivery, visits


# --------------------------------------------------------------------------
# summaries from observed scores
# --------------------------------------------------------------------------
def observed_summaries(t, ecm_o, caf_o, xl_o, ag_o, r_nm):
    out = {}
    res = compute_b_mab(t.A, ecm_o, xl_o, ag_o, t.vessel, r_nm=r_nm, **mab_kwargs())
    b = res["b_mab"]
    out["b_mab_core"] = float(np.nanmedian(b[t.core]))
    out["b_mab_reach_mean"] = float(np.nanmean(b[res["reachable"]])) if res["reachable"].any() else np.nan
    out["frac_core_clipped"] = float(np.mean(b[t.core] >= -np.log(1e-12) - 1e-6))
    fld = compute_b_cell_field(t.A, ecm_o, caf_o, t.vessel, **base.CELL)["b_cell_field"]
    out["b_cell_field_core"] = float(np.median(fld[t.core]))
    peri = (t.r > t.r_t) & (t.r <= t.r_t + 400)
    out["ecm_peritumoural_mean"] = float(ecm_o[peri].mean())
    out["ecm_global_mean"] = float(ecm_o.mean())
    out["vessel_boundary_distance"] = float(np.median(t.r[t.vessel] - t.r_t))
    out["_b_mab_field"] = b
    out["_reachable"] = res["reachable"]
    return out


# barrier-like direction: larger summary -> less delivery
DIRECTION = dict(b_mab_core=+1, b_mab_reach_mean=+1, b_cell_field_core=+1,
                 ecm_peritumoural_mean=+1, ecm_global_mean=+1, vessel_boundary_distance=+1)
KEYS = list(DIRECTION)


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return None
    gt = (pos[:, None] > neg[None, :]).mean()
    eq = (pos[:, None] == neg[None, :]).mean()
    return float(gt + 0.5 * eq)


def evaluate(rows, radius, T):
    """Tissue-level evaluation at one radius and one census time."""
    sub = [r for r in rows if r["radius"] == radius]
    d = {}
    for k in KEYS:
        v = np.array([r[k] for r in sub], float) * DIRECTION[k]
        lost = np.array([1.0 - r["delivery"][T] for r in sub])
        ok = np.isfinite(v)
        closed = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "closed"]
        gap05 = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "gap05"]
        band = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "band"]
        scat = [DIRECTION[k] * r[k] for r in sub if r["geom"] == "scattered"]
        d[k] = dict(spearman_vs_lost_delivery=float(spearmanr(v[ok], lost[ok])[0]) if ok.sum() > 3 else np.nan,
                    auc_closed_vs_gap05=auc(closed, gap05),
                    auc_closed_vs_band=auc(closed, band),
                    auc_closed_vs_scattered=auc(closed, scat))
    # ground-truth separability: does delivery itself separate the contrasts?
    dl = lambda g: [-r["delivery"][T] for r in sub if r["geom"] == g]
    d["_gt"] = dict(auc_closed_vs_gap05=auc(dl("closed"), dl("gap05")),
                    auc_closed_vs_band=auc(dl("closed"), dl("band")))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--W", type=int, default=6000)
    ap.add_argument("--quick", action="store_true", help="smoke test: 2 geoms, 1 rep, few particles")
    args = ap.parse_args()

    cfg = load_config(None)
    P = Paths(cfg)
    mab = mab_cfg()

    rng = np.random.default_rng(SEED)
    xy, A = base.hex_lattice()
    geoms = ["closed", "gap05"] if args.quick else GEOMS
    reps = 1 if args.quick else args.reps
    W = 600 if args.quick else args.W
    checkpoints = [100, 400, 1000] if args.quick else CHECKPOINTS
    radii = [0.5, 5.5] if args.quick else RADII
    do_sat = not args.quick

    print(f"lattice: {len(xy)} spots; geoms={geoms} reps={reps} W={W} radii={radii}", flush=True)
    rows = []
    node_rhos = {r: [] for r in radii}
    sat_rows = []

    def run_one(t, r_nm, saturate=False):
        mach = WalkMachinery(t, r_nm, mab)
        delivery, visits = particle_delivery(mach, t, rng, W=W, checkpoints=checkpoints,
                                             saturate=saturate)
        return mach, delivery, visits

    t0 = time.time()
    for rep in range(reps):
        for geom in geoms:
            t = plant_molecular(base.make_tissue(xy, A, geom, rng), rng)
            ecm_o = observe_field(t, t.ecm_true, rng, NOISE)
            caf_o = observe_field(t, t.caf_true, rng, NOISE)
            xl_o = observe_field(t, t.xl_true, rng, NOISE)
            ag_o = observe_field(t, t.ag_true, rng, NOISE)
            for r_nm in radii:
                mach, delivery, visits = run_one(t, r_nm)
                s = observed_summaries(t, ecm_o, caf_o, xl_o, ag_o, r_nm)
                field, reach = s.pop("_b_mab_field"), s.pop("_reachable")
                nv = ~np.isin(np.arange(len(xy)), t.vessel)
                m = reach & nv & np.isfinite(field)
                rho_node = float(spearmanr(field[m], -visits[m].astype(float))[0])
                node_rhos[r_nm].append(rho_node)
                rows.append(dict(geom=geom, rep=rep, radius=r_nm, noise=NOISE,
                                 delivery=delivery, node_rho=rho_node, **s))
                if do_sat and r_nm == 5.5:
                    _, d_sat, _ = run_one(t, r_nm, saturate=True)
                    sat_rows.append(dict(geom=geom, rep=rep, radius=r_nm,
                                         delivery_lin=delivery, delivery_sat=d_sat,
                                         b_mab_core=s["b_mab_core"]))
            print(f"  rep{rep} {geom:<9} ({time.time() - t0:.0f} s)", flush=True)

    res = dict(n_tissues=len(rows), geoms=geoms, reps=reps, W=W,
               radii=radii, checkpoints=checkpoints, noise=NOISE,
               sat_cap=SAT_CAP, seed=SEED,
               delivery_by_radius={}, gt_separability={}, per_radius={}, by_checkpoint={})

    for r in radii:
        sub = [x for x in rows if x["radius"] == r]
        res["delivery_by_radius"][f"{r:g}"] = {
            "mean": float(np.mean([x["delivery"][checkpoints[-1]] for x in sub])),
            "by_geom": {g: float(np.mean([x["delivery"][checkpoints[-1]] for x in sub if x["geom"] == g]))
                        for g in geoms}}
        res["per_radius"][f"{r:g}"] = evaluate(rows, r, checkpoints[-1])
        res["per_radius"][f"{r:g}"]["node_rho_median"] = float(np.median(node_rhos[r]))
        for T in checkpoints:
            res["by_checkpoint"].setdefault(f"{T}", {})[f"{r:g}"] = {
                k: v["spearman_vs_lost_delivery"] for k, v in evaluate(rows, r, T).items() if k != "_gt"}

    if sat_rows:
        lost_sat = [1.0 - x["delivery_sat"][checkpoints[-1]] for x in sat_rows]
        lost_lin = [1.0 - x["delivery_lin"][checkpoints[-1]] for x in sat_rows]
        bm = [x["b_mab_core"] for x in sat_rows]
        rel = [(a - b) / max(b, 1e-9) for x in sat_rows
               for a, b in [(x["delivery_sat"][checkpoints[-1]], x["delivery_lin"][checkpoints[-1]])] if b > 1e-9]
        res["saturation"] = dict(
            r_nm=5.5, cap=SAT_CAP, n=len(sat_rows),
            rho_mab_vs_lost_delivery_saturating=float(spearmanr(bm, lost_sat)[0]),
            rho_mab_vs_lost_delivery_linear=float(spearmanr(bm, lost_lin)[0]),
            delivery_spearman_lin_vs_sat=float(spearmanr(lost_lin, lost_sat)[0]),
            median_rel_delivery_change=float(np.median(rel)) if rel else None)

    res["clip_fraction_main"] = {
        f"{r:g}": dict(
            overall=float(np.mean([x["frac_core_clipped"] for x in rows if x["radius"] == r])),
            by_geom={g: float(np.mean([x["frac_core_clipped"] for x in rows
                                       if x["radius"] == r and x["geom"] == g]))
                     for g in geoms})
        for r in radii}

    res["meta"] = stamp_run(cfg, {"module": "M59-mab-ground-truth", "seed": SEED,
                                  "shared_local_law": True,
                                  "note": "particle process shares Eqs. (2)-(3) by construction; "
                                          "tests steady-state + linear-absorption approximations "
                                          "and noisy-input robustness, not the law itself"})
    save_json(P.validation("mab_ground_truth.json"), res)
    flat = []
    for x in rows:
        for T in checkpoints:
            flat.append(dict(geom=x["geom"], rep=x["rep"], radius=x["radius"], T=T,
                             delivery=x["delivery"][T], node_rho=x["node_rho"],
                             frac_core_clipped=x["frac_core_clipped"],
                             **{k: x[k] for k in KEYS}))
    with open(P.validation("mab_ground_truth_tissues.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(flat[0]))
        w.writeheader()
        w.writerows(flat)

    print(f"\n{len(rows)} tissue-radius runs done in {time.time() - t0:.0f} s")
    for r in radii:
        pr = res["per_radius"][f"{r:g}"]
        print(f"r={r:>3} nm  rho(mab_core, lost delivery)={pr['b_mab_core']['spearman_vs_lost_delivery']:+.3f}  "
              f"AUC gap05={pr['b_mab_core']['auc_closed_vs_gap05']:.3f} (gt {pr['_gt']['auc_closed_vs_gap05']:.3f})  "
              f"node rho={pr['node_rho_median']:+.3f}  clip={res['clip_fraction_main'][f'{r:g}']['overall']:.2f}")
    if "saturation" in res:
        s = res["saturation"]
        print(f"saturation: rho lin={s['rho_mab_vs_lost_delivery_linear']:+.3f} "
              f"sat={s['rho_mab_vs_lost_delivery_saturating']:+.3f} "
              f"median rel change={s['median_rel_delivery_change']:+.3f}")


if __name__ == "__main__":
    main()
