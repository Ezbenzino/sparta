#!/usr/bin/env python
"""
run_39_coupling_decomposition.py -- how much of the B_cell/B_mAb field coupling is built into the model?
======================================================================================================
Inputs : data/interim/{sid}.nodes.npz, {sid}.graph.npz
Outputs: results/validation/shared_input_spatial_null.json   (ECM-ablated association + its spatial null)
         results/validation/geometry_null.json               (construction / geometry nulls)

Three sources can make the two model fields co-vary after vessel-distance adjustment:
  (i)   anchoring geometry  -- both fields grow away from the same vessel/source set;
  (ii)  model construction  -- both operators read the same ECM score, and the
        CAF score that enters B_cell overlaps the ECM gene set;
  (iii) co-arrangement of the operator-specific inputs in tissue (CAF in B_cell;
        crosslinking and ligand absorption in B_mAb).
This script separates them with two input-level spectral nulls and one ablation:

  geometry null   B_mAb recomputed from graph-spectral surrogates of ALL its inputs
                  (ECM*, crosslink*, Ag*), real B_cell field kept;   -> (i)
  construction    B_mAb recomputed from the REAL ECM score but surrogate crosslink*
  null            and Ag*, real B_cell field kept;                    -> (i)+(ii)
  ECM ablation    b_ecm = 0, c_caf = 12 (configs/no_shared_ecm.yaml) and lam = 0;
                  association re-tested with the same 500-draw graph-spectral null
                  used for the primary analysis.
Input surrogates keep each input's graph power spectrum (sign randomisation of
|V^T x|) and are rank-normalised back to [0, 1] exactly like the real *_n scores.
The observed association minus the construction-null mean is the part that the
model structure alone does not produce.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from sparta.barrier import compute_b_cell_field, compute_b_mab, compute_b_meta  # noqa: E402
from sparta.io_ import (Paths, admitted_slides, load_config, load_graph, patient_map,  # noqa: E402
                        save_json, stamp_run)
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402
from sparta.spatial_stats import (control_projector, partial_spearman,  # noqa: E402
                                  partial_spearman_many, spectral_basis)

SEED = 20261005
PRIMARY = None  # filled from ledger


def rank01(x):
    return (rankdata(x) - 1.0) / (len(x) - 1.0)


def input_surrogate(x, V, rng):
    """Spectral sign surrogate of an input score, rank-normalised to [0, 1]."""
    x = np.asarray(x, float)
    if np.ptp(x) == 0:            # neutral-filled (missing) input stays constant
        return x.copy()
    fh = np.abs(V.T @ x)
    return rank01(V @ (rng.choice([-1.0, 1.0], size=len(fh)) * fh))


def bh(ps):
    ps = np.asarray(ps, float)
    m = len(ps)
    order = np.argsort(ps)
    adj = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = m - rank
        adj[idx] = min(prev, ps[idx] * m / (i + 1))
        prev = adj[idx]
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slides", nargs="*", default=None)
    ap.add_argument("--n-null", type=int, default=200, help="draws per input-level null")
    ap.add_argument("--n-sur", type=int, default=500, help="spatial surrogates for the ablation test")
    args = ap.parse_args()

    cfg = load_config(None)
    P = Paths(cfg)
    pmap = patient_map(P)
    bcfg = cfg["barrier"]
    cell_abl = dict(bcfg["b_cell"], b_ecm=0.0, c_caf=12.0)
    mab_abl = dict(bcfg["b_mab"], lam=0.0)
    slides = args.slides or admitted_slides(P)
    rng = np.random.default_rng(SEED)

    abl, geo = {}, {}
    for sid in slides:
        t0 = time.time()
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **bcfg["b_meta"])["d_vessel_um"]
        H = control_projector(dv)
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **bcfg["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **bcfg["b_mab"])["b_mab"]
        rho_obs = partial_spearman(bc, bm, dv)
        w, V = spectral_basis(A)
        n = len(bm)

        # ---------------- ECM ablation + spatial null -----------------
        bc0 = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cell_abl)["b_cell_field"]
        bm0 = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **mab_abl)["b_mab"]
        rho_abl = partial_spearman(bc0, bm0, dv)
        fh0 = np.abs(V.T @ bm0)
        null_abl = np.empty(args.n_sur)
        for k in range(args.n_sur):
            null_abl[k] = partial_spearman_many(bc0, V @ (rng.choice([-1.0, 1.0], size=n) * fh0), dv, H)[0]
        p_abl = (np.sum(null_abl >= rho_abl) + 1.0) / (args.n_sur + 1.0)
        abl[sid] = dict(patient=pmap.get(sid, sid), n_nodes=int(n), filled_keys=filled,
                        rho_partial_full_model=rho_obs, real_rho_partial=rho_abl,
                        retained_fraction=rho_abl / rho_obs if rho_obs != 0 else None,
                        null_mean=float(null_abl.mean()), null_std=float(null_abl.std()),
                        empirical_p_one_sided=float(p_abl), n_perm=args.n_sur)

        # ---------------- input-level nulls ---------------------------
        rho_geo = np.empty(args.n_null)
        rho_con = np.empty(args.n_null)
        for k in range(args.n_null):
            xl_s = input_surrogate(S["crosslink"], V, rng)
            ag_s = input_surrogate(S["ag_target"], V, rng)
            ecm_s = input_surrogate(S["ecm"], V, rng)
            bm_g = compute_b_mab(A, ecm_s, xl_s, ag_s, vessel, **bcfg["b_mab"])["b_mab"]
            bm_c = compute_b_mab(A, S["ecm"], xl_s, ag_s, vessel, **bcfg["b_mab"])["b_mab"]
            rho_geo[k] = partial_spearman_many(bc, bm_g, dv, H)[0]
            rho_con[k] = partial_spearman_many(bc, bm_c, dv, H)[0]
        p_geo = (np.sum(rho_geo >= rho_obs) + 1.0) / (args.n_null + 1.0)
        p_con = (np.sum(rho_con >= rho_obs) + 1.0) / (args.n_null + 1.0)
        cor = lambda a, b: float(spearmanr(S[a], S[b])[0]) if np.ptp(S[a]) > 0 and np.ptp(S[b]) > 0 else None  # noqa: E731
        geo[sid] = dict(patient=pmap.get(sid, sid), n_nodes=int(n), filled_keys=filled,
                        rho_obs=rho_obs,
                        geometry_null_mean=float(rho_geo.mean()), geometry_null_sd=float(rho_geo.std()),
                        geometry_null_ci95=[float(np.percentile(rho_geo, 2.5)), float(np.percentile(rho_geo, 97.5))],
                        p_vs_geometry_null=float(p_geo),
                        construction_null_mean=float(rho_con.mean()), construction_null_sd=float(rho_con.std()),
                        construction_null_ci95=[float(np.percentile(rho_con, 2.5)), float(np.percentile(rho_con, 97.5))],
                        p_vs_construction_null=float(p_con),
                        excess_over_geometry_null=float(rho_obs - rho_geo.mean()),
                        excess_over_construction_null=float(rho_obs - rho_con.mean()),
                        share_reproduced_by_construction=float(rho_con.mean() / rho_obs) if rho_obs > 0 else None,
                        input_correlations=dict(ecm_caf=cor("ecm", "caf"), ecm_crosslink=cor("ecm", "crosslink"),
                                                caf_crosslink=cor("caf", "crosslink"), ecm_ag=cor("ecm", "ag_target")),
                        n_null=args.n_null)
        print(f"{sid:<7} obs={rho_obs:+.3f} | ablated={rho_abl:+.3f} p={p_abl:.3f} | "
              f"geom null {rho_geo.mean():+.3f}±{rho_geo.std():.3f} (p={p_geo:.3f}) | "
              f"constr null {rho_con.mean():+.3f}±{rho_con.std():.3f} (p={p_con:.3f}) | "
              f"r(ECM,CAF)={geo[sid]['input_correlations']['ecm_caf']:.2f} ({time.time()-t0:.0f}s)", flush=True)
        save_json(P.validation("geometry_null.json"), dict(per_slide=geo, partial=True))
        save_json(P.validation("shared_input_spatial_null.json"), dict(per_slide=abl, partial=True))

    def summarise(d, primary):
        keys = [s for s in d if (s in primary)]
        return keys

    with open(P.root / "data" / "ledger.csv", encoding="utf-8") as f:
        import csv
        ext = {r["slide_id"] for r in csv.DictReader(f) if r.get("cancer_type") == "other"}
    prim = [s for s in abl if s not in ext]
    q_abl = bh([abl[s]["empirical_p_one_sided"] for s in prim])
    for s, q in zip(prim, q_abl):
        abl[s]["bh_padj"] = float(q)
    q_con = bh([geo[s]["p_vs_construction_null"] for s in prim])
    q_geo = bh([geo[s]["p_vs_geometry_null"] for s in prim])
    for s, a, b in zip(prim, q_con, q_geo):
        geo[s]["bh_q_construction"] = float(a)
        geo[s]["bh_q_geometry"] = float(b)
    med = lambda v: float(np.median(v))  # noqa: E731
    abl_sum = dict(n_slides=len(prim), median_full=med([abl[s]["rho_partial_full_model"] for s in prim]),
                   median_ablated=med([abl[s]["real_rho_partial"] for s in prim]),
                   n_positive=int(sum(abl[s]["real_rho_partial"] > 0 for s in prim)),
                   n_p_lt_05=int(sum(abl[s]["empirical_p_one_sided"] < 0.05 for s in prim)),
                   n_bh_lt_05=int(sum(abl[s]["bh_padj"] < 0.05 for s in prim)),
                   median_retained_fraction=med([abl[s]["retained_fraction"] for s in prim]),
                   config=dict(b_cell=cell_abl, b_mab=mab_abl))
    geo_sum = dict(n_slides=len(prim),
                   median_obs=med([geo[s]["rho_obs"] for s in prim]),
                   median_geometry_null_mean=med([geo[s]["geometry_null_mean"] for s in prim]),
                   median_construction_null_mean=med([geo[s]["construction_null_mean"] for s in prim]),
                   median_excess_over_construction=med([geo[s]["excess_over_construction_null"] for s in prim]),
                   n_excess_construction_positive=int(sum(geo[s]["excess_over_construction_null"] > 0 for s in prim)),
                   n_p_construction_lt_05=int(sum(geo[s]["p_vs_construction_null"] < 0.05 for s in prim)),
                   n_bh_construction_lt_05=int(sum(geo[s]["bh_q_construction"] < 0.05 for s in prim)),
                   n_p_geometry_lt_05=int(sum(geo[s]["p_vs_geometry_null"] < 0.05 for s in prim)),
                   n_bh_geometry_lt_05=int(sum(geo[s]["bh_q_geometry"] < 0.05 for s in prim)),
                   median_share_reproduced_by_construction=med([geo[s]["share_reproduced_by_construction"] for s in prim
                                                                if geo[s]["share_reproduced_by_construction"] is not None]),
                   median_r_ecm_caf=med([geo[s]["input_correlations"]["ecm_caf"] for s in prim]),
                   median_r_ecm_crosslink=med([geo[s]["input_correlations"]["ecm_crosslink"] for s in prim]),
                   external=[s for s in geo if s in ext])
    meta = stamp_run(cfg, {"module": "M39-coupling-decomposition", "seed": SEED})
    save_json(P.validation("shared_input_spatial_null.json"), dict(per_slide=abl, summary=abl_sum, meta=meta))
    save_json(P.validation("geometry_null.json"), dict(per_slide=geo, summary=geo_sum, meta=meta,
                                                       design=__doc__.split("This script separates")[1][:1500]))
    print("\nablation:", abl_sum)
    print("decomposition:", geo_sum)


if __name__ == "__main__":
    main()
