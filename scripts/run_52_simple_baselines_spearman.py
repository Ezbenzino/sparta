#!/usr/bin/env python
"""
run_52_simple_baselines_spearman.py -- SPARTA's two operators versus three simple spatial summaries
==================================================================================================
Inputs : data/interim/{sid}.graph.npz, {sid}.nodes.npz, {sid}.barrier.npz  (30 sections: 19 primary,
         3 external, 8 replication)
Outputs: results/validation/simple_baselines_spearman.json

Question answered: how much of what the two transport operators report could be read off three
simple summaries that need no transport model?

  stromal density        mean (ECM + CAF) rank score over a spot and its graph neighbours
                         (section level: over non-tumour spots adjacent to the tumour)
  distance to tumour     signed Euclidean distance to the tumour boundary (negative inside the tumour)
                         (section level: median distance from vessel spots to the tumour)
  niche enrichment       Squidpy-style neighbourhood enrichment of tumour-stroma contacts
                         (spot level: binomial z of the number of stromal neighbours; section level:
                         z of tumour-stroma contacts against 1,000 label permutations)

Two levels are reported:
  spot level     Spearman rho, within each section, between each operator field (B_cell migration-cost
                 field; B_mAb = -log penetration) and each baseline;
  section level  Spearman rho across sections between operator summaries (log relative cellular
                 barrier B_rel = F_open / F; median of each field over tumour-core spots) and the
                 section-level baselines.

Spots that cannot be reached from a source or vessel (other graph components) carry placeholder field
values and are excluded from every correlation. Compartments use the pipeline's own definitions (tumour = malignant score >= 70th percentile;
stroma versus immune among the remaining spots by the larger of the mean (ECM, CAF) and mean
(T/NK, myeloid, B/plasma) rank scores). No outcome is used here; the outcome-based comparison is
in run_38 (simulated ground truth) and run_47 (measured CD8+ T cells).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "is_figures"))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.stats import rankdata, spearmanr  # noqa: E402

from sparta.barrier import compute_b_cell, compute_b_cell_field  # noqa: E402
from sparta.cellgraph import contact_enrichment_z  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph, patient_map, save_json, stamp_run  # noqa: E402
from sparta.node_tables import coords_from_nodes, load_nodes, scores_from_nodes  # noqa: E402

SEED = 20261005
N_PERM = 1000
Q_MALIG = 0.70
BASELINES = ("stromal_density", "dist_tumour_boundary", "niche_z")
FIELDS = ("b_cell_field", "b_mab_field")


def rank01(v):
    v = np.asarray(v, float)
    return (rankdata(v) - 1) / max(len(v) - 1, 1)


def compartments(nodes):
    """tumour / stroma / immune labels (0 / 1 / 2) from the node-table scores."""
    def col(name):
        raw = nodes.get(f"obs__{name}_n")
        if raw is None and nodes.get(f"obs__{name}") is not None:
            raw = rank01(nodes[f"obs__{name}"])
        return None if raw is None else np.asarray(raw, float)
    mal = col("Malignant")
    tum = mal >= np.quantile(mal, Q_MALIG)
    stroma = np.nanmean(np.column_stack([c for c in (col("ECM_core"), col("CAF")) if c is not None]), axis=1)
    imm_cols = [c for c in (col("T_NK"), col("Myeloid"), col("B_Plasma")) if c is not None]
    immune = np.nanmean(np.column_stack(imm_cols), axis=1)
    lab = np.where(tum, 0, np.where(stroma >= immune, 1, 2))
    return lab, stroma


def signed_distance(xy, tum):
    d_in = cKDTree(xy[~tum]).query(xy, k=1)[0] if (~tum).any() else np.zeros(len(xy))
    d_out = cKDTree(xy[tum]).query(xy, k=1)[0] if tum.any() else np.zeros(len(xy))
    return np.where(tum, -d_in, d_out)


def spot_niche_z(A, is_stroma):
    """Binomial z of the number of stromal neighbours (spot-level neighbourhood enrichment)."""
    deg = np.asarray(A.sum(axis=1)).ravel()
    k = np.asarray(A @ is_stroma.astype(float)).ravel()
    p = float(is_stroma.mean())
    sd = np.sqrt(np.maximum(deg * p * (1 - p), 1e-12))
    return (k - deg * p) / sd


def section(sid, P, cfg):
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    A = sp.csr_matrix(A)
    nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
    S = scores_from_nodes(nodes)
    xy = coords_from_nodes(nodes)
    bz = np.load(P.barrier(sid), allow_pickle=True)
    cc = cfg["barrier"]["b_cell"]
    fc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cc)
    reach = np.asarray(fc["reachable"], bool) & np.asarray(bz["reachable"], bool)
    # unreachable spots carry a placeholder (largest finite value); they are excluded here
    f_cell = np.where(reach, fc["b_cell_field"], np.nan)
    f_mab = np.where(reach, np.asarray(bz["b_mab"], float), np.nan)
    lab, stroma_score = compartments(nodes)
    tum, is_str = lab == 0, lab == 1
    # ---- spot-level baselines ----
    Ah = A + sp.identity(A.shape[0], format="csr")
    dens = np.asarray(Ah @ (0.5 * (S["ecm"] + S["caf"]))).ravel() / np.asarray(Ah.sum(axis=1)).ravel()
    dist = signed_distance(xy, tum)
    niche = spot_niche_z(A, is_str)
    base = dict(stromal_density=dens, dist_tumour_boundary=dist, niche_z=niche)
    fields = dict(b_cell_field=f_cell, b_mab_field=f_mab)
    spot = {}
    for fn, fv in fields.items():
        ok = np.isfinite(fv)
        for bn, bv in base.items():
            m = ok & np.isfinite(bv)
            spot[f"{fn}__{bn}"] = float(spearmanr(fv[m], bv[m])[0])
    spot["b_cell_field__b_mab_field"] = float(spearmanr(f_cell[np.isfinite(f_cell) & np.isfinite(f_mab)],
                                                        f_mab[np.isfinite(f_cell) & np.isfinite(f_mab)])[0])
    # where the minimum cut sits on the density scale (is the cut just the densest matrix?)
    cut_nodes = np.asarray(bz["cut_nodes"], int)
    pct = (rankdata(dens) - 1) / max(len(dens) - 1, 1)
    spot["cut_density_pct_median"] = float(np.median(pct[cut_nodes])) if len(cut_nodes) else float("nan")
    spot["cut_share_top_density_decile"] = float(np.mean(pct[cut_nodes] >= 0.9)) if len(cut_nodes) else float("nan")
    spot["top_decile_share_on_cut"] = float(np.mean(np.isin(np.flatnonzero(pct >= 0.9), cut_nodes)))
    # ---- section-level summaries ----
    flow = compute_b_cell(A, S["ecm"], S["caf"], source, sink, **cc)["max_flow"]
    flow_open = compute_b_cell(A, np.zeros(A.shape[0]), np.zeros(A.shape[0]), source, sink, **cc)["max_flow"]
    peri = (~tum) & (np.asarray(A @ tum.astype(float)).ravel() > 0)
    ves_d = cKDTree(xy[tum]).query(xy[vessel], k=1)[0] if len(vessel) and tum.any() else np.array([np.nan])
    sec = dict(
        log_b_rel=float(np.log(flow_open / flow)),
        b_cell_field_core=float(np.nanmedian(f_cell[sink])),
        b_mab_field_core=float(np.nanmedian(f_mab[sink])),
        stromal_density=float(np.mean(0.5 * (S["ecm"] + S["caf"])[peri])) if peri.any() else float("nan"),
        dist_tumour_boundary=float(np.median(ves_d)),
        niche_z=float(contact_enrichment_z(A, tum, is_str, n_perm=N_PERM, seed=SEED)),
        frac_tumour=float(tum.mean()), frac_stroma=float(is_str.mean()), frac_immune=float((lab == 2).mean()),
        n_spots=int(A.shape[0]), n_reachable=int(reach.sum()))
    return dict(spot=spot, section=sec)


def main():
    from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, REPLICATION_ORDER
    cfg = load_config(None)
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    pmap = patient_map(P)
    cohorts = {**{s: "primary" for s in PRIMARY_ORDER}, **{s: "external" for s in EXTERNAL_ORDER},
               **{s: "replication" for s in REPLICATION_ORDER}}
    t0 = time.time()
    per = {}
    for sid, coh in cohorts.items():
        r = section(sid, P, cfg)
        r.update(cohort=coh, patient=pmap.get(sid, sid))
        per[sid] = r
        s = r["spot"]
        print(f"{sid:<14} {coh:<11} rho(Bcell, dens/dist/niche) = "
              f"{s['b_cell_field__stromal_density']:+.2f} {s['b_cell_field__dist_tumour_boundary']:+.2f} "
              f"{s['b_cell_field__niche_z']:+.2f} | rho(Bmab, ...) = {s['b_mab_field__stromal_density']:+.2f} "
              f"{s['b_mab_field__dist_tumour_boundary']:+.2f} {s['b_mab_field__niche_z']:+.2f}", flush=True)

    # spot level: distribution of within-section rho
    def dist_summary(keys):
        out = {}
        for k in keys:
            for coh in ("primary", "external", "replication", "all"):
                v = np.array([per[s]["spot"][k] for s in per if coh == "all" or per[s]["cohort"] == coh], float)
                v = v[np.isfinite(v)]
                if len(v):
                    out.setdefault(k, {})[coh] = dict(n=int(len(v)), median=float(np.median(v)),
                                                      q25=float(np.quantile(v, .25)), q75=float(np.quantile(v, .75)),
                                                      min=float(v.min()), max=float(v.max()),
                                                      median_abs=float(np.median(np.abs(v))))
        return out
    spot_keys = [f"{f}__{b}" for f in FIELDS for b in BASELINES] + ["b_cell_field__b_mab_field",
                                                                   "cut_density_pct_median",
                                                                   "cut_share_top_density_decile",
                                                                   "top_decile_share_on_cut"]
    spot_sum = dist_summary(spot_keys)
    # section level: Spearman across sections (tumour sections only; LN01 is a non-tumour control)
    tumour_sections = [s for s in per if s != "LN01"]
    sec_sum = {}
    for op in ("log_b_rel", "b_cell_field_core", "b_mab_field_core"):
        for b in BASELINES:
            x = np.array([per[s]["section"][op] for s in tumour_sections], float)
            y = np.array([per[s]["section"][b] for s in tumour_sections], float)
            m = np.isfinite(x) & np.isfinite(y)
            rho, p = spearmanr(x[m], y[m])
            # patient level: mean per patient
            pats = sorted({per[s]["patient"] for s in tumour_sections})
            xp = np.array([np.nanmean([per[s]["section"][op] for s in tumour_sections if per[s]["patient"] == q])
                           for q in pats])
            yp = np.array([np.nanmean([per[s]["section"][b] for s in tumour_sections if per[s]["patient"] == q])
                           for q in pats])
            mp = np.isfinite(xp) & np.isfinite(yp)
            rp, pp = spearmanr(xp[mp], yp[mp])
            sec_sum[f"{op}__{b}"] = dict(n_sections=int(m.sum()), rho=float(rho), p=float(p),
                                         n_patients=int(mp.sum()), rho_patient=float(rp), p_patient=float(pp))
    out = dict(per_section=per, spot_level=spot_sum, section_level=sec_sum,
               settings=dict(q_malig=Q_MALIG, n_perm=N_PERM, seed=SEED,
                             spot_density="closed 1-hop graph neighbourhood mean of (ECM + CAF)/2",
                             spot_distance="signed Euclidean distance to the tumour compartment boundary (um)",
                             spot_niche="binomial z of stromal-neighbour count (p = section stromal fraction)",
                             section_density="mean (ECM + CAF)/2 over non-tumour spots adjacent to tumour",
                             section_distance="median distance from vessel spots to the nearest tumour spot (um)",
                             section_niche="tumour-stroma contact z vs 1,000 label permutations (Squidpy-style)",
                             b_rel="max flow with uniform capacity sigma(a) / max flow with matrix capacities"),
               meta=stamp_run(cfg, {"module": "M52-simple-baselines-spearman", "seconds": round(time.time() - t0, 1)}))
    save_json(P.validation("simple_baselines_spearman.json"), out)
    print(json.dumps({k: v.get("all") for k, v in spot_sum.items()}, indent=1))
    print(json.dumps(sec_sum, indent=1))


if __name__ == "__main__":
    main()
