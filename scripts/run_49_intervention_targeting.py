#!/usr/bin/env python
"""
run_49_intervention_targeting.py -- how concentrated is each barrier, and does topology change the target?
=========================================================================================================
Inputs (read-only):
    results/intervention/{sid}.intervention.npz   per-spot single-ablation effects (run_44)
    data/interim/{sid}.graph.npz, {sid}.nodes.npz  graphs and model inputs
Output:
    results/validation/intervention_targeting.json

What this adds to run_44 (and what it does not claim)
----------------------------------------------------
run_44 ablates every spot once. Two of its summaries are consequences of max-flow/min-cut
duality rather than findings: raising capacities off a minimum cut cannot raise the max flow, so
single-spot effects on B_cell are zero away from the cut, and the top-ranked spots lie on the cut.
The informative questions are joint ones, asked here with the same in-model ablation
(ecm/caf/crosslink of the chosen spots set to the section's 5th percentile; absorption untouched):

1. Concentration: how many top-ranked spots must be ablated *together* to recover half of a full
   breach (B_cell: all cut nodes ablated; B_mAb: every spot ablated)? -> k50 (as a fraction of n)
2. Targeting: for the same budget k, what fraction do other selection rules achieve?
     density      the k spots with the highest local matrix score (no graph, no transport)
     random cut   k random spots of the minimum cut (B_cell only; 20 draws)
     random       k random spots (20 draws)
   and how much of the B_mAb limit does the *B_cell* ranking reach (cross-modality transfer)?
Everything here is a model-defined counterfactual; nothing is a measurement of permeability.

Usage:
    python scripts/run_49_intervention_targeting.py            # primary 19 + external 3 + replication 8
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "is_figures"))

N_DRAW = 20
SEED = 20261005
K_FIXED = 20


def main():
    import numpy as np

    from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, REPLICATION_ORDER
    from sparta.barrier import compute_b_cell, compute_b_mab
    from sparta.io_ import Paths, load_config, load_graph, patient_map, save_json, stamp_run
    from sparta.node_tables import load_nodes, scores_from_nodes

    cfg = load_config(ROOT / "configs" / "default.yaml")
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    pmap = patient_map(P)
    cc, cm = cfg["barrier"]["b_cell"], cfg["barrier"]["b_mab"]
    rng = np.random.default_rng(SEED)
    per = {}
    t0 = time.time()
    for sid in PRIMARY_ORDER + EXTERNAL_ORDER + REPLICATION_ORDER:
        z = np.load(P.results / "intervention" / f"{sid}.intervention.npz", allow_pickle=True)
        dcell, dmab, cut = z["delta_b_cell"], z["delta_b_mab"], z["cut_mask"].astype(bool)
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_nodes(load_nodes(P.interim / f"{sid}.nodes.npz"))
        E, F, X, G = S["ecm"], S["caf"], S["crosslink"], S["ag_target"]
        n = A.shape[0]
        qE, qF, qX = (float(np.quantile(v, 0.05)) for v in (E, F, X))
        b0 = compute_b_cell(A, E, F, source, sink, **cc)["b_cell"]
        base = compute_b_mab(A, E, X, G, vessel, **cm)
        reach = np.asarray(base["reachable"], bool)
        m0 = float(np.nanmean(base["b_mab"][reach]))

        def drop_cell(idx):
            e2, f2 = E.copy(), F.copy()
            e2[idx], f2[idx] = qE, qF
            return float(b0 - compute_b_cell(A, e2, f2, source, sink, **cc)["b_cell"])

        def drop_mab(idx):
            e2, x2 = E.copy(), X.copy()
            e2[idx], x2[idx] = qE, qX
            return float(m0 - np.nanmean(compute_b_mab(A, e2, x2, G, vessel, **cm)["b_mab"][reach]))

        cut_idx = np.flatnonzero(cut)
        anchor_c = drop_cell(cut_idx)
        anchor_m = drop_mab(np.arange(n))
        k5 = max(1, int(round(0.05 * n)))
        order_c = np.argsort(-dcell, kind="stable")
        order_m = np.argsort(-dmab, kind="stable")
        dens_c = np.argsort(-((8 * E + 4 * F) / 12.0), kind="stable")
        dens_m = np.argsort(-(0.5 * (E + X)), kind="stable")
        row = dict(patient=pmap.get(sid, sid), n=int(n), n_cut_nodes=int(cut.sum()),
                   cut_frac=float(cut.mean()), anchor_cell=anchor_c, anchor_mab=anchor_m,
                   frac_zero_cell=float(np.mean(dcell == 0)),
                   top5_on_cut=float(np.mean(cut[order_c[:k5]])))
        for k_name, k in (("k20", min(K_FIXED, n)), ("k5pct", k5)):
            r = {}
            r["cell_sparta"] = drop_cell(order_c[:k]) / anchor_c
            r["cell_density"] = drop_cell(dens_c[:k]) / anchor_c
            r["cell_random_cut"] = float(np.mean([drop_cell(rng.choice(cut_idx, size=min(k, len(cut_idx)),
                                                                       replace=False))
                                                  for _ in range(N_DRAW)])) / anchor_c
            r["cell_random"] = float(np.mean([drop_cell(rng.choice(n, size=k, replace=False))
                                              for _ in range(N_DRAW)])) / anchor_c
            r["mab_sparta"] = drop_mab(order_m[:k]) / anchor_m
            r["mab_density"] = drop_mab(dens_m[:k]) / anchor_m
            r["mab_random"] = float(np.mean([drop_mab(rng.choice(n, size=k, replace=False))
                                             for _ in range(N_DRAW)])) / anchor_m
            r["mab_with_cell_ranking"] = drop_mab(order_c[:k]) / anchor_m
            r["cell_with_mab_ranking"] = drop_cell(order_m[:k]) / anchor_c
            row[k_name] = r
        # concentration: smallest k reaching half of the anchor with the SPARTA ranking
        grid = sorted({1, 2, 5, 10, 20, 50, 100, 200, k5, int(round(0.1 * n)), int(round(0.2 * n)),
                       int(round(0.4 * n)), n})
        grid = [g for g in grid if 1 <= g <= n]
        curve_c = [(g, drop_cell(order_c[:g]) / anchor_c) for g in grid]
        curve_m = [(g, drop_mab(order_m[:g]) / anchor_m) for g in grid]
        k50c = next((g for g, f in curve_c if f >= 0.5), None)
        k50m = next((g for g, f in curve_m if f >= 0.5), None)
        row.update(curve_cell=curve_c, curve_mab=curve_m,
                   k50_cell=k50c, k50_cell_frac=(k50c / n) if k50c else None,
                   k50_mab=k50m, k50_mab_frac=(k50m / n) if k50m else None,
                   top5pct_jaccard=float(len(set(order_c[:k5]) & set(order_m[:k5]))
                                         / len(set(order_c[:k5]) | set(order_m[:k5]))))
        per[sid] = row
        print(f"{sid:<7} n={n:<5} cut={cut.mean():.2f} | k20 cell: sparta {row['k20']['cell_sparta']:.2f} "
              f"dens {row['k20']['cell_density']:.2f} rcut {row['k20']['cell_random_cut']:.2f} | "
              f"mab: sparta {row['k20']['mab_sparta']:.2f} dens {row['k20']['mab_density']:.2f} "
              f"rand {row['k20']['mab_random']:.3f} | k50 cell {row['k50_cell_frac']} mab {row['k50_mab_frac']} "
              f"({time.time() - t0:.0f}s)", flush=True)

    def med(sids, *path):
        vals = []
        for s in sids:
            v = per[s]
            for p in path:
                v = v[p]
            if v is not None:
                vals.append(v)
        return float(np.median(vals)) if vals else None

    def nfrac(sids, a, b, k="k20"):
        return int(sum(per[s][k][a] > per[s][k][b] for s in sids))

    summ = {}
    for name, sids in (("primary", PRIMARY_ORDER), ("external", EXTERNAL_ORDER), ("replication", REPLICATION_ORDER)):
        d = dict(n_sections=len(sids), median_cut_frac=med(sids, "cut_frac"),
                 median_frac_zero_cell=med(sids, "frac_zero_cell"), median_top5_on_cut=med(sids, "top5_on_cut"),
                 median_k50_cell_frac=med(sids, "k50_cell_frac"), median_k50_mab_frac=med(sids, "k50_mab_frac"),
                 n_k50_mab_reached=int(sum(per[s]["k50_mab"] is not None for s in sids)),
                 median_jaccard=med(sids, "top5pct_jaccard"))
        for k in ("k20", "k5pct"):
            for key in ("cell_sparta", "cell_density", "cell_random_cut", "cell_random", "mab_sparta",
                        "mab_density", "mab_random", "mab_with_cell_ranking", "cell_with_mab_ranking"):
                d[f"{k}_{key}"] = med(sids, k, key)
            d[f"{k}_n_sparta_gt_density_cell"] = nfrac(sids, "cell_sparta", "cell_density", k)
            d[f"{k}_n_sparta_gt_random_cut_cell"] = nfrac(sids, "cell_sparta", "cell_random_cut", k)
            d[f"{k}_n_sparta_gt_density_mab"] = nfrac(sids, "mab_sparta", "mab_density", k)
        summ[name] = d
    out = dict(per_section=per, summary=summ,
               settings=dict(k_fixed=K_FIXED, k_frac=0.05, n_draw=N_DRAW, seed=SEED, low_q=0.05),
               semantics=("In-model joint ablation (ecm/caf/crosslink of the selected spots set to the section's "
                          "5th percentile; absorption unchanged). Fractions are of a full minimum-cut breach "
                          "(B_cell) and of the all-spots ablation limit (B_mAb). 'density' ranks spots by their "
                          "local matrix score without any graph; 'random cut' draws spots from the minimum cut. "
                          "Zero single-spot B_cell effects away from the cut follow from max-flow/min-cut duality."),
               meta=stamp_run(cfg, {"module": "M49-intervention-targeting", "seconds": round(time.time() - t0, 1)}))
    save_json(P.validation("intervention_targeting.json"), out)
    for k, v in summ.items():
        print(k, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()})


if __name__ == "__main__":
    main()
