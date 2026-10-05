#!/usr/bin/env python
"""
run_43_nscore_spatial_null.py -- primary field association re-tested with rank-based spectral surrogates
=======================================================================================================
Inputs : data/interim/{sid}.nodes.npz, {sid}.graph.npz
Output : results/validation/spatial_null_nscore.json

Identical to run_27 (500 graph-spectral sign surrogates, one-sided p = (k+1)/501, BH within
the 19 primary sections and separately within the 3 external sections) except that the
surrogate amplitudes are taken from the spectrum of the NORMAL SCORES of the B_mAb field
instead of its raw values.  run_37 shows that this variant stays calibrated for heavily
skewed fields, where the raw-spectrum test becomes anti-conservative; for fields with the
marginal distribution of the real B_mAb both variants are calibrated.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.barrier import compute_b_cell_field, compute_b_mab, compute_b_meta  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph, patient_map, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402
from sparta.spatial_stats import (control_projector, normal_scores, partial_spearman,  # noqa: E402
                                  partial_spearman_many, spectral_basis)

N_PERM, SEED = 500, 20261008


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
    cfg = load_config(None)
    P = Paths(cfg)
    pmap = patient_map(P)
    with open(P.root / "data" / "ledger.csv", encoding="utf-8") as f:
        led = [r for r in csv.DictReader(f) if r.get("status") == "ingested"]
    rng = np.random.default_rng(SEED)
    rows = {}
    for r_ in led:
        sid = r_["slide_id"]
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg["barrier"]["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **cfg["barrier"]["b_mab"])["b_mab"]
        real = partial_spearman(bc, bm, dv)
        w, V = spectral_basis(A)
        H = control_projector(dv)
        fz = np.abs(V.T @ normal_scores(bm))
        n = len(bm)
        null = np.empty(N_PERM)
        for k in range(N_PERM):
            null[k] = partial_spearman_many(bc, V @ (rng.choice([-1.0, 1.0], size=n) * fz), dv, H)[0]
        p = (np.sum(null >= real) + 1.0) / (N_PERM + 1.0)
        rows[sid] = dict(patient=pmap.get(sid, sid), arm="external" if r_["cancer_type"] == "other" else "primary",
                         n_nodes=int(n), real_rho_partial=float(real), null_mean=float(null.mean()),
                         null_std=float(null.std()), empirical_p_one_sided=float(p), filled_keys=filled)
        print(f"{sid:<7} rho={real:+.3f} null={null.mean():+.3f}±{null.std():.3f} p={p:.4f}", flush=True)
    for arm in ("primary", "external"):
        ids = [s for s in rows if rows[s]["arm"] == arm]
        q = bh([rows[s]["empirical_p_one_sided"] for s in ids])
        for s, qq in zip(ids, q):
            rows[s]["bh_padj"] = float(qq)
    prim = [s for s in rows if rows[s]["arm"] == "primary"]
    summ = dict(n_primary=len(prim),
                n_p05=int(sum(rows[s]["empirical_p_one_sided"] < 0.05 for s in prim)),
                n_q05=int(sum(rows[s]["bh_padj"] < 0.05 for s in prim)),
                median_null_sd=float(np.median([rows[s]["null_std"] for s in prim])),
                external_q=[rows[s]["bh_padj"] for s in rows if rows[s]["arm"] == "external"])
    save_json(P.validation("spatial_null_nscore.json"),
              dict(per_slide=rows, summary=summ, n_perm=N_PERM, seed=SEED,
                   meta=stamp_run(cfg, {"module": "M43-nscore-spatial-null"})))
    print(summ)


if __name__ == "__main__":
    main()
