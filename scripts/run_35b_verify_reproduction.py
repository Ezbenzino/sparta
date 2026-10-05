#!/usr/bin/env python
"""
run_35b_verify_reproduction.py -- recompute the locked primary statistics from node tables
=========================================================================================
Inputs : data/interim/{sid}.nodes.npz + {sid}.graph.npz (+ barrier.npz for cross-checks)
         results/validation/spatial_null_check.json (locked values)
Output : results/validation/reproduction_check.json

Re-runs run_27 (field association + 500 graph-spectral surrogates, seed 20261003,
ledger order) using only numpy/scipy/networkx, and compares every per-section value
with the locked JSON.  Also checks the B_mAb field and d_vessel against barrier.npz.
A pass means the archived node tables + graphs are sufficient to regenerate the
paper's primary numbers on any machine.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.barrier import compute_b_cell_field, compute_b_mab, compute_b_meta  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph, load_json, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402
from sparta.spatial_stats import (control_projector, partial_spearman,  # noqa: E402
                                  partial_spearman_many, spectral_basis)

N_PERM, SEED = 500, 20261003


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    locked = load_json(P.validation("spatial_null_check.json"))["per_slide"]
    with open(P.root / "data" / "ledger.csv", encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f) if r.get("status") == "ingested"]
    rng = np.random.default_rng(SEED)
    rows, worst = {}, {}
    for sid in slides:
        if sid not in locked:          # external sections were not in the locked run
            continue
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
                            D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg["barrier"]["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                           **cfg["barrier"]["b_mab"])["b_mab"]
        z = np.load(P.barrier(sid), allow_pickle=True)
        d_bm = float(np.nanmax(np.abs(bm - z["b_mab"])))
        d_dv = float(np.nanmax(np.abs(dv - z["d_vessel_um"])))
        real = partial_spearman(bc, bm, dv)
        w, V = spectral_basis(A)
        fh = np.abs(V.T @ bm)
        n = len(bm)
        H = control_projector(dv)
        null = np.empty(N_PERM)
        # identical draw order to run_27: one rng.choice(size=n) per surrogate
        for k in range(N_PERM):
            s = rng.choice([-1.0, 1.0], size=n)
            null[k] = partial_spearman_many(bc, V @ (s * fh), dv, H)[0]
        p = (np.sum(null >= real) + 1) / (N_PERM + 1)
        L = locked[sid]
        rows[sid] = dict(
            real_rho=real, locked_real_rho=L["real_rho_partial"],
            null_std=float(null.std()), locked_null_std=L["null_std"],
            p=float(p), locked_p=L["empirical_p_one_sided"],
            filled=filled, locked_filled=L.get("filled_keys", []),
            max_abs_diff_bmab_vs_barrier_npz=d_bm, max_abs_diff_dvessel_um=d_dv)
        print(f"{sid:<7} rho {real:+.6f} vs {L['real_rho_partial']:+.6f} | "
              f"sd {null.std():.4f} vs {L['null_std']:.4f} | p {p:.4f} vs {L['empirical_p_one_sided']:.4f} "
              f"| dBmAb {d_bm:.1e} dDv {d_dv:.1e}", flush=True)
    diffs = dict(
        max_abs_rho=max(abs(r["real_rho"] - r["locked_real_rho"]) for r in rows.values()),
        max_abs_null_sd=max(abs(r["null_std"] - r["locked_null_std"]) for r in rows.values()),
        n_p_identical=sum(abs(r["p"] - r["locked_p"]) < 1e-12 for r in rows.values()),
        max_abs_bmab=max(r["max_abs_diff_bmab_vs_barrier_npz"] for r in rows.values()),
        max_abs_dvessel=max(r["max_abs_diff_dvessel_um"] for r in rows.values()),
    )
    diffs["n_sections"] = len(rows)
    save_json(P.validation("reproduction_check.json"),
              dict(per_slide=rows, summary=diffs,
                   meta=stamp_run(cfg, {"module": "M35b-reproduction-check"})))
    print("\nsummary:", diffs)


if __name__ == "__main__":
    main()
