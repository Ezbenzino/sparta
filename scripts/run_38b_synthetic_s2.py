#!/usr/bin/env python
"""
run_38b_synthetic_s2.py -- does the contiguous-gap counterfactual (S2) detect topology when it is known?
=====================================================================================================
Output: results/validation/synthetic_s2_thickness.json

Planted closed capsules of three thicknesses (160, 330, 500 um) on the benchmark lattice
(run_38).  For each capsule S2 removes 20 % of the cut nodes as one contiguous arc
(best of 8 candidate arcs) and compares it with scattered removals of the same number of
cut nodes, either unmatched (one random scattered set per draw, the original design)
or selection-matched (best of 8 scattered sets per draw, the design used for the real
sections).  Ratio > 1 means the contiguous gap leaves a lower model barrier.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import run_38_synthetic_benchmark as B  # noqa: E402
from sparta.counterfactual import s2_ring_breaking  # noqa: E402
from sparta.io_ import Paths, load_config, save_json, stamp_run  # noqa: E402

SEED = 20261007


def main(reps=10, widths=(160.0, 330.0, 500.0)):
    cfg = load_config(None)
    P = Paths(cfg)
    rng = np.random.default_rng(SEED)
    xy, A = B.hex_lattice()
    out = {}
    for w in widths:
        rows = []
        for rep in range(reps):
            t = B.make_tissue(xy, A, "closed", rng, ring_w=w)
            e, c = B.observe(t, rng, 0.15)
            row = dict(K=int(t.K))
            for matched in (False, True):
                o = s2_ring_breaking(A, e, c, t.vessel, t.core, B.CELL, coords=xy, k_frac=(0.20,),
                                     n_rand=40, seed=int(rng.integers(1_000_000_000)), low_q=0.05,
                                     match_selection=matched, n_match_groups=20)
                kk = next(iter(o["per_k"]))
                r = o["per_k"][kk]
                tag = "matched" if matched else "unmatched"
                row[f"ratio_{tag}"] = float(r["ratio_vs_in_cut"])
                row[f"p_{tag}"] = float(r["p_vs_in_cut"])
                row["n_cut_nodes"] = int(o["n_cut_nodes"])
            rows.append(row)
            print(f"w={w:.0f} rep={rep} unmatched={row['ratio_unmatched']:.2f} matched={row['ratio_matched']:.2f}",
                  flush=True)
        out[f"{int(w)}um"] = dict(
            per_rep=rows,
            median_ratio_unmatched=float(np.median([r["ratio_unmatched"] for r in rows])),
            median_ratio_matched=float(np.median([r["ratio_matched"] for r in rows])),
            n_matched_gt1=int(sum(r["ratio_matched"] > 1 for r in rows)),
            n_matched_p05=int(sum(r["p_matched"] < 0.05 for r in rows)),
            n=len(rows))
    save_json(P.validation("synthetic_s2_thickness.json"),
              dict(widths_um=list(widths), results=out,
                   meta=stamp_run(cfg, {"module": "M38b-synthetic-s2", "seed": SEED})))
    for k, v in out.items():
        print(k, {kk: vv for kk, vv in v.items() if kk != "per_rep"})


if __name__ == "__main__":
    main()
