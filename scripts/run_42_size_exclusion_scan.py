#!/usr/bin/env python
"""
run_42_size_exclusion_scan.py -- fraction of graph edges with complete steric exclusion over a beta grid
=======================================================================================================
Inputs : data/interim/{sid}.nodes.npz, {sid}.graph.npz
Output : results/validation/size_exclusion_scan.json

For every section and every beta, the edge crosslinking score x_e (mean of the two
endpoint ranks) is mapped to xi = xi0 * exp(-beta * x_e); an edge is completely excluded
for a molecule of radius r when r >= xi (Phi = 0 in Eq. 2).  Because x_e is built from
within-section ranks, the excluded fraction is close to a property of beta alone.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.barrier import edge_pairs  # noqa: E402
from sparta.io_ import Paths, admitted_slides, load_config, load_graph, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402

BETAS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    bm = cfg["barrier"]["b_mab"]
    xi0, r = float(bm["xi0_nm"]), float(bm["r_nm"])
    per = {}
    for sid in admitted_slides(P):
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, *_ = load_graph(P.graph(sid))
        S = scores_from_nodes(nodes)
        pr = edge_pairs(A)
        xe = 0.5 * (S["crosslink"][pr[:, 0]] + S["crosslink"][pr[:, 1]])
        per[sid] = dict(n_edges=int(len(pr)),
                        excluded_fraction=[float(np.mean(r >= xi0 * np.exp(-b * xe))) for b in BETAS])
    prim = [s for s in per if not s.startswith(("BRCA", "LN"))]
    k = BETAS.index(3.0)
    at3 = np.array([per[s]["excluded_fraction"][k] for s in prim])
    summ = dict(beta_default=3.0, threshold_x=float(np.log(xi0 / r) / 3.0),
                primary_min_pct=float(100 * at3.min()), primary_max_pct=float(100 * at3.max()),
                primary_median_pct=float(100 * np.median(at3)),
                median_curve_pct=[float(100 * np.median([per[s]["excluded_fraction"][i] for s in prim]))
                                  for i in range(len(BETAS))])
    save_json(P.validation("size_exclusion_scan.json"),
              dict(betas=BETAS, xi0_nm=xi0, r_nm=r, per_slide=per, summary=summ,
                   meta=stamp_run(cfg, {"module": "M42-size-exclusion-scan"})))
    print(summ)


if __name__ == "__main__":
    main()
