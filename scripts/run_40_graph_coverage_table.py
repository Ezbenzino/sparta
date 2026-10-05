#!/usr/bin/env python
"""
run_40_graph_coverage_table.py -- per-section graph, compartment and disconnection table (Online Resource)
========================================================================================================
Inputs : data/interim/{sid}.graph.npz, {sid}.barrier.npz, {sid}.nodes.npz, {sid}.mincut.json
         results/validation/metric_connectivity_sensitivity.json (largest-component re-estimate)
         data/ledger.csv
Outputs: results/validation/graph_coverage_table.json and .csv

One row per section: platform, patient, spots retained, edges, connected components,
largest-component share, vessel / source / sink counts, source and sink nodes outside
the largest component (and their share of the source+sink set), nodes unreachable from
every source (B_cell field) and from every vessel (B_mAb field), min-cut size, inputs
filled with the neutral value, and the change in the partial correlation when the
association is restricted to the largest component.  Codex readiness item 8.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

from sparta.barrier import compute_b_cell_field  # noqa: E402
from sparta.io_ import Paths, admitted_slides, load_config, load_graph, load_json, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes, scores_from_nodes  # noqa: E402


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    with open(P.root / "data" / "ledger.csv", encoding="utf-8") as f:
        led = {r["slide_id"]: r for r in csv.DictReader(f)}
    mc = load_json(P.validation("metric_connectivity_sensitivity.json"))["per_slide"]
    rows = []
    for sid in admitted_slides(P):
        A, D, source, sink, vessel, gmeta = load_graph(P.graph(sid))
        z = np.load(P.barrier(sid), allow_pickle=True)
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
        n = A.shape[0]
        ncomp, lab = connected_components(A, directed=False)
        sizes = np.bincount(lab)
        lcc = int(np.argmax(sizes))
        src_out = int(np.sum(lab[source] != lcc))
        snk_out = int(np.sum(lab[sink] != lcc))
        reach_src = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cfg["barrier"]["b_cell"])["reachable"]
        cut = json.load(open(P.mincut(sid), encoding="utf-8"))["cut_edges"]
        r = led[sid]
        row = dict(section=sid, patient=r["patient"], cancer_type=r["cancer_type"], platform=r["platform"],
                   accession=r["accession"], median_umi=float(r["median_umi"]) if r.get("median_umi") else None,
                   n_spots=int(n), n_edges=int(A.nnz // 2), n_components=int(ncomp),
                   lcc_share=float(sizes[lcc] / n),
                   n_vessel=int(len(vessel)), n_source=int(len(source)), n_sink=int(len(sink)),
                   source_outside_lcc=src_out, sink_outside_lcc=snk_out,
                   stranded_share_of_source_sink=float((src_out + snk_out) / (len(source) + len(sink))),
                   unreachable_from_source=int((~reach_src).sum()),
                   unreachable_from_vessel=int((~z["reachable"]).sum()),
                   n_cut_edges=int(len(cut)), n_cut_nodes=int(len(z["cut_nodes"])),
                   filled_inputs=";".join(filled))
        m = mc.get(sid)
        if m:
            row.update(rho_all_nodes=m["rho_all_nodes"], rho_lcc_only=m["rho_lcc_only"],
                       delta_lcc_minus_all=m["delta_lcc_minus_all"])
        rows.append(row)
        print(f"{sid:<7} n={n:>5} comp={ncomp:>3} lcc={row['lcc_share']:.3f} src_out={src_out:>3} "
              f"snk_out={snk_out:>3} ({100*row['stranded_share_of_source_sink']:.1f}%) cut={len(cut)}")
    prim = [r for r in rows if r["cancer_type"] != "other"]
    summ = dict(n_primary=len(prim),
                spots_range=[min(r["n_spots"] for r in prim), max(r["n_spots"] for r in prim)],
                components_range=[min(r["n_components"] for r in prim), max(r["n_components"] for r in prim)],
                lcc_share_range=[min(r["lcc_share"] for r in prim), max(r["lcc_share"] for r in prim)],
                n_sections_with_stranded_source_or_sink=int(sum((r["source_outside_lcc"] + r["sink_outside_lcc"]) > 0 for r in prim)),
                max_stranded_share=max(r["stranded_share_of_source_sink"] for r in prim),
                max_stranded_section=max(prim, key=lambda r: r["stranded_share_of_source_sink"])["section"])
    save_json(P.validation("graph_coverage_table.json"),
              dict(rows=rows, summary=summ, meta=stamp_run(cfg, {"module": "M40-graph-coverage"})))
    with open(P.validation("graph_coverage_table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) + [k for k in ("rho_all_nodes", "rho_lcc_only", "delta_lcc_minus_all") if k not in rows[0]])
        w.writeheader()
        w.writerows(rows)
    print(summ)


if __name__ == "__main__":
    main()
