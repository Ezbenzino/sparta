#!/usr/bin/env python
"""
run_50_exact_mincut_refresh.py -- recompute every stored minimum cut with exact (integer) capacities
====================================================================================================
Inputs : data/interim/{sid}.graph.npz, {sid}.nodes.npz, {sid}.mincut.json, {sid}.barrier.npz
Outputs: data/interim/{sid}.mincut.json        (cut edges of an exact minimum cut)
         data/interim/{sid}.barrier.npz        (cut_nodes / b_cell / max_flow fields refreshed; all
                                               other arrays copied unchanged)
         results/validation/mincut_exactness.json   old versus new cut, per section

Why (found 2026-10-05)
----------------------
Up to v2.1.0, ``compute_b_cell`` ran networkx's preflow-push on floating-point capacities. The
maximum-flow VALUE was correct (relative error < 1e-12 against an exact integer computation), but
the partition read off the residual graph was not reliable under floating-point round-off: in 12 of
22 sections the returned edge set had a total capacity 0.03-5.4% away from the maximum flow, i.e. it
was not exactly a minimum cut. ``barrier.exact_min_cut`` now solves the same problem with capacities
rounded to integer multiples of 2^-40, for which the cut capacity equals the maximum flow. B_cell is
unchanged; everything that uses the cut *geometry* (S2 gap experiments, domain-boundary comparison,
intervention maps, cut drawings) must be recomputed from the refreshed files.

Usage:
    python scripts/run_50_exact_mincut_refresh.py                  # refresh every stored cut + diagnose
    python scripts/run_50_exact_mincut_refresh.py --diagnose-only  # float-vs-exact gap on current graphs
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def cut_capacity(cut_edges, pairs, cap):
    """Total capacity of an edge set; edges outside ``pairs`` are counted separately."""
    d = {(int(u), int(v)): float(c) for (u, v), c in zip(pairs, cap)}
    tot, missing = 0.0, 0
    for u, v in cut_edges:
        k = (min(int(u), int(v)), max(int(u), int(v)))
        if k in d:
            tot += d[k]
        else:
            missing += 1
    return float(tot), missing


def float_partition_gap(pairs, cap, nodes, source, sink):
    """Gap of the partition that floating-point preflow-push returns on the same capacities
    (the implementation used up to v2.1.0): (cut capacity - max flow) / max flow."""
    import networkx as nx
    G = nx.DiGraph()
    G.add_nodes_from(int(u) for u in nodes)
    for (u, v), c in zip(pairs.tolist(), cap):
        G.add_edge(u, v, capacity=float(c))
        G.add_edge(v, u, capacity=float(c))
    big = float(sum(cap) * 10 + 1e6)
    for i in source:
        G.add_edge("s", int(i), capacity=big)
    for i in sink:
        G.add_edge(int(i), "t", capacity=big)
    f, (R, U) = nx.minimum_cut(G, "s", "t", capacity="capacity")
    d = {(int(u), int(v)): float(c) for (u, v), c in zip(pairs, cap)}
    cut = sum(d[(min(u, v), max(u, v))] for u in R if u != "s" for v in G.successors(u) if v in U and v != "t")
    return float((cut - f) / f), float(f)


def diagnose(P, cfg):
    """Float-versus-exact comparison on the current graphs (independent of the stored files)."""
    import numpy as np
    from scipy.sparse.csgraph import connected_components

    from sparta.barrier import compute_b_cell
    from sparta.io_ import load_graph, load_json, save_json
    from sparta.node_tables import load_nodes, scores_from_nodes
    out = load_json(P.validation("mincut_exactness.json"))
    for sid, row in out["per_section"].items():
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_nodes(load_nodes(P.interim / f"{sid}.nodes.npz"))
        r = compute_b_cell(A, S["ecm"], S["caf"], source, sink, return_capacity=True, **cfg["barrier"]["b_cell"])
        _, lab = connected_components(A, directed=False)
        big = np.bincount(lab).argmax()
        gap, f_float = float_partition_gap(r["edge_pairs"], r["edge_capacity"], np.flatnonzero(lab == big),
                                           source[lab[source] == big], sink[lab[sink] == big])
        row["float_partition_gap"] = gap
        row["float_flow_rel_diff"] = float((f_float - r["max_flow"]) / r["max_flow"])
        print(f"{sid:<14} float-partition gap {gap:+.3%}   float flow rel diff {row['float_flow_rel_diff']:+.1e}")
    g = {k: v["float_partition_gap"] for k, v in out["per_section"].items()}
    out["summary"].update(
        float_partition_n_not_min=int(sum(abs(x) > 1e-6 for x in g.values())),
        float_partition_max_abs_gap=float(max(abs(x) for x in g.values())),
        float_partition_max_abs_gap_section=max(g, key=lambda k: abs(g[k])),
        float_flow_max_abs_rel_diff=float(max(abs(v["float_flow_rel_diff"]) for v in out["per_section"].values())))
    save_json(P.validation("mincut_exactness.json"), out)
    print({k: v for k, v in out["summary"].items() if k.startswith("float")})


def main():
    import numpy as np

    from sparta.barrier import compute_b_cell
    from sparta.io_ import Paths, load_config, load_graph, load_json, save_json, stamp_run
    from sparta.node_tables import load_nodes, scores_from_nodes

    cfg = load_config(ROOT / "configs" / "default.yaml")
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    if "--diagnose-only" in sys.argv:
        return diagnose(P, cfg)
    sids = sorted(p.name.split(".")[0] for p in P.interim.glob("*.nodes.npz"))
    per = {}
    t0 = time.time()
    for sid in sids:
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        S = scores_from_nodes(load_nodes(P.interim / f"{sid}.nodes.npz"))
        new = compute_b_cell(A, S["ecm"], S["caf"], source, sink, return_capacity=True, **cfg["barrier"]["b_cell"])
        pairs, cap = new["edge_pairs"], new["edge_capacity"]
        old_edges = [tuple(e) for e in load_json(P.mincut(sid))["cut_edges"]]
        bz = dict(np.load(P.barrier(sid), allow_pickle=True))
        old_flow = float(np.asarray(bz.get("max_flow", np.nan)).ravel()[0]) if "max_flow" in bz else float("nan")
        # capacities of every edge of the graph (the stored cut may predate the largest-component rule)
        from sparta.barrier import _pair_mean, _sigmoid, edge_pairs
        allp = edge_pairs(A, upper_only=True)
        bc_ = cfg["barrier"]["b_cell"]
        allc = _sigmoid(bc_["a"] - bc_["b_ecm"] * _pair_mean(S["ecm"], allp) - bc_["c_caf"] * _pair_mean(S["caf"], allp))
        old_cap, old_missing = cut_capacity(old_edges, allp, allc) if old_edges else (float("nan"), 0)
        new_cap, _ = cut_capacity(new["cut_edges"], pairs, cap)
        lcc_set = {(int(u), int(v)) for u, v in pairs}
        old_outside_lcc = int(sum((min(u, v), max(u, v)) not in lcc_set for u, v in old_edges))
        so = {frozenset(e) for e in old_edges}
        sn = {frozenset(e) for e in new["cut_edges"]}
        per[sid] = dict(
            n_nodes=int(A.shape[0]), max_flow=new["max_flow"], max_flow_stored=old_flow,
            max_flow_rel_diff=float((new["max_flow"] - old_flow) / old_flow) if old_flow == old_flow else None,
            old_cut_capacity=old_cap, new_cut_capacity=new_cap,
            old_cut_gap=float((old_cap - new["max_flow"]) / new["max_flow"]),
            new_cut_gap=float((new_cap - new["max_flow"]) / new["max_flow"]),
            n_cut_edges_old=len(so), n_cut_edges_new=len(sn),
            jaccard_edges=float(len(so & sn) / len(so | sn)) if (so | sn) else 1.0,
            n_old_edges_outside_lcc=old_outside_lcc, n_old_edges_not_in_graph=old_missing,
            n_cut_nodes_old=int(len(np.asarray(bz.get("cut_nodes", []))) ),
            n_cut_nodes_new=int(len(new["cut_nodes"])),
        )
        # refresh the stored cut (all other arrays of barrier.npz are copied unchanged)
        save_json(P.mincut(sid), {"cut_edges": [[int(u), int(v)] for u, v in new["cut_edges"]]})
        bz["cut_nodes"] = np.asarray(new["cut_nodes"], int)
        bz["b_cell"] = np.asarray(new["b_cell"])
        bz["max_flow"] = np.asarray(new["max_flow"])
        np.savez_compressed(P.barrier(sid), **bz)
        r = per[sid]
        print(f"{sid:<14} flow {r['max_flow']:.6f} (stored rel diff {r['max_flow_rel_diff']:+.1e}) | "
              f"old cut gap {r['old_cut_gap']:+.3%} -> new {r['new_cut_gap']:+.1e} | edges {r['n_cut_edges_old']}"
              f" -> {r['n_cut_edges_new']} (Jaccard {r['jaccard_edges']:.2f})", flush=True)
    gaps = [abs(v["old_cut_gap"]) for v in per.values() if v["old_cut_gap"] == v["old_cut_gap"]]
    out = dict(per_section=per,
               summary=dict(n_sections=len(per),
                            n_old_cut_not_minimum=int(sum(g > 1e-6 for g in gaps)),
                            max_abs_old_cut_gap=float(max(gaps)),
                            max_abs_new_cut_gap=float(max(abs(v["new_cut_gap"]) for v in per.values())),
                            max_abs_flow_rel_diff=float(max(abs(v["max_flow_rel_diff"]) for v in per.values()
                                                            if v["max_flow_rel_diff"] is not None)),
                            sections_stored_flow_predates_lcc_rule=sorted(
                                k for k, v in per.items()
                                if v["max_flow_rel_diff"] is not None and abs(v["max_flow_rel_diff"]) > 1e-9),
                            max_abs_flow_rel_diff_other_sections=float(max(
                                abs(v["max_flow_rel_diff"]) for v in per.values()
                                if v["max_flow_rel_diff"] is not None and abs(v["max_flow_rel_diff"]) <= 1e-9)),
                            median_jaccard_edges=float(np.median([v["jaccard_edges"] for v in per.values()]))),
               semantics=("old = cut stored by run_04 with floating-point preflow-push; new = exact minimum cut "
                          "(barrier.exact_min_cut, capacities in units of 2^-40). gap = (cut capacity - max flow) / "
                          "max flow; a true minimum cut has gap 0. In the sections listed under "
                          "sections_stored_flow_predates_lcc_rule the stored max flow had been computed before the "
                          "largest-component rule (fragments with both sources and sinks added flow); these sections "
                          "have several graph components. B_cell scalars were not used for any reported statistic; "
                          "every statistic reads the field, the recomputed operator or the cut geometry."),
               meta=stamp_run(cfg, {"module": "M50-exact-mincut", "seconds": round(time.time() - t0, 1)}))
    save_json(P.validation("mincut_exactness.json"), out)
    print(json.dumps(out["summary"], indent=1))
    diagnose(P, cfg)


if __name__ == "__main__":
    main()
