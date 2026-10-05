#!/usr/bin/env python
"""Input-specific and disconnected-component sensitivity audit for CMPB.

The S2 counterfactual only uses ECM/CAF scores and the source-sink graph, so it
is reported on all 19 primary sections.  The B_cell/B_mAb field association
depends on Ag_target and is summarized separately for sections with that input.
The stored graph-spectral surrogate draws are section-specific; restricting the
test family therefore requires a new BH correction, not new surrogate draws.

For the spatial-field association, compare the current all-node statistic with
the same statistic restricted to each graph's largest connected component.
This measures the contribution of finite-value imputation for disconnected
nodes without changing any model inputs on the retained component.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from sparta.io_ import Paths, load_config, load_graph, load_json, save_json
from sparta.barrier import (compute_b_cell_field, compute_b_meta, compute_b_mab,
                            scores_from_adata)
from run_27_spatial_null import partial_spearman

AG_TARGET_MISSING = {"CSCC10", "CSCC14", "CSCC15", "CSCC16"}
EFFLUX_MISSING = {"MEL01", "CSCC14", "CSCC15", "CSCC16"}


def bh(pvalues):
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    q_sorted = np.minimum.accumulate((p[order] * len(p) /
                                      np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty(len(p), float)
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def stats(rows, value_key, p_key=None):
    vals = np.asarray([r[value_key] for r in rows], float)
    result = {
        "n": int(len(vals)),
        "median": float(np.median(vals)),
        "range": [float(np.min(vals)), float(np.max(vals))],
        "n_positive": int(np.sum(vals > 0)),
        "n_negative": int(np.sum(vals < 0)),
    }
    if p_key:
        qs = bh([r[p_key] for r in rows])
        result["n_bh_q_lt_0_05"] = int(np.sum(qs < .05))
        result["bh_q_values"] = [float(x) for x in qs]
    return result


def main():
    cfg = load_config(None)
    P = Paths(cfg)
    spatial = load_json(P.validation("spatial_null_check.json"))["per_slide"]
    matched = load_json(P.validation("s2_matched_selection.json"))["per_slide"]
    slides = sorted(set(spatial) & set(matched))
    rows = {}

    for sid in slides:
        adata = sc.read_h5ad(P.scored(sid))
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled = []
        S = scores_from_adata(adata, missing_out=filled)
        cfg_b = cfg["barrier"]
        dv = compute_b_meta(
            A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
            D=D, **cfg_b["b_meta"])["d_vessel_um"]
        bc = compute_b_cell_field(
            A, S["ecm"], S["caf"], source, **cfg_b["b_cell"])["b_cell_field"]
        bm = compute_b_mab(
            A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
            **cfg_b["b_mab"])["b_mab"]

        _, labels = connected_components(csr_matrix(A), directed=False)
        sizes = np.bincount(labels)
        lcc = labels == int(np.argmax(sizes))
        rho_all = partial_spearman(bc, bm, dv)
        rho_lcc = partial_spearman(bc[lcc], bm[lcc], dv[lcc])
        rows[sid] = {
            "patient": spatial[sid].get("patient"),
            "n_nodes": int(A.shape[0]),
            "lcc_nodes": int(lcc.sum()),
            "non_lcc_nodes": int((~lcc).sum()),
            "source_outside_lcc": int(np.sum(~lcc[np.asarray(source, int)])),
            "sink_outside_lcc": int(np.sum(~lcc[np.asarray(sink, int)])),
            "filled_keys": filled,
            "rho_all_nodes": float(rho_all),
            "rho_lcc_only": float(rho_lcc),
            "delta_lcc_minus_all": float(rho_lcc - rho_all),
            "spatial_null_p": float(spatial[sid]["empirical_p_one_sided"]),
            "s2_ratio": float(matched[sid]["ratio_vs_in_cut_matched"]),
            "s2_p": float(matched[sid]["p_vs_in_cut_matched"]),
        }

    all_rows = [rows[s] for s in slides]
    ag_ok_ids = [s for s in slides if s not in AG_TARGET_MISSING]
    ag_ok = [rows[s] for s in ag_ok_ids]
    s2_rows = [rows[s] for s in slides]
    spatial_q = bh([rows[s]["spatial_null_p"] for s in ag_ok_ids])
    s2_q = bh([r["s2_p"] for r in s2_rows])
    # Explicitly mark S2 q-values for reproducible manuscript rendering.
    for row, q in zip(s2_rows, s2_q):
        row["s2_bh_q"] = float(q)

    out = {
        "definitions": {
            "ag_target_missing": sorted(AG_TARGET_MISSING),
            "efflux_missing": sorted(EFFLUX_MISSING),
            "association_complete_case": "Exclude only sections missing Ag_target; efflux is not part of B_cell/B_mAb fields or graph-distance adjustment.",
            "spatial_null": "Existing 500 graph-spectral sign-randomized surrogates per section; the 15-section subset uses the same section-level empirical p-values and recomputes BH within that family.",
            "s2": "Selection-matched ratio from s2_matched_selection.json; S2 uses ECM/CAF and source/sink only, so all 19 sections are retained.",
            "lcc_sensitivity": "Recalculate partial Spearman on the largest connected component using the unchanged full-graph field values restricted to that component; component-local shortest-path and diffusion solves are independent of other components.",
        },
        "per_slide": rows,
        "summary": {
            "association_all_19": stats(all_rows, "rho_all_nodes"),
            "association_ag_target_available_15": stats(ag_ok, "rho_all_nodes"),
            "spatial_null_ag_target_available_15": {
                "n": len(ag_ok_ids),
                "n_positive": int(sum(rows[s]["rho_all_nodes"] > 0 for s in ag_ok_ids)),
                "n_empirical_p_lt_0_05": int(sum(rows[s]["spatial_null_p"] < .05 for s in ag_ok_ids)),
                "n_bh_q_lt_0_05": int(np.sum(spatial_q < .05)),
                "median_rho": float(np.median([rows[s]["rho_all_nodes"] for s in ag_ok_ids])),
                "bh_q_values_by_slide": {s: float(q) for s, q in zip(ag_ok_ids, spatial_q)},
            },
            "s2_matched_all_19": {
                "n": len(s2_rows),
                "median_ratio": float(np.median([r["s2_ratio"] for r in s2_rows])),
                "ratio_range": [float(min(r["s2_ratio"] for r in s2_rows)),
                                float(max(r["s2_ratio"] for r in s2_rows))],
                "n_ratio_gt_1": int(sum(r["s2_ratio"] > 1 for r in s2_rows)),
                "n_raw_p_lt_0_05": int(sum(r["s2_p"] < .05 for r in s2_rows)),
                "n_bh_q_lt_0_05": int(np.sum(s2_q < .05)),
            },
            "largest_component_association": {
                "all_19": stats(all_rows, "rho_lcc_only"),
                "ag_target_available_15": stats(ag_ok, "rho_lcc_only"),
                "n_sign_changed_all_19": int(sum(np.sign(r["rho_all_nodes"]) !=
                                                    np.sign(r["rho_lcc_only"])
                                                    for r in all_rows)),
                "median_abs_rho_change_all_19": float(np.median(
                    [abs(r["delta_lcc_minus_all"]) for r in all_rows])),
                "max_abs_rho_change_slide": max(slides,
                    key=lambda s: abs(rows[s]["delta_lcc_minus_all"])),
                "mel02": {k: rows["MEL02"][k] for k in (
                    "rho_all_nodes", "rho_lcc_only", "delta_lcc_minus_all",
                    "non_lcc_nodes", "source_outside_lcc", "sink_outside_lcc")},
            },
        },
    }
    save_json(P.validation("metric_connectivity_sensitivity.json"), out)
    print(json.dumps(out["summary"], indent=2, ensure_ascii=False))
    print(f"\nWrote {P.validation('metric_connectivity_sensitivity.json')}")


if __name__ == "__main__":
    main()
