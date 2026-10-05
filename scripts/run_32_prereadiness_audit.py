#!/usr/bin/env python
"""
run_32_prereadiness_audit.py —— 投稿前三个补强的定量输出
=========================================================
1) 条目6: 按指标实际依赖项计算主关联与匹配版 S2 敏感性
2) bulk ICB: 作为探索性分析报告 AUC(NR>R) 与跨队列方向
3) 条目8: 每张切片 source/sink 节点落在最大连通分量外的数量
"""
from __future__ import annotations
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from sparta.io_ import Paths, load_config, load_graph, load_json

DEGRADED_MISSING_AGTARGET = {"CSCC10", "CSCC14", "CSCC15", "CSCC16"}  # 缺 CD274/PDCD1LG2
DEGRADED_MISSING_EFFLUX   = {"MEL01", "CSCC14", "CSCC15", "CSCC16"}


def section_6_decoupling_robustness():
    cfg = load_config(None)
    P = Paths(cfg)
    d = load_json(P.validation("decoupling.json"))
    sn = load_json(P.validation("spatial_null_check.json"))["per_slide"]
    s2 = load_json(P.validation("s2_matched_selection.json"))["per_slide"]
    # B_mAb uses Ag_target. The graph-distance adjustment and S2 do not use
    # efflux, so input availability must follow the metric's dependencies.
    missing_ag = {"CSCC10", "CSCC14", "CSCC15", "CSCC16"}
    all_slides = sorted(set(d) & set(sn))
    ag_ok = [sid for sid in all_slides if sid not in missing_ag]
    out = {}
    for label, ids in [("all_19", all_slides), ("ag_target_available_15", ag_ok)]:
        rhos = [d[sid]["rho_partial"] for sid in ids]
        n_pos = sum(1 for x in rhos if x > 0)
        out[label] = dict(
            n=len(ids),
            median_rho=float(np.median(rhos)),
            n_positive=n_pos,
            n_negative=sum(1 for x in rhos if x < 0),
            range=[float(min(rhos)), float(max(rhos))],
        )
    ag_p = [sn[sid]["empirical_p_one_sided"] for sid in ag_ok]
    ag_q = _bh(ag_p)
    out["spatial_null_ag_target_available_15"] = dict(
        n=len(ag_ok),
        n_empirical_p_lt_0_05=int(sum(p < .05 for p in ag_p)),
        n_bh_q_lt_0_05=int(sum(q < .05 for q in ag_q)),
        note="Same 500 per-section graph-spectral surrogate tests; BH recalculated within the Ag_target-available family.",
    )
    s2_ids = sorted(s2)
    s2_p = [s2[sid]["p_vs_in_cut_matched"] for sid in s2_ids]
    s2_q = _bh(s2_p)
    s2_ratios = [s2[sid]["ratio_vs_in_cut_matched"] for sid in s2_ids]
    out["s2_matched_all_19"] = dict(
        n=len(s2_ids),
        n_sig_raw=int(sum(p < .05 for p in s2_p)),
        n_sig_bh=int(sum(q < .05 for q in s2_q)),
        median_ratio=float(np.median(s2_ratios)),
        note="S2 depends on ECM/CAF and source/sink only; no missing-input exclusion applied.",
    )
    out["input_missingness_by_metric"] = {
        "ag_target_missing": sorted(missing_ag),
        "efflux_missing": sorted(DEGRADED_MISSING_EFFLUX),
    }
    return out


def _bh(p_values):
    p = np.asarray(p_values, float)
    order = np.argsort(p)
    q_sorted = np.minimum.accumulate((p[order] * len(p) /
                                      np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty(len(p), float)
    q[order] = np.clip(q_sorted, 0, 1)
    return q


def section_bulk_icb():
    cfg = load_config(None)
    d = load_json(Paths(cfg).validation("icb_crosscohort_summary.json"))
    # Exploratory only: no time-stamped pre-outcome protocol has been identified.
    # AUC_R_vs_NR as stored = P(score_R > score_NR).
    # AUC(NR>R) = 1 - AUC_R_vs_NR.  >0.5 supports H1.
    rows = []
    for cohort, blk in d.items():
        b = blk["Barrier"]
        auc_rnr = b["AUC_R_vs_NR"]
        auc_nr_over_r = 1.0 - auc_rnr
        rows.append(dict(
            cohort=cohort,
            barrier_R_median=round(b["R_med"], 3),
            barrier_NR_median=round(b["NR_med"], 3),
            AUC_NR_over_R=round(auc_nr_over_r, 3),
            p_value=round(b["p"], 3),
            direction_consistent_with_H1=bool(auc_nr_over_r > 0.5),
        ))
    return dict(status="exploratory; pre-outcome protocol not documented",
                tested_direction="barrier score higher in non-responders than responders",
                cohorts=rows)


def section_8_unreachable_nodes():
    cfg = load_config(None)
    P = Paths(cfg)
    ledger = Path(__file__).resolve().parents[1] / "data" / "ledger.csv"
    import csv
    with open(ledger, encoding="utf-8") as f:
        slides = [r["slide_id"] for r in csv.DictReader(f) if r.get("status") == "ingested"]
    rows = []
    for sid in slides:
        gp = P.graph(sid)
        if not gp.exists():
            continue
        A, D, source, sink, vessel, _ = load_graph(gp)
        n = A.shape[0]
        n_comp, labels = connected_components(csr_matrix(A), directed=False)
        # largest component
        sizes = np.bincount(labels, minlength=n_comp)
        lcc = int(np.argmax(sizes))
        in_lcc = labels == lcc
        n_src_unreachable = int(sum(1 for s in source if not in_lcc[s]))
        n_snk_unreachable = int(sum(1 for s in sink if not in_lcc[s]))
        rows.append(dict(
            slide=sid,
            n_nodes=n,
            n_components=int(n_comp),
            lcc_size=int(sizes[lcc]),
            n_source=len(source),
            n_sink=len(sink),
            source_unreachable=n_src_unreachable,
            sink_unreachable=n_snk_unreachable,
            frac_unreachable=round((n_src_unreachable+n_snk_unreachable)/(len(source)+len(sink)), 3)
                                if (len(source)+len(sink)) else 0.0,
        ))
    return rows


def main():
    cfg = load_config(None)
    out = {
        "section_6_decoupling_robustness": section_6_decoupling_robustness(),
        "section_bulk_icb_exploratory": section_bulk_icb(),
        "section_8_unreachable_nodes": section_8_unreachable_nodes(),
    }
    P = Paths(cfg)
    op = P.validation("prereadiness_audit.json")
    Path(op).parent.mkdir(parents=True, exist_ok=True)
    with open(op, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"\nWrote {op}")


if __name__ == "__main__":
    main()
