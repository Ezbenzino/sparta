#!/usr/bin/env python
"""
run_48_replication_cohort.py -- independent melanoma replication cohort (Thrane et al. 2018)
============================================================================================
Inputs (read-only):
    data/external/thrane2018/ST-Melanoma-Datasets_1.zip   eight first-generation ST count tables
        (Thrane et al. 2018 Cancer Res 78:5970; lymph-node metastases, four stage III patients,
         two consecutive sections each; https://www.spatialresearch.org)
    configs/cscc_legacy_st.yaml                             first-generation ST settings (unchanged)
    docs/replication_melanoma_protocol.md                   pre-specified plan (SHA-256 recorded)
    data/raw/CSCC05|08|12 (optional, --validate)           raw tables of three primary sections
Outputs:
    data/interim/{sid}.admission.json / .graph.npz / .barrier.npz / .mincut.json /
        .barrier_meta.json / .nodes.npz                     same formats as the primary cohort
    data/ledger.csv                                         rows MEL_THR*: status "replication"
    results/validation/replication_melanoma.json            association, structural nulls,
                                                            patient-level inference
    results/validation/lite_scoring_equivalence.json        (--validate) Scanpy-free M1/M2 check

Why
---
The primary cohort has two melanoma patients among its four melanoma sections. This
cohort adds four patients, processed with the locked pipeline and parameters, as a
separate family that is never pooled into the primary numbers.
Scanpy is not needed: ``sparta/lite.py`` re-implements M1/M2 and ``--validate`` shows that it
reproduces the archived scores of primary first-generation sections.

Usage:
    python scripts/run_48_replication_cohort.py --validate      # equivalence check only
    python scripts/run_48_replication_cohort.py                 # build + statistics (~10 min)
    python scripts/run_48_replication_cohort.py --stats-only    # statistics from built files
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SECTIONS = {   # sid -> (member in the zip, patient, replicate)
    "MEL_THR1_rep1": ("ST_mel1_rep1_counts.tsv", "MEL_THR1", "rep1"),
    "MEL_THR1_rep2": ("ST_mel1_rep2_counts.tsv", "MEL_THR1", "rep2"),
    "MEL_THR2_rep1": ("ST_mel2_rep1_counts.tsv", "MEL_THR2", "rep1"),
    "MEL_THR2_rep2": ("ST_mel2_rep2_counts.tsv", "MEL_THR2", "rep2"),
    "MEL_THR3_rep1": ("ST_mel3_rep1_counts.tsv", "MEL_THR3", "rep1"),
    "MEL_THR3_rep2": ("ST_mel3_rep2_counts.tsv", "MEL_THR3", "rep2"),
    "MEL_THR4_rep1": ("ST_mel4_rep1_counts.tsv", "MEL_THR4", "rep1"),
    "MEL_THR4_rep2": ("ST_mel4_rep2_counts.tsv", "MEL_THR4", "rep2"),
}
COHORT = "mel_thrane2018"
VALIDATE = {  # primary first-generation sections used for the Scanpy-free equivalence check
    "CSCC05": ("GSM4284316_P2_ST_rep1_stdata.tsv.gz", "GSM4284316_spot_data-selection-P2_ST_rep1.tsv.gz"),
    "CSCC08": ("GSM4284319_P5_ST_rep1_stdata.tsv.gz", "GSM4284319_spot_data-selection-P5_ST_rep1.tsv.gz"),
    "CSCC12": ("GSM4284323_P9_ST_rep2_stdata.tsv.gz", "GSM4284323_spot_data-selection-P9_ST_rep2.tsv.gz"),
}
N_SUR = 500
N_NULL = 200
SEED = 20261005
ENDOTHELIAL = ["PECAM1", "VWF", "CDH5", "CLDN5", "ENG", "EGFL7"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate", action="store_true", help="only run the Scanpy-free equivalence check")
    ap.add_argument("--stats-only", action="store_true", help="skip building; recompute statistics")
    ap.add_argument("--raw-dir", default=str(ROOT / "data" / "raw"),
                    help="folder with CSCC05/08/12 raw tables (for --validate)")
    return ap.parse_args()


def _hypoxia(cfg):
    gp = ROOT / cfg["signatures"]["hypoxia_gmt"]
    with open(gp, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if "HYPOXIA" in parts[0].upper():
                return parts[2:], str(cfg["signatures"]["hypoxia_gmt"])
    raise SystemExit(f"no HYPOXIA row in {gp}")


# --------------------------------------------------------------------------
def validate(cfg, raw_dir: Path):
    import numpy as np
    from scipy.stats import spearmanr

    from sparta import lite
    from sparta.io_ import Paths, save_json, stamp_run
    from sparta.loaders import read_table_matrix
    from sparta.node_tables import load_nodes

    P = Paths(cfg)
    hyp, _ = _hypoxia(cfg)
    out = {}
    for sid, (cnt, sel) in VALIDATE.items():
        cf, sf = raw_dir / sid / cnt, raw_dir / sid / sel
        if not cf.exists():
            print(f"{sid}: raw table not found ({cf}); skipped")
            continue
        vals, obs, var, _ = read_table_matrix({"files": {"counts": str(cf), "coords": str(sf)}})
        var = lite.make_unique([str(v) for v in var])
        C, ks, _, varq = lite.qc_filter(vals, var, cfg["qc"]["min_counts"], cfg["qc"]["min_cells_per_gene"])
        X = lite.normalize_log1p(C)
        scores, _ = lite.score_all(X, varq, "cscc", hypoxia_genes=hyp,
                                   min_genes=cfg["signatures"]["min_genes_matched"],
                                   ctrl_size=cfg["signatures"]["ctrl_size"], seed=cfg["seed"])
        rn = lite.rank_normalized(scores)
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        same_obs = list(map(str, nodes["obs_names"])) == [o for o, k in zip(obs, ks) if k]
        per = {}
        for k, v in scores.items():
            key = f"obs__{k}"
            if key in nodes:
                per[k] = dict(max_abs_diff_raw=float(np.max(np.abs(v - nodes[key]))),
                              max_abs_diff_rank=float(np.max(np.abs(rn[f"{k}_n"] - nodes[f"obs__{k}_n"])))
                              if f"obs__{k}_n" in nodes and f"{k}_n" in rn else None,
                              spearman=float(spearmanr(v, nodes[key])[0]))
        out[sid] = dict(n_spots=int(C.shape[0]), n_genes=int(C.shape[1]),
                        archived=dict(n_obs=int(nodes["meta"]["n_obs"]), n_vars=int(nodes["meta"]["n_vars"])),
                        same_spot_order=bool(same_obs), signatures=per)
        mx = max(d["max_abs_diff_raw"] for d in per.values())
        mr = max((d["max_abs_diff_rank"] or 0.0) for d in per.values())
        print(f"{sid}: {C.shape[0]} spots x {C.shape[1]} genes; same order {same_obs}; "
              f"max |Δscore| {mx:.2e}; max |Δrank| {mr:.2e}")
    res = dict(sections=out,
               max_abs_diff_raw=max(d["max_abs_diff_raw"] for s in out.values() for d in s["signatures"].values()),
               max_abs_diff_rank=max((d["max_abs_diff_rank"] or 0.0) for s in out.values()
                                     for d in s["signatures"].values()),
               semantics="sparta/lite.py (no Scanpy) re-scores primary first-generation sections from their raw "
                         "tables; differences are float32 rounding of the normalised matrix.",
               meta=stamp_run(cfg, {"module": "M48-lite-equivalence"}))
    save_json(P.validation("lite_scoring_equivalence.json"), res)
    return res


# --------------------------------------------------------------------------
def build_section(sid, cfg, zf, ledger_rows):
    import numpy as np
    import scipy.sparse as sp
    from scipy.spatial import cKDTree

    from sparta import lite
    from sparta.barrier import compute_b_cell, compute_b_mab, compute_b_meta, save_barriers
    from sparta.graph import build_graph_radius, define_source_sink, graph_report
    from sparta.io_ import Paths, save_graph, save_json, stamp_run
    from sparta.loaders import parse_spot_coords

    P = Paths(cfg)
    member, patient, rep = SECTIONS[sid]
    import pandas as pd
    df = pd.read_csv(io.BytesIO(zf.read(member)), sep="\t", index_col=0)
    # genes in rows, spots ("AxB") in columns -> transpose (same rule as loaders._decide_orientation)
    if parse_spot_coords(df.columns[:5]) is not None:
        df = df.T
    obs = [str(x) for x in df.index]
    var_raw = [str(x) for x in df.columns]
    var, note = lite.clean_gene_names(var_raw)
    var = lite.make_unique(var)
    counts = df.values.astype(np.float32)
    xy_arr = parse_spot_coords(obs)

    # ---- admission (C1-C7), cohort override from the config ----
    adm = dict(cfg["admission"])
    ov = cfg["admission_overrides"][COHORT]
    applied = [k for k in ("min_spots_visium", "min_spots_legacy_st", "min_median_umi", "min_endothelial_spots")
               if k in ov]
    for k in applied:
        adm[k] = ov[k]
    tot = counts.sum(axis=1)
    pres = [g for g in ENDOTHELIAL if g in var]
    n_endo = int((counts[:, [var.index(g) for g in pres]].sum(axis=1) > 0).sum()) if pres else 0
    checks = {
        "C1": (True, f"{counts.shape[0]} spots x {counts.shape[1]} genes"),
        "C2": (xy_arr is not None, "coordinates in spot names"),
        "C3": (False, "no paired H&E in the public count release; waived at cohort level (see config)"),
        "C4": (counts.shape[0] >= adm["min_spots_legacy_st"],
               f"{counts.shape[0]} spots (threshold {adm['min_spots_legacy_st']}, 1,007-spot array)"),
        "C5": (float(np.median(tot)) >= adm["min_median_umi"],
               f"median UMI {np.median(tot):.0f} (threshold {adm['min_median_umi']})"),
        "C6": (n_endo >= adm["min_endothelial_spots"],
               f"{n_endo} spots with endothelial signal (threshold {adm['min_endothelial_spots']}; "
               f"{len(pres)}/{len(ENDOTHELIAL)} markers matched)"),
        "C7": (False, "treatment before sampling not reported; waived at cohort level (see config)"),
    }
    waived = set(ov.get("waived_criteria", []))
    failed = [k for k, (ok_, _) in checks.items() if not ok_ and k not in waived]
    record = dict(slide=sid, checks={k: {"pass": bool(v[0]), "detail": v[1]} for k, v in checks.items()},
                  failed=failed, waived=sorted(waived), admitted=not failed, forced=False, force_reason=None,
                  cohort=COHORT, gene_name_note=note,
                  thresholds_applied={k: adm[k] for k in ("min_spots_visium", "min_spots_legacy_st",
                                                          "min_median_umi", "min_endothelial_spots")},
                  thresholds_overridden=applied,
                  meta=stamp_run(cfg, {"module": "M48-admission", "slide": sid}))
    save_json(P.interim / f"{sid}.admission.json", record)
    if failed:
        print(f"{sid}: NOT ADMITTED ({failed}); recorded")
        return dict(sid=sid, admitted=False, failed=failed)

    # ---- M1 / M2 without Scanpy ----
    C, ks, kg, varq = lite.qc_filter(counts, var, cfg["qc"]["min_counts"], cfg["qc"]["min_cells_per_gene"])
    X = lite.normalize_log1p(C)
    obsq = [o for o, k in zip(obs, ks) if k]
    xy = xy_arr[ks]
    d_nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    scale = cfg["qc"]["spacing_um"] / max(float(np.median(d_nn)), 1e-9)
    xy_um = xy * scale
    hyp, hyp_src = _hypoxia(cfg)
    scores, skipped = lite.score_all(X, varq, "melanoma", hypoxia_genes=hyp,
                                     min_genes=cfg["signatures"]["min_genes_matched"],
                                     ctrl_size=cfg["signatures"]["ctrl_size"], seed=cfg["seed"])
    rn = lite.rank_normalized(scores)

    # ---- M3 graph ----
    A, D = build_graph_radius(xy_um, radius_um=cfg["graph"]["radius_um"])
    rep_g = graph_report(A)
    col = lambda k: rn.get(f"{k}_n", np.full(len(obsq), 0.5))  # noqa: E731
    ss = define_source_sink(A, col("Endothelial"), col("T_NK"), col("Malignant"),
                            **{k: cfg["source_sink"][k] for k in ("q_vessel", "q_immune_nbr", "q_malig", "q_core")},
                            fallback_border_coords=xy_um if cfg["source_sink"]["use_border_fallback"] else None)
    save_graph(P.graph(sid), A, D, ss["source"], ss["sink"], ss["vessel"],
               meta={**rep_g, "radius_um": cfg["graph"]["radius_um"], "used_fallback": ss["used_fallback"],
                     **stamp_run(cfg, {"module": "M48-graph", "slide": sid})})

    # ---- M4 barriers (operators and parameters unchanged) ----
    from sparta.barrier import _SCORE_KEYS
    filled = [k for k, c in _SCORE_KEYS.items() if f"{c}" not in rn]
    S = {k: np.asarray(rn[c], float) if c in rn else np.full(len(obsq), 0.5) for k, c in _SCORE_KEYS.items()}
    b = cfg["barrier"]
    bc = compute_b_cell(A, S["ecm"], S["caf"], ss["source"], ss["sink"], **b["b_cell"])
    bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], ss["vessel"], **b["b_mab"])
    bt = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], ss["vessel"], D=D, **b["b_meta"])
    save_barriers(P.barrier(sid), {**bc, **bm, **bt}, P.mincut(sid))
    save_json(P.interim / f"{sid}.barrier_meta.json",
              stamp_run(cfg, {"module": "M48-barrier", "slide": sid, "b_cell": bc["b_cell"],
                              "max_flow": bc["max_flow"], "r_nm": b["b_mab"]["r_nm"],
                              "n_unreachable": int((~bm["reachable"]).sum()),
                              "d_vessel_mode": bt["d_vessel_mode"], "filled_keys": filled,
                              "degraded": bool(filled)}))

    # ---- node table (same keys as run_35) ----
    from scripts_run35_genes import GENES  # noqa: E402  (module created at import time below)
    out = {f"obs__{k}": np.asarray(v, float) for k, v in scores.items()}
    out.update({f"obs__{k}": np.asarray(v, float) for k, v in rn.items()})
    out["obs__total_counts"] = C.sum(axis=1).astype(float)
    out["obs__n_counts"] = C.sum(axis=1).astype(float)
    out["obs__n_genes_by_counts"] = (C > 0).sum(axis=1).astype(float)
    out["obs_names"] = np.asarray(obsq)
    out["obsm__spatial_um"] = xy_um
    out["obsm__spatial"] = xy
    present = [g for g in GENES if g in varq]
    if present:
        idx = [varq.index(g) for g in present]
        out["genes"] = np.asarray(present)
        out["gene_logcp10k"] = X[:, idx].astype(np.float32)
        out["gene_counts"] = C[:, idx].astype(np.float32)
    out["total_umi"] = C.sum(axis=1).astype(float)
    meta = dict(slide=sid, n_obs=int(C.shape[0]), n_vars=int(C.shape[1]),
                numeric_obs_columns=[k[5:] for k in out if k.startswith("obs__")],
                genes_present=present, genes_missing=[g for g in GENES if g not in present],
                x_is_log1p_cp10k=True, source_file=f"{member} (ST-Melanoma-Datasets_1.zip)",
                m2_run=dict(module="M48-lite", tumor_type="melanoma", hypoxia_source=hyp_src,
                            hypoxia_n_genes=len(hyp), seed=cfg["seed"], scanpy_free=True,
                            skipped=[s[0] for s in skipped], gene_name_note=note, spacing_scale=scale))
    out["meta"] = np.array([json.dumps(meta, ensure_ascii=False)], dtype=object)
    np.savez_compressed(P.interim / f"{sid}.nodes.npz", **out)

    # ---- ledger row ----
    for r in ledger_rows:
        if r["slide_id"] == sid:
            r.update(patient=patient, replicate=rep, cancer_type="melanoma", platform="legacy_st",
                     source="Thrane et al. 2018 Cancer Res (spatialresearch.org)",
                     accession="10.1158/0008-5472.CAN-18-0747", n_spots=str(C.shape[0]),
                     n_genes=str(C.shape[1]), median_umi=f"{np.median(tot):.1f}", treatment="unknown",
                     site="lymph-node metastasis", status="replication",
                     notes=f"replication cohort (never pooled with the primary cohort); {note}; "
                           f"admission override {COHORT}")
    print(f"{sid}: {C.shape[0]} spots x {C.shape[1]} genes | graph {rep_g['n_components']} comp, "
          f"LCC {rep_g['largest_component_frac']:.2f} | vessel {ss['n_vessel']} source {ss['n_source']} "
          f"sink {ss['n_sink']} | B_cell {bc['b_cell']:.3f} | filled {filled} | skipped {[s[0] for s in skipped]}")
    return dict(sid=sid, admitted=True, filled=filled)


# --------------------------------------------------------------------------
def statistics(cfg, sids):
    import numpy as np
    from scipy.stats import spearmanr

    from sparta.barrier import compute_b_cell_field, compute_b_mab, compute_b_meta
    from sparta.io_ import Paths, load_graph, load_json, patient_map
    from sparta.node_tables import load_nodes, scores_from_nodes
    from sparta.spatial_stats import (control_projector, normal_scores, partial_spearman,
                                      partial_spearman_many, spectral_basis)
    sys.path.insert(0, str(ROOT / "scripts"))
    from run_36_patient_level import analyse  # noqa: E402
    from run_39_coupling_decomposition import bh, input_surrogate  # noqa: E402

    P = Paths(cfg)
    pmap = patient_map(P)
    bcfg = cfg["barrier"]
    cell_abl = dict(bcfg["b_cell"], b_ecm=0.0, c_caf=12.0)
    mab_abl = dict(bcfg["b_mab"], lam=0.0)
    rng = np.random.default_rng(SEED)
    per = {}
    for sid in sids:
        t0 = time.time()
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
        dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"], vessel, D=D,
                            **bcfg["b_meta"])["d_vessel_um"]
        H = control_projector(dv)
        bc = compute_b_cell_field(A, S["ecm"], S["caf"], source, **bcfg["b_cell"])["b_cell_field"]
        bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **bcfg["b_mab"])["b_mab"]
        rho = partial_spearman(bc, bm, dv)
        w, V = spectral_basis(A)
        n = len(bm)
        fh = np.abs(V.T @ bm)
        null = np.array([partial_spearman_many(bc, V @ (rng.choice([-1.0, 1.0], size=n) * fh), dv, H)[0]
                         for _ in range(N_SUR)])
        p = (np.sum(null >= rho) + 1.0) / (N_SUR + 1.0)
        # normal-score variant
        ns = normal_scores(bm)
        fhn = np.abs(V.T @ ns)
        null_ns = np.array([partial_spearman_many(bc, V @ (rng.choice([-1.0, 1.0], size=n) * fhn), dv, H)[0]
                            for _ in range(N_SUR)])
        p_ns = (np.sum(null_ns >= rho) + 1.0) / (N_SUR + 1.0)
        # structural nulls
        rho_geo = np.empty(N_NULL)
        rho_con = np.empty(N_NULL)
        for k in range(N_NULL):
            xl_s = input_surrogate(S["crosslink"], V, rng)
            ag_s = input_surrogate(S["ag_target"], V, rng)
            ecm_s = input_surrogate(S["ecm"], V, rng)
            bm_g = compute_b_mab(A, ecm_s, xl_s, ag_s, vessel, **bcfg["b_mab"])["b_mab"]
            bm_c = compute_b_mab(A, S["ecm"], xl_s, ag_s, vessel, **bcfg["b_mab"])["b_mab"]
            rho_geo[k] = partial_spearman_many(bc, bm_g, dv, H)[0]
            rho_con[k] = partial_spearman_many(bc, bm_c, dv, H)[0]
        # shared-ECM ablation with its own surrogate null
        bc0 = compute_b_cell_field(A, S["ecm"], S["caf"], source, **cell_abl)["b_cell_field"]
        bm0 = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **mab_abl)["b_mab"]
        rho_abl = partial_spearman(bc0, bm0, dv)
        fh0 = np.abs(V.T @ bm0)
        null_abl = np.array([partial_spearman_many(bc0, V @ (rng.choice([-1.0, 1.0], size=n) * fh0), dv, H)[0]
                             for _ in range(N_SUR)])
        # low B_cell (bottom quartile) and high B_mAb (top quartile); 6.25% under independence
        disc = float(np.mean((bc <= np.quantile(bc, 0.25)) & (bm >= np.quantile(bm, 0.75))))
        cor = lambda a, b: float(spearmanr(S[a], S[b])[0]) if np.ptp(S[a]) > 0 and np.ptp(S[b]) > 0 else None  # noqa: E731
        per[sid] = dict(
            patient=pmap.get(sid, sid), n_nodes=int(n), filled_keys=filled,
            n_source=int(len(source)), n_sink=int(len(sink)), n_vessel=int(len(vessel)),
            real_rho_partial=float(rho), null_mean=float(null.mean()), null_std=float(null.std()),
            null_ci95=[float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))],
            empirical_p_one_sided=float(p), empirical_p_nscore=float(p_ns),
            geometry_null_mean=float(rho_geo.mean()), geometry_null_sd=float(rho_geo.std()),
            geometry_null_ci95=[float(np.percentile(rho_geo, 2.5)), float(np.percentile(rho_geo, 97.5))],
            p_vs_geometry_null=float((np.sum(rho_geo >= rho) + 1.0) / (N_NULL + 1.0)),
            construction_null_mean=float(rho_con.mean()), construction_null_sd=float(rho_con.std()),
            construction_null_ci95=[float(np.percentile(rho_con, 2.5)), float(np.percentile(rho_con, 97.5))],
            p_vs_construction_null=float((np.sum(rho_con >= rho) + 1.0) / (N_NULL + 1.0)),
            excess_over_construction_null=float(rho - rho_con.mean()),
            share_reproduced_by_construction=float(rho_con.mean() / rho) if rho > 0 else None,
            ablated_rho=float(rho_abl), ablated_p=float((np.sum(null_abl >= rho_abl) + 1.0) / (N_SUR + 1.0)),
            ablated_null_sd=float(null_abl.std()),
            discordant_low_bcell_high_bmab=disc,
            input_correlations=dict(ecm_caf=cor("ecm", "caf"), ecm_crosslink=cor("ecm", "crosslink")),
        )
        print(f"{sid:<14} rho={rho:+.3f} null {null.mean():+.3f}±{null.std():.3f} p={p:.3f} (ns {p_ns:.3f}) | "
              f"constr {rho_con.mean():+.3f} (p={per[sid]['p_vs_construction_null']:.3f}) | "
              f"ablated {rho_abl:+.3f} | {time.time() - t0:.0f}s", flush=True)
    sids_ok = list(per)
    q_bh = bh([per[s]["empirical_p_one_sided"] for s in sids_ok])
    q_ns = bh([per[s]["empirical_p_nscore"] for s in sids_ok])
    q_con = bh([per[s]["p_vs_construction_null"] for s in sids_ok])
    q_abl = bh([per[s]["ablated_p"] for s in sids_ok])
    for s, a, b_, c, d in zip(sids_ok, q_bh, q_ns, q_con, q_abl):
        per[s].update(bh_padj=float(a), bh_padj_nscore=float(b_), bh_padj_construction=float(c),
                      bh_padj_ablated=float(d))
    vals = lambda key: {s: per[s][key] for s in sids_ok}  # noqa: E731
    sds = {s: per[s]["null_std"] for s in sids_ok}
    patient = dict(
        association=analyse(vals("real_rho_partial"), sds, pmap, "replication: partial rho"),
        excess=analyse(vals("excess_over_construction_null"),
                       {s: per[s]["construction_null_sd"] for s in sids_ok}, pmap,
                       "replication: excess over construction null"),
        ablated=analyse(vals("ablated_rho"), {s: per[s]["ablated_null_sd"] for s in sids_ok}, pmap,
                        "replication: shared ECM term ablated"),
    )
    # secondary: combined with the primary cohort in an explicitly labelled analysis
    prim = load_json(P.validation("spatial_null_check.json"))["per_slide"]
    primary_patients = {pmap[s] for s in prim}
    replication_patients = {per[s]["patient"] for s in sids_ok}
    pooled_vals = {**{s: r["real_rho_partial"] for s, r in prim.items()}, **vals("real_rho_partial")}
    pooled_sds = {**{s: r["null_std"] for s, r in prim.items()}, **sds}
    patient["pooled_primary_plus_replication"] = analyse(
        pooled_vals,
        pooled_sds,
        pmap,
        f"primary ({len(primary_patients)} patients) + replication ({len(replication_patients)} patients)",
    )
    rho = np.array([per[s]["real_rho_partial"] for s in sids_ok])
    summary = dict(
        n_sections=len(sids_ok), n_patients=len({per[s]["patient"] for s in sids_ok}),
        n_pos=int((rho > 0).sum()), median_rho=float(np.median(rho)),
        min_rho=float(rho.min()), max_rho=float(rho.max()),
        n_p05=int(sum(per[s]["empirical_p_one_sided"] < 0.05 for s in sids_ok)),
        n_q05=int(sum(per[s]["bh_padj"] < 0.05 for s in sids_ok)),
        n_q05_nscore=int(sum(per[s]["bh_padj_nscore"] < 0.05 for s in sids_ok)),
        median_construction_null=float(np.median([per[s]["construction_null_mean"] for s in sids_ok])),
        median_geometry_null=float(np.median([per[s]["geometry_null_mean"] for s in sids_ok])),
        median_share_construction=float(np.median([per[s]["share_reproduced_by_construction"]
                                                   for s in sids_ok
                                                   if per[s]["share_reproduced_by_construction"] is not None])),
        median_excess=float(np.median([per[s]["excess_over_construction_null"] for s in sids_ok])),
        n_excess_pos=int(sum(per[s]["excess_over_construction_null"] > 0 for s in sids_ok)),
        n_q05_construction=int(sum(per[s]["bh_padj_construction"] < 0.05 for s in sids_ok)),
        median_ablated=float(np.median([per[s]["ablated_rho"] for s in sids_ok])),
        n_ablated_pos=int(sum(per[s]["ablated_rho"] > 0 for s in sids_ok)),
        n_ablated_q05=int(sum(per[s]["bh_padj_ablated"] < 0.05 for s in sids_ok)),
        median_r_ecm_caf=float(np.median([per[s]["input_correlations"]["ecm_caf"] for s in sids_ok])),
        n_discordant_below_chance=int(sum(per[s]["discordant_low_bcell_high_bmab"] < 0.0625 for s in sids_ok)),
        replication_declared=bool(patient["association"]["nested_model"]["ci95"][0] > 0),
    )
    return per, patient, summary


def main():
    args = parse_args()
    from sparta.io_ import Paths, load_config, save_json, stamp_run

    cfg = load_config(ROOT / "configs" / "cscc_legacy_st.yaml")
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    if args.validate:
        validate(cfg, Path(args.raw_dir))
        return
    proto = ROOT / "docs" / "replication_melanoma_protocol.md"
    sha = hashlib.sha256(proto.read_bytes()).hexdigest()
    t0 = time.time()
    built = {}
    if not args.stats_only:
        # GENES of run_35 without importing scanpy-dependent code
        src = (ROOT / "scripts" / "run_35_export_node_tables.py").read_text(encoding="utf-8")
        start = src.index("GENES = [")
        ns = {}
        exec(src[start:src.index("]", start) + 1], ns)
        import types
        mod = types.ModuleType("scripts_run35_genes")
        mod.GENES = ns["GENES"]
        sys.modules["scripts_run35_genes"] = mod
        lp = ROOT / "data" / "ledger.csv"
        with open(lp, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames
            rows = list(reader)
        with zipfile.ZipFile(ROOT / "data" / "external" / "thrane2018" / "ST-Melanoma-Datasets_1.zip") as zf:
            for sid in SECTIONS:
                built[sid] = build_section(sid, cfg, zf, rows)
        with open(lp, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
    sids = [s for s in SECTIONS if (P.interim / f"{s}.nodes.npz").exists()
            and json.loads((P.interim / f"{s}.admission.json").read_text(encoding="utf-8")).get("admitted")]
    per, patient, summary = statistics(cfg, sids)
    out = dict(protocol=dict(path="docs/replication_melanoma_protocol.md", sha256=sha),
               data=dict(source="Thrane et al. 2018 Cancer Res 78:5970-5979", doi="10.1158/0008-5472.CAN-18-0747",
                         sections=list(SECTIONS)),
               per_slide=per, patient_level=patient, summary=summary,
               settings=dict(n_sur=N_SUR, n_null=N_NULL, seed=SEED, config="configs/cscc_legacy_st.yaml",
                             cohort_override=COHORT),
               semantics=("Independent replication family: never pooled into the primary-cohort numbers except "
                          "in the explicitly labelled 'pooled_primary_plus_replication' patient-level analysis."),
               meta=stamp_run(cfg, {"module": "M48-replication", "seconds": round(time.time() - t0, 1)}))
    save_json(P.validation("replication_melanoma.json"), out)
    s = summary
    pa = patient["association"]
    print(f"\nreplication: {s['n_pos']}/{s['n_sections']} positive, median {s['median_rho']:+.3f}, "
          f"q<0.05 {s['n_q05']}/{s['n_sections']}; patients {pa['sign_flip']['n_positive']}/{pa['n_patients']} "
          f"positive, sign-flip p={pa['sign_flip']['p_one_sided']:.4f}; pooled {pa['nested_model']['mean']:+.3f} "
          f"CI {pa['nested_model']['ci95']}")
    pp = patient["pooled_primary_plus_replication"]
    print(f"pooled {pp['n_patients']} patients: {pp['nested_model']['mean']:+.3f} CI {pp['nested_model']['ci95']} "
          f"sign-flip p={pp['sign_flip']['p_one_sided']:.5f} ({pp['sign_flip']['n_positive']}/{pp['n_patients']})")
    print(f"construction share median {s['median_share_construction']:.2f}; excess median "
          f"{s['median_excess']:+.3f}; replication declared: {s['replication_declared']}")


if __name__ == "__main__":
    main()
