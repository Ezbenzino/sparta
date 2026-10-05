#!/usr/bin/env python
"""
run_47_codex_validation.py -- SPARTA's cellular barrier against measured CD8+ T-cell positions
==============================================================================================
Inputs (read-only):
    data/external/codex_crc/CRC_clusters_neighborhoods_markers.csv
        Schürch et al. 2020 Cell 182:1341 (CODEX, colorectal cancer, 140 TMA cores, 35 patients;
        TCIA collection CRC_FFPE-CODEX_CELLNEIGHS, doi:10.7937/TCIA.2020.FQN0-0326, CC BY 4.0)
    docs/codex_validation_protocol.md      pre-specified protocol (its SHA-256 is recorded)
Outputs:
    results/validation/codex_validation.json   summary statistics (primary, comparators, secondary,
                                               sensitivity) and metadata
    results/validation/codex_cores.csv          one row per core and analysis variant

Why
---
The section analyses have no measured transport outcome. At single-cell resolution every T cell's
position is measured, so immune cells can be withheld from the graph and used as a held-out outcome:
the barrier is computed only from the tumour/stromal scaffold (Collagen IV and αSMA/vimentin
protein, unchanged Eq. 1 parameters), and the CD8+ T-cell density in the tumour core is compared
with it afterwards. See docs/codex_validation_protocol.md for every decision and its rationale.

Usage:
    python scripts/run_47_codex_validation.py               # full run (≈10–15 min, 2 workers)
    python scripts/run_47_codex_validation.py --jobs 4
    python scripts/run_47_codex_validation.py --cores 1_A 2_B --no-sensitivity   # smoke test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PX_UM = 0.37744                      # µm per pixel (CODEX, 20x/0.75, Schürch et al. 2020)
SCAFFOLD = ("tumor cells", "stroma", "smooth muscle", "vasculature", "lymphatics",
            "adipocytes", "nerves", "undefined")
TCELLS = ("CD8+ T cells", "CD4+ T cells CD45RO+", "CD4+ T cells", "CD4+ T cells GATA3+",
          "CD3+ T cells", "Tregs")
COL = dict(
    x="X:X", y="Y:Y", lab="ClusterName", spot="spots", patient="patients", group="groups",
    col4="Collagen IV - bas. memb.:Cyc_12_ch_2", asma="aSMA - smooth muscle:Cyc_11_ch_2",
    vim="Vimentin - cytoplasm:Cyc_8_ch_2",
)
A_CAP, B_ECM, C_CAF = 3.0, 8.0, 4.0  # Eq. (1), unchanged
VARIANTS = {   # name -> (max_len_um, f_mode, sink_mode, outcome, b, c)
    "primary":        (50.0, "asma_vim", "core", "cd8", B_ECM, C_CAF),
    "edge30":         (30.0, "asma_vim", "core", "cd8", B_ECM, C_CAF),
    "edge100":        (100.0, "asma_vim", "core", "cd8", B_ECM, C_CAF),
    "f_asma_only":    (50.0, "asma", "core", "cd8", B_ECM, C_CAF),
    "sinks_all_tumour": (50.0, "asma_vim", "all", "cd8", B_ECM, C_CAF),
    "outcome_all_t":  (50.0, "asma_vim", "core", "tcell", B_ECM, C_CAF),
    # post hoc (added 2026-10-05 after the first summary): Eq. (1) weights halved / increased by half,
    # and one input only with the total stromal weight kept (the convention of the ECM ablation)
    "eq1_half":       (50.0, "asma_vim", "core", "cd8", 4.0, 2.0),
    "eq1_x1.5":       (50.0, "asma_vim", "core", "cd8", 12.0, 6.0),
    "eq1_matrix_only": (50.0, "asma_vim", "core", "cd8", 12.0, 0.0),
    "eq1_fibroblast_only": (50.0, "asma_vim", "core", "cd8", 0.0, 12.0),
}
POST_HOC_VARIANTS = ("eq1_half", "eq1_x1.5", "eq1_matrix_only", "eq1_fibroblast_only")
INCL = dict(min_nodes=100, min_sources=5, min_tumour=25, min_sinks=10, min_cd8=10)
N_BOOT = 2000
SEED = 20261005


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--cores", nargs="*", default=None, help="subset of core ids (e.g. 1_A)")
    ap.add_argument("--no-sensitivity", action="store_true")
    ap.add_argument("--summary-only", action="store_true",
                    help="recompute the statistics from results/validation/codex_cores.csv")
    ap.add_argument("--csv", default=str(ROOT / "data/external/codex_crc/CRC_clusters_neighborhoods_markers.csv"))
    return ap.parse_args()


# --------------------------------------------------------------------------
# per-core analysis
# --------------------------------------------------------------------------
def analyse_core(payload):
    import numpy as np
    from scipy.spatial import cKDTree

    from sparta.barrier import compute_b_cell, compute_b_cell_field
    from sparta.cellgraph import (assign_nearest, contact_enrichment_z, core_sinks, cut_partition,
                                  delaunay_graph, largest_component, rank01, region_density_ratio,
                                  ripley_l, tissue_raster, uniform_max_flow)

    core, variant, d = payload
    max_len, f_mode, sink_mode, outcome, b_w, c_w = VARIANTS[variant]
    t0 = time.time()
    xy = np.column_stack([d["x"], d["y"]]) * PX_UM
    lab = d["lab"]
    ok = lab != "dirt"
    row = dict(core=core, patient=int(d["patient"]), group=int(d["group"]), variant=variant,
               n_cells=int(ok.sum()))
    scaf = np.isin(lab, SCAFFOLD) & ok
    S = np.flatnonzero(scaf)
    row["n_scaffold"] = int(len(S))
    if len(S) < INCL["min_nodes"]:
        row.update(included=False, reason="scaffold<100")
        return row
    A, D = delaunay_graph(xy[S], max_len)
    lcc = largest_component(A)
    keep = np.flatnonzero(lcc)
    A = A[keep][:, keep].tocsr()
    S = S[keep]
    xyL, labL = xy[S], lab[S]
    n = len(S)
    row.update(n_lcc=int(n), lcc_frac=float(n / scaf.sum()))
    tumour = labL == "tumor cells"
    vessel = labL == "vasculature"
    source = np.flatnonzero(vessel)
    sink = core_sinks(A, tumour) if sink_mode == "core" else np.flatnonzero(tumour)
    sink = np.setdiff1d(sink, source)
    tcell_types = ("CD8+ T cells",) if outcome == "cd8" else TCELLS
    tc = ok & np.isin(lab, tcell_types)
    row.update(n_sources=int(len(source)), n_tumour=int(tumour.sum()), n_sinks=int(len(sink)),
               n_outcome_cells=int(tc.sum()))
    fail = []
    if n < INCL["min_nodes"]:
        fail.append("lcc<100")
    if len(source) < INCL["min_sources"]:
        fail.append("sources<5")
    if tumour.sum() < INCL["min_tumour"]:
        fail.append("tumour<25")
    if len(sink) < INCL["min_sinks"]:
        fail.append("sinks<10")
    if tc.sum() < INCL["min_cd8"]:
        fail.append("outcome_cells<10")
    if fail:
        row.update(included=False, reason=";".join(fail))
        return row
    row["included"] = True

    # ---- inputs: protein ranks over LCC scaffold nodes ----
    E = rank01(d["col4"][S])
    if f_mode == "asma_vim":
        F = rank01(0.5 * (rank01(d["asma"][S]) + rank01(d["vim"][S])))
    else:
        F = rank01(d["asma"][S])

    # ---- SPARTA cellular barrier (unchanged operator) ----
    res = compute_b_cell(A, E, F, source, sink, a=A_CAP, b_ecm=b_w, c_caf=c_w,
                         return_capacity=True)
    pairs, cap = res["edge_pairs"], res["edge_capacity"]
    nodes = np.arange(n)
    f_sp, side_sp = cut_partition(pairs, cap, nodes, source, sink)
    f_open, side_open = uniform_max_flow(pairs, nodes, source, sink, 1.0 / (1.0 + np.exp(-A_CAP)))
    field = compute_b_cell_field(A, E, F, source, a=A_CAP, b_ecm=b_w, c_caf=c_w)["b_cell_field"]
    row.update(max_flow=float(res["max_flow"]), max_flow_check=float(f_sp), max_flow_open=float(f_open),
               b_cell=float(res["b_cell"]), log_b_cell=float(np.log(res["b_cell"])),
               log_b_rel=float(np.log(f_open / res["max_flow"])),
               n_cut_edges=int(len(res["cut_edges"])),
               field_core=float(np.median(field[sink])))

    # ---- held-out outcome ----
    P = tissue_raster(xy[ok], 5.0, 10.0)
    cell_area = 25.0                                      # µm² per raster point
    jr = assign_nearest(P, xyL)
    jc = assign_nearest(xy[tc], xyL)
    in_core = np.zeros(n, bool)
    in_core[sink] = True
    a_core = in_core[jr].sum() * cell_area
    a_rest = (~in_core[jr]).sum() * cell_area
    c_core = in_core[jc]
    row.update(area_tissue_mm2=float(len(P) * cell_area / 1e6), area_core_mm2=float(a_core / 1e6),
               n_cd8_core=int(c_core.sum()), n_cd8_rest=int((~c_core).sum()),
               ir=region_density_ratio(c_core, a_core, ~c_core, a_rest),
               cd8_core_density_mm2=float(c_core.sum() / max(a_core, 1e-9) * 1e6))
    in_tum = tumour.copy()
    c_tum = in_tum[jc]
    row["ir_tumour"] = region_density_ratio(c_tum, in_tum[jr].sum() * cell_area,
                                            ~c_tum, (~in_tum[jr]).sum() * cell_area)
    # blockade line: CD8 density beyond the cut (sink side) relative to the source side
    for name, side in (("sparta", side_sp[:n]), ("open", side_open[:n])):
        beyond = ~side
        cb = beyond[jc]
        row[f"cut_ratio_{name}"] = region_density_ratio(cb, beyond[jr].sum() * cell_area,
                                                        ~cb, (~beyond[jr]).sum() * cell_area)
        row[f"beyond_area_frac_{name}"] = float(beyond[jr].mean())
    row["cut_ratio_diff"] = row["cut_ratio_sparta"] - row["cut_ratio_open"]

    # ---- comparators (larger = more barrier) ----
    m = (B_ECM * E + C_CAF * F) / (B_ECM + C_CAF)
    stromal = np.isin(labL, ("stroma", "smooth muscle"))
    row["global_matrix"] = float(m.mean())
    nt = ~tumour
    dist_t, _ = cKDTree(xyL[tumour]).query(xyL[nt], k=1)
    peri = np.flatnonzero(nt)[dist_t <= 30.0]
    row["peritumoural_matrix"] = float(m[peri].mean()) if len(peri) else float("nan")
    row["stromal_fraction"] = float(stromal.mean())
    row["tumour_fraction"] = float(tumour.mean())
    row["contact_enrichment_z"] = contact_enrichment_z(A, tumour, stromal, n_perm=1000, seed=SEED)
    rich = E >= 0.7
    row["ripley_l_matrix"] = ripley_l(xyL[rich], 50.0, len(P) * cell_area)
    # post hoc geometry-only comparators (added 2026-10-05, after the outcome had been seen; reported
    # as post hoc): distance from vessels to the tumour boundary, and depth of the tumour core
    t_tree = cKDTree(xyL[tumour])
    row["vessel_tumour_distance"] = float(np.median(t_tree.query(xyL[source], k=1)[0]))
    row["core_depth"] = float(np.median(cKDTree(xyL[~tumour]).query(xyL[sink], k=1)[0])) if (~tumour).any() \
        else float("nan")
    row["seconds"] = float(time.time() - t0)
    return row


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def spearman(x, y):
    import numpy as np
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return float("nan"), float("nan")
    r, p = spearmanr(x[ok], y[ok])
    return float(r), float(p)


def cluster_boot(df, xcol, ycol, n_boot=N_BOOT, seed=SEED, fn=None):
    """Patient-cluster bootstrap of a core-level statistic (default: Spearman rho)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    pats = df["patient"].unique()
    groups = {p: df.index[df["patient"] == p].to_numpy() for p in pats}
    fn = fn or (lambda sub: spearman(sub[xcol], sub[ycol])[0])
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(pats, size=len(pats), replace=True)
        idx = np.concatenate([groups[p] for p in pick])
        out[b] = fn(df.loc[idx])
    out = out[np.isfinite(out)]
    return [float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))]


def patient_level(df, xcol, ycol, alternative="less"):
    """Spearman of patient means; one-sided Monte-Carlo permutation p (100,000 shuffles)."""
    import numpy as np
    from scipy.stats import rankdata
    g = df.groupby("patient")[[xcol, ycol]].mean().dropna()
    r, p2 = spearman(g[xcol], g[ycol])
    rng = np.random.default_rng(SEED)
    rx, ry = rankdata(g[xcol]), rankdata(g[ycol])
    rx = (rx - rx.mean()) / rx.std()
    ry = (ry - ry.mean()) / ry.std()
    obs = float(np.mean(rx * ry))
    perm = np.array([np.mean(rx * rng.permutation(ry)) for _ in range(100000)])
    if alternative == "less":
        p1 = float((np.sum(perm <= obs) + 1) / (len(perm) + 1))
    else:
        p1 = float((np.sum(perm >= obs) + 1) / (len(perm) + 1))
    return dict(rho=r, p_two_sided=p2, p_one_sided_perm=p1, n_patients=int(len(g)))


def partial_rank(df, xcol, ycol, covs):
    """Spearman of rank-regression residuals of x and y on the covariates."""
    import numpy as np
    from scipy.stats import rankdata
    sub = df[[xcol, ycol] + covs].dropna()
    Z = np.column_stack([np.ones(len(sub))] + [rankdata(sub[c]) for c in covs])
    H = Z @ np.linalg.pinv(Z)
    rx = rankdata(sub[xcol]) - H @ rankdata(sub[xcol])
    ry = rankdata(sub[ycol]) - H @ rankdata(sub[ycol])
    return spearman(rx, ry)[0]


def summarise(df, variant_rows):
    import numpy as np
    from scipy.stats import mannwhitneyu, wilcoxon
    out = {}
    prim = df[(df.variant == "primary") & (df.included)].copy().reset_index(drop=True)
    allp = df[df.variant == "primary"]
    out["counts"] = dict(
        n_cores_total=int(allp["core"].nunique()), n_cores_included=int(len(prim)),
        n_patients_total=int(allp["patient"].nunique()), n_patients_included=int(prim["patient"].nunique()),
        n_clr_patients=int(prim.loc[prim.group == 1, "patient"].nunique()),
        n_dii_patients=int(prim.loc[prim.group == 2, "patient"].nunique()),
        exclusion_reasons=allp.loc[~allp.included, "reason"].value_counts().to_dict(),
        max_flow_check_max_abs_diff=float(np.nanmax(np.abs(prim["max_flow"] - prim["max_flow_check"]))),
    )
    desc = {}
    for c in ("n_lcc", "lcc_frac", "n_sources", "n_sinks", "n_outcome_cells", "ir", "log_b_rel",
              "log_b_cell", "cd8_core_density_mm2", "n_cut_edges", "seconds"):
        v = prim[c].astype(float)
        desc[c] = dict(median=float(v.median()), q25=float(v.quantile(.25)), q75=float(v.quantile(.75)),
                       min=float(v.min()), max=float(v.max()))
    out["describe"] = desc

    preds = [("log_b_rel", "SPARTA relative barrier (min cut / matrix-free cut)"),
             ("log_b_cell", "SPARTA B_cell (min cut, absolute)"),
             ("field_core", "SPARTA B_cell field at the core"),
             ("peritumoural_matrix", "Peritumoural matrix mean"),
             ("global_matrix", "Global matrix mean"),
             ("stromal_fraction", "Stromal fraction"),
             ("contact_enrichment_z", "Tumour-stroma contact enrichment"),
             ("ripley_l_matrix", "Ripley's L of matrix-rich cells"),
             ("vessel_tumour_distance", "Vessel-to-tumour distance (post hoc)"),
             ("core_depth", "Tumour-core depth (post hoc)")]
    post_hoc_cols = {"vessel_tumour_distance", "core_depth"}
    comp = {}
    for col, label in preds:
        if float(prim[col].std()) < 1e-9:
            comp[col] = dict(label=label, degenerate=True,
                             note="constant across cores by construction (mean of within-core ranks)")
            continue
        if col not in prim or prim[col].isna().all():
            continue
        r_core, p_core = spearman(prim[col], prim["ir"])
        comp[col] = dict(label=label, rho_core=r_core, p_core_naive=p_core,
                         ci_core=cluster_boot(prim, col, "ir"),
                         patient=patient_level(prim, col, "ir", "less"), post_hoc=col in post_hoc_cols)
    out["primary_and_comparators"] = comp

    # SPARTA vs each comparator: patient-cluster bootstrap of |rho_Brel| - |rho_comp|
    rng = np.random.default_rng(SEED + 1)
    pats = prim["patient"].unique()
    groups = {p: prim.index[prim["patient"] == p].to_numpy() for p in pats}
    diffs = {c: [] for c, _ in preds if c != "log_b_rel" and c in comp and not comp[c].get("degenerate")}
    for _ in range(N_BOOT):
        pick = rng.choice(pats, size=len(pats), replace=True)
        sub = prim.loc[np.concatenate([groups[p] for p in pick])]
        r0 = spearman(sub["log_b_rel"], sub["ir"])[0]
        for c in diffs:
            diffs[c].append(-r0 - (-spearman(sub[c], sub["ir"])[0]))
    out["sparta_minus_comparator_neg_rho"] = {
        c: dict(median=float(np.nanmedian(v)), ci=[float(np.nanquantile(v, .025)), float(np.nanquantile(v, .975))],
                frac_gt0=float(np.nanmean(np.asarray(v) > 0))) for c, v in diffs.items()}

    covs = ["tumour_fraction", "stromal_fraction", "global_matrix", "peritumoural_matrix"]
    pr = partial_rank(prim, "log_b_rel", "ir", covs)
    ci = cluster_boot(prim, "log_b_rel", "ir", fn=lambda s: partial_rank(s, "log_b_rel", "ir", covs))
    out["added_value_partial"] = dict(covariates=covs, rho=pr, ci=ci)

    # ---- post hoc (not in the protocol; labelled as such everywhere) ----
    # the global-matrix covariate is constant (within-core ranks), so the pre-specified partial
    # reduces to tumour fraction, stromal fraction and peritumoural matrix; intermixing (contact
    # enrichment) turned out to be the strongest single correlate of the outcome, so the added value
    # of the barrier is also reported given it.
    ph = {}
    for name, cv in (("given_composition", ["tumour_fraction", "stromal_fraction", "peritumoural_matrix"]),
                     ("given_intermixing", ["contact_enrichment_z"]),
                     ("given_composition_and_intermixing",
                      ["tumour_fraction", "stromal_fraction", "peritumoural_matrix", "contact_enrichment_z"])):
        ph[f"b_rel_{name}"] = dict(covariates=cv, rho=partial_rank(prim, "log_b_rel", "ir", cv),
                                   ci=cluster_boot(prim, "log_b_rel", "ir",
                                                   fn=lambda s_, cv=cv: partial_rank(s_, "log_b_rel", "ir", cv)))
    for x in ("peritumoural_matrix", "contact_enrichment_z"):
        cv = [c for c in ("tumour_fraction", "stromal_fraction", "peritumoural_matrix",
                          "contact_enrichment_z", "log_b_rel") if c != x]
        ph[f"{x}_given_others"] = dict(covariates=cv, rho=partial_rank(prim, x, "ir", cv),
                                       ci=cluster_boot(prim, x, "ir",
                                                       fn=lambda s_, cv=cv, x=x: partial_rank(s_, x, "ir", cv)))
    # geometry-only distances (post hoc comparators): added value of the barrier given them
    for name, cv in (("given_distances", ["vessel_tumour_distance", "core_depth"]),
                     ("given_composition_and_distances",
                      ["tumour_fraction", "stromal_fraction", "peritumoural_matrix", "vessel_tumour_distance",
                       "core_depth"])):
        ph[f"b_rel_{name}"] = dict(covariates=cv, rho=partial_rank(prim, "log_b_rel", "ir", cv),
                                   ci=cluster_boot(prim, "log_b_rel", "ir",
                                                   fn=lambda s_, cv=cv: partial_rank(s_, "log_b_rel", "ir", cv)))
    # added value given each single simple summary (post hoc; answers "what does the cut add over X?")
    for x in ("peritumoural_matrix", "vessel_tumour_distance", "core_depth", "contact_enrichment_z"):
        ph[f"b_rel_given_{x}"] = dict(covariates=[x], rho=partial_rank(prim, "log_b_rel", "ir", [x]),
                                      ci=cluster_boot(prim, "log_b_rel", "ir",
                                                      fn=lambda s_, x=x: partial_rank(s_, "log_b_rel", "ir", [x])))
    allc = ["peritumoural_matrix", "vessel_tumour_distance", "core_depth", "contact_enrichment_z"]
    ph["b_rel_given_all_simple"] = dict(covariates=allc, rho=partial_rank(prim, "log_b_rel", "ir", allc),
                                        ci=cluster_boot(prim, "log_b_rel", "ir",
                                                        fn=lambda s_: partial_rank(s_, "log_b_rel", "ir", allc)))
    ph["spearman_matrix"] = prim[["log_b_rel", "log_b_cell", "field_core", "peritumoural_matrix",
                                  "contact_enrichment_z", "stromal_fraction", "tumour_fraction",
                                  "vessel_tumour_distance", "core_depth",
                                  "ir"]].corr(method="spearman").round(4).to_dict()
    out["post_hoc"] = ph

    # S1 absolute access
    out["s1_absolute"] = dict(
        rho_core=spearman(prim["log_b_cell"], prim["cd8_core_density_mm2"])[0],
        ci_core=cluster_boot(prim, "log_b_cell", "cd8_core_density_mm2"),
        patient=patient_level(prim, "log_b_cell", "cd8_core_density_mm2", "less"))
    # S2 blockade line
    g = prim.groupby("patient")[["cut_ratio_sparta", "cut_ratio_open", "cut_ratio_diff"]].mean()
    w = wilcoxon(g["cut_ratio_diff"], alternative="less")
    out["s2_blockade_line"] = dict(
        median_cut_ratio_sparta=float(prim["cut_ratio_sparta"].median()),
        median_cut_ratio_open=float(prim["cut_ratio_open"].median()),
        median_diff_core=float(prim["cut_ratio_diff"].median()),
        median_diff_patient=float(g["cut_ratio_diff"].median()),
        n_patients_diff_lt0=int((g["cut_ratio_diff"] < 0).sum()), n_patients=int(len(g)),
        wilcoxon_one_sided_p=float(w.pvalue),
        median_beyond_area_frac_sparta=float(prim["beyond_area_frac_sparta"].median()),
        median_beyond_area_frac_open=float(prim["beyond_area_frac_open"].median()))
    # S3 CLR vs DII (exploratory)
    gp = prim.groupby("patient").agg(log_b_rel=("log_b_rel", "mean"), ir=("ir", "mean"), group=("group", "first"))
    a, b = gp.loc[gp.group == 1, "log_b_rel"], gp.loc[gp.group == 2, "log_b_rel"]
    ai, bi = gp.loc[gp.group == 1, "ir"], gp.loc[gp.group == 2, "ir"]
    out["s3_clr_vs_dii"] = dict(
        median_log_b_rel_clr=float(a.median()), median_log_b_rel_dii=float(b.median()),
        mannwhitney_p_b_rel=float(mannwhitneyu(a, b, alternative="two-sided").pvalue),
        median_ir_clr=float(ai.median()), median_ir_dii=float(bi.median()),
        mannwhitney_p_ir=float(mannwhitneyu(ai, bi, alternative="two-sided").pvalue))
    # tumour-region outcome (secondary definition of the region)
    out["ir_tumour_region"] = dict(rho_core=spearman(prim["log_b_rel"], prim["ir_tumour"])[0],
                                   patient=patient_level(prim, "log_b_rel", "ir_tumour", "less"))
    # sensitivity variants
    sens = {}
    for v in variant_rows:
        if v == "primary":
            continue
        sub = df[(df.variant == v) & (df.included)].reset_index(drop=True)
        if len(sub) < 10:
            continue
        sens[v] = dict(post_hoc=bool(v in POST_HOC_VARIANTS),
                       n_cores=int(len(sub)), n_patients=int(sub["patient"].nunique()),
                       rho_core=spearman(sub["log_b_rel"], sub["ir"])[0],
                       patient=patient_level(sub, "log_b_rel", "ir", "less"),
                       rho_core_b_cell_density=spearman(sub["log_b_cell"], sub["cd8_core_density_mm2"])[0])
    out["sensitivity"] = sens
    return out


def main():
    args = parse_args()
    import numpy as np
    import pandas as pd

    from sparta.io_ import Paths, load_config, save_json, stamp_run

    cfg = load_config(ROOT / "configs" / "default.yaml")
    P = Paths(cfg)
    proto = ROOT / "docs" / "codex_validation_protocol.md"
    sha = hashlib.sha256(proto.read_bytes()).hexdigest()
    t0 = time.time()
    if args.summary_only:
        res = pd.read_csv(P.validation("codex_cores.csv"))
        res["included"] = res["included"].astype(bool)
        variants = sorted(res["variant"].unique(), key=lambda v: list(VARIANTS).index(v))
        return finish(args, cfg, P, sha, res, variants, t0)
    use = list(COL.values())
    df = pd.read_csv(args.csv, usecols=use)
    cores = sorted(df[COL["spot"]].unique(), key=lambda s: (int(s.split("_")[0]), s))
    if args.cores:
        cores = [c for c in cores if c in set(args.cores)]
    variants = ["primary"] if args.no_sensitivity else list(VARIANTS)
    jobs = []
    for c in cores:
        sub = df[df[COL["spot"]] == c]
        d = dict(x=sub[COL["x"]].to_numpy(float), y=sub[COL["y"]].to_numpy(float),
                 lab=sub[COL["lab"]].astype(str).to_numpy(), col4=sub[COL["col4"]].to_numpy(float),
                 asma=sub[COL["asma"]].to_numpy(float), vim=sub[COL["vim"]].to_numpy(float),
                 patient=int(sub[COL["patient"]].iloc[0]), group=int(sub[COL["group"]].iloc[0]))
        for v in variants:
            jobs.append((c, v, d))
    print(f"{len(cores)} cores x {len(variants)} variants = {len(jobs)} jobs; {args.jobs} workers", flush=True)
    rows = []
    if args.jobs > 1:
        from multiprocessing import Pool
        with Pool(args.jobs) as pool:
            for i, r in enumerate(pool.imap_unordered(analyse_core, jobs, chunksize=2)):
                rows.append(r)
                if (i + 1) % 50 == 0:
                    print(f"  {i + 1}/{len(jobs)}  {time.time() - t0:.0f}s", flush=True)
    else:
        for i, j in enumerate(jobs):
            rows.append(analyse_core(j))
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(jobs)}  {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(rows).sort_values(["variant", "patient", "core"]).reset_index(drop=True)
    res["included"] = res["included"].astype(bool)
    out_csv = P.validation("codex_cores.csv")
    res.to_csv(out_csv, index=False)
    return finish(args, cfg, P, sha, res, variants, t0)


def finish(args, cfg, P, sha, res, variants, t0):
    import numpy as np
    from sparta.io_ import save_json, stamp_run
    summary = summarise(res, variants)
    out = dict(
        protocol=dict(path="docs/codex_validation_protocol.md", sha256=sha),
        data=dict(source="Schürch et al. 2020 Cell 182:1341-1359 (CODEX, colorectal cancer)",
                  doi_paper="10.1016/j.cell.2020.07.005",
                  doi_single_cell_table="10.17632/mpjzbtfgfr.1 (Mendeley Data, CRC_clusters_neighborhoods_markers.csv)",
                  doi_images="10.7937/TCIA.2020.FQN0-0326 (TCIA, CC BY 4.0)", px_um=PX_UM),
        settings=dict(variants={k: dict(zip(("max_len_um", "f_mode", "sink_mode", "outcome", "b", "c"), v))
                                for k, v in VARIANTS.items()},
                      post_hoc_variants=list(POST_HOC_VARIANTS),
                      eq1=dict(a=A_CAP, b=B_ECM, c=C_CAF), inclusion=INCL, n_boot=N_BOOT, seed=SEED,
                      scaffold=list(SCAFFOLD), outcome_cells=dict(cd8=["CD8+ T cells"], tcell=list(TCELLS))),
        summary=summary,
        semantics=("Barrier computed from the non-immune scaffold only (Collagen IV, aSMA/vimentin protein; "
                   "Eq. 1 parameters unchanged); CD8+ T-cell positions withheld and used as the outcome. "
                   "IR = log2 relative CD8+ density in the tumour core (sink Voronoi region) vs the rest of the "
                   "tissue; negative = depleted. log_b_rel = log(max flow on matrix-free capacities / max flow). "
                   "Inference uses patients (35) as the unit; core-level CIs are patient-cluster bootstraps."),
        meta=stamp_run(cfg, dict(module="M47-codex-validation", seconds=round(time.time() - t0, 1))),
    )
    save_json(P.validation("codex_validation.json"), out)
    s = summary
    pc = s["primary_and_comparators"]["log_b_rel"]
    print(f"\nincluded cores {s['counts']['n_cores_included']}/{s['counts']['n_cores_total']}, "
          f"patients {s['counts']['n_patients_included']}")
    print(f"PRIMARY log B_rel vs IR: core rho={pc['rho_core']:.3f} CI {pc['ci_core']}, "
          f"patient rho={pc['patient']['rho']:.3f} p1={pc['patient']['p_one_sided_perm']:.4f}")
    for c, v in s["primary_and_comparators"].items():
        if v.get("degenerate"):
            print(f"  {c:24s} degenerate ({v['note']})")
            continue
        print(f"  {c:24s} core {v['rho_core']:+.3f}  patient {v['patient']['rho']:+.3f} "
              f"(p1 {v['patient']['p_one_sided_perm']:.4f})")
    for k, v in s["post_hoc"].items():
        if k != "spearman_matrix":
            print(f"  post hoc {k}: {v['rho']:+.3f} CI {np.round(v['ci'], 3)}")
    print("partial:", s["added_value_partial"])
    print("S1:", s["s1_absolute"]["rho_core"], s["s1_absolute"]["patient"])
    print("S2:", s["s2_blockade_line"])
    print("S3:", s["s3_clr_vs_dii"])
    for k, v in s["sensitivity"].items():
        print(f"  sens {k:18s} n={v['n_cores']} core {v['rho_core']:+.3f} patient {v['patient']['rho']:+.3f} "
              f"p1 {v['patient']['p_one_sided_perm']:.4f}")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
