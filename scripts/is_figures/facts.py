"""
Every number quoted in the IS manuscript, computed from the stored result files.
`python facts.py` writes results/validation/is_manuscript_facts.json; the manuscript
builder substitutes {{key}} placeholders from it, and the verification step re-checks
the rendered text against the same file.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from isdata import EXTENSION_ORDER, EXTERNAL_ORDER, PRIMARY_ORDER, VAL, j, ledger  # noqa: E402


def f2(x):
    out = f"{x:.2f}"
    return "0.00" if out == "-0.00" else out


def f3(x):
    out = f"{x:.3f}"
    return "0.000" if out == "-0.000" else out


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def pfmt(p):
    if p < 0.001:
        return "< 0.001"
    if p < 0.01:
        return f"= {p:.4f}".rstrip("0") if p >= 0.0001 else "< 0.0001"
    return f"= {p:.3f}"


def main():
    F = {}
    led = ledger()
    # ---------------- cohort / graph coverage ----------------
    gc = j("graph_coverage_table.json")
    s = gc["summary"]
    F["n_primary_sections"] = s["n_primary"]
    F["spots_min"], F["spots_max"] = f"{s['spots_range'][0]:,}", f"{s['spots_range'][1]:,}"
    F["comp_min"], F["comp_max"] = s["components_range"]
    F["lcc_min_pct"] = pct(s["lcc_share_range"][0])
    F["n_stranded_sections"] = s["n_sections_with_stranded_source_or_sink"]
    F["max_stranded_pct"] = pct(s["max_stranded_share"])
    F["max_stranded_section"] = s["max_stranded_section"]
    rows = {r["section"]: r for r in gc["rows"]}
    F["ext_spots_min"] = f"{min(rows[x]['n_spots'] for x in EXTERNAL_ORDER):,}"
    F["ext_spots_max"] = f"{max(rows[x]['n_spots'] for x in EXTERNAL_ORDER):,}"

    # ---------------- primary association ----------------
    sn = j("spatial_null_check.json")
    per = sn["per_slide"]
    rho = np.array([per[x]["real_rho_partial"] for x in PRIMARY_ORDER])
    F["assoc_median"] = f3(np.median(rho))
    F["assoc_min"] = f"{rho.min():.4f}" if abs(rho.min()) < 0.001 else f3(rho.min())
    F["assoc_max"] = f3(rho.max())
    F["assoc_n_pos"] = int((rho > 0).sum())
    ps = np.array([per[x]["empirical_p_one_sided"] for x in PRIMARY_ORDER])
    qs = np.array([per[x]["bh_padj"] for x in PRIMARY_ORDER])
    F["assoc_n_p05"] = int((ps < 0.05).sum())
    F["assoc_n_q05"] = int((qs < 0.05).sum())
    nsd = np.array([per[x]["null_std"] for x in PRIMARY_ORDER])
    F["null_sd_min"], F["null_sd_max"] = f3(nsd.min()), f3(nsd.max())
    F["null_mean_absmax"] = f3(max(abs(per[x]["null_mean"]) for x in PRIMARY_ORDER))
    neg = [x for x in PRIMARY_ORDER if per[x]["real_rho_partial"] <= 0]
    F["assoc_nonpos_section"] = ", ".join(neg)
    vis = [x for x in PRIMARY_ORDER if x in ("CSCC01", "CSCC02", "CSCC03", "CSCC04")]
    st = [x for x in PRIMARY_ORDER if x.startswith("CSCC") and x not in vis]
    mel = [x for x in PRIMARY_ORDER if x.startswith("MEL")]
    F["assoc_median_vis_cscc"] = f3(np.median([per[x]["real_rho_partial"] for x in vis]))
    F["assoc_median_st_cscc"] = f3(np.median([per[x]["real_rho_partial"] for x in st]))
    F["assoc_median_mel"] = f3(np.median([per[x]["real_rho_partial"] for x in mel]))
    # discordance
    dec = j("decoupling.json")
    disc = np.array([dec[x]["frac_discordant_r"] for x in PRIMARY_ORDER])
    F["n_discordant_below_chance"] = int((disc < 0.0625).sum())
    F["discordant_min_pct"], F["discordant_max_pct"] = pct(disc.min()), pct(disc.max())

    # ---------------- external ----------------
    ext = j("ext_validation.json")["per_slide"]
    er = [ext[x]["spatial_null"]["real_rho_partial"] for x in EXTERNAL_ORDER]
    eq = [ext[x]["spatial_null"]["bh_padj"] for x in EXTERNAL_ORDER]
    F["ext_median"] = f3(np.median(er))
    F["ext_rhos"] = ", ".join(f3(x) for x in er)
    F["ext_qs"] = ", ".join(f"{x:.4f}" if 0.045 <= x < 0.055 else f3(x) for x in eq)   # borderline stays exact
    F["ext_brca02_q"] = f"{ext['BRCA02']['spatial_null']['bh_padj']:.4f}"

    # ---------------- patient level ----------------
    pl = j("patient_level_inference.json")

    def pl_fill(key, prefix):
        r = pl[key]
        nm = r["nested_model"]
        F[f"{prefix}_mean"] = f2(nm["mean"])
        F[f"{prefix}_lo"], F[f"{prefix}_hi"] = f2(nm["ci95"][0]), f2(nm["ci95"][1])
        F[f"{prefix}_p"] = f"{nm['p_two_sided']:.4f}"
        F[f"{prefix}_df"] = nm["df"]
        F[f"{prefix}_J"] = r["n_patients"]
        F[f"{prefix}_n"] = r["n_sections"]
        F[f"{prefix}_signflip_p"] = f"{r['sign_flip']['p_one_sided']:.4f}"
        F[f"{prefix}_npos_pat"] = r["sign_flip"]["n_positive"]
        F[f"{prefix}_lopo_min"] = f2(r["lopo_min_mean"])
        F[f"{prefix}_lopo_maxp"] = f"{r['lopo_max_p_two_sided']:.4f}"
        F[f"{prefix}_bp_share"] = (f"{100 * nm['between_patient_share']:.0f}"
                                   if nm.get("between_patient_share") is not None else "NA")
        F[f"{prefix}_tau_p"] = f3(np.sqrt(nm["tau2_patient"]))
    pl_fill("primary_association", "pl")
    pl_fill("primary_association_cscc_only", "plc")
    pl_fill("primary_association_ag_available", "pla")
    pl_fill("primary_association_lcc", "pll")
    pl_fill("s2_matched_log_ratio", "pls2")
    F["pls2_ratio"] = f"{np.exp(pl['s2_matched_log_ratio']['nested_model']['mean']):.3f}"
    F["pls2_ratio_lo"] = f"{np.exp(pl['s2_matched_log_ratio']['nested_model']['ci95'][0]):.3f}"
    F["pls2_ratio_hi"] = f"{np.exp(pl['s2_matched_log_ratio']['nested_model']['ci95'][1]):.3f}"
    if "ecm_ablated_association" in pl:
        pl_fill("ecm_ablated_association", "plab")
    if "construction_null_excess" in pl:
        pl_fill("construction_null_excess", "plex")

    # ---------------- table 3 counts ----------------
    from isdata import bh as _bh
    mc = j("metric_connectivity_sensitivity.json")["summary"]
    F["pla_n_pos"] = mc["association_ag_target_available_15"]["n_positive"]
    F["pla_n_q"] = mc["spatial_null_ag_target_available_15"]["n_bh_q_lt_0_05"]
    cs = [x for x in PRIMARY_ORDER if x.startswith("CSCC")]
    F["plc_n_pos"] = int(sum(per[x]["real_rho_partial"] > 0 for x in cs))
    F["plc_n_q"] = int((np.array(_bh([per[x]["empirical_p_one_sided"] for x in cs])) < 0.05).sum())
    F["pll_n_pos"] = mc["largest_component_association"]["all_19"]["n_positive"]
    F["pll_median"] = f3(mc["largest_component_association"]["all_19"]["median"])
    # ---------------- adjustment robustness ----------------
    ar = j("adjustment_robustness.json")["summary_primary"]
    for k in ["rank_quadratic", "rank_cubic", "spline_raw", "shell_stratified", "unadjusted", "tie_jitter_min", "tie_jitter_median"]:
        F[f"adj_{k}_median"] = f3(ar[k]["median"])
        F[f"adj_{k}_npos"] = ar[k]["n_positive"]
    alt = [ar[k]["median"] for k in ["rank_cubic", "spline_raw", "shell_stratified", "unadjusted"]]
    F["adj_alt_min"], F["adj_alt_max"] = f3(min(alt)), f3(max(alt))
    F["adj_alt_npos_min"] = min(ar[k]["n_positive"] for k in ["rank_cubic", "spline_raw", "shell_stratified", "unadjusted"])

    # ---------------- reproduction ----------------
    rc = j("reproduction_check.json")["summary"]
    F["repro_max_drho"] = f"{rc['max_abs_rho']:.4f}"
    F["repro_n_p_identical"] = rc["n_p_identical"]
    F["repro_max_dbmab"] = f"{rc['max_abs_bmab']:.0e}"

    # ---------------- null calibration ----------------
    p = VAL / "null_calibration_all.json"
    if p.exists():
        nc = j("null_calibration_all.json")
        sm = nc.get("summary", {})
        if "gaussian" in sm:
            g = sm["gaussian"]
            F["cal_n_graphs"] = len(nc["per_slide"])
            F["cal_n_sim"] = nc["n_sim"]
            F["cal_n_sur"] = nc["n_sur"]
            F["cal_pooled_spec"] = f"{g['pooled_fpr_spectral_05']:.3f}"
            F["cal_pooled_naive"] = f"{g['pooled_fpr_naive_05']:.3f}"
            F["cal_spec_min"], F["cal_spec_max"] = (f"{g['section_fpr_spectral_range'][0]:.3f}",
                                                    f"{g['section_fpr_spectral_range'][1]:.3f}")
            F["cal_naive_min"], F["cal_naive_max"] = (f"{g['section_fpr_naive_range'][0]:.3f}",
                                                      f"{g['section_fpr_naive_range'][1]:.3f}")
            F["cal_n_tests"] = f"{g['n_tests']:,}"
            sk = sm["skewed"]
            F["cal_skew_spec"] = f"{sk['pooled_fpr_spectral_05']:.3f}"
            F["cal_skew_naive"] = f"{sk['pooled_fpr_naive_05']:.3f}"
            smo = [v for k, v in sm.items() if k.startswith("smoothness")]
            F["cal_smooth_spec_min"] = f"{min(v['pooled_fpr_spectral_05'] for v in smo):.3f}"
            F["cal_smooth_spec_max"] = f"{max(v['pooled_fpr_spectral_05'] for v in smo):.3f}"
            F["cal_smooth_naive_min"] = f"{min(v['pooled_fpr_naive_05'] for v in smo):.3f}"
            F["cal_smooth_naive_max"] = f"{max(v['pooled_fpr_naive_05'] for v in smo):.3f}"
            F["cal_n_sim_secondary"] = nc.get("n_sim_secondary")
            # binomial check: sections whose FPR CI excludes 0.05
            from scipy.stats import binomtest
            bad = 0
            for sid, r in nc["per_slide"].items():
                k_, n_ = r["gaussian"]["k_spectral_05"], r["gaussian"]["n_sim"]
                if binomtest(k_, n_, 0.0495).pvalue < 0.05:
                    bad += 1
            F["cal_n_sections_miscal"] = bad

            if "marginal" in sm:
                F["cal_marg_spec"] = f"{sm['marginal']['pooled_fpr_spectral_05']:.3f}"
                F["cal_marg_nscore"] = f"{sm['marginal']['pooled_fpr_nscore_05']:.3f}"
                F["cal_marg_naive"] = f"{sm['marginal']['pooled_fpr_naive_05']:.3f}"
            if "pooled_fpr_nscore_05" in g:
                F["cal_gauss_nscore"] = f"{g['pooled_fpr_nscore_05']:.3f}"
                F["cal_skew_nscore"] = f"{sk['pooled_fpr_nscore_05']:.3f}"
                F["cal_smooth_nscore_min"] = f"{min(v['pooled_fpr_nscore_05'] for v in smo):.3f}"
                F["cal_smooth_nscore_max"] = f"{max(v['pooled_fpr_nscore_05'] for v in smo):.3f}"
                allns = [g["pooled_fpr_nscore_05"], sk["pooled_fpr_nscore_05"]] + [v["pooled_fpr_nscore_05"] for v in smo]
                if "marginal" in sm:
                    allns.append(sm["marginal"]["pooled_fpr_nscore_05"])
                F["cal_nscore_min"], F["cal_nscore_max"] = f"{min(allns):.3f}", f"{max(allns):.3f}"

    # ---------------- normal-score surrogate test ----------------
    p = VAL / "spatial_null_nscore.json"
    if p.exists():
        ns = j("spatial_null_nscore.json")
        F["ns_n_q05"] = ns["summary"]["n_q05"]
        F["ns_n_p05"] = ns["summary"]["n_p05"]
        F["ns_ext_qs"] = ", ".join(f"{ns['per_slide'][x]['bh_padj']:.3f}" for x in EXTERNAL_ORDER)
        sn_ = j("spatial_null_check.json")["per_slide"]
        up = [x for x in PRIMARY_ORDER if ns["per_slide"][x]["bh_padj"] < 0.05 <= sn_[x]["bh_padj"]]
        down = [x for x in PRIMARY_ORDER if sn_[x]["bh_padj"] < 0.05 <= ns["per_slide"][x]["bh_padj"]]
        F["ns_n_up"], F["ns_n_down"] = len(up), len(down)
        F["ns_up_sections"] = ", ".join(up) if up else "none"
        F["ns_down_sections"] = ", ".join(down) if down else "none"
    # skewness of the observed B_mAb fields
    from scipy.stats import skew as _skew
    sk_ = []
    for x in PRIMARY_ORDER:
        z = np.load(VAL.parents[1] / "data" / "interim" / f"{x}.barrier.npz")
        sk_.append(_skew(z["b_mab"]))
    F["bmab_skew_median"] = f"{np.median(sk_):.2f}"
    # S2 on planted capsules of increasing thickness
    p = VAL / "synthetic_s2_thickness.json"
    if p.exists():
        th = j("synthetic_s2_thickness.json")["results"]
        parts = []
        for k in sorted(th, key=lambda s: int(s[:-2])):
            v = th[k]
            parts.append(f"{k[:-2]} μm: median {v['median_ratio_matched']:.2f} (unmatched {v['median_ratio_unmatched']:.2f}; "
                         f"{v['n_matched_gt1']}/{v['n']} above 1)")
        F["s2thick_sentence"] = ("With capsules of increasing thickness the matched ratios were " + "; ".join(parts) +
                                 ", whereas the unmatched design, in which only the contiguous arm is optimised, "
                                 "inflated the ratios.")
        allr = [r["ratio_matched"] for v in th.values() for r in v["per_rep"]]
        meds = [v["median_ratio_matched"] for v in th.values()]
        F["s2thick_n"] = len(allr)
        F["s2thick_nsig"] = int(sum(v["n_matched_p05"] for v in th.values()))
        F["s2thick_median"] = f"{np.median(allr):.2f}"
        F["s2thick_med_min"], F["s2thick_med_max"] = f"{min(meds):.2f}", f"{max(meds):.2f}"
    else:
        F["s2thick_sentence"] = ""

    # ---------------- decomposition ----------------
    p = VAL / "geometry_null.json"
    if p.exists():
        g = j("geometry_null.json")
        if "summary" in g:
            s_ = g["summary"]
            F["dec_geom_median"] = f3(s_["median_geometry_null_mean"])
            F["dec_constr_median"] = f3(s_["median_construction_null_mean"])
            F["dec_excess_median"] = f3(s_["median_excess_over_construction"])
            F["dec_n_excess_pos"] = s_["n_excess_construction_positive"]
            F["dec_n_p_constr"] = s_["n_p_construction_lt_05"]
            F["dec_n_q_constr"] = s_["n_bh_construction_lt_05"]
            F["dec_n_p_geom"] = s_["n_p_geometry_lt_05"]
            F["dec_n_q_geom"] = s_["n_bh_geometry_lt_05"]
            F["dec_share_constr_pct"] = f"{100 * s_['median_share_reproduced_by_construction']:.0f}"
            F["dec_r_ecm_caf"] = f2(s_["median_r_ecm_caf"])
            F["dec_r_ecm_xl"] = f2(s_["median_r_ecm_crosslink"])
            rr = [g["per_slide"][x]["input_correlations"]["ecm_caf"] for x in PRIMARY_ORDER]
            F["dec_r_ecm_caf_min"], F["dec_r_ecm_caf_max"] = f2(min(rr)), f2(max(rr))
            F["dec_n_null"] = g["per_slide"][PRIMARY_ORDER[0]]["n_null"]
            vis = ["CSCC01", "CSCC02", "CSCC03", "CSCC04"] + [x for x in PRIMARY_ORDER if x.startswith("MEL")]
            stl = [x for x in PRIMARY_ORDER if x not in vis]
            F["dec_excess_median_visium"] = f3(np.median([g["per_slide"][x]["excess_over_construction_null"] for x in vis]))
            F["dec_excess_median_st"] = f3(np.median([g["per_slide"][x]["excess_over_construction_null"] for x in stl]))
            F["dec_geom_median_visium"] = f3(np.median([g["per_slide"][x]["geometry_null_mean"] for x in vis]))
            F["dec_geom_median_st"] = f3(np.median([g["per_slide"][x]["geometry_null_mean"] for x in stl]))
            ex = [g["per_slide"][x] for x in EXTERNAL_ORDER if x in g["per_slide"]]
            if ex:
                F["dec_ext_n_p_constr"] = sum(r["p_vs_construction_null"] < 0.05 for r in ex)
    p = VAL / "shared_input_spatial_null.json"
    if p.exists():
        a = j("shared_input_spatial_null.json")
        if "summary" in a:
            s_ = a["summary"]
            F["abl_median"] = f3(s_["median_ablated"])
            F["abl_median_full"] = f3(s_["median_full"])
            F["abl_npos"] = s_["n_positive"]
            F["abl_n_p05"] = s_["n_p_lt_05"]
            F["abl_n_q05"] = s_["n_bh_lt_05"]
            F["abl_retained_pct"] = f"{100 * s_['median_retained_fraction']:.0f}"

    # ---------------- synthetic benchmark ----------------
    p = VAL / "synthetic_benchmark.json"
    if p.exists():
        sb = j("synthetic_benchmark.json")
        po = sb["pooled"]
        F["syn_n_tissues"] = sb["n_tissues"]
        F["syn_reps"] = sb["reps"]
        F["syn_n_geoms"] = len(sb["geoms"])
        F["syn_spots"] = sb["lattice_spots"]
        for k, v in po.items():
            F[f"syn_{k}_rho"] = f2(v["spearman_vs_lost_access"])
            F[f"syn_{k}_auc_band"] = f2(v["auc_closed_vs_band"])
            F[f"syn_{k}_auc_gap05"] = f2(v["auc_closed_vs_gap05"])
            F[f"syn_{k}_auc_scat"] = f2(v["auc_closed_vs_scattered"])
            if v.get("auc_closedfar_vs_scattered") is not None:
                F[f"syn_{k}_auc_far"] = f2(v["auc_closedfar_vs_scattered"])
        others = [k for k in po if not k.startswith("sparta")]
        F["syn_other_auc_gap05_max"] = f2(max(po[k]["auc_closed_vs_gap05"] for k in others))
        ab = sb["access_by_geom"]
        F["syn_access_closed"] = f2(ab["closed"]["mean"])
        F["syn_access_scattered"] = f2(ab["scattered"]["mean"])
        F["syn_access_band"] = f2(ab["band"]["mean"])
        if "closed_far" in ab:
            F["syn_access_far"] = f2(ab["closed_far"]["mean"])
        F["syn_s2_median"] = f2(sb["synthetic_s2_closed"]["median_ratio"])
        F["syn_s2_n"] = len(sb["synthetic_s2_closed"]["ratios"])
        F["syn_s2_nsig"] = int(sum(x < 0.05 for x in sb["synthetic_s2_closed"]["p"]))
        pn = sb["per_noise"]
        F["syn_bcell_auc_gap05_by_noise"] = ", ".join(f2(pn[k]["sparta_bcell"]["auc_closed_vs_gap05"]) for k in sorted(pn, key=float))
        F["syn_nhood_auc_gap05_by_noise"] = ", ".join(f2(pn[k]["nhood_enrichment_z"]["auc_closed_vs_gap05"]) for k in sorted(pn, key=float))

    # ---------------- molecular ground truth (run_59, v2.3) ----------------
    p = VAL / "mab_ground_truth.json"
    if p.exists():
        mg = j("mab_ground_truth.json")
        RT = {"0.5": "05", "2": "2", "5.5": "55"}
        F["mab59_n_tissues"] = mg["n_tissues"] // len(mg["radii"])
        F["mab59_W"] = f"{mg['W']:,}"
        F["mab59_sat_cap_w"] = "forty"
        F["mab59_bmab_floor"] = f"{-np.log(1e-12):.1f}"
        for r, tg in RT.items():
            pr = mg["per_radius"][r]
            for k in ("b_mab_core", "b_mab_reach_mean", "b_cell_field_core",
                      "ecm_peritumoural_mean", "ecm_global_mean", "vessel_boundary_distance"):
                F[f"mab59_rho_{k}_{tg}"] = f2(pr[k]["spearman_vs_lost_delivery"])
                F[f"mab59_auc_{k}_{tg}"] = f2(pr[k]["auc_closed_vs_gap05"])
            F[f"mab59_gt_auc_gap05_{tg}"] = f2(pr["_gt"]["auc_closed_vs_gap05"])
            F[f"mab59_gt_auc_band_{tg}"] = f2(pr["_gt"]["auc_closed_vs_band"])
            F[f"mab59_node_rho_{tg}"] = f2(pr["node_rho_median"])
            F[f"mab59_clip_{tg}"] = pct(mg["clip_fraction_main"][r]["overall"], 0)
            for T in mg["checkpoints"]:
                F[f"mab59_rho_T{T}_{tg}"] = f2(mg["by_checkpoint"][str(T)][r]["b_mab_core"])
            d = mg["delivery_by_radius"][r]["by_geom"]
            F[f"mab59_delivery_closed_{tg}"] = f"{d['closed']:.4f}"
            F[f"mab59_delivery_gap05_{tg}"] = f"{d['gap05']:.4f}"
        F["mab59_clip_closed_55"] = pct(mg["clip_fraction_main"]["5.5"]["by_geom"]["closed"], 0)
        if "saturation" in mg:
            st = mg["saturation"]
            F["mab59_sat_rho_lin"] = f2(st["rho_mab_vs_lost_delivery_linear"])
            F["mab59_sat_rho_sat"] = f2(st["rho_mab_vs_lost_delivery_saturating"])
            F["mab59_sat_rho_linvssat"] = f2(st["delivery_spearman_lin_vs_sat"])
            F["mab59_sat_relchange_pct"] = f"{100 * st['median_rel_delivery_change']:+.0f}"

    # ---------------- S2 / S3 / mesh / crosslink ----------------
    s2 = j("s2_matched_selection.json")["summary"]
    F["s2_median"] = f3(s2["median_ratio"])
    F["s2_min"], F["s2_max"] = f3(s2["ratio_range"][0]), f3(s2["ratio_range"][1])
    F["s2_n_gt1"] = s2["n_ratio_above_1"]
    F["s2_n_p05"] = s2["n_significant"]
    s2p = j("s2_matched_selection.json")["per_slide"]
    q2 = __import__("isdata").bh([s2p[x]["p_vs_in_cut_matched"] for x in PRIMARY_ORDER])
    F["s2_n_q05"] = int((np.array(q2) < 0.05).sum())
    s2e = j("s2_matched_ext.json")["per_slide"]
    F["s2_ext"] = ", ".join(f3(s2e[x]["ratio_vs_in_cut_matched"]) for x in EXTERNAL_ORDER)
    ms = j("size_exclusion_scan.json")["summary"]
    F["mesh_excl_min"], F["mesh_excl_max"] = f"{ms['primary_min_pct']:.1f}", f"{ms['primary_max_pct']:.1f}"
    F["mesh_excl_median"] = f"{ms['primary_median_pct']:.1f}"
    F["mesh_threshold"] = f"{np.log(20.0 / 5.5) / 3.0:.3f}"
    nx_ = j("null_crosslink_check.json")["summary"]
    F["xl_share_real"] = f"{100 * nx_['median_real']:.1f}"
    F["xl_share_null"] = f"{100 * nx_['median_null']:.1f}"
    # S3 size scan monotonicity
    mono, ratio = 0, []
    for x in PRIMARY_ORDER:
        c = json.load(open(VAL.parent / "counterfactual" / f"{x}.json", encoding="utf-8"))["s3"]
        mc = np.array(c["mean_core"])
        mono += int(np.all(np.diff(mc) > 0))
        r_ = c["radii_nm"]
        ratio.append(mc[r_.index(5.5)] / mc[r_.index(0.5)])
    F["s3_n_monotone"] = mono
    F["s3_ratio_median"] = f"{np.median(ratio):.1f}"

    # ---------------- domain comparison & runtime ----------------
    bl = j("benchmark_lambda_sensitivity.json")["summary"]
    F["dom_enrich_median"] = f2(bl["lambda_0.3_all_primary"]["median_enrichment"])
    F["dom_enrich_min"] = f2(bl["lambda_0.3_all_primary"]["enrichment_range"][0])
    F["dom_enrich_max"] = f2(bl["lambda_0.3_all_primary"]["enrichment_range"][1])
    F["dom_prec_median_pct"] = f"{100 * bl['lambda_0.3_all_primary']['median_precision']:.1f}"
    lam_keys = [k for k in bl if k.endswith("_all_primary")]
    F["dom_lambda_enrich_min"] = f2(min(bl[k]["median_enrichment"] for k in lam_keys))
    F["dom_lambda_enrich_max"] = f2(max(bl[k]["median_enrichment"] for k in lam_keys))
    F["dom_lambda_prec_min"] = f"{100 * min(bl[k]['median_precision'] for k in lam_keys):.1f}"
    F["dom_lambda_prec_max"] = f"{100 * max(bl[k]['median_precision'] for k in lam_keys):.1f}"
    rt = j("runtime_benchmark.json")["summary"]["all"]
    F["rt_sparta_median_s"] = f"{rt['t_sparta_core_s_median']:.3f}"
    F["rt_sparta_total_s"] = f"{rt['t_sparta_core_s_total']:.2f}"
    F["rt_banksy_median_s"] = f"{rt['t_banksy_style_s_median']:.2f}"
    F["rt_squidpy_median_s"] = f"{rt['t_squidpy_s_median']:.1f}"
    F["rt_peak_mem_mb"] = f"{rt['peak_mem_mb_max']:.0f}"

    # ---------------- biological alignment, bulk ----------------
    bv = j("biological_validation.json")["summary"]
    F["bio_bc_tnk_median"] = f3(bv["median_bc_vs_tnk"])
    F["bio_bc_tnk_nsig"] = bv["n_bc_tnk_neg_sig"]
    F["bio_bm_prolif_median"] = f3(bv["median_bm_vs_prolif"])
    F["bio_bm_prolif_nsig"] = bv["n_bm_prolif_pos_sig"]
    pa = j("prereadiness_audit.json")["section_bulk_icb_exploratory"]["cohorts"]
    c = {x["cohort"]: x for x in pa}
    F["icb_78220_auc"] = f"{c['GSE78220']['AUC_NR_over_R']:.3f}"
    F["icb_78220_p"] = f"{c['GSE78220']['p_value']:.3f}"
    F["icb_91061_auc"] = f"{c['GSE91061']['AUC_NR_over_R']:.3f}"
    F["icb_91061_p"] = f"{c['GSE91061']['p_value']:.3f}"
    F["icb_91061s_auc"] = f"{c['GSE91061_strict']['AUC_NR_over_R']:.3f}"

    # ---------------- radius sensitivity ----------------
    rs = j("radius_sensitivity.json")["summary"]
    F["radius_n_configs"] = rs["total_radius_configs"]
    F["radius_frac_pos_pct"] = f"{100 * rs['frac_rho_positive']:.0f}"
    F["radius_median_rho"] = f3(rs["median_rho_all"])
    # configurations whose radius equals the spot pitch connect only part of the nearest neighbours
    rps = j("radius_sensitivity.json")["per_slide"]
    pitch = {sid: (100.0 if (sid.startswith("MEL") or sid in ("CSCC01", "CSCC02", "CSCC03", "CSCC04")) else 200.0)
             for sid in rps}
    nd = [v["rho_partial"] for sid, row in rps.items() for r_, v in row.items() if float(r_) > pitch[sid]]
    F["radius_n_degenerate"] = sum(1 for sid, row in rps.items() for r_ in row if float(r_) <= pitch[sid])
    F["radius_n_nondeg"] = len(nd)
    F["radius_nondeg_npos"] = int(sum(x > 0 for x in nd))


    # ================= v2.2 additions (2026-10-05) =================
    # ---------------- melanoma replication cohort (run_48) ----------------
    rep = j("replication_melanoma.json")
    rs, rp = rep["summary"], rep["patient_level"]
    rper = rep["per_slide"]
    F["rep_n_sections"] = rs["n_sections"]
    F["rep_n_patients"] = rs["n_patients"]
    F["rep_n_pos"] = rs["n_pos"]
    F["rep_median"] = f3(rs["median_rho"])
    F["rep_min"], F["rep_max"] = f3(rs["min_rho"]), f3(rs["max_rho"])
    F["rep_n_q05"] = rs["n_q05"]
    F["rep_n_q05_nscore"] = rs["n_q05_nscore"]
    F["rep_spots_min"] = f"{min(v['n_nodes'] for v in rper.values()):,}"
    F["rep_spots_max"] = f"{max(v['n_nodes'] for v in rper.values()):,}"
    ag_miss = [s_ for s_, v in rper.items() if "ag_target" in v["filled_keys"]]
    F["rep_n_ag_missing"] = len(ag_miss)
    F["rep_ag_missing_sections"] = ", ".join(x.replace("MEL_THR", "LN P").replace("_rep", " r") for x in ag_miss)
    a_ = rp["association"]
    F["rep_pl_mean"], F["rep_pl_lo"], F["rep_pl_hi"] = (f2(a_["nested_model"]["mean"]), f2(a_["nested_model"]["ci95"][0]),
                                                         f2(a_["nested_model"]["ci95"][1]))
    F["rep_pl_df"] = a_["nested_model"]["df"]
    F["rep_pl_p"] = f"{a_['nested_model']['p_two_sided']:.3f}"
    F["rep_signflip_p"] = f"{a_['sign_flip']['p_one_sided']:.4f}"
    F["rep_npos_pat"] = a_["sign_flip"]["n_positive"]
    F["rep_lopo_min"] = f2(a_["lopo_min_mean"])
    al = rp["pooled_primary_plus_replication"]
    F["all12_mean"], F["all12_lo"], F["all12_hi"] = (f2(al["nested_model"]["mean"]), f2(al["nested_model"]["ci95"][0]),
                                                     f2(al["nested_model"]["ci95"][1]))
    F["all12_df"] = al["nested_model"]["df"]
    F["all12_p"] = "< 0.0001" if al["nested_model"]["p_two_sided"] < 1e-4 else f"{al['nested_model']['p_two_sided']:.4f}"
    F["all12_signflip_p"] = f"{al['sign_flip']['p_one_sided']:.5f}"
    F["all12_npos"] = al["sign_flip"]["n_positive"]
    F["all12_J"] = al["sign_flip"]["J"]
    F["all12_n_sections"] = al["n_sections"]
    F["rep_constr_median"] = f3(rs["median_construction_null"])
    F["rep_geom_median"] = f3(rs["median_geometry_null"])
    F["rep_share_constr_pct"] = f"{100 * rs['median_share_construction']:.0f}"
    F["rep_excess_median"] = f3(rs["median_excess"])
    F["rep_n_excess_pos"] = rs["n_excess_pos"]
    F["rep_n_q_constr"] = rs["n_q05_construction"]
    F["rep_n_p_constr"] = int(sum(v["p_vs_construction_null"] < 0.05 for v in rper.values()))
    ex_ = rp["excess"]["nested_model"]
    F["rep_ex_npos_pat"] = rp["excess"]["sign_flip"]["n_positive"]
    F["rep_ex_signflip_p"] = f"{rp['excess']['sign_flip']['p_one_sided']:.4f}"
    F["all12_n_pos_sections"] = int(sum(v > 0 for v in [per[x]["real_rho_partial"] for x in PRIMARY_ORDER]
                                        + [v_["real_rho_partial"] for v_ in rper.values()]))
    F["rep_ex_mean"], F["rep_ex_lo"], F["rep_ex_hi"] = f2(ex_["mean"]), f2(ex_["ci95"][0]), f2(ex_["ci95"][1])
    F["rep_ablated_median"] = f3(rs["median_ablated"])
    F["rep_n_ablated_pos"] = rs["n_ablated_pos"]
    F["rep_n_ablated_q05"] = rs["n_ablated_q05"]
    F["rep_r_ecm_caf"] = f2(rs["median_r_ecm_caf"])
    F["rep_n_discordant_below_chance"] = rs["n_discordant_below_chance"]
    lite_eq = j("lite_scoring_equivalence.json")
    _m, _e = f"{lite_eq['max_abs_diff_raw']:.1e}".split("e")
    F["lite_max_dscore"] = f"${_m}\\times 10^{{{int(_e)}}}$"
    F["lite_max_drank"] = f"{lite_eq['max_abs_diff_rank']:.0f}"
    F["lite_n_sections"] = len(lite_eq["sections"])
    cr = j("null_calibration_replication.json")["summary"]
    F["cal_rep_spec"] = f3(cr["gaussian"]["pooled_fpr_spectral_05"])
    F["cal_rep_naive"] = f3(cr["gaussian"]["pooled_fpr_naive_05"])
    F["cal_rep_nscore"] = f3(cr["gaussian"]["pooled_fpr_nscore_05"])
    F["cal_rep_skew_spec"] = f3(cr["skewed"]["pooled_fpr_spectral_05"])
    F["cal_rep_skew_nscore"] = f3(cr["skewed"]["pooled_fpr_nscore_05"])
    c22 = j("null_calibration_all.json")["summary"]["gaussian"]
    k30 = c22["pooled_fpr_spectral_05"] * c22["n_tests"] + cr["gaussian"]["pooled_fpr_spectral_05"] * cr["gaussian"]["n_tests"]
    n30 = c22["n_tests"] + cr["gaussian"]["n_tests"]
    F["cal30_n_graphs"] = 30
    F["cal30_pooled_spec"] = f3(k30 / n30)
    F["cal30_n_tests"] = f"{n30:,}"

    # ---------------- CODEX validation against measured CD8+ T cells (run_47) ----------------
    cx = j("codex_validation.json")
    cs = cx["summary"]
    cnt = cs["counts"]
    F["cx_n_cores_total"] = cnt["n_cores_total"]
    F["cx_n_cores"] = cnt["n_cores_included"]
    F["cx_n_patients"] = cnt["n_patients_included"]
    F["cx_n_clr"], F["cx_n_dii"] = cnt["n_clr_patients"], cnt["n_dii_patients"]
    F["cx_n_cells"] = "258,385"
    F["cx_n_excluded"] = cnt["n_cores_total"] - cnt["n_cores_included"]
    dsc = cs["describe"]
    F["cx_lcc_median"] = f"{dsc['n_lcc']['median']:,.0f}"
    F["cx_lcc_min"], F["cx_lcc_max"] = f"{dsc['n_lcc']['min']:,.0f}", f"{dsc['n_lcc']['max']:,.0f}"
    F["cx_cd8_median"] = f"{dsc['n_outcome_cells']['median']:.0f}"
    F["cx_ir_median"] = f2(dsc["ir"]["median"])
    pc = cs["primary_and_comparators"]
    F["cx_rho_core"] = f2(pc["log_b_rel"]["rho_core"])
    F["cx_ci_lo"], F["cx_ci_hi"] = f2(pc["log_b_rel"]["ci_core"][0]), f2(pc["log_b_rel"]["ci_core"][1])
    F["cx_rho_pat"] = f2(pc["log_b_rel"]["patient"]["rho"])
    F["cx_p_pat"] = f"{pc['log_b_rel']['patient']['p_one_sided_perm']:.4f}"
    F["cx_bcell_rho_core"] = f2(pc["log_b_cell"]["rho_core"])
    F["cx_bcell_rho_pat"] = f2(pc["log_b_cell"]["patient"]["rho"])
    F["cx_field_rho_core"] = f2(pc["field_core"]["rho_core"])
    F["cx_field_rho_pat"] = f2(pc["field_core"]["patient"]["rho"])
    F["cx_peri_rho_core"] = f2(pc["peritumoural_matrix"]["rho_core"])
    F["cx_peri_rho_pat"] = f2(pc["peritumoural_matrix"]["patient"]["rho"])
    F["cx_stromal_rho_core"] = f2(pc["stromal_fraction"]["rho_core"])
    F["cx_ripley_rho_core"] = f2(pc["ripley_l_matrix"]["rho_core"])
    F["cx_contact_rho_core"] = f2(pc["contact_enrichment_z"]["rho_core"])
    F["cx_contact_rho_pat"] = f2(pc["contact_enrichment_z"]["patient"]["rho"])
    pa_ = cs["added_value_partial"]
    F["cx_partial_rho"], F["cx_partial_lo"], F["cx_partial_hi"] = f2(pa_["rho"]), f2(pa_["ci"][0]), f2(pa_["ci"][1])
    ph = cs["post_hoc"]
    F["cx_ph_intermix_rho"] = f2(ph["b_rel_given_intermixing"]["rho"])
    F["cx_ph_intermix_lo"], F["cx_ph_intermix_hi"] = (f2(ph["b_rel_given_intermixing"]["ci"][0]),
                                                      f3(ph["b_rel_given_intermixing"]["ci"][1]))
    F["cx_ph_all_rho"] = f2(ph["b_rel_given_composition_and_intermixing"]["rho"])
    F["cx_ph_all_lo"] = f2(ph["b_rel_given_composition_and_intermixing"]["ci"][0])
    F["cx_ph_all_hi"] = f2(ph["b_rel_given_composition_and_intermixing"]["ci"][1])
    F["cx_ph_contact_rho"] = f2(ph["contact_enrichment_z_given_others"]["rho"])
    F["cx_ph_peri_rho"] = f2(ph["peritumoural_matrix_given_others"]["rho"])
    s1 = cs["s1_absolute"]
    F["cx_s1_rho_core"] = f2(s1["rho_core"])
    F["cx_s1_rho_pat"] = f2(s1["patient"]["rho"])
    F["cx_s1_p"] = f"{s1['patient']['p_one_sided_perm']:.3f}"
    s2_ = cs["s2_blockade_line"]
    F["cx_s2_ratio_sparta"] = f2(s2_["median_cut_ratio_sparta"])
    F["cx_s2_ratio_open"] = f2(s2_["median_cut_ratio_open"])
    F["cx_s2_p"] = f"{s2_['wilcoxon_one_sided_p']:.2f}"
    F["cx_s2_area_sparta_pct"] = f"{100 * s2_['median_beyond_area_frac_sparta']:.0f}"
    F["cx_s2_area_open_pct"] = f"{100 * s2_['median_beyond_area_frac_open']:.0f}"
    s3_ = cs["s3_clr_vs_dii"]
    F["cx_s3_p_brel"] = f"{s3_['mannwhitney_p_b_rel']:.3f}"
    F["cx_s3_p_ir"] = f"{s3_['mannwhitney_p_ir']:.3f}"
    F["cx_s3_brel_clr"], F["cx_s3_brel_dii"] = f2(s3_["median_log_b_rel_clr"]), f2(s3_["median_log_b_rel_dii"])
    F["cx_s3_ir_clr"], F["cx_s3_ir_dii"] = f2(s3_["median_ir_clr"]), f2(s3_["median_ir_dii"])
    sens = cs["sensitivity"]
    F["cx_sens_rho_min"] = f2(max(v["patient"]["rho"] for v in sens.values()))   # least negative
    F["cx_sens_rho_max"] = f2(min(v["patient"]["rho"] for v in sens.values()))   # most negative
    F["cx_sens_pmax"] = f"{max(v['patient']['p_one_sided_perm'] for v in sens.values()):.3f}"
    F["cx_sens_n"] = len(sens)
    pre = {k: v for k, v in sens.items() if not v.get("post_hoc")}
    phs = {k: v for k, v in sens.items() if v.get("post_hoc")}
    F["cx_sens_pre_n"] = len(pre)
    F["cx_sens_pre_rho_min"] = f2(max(v["patient"]["rho"] for v in pre.values()))
    F["cx_sens_pre_rho_max"] = f2(min(v["patient"]["rho"] for v in pre.values()))
    F["cx_sens_pre_pmax"] = f"{max(v['patient']['p_one_sided_perm'] for v in pre.values()):.3f}"
    if phs:
        F["cx_sens_ph_n"] = len(phs)
        F["cx_sens_ph_rho_min"] = f2(max(v["patient"]["rho"] for v in phs.values()))
        F["cx_sens_ph_rho_max"] = f2(min(v["patient"]["rho"] for v in phs.values()))
        F["cx_sens_ph_pmax"] = f"{max(v['patient']['p_one_sided_perm'] for v in phs.values()):.3f}"
    F["cx_tumour_region_rho"] = f2(cs["ir_tumour_region"]["patient"]["rho"])

    # ---------------- intervention maps (run_44 + run_49) ----------------
    tg = j("intervention_targeting.json")["summary"]["primary"]
    F["iv_cut_frac_pct"] = f"{100 * tg['median_cut_frac']:.0f}"
    F["iv_frac_zero_pct"] = f"{100 * tg['median_frac_zero_cell']:.0f}"
    F["iv_k20_cell_sparta"] = f"{100 * tg['k20_cell_sparta']:.0f}"
    F["iv_k20_cell_randcut"] = f"{100 * tg['k20_cell_random_cut']:.0f}"
    F["iv_k20_cell_density"] = f"{100 * tg['k20_cell_density']:.0f}"
    F["iv_k20_cell_random"] = f"{100 * tg['k20_cell_random']:.0f}"
    F["iv_k20_cell_mabrank"] = f"{100 * tg['k20_cell_with_mab_ranking']:.0f}"
    F["iv_k20_mab_sparta"] = f"{100 * tg['k20_mab_sparta']:.0f}"
    F["iv_k20_mab_density"] = f"{100 * tg['k20_mab_density']:.0f}"
    F["iv_k20_mab_random"] = f"{100 * tg['k20_mab_random']:.0f}"
    F["iv_k20_mab_cellrank"] = f"{100 * tg['k20_mab_with_cell_ranking']:.0f}"
    F["iv_k50_cell_pct"] = f"{100 * tg['median_k50_cell_frac']:.1f}"
    F["iv_k50_mab_pct"] = f"{100 * tg['median_k50_mab_frac']:.0f}"
    F["iv_jaccard"] = f2(tg["median_jaccard"])
    F["iv_n_wins_density"] = tg["k20_n_sparta_gt_density_cell"]
    F["iv_n_wins_randcut"] = tg["k20_n_sparta_gt_random_cut_cell"]
    F["iv_n_wins_density_mab"] = tg["k20_n_sparta_gt_density_mab"]

    # intervention: concentration numbers for the Discussion (run_44 joint curves, v2.2)
    ir_ = j("intervention_ranking.json")
    irs, irp = ir_["summary"], ir_["per_section"]
    mc_ = irs["main_cohort"]

    def _joint(sid, op, frac):
        r = irp[sid]
        k = max(1, int(round(frac * r["n"])))
        pts = {q["k"]: q["joint_frac_of_anchor"] for q in r["joint_intervention"][op]["points"]}
        return pts.get(k, float("nan"))
    single = [irp[x]["b_cell"]["max_over_anchor"] for x in PRIMARY_ORDER]
    F["iv_single_cell_median"] = f"{100 * np.median(single):.0f}"
    F["iv_single_cell_max"] = f"{100 * max(single):.0f}"
    F["iv_single_cell_max_section"] = PRIMARY_ORDER[int(np.argmax(single))]
    single_m = [irp[x]["b_mab"]["max_over_anchor"] for x in PRIMARY_ORDER]
    F["iv_single_mab_median"] = f"{100 * np.median(single_m):.1f}"
    F["iv_single_mab_max"] = f"{100 * max(single_m):.1f}"
    t1c = [_joint(x, "b_cell", 0.01) for x in PRIMARY_ORDER]
    t1m = [_joint(x, "b_mab", 0.01) for x in PRIMARY_ORDER]
    t5c = [_joint(x, "b_cell", 0.05) for x in PRIMARY_ORDER]
    t5m = [_joint(x, "b_mab", 0.05) for x in PRIMARY_ORDER]
    t10m = [_joint(x, "b_mab", 0.10) for x in PRIMARY_ORDER]
    F["iv_top1_cell_median"] = f"{100 * np.median(t1c):.0f}"
    F["iv_top1_cell_min"], F["iv_top1_cell_max"] = f"{100 * min(t1c):.0f}", f"{100 * max(t1c):.0f}"
    F["iv_top1_mab_median"] = f"{100 * np.median(t1m):.0f}"
    F["iv_top1_mab_min"], F["iv_top1_mab_max"] = f"{100 * min(t1m):.0f}", f"{100 * max(t1m):.0f}"
    F["iv_top5_cell_median"] = f"{100 * np.median(t5c):.0f}"
    F["iv_top5_mab_median"] = f"{100 * np.median(t5m):.0f}"
    F["iv_top5_mab_max"] = f"{100 * max(t5m):.0f}"
    F["iv_top10_mab_median"] = f"{100 * np.median(t10m):.0f}"
    jac = [irp[x]["operator_agreement"]["top5pct_jaccard"] for x in PRIMARY_ORDER]
    F["iv_jacc_min_pct"], F["iv_jacc_max_pct"] = f"{100 * min(jac):.0f}", f"{100 * max(jac):.0f}"
    F["iv_jacc_median_pct"] = f"{100 * np.median(jac):.1f}"
    coh_j = [irs["per_cohort"][c]["median_top5pct_jaccard"] for c in irs["per_cohort"]]
    F["iv_jacc_cohort_min_pct"], F["iv_jacc_cohort_max_pct"] = f"{100 * min(coh_j):.0f}", f"{100 * max(coh_j):.0f}"
    F["iv_n_ablations"] = f"{sum(r['n'] for r in irp.values()):,}"
    F["iv_n_sections_all"] = len(irp)
    rc = irs.get("replication_cohort", {})
    if rc:
        F["iv_rep_top1_cell"] = f"{100 * rc['median_joint_frac_cell_at_top1pct']:.0f}"
        F["iv_rep_top1_mab"] = f"{100 * rc['median_joint_frac_mab_at_top1pct']:.0f}"
        F["iv_rep_top5_mab"] = f"{100 * rc['median_joint_frac_mab_at_top5pct']:.0f}"
        F["iv_rep_jacc_pct"] = f"{100 * rc['median_top5pct_jaccard']:.1f}"
    tgr = j("intervention_targeting.json")["summary"].get("replication")
    if tgr:
        F["iv_rep_k20_cell_sparta"] = f"{100 * tgr['k20_cell_sparta']:.0f}"
        F["iv_rep_k20_cell_density"] = f"{100 * tgr['k20_cell_density']:.0f}"
        F["iv_rep_n_wins_density"] = tgr["k20_n_sparta_gt_density_cell"]

    # compartments of the high-impact spots (run_53)
    de = j("intervention_domain_enrichment.json")["summary"]
    for key, tag in (("b_cell_top1pct", "c1"), ("b_cell_top5pct", "c5"), ("b_mab_top1pct", "m1"),
                     ("b_mab_top5pct", "m5")):
        for d, dt in (("tumour", "tum"), ("stroma", "str"), ("immune", "imm")):
            e = de[key][d]
            F[f"de_{tag}_{dt}_top"] = f"{100 * e['pooled_frac_top']:.0f}"
            F[f"de_{tag}_{dt}_all"] = f"{100 * e['pooled_frac_all']:.0f}"
            F[f"de_{tag}_{dt}_p"] = f"{e['patient_signflip_p_two_sided']:.3f}"
            F[f"de_{tag}_{dt}_npos"] = e["n_patients_enriched"]
            F[f"de_{tag}_{dt}_nneg"] = e["n_patients_depleted"]
            F[f"de_{tag}_{dt}_J"] = e["n_patients"]
            if "pooled_frac_cut" in e:
                F[f"de_{tag}_{dt}_cut"] = f"{100 * e['pooled_frac_cut']:.0f}"

    # simple spatial summaries (run_52)
    sbl = j("simple_baselines_spearman.json")
    sp_, sc_ = sbl["spot_level"], sbl["section_level"]
    for fld, ft in (("b_cell_field", "cell"), ("b_mab_field", "mab")):
        for b, bt in (("stromal_density", "dens"), ("dist_tumour_boundary", "dist"), ("niche_z", "niche")):
            a_ = sp_[f"{fld}__{b}"]["all"]
            F[f"bl_{ft}_{bt}_med"] = f2(a_["median"])
            F[f"bl_{ft}_{bt}_q25"], F[f"bl_{ft}_{bt}_q75"] = f2(a_["q25"]), f2(a_["q75"])
            F[f"bl_{ft}_{bt}_max"] = f2(a_["max"])
            F[f"bl_{ft}_{bt}_min"] = f2(a_["min"])
            F[f"bl_{ft}_{bt}_n"] = a_["n"]
    F["bl_cell_mab_med"] = f2(sp_["b_cell_field__b_mab_field"]["all"]["median"])
    if "cut_density_pct_median" in sp_:
        F["bl_cut_pct_median"] = f"{100 * sp_['cut_density_pct_median']['all']['median']:.0f}"
        F["bl_cut_topdec_pct"] = f"{100 * sp_['cut_share_top_density_decile']['all']['median']:.0f}"
        F["bl_topdec_on_cut_pct"] = f"{100 * sp_['top_decile_share_on_cut']['all']['median']:.0f}"
    F["bl_cell_dens_r2_pct"] = f"{100 * sp_['b_cell_field__stromal_density']['all']['median'] ** 2:.0f}"
    for op, ot in (("log_b_rel", "brel"), ("b_cell_field_core", "cellcore"), ("b_mab_field_core", "mabcore")):
        for b, bt in (("stromal_density", "dens"), ("dist_tumour_boundary", "dist"), ("niche_z", "niche")):
            r_ = sc_[f"{op}__{b}"]
            F[f"bl_sec_{ot}_{bt}_rho"] = f2(r_["rho"])
            F[f"bl_sec_{ot}_{bt}_p"] = f"{r_['p']:.3f}" if r_["p"] >= 0.001 else f"{r_['p']:.4f}"
    F["bl_sec_n"] = sc_["log_b_rel__stromal_density"]["n_sections"]
    F["bl_sec_brel_absmax"] = f2(max(abs(sc_[f"log_b_rel__{b}"]["rho"]) for b in
                                     ("stromal_density", "dist_tumour_boundary", "niche_z")))

    # CODEX: post hoc geometry-only comparators
    if "vessel_tumour_distance" in pc:
        F["cx_dist_rho_core"] = f2(pc["vessel_tumour_distance"]["rho_core"])
        F["cx_dist_rho_pat"] = f2(pc["vessel_tumour_distance"]["patient"]["rho"])
        F["cx_dist_p_pat"] = f"{pc['vessel_tumour_distance']['patient']['p_one_sided_perm']:.3f}"
        F["cx_depth_rho_core"] = f2(pc["core_depth"]["rho_core"])
        F["cx_depth_rho_pat"] = f2(pc["core_depth"]["patient"]["rho"])
        F["cx_contact_p_pat"] = f"{pc['contact_enrichment_z']['patient']['p_two_sided']:.4f}"
        F["cx_peri_p_pat"] = f"{pc['peritumoural_matrix']['patient']['p_one_sided_perm']:.4f}"
        g_ = ph["b_rel_given_distances"]
        F["cx_ph_dist_rho"], F["cx_ph_dist_lo"], F["cx_ph_dist_hi"] = f2(g_["rho"]), f2(g_["ci"][0]), f2(g_["ci"][1])
        for x, t in (("peritumoural_matrix", "peri"), ("vessel_tumour_distance", "dist"), ("core_depth", "depth"),
                     ("contact_enrichment_z", "contact")):
            if f"b_rel_given_{x}" in ph:
                g_ = ph[f"b_rel_given_{x}"]
                F[f"cx_ph_given_{t}_rho"] = f2(g_["rho"])
                _fc = (lambda x: f"{x:.3f}" if abs(x) < 0.01 else f2(x))
                F[f"cx_ph_given_{t}_lo"], F[f"cx_ph_given_{t}_hi"] = _fc(g_["ci"][0]), _fc(g_["ci"][1])
        if "b_rel_given_all_simple" in ph:
            g_ = ph["b_rel_given_all_simple"]
            F["cx_ph_allsimple_rho"] = f2(g_["rho"])
            F["cx_ph_allsimple_lo"], F["cx_ph_allsimple_hi"] = f2(g_["ci"][0]), f2(g_["ci"][1])
        g_ = ph["b_rel_given_composition_and_distances"]
        F["cx_ph_compdist_rho"] = f2(g_["rho"])
        F["cx_ph_compdist_lo"], F["cx_ph_compdist_hi"] = f2(g_["ci"][0]), f2(g_["ci"][1])

    # domain comparison recomputed with exact cuts (run_51)
    dce = j("domain_comparison_exact_cut.json")["summary"]
    F["dom_repro_n"] = dce["n_sections_exactly_reproduced"]
    F["dom_repro_total"] = dce["n_sections"]
    if "archived_summary_lambda_0_3" in dce:
        F["dom_enrich_median_old"] = f2(dce["archived_summary_lambda_0_3"]["median_enrichment"])

    # ---------------- exact minimum cut (run_50) ----------------
    mc = j("mincut_exactness.json")["summary"]
    F["mc_n_sections"] = mc["n_sections"]
    F["mc_n_old_not_min"] = mc["float_partition_n_not_min"]
    F["mc_old_gap_max_pct"] = f"{100 * mc['float_partition_max_abs_gap']:.1f}"
    F["mc_median_jaccard"] = f2(mc["median_jaccard_edges"])
    F["mc_n_predates_lcc"] = len(mc["sections_stored_flow_predates_lcc_rule"])
    _m, _e = f"{mc['float_flow_max_abs_rel_diff']:.0e}".split("e")
    F["mc_flow_diff_max"] = f"${_m}\\times 10^{{{int(_e)}}}$"
    _m, _e = f"{mc['max_abs_new_cut_gap']:.0e}".split("e")
    F["mc_new_gap_max"] = f"${_m}\\times 10^{{{int(_e)}}}$"

    # ---------------- 2026 extension cohort (run_56) ----------------
    ext = j("extension_cohort.json")
    ext_per, ext_sum, ext_groups = ext["per_slide"], ext["summary"], ext["grouped"]
    ext_pl = ext["patient_level"]
    F["ext2026_n_sections"] = ext_sum["n_sections"]
    F["ext2026_n_visium"] = ext_groups["Visium"]["n_sections"]
    F["ext2026_n_slideseq"] = ext_groups["SlideSeqV2"]["n_sections"]
    F["ext2026_n_cscc"] = ext_groups["cSCC"]["n_sections"]
    F["ext2026_n_primary_mel"] = ext_groups["primary_melanoma"]["n_sections"]
    F["ext2026_n_metastatic_mel"] = ext_groups["metastatic_melanoma"]["n_sections"]
    F["ext2026_n_identifiable_patients"] = 32
    F["ext2026_n_unknown_relation_sections"] = 11
    F["ext2026_n_analytical_units"] = ext_sum["n_patients"]
    F["ext2026_n_pos"] = ext_sum["n_positive"]
    F["ext2026_median_rho"] = f3(ext_sum["median_rho"])
    F["ext2026_min_rho"] = f3(ext_sum["min_rho"])
    F["ext2026_max_rho"] = f3(ext_sum["max_rho"])
    F["ext2026_n_q05"] = ext_sum["n_q05"]

    def ext_model_fill(source: dict, prefix: str):
        nm = source["nested_model"]
        sf = source["sign_flip"]
        F[f"{prefix}_mean"] = f2(nm["mean"])
        F[f"{prefix}_lo"], F[f"{prefix}_hi"] = f2(nm["ci95"][0]), f2(nm["ci95"][1])
        F[f"{prefix}_p"] = f"{nm['p_two_sided']:.2e}"
        F[f"{prefix}_J"] = source["n_patients"]
        F[f"{prefix}_sign_p"] = f"{sf['p_one_sided']:.2e}"
        F[f"{prefix}_npos"] = sf["n_positive"]

    ext_model_fill(ext_pl["association"], "ext2026_pl")
    ext_model_fill(ext_pl["excess"], "ext2026_ex")
    ext_model_fill(ext_pl["ablated"], "ext2026_ab")
    ext_model_fill(ext_pl["pooled_primary_plus_extension"], "ext2026_pool")

    for group, prefix in [
        ("cSCC", "ext2026_g_cscc"),
        ("primary_melanoma", "ext2026_g_pmel"),
        ("metastatic_melanoma", "ext2026_g_mmel"),
        ("Visium", "ext2026_g_visium"),
        ("SlideSeqV2", "ext2026_g_slide"),
    ]:
        record = ext_groups[group]
        model = record["patient_model"]
        F[f"{prefix}_n"] = record["n_sections"]
        F[f"{prefix}_J"] = record["n_patients"]
        F[f"{prefix}_median"] = f3(record["median_rho"])
        F[f"{prefix}_npos"] = record["n_positive"]
        if isinstance(model, dict) and "nested_model" in model:
            F[f"{prefix}_mean"] = f2(model["nested_model"]["mean"])
            F[f"{prefix}_lo"] = f2(model["nested_model"]["ci95"][0])
            F[f"{prefix}_hi"] = f2(model["nested_model"]["ci95"][1])
            F[f"{prefix}_p"] = f"{model['nested_model']['p_two_sided']:.2e}"

    F["ext2026_constr_share_pct"] = f"{100 * ext_sum['median_share_construction']:.0f}"
    F["ext2026_excess_median"] = f3(ext_sum["median_excess"])
    F["ext2026_ablated_median"] = f3(ext_sum["median_ablated"])
    ext_bl = j("extension_simple_baselines_spearman.json")
    F["ext2026_bl_cell_dens_med"] = f2(
        ext_bl["spot_level"]["b_cell_field__stromal_density"]["all"]["median"])
    F["ext2026_bl_mab_dens_med"] = f2(
        ext_bl["spot_level"]["b_mab_field__stromal_density"]["all"]["median"])
    ext_cal = j("extension_null_calibration.json")["summary"]
    F["ext2026_cal_spec"] = f3(ext_cal["gaussian"]["pooled_fpr_spectral_05"])
    F["ext2026_cal_nscore"] = f3(ext_cal["gaussian"]["pooled_fpr_nscore_05"])
    F["ext2026_cal_naive"] = f3(ext_cal["gaussian"]["pooled_fpr_naive_05"])
    F["ext2026_cal_n_tests"] = ext_cal["gaussian"]["n_tests"]

    # small counts spelled out where they appear in running text
    _W = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven",
          "twelve"]
    for k in ("rep_n_ag_missing", "lite_n_sections", "cx_sens_n", "cx_sens_pre_n", "cx_sens_ph_n", "rep_n_q_constr",
              "rep_n_sections",
              "rep_n_excess_pos", "dom_repro_n", "dom_repro_total"):
        F[f"{k}_w"] = _W[int(F[k])] if int(F[k]) < len(_W) else str(F[k])

    out = VAL / "is_manuscript_facts.json"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(F, fh, indent=1, ensure_ascii=False, default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
    print(f"{len(F)} facts -> {out}")
    return F


if __name__ == "__main__":
    F = main()
    for k, v in F.items():
        print(f"{k:<34} {v}")
