#!/usr/bin/env python
"""
run_36_patient_level.py -- patient-level (hierarchical) inference for the section-level statistics
=================================================================================================
Inputs (read-only; nothing is recomputed from expression data):
    results/validation/spatial_null_check.json            primary field association, 19 sections
    results/validation/metric_connectivity_sensitivity.json largest-component re-estimate
    results/validation/s2_matched_selection.json           selection-matched S2 ratios
    results/validation/shared_input_spatial_null.json      (optional) ECM-ablated association
    results/validation/geometry_null.json                  (optional) construction-null excess
    data/ledger.csv                                        section -> patient map
Output:
    results/validation/patient_level_inference.json

Why (CMPB readiness review, item 4; IS submission)
--------------------------------------------------
Nineteen sections come from eight people (four melanoma sections from two people).
Section-level BH tests answer "is there spatial alignment in this section?", not
"is the association positive across patients?".  This script answers the second
question three ways, from weakest to strongest assumptions:

  1. Exact patient-level sign-flip test on patient means (2^J relabellings).
  2. Nested random-effects model (sections within patients) on the partial rho,
     using the graph-spectral surrogate SD of each section as its known sampling
     SD; between-patient and between-section variances by REML; inference on the
     pooled mean with a t reference on J-1 degrees of freedom (conservative).
  3. Leave-one-patient-out refits.

All three use the patient, never the section, as the unit of replication.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import optimize, stats  # noqa: E402

from sparta.io_ import Paths, load_config, load_json, patient_map, save_json, stamp_run  # noqa: E402


# --------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------
def sign_flip_exact(patient_means: np.ndarray) -> dict:
    """Exact sign-flip test of patient-level effects.

    For small J, enumerate all 2^J magnitude-weighted sign assignments. For
    larger J, use the exact binomial distribution of the number of positive
    patient effects; this avoids materialising 2^J assignments.
    """
    m = np.asarray(patient_means, float)
    J = len(m)
    obs = m.mean()
    n_positive = int((m > 0).sum())
    if J <= 20:
        flips = np.array(list(itertools.product([-1.0, 1.0], repeat=J)))
        null = (flips * np.abs(m)).mean(axis=1)
        p_one = float(np.mean(null >= obs - 1e-15))
        p_two = float(np.mean(np.abs(null) >= abs(obs) - 1e-15))
        method = "exact enumeration of magnitude-weighted sign flips"
    else:
        p_one = float(stats.binomtest(n_positive, J, alternative="greater").pvalue)
        p_two = float(stats.binomtest(n_positive, J).pvalue)
        method = "exact binomial sign test on number of positive patient effects"
    return dict(J=J, mean_of_patient_means=float(obs), n_positive=n_positive,
                p_one_sided=p_one, p_two_sided=p_two, method=method,
                min_attainable_p=float(2.0 ** -J))


def _nested_reml(y, groups, v_known, X=None, fit_section_var=True):
    """REML fit of y_i = X b + u_g(i) + w_i + e_i.

    u ~ N(0, tau_p^2) per patient, w ~ N(0, tau_s^2) per section (optional),
    e ~ N(0, v_i) with v_i known (pass zeros and fit_section_var=True to get an
    ordinary random-intercept LMM with an estimated residual variance).
    Returns dict with beta, its covariance, variance components and log-lik.
    """
    y = np.asarray(y, float)
    n = len(y)
    X = np.ones((n, 1)) if X is None else np.asarray(X, float)
    labels, g = np.unique(groups, return_inverse=True)
    Z = np.zeros((n, len(labels)))
    Z[np.arange(n), g] = 1.0
    ZZ = Z @ Z.T
    v_known = np.asarray(v_known, float)

    def build_V(theta):
        tp2 = theta[0]
        ts2 = theta[1] if fit_section_var else 0.0
        return np.diag(v_known + ts2) + tp2 * ZZ

    def neg_reml(log_theta):
        theta = np.exp(log_theta)
        V = build_V(theta)
        try:
            L = np.linalg.cholesky(V)
        except np.linalg.LinAlgError:
            return 1e10
        Vi = np.linalg.inv(V)
        XtViX = X.T @ Vi @ X
        b = np.linalg.solve(XtViX, X.T @ Vi @ y)
        r = y - X @ b
        logdetV = 2.0 * np.log(np.diag(L)).sum()
        sign, logdetX = np.linalg.slogdet(XtViX)
        return 0.5 * (logdetV + logdetX + r @ Vi @ r)

    k = 2 if fit_section_var else 1
    scale = max(np.var(y), 1e-6)
    best = None
    # multi-start on a log grid; variance components may sit at the boundary (0)
    for s0 in itertools.product([-12.0, np.log(scale * 0.1), np.log(scale)], repeat=k):
        res = optimize.minimize(neg_reml, np.array(s0), method="L-BFGS-B",
                                bounds=[(-25.0, np.log(scale * 50.0))] * k)
        if best is None or res.fun < best.fun:
            best = res
    theta = np.exp(best.x)
    V = build_V(theta)
    Vi = np.linalg.inv(V)
    XtViX = X.T @ Vi @ X
    cov_b = np.linalg.inv(XtViX)
    b = cov_b @ X.T @ Vi @ y
    out = dict(beta=b.tolist(), cov_beta=cov_b.tolist(), tau2_patient=float(theta[0]),
               tau2_section=float(theta[1]) if fit_section_var else None,
               neg_reml=float(best.fun), n=int(n), J=int(len(labels)))
    # patient BLUPs (conditional modes)
    r = y - X @ b
    out["patient_blup"] = {str(lab): float(theta[0] * (Z[:, i] @ Vi @ r))
                           for i, lab in enumerate(labels)}
    return out


def nested_mean_test(y, groups, v_known, fit_section_var=True, label=""):
    """Pooled-mean inference with t(J-1) reference (patients are the replication unit)."""
    fit = _nested_reml(y, groups, v_known, fit_section_var=fit_section_var)
    b = fit["beta"][0]
    se = float(np.sqrt(fit["cov_beta"][0][0]))
    df = fit["J"] - 1
    t = b / se
    tcrit = stats.t.ppf(0.975, df)
    fit.update(dict(label=label, mean=float(b), se=se, df=int(df), t=float(t),
                    ci95=[float(b - tcrit * se), float(b + tcrit * se)],
                    p_two_sided=float(2 * stats.t.sf(abs(t), df)),
                    p_one_sided_positive=float(stats.t.sf(t, df))))
    # share of total between-unit variance that is between patients (ICC-like)
    tp, ts = fit["tau2_patient"], (fit["tau2_section"] or 0.0)
    fit["between_patient_share"] = float(tp / (tp + ts)) if (tp + ts) > 0 else None
    return fit


def patient_means(values: dict, pmap: dict) -> dict:
    by = {}
    for sid, v in values.items():
        by.setdefault(pmap[sid], []).append(v)
    return {p: float(np.mean(v)) for p, v in sorted(by.items())}


def analyse(values: dict, sds: dict | None, pmap: dict, label: str,
            fit_section_var=True) -> dict:
    sids = list(values)
    y = np.array([values[s] for s in sids], float)
    groups = np.array([pmap[s] for s in sids])
    v = np.array([sds[s] ** 2 for s in sids], float) if sds else np.zeros(len(sids))
    pm = patient_means(values, pmap)
    res = dict(label=label, n_sections=len(sids), n_patients=len(pm),
               sections=sids, patient_means=pm,
               section_n_positive=int((y > 0).sum()),
               sign_flip=sign_flip_exact(np.array(list(pm.values()))))
    res["nested_model"] = nested_mean_test(y, groups, v, fit_section_var=fit_section_var,
                                           label=label)
    # leave-one-patient-out
    lopo = {}
    for p in sorted(set(groups)):
        keep = groups != p
        if len(set(groups[keep])) < 2:
            continue
        f = nested_mean_test(y[keep], groups[keep], v[keep],
                             fit_section_var=fit_section_var, label=f"{label}-minus-{p}")
        lopo[p] = dict(mean=f["mean"], ci95=f["ci95"], p_two_sided=f["p_two_sided"],
                       n_sections=int(keep.sum()), n_patients=int(len(set(groups[keep]))))
    res["leave_one_patient_out"] = lopo
    res["lopo_min_mean"] = float(min(d["mean"] for d in lopo.values())) if lopo else None
    res["lopo_max_p_two_sided"] = float(max(d["p_two_sided"] for d in lopo.values())) if lopo else None
    return res


# --------------------------------------------------------------------------
def main():
    cfg = load_config(None)
    P = Paths(cfg)
    pmap = patient_map(P)
    out = dict(meta=stamp_run(cfg, {"module": "M36-patient-level"}))

    sn = load_json(P.validation("spatial_null_check.json"))["per_slide"]
    primary = {s: r["real_rho_partial"] for s, r in sn.items()}
    sd = {s: r["null_std"] for s, r in sn.items()}
    out["primary_association"] = analyse(primary, sd, pmap,
                                         "partial rho, all 19 primary sections")

    filled = {s: r.get("filled_keys", []) for s, r in sn.items()}
    ag_ok = {s: v for s, v in primary.items() if "ag_target" not in filled[s]}
    out["primary_association_ag_available"] = analyse(
        ag_ok, {s: sd[s] for s in ag_ok}, pmap, "partial rho, Ag_target-available sections")

    cscc = {s: v for s, v in primary.items() if s.startswith("CSCC")}
    out["primary_association_cscc_only"] = analyse(
        cscc, {s: sd[s] for s in cscc}, pmap, "partial rho, cSCC sections only")

    # largest-component re-estimate (same SDs reused: the component restriction
    # changes at most a few per cent of nodes)
    mc = load_json(P.validation("metric_connectivity_sensitivity.json"))
    lcc = {}
    for s, r in mc["per_slide"].items():
        for key in ("rho_lcc_only", "rho_partial_lcc", "lcc_rho_partial"):
            if isinstance(r, dict) and key in r:
                lcc[s] = r[key]
                break
        else:
            sub = r.get("largest_component_association") if isinstance(r, dict) else None
            if isinstance(sub, dict) and "rho_partial" in sub:
                lcc[s] = sub["rho_partial"]
    if len(lcc) == len(primary):
        out["primary_association_lcc"] = analyse(lcc, {s: sd[s] for s in lcc}, pmap,
                                                 "partial rho, largest component only")

    # S2 selection-matched counterfactual: log ratio, residual variance estimated
    s2 = load_json(P.validation("s2_matched_selection.json"))["per_slide"]
    s2log = {s: float(np.log(r["ratio_vs_in_cut_matched"])) for s, r in s2.items()}
    out["s2_matched_log_ratio"] = analyse(s2log, None, pmap,
                                          "log(ratio) selection-matched S2, 19 sections",
                                          fit_section_var=True)

    # optional inputs produced later in the pipeline
    for fname, key, field, sdfield in [
        ("shared_input_spatial_null.json", "ecm_ablated_association", "real_rho_partial", "null_std"),
        ("geometry_null.json", "construction_null_excess", "excess_over_construction_null", "construction_null_sd"),
    ]:
        p = P.validation(fname)
        if p.exists():
            d = load_json(p)["per_slide"]
            d = {s: r for s, r in d.items() if s in primary}
            vals = {s: r[field] for s, r in d.items() if r.get(field) is not None}
            sds = {s: r[sdfield] for s, r in d.items() if r.get(sdfield) is not None}
            if len(vals) >= 3:
                out[key] = analyse(vals, sds if len(sds) == len(vals) else None, pmap,
                                   f"{field} ({fname})")

    save_json(P.validation("patient_level_inference.json"), out)

    def line(r):
        nm = r["nested_model"]
        sf = r["sign_flip"]
        return (f"{r['label']:<52} J={r['n_patients']} n={r['n_sections']:>2}  "
                f"mean={nm['mean']:+.3f} [{nm['ci95'][0]:+.3f},{nm['ci95'][1]:+.3f}] "
                f"p2={nm['p_two_sided']:.4f}  sign-flip p1={sf['p_one_sided']:.4f} "
                f"({sf['n_positive']}/{sf['J']} patients >0)  "
                f"LOPO min={r['lopo_min_mean']:+.3f} max p={r['lopo_max_p_two_sided']:.4f}")
    for k, r in out.items():
        if isinstance(r, dict) and "nested_model" in r:
            print(line(r))
    print("\nWrote", P.validation("patient_level_inference.json"))


if __name__ == "__main__":
    main()
