#!/usr/bin/env python
"""
run_37_null_calibration.py -- type-I error of the graph-spectral surrogate test on every real graph
===================================================================================================
Inputs : data/interim/{sid}.graph.npz, {sid}.barrier.npz  (22 sections: 19 primary + 3 external)
Output : results/validation/null_calibration_all.json

Question (CMPB readiness review, item 5)
----------------------------------------
The per-section p-values of the field association come from 500 graph-spectral
sign-randomised surrogates of B_mAb.  The earlier calibration (null_calibration.json)
used one graph only.  Here we estimate the false-positive rate on all 22 real
graphs, under fields whose spatial autocorrelation matches the real B_mAb field.

Design (H0 true by construction)
--------------------------------
For each section:
  x  = the real B_mAb field of that section (a smooth, vessel-anchored field
       standing in for "a fixed spatially structured field");
  c  = the real weighted graph distance to the nearest vessel (the covariate);
  y* = an independent random field with the SAME expected graph power spectrum
       as the real B_mAb field:  y* = V (|V^T b_mAb| * e),  e ~ iid N(0, 1)
       ("gaussian" scenario), or a skewed monotone transform exp(2 z(y*))
       ("skewed" scenario), or heat-kernel fields over a smoothness grid
       ("smoothness" scenario).
Then the identical test used in the paper is applied to (x, y*, c): partial
Spearman after quadratic rank-space adjustment, one-sided p from n_sur graph-
spectral sign surrogates of y*, p = (k+1)/(n_sur+1).  The naive point-level
p-value (spots treated as independent) is computed alongside.
Because y* is independent of x, the fraction of p < 0.05 estimates the
type-I error of each test on that graph.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from sparta.io_ import Paths, load_config, load_graph, patient_map, save_json, stamp_run  # noqa: E402
from sparta.spatial_stats import (control_projector, morans_i, naive_spearman_p,  # noqa: E402
                                  normal_scores, partial_spearman_many, spectral_basis)

SEED = 20261004


def one_sided_p(stat, null):
    return (np.sum(null >= stat - 1e-15) + 1.0) / (len(null) + 1.0)


def run_section(sid, P, rng, n_sim, n_sur, scenarios, taus, n_sim_secondary=None):
    A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
    z = np.load(P.barrier(sid), allow_pickle=True)
    x = np.asarray(z["b_mab"], float)
    c = np.asarray(z["d_vessel_um"], float)
    n = len(x)
    t0 = time.time()
    w, V = spectral_basis(A)
    H = control_projector(c)
    amp = np.abs(V.T @ x)
    out = dict(n_nodes=int(n), moran_I_bmab=morans_i(x, A), t_eig_s=round(time.time() - t0, 2))

    def evaluate(make_field, n_sim=n_sim):
        p_spec, p_ns, p_naive, rho = [], [], [], []
        for _ in range(n_sim):
            y = make_field()
            r = float(partial_spearman_many(x, y, c, H)[0])
            # (i) surrogates from the raw field's spectrum (the test used for the locked results)
            fh = np.abs(V.T @ y)
            S = rng.choice([-1.0, 1.0], size=(n, n_sur)) * fh[:, None]
            null = partial_spearman_many(x, V @ S, c, H)
            p_spec.append(one_sided_p(r, null))
            # (ii) surrogates from the spectrum of the normal scores of the field (rank-based variant)
            fz = np.abs(V.T @ normal_scores(y))
            Sz = rng.choice([-1.0, 1.0], size=(n, n_sur)) * fz[:, None]
            null_z = partial_spearman_many(x, V @ Sz, c, H)
            p_ns.append(one_sided_p(r, null_z))
            p_naive.append(float(naive_spearman_p(r, n)) / 2.0 if r > 0 else 1 - float(naive_spearman_p(r, n)) / 2.0)
            rho.append(r)
        p_spec, p_ns, p_naive = np.array(p_spec), np.array(p_ns), np.array(p_naive)
        k_s, k_z, k_n = int((p_spec < 0.05).sum()), int((p_ns < 0.05).sum()), int((p_naive < 0.05).sum())
        return dict(n_sim=n_sim, n_sur=n_sur,
                    fpr_spectral_05=k_s / n_sim, k_spectral_05=k_s,
                    fpr_nscore_05=k_z / n_sim, k_nscore_05=k_z,
                    fpr_spectral_10=float((p_spec < 0.10).mean()),
                    fpr_naive_05=k_n / n_sim, k_naive_05=k_n,
                    mean_p_spectral=float(p_spec.mean()),
                    rho_sd=float(np.std(rho)), rho_mean=float(np.mean(rho)),
                    p_spectral_deciles=np.histogram(p_spec, bins=np.linspace(0, 1, 11))[0].tolist(),
                    p_nscore_deciles=np.histogram(p_ns, bins=np.linspace(0, 1, 11))[0].tolist())

    if "gaussian" in scenarios:
        out["gaussian"] = evaluate(lambda: V @ (amp * rng.standard_normal(n)))
    if "marginal" in scenarios:
        srt = np.sort(x)
        def marg():
            g = V @ (amp * rng.standard_normal(n))
            return srt[np.argsort(np.argsort(g))]      # the real B_mAb marginal imposed on a spectral field
        out["marginal"] = evaluate(marg)
    if "skewed" in scenarios:
        def skew():
            g = V @ (amp * rng.standard_normal(n))
            g = (g - g.mean()) / (g.std() + 1e-12)
            return np.exp(2.0 * g)
        out["skewed"] = evaluate(skew, n_sim_secondary or n_sim)
    if "smoothness" in scenarios:
        sm = {}
        for tau in taus:
            h = np.exp(-tau * np.clip(w, 0, None))
            f = lambda h=h: V @ (h * rng.standard_normal(n))  # noqa: E731
            y0 = f()
            res = evaluate(f, n_sim_secondary or n_sim)
            res["moran_I_example"] = morans_i(y0, A)
            sm[f"tau_{tau:g}"] = res
        out["smoothness"] = sm
    out["t_total_s"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slides", nargs="*", default=None)
    ap.add_argument("--n-sim", type=int, default=200)
    ap.add_argument("--n-sur", type=int, default=100)
    ap.add_argument("--n-sim-secondary", type=int, default=60,
                    help="simulations for the skewed and smoothness scenarios")
    ap.add_argument("--scenarios", nargs="*", default=["gaussian", "marginal", "skewed", "smoothness"])
    ap.add_argument("--taus", nargs="*", type=float, default=[2.0, 10.0, 50.0])
    ap.add_argument("--out", default="null_calibration_all.json")
    args = ap.parse_args()

    cfg = load_config(None)
    P = Paths(cfg)
    pmap = patient_map(P)
    from sparta.io_ import admitted_slides
    slides = args.slides or admitted_slides(P)
    rng = np.random.default_rng(SEED)
    per = {}
    for sid in slides:
        per[sid] = run_section(sid, P, rng, args.n_sim, args.n_sur, args.scenarios, args.taus,
                               args.n_sim_secondary)
        per[sid]["patient"] = pmap.get(sid, sid)
        g = per[sid].get("gaussian", {})
        print(f"{sid:<7} n={per[sid]['n_nodes']:>5}  I={per[sid]['moran_I_bmab']:.3f}  "
              f"FPR spectral={g.get('fpr_spectral_05', float('nan')):.3f}  "
              f"naive={g.get('fpr_naive_05', float('nan')):.3f}  ({per[sid]['t_total_s']} s)",
              flush=True)
        save_json(P.validation(args.out), dict(per_slide=per, partial=True))

    summ = {}
    for scen in args.scenarios:
        if scen == "smoothness":
            for tau in args.taus:
                key = f"tau_{tau:g}"
                ks = sum(per[s]["smoothness"][key]["k_spectral_05"] for s in per)
                kz = sum(per[s]["smoothness"][key]["k_nscore_05"] for s in per)
                kn = sum(per[s]["smoothness"][key]["k_naive_05"] for s in per)
                N = sum(per[s]["smoothness"][key]["n_sim"] for s in per)
                summ[f"smoothness_{key}"] = dict(pooled_fpr_spectral_05=ks / N, pooled_fpr_naive_05=kn / N,
                                                 pooled_fpr_nscore_05=kz / N,
                                                 n_tests=N,
                                                 median_moran_I=float(np.median([per[s]["smoothness"][key]["moran_I_example"] for s in per])))
        else:
            ks = sum(per[s][scen]["k_spectral_05"] for s in per)
            kz = sum(per[s][scen]["k_nscore_05"] for s in per)
            kn = sum(per[s][scen]["k_naive_05"] for s in per)
            N = sum(per[s][scen]["n_sim"] for s in per)
            fpr = [per[s][scen]["fpr_spectral_05"] for s in per]
            summ[scen] = dict(pooled_fpr_spectral_05=ks / N, pooled_fpr_naive_05=kn / N, n_tests=N,
                              pooled_fpr_nscore_05=kz / N,
                              section_fpr_nscore_range=[float(min(per[s][scen]["fpr_nscore_05"] for s in per)),
                                                        float(max(per[s][scen]["fpr_nscore_05"] for s in per))],
                              section_fpr_spectral_range=[float(min(fpr)), float(max(fpr))],
                              section_fpr_naive_range=[float(min(per[s][scen]["fpr_naive_05"] for s in per)),
                                                       float(max(per[s][scen]["fpr_naive_05"] for s in per))],
                              pooled_p_deciles=np.sum([per[s][scen]["p_spectral_deciles"] for s in per], axis=0).tolist(),
                              pooled_p_nscore_deciles=np.sum([per[s][scen]["p_nscore_deciles"] for s in per], axis=0).tolist())
    out = dict(per_slide=per, summary=summ, seed=SEED, n_sim=args.n_sim, n_sur=args.n_sur,
               n_sim_secondary=args.n_sim_secondary,
               design=__doc__.split("Design (H0 true by construction)")[1].strip()[:1200],
               meta=stamp_run(cfg, {"module": "M37-null-calibration"}))
    save_json(P.validation(args.out), out)
    for k, v in summ.items():
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "pooled_p_deciles"})


if __name__ == "__main__":
    main()
