#!/usr/bin/env python
"""
run_51_domain_comparison_exact_cut.py -- BANKSY-style domain boundaries versus the exact minimum cut
===================================================================================================
Inputs : data/interim/{sid}.scored.h5ad (log-normalised expression), {sid}.graph.npz, {sid}.mincut.json
Outputs: results/validation/benchmark_lambda_sensitivity.json   (same schema as run_34)
         results/validation/ndomains_sensitivity.json           (same schema as run_13b)
         results/validation/domain_comparison_exact_cut.json    (reproduction check + old/new cut)

Why this script exists (2026-10-05)
-----------------------------------
run_13 / run_13b / run_34 compare minimum-cut edges with boundaries between BANKSY-style spatial
domains. Their archived outputs were computed with the minimum cuts stored before v2.2, which were
not always exact minimum cuts (see run_50). This script recomputes the same comparison with the
exact cuts. It runs with or without Scanpy:

* with Scanpy, the domains come from ``run_13_benchmark_tools.banksy_style_domains`` exactly as
  before;
* without Scanpy, ``sparta.lite`` re-implements the same recipe (Seurat-flavour HVG, scale with
  max_value 10, ARPACK PCA, k-means on own and neighbour-mean PCs) and ``_h5ad_reader`` reads the
  .h5ad through libhdf5.

Before anything is overwritten, ``--old-cuts DIR`` recomputes the comparison with the previous
cuts and checks it against the archived run_34 numbers; the check is stored in
``domain_comparison_exact_cut.json`` so the equivalence is measured, not assumed.

Usage:
    python scripts/run_51_domain_comparison_exact_cut.py [--h5ad-dir DIR] [--old-cuts DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

from sparta.barrier import edge_pairs  # noqa: E402
from sparta.io_ import Paths, load_config, load_graph, load_json, save_json, stamp_run  # noqa: E402

LAM_GRID = (0.1, 0.3, 0.5, 0.7)        # run_34
N_DOMAINS = 8
N_PCS = 20
K_GRID = [4, 6, 8, 10, 12]             # run_13b
SCALE_DENOM, SCALE_MIN, SCALE_MAX = 200, 4, 16
VISIUM = {"MEL01", "MEL02", "MEL03", "MEL04", "CSCC01", "CSCC02", "CSCC03", "CSCC04"}


def pcs_for(sid, P, h5ad_dir, seed):
    """X_pca as banksy_style_domains computes it (Scanpy if present, else the lite path)."""
    path = (Path(h5ad_dir) / f"{sid}.scored.h5ad") if h5ad_dir else P.scored(sid)
    try:
        import scanpy as sc
        adata = sc.read_h5ad(path)
        if "X_pca" in adata.obsm:
            return np.asarray(adata.obsm["X_pca"]), "obsm", adata.n_obs
        ad = adata.copy()
        sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat")
        ad = ad[:, ad.var.highly_variable].copy()
        sc.pp.scale(ad, max_value=10)
        sc.tl.pca(ad, n_comps=N_PCS, random_state=seed)
        return np.asarray(ad.obsm["X_pca"]), "scanpy", adata.n_obs
    except ImportError:
        from _h5ad_reader import read_h5ad

        from sparta import lite
        d = read_h5ad(path)
        if "X_pca" in d["obsm"]:
            return np.asarray(d["obsm"]["X_pca"]), "obsm", d["X"].shape[0]
        X = d["X"]
        if hasattr(X, "toarray"):
            X = X.toarray()
        return lite.banksy_pcs(np.asarray(X), N_PCS, seed, d["uns_log1p_base"]), "lite", X.shape[0]


def labels_for(X_pca, A, k, lam, seed):
    from sparta import lite
    return lite.banksy_style_domains(X_pca, A, k, lam, N_PCS, seed)


def compare(labels, pairs, cut_idx):
    boundary = labels[pairs[:, 0]] != labels[pairs[:, 1]]
    fb = float(boundary.mean())
    if len(cut_idx):
        fc = float(boundary[cut_idx].mean())
        enr = fc / max(fb, 1e-12)
        prec = float(boundary[cut_idx].sum() / max(boundary.sum(), 1))
    else:
        fc = enr = prec = float("nan")
    return dict(n_domains=int(len(np.unique(labels))), frac_boundary_edges=fb, frac_boundary_within_cut=fc,
                enrichment=float(enr), precision_of_boundary_for_cut=float(prec))


def cut_index(cut, pairs):
    edge_id = {(int(min(u, v)), int(max(u, v))): i for i, (u, v) in enumerate(pairs)}
    return np.array([edge_id[(min(u, v), max(u, v))] for u, v in cut if (min(u, v), max(u, v)) in edge_id], int)


def lambda_summary(rows):
    summary = {}
    for lam in LAM_GRID:
        key = str(lam)
        for cohort, ids in (("all_primary", list(rows)), ("visium", [s for s in rows if s in VISIUM]),
                            ("legacy_st", [s for s in rows if s.startswith("CSCC") and s not in VISIUM])):
            enr = [rows[s]["runs"][key]["enrichment"] for s in ids if np.isfinite(rows[s]["runs"][key]["enrichment"])]
            prc = [rows[s]["runs"][key]["precision_of_boundary_for_cut"] for s in ids
                   if np.isfinite(rows[s]["runs"][key]["precision_of_boundary_for_cut"])]
            summary[f"lambda_{key}_{cohort}"] = dict(n=len(enr), median_enrichment=float(np.median(enr)),
                                                     enrichment_range=[float(np.min(enr)), float(np.max(enr))],
                                                     median_precision=float(np.median(prc)))
    return summary


def ndomain_summary(per):
    legacy = [s for s in per if s not in VISIUM]
    vis = [s for s in per if s in VISIUM]
    summ = {}
    for k in K_GRID:
        for name, grp in (("visium", vis), ("legacy_st", legacy)):
            es = [per[s]["runs"][str(k)]["enrichment"] for s in grp if np.isfinite(per[s]["runs"][str(k)]["enrichment"])]
            ps = [per[s]["runs"][str(k)]["precision_of_boundary_for_cut"] for s in grp
                  if np.isfinite(per[s]["runs"][str(k)]["precision_of_boundary_for_cut"])]
            if es:
                summ[f"fixed_k{k}_{name}"] = dict(median_enrichment=float(np.median(es)),
                                                  median_precision=float(np.median(ps)), n=len(es))
    for name, grp in (("visium", vis), ("legacy_st", legacy)):
        es, ps = [], []
        for s in grp:
            r = per[s]["runs"][str(per[s]["scaled_n_domains"])]
            if np.isfinite(r["enrichment"]):
                es.append(r["enrichment"])
                ps.append(r["precision_of_boundary_for_cut"])
        if es:
            summ[f"scaled_{name}"] = dict(median_enrichment=float(np.median(es)),
                                          median_precision=float(np.median(ps)), n=len(es))
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad-dir", default=None, help="directory holding {sid}.scored.h5ad (default data/interim)")
    ap.add_argument("--old-cuts", default=None, help="directory with the pre-v2.2 {sid}.mincut.json files")
    ap.add_argument("--slides", nargs="*", default=None)
    args = ap.parse_args()

    cfg = load_config(None)
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    seed = int(cfg["seed"])
    s2 = load_json(P.validation("s2_matched_selection.json"))["per_slide"]
    sids = args.slides or sorted(s2)
    old_lam = load_json(P.validation("benchmark_lambda_sensitivity.json"))
    old_nd = load_json(P.validation("ndomains_sensitivity.json"))

    lam_rows, nd_rows, check = {}, {}, {}
    t0 = time.time()
    for sid in sids:
        A, _, _, _, _, _ = load_graph(P.graph(sid))
        X_pca, how, n_obs = pcs_for(sid, P, args.h5ad_dir, seed)
        if n_obs != A.shape[0]:
            raise SystemExit(f"{sid}: {n_obs} spots in the .h5ad but {A.shape[0]} graph nodes")
        pairs = edge_pairs(A, upper_only=True)
        cut_new = cut_index(load_json(P.mincut(sid))["cut_edges"], pairs)
        cut_old = None
        if args.old_cuts:
            cut_old = cut_index(json.load(open(Path(args.old_cuts) / f"{sid}.mincut.json"))["cut_edges"], pairs)
        n_nodes = int(A.shape[0])
        k_scaled = int(np.clip(round(n_nodes / SCALE_DENOM), SCALE_MIN, SCALE_MAX))
        lam_rec = {"n_nodes": n_nodes, "runs": {}}
        nd_rec = {"n_nodes": n_nodes, "scaled_n_domains": k_scaled, "runs": {}}
        chk = {"pcs_from": how, "lambda": {}, "n_domains": {}}
        for lam in LAM_GRID:
            lab = labels_for(X_pca, A, N_DOMAINS, lam, seed)
            lam_rec["runs"][str(lam)] = compare(lab, pairs, cut_new)
            if cut_old is not None:
                o = compare(lab, pairs, cut_old)
                ref = old_lam["per_slide"].get(sid, {}).get("runs", {}).get(str(lam))
                chk["lambda"][str(lam)] = dict(
                    recomputed_old_cut=o, archived=ref,
                    abs_diff_frac_boundary=None if ref is None else abs(o["frac_boundary_edges"] - ref["frac_boundary_edges"]),
                    abs_diff_enrichment=None if ref is None else abs(o["enrichment"] - ref["enrichment"]),
                    enrichment_new_cut=lam_rec["runs"][str(lam)]["enrichment"])
        for k in sorted(set(K_GRID) | {k_scaled}):
            lab = labels_for(X_pca, A, k, 0.3, seed)
            nd_rec["runs"][str(k)] = compare(lab, pairs, cut_new)
            ref = old_nd["per_slide"].get(sid, {}).get("runs", {}).get(str(k))
            if ref is not None:
                chk["n_domains"][str(k)] = dict(
                    abs_diff_frac_boundary=abs(nd_rec["runs"][str(k)]["frac_boundary_edges"] - ref["frac_boundary_edges"]),
                    enrichment_archived=ref["enrichment"], enrichment_new_cut=nd_rec["runs"][str(k)]["enrichment"])
        lam_rows[sid], nd_rows[sid], check[sid] = lam_rec, nd_rec, chk
        e = [lam_rec["runs"][str(l)]["enrichment"] for l in LAM_GRID]
        fbd = [c["abs_diff_frac_boundary"] for c in chk["lambda"].values() if c["abs_diff_frac_boundary"] is not None]
        print(f"{sid:<8} pcs={how:<6} new-cut enrichment {min(e):.3f}-{max(e):.3f}"
              + (f" | reproduction max |d frac_boundary| = {max(fbd):.2e}" if fbd else ""), flush=True)

    meta = stamp_run(cfg, {"module": "M51-domain-comparison-exact-cut", "seconds": round(time.time() - t0, 1),
                           "pcs_from": sorted({c["pcs_from"] for c in check.values()})})
    lam_out = dict(params=dict(implementation="BANKSY-style augmented PCA/neighborhood-mean k-means; not official BANKSY",
                               n_domains=N_DOMAINS, lambda_grid=list(LAM_GRID), n_pcs=N_PCS, seed=seed,
                               fixed_graph=True, comparison="minimum-cut edge overlap with domain-boundary edges",
                               cut="exact minimum cut (v2.2, barrier.exact_min_cut)"),
                   per_slide=lam_rows, summary=lambda_summary(lam_rows), meta=meta)
    nd_out = dict(params=dict(grid=K_GRID, lam=0.3, n_pcs=N_PCS,
                              scale_rule=f"clip(round(n_nodes/{SCALE_DENOM}), {SCALE_MIN}, {SCALE_MAX})", seed=seed,
                              cut="exact minimum cut (v2.2, barrier.exact_min_cut)"),
                  per_slide=nd_rows, summary=ndomain_summary(nd_rows), meta=meta)

    # reproduction summary (only meaningful with --old-cuts)
    rep = {}
    if args.old_cuts:
        fb = [c["abs_diff_frac_boundary"] for s in check.values() for c in s["lambda"].values()
              if c["abs_diff_frac_boundary"] is not None]
        en = [c["abs_diff_enrichment"] for s in check.values() for c in s["lambda"].values()
              if c["abs_diff_enrichment"] is not None]
        exact = [s for s, c in check.items()
                 if all(v["abs_diff_frac_boundary"] is not None and v["abs_diff_frac_boundary"] < 1e-12
                        and v["abs_diff_enrichment"] < 1e-9 for v in c["lambda"].values())]
        old_summary = lambda_summary({s: {"runs": {k: v["recomputed_old_cut"] for k, v in c["lambda"].items()}}
                                      for s, c in check.items()})
        rep = dict(n_sections=len(check), n_sections_exactly_reproduced=len(exact),
                   sections_exactly_reproduced=exact,
                   max_abs_diff_frac_boundary=float(max(fb)) if fb else None,
                   max_abs_diff_enrichment=float(max(en)) if en else None,
                   archived_summary_lambda_0_3=old_lam["summary"]["lambda_0.3_all_primary"],
                   recomputed_old_cut_summary_lambda_0_3=old_summary["lambda_0.3_all_primary"],
                   new_cut_summary_lambda_0_3=lam_out["summary"]["lambda_0.3_all_primary"])
    nd_fb = [c["abs_diff_frac_boundary"] for s in check.values() for c in s["n_domains"].values()]
    rep["ndomains_max_abs_diff_frac_boundary_vs_archived"] = float(max(nd_fb)) if nd_fb else None
    save_json(P.validation("domain_comparison_exact_cut.json"),
              dict(per_section=check, summary=rep, meta=meta,
                   semantics=("Reproduction check: domains recomputed with the path in pcs_from and compared "
                              "with the cuts stored before v2.2 (--old-cuts) must give the archived run_34 "
                              "numbers; frac_boundary_edges depends on the domains only, enrichment also on "
                              "the cut. New-cut numbers replace benchmark_lambda_sensitivity.json and "
                              "ndomains_sensitivity.json.")))
    save_json(P.validation("benchmark_lambda_sensitivity.json"), lam_out)
    save_json(P.validation("ndomains_sensitivity.json"), nd_out)
    print(json.dumps(rep, indent=1)[:3000])


if __name__ == "__main__":
    main()
