#!/usr/bin/env python
"""
run_53_intervention_domain_enrichment.py -- where do the high-impact intervention spots sit?
==========================================================================================
Inputs : results/intervention/{sid}.intervention.npz (run_44: per-spot single-ablation effects and
         the exact minimum cut), data/interim/{sid}.nodes.npz, {sid}.graph.npz
Output : results/validation/intervention_domain_enrichment.json

For each operator (B_cell, B_mAb) and budget (top 1 % and top 5 % of spots by single-spot effect),
the share of top spots in the tumour, stromal and immune compartments is compared with the share
of all spots (random placement). Compartments are those of run_52 (tumour = malignant score >= 70th
percentile; other spots stromal or immune by the larger of the mean (ECM, CAF) and mean (T/NK,
myeloid, B/plasma) rank scores).

Inference
  within section   hypergeometric tail probability (descriptive: spots are spatially
                   autocorrelated, so these p-values are optimistic);
  pooled           exact null distribution of the total count over sections (convolution of the
                   per-section hypergeometric distributions; random placement within each section);
  across patients  difference between the top-spot share and the section share, averaged within
                   patient, then an exact sign-flip test across patients -- robust to the spatial
                   autocorrelation of spots and to the small budgets of small sections;
  cut-band check   B_cell effects are non-zero only on the minimum cut (max-flow/min-cut duality), so
                   B_cell top spots are also compared with the composition of the cut band itself.
LN01 (lymph node, no tumour) is reported but excluded from the pooled tests.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts" / "is_figures"))

import numpy as np  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

from sparta.io_ import Paths, load_config, patient_map, save_json, stamp_run  # noqa: E402
from sparta.node_tables import load_nodes  # noqa: E402

DOMAINS = ("tumour", "stroma", "immune")
FRACS = (0.01, 0.05)


def signflip_p(x):
    """Two-sided exact sign-flip p for the mean of x (enumerates all sign patterns up to 16)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    J = len(x)
    if J == 0:
        return float("nan"), float("nan")
    obs = abs(x.mean())
    if J <= 16:
        signs = np.array(list(itertools.product((-1, 1), repeat=J)))
    else:
        rng = np.random.default_rng(0)
        signs = rng.choice((-1, 1), size=(100000, J))
    null = np.abs((signs * x).mean(axis=1))
    p2 = float((np.sum(null >= obs - 1e-12)) / len(null))
    pos = (signs * x).mean(axis=1)
    p_hi = float(np.mean(pos >= x.mean() - 1e-12))
    return p2, p_hi


def pooled_hypergeom_p(rows):
    """One-sided p-values (over, under) for the total count sum_s x_s under independent
    hypergeometric draws per section; rows = [(x, n, K, k), ...]."""
    pmf = np.array([1.0])
    for _, n, K, k in rows:
        lo, hi = max(0, k - (n - K)), min(k, K)
        support = np.arange(0, hi + 1)
        p = np.zeros(hi + 1)
        p[lo:] = hypergeom.pmf(support[lo:], n, K, k)
        pmf = np.convolve(pmf, p)
    X = int(sum(r[0] for r in rows))
    cdf = np.cumsum(pmf)
    p_over = float(1.0 - (cdf[X - 1] if X > 0 else 0.0))
    p_under = float(cdf[min(X, len(cdf) - 1)])
    return min(1.0, max(p_over, 0.0)), min(1.0, p_under)


def main():
    from isdata import EXTERNAL_ORDER, PRIMARY_ORDER, REPLICATION_ORDER
    from run_52_simple_baselines_spearman import compartments
    cfg = load_config(None)
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    pmap = patient_map(P)
    cohort = {**{s: "primary" for s in PRIMARY_ORDER}, **{s: "external" for s in EXTERNAL_ORDER},
              **{s: "replication" for s in REPLICATION_ORDER}}
    t0 = time.time()
    per = {}
    for sid, coh in cohort.items():
        z = np.load(Path(P.results) / "intervention" / f"{sid}.intervention.npz", allow_pickle=True)
        lab, _ = compartments(load_nodes(P.interim / f"{sid}.nodes.npz"))
        n = len(lab)
        if len(z["delta_b_cell"]) != n:
            raise SystemExit(f"{sid}: intervention array has {len(z['delta_b_cell'])} spots, node table {n}")
        cut = np.asarray(z["cut_mask"], bool)
        base = {d: float(np.mean(lab == i)) for i, d in enumerate(DOMAINS)}
        rec = dict(cohort=coh, patient=pmap.get(sid, sid), n=n, base_frac=base,
                   cut_frac={d: float(np.mean(lab[cut] == i)) if cut.any() else float("nan")
                             for i, d in enumerate(DOMAINS)}, n_cut=int(cut.sum()), top={})
        for op in ("b_cell", "b_mab"):
            delta = np.asarray(z[f"delta_{op}"], float)
            order = np.argsort(-delta, kind="stable")
            pos = int(np.sum(delta > 0))
            for f in FRACS:
                k = max(1, int(round(f * n)))
                k_eff = min(k, pos) if pos > 0 else k
                top = order[:k_eff]
                r = dict(k=k, k_used=int(k_eff))
                for i, d in enumerate(DOMAINS):
                    K = int(np.sum(lab == i))
                    x = int(np.sum(lab[top] == i))
                    frac = x / k_eff
                    r[d] = dict(count=x, frac=float(frac),
                                enrichment=float(frac / base[d]) if base[d] > 0 else float("nan"),
                                p_over=float(hypergeom.sf(x - 1, n, K, k_eff)),
                                p_under=float(hypergeom.cdf(x, n, K, k_eff)))
                    if op == "b_cell" and cut.any():
                        Kc = int(np.sum(lab[cut] == i))
                        top_in_cut = top[cut[top]]
                        xc = int(np.sum(lab[top_in_cut] == i))
                        r[d]["enrichment_vs_cut"] = (float((xc / max(len(top_in_cut), 1)) / (Kc / cut.sum()))
                                                     if Kc > 0 else float("nan"))
                r["frac_in_cut"] = float(np.mean(cut[top])) if cut.any() else float("nan")
                rec["top"][f"{op}_top{int(f * 100)}pct"] = r
        per[sid] = rec
        t1 = rec["top"]["b_cell_top1pct"]
        t2 = rec["top"]["b_mab_top5pct"]
        print(f"{sid:<14} base T/S/I = {base['tumour']:.2f}/{base['stroma']:.2f}/{base['immune']:.2f} | "
              f"B_cell top1% {t1['tumour']['frac']:.2f}/{t1['stroma']['frac']:.2f}/{t1['immune']['frac']:.2f} | "
              f"B_mAb top5% {t2['tumour']['frac']:.2f}/{t2['stroma']['frac']:.2f}/{t2['immune']['frac']:.2f}", flush=True)

    # ---- pooled summaries (tumour sections only) ----
    tumour_sids = [s for s in per if s != "LN01"]
    summ = {}
    for key in per[tumour_sids[0]]["top"]:
        for i, d in enumerate(DOMAINS):
            diff = {s: per[s]["top"][key][d]["frac"] - per[s]["base_frac"][d] for s in tumour_sids}
            pats = sorted({per[s]["patient"] for s in tumour_sids})
            pm = np.array([np.mean([diff[s] for s in tumour_sids if per[s]["patient"] == q]) for q in pats])
            p2, p_hi = signflip_p(pm)
            rows = [(per[s]["top"][key][d]["count"], per[s]["n"], int(round(per[s]["base_frac"][d] * per[s]["n"])),
                     per[s]["top"][key]["k_used"]) for s in tumour_sids]
            p_over, p_under = pooled_hypergeom_p(rows)
            e = [per[s]["top"][key][d]["enrichment"] for s in tumour_sids]
            ent = dict(n_sections=len(tumour_sids), n_patients=len(pats),
                       pooled_frac_top=float(sum(r[0] for r in rows) / sum(r[3] for r in rows)),
                       pooled_frac_all=float(sum(r[2] for r in rows) / sum(r[1] for r in rows)),
                       pooled_count=int(sum(r[0] for r in rows)), pooled_k=int(sum(r[3] for r in rows)),
                       pooled_p_over=p_over, pooled_p_under=p_under,
                       median_enrichment=float(np.median(e)),
                       patient_mean_diff=float(pm.mean()),
                       patient_signflip_p_two_sided=p2, patient_signflip_p_enriched=p_hi,
                       n_patients_enriched=int(np.sum(pm > 0)), n_patients_depleted=int(np.sum(pm < 0)),
                       n_sections_p_over_lt_005=int(np.sum([per[s]["top"][key][d]["p_over"] < 0.05
                                                            for s in tumour_sids])),
                       n_sections_p_under_lt_005=int(np.sum([per[s]["top"][key][d]["p_under"] < 0.05
                                                             for s in tumour_sids])))
            if key.startswith("b_cell"):
                cut_rows = []
                dc = {}
                for s in tumour_sids:
                    z = np.load(Path(P.results) / "intervention" / f"{s}.intervention.npz", allow_pickle=True)
                    cut = np.asarray(z["cut_mask"], bool)
                    if not cut.any():
                        continue
                    k_used = per[s]["top"][key]["k_used"]
                    frac_cut = per[s]["cut_frac"][d]
                    dc[s] = per[s]["top"][key][d]["frac"] - frac_cut
                    cut_rows.append((per[s]["top"][key][d]["count"], int(cut.sum()),
                                     int(round(frac_cut * cut.sum())), min(k_used, int(cut.sum()))))
                pmc = np.array([np.mean([dc[s] for s in dc if per[s]["patient"] == q]) for q in pats
                                if any(per[s]["patient"] == q for s in dc)])
                p2c, _ = signflip_p(pmc)
                po_c, pu_c = pooled_hypergeom_p(cut_rows)
                ent.update(pooled_frac_cut=float(sum(r[2] for r in cut_rows) / sum(r[1] for r in cut_rows)),
                           pooled_p_over_vs_cut=po_c, pooled_p_under_vs_cut=pu_c,
                           patient_mean_diff_vs_cut=float(pmc.mean()), patient_signflip_p_vs_cut=p2c)
            summ.setdefault(key, {})[d] = ent
    out = dict(per_section=per, summary=summ,
               settings=dict(domains=list(DOMAINS), fracs=list(FRACS),
                             compartments="run_52.compartments (tumour = malignant >= q0.70; stroma vs immune by larger mean rank score)",
                             pooled_sections="all sections except LN01 (non-tumour)",
                             test=("pooled: exact convolution of per-section hypergeometric nulls; "
                                   "patient level: exact sign-flip on the mean difference top share - section share")),
               meta=stamp_run(cfg, {"module": "M53-intervention-domain-enrichment",
                                    "seconds": round(time.time() - t0, 1)}))
    save_json(P.validation("intervention_domain_enrichment.json"), out)
    for key, v in summ.items():
        print(key)
        for d, x in v.items():
            print(f"   {d:<7} top {x['pooled_frac_top']:.3f} vs all {x['pooled_frac_all']:.3f}"
                  + (f" (cut {x['pooled_frac_cut']:.3f}, p_vs_cut over/under {x['pooled_p_over_vs_cut']:.3g}/{x['pooled_p_under_vs_cut']:.3g})" if 'pooled_frac_cut' in x else "")
                  + f" | pooled p over/under {x['pooled_p_over']:.3g}/{x['pooled_p_under']:.3g}"
                  f" | patients +/-: {x['n_patients_enriched']}/{x['n_patients_depleted']} sign-flip p {x['patient_signflip_p_two_sided']:.3f}")


if __name__ == "__main__":
    main()
