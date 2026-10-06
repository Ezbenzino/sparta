#!/usr/bin/env python
"""
run_56_extension_cohort.py
==========================

Run the locked SPARTA spatial pipeline and patient-level inference for the 2026
extension cohort.

Inputs:
    data/interim/extension_2026/{slide}.h5ad  QC-filtered raw counts
    data/ledger.csv                            extension metadata
Outputs:
    data/interim/{slide}.graph.npz
    data/interim/{slide}.barrier.npz
    data/interim/{slide}.nodes.npz
    results/validation/extension_cohort.json

The extension cohort is reported as a separate family. It is not silently pooled
into the primary cohort; a combined primary+extension analysis is explicitly
labelled.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from sparta import lite  # noqa: E402
from sparta.barrier import (  # noqa: E401
    _SCORE_KEYS, compute_b_cell, compute_b_cell_field, compute_b_mab,
    compute_b_meta, save_barriers,
)
from sparta.graph import build_graph_radius, define_source_sink, graph_report  # noqa: E402
from sparta.io_ import (  # noqa: E402
    Paths, load_config, load_graph, load_json, patient_map, save_graph,
    save_json, stamp_run,
)

N_SUR = 500
N_NULL = 200
SEED = 20261005
VISIUM_RADIUS_UM = 150.0
SLIDESEQ_BIN_UM = 50.0
SLIDESEQ_RADIUS_UM = 100.0


@dataclass
class ExtensionSection:
    """Ledger identity plus platform-specific analysis settings."""

    slide_id: str
    patient_id: str
    cancer_type: str
    platform: str
    project_id: str
    site: str
    treatment: str
    replicate: str
    h5ad_path: Path
    graph_radius_um: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats-only", action="store_true",
                    help="skip graph/barrier/node construction and recompute statistics")
    ap.add_argument("--slides", nargs="*", default=None,
                    help="optional subset of extension slide IDs")
    return ap.parse_args()


def load_extension_sections(slide_subset: set[str] | None = None) -> list[ExtensionSection]:
    """Read extension rows from the ledger and pair them with h5ad files."""
    with (ROOT / "data" / "ledger.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    sections: list[ExtensionSection] = []
    for row in rows:
        if row.get("status") != "extension":
            continue
        slide_id = row["slide_id"]
        if slide_subset is not None and slide_id not in slide_subset:
            continue
        platform_label = row["platform"]
        platform = "Slide-seqV2" if platform_label == "slideseqv2" else "Visium"
        radius = SLIDESEQ_RADIUS_UM if platform == "Slide-seqV2" else VISIUM_RADIUS_UM
        h5ad = ROOT / "data" / "interim" / "extension_2026" / f"{slide_id}.h5ad"
        if not h5ad.exists():
            raise FileNotFoundError(f"Missing h5ad for {slide_id}: {h5ad}")
        sections.append(ExtensionSection(
            slide_id=slide_id,
            patient_id=row["patient"],
            cancer_type=row["cancer_type"],
            platform=platform,
            project_id=row["source"].replace("GEO ", ""),
            site=row["site"],
            treatment=row["treatment"],
            replicate=row["replicate"],
            h5ad_path=h5ad,
            graph_radius_um=radius,
        ))
    return sorted(sections, key=lambda item: (item.project_id, item.slide_id))


def hypoxia_genes(cfg) -> tuple[list[str], str]:
    """Load HALLMARK_HYPOXIA genes from the configured GMT."""
    path = ROOT / cfg["signatures"]["hypoxia_gmt"]
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if "HYPOXIA" in fields[0].upper():
                return fields[2:], str(cfg["signatures"]["hypoxia_gmt"])
    raise SystemExit(f"No HYPOXIA row in {path}")


def rank_column(rn: dict, key: str, n: int) -> np.ndarray:
    """Return a rank score or a neutral fill."""
    value = rn.get(key)
    return np.asarray(value, float) if value is not None else np.full(n, 0.5)


def aggregate_slideseq_bins(adata, bin_um: float = SLIDESEQ_BIN_UM):
    """Aggregate QC-passed Slide-seq beads into fixed micrometre grid bins."""
    import scipy.sparse as sp
    from anndata import AnnData

    xy = np.asarray(adata.obsm["spatial_um"], dtype=float)
    origin = np.floor(xy.min(axis=0))
    bin_index = np.floor((xy - origin) / bin_um).astype(int)
    bin_keys = np.char.add(
        bin_index[:, 0].astype(str),
        np.char.add("|", bin_index[:, 1].astype(str)),
    )
    unique_bins, inverse = np.unique(bin_keys, return_inverse=True)
    n_beads = adata.n_obs
    n_bins = len(unique_bins)
    row_sum = sp.csr_matrix(
        (np.ones(n_beads), (inverse, np.arange(n_beads))),
        shape=(n_bins, n_beads),
    )
    counts_binned = (row_sum @ adata.layers["counts"].tocsr()).tocsr()
    counts_per_bin = np.bincount(inverse)
    centers = np.column_stack([
        np.bincount(inverse, weights=xy[:, 0]) / counts_per_bin,
        np.bincount(inverse, weights=xy[:, 1]) / counts_per_bin,
    ])
    representative_positions = [
        int(np.flatnonzero(inverse == j)[0]) for j in range(n_bins)
    ]
    representative_barcodes = adata.obs["barcode"].to_numpy()[representative_positions]
    grouped_barcodes = [
        adata.obs["barcode"].to_numpy()[inverse == j].tolist()
        for j in range(n_bins)
    ]
    obs = adata.obs.iloc[representative_positions].copy()
    obs["barcode"] = representative_barcodes
    obs["n_beads_aggregated"] = counts_per_bin
    obs["group_barcodes_json"] = [
        json.dumps(items, ensure_ascii=False) for items in grouped_barcodes
    ]
    binned = AnnData(
        X=counts_binned,
        obs=obs.set_index(adata.obs.index[representative_positions], drop=True),
        var=adata.var.copy(),
        obsm={"spatial_um": centers},
        layers={"counts": counts_binned},
        uns=dict(adata.uns),
    )
    binned.uns["slideseq_grid_aggregation"] = {
        "bin_size_um": bin_um,
        "n_beads_before": int(n_beads),
        "n_bins_after": int(n_bins),
        "counts_summed": True,
    }
    return binned


def build_section(section: ExtensionSection, cfg, P, hyp: list[str], hyp_src: str) -> None:
    """Build scores, graph, barriers and node table for one extension section."""
    import anndata

    adata = anndata.read_h5ad(section.h5ad_path)
    if section.platform == "Slide-seqV2":
        adata = aggregate_slideseq_bins(adata)
    counts = adata.layers["counts"].tocsr()
    gene_names = adata.var["gene_name"].astype(str).tolist()
    X = lite.normalize_log1p_sparse(counts)
    scores, skipped = lite.score_all_sparse(
        X,
        gene_names,
        section.cancer_type,
        hypoxia_genes=hyp,
        min_genes=cfg["signatures"]["min_genes_matched"],
        ctrl_size=cfg["signatures"]["ctrl_size"],
        seed=cfg["seed"],
    )
    rn = lite.rank_normalized(scores)
    coords_um = np.asarray(adata.obsm["spatial_um"], float)
    n = counts.shape[0]

    A, D = build_graph_radius(coords_um, radius_um=section.graph_radius_um)
    rep = graph_report(A)
    source_sink = define_source_sink(
        A,
        rank_column(rn, "Endothelial_n", n),
        rank_column(rn, "T_NK_n", n),
        rank_column(rn, "Malignant_n", n),
        **{key: cfg["source_sink"][key]
           for key in ("q_vessel", "q_immune_nbr", "q_malig", "q_core")},
        fallback_border_coords=(
            coords_um if cfg["source_sink"]["use_border_fallback"] else None
        ),
    )
    save_graph(
        P.graph(section.slide_id), A, D,
        source_sink["source"], source_sink["sink"], source_sink["vessel"],
        meta={
            **rep,
            "radius_um": section.graph_radius_um,
            "used_fallback": source_sink["used_fallback"],
            "extension_platform": section.platform,
            **stamp_run(cfg, {"module": "M56-extension-graph",
                              "slide": section.slide_id}),
        },
    )

    S = {
        key: rank_column(rn, column, n)
        for key, column in _SCORE_KEYS.items()
    }
    barrier_cfg = cfg["barrier"]
    bc = compute_b_cell(
        A, S["ecm"], S["caf"], source_sink["source"], source_sink["sink"],
        **barrier_cfg["b_cell"],
    )
    bm = compute_b_mab(
        A, S["ecm"], S["crosslink"], S["ag_target"], source_sink["vessel"],
        **barrier_cfg["b_mab"],
    )
    bt = compute_b_meta(
        A, S["hypoxia"], S["proliferation"], S["efflux"],
        source_sink["vessel"], D=D, **barrier_cfg["b_meta"],
    )
    save_barriers(
        P.barrier(section.slide_id), {**bc, **bm, **bt},
        P.mincut(section.slide_id),
    )

    total_counts = np.asarray(counts.sum(axis=1)).ravel()
    node_table: dict[str, Any] = {
        f"obs__{key}": np.asarray(value, float)
        for key, value in scores.items()
    }
    node_table.update({
        f"obs__{key}": np.asarray(value, float)
        for key, value in rn.items()
    })
    node_table.update({
        "obs_names": np.asarray(adata.obs["barcode"].astype(str)),
        "obsm__spatial_um": coords_um,
        "obs__total_counts": total_counts,
        "obs__n_counts": total_counts,
        "obs__n_genes_by_counts": np.asarray((counts > 0).sum(axis=1)).ravel(),
    })
    meta = {
        "slide": section.slide_id,
        "n_obs": int(counts.shape[0]),
        "n_vars": int(counts.shape[1]),
        "source_h5ad": str(section.h5ad_path),
        "platform": section.platform,
        "graph_radius_um": section.graph_radius_um,
        "slideseq_bin_size_um": SLIDESEQ_BIN_UM if section.platform == "Slide-seqV2" else None,
        "skipped_signatures": [item[0] for item in skipped],
        "used_fallback": bool(source_sink["used_fallback"]),
        "m2_run": {
            "module": "M56-extension",
            "tumor_type": section.cancer_type,
            "hypoxia_source": hyp_src,
            "seed": cfg["seed"],
            "scanpy_free": True,
        },
    }
    node_table["meta"] = np.asarray(
        [json.dumps(meta, ensure_ascii=False)], dtype=object
    )
    np.savez_compressed(
        P.interim / f"{section.slide_id}.nodes.npz", **node_table
    )
    print(
        f"{section.slide_id:<22} {section.platform:<10} n={n:>5} "
        f"graph {rep['n_components']} comp LCC {rep['largest_component_frac']:.2f} "
        f"vessel {source_sink['n_vessel']} source {source_sink['n_source']} "
        f"sink {source_sink['n_sink']} Bcell {bc['b_cell']:.3f} "
        f"fallback={source_sink['used_fallback']}",
        flush=True,
    )


def statistics(sections: list[ExtensionSection], cfg, P, pmap) -> dict[str, Any]:
    """Run spatial association, structural nulls, ablation and patient inference."""
    from scipy.stats import spearmanr

    from sparta.node_tables import load_nodes, scores_from_nodes
    from sparta.spatial_stats import (
        control_projector, normal_scores, partial_spearman,
        partial_spearman_many, spectral_basis,
    )
    from run_36_patient_level import analyse
    from run_39_coupling_decomposition import bh, input_surrogate

    barrier_cfg = cfg["barrier"]
    cell_abl = dict(barrier_cfg["b_cell"], b_ecm=0.0, c_caf=12.0)
    mab_abl = dict(barrier_cfg["b_mab"], lam=0.0)
    per: dict[str, dict[str, Any]] = {}
    checkpoint_path = P.validation("extension_cohort_checkpoint.json")
    checkpoint_settings = {
        "n_surrogates": N_SUR,
        "n_structural_nulls": N_NULL,
        "seed": SEED,
        "visium_radius_um": VISIUM_RADIUS_UM,
        "slideseqv2_bin_size_um": SLIDESEQ_BIN_UM,
        "slideseqv2_radius_um": SLIDESEQ_RADIUS_UM,
    }
    if checkpoint_path.exists():
        checkpoint = load_json(checkpoint_path)
        if checkpoint.get("settings") == checkpoint_settings:
            per = checkpoint.get("per_slide", {})
            print(f"resumed {len(per)} completed slides from checkpoint", flush=True)

    for index, section in enumerate(sections):
        sid = section.slide_id
        if sid in per:
            print(f"{sid:<22} checkpoint present, skipped", flush=True)
            continue
        t0 = time.time()
        nodes = load_nodes(P.interim / f"{sid}.nodes.npz")
        A, D, source, sink, vessel, _ = load_graph(P.graph(sid))
        filled: list[str] = []
        S = scores_from_nodes(nodes, missing_out=filled)
        dv = compute_b_meta(
            A, S["hypoxia"], S["proliferation"], S["efflux"], vessel,
            D=D, **barrier_cfg["b_meta"],
        )["d_vessel_um"]
        H = control_projector(dv)
        bc = compute_b_cell_field(
            A, S["ecm"], S["caf"], source, **barrier_cfg["b_cell"],
        )["b_cell_field"]
        bm = compute_b_mab(
            A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
            **barrier_cfg["b_mab"],
        )["b_mab"]
        rho = partial_spearman(bc, bm, dv)
        w, V = spectral_basis(A)
        n = len(bm)
        rng = np.random.default_rng(SEED + index)

        fh = np.abs(V.T @ bm)
        signs = rng.choice([-1.0, 1.0], size=(n, N_SUR))
        null = partial_spearman_many(bc, V @ (signs * fh[:, None]), dv, H)
        p = (np.sum(null >= rho) + 1.0) / (N_SUR + 1.0)

        ns = normal_scores(bm)
        fhn = np.abs(V.T @ ns)
        signs_ns = rng.choice([-1.0, 1.0], size=(n, N_SUR))
        null_ns = partial_spearman_many(
            bc, V @ (signs_ns * fhn[:, None]), dv, H
        )
        p_ns = (np.sum(null_ns >= rho) + 1.0) / (N_SUR + 1.0)

        rho_geo = np.empty(N_NULL)
        rho_con = np.empty(N_NULL)
        for k in range(N_NULL):
            ecm_s = input_surrogate(S["ecm"], V, rng)
            xl_s = input_surrogate(S["crosslink"], V, rng)
            ag_s = input_surrogate(S["ag_target"], V, rng)
            bm_geo = compute_b_mab(
                A, ecm_s, xl_s, ag_s, vessel, **barrier_cfg["b_mab"],
            )["b_mab"]
            bm_con = compute_b_mab(
                A, S["ecm"], xl_s, ag_s, vessel, **barrier_cfg["b_mab"],
            )["b_mab"]
            rho_geo[k] = partial_spearman_many(bc, bm_geo, dv, H)[0]
            rho_con[k] = partial_spearman_many(bc, bm_con, dv, H)[0]

        bc0 = compute_b_cell_field(
            A, S["ecm"], S["caf"], source, **cell_abl,
        )["b_cell_field"]
        bm0 = compute_b_mab(
            A, S["ecm"], S["crosslink"], S["ag_target"], vessel, **mab_abl,
        )["b_mab"]
        rho_abl = partial_spearman(bc0, bm0, dv)
        fh0 = np.abs(V.T @ bm0)
        signs_abl = rng.choice([-1.0, 1.0], size=(n, N_SUR))
        null_abl = partial_spearman_many(
            bc0, V @ (signs_abl * fh0[:, None]), dv, H
        )
        disc = float(np.mean(
            (bc <= np.quantile(bc, 0.25)) & (bm >= np.quantile(bm, 0.75))
        ))

        def input_cor(a: str, b: str):
            if np.ptp(S[a]) == 0 or np.ptp(S[b]) == 0:
                return None
            return float(spearmanr(S[a], S[b])[0])

        per[sid] = {
            "patient": pmap.get(sid, section.patient_id),
            "platform": section.platform,
            "project": section.project_id,
            "cancer_type": section.cancer_type,
            "site": section.site,
            "n_nodes": int(n),
            "filled_keys": filled,
            "n_source": int(len(source)),
            "n_sink": int(len(sink)),
            "n_vessel": int(len(vessel)),
            "real_rho_partial": float(rho),
            "null_mean": float(null.mean()),
            "null_std": float(null.std()),
            "null_ci95": [float(np.percentile(null, 2.5)),
                          float(np.percentile(null, 97.5))],
            "empirical_p_one_sided": float(p),
            "empirical_p_nscore": float(p_ns),
            "geometry_null_mean": float(rho_geo.mean()),
            "geometry_null_sd": float(rho_geo.std()),
            "geometry_null_ci95": [float(np.percentile(rho_geo, 2.5)),
                                   float(np.percentile(rho_geo, 97.5))],
            "p_vs_geometry_null": float((np.sum(rho_geo >= rho) + 1) / (N_NULL + 1)),
            "construction_null_mean": float(rho_con.mean()),
            "construction_null_sd": float(rho_con.std()),
            "construction_null_ci95": [float(np.percentile(rho_con, 2.5)),
                                       float(np.percentile(rho_con, 97.5))],
            "p_vs_construction_null": float((np.sum(rho_con >= rho) + 1) / (N_NULL + 1)),
            "excess_over_construction_null": float(rho - rho_con.mean()),
            "share_reproduced_by_construction": (
                float(rho_con.mean() / rho) if rho > 0 else None
            ),
            "ablated_rho": float(rho_abl),
            "ablated_p": float((np.sum(null_abl >= rho_abl) + 1) / (N_SUR + 1)),
            "ablated_null_sd": float(null_abl.std()),
            "discordant_low_bcell_high_bmab": disc,
            "input_correlations": {
                "ecm_caf": input_cor("ecm", "caf"),
                "ecm_crosslink": input_cor("ecm", "crosslink"),
            },
        }
        save_json(checkpoint_path, {
            "settings": checkpoint_settings,
            "per_slide": per,
        })
        print(
            f"{sid:<22} rho={rho:+.3f} null {null.mean():+.3f}±{null.std():.3f} "
            f"p={p:.3f} ns={p_ns:.3f} constr {rho_con.mean():+.3f} "
            f"abl {rho_abl:+.3f} {time.time() - t0:.0f}s",
            flush=True,
        )
        gc.collect()

    slide_ids = [item.slide_id for item in sections]
    q_bh = bh([per[s]["empirical_p_one_sided"] for s in slide_ids])
    q_ns = bh([per[s]["empirical_p_nscore"] for s in slide_ids])
    q_con = bh([per[s]["p_vs_construction_null"] for s in slide_ids])
    q_abl = bh([per[s]["ablated_p"] for s in slide_ids])
    for sid, a, b, c, d in zip(slide_ids, q_bh, q_ns, q_con, q_abl):
        per[sid].update({
            "bh_padj": float(a),
            "bh_padj_nscore": float(b),
            "bh_padj_construction": float(c),
            "bh_padj_ablated": float(d),
        })

    values = lambda key: {s: per[s][key] for s in slide_ids}
    sds = {s: per[s]["null_std"] for s in slide_ids}
    patient_level = {
        "association": analyse(
            values("real_rho_partial"), sds, pmap,
            "extension: partial rho",
        ),
        "excess": analyse(
            values("excess_over_construction_null"),
            {s: per[s]["construction_null_sd"] for s in slide_ids},
            pmap, "extension: excess over construction null",
        ),
        "ablated": analyse(
            values("ablated_rho"),
            {s: per[s]["ablated_null_sd"] for s in slide_ids},
            pmap, "extension: shared ECM term ablated",
        ),
    }

    grouped = grouped_inference(sections, per, pmap, values, sds)
    primary_sections = [item for item in sections if item.project_id in {"GSE144239", "GSE250636"}]
    primary_ids = {item.slide_id for item in primary_sections}
    pooled_values = {
        s: record["real_rho_partial"]
        for s, record in load_json(
            P.validation("spatial_null_check.json")
        )["per_slide"].items()
    }
    pooled_values.update(values("real_rho_partial"))
    pooled_sds = {
        s: record["null_std"]
        for s, record in load_json(
            P.validation("spatial_null_check.json")
        )["per_slide"].items()
    }
    pooled_sds.update(sds)
    patient_level["pooled_primary_plus_extension"] = analyse(
        pooled_values, pooled_sds, pmap,
        "primary + extension (explicit combined analysis)",
    )

    rho = np.asarray([per[s]["real_rho_partial"] for s in slide_ids])
    summary = {
        "n_sections": len(slide_ids),
        "n_patients": len({per[s]["patient"] for s in slide_ids}),
        "n_positive": int((rho > 0).sum()),
        "median_rho": float(np.median(rho)),
        "min_rho": float(rho.min()),
        "max_rho": float(rho.max()),
        "n_p05": int(sum(per[s]["empirical_p_one_sided"] < 0.05 for s in slide_ids)),
        "n_q05": int(sum(per[s]["bh_padj"] < 0.05 for s in slide_ids)),
        "n_q05_nscore": int(sum(per[s]["bh_padj_nscore"] < 0.05 for s in slide_ids)),
        "median_geometry_null": float(np.median([per[s]["geometry_null_mean"] for s in slide_ids])),
        "median_construction_null": float(np.median([per[s]["construction_null_mean"] for s in slide_ids])),
        "median_share_construction": float(np.median([
            per[s]["share_reproduced_by_construction"]
            for s in slide_ids
            if per[s]["share_reproduced_by_construction"] is not None
        ])),
        "median_excess": float(np.median([per[s]["excess_over_construction_null"] for s in slide_ids])),
        "n_excess_positive": int(sum(per[s]["excess_over_construction_null"] > 0 for s in slide_ids)),
        "n_construction_q05": int(sum(per[s]["bh_padj_construction"] < 0.05 for s in slide_ids)),
        "median_ablated": float(np.median([per[s]["ablated_rho"] for s in slide_ids])),
        "n_ablated_positive": int(sum(per[s]["ablated_rho"] > 0 for s in slide_ids)),
    }
    return {
        "per_slide": per,
        "patient_level": patient_level,
        "grouped": grouped,
        "summary": summary,
    }


def grouped_inference(sections, per, pmap, values, sds) -> dict[str, Any]:
    """Run disease-strand and platform-strand patient inference."""
    from run_36_patient_level import analyse

    groups: dict[str, list[str]] = {
        "cSCC": [],
        "primary_melanoma": [],
        "metastatic_melanoma": [],
        "Visium": [],
        "SlideSeqV2": [],
    }
    for section in sections:
        sid = section.slide_id
        if section.cancer_type == "cscc":
            groups["cSCC"].append(sid)
        if section.project_id in {"GSE300445", "GSE316760"}:
            groups["primary_melanoma"].append(sid)
        if section.project_id in {"GSE200278", "GSE320041"}:
            groups["metastatic_melanoma"].append(sid)
        if section.platform == "Visium":
            groups["Visium"].append(sid)
        else:
            groups["SlideSeqV2"].append(sid)

    out: dict[str, Any] = {}
    for group, slide_ids in groups.items():
        if not slide_ids:
            continue
        group_values = {s: values("real_rho_partial")[s] for s in slide_ids}
        group_sds = {s: sds[s] for s in slide_ids}
        n_patients = len({per[s]["patient"] for s in slide_ids})
        model = (
            analyse(group_values, group_sds, pmap, f"extension group: {group}")
            if len(slide_ids) >= 2 and n_patients >= 2
            else {"skipped": "fewer than two sections or patients"}
        )
        rho = np.asarray([per[s]["real_rho_partial"] for s in slide_ids])
        out[group] = {
            "n_sections": len(slide_ids),
            "n_patients": n_patients,
            "median_rho": float(np.median(rho)),
            "n_positive": int((rho > 0).sum()),
            "patient_model": model,
        }

    # Primary-cohort threshold sensitivity: exclude small Visium deposits.
    large_visium = [
        item.slide_id for item in sections
        if item.platform == "Visium" and per[item.slide_id]["n_nodes"] >= 1000
    ]
    group_values = {s: values("real_rho_partial")[s] for s in large_visium}
    group_sds = {s: sds[s] for s in large_visium}
    n_large_patients = len({per[s]["patient"] for s in large_visium})
    large_model = (
        analyse(
            group_values, group_sds, pmap,
            "extension Visium excluding <1000 spots",
        )
        if len(large_visium) >= 2 and n_large_patients >= 2
        else {"skipped": "fewer than two sections or patients"}
    )
    out["Visium_excluding_sections_below_1000_spots"] = {
        "n_sections": len(large_visium),
        "n_patients": n_large_patients,
        "patient_model": large_model,
        "excluded_slides": [
            item.slide_id for item in sections
            if item.platform == "Visium" and item.slide_id not in large_visium
        ],
    }
    return out


def main() -> None:
    """Build and analyse the extension cohort end to end."""
    args = parse_args()
    t0 = time.time()
    cfg = load_config(ROOT / "configs" / "default.yaml")
    cfg["paths"]["root"] = str(ROOT)
    P = Paths(cfg)
    subset = set(args.slides) if args.slides else None
    sections = load_extension_sections(subset)
    if not sections:
        raise SystemExit("No extension sections selected")
    pmap = patient_map(P)
    hyp, hyp_src = hypoxia_genes(cfg)

    if not args.stats_only:
        for section in sections:
            build_section(section, cfg, P, hyp, hyp_src)

    result = statistics(sections, cfg, P, pmap)
    out = {
        "data": {
            "source": "GEO extension datasets",
            "projects": sorted({item.project_id for item in sections}),
        },
        **result,
        "settings": {
            "n_surrogates": N_SUR,
            "n_structural_nulls": N_NULL,
            "seed": SEED,
            "visium_radius_um": VISIUM_RADIUS_UM,
            "slideseqv2_bin_size_um": SLIDESEQ_BIN_UM,
            "slideseqv2_radius_um": SLIDESEQ_RADIUS_UM,
            "config": "configs/default.yaml",
        },
        "semantics": (
            "Extension family: not pooled into primary-cohort results except in "
            "the explicitly labelled pooled analysis."
        ),
        "meta": stamp_run(cfg, {
            "module": "M56-extension-cohort",
            "seconds": round(time.time() - t0, 1),
        }),
    }
    save_json(P.validation("extension_cohort.json"), out)
    s = result["summary"]
    print(
        f"\nextension: {s['n_positive']}/{s['n_sections']} positive; "
        f"median rho {s['median_rho']:+.3f}; q<0.05 {s['n_q05']}/{s['n_sections']}; "
        f"patients {result['patient_level']['association']['sign_flip']['n_positive']}"
        f"/{s['n_patients']}; construction share {s['median_share_construction']:.2f}"
    )


if __name__ == "__main__":
    main()
