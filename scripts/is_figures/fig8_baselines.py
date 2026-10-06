"""Fig. 8 -- cross-cohort robustness, subgroups, baselines and perturbations."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import j  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/figures/is"


def jitter_points(ax, values, x, rng, color=None, marker="o", width=0.22):
    dy = rng.permutation(np.linspace(-width, width, len(values)))
    ax.scatter(np.full(len(values), x) + dy, values, s=7, facecolors="white",
               edgecolors=S.MUTED, linewidths=0.5, marker=marker, zorder=2)


def main():
    S.apply()
    rng = np.random.default_rng(20261005)
    primary_null = j("spatial_null_check.json")["per_slide"]
    geometry = j("geometry_null.json")["per_slide"]
    shared_summary = j("shared_input_spatial_null.json")["summary"]
    replication = j("replication_melanoma.json")["per_slide"]
    extension = j("extension_cohort.json")["per_slide"]
    primary_baseline = j("simple_baselines_spearman.json")["spot_level"]
    extension_baseline = j("extension_simple_baselines_spearman.json")["spot_level"]
    radius = j("radius_sensitivity.json")["per_slide"]
    primary_cal = j("null_calibration_all.json")["summary"]["gaussian"]
    extension_cal = j("extension_null_calibration.json")["summary"]["gaussian"]

    fig = plt.figure(figsize=(S.FULL_W, 132 * S.MM))
    gs = fig.add_gridspec(2, 3, wspace=0.58, hspace=0.62, left=0.075,
                          right=0.985, top=0.95, bottom=0.10)

    # ---------------- a: cross-cohort observed / construction / ablated ----------------
    ax = fig.add_subplot(gs[0, 0])
    cohorts = ["Primary", "Replication", "Extension"]
    observed = [
        np.median([v["real_rho_partial"] for v in primary_null.values()]),
        np.median([v["real_rho_partial"] for v in replication.values()]),
        np.median([v["real_rho_partial"] for v in extension.values()]),
    ]
    construction = [
        np.median([v["construction_null_mean"] for v in geometry.values()]),
        np.median([v["construction_null_mean"] for v in replication.values()]),
        np.median([v["construction_null_mean"] for v in extension.values()]),
    ]
    ablated = [
        shared_summary["median_ablated"],
        np.median([v["ablated_rho"] for v in replication.values()]),
        np.median([v["ablated_rho"] for v in extension.values()]),
    ]
    x = np.arange(3)
    for vals, marker, face, label in [
        (observed, "o", S.INK, "Observed"),
        (construction, "s", "white", "Construction null"),
        (ablated, "D", "white", "Shared ECM ablated"),
    ]:
        ax.plot(x, vals, marker=marker, mfc=face, mec=S.INK, mew=0.8,
                ms=4.6, ls="-", lw=0.8, color=S.MUTED, label=label)
    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts, fontsize=6.8)
    ax.set_ylabel("Median partial Spearman ρ")
    ax.set_ylim(0.0, 0.36)
    ax.legend(fontsize=5.8, handletextpad=0.3, loc="upper right")
    S.panel(ax, "a", x=-0.25, y=1.03)

    # ---------------- b: extension disease subgroups ----------------
    ax = fig.add_subplot(gs[0, 1])
    subgroup_defs = [
        ("cSCC", [v for v in extension.values() if v["cancer_type"] == "cscc"]),
        ("Primary\nmelanoma", [
            v for v in extension.values()
            if v["cancer_type"] == "melanoma" and v["site"] in ("skin", "primary skin")
        ]),
        ("Metastatic\nmelanoma", [
            v for v in extension.values()
            if v["cancer_type"] == "melanoma" and v["site"] not in ("skin", "primary skin")
        ]),
    ]
    for i, (label, records) in enumerate(subgroup_defs):
        vals = [v["real_rho_partial"] for v in records]
        jitter_points(ax, vals, i, rng)
        ax.plot(i, np.median(vals), marker="o", mfc=S.INK, mec=S.INK, ms=4.8, zorder=4)
    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.set_xticks(range(3))
    ax.set_xticklabels([item[0] for item in subgroup_defs], fontsize=6.5)
    ax.set_ylabel("Section partial Spearman ρ")
    ax.set_ylim(-0.12, 0.68)
    S.panel(ax, "b", x=-0.28, y=1.03)

    # ---------------- c: platform subgroups ----------------
    ax = fig.add_subplot(gs[0, 2])
    platform_defs = [
        ("Visium", [v for v in extension.values() if v["platform"] == "Visium"]),
        ("Slide-seqV2", [v for v in extension.values() if v["platform"] == "Slide-seqV2"]),
    ]
    for i, (label, records) in enumerate(platform_defs):
        vals = [v["real_rho_partial"] for v in records]
        jitter_points(ax, vals, i, rng, width=0.28)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.plot([i, i], [q1, q3], color=S.LIGHT, lw=5, zorder=3)
        ax.plot(i, med, "o", mfc=S.INK, mec=S.INK, ms=4.6, zorder=4)
    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.set_xticks(range(2))
    ax.set_xticklabels([item[0] for item in platform_defs], fontsize=6.8)
    ax.set_ylabel("Section partial Spearman ρ")
    ax.set_ylim(-0.12, 0.68)
    S.panel(ax, "c", x=-0.25, y=1.03)

    # ---------------- d: simple baseline comparison ----------------
    ax = fig.add_subplot(gs[1, 0])
    labels = ["cell\ndens", "cell\ndist", "cell\nniche", "mAb\ndens", "mAb\ndist", "mAb\nniche"]
    keys = [
        ("b_cell_field", "stromal_density"),
        ("b_cell_field", "dist_tumour_boundary"),
        ("b_cell_field", "niche_z"),
        ("b_mab_field", "stromal_density"),
        ("b_mab_field", "dist_tumour_boundary"),
        ("b_mab_field", "niche_z"),
    ]
    primary_vals = [primary_baseline[f"{a}__{b}"]["all"]["median"] for a, b in keys]
    extension_vals = [extension_baseline[f"{a}__{b}"]["all"]["median"] for a, b in keys]
    x = np.arange(len(keys))
    ax.plot(x, primary_vals, "o", mfc="white", mec=S.INK, mew=0.8, ms=4.2,
            color=S.MUTED, lw=0.8, label="Primary/replication")
    ax.plot(x, extension_vals, "D", mfc=S.INK, mec=S.INK, ms=4.0,
            color=S.MUTED, lw=0.8, label="Extension")
    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6.2)
    ax.set_ylabel("Median within-section ρ")
    ax.set_ylim(-0.28, 0.62)
    ax.legend(fontsize=5.8, loc="upper right")
    S.panel(ax, "d", x=-0.24, y=1.03)

    # ---------------- e: graph radius perturbation ----------------
    ax = fig.add_subplot(gs[1, 1])
    ratio_bins = np.array([2 / 3, 1.0, 4 / 3, 5 / 3, 2.0])
    binned = {r: [] for r in ratio_bins}
    for sid, record in radius.items():
        default_radius = 150.0 if sid.startswith(("MEL", "CSCC01", "CSCC02", "CSCC03", "CSCC04")) else 200.0
        for radius_text, vals in record.items():
            ratio = float(radius_text) / default_radius
            nearest = ratio_bins[np.argmin(np.abs(ratio_bins - ratio))]
            binned[nearest].append(vals["rho_partial"])
    medians = [np.median(binned[r]) for r in ratio_bins]
    ax.plot(ratio_bins, medians, "-o", color=S.INK, ms=4.0, lw=1.0)
    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.set_xlabel("Graph radius / default radius")
    ax.set_ylabel("Median primary-cohort ρ")
    ax.set_ylim(0.0, 0.34)
    S.panel(ax, "e", x=-0.22, y=1.03)

    # ---------------- f: calibration across cohorts ----------------
    ax = fig.add_subplot(gs[1, 2])
    methods = ["Spectral", "Normal score", "Point level"]
    primary_vals = [
        primary_cal["pooled_fpr_spectral_05"],
        primary_cal["pooled_fpr_nscore_05"],
        primary_cal["pooled_fpr_naive_05"],
    ]
    extension_vals = [
        extension_cal["pooled_fpr_spectral_05"],
        extension_cal["pooled_fpr_nscore_05"],
        extension_cal["pooled_fpr_naive_05"],
    ]
    x = np.arange(3)
    ax.plot(x, primary_vals, "o", mfc="white", mec=S.INK, mew=0.8, ms=4.2,
            color=S.MUTED, lw=0.8, label="Primary/replication")
    ax.plot(x, extension_vals, "D", mfc=S.INK, mec=S.INK, ms=4.0,
            color=S.MUTED, lw=0.8, label="Extension")
    ax.axhline(0.05, color=S.MUTED, lw=0.7, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=6.4)
    ax.set_ylabel("Pooled false-positive rate")
    ax.set_ylim(0, 0.72)
    ax.legend(fontsize=5.8, loc="upper left")
    S.panel(ax, "f", x=-0.24, y=1.03)

    S.save(fig, OUT, "Fig8_baselines", tiff=True, eps=True)


if __name__ == "__main__":
    main()
