"""Fig. S9 -- technical-summary PCA for extension cohort batch structure."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E409
from sklearn.preprocessing import StandardScaler  # noqa: E402

import figstyle as S  # noqa: E402
from isdata import j  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/figures/is/supplement"


def safe(value, default=np.nan):
    try:
        if value is None:
            return default
        return float(value)
    except Exception:
        return default


def main():
    S.apply()
    cohort = j("extension_cohort.json")["per_slide"]
    qc_path = ROOT / "results/qc/extension_2026_qc_summary.json"
    qc = __import__("json").loads(qc_path.read_text(encoding="utf-8"))
    slides = sorted(set(cohort) & set(qc))

    feature_names = [
        "log raw spots", "raw genes", "log retained spots", "retained genes",
        "log raw median UMI", "log retained median UMI", "retained detected genes",
        "median MT fraction", "log n nodes", "log n source", "log n sink",
        "log n vessel", "rho", "null SD", "geometry null", "construction null",
        "construction share", "ablated rho",
    ]
    X = []
    for sid in slides:
        c, q = cohort[sid], qc[sid]
        n_edges = c.get("n_edges", np.nan)
        if not np.isfinite(safe(n_edges)):
            # Not stored directly; graph artifacts could be read, but degree is not needed.
            n_edges = np.nan
        X.append([
            np.log10(q["raw_spots"]), q["raw_genes"], np.log10(q["retained_spots"]),
            q["retained_genes"], np.log10(max(q["raw_median_umi"], 1)),
            np.log10(max(q["retained_median_umi"], 1)), q["retained_median_detected_genes"],
            q.get("median_mito_fraction", -1) if q.get("mitochondria_assessable") else -1,
            np.log10(c["n_nodes"]), np.log10(c["n_source"]), np.log10(c["n_sink"]),
            np.log10(c["n_vessel"]), c["real_rho_partial"], c["null_std"],
            c["geometry_null_mean"], c["construction_null_mean"],
            c["share_reproduced_by_construction"], c["ablated_rho"],
        ])
    X = np.asarray(X, float)
    # Impute rare missing values using feature medians.
    for j_col in range(X.shape[1]):
        col = X[:, j_col]
        col[~np.isfinite(col)] = np.nanmedian(col[np.isfinite(col)])

    Z = StandardScaler().fit_transform(X)
    P = PCA(n_components=2, random_state=20261005).fit_transform(Z)

    projects = sorted({cohort[s]["project"] for s in slides})
    cmap = S.mpl.colormaps["tab10"].resampled(len(projects))
    colors = {p: cmap(i) for i, p in enumerate(projects)}
    markers = {"Visium": "o", "Slide-seqV2": "s"}

    fig, ax = plt.subplots(figsize=(S.FULL_W, 78 * S.MM))
    for project in projects:
        for platform in ("Visium", "Slide-seqV2"):
            idx = [i for i, s in enumerate(slides)
                   if cohort[s]["project"] == project and cohort[s]["platform"] == platform]
            if not idx:
                continue
            ax.scatter(P[idx, 0], P[idx, 1], s=24, marker=markers[platform],
                       facecolors=colors[project], edgecolors=S.INK, linewidths=0.45,
                       label=project if platform == "Visium" else None, zorder=3)
        idx = [i for i, s in enumerate(slides) if cohort[s]["project"] == project]
        ax.text(P[idx, 0].mean(), P[idx, 1].mean(), project.replace("GSE", ""),
                fontsize=6.5, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7))

    ax.axhline(0, color=S.HAIR, lw=0.6)
    ax.axvline(0, color=S.HAIR, lw=0.6)
    ax.set_xlabel("PC1 of technical and structural summaries")
    ax.set_ylabel("PC2")
    project_handles = [
        __import__("matplotlib.lines", fromlist=["Line2D"]).Line2D(
            [0], [0], marker="o", ls="none", ms=5, mfc=colors[p], mec=S.INK, mew=0.5)
        for p in projects
    ]
    platform_handles = [
        __import__("matplotlib.lines", fromlist=["Line2D"]).Line2D(
            [0], [0], marker=markers[p], ls="none", ms=5, mfc="white", mec=S.INK, mew=0.7)
        for p in ("Visium", "Slide-seqV2")
    ]
    ax.legend(project_handles + platform_handles,
              projects + ["Visium", "Slide-seqV2"], fontsize=6, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)
    fig.tight_layout()
    S.save(fig, OUT, "FigS9_batch_effect", eps=False)


if __name__ == "__main__":
    main()
