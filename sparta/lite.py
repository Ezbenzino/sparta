"""
Scanpy-free M1/M2 path for table-format spatial data (first-generation ST)
=========================================================================
Input : a spots x genes count table (``loaders.read_table_matrix`` output)
Output: the per-spot arrays that ``run_35_export_node_tables`` writes to ``{sid}.nodes.npz``
        (raw and rank-normalised signature scores, µm coordinates, marker panel)
Upstream  : ``sparta.loaders``
Downstream: ``scripts/run_48_replication_cohort.py`` -> graph / barrier / node tables

Why this exists
---------------
The section pipeline uses Scanpy for quality control and signature scoring (M1, M2). The
melanoma replication cohort had to be processed on a machine without Scanpy, so this module
re-implements exactly the steps the pipeline uses, with NumPy / pandas only:

* ``sc.pp.filter_cells(min_counts)``, ``sc.pp.filter_genes(min_cells)``,
  ``sc.pp.normalize_total(target_sum=1e4)``, ``sc.pp.log1p``;
* ``sc.tl.score_genes`` (Seurat-style: genes binned into 25 expression bins, 50 control genes
  drawn per bin of the gene list, score = mean(list) - mean(controls)), with the legacy
  global-seed sampling of Scanpy <= 1.11;
* within-section rank normalisation ``(rank - 1) / (n - 1)``;
* (added 2026-10-05) the BANKSY-style domain recipe of ``run_13_benchmark_tools``:
  ``sc.pp.highly_variable_genes(flavor="seurat", n_top_genes=2000)``, ``sc.pp.scale(max_value=10)``,
  ``sc.tl.pca(n_comps=20, random_state=0)`` (sklearn PCA, ARPACK solver) and k-means on own and
  neighbour-averaged PCs. ``scripts/run_51_domain_comparison_exact_cut.py`` checks that this path
  reproduces the archived Scanpy results before using it.

The control-gene draw is the only step whose random stream may differ between Scanpy
versions. ``scripts/run_48_replication_cohort.py --validate`` re-scores primary first-generation
sections from their raw tables and compares the result with the archived node tables, so
the equivalence is measured, not assumed.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .graph import rank_normalize_array
from .signatures import SCORE_COLUMNS, get_gene_sets

__all__ = ["clean_gene_names", "make_unique", "qc_filter", "normalize_log1p",
           "score_genes", "score_all", "rank_normalized",
           "highly_variable_seurat", "scale_dense", "pca_arpack", "banksy_style_domains"]


def clean_gene_names(var: list[str]) -> tuple[list[str], str]:
    """Strip an Ensembl suffix from names written as "SYMBOL ENSG..." (Thrane et al. 2018).

    Returns the cleaned names and a short note. Names without the pattern are unchanged.
    """
    var = [str(v) for v in var]
    pat = [v.split()[0] for v in var if len(v.split()) == 2 and v.split()[1].startswith("ENSG")]
    if len(pat) >= 0.5 * len(var):
        out = [v.split()[0] if len(v.split()) >= 1 else v for v in var]
        return out, f"gene names 'SYMBOL ENSG...' -> SYMBOL ({len(pat)}/{len(var)})"
    return var, "gene names unchanged"


def make_unique(names: list[str], join: str = "-") -> list[str]:
    """anndata's var_names_make_unique: later duplicates get -1, -2, ..."""
    seen: dict[str, int] = {}
    out = []
    taken = set(names)
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
            continue
        k = seen[n]
        while True:
            k += 1
            cand = f"{n}{join}{k}"
            if cand not in taken:
                break
        seen[n] = k
        taken.add(cand)
        out.append(cand)
    return out


def qc_filter(counts: np.ndarray, var: list[str], min_counts: float = 500, min_cells: int = 3):
    """filter_cells(min_counts) then filter_genes(min_cells), as in run_01_qc.

    Returns (counts_kept, spot_mask, gene_mask, var_kept).
    """
    C = np.asarray(counts, dtype=np.float32)
    keep_s = C.sum(axis=1) >= min_counts
    C = C[keep_s]
    keep_g = (C > 0).sum(axis=0) >= min_cells
    C = C[:, keep_g]
    return C, keep_s, keep_g, [v for v, k in zip(var, keep_g) if k]


def normalize_log1p(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """normalize_total(target_sum) followed by natural log1p (float32, as Scanpy)."""
    C = np.asarray(counts, dtype=np.float32)
    tot = C.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    X = C / tot * np.float32(target_sum)
    return np.log1p(X).astype(np.float32)


def score_genes(X: np.ndarray, var_names, gene_list, *, ctrl_size: int = 50, n_bins: int = 25,
                random_state: int = 0) -> np.ndarray:
    """Re-implementation of ``scanpy.tl.score_genes`` (ctrl_as_ref=True, legacy global seed)."""
    np.random.seed(random_state)
    var_names = pd.Index([str(v) for v in var_names])
    gl = pd.Index(list(gene_list))
    gl = gl.intersection(var_names)
    if len(gl) == 0:
        raise ValueError("No valid genes were passed for scoring.")
    pool = var_names.astype("string")
    obs_avg = pd.Series(np.nanmean(X, axis=0, dtype=np.float64), index=pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    control = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gl]):
        r = obs_cut[obs_cut == cut].index
        if ctrl_size < len(r):
            r = r.to_series().sample(ctrl_size).index
        r = r.difference(gl)
        control = control.union(r)
    idx_l = var_names.get_indexer(gl)
    idx_c = var_names.get_indexer(control)
    return (np.nanmean(X[:, idx_l], axis=1, dtype=np.float64)
            - np.nanmean(X[:, idx_c], axis=1, dtype=np.float64))


def score_all(X: np.ndarray, var_names, tumor_type: str, hypoxia_genes=None,
              min_genes: int = 2, ctrl_size: int = 50, seed: int = 0) -> tuple[dict, list]:
    """Same contract as ``signatures.score_all``: raw scores of every signature with enough genes."""
    sets = get_gene_sets(tumor_type, hypoxia_genes)
    vs = set(str(v) for v in var_names)
    out, skipped = {}, []
    for name, genes in sets.items():
        present = [g for g in genes if g in vs]
        if len(present) < min_genes:
            skipped.append((name, len(present), len(genes)))
            continue
        out[name] = score_genes(X, var_names, present, ctrl_size=ctrl_size, random_state=seed)
    return out, skipped


def rank_normalized(scores: dict) -> dict:
    """``<key>_n`` columns for the model signatures, as ``signatures.rank_normalize``."""
    return {f"{k}_n": rank_normalize_array(np.asarray(scores[k], float))
            for k in SCORE_COLUMNS if k in scores}


# ----------------------------------------------------------------------------------------------
# BANKSY-style domains (run_13 / run_13b / run_34 recipe) without Scanpy
# ----------------------------------------------------------------------------------------------
def _mean_var(X: np.ndarray):
    """Column mean and unbiased variance as ``scanpy.pp._utils._get_mean_var`` (float64 result)."""
    mean = X.mean(axis=0, dtype=np.float64)
    mean_sq = np.multiply(X, X).mean(axis=0, dtype=np.float64)
    var = mean_sq - mean ** 2
    n = X.shape[0]
    if n > 1:
        var *= n / (n - 1)
    return mean, var


def highly_variable_seurat(X: np.ndarray, n_top_genes: int = 2000, n_bins: int = 20,
                           log1p_base=None) -> np.ndarray:
    """Boolean mask of ``sc.pp.highly_variable_genes(flavor="seurat", n_top_genes=...)``.

    X is the log1p-normalised (cells x genes) matrix. Dispersions are computed on expm1(X),
    log-transformed, normalised within 20 equal-width bins of log1p(mean) and the top
    ``n_top_genes`` normalised dispersions are kept (ties at the cut-off are all kept)."""
    X = np.array(X, copy=True)
    if log1p_base is not None:
        X *= np.log(log1p_base)
    np.expm1(X, out=X)
    mean, var = _mean_var(X)
    mean[mean == 0] = 1e-12
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)
    df = pd.DataFrame({"means": mean, "dispersions": dispersion})
    df["mean_bin"] = pd.cut(df["means"], bins=n_bins)
    stats = df.groupby("mean_bin", observed=True)["dispersions"].agg(avg="mean", dev="std")
    one = stats["dev"].isnull()
    stats.loc[one, "dev"] = stats.loc[one, "avg"]
    stats.loc[one, "avg"] = 0
    stats = stats.loc[df["mean_bin"]].set_index(df.index)
    norm = ((df["dispersions"] - stats["avg"]) / stats["dev"]).to_numpy()
    x = norm[~np.isnan(norm)].copy()
    x[::-1].sort()
    cut = x[min(n_top_genes, x.size) - 1]
    return np.nan_to_num(norm, nan=-np.inf) >= cut


def scale_dense(X: np.ndarray, max_value: float = 10.0) -> np.ndarray:
    """``sc.pp.scale(zero_center=True, max_value)`` on a dense array (dtype preserved)."""
    X = np.array(X, copy=True)
    mean, var = _mean_var(X)
    std = np.sqrt(var)
    std[std == 0] = 1
    X -= mean
    X /= std
    np.clip(X, -max_value, max_value, out=X)
    return X


def pca_arpack(X: np.ndarray, n_comps: int = 20, random_state: int = 0) -> np.ndarray:
    """``sc.tl.pca`` default path for dense input: sklearn PCA with the ARPACK solver."""
    from sklearn.decomposition import PCA
    return PCA(n_components=n_comps, svd_solver="arpack", random_state=random_state).fit_transform(X)


def banksy_pcs(X_log: np.ndarray, n_pcs: int = 20, seed: int = 0, log1p_base=None) -> np.ndarray:
    """PCs used by ``run_13_benchmark_tools.banksy_style_domains`` when X_pca is absent."""
    hv = highly_variable_seurat(X_log, n_top_genes=2000, log1p_base=log1p_base)
    Z = scale_dense(np.asarray(X_log)[:, hv], max_value=10)
    return pca_arpack(Z, n_comps=n_pcs, random_state=seed)


def banksy_style_domains(X_pca: np.ndarray, A, n_domains: int, lam: float, n_pcs: int = 20,
                         seed: int = 0) -> np.ndarray:
    """k-means on [sqrt(1-lam) * PCs, sqrt(lam) * neighbour-mean PCs] (same as run_13)."""
    from sklearn.cluster import KMeans
    X = np.asarray(X_pca[:, :n_pcs], float)
    deg = np.maximum(np.asarray(A.sum(axis=1)).ravel(), 1.0)
    Xn = (A @ X) / deg[:, None]
    Z = np.hstack([np.sqrt(1 - lam) * X, np.sqrt(lam) * Xn])
    return KMeans(n_clusters=n_domains, n_init=10, random_state=seed).fit_predict(Z)
