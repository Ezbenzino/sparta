"""
Spatial statistics shared by the IS-revision scripts (run_36 ... run_40)
========================================================================
Only numpy / scipy.  Nothing in here imports scanpy, anndata or h5py, so every
statistic can be recomputed from the archived node tables and graph files.

Functions mirror scripts/run_27_spatial_null.py exactly (same rank transform,
same quadratic rank-space adjustment, same Spearman on residuals), but are
vectorised so that thousands of surrogate or simulated fields can be scored at
once.  ``partial_spearman`` reproduces run_27's statistic to machine precision.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

__all__ = [
    "rankz", "control_projector", "partial_spearman", "partial_spearman_many",
    "normalized_laplacian", "spectral_basis", "spectral_sign_null",
    "morans_i", "smooth_random_fields", "naive_spearman_p", "normal_scores",
]


def rankz(v: np.ndarray, axis: int = 0) -> np.ndarray:
    """Average-tie ranks, centred and scaled by the population SD (ddof=0)."""
    r = rankdata(v, axis=axis)
    r = r - r.mean(axis=axis, keepdims=True)
    sd = r.std(axis=axis, keepdims=True)
    return r / np.where(sd > 0, sd, 1.0)


def control_projector(control: np.ndarray) -> np.ndarray:
    """Hat matrix of the quadratic rank-space adjustment [1, z_c, z_c^2]."""
    zc = rankz(np.asarray(control, float))
    X = np.column_stack([np.ones_like(zc), zc, zc ** 2])
    return X @ np.linalg.pinv(X)


def _resid(Zv: np.ndarray, H: np.ndarray) -> np.ndarray:
    return Zv - H @ Zv


def partial_spearman(x, y, control) -> float:
    """Spearman correlation of rank-space residuals (identical to run_27)."""
    x, y, c = (np.asarray(a, float) for a in (x, y, control))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    x, y, c = x[ok], y[ok], c[ok]
    H = control_projector(c)
    rx = _resid(rankz(x), H)
    ry = _resid(rankz(y), H)
    a, b = rankz(rx), rankz(ry)
    return float(np.mean(a * b))


def partial_spearman_many(x, Y, control, H=None) -> np.ndarray:
    """Partial Spearman of one fixed field x against every column of Y (n x S)."""
    x = np.asarray(x, float)
    Y = np.asarray(Y, float)
    if Y.ndim == 1:
        Y = Y[:, None]
    H = control_projector(control) if H is None else H
    ax = rankz(_resid(rankz(x), H))
    RY = _resid(rankz(Y, axis=0), H)
    BY = rankz(RY, axis=0)
    return (ax[:, None] * BY).mean(axis=0)


def naive_spearman_p(rho: np.ndarray, n: int) -> np.ndarray:
    """Two-sided p of a Spearman rho treating the n spots as independent (t approx.)."""
    from scipy.stats import t as tdist
    rho = np.clip(np.asarray(rho, float), -0.999999, 0.999999)
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    return 2 * tdist.sf(np.abs(t), n - 2)


def normalized_laplacian(A: sp.spmatrix) -> sp.spmatrix:
    """L = I - D^-1/2 A D^-1/2 (isolated nodes get degree 1, as in run_27)."""
    A = sp.csr_matrix(A, dtype=float)
    deg = np.asarray(A.sum(axis=1)).ravel()
    dis = 1.0 / np.sqrt(np.where(deg > 0, deg, 1.0))
    Dm = sp.diags(dis)
    return sp.eye(A.shape[0]) - Dm @ A @ Dm


def spectral_basis(A: sp.spmatrix):
    """Eigen-decomposition of the normalised Laplacian (dense; n <= ~5000)."""
    L = normalized_laplacian(A).toarray()
    w, V = np.linalg.eigh(L)
    return w, V


def spectral_sign_null(field: np.ndarray, V: np.ndarray, n_draw: int, rng) -> np.ndarray:
    """Graph-spectral sign randomisation of ``field`` (n x n_draw surrogates).

    The absolute spectral coefficients |V^T f| are kept, signs are drawn
    independently and uniformly from {-1, +1}.  This keeps the graph power
    spectrum of the field and randomises its alignment with any other field.
    Draw order matches run_27 (one rng.choice(size=n) per surrogate).
    """
    fh = np.abs(V.T @ np.asarray(field, float))
    n = len(fh)
    S = np.empty((n, n_draw))
    for k in range(n_draw):
        S[:, k] = rng.choice([-1.0, 1.0], size=n) * fh
    return V @ S


def morans_i(x: np.ndarray, A: sp.spmatrix) -> float:
    """Global Moran's I with binary symmetric weights given by A's pattern."""
    W = sp.csr_matrix(A, dtype=float)
    W.data[:] = 1.0
    z = np.asarray(x, float) - np.mean(x)
    num = z @ (W @ z)
    return float(len(z) / W.sum() * num / (z @ z))


def smooth_random_fields(V: np.ndarray, w: np.ndarray, tau: float, n_fields: int, rng,
                         rank_normalise: bool = True) -> np.ndarray:
    """Independent Gaussian fields with a heat-kernel spectrum exp(-tau * lambda).

    Larger tau -> smoother field (higher Moran's I).  When rank_normalise is True
    each column is mapped to (0, 1] ranks / n, matching how signature scores are
    rank-normalised within a section before entering the operators.
    """
    h = np.exp(-tau * np.clip(w, 0, None))
    E = rng.standard_normal((len(w), n_fields))
    F = V @ (h[:, None] * E)
    if rank_normalise:
        F = rankdata(F, axis=0) / F.shape[0]
    return F


def normal_scores(v: np.ndarray) -> np.ndarray:
    """Rank-based inverse-normal (van der Waerden) scores, (rank - 0.5) / n -> Phi^-1.

    Building graph-spectral surrogates from the normal scores of a field, rather than
    from its raw values, matches the surrogate spectrum to the rank structure that the
    Spearman-type statistic actually uses; this keeps the test calibrated for strongly
    skewed fields (run_37, 'skewed' scenario).
    """
    from scipy.stats import norm
    v = np.asarray(v, float)
    return norm.ppf((rankdata(v) - 0.5) / len(v))
