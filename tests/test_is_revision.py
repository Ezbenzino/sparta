"""
Tests for the statistics added in the IS revision (spatial_stats, node_tables, run_36, run_38).

Dual-mode like the other test files: `pytest tests/` or `python tests/test_is_revision.py`.
None of these tests needs scanpy, anndata or h5py.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from sparta.spatial_stats import (control_projector, normal_scores, partial_spearman,  # noqa: E402
                                  partial_spearman_many, spectral_basis, spectral_sign_null)


def _grid(n=12):
    idx = np.arange(n * n).reshape(n, n)
    rows, cols = [], []
    for i in range(n):
        for j in range(n):
            for di, dj in ((0, 1), (1, 0)):
                a, b = i + di, j + dj
                if a < n and b < n:
                    rows += [idx[i, j], idx[a, b]]
                    cols += [idx[a, b], idx[i, j]]
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n * n, n * n)).tocsr()
    return A


def _reference_partial(x, y, c):
    """The formula of scripts/run_27_spatial_null.py, written out independently."""
    def rz(v):
        r = rankdata(v)
        return (r - r.mean()) / r.std()

    def res(v):
        zv, zc = rz(v), rz(c)
        X = np.column_stack([np.ones_like(zc), zc, zc ** 2])
        b, *_ = np.linalg.lstsq(X, zv, rcond=None)
        return zv - X @ b
    return spearmanr(res(x), res(y))[0]


def test_partial_spearman_matches_reference_with_ties():
    rng = np.random.default_rng(0)
    n = 300
    c = np.round(rng.gamma(2, 50, n), -1)          # heavily tied covariate, like graph distances
    x = c / 100 + rng.normal(size=n)
    Y = 0.4 * x[:, None] + rng.normal(size=(n, 4))
    Y[:, 0] = np.round(Y[:, 0], 1)
    ref = [_reference_partial(x, Y[:, k], c) for k in range(4)]
    got = partial_spearman_many(x, Y, c)
    assert np.allclose(ref, got, atol=1e-12)
    assert abs(partial_spearman(x, Y[:, 1], c) - ref[1]) < 1e-12


def test_spectral_surrogates_keep_power_spectrum():
    A = _grid(8)
    w, V = spectral_basis(A)
    rng = np.random.default_rng(1)
    y = rng.normal(size=A.shape[0])
    S = spectral_sign_null(y, V, 5, rng)
    for k in range(5):
        assert np.allclose(np.abs(V.T @ S[:, k]), np.abs(V.T @ y), atol=1e-9)


def test_normal_scores_are_rank_based():
    rng = np.random.default_rng(2)
    v = rng.lognormal(0, 2, 500)
    z = normal_scores(v)
    assert np.all(np.argsort(z) == np.argsort(v))
    assert abs(z.mean()) < 1e-9 and abs(np.median(z)) < 0.01
    assert np.allclose(normal_scores(np.exp(v)), z)            # invariant to monotone transforms


def test_surrogate_test_is_calibrated_on_a_grid():
    """Independent smooth fields on a grid: rejection rate at 0.05 stays near nominal."""
    A = _grid(12)
    w, V = spectral_basis(A)
    rng = np.random.default_rng(3)
    n = A.shape[0]
    h = np.exp(-8 * np.clip(w, 0, None))
    x = V @ (h * rng.normal(size=n))
    c = np.arange(n, dtype=float) % 12                         # a smooth, tied covariate
    H = control_projector(c)
    rej = 0
    n_sim, n_sur = 120, 60
    for _ in range(n_sim):
        y = V @ (h * rng.normal(size=n))
        r = partial_spearman_many(x, y, c, H)[0]
        fz = np.abs(V.T @ normal_scores(y))
        null = partial_spearman_many(x, V @ (rng.choice([-1.0, 1.0], size=(n, n_sur)) * fz[:, None]), c, H)
        rej += (np.sum(null >= r) + 1) / (n_sur + 1) < 0.05
    assert rej / n_sim < 0.13, rej / n_sim


def test_sign_flip_exact_minimum_p():
    from run_36_patient_level import sign_flip_exact
    r = sign_flip_exact(np.array([0.1, 0.2, 0.15, 0.3, 0.05, 0.25, 0.12]))
    assert r["n_positive"] == 7
    assert abs(r["p_one_sided"] - 1 / 128) < 1e-12
    r2 = sign_flip_exact(np.array([0.1, -0.2, 0.15, -0.3, 0.05, -0.25, 0.12]))
    assert r2["p_one_sided"] > 0.4


def test_nested_reml_recovers_pooled_mean():
    from run_36_patient_level import nested_mean_test
    rng = np.random.default_rng(4)
    est = []
    for _ in range(60):
        ns = [4, 2, 2, 3, 3, 2, 3]
        g = np.repeat(np.arange(7), ns)
        v = rng.uniform(0.04, 0.1, len(g)) ** 2
        y = 0.2 + rng.normal(0, 0.05, 7)[g] + rng.normal(0, 0.03, len(g)) + rng.normal(0, np.sqrt(v))
        f = nested_mean_test(y, g, v)
        est.append(f["mean"])
        assert f["df"] == 6
    assert abs(np.mean(est) - 0.2) < 0.02


def test_node_table_roundtrip_and_neutral_fill():
    from sparta.node_tables import load_nodes, scores_from_nodes
    n = 10
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "X.nodes.npz"
        np.savez_compressed(p, obs__ECM_core_n=np.linspace(0, 1, n), obs__CAF_n=np.linspace(1, 0, n),
                            obs_names=np.array([f"s{i}" for i in range(n)]),
                            obsm__spatial_um=np.zeros((n, 2)),
                            meta=np.array(['{"slide": "X"}'], dtype=object))
        nodes = load_nodes(p)
        filled = []
        S = scores_from_nodes(nodes, missing_out=filled)
    assert np.allclose(S["ecm"], np.linspace(0, 1, n))
    assert "ag_target" in filled and np.all(S["ag_target"] == 0.5)
    assert nodes["meta"]["slide"] == "X"


def test_min_cut_separates_closed_capsule_from_gap():
    import run_38_synthetic_benchmark as B
    from sparta.barrier import compute_b_cell
    xy, A = B.hex_lattice(R=1200)
    rng = np.random.default_rng(5)
    t_closed = B.make_tissue(xy, A, "closed", rng)
    t_gap = B.make_tissue(xy, A, "gap20", rng)
    b_closed = compute_b_cell(A, t_closed.ecm_true, t_closed.caf_true, t_closed.vessel, t_closed.core, **B.CELL)["b_cell"]
    b_gap = compute_b_cell(A, t_gap.ecm_true, t_gap.caf_true, t_gap.vessel, t_gap.core, **B.CELL)["b_cell"]
    assert b_closed > 3 * b_gap


if __name__ == "__main__":
    fns = [v for k, v in dict(globals()).items() if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(fns)}/{len(fns)} passed")
