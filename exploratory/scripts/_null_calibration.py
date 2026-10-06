"""空模型校准：在已知无关联的数据上测假阳性率。"""
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import scanpy as sc
from scipy.stats import spearmanr

from sparta.io_ import Paths, load_config, load_graph
from sparta.barrier import (compute_b_cell_field, compute_b_mab, compute_b_meta,
                             scores_from_adata)
import scipy.sparse as sp


def _residualize(x, c):
    """OLS residual of x on c (1D)."""
    X = np.column_stack([np.ones_like(c), c])
    beta, *_ = np.linalg.lstsq(X, x, rcond=None)
    return x - X @ beta


def partial_spearman(x, y, control):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(control)
    x, y, c = x[ok], y[ok], control[ok]
    rx, ry = _residualize(x, c), _residualize(y, c)
    rho, _ = spearmanr(rx, ry)
    return float(rho)


def normalized_laplacian(A):
    A = A.tocsr()
    n = A.shape[0]
    deg = np.asarray(A.sum(axis=1)).flatten()
    deg_safe = np.where(deg > 0, deg, 1.0)
    d_inv_sqrt = 1.0 / np.sqrt(deg_safe)
    D_inv_sqrt = sp.diags(d_inv_sqrt)
    return sp.eye(n) - D_inv_sqrt @ A @ D_inv_sqrt


cfg = load_config(None)
P = Paths(cfg)

sid = "CSCC01"
adata = sc.read_h5ad(P.scored(sid))
A, D, source, sink, vessel, _meta = load_graph(P.graph(sid))

filled = []
S = scores_from_adata(adata, missing_out=filled)
dv = compute_b_meta(A, S["hypoxia"], S["proliferation"], S["efflux"],
                    vessel, D=D, **cfg["barrier"]["b_meta"])["d_vessel_um"]
bc = compute_b_cell_field(A, S["ecm"], S["caf"], source,
                          **cfg["barrier"]["b_cell"])["b_cell_field"]
bm = compute_b_mab(A, S["ecm"], S["crosslink"], S["ag_target"], vessel,
                   **cfg["barrier"]["b_mab"])["b_mab"]

n = A.shape[0]
print(f"{sid}: {n} nodes", flush=True)

L = normalized_laplacian(A)
eigenvalues, V = np.linalg.eigh(L.toarray())

f_hat_bm_abs = np.abs(V.T @ bm)
f_hat_bc_abs = np.abs(V.T @ bc)

rng = np.random.default_rng(42)
N_CALIB = 200
N_PERM = 100
ps = []

for calib_i in range(N_CALIB):
    signs_bc = rng.choice([-1.0, 1.0], size=n)
    bc_null = V @ (signs_bc * f_hat_bc_abs)

    null_rhos = []
    for _ in range(N_PERM):
        signs_bm = rng.choice([-1.0, 1.0], size=n)
        bm_null = V @ (signs_bm * f_hat_bm_abs)
        null_rhos.append(partial_spearman(bc_null, bm_null, dv))
    null_rhos = np.array(null_rhos)

    obs = partial_spearman(bc_null, bm, dv)
    k_ge = int((null_rhos >= obs).sum())
    ps.append((k_ge + 1) / (N_PERM + 1))
    if (calib_i + 1) % 25 == 0:
        print(f"  calib {calib_i+1}/{N_CALIB}: FPR={np.mean(np.array(ps)<0.05):.3f}", flush=True)

ps = np.array(ps)
print(f"\n=== Null calibration on {sid} ===")
print(f"N calib={N_CALIB}, N perm={N_PERM}")
print(f"FPR@0.05: {np.mean(ps<0.05):.3f} (expect ~0.05)")
print(f"FPR@0.01: {np.mean(ps<0.01):.3f} (expect ~0.01)")
print(f"Mean p: {ps.mean():.3f} (expect ~0.5)")

import json
out = {"slide": sid, "n_calib": N_CALIB, "n_perm": N_PERM,
       "fpr_05": float(np.mean(ps<0.05)), "fpr_01": float(np.mean(ps<0.01)),
       "mean_p": float(ps.mean())}
Path(r"D:\sparta\results\validation\null_calibration.json").write_text(
    json.dumps(out, indent=2), encoding="utf-8")
print("Wrote null_calibration.json")
