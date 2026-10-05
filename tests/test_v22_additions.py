"""
Tests for the v2.2 additions (2026-10-05): exact minimum cut, cell-resolution graphs, Scanpy-free scoring.

    python tests/test_v22_additions.py        # no pytest needed
    pytest tests/test_v22_additions.py -v

None of these tests needs Scanpy or real data.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402


def _random_geometric(n=500, r=0.075, seed=0):
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    xy = rng.random((n, 2))
    pairs = np.array(sorted(cKDTree(xy).query_pairs(r)))
    A = sp.coo_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])),
                      shape=(n, n)).tocsr()
    return xy, A, rng


def test_exact_min_cut_capacity_equals_max_flow():
    """The returned cut edges must add up to the maximum flow (a true minimum cut).

    With floating-point preflow-push the residual-graph partition can be off by several per cent;
    the integer-capacity solver must close that gap on graphs with many near-ties."""
    from sparta.barrier import compute_b_cell, edge_pairs
    worst = 0.0
    for seed in range(6):
        xy, A, rng = _random_geometric(seed=seed)
        n = A.shape[0]
        ecm = rng.random(n)
        caf = 0.7 * ecm + 0.3 * rng.random(n)
        source = np.flatnonzero(xy[:, 0] < 0.15)
        sink = np.flatnonzero(np.linalg.norm(xy - 0.6, axis=1) < 0.12)
        out = compute_b_cell(A, ecm, caf, source, sink, return_capacity=True)
        cap = {(int(u), int(v)): c for (u, v), c in zip(out["edge_pairs"], out["edge_capacity"])}
        cut = sum(cap[(min(u, v), max(u, v))] for u, v in out["cut_edges"])
        worst = max(worst, abs(cut - out["max_flow"]) / out["max_flow"])
    assert worst < 1e-9, f"cut capacity differs from the max flow by {worst:.2e}"


def test_exact_min_cut_matches_float_flow_value():
    """The integer solver must not change the maximum-flow value (B_cell) beyond round-off."""
    import networkx as nx
    from sparta.barrier import exact_min_cut
    xy, A, rng = _random_geometric(n=300, seed=11)
    pairs = np.array(sp.triu(A, 1).nonzero()).T
    cap = rng.random(len(pairs)) ** 3
    source, sink = np.arange(0, 10), np.arange(290, 300)
    f_exact, *_ = exact_min_cut(range(300), pairs, cap, source, sink)
    G = nx.DiGraph()
    for (u, v), c in zip(pairs.tolist(), cap):
        G.add_edge(u, v, capacity=float(c))
        G.add_edge(v, u, capacity=float(c))
    for i in source:
        G.add_edge("s", int(i), capacity=1e9)
    for i in sink:
        G.add_edge(int(i), "t", capacity=1e9)
    f_float = nx.maximum_flow_value(G, "s", "t")
    assert abs(f_exact - f_float) <= 1e-9 * max(1.0, f_float)


def test_exact_min_cut_on_archived_section():
    """Regression case: on CSCC08 the floating-point partition was 5.4% away from the max flow.

    Skipped when the archived node table is not present."""
    root = Path(__file__).resolve().parents[1] / "data" / "interim"
    if not (root / "CSCC08.nodes.npz").exists():
        print("  (skipped: CSCC08 node table not present)")
        return
    from sparta.barrier import compute_b_cell
    from sparta.io_ import load_graph
    from sparta.node_tables import load_nodes, scores_from_nodes
    A, D, source, sink, vessel, _ = load_graph(root / "CSCC08.graph.npz")
    S = scores_from_nodes(load_nodes(root / "CSCC08.nodes.npz"))
    out = compute_b_cell(A, S["ecm"], S["caf"], source, sink, return_capacity=True)
    cap = {(int(u), int(v)): c for (u, v), c in zip(out["edge_pairs"], out["edge_capacity"])}
    cut = sum(cap[(min(u, v), max(u, v))] for u, v in out["cut_edges"])
    assert abs(cut - out["max_flow"]) / out["max_flow"] < 1e-9
    assert abs(out["max_flow"] - 30.675315) < 1e-5


def test_delaunay_graph_and_regions():
    """Delaunay edges respect the length cut-off; region density ratio is 0 for equal densities."""
    from sparta.cellgraph import delaunay_graph, largest_component, rank01, region_density_ratio
    rng = np.random.default_rng(2)
    xy = rng.random((400, 2)) * 200.0
    A, D = delaunay_graph(xy, max_len_um=25.0)
    assert D.max() <= 25.0 + 1e-9
    assert A.nnz == D.nnz and (A != A.T).nnz == 0
    assert largest_component(A).sum() > 0.9 * 400
    r = rank01(np.array([3.0, 1.0, 2.0, 2.0]))
    assert np.allclose(r, [1.0, 0.0, 0.5, 0.5])
    inside = np.r_[np.ones(10, bool), np.zeros(20, bool)]
    assert abs(region_density_ratio(inside, 1.0, ~inside, 2.0, pseudo=0.0)) < 1e-12


def test_lite_gene_names_and_scoring():
    """Ensembl suffixes are stripped, duplicates made unique, and score_genes = list mean - control mean."""
    from sparta import lite
    names, note = lite.clean_gene_names(["ANXA2 ENSG00000182718", "VWF ENSG00000110799", "VWF ENSG00000999999"])
    assert names == ["ANXA2", "VWF", "VWF"] and "SYMBOL" in note
    assert lite.make_unique(names) == ["ANXA2", "VWF", "VWF-1"]
    rng = np.random.default_rng(0)
    X = np.log1p(rng.poisson(3.0, size=(50, 400)).astype(np.float32))
    var = [f"G{i}" for i in range(400)]
    sc = lite.score_genes(X, var, ["G1", "G2", "G3"], ctrl_size=50, random_state=0)
    assert sc.shape == (50,) and np.all(np.isfinite(sc))
    # same seed, same result (the control draw is deterministic)
    assert np.allclose(sc, lite.score_genes(X, var, ["G1", "G2", "G3"], ctrl_size=50, random_state=0))
    # removing a constant from every gene leaves the score unchanged
    assert np.allclose(sc, lite.score_genes(X + 1.0, var, ["G1", "G2", "G3"], ctrl_size=50, random_state=0))


def test_lite_banksy_preprocessing():
    """Seurat-flavour HVG, scaling and ARPACK PCA: shapes, clipping, centring and determinism."""
    from sparta import lite
    rng = np.random.default_rng(4)
    lam = rng.gamma(0.5, 2.0, size=500)
    X = np.log1p(rng.poisson(lam, size=(200, 500)).astype(np.float32))
    hv = lite.highly_variable_seurat(X, n_top_genes=100)
    assert hv.dtype == bool and hv.sum() >= 100
    Z = lite.scale_dense(X[:, hv], max_value=10)
    assert Z.dtype == np.float32 and np.abs(Z).max() <= 10 + 1e-6
    assert np.allclose(Z.mean(axis=0)[np.abs(Z).max(axis=0) < 10], 0, atol=1e-4)
    p1 = lite.pca_arpack(Z, n_comps=5, random_state=0)
    p2 = lite.pca_arpack(Z, n_comps=5, random_state=0)
    assert p1.shape == (200, 5) and np.allclose(p1, p2)


def test_pooled_hypergeometric_null():
    """With one section the pooled tail probability is the hypergeometric tail itself."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from scipy.stats import hypergeom
    from run_53_intervention_domain_enrichment import pooled_hypergeom_p
    p_over, p_under = pooled_hypergeom_p([(7, 100, 30, 15)])
    assert abs(p_over - hypergeom.sf(6, 100, 30, 15)) < 1e-12
    assert abs(p_under - hypergeom.cdf(7, 100, 30, 15)) < 1e-12


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    ok = 0
    for t in tests:
        try:
            t()
            print(f"[PASS] {t.__name__}")
            ok += 1
        except AssertionError as e:
            print(f"[FAIL] {t.__name__}: {e}")
    print("=" * 60)
    print(f"{ok} / {len(tests)} passed")
    sys.exit(0 if ok == len(tests) else 1)
