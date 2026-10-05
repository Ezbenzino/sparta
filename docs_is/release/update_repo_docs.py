"""
One-off repository update for release v2.1.0 (IS submission). Standard library only.

    python update_repo_docs.py <repo-root>

Edits README.md, pyproject.toml, sparta/__init__.py and .gitignore in place (read-modify-write,
each step idempotent), and writes CITATION.cff and .zenodo.json from the copies next to this file.
README numbers are read from results/validation/is_manuscript_facts.json.
"""
from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

IS_SECTION = """## v2.1.0: Interdisciplinary Sciences submission

Release v2.1.0 accompanies the manuscript submitted to *Interdisciplinary Sciences:
Computational Life Sciences*. The v2.0.0 values above remain correct as raw summaries, but the
analyses below change how they should be read. Every number quoted in the manuscript is written
by `scripts/is_figures/facts.py` to `results/validation/is_manuscript_facts.json`.

- **Simulation benchmark with ground truth** (`run_38`, `run_38b`): {syn_n_tissues} simulated
  tissues in {syn_n_geoms} planted-barrier geometries, with access measured by agent-based random
  walks. The minimum cut tracked lost access (Spearman ρ = {syn_sparta_bcell_rho}) and was the only
  summary that separated closed capsules from capsules with a 5% gap (AUC {syn_sparta_bcell_auc_gap05};
  best alternative {syn_other_auc_gap05_max}). The selection-matched gap experiment (S2) gave ratios
  near 1 on planted closed capsules, so it is reported as descriptive only.
- **Calibration on every real graph** (`run_37`): the graph-spectral surrogate test held the
  false-positive rate near 5% on all {cal_n_graphs} section graphs (pooled {cal_pooled_spec}), whereas a
  spot-level test reached {cal_pooled_naive}. For strongly skewed fields use the normal-score variant
  (`sparta.spatial_stats.normal_scores`; pooled {cal_skew_nscore} against {cal_skew_spec} for the raw spectrum).
- **Patient-level inference** (`run_36`): nested random-effects model (sections within patients),
  exact patient sign-flip test and leave-one-patient-out. Pooled ρ = {pl_mean} (95% CI {pl_lo}–{pl_hi});
  positive in {pl_npos_pat}/{pl_J} patients.
- **Structural nulls** (`run_39`): with the real ECM score kept and the crosslinking and ligand
  inputs replaced by graph-spectral surrogates, the model reproduces a median {dec_share_constr_pct}% of
  the observed coupling. The tissue-specific excess is small (pooled {plex_mean}, 95% CI
  {plex_lo}–{plex_hi}).
- **Robustness** (`run_41`, `run_42`, `run_43`): alternative vessel-distance adjustments (median ρ
  {adj_alt_min}–{adj_alt_max}), a size-exclusion scan over β and a normal-score re-test of the primary
  association ({ns_n_q05}/19 sections BH q < 0.05).
- **Reproduction without Scanpy** (`run_35`, `run_35b`, `sparta/node_tables.py`): the per-spot node
  tables in `data/interim/*.nodes.npz` reproduce the stored field associations (max |Δρ| =
  {repro_max_drho}) without Scanpy, AnnData or h5py.

New statistics live in `sparta/spatial_stats.py` and are covered by `tests/test_is_revision.py`.
The two fields co-vary in {assoc_n_pos}/19 sections and in every patient, but most of that coupling is
produced by shared inputs and graph geometry. Outputs are model-defined; they are not measurements
of antibody exposure, T-cell passage or treatment response.

| Script | Purpose |
|---|---|
| `run_35_export_node_tables.py` | export Scanpy-free per-spot node tables (run once on the machine that holds the `.h5ad` files) |
| `run_35b_verify_reproduction.py` | recompute the field associations from the node tables and compare with the stored results |
| `run_36_patient_level.py` | nested random-effects model, patient sign-flip test, leave-one-patient-out |
| `run_37_null_calibration.py` | false-positive rate of the surrogate tests on every real section graph |
| `run_38_synthetic_benchmark.py`, `run_38b_synthetic_s2.py` | planted-barrier simulation benchmark with agent-based access |
| `run_39_coupling_decomposition.py` | geometry and construction null models for the field coupling |
| `run_40_graph_coverage_table.py` | per-section graph coverage, components and stranded spots |
| `run_41_adjustment_robustness.py` | alternative adjustments for vessel distance |
| `run_42_size_exclusion_scan.py` | fraction of IgG-excluding edges as a function of β |
| `run_43_nscore_spatial_null.py` | primary association re-tested with normal-score surrogates |
| `scripts/is_figures/*.py` | all IS figures (`fig1`–`fig6`, `figS_supplement`) and `facts.py` |

"""

TEST_LINE = "python tests/test_is_revision.py        # 8 tests (IS revision statistics; no Scanpy)\n"

DATA_OLD = "**This repository contains no data.** `data/` and `results/` are gitignored."
DATA_NEW = ("**This repository contains no raw data.** From v2.1.0 the derived per-spot node tables, spot\n"
            "graphs, model fields and minimum-cut records of the 22 analysed sections (`data/interim/*.nodes.npz`,\n"
            "`*.graph.npz`, `*.barrier.npz`, `*.mincut.json`), every result file (`results/validation/`,\n"
            "`results/counterfactual/`) and the IS figures are versioned, so every manuscript number can be\n"
            "regenerated without the raw expression matrices.")

CITE_OLD_START = "The exact code-and-results release archive is being prepared for public deposit."
CITE_NEW = ("Release v2.1.0 (code, node tables, result files and figure scripts) accompanies the manuscript\n"
            "submitted to *Interdisciplinary Sciences: Computational Life Sciences*. Its Zenodo DOI will be added\n"
            "here and in [`CITATION.cff`](CITATION.cff) once the deposit and landing page have been verified.")

GITIGNORE_INTERIM = """!data/interim/*.nodes.npz
!data/interim/*.graph.npz
!data/interim/*.barrier.npz
!data/interim/*.mincut.json
!data/interim/*.barrier_meta.json
!data/interim/*.admission.json
"""
GITIGNORE_RESULTS = """!results/validation/
!results/counterfactual/
!results/figures/
results/figures/*
!results/figures/is/
results/figures/is/**/*.tif
results/figures/is/**/*.eps
"""
GITIGNORE_TAIL = """
# IS submission: manuscript renders and the submission package stay local
docs_is/build/
submission_IS/
"""


def fmt(template: str, facts: dict) -> str:
    def rep(m):
        k = m.group(1)
        if k not in facts:
            raise KeyError(f"fact '{k}' missing from is_manuscript_facts.json")
        return str(facts[k])
    return re.sub(r"\{([A-Za-z0-9_]+)\}", rep, template)


def _read(p: Path):
    """Read without newline translation; return text and the file's dominant line ending."""
    with open(p, encoding="utf-8", newline="") as f:
        s = f.read()
    nl = "\r\n" if s.count("\r\n") > s.count("\n") / 2 else "\n"
    return s, nl


def _write(p: Path, s: str) -> None:
    with open(p, "w", encoding="utf-8", newline="") as f:
        f.write(s)


def _nl(block: str, nl: str) -> str:
    return block.replace("\r\n", "\n").replace("\n", nl)


def update_readme(root: Path, facts: dict) -> None:
    p = root / "README.md"
    s, nl = _read(p)
    if "## v2.1.0: Interdisciplinary Sciences submission" not in s:
        anchor = s.index("## Results reported in the manuscript")
        m = re.compile(r"\r?\n---\r?\n").search(s, anchor)
        cut = m.start() + len(m.group(0)) - len(m.group(0).lstrip("\r\n"))   # start of the '---' line
        s = s[:cut] + nl + _nl(fmt(IS_SECTION, facts), nl) + s[cut:]
    if "python tests/test_is_revision.py" not in s:
        s = s.replace("python tests/test_new_experiments.py", _nl(TEST_LINE, nl) + "python tests/test_new_experiments.py", 1)
        s = re.sub(r"All five scripts run without pytest\.\s+The four synthetic/core scripts contain\s+22 checks;",
                   "All six scripts run without pytest. The five synthetic/core scripts contain" + nl + "30 checks;", s)
    if DATA_OLD in s:
        s = s.replace(DATA_OLD, _nl(DATA_NEW, nl))
    s = s.replace("permutation check for rank-normalised crosslinking (6 representative sections)",
                  "permutation check for rank-normalised crosslinking (all 19 primary sections)")
    if CITE_OLD_START in s:
        i = s.index(CITE_OLD_START)
        j = re.compile(r"\r?\n\r?\n").search(s, i).start()
        s = s[:i] + _nl(CITE_NEW, nl) + s[j:]
    _write(p, s)


def bump_versions(root: Path) -> None:
    for rel, pat, new in [("pyproject.toml", r'(?m)^version = "2\.0\.0"', 'version = "2.1.0"'),
                          ("sparta/__init__.py", r'(?m)^__version__ = "2\.0\.0"', '__version__ = "2.1.0"')]:
        p = root / rel
        s, _ = _read(p)
        _write(p, re.sub(pat, new, s))


def update_gitignore(root: Path) -> None:
    p = root / ".gitignore"
    s, nl = _read(p)
    if "!data/interim/*.nodes.npz" not in s:
        s = re.sub(r"(data/external/\*\r?\n)", lambda m: m.group(1) + _nl(GITIGNORE_INTERIM, nl), s, count=1)
    if "!results/validation/" not in s:
        s = re.sub(r"(!results/\.gitkeep\r?\n)", lambda m: m.group(1) + _nl(GITIGNORE_RESULTS, nl), s, count=1)
    if "submission_IS/" not in s:
        s = s.rstrip("\r\n") + nl + _nl(GITIGNORE_TAIL, nl)
    _write(p, s)


def main(root: Path) -> None:
    facts = json.loads((root / "results/validation/is_manuscript_facts.json").read_text(encoding="utf-8"))
    update_readme(root, facts)
    bump_versions(root)
    update_gitignore(root)
    shutil.copyfile(HERE / "CITATION.cff", root / "CITATION.cff")
    shutil.copyfile(HERE / "zenodo.json", root / ".zenodo.json")
    print("README.md, pyproject.toml, sparta/__init__.py, .gitignore, CITATION.cff, .zenodo.json updated")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
