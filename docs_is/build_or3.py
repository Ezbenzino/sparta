"""
Build Online Resource 3: the anonymised code-and-results archive given to reviewers.

Standard library only, so it runs with any Python >= 3.9 (no numpy needed):

    python docs_is/build_or3.py
    python docs_is/build_or3.py --out submission_IS/Online_Resource_3_code_and_results.zip

What it does
  1. collects the package, analysis scripts, configs, tests, ledger, derived per-section files,
     result files and IS figures (raw expression matrices are never included);
  2. replaces the copyright holder in LICENSE with an anonymous placeholder and adds
     README_REVIEWERS.md and MANIFEST.sha256;
  3. scans every included file (inside .npz members too) for identifying strings and refuses to
     write the archive if any are found.

The archive is written next to, never over, existing files: an existing output is renamed *.prev.zip.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOP = "SPARTA_code_and_results"

# (directory, glob) pairs, relative to the repository root
INCLUDE = [
    ("sparta", "*.py"),
    ("scripts", "*.py"),
    ("scripts/is_figures", "*.py"),
    ("configs", "*.yaml"),
    ("tests", "*.py"),
    ("data", "ledger.csv"),
    ("data", "admission_audit.csv"),
    ("data/interim", "*.nodes.npz"),
    ("data/interim", "*.graph.npz"),
    ("data/interim", "*.barrier.npz"),
    ("data/interim", "*.mincut.json"),
    ("data/interim", "*.barrier_meta.json"),
    ("data/interim", "*.admission.json"),
    ("results/validation", "*.json"),
    ("results/validation", "*.csv"),
    ("results/intervention", "*.npz"),
    ("docs", "codex_validation_protocol.md"),
    ("docs", "replication_melanoma_protocol.md"),
    ("results/counterfactual", "*.json"),
    ("results/figures/is", "*.pdf"),
    ("results/figures/is/supplement", "*.pdf"),
    (".", "requirements.txt"),
    (".", "pyproject.toml"),
    (".", "environment.yml"),
]
# builders for earlier manuscript versions carry author details and are not needed for review
EXCLUDE_NAMES = {"build_cmpb_docx.py", "build_cmpb_release_candidate.py", "generate_cmpb_figures.py"}

IDENTITY = re.compile(rb"(?i)yize|ezbenzino|lllyz|630031258|hangzhou|gmail\.com|github\.com/ezbenzino")
LICENSE_HOLDER = re.compile(r"(?im)^copyright \(c\) .*$")

README = """# SPARTA: anonymised code and results archive (Online Resource 3)

This archive accompanies the manuscript "SPARTA: graph transport operators and structural null
models for immune-cell and IgG-transport barriers in spatial omics", submitted to
*Interdisciplinary Sciences: Computational Life Sciences*. Author names, affiliations and
repository links have been removed for double-blind review. The software is released under the
MIT licence (LICENSE). Version 2.2.0.

Raw data are not redistributed. They are public: GEO GSE144239 (cSCC), GSE250636 (melanoma), the
10x Genomics Visium datasets Human Breast Cancer Block A Sections 1 and 2 and Human Lymph Node
(Space Ranger 1.1.0), the melanoma lymph-node ST data of Thrane et al. (2018, Cancer Research
78:5970-5979; distributed by the original authors) and the CODEX colorectal-cancer data of Schurch
et al. (2020, Cell 182:1341-1359; Mendeley Data doi:10.17632/mpjzbtfgfr.1).

## Contents

| Path | What it holds |
|---|---|
| `sparta/` | the package; the operators (`barrier.py`, `graph.py`, `counterfactual.py`, `synthetic.py`) use NumPy, SciPy and NetworkX only |
| `sparta/barrier.py` | `exact_min_cut`: minimum cut with integer-scaled capacities (Section S13 of Online Resource 1) |
| `sparta/spatial_stats.py` | partial Spearman, graph-spectral and normal-score surrogates |
| `sparta/node_tables.py` | reader for the derived node tables (no Scanpy or h5py needed) |
| `sparta/cellgraph.py` | cell-resolution graphs for the CODEX analysis (Delaunay graph, core sinks, outcome raster) |
| `sparta/lite.py` | Scanpy-free quality control, scoring and BANKSY-style preprocessing, checked against Scanpy outputs |
| `scripts/run_00 ... run_34` | ingestion, QC, scoring, graph construction, operators and the earlier analyses |
| `scripts/run_35 ... run_53` | analyses added for this submission (see the map below) |
| `scripts/_h5ad_reader.py` | reads the stored .h5ad files through libhdf5 when h5py is absent (used by run_51) |
| `scripts/is_figures/` | every manuscript and supplementary figure, and `facts.py`, which writes every number quoted in the text |
| `configs/default.yaml` | every model parameter (`cscc_legacy_st.yaml` for first-generation ST) |
| `docs/codex_validation_protocol.md`, `docs/replication_melanoma_protocol.md` | the two analysis plans written before the corresponding results, with their amendments |
| `tests/` | unit tests; each file runs with `pytest` or as `python tests/<file>.py` |
| `data/ledger.csv` | section-to-patient mapping (status `ingested` = primary or external, `replication` = replication cohort) |
| `data/interim/*.nodes.npz` | per-spot node tables for the 30 analysed sections: coordinates, signature scores, QC columns and log-normalised expression of the model genes |
| `data/interim/*.graph.npz`, `*.barrier.npz`, `*.mincut.json` | spot graphs, the two model fields and the exact minimum cuts |
| `results/validation/`, `results/counterfactual/` | every per-section result file; `is_manuscript_facts.json` holds every number in the manuscript |
| `results/intervention/` | per-spot single-ablation effects of both operators for all 30 sections |
| `results/figures/is/` | figures as PDF |
| `MANIFEST.sha256` | SHA-256 of every file in this archive |

## Quick check (a few minutes on a laptop CPU)

```bash
python -m venv .venv
.venv/bin/python -m pip install numpy scipy networkx matplotlib pyyaml pandas scikit-learn
.venv/bin/python tests/test_is_revision.py              # statistics, node tables, min-cut
.venv/bin/python tests/test_v22_additions.py            # exact minimum cut, cell graphs, Scanpy-free scoring
.venv/bin/python tests/test_barrier.py                  # operator tests
.venv/bin/python scripts/run_35b_verify_reproduction.py # recomputes the field associations from the node tables
.venv/bin/python scripts/is_figures/facts.py            # regenerates every manuscript number from the result files
```

On Windows use `.venv\\Scripts\\python.exe` in place of `.venv/bin/python`. Scripts that start from
the raw `.h5ad` files (`run_00b` to `run_04`) additionally need Scanpy (`requirements.txt` lists the
exact versions used); none of the analyses added for this submission does.

## Where each result comes from

| Manuscript item | Script | Output in `results/validation/` |
|---|---|---|
| Table 1, Fig. S1 | `run_40_graph_coverage_table.py` | `graph_coverage_table.json` |
| Fig. 2, Table 2 | `run_38_synthetic_benchmark.py`, `run_38b_synthetic_s2.py` | `synthetic_benchmark.json`, `synthetic_s2_thickness.json` |
| Fig. 3, Fig. S2, Table S4 | `run_37_null_calibration.py` | `null_calibration_all.json`, `null_calibration_replication.json` |
| Fig. 4a | `run_27_spatial_null.py`, `run_43_nscore_spatial_null.py` | `spatial_null_check.json`, `spatial_null_nscore.json` |
| Fig. 4b, Table 3 | `run_36_patient_level.py` | `patient_level_inference.json` |
| Fig. 4c | `run_41_adjustment_robustness.py` | `adjustment_robustness.json` |
| Fig. 5 | `run_39_coupling_decomposition.py`, `run_14_shared_ecm_check.py` | `geometry_null.json`, `shared_input_spatial_null.json`, `shared_ecm_check.json` |
| Fig. 6 | `run_42_size_exclusion_scan.py`, `run_28_s2_matched.py`, `run_05_counterfactual.py` | `size_exclusion_scan.json`, `s2_matched_selection.json`, `results/counterfactual/` |
| Replication cohort, Table S7 | `run_48_replication_cohort.py` | `replication_melanoma.json`, `lite_scoring_equivalence.json` |
| Fig. 7, Table S8 | `run_47_codex_validation.py` | `codex_validation.json`, `codex_cores.csv` |
| Fig. 8, Table S10 | `run_52_simple_baselines_spearman.py` (with `run_38`, `run_47`) | `simple_baselines_spearman.json` |
| Fig. S7, Table S9 | `run_44_intervention_ranking.py`, `run_49_intervention_targeting.py`, `run_53_intervention_domain_enrichment.py` | `intervention_ranking.json`, `intervention_targeting.json`, `intervention_domain_enrichment.json` |
| Domain comparison, Fig. S3 | `run_51_domain_comparison_exact_cut.py` | `benchmark_lambda_sensitivity.json`, `ndomains_sensitivity.json`, `domain_comparison_exact_cut.json` |
| Exact minimum cuts (Section S13) | `run_50_exact_mincut_refresh.py` | `mincut_exactness.json` |
| Fig. S5 | `run_35b_verify_reproduction.py` | `reproduction_check.json` |
| External sections | `run_29_ext_validation.py`, `run_30_s2_matched_ext.py` | `ext_validation.json`, `s2_matched_ext.json` |

All random seeds are fixed in the scripts or in `configs/default.yaml` and are recorded in the
`meta` block of each result file. Each figure script reads only files in this archive:

```bash
python scripts/is_figures/fig2_synthetic.py      # likewise fig1 ... fig8, figS_supplement.py, figS7_intervention.py
```

Recomputing the calibration (`run_37`, about 45 min on one CPU core), the simulation benchmark
(`run_38`, about 10 min) and the CODEX analysis (`run_47`, about 7 min; needs the public CODEX table)
is optional; their stored outputs are what the figures use.
"""


def collect(root: Path) -> list[Path]:
    files: list[Path] = []
    for d, pat in INCLUDE:
        base = root / d
        if not base.exists():
            continue
        for p in sorted(base.glob(pat)):
            if p.is_file() and p.name not in EXCLUDE_NAMES and "__pycache__" not in p.parts:
                files.append(p)
    return files


def scan(name: str, data: bytes) -> list[str]:
    hits = [m.group(0).decode("latin-1") for m in IDENTITY.finditer(data)]
    if name.endswith(".npz"):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for member in z.namelist():
                hits += [f"{member}:{m.group(0).decode('latin-1')}" for m in IDENTITY.finditer(z.read(member))]
    return hits


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--out", default=str(ROOT / "submission_IS" / "Online_Resource_3_code_and_results.zip"))
    a = ap.parse_args(argv)
    root, out = Path(a.root), Path(a.out)
    files = collect(root)
    payload: dict[str, bytes] = {}
    problems: list[str] = []
    for p in files:
        rel = p.relative_to(root).as_posix()
        data = p.read_bytes()
        hits = scan(rel, data)
        if hits:
            problems.append(f"{rel}: {sorted(set(hits))[:5]}")
        payload[rel] = data
    lic = root / "LICENSE"
    if lic.exists():
        txt = LICENSE_HOLDER.sub("Copyright (c) 2026 The authors (names withheld for double-blind review)",
                                 lic.read_text(encoding="utf-8"))
        payload["LICENSE"] = txt.encode("utf-8")
    payload["README_REVIEWERS.md"] = README.encode("utf-8")
    for rel in ("LICENSE", "README_REVIEWERS.md"):
        hits = scan(rel, payload[rel])
        if hits:
            problems.append(f"{rel}: {hits}")
    if problems:
        print("Identifying strings found; archive NOT written:")
        for line in problems:
            print("  " + line)
        return 1
    manifest = "".join(f"{hashlib.sha256(b).hexdigest()}  {rel}\n" for rel, b in sorted(payload.items()))
    payload["MANIFEST.sha256"] = manifest.encode("utf-8")
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.replace(out.with_suffix(".prev.zip"))
    stamp = datetime.now(timezone.utc).timetuple()[:6]
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for rel, data in sorted(payload.items()):
            info = zipfile.ZipInfo(f"{TOP}/{rel}", date_time=stamp)
            info.compress_type = zipfile.ZIP_STORED if rel.endswith(".npz") else zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            z.writestr(info, data)
    size = out.stat().st_size / 1e6
    n_npz = sum(rel.endswith(".npz") for rel in payload)
    print(f"wrote {out}  ({len(payload)} files, {n_npz} .npz, {size:.1f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
