# SPARTA

SPARTA is a graph-based methodological framework for characterizing spatial tissue architecture and local microenvironmental structure in spatial omics data. It constructs a spatial graph for each tissue section and applies graph structural operators, calibrated spatial null models and perturbation analyses to describe how local tissue fields are organized.

> Chinese guide: [`README.zh.md`](README.zh.md). The English README is the canonical project overview for the current submission.

## Positioning statement

SPARTA is **not** presented as a treatment-response, prognostic, therapeutic-efficacy or clinical-utility method. The IgG-sized molecular field is an in silico structural benchmark; no measured drug exposure is used or estimated. The framework characterizes graph-defined tissue structure and generates analytical hypotheses.

## Repository layout

| Path | Contents |
|---|---|
| `sparta/` | Importable production package |
| `scripts/` | Versioned pipeline and figure-generation scripts |
| `configs/` | Default parameters and cohort admission rules |
| `docs_is/` | Manuscript, cover letter, references and supplement sources |
| `docs/` | Scientific protocol notes, cohort correction record, gene-set sources and manuscript tables |
| `results/` | Validation outputs and manuscript figures |
| `exploratory/` | One-off inspection and drafting work; not part of the clean pipeline |
| `legacy/` | Obsolete scripts retained for traceability; do not run |
| `tests/` | Reproducibility and correctness checks |

## Cohorts

| Cohort | Sections | Patient/analytical units | Platforms and tissue types |
|---|---:|---:|---|
| Primary | 19 | 8 patients | 15 cSCC sections from 6 patients; 4 melanoma sections from 2 patients |
| Replication | 8 | 4 patients | First-generation ST melanoma lymph-node metastases |
| 2026 extension | 50 | 32 identifiable patients plus 11 relationship-unavailable cSCC sections | 34 Visium and 16 Slide-seqV2 sections |
| Combined transcriptomics | 77 | 44 identifiable patients plus 11 conservative section-level units | cSCC, primary melanoma and metastatic melanoma |

Extension disease strata: 13 cSCC, 6 primary melanoma and 31 metastatic melanoma sections. The 11 GSE289745 cSCC sections do not have patient identifiers and are not counted as distinct patients.

## Main extension datasets

| Accession | Platform | Included sections | Stratum |
|---|---|---:|---|
| GSE289745 | Visium | 11 | cSCC; patient relationship unavailable |
| GSE321832 | Visium | 2 | Cutaneous SCC |
| GSE300445 | Visium | 4 | Primary melanoma |
| GSE316760 | Visium | 2 | Primary melanoma |
| GSE320041 | Visium | 15 | Metastatic melanoma |
| GSE200278 | Slide-seqV2 | 16 | Metastatic melanoma |

Exclusions: WU1457 and WU2109 in GSE320041 were excluded because they are Visium HD 16 µm samples; oral1, oral2 and lung1 in GSE321832 were excluded because they are not cutaneous SCC.

## Core outputs

- Extension association: 48/50 sections positive; median partial Spearman ρ = 0.273; patient/analytical-unit pooled estimate 0.28 (95% CI 0.23–0.32).
- Disease-stratified extension pooled estimates: cSCC 0.37; primary melanoma 0.15; metastatic melanoma 0.26.
- Platform medians: Visium 0.320; Slide-seqV2 0.193.
- Extension calibration false-positive rates: graph spectral 0.042; normal score 0.055; point level 0.314.
- Construction and shared-input nulls show that part, but not all, of the association is reproduced by graph construction and shared matrix structure.

These are structural validation results, not clinical predictions.

## Environment

Use the project virtual environment:

```powershell
D:\sparta\.venv\Scripts\python.exe -m pip install -r requirements.txt
D:\sparta\.venv\Scripts\python.exe -m pip install -e .
```

The system Python may not contain project dependencies. On this machine, use:

```powershell
D:\sparta\.venv\Scripts\python.exe
```

## Quality control and platform settings

Primary and public-data QC filters are recorded in `configs/default.yaml`, `results/qc/extension_2026_qc_summary.csv` and the extension README.

Final extension rules:

- Visium spots: total UMI ≥ 500.
- Slide-seqV2 beads: total UMI ≥ 100.
- Genes: detected in at least three retained locations.
- Mitochondrial fraction: ≤ 0.20 when mitochondrial genes are measurable; otherwise marked not assessable.
- Extension Visium sections: ≥ 300 retained spots and retained median UMI ≥ 1500.
- Extension Slide-seqV2 sections: ≥ 1000 retained beads and retained median UMI ≥ 100.

Slide-seqV2 is QC'd at native bead level, then counts are summed into fixed 50 µm grid bins for structural analysis. The Slide-seqV2 graph radius is 100 µm and is not reused from the 150 µm Visium setting.

## Reproducing the work

After data are placed according to the dataset README files, the principal extension modules are:

```powershell
# Primary melanoma mapping correction
D:\sparta\.venv\Scripts\python.exe scripts\run_54_correct_melanoma_mapping.py

# Extension ingestion and QC
D:\sparta\.venv\Scripts\python.exe scripts\run_55_extension_ingest.py

# Extension structural analysis
D:\sparta\.venv\Scripts\python.exe scripts\run_56_extension_cohort.py

# Extension supplement tables
D:\sparta\.venv\Scripts\python.exe scripts\run_57_extension_supplement.py

# Manuscript facts and figures
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\facts.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\fig4_association.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\fig8_baselines.py
D:\sparta\.venv\Scripts\python.exe scripts\is_figures\figS1_parameter_heatmaps.py
```

Random seeds used in the final extension work include `20261005`; the extension pipeline uses fixed surrogate and null counts and writes checkpoints for recovery.

## Tests

```powershell
D:\sparta\.venv\Scripts\python.exe tests\test_barrier.py
D:\sparta\.venv\Scripts\python.exe tests\test_statistics.py
D:\sparta\.venv\Scripts\python.exe tests\test_is_revision.py
D:\sparta\.venv\Scripts\python.exe tests\test_v22_additions.py
```

## Cohort change log

- **2026-10-05, primary melanoma mapping correction.** GSE250636 sample characteristics show MEL01 and MEL04 belong to patient A; MEL02 and MEL03 belong to patient B. The primary cohort is 19 sections from 8 patients. See `docs/melanoma_mapping_correction.md`.
- **2026 extension cohort.** Fifty public sections from six GEO series were ingested, QC'd and analysed separately. See `data/external/extension_2026/README.md`.
- **Manuscript scope revision.** Treatment-response, prognosis, therapeutic-efficacy and clinical-utility statements were removed, and limitations were moved to the beginning of the Discussion.
- **Parameter sensitivity placement.** ξ₀/β/λ sensitivity heatmaps are Section S1, Table S1 and Fig. S1 in Online Resource 1.

## Data availability

Primary transcriptomic data come from public spatial-omics datasets, including GSE144239 and GSE250636. The 2026 extension uses GSE200278, GSE289745, GSE300445, GSE316760, GSE320041 and GSE321832. Raw public data are not redistributed in this repository beyond the terms of their source archives.

## License

MIT. See [`LICENSE`](LICENSE).
