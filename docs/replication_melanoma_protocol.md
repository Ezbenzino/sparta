# Pre-specified plan — melanoma replication cohort (Thrane et al. 2018)

Written 2026-10-05 01:20 (Asia/Shanghai) **before any graph, barrier or association was computed
for these sections**. Earlier the same day the eight count tables had been ingested on the author's
machine, but every section failed admission criterion C6 because the gene names were stored as
"SYMBOL ENSG…" strings and no endothelial marker matched; nothing downstream had been run.

## Data

Thrane K, Eriksson H, Maaskola J, Hansson J, Lundeberg J (2018) Spatially resolved transcriptomics
enables dissection of genetic heterogeneity in stage III cutaneous malignant melanoma. *Cancer Res*
78:5970–5979 (doi:10.1158/0008-5472.CAN-18-0747). Lymph-node metastases from four patients with
stage III cutaneous melanoma, fresh frozen; two consecutive sections per metastasis (eight sections);
first-generation ST arrays with 1,007 spots of 100 µm diameter at 200 µm centre-to-centre spacing.
Treatment before sampling is not reported.

## Fixed pipeline (no parameter is chosen for this cohort)

- Gene names: the Ensembl suffix is stripped ("ANXA2 ENSG00000182718" → "ANXA2"); duplicates are
  made unique as AnnData does.
- QC, normalisation, signature scoring and rank normalisation exactly as for the primary cohort
  (`sparta/lite.py`, a Scanpy-free re-implementation verified to reproduce the archived scores of
  three primary first-generation sections to float32 precision; melanoma marker set; HALLMARK_HYPOXIA).
- Coordinates, graph and compartments as for the primary first-generation ST sections
  (`configs/cscc_legacy_st.yaml`: 200 µm pitch, 300 µm radius; the same quantile thresholds).
- Operators and parameters unchanged (`configs/default.yaml`).

## Admission (decided now, recorded in `configs/default.yaml` as `mel_thrane2018`)

- Minimum spot number 200 instead of 300: the 2016 array has 1,007 spots (the primary first-generation
  sections used a 1,933-spot array), and the published sections average 286 tissue spots.
- C3 (paired H&E) and C7 (treatment status) are waived at cohort level and disclosed; neither
  enters any computation.
- All other criteria at their defaults (median UMI ≥ 1,500; ≥ 20 spots with endothelial signal).

## Analyses (the primary-cohort analyses, unchanged, on this cohort as a separate family)

1. Vessel-distance-adjusted partial Spearman association of the $B_{\mathrm{cell}}$ and
   $B_{\mathrm{mAb}}$ fields; graph-spectral surrogate test (500 surrogates, one-sided); BH across
   the admitted replication sections.
2. Geometry and construction nulls (200 draws each), excess over the construction null.
3. Patient-level inference: exact sign-flip test over the four patients (smallest attainable
   one-sided p = 1/16) and the nested random-effects model; and, as a secondary analysis, the
   same model pooled over the 7 primary and 4 replication patients (11 patients).

Replication is declared if the pooled replication-cohort association is positive with a 95% CI
excluding zero; agreement with the primary cohort's construction-null share is reported
descriptively. Every number is reported whatever its direction.

## Post-plan note

The pooled primary-plus-replication secondary analysis in item 3 was specified as 11 patients
(7 primary + 4 replication) when this plan was written. During the study the primary-cohort
melanoma ledger was corrected from one patient to two (MEL01/MEL04 patient A; MEL02/MEL03
patient B; evidence in `docs/melanoma_mapping_correction.md`), so the primary cohort contains
eight patients and the pooled analysis reported in the manuscript covers 12 patients. No
section-level measurement changed, and the pre-specified replication criterion above, which
concerns the replication cohort alone, is unaffected.
