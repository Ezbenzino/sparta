# CMPB submission readiness review

**Review date:** 4 October 2026  
**Target journal:** *Computer Methods and Programs in Biomedicine* (CMPB)  
**Manuscript:** `docs/manuscript_cmpb_draft.md`

## Editorial position

The work can be presented as a computational method demonstration with explicitly model-defined outputs. The current evidence does not support claims of measured tissue permeability, functional antibody penetration, clinical prediction, or patient-level generalisation. The manuscript now states those limits and treats the external sections and bulk response analysis as exploratory.

This is a content-ready draft, not a submission-ready package yet. The main scientific wording and sensitivity calculations have been reconciled. A versioned, immutable code-and-results archive with a real DOI is still required to close reproducibility and data-availability statements. The manuscript author must also confirm the declarations and final scientific interpretation.

## Claims that are now aligned to the saved analyses

- The primary field association is positive in 18/19 sections; 14/19 pass BH correction under the 19-section graph-spectral surrogate family.
- The `Ag_target`-available subset contains 15 sections. It is positive in 15/15, with 13/15 below BH q<0.05 after correction within this subset. The same stored per-section 500-draw surrogate tests are used; only the multiplicity family changes.
- S2 is the selection-matched analysis from `s2_matched_selection.json`: median ratio 1.053, with 8/19 sections below BH q<0.05. All 19 sections are included because S2 uses ECM/CAF scores and the source–sink graph, not `Ag_target` or `efflux`.
- The old clean14 S2 statistic (1.157; 12/14) came from `paper_stats.json` and the unmatched legacy selection procedure. It has been removed as evidence for the matched analysis.
- Restricting the field association to each largest connected component did not change signs (18/19 positive); the median partial correlation changed from 0.217 to 0.201. This is an input-field sensitivity check, not a new spatial-null test.
- The bulk ICB check is exploratory. No dated pre-outcome protocol has been identified, the two cohorts disagree in direction, and the results do not support a spatial-versus-bulk performance claim.
- The domain comparison uses a BANKSY-style implementation, not the official BANKSY implementation. The Squidpy structure check shares the graph used in label smoothing. Neither is described as an independent head-to-head validation.
- Domain-count and neighborhood-weight scans are now available in the stored results. They show selected parameter sensitivity, not that these settings are optimal or that SPARTA outperforms domain methods.

## Remaining scientific limitations

1. **No functional calibration.** There are no direct measurements of antibody exposure, molecular mesh size, cell migration, or barrier opening. The model's rank-based mesh law and `Ag_target` absorption term are qualitative assumptions.
2. **Few independent patients.** Nineteen sections come from seven people; all four melanoma sections come from one person. Section-level BH tests do not substitute for a patient-clustered model.
3. **Small external arm.** Three sections represent two people; two breast sections are adjacent sections from one patient and the lymph-node sample is non-tumour. The external results are exploratory computational replication only.
4. **Graph and compartment assumptions.** The source, sink, radius graph and largest-component rule shape the result. The MEL02 disconnection check suggests that finite-value imputation has a small effect on the field correlation, but it does not test alternative segmentation, graph radius, or pathology-defined compartments.
5. **Spatial surrogate scope.** The graph-spectral sign-randomisation preserves spectral power, not all spatial statistics. Its per-section p-values are conditional on this surrogate and are not patient-level evidence.
6. **Benchmark scope.** The BANKSY-style / Squidpy analyses are implementation-specific descriptions of complementary outputs. Runtime measurements are machine-dependent. Do not claim superiority, accuracy, or general efficiency from these comparisons.

## CMPB package checks before submission

- Re-open the live CMPB Guide for Authors and Editorial Manager fields immediately before upload. The publisher page was not accessible from this environment during this review, so word limits, required declarations, file types and figure specifications need a final live check.
- Ensure the submission has an editable manuscript, a separate editable Highlights file if required, figure files at the requested resolution, and any supplementary tables/data referenced in the manuscript.
- Verify that the structured abstract, keyword count, title page, corresponding-author details, figure and table citations, references, and legends match current journal instructions.
- Confirm the final author list, affiliations, author contributions, funding, competing interests, ethics wording, data statement, and generative-AI declaration with the author.
- Upload only after the archived analysis version and DOI resolve, and the exact manuscript/figure source corresponds to that archived version.
- Check the cover letter explains CMPB's computational-method contribution and its biomedical scope without claiming clinical readiness or superiority over spatial-domain methods.

## Recommended next scientific investment

The most valuable new evidence would be an independent patient cohort with pathology-confirmed spatial compartments and a direct readout of antibody/tracer penetration or cell migration. More sensitivity plots generated from the same assumed transport law would add less evidence. Until functional validation is available, retain the current model-demonstration framing.

## Source outputs checked

- `results/validation/s2_matched_selection.json`
- `results/validation/spatial_null_check.json`
- `results/validation/metric_connectivity_sensitivity.json`
- `results/validation/benchmark_lambda_sensitivity.json`
- `results/validation/prereadiness_audit.json`
- `scripts/run_28_s2_matched.py`, `scripts/run_27_spatial_null.py`, and `scripts/run_33_input_connectivity_sensitivity.py`
