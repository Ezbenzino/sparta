# Melanoma patient-mapping correction

## Summary

The initial primary-cohort ledger assigned all four extracranial melanoma
sections to one patient (`MEL_PtB`). The GSE250636 sample characteristics
identify two patients among those sections: patient A contributed MEL01 and
MEL04, and patient B contributed MEL02 and MEL03.

The primary cohort therefore contains **19 sections from 8 patients**, not 7.
Section-level measurements did not change; patient-level summaries, nested
models, confidence intervals, sign-flip tests and pooled analyses with the
Thrane replication cohort were recomputed.

## Evidence and corrected mapping

| Slide | GEO sample | GEO patient ID | Old ledger patient | Corrected patient |
|---|---|---|---|---|
| MEL01 | GSM7983359 | A | MEL_PtB | MEL_PtA |
| MEL02 | GSM7983364 | B | MEL_PtB | MEL_PtB |
| MEL03 | GSM7983365 | B | MEL_PtB | MEL_PtB |
| MEL04 | GSM7983366 | A | MEL_PtB | MEL_PtA |

## Source

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250636
- Sample characteristics field: `patient id`

## Reproduction

1. Run `python scripts/run_54_correct_melanoma_mapping.py`.
2. Re-run patient-level and table-generation scripts:
   `python scripts/run_40_graph_coverage_table.py`;
   `python scripts/run_36_patient_level.py`;
   `python scripts/run_48_replication_cohort.py --stats-only`;
   `python scripts/run_18_table1.py`;
   `python scripts/is_figures/facts.py`.

Git supplies the exact before/after version of the ledger.
