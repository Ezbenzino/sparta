# Table 1 — Sections analysed

由 `scripts/run_18_table1.py` 从 `data/ledger.csv`、`data/interim/*.graph.npz`
与 `data/interim/*.admission.json` 生成。**不要手工编辑**——队列一变就会不一致，
而且不会报错。重新生成：`python scripts/run_18_table1.py`

生成时间戳与环境见文件末尾。

| Section | Patient | Replicate / site | Cohort | Platform | GEO | Spots (raw) | Median UMI | Nodes after QC | Mean degree | Components | Largest comp. | Source / Sink / Vessel | Admission |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MEL01 | MEL_PtB | sternum | Melanoma | Visium | GSM7983359 | 1,667 | 5,514 | 1,660 | 5.74 | 1 | 100.0 % | 133 / 466 / 332 | admitted |
| MEL02 | MEL_PtB | cecal-nodule | Melanoma | Visium | GSM7983364 | 840 | 3,098 | 801 | 5.35 | 8 | 91.0 % | 65 / 232 / 161 | admitted |
| MEL03 | MEL_PtB | chest-wall | Melanoma | Visium | GSM7983365 | 1,789 | 3,577 | 1,724 | 5.40 | 10 | 96.4 % | 138 / 486 / 345 | admitted |
| MEL04 | MEL_PtB | ribcage | Melanoma | Visium | GSM7983366 | 3,263 | 1,488 | 2,673 | 5.41 | 15 | 98.9 % | 214 / 762 / 535 | admitted |
| CSCC01 | CSCC_P4 | rep1 | cSCC | Visium | GSM4565823 | 744 | 15,714 | 701 | 5.39 | 3 | 99.4 % | 57 / 201 / 141 | admitted |
| CSCC02 | CSCC_P4 | rep2 | cSCC | Visium | GSM4565824 | 696 | 16,686 | 652 | 5.44 | 5 | 98.8 % | 53 / 182 / 131 | admitted |
| CSCC03 | CSCC_P6 | rep1 | cSCC | Visium | GSM4565825 | 3,650 | 982 | 2,645 | 5.29 | 50 | 94.2 % | 212 / 481 / 529 | admitted |
| CSCC04 | CSCC_P6 | rep2 | cSCC | Visium | GSM4565826 | 3,838 | 636 | 2,328 | 4.92 | 54 | 90.0 % | 187 / 358 / 466 | admitted |
| CSCC05 | CSCC_P2 | rep1 | cSCC | 1st-gen ST | GSM4284316 | 666 | 5,388 | 664 | 7.50 | 1 | 100.0 % | 53 / 197 / 133 | admitted |
| CSCC06 | CSCC_P2 | rep2 | cSCC | 1st-gen ST | GSM4284317 | 646 | 5,737 | 643 | 7.46 | 1 | 100.0 % | 52 / 190 / 129 | admitted |
| CSCC07 | CSCC_P2 | rep3 | cSCC | 1st-gen ST | GSM4284318 | 638 | 7,219 | 637 | 7.46 | 1 | 100.0 % | 51 / 187 / 128 | admitted |
| CSCC08 | CSCC_P5 | rep1 | cSCC | 1st-gen ST | GSM4284319 | 590 | 1,956 | 566 | 7.22 | 2 | 99.6 % | 46 / 156 / 114 | admitted |
| CSCC09 | CSCC_P5 | rep2 | cSCC | 1st-gen ST | GSM4284320 | 521 | 2,846 | 511 | 7.20 | 1 | 100.0 % | 41 / 142 / 103 | admitted |
| CSCC10 | CSCC_P5 | rep3 | cSCC | 1st-gen ST | GSM4284321 | 521 | 1,648 | 448 | 6.87 | 1 | 100.0 % | 36 / 128 / 90 | admitted |
| CSCC11 | CSCC_P9 | rep1 | cSCC | 1st-gen ST | GSM4284322 | 1,145 | 602 | 672 | 6.25 | 15 | 94.5 % | 54 / 185 / 135 | admitted |
| CSCC12 | CSCC_P9 | rep2 | cSCC | 1st-gen ST | GSM4284323 | 1,071 | 567 | 595 | 6.30 | 16 | 95.1 % | 48 / 161 / 119 | admitted |
| CSCC13 | CSCC_P9 | rep3 | cSCC | 1st-gen ST | GSM4284324 | 1,182 | 290 | — | — | — | — | — | **rejected — C5** |
| CSCC14 | CSCC_P10 | rep1 | cSCC | 1st-gen ST | GSM4284325 | 608 | 675 | 370 | 5.94 | 3 | 98.6 % | 30 / 106 / 74 | admitted |
| CSCC15 | CSCC_P10 | rep2 | cSCC | 1st-gen ST | GSM4284326 | 621 | 1,817 | 597 | 7.24 | 1 | 100.0 % | 48 / 166 / 120 | admitted |
| CSCC16 | CSCC_P10 | rep3 | cSCC | 1st-gen ST | GSM4284327 | 462 | 1,442 | 380 | 6.75 | 2 | 99.7 % | 31 / 103 / 76 | admitted |

**19 sections from 7 patients admitted; 1 rejected at admission.**

- **Visium** — 8 sections / 3 patients; mean degree 4.92–5.74; largest component 90.0–100.0 % (single component in 1/8)
- **1st-gen ST** — 11 sections / 4 patients; mean degree 5.94–7.50; largest component 94.5–100.0 % (single component in 6/11)
