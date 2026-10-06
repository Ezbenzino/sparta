# Table 1 — Sections analysed

由 `scripts/run_18_table1.py` 从 `data/ledger.csv`、`data/interim/*.graph.npz`
与 `data/interim/*.admission.json` 生成。**不要手工编辑**——队列一变就会不一致，
而且不会报错。重新生成：`python scripts/run_18_table1.py`

| Section | Patient | Replicate / site | Cohort | Platform | GEO | Spots (raw) | Median UMI | Nodes after QC | Mean degree | Components | Largest comp. | Source / Sink / Vessel | Admission |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MEL01 | MEL_PtA | sternum | Melanoma | Visium | GSM7983359 | 1,667 | 5,514 | 1,660 | 5.74 | 1 | 100.0 % | 133 / 466 / 332 | admitted**（forced）** |
| MEL02 | MEL_PtB | cecal-nodule | Melanoma | Visium | GSM7983364 | 840 | 3,098 | 801 | 5.35 | 8 | 91.0 % | 65 / 232 / 161 | admitted**（forced）** |
| MEL03 | MEL_PtB | chest-wall | Melanoma | Visium | GSM7983365 | 1,789 | 3,577 | 1,724 | 5.40 | 10 | 96.4 % | 138 / 486 / 345 | admitted**（forced）** |
| MEL04 | MEL_PtA | ribcage | Melanoma | Visium | GSM7983366 | 3,263 | 1,488 | 2,673 | 5.41 | 15 | 98.9 % | 214 / 762 / 535 | admitted**（forced）** |
| MEL300_0019 | MEL300_0019 | AdjIII-0019 | Melanoma | Visium | GSM9060732 | 2,534 | 34,358 | 2,534 | 11.33 | 3 | 99.9 % | 203 / 754 / 507 | 2026 extension cohort — analysed |
| MEL300_0022 | MEL300_0022 | AdjIII-0022 | Melanoma | Visium | GSM9060733 | 2,762 | 27,524 | 2,762 | 11.54 | 1 | 100.0 % | 221 / 810 / 553 | 2026 extension cohort — analysed |
| MEL300_0113 | MEL300_0113 | AdjIII-0113 | Melanoma | Visium | GSM9060734 | 2,132 | 6,802 | 2,132 | 10.91 | 21 | 97.2 % | 171 / 613 / 427 | 2026 extension cohort — analysed |
| MEL300_0133 | MEL300_0133 | AdjIII-0133 | Melanoma | Visium | GSM9060735 | 4,080 | 37,024 | 4,080 | 11.64 | 2 | 100.0 % | 327 / 1147 / 816 | 2026 extension cohort — analysed |
| MEL316_mel2 | MEL316_mel2 | mel2 | Melanoma | Visium | GSM9459774 | 4,082 | 6,694 | 4,082 | 11.59 | 1 | 100.0 % | 327 / 1153 / 817 | 2026 extension cohort — analysed |
| MEL316_mel3 | MEL316_mel3 | mel3 | Melanoma | Visium | GSM9459775 | 4,933 | 9,652 | 4,933 | 11.65 | 1 | 100.0 % | 395 / 1416 / 987 | 2026 extension cohort — analysed |
| MEL320_WU1340 | MEL320_WU1340 | WU1340 | Melanoma | Visium | GSM9532655 | 4,989 | 28,288 | 4,989 | 11.79 | 1 | 100.0 % | 399 / 1489 / 998 | 2026 extension cohort — analysed |
| MEL320_WU1373 | MEL320_WU1373 | WU1373 | Melanoma | Visium | GSM9532656 | 4,651 | 37,893 | 4,651 | 11.59 | 3 | 96.8 % | 373 / 1381 / 931 | 2026 extension cohort — analysed |
| MEL320_WU1384_1 | MEL320_WU1384 | WU1384_1 | Melanoma | Visium | GSM9532657 | 4,782 | 52,328 | 4,782 | 11.79 | 1 | 100.0 % | 383 / 1428 / 957 | 2026 extension cohort — analysed |
| MEL320_WU1384_2 | MEL320_WU1384 | WU1384_2 | Melanoma | Visium | GSM9532658 | 4,571 | 37,535 | 4,571 | 11.70 | 1 | 100.0 % | 366 / 1345 / 915 | 2026 extension cohort — analysed |
| MEL320_WU1609 | MEL320_WU1609 | WU1609 | Melanoma | Visium | GSM9532660 | 4,989 | 63,900 | 4,989 | 11.80 | 1 | 100.0 % | 399 / 1400 / 998 | 2026 extension cohort — analysed |
| MEL320_WU2130_1 | MEL320_WU2130 | WU2130_1 | Melanoma | Visium | GSM9532662 | 4,251 | 3,564 | 4,251 | 11.55 | 3 | 100.0 % | 341 / 1248 / 851 | 2026 extension cohort — analysed |
| MEL320_WU2130_2 | MEL320_WU2130 | WU2130_2 | Melanoma | Visium | GSM9532663 | 3,771 | 42,248 | 3,771 | 10.73 | 8 | 99.8 % | 302 / 1084 / 755 | 2026 extension cohort — analysed |
| MEL320_WU2415 | MEL320_WU2415 | WU2415 | Melanoma | Visium | GSM9532664 | 4,807 | 60,766 | 4,807 | 11.61 | 1 | 100.0 % | 385 / 1385 / 962 | 2026 extension cohort — analysed |
| MEL320_WU3049 | MEL320_WU3049 | WU3049 | Melanoma | Visium | GSM9532665 | 3,897 | 51,182 | 3,897 | 11.68 | 4 | 99.9 % | 312 / 1168 / 780 | 2026 extension cohort — analysed |
| MEL320_WU3244 | MEL320_WU3244 | WU3244 | Melanoma | Visium | GSM9532666 | 4,867 | 51,462 | 4,867 | 11.70 | 1 | 100.0 % | 390 / 1372 / 974 | 2026 extension cohort — analysed |
| MEL320_YUADD | MEL320_YUADD | YUADD | Melanoma | Visium | GSM9532667 | 482 | 10,758 | 482 | 10.78 | 3 | 99.4 % | 39 / 136 / 97 | 2026 extension cohort — analysed |
| MEL320_YUALT | MEL320_YUALT | YUALT | Melanoma | Visium | GSM9532668 | 2,393 | 6,460 | 2,393 | 9.87 | 12 | 86.9 % | 192 / 712 / 479 | 2026 extension cohort — analysed |
| MEL320_YUBOISE | MEL320_YUBOISE | YUBOISE | Melanoma | Visium | GSM9532669 | 369 | 56,178 | 369 | 9.41 | 12 | 74.8 % | 30 / 111 / 74 | 2026 extension cohort — analysed |
| MEL320_YUMAZO | MEL320_YUMAZO | YUMAZO | Melanoma | Visium | GSM9532670 | 2,355 | 23,859 | 2,355 | 8.72 | 58 | 92.9 % | 189 / 640 / 471 | 2026 extension cohort — analysed |
| MEL320_YUSTE | MEL320_YUSTE | YUSTE | Melanoma | Visium | GSM9532671 | 1,458 | 15,102 | 1,458 | 10.92 | 4 | 99.5 % | 117 / 431 / 292 | 2026 extension cohort — analysed |
| ECM01_rep1 | MPM01 | ECM01_rep1 | Melanoma | Slide-seqV2 | GSM6025946 | 27,325 | 556 | 7,100 | 10.51 | 1 | 100.0 % | 568 / 1988 / 1420 | 2026 extension cohort — analysed |
| ECM01_rep2 | MPM01 | ECM01_rep2 | Melanoma | Slide-seqV2 | GSM6025947 | 24,150 | 581 | 7,088 | 10.57 | 4 | 99.9 % | 567 / 1981 / 1418 | 2026 extension cohort — analysed |
| ECM06 | MPM06 | ECM06 | Melanoma | Slide-seqV2 | GSM6025948 | 25,461 | 284 | 6,152 | 10.06 | 20 | 99.5 % | 493 / 1717 / 1231 | 2026 extension cohort — analysed |
| ECM08 | MPM08 | ECM08 | Melanoma | Slide-seqV2 | GSM6025949 | 30,015 | 338 | 5,984 | 9.82 | 50 | 98.9 % | 479 / 1718 / 1197 | 2026 extension cohort — analysed |
| ECM10 | MPM10 | ECM10 | Melanoma | Slide-seqV2 | GSM6025950 | 37,345 | 412 | 6,112 | 10.14 | 16 | 99.7 % | 489 / 1745 / 1223 | 2026 extension cohort — analysed |
| MBM05_rep1 | MBM05 | MBM05_rep1 | Melanoma | Slide-seqV2 | GSM6025935 | 29,526 | 344 | 6,082 | 10.30 | 13 | 99.8 % | 487 / 1695 / 1217 | 2026 extension cohort — analysed |
| MBM05_rep2 | MBM05 | MBM05_rep2 | Melanoma | Slide-seqV2 | GSM6025936 | 32,146 | 358 | 6,422 | 9.94 | 47 | 98.2 % | 514 / 1846 / 1285 | 2026 extension cohort — analysed |
| MBM05_rep3 | MBM05 | MBM05_rep3 | Melanoma | Slide-seqV2 | GSM6025937 | 5,999 | 433 | 4,066 | 6.39 | 20 | 99.2 % | 326 / 1136 / 815 | 2026 extension cohort — analysed |
| MBM06 | MBM06 | MBM06 | Melanoma | Slide-seqV2 | GSM6025938 | 27,026 | 470 | 4,803 | 9.69 | 61 | 98.3 % | 385 / 1355 / 961 | 2026 extension cohort — analysed |
| MBM07 | MBM07 | MBM07 | Melanoma | Slide-seqV2 | GSM6025939 | 38,460 | 472 | 6,355 | 9.92 | 59 | 97.9 % | 617 / 1739 / 1542 | 2026 extension cohort — analysed |
| MBM08 | MBM08 | MBM08 | Melanoma | Slide-seqV2 | GSM6025940 | 35,054 | 352 | 6,888 | 9.77 | 49 | 98.4 % | 551 / 1980 / 1378 | 2026 extension cohort — analysed |
| MBM11_rep1 | MBM11 | MBM11_rep1 | Melanoma | Slide-seqV2 | GSM6025941 | 27,470 | 280 | 6,325 | 10.13 | 33 | 99.4 % | 506 / 1751 / 1265 | 2026 extension cohort — analysed |
| MBM11_rep2 | MBM11 | MBM11_rep2 | Melanoma | Slide-seqV2 | GSM6025942 | 39,044 | 426 | 6,296 | 10.13 | 31 | 99.4 % | 504 / 1784 / 1260 | 2026 extension cohort — analysed |
| MBM11_rep3 | MBM11 | MBM11_rep3 | Melanoma | Slide-seqV2 | GSM6025943 | 9,366 | 256 | 5,220 | 8.36 | 21 | 98.9 % | 418 / 1434 / 1044 | 2026 extension cohort — analysed |
| MBM13 | MBM13 | MBM13 | Melanoma | Slide-seqV2 | GSM6025944 | 32,656 | 337 | 6,351 | 10.02 | 28 | 99.2 % | 509 / 1769 / 1271 | 2026 extension cohort — analysed |
| MBM18 | MBM18 | MBM18 | Melanoma | Slide-seqV2 | GSM6025945 | 38,418 | 396 | 6,689 | 10.09 | 26 | 99.4 % | 535 / 1927 / 1338 | 2026 extension cohort — analysed |
| MEL_THR1_rep1 | MEL_THR1 | rep1 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 271 | 3,662 | 271 | 7.00 | 2 | 99.6 % | 22 / 79 / 55 | replication cohort — registered override |
| MEL_THR1_rep2 | MEL_THR1 | rep2 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 292 | 4,782 | 292 | 7.14 | 1 | 100.0 % | 24 / 87 / 59 | replication cohort — registered override |
| MEL_THR2_rep1 | MEL_THR2 | rep1 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 382 | 2,722 | 382 | 7.34 | 1 | 100.0 % | 31 / 115 / 77 | replication cohort — registered override |
| MEL_THR2_rep2 | MEL_THR2 | rep2 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 375 | 2,644 | 375 | 7.34 | 1 | 100.0 % | 30 / 113 / 75 | replication cohort — registered override |
| MEL_THR3_rep1 | MEL_THR3 | rep1 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 255 | 4,078 | 255 | 7.15 | 1 | 100.0 % | 21 / 71 / 51 | replication cohort — registered override |
| MEL_THR3_rep2 | MEL_THR3 | rep2 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 282 | 3,952 | 282 | 7.01 | 1 | 100.0 % | 23 / 79 / 57 | replication cohort — registered override |
| MEL_THR4_rep1 | MEL_THR4 | rep1 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 212 | 3,371 | 212 | 7.08 | 1 | 100.0 % | 17 / 64 / 43 | replication cohort — registered override |
| MEL_THR4_rep2 | MEL_THR4 | rep2 | Melanoma | 1st-gen ST | 10.1158/0008-5472.CAN-18-0747 | 248 | 9,326 | 248 | 7.15 | 1 | 100.0 % | 20 / 73 / 50 | replication cohort — registered override |
| BRCA01 | BRCA_P1 | s1 | other | Visium | V1_Breast_Cancer_Block_A_Section_1 | 3,798 | 20,762 | 3,798 | 5.81 | 1 | 100.0 % | 304 / 1092 / 760 | admitted**（forced）** |
| BRCA02 | BRCA_P1 | s2 | other | Visium | V1_Breast_Cancer_Block_A_Section_2 | 3,987 | 18,828 | 3,985 | 5.82 | 1 | 100.0 % | 319 / 1143 / 797 | admitted**（forced）** |
| CSCC01 | CSCC_P4 | rep1 | cSCC | Visium | GSM4565823 | 744 | 15,714 | 701 | 5.39 | 3 | 99.4 % | 57 / 201 / 141 | admitted |
| CSCC02 | CSCC_P4 | rep2 | cSCC | Visium | GSM4565824 | 696 | 16,686 | 652 | 5.44 | 5 | 98.8 % | 53 / 182 / 131 | admitted |
| CSCC03 | CSCC_P6 | rep1 | cSCC | Visium | GSM4565825 | 3,650 | 982 | 2,645 | 5.29 | 50 | 94.2 % | 212 / 481 / 529 | admitted |
| CSCC04 | CSCC_P6 | rep2 | cSCC | Visium | GSM4565826 | 3,838 | 636 | 2,328 | 4.92 | 54 | 90.0 % | 187 / 358 / 466 | admitted |
| CSCC289_S1 | CSCC289_S1 | S1 | cSCC | Visium | GSM8797973 | 2,860 | 5,937 | 2,860 | 11.05 | 1 | 100.0 % | 229 / 846 / 572 | 2026 extension cohort — analysed |
| CSCC289_S10 | CSCC289_S10 | S10 | cSCC | Visium | GSM8797974 | 1,146 | 18,764 | 1,146 | 10.03 | 5 | 54.8 % | 92 / 183 / 230 | 2026 extension cohort — analysed |
| CSCC289_S11 | CSCC289_S11 | S11 | cSCC | Visium | GSM8797975 | 1,483 | 3,378 | 1,483 | 11.04 | 3 | 99.9 % | 119 / 431 / 297 | 2026 extension cohort — analysed |
| CSCC289_S15 | CSCC289_S15 | S15 | cSCC | Visium | GSM8797976 | 2,029 | 7,495 | 2,029 | 11.05 | 10 | 98.6 % | 163 / 601 / 406 | 2026 extension cohort — analysed |
| CSCC289_S3 | CSCC289_S3 | S3 | cSCC | Visium | GSM8797977 | 2,858 | 15,935 | 2,858 | 11.37 | 1 | 100.0 % | 229 / 441 / 572 | 2026 extension cohort — analysed |
| CSCC289_S4 | CSCC289_S4 | S4 | cSCC | Visium | GSM8797978 | 3,195 | 15,021 | 3,195 | 11.43 | 6 | 99.1 % | 256 / 953 / 639 | 2026 extension cohort — analysed |
| CSCC289_S5 | CSCC289_S5 | S5 | cSCC | Visium | GSM8797979 | 2,729 | 7,218 | 2,729 | 11.17 | 30 | 97.6 % | 219 / 810 / 546 | 2026 extension cohort — analysed |
| CSCC289_S6 | CSCC289_S6 | S6 | cSCC | Visium | GSM8797980 | 2,059 | 5,758 | 2,059 | 11.13 | 10 | 99.1 % | 165 / 595 / 412 | 2026 extension cohort — analysed |
| CSCC289_S7 | CSCC289_S7 | S7 | cSCC | Visium | GSM8797981 | 2,746 | 27,858 | 2,746 | 11.47 | 1 | 100.0 % | 220 / 824 / 550 | 2026 extension cohort — analysed |
| CSCC289_S8 | CSCC289_S8 | S8 | cSCC | Visium | GSM8797982 | 2,156 | 17,300 | 2,156 | 11.03 | 14 | 95.2 % | 173 / 647 / 432 | 2026 extension cohort — analysed |
| CSCC289_S9 | CSCC289_S9 | S9 | cSCC | Visium | GSM8797983 | 4,261 | 10,234 | 4,261 | 11.59 | 1 | 100.0 % | 341 / 1209 / 853 | 2026 extension cohort — analysed |
| CSCC321_cut1 | CSCC321_cut1 | cut1 | cSCC | Visium | GSM9550336 | 4,258 | 4,744 | 4,258 | 11.19 | 3 | 99.8 % | 341 / 1231 / 852 | 2026 extension cohort — analysed |
| CSCC321_cut2 | CSCC321_cut2 | cut2 | cSCC | Visium | GSM9550337 | 4,195 | 13,472 | 4,195 | 11.59 | 1 | 100.0 % | 336 / 1232 / 839 | 2026 extension cohort — analysed |
| LN01 | LN_P1 | s1 | other | Visium | V1_Human_Lymph_Node | 4,035 | 20,239 | 4,025 | 5.81 | 7 | 99.7 % | 322 / 1096 / 805 | admitted**（forced）** |
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

**22 ingested sections from 10 patient groups (19 current primary sections from 8 patients plus 3 legacy external-validation sections); 8 replication sections from 4 patients analysed under a registered override; 50 extension sections from 32 identifiable patients (plus 11 sections with unavailable patient relationships) analysed; 1 rejected at admission.**

- **Primary / 1st-gen ST** — 11 sections / 4 patients; mean degree 5.94–7.50; largest component 94.5–100.0 % (single component in 6/11)
- **Primary / Visium** — 11 sections / 6 patients; mean degree 4.92–5.82; largest component 90.0–100.0 % (single component in 3/11)
- **Replication / 1st-gen ST** — 8 sections / 4 patients; mean degree 7.00–7.34; largest component 99.6–100.0 % (single component in 7/8)
- **Extension / Slide-seqV2** — 16 sections / 11 patients; mean degree 6.39–10.57; largest component 97.9–100.0 % (single component in 1/16)
- **Extension / Visium** — 34 sections / 21 identifiable patients plus 11 patient-relationship-unavailable sections; mean degree 8.72–11.80; largest component 54.8–100.0 % (single component in 14/34)

---

生成于 2026-10-06T03:02:23+00:00　seed=0　numpy 2.5.2　scipy 1.18.1
