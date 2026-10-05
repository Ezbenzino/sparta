"""
Loads every stored result used by the IS manuscript into plain Python objects.
Figures and the manuscript text both read from here, so a number printed in the
paper and the number drawn in a figure come from the same JSON field.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "results" / "validation"


def j(name):
    with open(VAL / name, encoding="utf-8") as f:
        return json.load(f)


def ledger():
    with open(ROOT / "data" / "ledger.csv", encoding="utf-8") as f:
        return {r["slide_id"]: r for r in csv.DictReader(f)}


PRIMARY_ORDER = ["CSCC01", "CSCC02", "CSCC03", "CSCC04",          # Visium cSCC (P4, P6)
                 "CSCC05", "CSCC06", "CSCC07",                    # ST P2
                 "CSCC08", "CSCC09", "CSCC10",                    # ST P5
                 "CSCC11", "CSCC12",                              # ST P9
                 "CSCC14", "CSCC15", "CSCC16",                    # ST P10
                 "MEL01", "MEL02", "MEL03", "MEL04"]              # melanoma PtB
EXTERNAL_ORDER = ["BRCA01", "BRCA02", "LN01"]
PATIENT_ORDER = ["CSCC_P4", "CSCC_P6", "CSCC_P2", "CSCC_P5", "CSCC_P9", "CSCC_P10", "MEL_PtB"]
# independent melanoma replication cohort (Thrane et al. 2018; run_48) -- a separate family, never
# pooled into the primary numbers except in the explicitly labelled 11-patient analysis
REPLICATION_ORDER = ["MEL_THR1_rep1", "MEL_THR1_rep2", "MEL_THR2_rep1", "MEL_THR2_rep2",
                     "MEL_THR3_rep1", "MEL_THR3_rep2", "MEL_THR4_rep1", "MEL_THR4_rep2"]
REPLICATION_PATIENTS = ["MEL_THR1", "MEL_THR2", "MEL_THR3", "MEL_THR4"]


def bh(ps):
    import numpy as np
    ps = np.asarray(ps, float)
    m = len(ps)
    order = np.argsort(ps)
    adj = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        i = m - rank
        adj[idx] = min(prev, ps[idx] * m / (i + 1))
        prev = adj[idx]
    return adj
