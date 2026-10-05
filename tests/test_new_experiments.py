"""Smoke checks for manuscript result JSONs (skipped before real-data runs)."""
from __future__ import annotations

import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "results" / "validation"
REQUIRED = (
    "null_crosslink_check.json",
    "stromal_intervention.json",
    "naive_baseline.json",
)


def _load_outputs() -> dict[str, dict]:
    missing = [name for name in REQUIRED if not (VALIDATION / name).is_file()]
    if missing:
        raise unittest.SkipTest(
            "real-data validation outputs are not present: " + ", ".join(missing)
        )
    return {
        name: json.loads((VALIDATION / name).read_text(encoding="utf-8"))
        for name in REQUIRED
    }


def test_null_crosslink_records_the_permutation_design():
    outputs = _load_outputs()
    result = outputs["null_crosslink_check.json"]
    assert result["seed"] == 20261001
    assert result["n_perm"] == 50
    # rerun on 2026-10-03 across all 19 primary sections (originally 6 representative sections)
    assert result["summary"]["n_slides"] == 19


def test_stromal_intervention_covers_all_admitted_sections():
    outputs = _load_outputs()
    result = outputs["stromal_intervention.json"]
    with (ROOT / "data" / "ledger.csv").open(encoding="utf-8", newline="") as f:
        import csv
        # the three external sections (cancer_type "other") are ingested too but analysed separately
        admitted = {
            row["slide_id"] for row in csv.DictReader(f)
            if row.get("status") == "ingested" and row.get("cancer_type") != "other"
        }
    assert len(admitted) == 19
    assert set(result["per_slide"]) == admitted


def test_stromal_intervention_and_naive_baseline_smoke_summaries():
    outputs = _load_outputs()
    intervention = outputs["stromal_intervention.json"]["summary"]
    naive = outputs["naive_baseline.json"]["summary"]
    assert 50.0 <= intervention["median_pct_drop_b_cell_30"] <= 70.0
    assert naive["n_slides"] == 19


if __name__ == "__main__":
    tests = [value for name, value in sorted(globals().items())
             if name.startswith("test_") and callable(value)]
    failed = 0
    skipped = 0
    passed = 0
    for test in tests:
        print(f"\n[RUN ] {test.__name__}")
        try:
            test()
            passed += 1
            print(f"[PASS] {test.__name__}")
        except unittest.SkipTest as exc:
            skipped += 1
            print(f"[SKIP] {test.__name__}: {exc}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"[FAIL] {test.__name__}: {type(exc).__name__}: {exc}")
    print(f"\n{passed} passed, {skipped} skipped, {failed} failed")
    raise SystemExit(1 if failed else 0)
