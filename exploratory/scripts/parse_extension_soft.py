"""Parse GEO family SOFT files and summarize sample-level metadata.

Input:  downloads/GSE*_family.soft.gz
Output: plain-text report to stdout (redirected by the caller when desired)

This is a diagnostic script for the 2026 extension-cohort ingestion. It does
not modify raw data or analysis outputs.
"""
from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOWNLOADS = ROOT / "downloads"
GSE_IDS = [
    "GSE200278",
    "GSE289745",
    "GSE300445",
    "GSE316760",
    "GSE320041",
    "GSE321832",
]


def parse_soft(path: Path) -> list[dict[str, object]]:
    """Return sample records from a GEO family SOFT file."""
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"accession": line.split("=", 1)[1].strip(), "fields": {}}
                continue
            if current is None or not line.startswith("!Sample_"):
                continue
            key, _, value = line.partition(" = ")
            fields = current["fields"]
            assert isinstance(fields, dict)
            fields.setdefault(key, []).append(value)
    if current is not None:
        samples.append(current)
    return samples


def main() -> None:
    for gse in GSE_IDS:
        path = DOWNLOADS / f"{gse}_family.soft.gz"
        print(f"\n===== {gse} =====")
        for record in parse_soft(path):
            fields = record["fields"]
            print(f"\n{record['accession']}")
            for key in (
                "!Sample_title",
                "!Sample_source_name_ch1",
                "!Sample_organism_ch1",
                "!Sample_characteristics_ch1",
                "!Sample_library_strategy",
                "!Sample_instrument_model",
            ):
                values = fields.get(key, [])
                if values:
                    print(f"  {key}: {' | '.join(values)}")
            supp = fields.get("!Sample_supplementary_file", [])
            if supp:
                print(f"  supplementary files: {len(supp)}")
                for value in supp:
                    print(f"    - {value}")


if __name__ == "__main__":
    main()
