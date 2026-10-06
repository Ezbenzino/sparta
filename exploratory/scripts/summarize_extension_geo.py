"""Summarize series-level and sample-level GEO metadata for extension datasets."""
from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOWNLOADS = ROOT / "downloads"
GSES = ["GSE200278", "GSE289745", "GSE300445", "GSE316760", "GSE320041", "GSE321832"]
SERIES_KEYS = ("!Series_title", "!Series_summary", "!Series_overall_design")
SAMPLE_KEYS = ("!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1")


def main() -> None:
    for gse in GSES:
        path = DOWNLOADS / f"{gse}_family.soft.gz"
        text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
        print(f"\n===== {gse} =====")
        for line in text.splitlines():
            if line.startswith(SERIES_KEYS):
                print(line)
        for block in text.split("^SAMPLE = ")[1:]:
            lines = block.splitlines()
            print(f"\n{lines[0]}")
            for line in lines:
                if line.startswith(SAMPLE_KEYS):
                    print(" ", line)


if __name__ == "__main__":
    main()
