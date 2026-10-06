import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATTERN = re.compile(
    r"treatment response|therapeutic efficacy|clinical utility|prognosis|prognostic|predict|predicted|"
    r"survival|clinical outcome|response to|efficacy|therapeutic benefit|patient outcome|drug delivery|"
    r"drug exposure|treatment|therapeutic|therapy|response|clinical",
    re.IGNORECASE,
)
files = [
    ROOT / "docs_is/manuscript_IS.md",
    ROOT / "docs_is/online_resource_1.md",
    ROOT / "docs_is/cover_letter.md",
    ROOT / "docs_is/references.md",
]
out = []
for path in files:
    out.append(f"# {path.relative_to(ROOT)}")
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if PATTERN.search(line):
            out.append(f"L{i}: {line}")
    out.append("")
target = ROOT / "results/logs/clinical_term_scan.txt"
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text("\n".join(out), encoding="utf-8")
print(target)
