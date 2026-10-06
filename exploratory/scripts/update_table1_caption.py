from pathlib import Path

path = Path(__file__).resolve().parents[2] / "docs_is/manuscript_IS.md"
lines = path.read_text(encoding="utf-8").splitlines()
new_caption = (
    "**Table 1** Transcriptomics cohorts used in the study. The primary and replication cohorts were "
    "analysed as separate families; the 2026 extension cohort was used for external validation and "
    "platform comparisons. Patient relationship was unavailable for 11 cSCC sections in GSE289745, "
    "which were therefore treated as conservative section-level sampling units rather than counted as "
    "distinct patients. The full section-level list, including raw and retained spot counts, is given in "
    "Online Resource 1 (Table S1)."
)
for i, line in enumerate(lines):
    if line.startswith("**Table 1**"):
        lines[i] = new_caption
        break
else:
    raise SystemExit("Table 1 caption not found")
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
