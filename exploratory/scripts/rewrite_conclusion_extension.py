from pathlib import Path

path = Path(__file__).resolve().parents[2] / "docs_is/manuscript_IS.md"
lines = path.read_text(encoding="utf-8").splitlines()
conclusion = (
    "SPARTA provides two complementary graph structural operators for spatial omics, a minimum cut for graph-defined cellular separation and a screened diffusion–absorption field for a molecular-sized in silico probe. Across 77 transcriptomics sections from primary, replication and heterogeneous public extension cohorts, the operators showed reproducible spatial associations that were partly explained by construction and partly reduced by removal of shared matrix inputs. The minimum cut also aligned with withheld CD8+ T-cell positions in independent CODEX images. These results support SPARTA as a reproducible framework for characterizing spatial tissue architecture and comparing structural hypotheses; they do not establish treatment response, prognosis, therapeutic efficacy or clinical utility."
)
for i, line in enumerate(lines):
    if line.startswith("SPARTA provides two complementary graph structural operators"):
        lines[i] = conclusion
        break
else:
    raise SystemExit("conclusion not found")
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
