from pathlib import Path

path = Path(__file__).resolve().parents[2] / "docs_is/manuscript_IS.md"
lines = path.read_text(encoding="utf-8").splitlines()
abstract = (
    "Spatial omics maps stromal and immune programmes within tissue, but simple composition summaries do not describe how those programmes are connected. We present SPARTA, a graph-based methodological framework with two structural operators: a minimum-cut operator for graph-defined separation between vessel-proximal immune sources and tumour-nest sinks, and a screened diffusion–absorption field for an in silico molecular-sized probe. Both operators are accompanied by graph-spectral surrogate tests, construction and geometry null models, parameter perturbations and held-out structural comparisons. We analysed 19 primary spatial transcriptomics sections from eight patients, eight independent melanoma sections from four replication patients, and 50 public extension sections from six GEO series, spanning Visium, first-generation Spatial Transcriptomics and Slide-seqV2; patient relationships were unavailable for 11 extension cSCC sections, which were treated as conservative section-level units. The two fields were positively associated in 18 primary sections and 48 extension sections; the extension cohort had a pooled patient-level estimate of {{ext2026_pl_mean}} (95% CI {{ext2026_pl_lo}}–{{ext2026_pl_hi}}). Structural nulls reproduced part of the association, whereas ablating shared matrix inputs reduced but did not eliminate it. In 114 CODEX colorectal-cancer cores from 35 patients, the minimum-cut field was associated with withheld CD8+ T-cell positions. SPARTA characterizes spatial tissue architecture and operator behaviour across platforms; it does not establish treatment response, prognosis, therapeutic efficacy or clinical utility."
)
for i, line in enumerate(lines):
    if line.startswith("Spatial omics maps stromal and immune programmes"):
        lines[i] = abstract
        break
else:
    raise SystemExit("abstract paragraph not found")
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
