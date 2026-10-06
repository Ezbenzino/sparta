"""Insert parameter sensitivity as Online Resource 1 Section S1 and renumber old items."""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OR = ROOT / "docs_is/online_resource_1.md"
MS = ROOT / "docs_is/manuscript_IS.md"
BUILD_OR = ROOT / "docs_is/build_or1.py"
BUILD_IS = ROOT / "docs_is/build_is.py"


def bump_number(match: re.Match) -> str:
    return f"{match.group(1)}{int(match.group(2)) + 1}{match.group(3) or ''}"


def bump_text(text: str) -> str:
    patterns = [
        r"(# S)(\d+)()",
        r"(\*\*Table S)(\d+)()(\*\*)",
        r"(Table S)(\d+)()",
        r"(Fig\. S)(\d+)([a-z]?)",
        r"(Section S)(\d+)()",
        r"(Sections S)(\d+)()",
        r"(<<FigS)(\d+)([a-z]?)",
    ]
    for pattern in patterns:
        text = re.sub(pattern, bump_number, text)
    return text


or_text = OR.read_text(encoding="utf-8")
or_text = bump_text(or_text)
# Explicit prose ranges and references.
or_text = or_text.replace("Figs. S2 and S4–S7 follow; Fig. S3 accompanies Section S6.",
                          "Figs. S2 and S4–S8 follow; Fig. S3 accompanies Section S6.")
new_section = """# S1 Parameter sensitivity (ξ₀/β/λ)

We examined parameter dependence before interpreting the structural fields. Table S1 lists the grids and evaluation metrics; Fig. S1 shows median values across five representative primary sections (MEL01, MEL03, CSCC01, CSCC03 and CSCC05). The default setting is marked on each heatmap. The qualitative pattern is not confined to a narrow parameter ridge: high crosslink contribution persists over the default β range, while dissociation potential changes smoothly with λ, ξ₀ and β. These parameters remain qualitative effective parameters rather than calibrated physical constants in tissue.

**Table S1** Parameter sensitivity grids and evaluation metrics. Crosslink contribution is the percentage of model-derived edge fraction attributed to crosslinking; dissociation potential is the molecular barrier output before rank standardisation. Heatmap values are medians across the five representative sections.

| Parameter | Symbol | Grid | Default | Assessment |
|---|---:|---:|---:|---|
| ECM decay | λ | 0.5, 1, 2, 3, 5, 8, 12, 20 | 3 | Crosslink contribution; log10 potential |
| Mesh contraction | β | 0.25, 0.5, 1, 1.5, 2, 3, 4.5, 6, 9, 12 | 3 | Crosslink contribution; log10 potential |
| Baseline mesh size (nm) | ξ₀ | 8, 12, 16, 20, 30, 40, 60, 80 | 20 | Crosslink contribution; log10 potential |

![**Fig. S1** Parameter sensitivity heatmaps. **a** Crosslink contribution over λ and β; **b** log10 dissociation potential over λ and β; **c** crosslink contribution over ξ₀ and β; **d** log10 dissociation potential over ξ₀ and β. Black rectangles mark the default setting](<<FigS1_parameter_sensitivity>>){width=100%}

"""
marker = "# S2 Section admission and quality control"
or_text = or_text.replace(marker, new_section + marker)
OR.write_text(or_text, encoding="utf-8")

ms_text = bump_text(MS.read_text(encoding="utf-8"))
ms_text = ms_text.replace(
    "Parameter dependence of the size-exclusion law was examined over grids of $\\beta$, $\\lambda$ and $\\xi_0$, and with a permutation control that shuffled the crosslinking score within sections.",
    "Parameter dependence of the size-exclusion law was examined over grids of $\\beta$, $\\lambda$ and $\\xi_0$ (Online Resource 1, Section S1, Table S1 and Fig. S1), and with a permutation control that shuffled the crosslinking score within sections."
)
ms_text = ms_text.replace(
    "The full section-level list, including raw and retained spot counts, is given in Online Resource 1 (Table S2).",
    "The full section-level list, including raw and retained spot counts, is given in Online Resource 1 (Tables S2 and S12–S13)."
)
MS.write_text(ms_text, encoding="utf-8")

build_or = BUILD_OR.read_text(encoding="utf-8")
build_or = bump_text(build_or)
build_or = build_or.replace(
    'for tag in ("FigS2a_maps", "FigS2b_maps", "FigS3_calibration", "FigS4_domain_scans", "FigS5_radius",\n                "FigS6_reproduction", "FigS7_runtime", "FigS8_intervention"):',
    'for tag in ("FigS1_parameter_sensitivity", "FigS2a_maps", "FigS2b_maps", "FigS3_calibration",\n                "FigS4_domain_scans", "FigS5_radius", "FigS6_reproduction", "FigS7_runtime",\n                "FigS8_intervention"):'
)
BUILD_OR.write_text(build_or, encoding="utf-8")

build_is = BUILD_IS.read_text(encoding="utf-8")
build_is = build_is.replace(
    "supplementary methods, Tables S1–S9 and Figs. S1–S6.",
    "supplementary methods, Tables S1–S13 and Figs. S1–S8."
)
BUILD_IS.write_text(build_is, encoding="utf-8")
