#!/usr/bin/env python
"""
Build the IS submission manuscript (anonymised) from docs_is/manuscript_IS.md.

    python docs_is/build_is.py            # facts -> filled markdown -> DOCX + LaTeX PDF

Every {{key}} placeholder is replaced from results/validation/is_manuscript_facts.json;
an unknown key aborts the build, so no number can be typed by hand.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs_is"
sys.path.insert(0, str(DOCS))
from ooxml_fix import fix_docx  # noqa: E402

BUILD = DOCS / "build"
FIG = ROOT / "results" / "figures" / "is"
BUILD.mkdir(parents=True, exist_ok=True)

CAPTIONS = {
    "Fig1_overview": (
        "**Fig. 1** SPARTA overview. **a** Workflow. **b** Cellular barrier on a toy hexagonal lattice: edge "
        "capacity (width and shade) falls with the matrix (ECM) and fibroblast (CAF) scores (Eq. 1); the minimum "
        "cut (black bars) separates immune-entry sources (blue) from tumour-core sinks (black squares) and runs "
        "through the gap in a matrix-rich band (grey nodes). **c** IgG-sized transport on the same lattice: steady "
        "diffusion from vessels ($\\varphi=1$) along edges whose conductance falls with matrix density and steric "
        "exclusion (Eq. 2), with absorption at ligand-expressing spots (rings; Eq. 3); fill shows "
        "$B_{\\mathrm{mAb}}=-\\log\\varphi$. **d**–**g** The primary-cohort section whose field association equals "
        "the cohort median (CSCC04; cSCC, Visium): compartments (**d**), ECM score (**e**), $B_{\\mathrm{cell}}$ "
        "field with minimum-cut edges (black bars; **f**) and $B_{\\mathrm{mAb}}$ field (**g**), shown as "
        "within-section ranks. Scale bars, 1 mm"),
    "Fig2_synthetic": (
        "**Fig. 2** Simulation benchmark with planted barriers. **a** Four of the nine geometries, with vessel "
        "sources (blue), tumour-core sinks (black squares), matrix-rich spots (dark grey), the tumour nest (light "
        "grey) and the minimum cut computed from noisy observed scores (black bars); access is the fraction of "
        "agents reaching the core. **b** Ground-truth access for every simulated tissue (20 replicates at each of "
        "three noise levels per geometry); bars, means. **c** Spearman correlation of each summary with lost "
        "access over all tissues (dots) and within each noise level (ticks); SPARTA summaries in blue. "
        "**d** AUC for four contrasts, rows as in **c**; 0.5 is chance. Composition is matched across "
        "geometries except for the distant capsule, which needs more spots to close"),
    "Fig3_calibration": (
        "**Fig. 3** Calibration of the surrogate test on the {{cal_n_graphs}} real section graphs. **a** False-"
        "positive rate at α = 0.05 for fields that share the expected graph power spectrum of the section's "
        "$B_{\\mathrm{mAb}}$ field but are independent of it by construction ({{cal_n_sim}} simulations per graph): "
        "graph-spectral surrogate test on the raw spectrum (filled circles) and on normal scores (open squares), and a "
        "point-level test that treats spots as independent (open triangles); the grey band is the 95% binomial "
        "range around 0.05. **b** Pooled false-positive rate by scenario for the "
        "raw-spectrum surrogate test, its normal-score variant and the point-level test. **c** Distribution of "
        "pooled surrogate-test p-values in the spectrum-matched scenario; the line marks the uniform expectation"),
    "Fig4_association": (
        "**Fig. 4** Field association within sections and patient-level cohorts. **a** Vessel-distance-adjusted "
        "partial Spearman correlations for all 77 transcriptomics sections, grouped by cohort and tumour stratum; "
        "marker shape encodes platform, filled markers indicate BH q < 0.05 within that section's analysis family, "
        "and grey bars show the stratum interquartile range with the median marked. The 11 GSE289745 sections with "
        "unavailable patient relationships are included as conservative section-level units. **b** Pooled estimates "
        "with 95% CIs from nested random-effects models for the primary, replication and extension cohorts and for "
        "the explicitly labelled combined analyses; text gives the number of positive patient or analytical units. "
        "**c** The primary-cohort association under alternative adjustments for vessel distance (*primary analysis)"
    ),
    "Fig5_structure": (
        "**Fig. 5** How much of the coupling does the model produce by itself? **a** Observed association "
        "(markers) against the 95% ranges of the geometry null (light grey; all $B_{\\mathrm{mAb}}$ inputs replaced "
        "by graph-spectral surrogates) and the construction null (dark grey; real ECM score retained, crosslinking "
        "and ligand inputs replaced); filled markers, observed association above the construction null at "
        "p < 0.05. **b** Excess of the observed association over the construction-null mean per patient (dots, "
        "means; ticks, sections) and pooled across patients (diamond with 95% CI, nested random-effects model). "
        "**c** Association before and after ablating the shared ECM term (grey lines, sections; black, medians). "
        "**d** Within-section Spearman correlations between model inputs (XL, crosslinking score; lig., ligand "
        "score)"),
    "Fig6_parameters": (
        "**Fig. 6** Parameter dependence and in-model experiments. **a** Effective pore radius as a function of "
        "the within-section crosslinking rank for three values of β; for an IgG-sized molecule, edges right of "
        "$X^{*}$ are completely excluded at the default β = 3 (shaded). **b** Fraction of completely excluded "
        "edges as a function of β in the 19 primary sections (grey) and their median (orange). **c** Mean "
        "tumour-core $B_{\\mathrm{mAb}}$ as the specified molecular radius increases with all inputs held fixed. "
        "**d** Selection-matched gap experiment: residual $B_{\\mathrm{cell}}$ after scattered removal divided by "
        "that after contiguous removal of 20% of the cut nodes, for each primary section (filled, BH q < 0.05) "
        "and for simulated closed capsules of three thicknesses (open circles, 10 replicates each; arrowheads, "
        "values beyond the axis, with their number); ratios above 1 favour the contiguous gap"),
    "Fig7_codex": (
        "**Fig. 7** The minimum cut against measured CD8$^+$ T-cell positions in CODEX images of colorectal cancer. "
        "**a** Two cores with a high (left) and a low (right) geometry-normalised barrier $B_{\\mathrm{rel}}$: "
        "scaffold cells (stroma shaded by collagen IV, other tumour cells light grey, tumour-core sinks dark squares, "
        "vascular sources blue), minimum-cut edges (black) and the withheld CD8$^+$ T cells (green crosses); IR, "
        "log$_2$ infiltration ratio of the tumour core. **b** Patient means of $\\log B_{\\mathrm{rel}}$ and IR (open "
        "circles, Crohn's-like reaction, CLR; filled squares, diffuse inflammatory infiltration, DII). **c** Spearman "
        "correlation of each summary with IR: dots, core level with patient-cluster bootstrap 95% CI (bars); ticks, "
        "patient level; SPARTA summaries in blue; the two geometry-only summaries were added post hoc. **d** "
        "Patient-mean $\\log B_{\\mathrm{rel}}$ by group; lines, medians"),
    "Fig8_baselines": (
        "**Fig. 8** Cross-cohort robustness and incremental evidence. **a** Median observed field association, "
        "construction-null association and association after shared-ECM ablation in the primary, replication and "
        "extension cohorts. **b** Extension sections grouped by disease stratum; black dots, medians. **c** Extension "
        "sections grouped by platform; vertical grey bars show interquartile ranges and black dots medians. "
        "**d** Median within-section correlations of the two fields with stromal density, tumour distance and "
        "neighbourhood enrichment, comparing the primary/replication and extension cohorts. **e** Primary-cohort "
        "association as the graph radius is perturbed relative to the platform default. **f** Pooled false-positive "
        "rates for graph-spectral, normal-score and point-level tests in the primary/replication and extension "
        "calibrations; dashed line marks α = 0.05. The figure supports robustness of the structural characterization; "
        "it is not evidence of treatment response, prognosis or clinical utility"
    ),
}

# figure -> insert before this exact line of the markdown
PLACE = {
    "Fig1_overview": "# 2 Materials and Methods",
    "Fig2_synthetic": "**Table 2** Simulation benchmark.",
    "Fig3_calibration": "## 3.3 The two model fields co-vary within sections and across patients",
    "Fig4_association": "**Table 3** Section-level and patient-level inference.",
    "Fig5_structure": "## 3.5 The size-exclusion term is parameter-driven",
    "Fig6_parameters": "## 3.7 The minimum cut tracks measured CD8^+^ T-cell exclusion in multiplexed images",
    "Fig7_codex": "## 3.8 Simple spatial summaries reproduce only part of either operator",
    "Fig8_baselines": "## 3.9 Modality-specific in-model perturbation maps",
}

DECLARATIONS = """
# Statements and declarations

Author names, affiliations, funding and acknowledgements are given on the separate title page to preserve
anonymity during double-blind review.

**Competing interests** None declared.

**Ethics approval and consent** Not applicable. This study is a secondary analysis of publicly available,
de-identified data; no new samples or participants were involved.

**Data availability** Primary data are available from the Gene Expression Omnibus under accessions GSE144239 and
GSE250636; external sections are 10x Genomics public datasets [27, 28]. The replication cohort is available from the
publishers of the original study [@thrane2018] (https://www.spatialresearch.org), and the CODEX single-cell table from
Mendeley Data [@schurch2020data] (images: The Cancer Imaging Archive, https://doi.org/10.7937/TCIA.2020.FQN0-0326).
Derived node tables, graph files and all per-section and per-core result files are included in the anonymised archive
provided to reviewers (Online Resource 3).

**Code availability** The complete source code, configuration files, analysis scripts and tests are provided to
reviewers as an anonymised archive (Online Resource 3). They will be released under the MIT licence in a public
repository with a versioned archive DOI upon acceptance (identifiers withheld for double-blind review).

**Online Resources** Online Resource 1 (PDF): supplementary methods, Tables S1–S13 and Figs. S1–S9. Online Resource 2
(XLSX): per-section results. Online Resource 3 (ZIP): anonymised code and result archive.
"""


def load_facts():
    p = ROOT / "results" / "validation" / "is_manuscript_facts.json"
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def fill(text, facts):
    missing = []

    def rep(m):
        k = m.group(1).strip()
        if k not in facts:
            missing.append(k)
            return "{{" + k + "}}"
        # typographic minus for negative numbers (also inside comma-separated lists)
        return re.sub(r"(?<![\w.])-(?=\d)", "\u2212", str(facts[k]))
    out = re.sub(r"\{\{\s*([A-Za-z0-9_]+)\s*\}\}", rep, text)
    if missing:
        import os
        if os.environ.get("IS_DRAFT"):
            print(f"[draft] unfilled placeholders: {sorted(set(missing))}")
        else:
            raise SystemExit(f"unfilled placeholders: {sorted(set(missing))}")
    return out


# keyed citations used by additions after v2.1.0 -> their entry number in references.md
REF_KEYS = {"thrane2018": 43, "schurch2020": 44, "schurch2020data": 45}
_CITE = re.compile(r"\[(\d+(?:\s*[,\u2013-]\s*\d+)*)\]")


def _expand(tok):
    out = []
    for part in tok.split(","):
        part = part.strip()
        m = re.match(r"^(\d+)\s*[\u2013-]\s*(\d+)$", part)
        if m:
            out += list(range(int(m.group(1)), int(m.group(2)) + 1))
        elif part:
            out.append(int(part))
    return out


def _compress(nums):
    nums = sorted(set(nums))
    runs, start = [], None
    for i, n in enumerate(nums):
        if start is None:
            start = prev = n
            continue
        if n == prev + 1:
            prev = n
            continue
        runs.append((start, prev))
        start = prev = n
    if start is not None:
        runs.append((start, prev))
    txt = []
    for a, b in runs:
        if b - a >= 2:
            txt.append(f"{a}\u2013{b}")
        elif b == a + 1:
            txt += [str(a), str(b)]
        else:
            txt.append(str(a))
    return "[" + ", ".join(txt) + "]"


def resolve_and_renumber(md, refs_text):
    """[@key] -> numbers, then renumber every citation by order of first appearance (math is skipped)
    and reorder the reference list accordingly. Aborts on unknown keys or uncited references."""
    def key_rep(m):
        keys = [k.strip().lstrip("@") for k in m.group(1).split(";")]
        missing = [k for k in keys if k not in REF_KEYS]
        if missing:
            raise SystemExit(f"unknown citation keys {missing}")
        return "[" + ", ".join(str(REF_KEYS[k]) for k in keys) + "]"
    md = re.sub(r"\[(@[A-Za-z0-9_]+(?:\s*;\s*@[A-Za-z0-9_]+)*)\]", key_rep, md)
    parts = re.split(r"(\$[^$]*\$)", md)
    order = []
    for i, part in enumerate(parts):
        if i % 2:
            continue
        for m in _CITE.finditer(part):
            for n in _expand(m.group(1)):
                if n not in order:
                    order.append(n)
    new_of = {old: k + 1 for k, old in enumerate(order)}
    for i, part in enumerate(parts):
        if i % 2:
            continue
        parts[i] = _CITE.sub(lambda m: _compress([new_of[n] for n in _expand(m.group(1))]), part)
    entries = {}
    for line in refs_text.strip().splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m:
            entries[int(m.group(1))] = m.group(2)
    uncited = sorted(set(entries) - set(order))
    if uncited:
        raise SystemExit(f"references never cited: {uncited}")
    missing = [n for n in order if n not in entries]
    if missing:
        raise SystemExit(f"citations without a reference entry: {missing}")
    refs_new = "\n".join(f"{k + 1}. {entries[old]}" for k, old in enumerate(order))
    return "".join(parts), refs_new


def assemble(facts, target):
    md = (DOCS / "manuscript_IS.md").read_text(encoding="utf-8")
    md = fill(md, facts)
    for name, anchor in PLACE.items():
        ext = "png" if target == "docx" else "pdf"
        path = FIG / f"{name}.{ext}"
        if not path.exists():
            print(f"[warn] missing figure {path.name}; skipped")
            continue
        cap = fill(CAPTIONS[name], facts)
        width = "16cm" if target == "docx" else "100%"
        block = f"![{cap}]({path.as_posix()}){{width={width}}}\n\n"
        if anchor not in md:
            raise SystemExit(f"anchor not found for {name}: {anchor}")
        md = md.replace(anchor, block + anchor, 1)
    refs = (DOCS / "references.md").read_text(encoding="utf-8")
    md = md.rstrip() + "\n\n" + DECLARATIONS.strip()
    md, refs = resolve_and_renumber(md, refs)
    md = md + "\n\n# References\n\n" + refs + "\n"
    return md


def build_docx(facts):
    md = assemble(facts, "docx")
    src = BUILD / "manuscript_IS_docx.md"
    src.write_text(md, encoding="utf-8")
    out = BUILD / "Manuscript_IS_anonymised.docx"
    subprocess.run(["pandoc", str(src), "-f", "markdown+tex_math_dollars+pipe_tables+implicit_figures",
                    "-t", "docx", "--reference-doc", str(BUILD / "reference_IS.docx"), "-o", str(out)],
                   check=True)
    postprocess_docx(out)
    fix_docx(out)
    return out


def postprocess_docx(path):
    """Line numbers, page numbers, table borders; keep everything as plain Word content."""
    from docx import Document
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Pt

    d = Document(str(path))
    sec = d.sections[0]
    sp = sec._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "283")
    sp.append(ln)
    # page number footer (PAGE field is the automatic numbering Springer asks for)
    footer = sec.footer
    p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
    p.alignment = 1
    p._p.get_or_add_pPr().append(OxmlElement("w:suppressLineNumbers"))   # footer is not a numbered line
    r = p.add_run()
    for tag, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            fc = OxmlElement("w:fldChar")
            fc.set(qn("w:fldCharType"), tag)
            r._r.append(fc)
        else:
            it = OxmlElement("w:instrText")
            it.set(qn("xml:space"), "preserve")
            it.text = txt
            r._r.append(it)
    r.font.size = Pt(9)
    # table captions ("Table N ...") stay on the same page as their table
    for para in d.paragraphs:
        if re.match(r"Table \d+\b", para.text or ""):
            para.paragraph_format.keep_with_next = True
    # simple three-line tables: top/bottom rule + header rule, 9-pt text
    for t in d.tables:
        tblPr = t._tbl.tblPr
        borders = OxmlElement("w:tblBorders")
        for edge, val in (("top", "single"), ("bottom", "single"), ("left", "nil"), ("right", "nil"),
                          ("insideH", "nil"), ("insideV", "nil")):
            el = OxmlElement(f"w:{edge}")
            el.set(qn("w:val"), val)
            if val == "single":
                el.set(qn("w:sz"), "6")
                el.set(qn("w:color"), "000000")
            borders.append(el)
        tblPr.append(borders)
        if t.rows:
            for cell in t.rows[0].cells:
                tcPr = cell._tc.get_or_add_tcPr()
                tb = OxmlElement("w:tcBorders")
                b = OxmlElement("w:bottom")
                b.set(qn("w:val"), "single")
                b.set(qn("w:sz"), "4")
                b.set(qn("w:color"), "000000")
                tb.append(b)
                tcPr.append(tb)
        for row in t.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    para.paragraph_format.line_spacing = 1.0
                    para.paragraph_format.space_after = Pt(0)
                    for run in para.runs:
                        run.font.size = Pt(8.5)
    d.save(str(path))


LATEX_PREAMBLE = r"""
\usepackage[a4paper,margin=2.5cm]{geometry}
\DeclareUnicodeCharacter{03B1}{\ensuremath{\alpha}}
\DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
\DeclareUnicodeCharacter{03BB}{\ensuremath{\lambda}}
\DeclareUnicodeCharacter{03BC}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{03C1}{\ensuremath{\rho}}
\DeclareUnicodeCharacter{03C4}{\ensuremath{\tau}}
\DeclareUnicodeCharacter{03BE}{\ensuremath{\xi}}
\DeclareUnicodeCharacter{03C6}{\ensuremath{\varphi}}
\DeclareUnicodeCharacter{03BA}{\ensuremath{\kappa}}
\DeclareUnicodeCharacter{03A6}{\ensuremath{\Phi}}
\DeclareUnicodeCharacter{03C3}{\ensuremath{\sigma}}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{00B7}{\ensuremath{\cdot}}
\DeclareUnicodeCharacter{00D7}{\ensuremath{\times}}
\usepackage{mathptmx}
\usepackage[scaled=0.92]{helvet}
\usepackage{setspace}\onehalfspacing
\usepackage[running]{lineno}\linenumbers
\renewcommand\linenumberfont{\normalfont\tiny\sffamily\color{gray}}
\usepackage{xcolor}
\usepackage{booktabs,longtable,array}
\setlength{\LTcapwidth}{\textwidth}
\usepackage{needspace}
\usepackage{etoolbox}
\usepackage{caption}
\captionsetup{labelformat=empty,font=small,justification=justified,singlelinecheck=false}
\usepackage{float}
\makeatletter\def\fps@figure{!tbp}\makeatother
\usepackage{xurl}
\usepackage[hidelinks]{hyperref}
\usepackage{titlesec}
\titleformat{\section}{\normalfont\large\bfseries}{}{0pt}{}
\titleformat{\subsection}{\normalfont\normalsize\bfseries}{}{0pt}{}
\titlespacing*{\section}{0pt}{12pt}{4pt}
\titlespacing*{\subsection}{0pt}{8pt}{3pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\usepackage{fancyhdr}\pagestyle{fancy}\fancyhf{}\cfoot{\small\thepage}\renewcommand{\headrulewidth}{0pt}
\sloppy
"""


def keep_tables_together(t):
    """Move each "Table N" / "Table SN" caption paragraph into its longtable as an unnumbered caption,
    so the caption is part of the table's first head and can never be separated from it."""
    pat = re.compile(r"\n(\\textbf\{Table S?\d+\}(?:[^\n]+\n)+)\n(\\begin\{longtable\}\[\]\{@\{\}.*?@\{\}\}\n)",
                     re.S)

    def rep(m):
        cap = " ".join(m.group(1).split())
        return "\n" + m.group(2) + "\\caption*{" + cap + "}\\tabularnewline\n"
    t = pat.sub(rep, t)
    # the caption belongs to the first head only; continuation pages repeat just the column header
    head = re.compile(r"(\\caption\*\{[^\n]*\}\\tabularnewline\n)(.*?)\\endhead\n", re.S)
    return head.sub(lambda m: m.group(1) + m.group(2) + "\\endfirsthead\n" + m.group(2) + "\\endhead\n", t)


def build_pdf(facts):
    md = assemble(facts, "pdf")
    src = BUILD / "manuscript_IS_tex.md"
    src.write_text(md, encoding="utf-8")
    hdr = BUILD / "preamble.tex"
    hdr.write_text(LATEX_PREAMBLE, encoding="utf-8")
    tex = BUILD / "Manuscript_IS_anonymised.tex"
    subprocess.run(["pandoc", str(src), "-f", "markdown+tex_math_dollars+pipe_tables+implicit_figures",
                    "-t", "latex", "-s", "-H", str(hdr), "-V", "documentclass=article", "-V", "fontsize=10pt",
                    "-o", str(tex)], check=True)
    t = tex.read_text(encoding="utf-8")
    # pandoc numbers sections when asked; we carry explicit decimal numbers in the headings
    t = t.replace("\\usepackage{lmodern}", "")
    t = keep_tables_together(t)
    tex.write_text(t, encoding="utf-8")
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name], cwd=BUILD,
                           capture_output=True, text=True, errors="replace")
    if r.returncode != 0:
        print(r.stdout[-3000:])
        raise SystemExit("pdflatex failed")
    return BUILD / "Manuscript_IS_anonymised.pdf"


def word_counts(facts):
    md = fill((DOCS / "manuscript_IS.md").read_text(encoding="utf-8"), facts)
    abstract = md.split("# Abstract", 1)[1].split("**Keywords**", 1)[0]
    body = md.split("# 1 Introduction", 1)[1]
    keep = []
    for line in body.splitlines():
        if line.startswith("|") or line.startswith("**Table") or line.startswith("!["):
            continue
        keep.append(line)
    body = "\n".join(keep)
    count = lambda t: len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'’.\-–/]*", re.sub(r"\$[^$]*\$", "x", t)))  # noqa: E731
    return count(abstract), count(body)


def build_side_docs(facts):
    import datetime
    a, m = word_counts(facts)
    f2 = dict(facts, abstract_words=a, main_words=f"{round(m, -1):,}",
              date=datetime.date.today().strftime("%d %B %Y").lstrip("0"))
    outs = []
    for stem, title in (("title_page", "Title_page_IS"), ("cover_letter", "Cover_letter_IS")):
        md = fill((DOCS / f"{stem}.md").read_text(encoding="utf-8"), f2)
        src = BUILD / f"{stem}_filled.md"
        src.write_text(md, encoding="utf-8")
        out = BUILD / f"{title}.docx"
        subprocess.run(["pandoc", str(src), "-f", "markdown+hard_line_breaks" if stem == "cover_letter" else "markdown",
                        "-t", "docx", "--reference-doc", str(BUILD / "reference_IS.docx"), "-o", str(out)], check=True)
        fix_docx(out)
        hdr = BUILD / "preamble_side.tex"
        hdr.write_text(LATEX_PREAMBLE.replace("\\usepackage[running]{lineno}\\linenumbers", "")
                       .replace("\\renewcommand\\linenumberfont{\\normalfont\\tiny\\sffamily\\color{gray}}", "")
                       .replace("\\usepackage{setspace}\\onehalfspacing", "\\usepackage{setspace}\\setstretch{1.15}"),
                       encoding="utf-8")
        tex = BUILD / f"{title}.tex"
        subprocess.run(["pandoc", str(src), "-f", "markdown+hard_line_breaks" if stem == "cover_letter" else "markdown",
                        "-t", "latex", "-s", "-H", str(hdr), "-V", "fontsize=11pt", "-o", str(tex)], check=True)
        t = tex.read_text(encoding="utf-8").replace("\\usepackage{lmodern}", "")
        if stem == "cover_letter":            # one page: slightly tighter margins and spacing
            t = t.replace("\\maketitle", "")
            t = t.replace("\\begin{document}", "\\newgeometry{margin=2cm}\n\\setstretch{1.05}\n\\begin{document}", 1)
        else:                                  # keep the title page on one page
            t = t.replace("\\begin{document}", "\\usepackage{titling}\n\\setlength{\\droptitle}{-6em}\n"
                          "\\posttitle{\\par\\end{center}\\vspace{-2.5em}}\n\\begin{document}", 1)
        tex.write_text(t, encoding="utf-8")
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name], cwd=BUILD,
                               capture_output=True, text=True, errors="replace")
        if r.returncode != 0:
            print(r.stdout[-2000:])
            raise SystemExit(f"pdflatex failed ({stem})")
        outs += [out, BUILD / f"{title}.pdf"]
    print(f"abstract {a} words; main text {m} words")
    return outs


if __name__ == "__main__":
    facts = load_facts()
    which = sys.argv[1:] or ["docx", "pdf", "side"]
    if "side" in which:
        print("side:", build_side_docs(facts))
    if "docx" in which:
        print("DOCX:", build_docx(facts))
    if "pdf" in which:
        print("PDF:", build_pdf(facts))
