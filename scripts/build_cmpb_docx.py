#!/usr/bin/env python
"""Build the editable CMPB submission manuscript from its Markdown source."""
from __future__ import annotations

import re
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "manuscript_cmpb_draft.md"
OUT = ROOT / "docs" / "manuscript_cmpb_draft.docx"
FIGURES = ROOT / "results" / "figures" / "cmpb"
FIGURE_NAMES = {
    "figure_1_framework": "Figure 1. Framework schematic; illustrative graph only, not patient tissue.",
    "figure_2_coupling": "Figure 2. Section-level adjusted field association and BH q-values.",
    "figure_3_parameter_sensitivity": "Figure 3. Parameterized size-exclusion law and edge-exclusion sensitivity.",
    "figure_4_s2_matched": "Figure 4. Selection-matched in-model counterfactual ratios.",
    "figure_5_size_scan": "Figure 5. Fixed-input molecular-radius sensitivity scan.",
}


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)


def add_field(paragraph, name: str) -> None:
    run = paragraph.add_run()
    fld_char1 = OxmlElement("w:fldChar")
    fld_char1.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = name
    fld_char2 = OxmlElement("w:fldChar")
    fld_char2.set(qn("w:fldCharType"), "end")
    run._r.append(fld_char1)
    run._r.append(instr)
    run._r.append(fld_char2)


def format_run(run, *, bold=False, italic=False, code=False, size=None, color=None):
    run.bold = bold
    run.italic = italic
    if code:
        run.font.name = "Consolas"
    if size:
        run.font.size = Pt(size)
    if color:
        run.font.color.rgb = RGBColor.from_string(color)


def add_inline(paragraph, text: str, size=None):
    # Small Markdown subset: bold, italic, and inline code.
    pattern = re.compile(r"(\*\*[^*]+\*\*|(?<!\*)\*[^*]+\*(?!\*)|`[^`]+`)")
    pos = 0
    for match in pattern.finditer(text):
        if match.start() > pos:
            run = paragraph.add_run(text[pos:match.start()])
            format_run(run, size=size)
        token = match.group(0)
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2]); format_run(run, bold=True, size=size)
        elif token.startswith("*"):
            run = paragraph.add_run(token[1:-1]); format_run(run, italic=True, size=size)
        else:
            run = paragraph.add_run(token[1:-1]); format_run(run, code=True, size=size)
        pos = match.end()
    if pos < len(text):
        run = paragraph.add_run(text[pos:]); format_run(run, size=size)


def paragraph_defaults(p, *, alignment=None, spacing=2.0, after=Pt(7), before=Pt(0)):
    p.paragraph_format.line_spacing = spacing
    p.paragraph_format.space_after = after
    p.paragraph_format.space_before = before
    if alignment is not None:
        p.alignment = alignment


def add_table(doc, lines):
    rows = []
    for line in lines:
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and all(re.fullmatch(r":?-{3,}:?", c or "-") for c in cells):
            continue
        rows.append(cells)
    if not rows:
        return
    ncols = max(map(len, rows))
    table = doc.add_table(rows=0, cols=ncols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    table.autofit = False
    widths = [1.15, 2.35, .65, .95, 1.4] if ncols == 5 else [6.5 / ncols] * ncols
    for ci, width in enumerate(widths):
        table.columns[ci].width = Inches(width)
    for ri, values in enumerate(rows):
        row = table.add_row()
        if ri == 0:
            set_repeat_table_header(row)
        for ci in range(ncols):
            cell = row.cells[ci]
            cell.width = Inches(widths[ci])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.line_spacing = 1.05
            p.paragraph_format.space_after = Pt(2)
            add_inline(p, values[ci] if ci < len(values) else "", size=8.5)
            if ri == 0:
                set_cell_shading(cell, "DCE6EF")
                for run in p.runs:
                    run.bold = True
    doc.add_paragraph().paragraph_format.space_after = Pt(2)


def configure(doc):
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(1)
    sec.bottom_margin = Inches(1)
    sec.left_margin = Inches(1)
    sec.right_margin = Inches(1)
    sec.header_distance = Inches(.35)
    sec.footer_distance = Inches(.35)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(11)
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    normal.paragraph_format.line_spacing = 2.0
    normal.paragraph_format.space_after = Pt(7)

    for style_name, size in (("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 11)):
        style = doc.styles[style_name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        style.paragraph_format.keep_with_next = True
        style.paragraph_format.line_spacing = 1.15
        style.paragraph_format.space_before = Pt(12 if style_name == "Heading 1" else 8)
        style.paragraph_format.space_after = Pt(5)

    title_style = doc.styles["Title"]
    title_style.font.name = "Times New Roman"
    title_style.font.size = Pt(18)
    title_style.font.bold = True
    title_style.font.color.rgb = RGBColor(0, 0, 0)
    title_style._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    title_style.paragraph_format.keep_with_next = True
    title_style.paragraph_format.line_spacing = 1.15
    title_style.paragraph_format.space_after = Pt(18)
    # Word's built-in Title style carries a blue bottom rule in its default
    # template. Remove it explicitly for a plain submission manuscript title.
    title_ppr = title_style._element.get_or_add_pPr()
    for border in title_ppr.findall(qn("w:pBdr")):
        title_ppr.remove(border)

    # Continuous line numbers for editorial review.
    sect_pr = sec._sectPr
    line_num = OxmlElement("w:lnNumType")
    line_num.set(qn("w:countBy"), "1")
    line_num.set(qn("w:restart"), "continuous")
    line_num.set(qn("w:distance"), "360")
    sect_pr.append(line_num)

    footer = sec.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    footer.paragraph_format.line_spacing = 1
    add_field(footer, "PAGE")


def add_image(doc, name):
    img = FIGURES / f"{name}.png"
    if not img.exists():
        raise FileNotFoundError(img)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.line_spacing = 1.0
    p.paragraph_format.space_after = Pt(7)
    p.paragraph_format.keep_with_next = True
    shape = p.add_run().add_picture(str(img), width=Inches(6.38))
    shape._inline.docPr.set("descr", FIGURE_NAMES.get(name, name))


def main():
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    doc = Document()
    configure(doc)
    i = 0
    first_figure = True
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        if line.startswith("[[FIGURE:") and line.endswith("]]" ):
            name = line[len("[[FIGURE:"):-2]
            if first_figure:
                first_figure = False
            else:
                doc.add_page_break()
            add_image(doc, name)
            i += 1
            continue

        if line.startswith("|"):
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i].strip()); i += 1
            add_table(doc, block)
            continue

        if line.startswith("# "):
            p = doc.add_paragraph(style="Title")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            paragraph_defaults(p, alignment=WD_ALIGN_PARAGRAPH.CENTER, spacing=1.15,
                               after=Pt(18), before=Pt(18))
            add_inline(p, line[2:], size=18)
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if line.startswith("## "):
            title = line[3:]
            p = doc.add_paragraph(style="Heading 1")
            if title == "Figure legends":
                p.paragraph_format.page_break_before = True
            add_inline(p, title)
            i += 1
            continue

        if line.startswith("### "):
            p = doc.add_paragraph(style="Heading 2")
            add_inline(p, line[4:])
            i += 1
            continue

        if line.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            p.paragraph_format.line_spacing = 1.15
            p.paragraph_format.space_after = Pt(2)
            add_inline(p, line[2:])
            i += 1
            continue

        # Figure captions follow their image and should remain together.
        is_figure_caption = line.startswith("**Figure ")
        p = doc.add_paragraph()
        if line == "**Yize Li**" or line.startswith("Hangzhou Medical College,") or line.startswith("Corresponding author:"):
            paragraph_defaults(p, alignment=WD_ALIGN_PARAGRAPH.CENTER, spacing=1.1,
                               after=Pt(2), before=Pt(0))
            add_inline(p, line, size=10.5)
        elif is_figure_caption:
            paragraph_defaults(p, spacing=1.1, after=Pt(14))
            p.paragraph_format.keep_together = True
            add_inline(p, line, size=9)
        elif line.startswith("[1]") or re.match(r"^\[\d+\]", line):
            paragraph_defaults(p, spacing=1.15, after=Pt(6))
            p.paragraph_format.left_indent = Inches(.28)
            p.paragraph_format.first_line_indent = Inches(-.28)
            add_inline(p, line, size=9.5)
        else:
            paragraph_defaults(p)
            add_inline(p, line)
        i += 1

    # Avoid orphaning the author block: modestly tighten the title header.
    doc.core_properties.title = "SPARTA: dual graph operators for modelling immune-cell migration and IgG-sized transport barriers from spatial transcriptomics"
    doc.core_properties.subject = "Submission-style manuscript draft for Computer Methods and Programs in Biomedicine"
    doc.core_properties.author = "Yize Li"
    doc.core_properties.keywords = "spatial transcriptomics; graph theory; minimum cut; diffusion–absorption"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUT)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
