"""Builds the Word reference document (styles) used by pandoc for the IS manuscript.

Springer guidance: plain 10-pt Times Roman, automatic page numbers, decimal headings
(max three levels), no field functions in the text.  We add 1.5 line spacing and
continuous line numbers for reviewers.
"""
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor, Cm

HERE = Path(__file__).resolve().parent
src = HERE / "build" / "ref_default.docx"
dst = HERE / "build" / "reference_IS.docx"

d = Document(str(src))
INK = RGBColor(0, 0, 0)


def font(style, size=10, bold=None, italic=None, name="Times New Roman"):
    f = style.font
    f.name = name
    f.size = Pt(size)
    f.color.rgb = INK
    if bold is not None:
        f.bold = bold
    if italic is not None:
        f.italic = italic
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = rpr.makeelement(qn("w:rFonts"), {})
        rpr.append(rfonts)
    for k in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rfonts.set(qn(k), name)
    for k in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rfonts.get(qn(k)) is not None:
            del rfonts.attrib[qn(k)]


def para(style, before=0, after=6, spacing=1.5, align=None, keep=False):
    pf = style.paragraph_format
    pf.space_before = Pt(before)
    pf.space_after = Pt(after)
    pf.line_spacing = spacing
    if align is not None:
        pf.alignment = align
    if keep:
        pf.keep_with_next = True


class _Styles:
    def __init__(self, styles):
        self._by = {st.name: st for st in styles}

    def __getitem__(self, k):
        return self._by[k]


S = _Styles(d.styles)
for nm in ("Normal", "Body Text", "First Paragraph", "Compact", "Block Text", "Bibliography", "Abstract"):
    font(S[nm], 10)
    para(S[nm], 0, 6, 1.5, WD_ALIGN_PARAGRAPH.JUSTIFY)
para(S["Compact"], 0, 2, 1.0, WD_ALIGN_PARAGRAPH.LEFT)
font(S["Title"], 14, bold=True)
para(S["Title"], 0, 12, 1.15, WD_ALIGN_PARAGRAPH.LEFT)
font(S["Heading 1"], 12, bold=True, italic=False)
para(S["Heading 1"], 14, 6, 1.15, WD_ALIGN_PARAGRAPH.LEFT, keep=True)
font(S["Heading 2"], 10, bold=True, italic=False)
para(S["Heading 2"], 10, 4, 1.15, WD_ALIGN_PARAGRAPH.LEFT, keep=True)
font(S["Heading 3"], 10, bold=False, italic=True)
para(S["Heading 3"], 8, 4, 1.15, WD_ALIGN_PARAGRAPH.LEFT, keep=True)
for nm in ("Caption", "Table Caption", "Image Caption"):
    font(S[nm], 9, italic=False)
    para(S[nm], 4, 10, 1.15, WD_ALIGN_PARAGRAPH.JUSTIFY)
para(S["Table Caption"], 10, 4, 1.15, WD_ALIGN_PARAGRAPH.JUSTIFY, keep=True)
font(S["Footnote Text"], 9)

# page setup: A4, 2.5 cm margins, continuous line numbers
sec = d.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
    setattr(sec, side, Cm(2.5))
d.save(str(dst))
print("wrote", dst)
