"""
Schema-order repair for pandoc/python-docx output (ECMA-376 transitional).

Word opens pandoc's .docx files, but strict validators flag a handful of ordering quirks that
come from pandoc's default reference styles: child order inside w:settings, w:style, w:tcPr and
m:dPr, two-digit list nsid values and w:pgMar without header/footer/gutter. Editorial systems
sometimes run such validators, so the submission files are repaired in place.

    python docs_is/ooxml_fix.py file.docx [...]
"""
from __future__ import annotations

import re
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"


def q(ns, tag):
    return f"{{{ns}}}{tag}"


SETTINGS_ORDER = """writeProtection view zoom removePersonalInformation removeDateAndTime doNotDisplayPageBoundaries
displayBackgroundShape printPostScriptOverText printFractionalCharacterWidth printFormsData embedTrueTypeFonts
embedSystemFonts saveSubsetFonts saveFormsData mirrorMargins alignBordersAndEdges bordersDoNotSurroundHeader
bordersDoNotSurroundFooter gutterAtTop hideSpellingErrors hideGrammaticalErrors activeWritingStyle proofState
formsDesign attachedTemplate linkStyles stylePaneFormatFilter stylePaneSortMethod documentType mailMerge
revisionView trackRevisions doNotTrackMoves doNotTrackFormatting documentProtection autoFormatOverride
styleLockTheme styleLockQFSet defaultTabStop autoHyphenation consecutiveHyphenLimit hyphenationZone
doNotHyphenateCaps showEnvelope summaryLength clickAndTypeStyle defaultTableStyle evenAndOddHeaders
bookFoldRevPrinting bookFoldPrinting bookFoldPrintingSheets drawingGridHorizontalSpacing drawingGridVerticalSpacing
displayHorizontalDrawingGridEvery displayVerticalDrawingGridEvery doNotUseMarginsForDrawingGridOrigin
drawingGridHorizontalOrigin drawingGridVerticalOrigin doNotShadeFormData noPunctuationKerning
characterSpacingControl printTwoOnOne strictFirstAndLastChars noLineBreaksAfter noLineBreaksBefore
savePreviewPicture doNotValidateAgainstSchema saveInvalidXml ignoreMixedContent alwaysShowPlaceholderText
doNotDemarcateInvalidXml saveXmlDataOnly useXSLTWhenSaving saveThroughXslt showXMLTags alwaysMergeEmptyNamespace
updateFields hdrShapeDefaults footnotePr endnotePr compat docVars rsids mathPr attachedSchema themeFontLang
clrSchemeMapping doNotIncludeSubdocsInStats doNotAutoCompressPictures forceUpgrade captions readModeInkLockDown
smartTagType schemaLibrary shapeDefaults doNotEmbedSmartTags decimalSymbol listSeparator""".split()

STYLE_ORDER = """name aliases basedOn next link autoRedefine hidden uiPriority semiHidden unhideWhenUsed qFormat locked
personal personalCompose personalReply rsid pPr rPr tblPr trPr tcPr tblStylePr""".split()

TCPR_ORDER = """cnfStyle tcW gridSpan hMerge vMerge tcBorders shd noWrap tcMar textDirection tcFitText vAlign hideMark
headers cellIns cellDel cellMerge tcPrChange""".split()

DPR_ORDER = "begChr sepChr endChr grow shp ctrlPr".split()

PPR_ORDER = """pStyle keepNext keepLines pageBreakBefore framePr widowControl numPr suppressLineNumbers pBdr shd tabs
suppressAutoHyphens kinsoku wordWrap overflowPunct topLinePunct autoSpaceDE autoSpaceDN bidi adjustRightInd
snapToGrid spacing ind contextualSpacing mirrorIndents suppressOverlap jc textDirection textAlignment
textboxTightWrap outlineLvl divId cnfStyle rPr sectPr pPrChange""".split()

RPR_ORDER = """rStyle rFonts b bCs i iCs caps smallCaps strike dstrike outline shadow emboss imprint noProof snapToGrid
vanish webHidden color spacing w kern position sz szCs highlight u effect bdr shd fitText vertAlign rtl cs em lang
eastAsianLayout specVanish oMath""".split()

TBLPR_ORDER = """tblStyle tblpPr tblOverlap bidiVisual tblStyleRowBandSize tblStyleColBandSize tblW jc tblCellSpacing
tblInd tblBorders shd tblLayout tblCellMar tblLook tblCaption tblDescription tblPrChange""".split()

SECTPR_ORDER = """headerReference footerReference footnotePr endnotePr type pgSz pgMar paperSrc pgBorders lnNumType
pgNumType cols formProt vAlign noEndnote titlePg textDirection bidi rtlGutter docGrid printerSettings
sectPrChange""".split()

MRPR_ORDER = "lit nor scr sty brk aln".split()

BORDER_ORDER = "top start left bottom end right between bar insideH insideV tl2br tr2bl".split()


def fix_common(root):
    """Property-element child order shared by styles.xml and document.xml."""
    for el in root.iter(q(W, "pPr")):
        reorder(el, PPR_ORDER)
    for el in root.iter(q(W, "rPr")):
        reorder(el, RPR_ORDER)
    for el in root.iter(q(W, "tblPr")):
        reorder(el, TBLPR_ORDER)
    for el in root.iter(q(W, "tcPr")):
        reorder(el, TCPR_ORDER)
    for el in root.iter(q(W, "sectPr")):
        reorder(el, SECTPR_ORDER)
    for el in root.iter(q(M, "rPr")):
        # m:nor (normal text) and m:sty are alternatives in the schema; keep m:nor
        if el.find(q(M, "nor")) is not None:
            for sty in el.findall(q(M, "sty")):
                el.remove(sty)
        reorder(el, MRPR_ORDER)
    for tag in ("tblBorders", "tcBorders", "pBdr"):
        for el in root.iter(q(W, tag)):
            reorder(el, BORDER_ORDER)
    for el in root.iter(q(M, "dPr")):
        reorder(el, DPR_ORDER)


def reorder(parent, order):
    rank = {name: i for i, name in enumerate(order)}
    kids = list(parent)
    if not kids:
        return

    def key(item):
        idx, el = item
        if not isinstance(el.tag, str):          # comments / processing instructions stay in place order
            return (len(order) + 1, idx)
        local = etree.QName(el).localname
        return (rank.get(local, len(order)), idx)
    ordered = [el for _, el in sorted(enumerate(kids), key=key)]
    if ordered != kids:
        for el in kids:
            parent.remove(el)
        for el in ordered:
            parent.append(el)


def fix_settings(root):
    reorder(root, SETTINGS_ORDER)


def fix_styles(root):
    for st in root.iter(q(W, "style")):
        reorder(st, STYLE_ORDER)
    fix_common(root)
    for rpr in root.iter(q(W, "rPr")):                # stray character content inside property elements
        rpr.text = None
        for ch in rpr:
            ch.tail = None


def fix_numbering(root):
    for el in root.iter(q(W, "nsid")):
        v = el.get(q(W, "val"))
        if v and len(v) < 8:
            el.set(q(W, "val"), v.upper().zfill(8))


def fix_document(root):
    for pg in root.iter(q(W, "pgMar")):
        for att, default in (("header", "709"), ("footer", "709"), ("gutter", "0")):
            if pg.get(q(W, att)) is None:
                pg.set(q(W, att), default)
    fix_common(root)


FIXERS = {"word/settings.xml": fix_settings, "word/styles.xml": fix_styles,
          "word/numbering.xml": fix_numbering, "word/document.xml": fix_document}


def fix_docx(path):
    path = Path(path)
    tmp = Path(tempfile.mkdtemp())
    out = tmp / path.name
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            fx = FIXERS.get(info.filename)
            if fx is None and re.fullmatch(r"word/(header|footer)\d*\.xml", info.filename):
                fx = fix_common
            if fx is not None:
                root = etree.fromstring(data)
                fx(root)
                data = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
            zout.writestr(info, data)
    shutil.move(str(out), str(path))
    shutil.rmtree(tmp, ignore_errors=True)
    return path


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print("fixed", fix_docx(p))
