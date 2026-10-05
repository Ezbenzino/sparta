#!/usr/bin/env python
"""
Assemble the Editorial Manager upload folder ``submission_IS/`` from the built documents and figures.

    python docs_is/build_is.py            # manuscript, title page, cover letter
    python docs_is/build_or1.py           # Online Resources 1 and 2
    python docs_is/build_or3.py           # Online Resource 3 (anonymised archive)
    python docs_is/package_submission.py  # this script

Figures: combination figures (raster maps or rasterised scatter) go in as 600-dpi TIFF; line art goes in
as EPS with 1200-dpi TIFF fall-backs; vector PDFs of every figure are kept for production. Existing files
in ``submission_IS/`` are replaced; files that no longer belong to the submission are moved to
``submission_IS/_superseded/`` rather than deleted.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "docs_is" / "build"
FIG = ROOT / "results" / "figures" / "is"
OUT = ROOT / "submission_IS"

COMBINATION = ["Fig1_overview", "Fig7_codex", "Fig8_baselines"]
LINE_ART = ["Fig2_synthetic", "Fig3_calibration", "Fig4_association", "Fig5_structure", "Fig6_parameters"]
DOCS = {
    "Manuscript_IS_anonymised.docx": BUILD / "Manuscript_IS_anonymised.docx",
    "Manuscript_IS_anonymised.pdf": BUILD / "Manuscript_IS_anonymised.pdf",
    "Title_page_IS.docx": BUILD / "Title_page_IS.docx",
    "Title_page_IS.pdf": BUILD / "Title_page_IS.pdf",
    "Cover_letter_IS.docx": BUILD / "Cover_letter_IS.docx",
    "Cover_letter_IS.pdf": BUILD / "Cover_letter_IS.pdf",
    "ESM_1.pdf": BUILD / "Online_Resource_1.pdf",
    "ESM_2.xlsx": BUILD / "Online_Resource_2_per_section_results.xlsx",
    "ESM_3.zip": BUILD / "Online_Resource_3_code_and_results.zip",
    "suggested_reviewers.md": ROOT / "docs_is" / "suggested_reviewers.md",
    "SUBMISSION_CHECKLIST_zh.md": ROOT / "docs_is" / "SUBMISSION_CHECKLIST_zh.md",
}


def line_art_tiff(pdf: Path, out: Path, dpi: int = 1200) -> None:
    """Render a vector PDF to an LZW-compressed RGB TIFF (Springer line-art fall-back)."""
    stem = out.with_suffix("")
    subprocess.run(["pdftoppm", "-r", str(dpi), "-tiff", "-tiffcompression", "lzw", "-singlefile", str(pdf),
                    str(stem)], check=True)
    from PIL import Image
    im = Image.open(out)
    if im.mode != "RGB":
        im = im.convert("RGB")
    im.save(out, compression="tiff_lzw", dpi=(dpi, dpi))


def main() -> int:
    # build_or3.py writes its archive into submission_IS/ by default; keep it in the build folder instead
    for name in ("Online_Resource_3_code_and_results.zip", "Online_Resource_3_code_and_results.prev.zip"):
        if (OUT / name).exists():
            shutil.move(str(OUT / name), str(BUILD / name))
    missing = [k for k, p in DOCS.items() if not p.exists()]
    missing += [f"{n}.pdf" for n in COMBINATION + LINE_ART if not (FIG / f"{n}.pdf").exists()]
    if missing:
        print("missing inputs:", missing)
        return 1
    OUT.mkdir(exist_ok=True)
    keep = set(DOCS) | {"Figures", "_superseded"}
    sup = OUT / "_superseded"
    for p in OUT.iterdir():
        if p.name not in keep:
            sup.mkdir(exist_ok=True)
            shutil.move(str(p), str(sup / p.name))
    for name, src in DOCS.items():
        if src.resolve() != (OUT / name).resolve():
            shutil.copy2(src, OUT / name)
    figs = OUT / "Figures"
    for sub in ("eps", "pdf", "line_art_1200dpi"):
        (figs / sub).mkdir(parents=True, exist_ok=True)
    expected = ({f"{n}.tif" for n in COMBINATION} | {f"eps/{n}.eps" for n in LINE_ART}
                | {f"pdf/{n}.pdf" for n in COMBINATION + LINE_ART} | {f"line_art_1200dpi/{n}.tif" for n in LINE_ART})
    for p in list(figs.rglob("*")):
        if p.is_file() and p.relative_to(figs).as_posix() not in expected:
            sup.mkdir(exist_ok=True)
            shutil.move(str(p), str(sup / f"Figures_{p.relative_to(figs).as_posix().replace('/', '_')}"))
    for n in COMBINATION:
        shutil.copy2(FIG / f"{n}.tif", figs / f"{n}.tif")
    for n in COMBINATION + LINE_ART:
        shutil.copy2(FIG / f"{n}.pdf", figs / "pdf" / f"{n}.pdf")
    for n in LINE_ART:
        shutil.copy2(FIG / f"{n}.eps", figs / "eps" / f"{n}.eps")
        line_art_tiff(FIG / f"{n}.pdf", figs / "line_art_1200dpi" / f"{n}.tif")
    print(f"submission folder ready: {OUT} ({date.today().isoformat()})")
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and "_superseded" not in p.parts:
            print(f"  {p.relative_to(OUT)}  {p.stat().st_size / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
