"""
Shared figure style for the IS submission (Springer artwork guidelines)
=======================================================================
* Widths: 84 mm (one column) or 174 mm (full width); height <= 234 mm.
* Lettering: sans-serif (Arial metrics; Liberation Sans here), 8-9 pt at final size,
  panel labels 10 pt bold lower-case (a, b, c ...), no titles inside the artwork
  (captions live in the manuscript text).
* Lines >= 0.3 pt.  No transparency in vector output (EPS cannot carry alpha).
* Colour carries ONE meaning across the paper: which barrier operator.
      cell / B_cell     #2a78d6  (blue)
      antibody / B_mAb  #eb6834  (orange)
  Cohorts / platforms are encoded by marker SHAPE and position, never by colour.
  Statistical significance is encoded by FILLED vs OPEN markers.
  Null / reference distributions are grey.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MM = 1 / 25.4
FULL_W = 174 * MM
COL_W = 84 * MM
MAX_H = 234 * MM

CELL = "#2a78d6"
MAB = "#eb6834"
AQUA = "#1baf7a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
LIGHT = "#c3c2b7"
HAIR = "#e1e0d9"
WASH_CELL = "#cde2fb"   # blue step 100
WASH_MAB = "#fbe0d3"    # light orange wash (opaque)
WASH_GREY = "#f0efec"   # diverging midpoint grey

BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

# marker shape per platform / arm (shape, not colour, carries the cohort)
SHAPE = {"visium_cscc": "o", "legacy_st": "s", "visium_mel": "^", "external": "D", "legacy_mel": "v"}


def platform_key(sid: str) -> str:
    if sid.startswith("MEL_THR"):
        return "legacy_mel"
    if sid.startswith("MEL"):
        return "visium_mel"
    if sid in ("CSCC01", "CSCC02", "CSCC03", "CSCC04"):
        return "visium_cscc"
    if sid.startswith("CSCC"):
        return "legacy_st"
    return "external"


SHAPE_LABEL = {"visium_cscc": "cSCC, Visium", "legacy_st": "cSCC, first-generation ST",
               "visium_mel": "Melanoma, Visium", "external": "External, Visium",
               "legacy_mel": "Melanoma LN, first-generation ST"}

PATIENT_LABEL = {"CSCC_P2": "cSCC P2", "CSCC_P4": "cSCC P4", "CSCC_P5": "cSCC P5",
                 "CSCC_P6": "cSCC P6", "CSCC_P9": "cSCC P9", "CSCC_P10": "cSCC P10",
                 "MEL_PtB": "Melanoma PtB", "BRCA_P1": "Breast P1", "LN_P1": "Lymph node",
                 "MEL_THR1": "Melanoma LN P1", "MEL_THR2": "Melanoma LN P2", "MEL_THR3": "Melanoma LN P3",
                 "MEL_THR4": "Melanoma LN P4"}


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.6,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": HAIR,
        "grid.linewidth": 0.5,
        "legend.frameon": False,
        "legend.handlelength": 1.2,
        "legend.borderaxespad": 0.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
        "figure.dpi": 150,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
    })


def panel(ax, letter, x=-0.12, y=1.04, fig=None):
    """Lower-case bold panel label in the upper-left corner (Springer style)."""
    t = ax.transAxes if fig is None else fig.transFigure
    (ax if fig is None else fig).text(x, y, letter, transform=t, fontsize=10, fontweight="bold",
                                      va="bottom", ha="left", color=INK)


def save(fig, out_dir: Path, name: str, tiff=True, eps=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", pad_inches=0.01)
    if eps:
        fig.savefig(out_dir / f"{name}.eps", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(out_dir / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.01)
    if tiff:
        fig.savefig(out_dir / f"{name}.tif", dpi=600, bbox_inches="tight", pad_inches=0.01,
                    pil_kwargs={"compression": "tiff_lzw"})
        flatten_tiff(out_dir / f"{name}.tif")
    plt.close(fig)


def flatten_tiff(path: Path, dpi: int = 600):
    """Springer artwork: RGB without an alpha channel, LZW-compressed, resolution tag kept."""
    from PIL import Image
    im = Image.open(path)
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        bg.save(path, compression="tiff_lzw", dpi=(dpi, dpi))
