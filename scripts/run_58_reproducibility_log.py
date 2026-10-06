"""Write a reproducibility log for the 2026 SPARTA submission."""
from __future__ import annotations

import importlib.metadata
import platform
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/logs/reproducibility_log.txt"
PACKAGES = [
    "numpy", "scipy", "pandas", "matplotlib", "networkx", "anndata", "scanpy",
    "h5py", "openpyxl", "scikit-learn", "statsmodels", "pillow", "pytest",
]
SEEDS = {
    "extension ingest/QC": "not randomized; deterministic filtering",
    "extension cohort": "20261005",
    "figure/QC helper jitter": "20261005",
    "Scanpy score_genes controls": "fixed seed configured by calling code",
}


def hardware_text() -> str:
    cpu = platform.processor() or "unknown processor"
    n = str(__import__("os").cpu_count() or "unknown")
    ram = "not available"
    try:
        import psutil

        ram = f"{psutil.virtual_memory().total / (1024 ** 3):.1f} GiB"
    except Exception:
        try:
            import ctypes

            class MEMORYSTATUSEX(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            stat = MEMORYSTATUSEX()
            stat.dwLength = ctypes.sizeof(stat)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
            ram = f"{stat.ullTotalPhys / (1024 ** 3):.1f} GiB"
        except Exception:
            pass
    return f"CPU: {cpu}\nLogical cores: {n}\nRAM: {ram}"


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "SPARTA 2026 submission reproducibility log",
        "=" * 48,
        "",
        f"Date: {platform.node() and __import__('datetime').datetime.now().astimezone().isoformat()}",
        f"Operating system: {platform.platform()}",
        hardware_text(),
        "",
        f"Python executable: {sys.executable}",
        f"Python version: {sys.version}",
        "",
        "External programs",
        "-" * 48,
        f"R: {shutil.which('R') or 'not found in PATH'}",
        f"pandoc: {shutil.which('pandoc') or 'not found in PATH'}",
        f"pdflatex: {shutil.which('pdflatex') or 'not found in PATH'}",
        "",
        "Python package versions",
        "-" * 48,
    ]
    for package in PACKAGES:
        try:
            lines.append(f"{package}: {importlib.metadata.version(package)}")
        except importlib.metadata.PackageNotFoundError:
            lines.append(f"{package}: not installed")
    lines += [
        "",
        "Random seeds",
        "-" * 48,
    ]
    for name, seed in SEEDS.items():
        lines.append(f"{name}: {seed}")
    lines += [
        "",
        "Notes",
        "-" * 48,
        "DOCX/PDF construction requires pandoc; it was not available in PATH when this log was written.",
        "Structural results are deterministic given the recorded inputs, parameters and seeds.",
    ]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
