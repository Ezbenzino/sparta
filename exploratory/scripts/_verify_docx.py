# -*- coding: utf-8 -*-
import zipfile, os
z = zipfile.ZipFile(r"D:\sparta\docs\manuscript_cbc_draft.docx")
xml = z.read("word/document.xml").decode("utf-8")
checks = [
    ("stromal co-targeting", "In silico stromal co-targeting"),
    ("null permutation", "Permutation control for the rank-normalised"),
    ("why a graph", "Why a graph at all"),
    ("single-author", "single-author work"),
    ("30% reduction", "30% lowers both barriers"),
    ("naive baseline", "0.5"),
    ("positioning", "methodological contribution on public data"),
    ("Yize Li", "Yize Li"),
    ("GitHub URL", "github.com/Ezbenzino/sparta"),
    ("Zenodo DOI", "10.5281/zenodo.23086431"),
    ("no funding", "no specific funding"),
    ("97.6 null", "97.6%"),
]
ok = 0
for label, kw in checks:
    found = kw in xml
    if found: ok += 1
    print(f"  [{'OK' if found else 'MISSING'}] {label}")
print(f"\n{ok}/{len(checks)} checks passed")
print(f"docx size: {os.path.getsize(r'D:\sparta\docs\manuscript_cbc_draft.docx')/1024:.0f} KB")
