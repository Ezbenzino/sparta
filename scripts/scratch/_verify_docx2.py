# -*- coding: utf-8 -*-
import zipfile, os
z = zipfile.ZipFile(r"D:\sparta\docs\manuscript_cbc_draft.docx")
xml = z.read("word/document.xml").decode("utf-8")
checks = [
    ("radius sensitivity", "Graph-radius sensitivity"),
    ("65 configurations", "65 configurations"),
    ("+0.225", "+0.225"),
    ("biological alignment", "Alignment with measured cell distributions"),
    ("-0.022", "0.022"),
    ("+0.046", "0.046"),
    ("T/NK", "T/NK"),
    ("96.9%", "96.9%"),
    ("89.2%", "89.2%"),
]
ok = 0
for label, kw in checks:
    found = kw in xml
    if found: ok += 1
    print(f"  [{'OK' if found else 'MISSING'}] {label}")
print(f"\n{ok}/{len(checks)} passed, docx {os.path.getsize(r'D:\sparta\docs\manuscript_cbc_draft.docx')/1024:.0f} KB")
