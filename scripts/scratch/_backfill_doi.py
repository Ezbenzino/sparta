# -*- coding: utf-8 -*-
import io

# 1) 稿件 Data availability
p = r"D:\sparta\docs\manuscript_cbc_draft.md"
s = io.open(p, encoding="utf-8").read()
before = s.count("待回填")
s = s.replace("[Zenodo DOI — 待回填]",
              "https://doi.org/10.5281/zenodo.23086431")
io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("manuscript: 待回填 count", before, "->", s.count("待回填"))

# 2) CITATION.cff 加 repository-artifact (Zenodo)
p2 = r"D:\sparta\CITATION.cff"
c = io.open(p2, encoding="utf-8").read()
if "repository-artifact" not in c:
    c = c.replace('repository-code: "https://github.com/Ezbenzino/sparta"',
                  'repository-code: "https://github.com/Ezbenzino/sparta"\nrepository-artifact: "https://doi.org/10.5281/zenodo.23086431"')
io.open(p2, "w", encoding="utf-8", newline="\n").write(c)
print("CITATION.cff updated")
print("---- CITATION tail ----")
print("\n".join(c.splitlines()[-8:]))
