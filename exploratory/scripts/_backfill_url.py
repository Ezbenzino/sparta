# -*- coding: utf-8 -*-
import io
# 1) CITATION.cff
p1 = r"D:\sparta\CITATION.cff"
s = io.open(p1, encoding="utf-8").read()
s = s.replace('https://github.com/PLACEHOLDER_USER/sparta',
              'https://github.com/Ezbenzino/sparta')
io.open(p1, "w", encoding="utf-8", newline="\n").write(s)
print("CITATION.cff updated")

# 2) 稿件 Data availability 段回填 GitHub URL
p2 = r"D:\sparta\docs\manuscript_cbc_draft.md"
t = io.open(p2, encoding="utf-8").read()
t = t.replace("[GitHub URL — 待 git 首次提交后回填]",
             "https://github.com/Ezbenzino/sparta")
io.open(p2, "w", encoding="utf-8", newline="\n").write(t)
print("manuscript Data availability updated; '待 git' left?",
      "待 git" in t, "| Zenodo DOI left?", "待回填" in t)
