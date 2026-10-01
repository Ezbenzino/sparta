# -*- coding: utf-8 -*-
import io
p = r"D:\sparta\docs\manuscript_cbc_draft.md"
s = io.open(p, encoding="utf-8").read()

# CRediT 单作者
s = s.replace(
"""**[作者姓名]**: Conceptualisation, Methodology, Software, Validation, Formal
analysis, Investigation, Data curation, Writing — original draft,
Visualisation. *(单作者声明；若导师或其他贡献者列名，按实际 CRediT 角色补)*""",
"""**Yize Li**: Conceptualisation, Methodology, Software, Validation, Formal
analysis, Investigation, Data curation, Writing — original draft,
Visualisation.""")

# 致谢：无基金
s = s.replace("**[基金号/致谢待补]**",
              "The author received no specific funding for this work.")

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("placeholders left:", "[作者" in s or "待补" in s or "待 git" in s)
print("CRediT has Yize Li:", "Yize Li" in s)
