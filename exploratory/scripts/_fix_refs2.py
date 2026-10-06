# -*- coding: utf-8 -*-
import io
p = r"D:\sparta\docs\manuscript_cbc_draft.md"
lines = io.open(p, encoding="utf-8").read().split("\n")

out = []
i = 0
while i < len(lines):
    ln = lines[i]
    # [14] 块：从 "[14] ⚠" 行起，替换三行
    if ln.startswith("[14] ⚠"):
        out.append('[14] R.K. Jain, Delivery of molecular and cellular medicine to tumors,')
        out.append('     Nature Reviews Drug Discovery 4 (2005) 619-632.')
        i += 3  # 跳过原来的 3 行
        continue
    # 删除已过时的编辑注第 3 条
    if ln.strip().startswith("> 3. 参考文献中标"):
        i += 1
        continue
    out.append(ln)
    i += 1

io.open(p, "w", encoding="utf-8", newline="\n").write("\n".join(out))
# 复查
s = io.open(p, encoding="utf-8").read()
left = [l for l in s.splitlines() if "⚠" in l]
print("leftover ⚠ lines:", len(left))
for l in left:
    print("  ", l)
