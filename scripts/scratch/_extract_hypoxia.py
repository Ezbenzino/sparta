"""提取 HALLMARK_HYPOXIA 基因集为 gmt 文件。"""
import json, zipfile
zp = r"D:\sparta\data\external\genesets.json.zip"
z = zipfile.ZipFile(zp)
name = [n for n in z.namelist() if n.endswith(".json")][0]
df = json.load(z.open(name))
h = df.get("HALLMARK_HYPOXIA", {})
genes = h.get("symbols", [])
print(f"HALLMARK_HYPOXIA: {len(genes)} genes (first 10: {genes[:10]})")
# 写 gmt（一行：name, desc, genes...）
gmt = "HALLMARK_HYPOXIA\tHALLMARK_HYPOXIA\t" + "\t".join(genes) + "\n"
out = r"D:\sparta\data\external\h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt"
open(out, "w", encoding="utf-8").write(gmt)
print("saved", out)
