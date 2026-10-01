"""维度③概念验证（1/2）：解析 GSE91061 响应标签 + mygene Entrez 映射缓存。"""
import gzip, json, re
import requests

# ---- 1. 解析 series matrix：Sample_title / GSM / response ----
raw = gzip.open(r"D:\sparta\data\external\GSE91061_series_matrix.txt.gz",
                "rt", encoding="utf-8", errors="replace").read()
lines = raw.splitlines()
def col(line):
    parts = line.split("\t")
    return [p.strip().strip('"') for p in parts[1:]]
titles = None; gsms = None; resp = None
for l in lines:
    if l.startswith("!Sample_title"): titles = col(l)
    elif l.startswith("!Sample_geo_accession"): gsms = col(l)
    elif l.startswith('!Sample_characteristics') and 'response' in l: resp = col(l)
assert titles and gsms and resp, "missing columns"
print(f"n_samples={len(gsms)}")
meta = {}
for t, g, r in zip(titles, gsms, resp):
    meta[t] = {"gsm": g, "response": r.split(": ")[-1]}
json.dump(meta, open(r"D:\sparta\data\external\GSE91061_meta.json", "w"), indent=1)
from collections import Counter
print("response counts:", Counter(m["response"] for m in meta.values()))

# ---- 2. mygene 查询：signatures 全部基因 symbol -> entrez ----
import sys; sys.path.insert(0, r"D:\sparta")
from sparta.signatures import GENE_SETS
genes = set()
for tt in GENE_SETS:
    for gs in GENE_SETS[tt].values():
        genes |= set(gs)
# 加 HALLMARK_HYPOXIA
for line in open(r"D:\sparta\data\external\h.all.v7.1.HALLMARK_HYPOXIA.symbols.gmt", encoding="utf-8"):
    genes |= set(line.rstrip("\n").split("\t")[2:])
genes = sorted(genes)
print(f"total unique signature genes: {len(genes)}")

# mygene POST 批量查询（每次 1000 内）
sym2entrez = {}
url = "https://mygene.info/v3/query"
for i in range(0, len(genes), 1000):
    chunk = genes[i:i+1000]
    r = requests.post(url, json={"q": chunk, "scopes": "symbol", "fields": "entrezgene", "species": "human"}, timeout=120)
    for hit in r.json():
        q = hit.get("query")
        ez = hit.get("entrezgene")
        if q and ez:
            sym2entrez[q] = int(ez)
print(f"mapped {len(sym2entrez)}/{len(genes)}")
json.dump(sym2entrez, open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", "w"))
print("saved meta + sym2entrez")
