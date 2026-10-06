import gzip, json
raw = gzip.open(r"D:\sparta\data\external\GSE78220_matrix.txt.gz", "rt", encoding="utf-8", errors="replace").read()
lines = raw.splitlines()
gsms = resp = titles = None
for l in lines:
    if l.startswith("!Sample_geo_accession"):
        gsms = [p.strip().strip('"') for p in l.split("\t")[1:]]
    if l.startswith("!Sample_characteristics") and "response" in l:
        resp = [p.split(": ")[-1].strip().strip('"') for p in l.split("\t")[1:]]
    if l.startswith("!Sample_title"):
        titles = [p.strip().strip('"') for p in l.split("\t")[1:]]
print("n:", len(gsms))
meta = {}
for t, g, r in zip(titles, gsms, resp):
    meta[t] = {"gsm": g, "response": r}
    print(t, g, r)
json.dump(meta, open(r"D:\sparta\data\external\GSE78220_meta.json", "w"), indent=1)
from collections import Counter
print(Counter(meta[m]["response"] for m in meta))
