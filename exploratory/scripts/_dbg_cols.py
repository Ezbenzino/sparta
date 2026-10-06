import pandas as pd, json
df = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df.columns = [c.strip().strip('"') for c in df.columns]
meta = json.load(open(r"D:\sparta\data\external\GSE91061_meta.json", encoding="utf-8"))
pre = {t for t, m in meta.items() if "_Pre" in t and m["response"] != "UNK"}
cols = set(df.columns)
inter = cols & pre
print("rld cols:", len(cols), " pre keys:", len(pre), " intersect:", len(inter))
print("pre not in cols:", sorted(pre - cols)[:5])
print("cols not in pre:", sorted(cols - pre)[:5])
