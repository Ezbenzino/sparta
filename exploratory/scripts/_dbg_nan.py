import gzip, json, sys
import numpy as np
import pandas as pd
sys.path.insert(0, r"D:\sparta")
S2E = json.load(open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", encoding="utf-8"))
df = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df.columns = [c.strip().strip('"') for c in df.columns]
df.index = df.index.astype(str).str.strip().str.strip('"')
print("df head index:", list(df.index[:3]))
zdf = (df - df.mean(axis=1)) / df.std(axis=1)
print("zdf any NaN:", zdf.isna().sum().sum(), "of", zdf.size)
# ECM_crosslink
syms = ["LOX", "LOXL1", "LOXL2", "LOXL3", "PLOD1", "PLOD2", "PLOD3"]
ezs = [S2E[s] for s in syms if s in S2E and str(S2E[s]) in zdf.index]
print("ezs:", ezs)
print("match check:", [str(e) in zdf.index for e in ezs])
sub = zdf.loc[[str(e) for e in ezs]]
print("sub shape:", sub.shape, "sub NaN:", sub.isna().sum().sum())
s = sub.mean(axis=0)
print("score NaN:", s.isna().sum(), "first values:", s.values[:5])
