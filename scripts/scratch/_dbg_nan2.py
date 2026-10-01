import gzip, json, sys
import numpy as np
import pandas as pd
sys.path.insert(0, r"D:\sparta")
S2E = json.load(open(r"D:\sparta\data\external\GSE91061_sym2entrez.json", encoding="utf-8"))
df = pd.read_csv(r"D:\sparta\data\external\GSE91061_rld.csv.gz", index_col=0)
df.columns = [c.strip().strip('"') for c in df.columns]
df.index = df.index.astype(str).str.strip().str.strip('"')
print("df shape:", df.shape)
print("dup index:", df.index.duplicated().sum())
print("dup cols:", df.columns.duplicated().sum())
zdf = (df - df.mean(axis=1)) / df.std(axis=1)
print("zdf shape:", zdf.shape)
ezs = [str(S2E["LOX"]), str(S2E["LOXL1"])]
print("LOX ez:", ezs)
try:
    sub = zdf.loc[ezs]
    print("sub shape:", sub.shape)
except Exception as e:
    print("loc err:", e)
# 用 iloc 测试
print("zdf iloc[:2,:3]:")
print(zdf.iloc[:2, :3])
