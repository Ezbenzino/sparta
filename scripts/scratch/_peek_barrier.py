"""查看 barrier npz 与 mincut.json 结构，为割带特征分析做准备。"""
import numpy as np, json
for sid in ["MEL01","MEL02","CSCC03","CSCC04"]:
    npz = np.load(rf"D:\sparta\data\interim\{sid}.barrier.npz", allow_pickle=True)
    print("="*50, sid)
    print("npz keys:", list(npz.keys()))
    for k in npz.keys():
        a = npz[k]
        print(f"  {k}: shape={getattr(a,'shape',None)} dtype={getattr(a,'dtype',None)}")
    try:
        m = json.load(open(rf"D:\sparta\data\interim\{sid}.mincut.json", encoding="utf-8"))
        print("mincut.json keys:", list(m.keys()))
        for k, v in m.items():
            if isinstance(v, list):
                print(f"  {k}: list[{len(v)}] 前5个={v[:5]}")
            else:
                print(f"  {k}: {v}")
    except Exception as e:
        print("mincut err:", e)
    break
