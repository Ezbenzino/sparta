import json
for f in [r"D:\sparta\results\counterfactual\MEL01.json",
          r"D:\sparta\results\validation\decoupling.json"]:
    d = json.load(open(f, encoding="utf-8"))
    print("="*40, f.split("\\")[-1])
    print("top keys:", list(d.keys()))
    if "s1" in d:
        s1 = d["s1"]
        print("s1 type:", type(s1).__name__)
        if isinstance(s1, dict):
            print("s1 keys:", list(s1.keys()))
            for k, v in s1.items():
                if isinstance(v, dict):
                    print("  ", k, "->", {kk: (round(vv,3) if isinstance(vv,float) else vv) for kk,vv in v.items() if kk in ('mode','z','p_emp','b_real','n_perm')})
