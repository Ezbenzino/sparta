import json
d = json.load(open(r"D:\sparta\results\validation\screen_decision.json", encoding="utf-8"))
print("verdict:", d["verdict"], " median:", round(d["median_potential"],2))
for sid, r in d["per_slide"].items():
    print(f"{sid}: potential={r['dissociation_potential']:.2f}  ECM={r['frac_shared_ecm']*100:.1f}%  "
          f"antigen={r['frac_antigen']*100:.1f}%  crosslink={r['frac_crosslink']*100:.1f}%  "
          f"col={r['crosslink_collinear']}")
