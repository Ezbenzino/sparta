import json
d = json.load(open(r"D:\sparta\results\validation\screen_decision.json", encoding="utf-8"))
print("verdict:", d["verdict"], " median:", d["median_potential"])
for sid, r in d["per_slide"].items():
    print(f"{sid}: potential={r['dissociation_potential']:.2f} "
          f"ECM%={r['frac_shared_ecm']*100:.1f} antigen%={r['frac_antigen']*100:.1f} "
          f"crosslink%={r['frac_crosslink']*100:.1f} collinear={r['crosslink_collinear']} "
          f"n_spots={r['n_spots']}")
