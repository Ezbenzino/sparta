import json
o = json.load(open(r'D:\sparta\results\validation\s2_matched_ext.json', encoding='utf-8'))
print(json.dumps(o['summary'], indent=2))
print()
for s, r in o['per_slide'].items():
    print(f"{s:<8} ratio={r['ratio_vs_in_cut_matched']:.3f}x  "
          f"p={r['p_vs_in_cut_matched']:.3f}  k={r['k']}")
