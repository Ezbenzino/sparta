import json
o = json.load(open(r'D:\sparta\results\validation\ext_validation.json', encoding='utf-8'))
for sid, r in o['per_slide'].items():
    d = r['decoupling']; s = r['spatial_null']
    print(f"{sid:<8} rho_partial={d['rho_partial']:+.3f} p={d['p_partial']:.2e}  "
          f"null_real={s['real_rho_partial']:+.3f} null_p={s['empirical_p_one_sided']:.3f}  "
          f"filled={r['filled_keys']}")
