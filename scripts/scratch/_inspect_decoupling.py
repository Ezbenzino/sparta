import json
o = json.load(open(r'D:\sparta\results\validation\decoupling.json', encoding='utf-8'))
print(json.dumps(o['MEL01'], indent=2, ensure_ascii=False)[:1500])
print('---')
print(json.dumps(o['CSCC14'], indent=2, ensure_ascii=False)[:1500])
