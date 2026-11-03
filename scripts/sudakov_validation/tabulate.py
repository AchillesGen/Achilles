"""Tabulate compare_achilles.py output as Markdown. Usage: tabulate.py <json> <interaction label>"""
import json, sys
r = json.load(open(sys.argv[1])); lab = sys.argv[2]
algs = ["Veto", "B_k1_Atom", "Continuous", "B_k1_Box", "B_eik_Gauss", "A_rho", "A_fold_finiteA"]
for mode, key, fmt in (("CrossSection", "sigmaR", "{:.0f} ± {:.0f}"), ("Transparency", "T", "{:.3f} ± {:.3f}")):
    print(f"\n**{lab}, {'sigma_R [mb]' if key == 'sigmaR' else 'transparency'}**\n")
    print("| p [MeV/c] | " + " | ".join(algs) + " |")
    print("|---" * (len(algs) + 1) + "|")
    first = next(v for k, v in r.items() if mode in k and "result" in v)
    for i in range(len(first["result"])):
        row = [f"{first['result'][i]['p']:.0f}"]
        for a in algs:
            x = r.get(f"{sys.argv[3] if len(sys.argv) > 3 else lab}_{mode}_{a}", {})
            row.append("err" if "result" not in x else fmt.format(x["result"][i][key], x["result"][i]["err"]))
        print("| " + " | ".join(row) + " |")
