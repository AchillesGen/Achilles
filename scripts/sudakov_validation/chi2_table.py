"""chi^2/N of p + 12C sigma_R (Carlson, digitised) and 12C transparency against data_runs.py output.

Usage: python3 chi2_table.py <carlson.csv> <transparency.csv> <runs.json> [<runs.json> ...]
Prints, per algorithm: chi^2/N and mean model/data for sigma_R below and above 100 MeV,
the same with the suspect 61.0 MeV point and the < 20 MeV points dropped, and chi^2/N for
transparency (any-interaction definition and, if available, the acceptance definition).
Use --json <file> to also write the numbers."""
import json, math, sys
import numpy as np

MN = 938.92
EL = [8, 12, 16, 20, 25, 30, 36, 43, 52, 62, 75, 90, 105]
EH = [105, 130, 160, 200, 250, 310, 380, 460, 560, 680, 820, 980]
ET = [140, 220, 300, 380, 460, 560, 680, 820, 980, 1150, 1350, 1550, 1800, 1950]


def T_of_p(p):
    return np.sqrt(p**2 + MN**2) - MN


def regroup(e, edges, mode):
    lo, hi = np.array(e["lo"]), np.array(e["hi"])
    h, n = np.array(e["hits"]), np.array(e["nohits"])
    acc = np.array(e["accepted"]) if mode == "a" else None
    esc = np.array(e["escape"]) if mode == "e" else None
    area = 10 * math.pi * e.get("radius", 7.0) ** 2  # mb
    Tc = T_of_p(0.5 * (lo + hi))
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (Tc >= a) & (Tc < b)
        H, N = h[s].sum(), n[s].sum()
        t = H + N
        if t < 50:
            continue
        Tm = np.average(Tc[s], weights=(h + n)[s])
        if mode == "a":
            f = acc[s].sum() / t
            out.append((Tm, f, math.sqrt(f * (1 - f) / t)))
            continue
        if mode == "e":
            f = esc[s].sum() / t
            out.append((Tm, area * f, area * math.sqrt(f * (1 - f) / t)))
            continue
        f = H / t
        er = math.sqrt(f * (1 - f) / t)
        out.append((Tm, area * f, area * er) if mode == "s" else (Tm, 1 - f, er))
    return np.array(out)


def load(fn):
    return np.array([[float(x) for x in l.split(",")] for l in open(fn)
                     if not (l.startswith("#") or l.startswith("Tp"))])


def chi2(data, curve, sel, logx=True):
    x = np.log(data[:, 0]) if logx else data[:, 0]
    cx = np.log(curve[:, 0]) if logx else curve[:, 0]
    m = np.interp(x, cx, curve[:, 1], left=np.nan, right=np.nan)
    me = np.interp(x, cx, curve[:, 2])
    k = sel & ~np.isnan(m)
    c2 = np.sum((data[k, 1] - m[k]) ** 2 / (data[k, 2] ** 2 + me[k] ** 2)) / k.sum()
    return float(c2), float(np.mean(m[k] / data[k, 1])), int(k.sum())


def main():
    args = sys.argv[1:]
    jout = None
    if "--json" in args:
        i = args.index("--json")
        jout = args[i + 1]
        del args[i:i + 2]
    car, tra = load(args[0]), load(args[1])
    low, high = car[:, 0] < 100, car[:, 0] >= 100
    robust = low & (car[:, 0] >= 20) & ~((car[:, 0] > 60.9) & (car[:, 0] < 61.1))
    table = {}
    for fn in args[2:]:
        r = json.load(open(fn))
        for a in sorted(set(v["alg"] for v in r.values() if "alg" in v)):
            row = {}
            if "lo" in r.get(f"sigR_low_{a}", {}) and "lo" in r.get(f"sigR_high_{a}", {}):
                g = np.vstack([regroup(r[f"sigR_low_{a}"], EL, "s"), regroup(r[f"sigR_high_{a}"], EH, "s")])
                row["low"] = chi2(car, g, low)
                row["low_robust"] = chi2(car, g, robust)
                row["high"] = chi2(car, g, high)
                row["sigR"] = [(round(t, 1), round(s, 1), round(e, 1)) for t, s, e in g]
                if "escape" in r[f"sigR_low_{a}"] and "escape" in r[f"sigR_high_{a}"]:
                    ge = np.vstack([regroup(r[f"sigR_low_{a}"], EL, "e"), regroup(r[f"sigR_high_{a}"], EH, "e")])
                    row["low_esc"] = chi2(car, ge, low)
                    row["high_esc"] = chi2(car, ge, high)
                    row["sigR_esc"] = [(round(t, 1), round(s, 1), round(e, 1)) for t, s, e in ge]
            if "lo" in r.get(f"transp_{a}", {}):
                tr = r[f"transp_{a}"]
                g = regroup(tr, ET, "t")
                row["T"] = chi2(tra, g, np.ones(len(tra), bool), logx=False)
                row["Tcurve"] = [(round(t, 1), round(s, 4), round(e, 4)) for t, s, e in g]
                if "accepted" in tr:
                    g = regroup(tr, ET, "a")
                    row["Tacc"] = chi2(tra, g, np.ones(len(tra), bool), logx=False)
                    row["Tacc_curve"] = [(round(t, 1), round(s, 4), round(e, 4)) for t, s, e in g]
            table[f"{fn}:{a}"] = row
            fmt = lambda c: f"{c[0]:6.1f} ({c[1]:.2f})" if c else "      -     "
            print(f"{a:28s} sR<100 {fmt(row.get('low'))}  robust {fmt(row.get('low_robust'))}  "
                  f"sR>=100 {fmt(row.get('high'))}  T {fmt(row.get('T'))}  Tacc {fmt(row.get('Tacc'))}"
                  + (f"  | escape-def sR<100 {fmt(row.get('low_esc'))} sR>=100 {fmt(row.get('high_esc'))}" if 'low_esc' in row else "")
                  + f"   [{fn.split('/')[-1]}]")
    if jout:
        json.dump(table, open(jout, "w"), indent=1)


if __name__ == "__main__":
    main()
