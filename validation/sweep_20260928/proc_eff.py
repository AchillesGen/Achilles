"""Per-process 1% TailFraction efficiency, as Achilles unweights (one cap per process)."""
import sys, numpy as np, pyhepmc

def cap_of(w, eps=0.01):
    s = np.sort(w); suffix = np.cumsum(s[::-1])[::-1]
    ok = suffix <= eps * s.sum()
    i = int(np.argmax(ok)) if ok.any() else len(s)
    return s[i - 1] if i > 0 else s[0]

def per_process(path):
    procs = {}
    with pyhepmc.open(path) as f:
        for e in f:
            prim = [v for v in e.vertices if v.status == 1][0]
            key = tuple(sorted(p.pid for p in prim.particles_out if abs(p.pid) != 13))
            procs.setdefault(key, []).append(abs(e.weight("CV")))
    num = den = 0.0; parts = {}
    for k, w in procs.items():
        w = np.array(w); c = cap_of(w)
        parts[k] = (len(w), np.minimum(w, c).mean() / c)
        num += np.minimum(w, c).mean(); den += c
    return num / den, parts

if __name__ == "__main__":
    for path in sys.argv[1:]:
        eff, parts = per_process(path)
        print(f"{path.split('/')[-1]:24s} eff1={eff*100:6.2f}%  " +
              "  ".join(f"{k}:n={n},eff={e*100:.1f}%" for k, (n, e) in sorted(parts.items())))
