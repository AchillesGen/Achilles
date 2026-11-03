"""Step 4(h) runs: p + 12C sigma_R(T) and 12C transparency T(T_p) for comparison with data.
Writes per-histogram-bin hit / no-hit counts so the plotting script can regroup in kinetic
energy.  Usage: python3 data_runs.py <build_dir> <out_dir> [nevents] [algs] [json] [NN|Geant|NASA]
(env RANGES=transp,... restricts the ranges)"""
import json, os, pathlib, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import compare_achilles as ca

INTER = "NN"  # or "Geant" (6th command-line argument)
# beam radius (fm) for CrossSection runs; env BEAM_RADIUS. Achilles Base needs >= 14 fm below
# ~40 MeV (its Gaussian pair probability follows the divergent low-energy NN cross section).
RADIUS = float(os.environ.get("BEAM_RADIUS", ca.BEAM_RADIUS))
ALGS = ["Continuous", "B_k1_Box", "B_eik_Gauss", "A_fold_finiteA", "A_rho"]
RANGES = {  # (mode, pmin, pmax) in MeV/c; uniform in p, so split the low-energy region
    "sigR_low": ("CrossSection", 130, 470),
    "sigR_high": ("CrossSection", 400, 1830),
    "transp": ("Transparency", 550, 2700),
}


def counts_per_bin(name):
    lo, hi, hv, he = ca.read_hist(f"{name}_hits.txt")
    _, _, nv, ne = ca.read_hist(f"{name}_nohits.txt")
    return lo.tolist(), hi.tolist(), ca.counts(hv, he).tolist(), ca.counts(nv, ne).tolist()


def run(build, outdir, tag, mode, pmin, pmax, alg_key, nevents):
    alg, sud = ca.ALGOS[alg_key.replace("_pot", "").replace("_cap1000", "")]
    name = outdir / tag
    card = outdir / f"{tag}.yml"
    if INTER in ("Geant", "NASA"):
        inter = ca.INTERACTIONS[INTER]
    else:
        lines = (build / "data" / "default" / "VirtResInteractions.yml").read_text().splitlines()
        inter = "\n".join("    " + l for l in lines if l.strip())
    text = (ca.RUNCARD.format(seed=2468, mode=mode, interactions=inter, alg=alg,
                                      sudakov=sud.replace("{cap}", str(ca.CAP[INTER])),
                                      radius=RADIUS, pmin=pmin, pmax=pmax,
                                      nevents=nevents, name=name))
    if "_cyl" in alg_key:
        text = text.replace("Probability: Gaussian", "Probability: Cylinder")
    if "_med" in alg_key:
        text = text.replace("InMedium: None", "InMedium: NonRelativistic")
    if "{ipm}" in text:
        import ipm_veto
        d = np.loadtxt(build / "data" / "densities" / "c12.prova.txt", comments="#")
        qs, vs, _, _ = ipm_veto.table(d[:, 0], d[:, 1], 600.0, 10.0)
        text = text.replace("{ipm}", "    CoherentTable:\n      q: [" + ", ".join(f"{q:g}" for q in qs)
                            + "]\n      v: [" + ", ".join(f"{v:.4f}" for v in vs) + "]\n")
    if "_pot" in alg_key:  # mean-field propagation (Base/MFP/Veto/Continuous only)
        text = text.replace("PotentialProp: False", "PotentialProp: True")
    if "_cap1000" in alg_key:
        text = text.replace(f"SigmaCap: {ca.CAP[INTER]}", "SigmaCap: 1000")
    card.write_text(text)
    t0 = time.time()
    res = subprocess.run([str(build / "bin" / "achilles-cascade"), str(card)], cwd=build,
                         capture_output=True, text=True)
    for f in outdir.glob(f"{tag}.hepmc*"):
        f.unlink()
    if res.returncode != 0:
        return tag, {"error": res.stderr[-1500:]}
    lo, hi, h, n = counts_per_bin(name)
    out = {"mode": mode, "alg": alg_key, "seconds": time.time() - t0,
           "lo": lo, "hi": hi, "hits": h, "nohits": n}
    out["radius"] = RADIUS
    esc = pathlib.Path(f"{name}_escape.txt")
    if esc.exists():  # sigma_R with the 2020-code hit definition (CrossSection mode)
        _, _, ev, ee = ca.read_hist(esc)
        out["escape"] = ca.counts(ev, ee).tolist()
    acc = pathlib.Path(f"{name}_accepted.txt")
    if acc.exists():  # acceptance-defined transparency (Transparency mode)
        _, _, av, ae = ca.read_hist(acc)
        out["accepted"] = ca.counts(av, ae).tolist()
    return tag, out


def main():
    build = pathlib.Path(sys.argv[1]).resolve()
    outdir = pathlib.Path(sys.argv[2]).resolve()
    nevents = int(sys.argv[3]) if len(sys.argv) > 3 else 100000
    algs = sys.argv[4].split(",") if len(sys.argv) > 4 else ALGS
    jsonname = sys.argv[5] if len(sys.argv) > 5 else "data_runs.json"
    global INTER
    INTER = sys.argv[6] if len(sys.argv) > 6 else "NN"
    outdir.mkdir(parents=True, exist_ok=True)
    ranges = os.environ.get("RANGES", ",".join(RANGES)).split(",")
    jobs = [(f"{rk}_{a}", *RANGES[rk], a) for rk in ranges for a in algs]
    results = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for tag, res in pool.map(lambda j: run(build, outdir, j[0], j[1], j[2], j[3], j[4], nevents), jobs):
            results[tag] = res
            print(tag, res.get("seconds") or res.get("error", "")[:300], flush=True)
            (outdir / jsonname).write_text(json.dumps(results))


if __name__ == "__main__":
    main()
