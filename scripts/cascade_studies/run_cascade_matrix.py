"""Run achilles-cascade over a matrix of algorithms and nuclear configurations.

For each run the hit / no-hit histograms (in probe momentum) give
  Transparency mode: T(p) = nohits / (hits + nohits)       (kicked nucleon from a vertex)
  CrossSection mode: sigma_R(p) = pi R_beam^2 hits / total (external proton beam)

Run from the build directory (it needs data/ and bin/achilles-cascade):
  python3 <repo>/scripts/cascade_studies/run_cascade_matrix.py <build_dir> <out_dir> [study]
"""

import json
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

RUNCARD = """Initialize:
  seed: {seed}
Nucleus:
  Name: 12C
  Binding: 8.6
  Fermi Momentum: 225
  Density:
    ProtonFile: data/densities/c12.prova.txt
    NeutronFile: data/densities/c12.prova.txt
    Configs: {configs}
  FermiGas:
    Type: Local
    Correlated: False
    SRCfraction: 0.2
    LambdaSRC: 2.75
    Params: []
  Potential:
    Name: Schroedinger
    r0: 0.16
    Mode: 3
Cascade:
  Mode: {mode}
  Interactions:
    - Name: GeantInteraction
      Options:
        GeantData: data/GeantData.hdf5
  Step: {step}
  Probability: {prob}
  InMedium: None
  PotentialProp: False
  Algorithm: {alg}
  Params:
    radius: {radius}
    external: 0
KickMomentum: [{pmin}, {pmax}]
NEvents: {nevents}
PID: 2212
Output:
  Format: Achilles
  Name: {name}.hepmc
  Zipped: True
"""

BEAM_RADIUS = 7.0  # fm, CrossSection mode
PBINS = [200, 400, 600, 800, 1000, 1200, 1400]  # MeV, probe momentum


def read_hist(path):
    rows = np.loadtxt(path, skiprows=2)
    return rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]


def counts(value, error):
    """Event counts per bin: all events share one weight w, so N = (N w)^2 / (N w^2)."""
    return np.where(error > 0, value**2 / np.where(error > 0, error**2, 1), 0.0)


def summarize(name, mode):
    lo, hi, hv, he = read_hist(f"{name}_hits.txt")
    _, _, nv, ne = read_hist(f"{name}_nohits.txt")
    hits, nohits = counts(hv, he), counts(nv, ne)
    out = []
    for a, b in zip(PBINS[:-1], PBINS[1:]):
        sel = (lo >= a) & (hi <= b)
        h, n = hits[sel].sum(), nohits[sel].sum()
        tot = h + n
        frac, err = h / tot, np.sqrt(h * n / tot) / tot
        if mode == "Transparency":
            out.append({"p": 0.5 * (a + b), "T": 1 - frac, "err": err})
        else:
            area = np.pi * BEAM_RADIUS**2 * 10  # mb
            out.append({"p": 0.5 * (a + b), "sigmaR": area * frac, "err": area * err})
    return out


def run_one(build, outdir, spec):
    name = outdir / spec["tag"]
    card = outdir / f"{spec['tag']}.yml"
    card.write_text(RUNCARD.format(name=name, pmin=PBINS[0], pmax=PBINS[-1], radius=BEAM_RADIUS, **spec["card"]))
    t0 = time.time()
    res = subprocess.run([str(build / "bin" / "achilles-cascade"), str(card)], cwd=build,
                         capture_output=True, text=True)
    if res.returncode != 0:
        return spec["tag"], {"error": res.stderr[-2000:]}
    # the event file is not needed for these observables
    for f in outdir.glob(f"{spec['tag']}.hepmc*"):
        f.unlink()
    return spec["tag"], {"spec": spec, "seconds": time.time() - t0,
                         "result": summarize(name, spec["card"]["mode"])}


def matrix(study, cfgdir):
    qmc = "data/configurations/QMC_configs.out.gz"
    specs = []

    def add(tag, **card):
        base = {"seed": 12345, "step": 0.04, "prob": "Gaussian", "alg": "Base",
                "configs": qmc, "nevents": 60000}
        base.update(card)
        specs.append({"tag": tag, "card": base})

    if study in ("steps", "all"):
        for mode in ("Transparency", "CrossSection"):
            for step in (0.02, 0.04, 0.1, 0.25, 0.5):
                add(f"steps_{mode}_Base_{step}", mode=mode, step=step)
            add(f"steps_{mode}_Veto", mode=mode, alg="Veto")
    if study in ("configs", "all"):
        ens = {"QMC": qmc, "MF": "data/configurations/MF_configs.out.gz",
               "Mixed": f"{cfgdir}/12C_mixed_configs.out.gz",
               "Jastrow": f"{cfgdir}/12C_jastrow_configs.out.gz"}
        for mode in ("Transparency", "CrossSection"):
            for prob in ("Gaussian", "Cylinder"):
                for e, path in ens.items():
                    for alg in ("Veto", "Continuous"):
                        add(f"configs_{mode}_{prob}_{e}_{alg}", mode=mode, prob=prob, alg=alg,
                            configs=path, nevents=40000)
    return specs


def main():
    build = pathlib.Path(sys.argv[1]).resolve()
    outdir = pathlib.Path(sys.argv[2]).resolve()
    study = sys.argv[3] if len(sys.argv) > 3 else "all"
    cfgdir = sys.argv[4] if len(sys.argv) > 4 else str(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    specs = matrix(study, cfgdir)
    results = {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for tag, res in pool.map(lambda s: run_one(build, outdir, s), specs):
            results[tag] = res
            print(tag, json.dumps(res.get("result", res))[:300], flush=True)
    (outdir / f"matrix_{study}.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
