"""Step 4(f): full-kinematics comparison of the Sudakov tracks with Achilles Veto/Continuous.

Runs achilles-cascade (CrossSection and Transparency modes, p + 12C, QMC configurations) for
each algorithm and summarises sigma_R(p) and T(p) in momentum bins.
Run from anywhere:  python3 compare_achilles.py <build_dir> <out_dir> [interaction] [nevents]
interaction: Geant (default; the Continuous study harness) or NN (Achilles NucleonNucleon)."""
import json, pathlib, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np

INTERACTIONS = {
    "Geant": """    - Name: GeantInteraction
      Options:
        GeantData: data/GeantData.hdf5""",
    # NASA parametrisation ("Tot" curves of Isaacson et al., PRC 103, 015502)
    "NASA": """    - Name: NasaInteraction
      Options: {}""",
    # Achilles default interaction set (data/default/VirtResInteractions.yml), filled in main()
    "NN": None,
}

RUNCARD = """Initialize:
  seed: {seed}
Nucleus:
  Name: 12C
  Binding: 8.6
  Fermi Momentum: 225
  Density:
    ProtonFile: data/densities/c12.prova.txt
    NeutronFile: data/densities/c12.prova.txt
    Configs: data/configurations/QMC_configs.out.gz
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
{interactions}
  Step: 0.04
  Probability: Gaussian
  InMedium: None
  PotentialProp: False
  Algorithm: {alg}
{sudakov}  Params:
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
BEAM_RADIUS = 7.0
# NN cross-section cap for the Sudakov algorithms (mb; -1: none). See AchillesCollision.
CAP = {"Geant": -1, "NN": 500, "NASA": 500}
PBINS = [200, 400, 600, 800, 1000, 1200, 1400]

_EIK = "    Law: Eikonal\n    X0: 2.6\n    Profile: Gauss\n    Headroom: 1.3\n"
_AFF = "    Shells: 16\n    FoldWidth: 0.8\n    FiniteA: true\n"
_CAP = "    SigmaCap: {cap}\n"
_COH = "    CoherentVeto: true\n"
# Keys: a "_cyl" in the key sets Probability: Cylinder (Base only) and a "_med" sets
# InMedium: NonRelativistic (data_runs.py). PauliAttempt defaults to Single in the engine;
# the cards that predate that default say Retry explicitly so earlier results reproduce.
ALGOS = {
    "Base": ("Base", ""),
    "Base_cyl": ("Base", ""),
    "Base_med": ("Base", ""),
    "Base_cyl_med": ("Base", ""),
    "Veto": ("Veto", ""),
    "Continuous": ("Continuous", ""),
    # identity settings: same pair law (kappa = 1), same time profile, hard Pauli cut
    "B_k1_Atom": ("SudakovB", "  Sudakov:\n    Law: Kappa1\n    Profile: Atom\n    Pauli: Step\n    PauliAttempt: Retry\n" + _CAP),
    "B_k1_Box": ("SudakovB", "  Sudakov:\n    Law: Kappa1\n    Profile: Box\n    Pauli: Step\n    PauliAttempt: Retry\n" + _CAP),
    # proposed Track B: eikonal-exponent law, isotropic Gaussian profile
    "B_eik_Gauss": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: Step\n    PauliAttempt: Retry\n" + _CAP),
    # Track B with Achilles' Pauli prescription (single attempt per pair; now the default)
    "B_k1_Box_single": ("SudakovB", "  Sudakov:\n    Law: Kappa1\n    Profile: Box\n    Pauli: Step\n    PauliAttempt: Single\n" + _CAP),
    "B_eik_Gauss_single": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: Step\n    PauliAttempt: Single\n" + _CAP),
    # Pauli blocking off (Fermi motion kept)
    "B_eik_Gauss_nopauli": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: None\n" + _CAP),
    # coherent F^2(q) veto on elastic collisions (single attempt)
    "B_eik_Gauss_single_coh": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: Step\n" + _COH + _CAP),
    "B_eik_Gauss_single_coh_med": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: Step\n" + _COH + _CAP),
    "B_k1_Box_single_coh": ("SudakovB", "  Sudakov:\n    Law: Kappa1\n    Profile: Box\n    Pauli: Step\n" + _COH + _CAP),
    # Track A: density backbone, plain and with the two-body corrections that matter (2.E)
    "A_rho": ("SudakovA", "  Sudakov:\n    Shells: 16\n    Pauli: Step\n" + _CAP),
    "A_fold_finiteA": ("SudakovA", "  Sudakov:\n" + _AFF + "    Pauli: Step\n" + _CAP),
    "A_fold_finiteA_nopauli": ("SudakovA", "  Sudakov:\n" + _AFF + "    Pauli: None\n" + _CAP),
    "A_fold_finiteA_coh": ("SudakovA", "  Sudakov:\n" + _AFF + "    Pauli: Step\n" + _COH + _CAP),
    "A_fold_finiteA_coh_med": ("SudakovA", "  Sudakov:\n" + _AFF + "    Pauli: Step\n" + _COH + _CAP),
    # IPM-corrected coherent veto: v(q) = 1 - S_inc/L_LFG ({ipm} is filled by data_runs.py
    # from ipm_veto.py with the run's density file)
    "B_eik_Gauss_single_ipm": ("SudakovB", "  Sudakov:\n" + _EIK + "    Pauli: Step\n{ipm}" + _CAP),
    "B_k1_Box_single_ipm": ("SudakovB", "  Sudakov:\n    Law: Kappa1\n    Profile: Box\n    Pauli: Step\n{ipm}" + _CAP),
    "A_fold_finiteA_ipm": ("SudakovA", "  Sudakov:\n" + _AFF + "    Pauli: Step\n{ipm}" + _CAP),
}


def read_hist(path):
    rows = np.loadtxt(path, skiprows=2)
    return rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]


def counts(value, error):
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
            out.append({"p": 0.5 * (a + b), "T": 1 - frac, "err": err, "n": int(tot)})
        else:
            area = np.pi * BEAM_RADIUS**2 * 10
            out.append({"p": 0.5 * (a + b), "sigmaR": area * frac, "err": area * err, "n": int(tot)})
    return out


def run_one(build, outdir, tag, mode, alg, sud, inter, nevents):
    name = outdir / tag
    card = outdir / f"{tag}.yml"
    card.write_text(RUNCARD.format(seed=12345, mode=mode, interactions=INTERACTIONS[inter], alg=alg,
                                   sudakov=sud.replace("{cap}", str(CAP[inter])), radius=BEAM_RADIUS, pmin=PBINS[0], pmax=PBINS[-1],
                                   nevents=nevents, name=name))
    t0 = time.time()
    res = subprocess.run([str(build / "bin" / "achilles-cascade"), str(card)], cwd=build,
                         capture_output=True, text=True)
    for f in outdir.glob(f"{tag}.hepmc*"):
        f.unlink()
    if res.returncode != 0:
        return tag, {"error": res.stderr[-1500:]}
    return tag, {"seconds": time.time() - t0, "result": summarize(name, mode)}


def main():
    build = pathlib.Path(sys.argv[1]).resolve()
    outdir = pathlib.Path(sys.argv[2]).resolve()
    inter = sys.argv[3] if len(sys.argv) > 3 else "Geant"
    nevents = int(sys.argv[4]) if len(sys.argv) > 4 else 40000
    outdir.mkdir(parents=True, exist_ok=True)
    if INTERACTIONS["NN"] is None:
        lines = (build / "data" / "default" / "VirtResInteractions.yml").read_text().splitlines()
        INTERACTIONS["NN"] = "\n".join("    " + l for l in lines if l.strip())
    jobs = [(f"{inter}_{mode}_{key}", mode, alg, sud) for mode in ("CrossSection", "Transparency")
            for key, (alg, sud) in ALGOS.items()]
    results = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for tag, res in pool.map(lambda j: run_one(build, outdir, j[0], j[1], j[2], j[3], inter, nevents), jobs):
            results[tag] = res
            print(tag, json.dumps(res), flush=True)
    (outdir / f"compare_{inter}.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
