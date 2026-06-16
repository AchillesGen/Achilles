#!/usr/bin/env python3
"""Cross-section consistency validation across integrator backends.

The total cross section is a physical quantity: every integrator (Vegas-based
``MultiChannel`` or any normalizing-flow backend) is an *unbiased* Monte Carlo estimator
of the same integral, so they must all agree within their statistical errors -- regardless
of how well a flow is trained. This script runs the same physics (the default run card)
with each integrator and checks that:

  * every backend's cross section agrees with the C++ ``MultiChannel`` baseline within
    ``--nsigma`` combined statistical errors, and
  * the standalone C++ binary and the Python entry point give the same MultiChannel result
    ("consistent with the C++-only setup").

It is deliberately fast: the cross section is computed from the integrator's fixed-size
maximum-weight sampling pass, so a small ``--nevents`` suffices and flow training can be
short (an under-trained flow just has larger error bars, not a biased mean).

Usage (after ``pip install -e .[flow]`` from the repo root)::

    python validation/validate_integrators.py                 # all available backends
    python validation/validate_integrators.py --only uniform normflow
    python validation/validate_integrators.py --nevents 2000 --nsigma 5

Each backend runs in its own subprocess (isolating RNG/interpreter state); generated run
cards, logs and output land under ``validation/{cards,logs,out}/``.
"""

import argparse
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
VALID_DIR = REPO_ROOT / "validation"
CARDS_DIR = VALID_DIR / "cards"
LOGS_DIR = VALID_DIR / "logs"
OUT_DIR = VALID_DIR / "out"

# Physics copied verbatim from data/default/run.yml, with Options inlined so the integrator
# can be varied per card. Placeholders: {nevents}, {output}, {seed}, {integrator}.
CARD_TEMPLATE = """\
Main:
  NEvents: {nevents}
  HardCuts: true
  EventCuts: false
  DoRotate: false
  RunDecays: False
  Output:
      Format: NuHepMC
      Name: {output}
      Zipped: True

SherpaOptions: !include "data/default/SherpaOptions.yml"

Processes:
  - Leptons: [11, [11]]

Beams:
  - Beam:
      PID: 11
      Beam Params:
        Type: Monochromatic
        Energy: 1108

Cascade:
  Run: False
  Interactions: !include "data/default/VirtResInteractions.yml"
  Step: 0.04
  Probability: Cylinder
  InMedium: None
  PotentialProp: False
  Algorithm: Base

NuclearModels:
- NuclearModel:
   Model: QESpectral
   ConfigFile: data/info_C12_pke.data
   SpectralP: data/Spectral_Functions/pke12p_tot.data
   SpectralN: data/Spectral_Functions/pke12n_tot.data
   FormFactorFile: "FormFactors.yml"
   Ward: None

Nuclei:
  - Nucleus: !include "data/default/12C.yml"

HardCuts:
  - Type: AngleTheta
    PIDs: 11
    range: [36.5, 38.5]

Options:
  Initialize:
    Seed: {seed}
  Integrator:
    BatchSize: 1024
{integrator}
  Unweighting:
    Name: Percentile
    percentile: 99

Backend:
  Name: Default
  Options: []

# Each validation run computes from scratch (no cross-run cache reuse).
Cache:
  Load: false
  Save: false
"""

# Flow training profiles. "fast" just checks the unbiased mean cheaply; "quality" trains
# long enough to test whether weak agreement was merely under-training (the mean should
# tighten toward the baseline and the error shrink).
TRAIN_PROFILES = {
    "fast": {"Epochs": 15, "NIterations": 3, "NCalls": 2000},
    "quality": {"Epochs": 300, "NIterations": 20, "NCalls": 10000},
}

# name -> dict(kind=MultiChannel|Flow, params={...}, needs=[modules], runner=python|cpp).
# `params` for a Flow are the constructor kwargs (Class + architecture); the training
# profile (Epochs/NIterations/NCalls) is merged in per run.
CONFIGS = {
    "multichannel": {"kind": "MultiChannel", "params": {"Accuracy": "1.0e-2"},
                     "needs": [], "runner": "python"},
    "uniform": {"kind": "Flow", "params": {"Class": "UniformFlow"},
                "needs": [], "runner": "python"},
    "normflow": {"kind": "Flow", "params": {"Class": "NormFlow", "hidden": 32, "blocks": 4},
                 "needs": ["torch", "normflows"], "runner": "python"},
    "spline": {"kind": "Flow",
               "params": {"Class": "SplineFlow", "hidden": 32, "blocks": 4, "bins": 8},
               "needs": ["torch", "normflows"], "runner": "python"},
    "zuko_nsf": {"kind": "Flow",
                 "params": {"Class": "ZukoFlow", "flow": "NSF", "hidden": 32, "transforms": 3},
                 "needs": ["torch", "zuko"], "runner": "python"},
    "zuko_cnf": {"kind": "Flow", "params": {"Class": "ZukoFlow", "flow": "CNF", "hidden": 32},
                 "needs": ["torch", "zuko"], "runner": "python"},
    "flowmatching": {"kind": "Flow",
                     "params": {"Class": "FlowMatchingFlow", "hidden": 32, "step_size": 0.1},
                     "needs": ["torch", "flow_matching"], "runner": "python"},
}


def integrator_yaml(cfg, profile):
    """Render the Options/Integrator body for a config and training profile."""
    params = dict(cfg["params"])
    if cfg["kind"] == "Flow":
        params.update(profile)  # training profile only applies to flows
    lines = [f"    Type: {cfg['kind']}", f"    {cfg['kind']}:"]
    lines += [f"      {k}: {v}" for k, v in params.items()]
    return "\n".join(lines) + "\n"

XSEC_RE = re.compile(r"Total xsec:\s*([0-9eE.+-]+)\s*\+/-\s*([0-9eE.+-]+)")


def have_modules(mods):
    return all(importlib.util.find_spec(m) is not None for m in mods)


def write_card(name, cfg, nevents, seed, profile):
    CARDS_DIR.mkdir(parents=True, exist_ok=True)
    card = CARDS_DIR / f"run_{name}.yml"
    card.write_text(
        CARD_TEMPLATE.format(
            nevents=nevents,
            output=f"validation/out/{name}.hepmc",  # relative to the repo root (the run cwd)
            seed=seed,
            integrator=integrator_yaml(cfg, profile).rstrip("\n"),
        )
    )
    return card


def parse_xsec(logfile):
    """Sum the per-process-group 'Total xsec' lines (errors added in quadrature)."""
    mean, var = 0.0, 0.0
    found = False
    for m in XSEC_RE.finditer(Path(logfile).read_text(errors="replace")):
        found = True
        mean += float(m.group(1))
        var += float(m.group(2)) ** 2
    return (mean, var**0.5) if found else (None, None)


def run_python(card, logfile):
    code = (
        "import Achilles;"
        f"Achilles.initialize_logging(3, 2, {str(logfile)!r});"
        f"Achilles.generate_events({str(card)!r}, batch_mode=True)"
    )
    return subprocess.run(
        [sys.executable, "-c", code], cwd=REPO_ROOT, capture_output=True, text=True
    )


def run_cpp(card, logfile):
    binary = REPO_ROOT / "build" / "bin" / "achilles"
    if not binary.exists():
        return None
    return subprocess.run(
        [str(binary), str(card), "-b", "--logfile", str(logfile)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def run_one(name, cfg, nevents, seed, profile):
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    card = write_card(name, cfg, nevents, seed, profile)
    logfile = LOGS_DIR / f"{name}.log"
    if logfile.exists():
        logfile.unlink()

    runner = run_cpp if cfg.get("runner") == "cpp" else run_python
    proc = runner(card, logfile)
    if proc is None:
        return {"status": "SKIP", "reason": "binary not found"}
    if proc.returncode != 0:
        output = (proc.stderr or "") + (proc.stdout or "")
        # A missing optional flow dependency is not a validation failure: skip it.
        # (find_spec can miss a package whose submodules fail to import, so we also
        # catch the ImportError raised at run time here.)
        if "ImportError" in output or "ModuleNotFoundError" in output:
            return {"status": "SKIP", "reason": "optional dependency not installed"}
        tail = output.strip().splitlines()[-5:]
        return {"status": "ERROR", "reason": "\n        ".join(tail) or "nonzero exit"}

    mean, err = parse_xsec(logfile)
    if mean is None:
        return {"status": "ERROR", "reason": f"no 'Total xsec' in {logfile}"}
    return {"status": "OK", "mean": mean, "err": err}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nevents", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=123456)
    ap.add_argument("--nsigma", type=float, default=5.0)
    ap.add_argument("--only", nargs="*", help="restrict to these config names")
    ap.add_argument("--also-cpp-binary", action="store_true",
                    help="also run the standalone C++ binary for MultiChannel as an "
                         "extra C++-only cross-check")
    ap.add_argument("--quality", action="store_true",
                    help="train the flows much longer (the 'quality' profile) to check "
                         "whether weak agreement was just under-training")
    args = ap.parse_args()

    profile = TRAIN_PROFILES["quality" if args.quality else "fast"]
    print(f"Flow training profile: {'quality' if args.quality else 'fast'} {profile}")

    if importlib.util.find_spec("Achilles") is None:
        sys.exit("Achilles is not importable. Run `pip install -e .[flow]` from the repo root.")

    names = args.only or list(CONFIGS)
    # The MultiChannel baseline is always needed to judge consistency.
    if "multichannel" not in names:
        names = ["multichannel"] + list(names)
    if args.also_cpp_binary:
        CONFIGS["multichannel_cpp"] = dict(CONFIGS["multichannel"], runner="cpp")
        if not args.only:
            names.append("multichannel_cpp")

    results = {}
    for name in names:
        cfg = CONFIGS[name]
        if not have_modules(cfg["needs"]):
            results[name] = {"status": "SKIP",
                             "reason": f"missing {', '.join(cfg['needs'])}"}
            print(f"[skip] {name}: missing {', '.join(cfg['needs'])}")
            continue
        print(f"[run ] {name} ...", flush=True)
        results[name] = run_one(name, cfg, args.nevents, args.seed, profile)

    # Baseline: the C++ MultiChannel integrator (run through the Python entry point, i.e.
    # the same C++ engine). The optional standalone-binary result is checked against it too.
    base = results.get("multichannel")
    base_ok = base and base["status"] == "OK"

    print("\n=== Cross-section consistency ===")
    print(f"{'backend':<18}{'xsec':>14}{'error':>12}{'Δ/σ vs MC':>12}  status")
    overall_ok = True
    for name in results:
        r = results[name]
        if r["status"] != "OK":
            overall_ok = overall_ok and r["status"] == "SKIP"
            print(f"{name:<18}{'-':>14}{'-':>12}{'-':>12}  {r['status']}: {r.get('reason','')}")
            continue
        if base_ok and name != "multichannel":
            sigma = (r["err"] ** 2 + base["err"] ** 2) ** 0.5
            nsig = abs(r["mean"] - base["mean"]) / sigma if sigma > 0 else float("inf")
            ok = nsig <= args.nsigma
            overall_ok = overall_ok and ok
            verdict = "PASS" if ok else "FAIL"
            print(f"{name:<18}{r['mean']:>14.6g}{r['err']:>12.3g}{nsig:>12.2f}  {verdict}")
        else:
            print(f"{name:<18}{r['mean']:>14.6g}{r['err']:>12.3g}{'(baseline)':>12}  OK")

    if not base_ok:
        print("\nWARNING: MultiChannel baseline did not produce a cross section; "
              "cannot judge consistency.")
        overall_ok = False

    print(f"\nOverall: {'PASS' if overall_ok else 'FAIL'}")
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
