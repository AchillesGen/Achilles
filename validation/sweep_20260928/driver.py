"""Run one flow config through Achilles and report 1% TailFraction efficiency + timings.

usage: driver.py TAG MODE EPOCHS NCALLS 'json-of-flow-params' [NEVENTS]
"""
import json, re, sys, os
from pathlib import Path
import numpy as np
REPO = Path("/home/isaacson/Documents/Projects/NeutrinoGenerator/Achilles/138-feature-pip-installable")
sys.path.insert(0, str(REPO / "validation"))
import validate_integrators as v

tag, mode, epochs, ncalls, params = sys.argv[1:6]
nevents = int(sys.argv[6]) if len(sys.argv) > 6 else 100000
params = json.loads(params)
kind = params.pop("Kind", "Flow")  # "MultiChannel" for the VEGAS multichannel baseline
params.setdefault("Accuracy", "1.0e-4")
if kind == "Flow":
    params["WeightsPath"] = f"validation/out/{tag}.weights"
cfg = {"kind": kind, "params": params, "needs": [], "runner": "python"}
profile = {"Epochs": int(epochs), "NIterations": 3, "NCalls": int(ncalls)}
res = v.run_one(tag, cfg, nevents, 123456, profile, "weighted", mode, "monochromatic")
out = {"tag": tag, "mode": mode, **{k: res.get(k) for k in ("status", "mean", "err", "reason")}}
if res["status"] == "OK":
    import pyhepmc
    hep = REPO / "validation" / "out" / f"{tag}.hepmc"
    w = []
    with pyhepmc.open(str(hep)) as f:
        for e in f: w.append(e.weight("CV"))
    if not os.environ.get("KEEP"): hep.unlink()
    w = np.abs(np.array(w))
    s = np.sort(w); total = s.sum()
    suffix = np.cumsum(s[::-1])[::-1]
    i = int(np.argmax(suffix <= 0.01 * total)) if np.any(suffix <= 0.01 * total) else len(s)
    cap = s[i - 1] if i > 0 else s[0]   # TailFractionUnweighter::ComputeCap
    out["eff1"] = float(np.minimum(w, cap).sum() / cap / len(w))
    out["relstd"] = float(w.std() / w.mean())
    log = (REPO / "validation" / "logs" / f"{tag}.log").read_text()
    def secs(name):
        m = re.search(rf"{name}\s+\d+\s+((?:\d+h )?(?:\d+m )?[\d.]+ ?m?s)", log)
        if not m: return None
        t, tot = m.group(1), 0.0
        for val, unit in re.findall(r"([\d.]+) ?(h|ms|m|s)", t):
            tot += float(val) * {"h": 3600, "m": 60, "s": 1, "ms": 1e-3}[unit]
        return tot
    out["t_train"], out["t_scan"], out["t_gen"] = secs("Integrator Optimization") or 0.0, secs("Max Weight Scan") or 0.0, secs("Event Generation")
    errs = [float(x) for x in re.findall(r"^\[info\] +\d+ +[0-9.e+-]+ \+/- +[0-9.e+-]+ nb \( *([0-9.]+) +%\)", log, re.M)]
    out["iter_err_med_last20"] = float(np.median(errs[-20:])) if errs else None
    out["iters"] = len(errs)
    out["us_per_trial"] = out["t_gen"] / len(w) * 1e6
    out["hours_100M"] = (out["t_train"] + out["t_scan"] + 1e8 / out["eff1"] * out["t_gen"] / len(w)) / 3600
print("RESULT " + json.dumps(out), flush=True)
