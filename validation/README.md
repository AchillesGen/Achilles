# Integrator cross-section validation

End-to-end check that every integrator backend (Vegas `MultiChannel` and the
normalizing-flow backends) produces a **consistent total cross section** for the same
physics, and that the flow path agrees with the C++-only `MultiChannel` setup.

Why this works: the cross section is a physical quantity, and every integrator is an
*unbiased* Monte Carlo estimator of the same integral. So they must agree within their
statistical errors no matter how well a flow is trained — an under-trained flow just has
larger error bars, not a biased mean. (The reported error comes from the integrator's
fixed maximum-weight sampling pass, so a small `--nevents` is enough and flow training can
be short.)

## Run

From the repository root, after building the Python package with the flow extra:

```bash
pip install -e .[flow]
python validation/validate_integrators.py
```

Useful flags:

```bash
python validation/validate_integrators.py --only multichannel uniform normflow
python validation/validate_integrators.py --nevents 2000 --nsigma 5
python validation/validate_integrators.py --also-cpp-binary   # extra C++-only cross-check
```

Backends whose Python dependencies are missing are skipped (e.g. `zuko_cnf` without
`zuko`). `uniform` needs no ML dependency, so it always exercises the full
C++↔Python flow pipeline.

## What it does

- `cards/run_<name>.yml` — one run card per backend, generated from the default run card
  with only the `Options/Integrator` block changed (regenerated on each run).
- Runs each backend in its own subprocess; cross sections are parsed from
  `logs/<name>.log`; events go to `out/<name>.hepmc`.
- Compares every backend to the `multichannel` baseline within `--nsigma` combined errors
  and prints a table. Exit code is non-zero if any backend is inconsistent.
