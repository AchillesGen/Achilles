# Achilles physics validation (physval)

Compares an Achilles branch against the stored `main` baseline and against experimental
data, using NUISANCE3 (legacy NUISANCE2 record) as the comparison engine, and posts the
result as a PR comment. The CI wiring is `.github/workflows/physval.yml`.

| File | Role |
|---|---|
| `physval.py` | Driver: generate → bin → compare → report; `--make-baseline`, `--merge` |
| `adapters.py` | Achilles + NUISANCE3 boundary (`Nuisance3Adapter`) and the synthetic stand-in for `--dry-run` |
| `stats.py` | Trial-count covariance, compatibility / goodness-of-fit χ², BH |
| `report.py` | `comment.md` and `summary.json` |
| `plots.py` | Data / main / branch overlay with a ratio panel, one PNG per measurement |
| `ci_scope.py` | Decides what a workflow run does (mode, setups, whether the baseline is current) |
| `test_physval.py` | The tests; needs neither Achilles nor NUISANCE |
| `measurements.yml` | Setups (one run card each) and their measurements; per-sample overrides |

## Running locally

```bash
pip install numpy scipy pyyaml matplotlib pytest
python3 -m pytest

# Synthetic end to end:
python3 physval.py --dry-run --make-baseline --key dev --out-dir /tmp/pv/baselines/dev
python3 physval.py --dry-run --key dev --baseline-dir /tmp/pv --out-dir /tmp/pv/out
```

The real path needs `achilles` on PATH and pyNUISANCE, i.e. the
`ghcr.io/achillesgen/achilles-physval` image.

## CI

| trigger | what runs |
|---|---|
| `!physval` in a commit subject on a branch | compare, whole suite |
| `!physval(A,B)` | compare, setups A and B only |
| `!physval(dry-run[,A])` | the same through the synthetic adapter |
| nightly | baseline on `main`, skipped if the stored one is current |
| dispatch | either mode, with `events` / `seed` / `dry_run` / `only_experiment` |

Only commit subjects are read, so a body that documents the syntax triggers nothing.
Markers are ignored on `main`, including those a merge brings in.
An unknown setup name fails the run with the list of valid ones.

Both modes build Achilles once in the physval image and run one job per setup, since
generation is the expensive step and every measurement in a setup reuses its events.
A setup whose generation fails is still reported, with the reason, and the run goes red.

### The `physval-artifacts` repo

Results live in [AchillesGen/physval-artifacts](https://github.com/AchillesGen/physval-artifacts),
so they never add to a clone of this repo. Its `main` branch is served by GitHub Pages:

```
baselines/<image-key>/<setup>.json   main's predictions + MC covariance, and the main sha
plots/<sha>/*.png                    overlays linked from the comment
plots/<sha>/.published               publish time, which sets the plots' expiry
```

The workflow pushes with the deploy key in the `PHYSVAL_ARTIFACTS_KEY` secret. Each publish
replaces the repo's single commit, so it keeps no history, and leaves out plots older than
60 days. Expiry therefore happens at the next publish, not on the day.

The key is the image's manifest digest, so rebuilding the toolchain starts a fresh
baseline and older keys are dropped. The nightly run regenerates every setup's baseline
when the image changed, a setup has none, or `main` touched physics paths
(`ci_scope.PHYSICS_PATHS`, which includes `measurements.yml`) since the
stored sha. Other changes under `validation/physval` don't count. If one moves predictions,
run the workflow by hand on `main` with mode `baseline`, which always regenerates.
A baseline is only used if it came from the same adapter and has the
same binning; otherwise that measurement is compared with data only and marked ⚪.
Dry-run results are never published.

## Reading the comment

* **p_cmp**: main vs branch, correlated χ² with both MC covariances. **q**: the same after
  a Benjamini-Hochberg correction over the suite; a row is flagged when q < 0.05.
* **Δχ²**: the change in agreement with data. It only labels a flag: 🚩 worse, ⭐ better,
  🔀 |Δχ²| < 1.
* **p_data**: the branch's goodness of fit to data, for context.
* `summary.json` also records `selected_events` per measurement. Samples that share a
  selection should agree, so an outlier points at a selection bug.

## Statistics

The MC covariance comes from the generator's trial counters, which NuHepMC writes in
every event's `GenCrossSection`:

    C_jk = δ_jk Σ_{i∈k} w_i² − S_j S_k / N_nonzero + (δσ/σ)² S_j S_k

The multinomial is over all N non-zero trials, so the selection's acceptance is free to
fluctuate. Resampling only the selected events holds their count fixed and under-covers
the normalisation (~3× at 19% acceptance). This error once made two runs of one
configuration look incompatible at p ~ 1e-150. `test_physval.py` checks that the
p-values of two such runs are uniform.

## Adapter notes

* The adapter divides the per-atom `fatx` by each setup's `xsec_divisor` itself.
  NUISANCE's PerNucleon divides CH by 12, not 13.
* The legacy record never calls a sample's `ConvertEventRates`, so the adapter calls it
  on unit vectors through the `convert_event_rates` patch the image carries, and applies
  the resulting matrix as R S, R C Rᵀ. This covers bin widths and any smearing. The
  exceptions are in `measurements.yml`: `solid_angle` for per-steradian electron data,
  whose own conversion has no normalisation, and `data_scale` for a table shipped in the
  wrong units.
* A measurement that selects no events gets no χ² or plot. It is reported as "no events
  selected" and fails the run, unless `may_be_empty` in `measurements.yml` gives a reason
  it could happen by chance, which makes it a warning. A baseline leaves such a sample out.
* `MicroBooNE_NCpi0_*` and `ElectronData_*` need nuisance#115 and #116, which the image's
  pinned sha includes, and `MicroBooNE_CC1Mu0pNp_*` needs a patch the image carries.
  Without them they come back empty rather than red, so check `versions.json`'s
  `nuisance2`.
* Left out on purpose: `*_XSec_1DEnu*` and `MicroBooNE_CC1Mu0pNp_XSec_ENu_nu`
  (flux-integrated scaling not implemented),
  MiniBooNE CCQE (per-neutron convention), bubble-chamber deuterium, Fe/Pb ratios,
  coherent, DIS, and the 3D `multidif` samples.
* The image's `LD_LIBRARY_PATH` puts NUISANCE2's spdlog/fmt ahead of Achilles' own, and
  mixing them segfaults silently. The workflow puts Achilles' lib first; do the same when
  running locally.
