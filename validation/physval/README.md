# Achilles physics validation (physval)

Compares an Achilles branch against the stored `main` baseline and against experimental
data, using NUISANCE3 (legacy NUISANCE2 record) as the comparison engine, and posts the
result as a PR comment. The CI wiring is `.github/workflows/physval.yml`.

| File | Role |
|---|---|
| `physval.py` | Driver: generate → bin → compare → report; `--make-baseline`, `--merge`, `--unweighter-scan` |
| `adapters.py` | Achilles + NUISANCE3 boundary (`Nuisance3Adapter`) and the synthetic stand-in for `--dry-run` |
| `stats.py` | Trial-count covariance, compatibility / goodness-of-fit χ², BH and Bonferroni |
| `report.py` | `comment.md` and `summary.json` |
| `plots.py` | Data / main / branch overlay with a ratio panel, one PNG per measurement |
| `ci_scope.py` | Decides what a workflow run does (mode, setups, whether the baseline is current) |
| `measurements.yml` | Setups (one run card each) and their measurements; per-sample overrides; scan variants |
| `smearing/` | Wiener-SVD A_C matrices |

## Running locally

```bash
pip install numpy scipy pyyaml matplotlib
for f in stats report adapters plots physval ci_scope; do python3 $f.py --selftest; done

# Synthetic end to end:
python3 physval.py --dry-run --make-baseline --key dev --out-dir /tmp/pv/baselines/dev
python3 physval.py --dry-run --key dev --baseline-dir /tmp/pv --out-dir /tmp/pv/out
```

The real path needs `achilles` on PATH and pyNUISANCE, i.e. the
`ghcr.io/achillesgen/achilles-physval` image.

## CI

| trigger | what runs |
|---|---|
| `!physval` in a commit subject | compare, whole suite |
| `!physval(A,B)` | compare, setups A and B only |
| `!physval(dry-run[,A])` | the same through the synthetic adapter |
| nightly | baseline on `main`, skipped if the stored one is current |
| dispatch | either mode, with `events` / `seed` / `dry_run` / `only_experiment` |

Only commit subjects are read, so a body that documents the syntax triggers nothing.
An unknown setup name fails the run with the list of valid ones.

Both modes build Achilles once in the physval image and run one job per setup, since
generation is the expensive step and every measurement in a setup reuses its events.
A setup whose generation fails is still reported, with the reason, and the run goes red.

### The `physval-baselines` branch

```
baselines/<image-key>/<setup>.json   main's predictions + MC covariance, and the main sha
plots/<sha>/*.png                    overlays linked from the comment
```

The key is the image's manifest digest, so rebuilding the toolchain starts a fresh
baseline and older keys are dropped. The nightly run regenerates a setup's baseline only
when the image changed or `main` touched physics paths (`ci_scope.PHYSICS_PATHS`) since
the stored sha. A baseline is only used if it came from the same adapter and has the
same binning; otherwise that measurement is compared with data only and marked ⚪.
`physval-prune.yml` drops plot folders older than 45 days. Dry-run results are never
published.

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
configuration look incompatible at p ~ 1e-150. `stats.py --selftest` checks that the
p-values of two such runs are uniform.

## Adapter notes

* The adapter divides the per-atom `fatx` by `target_nucleons` itself. NUISANCE's
  PerNucleon divides CH by 12, not 13.
* The legacy record never calls a sample's `ConvertEventRates`, so `measurements.yml`
  supplies what it would have done: `bin_edges` / `bin_widths` for bin-number
  histograms, `solid_angle` for per-steradian electron data, `smearing` for Wiener-SVD
  samples (applied to the prediction as R S, R C Rᵀ), and `data_scale` for a table
  shipped in the wrong units.
* `MicroBooNE_NCpi0_*`, `MicroBooNE_CC1Mu0pNp_*` and `ElectronData_*` need NUISANCE2
  patches (nuisance#115, #116, and the CC1Mu0pNp one) that the image carries; without
  them they come back empty rather than red. Check `versions.json`'s
  `nuisance2.patches`.
* Left out on purpose: `*_XSec_1DEnu*` (flux-integrated scaling not implemented),
  MiniBooNE CCQE (per-neutron convention), bubble-chamber deuterium, Fe/Pb ratios,
  coherent, DIS, and the 3D `multidif` samples.
* The image's `LD_LIBRARY_PATH` puts NUISANCE2's spdlog/fmt ahead of Achilles' own, and
  mixing them segfaults silently. The adapter and the build job put Achilles' lib first.

## Unweighting scan

`--unweighter-scan` generates each setup once per `unweighting.variants` entry, with the
same seed and the same accepted-event count, and compares every variant with
`unweighting.reference`. Unweighting is variance reduction, so a flag means a scheme
changed the physics. The table also reports what each scheme costs and buys: ESS per
event, MC error, Δnorm, and wall time. `null-control` is the reference reseeded; the
flagging thresholds are floored at its score. `adapters.unweighting_cap` mirrors the C++
cap rules, so `--dry-run` exercises them too.

```bash
python3 physval.py --unweighter-scan --only-experiment MINERvA_numu_CH --events 200000 --out-dir out
python3 physval.py --unweighter-scan --merge out/*/summary.json --out-dir final
```
