# Achilles distribution-based physics validation (physval)

A physics-validation pipeline that compares an Achilles feature branch against
**experimental data** and against the **`main` baseline**, using NUISANCE3 (with
the legacy NUISANCE interface enabled) as the comparison engine. It reports χ² and
p-values per measurement and posts a Sherpa-style summary table as a PR comment,
flagging any measurement whose `main`-vs-feature **compatibility p-value drops below
0.05**.

This directory holds the driver and its statistics; the CI wiring lives in
`.github/workflows/physval*.yml`. It is intended to be relocated into the new
Achilles Python package once that lands — nothing here imports Achilles Python, so
the move is a path change.

## Layout

| File | Role |
|------|------|
| `stats.py` | Bootstrap covariance, correlated-χ² compatibility (`p_compat`), goodness-of-fit vs data, Bonferroni overall p. Hartlap-corrected. **No external deps beyond numpy/scipy.** |
| `report.py` | Renders the Sherpa-style `comment.md` table + `summary.json`. |
| `adapters.py` | The boundary to Achilles/NUISANCE3: `generate` (once per experimental setup) + `histogram` (per measurement). `Nuisance3Adapter` is the real path; `SyntheticAdapter` backs `--dry-run` and the self-tests. |
| `plots.py` | Publication-style data/main/branch overlay + ratio panel, one PNG per measurement. |
| `physval.py` | Driver: config → generate → stats → report; plus `--make-baseline`. |
| `measurements.yml` | Experiments (run cards) each grouping the measurements that reuse their events. |

## Run it locally

```bash
pip install numpy scipy pyyaml

# Unit self-tests (no external tools):
python3 stats.py  --selftest      # bootstrap calibration + power + Bonferroni
python3 report.py --selftest      # table rendering / flagging / sorting
python3 physval.py --selftest     # end-to-end w/ an injected regression

# Full synthetic dry-run from the config:
python3 physval.py --config measurements.yml --dry-run --make-baseline \
    --key dev --baseline-dir /tmp/pv
python3 physval.py --config measurements.yml --dry-run \
    --key dev --baseline-dir /tmp/pv --out-dir /tmp/pv/out --feature-sha $(git rev-parse HEAD)
# -> /tmp/pv/out/comment.md  and  /tmp/pv/out/summary.json
```

## The statistics

- **Bootstrap covariance** (`bootstrap_covariance`): the feature (and stored `main`)
  prediction's MC uncertainty is estimated by resampling its **weighted events with
  replacement**. No re-generation, works with weighted / negative-weight events.
- **Compatibility** (`compatibility`, drives the flag): a correlated χ²
  `Δᵀ (C_main + C_feature)⁻¹ Δ`, `Δ = h_feature − h_main`; `p_compat` from the χ²
  survival function. **Flag when `p_compat < 0.05`.**
- **Hartlap correction** (`hartlap_factor`): the inverse of a bootstrap covariance
  is biased high, inflating χ². The Hartlap factor `(N−p−2)/(N−1)` debiases it.
  **Requirement: `n_boot ≫ n_bins`** (need `n_boot > n_bins + 2` at minimum). The
  self-tests confirm the false-flag rate sits at ~0.05 once this holds.
- **Goodness-of-fit** (`goodness_of_fit`): χ² of each prediction vs data, reported
  per row for context; does **not** drive the flag.
- **Bonferroni** (`bonferroni`): overall family-wise p `min(1, N·min_i p_i)` in the
  comment header; overall concern when `< 0.05`. Standard particle-physics
  multiple-testing control.

Non-Gaussian fallback: `empirical_pvalue` reads `p_compat` straight off a bootstrap
Δχ² null distribution if a measurement's statistic is far from χ²-distributed.

## Adapter boundary

The CI runs inside the public `ghcr.io/achillesgen/achilles-physval` container
(NUISANCE3 + ROOT + HepMC3 + ProSelecta + Achilles build deps), builds the current
Achilles checkout there, and drives the adapter with `achilles` and NUISANCE3 on
PATH. The adapter has two stages so the expensive step runs once per setup:

* `generate(experiment, …)` renders the run card with the pinned seed / event count /
  output path (`!include`s are round-tripped, and the `Options` include is expanded so
  `Initialize.Seed` can be set), runs `achilles <card>` with `cwd` = repo root, and
  returns the NuHepMC file. It is deleted once the setup's last measurement is binned.
* `histogram(generated, measurement)` opens the file with `pn.EventSource`, lets
  `IAnalysis.add_to_framegen` register the sample's selection + projections, and walks
  the frame to get each selected event's bin (`Binning.find_bin`) and weight.

Normalisation is **not** reimplemented here: `IAnalysis.process` yields the
cross-section-scaled, bin-width-divided prediction, and the per-event weights are
rescaled so they sum to it. The bootstrap therefore measures MC uncertainty directly
in the data's units. `data_table` takes the published values and the *full*
`get_covariance_matrix()`, falling back to diagonal errors only if that is absent.

Note the legacy NUISANCE2 record resolves analyses lazily — `get_analyses()` returns
an empty list, but `record.analysis("<sample>")` works.

### Achilles' libraries must precede the image's

The image exports `LD_LIBRARY_PATH=/opt/root/lib:/opt/nuisance3/lib:/opt/nuisance2/lib:…`,
which outranks the binary's `RUNPATH`. NUISANCE2 ships `libspdlog` built against
`fmt v10` while Achilles bundles `fmt v11`; loading both leaves Achilles calling a
spdlog with a mismatched fmt and it **segfaults inside `InitializeLogging`, printing
nothing at all** (the splash is lost to stdout buffering, so it looks like an instant
silent crash). Putting Achilles' own `lib` first fixes it. `Nuisance3Adapter` does this
itself when it spawns achilles, and the CI build job smoke-tests `achilles --version`
so a regression fails early instead of mid-generation.

## `physval-baselines` branch

Durable store for the `main` baseline and the run plots (created on demand):

```
baselines/<key>.json        # main central histograms + bootstrap covariance + data
plots/<feature-sha>/*.png    # per-run new/old/data overlays (inline in the comment)
```

`<key>` is the short physval-image manifest digest, so any toolchain rebuild
(NUISANCE2/3, ROOT, HepMC3, ProSelecta) invalidates the baseline. The baseline JSON
is tiny and kept indefinitely; the plots are pruned on a time window by
`physval-prune.yml`. Where a run has no matching baseline, the driver recomputes
`main` inline and notes it in the comment header.

## CI workflows

- `physval.yml` — triggered by a `!physval` commit-message marker,
  `workflow_dispatch`, or nightly `schedule`. Builds Achilles once inside the
  container, shards event generation across a matrix (one job per **experimental
  setup**; measurements in a setup reuse its events), aggregates, commits plots to
  `physval-baselines/plots/<sha>/`, and posts/updates the PR comment.
- `physval-baseline-refresh.yml` — on push to `main`, recompute and store the
  baseline (predictions + bootstrap covariance + data) on `physval-baselines`,
  building Achilles in the same container.
- `physval-prune.yml` — scheduled pruning of old `plots/<sha>/` folders.

See the plan for the full rationale (bootstrap vs same-seed, merge-base attribution,
public-runner compute budget).
