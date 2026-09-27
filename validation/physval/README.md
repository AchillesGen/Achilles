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
| `adapters.py` | The boundary to Achilles/NUISANCE3: `generate` (once per setup) + `histogram_many` (all of a setup's measurements in one pass). `Nuisance3Adapter` is the real path; `SyntheticAdapter` backs `--dry-run` and the self-tests. |
| `plots.py` | Publication-style data/main/branch overlay + ratio panel, one PNG per measurement. |
| `physval.py` | Driver: config → generate → stats → report; plus `--make-baseline` and `--unweighter-scan`. |
| `measurements.yml` | Experiments (run cards) each grouping the measurements that reuse their events, the per-sample override maps (`data_scale`, `bin_edges`, `solid_angle`, `smearing`), plus the `unweighting:` variant list. |
| `smearing/` | Regularisation matrices (A_C) for the Wiener-SVD unfolded samples, as csv. |

## Run it locally

```bash
pip install numpy scipy pyyaml

# Unit self-tests (no external tools):
python3 stats.py    --selftest    # bootstrap calibration + power + Bonferroni
python3 report.py   --selftest    # comment rendering / grouping / flagging / sorting
python3 adapters.py --selftest    # the batched frame walk, against a fake pyNUISANCE
python3 physval.py  --selftest    # end-to-end w/ an injected regression

# Full synthetic dry-run from the config:
python3 physval.py --config measurements.yml --dry-run --make-baseline \
    --key dev --baseline-dir /tmp/pv
python3 physval.py --config measurements.yml --dry-run \
    --key dev --baseline-dir /tmp/pv --out-dir /tmp/pv/out --feature-sha $(git rev-parse HEAD)
# -> /tmp/pv/out/comment.md  and  /tmp/pv/out/summary.json
```

## Unweighting scan (`--unweighter-scan`)

A second comparison axis over the same machinery: instead of *this branch vs `main`*,
it runs *one unweighting scheme vs another* on a single checkout. Achilles registers
several `Options/Unweighting` schemes — `None` (keep the weights), `Percentile`
(cap at the p-th percentile of |w|), `Excess` (smallest cap `C` with
`Σ max(|w|−C,0) ≤ ε·Σ|w|`) and `TailFraction` (smallest `C` with
`Σ_{|w|>C} |w| ≤ ε·Σ|w|`) — and all of them are supposed to be *variance reduction*,
not physics: an event above the cap is kept carrying its excess weight, so no scheme
may move a distribution.

```bash
python3 physval.py --config measurements.yml --unweighter-scan \
    --only-experiment MINERvA_CC0pinp_STV_XSec --events 200000 --out-dir out
# add --dry-run to exercise the plumbing without Achilles/NUISANCE
# add --only-variant <name> (repeatable) to shorten the scan; the reference is kept

# One job per setup, then merge — same shard axis as the branch comparison:
python3 physval.py --unweighter-scan --merge out/*/summary.json --out-dir final
```

The variants live under `unweighting:` in `measurements.yml`. Each setup is generated
once per variant, at the **same seed and the same accepted-event count**, so the only
difference between two runs is the scheme. Every variant is then compared to
`unweighting.reference` — the loosest cap in the list, which leaves the fewest events
overweighted — with the same correlated χ² the branch comparison uses.

`Name: None` is **not** usable as that reference today: its events arrive normalised
low by a large factor, because `NoUnweighter::AcceptEvent` returns the raw weight
while every capped scheme returns `weight/cap`, and `EventGen::GenerateSingleEvent`
multiplies by the summed caps regardless. Fixing that would make `None` the natural
reference, since it is unbiased by construction.

### The null control is not optional

The variant list carries a `null-control`: the reference's own options at
`seed_offset: 1`, i.e. two runs of one configuration that differ only by seed. Its
rows are draws from the null hypothesis and must be uniform. The flagging thresholds
are floored at whatever it scores (`min(alpha, null)`, so a well-behaved control
leaves α alone and never tightens it), and with a correct covariance that floor is
inert — it is a guard, not a calibration.

It is worth knowing what this control caught, because the same trap is easy to walk
back into. `p_compat` used to come out at ~1e-150 for the control — two runs of the
identical configuration, judged incompatible — while every individual bin agreed
within its error (empirical run-to-run scatter / bootstrap σ = 0.9 per bin). The cause
was the resampling, not the physics: `bootstrap_covariance` drew a fixed number of the
**selected** events, which holds their count exact and so denies the acceptance
fluctuation. On `MINERvA_CC0pinp_STV_XSec` (19% acceptance) it claimed the
normalisation was known to 0.11% where five independent runs actually scatter by
0.41%, and a correlated χ² divides that real difference by a near-null eigenvalue.
Measured over 5 runs x 8 measurements = 80 seed pairs: **67/80 flagged before, 5/80
after** (median p 0.57). `e12C_1108` never showed it because its 72% acceptance leaves
the selected count nearly deterministic, which is the regime where fixed-count
resampling is accidentally right.

Two things follow. `stats.py --selftest` carries a `[trial-null]`/`[trial-cover]`
regression test that simulates runs the way a generator delivers them and asserts the
flag rate sits at α and that the reported normalisation uncertainty matches the runs'
actual scatter — the second is what the bootstrap failed. And the **branch
comparison** is covered by the same fix: it also compares two independent runs
(`main` at `seed`, feature at `seed + 1`) whenever there is no stored baseline.

`p (shape only)` — the same χ² with the normalisation fitted out and one dof given up
for it — is now a diagnostic rather than a crutch: a scheme that moved a distribution
fails both columns, two runs that merely disagree on the total cross section fail only
the first. The scale is fitted in the covariance metric, not as a ratio of plain bin
sums, because the weakly-constrained direction is the bin-content one and the two
differ once bin widths are unequal (that mismatch is why the old control reported a
broken *shape* for `dphit`, width CV 1.4, while `thmu`, width CV 0.2, looked fine).

A stored baseline records which estimator built its covariances; a baseline from a
different one is treated as stale and `main` is recomputed inline, since mixing the
two would compare a prediction against an uncertainty that was never meant for it.

Reading the summary table:

- **p (vs reference)** — Bonferroni over the setup's measurements, judged against the
  null control. A flag here means the scheme *changed the physics*, which is a bug
  rather than a trade-off.
- **p (shape only)** — as above with the overall normalisation fitted out.
- **ESS/event** — Kish effective sample size `(Σ|w|)²/(N Σw²)` on the raw generator
  weights, i.e. the statistical power the scheme delivers per accepted event
  (`1.0` = perfect unit weights). Computed pre-normalisation, so the bin-width
  division does not masquerade as weight spread.
- **MC error** — mean MC σ relative to the reference's. At a fixed accepted-
  event count, a harsher cap buys precision here and pays for it in **wall**.
- **Δnorm** — change in the integrated cross section, which no scheme should move.
- **max/mean w** — the heaviest surviving overweight, the tail the cap left behind.

`adapters.unweighting_cap` reimplements the three cap rules in numpy so `--dry-run`
is a real test of the scan (and a second opinion on the C++).

## The statistics

- **Trial-count covariance** (`trial_covariance`, what the table is built from): the
  MC uncertainty is written down analytically from the generator's own counters, which
  NuHepMC carries in each event's `GenCrossSection` (total xsec, its uncertainty,
  non-zero trials, total trials — the last event's copy has seen the whole file):

      C_jk = δ_jk Σ_{i∈k} w_i²  −  S_j S_k / N_nonzero  +  (δσ/σ)² S_j S_k

  The first two terms are the multinomial structure over **all** `N_nonzero` events
  the generator produced, so a selection's acceptance is free to fluctuate; the third
  is the cross-section uncertainty, which scales every bin together. Exact for
  weighted and negative-weight events, exact under a response matrix (`R C Rᵀ`),
  deterministic, and free.
- **Bootstrap covariance** (`bootstrap_covariance`): the same quantity by resampling
  the weighted events with replacement. Kept as a cross-check and as the fallback for
  events that carry no trial counters. It must **not** be used on a selected
  subsample: resampling a fixed number of selected events asserts that the number
  passing the selection is exact, which under-covers the normalisation by ~3x at a
  19% acceptance.
- **Compatibility** (`compatibility`, drives the flag): a correlated χ²
  `Δᵀ (C_main + C_feature)⁻¹ Δ`, `Δ = h_feature − h_main`; `p_compat` from the χ²
  survival function. **Flag when `p_compat < 0.05`.**
- **Hartlap correction** (`hartlap_factor`): applies only to the bootstrap fallback,
  whose inverse covariance is biased high (`(N−p−2)/(N−1)` debiases it, requiring
  `n_boot ≫ n_bins`). `trial_covariance` is not a finite-sample estimate, carries
  `n_boot = None`, and is left uncorrected — so `N_BOOT` no longer constrains how many
  bins a measurement may have.
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
* `histogram_many(generated, measurements)` opens the file once with `pn.EventSource`
  and puts **every** sample of the setup on one `EventFrameGen` — one `add_int_column`
  for each selection, one `add_double_column` per projection — then walks the frame in
  blocks, reading each sample's own columns to get its selected events' bins
  (`Binning.find_bin`) and weights. This is the notebook's "lots of projections"
  pattern (`nuisance3/notebooks/nuisance2.ipynb`), and it is why a setup with eighteen
  measurements costs one pass over its event file instead of eighteen: the file is
  opened, parsed and walked once. Blocks (250k events) keep peak memory flat in the
  run length, and the `fatx_per_sumw` estimate is taken from the last block, where it
  has seen the whole file. `histogram(generated, measurement)` is a one-measurement
  wrapper over the same code.

Normalisation is **not** reimplemented here: `IAnalysis.process` yields the
cross-section-scaled, bin-width-divided prediction, and the per-event weights are
rescaled so they sum to it. The bootstrap therefore measures MC uncertainty directly
in the data's units. `data_table` takes the published values and the *full*
`get_covariance_matrix()`, falling back to diagonal errors only if that is absent.

Note the legacy NUISANCE2 record resolves analyses lazily — `get_analyses()` returns
an empty list, but `record.analysis("<sample>")` works.

### What a sample needs beyond the record

Four things the legacy record does not do for us are declared per sample in
`measurements.yml` and applied in the adapter, because the record calls neither the
sample's `ConvertEventRates` nor its normalisation:

* `bin_edges` — the sample's histogram is indexed by bin *number* (the WireCell
  analyses stack several blocks), so the widths from the NUISANCE binning are 1 and the
  prediction has to be divided by the real ones.
* `bin_widths` — the same thing for an axis that stacks channels (a 0p block followed by
  an Np block). There is no monotonic edge list to give, but the per-bin widths are
  still defined, so they are listed directly and the plot runs against bin number.
* `solid_angle` — the Durham electron data is published per steradian while the
  selection integrates over its ±4° acceptance window.
* `smearing` — a Wiener-SVD unfolded measurement is only comparable to `A_C ×
  prediction`. The matrix is applied inside `bootstrap_covariance`, per replica, so it
  reaches the MC covariance as `A C Aᵀ` and not just the central values. It is never
  applied to the data.
* `data_scale` — one shipped table (MiniBooNE dσ/dQ²) is written 1e6 low; the factor
  multiplies the data and its covariance, not the prediction.

### Samples that need patched NUISANCE2

`MicroBooNE_NCpi0_*`, `MicroBooNE_CC1Mu0pNp_*` and `ElectronData_*` do not work against
NUISANCE2 as shipped: the two WireCell families set their bin index only in
`FillHistograms`, which the legacy record never calls, and the Durham samples cut on
variables the record fills *after* `isSignal`. All of them fail silently — an empty
prediction, not an error. The fixes are upstream as NUISANCEMC/nuisance#115, #116 and
the CC1Mu0pNp one, and are carried as patches in the physval image on top of its pinned
NUISANCE2 sha; `versions.json`'s `nuisance2.patches` records them. The δp_n unit mixing
that used to make `MicroBooNE_CC1Mu1p_XSec_1DDeltaPn_nu` incomparable to data is fixed
upstream (#114) and needs NUISANCE2 ≥ `86c64b44` in the image.

### Samples deliberately left out

* anything named `*_XSec_1DEnu*` — NUISANCE scales those **flux-integrated, not
  flux-averaged** (`fIsEnu1D`, `Measurement1D.cxx`), which the adapter does not
  implement, so they would come out silently mis-normalised;
* MiniBooNE `CCQE`/`CCQELike` — the published denominator is per *neutron* via a
  14.08/6.0 factor in the sample's own scaling, which is a convention call rather than
  something to guess, and getting it wrong is a factor of six;
* the bubble-chamber sets (ANL, BNL, BEBC, FNAL, GGM) — they need a deuteron initial
  state;
* Fe/Pb nuclear-target ratios (no spectral function), coherent pion and DIS/inclusive
  samples (not in the model), and the 3D/`multidif` samples, whose global-bin-number
  binning the record's own docs call out as snowflakes.

### Achilles' libraries must precede the image's

The image exports `LD_LIBRARY_PATH=/opt/root/lib:/opt/nuisance3/lib:/opt/nuisance2/lib:…`,
which outranks the binary's `RUNPATH`. NUISANCE2 ships `libspdlog` built against
`fmt v10` while Achilles bundles `fmt v11`; loading both leaves Achilles calling a
spdlog with a mismatched fmt and it **segfaults inside `InitializeLogging`, printing
nothing at all** (the splash is lost to stdout buffering, so it looks like an instant
silent crash). Putting Achilles' own `lib` first fixes it. `Nuisance3Adapter` does this
itself when it spawns achilles, and the CI build job smoke-tests `achilles --version`
so a regression fails early instead of mid-generation.

## Triggering a run

| trigger | effect |
|---|---|
| `!physval` in a pushed commit message | the whole suite, real NUISANCE3 path |
| `!physval(dry-run)` | the whole suite through the synthetic adapter |
| `!physval(<setup>[,<setup>...])` | just those setups, e.g. `!physval(MiniBooNE_CC1pi)` or `!physval(MiniBooNE_CC1pi,T2K_CC)` — one generation each instead of twenty |
| `!physval(dry-run,<setup>)` | both; order does not matter, and duplicates collapse |
| `workflow_dispatch` | same, via the `events` / `seed` / `dry_run` / `only_experiment` inputs (`only_experiment` takes a comma-separated list) |
| nightly `schedule` | the whole suite, for real |

Markers are read from commit **subject lines only**. A body that documents the syntax —
this repo's own history does — neither scopes a run nor triggers one; if a push mentions
`!physval` only in a body, the setup job says so and the rest of the workflow is skipped.
Scopes from every marker in the push are unioned, and shards run in config order
whatever order they were typed. Any name that is not a setup fails the `setup` job with
the list of valid names, rather than quietly running all twenty. The scope reaches the
aggregate too, so a scoped run expects only the setups it asked for and does not report
the rest as missing.

## The PR comment

`report.py` renders one comment per run, updated in place via `COMMENT_MARKER`. It is
built to stay readable as the suite grows:

* a verdict line first — Bonferroni p, how many measurements are flagged, how many
  setups they are spread over;
* **Needs attention**: every flagged measurement, whatever setup it came from, with the
  setup named per row, and thumbnails for the first `MAX_INLINE_PLOTS` of them. A
  comment with fifty embedded PNGs is unreadable, so the rest are one click away;
* then one `<details>` per experimental setup, worst setup first. The `<summary>` line
  carries the setup's status, its measurement count, how many are flagged and its
  lowest `p_compat`, so a setup can be judged without expanding it;
* plot links are reference-style (`[name][p12]`, definitions at the end) and names are
  stripped of the prefix their setup shares, which keeps rows short in the raw
  markdown; the untruncated name is the link's title;
* past `MAX_COMMENT_CHARS` the tables of setups with nothing flagged collapse to a
  single line each, so the comment cannot exceed GitHub's 65536-character limit. The
  full numbers are always in `summary.json`;
* `summary.json` carries `selected_events` per measurement — the count each sample
  actually binned. Samples that share a selection must agree; when one of them comes
  out empty, that column says so immediately (this is what identified the beam-particle
  bug, where `Q2` binned 8 events against `Tpi`'s 2547 from the same selection);
* a setup whose job crashed is reported, not dropped. `aggregate` runs with `always()`,
  the merge takes `--config` to learn what was expected, and anything no shard reported
  is listed under **Did not run** with the verdict forced to *incomplete* — then the job
  fails, so the run still goes red. `summary.json` carries `did_not_run` and
  `failed_setups`.

Plot files are named by `plot_basename()`, not by the measurement: `upload-artifact`
rejects a path containing `:`, and the Durham samples carry their reference in the name
(`..._Barreau:1983ht`). Only the file name is sanitised; the measurement keeps its
published name everywhere it is displayed or looked up, and the plot URL is built from
the stored basename.

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
