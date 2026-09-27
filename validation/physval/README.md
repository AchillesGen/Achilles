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

`p_compat` is calibrated in `stats.py --selftest` on **i.i.d.** synthetic events,
where the bootstrap is exact. Real Achilles runs are not i.i.d.: the events of a run
share one adapted VEGAS grid and one weight cap, and the whole histogram is scaled by
that run's flux-averaged cross section. The bootstrap resamples *within* a run, so it
carries none of the run-to-run scatter in that overall scale — measured at up to ±5%
between independent runs of the same setup, against per-bin bootstrap errors of ~1%.
Two independent runs of the *identical* configuration therefore come out with
p_compat ≈ 0, and every variant looks "biased".

So the variant list carries a `null-control`: the reference's own options at
`seed_offset: 1`. Its rows are draws from the null hypothesis, and the flagging
thresholds are floored at whatever it scores (`min(alpha, null)`, so a well-behaved
control leaves α alone and never tightens it). Nothing is called out for doing as
well as an identical rerun. `p (shape only)` — the same χ² with the normalisation
divided out and one dof given up for it — separates the two failure modes: a scheme
that moved a distribution fails both columns, two runs that merely disagree on the
total cross section fail only the first.

The same caveat applies to the **branch comparison**, which also compares two
independent runs (`main` at `seed`, feature at `seed + 1`) whenever there is no
stored baseline. Its `p_compat` is anticonservative for the same reason.

Reading the summary table:

- **p (vs reference)** — Bonferroni over the setup's measurements, judged against the
  null control. A flag here means the scheme *changed the physics*, which is a bug
  rather than a trade-off.
- **p (shape only)** — as above with the overall normalisation fitted out.
- **ESS/event** — Kish effective sample size `(Σ|w|)²/(N Σw²)` on the raw generator
  weights, i.e. the statistical power the scheme delivers per accepted event
  (`1.0` = perfect unit weights). Computed pre-normalisation, so the bin-width
  division does not masquerade as weight spread.
- **MC error** — mean bootstrap σ relative to the reference's. At a fixed accepted-
  event count, a harsher cap buys precision here and pays for it in **wall**.
- **Δnorm** — change in the integrated cross section, which no scheme should move.
- **max/mean w** — the heaviest surviving overweight, the tail the cap left behind.

`adapters.unweighting_cap` reimplements the three cap rules in numpy so `--dry-run`
is a real test of the scan (and a second opinion on the C++).

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

* `bin_edges` — the sample's histogram is indexed by bin *number* (NCπ⁰ stacks several
  blocks), so the widths from the NUISANCE binning are 1 and the prediction has to be
  divided by the real ones.
* `solid_angle` — the Durham electron data is published per steradian while the
  selection integrates over its ±4° acceptance window.
* `smearing` — a Wiener-SVD unfolded measurement is only comparable to `A_C ×
  prediction`. The matrix is applied inside `bootstrap_covariance`, per replica, so it
  reaches the MC covariance as `A C Aᵀ` and not just the central values. It is never
  applied to the data.
* `data_scale` — one shipped table (MiniBooNE dσ/dQ²) is written 1e6 low; the factor
  multiplies the data and its covariance, not the prediction.

### Samples that need patched NUISANCE2

`MicroBooNE_NCpi0_*` and `ElectronData_*` do not work against NUISANCE2 as shipped:
NCπ⁰ sets its bin index only in `FillHistograms`, which the legacy record never calls,
and the Durham samples cut on variables the record fills *after* `isSignal`. Both fail
silently — an empty prediction, not an error. The fixes are upstream as
NUISANCEMC/nuisance#115 and #116 and are carried as patches in the physval image on top
of its pinned NUISANCE2 sha; `versions.json`'s `nuisance2.patches` records them. The
δp_n unit mixing that used to make `MicroBooNE_CC1Mu1p_XSec_1DDeltaPn_nu` incomparable
to data is fixed upstream (#114) and needs NUISANCE2 ≥ `86c64b44` in the image.

### Achilles' libraries must precede the image's

The image exports `LD_LIBRARY_PATH=/opt/root/lib:/opt/nuisance3/lib:/opt/nuisance2/lib:…`,
which outranks the binary's `RUNPATH`. NUISANCE2 ships `libspdlog` built against
`fmt v10` while Achilles bundles `fmt v11`; loading both leaves Achilles calling a
spdlog with a mismatched fmt and it **segfaults inside `InitializeLogging`, printing
nothing at all** (the splash is lost to stdout buffering, so it looks like an instant
silent crash). Putting Achilles' own `lib` first fixes it. `Nuisance3Adapter` does this
itself when it spawns achilles, and the CI build job smoke-tests `achilles --version`
so a regression fails early instead of mid-generation.

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
  full numbers are always in `summary.json`.

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
