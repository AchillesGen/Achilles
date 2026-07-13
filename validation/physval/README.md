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
| `adapters.py` | The boundary to Achilles/NUISANCE3: `SyntheticAdapter` (tests/dry-run) and `Nuisance3Adapter` (real path — skeleton). |
| `physval.py` | Driver: config → generate → stats → report; plus `--make-baseline`. |
| `measurements.yml` | The measurement list (starter set; sample names are placeholders). |

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

## Adapter boundary — external hand-offs (TODO)

`adapters.Nuisance3Adapter` is a documented skeleton. To make the real path live:

1. **NUISANCE3 container image** — pin an image with NUISANCE3 + ROOT + HepMC3 +
   ProSelecta and the legacy interface enabled. Set `NUISANCE_IMAGE` /
   `--nuisance-image`.
2. **Per-measurement mapping** in `measurements.yml`: replace the `TODO_*`
   `nuisance_sample` names with real NUISANCE3 sample names and `achilles_run` with
   real Achilles run cards.
3. Implement `Nuisance3Adapter.generate` / `.data_table`:
   - render the Achilles run card with the pinned `seed` + `n_events`,
     `achilles <card>` → NuHepMC-v1.0 event file (stream / gzip; delete after);
   - run the NUISANCE3 comparison for the sample over that file; return the
     per-event `(bin_index, weight)` arrays and the sample's published
     `DataTable(values, covariance)`.

Everything downstream (bootstrap, `p_compat`, Bonferroni, table, plots) is already
implemented and tested against `SyntheticAdapter`.

## `physval-baselines` branch

Durable store for the `main` baseline and the run plots (created on demand):

```
baselines/<key>.json        # main central histograms + bootstrap covariance + data
plots/<feature-sha>/*.png    # per-run new/old/data overlays (inline in the comment)
```

`<key> = {nuisance_version}-{data_release}-{config_hash}-{sha}`. The baseline JSON is
tiny and kept indefinitely; the plots are pruned on a time window by
`physval-prune.yml`. If the feature run's NUISANCE3 / data version doesn't match the
stored baseline's, recompute `main` inline (the driver already falls back and notes
it in the comment header).

## CI workflows

- `physval.yml` — triggered by a `!physval` commit-message marker,
  `workflow_dispatch`, or nightly `schedule`. Builds Achilles, shards event
  generation across a matrix (one job per measurement), aggregates, commits plots to
  `physval-baselines/plots/<sha>/`, and posts/updates the PR comment.
- `physval-baseline-refresh.yml` — on push to `main`, recompute and store the
  baseline (predictions + bootstrap covariance + data) on `physval-baselines`.
- `physval-prune.yml` — scheduled pruning of old `plots/<sha>/` folders.

See the plan for the full rationale (bootstrap vs same-seed, merge-base attribution,
public-runner compute budget).
