#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles distribution-based physics-validation driver.

Per run: for each experimental setup, generate events once, histogram each of its
measurements, bootstrap a covariance, compare against the stored 'main' baseline
(compatibility p_compat) and the data (goodness-of-fit), and render a Sherpa-style
PR comment + summary.json.

Modes: --make-baseline stores the 'main' baseline; --dry-run swaps in the synthetic
adapter (no Achilles/NUISANCE); --selftest runs an injected-regression check.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Dict, List, Optional

import numpy as np
import yaml

from adapters import (DataTable, GeneratedEvents, Nuisance3Adapter,
                      SyntheticAdapter)
from report import (ALPHA, MeasurementResult, MissingMeasurement, Report,
                    ScanReport, plot_basename,
                    VariantMeasurement, VariantSummary)
from stats import (Prediction, bonferroni, bootstrap_covariance, compatibility,
                   goodness_of_fit, shape_compatibility, trial_covariance)


# ---------------------------------------------------------------------------
# Baseline (stored 'main') I/O
# ---------------------------------------------------------------------------

def baseline_path(baseline_dir: str, key: str) -> str:
    return os.path.join(baseline_dir, "baselines", f"{key}.json")


def load_baseline(path: str) -> Optional[dict]:
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def _measurement_baseline(baseline: dict, name: str):
    """Return (main Prediction, DataTable) for a measurement from a baseline dict."""
    entry = baseline["measurements"][name]
    main = Prediction.from_dict(entry["prediction"])
    data = DataTable(values=np.asarray(entry["data"]["values"], dtype=float),
                     covariance=np.asarray(entry["data"]["covariance"], dtype=float))
    return main, data


# Name of the estimator the covariances are built with, recorded in a stored
# baseline: a baseline built with a different one cannot be mixed with fresh feature
# predictions, so `run` treats it as stale and recomputes main inline.
ESTIMATOR = "trial"


def _predict(adapter, generated: GeneratedEvents, measurement: dict, sample,
             n_boot: int, rng: np.random.Generator) -> Prediction:
    """One measurement's prediction and MC covariance.

    Uses the generator's trial counts and cross-section uncertainty when the events
    carry them (`stats.trial_covariance`, which says why that is the correct
    normalisation), and falls back to resampling when they do not -- an event sample
    from some other source, or one written before the counters were recorded.
    """
    response = adapter.response_matrix(measurement, sample.nbins)
    if generated.n_nonzero_trials:
        return trial_covariance(sample.bin_index, sample.weights, sample.nbins,
                                generated.n_nonzero_trials,
                                rel_xsec_err=generated.rel_xsec_err or 0.0,
                                response=response)
    return bootstrap_covariance(sample.bin_index, sample.weights, sample.nbins,
                                n_boot=n_boot, rng=rng, response=response)


def predict_all(adapter, generated: GeneratedEvents, measurements,
                n_boot: int, rng: np.random.Generator) -> Dict[str, Prediction]:
    """Bin a setup's shared events onto all of its measurements, then add errors.

    One call into the adapter for the whole setup: the real one walks the event file
    once with every sample's columns on the same frame, so this is where the batching
    pays off.
    """
    samples = adapter.histogram_many(generated, measurements)
    out, counts = {}, {}
    for m in measurements:
        sample = samples[m["name"]]
        out[m["name"]] = _predict(adapter, generated, m, sample, n_boot, rng)
        # Two samples that share a selection must bin the same number of events; when
        # one of them comes out empty, this is the column that says so.
        counts[m["name"]] = int(sample.bin_index.size)
    return out, counts


def predict(adapter, generated: GeneratedEvents, measurement: dict,
            n_boot: int, rng: np.random.Generator) -> Prediction:
    """One measurement, for callers that do not have the whole setup in hand."""
    return predict_all(adapter, generated, [measurement], n_boot, rng)[0][
        measurement["name"]]


# ---------------------------------------------------------------------------
# make-baseline: compute and store 'main' predictions + data + covariance
# ---------------------------------------------------------------------------

def make_baseline(adapter, config: dict, key: str, seed: int,
                  n_events: int, n_boot: int, out_dir: str) -> str:
    rng = np.random.default_rng(seed)
    measurements: Dict[str, dict] = {}
    for exp in config["experiments"]:
        gen_main = adapter.generate(exp, "main", seed, n_events)  # once per setup
        mains, _ = predict_all(adapter, gen_main, exp["measurements"], n_boot, rng)
        for m in exp["measurements"]:
            name = m["name"]
            main = mains[name]
            data = adapter.data_table(m)
            measurements[name] = {
                "prediction": main.to_dict(),
                "data": {"values": data.values.tolist(),
                         "covariance": data.covariance.tolist()},
            }
        gen_main.cleanup()  # the event file can be many GB; drop it once binned
    baseline = {
        "key": key,
        "seed": seed,
        "events_per_measurement": n_events,
        "n_boot": n_boot,
        "estimator": ESTIMATOR,
        "measurements": measurements,
    }
    path = baseline_path(out_dir, key)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(baseline, fh, indent=2)
    return path


# ---------------------------------------------------------------------------
# run: compare feature branch against stored baseline, emit report
# ---------------------------------------------------------------------------

def run(adapter, config: dict, *, seed: int, n_events: int, n_boot: int,
        baseline: Optional[dict], repo: str, feature_sha: str,
        nuisance_version: str,
        out_dir: str) -> Report:
    # Imported here, not at module scope: only this path renders anything, so the
    # aggregate job (which just merges shard summaries) needs no matplotlib.
    from plots import plot_measurement

    rng = np.random.default_rng(seed + 101)
    results = []
    warnings = []

    for exp in config["experiments"]:
        # Feature events: once per setup, reused by every measurement. Main events
        # are generated (also once) only for measurements with no stored baseline.
        gen_feature = adapter.generate(exp, "feature", seed, n_events)
        features, n_selected = predict_all(adapter, gen_feature, exp["measurements"],
                                           n_boot, rng)

        # Measurements with no stored baseline need main computed inline. Collect them
        # first so that generation, and the pass over its events, happens once.
        # A baseline whose covariances came from a different estimator cannot be
        # compared against fresh predictions -- the two would carry different
        # uncertainties for the same thing -- so recompute main inline instead.
        usable = baseline is not None and baseline.get("estimator") == ESTIMATOR
        stale = [m for m in exp["measurements"]
                 if not usable
                 or m["name"] not in baseline.get("measurements", {})]
        mains = {}
        gen_main = None
        if stale:
            warnings.extend(m["name"] for m in stale)
            gen_main = adapter.generate(exp, "main", seed, n_events)
            mains, _ = predict_all(adapter, gen_main, stale, n_boot, rng)

        for m in exp["measurements"]:
            name = m["name"]
            feature = features[name]

            if name in mains:
                main, data = mains[name], adapter.data_table(m)
            else:
                main, data = _measurement_baseline(baseline, name)

            compat = compatibility(main, feature)
            gof_pr = goodness_of_fit(feature, data.values, data.covariance)
            gof_main = goodness_of_fit(main, data.values, data.covariance)

            os.makedirs(out_dir, exist_ok=True)
            plot_measurement(
                os.path.join(out_dir, plot_basename(name)), name,
                data=data.values, data_cov=data.covariance,
                main=main.values, main_cov=main.covariance,
                feature=feature.values, feature_cov=feature.covariance,
                edges=data.edges, xlabel=data.xlabel, ylabel=data.ylabel,
                subtitle=f"NUISANCE3 {nuisance_version}  |  seed {seed}  |  "
                         f"{n_events:,} events  |  {feature_sha[:8]}")

            results.append(MeasurementResult(
                name=name,
                ndof=compat.ndof,
                chi2_ndof_main=gof_main.chi2_per_ndof,
                chi2_ndof_pr=gof_pr.chi2_per_ndof,
                delta_chi2=gof_pr.chi2 - gof_main.chi2,
                p_compat=compat.pvalue,
                p_data=gof_pr.pvalue,
                plot=plot_basename(name),
                experiment=exp["name"],
                selected_events=n_selected.get(name, 0),
            ))

        gen_feature.cleanup()  # event files can be many GB; drop once binned
        if gen_main is not None:
            gen_main.cleanup()

    extra_header = []
    if warnings:
        extra_header.append(
            f"> ⚠️ No stored baseline for {len(warnings)} measurement(s) "
            f"({', '.join(warnings)}); main was computed inline for this run.")

    report = Report(results=results, repo=repo, feature_sha=feature_sha,
                    nuisance_version=nuisance_version,
                    seed=seed, events_per_measurement=n_events,
                    extra_header=extra_header)

    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# Unweighting scan: the same setup generated once per unweighting scheme
# ---------------------------------------------------------------------------

def _mean_error(pred: Prediction) -> float:
    """Mean per-bin bootstrap 1-sigma — the size of the prediction's error band."""
    return float(np.mean(np.sqrt(np.clip(np.diag(pred.covariance), 0.0, None))))


def run_unweighting_scan(adapter, config: dict, *, seed: int, n_events: int,
                         n_boot: int, repo: str, feature_sha: str,
                         nuisance_version: str, out_dir: str) -> ScanReport:
    """Generate each setup once per unweighting variant and compare them.

    Unweighting is variance reduction, so every variant must reproduce the reference
    variant's distributions; what legitimately differs is the MC noise per event and
    the wall time. Variants share a seed and an event count, so the only difference
    between two runs of a setup is the scheme itself.
    """
    from plots import plot_variants

    spec = config.get("unweighting")
    if not spec:
        raise SystemExit("--unweighter-scan needs an 'unweighting:' block in the config")
    variants = spec["variants"]
    reference = spec["reference"]
    names = [v["name"] for v in variants]
    if reference not in names:
        raise SystemExit(f"unweighting.reference {reference!r} is not one of {names}")

    rng = np.random.default_rng(seed + 202)
    rows: List[VariantMeasurement] = []
    # variant -> per-measurement pieces, rolled up once every setup has been seen.
    per_variant: Dict[str, dict] = {v["name"]: {"pvalues": [], "ess": [], "err": [],
                                                "norm": [], "seconds": 0.0,
                                                "eff": []} for v in variants}

    for exp in config["experiments"]:
        data = {m["name"]: adapter.data_table(m) for m in exp["measurements"]}
        preds: Dict[str, Dict[str, Prediction]] = {}
        samples: Dict[str, Dict[str, object]] = {}

        # One generation per variant, immediately binned into every measurement of
        # the setup so only a single event file is on disk at a time.
        for variant in variants:
            name = variant["name"]
            # seed_offset is normally 0 so every variant shares a stream and only
            # the scheme differs. A null-control variant repeats the reference's
            # options at a different offset, which calibrates the test: its p-values
            # are drawn from the null and should be uniform.
            gen = adapter.generate(exp, name, seed, n_events,
                                   unweighting=variant["options"],
                                   seed_offset=int(variant.get("seed_offset", 0)))
            preds[name], samples[name] = {}, {}
            binned = adapter.histogram_many(gen, exp["measurements"])
            for m in exp["measurements"]:
                sample = binned[m["name"]]
                preds[name][m["name"]] = _predict(adapter, gen, m, sample,
                                                  n_boot, rng)
                samples[name][m["name"]] = sample
            if gen.run is not None:
                if gen.run.seconds:
                    per_variant[name]["seconds"] += gen.run.seconds
                if gen.run.unweight_eff:
                    per_variant[name]["eff"].append(gen.run.unweight_eff)
            gen.cleanup()

        os.makedirs(out_dir, exist_ok=True)
        for m in exp["measurements"]:
            mname = m["name"]
            ref_pred = preds[reference][mname]
            ref_total = float(np.sum(ref_pred.values))
            ref_err = _mean_error(ref_pred)

            for variant in variants:
                vname = variant["name"]
                pred = preds[vname][mname]
                sample = samples[vname][mname]
                compat = compatibility(ref_pred, pred)
                shape = shape_compatibility(ref_pred, pred)
                gof = goodness_of_fit(pred, data[mname].values, data[mname].covariance)
                total = float(np.sum(pred.values))
                row = VariantMeasurement(
                    measurement=mname, variant=vname, ndof=compat.ndof,
                    chi2_ndof_data=gof.chi2_per_ndof,
                    # The reference has nothing to be compared against; NaN keeps it
                    # out of the flagged rows and renders as a dash.
                    p_compat=float("nan") if vname == reference else compat.pvalue,
                    p_shape=float("nan") if vname == reference else shape.pvalue,
                    p_data=gof.pvalue,
                    norm_shift=(total - ref_total) / ref_total if ref_total else float("nan"),
                    ess_fraction=sample.ess_fraction(),
                    max_over_mean=sample.max_over_mean(),
                    mc_error_ratio=_mean_error(pred) / ref_err if ref_err else float("nan"))
                rows.append(row)
                agg = per_variant[vname]
                agg["ess"].append(row.ess_fraction)
                agg["err"].append(row.mc_error_ratio)
                agg["norm"].append(row.norm_shift)
                if vname != reference:
                    agg["pvalues"].append(compat.pvalue)

            plot_variants(
                os.path.join(out_dir, f"{mname}.unweighting.png"), mname,
                data=data[mname].values, data_cov=data[mname].covariance,
                variants={v["name"]: (preds[v["name"]][mname].values,
                                      preds[v["name"]][mname].covariance)
                          for v in variants},
                reference=reference,
                edges=data[mname].edges, xlabel=data[mname].xlabel,
                ylabel=data[mname].ylabel,
                subtitle=f"NUISANCE3 {nuisance_version}  |  seed {seed}  |  "
                         f"{n_events:,} events/variant  |  {feature_sha[:8]}")

    report = ScanReport(
        summaries=_scan_summaries(
            [(v["name"], v["options"]) for v in variants], rows, reference,
            {name: (agg["seconds"], agg["eff"]) for name, agg in per_variant.items()}),
        rows=rows, reference=reference, repo=repo, feature_sha=feature_sha,
        nuisance_version=nuisance_version, seed=seed,
        events_per_measurement=n_events)
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# The variant that repeats the reference at a different seed. Its p-values come
# from the null hypothesis, so they set the floor any real scheme is judged against.
NULL_CONTROL = "null-control"


def _scan_summaries(variants, rows: List[VariantMeasurement], reference: str,
                    runtime: Dict[str, tuple]) -> List[VariantSummary]:
    """Roll the per-measurement scan rows up into one line per variant.

    ``variants`` is an ordered ``(name, options)`` sequence; ``runtime`` maps a
    variant to ``(total_seconds, [per-setup efficiencies])``. Shared by the direct
    run and by ``merge_scan_shards`` so a sharded scan reports identically.
    """
    by_variant: Dict[str, List[VariantMeasurement]] = {}
    for r in rows:
        by_variant.setdefault(r.variant, []).append(r)

    summaries = []
    for name, options in variants:
        mine = by_variant.get(name, [])
        pv = [r.p_compat for r in mine if r.p_compat == r.p_compat]  # drops the ref's NaN
        ps = [r.p_shape for r in mine if r.p_shape == r.p_shape]
        seconds, effs = runtime.get(name, (0.0, []))
        summaries.append(VariantSummary(
            variant=name, options=dict(options), n_measurements=len(mine),
            p_worst=float(np.min(pv)) if pv else float("nan"),
            p_overall=bonferroni(pv) if pv else float("nan"),
            p_shape_worst=float(np.min(ps)) if ps else float("nan"),
            p_shape_overall=bonferroni(ps) if ps else float("nan"),
            n_flagged=sum(1 for p in pv if p < ALPHA),
            ess_fraction=float(np.nanmedian([r.ess_fraction for r in mine]))
            if mine else float("nan"),
            mc_error_ratio=float(np.nanmedian([r.mc_error_ratio for r in mine]))
            if mine else float("nan"),
            max_norm_shift=max((r.norm_shift for r in mine), key=abs)
            if mine else float("nan"),
            seconds=seconds or None,
            unweight_eff=float(np.min(effs)) if effs else None,
            is_reference=name == reference,
            is_null_control=name == NULL_CONTROL))
    return summaries


def merge_scan_shards(shard_paths, out_dir: str) -> ScanReport:
    """Combine per-experiment unweighting-scan ``summary.json`` files into one report.

    The scan shards on the same axis as the branch comparison — one experimental
    setup per job — so the rows just concatenate; only the per-variant rollup has to
    be recomputed across the whole family.
    """
    shards = []
    for p in shard_paths:
        with open(p) as fh:
            shards.append(json.load(fh))
    if not shards:
        raise SystemExit("merge_scan_shards: no shard summaries given")
    if any(s.get("kind") != "unweighting-scan" for s in shards):
        raise SystemExit("merge_scan_shards: not every shard is an unweighting scan")

    head = shards[0]
    rows = [VariantMeasurement(**r) for s in shards for r in s["rows"]]
    # Variant order and options come from the first shard; every shard runs the
    # same variant list, so this is just the display order.
    variants = [(v["variant"], v["options"]) for v in head["variants"]]
    runtime: Dict[str, tuple] = {name: (0.0, []) for name, _ in variants}
    for s in shards:
        for v in s["variants"]:
            seconds, effs = runtime[v["variant"]]
            runtime[v["variant"]] = (seconds + (v["seconds"] or 0.0),
                                     effs + ([v["unweight_eff"]]
                                             if v["unweight_eff"] else []))

    report = ScanReport(
        summaries=_scan_summaries(variants, rows, head["reference"], runtime),
        rows=rows, reference=head["reference"], repo=head["repo"],
        feature_sha=head["feature_sha"],
        nuisance_version=head["nuisance_version"], seed=head["seed"],
        events_per_measurement=head["events_per_measurement"])
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

# Per-sample overrides declared as top-level maps in the config, folded onto the
# measurement they name so the adapter sees them next to the sample.
_MEASUREMENT_KEYS = ("data_scale", "bin_edges", "bin_widths", "solid_angle",
                     "smearing")


def _load_config(path: str) -> dict:
    with open(path) as fh:
        config = yaml.safe_load(fh)
    # A measurement is normally just its NUISANCE3 sample name; allow a mapping too
    # for the odd entry that needs extra keys (e.g. a dryrun nbins override).
    for exp in config["experiments"]:
        exp["measurements"] = [m if isinstance(m, dict) else {"name": m}
                               for m in exp["measurements"]]

    known = {m["name"] for exp in config["experiments"] for m in exp["measurements"]}
    for key in _MEASUREMENT_KEYS:
        entries = config.get(key) or {}
        unknown = set(entries) - known
        if unknown:
            raise SystemExit(f"{key}: no such measurement(s) {sorted(unknown)}")
        for exp in config["experiments"]:
            for m in exp["measurements"]:
                if m["name"] in entries:
                    value = entries[m["name"]]
                    # smearing names a csv, relative to the config it is declared in.
                    if key == "smearing":
                        value = os.path.join(os.path.dirname(os.path.abspath(path)),
                                             value)
                    m[key] = value
    return config


def _filter_config(config: dict, only_experiments, only_measurements,
                   only_variants=None) -> dict:
    """Shard the config by experiment (primary) and/or by measurement name.

    ``--only-experiment`` selects whole experimental setups (the CI shard axis);
    ``--only`` further narrows to individual measurements within them;
    ``--only-variant`` narrows the unweighting scan's variant list.
    """
    experiments = config["experiments"]

    if only_experiments:
        wanted = set(only_experiments)
        experiments = [e for e in experiments if e["name"] in wanted]
        missing = wanted - {e["name"] for e in experiments}
        if missing:
            raise SystemExit(
                f"--only-experiment names not in config: {sorted(missing)}")

    if only_measurements:
        wanted = set(only_measurements)
        kept = []
        for e in experiments:
            ms = [m for m in e["measurements"] if m["name"] in wanted]
            if ms:
                kept.append({**e, "measurements": ms})
        found = {m["name"] for e in kept for m in e["measurements"]}
        missing = wanted - found
        if missing:
            raise SystemExit(f"--only names not in config: {sorted(missing)}")
        experiments = kept

    config = {**config, "experiments": experiments}

    if only_variants:
        spec = config.get("unweighting")
        if not spec:
            raise SystemExit("--only-variant needs an 'unweighting:' block in the config")
        wanted = set(only_variants) | {spec["reference"]}  # the reference is required
        kept_variants = [v for v in spec["variants"] if v["name"] in wanted]
        missing = wanted - {v["name"] for v in kept_variants}
        if missing:
            raise SystemExit(f"--only-variant names not in config: {sorted(missing)}")
        config["unweighting"] = {**spec, "variants": kept_variants}

    return config


def merge_shards(shard_paths, out_dir: str,
                 config: Optional[dict] = None) -> Report:
    """Combine per-shard ``summary.json`` files into one Report.

    A shard whose job failed leaves no summary behind, and the aggregate still has to
    run: reporting the setups that did report, with the rest marked as not run, beats
    producing nothing at all. Paths that do not exist are therefore skipped (the caller
    passes a glob), and ``config`` — when given — says which measurements were expected,
    so the ones nobody reported can be named instead of silently vanishing.
    """
    shards = []
    for p in shard_paths:
        if not os.path.exists(p):
            continue  # that shard's job failed, or the glob matched nothing
        with open(p) as fh:
            shards.append(json.load(fh))

    expected = ([(m["name"], exp["name"])
                 for exp in config["experiments"] for m in exp["measurements"]]
                if config else [])
    if not shards:
        if not expected:
            raise SystemExit("merge_shards: no shard summaries and no config to "
                             "say what was expected")
        # Every shard failed. Still emit a comment, so the run says what happened.
        head = {"repo": os.environ.get("GITHUB_REPOSITORY", "AchillesGen/Achilles"),
                "feature_sha": os.environ.get("GITHUB_SHA", "unknown"),
                "nuisance_version": "unknown", "seed": 0,
                "events_per_measurement": 0}
        shards = []
    else:
        head = shards[0]

    results = []
    for s in shards:
        for m in s["measurements"]:
            results.append(MeasurementResult(
                name=m["name"], ndof=m["ndof"],
                chi2_ndof_main=m["chi2_ndof_main"], chi2_ndof_pr=m["chi2_ndof_pr"],
                delta_chi2=m["delta_chi2"], p_compat=m["p_compat"],
                p_data=m["p_data"], plot=m.get("plot"),
                experiment=m.get("experiment", ""),
                selected_events=m.get("selected_events", 0)))

    reported = {r.name for r in results}
    did_not_run = [MissingMeasurement(name=name, experiment=exp)
                   for name, exp in expected if name not in reported]

    report = Report(results=results, repo=head["repo"],
                    feature_sha=head["feature_sha"],
                    nuisance_version=head["nuisance_version"],
                    seed=head["seed"],
                    events_per_measurement=head["events_per_measurement"],
                    did_not_run=did_not_run)
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="measurements.yml",
                   help="measurement-list YAML")
    p.add_argument("--workdir", default="physval-work")
    p.add_argument("--out-dir", default="physval-out")
    p.add_argument("--baseline-dir", default="physval-baselines",
                   help="checkout of the physval-baselines branch")
    p.add_argument("--key", default="dev",
                   help="baseline key {nuisance}-{data}-{confighash}-{sha}")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--events", type=int, default=200000)
    p.add_argument("--n-boot", type=int, default=200)
    p.add_argument("--repo", default="AchillesGen/Achilles")
    p.add_argument("--feature-sha", default=os.environ.get("GITHUB_SHA", "local"))
    p.add_argument("--nuisance-version", default=os.environ.get("NUISANCE_VERSION",
                                                                "unknown"))
    p.add_argument("--only-experiment", action="append", default=[],
                   dest="only_experiments",
                   help="run only this experimental setup (repeatable); the CI "
                        "shard axis — its events are generated once and reused")
    p.add_argument("--only", action="append", default=[],
                   dest="only_measurements",
                   help="run only this measurement (repeatable); narrows within "
                        "the selected experiment(s)")
    p.add_argument("--only-variant", action="append", default=[],
                   dest="only_variants",
                   help="run only this unweighting variant (repeatable); the "
                        "reference variant is always kept")
    p.add_argument("--unweighter-scan", action="store_true",
                   help="generate each setup once per Options/Unweighting variant "
                        "and compare them against the reference variant")
    p.add_argument("--merge", nargs="+", default=None,
                   help="merge these shard summary.json files into one report")
    p.add_argument("--dry-run", action="store_true",
                   help="use the synthetic adapter (no Achilles/NUISANCE)")
    p.add_argument("--make-baseline", action="store_true",
                   help="produce the stored main baseline instead of a comparison")
    p.add_argument("--selftest", action="store_true")
    return p


def _adapter_from_args(args):
    if args.dry_run:
        return SyntheticAdapter(base_seed=args.seed)
    return Nuisance3Adapter(workdir=args.workdir)


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    if args.selftest:
        return _selftest()

    if args.merge:
        if args.unweighter_scan:
            scan = merge_scan_shards(args.merge, out_dir=args.out_dir)
            print(f"merged {len(args.merge)} scan shard(s): rows={len(scan.rows)} "
                  f"biased={[s.variant for s in scan.biased()] or 'none'}")
            print(f"wrote {args.out_dir}/comment.md and {args.out_dir}/summary.json")
            return 0
        # The config is what says which measurements were expected, so a shard whose
        # job failed can be named in the comment rather than quietly dropped. A scoped
        # run narrows that expectation: the setups nobody asked for are not missing.
        merge_config = (_load_config(args.config)
                        if os.path.exists(args.config) else None)
        if merge_config and (args.only_experiments or args.only_measurements):
            merge_config = _filter_config(merge_config, args.only_experiments,
                                          args.only_measurements)
        report = merge_shards(args.merge, out_dir=args.out_dir, config=merge_config)
        print(f"merged {len(report.results)} measurement(s): "
              f"p_overall={report.p_overall():.4g} "
              f"flagged={report.n_flagged()}/{len(report.results)}")
        if report.did_not_run:
            print(f"did not report: {len(report.did_not_run)} measurement(s) in "
                  f"{len(report.failed_setups())} setup(s): "
                  f"{', '.join(report.failed_setups())}")
        print(f"wrote {args.out_dir}/comment.md and {args.out_dir}/summary.json")
        return 0

    config = _filter_config(_load_config(args.config), args.only_experiments,
                            args.only_measurements, args.only_variants)
    adapter = _adapter_from_args(args)

    if args.unweighter_scan:
        report = run_unweighting_scan(
            adapter, config, seed=args.seed, n_events=args.events,
            n_boot=args.n_boot, repo=args.repo, feature_sha=args.feature_sha,
            nuisance_version=args.nuisance_version, out_dir=args.out_dir)
        biased = [s.variant for s in report.biased()]
        print(f"reference={report.reference} "
              f"variants={len(report.summaries)} rows={len(report.rows)} "
              f"biased={biased or 'none'}")
        print(f"wrote {args.out_dir}/comment.md and {args.out_dir}/summary.json")
        return 0

    if args.make_baseline:
        path = make_baseline(adapter, config, key=args.key, seed=args.seed,
                             n_events=args.events, n_boot=args.n_boot,
                             out_dir=args.baseline_dir)
        print(f"wrote baseline: {path}")
        return 0

    baseline = load_baseline(baseline_path(args.baseline_dir, args.key))
    report = run(adapter, config, seed=args.seed, n_events=args.events,
                 n_boot=args.n_boot, baseline=baseline, repo=args.repo,
                 feature_sha=args.feature_sha,
                 nuisance_version=args.nuisance_version,
                 out_dir=args.out_dir)

    po = report.p_overall()
    print(f"p_overall={po:.4g} flagged={report.n_flagged()}/{len(report.results)} "
          f"overall_ok={report.overall_ok()}")
    print(f"wrote {args.out_dir}/comment.md and {args.out_dir}/summary.json")
    # Advisory only: exit 0 regardless so the build stays green (per plan).
    return 0


# ---------------------------------------------------------------------------
# End-to-end self-test: baseline -> compatible run, then injected-regression run
# ---------------------------------------------------------------------------

def _selftest() -> int:
    import tempfile

    # Two experimental setups: one whose feature generation is unchanged and one
    # with an injected physics shift; each carries a single measurement.
    config = {"experiments": [
        {"name": "SYNTH_stable_exp", "dryrun": {"feature_shift": 0.0},
         "measurements": [{"name": "SYNTH_stable", "dryrun": {"nbins": 12}}]},
        {"name": "SYNTH_regressed_exp", "dryrun": {"feature_shift": 0.05},
         "measurements": [{"name": "SYNTH_regressed", "dryrun": {"nbins": 12}}]},
    ]}
    with tempfile.TemporaryDirectory() as tmp:
        adapter = SyntheticAdapter(base_seed=7)
        bpath = make_baseline(adapter, config, key="selftest", seed=7,
                              n_events=60000, n_boot=150, out_dir=tmp)
        baseline = load_baseline(bpath)
        report = run(adapter, config, seed=7, n_events=60000, n_boot=150,
                     baseline=baseline, repo="AchillesGen/Achilles",
                     feature_sha="deadbeefcafef00d", nuisance_version="selftest",
                     out_dir=tmp)
        summary = report.to_summary_dict()
        by_name = {m["name"]: m for m in summary["measurements"]}

        checks = {
            "baseline written": os.path.exists(bpath),
            "stable is compatible": by_name["SYNTH_stable"]["status"] == "compatible",
            "regressed is flagged": by_name["SYNTH_regressed"]["status"] != "compatible",
            "regressed p_compat < 0.05": by_name["SYNTH_regressed"]["p_compat"] < 0.05,
            "comment written": os.path.exists(os.path.join(tmp, "comment.md")),
            "summary written": os.path.exists(os.path.join(tmp, "summary.json")),
        }
        checks.update(_selftest_scan(tmp))
        for name, ok in checks.items():
            print(f"[{'ok' if ok else 'FAIL'}] {name}")
        ok = all(checks.values())
        print("SELFTEST:", "PASS" if ok else "FAIL")
        return 0 if ok else 1


def _selftest_scan(tmp: str) -> dict:
    """Unweighting scan on synthetic events: unbiased, sharper per event, slower.

    The synthetic adapter applies the real cap rules (``adapters.unweighting_cap``),
    so this checks both the scan plumbing and the trade the schemes are supposed to
    make: at a fixed accepted-event count, accept-with-excess unweighting leaves the
    distribution alone and buys precision, paid for in extra trials.
    """
    config = {
        "unweighting": {
            "reference": "weighted",
            "variants": [
                {"name": "weighted", "options": {"Name": "None"}},
                {"name": "percentile-99", "options": {"Name": "Percentile",
                                                      "percentile": 99}},
                {"name": "excess-1e-2", "options": {"Name": "Excess",
                                                    "epsilon": 0.01}},
                {"name": "tailfrac-1e-2", "options": {"Name": "TailFraction",
                                                      "epsilon": 0.01}},
            ],
        },
        "experiments": [
            {"name": "SYNTH_exp", "dryrun": {"feature_shift": 0.0},
             "measurements": [{"name": "SYNTH_scan", "dryrun": {"nbins": 12}}]},
        ],
    }
    out = os.path.join(tmp, "scan")
    report = run_unweighting_scan(SyntheticAdapter(base_seed=3), config, seed=3,
                                 n_events=40000, n_boot=200,
                                 repo="AchillesGen/Achilles",
                                 feature_sha="deadbeefcafef00d",
                                 nuisance_version="selftest", out_dir=out)
    by_name = {s.variant: s for s in report.summaries}
    ref, cut = by_name["weighted"], by_name["percentile-99"]
    return {
        "scan: every variant reported": len(report.summaries) == 4,
        "scan: no variant biases the distribution": report.biased() == [],
        "scan: unweighting raises ESS/event": cut.ess_fraction > ref.ess_fraction,
        "scan: sharper per accepted event": cut.mc_error_ratio < 1.0,
        "scan: paid for in wall time": cut.seconds > ref.seconds,
        "scan: normalisation preserved": abs(cut.max_norm_shift) < 0.05,
        "scan: comment written": os.path.exists(os.path.join(out, "comment.md")),
        "scan: overlay plot written": os.path.exists(
            os.path.join(out, "SYNTH_scan.unweighting.png")),
    }


if __name__ == "__main__":
    raise SystemExit(main())
