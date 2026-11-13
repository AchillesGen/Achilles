#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles distribution-based physics-validation driver.

For each experimental setup: generate events once, bin every measurement, and compare
the PR against the stored main baseline (``p_compat``) and the data (``p_data``).

  --make-baseline   write main's per-setup predictions to --out-dir
  (default)         compare against <baseline-dir>/baselines/<key>/<setup>.json
  --merge           combine shard summaries into one comment
  --dry-run         synthetic adapter, no Achilles/NUISANCE
"""

from __future__ import annotations

import argparse
import json
import os
from typing import List, Optional

import yaml

from adapters import GeneratedEvents, Nuisance3Adapter, SyntheticAdapter
from report import MeasurementResult, MissingMeasurement, Report, plot_basename
from stats import Prediction, compatibility, goodness_of_fit, trial_covariance

# Bump when a stored baseline stops being comparable with fresh predictions.
BASELINE_SCHEMA = 1


# ---------------------------------------------------------------------------
# Predictions and baselines
# ---------------------------------------------------------------------------

def predict_all(adapter, generated: GeneratedEvents, measurements) -> dict:
    """Predictions by name, leaving out the measurements that selected no events."""
    if not generated.n_nonzero_trials:
        raise RuntimeError("event file carries no trial counts; the MC covariance "
                           "cannot be built")
    samples = adapter.histogram_many(generated, measurements)
    preds = {}
    for m in measurements:
        s = samples[m["name"]]
        if s.bin_index.size:
            preds[m["name"]] = trial_covariance(
                s.bin_index, s.weights, s.nbins, generated.n_nonzero_trials,
                rel_xsec_err=generated.rel_xsec_err or 0.0,
                response=adapter.response_matrix(m, s.nbins))
    return preds


def no_events(exp: dict, preds: dict) -> List[MissingMeasurement]:
    """The setup's measurements without a prediction, announced to the CI log."""
    empty = [MissingMeasurement(m["name"], exp["name"], m.get("may_be_empty", ""))
             for m in exp["measurements"] if m["name"] not in preds]
    for m in empty:
        print(f"::{'warning' if m.reason else 'error'}::{m.name}: no events selected")
    return empty


def baseline_path(baseline_dir: str, key: str, experiment: str) -> str:
    return os.path.join(baseline_dir, "baselines", key, f"{experiment}.json")


def make_baseline(adapter, config: dict, *, seed: int, n_events: int, main_sha: str,
                  out_dir: str) -> List[MissingMeasurement]:
    """Write one ``<setup>.json`` of main predictions per setup; return the empty ones."""
    os.makedirs(out_dir, exist_ok=True)
    empty = []
    for exp in config["experiments"]:
        gen = adapter.generate(exp, "main", seed, n_events)
        preds = predict_all(adapter, gen, exp["measurements"])
        gen.cleanup()
        empty += no_events(exp, preds)
        path = os.path.join(out_dir, f"{exp['name']}.json")
        with open(path, "w") as fh:
            json.dump({"schema": BASELINE_SCHEMA, "adapter": adapter.name,
                       "main_sha": main_sha, "seed": seed, "events": n_events,
                       "measurements": {k: p.to_dict() for k, p in preds.items()}},
                      fh)
        print(f"wrote baseline: {path}")
    return empty


def load_baseline(path: str, adapter_name: str) -> Optional[dict]:
    """The stored baseline at ``path``, or None if absent or not comparable."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        b = json.load(fh)
    if b.get("schema") != BASELINE_SCHEMA or b.get("adapter") != adapter_name:
        return None
    return b


# ---------------------------------------------------------------------------
# Comparison run
# ---------------------------------------------------------------------------

def run(adapter, config: dict, *, seed: int, n_events: int, baseline_dir: str,
        key: str, feature_sha: str, nuisance_version: str, out_dir: str) -> Report:
    from plots import plot_measurement  # only this path needs matplotlib

    os.makedirs(out_dir, exist_ok=True)
    results, missing, empty, main_shas = [], [], [], set()
    for exp in config["experiments"]:
        try:
            gen = adapter.generate(exp, "feature", seed, n_events)
            features = predict_all(adapter, gen, exp["measurements"])
            gen.cleanup()
        except Exception as err:  # report the setup as failed, keep the rest
            print(f"::error::{exp['name']}: {err}")
            missing += [MissingMeasurement(m["name"], exp["name"], str(err))
                        for m in exp["measurements"]]
            continue
        empty += no_events(exp, features)

        base = load_baseline(baseline_path(baseline_dir, key, exp["name"]), adapter.name)
        stored = base["measurements"] if base else {}
        if base:
            main_shas.add(base["main_sha"])

        for m in exp["measurements"]:
            name = m["name"]
            if name not in features:
                continue
            feature = features[name]
            data = adapter.data_table(m)
            gof_pr = goodness_of_fit(feature, data.values, data.covariance)
            main = Prediction.from_dict(stored[name]) if name in stored else None
            if main is not None and main.nbins != feature.nbins:
                main = None  # binning changed since the baseline
            r = MeasurementResult(name=name, experiment=exp["name"], ndof=gof_pr.ndof,
                                  chi2_ndof_pr=gof_pr.chi2_per_ndof, p_data=gof_pr.pvalue)
            if main is not None:
                gof_main = goodness_of_fit(main, data.values, data.covariance)
                r.chi2_ndof_main = gof_main.chi2_per_ndof
                r.delta_chi2 = gof_pr.chi2 - gof_main.chi2
                r.p_compat = compatibility(main, feature).pvalue
            results.append(r)
            plot_measurement(
                os.path.join(out_dir, plot_basename(name)), name,
                data=data.values, data_cov=data.covariance,
                feature=feature.values, feature_cov=feature.covariance,
                main=None if main is None else main.values,
                main_cov=None if main is None else main.covariance,
                edges=data.edges,
                subtitle=f"NUISANCE3 {nuisance_version}  |  seed {seed}  |  "
                         f"{n_events:,} events  |  {feature_sha[:8]}")

    report = Report(results=results, feature_sha=feature_sha,
                    nuisance_version=nuisance_version, seed=seed,
                    events_per_measurement=n_events, did_not_run=missing, no_events=empty,
                    baseline_shas=sorted(main_shas), image_key=key)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


def _run_url() -> str:
    env = os.environ
    if not env.get("GITHUB_RUN_ID"):
        return ""
    return (f"{env.get('GITHUB_SERVER_URL', 'https://github.com')}/"
            f"{env.get('GITHUB_REPOSITORY')}/actions/runs/{env['GITHUB_RUN_ID']}")


def merge_shards(shard_paths, out_dir: str, config: dict) -> Report:
    """Combine shard summaries; ``config`` names the measurements no shard reported."""
    shards = []
    for p in shard_paths:
        if os.path.exists(p):
            with open(p) as fh:
                shards.append(json.load(fh))

    results = [MeasurementResult(**m) for s in shards for m in s["results"]]
    missing = [MissingMeasurement(**m) for s in shards for m in s["did_not_run"]]
    empty = [MissingMeasurement(**m) for s in shards for m in s["no_events"]]
    seen = {m.name for m in results + missing + empty}
    for exp in config["experiments"]:
        missing += [MissingMeasurement(m["name"], exp["name"], "the shard's job failed")
                    for m in exp["measurements"] if m["name"] not in seen]

    head = shards[0] if shards else {}
    meta = {k: head[k] for k in ("feature_sha", "nuisance_version", "seed",
                                 "events_per_measurement", "image_key") if k in head}
    report = Report(results=results, did_not_run=missing, no_events=empty,
                    baseline_shas=sorted({sha for s in shards for sha in s["baseline_shas"]}),
                    run_url=_run_url(), **meta)
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# Config and CLI
# ---------------------------------------------------------------------------

# Top-level per-sample maps in the config, folded onto the measurement they name.
_MEASUREMENT_KEYS = ("data_scale", "solid_angle", "may_be_empty")


def load_config(path: str, experiments=()) -> dict:
    """The config, narrowed to the named setups if any are given."""
    with open(path) as fh:
        config = yaml.safe_load(fh)
    for exp in config["experiments"]:
        exp["measurements"] = [{"name": m} for m in exp["measurements"]]
    by_name = {m["name"]: m for exp in config["experiments"] for m in exp["measurements"]}
    for key in _MEASUREMENT_KEYS:
        entries = config.get(key) or {}
        unknown = set(entries) - set(by_name)
        if unknown:
            raise SystemExit(f"{key}: no such measurement(s) {sorted(unknown)}")
        for name, value in entries.items():
            by_name[name][key] = value
    if experiments:
        unknown = set(experiments) - {e["name"] for e in config["experiments"]}
        if unknown:
            raise SystemExit(f"setups not in config: {sorted(unknown)}")
        config["experiments"] = [e for e in config["experiments"] if e["name"] in experiments]
    return config


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="measurements.yml")
    p.add_argument("--workdir", default="physval-work")
    p.add_argument("--out-dir", default="physval-out")
    p.add_argument("--baseline-dir", default="physval-artifacts",
                   help="checkout of the physval-artifacts repo")
    p.add_argument("--key", default="dev", help="baseline key (the image digest in CI)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--events", type=int, default=200000)
    p.add_argument("--feature-sha", default=os.environ.get("GITHUB_SHA", "local"))
    p.add_argument("--nuisance-version",
                   default=os.environ.get("NUISANCE_VERSION", "unknown"))
    p.add_argument("--only-experiment", action="append", default=[],
                   dest="only_experiments", help="run only this setup (repeatable)")
    p.add_argument("--merge", nargs="+", default=None,
                   help="merge these shard summary.json files")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--make-baseline", action="store_true")
    args = p.parse_args(argv)

    config = load_config(args.config, args.only_experiments)
    if args.merge:
        report = merge_shards(args.merge, args.out_dir, config)
        print(f"merged {len(report.results)} measurement(s); "
              f"flagged={len(report.flagged())}; did not report="
              f"{', '.join(report.failed_setups()) or 'none'}; "
              f"no events={len(report.no_events)}")
        return 0

    adapter = (SyntheticAdapter(base_seed=args.seed) if args.dry_run
               else Nuisance3Adapter(workdir=args.workdir))
    common = dict(seed=args.seed, n_events=args.events, out_dir=args.out_dir)
    if args.make_baseline:
        empty = make_baseline(adapter, config, main_sha=args.feature_sha, **common)
        return 1 if any(not m.reason for m in empty) else 0

    report = run(adapter, config, baseline_dir=args.baseline_dir, key=args.key,
                 feature_sha=args.feature_sha, nuisance_version=args.nuisance_version,
                 **common)
    print(f"flagged={len(report.flagged())}/{len(report.results)} "
          f"without baseline={len(report.without_baseline())} "
          f"no events={len(report.no_events)}")
    # A failed setup still wrote its summary, but the job must go red.
    return 1 if report.failed() else 0


if __name__ == "__main__":
    raise SystemExit(main())
