#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles distribution-based physics-validation driver.

Ties the pieces together for one physval run:

    config (measurements.yml)
        -> adapter.generate(feature)  -> bootstrap covariance  -> feature Prediction
        -> baseline (stored main Prediction + data)            -> main Prediction
        -> compatibility p_compat + goodness-of-fit vs data
        -> Sherpa-style PR comment (comment.md) + summary.json

Modes:
    --dry-run           use the SyntheticAdapter (no Achilles/NUISANCE) end to end
    --make-baseline     produce the stored 'main' baseline JSON (predictions + data
                        + bootstrap covariance) for the physval-baselines branch
    --selftest          dry-run with an injected regression; assert it is flagged

The real path (NUISANCE3) is selected with --adapter nuisance3 once the container
and sample mapping are wired up in adapters.Nuisance3Adapter.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Dict, Optional

import numpy as np
import yaml

from adapters import Adapter, DataTable, make_adapter
from report import MeasurementResult, Report
from stats import (Prediction, bootstrap_covariance, compatibility,
                   goodness_of_fit)


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


# ---------------------------------------------------------------------------
# Core: build a Prediction for one branch/measurement
# ---------------------------------------------------------------------------

def predict(adapter: Adapter, measurement: dict, branch: str, seed: int,
            n_events: int, n_boot: int,
            rng: np.random.Generator) -> Prediction:
    sample = adapter.generate(measurement, branch=branch, seed=seed,
                              n_events=n_events)
    return bootstrap_covariance(sample.bin_index, sample.weights, sample.nbins,
                                n_boot=n_boot, rng=rng)


# ---------------------------------------------------------------------------
# make-baseline: compute and store 'main' predictions + data + covariance
# ---------------------------------------------------------------------------

def make_baseline(adapter: Adapter, config: dict, key: str, seed: int,
                  n_events: int, n_boot: int, out_dir: str) -> str:
    rng = np.random.default_rng(seed)
    measurements: Dict[str, dict] = {}
    for m in config["measurements"]:
        name = m["name"]
        main = predict(adapter, m, "main", seed, n_events, n_boot, rng)
        data = adapter.data_table(m)
        measurements[name] = {
            "prediction": main.to_dict(),
            "data": {"values": data.values.tolist(),
                     "covariance": data.covariance.tolist()},
        }
    baseline = {
        "key": key,
        "seed": seed,
        "events_per_measurement": n_events,
        "n_boot": n_boot,
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

def run(adapter: Adapter, config: dict, *, seed: int, n_events: int, n_boot: int,
        baseline: Optional[dict], repo: str, feature_sha: str,
        nuisance_version: str, data_release: str,
        out_dir: str) -> Report:
    rng = np.random.default_rng(seed + 101)
    results = []
    warnings = []

    for m in config["measurements"]:
        name = m["name"]
        feature = predict(adapter, m, "feature", seed, n_events, n_boot, rng)

        if baseline is not None and name in baseline.get("measurements", {}):
            main, data = _measurement_baseline(baseline, name)
        else:
            # No stored baseline (e.g. dry-run, or a new measurement): compute main
            # inline so the run still produces a complete table.
            warnings.append(name)
            main = predict(adapter, m, "main", seed, n_events, n_boot, rng)
            data = adapter.data_table(m)

        compat = compatibility(main, feature)
        gof_pr = goodness_of_fit(feature, data.values, data.covariance)
        gof_main = goodness_of_fit(main, data.values, data.covariance)

        results.append(MeasurementResult(
            name=name,
            ndof=compat.ndof,
            chi2_ndof_main=gof_main.chi2_per_ndof,
            chi2_ndof_pr=gof_pr.chi2_per_ndof,
            delta_chi2=gof_pr.chi2 - gof_main.chi2,
            p_compat=compat.pvalue,
            p_data=gof_pr.pvalue,
            plot=f"{name}.png",
        ))

    extra_header = []
    if warnings:
        extra_header.append(
            f"> ⚠️ No stored baseline for {len(warnings)} measurement(s) "
            f"({', '.join(warnings)}); main was computed inline for this run.")

    report = Report(results=results, repo=repo, feature_sha=feature_sha,
                    nuisance_version=nuisance_version, data_release=data_release,
                    seed=seed, events_per_measurement=n_events,
                    extra_header=extra_header)

    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_config(path: str) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def _filter_measurements(config: dict, only) -> dict:
    if not only:
        return config
    wanted = set(only)
    kept = [m for m in config["measurements"] if m["name"] in wanted]
    missing = wanted - {m["name"] for m in kept}
    if missing:
        raise SystemExit(f"--only names not in config: {sorted(missing)}")
    return {**config, "measurements": kept}


def merge_shards(shard_paths, out_dir: str) -> Report:
    """Combine per-measurement shard ``summary.json`` files into one Report.

    Each shard is a summary.json emitted by a sharded ``run`` (typically one
    measurement). Run-level metadata is taken from the first shard.
    """
    shards = []
    for p in shard_paths:
        with open(p) as fh:
            shards.append(json.load(fh))
    if not shards:
        raise SystemExit("merge_shards: no shard summaries given")

    head = shards[0]
    results = []
    for s in shards:
        for m in s["measurements"]:
            results.append(MeasurementResult(
                name=m["name"], ndof=m["ndof"],
                chi2_ndof_main=m["chi2_ndof_main"], chi2_ndof_pr=m["chi2_ndof_pr"],
                delta_chi2=m["delta_chi2"], p_compat=m["p_compat"],
                p_data=m["p_data"], plot=m.get("plot")))

    report = Report(results=results, repo=head["repo"],
                    feature_sha=head["feature_sha"],
                    nuisance_version=head["nuisance_version"],
                    data_release=head["data_release"], seed=head["seed"],
                    events_per_measurement=head["events_per_measurement"])
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="measurements.yml",
                   help="measurement-list YAML")
    p.add_argument("--adapter", default="synthetic",
                   choices=["synthetic", "nuisance3"])
    p.add_argument("--nuisance-image", default=os.environ.get("NUISANCE_IMAGE", ""),
                   help="NUISANCE3 container image (nuisance3 adapter)")
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
    p.add_argument("--data-release", default=os.environ.get("PHYSVAL_DATA_RELEASE",
                                                            "unknown"))
    p.add_argument("--only", action="append", default=[],
                   help="run only this measurement (repeatable); used for sharding")
    p.add_argument("--merge", nargs="+", default=None,
                   help="merge these shard summary.json files into one report")
    p.add_argument("--dry-run", action="store_true",
                   help="force the synthetic adapter")
    p.add_argument("--make-baseline", action="store_true",
                   help="produce the stored main baseline instead of a comparison")
    p.add_argument("--selftest", action="store_true")
    return p


def _adapter_from_args(args) -> Adapter:
    kind = "synthetic" if args.dry_run else args.adapter
    return make_adapter(kind, base_seed=args.seed, workdir=args.workdir,
                        nuisance_image=args.nuisance_image)


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    if args.selftest:
        return _selftest()

    if args.merge:
        report = merge_shards(args.merge, out_dir=args.out_dir)
        print(f"merged {len(args.merge)} shard(s): p_overall={report.p_overall():.4g} "
              f"flagged={report.n_flagged()}/{len(report.results)}")
        print(f"wrote {args.out_dir}/comment.md and {args.out_dir}/summary.json")
        return 0

    config = _filter_measurements(_load_config(args.config), args.only)
    adapter = _adapter_from_args(args)

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
                 data_release=args.data_release, out_dir=args.out_dir)

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

    config = {"measurements": [
        {"name": "SYNTH_stable", "dryrun": {"nbins": 12, "feature_shift": 0.0}},
        {"name": "SYNTH_regressed", "dryrun": {"nbins": 12, "feature_shift": 0.05}},
    ]}
    with tempfile.TemporaryDirectory() as tmp:
        adapter = make_adapter("synthetic", base_seed=7)
        bpath = make_baseline(adapter, config, key="selftest", seed=7,
                              n_events=60000, n_boot=150, out_dir=tmp)
        baseline = load_baseline(bpath)
        report = run(adapter, config, seed=7, n_events=60000, n_boot=150,
                     baseline=baseline, repo="AchillesGen/Achilles",
                     feature_sha="deadbeefcafef00d", nuisance_version="selftest",
                     data_release="selftest", out_dir=tmp)
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
        for name, ok in checks.items():
            print(f"[{'ok' if ok else 'FAIL'}] {name}")
        ok = all(checks.values())
        print("SELFTEST:", "PASS" if ok else "FAIL")
        return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
