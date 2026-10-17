#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles distribution-based physics-validation driver.

For each experimental setup: generate events once, bin every measurement, and compare
the PR against the stored main baseline (``p_compat``) and the data (``p_data``).

  --make-baseline   write main's per-setup predictions to --out-dir
  (default)         compare against <baseline-dir>/baselines/<key>/<setup>.json
  --merge           combine shard summaries into one comment
  --unweighter-scan compare Options/Unweighting variants against a reference
  --dry-run         synthetic adapter, no Achilles/NUISANCE
  --selftest        end-to-end check on the synthetic adapter
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Dict, List, Optional

import numpy as np
import yaml

from adapters import GeneratedEvents, Nuisance3Adapter, SyntheticAdapter
from report import (ALPHA, MeasurementResult, MissingMeasurement, Report, ScanReport,
                    VariantMeasurement, VariantSummary, plot_basename)
from stats import (Prediction, bonferroni, compatibility, goodness_of_fit,
                   shape_compatibility, trial_covariance)

# Bump when a stored baseline stops being comparable with fresh predictions.
BASELINE_SCHEMA = 1


# ---------------------------------------------------------------------------
# Predictions and baselines
# ---------------------------------------------------------------------------

def predict_all(adapter, generated: GeneratedEvents, measurements):
    """Bin a setup's events onto all its measurements: (predictions, samples) by name."""
    if not generated.n_nonzero_trials:
        raise RuntimeError("event file carries no trial counts; the MC covariance "
                           "cannot be built")
    samples = adapter.histogram_many(generated, measurements)
    preds = {}
    for m in measurements:
        s = samples[m["name"]]
        preds[m["name"]] = trial_covariance(
            s.bin_index, s.weights, s.nbins, generated.n_nonzero_trials,
            rel_xsec_err=generated.rel_xsec_err or 0.0,
            response=adapter.response_matrix(m, s.nbins))
    return preds, samples


def baseline_path(baseline_dir: str, key: str, experiment: str) -> str:
    return os.path.join(baseline_dir, "baselines", key, f"{experiment}.json")


def make_baseline(adapter, config: dict, *, seed: int, n_events: int, main_sha: str,
                  out_dir: str) -> List[str]:
    """Write one ``<setup>.json`` of main predictions per setup; return the paths."""
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for exp in config["experiments"]:
        gen = adapter.generate(exp, "main", seed, n_events)
        preds, _ = predict_all(adapter, gen, exp["measurements"])
        gen.cleanup()
        path = os.path.join(out_dir, f"{exp['name']}.json")
        with open(path, "w") as fh:
            json.dump({"schema": BASELINE_SCHEMA, "adapter": adapter.name,
                       "main_sha": main_sha, "seed": seed, "events": n_events,
                       "measurements": {k: p.to_dict() for k, p in preds.items()}},
                      fh)
        paths.append(path)
    return paths


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
        key: str, repo: str, feature_sha: str, nuisance_version: str,
        out_dir: str) -> Report:
    from plots import plot_measurement  # only this path needs matplotlib

    os.makedirs(out_dir, exist_ok=True)
    results, missing, main_shas = [], [], set()
    for exp in config["experiments"]:
        try:
            gen = adapter.generate(exp, "feature", seed, n_events)
            features, samples = predict_all(adapter, gen, exp["measurements"])
            gen.cleanup()
        except Exception as err:  # report the setup as failed, keep the rest
            print(f"::error::{exp['name']}: {err}")
            missing += [MissingMeasurement(m["name"], exp["name"], str(err))
                        for m in exp["measurements"]]
            continue

        base = load_baseline(baseline_path(baseline_dir, key, exp["name"]), adapter.name)
        stored = base["measurements"] if base else {}
        if base:
            main_shas.add(base["main_sha"])

        for m in exp["measurements"]:
            name, feature = m["name"], features[m["name"]]
            data = adapter.data_table(m)
            gof_pr = goodness_of_fit(feature, data.values, data.covariance)
            main = Prediction.from_dict(stored[name]) if name in stored else None
            if main is not None and main.nbins != feature.nbins:
                main = None  # binning changed since the baseline
            r = MeasurementResult(name=name, experiment=exp["name"], ndof=gof_pr.ndof,
                                  chi2_ndof_pr=gof_pr.chi2_per_ndof, p_data=gof_pr.pvalue,
                                  plot=plot_basename(name),
                                  selected_events=int(samples[name].bin_index.size))
            if main is not None:
                gof_main = goodness_of_fit(main, data.values, data.covariance)
                r.chi2_ndof_main = gof_main.chi2_per_ndof
                r.delta_chi2 = gof_pr.chi2 - gof_main.chi2
                r.p_compat = compatibility(main, feature).pvalue
            results.append(r)
            plot_measurement(
                os.path.join(out_dir, r.plot), name,
                data=data.values, data_cov=data.covariance,
                main=None if main is None else main.values,
                main_cov=None if main is None else main.covariance,
                feature=feature.values, feature_cov=feature.covariance,
                edges=data.edges, xlabel=data.xlabel, ylabel=data.ylabel,
                subtitle=f"NUISANCE3 {nuisance_version}  |  seed {seed}  |  "
                         f"{n_events:,} events  |  {feature_sha[:8]}")

    report = Report(results=results, repo=repo, feature_sha=feature_sha,
                    nuisance_version=nuisance_version, seed=seed,
                    events_per_measurement=n_events, did_not_run=missing,
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


def merge_shards(shard_paths, out_dir: str, config: Optional[dict] = None) -> Report:
    """Combine shard summaries; ``config`` names the measurements no shard reported."""
    shards = []
    for p in shard_paths:
        if os.path.exists(p):
            with open(p) as fh:
                shards.append(json.load(fh))
    if not shards and not config:
        raise SystemExit("merge: no shard summaries and no config")

    results = [MeasurementResult.from_dict(m) for s in shards for m in s["measurements"]]
    missing = [MissingMeasurement(**m) for s in shards for m in s.get("did_not_run", [])]
    seen = {r.name for r in results} | {m.name for m in missing}
    for exp in (config or {}).get("experiments", []):
        missing += [MissingMeasurement(m["name"], exp["name"], "the shard's job failed")
                    for m in exp["measurements"] if m["name"] not in seen]

    head = shards[0] if shards else {}
    report = Report(results=results,
                    repo=head.get("repo", os.environ.get("GITHUB_REPOSITORY",
                                                         "AchillesGen/Achilles")),
                    feature_sha=head.get("feature_sha",
                                         os.environ.get("GITHUB_SHA", "unknown")),
                    nuisance_version=head.get("nuisance_version", "unknown"),
                    seed=head.get("seed", 0),
                    events_per_measurement=head.get("events_per_measurement", 0),
                    did_not_run=missing,
                    baseline_shas=sorted({sha for s in shards
                                          for sha in s.get("baseline_shas", [])}),
                    image_key=head.get("image_key", ""), run_url=_run_url())
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# Unweighting scan: the same setup generated once per unweighting scheme
# ---------------------------------------------------------------------------

NULL_CONTROL = "null-control"  # the reference reseeded; calibrates the scan


def _mean_error(pred: Prediction) -> float:
    return float(np.mean(np.sqrt(np.clip(np.diag(pred.covariance), 0.0, None))))


def run_unweighting_scan(adapter, config: dict, *, seed: int, n_events: int,
                         repo: str, feature_sha: str, nuisance_version: str,
                         out_dir: str) -> ScanReport:
    """Generate each setup once per variant; every variant must match the reference."""
    from plots import plot_variants

    spec = config.get("unweighting")
    if not spec:
        raise SystemExit("--unweighter-scan needs an 'unweighting:' block in the config")
    variants, reference = spec["variants"], spec["reference"]
    if reference not in [v["name"] for v in variants]:
        raise SystemExit(f"unweighting.reference {reference!r} is not a variant")

    rows: List[VariantMeasurement] = []
    runtime: Dict[str, tuple] = {v["name"]: (0.0, []) for v in variants}
    os.makedirs(out_dir, exist_ok=True)
    for exp in config["experiments"]:
        data = {m["name"]: adapter.data_table(m) for m in exp["measurements"]}
        preds, samples = {}, {}
        for variant in variants:
            name = variant["name"]
            gen = adapter.generate(exp, name, seed, n_events,
                                   unweighting=variant["options"],
                                   seed_offset=int(variant.get("seed_offset", 0)))
            preds[name], samples[name] = predict_all(adapter, gen, exp["measurements"])
            if gen.run is not None:
                secs, effs = runtime[name]
                runtime[name] = (secs + (gen.run.seconds or 0.0),
                                 effs + ([gen.run.unweight_eff]
                                         if gen.run.unweight_eff else []))
            gen.cleanup()

        for m in exp["measurements"]:
            mname = m["name"]
            ref = preds[reference][mname]
            ref_total, ref_err = float(np.sum(ref.values)), _mean_error(ref)
            for variant in variants:
                vname = variant["name"]
                pred, sample = preds[vname][mname], samples[vname][mname]
                is_ref = vname == reference
                gof = goodness_of_fit(pred, data[mname].values, data[mname].covariance)
                total = float(np.sum(pred.values))
                rows.append(VariantMeasurement(
                    measurement=mname, variant=vname, ndof=pred.nbins,
                    chi2_ndof_data=gof.chi2_per_ndof,
                    p_compat=float("nan") if is_ref else compatibility(ref, pred).pvalue,
                    p_shape=(float("nan") if is_ref
                             else shape_compatibility(ref, pred).pvalue),
                    p_data=gof.pvalue,
                    norm_shift=(total - ref_total) / ref_total if ref_total else float("nan"),
                    ess_fraction=sample.ess_fraction(),
                    max_over_mean=sample.max_over_mean(),
                    mc_error_ratio=_mean_error(pred) / ref_err if ref_err else float("nan")))
            plot_variants(
                os.path.join(out_dir, f"{mname}.unweighting.png"), mname,
                data=data[mname].values, data_cov=data[mname].covariance,
                variants={v["name"]: (preds[v["name"]][mname].values,
                                      preds[v["name"]][mname].covariance)
                          for v in variants},
                reference=reference, edges=data[mname].edges,
                xlabel=data[mname].xlabel, ylabel=data[mname].ylabel,
                subtitle=f"NUISANCE3 {nuisance_version}  |  seed {seed}  |  "
                         f"{n_events:,} events/variant  |  {feature_sha[:8]}")

    report = ScanReport(
        summaries=_scan_summaries([(v["name"], v["options"]) for v in variants], rows,
                                  reference, runtime),
        rows=rows, reference=reference, repo=repo, feature_sha=feature_sha,
        nuisance_version=nuisance_version, seed=seed, events_per_measurement=n_events)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


def _scan_summaries(variants, rows: List[VariantMeasurement], reference: str,
                    runtime: Dict[str, tuple]) -> List[VariantSummary]:
    """One line per variant; ``runtime`` maps a variant to (seconds, [efficiencies])."""
    summaries = []
    for name, options in variants:
        mine = [r for r in rows if r.variant == name]
        pv = [r.p_compat for r in mine if r.p_compat == r.p_compat]
        ps = [r.p_shape for r in mine if r.p_shape == r.p_shape]
        seconds, effs = runtime.get(name, (0.0, []))
        nan = float("nan")
        summaries.append(VariantSummary(
            variant=name, options=dict(options), n_measurements=len(mine),
            p_worst=min(pv, default=nan), p_overall=bonferroni(pv),
            p_shape_worst=min(ps, default=nan), p_shape_overall=bonferroni(ps),
            n_flagged=sum(p < ALPHA for p in pv),
            ess_fraction=float(np.nanmedian([r.ess_fraction for r in mine])) if mine else nan,
            mc_error_ratio=(float(np.nanmedian([r.mc_error_ratio for r in mine]))
                            if mine else nan),
            max_norm_shift=max((r.norm_shift for r in mine), key=abs, default=nan),
            seconds=seconds or None, unweight_eff=min(effs, default=None),
            is_reference=name == reference, is_null_control=name == NULL_CONTROL))
    return summaries


def merge_scan_shards(shard_paths, out_dir: str) -> ScanReport:
    shards = []
    for p in shard_paths:
        with open(p) as fh:
            shards.append(json.load(fh))
    if not shards or any(s.get("kind") != "unweighting-scan" for s in shards):
        raise SystemExit("merge: expected unweighting-scan shard summaries")
    head = shards[0]
    rows = [VariantMeasurement(**r) for s in shards for r in s["rows"]]
    variants = [(v["variant"], v["options"]) for v in head["variants"]]
    runtime: Dict[str, tuple] = {name: (0.0, []) for name, _ in variants}
    for s in shards:
        for v in s["variants"]:
            secs, effs = runtime[v["variant"]]
            runtime[v["variant"]] = (secs + (v["seconds"] or 0.0),
                                     effs + ([v["unweight_eff"]] if v["unweight_eff"] else []))
    report = ScanReport(
        summaries=_scan_summaries(variants, rows, head["reference"], runtime),
        rows=rows, reference=head["reference"], repo=head["repo"],
        feature_sha=head["feature_sha"], nuisance_version=head["nuisance_version"],
        seed=head["seed"], events_per_measurement=head["events_per_measurement"])
    os.makedirs(out_dir, exist_ok=True)
    report.write(os.path.join(out_dir, "comment.md"),
                 os.path.join(out_dir, "summary.json"))
    return report


# ---------------------------------------------------------------------------
# Config and CLI
# ---------------------------------------------------------------------------

# Top-level per-sample maps in the config, folded onto the measurement they name.
_MEASUREMENT_KEYS = ("data_scale", "bin_edges", "bin_widths", "solid_angle",
                     "smearing")


def load_config(path: str) -> dict:
    with open(path) as fh:
        config = yaml.safe_load(fh)
    for exp in config["experiments"]:
        exp["measurements"] = [m if isinstance(m, dict) else {"name": m}
                               for m in exp["measurements"]]
    by_name = {m["name"]: m for exp in config["experiments"] for m in exp["measurements"]}
    for key in _MEASUREMENT_KEYS:
        entries = config.get(key) or {}
        unknown = set(entries) - set(by_name)
        if unknown:
            raise SystemExit(f"{key}: no such measurement(s) {sorted(unknown)}")
        for name, value in entries.items():
            if key == "smearing":  # a csv relative to the config
                value = os.path.join(os.path.dirname(os.path.abspath(path)), value)
            by_name[name][key] = value
    return config


def filter_config(config: dict, experiments=(), measurements=(), variants=()) -> dict:
    """Narrow the config to named setups, measurements and scan variants."""
    exps = config["experiments"]
    if experiments:
        exps = [e for e in exps if e["name"] in experiments]
        missing = set(experiments) - {e["name"] for e in exps}
        if missing:
            raise SystemExit(f"setups not in config: {sorted(missing)}")
    if measurements:
        exps = [{**e, "measurements": [m for m in e["measurements"]
                                       if m["name"] in measurements]} for e in exps]
        exps = [e for e in exps if e["measurements"]]
        missing = set(measurements) - {m["name"] for e in exps for m in e["measurements"]}
        if missing:
            raise SystemExit(f"measurements not in config: {sorted(missing)}")
    config = {**config, "experiments": exps}
    if variants:
        spec = config.get("unweighting") or {}
        wanted = set(variants) | {spec.get("reference")}
        kept = [v for v in spec.get("variants", []) if v["name"] in wanted]
        missing = wanted - {v["name"] for v in kept}
        if missing:
            raise SystemExit(f"variants not in config: {sorted(missing)}")
        config["unweighting"] = {**spec, "variants": kept}
    return config


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="measurements.yml")
    p.add_argument("--workdir", default="physval-work")
    p.add_argument("--out-dir", default="physval-out")
    p.add_argument("--baseline-dir", default="physval-baselines",
                   help="checkout of the physval-baselines branch")
    p.add_argument("--key", default="dev", help="baseline key (the image digest in CI)")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--events", type=int, default=200000)
    p.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY",
                                                    "AchillesGen/Achilles"))
    p.add_argument("--feature-sha", default=os.environ.get("GITHUB_SHA", "local"))
    p.add_argument("--nuisance-version",
                   default=os.environ.get("NUISANCE_VERSION", "unknown"))
    p.add_argument("--only-experiment", action="append", default=[],
                   dest="only_experiments", help="run only this setup (repeatable)")
    p.add_argument("--only", action="append", default=[], dest="only_measurements",
                   help="run only this measurement (repeatable)")
    p.add_argument("--only-variant", action="append", default=[], dest="only_variants",
                   help="run only this scan variant (repeatable)")
    p.add_argument("--unweighter-scan", action="store_true")
    p.add_argument("--merge", nargs="+", default=None,
                   help="merge these shard summary.json files")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--make-baseline", action="store_true")
    p.add_argument("--selftest", action="store_true")
    return p


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)
    if args.selftest:
        return _selftest()

    config = filter_config(load_config(args.config), args.only_experiments,
                           args.only_measurements, args.only_variants)
    if args.merge:
        if args.unweighter_scan:
            scan = merge_scan_shards(args.merge, out_dir=args.out_dir)
            print(f"merged {len(args.merge)} scan shard(s); "
                  f"biased={[s.variant for s in scan.biased()] or 'none'}")
            return 0
        report = merge_shards(args.merge, out_dir=args.out_dir, config=config)
        print(f"merged {len(report.results)} measurement(s); "
              f"flagged={len(report.flagged())}; did not report="
              f"{', '.join(report.failed_setups()) or 'none'}")
        return 0

    adapter = (SyntheticAdapter(base_seed=args.seed) if args.dry_run
               else Nuisance3Adapter(workdir=args.workdir))
    common = dict(seed=args.seed, n_events=args.events, out_dir=args.out_dir)
    if args.unweighter_scan:
        scan = run_unweighting_scan(adapter, config, repo=args.repo,
                                    feature_sha=args.feature_sha,
                                    nuisance_version=args.nuisance_version, **common)
        print(f"reference={scan.reference} "
              f"biased={[s.variant for s in scan.biased()] or 'none'}")
        return 0
    if args.make_baseline:
        for path in make_baseline(adapter, config, main_sha=args.feature_sha, **common):
            print(f"wrote baseline: {path}")
        return 0

    report = run(adapter, config, baseline_dir=args.baseline_dir, key=args.key,
                 repo=args.repo, feature_sha=args.feature_sha,
                 nuisance_version=args.nuisance_version, **common)
    print(f"flagged={len(report.flagged())}/{len(report.results)} "
          f"without baseline={len(report.without_baseline())}")
    # A failed setup still wrote its summary, but the job must go red.
    return 1 if report.did_not_run else 0


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

class _FailingAdapter(SyntheticAdapter):
    def generate(self, experiment, *a, **kw):
        if experiment["name"] == "SYNTH_crash":
            raise RuntimeError("achilles failed for SYNTH_crash (signal 11)")
        return super().generate(experiment, *a, **kw)


def _selftest() -> int:
    import tempfile

    config = {"experiments": [
        {"name": "SYNTH_stable", "dryrun": {"feature_shift": 0.0},
         "measurements": [{"name": "SYNTH_stable_1Dx", "dryrun": {"nbins": 12}}]},
        {"name": "SYNTH_regressed", "dryrun": {"feature_shift": 0.05},
         "measurements": [{"name": "SYNTH_regressed_1Dx", "dryrun": {"nbins": 12}}]},
    ]}
    crash = {"name": "SYNTH_crash", "measurements": [{"name": "SYNTH_crash_1Dx"}]}
    with tempfile.TemporaryDirectory() as tmp:
        kw = dict(seed=7, n_events=60000, repo="AchillesGen/Achilles",
                  feature_sha="deadbeefcafef00d", nuisance_version="selftest")
        make_baseline(SyntheticAdapter(base_seed=7), config, seed=7, n_events=60000,
                      main_sha="0123456789ab",
                      out_dir=os.path.join(tmp, "baselines", "k"))
        rep = run(_FailingAdapter(base_seed=7),
                  {"experiments": config["experiments"] + [crash]},
                  baseline_dir=tmp, key="k", out_dir=os.path.join(tmp, "a"), **kw)
        status = {r.name: r.status() for r in rep.results}

        # A baseline from another adapter must not be used.
        with open(baseline_path(tmp, "k", "SYNTH_stable")) as fh:
            b = json.load(fh)
        b["adapter"] = "nuisance3"
        with open(baseline_path(tmp, "k", "SYNTH_stable"), "w") as fh:
            json.dump(b, fh)
        other = run(SyntheticAdapter(base_seed=7), config, baseline_dir=tmp, key="k",
                    out_dir=os.path.join(tmp, "b"), **kw)

        merged = merge_shards([os.path.join(tmp, "a", "summary.json"),
                               os.path.join(tmp, "nope", "summary.json")],
                              os.path.join(tmp, "m"),
                              config={"experiments": config["experiments"] + [crash]})
        checks = {
            "stable is compatible": status["SYNTH_stable_1Dx"] == "compatible",
            "regression is flagged": status["SYNTH_regressed_1Dx"] in
                                     ("regression", "changed", "improvement"),
            "baseline sha reported": rep.baseline_shas == ["0123456789ab"],
            "crash reported with reason": [m.reason for m in rep.did_not_run] ==
                                          ["achilles failed for SYNTH_crash (signal 11)"],
            "other adapter's baseline ignored":
                {r.name: r.status() for r in other.results}["SYNTH_stable_1Dx"]
                == "no-baseline",
            "merge keeps the reason": "(signal 11)" in merged.to_markdown(),
            "merge counts once": len(merged.did_not_run) == 1,
            "plot written": os.path.exists(os.path.join(tmp, "a",
                                                        "SYNTH_stable_1Dx.png")),
        }
        checks.update(_selftest_scan(tmp))
    for name, ok in checks.items():
        print(f"[{'ok' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def _selftest_scan(tmp: str) -> dict:
    """Unweighting on synthetic events: unbiased, sharper per event, slower."""
    config = {
        "unweighting": {"reference": "weighted", "variants": [
            {"name": "weighted", "options": {"Name": "None"}},
            {"name": "percentile-99", "options": {"Name": "Percentile", "percentile": 99}},
            {"name": "excess-1e-2", "options": {"Name": "Excess", "epsilon": 0.01}},
            {"name": "tailfrac-1e-2", "options": {"Name": "TailFraction", "epsilon": 0.01}},
        ]},
        "experiments": [{"name": "SYNTH_exp", "dryrun": {"feature_shift": 0.0},
                         "measurements": [{"name": "SYNTH_scan", "dryrun": {"nbins": 12}}]}],
    }
    out = os.path.join(tmp, "scan")
    report = run_unweighting_scan(SyntheticAdapter(base_seed=3), config, seed=3,
                                  n_events=40000, repo="AchillesGen/Achilles",
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
        "scan: overlay plot written": os.path.exists(
            os.path.join(out, "SYNTH_scan.unweighting.png")),
    }


if __name__ == "__main__":
    raise SystemExit(main())
