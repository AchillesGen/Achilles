# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the physval pipeline; run with pytest. No Achilles or NUISANCE needed."""

import json
import os
import types

import numpy as np
import pytest
from scipy import stats as scipy_stats

import ci_scope
import physval
from adapters import EventSample, GeneratedEvents, Nuisance3Adapter, SyntheticAdapter
from plots import plot_measurement
from report import MAX_COMMENT_CHARS, MeasurementResult, MissingMeasurement, Report, plot_basename
from stats import benjamini_hochberg, compatibility, trial_covariance

HERE = os.path.dirname(os.path.abspath(__file__))


# -- stats --------------------------------------------------------------------

def _runs(rng, n_runs, shift=0.0, nbins=25, n_nonzero=50_000, acceptance=0.2,
          rel_xsec_err=0.0015):
    """Runs of one configuration as a generator delivers them: N non-zero trials fixed,
    Binomial(N, acceptance) selected, a per-run xsec error and uneven bin widths."""
    shape = np.linspace(2.0, 1.0, nbins) + shift * np.linspace(-1.0, 1.0, nbins)
    probs = np.append(shape / shape.sum() * acceptance, 1.0 - acceptance)
    widths = np.where(np.arange(nbins) < nbins // 2, 0.1, 1.0)
    for _ in range(n_runs):
        counts = rng.multinomial(n_nonzero, probs)[:nbins]
        bins = np.repeat(np.arange(nbins), counts)
        weights = (1.0 + rel_xsec_err * rng.normal()) / n_nonzero / widths[bins]
        yield trial_covariance(bins, weights, nbins, n_nonzero, rel_xsec_err=rel_xsec_err)


def test_two_runs_of_one_configuration_give_uniform_p_values():
    rng = np.random.default_rng(98765)
    runs = list(_runs(rng, 800))
    pvals = np.array([compatibility(runs[2 * i], runs[2 * i + 1]).pvalue for i in range(400)])
    assert abs(np.mean(pvals < 0.05) - 0.05) <= 3 * np.sqrt(0.05 * 0.95 / 400)
    assert scipy_stats.kstest(pvals, "uniform").pvalue > 0.01


def test_a_shape_change_is_detected():
    rng = np.random.default_rng(1)
    assert compatibility(next(_runs(rng, 1)), next(_runs(rng, 1, shift=0.5))).pvalue < 0.05


def test_response_is_applied_to_values_and_covariance():
    bins, weights, response = np.array([0, 0, 1]), np.ones(3), np.array([[0.0, 2.0], [1.0, 0.0]])
    plain = trial_covariance(bins, weights, 2, 10)
    smeared = trial_covariance(bins, weights, 2, 10, response=response)
    assert np.allclose(smeared.values, response @ plain.values)
    assert np.allclose(smeared.covariance, response @ plain.covariance @ response.T)


def test_benjamini_hochberg():
    assert np.allclose(benjamini_hochberg([0.01, 0.04, 0.03, 0.5]),
                       [0.04, 0.16 / 3, 0.16 / 3, 0.5])


# -- report -------------------------------------------------------------------

def _row(name, exp, p, delta=2.0):
    return MeasurementResult(name, exp, ndof=12, chi2_ndof_pr=1.2, p_data=0.4,
                             chi2_ndof_main=1.1, delta_chi2=delta, p_compat=p)


def test_statuses_and_flags():
    rep = Report([
        _row("quiet", "A", 0.6),
        _row("nominal_only", "A", 0.04),  # p < 0.05, but not after BH
        _row("worse", "B", 1e-4, delta=6.1),
        _row("better", "B", 2e-4, delta=-11.0),
        _row("tiny", "B", 3e-4, delta=0.2),
        MeasurementResult("new", "C", 8, 1.5, 0.3),
    ])
    assert {r.name: r.status() for r in rep.results} == {
        "quiet": "compatible", "nominal_only": "compatible", "worse": "regression",
        "better": "improvement", "tiny": "changed", "new": "no-baseline"}
    assert len(rep.flagged()) == 3 and not rep.failed()
    assert "compared with data only" in rep.to_markdown()
    md = rep.to_markdown()
    assert md.index("<b>B</b>") < md.index("<b>A</b>")  # setups needing attention first


def test_a_crash_shows_only_the_first_line_of_its_reason():
    rep = Report([], did_not_run=[MissingMeasurement("x", "A", "signal 11\nlog tail")])
    md = rep.to_markdown()
    assert rep.failed() and "incomplete run" in md
    assert "signal 11" in md and "log tail" not in md


def test_an_unexpected_empty_sample_fails_and_a_flagged_one_warns():
    bad = MissingMeasurement("bad", "A")
    rare = MissingMeasurement("rare", "A", "tail of the spectrum")
    md = Report([], no_events=[rare, bad]).to_markdown()
    assert Report([], no_events=[bad]).failed() and "incomplete run" in md
    assert "No events selected (2)" in md and "**1** with no events selected" in md
    assert "may be empty: tail of the spectrum" in md
    assert not Report([], no_events=[rare]).failed()


def test_a_huge_suite_is_trimmed_to_fit_a_comment():
    rows = [_row(f"Exp{s}_XSec_1DVar{i}_nu", f"Exp{s}", 0.5) for s in range(40) for i in range(15)]
    md = Report(rows, feature_sha="abc").to_markdown()
    assert len(md) <= MAX_COMMENT_CHARS and "Trimmed" in md


def test_plot_names_are_safe_for_artifacts():
    assert plot_basename("A_Barreau:1983ht") == "A_Barreau_1983ht.png"


# -- plots --------------------------------------------------------------------

def test_plots_survive_zero_data_bins_and_no_baseline(tmp_path):
    values, cov = np.array([1.0, 0.0, 2.0]), np.diag([0.1, 0.1, 0.1])
    for name, extra in (("with_main", dict(main=values, main_cov=cov, edges=[0, 1, 2, 4])),
                        ("no_main", {})):
        path = tmp_path / f"{name}.png"
        plot_measurement(str(path), name, data=values, data_cov=cov, feature=values * 1.1,
                         feature_cov=cov, subtitle="seed 1", **extra)
        assert path.read_bytes().startswith(b"\x89PNG")


# -- adapter ------------------------------------------------------------------

def _fake_nuisance_adapter(blocks):
    """A Nuisance3Adapter on a fake pyNUISANCE serving ``blocks``, whose columns are
    [weight, fatx, then (selection, projection) per sample]; bin edges are 0, 1, 2, 3."""
    class Binning:
        def find_bins(self, frame, columns):
            values = frame.table[:, frame.column_names.index(columns[0])]
            bins = np.digitize(values, [0.0, 1.0, 2.0, 3.0]) - 1
            return np.where((bins >= 0) & (bins < 3), bins, 2 ** 31 - 1)

        def bin_sizes(self):
            return [1.0, 1.0, 1.0]

    class FrameGen:
        def __init__(self):
            self.names, self.i = ["weight.cv", "fatx_per_sumw.pb_per_target.estimate"], 0
        def add_int_column(self, name, op): self.names.append(name)
        add_double_column = add_int_column
        def first(self): return self.next()
        def next(self):
            table = blocks[self.i] if self.i < len(blocks) else np.empty((0, 6))
            self.i += 1
            return types.SimpleNamespace(table=table, column_names=self.names)

    op = types.SimpleNamespace(op=None)
    binned = types.SimpleNamespace(binning=Binning(), values=np.zeros(3))
    analysis = types.SimpleNamespace(get_selection=lambda: op, get_projections=lambda: [op],
                                     get_data=lambda: [binned],
                                     convert_event_rates=lambda rates, scale: scale * rates[::-1])
    adapter = Nuisance3Adapter.__new__(Nuisance3Adapter)
    adapter._pn = lambda: types.SimpleNamespace(EventSource=lambda path: object(),
                                                EventFrameGen=lambda evs, n: FrameGen())
    adapter._analysis = lambda name: analysis
    return adapter


def test_each_sample_takes_its_own_rows_scaled_by_the_last_fatx():
    blocks = [np.array([[2.0, 10.0, 1, 0.5, 1, 2.5],   # both select, different bins
                        [3.0, 10.0, 0, 9.9, 1, 1.5],   # only B selects
                        [1.0, 10.0, 1, 7.0, 0, 0.5]]),  # A selects, out of range
              np.array([[4.0, 20.0, 1, 2.5, 0, 0.5]])]  # fatx has moved
    got = _fake_nuisance_adapter(blocks).histogram_many(
        GeneratedEvents(path="x", xsec_divisor=2.0), [{"name": "A"}, {"name": "B"}])
    scale = 20.0 * Nuisance3Adapter._PB_TO_CM2 / 2.0
    assert list(got["A"].bin_index) == [0, 2] and list(got["B"].bin_index) == [2, 1]
    assert np.allclose(got["A"].weights, np.array([2.0, 4.0]) * scale)
    assert np.allclose(got["B"].weights, np.array([2.0, 3.0]) * scale)


def test_response_probes_the_conversion_unless_a_solid_angle_replaces_it():
    adapter = _fake_nuisance_adapter([])
    assert np.allclose(adapter.response_matrix({"name": "A"}, 3), np.eye(3)[::-1])
    assert np.allclose(adapter.response_matrix({"name": "A", "solid_angle": 4.0}, 3),
                       np.eye(3) / 4.0)


# -- driver -------------------------------------------------------------------

class _Adapter(SyntheticAdapter):
    """Synthetic, except for the setup named ``crash`` and samples named ``empty*``."""

    def generate(self, experiment, *args):
        if experiment["name"] == "crash":
            raise RuntimeError("achilles failed (signal 11)")
        return super().generate(experiment, *args)

    def histogram_many(self, generated, measurements):
        samples = super().histogram_many(generated, measurements)
        for name, s in samples.items():
            if name.startswith("empty"):
                samples[name] = EventSample(s.bin_index[:0], s.weights[:0], s.nbins)
        return samples


def _setup(name, *samples, **kw):
    return {"name": name, "measurements": [{"name": s} for s in samples], **kw}


COMPARE = dict(seed=7, n_events=60000, feature_sha="deadbeef", nuisance_version="test", key="k")


def _baseline(config, tmp_path):
    return physval.make_baseline(_Adapter(7), config, seed=7, n_events=60000, main_sha="0123abcd",
                                 out_dir=str(tmp_path / "baselines" / "k"))


def _compare(config, tmp_path, adapter=None, out="out"):
    return physval.run(adapter or _Adapter(7), config, baseline_dir=str(tmp_path),
                       out_dir=str(tmp_path / out), **COMPARE)


def test_compare_flags_a_real_shift_and_writes_plots(tmp_path):
    config = {"experiments": [_setup("stable", "stable_x"),
                              _setup("shifted", "shifted_x", feature_shift=0.05)]}
    _baseline(config, tmp_path)
    rep = _compare(config, tmp_path)
    status = {r.name: r.status() for r in rep.results}
    assert status["stable_x"] == "compatible"
    assert status["shifted_x"] in ("regression", "changed", "improvement")
    assert rep.baseline_shas == ["0123abcd"]
    assert (tmp_path / "out" / "stable_x.png").exists()


def test_another_adapters_baseline_is_ignored(tmp_path):
    config = {"experiments": [_setup("stable", "stable_x")]}
    _baseline(config, tmp_path)
    path = physval.baseline_path(str(tmp_path), "k", "stable")
    stored = json.load(open(path))
    json.dump({**stored, "adapter": "nuisance3"}, open(path, "w"))
    assert _compare(config, tmp_path).results[0].status() == "no-baseline"


def test_no_events_is_reported_without_a_chi2_or_plot(tmp_path):
    exp = _setup("A", "fine", "empty_bad", "empty_rare")
    exp["measurements"][2]["may_be_empty"] = "tail of the spectrum"
    rep = _compare({"experiments": [exp]}, tmp_path)
    assert [r.name for r in rep.results] == ["fine"]
    assert {m.name: m.reason for m in rep.no_events} == {
        "empty_bad": "", "empty_rare": "tail of the spectrum"}
    assert rep.failed() and not (tmp_path / "out" / "empty_bad.png").exists()


def test_baseline_leaves_out_empty_samples_and_reports_them(tmp_path):
    empty = _baseline({"experiments": [_setup("A", "fine", "empty_bad")]}, tmp_path)
    assert [m.name for m in empty] == ["empty_bad"]
    stored = json.load(open(physval.baseline_path(str(tmp_path), "k", "A")))
    assert list(stored["measurements"]) == ["fine"]


def test_merge_keeps_every_outcome_once(tmp_path):
    config = {"experiments": [_setup("A", "fine", "empty_bad"), _setup("crash", "crash_x"),
                              _setup("lost", "lost_x")]}
    shard = {"experiments": config["experiments"][:2]}
    _compare(shard, tmp_path)
    merged = physval.merge_shards([str(tmp_path / "out" / "summary.json"), "nope/summary.json"],
                                  str(tmp_path / "merged"), config)
    assert [r.name for r in merged.results] == ["fine"]
    assert [m.name for m in merged.no_events] == ["empty_bad"]
    assert {m.name: m.reason for m in merged.did_not_run} == {
        "crash_x": "achilles failed (signal 11)", "lost_x": "the shard's job failed"}
    assert merged.feature_sha == "deadbeef"


def test_the_real_config_is_consistent():
    config = physval.load_config(os.path.join(HERE, "measurements.yml"))
    for exp in config["experiments"]:
        assert float(exp["xsec_divisor"]) > 0, exp["name"]
        assert os.path.isfile(os.path.join(HERE, "..", "..", exp["achilles_run"])), exp["name"]
    assert physval.load_config(os.path.join(HERE, "measurements.yml"), ["T2K_CC"])[
        "experiments"][0]["name"] == "T2K_CC"
    with pytest.raises(SystemExit):
        physval.load_config(os.path.join(HERE, "measurements.yml"), ["nope"])


# -- CI scope -----------------------------------------------------------------

@pytest.mark.parametrize("messages, expected", [
    (["fix: x !physval"], (True, False, [])),
    (["fix: x\n\n!physval(A)"], (False, False, [])),
    (["x !physval(A, dry-run)", "y !physval(B,A)"], (True, True, ["A", "B"])),
])
def test_commit_markers(messages, expected):
    assert ci_scope.parse_markers(messages) == expected
