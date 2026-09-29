#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Render the physval summary as a PR comment (markdown) and summary.json.

A row is flagged when its main-vs-PR compatibility survives a Benjamini-Hochberg
correction across the suite; the change in chi-square against data only labels the
direction. No plotting or NUISANCE dependency, so it is testable on synthetic rows.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields
from typing import List, Optional

from stats import benjamini_hochberg, bonferroni

ALPHA = 0.05
# Below this |Δχ²| a flagged row is called "changed" rather than better or worse.
DELTA_CHI2_MIN = 1.0
COMMENT_MARKER = "<!-- achilles-physval-summary -->"
RAW_URL_TEMPLATE = ("https://raw.githubusercontent.com/{repo}/physval-baselines/"
                    "plots/{sha}/{plot}")
MAX_INLINE_PLOTS = 6
# GitHub rejects comments over 65536 characters; past this, quiet setups collapse.
MAX_COMMENT_CHARS = 60000

# Characters upload-artifact and git trees reject; Durham sample names contain ':'.
_UNSAFE_IN_FILENAMES = '":<>|*?\r\n'


def plot_basename(measurement: str) -> str:
    return "".join("_" if c in _UNSAFE_IN_FILENAMES else c for c in measurement) + ".png"


def _fmt(value: Optional[float], spec: str = ".3g") -> str:
    return "—" if value is None or value != value else format(value, spec)


def _plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


@dataclass
class MeasurementResult:
    name: str
    experiment: str
    ndof: int
    chi2_ndof_pr: float                    # PR vs data
    p_data: float                          # PR vs data
    chi2_ndof_main: Optional[float] = None  # None: no usable baseline
    delta_chi2: Optional[float] = None      # chi2_pr - chi2_main vs data; + = worse
    p_compat: Optional[float] = None        # main vs PR
    q_compat: Optional[float] = None        # BH-adjusted p_compat, set by Report
    plot: Optional[str] = None
    selected_events: int = 0

    def status(self, alpha: float = ALPHA) -> str:
        if self.p_compat is None:
            return "no-baseline"
        if self.q_compat is None or self.q_compat >= alpha:
            return "compatible"
        if abs(self.delta_chi2) < DELTA_CHI2_MIN:
            return "changed"
        return "regression" if self.delta_chi2 > 0 else "improvement"

    @classmethod
    def from_dict(cls, d: dict) -> "MeasurementResult":
        names = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in names})


_EMOJI = {"regression": "🚩", "changed": "🔀", "improvement": "⭐",
          "compatible": "✅", "no-baseline": "⚪"}
_RANK = {s: i for i, s in enumerate(_EMOJI)}
_FLAGGED = ("regression", "changed", "improvement")


@dataclass
class MissingMeasurement:
    """A measurement the config asked for that no shard reported."""

    name: str
    experiment: str = ""
    reason: str = ""


@dataclass
class Report:
    results: List[MeasurementResult]
    repo: str = "AchillesGen/Achilles"
    feature_sha: str = "unknown"
    nuisance_version: str = "unknown"
    seed: int = 0
    events_per_measurement: int = 0
    alpha: float = ALPHA
    did_not_run: List[MissingMeasurement] = field(default_factory=list)
    baseline_shas: List[str] = field(default_factory=list)  # main commits compared to
    image_key: str = ""
    run_url: str = ""

    def __post_init__(self):
        compared = [r for r in self.results if r.p_compat is not None]
        for r, q in zip(compared, benjamini_hochberg([r.p_compat for r in compared])):
            r.q_compat = float(q)

    def flagged(self) -> List[MeasurementResult]:
        return self._ordered([r for r in self.results if r.status(self.alpha) in _FLAGGED])

    def without_baseline(self) -> List[MeasurementResult]:
        return [r for r in self.results if r.p_compat is None]

    def overall_ok(self) -> bool:
        return not (self.did_not_run or self.flagged())

    def failed_setups(self) -> List[str]:
        return sorted({m.experiment or "ungrouped" for m in self.did_not_run})

    # -- rendering ------------------------------------------------------------

    def _plot_url(self, r: MeasurementResult) -> Optional[str]:
        return (RAW_URL_TEMPLATE.format(repo=self.repo, sha=self.feature_sha, plot=r.plot)
                if r.plot else None)

    def _link(self, r: MeasurementResult, text: str) -> str:
        """Reference-style link; the definitions go once at the end of the comment."""
        url = self._plot_url(r)
        if not url:
            return text
        if url not in self._ref_ids:
            self._ref_ids[url] = f"p{len(self._refs) + 1}"
            self._refs.append((self._ref_ids[url], url, r.name))
        return f"[{text}][{self._ref_ids[url]}]"

    @staticmethod
    def _common_prefix(names: List[str]) -> str:
        """The shared leading text of a setup's sample names, cut after an '_'."""
        if len(names) < 2:
            return ""
        first, last = min(names), max(names)
        n = next((i for i, (a, b) in enumerate(zip(first, last)) if a != b), len(first))
        return first[:first.rfind("_", 0, n) + 1]

    @staticmethod
    def _short(name: str, prefix: str) -> str:
        short = name[len(prefix):] if prefix else name
        return (short[:-3] if short.endswith("_nu") else short) or name

    def _row(self, r: MeasurementResult, *, setup: bool = False, prefix: str = "") -> str:
        cells = [_EMOJI[r.status(self.alpha)], self._link(r, self._short(r.name, prefix))]
        if setup:
            cells.append(f"`{r.experiment}`")
        cells += [str(r.ndof),
                  f"{_fmt(r.chi2_ndof_main, '.2f')} → {r.chi2_ndof_pr:.2f}",
                  _fmt(r.delta_chi2, "+.1f"), _fmt(r.p_compat), _fmt(r.q_compat),
                  _fmt(r.p_data)]
        return "| " + " | ".join(cells) + " |"

    @staticmethod
    def _table_header(setup: bool = False) -> str:
        cols = ["", "Measurement"] + (["Setup"] if setup else []) + [
            "ndof", "χ²/ndof main → PR", "Δχ²", "p_cmp", "q", "p_data"]
        align = [":-:", "---"] + (["---"] if setup else []) + [
            "--:", ":-:", "--:", "--:", "--:", "--:"]
        return f"| {' | '.join(cols)} |\n|{'|'.join(align)}|"

    def _ordered(self, rows: List[MeasurementResult]) -> List[MeasurementResult]:
        return sorted(rows, key=lambda r: (_RANK[r.status(self.alpha)],
                                           r.p_compat if r.p_compat is not None else 1.0))

    def by_experiment(self) -> List[tuple]:
        """(setup, rows), the setups that need attention first."""
        groups: dict = {}
        for r in self.results:
            groups.setdefault(r.experiment or "ungrouped", []).append(r)
        ordered = {name: self._ordered(rows) for name, rows in groups.items()}
        return sorted(ordered.items(),
                      key=lambda kv: (_RANK[kv[1][0].status(self.alpha)], kv[0]))

    def _setup_summary(self, name: str, rows: List[MeasurementResult]) -> str:
        worst = rows[0].status(self.alpha)
        n_flagged = sum(r.status(self.alpha) in _FLAGGED for r in rows)
        note = (f"{n_flagged} flagged" if n_flagged else
                "no baseline" if worst == "no-baseline" else "all compatible")
        return (f"{_EMOJI[worst]} <b>{name}</b> — "
                f"{_plural(len(rows), 'measurement')}, {note}")

    def to_markdown(self) -> str:
        full = self._render(compact=False)
        return full if len(full) <= MAX_COMMENT_CHARS else self._render(compact=True)

    def _render(self, compact: bool) -> str:
        self._refs: List[tuple] = []
        self._ref_ids: dict = {}
        flagged = self.flagged()
        groups = self.by_experiment()
        missing_base = self.without_baseline()
        n = len(self.results)

        if self.did_not_run:
            verdict = "⚠️ incomplete run"
        elif n and len(missing_base) == n:
            verdict = "⚪ no baseline to compare against"
        elif flagged:
            verdict = "⚠️ significant change"
        else:
            verdict = "✅ no significant change"
        summary = (f"**{verdict}** — **{len(flagged)}** flagged of **{n}** measurements "
                   f"in **{len(groups)}** setups")
        if self.did_not_run:
            summary += (f" · ❌ **{len(self.did_not_run)}** in "
                        f"**{_plural(len(self.failed_setups()), 'setup')}** did not report")

        base = ", ".join(f"`{s[:8]}`" for s in self.baseline_shas) or "none"
        meta = (f"NUISANCE3 `{self.nuisance_version}` · seed `{self.seed}` · "
                f"{self.events_per_measurement:,} events/setup · "
                f"PR `{self.feature_sha[:8]}` vs main {base}")
        if self.run_url:
            meta += f" · [run]({self.run_url})"
        lines = [COMMENT_MARKER, "## 🔬 Physics validation (NUISANCE3)", "",
                 summary, "", meta, ""]

        if missing_base and len(missing_base) < n:
            k = len(missing_base)
            lines += [f"> ⚪ {_plural(k, 'measurement')} {'has' if k == 1 else 'have'} "
                      f"no stored baseline for image `{self.image_key}` and are compared with data "
                      "only. The nightly baseline run fills them in.", ""]
        elif missing_base:
            lines += [f"> ⚪ No stored baseline for image `{self.image_key}`; every "
                      "measurement is compared with data only. Run *Physics Validation* "
                      "with mode `baseline` on main, or wait for the nightly run.", ""]

        if self.did_not_run:
            lines += [f"### ❌ Did not run ({len(self.did_not_run)})", ""]
            by_setup: dict = {}
            for m in self.did_not_run:
                by_setup.setdefault(m.experiment or "ungrouped", []).append(m)
            for setup, ms in sorted(by_setup.items()):
                reason = next((m.reason for m in ms if m.reason), "")
                why = f": `{reason.splitlines()[0][:200]}`" if reason else ""
                lines.append(f"* **{setup}** — {_plural(len(ms), 'measurement')}{why}")
            lines.append("")

        if flagged:
            lines += [f"### Needs attention ({len(flagged)})", "",
                      self._table_header(setup=True)]
            lines += [self._row(r, setup=True) for r in flagged]
            lines.append("")
            with_plots = [r for r in flagged if r.plot]
            for r in with_plots[:MAX_INLINE_PLOTS]:
                lines += [f"<b>{r.name}</b> — {r.status(self.alpha)}<br>",
                          f'<img src="{self._plot_url(r)}" width="420">', ""]
            if len(with_plots) > MAX_INLINE_PLOTS:
                lines += [f"_{len(with_plots) - MAX_INLINE_PLOTS} further flagged plots "
                          "are linked from the table._", ""]

        lines += [f"### All {n} measurements by setup", ""]
        for name, rows in groups:
            if compact and not any(r.status(self.alpha) in _FLAGGED for r in rows):
                lines.append(f"* {self._setup_summary(name, rows)}")
                continue
            prefix = self._common_prefix([r.name for r in rows])
            lines += ["<details>", f"<summary>{self._setup_summary(name, rows)}</summary>",
                      "", self._table_header()]
            lines += [self._row(r, prefix=prefix) for r in rows]
            lines += ["", "</details>"]
        lines.append("")
        if compact:
            lines += ["_Trimmed to fit a PR comment; every row is in `summary.json`._", ""]

        lines += [
            "<details><summary>How to read this</summary>", "",
            f"`p_cmp` is the main-vs-PR compatibility; `q` is it after a Benjamini-Hochberg "
            f"correction over the {n} measurements, and a row is flagged when q < "
            f"{self.alpha} (so ~{self.alpha:.0%} of flags are expected to be noise). Δχ² is "
            "the change in agreement with data and only labels the flag: "
            f"{_EMOJI['regression']} worse, {_EMOJI['improvement']} better, "
            f"{_EMOJI['changed']} |Δχ²| < {DELTA_CHI2_MIN:g}. `p_data` is the PR's fit to "
            "data, for context.", "", "</details>"]
        if self._refs:
            lines.append("")
            lines += [f'[{ref}]: {url} "{name}"' for ref, url, name in self._refs]
        return "\n".join(lines)

    def to_summary_dict(self) -> dict:
        return {
            "repo": self.repo, "feature_sha": self.feature_sha,
            "nuisance_version": self.nuisance_version, "seed": self.seed,
            "events_per_measurement": self.events_per_measurement, "alpha": self.alpha,
            "baseline_shas": self.baseline_shas, "image_key": self.image_key,
            "overall_ok": self.overall_ok(), "n_flagged": len(self.flagged()),
            "did_not_run": [asdict(m) for m in self.did_not_run],
            "failed_setups": self.failed_setups(),
            "measurements": [{**asdict(r), "status": r.status(self.alpha)}
                             for r in self.results],
        }

    def write(self, comment_path: str, summary_path: str) -> None:
        with open(comment_path, "w") as fh:
            fh.write(self.to_markdown())
        with open(summary_path, "w") as fh:
            json.dump(self.to_summary_dict(), fh, indent=2)


# ---------------------------------------------------------------------------
# Unweighting scan
# ---------------------------------------------------------------------------

@dataclass
class VariantMeasurement:
    """One unweighting variant's prediction for one measurement, vs the reference."""

    measurement: str
    variant: str
    ndof: int
    chi2_ndof_data: float
    p_compat: float
    p_shape: float           # p_compat with the normalisation fitted out
    p_data: float
    norm_shift: float        # (Σ variant − Σ reference) / Σ reference
    ess_fraction: float
    max_over_mean: float
    mc_error_ratio: float    # mean sigma_i(variant) / mean sigma_i(reference)


@dataclass
class VariantSummary:
    variant: str
    options: dict
    n_measurements: int
    p_worst: float
    p_overall: float          # Bonferroni over this variant's measurements
    p_shape_worst: float
    p_shape_overall: float
    n_flagged: int
    ess_fraction: float       # medians across measurements
    mc_error_ratio: float
    max_norm_shift: float
    seconds: Optional[float] = None
    unweight_eff: Optional[float] = None
    is_reference: bool = False
    is_null_control: bool = False


@dataclass
class ScanReport:
    summaries: List[VariantSummary]
    rows: List[VariantMeasurement]
    reference: str
    repo: str = "AchillesGen/Achilles"
    feature_sha: str = "unknown"
    nuisance_version: str = "unknown"
    seed: int = 0
    events_per_measurement: int = 0
    alpha: float = ALPHA

    def null_control(self) -> Optional[VariantSummary]:
        return next((s for s in self.summaries if s.is_null_control), None)

    def calibrated(self) -> bool:
        """The null control (the reference reseeded) must itself pass the test."""
        null = self.null_control()
        return null is None or (null.p_overall >= self.alpha
                                and null.p_shape_overall >= self.alpha)

    def thresholds(self) -> tuple:
        """(p, p_shape) thresholds, lowered to whatever the null control scored."""
        null = self.null_control()
        if null is None:
            return (self.alpha, self.alpha)
        return tuple(min(self.alpha, p) if p == p else self.alpha
                     for p in (null.p_overall, null.p_shape_overall))

    def _mark(self, s: VariantSummary) -> str:
        if s.is_reference:
            return "🎯"
        if s.is_null_control:
            return "🧪"
        if not self.calibrated():
            return "❔"
        p_thr, s_thr = self.thresholds()
        return "🚩" if (s.p_overall < p_thr or s.p_shape_overall < s_thr) else "✅"

    def biased(self) -> List[VariantSummary]:
        return [s for s in self.summaries if self._mark(s) == "🚩"]

    def _summary_table(self) -> List[str]:
        lines = ["| Unweighting | p (vs reference) | p (shape only) | ESS/event | "
                 "MC error | max Δnorm | Achilles eff | wall | |",
                 "|---|---|---|---|---|---|---|---|---|"]
        for s in self.summaries:
            p, ps = (("_reference_",) * 2 if s.is_reference
                     else (_fmt(s.p_overall), _fmt(s.p_shape_overall)))
            wall = f"{s.seconds / 60:.1f} min" if s.seconds else "—"
            lines.append(f"| `{s.variant}` | {p} | {ps} | {_fmt(s.ess_fraction, '.3f')} | "
                         f"×{_fmt(s.mc_error_ratio, '.2f')} | "
                         f"{_fmt(s.max_norm_shift, '+.2%')} | {_fmt(s.unweight_eff)} | "
                         f"{wall} | {self._mark(s)} |")
        return lines

    @staticmethod
    def _detail_table(rows: List[VariantMeasurement]) -> List[str]:
        lines = ["| Measurement | Unweighting | ndof | χ²/ndof (data) | p (vs ref) | "
                 "p (shape) | Δnorm | ESS/event | max/mean w |",
                 "|---|---|---|---|---|---|---|---|---|"]
        lines += [f"| {r.measurement} | `{r.variant}` | {r.ndof} | "
                  f"{_fmt(r.chi2_ndof_data, '.2f')} | {_fmt(r.p_compat)} | "
                  f"{_fmt(r.p_shape)} | {_fmt(r.norm_shift, '+.2%')} | "
                  f"{_fmt(r.ess_fraction, '.3f')} | {_fmt(r.max_over_mean, '.1f')} |"
                  for r in rows]
        return lines

    def to_markdown(self) -> str:
        biased = self.biased()
        n_var = sum(not (s.is_reference or s.is_null_control) for s in self.summaries)
        if not self.calibrated():
            verdict = "❔ uncalibrated — the null control fails its own test"
        elif biased:
            verdict = f"⚠️ {len(biased)} of {n_var} scheme(s) differ from the reference"
        else:
            verdict = "✅ every scheme reproduces the reference"
        lines = [COMMENT_MARKER, "## 🎚️ Unweighting scan", "",
                 f"**Reference `{self.reference}` · {verdict}**", "",
                 f"NUISANCE3 `{self.nuisance_version}` · seed `{self.seed}` · "
                 f"{self.events_per_measurement:,} events/variant/setup · "
                 f"`{self.feature_sha[:8]}`", ""]
        lines += self._summary_table()
        lines.append("")
        if self.calibrated():
            flagged = sorted((r for r in self.rows if r.p_compat < self.alpha),
                             key=lambda r: r.p_compat)
            if flagged:
                lines += ["### Measurements differing from the reference"]
                lines += self._detail_table(flagged)
                lines.append("")
        lines += [f"<details><summary>All {len(self.rows)} variant × measurement "
                  "rows</summary>\n"] + self._detail_table(self.rows) + ["\n</details>", ""]
        null = self.null_control()
        if null is not None:
            if self.calibrated():
                p_thr, s_thr = self.thresholds()
                lines.append(f"🧪 The null control (reference reseeded) scores p = "
                             f"{_fmt(null.p_overall)} (shape {_fmt(null.p_shape_overall)}); "
                             f"thresholds are floored there: {_fmt(p_thr)} and "
                             f"{_fmt(s_thr)}.")
            else:
                lines.append(f"❔ The null control scores p = {_fmt(null.p_overall)} "
                             f"(shape {_fmt(null.p_shape_overall)}), so the p columns are "
                             "not usable this run; Δnorm, ESS, MC error and wall time "
                             "still are.")
            lines.append("")
        lines.append("**p (shape)** fits the normalisation out; **ESS/event** is the Kish "
                     "effective sample size per selected event; **MC error** is relative "
                     "to the reference; **Δnorm** is the change in integrated cross "
                     "section.")
        return "\n".join(lines)

    def to_summary_dict(self) -> dict:
        return {
            "kind": "unweighting-scan", "reference": self.reference, "repo": self.repo,
            "feature_sha": self.feature_sha, "nuisance_version": self.nuisance_version,
            "seed": self.seed, "events_per_measurement": self.events_per_measurement,
            "alpha": self.alpha, "biased_variants": [s.variant for s in self.biased()],
            "variants": [asdict(s) for s in self.summaries],
            "rows": [asdict(r) for r in self.rows],
        }

    def write(self, comment_path: str, summary_path: str) -> None:
        with open(comment_path, "w") as fh:
            fh.write(self.to_markdown())
        with open(summary_path, "w") as fh:
            json.dump(self.to_summary_dict(), fh, indent=2)


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _row(name, exp, p, delta=2.0, **kw):
    return MeasurementResult(name, exp, ndof=12, chi2_ndof_pr=1.2, p_data=0.4,
                             chi2_ndof_main=1.1, delta_chi2=delta, p_compat=p,
                             plot=plot_basename(name), **kw)


def _selftest() -> int:
    rows = [
        _row("MINERvA_CC0pi_XSec_1DTp_nu", "MINERvA_CH", 0.62),
        _row("MINERvA_CC0pi_XSec_1Dpmu_nu", "MINERvA_CH", 0.55),
        _row("MiniBooNE_CC1pip_XSec_1DQ2_nu", "MiniBooNE", 1e-4, delta=6.1),
        _row("T2K_CC0pi_XSec_1Dcos_nu", "T2K", 2e-4, delta=-11.0),
        _row("T2K_CC0pi_XSec_1Dp_nu", "T2K", 3e-4, delta=0.2),
        _row("T2K_CC0pi_XSec_1Dq_nu", "T2K", 0.04),   # nominal p < 0.05, not after BH
        MeasurementResult("e12C_new", "e12C", 8, 1.5, 0.3),  # no baseline
    ]
    rep = Report(rows, feature_sha="abcdef1234567890", seed=42,
                 events_per_measurement=500000, baseline_shas=["1234567890ab"],
                 image_key="deadbeef", run_url="https://example/run/1",
                 did_not_run=[MissingMeasurement(
                     "ElectronData_6_12_1.108", "e12C_1108",
                     "achilles failed for e12C_1108 (signal 11)\nlog-tail-line")])
    md = rep.to_markdown()
    summary = rep.to_summary_dict()
    status = {m["name"]: m["status"] for m in summary["measurements"]}

    big = Report([_row(f"Exp{s}_XSec_1DVar{i}_nu", f"Exp{s}", 0.5)
                  for s in range(40) for i in range(15)], feature_sha="abc").to_markdown()
    checks = {
        "marker present": COMMENT_MARKER in md,
        "regression": status["MiniBooNE_CC1pip_XSec_1DQ2_nu"] == "regression",
        "improvement": status["T2K_CC0pi_XSec_1Dcos_nu"] == "improvement",
        "tiny Δχ² is only 'changed'": status["T2K_CC0pi_XSec_1Dp_nu"] == "changed",
        "BH keeps a nominal 0.04 quiet": status["T2K_CC0pi_XSec_1Dq_nu"] == "compatible",
        "no baseline is not flagged": status["e12C_new"] == "no-baseline",
        "no-baseline note": "compared with data only" in md,
        "three flagged": summary["n_flagged"] == 3,
        "incomplete is not ok": not summary["overall_ok"] and "incomplete run" in md,
        "failure reason shown": "(signal 11)`" in md and "log-tail-line" not in md,
        "baseline provenance": "vs main `12345678`" in md,
        "run link": "[run](https://example/run/1)" in md,
        "names cut at an underscore": "| [1DTp][" in md and "| [1Dpmu][" in md,
        "flagged table keeps full names": "[MiniBooNE_CC1pip_XSec_1DQ2][" in md,
        "reference links": "\n[p1]: https://raw.githubusercontent.com" in md,
        "flagged setup first": md.index("<b>T2K</b>") < md.index("<b>MINERvA_CH</b>"),
        "guard trims a huge suite": len(big) <= MAX_COMMENT_CHARS and "Trimmed" in big,
        "plot name sanitised": plot_basename("A_Barreau:1983ht") == "A_Barreau_1983ht.png",
        "round trip": MeasurementResult.from_dict(summary["measurements"][0]).name
                      == rows[0].name,
    }
    for name, ok in checks.items():
        print(f"[{'ok' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("SELFTEST:", "PASS" if ok else "FAIL")
    if "--print" in __import__("sys").argv:
        print("\n" + md)
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    print(__doc__)
