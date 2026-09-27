#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Render the Achilles physval summary as a PR comment and a machine-readable JSON.

The comment is modeled on the Sherpa physval summary table: one row per
measurement, chi2/ndof for main vs the pull request side by side, flagged when the
main-vs-feature compatibility p-value (``p_compat``) drops below ``alpha``.  The
flag is driven *solely* by ``p_compat``; the sign of the change in agreement with
data (``delta_chi2``) only labels a flagged row as a regression or an improvement.

No plotting or NUISANCE dependency here, so it is unit-testable on synthetic rows.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from typing import List, Optional

from stats import bonferroni, bonferroni_threshold


ALPHA = 0.05

# Marker used to find-and-update the single PR comment instead of spamming.
COMMENT_MARKER = "<!-- achilles-physval-summary -->"

RAW_URL_TEMPLATE = (
    "https://raw.githubusercontent.com/{repo}/physval-baselines/"
    "plots/{sha}/{plot}"
)

# actions/upload-artifact rejects these outright in a *path* inside the artifact, and
# they are no better in a git tree or a Windows checkout. Sample names are not free of
# them: the Durham electron samples carry their reference as
# ElectronData_6_12_0.560_36.000_Barreau:1983ht.
_UNSAFE_IN_FILENAMES = '":<>|*?\r\n'


def plot_basename(measurement: str) -> str:
    """The on-disk name of a measurement's overlay plot.

    Only the file name is sanitised: the measurement keeps its published name
    everywhere it is displayed, compared or looked up.
    """
    return "".join("_" if c in _UNSAFE_IN_FILENAMES else c
                   for c in measurement) + ".png"

# Thumbnails are only inlined for flagged rows, and only this many: a comment with
# fifty embedded PNGs is unreadable and slow to load. Everything else is one click
# away behind its measurement name.
MAX_INLINE_PLOTS = 6

# GitHub rejects a comment body over 65536 characters. Past this, the renderer drops
# the tables of setups that have nothing flagged.
MAX_COMMENT_CHARS = 60000


@dataclass
class MeasurementResult:
    name: str
    ndof: int
    chi2_ndof_main: float      # main prediction vs data
    chi2_ndof_pr: float        # PR prediction vs data
    delta_chi2: float          # chi2_pr(vs data) - chi2_main(vs data); + = worse
    p_compat: float            # main vs feature compatibility p-value (drives flag)
    p_data: float              # PR vs data goodness-of-fit p-value (context only)
    plot: Optional[str] = None  # basename of the overlay plot, e.g. "<name>.png"
    experiment: str = ""        # the setup it was generated with; groups the comment

    def status(self, alpha: float = ALPHA) -> str:
        """One of 'regression', 'improvement', 'compatible'."""
        if self.p_compat >= alpha:
            return "compatible"
        return "regression" if self.delta_chi2 > 0 else "improvement"


_EMOJI = {"regression": "🚩", "improvement": "⭐", "compatible": "✅"}
# Flagged rows first (regression, then improvement), compatible last.
_SORT_RANK = {"regression": 0, "improvement": 1, "compatible": 2}


@dataclass
class MissingMeasurement:
    """A measurement the config asked for that no shard reported.

    Its setup's job crashed or timed out. Kept in the report rather than dropped: a
    comment that silently omits a third of the suite reads like a pass.
    """

    name: str
    experiment: str = ""


@dataclass
class Report:
    results: List[MeasurementResult]
    repo: str = "AchillesGen/Achilles"
    feature_sha: str = "unknown"
    nuisance_version: str = "unknown"
    seed: int = 0
    events_per_measurement: int = 0
    alpha: float = ALPHA
    extra_header: List[str] = field(default_factory=list)
    did_not_run: List[MissingMeasurement] = field(default_factory=list)

    # -- derived quantities ---------------------------------------------------

    def p_overall(self) -> float:
        return bonferroni([r.p_compat for r in self.results])

    def n_flagged(self) -> int:
        return sum(1 for r in self.results if r.status(self.alpha) != "compatible")

    def overall_ok(self) -> bool:
        if self.did_not_run:
            return False  # part of the suite never reported; nothing to be ok about
        po = self.p_overall()
        return not (po == po and po < self.alpha)  # NaN-safe: ok if not < alpha

    def failed_setups(self) -> "List[str]":
        """Setups with at least one measurement that never reported."""
        return sorted({m.experiment or "ungrouped" for m in self.did_not_run})

    def _plot_url(self, r: MeasurementResult) -> Optional[str]:
        if not r.plot:
            return None
        # Built from the stored basename, not the measurement name: the two differ
        # wherever the name carries a character a file path cannot.
        return RAW_URL_TEMPLATE.format(repo=self.repo, sha=self.feature_sha,
                                       plot=r.plot)

    # -- rendering ------------------------------------------------------------

    def _link(self, r: MeasurementResult, text: str) -> str:
        """A reference-style link to the overlay plot, so table rows stay short.

        The definitions are emitted once at the end of the comment. Inline URLs make
        every row ~180 characters of raw markdown, which is unreadable in a diff or in
        the edit box once there are fifty of them.
        """
        url = self._plot_url(r)
        if not url:
            return text
        ref = self._ref_ids.get(url)
        if ref is None:
            ref = f"p{len(self._refs) + 1}"
            self._ref_ids[url] = ref
            self._refs.append((ref, url, r.name))
        return f"[{text}][{ref}]"

    @staticmethod
    def _short(name: str, prefix: str) -> str:
        """Drop the part of a name every row in its setup shares, plus the `_nu` tail."""
        short = name[len(prefix):] if prefix and name.startswith(prefix) else name
        short = short[:-3] if short.endswith("_nu") else short
        return short.strip("_") or name

    @staticmethod
    def _common_prefix(names: "List[str]") -> str:
        """The shared leading text of a setup's sample names, cut at a separator."""
        if len(names) < 2:
            return ""
        prefix = os.path.commonprefix(names)
        cut = max(prefix.rfind("_"), prefix.rfind("1D"))
        return prefix[:cut + 1] if cut > 0 else ""

    def _row(self, r: MeasurementResult, *, setup: bool = False,
             prefix: str = "") -> str:
        """One table row. ``χ²/ndof`` is written main → PR to keep the table narrow."""
        st = r.status(self.alpha)
        name_cell = self._link(r, self._short(r.name, prefix))
        cells = [_EMOJI[st], name_cell]
        if setup:
            cells.append(f"`{r.experiment}`" if r.experiment else "")
        cells += [str(r.ndof),
                  f"{r.chi2_ndof_main:.2f} → {r.chi2_ndof_pr:.2f}",
                  f"{r.delta_chi2:+.1f}",
                  f"{r.p_compat:.3g}",
                  f"{r.p_data:.3g}"]
        return "| " + " | ".join(cells) + " |"

    @staticmethod
    def _table_header(setup: bool = False) -> str:
        cols = ["", "Measurement"] + (["Setup"] if setup else []) + [
            "ndof", "χ²/ndof main → PR", "Δχ²", "p_cmp", "p_data"]
        return ("| " + " | ".join(cols) + " |\n"
                "|" + "|".join([":-:", "---"] + (["---"] if setup else []) +
                               ["--:", ":-:", "--:", "--:", "--:"]) + "|")

    def _ordered(self, results: List[MeasurementResult]) -> List[MeasurementResult]:
        """Flagged first (regression, then improvement), then by how big the move was."""
        return sorted(results, key=lambda r: (_SORT_RANK[r.status(self.alpha)],
                                              -abs(r.delta_chi2)))

    def by_experiment(self) -> "List[tuple]":
        """(setup, its results) with the setups that need attention first."""
        groups: dict = {}
        for r in self.results:
            groups.setdefault(r.experiment or "ungrouped", []).append(r)

        def rank(item):
            name, rows = item
            worst = min((_SORT_RANK[r.status(self.alpha)] for r in rows), default=2)
            return (worst, name)

        return [(name, self._ordered(rows)) for name, rows in sorted(groups.items(),
                                                                     key=rank)]

    def _setup_summary(self, name: str, rows: List[MeasurementResult]) -> str:
        """The one line you read without expanding a setup."""
        flagged = [r for r in rows if r.status(self.alpha) != "compatible"]
        emoji = _EMOJI["compatible"]
        if flagged:
            emoji = _EMOJI[min((r.status(self.alpha) for r in flagged),
                               key=lambda st: _SORT_RANK[st])]
        worst_p = min((r.p_compat for r in rows), default=float("nan"))
        note = (f"{len(flagged)} flagged" if flagged else "all compatible")
        return (f"{emoji} <b>{name}</b> — {len(rows)} measurement"
                f"{'s' if len(rows) != 1 else ''}, {note} "
                f"<i>(lowest p_compat {worst_p:.3g})</i>")

    def to_markdown(self) -> str:
        """The comment, trimmed if the full one would not fit in a PR comment."""
        full = self._render(compact=False)
        if len(full) <= MAX_COMMENT_CHARS:
            return full
        return self._render(compact=True)

    def _render(self, compact: bool) -> str:
        self._refs: List[tuple] = []
        self._ref_ids: dict = {}
        po = self.p_overall()
        n = len(self.results)
        groups = self.by_experiment()
        flagged = [r for r in self._ordered(self.results)
                   if r.status(self.alpha) != "compatible"]

        if self.did_not_run:
            verdict = "⚠️ incomplete run"
        elif self.overall_ok():
            verdict = "✅ no significant change"
        else:
            verdict = "⚠️ significant change"
        lines: List[str] = [COMMENT_MARKER, "## 🔬 Physics validation (NUISANCE3)", ""]
        summary = (f"**{verdict}** — Bonferroni p = `{po:.3g}` · "
                   f"**{len(flagged)}** flagged of **{n}** measurements in "
                   f"**{len(groups)}** setups")
        if self.did_not_run:
            setups = self.failed_setups()
            summary += (f" · ❌ **{len(self.did_not_run)}** measurements in "
                        f"**{len(setups)}** setup{'s' if len(setups) != 1 else ''} "
                        f"did not report")
        lines.append(summary)
        lines.append("")
        meta = (f"NUISANCE3 `{self.nuisance_version}` · seed `{self.seed}` · "
                f"{self.events_per_measurement:,} events/setup · "
                f"feature `{self.feature_sha[:8]}`")
        lines.append(meta)
        for extra in self.extra_header:
            lines.append(extra)
        lines.append("")

        # Setups whose job died. Their numbers are absent, not compatible, so they get
        # their own block above everything else.
        if self.did_not_run:
            by_setup: dict = {}
            for m in self.did_not_run:
                by_setup.setdefault(m.experiment or "ungrouped", []).append(m.name)
            lines.append(f"### ❌ Did not run ({len(self.did_not_run)})")
            lines.append("")
            lines.append("These setups' jobs failed, so the suite is incomplete and the "
                         "verdict above covers only what reported. Check the run's job "
                         "logs.")
            lines.append("")
            for setup, names in sorted(by_setup.items()):
                lines.append(f"* **{setup}** — {len(names)} measurement"
                             f"{'s' if len(names) != 1 else ''}")
            lines.append("")

        # What needs attention, across every setup, with the setup named per row.
        if flagged:
            lines.append(f"### Needs attention ({len(flagged)})")
            lines.append("")
            lines.append(self._table_header(setup=True))
            for r in flagged:
                lines.append(self._row(r, setup=True))
            lines.append("")
            shown = [r for r in flagged if self._plot_url(r)][:MAX_INLINE_PLOTS]
            if shown:
                for r in shown:
                    lines.append(f"<b>{r.name}</b> — {r.status(self.alpha)}<br>")
                    lines.append(f'<img src="{self._plot_url(r)}" width="420">')
                    lines.append("")
                hidden = len([r for r in flagged if self._plot_url(r)]) - len(shown)
                if hidden > 0:
                    lines.append(f"_{hidden} further flagged plot"
                                 f"{'s' if hidden != 1 else ''} are linked from the "
                                 f"table above._")
                    lines.append("")
        elif not self.did_not_run:
            lines.append("Every measurement is compatible with the stored `main` "
                         "baseline.")
            lines.append("")
        else:
            lines.append("Nothing that reported is flagged.")
            lines.append("")

        # One collapsible per setup: the comment stays the same length whether the
        # suite has three setups or thirty.
        lines.append(f"### All {n} measurements by setup")
        lines.append("")
        for name, rows in groups:
            quiet = all(r.status(self.alpha) == "compatible" for r in rows)
            if compact and quiet:
                # Too many measurements to print every table: a quiet setup is one
                # line, and its numbers stay in summary.json.
                lines.append(f"* {self._setup_summary(name, rows)}")
                continue
            lines.append("<details>")
            lines.append(f"<summary>{self._setup_summary(name, rows)}</summary>")
            lines.append("")
            lines.append(self._table_header())
            prefix = self._common_prefix([r.name for r in rows])
            for r in rows:
                lines.append(self._row(r, prefix=prefix))
            lines.append("")
            lines.append("</details>")
        lines.append("")
        if compact:
            lines.append("_Trimmed to fit a PR comment: setups with nothing flagged "
                         "are summarised rather than tabulated. Every measurement is "
                         "in `summary.json`._")
            lines.append("")

        # Legend and the multiple-comparison arithmetic, out of the way.
        thr = bonferroni_threshold(n, self.alpha)
        expected_false = self.alpha * n
        lines.append("<details><summary>How to read this</summary>")
        lines.append("")
        lines.append(f"* {_EMOJI['regression']} **regression** — p_compat < "
                     f"{self.alpha} and Δχ² > 0 (agreement with data got worse)")
        lines.append(f"* {_EMOJI['improvement']} **improvement** — p_compat < "
                     f"{self.alpha} and Δχ² < 0")
        lines.append(f"* {_EMOJI['compatible']} **compatible** — p_compat ≥ "
                     f"{self.alpha}")
        lines.append("")
        lines.append("The flag is driven *only* by `p_compat`, the main-vs-PR "
                     "compatibility; the sign of Δχ² only labels its direction. "
                     "`p (data)` is the PR's goodness of fit to the published data, "
                     "for context — a sample can disagree with data and still be "
                     "perfectly compatible with `main`.")
        lines.append("")
        lines.append(f"At uncorrected α={self.alpha} across N={n} measurements, "
                     f"~{expected_false:.1f} false flags are expected by chance; the "
                     f"Bonferroni per-measurement threshold is α/N = {thr:.4g}.")
        lines.append("")
        lines.append("</details>")

        if self._refs:
            lines.append("")
            for ref, url, name in self._refs:
                lines.append(f"[{ref}]: {url} \"{name}\"")
        return "\n".join(lines)

    def to_summary_dict(self) -> dict:
        return {
            "repo": self.repo,
            "feature_sha": self.feature_sha,
            "nuisance_version": self.nuisance_version,
            "seed": self.seed,
            "events_per_measurement": self.events_per_measurement,
            "alpha": self.alpha,
            "p_overall": self.p_overall(),
            "overall_ok": self.overall_ok(),
            "n_flagged": self.n_flagged(),
            "did_not_run": [asdict(m) for m in self.did_not_run],
            "failed_setups": self.failed_setups(),
            "measurements": [
                {**asdict(r), "status": r.status(self.alpha)} for r in self.results
            ],
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
    """One unweighting variant's prediction for one measurement.

    ``p_compat`` is against the *reference* variant (the fully weighted sample), so a
    small value means the scheme is biasing the distribution — a correctness failure,
    not a cost. ``ess_fraction`` is the statistical power the scheme delivered per
    selected event, and ``mc_error_ratio`` the resulting error band relative to the
    reference's.
    """

    measurement: str
    variant: str
    ndof: int
    chi2_ndof_data: float
    p_compat: float
    p_shape: float          # as p_compat, with the overall normalisation divided out
    p_data: float
    norm_shift: float        # (Σ variant − Σ reference) / Σ reference
    ess_fraction: float
    max_over_mean: float
    mc_error_ratio: float    # mean sigma_i(variant) / mean sigma_i(reference)


@dataclass
class VariantSummary:
    """Per-variant rollup across every measurement in the scan."""

    variant: str
    options: dict
    n_measurements: int
    p_worst: float
    p_overall: float          # Bonferroni over this variant's measurements
    p_shape_worst: float
    p_shape_overall: float
    n_flagged: int
    ess_fraction: float       # median across measurements
    mc_error_ratio: float     # median across measurements
    max_norm_shift: float
    seconds: Optional[float] = None
    unweight_eff: Optional[float] = None
    is_reference: bool = False
    is_null_control: bool = False


@dataclass
class ScanReport:
    """The unweighting-scan comment: a per-variant verdict plus the detail rows."""

    summaries: List[VariantSummary]
    rows: List[VariantMeasurement]
    reference: str
    repo: str = "AchillesGen/Achilles"
    feature_sha: str = "unknown"
    nuisance_version: str = "unknown"
    seed: int = 0
    events_per_measurement: int = 0
    alpha: float = ALPHA
    extra_header: List[str] = field(default_factory=list)
    did_not_run: List[MissingMeasurement] = field(default_factory=list)

    def biased(self) -> List[VariantSummary]:
        """Variants whose distributions differ from the reference beyond MC noise."""
        return [s for s in self.summaries
                if not (s.is_reference or s.is_null_control)
                and self._mark(s) == "🚩"]

    @staticmethod
    def _fmt(value: Optional[float], spec: str = ".3g", dash: str = "—") -> str:
        if value is None or value != value:
            return dash
        return format(value, spec)

    def _summary_table(self) -> List[str]:
        lines = [
            "| Unweighting | p (vs reference) | p (shape only) | ESS/event | "
            "MC error | max Δnorm | Achilles eff | wall | |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for s in self.summaries:
            if s.is_reference:
                mark, pcell, scell = "🎯", "_reference_", "_reference_"
            else:
                mark = self._mark(s)
                pcell, scell = (self._fmt(s.p_overall),
                                self._fmt(s.p_shape_overall))
            wall = (f"{s.seconds / 60:.1f} min" if s.seconds else "—")
            lines.append(
                f"| `{s.variant}` | {pcell} | {scell} | "
                f"{self._fmt(s.ess_fraction, '.3f')} | "
                f"×{self._fmt(s.mc_error_ratio, '.2f')} | "
                f"{self._fmt(s.max_norm_shift, '+.2%')} | "
                f"{self._fmt(s.unweight_eff, '.3g')} | {wall} | {mark} |")
        return lines

    def _mark(self, s: "VariantSummary") -> str:
        """🎯 reference · 🧪 null control · ❔ uncalibrated · 🚩 flagged · ✅ ok."""
        if s.is_null_control:
            return "🧪"
        if not self.calibrated():
            return "❔"
        p_thr, s_thr = self.thresholds()
        return "🚩" if (s.p_overall < p_thr or s.p_shape_overall < s_thr) else "✅"

    def null_control(self) -> Optional["VariantSummary"]:
        return next((s for s in self.summaries if s.is_null_control), None)

    def calibrated(self) -> bool:
        """Whether the p-values mean anything for this run.

        The null control is the reference's own configuration at a different seed, so
        it is a draw from the null hypothesis and ought to sit above ``alpha``. When
        it does not, the covariance is missing variance that is present between any
        two runs, every p-value is compressed against zero, and ranking the variants
        by p would be reading noise — several of them underflow to 0.0 outright.
        The report says so instead of naming a culprit.
        """
        null = self.null_control()
        if null is None:
            return True  # no control was run; fall back to the nominal alpha
        return null.p_overall >= self.alpha and null.p_shape_overall >= self.alpha

    def thresholds(self) -> tuple:
        """Flagging thresholds for (p_compat, p_shape), floored by the null control.

        Whatever the control scores is what an *identical* configuration costs, and
        no real scheme should be called out for doing at least as well. The control
        can therefore only make the test more conservative — ``min`` with ``alpha`` —
        never less, so a control that lands at p ≈ 1 leaves the threshold untouched.
        """
        null = self.null_control()
        if null is None:
            return (self.alpha, self.alpha)
        p, q = null.p_overall, null.p_shape_overall
        return (min(self.alpha, p if p == p else self.alpha),
                min(self.alpha, q if q == q else self.alpha))

    def _detail_table(self, rows: List[VariantMeasurement]) -> List[str]:
        lines = [
            "| Measurement | Unweighting | ndof | χ²/ndof (data) | p (vs ref) | "
            "p (shape) | Δnorm | ESS/event | max/mean w |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for r in rows:
            lines.append(
                f"| {r.measurement} | `{r.variant}` | {r.ndof} | "
                f"{self._fmt(r.chi2_ndof_data, '.2f')} | {self._fmt(r.p_compat)} | "
                f"{self._fmt(r.p_shape)} | {self._fmt(r.norm_shift, '+.2%')} | "
                f"{self._fmt(r.ess_fraction, '.3f')} | "
                f"{self._fmt(r.max_over_mean, '.1f')} |")
        return lines

    def to_markdown(self) -> str:
        biased = self.biased()
        n_var = len([s for s in self.summaries
                     if not (s.is_reference or s.is_null_control)])
        if not self.calibrated():
            verdict = ("❔ uncalibrated — the null control fails its own test, so "
                       "no scheme can be judged on p")
        elif biased:
            verdict = f"⚠️ {len(biased)} of {n_var} scheme(s) differ from the reference"
        else:
            verdict = "✅ every scheme reproduces the reference"

        lines: List[str] = [
            COMMENT_MARKER, "## 🎚️ Unweighting scan", "",
            f"**Reference `{self.reference}` · {verdict}**", "",
            f"NUISANCE3 `{self.nuisance_version}` · seed `{self.seed}` · "
            f"{self.events_per_measurement:,} events/variant/setup "
            f"· `{self.feature_sha[:8]}`",
        ]
        lines.extend(self.extra_header)
        lines.append("")
        lines.extend(self._summary_table())
        lines.append("")

        flagged_rows = ([] if not self.calibrated()
                        else [r for r in self.rows if r.p_compat < self.alpha])
        if flagged_rows:
            lines.append("### Measurements differing from the reference")
            lines.extend(self._detail_table(
                sorted(flagged_rows, key=lambda r: r.p_compat)))
            lines.append("")

        lines.append(f"<details><summary>All {len(self.rows)} variant × "
                     "measurement rows</summary>\n")
        lines.extend(self._detail_table(self.rows))
        lines.append("\n</details>")
        lines.append("")
        null = self.null_control()
        if null is not None and not self.calibrated():
            lines.append(
                f"❔ **The p columns above are not usable for this run.** 🧪 "
                f"`null-control` is the reference's own configuration at a different "
                f"seed — a draw from the null hypothesis, which should sit above "
                f"α={self.alpha}. It scores p = {self._fmt(null.p_overall)} "
                f"(shape {self._fmt(null.p_shape_overall)}). Two *identical* "
                f"configurations are therefore \"incompatible\", so the covariance is "
                f"missing variance that is present between any two runs: the bootstrap "
                f"is estimated inside a single run and carries neither the run-to-run "
                f"scatter of the flux-averaged cross section nor the spread from each "
                f"run's own adapted integration grid. Ranking variants by p here would "
                f"be reading noise — some p-values underflow to 0.0 outright.")
            lines.append("")
            lines.append(
                "The **effect sizes are still valid**: Δnorm, ESS/event, MC error and "
                "wall time are direct measurements, not test statistics. Compare each "
                "scheme's Δnorm against the null control's — a scheme inside that is "
                "indistinguishable from rerunning the reference.")
        elif null is not None:
            p_thr, s_thr = self.thresholds()
            lines.append(
                f"🧪 **null control** — the reference rerun at a different seed, so its "
                f"rows are drawn from the null. It scores p = "
                f"{self._fmt(null.p_overall)} (shape "
                f"{self._fmt(null.p_shape_overall)}); thresholds are floored there "
                f"— {self._fmt(p_thr)} and {self._fmt(s_thr)} — so no scheme is "
                f"flagged for doing as well as an identical rerun.")
            lines.append("")
        lines.append(
            "Legend: **p (vs ref)** — correlated χ² against the reference variant's "
            "histogram using both bootstrap covariances. **p (shape)** — the same "
            "with the overall normalisation divided out and one dof given up for it; "
            "the pair separates \"this scheme moved the distribution\" from "
            "\"these two runs disagree on the total cross section\". "
            "**ESS/event** — Kish effective sample size per selected event "
            "(1.0 = unit weights); **MC error** — mean bootstrap error relative to "
            "the reference; **Δnorm** — change in the integrated cross section; "
            "**max/mean w** — heaviest surviving overweight."
        )
        return "\n".join(lines)

    def to_summary_dict(self) -> dict:
        return {
            "kind": "unweighting-scan",
            "reference": self.reference,
            "repo": self.repo,
            "feature_sha": self.feature_sha,
            "nuisance_version": self.nuisance_version,
            "seed": self.seed,
            "events_per_measurement": self.events_per_measurement,
            "alpha": self.alpha,
            "biased_variants": [s.variant for s in self.biased()],
            "variants": [asdict(s) for s in self.summaries],
            "rows": [asdict(r) for r in self.rows],
        }

    def write(self, comment_path: str, summary_path: str) -> None:
        with open(comment_path, "w") as fh:
            fh.write(self.to_markdown())
        with open(summary_path, "w") as fh:
            json.dump(self.to_summary_dict(), fh, indent=2)


def _selftest_incomplete() -> bool:
    """A report missing a setup says so, and cannot come out green."""
    rows = [MeasurementResult("A_XSec_1DVar_nu", 12, 1.1, 1.1, 0.0, 0.8, 0.5,
                              plot="A_XSec_1DVar_nu.png", experiment="ExpA")]
    rep = Report(rows, feature_sha="abc",
                 did_not_run=[MissingMeasurement("B_XSec_1DVar_nu", "ExpB"),
                              MissingMeasurement("B_XSec_1DOther_nu", "ExpB")])
    md = rep.to_markdown()
    summary = rep.to_summary_dict()
    return (not rep.overall_ok()                     # cannot be green
            and "incomplete run" in md               # verdict says so
            and "Did not run (2)" in md              # and lists them
            and "**ExpB** — 2 measurements" in md
            and rep.failed_setups() == ["ExpB"]
            and len(summary["did_not_run"]) == 2
            and summary["failed_setups"] == ["ExpB"])


def _selftest_guard() -> bool:
    """A suite far larger than today's must still fit in one PR comment."""
    rows = [MeasurementResult(f"Exp{s}_XSec_1DVar{i}_nu", 12, 1.1, 1.2, 0.4, 0.5, 0.4,
                              plot=f"Exp{s}_XSec_1DVar{i}_nu.png",
                              experiment=f"Exp{s}")
            for s in range(40) for i in range(15)]
    md = Report(rows, feature_sha="abc123def456", seed=1,
                events_per_measurement=500000).to_markdown()
    return len(md) <= MAX_COMMENT_CHARS and "Trimmed to fit" in md


def _selftest() -> int:
    results = [
        MeasurementResult("MINERvA_CC0pi_Tp", 38, 1.11, 1.10, -0.4, 0.62, 0.31,
                          plot="MINERvA_CC0pi_Tp.png", experiment="MINERvA_CC"),
        MeasurementResult("MINERvA_CC0pi_pmu", 21, 1.02, 1.04, +0.2, 0.55, 0.44,
                          plot="MINERvA_CC0pi_pmu.png", experiment="MINERvA_CC"),
        MeasurementResult("MiniBooNE_CC1pip_Q2", 40, 1.38, 1.53, +6.1, 0.013, 0.04,
                          plot="MiniBooNE_CC1pip_Q2.png", experiment="MiniBooNE_CC1pi"),
        MeasurementResult("T2K_CC0pi_cosTheta", 58, 0.98, 0.79, -11.0, 0.008, 0.71,
                          plot="T2K_CC0pi_cosTheta.png", experiment="T2K_CC"),
    ]
    rep = Report(results, feature_sha="abcdef1234567890", nuisance_version="v3.0.1",
                 seed=42, events_per_measurement=500000)
    md = rep.to_markdown()
    summary = rep.to_summary_dict()

    checks = {
        "marker present": COMMENT_MARKER in md,
        "bonferroni p = min(1, 4*0.008)=0.032": abs(summary["p_overall"] - 0.032) < 1e-9,
        "overall flagged": summary["overall_ok"] is False,
        "two flagged rows": summary["n_flagged"] == 2,
        "regression labelled": next(m for m in summary["measurements"]
                                    if m["name"].startswith("MiniBooNE"))["status"] == "regression",
        "improvement labelled": next(m for m in summary["measurements"]
                                     if m["name"].startswith("T2K"))["status"] == "improvement",
        "compatible collapsed": "<details>" in md,
        "raw url embedded": "raw.githubusercontent.com" in md,
        # grouping: one collapsible per setup, worst setup first
        "one details per setup": md.count("<summary>") == 4,  # 3 setups + the legend
        "setups named in summaries": all(f"<b>{s}</b>" in md for s in
                                         ("MINERvA_CC", "MiniBooNE_CC1pi", "T2K_CC")),
        "flagged setup listed first": (md.index("<b>MiniBooNE_CC1pi</b>")
                                       < md.index("<b>MINERvA_CC</b>")),
        "quiet setup says so": "all compatible" in md,
        # rows are reference links, with the definitions emitted once each
        "reference links used": "][p1]" in md and "\n[p1]: http" in md,
        "no inline urls in rows": "| [" in md and "](http" not in md,
        # inside a setup the shared prefix goes; the full name stays in the link title
        "names shortened per setup": "| [Tp][" in md and "| [pmu][" in md,
        "full name kept in the link": '"MINERvA_CC0pi_Tp"' in md,
        "flagged table keeps full names": "[MiniBooNE_CC1pip_Q2][" in md,
        # the size guard trims instead of producing an over-long comment
        "guard trims a huge suite": _selftest_guard(),
        # a name that cannot be a file path still gets a usable plot file + url
        "colon stripped from plot name":
            plot_basename("ElectronData_6_12_0.560_36.000_Barreau:1983ht")
            == "ElectronData_6_12_0.560_36.000_Barreau_1983ht.png",
        "plot name keeps the rest verbatim":
            plot_basename("MINERvA_CC0pi_Tp") == "MINERvA_CC0pi_Tp.png",
        # an incomplete run is reported, not silently shrunk
        "incomplete run marked": _selftest_incomplete(),
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
