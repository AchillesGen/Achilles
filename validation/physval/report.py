#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Render the physval summary as a PR comment (markdown) and summary.json.

A row is flagged when its main-vs-PR compatibility survives a Benjamini-Hochberg
correction across the suite; the change in chi-square against data only labels the
direction.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import List, Optional

from stats import benjamini_hochberg

ALPHA = 0.05
# Below this |Δχ²| a flagged row is called "changed" rather than better or worse.
DELTA_CHI2_MIN = 1.0
COMMENT_MARKER = "<!-- achilles-physval-summary -->"
PLOT_URL_TEMPLATE = "https://achillesgen.github.io/physval-artifacts/plots/{sha}/{plot}"
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

    def status(self) -> str:
        if self.p_compat is None:
            return "no-baseline"
        if self.q_compat is None or self.q_compat >= ALPHA:
            return "compatible"
        if abs(self.delta_chi2) < DELTA_CHI2_MIN:
            return "changed"
        return "regression" if self.delta_chi2 > 0 else "improvement"


_EMOJI = {"regression": "🚩", "changed": "🔀", "improvement": "⭐",
          "compatible": "✅", "no-baseline": "⚪"}
_RANK = {s: i for i, s in enumerate(_EMOJI)}
_FLAGGED = ("regression", "changed", "improvement")


@dataclass
class MissingMeasurement:
    """A measurement with no prediction. ``reason`` is why, or for one that selected no
    events, why that may be acceptable (empty: it is not)."""

    name: str
    experiment: str
    reason: str = ""


@dataclass
class Report:
    results: List[MeasurementResult]
    feature_sha: str = "unknown"
    nuisance_version: str = "unknown"
    seed: int = 0
    events_per_measurement: int = 0
    did_not_run: List[MissingMeasurement] = field(default_factory=list)
    no_events: List[MissingMeasurement] = field(default_factory=list)
    baseline_shas: List[str] = field(default_factory=list)  # main commits compared to
    image_key: str = ""
    run_url: str = ""

    def __post_init__(self):
        compared = [r for r in self.results if r.p_compat is not None]
        for r, q in zip(compared, benjamini_hochberg([r.p_compat for r in compared])):
            r.q_compat = float(q)

    def flagged(self) -> List[MeasurementResult]:
        return self._ordered([r for r in self.results if r.status() in _FLAGGED])

    def without_baseline(self) -> List[MeasurementResult]:
        return [r for r in self.results if r.p_compat is None]

    def unexpected_empty(self) -> List[MissingMeasurement]:
        return [m for m in self.no_events if not m.reason]

    def failed(self) -> bool:
        return bool(self.did_not_run or self.unexpected_empty())

    def failed_setups(self) -> List[str]:
        return sorted({m.experiment for m in self.did_not_run})

    # -- rendering ------------------------------------------------------------

    def _plot_url(self, r: MeasurementResult) -> str:
        return PLOT_URL_TEMPLATE.format(sha=self.feature_sha, plot=plot_basename(r.name))

    def _link(self, r: MeasurementResult, text: str) -> str:
        """Reference-style link; the definitions go once at the end of the comment."""
        url = self._plot_url(r)
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
        cells = [_EMOJI[r.status()], self._link(r, self._short(r.name, prefix))]
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
        return sorted(rows, key=lambda r: (_RANK[r.status()],
                                           r.p_compat if r.p_compat is not None else 1.0))

    def by_experiment(self) -> List[tuple]:
        """(setup, rows), the setups that need attention first."""
        groups: dict = {}
        for r in self.results:
            groups.setdefault(r.experiment, []).append(r)
        ordered = {name: self._ordered(rows) for name, rows in groups.items()}
        return sorted(ordered.items(), key=lambda kv: (_RANK[kv[1][0].status()], kv[0]))

    def _setup_summary(self, name: str, rows: List[MeasurementResult]) -> str:
        worst = rows[0].status()
        n_flagged = sum(r.status() in _FLAGGED for r in rows)
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
        empty = self.unexpected_empty()
        n = len(self.results)

        if self.failed():
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
        if empty:
            summary += f" · ❌ **{len(empty)}** with no events selected"

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
                by_setup.setdefault(m.experiment, []).append(m)
            for setup, ms in sorted(by_setup.items()):
                reason = next((m.reason for m in ms if m.reason), "")
                why = f": `{reason.splitlines()[0][:200]}`" if reason else ""
                lines.append(f"* **{setup}** — {_plural(len(ms), 'measurement')}{why}")
            lines.append("")

        if self.no_events:
            lines += [f"### No events selected ({len(self.no_events)})", ""]
            for m in sorted(self.no_events, key=lambda m: bool(m.reason)):
                note = f" — ⚠️ may be empty: {m.reason}" if m.reason else ""
                lines.append(f"* {'' if m.reason else '❌ '}`{m.name}` (`{m.experiment}`){note}")
            lines.append("")

        if flagged:
            lines += [f"### Needs attention ({len(flagged)})", "",
                      self._table_header(setup=True)]
            lines += [self._row(r, setup=True) for r in flagged]
            lines.append("")
            for r in flagged[:MAX_INLINE_PLOTS]:
                lines += [f"<b>{r.name}</b> — {r.status()}<br>",
                          f'<img src="{self._plot_url(r)}" width="420">', ""]
            if len(flagged) > MAX_INLINE_PLOTS:
                lines += [f"_{len(flagged) - MAX_INLINE_PLOTS} further flagged plots "
                          "are linked from the table._", ""]

        lines += [f"### All {n} measurements by setup", ""]
        for name, rows in groups:
            if compact and not any(r.status() in _FLAGGED for r in rows):
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
            f"{ALPHA} (so ~{ALPHA:.0%} of flags are expected to be noise). Δχ² is "
            "the change in agreement with data and only labels the flag: "
            f"{_EMOJI['regression']} worse, {_EMOJI['improvement']} better, "
            f"{_EMOJI['changed']} |Δχ²| < {DELTA_CHI2_MIN:g}. `p_data` is the PR's fit to "
            "data, for context.", "", "</details>"]
        if self._refs:
            lines.append("")
            lines += [f'[{ref}]: {url} "{name}"' for ref, url, name in self._refs]
        return "\n".join(lines)

    def write(self, comment_path: str, summary_path: str) -> None:
        with open(comment_path, "w") as fh:
            fh.write(self.to_markdown())
        with open(summary_path, "w") as fh:
            json.dump(asdict(self), fh, indent=2)
