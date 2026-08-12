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
from dataclasses import asdict, dataclass, field
from typing import List, Optional

from stats import bonferroni, bonferroni_threshold


ALPHA = 0.05

# Marker used to find-and-update the single PR comment instead of spamming.
COMMENT_MARKER = "<!-- achilles-physval-summary -->"

RAW_URL_TEMPLATE = (
    "https://raw.githubusercontent.com/{repo}/physval-baselines/"
    "plots/{sha}/{measurement}.png"
)


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

    def status(self, alpha: float = ALPHA) -> str:
        """One of 'regression', 'improvement', 'compatible'."""
        if self.p_compat >= alpha:
            return "compatible"
        return "regression" if self.delta_chi2 > 0 else "improvement"


_EMOJI = {"regression": "🚩", "improvement": "⭐", "compatible": "✅"}
# Flagged rows first (regression, then improvement), compatible last.
_SORT_RANK = {"regression": 0, "improvement": 1, "compatible": 2}


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

    # -- derived quantities ---------------------------------------------------

    def p_overall(self) -> float:
        return bonferroni([r.p_compat for r in self.results])

    def n_flagged(self) -> int:
        return sum(1 for r in self.results if r.status(self.alpha) != "compatible")

    def overall_ok(self) -> bool:
        po = self.p_overall()
        return not (po == po and po < self.alpha)  # NaN-safe: ok if not < alpha

    def _plot_url(self, r: MeasurementResult) -> Optional[str]:
        if not r.plot:
            return None
        return RAW_URL_TEMPLATE.format(repo=self.repo, sha=self.feature_sha,
                                       measurement=r.name)

    # -- rendering ------------------------------------------------------------

    def _row(self, r: MeasurementResult) -> str:
        st = r.status(self.alpha)
        emoji = _EMOJI[st]
        url = self._plot_url(r)
        name_cell = f"[{r.name}]({url})" if url else r.name
        return (f"| {name_cell} | {r.ndof} | {r.chi2_ndof_main:.2f} | "
                f"{r.chi2_ndof_pr:.2f} | {r.delta_chi2:+.1f} | {r.p_compat:.3g} | "
                f"{r.p_data:.3g} | {emoji} |")

    @staticmethod
    def _table_header() -> str:
        return ("| Measurement | ndof | χ²/ndof (main) | χ²/ndof (PR) | Δχ² | "
                "p_compat | p (PR vs data) | |\n"
                "|---|---|---|---|---|---|---|---|")

    def to_markdown(self) -> str:
        po = self.p_overall()
        n = len(self.results)
        ordered = sorted(self.results,
                         key=lambda r: (_SORT_RANK[r.status(self.alpha)], -abs(r.delta_chi2)))
        flagged = [r for r in ordered if r.status(self.alpha) != "compatible"]
        compatible = [r for r in ordered if r.status(self.alpha) == "compatible"]

        verdict = ("✅ no significant change" if self.overall_ok()
                   else "⚠️ significant change")
        lines: List[str] = [COMMENT_MARKER, "## 🔬 Physics validation (NUISANCE3)", ""]
        lines.append(
            f"**Overall compatibility (Bonferroni, N={n}): p = {po:.3g} → {verdict}**")
        lines.append("")
        lines.append(
            f"NUISANCE3 `{self.nuisance_version}` · seed `{self.seed}` · "
            f"{self.events_per_measurement:,} events/measurement "
            f"· feature `{self.feature_sha[:8]}`")
        for extra in self.extra_header:
            lines.append(extra)
        lines.append("")

        # Flagged rows (with inline plot thumbnails) shown up front.
        lines.append(self._table_header())
        for r in flagged:
            lines.append(self._row(r))
        if not flagged:
            lines.append("| _all measurements compatible_ | | | | | | | ✅ |")
        lines.append("")

        if flagged:
            lines.append("### Flagged distributions")
            for r in flagged:
                url = self._plot_url(r)
                if url:
                    lines.append(f"**{r.name}** — {_EMOJI[r.status(self.alpha)]} "
                                 f"{r.status(self.alpha)}")
                    lines.append(f"![{r.name}]({url})")
            lines.append("")

        # Compatible rows collapsed to keep the comment scannable.
        if compatible:
            lines.append("<details><summary>"
                         f"{len(compatible)} compatible measurements</summary>\n")
            lines.append(self._table_header())
            for r in compatible:
                lines.append(self._row(r))
            lines.append("\n</details>")
            lines.append("")

        # Legend + multiple-comparison footer.
        thr = bonferroni_threshold(n, self.alpha)
        expected_false = self.alpha * n
        lines.append(
            "Legend: 🚩 regression (p_compat < {a}, Δχ² > 0) · "
            "⭐ significant improvement (p_compat < {a}, Δχ² < 0) · "
            "✅ compatible (p_compat ≥ {a}). Flag driven only by p_compat; Δχ² sign "
            "labels direction.".format(a=self.alpha))
        lines.append("")
        lines.append(
            f"_At uncorrected α={self.alpha} across N={n} measurements, "
            f"~{expected_false:.1f} false flags are expected by chance; the "
            f"Bonferroni per-measurement threshold is α/N = {thr:.4g}._")
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
        """🎯 reference · 🧪 null control · 🚩 worse than the null · ✅ otherwise."""
        if s.is_null_control:
            return "🧪"
        p_thr, s_thr = self.thresholds()
        return "🚩" if (s.p_overall < p_thr or s.p_shape_overall < s_thr) else "✅"

    def thresholds(self) -> tuple:
        """Flagging thresholds for (p_compat, p_shape), floored by the null control.

        The null control is the reference scheme rerun at a different seed, so its
        p-values are drawn from the null. Whatever it scores is what an *identical*
        configuration costs, and no real scheme should be called out for doing at
        least as well. The control can therefore only make the test more
        conservative — ``min`` with ``alpha`` — never less, so a control that happens
        to land at p ≈ 1 leaves the nominal threshold untouched.
        """
        for s in self.summaries:
            if s.is_null_control:
                p, q = s.p_overall, s.p_shape_overall
                return (min(self.alpha, p if p == p else self.alpha),
                        min(self.alpha, q if q == q else self.alpha))
        return (self.alpha, self.alpha)

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
        n_var = len([s for s in self.summaries if not s.is_reference])
        verdict = ("✅ every scheme reproduces the reference"
                   if not biased else
                   f"⚠️ {len(biased)} of {n_var} scheme(s) differ from the reference")

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

        flagged_rows = [r for r in self.rows if r.p_compat < self.alpha]
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
        null = next((x for x in self.summaries if x.is_null_control), None)
        if null is not None:
            p_thr, s_thr = self.thresholds()
            lines.append(
                f"🧪 **null control** — the reference scheme rerun at a different "
                f"seed, so its rows are drawn from the null hypothesis. It scores "
                f"p = {self._fmt(null.p_overall)} (shape "
                f"{self._fmt(null.p_shape_overall)}), which is what two *identical* "
                f"configurations cost: the bootstrap covariance is estimated inside a "
                f"single run and so does not carry the generator's run-to-run scatter "
                f"in the overall normalisation. Thresholds are floored at that value "
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


def _selftest() -> int:
    results = [
        MeasurementResult("MINERvA_CC0pi_Tp", 38, 1.11, 1.10, -0.4, 0.62, 0.31,
                          plot="MINERvA_CC0pi_Tp.png"),
        MeasurementResult("MiniBooNE_CC1pip_Q2", 40, 1.38, 1.53, +6.1, 0.013, 0.04,
                          plot="MiniBooNE_CC1pip_Q2.png"),
        MeasurementResult("T2K_CC0pi_cosTheta", 58, 0.98, 0.79, -11.0, 0.008, 0.71,
                          plot="T2K_CC0pi_cosTheta.png"),
    ]
    rep = Report(results, feature_sha="abcdef1234567890", nuisance_version="v3.0.1",
                 seed=42, events_per_measurement=500000)
    md = rep.to_markdown()
    summary = rep.to_summary_dict()

    checks = {
        "marker present": COMMENT_MARKER in md,
        "bonferroni p = min(1, 3*0.008)=0.024": abs(summary["p_overall"] - 0.024) < 1e-9,
        "overall flagged": summary["overall_ok"] is False,
        "two flagged rows": summary["n_flagged"] == 2,
        "regression labelled": next(m for m in summary["measurements"]
                                    if m["name"].startswith("MiniBooNE"))["status"] == "regression",
        "improvement labelled": next(m for m in summary["measurements"]
                                     if m["name"].startswith("T2K"))["status"] == "improvement",
        "compatible collapsed": "<details>" in md,
        "raw url embedded": "raw.githubusercontent.com" in md,
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
