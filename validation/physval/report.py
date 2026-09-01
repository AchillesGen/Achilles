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
