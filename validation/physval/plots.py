#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Publication-style overlay plots for physval.

One PNG per measurement: data, main and this branch with their 1-sigma bands, plus a
ratio-to-data panel.
"""

from __future__ import annotations

from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")  # headless CI
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito (colour-blind safe).
C_DATA = "#333333"
C_MAIN = "#0072B2"
C_FEATURE = "#D55E00"

_STYLE = {
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "semibold",
    "axes.labelsize": 10,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#444444",
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "lines.linewidth": 1.4,
}


def _errors(covariance: np.ndarray) -> np.ndarray:
    """Per-bin 1-sigma from a covariance diagonal."""
    return np.sqrt(np.clip(np.diag(np.asarray(covariance, dtype=float)), 0.0, None))


def _band(ax, edges: np.ndarray, values: np.ndarray, err: np.ndarray, color: str):
    if np.any(err):
        ax.stairs(values + err, edges, baseline=values - err, fill=True,
                  color=color, alpha=0.18, linewidth=0)


def plot_measurement(path: str, name: str, *,
                     data: np.ndarray,
                     data_cov: np.ndarray,
                     feature: np.ndarray,
                     feature_cov: np.ndarray,
                     main: Optional[np.ndarray] = None,
                     main_cov: Optional[np.ndarray] = None,
                     edges: Optional[Sequence[float]] = None,
                     subtitle: str = "") -> None:
    """Render the data/main/feature overlay for one measurement to ``path``."""
    data = np.asarray(data, dtype=float)
    nbins = data.shape[0]
    edges = (np.asarray(edges, dtype=float) if edges is not None
             else np.arange(nbins + 1, dtype=float))
    centers = 0.5 * (edges[:-1] + edges[1:])
    d_err = _errors(data_cov)
    series = [(np.asarray(feature, dtype=float), _errors(feature_cov), C_FEATURE, "this branch")]
    if main is not None:
        series.insert(0, (np.asarray(main, dtype=float), _errors(main_cov), C_MAIN, "main"))

    with plt.rc_context(_STYLE):
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(6.0, 5.2), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07})

        ax.errorbar(centers, data, yerr=d_err, fmt="o", color=C_DATA,
                    markersize=3.4, elinewidth=1.0, capsize=0, label="Data",
                    zorder=5)
        for values, err, color, label in series:
            _band(ax, edges, values, err, color)
            ax.stairs(values, edges, color=color, label=label, linewidth=1.4)

        ax.set_ylabel("d$\\sigma$/dx")
        # Title sits above the subtitle line; pad leaves room for both.
        ax.set_title(name, loc="left", pad=18 if subtitle else 6)
        if subtitle:
            ax.text(0.0, 1.012, subtitle, transform=ax.transAxes, fontsize=8,
                    color="#666666", va="bottom")
        # Keep the ~1e-38 offset text clear of the subtitle.
        offset = ax.yaxis.get_offset_text()
        offset.set_horizontalalignment("right")
        offset.set_position((1.0, 1.0))
        ax.legend(loc="best")
        ax.margins(x=0.01)

        # Ratio to data; bins with zero data are left blank rather than inf.
        safe = np.where(data != 0.0, data, np.nan)
        rax.axhline(1.0, color=C_DATA, linewidth=1.0, zorder=4)
        with np.errstate(invalid="ignore", divide="ignore"):
            _band(rax, edges, np.ones(nbins), d_err / np.abs(safe), C_DATA)
            for values, err, color, _ in series:
                ratio = values / safe
                _band(rax, edges, ratio, err / np.abs(safe), color)
                rax.stairs(ratio, edges, color=color, linewidth=1.4)

        rax.set_ylabel("MC / data")
        rax.set_xlabel("x")
        rax.set_ylim(0.5, 1.5)
        rax.margins(x=0.01)

        fig.savefig(path)
        plt.close(fig)
