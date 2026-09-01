#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Publication-style overlay plots for physval.

One PNG per measurement: data with its uncertainty band, the stored `main`
prediction and this branch's prediction on a shared axis, plus a ratio-to-data
panel underneath. Prediction bands are the bootstrap MC uncertainty (the sqrt of
the covariance diagonal).
"""

from __future__ import annotations

from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")  # headless CI
import matplotlib.pyplot as plt
import numpy as np

# Colour-blind-safe (Okabe-Ito); data stays neutral so the two predictions read as
# the signal. Kept in one place so every measurement's plot looks identical.
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


def _errors(covariance: Optional[np.ndarray], n: int) -> np.ndarray:
    """Per-bin 1-sigma from a covariance diagonal (zeros if absent)."""
    if covariance is None:
        return np.zeros(n)
    return np.sqrt(np.clip(np.diag(np.asarray(covariance, dtype=float)), 0.0, None))


def _step(ax, edges: np.ndarray, values: np.ndarray, color: str, label: str):
    """Histogram outline plus a shaded 1-sigma band, drawn bin-edge to bin-edge."""
    ax.stairs(values, edges, color=color, label=label, linewidth=1.4)


def _band(ax, edges: np.ndarray, values: np.ndarray, err: np.ndarray, color: str):
    if not np.any(err):
        return
    ax.stairs(values + err, edges, baseline=values - err, fill=True,
              color=color, alpha=0.18, linewidth=0)


def plot_measurement(path: str, name: str, *,
                     data: np.ndarray,
                     data_cov: Optional[np.ndarray] = None,
                     main: Optional[np.ndarray] = None,
                     main_cov: Optional[np.ndarray] = None,
                     feature: Optional[np.ndarray] = None,
                     feature_cov: Optional[np.ndarray] = None,
                     edges: Optional[Sequence[float]] = None,
                     xlabel: str = "bin",
                     ylabel: str = "d$\\sigma$/dx",
                     subtitle: str = "") -> str:
    """Render the data/main/feature overlay for one measurement to ``path``."""
    data = np.asarray(data, dtype=float)
    nbins = data.shape[0]
    edges = (np.asarray(edges, dtype=float) if edges is not None
             else np.arange(nbins + 1, dtype=float))
    centers = 0.5 * (edges[:-1] + edges[1:])

    d_err = _errors(data_cov, nbins)

    with plt.rc_context(_STYLE):
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(6.0, 5.2), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.07})

        # -- main panel ------------------------------------------------------
        ax.errorbar(centers, data, yerr=d_err, fmt="o", color=C_DATA,
                    markersize=3.4, elinewidth=1.0, capsize=0, label="Data",
                    zorder=5)

        for values, cov, color, label in (
                (main, main_cov, C_MAIN, "main"),
                (feature, feature_cov, C_FEATURE, "this branch")):
            if values is None:
                continue
            values = np.asarray(values, dtype=float)
            _band(ax, edges, values, _errors(cov, nbins), color)
            _step(ax, edges, values, color, label)

        ax.set_ylabel(ylabel)
        # Title sits above the subtitle line; pad leaves room for both.
        ax.set_title(name, loc="left", pad=18 if subtitle else 6)
        if subtitle:
            ax.text(0.0, 1.012, subtitle, transform=ax.transAxes, fontsize=8,
                    color="#666666", va="bottom")
        # Cross sections are ~1e-38, so matplotlib parks a shared exponent at the top
        # left -- exactly where the subtitle goes. Move it to the right instead.
        offset = ax.yaxis.get_offset_text()
        offset.set_horizontalalignment("right")
        offset.set_position((1.0, 1.0))
        ax.legend(loc="best")
        ax.margins(x=0.01)

        # -- ratio panel -----------------------------------------------------
        # Ratio to data; bins with zero data are left blank rather than inf.
        safe = np.where(data != 0.0, data, np.nan)
        rax.axhline(1.0, color=C_DATA, linewidth=1.0, zorder=4)
        with np.errstate(invalid="ignore", divide="ignore"):
            _band(rax, edges, np.ones(nbins), d_err / np.abs(safe), C_DATA)
            for values, cov, color in ((main, main_cov, C_MAIN),
                                       (feature, feature_cov, C_FEATURE)):
                if values is None:
                    continue
                values = np.asarray(values, dtype=float)
                ratio = values / safe
                err = _errors(cov, nbins) / np.abs(safe)
                _band(rax, edges, ratio, err, color)
                _step(rax, edges, ratio, color, "")

        rax.set_ylabel("MC / data")
        rax.set_xlabel(xlabel)
        rax.set_ylim(0.5, 1.5)
        rax.margins(x=0.01)

        fig.savefig(path)
        plt.close(fig)
    return path


def _selftest() -> int:
    """Render a plot to a temp dir and assert a non-trivial PNG came out."""
    import os
    import tempfile

    rng = np.random.default_rng(0)
    nbins = 14
    edges = np.linspace(0.0, 2.0, nbins + 1)
    truth = 40 * np.exp(-0.5 * ((0.5 * (edges[:-1] + edges[1:]) - 0.8) / 0.35) ** 2)
    data = truth * rng.normal(1.0, 0.04, nbins)
    dcov = np.diag((0.06 * truth) ** 2)
    main = truth * rng.normal(1.0, 0.02, nbins)
    feat = truth * 1.08 * rng.normal(1.0, 0.02, nbins)
    mcov = np.diag((0.03 * truth) ** 2)

    checks = {}
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "SELFTEST.png")
        plot_measurement(p, "SELFTEST_sample_1Dx_nu", data=data, data_cov=dcov,
                         main=main, main_cov=mcov, feature=feat, feature_cov=mcov,
                         edges=edges, xlabel="x [GeV]", subtitle="seed 42")
        checks["png written"] = os.path.exists(p)
        checks["png non-trivial"] = os.path.getsize(p) > 5000
        with open(p, "rb") as fh:
            checks["png magic"] = fh.read(8) == b"\x89PNG\r\n\x1a\n"

        # Degenerate inputs must not raise: no edges, no covariances, zero-data bins.
        p2 = os.path.join(tmp, "MINIMAL.png")
        z = data.copy()
        z[3] = 0.0
        plot_measurement(p2, "MINIMAL", data=z, main=main, feature=feat)
        checks["minimal png written"] = os.path.exists(p2)

    for k, ok in checks.items():
        print(f"[{'ok' if ok else 'FAIL'}] {k}")
    ok = all(checks.values())
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()
    raise SystemExit(_selftest() if args.selftest else 0)
