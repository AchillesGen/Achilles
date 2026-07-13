#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Statistics for the Achilles distribution-based physics validation.

This module is deliberately free of any NUISANCE / Achilles dependency so that it
can be unit-tested in isolation (see ``--selftest``).  It implements the three
quantities the physval summary table is built from:

* ``bootstrap_covariance`` - the bin-by-bin MC covariance of a prediction, obtained
  by resampling its weighted events with replacement.  The stored ``main`` baseline
  keeps this covariance so a pull request never has to regenerate ``main``.
* ``compatibility`` - the ``main`` vs feature-branch compatibility p-value
  (``p_compat``) from a correlated chi-square using ``C_main + C_feature``.  This is
  the statistic that drives the flag (flag when ``p_compat < 0.05``).
* ``goodness_of_fit`` - the chi-square / p-value of a prediction against the
  experimental data, reported per row for context only.
* ``bonferroni`` - the family-wise overall p-value across all measurements.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# Histograms and bootstrap covariance
# ---------------------------------------------------------------------------

def weighted_histogram(bin_index: np.ndarray, weights: np.ndarray,
                       nbins: int) -> np.ndarray:
    """Sum event ``weights`` into ``nbins`` bins indexed by ``bin_index``."""
    return np.bincount(bin_index, weights=weights, minlength=nbins)[:nbins]


@dataclass
class Prediction:
    """A histogrammed prediction together with its bootstrap MC covariance."""

    values: np.ndarray            # central histogram, shape (nbins,)
    covariance: np.ndarray        # MC covariance, shape (nbins, nbins)
    n_boot: Optional[int] = None  # bootstrap replicas used (for Hartlap correction)

    @property
    def nbins(self) -> int:
        return self.values.shape[0]

    def to_dict(self) -> dict:
        return {
            "values": self.values.tolist(),
            "covariance": self.covariance.tolist(),
            "n_boot": self.n_boot,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Prediction":
        return cls(values=np.asarray(d["values"], dtype=float),
                   covariance=np.asarray(d["covariance"], dtype=float),
                   n_boot=d.get("n_boot"))


def bootstrap_covariance(bin_index: np.ndarray, weights: np.ndarray, nbins: int,
                         n_boot: int = 200,
                         rng: Optional[np.random.Generator] = None) -> Prediction:
    """Central histogram and bootstrap covariance for a weighted event sample.

    The events (``bin_index``/``weights`` pairs) are resampled with replacement
    ``n_boot`` times; the covariance is estimated across the resulting ensemble of
    histograms.  This captures the prediction's MC statistical uncertainty without
    any re-generation, and works for weighted (including negative-weight) events.
    """
    if rng is None:
        rng = np.random.default_rng()
    bin_index = np.asarray(bin_index)
    weights = np.asarray(weights, dtype=float)
    n_events = bin_index.shape[0]

    central = weighted_histogram(bin_index, weights, nbins)

    ensemble = np.empty((n_boot, nbins), dtype=float)
    for b in range(n_boot):
        pick = rng.integers(0, n_events, size=n_events)
        ensemble[b] = weighted_histogram(bin_index[pick], weights[pick], nbins)

    # rowvar=False: each column is a bin, each row a bootstrap replica.
    cov = np.cov(ensemble, rowvar=False)
    cov = np.atleast_2d(cov)
    return Prediction(values=central, covariance=cov, n_boot=n_boot)


# ---------------------------------------------------------------------------
# Chi-square helpers
# ---------------------------------------------------------------------------

def hartlap_factor(n_samples: int, n_bins: int) -> float:
    """Debias the inverse of a covariance estimated from ``n_samples`` draws.

    The inverse of a sample covariance is biased high (inverse-Wishart), which
    inflates chi-square. Hartlap et al. (2007) correct it by
    ``(n_samples - n_bins - 2) / (n_samples - 1)``. Requires ``n_samples > n_bins + 2``;
    returns 1.0 (no correction) when the sample size / dimension are unknown or the
    correction would be non-positive (caller should also enforce n_boot >> n_bins).
    """
    if n_samples is None or n_samples <= n_bins + 2:
        return 1.0
    return (n_samples - n_bins - 2) / (n_samples - 1)


def _chi2_quadratic_form(delta: np.ndarray, cov: np.ndarray,
                         n_samples: Optional[int] = None) -> float:
    """delta^T cov^{-1} delta, robust to (near-)singular covariance.

    When ``n_samples`` (the effective number of bootstrap replicas behind ``cov``)
    is given, the inverse is Hartlap-corrected to remove the finite-sample bias.
    """
    # A pseudo-inverse degrades gracefully when a bin has zero variance (e.g. an
    # empty bin in both predictions) instead of raising / returning inf.
    inv = np.linalg.pinv(cov)
    inv = inv * hartlap_factor(n_samples, delta.shape[0])
    return float(delta @ inv @ delta)


def _combine_n_boot(a: Optional[int], b: Optional[int]) -> Optional[int]:
    """Effective replica count for a sum of two bootstrap covariances."""
    vals = [x for x in (a, b) if x is not None]
    return min(vals) if vals else None


@dataclass
class ChiSquareResult:
    chi2: float
    ndof: int
    pvalue: float

    @property
    def chi2_per_ndof(self) -> float:
        return self.chi2 / self.ndof if self.ndof > 0 else float("nan")


def compatibility(main: Prediction, feature: Prediction) -> ChiSquareResult:
    """Compatibility of two predictions via a correlated chi-square.

    ``chi2 = (h_feature - h_main)^T (C_main + C_feature)^{-1} (h_feature - h_main)``

    A large ``p_compat`` means the two predictions are statistically the same given
    their MC uncertainties; a small ``p_compat`` (< 0.05) flags a significant change.
    """
    delta = feature.values - main.values
    cov = main.covariance + feature.covariance
    # Both terms are bootstrap sample covariances; the sum's effective replica
    # count is bounded below by the smaller ensemble. Use it for a conservative
    # Hartlap debias of the inverse.
    n_eff = _combine_n_boot(main.n_boot, feature.n_boot)
    chi2 = _chi2_quadratic_form(delta, cov, n_samples=n_eff)
    ndof = delta.shape[0]
    pvalue = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else float("nan")
    return ChiSquareResult(chi2=chi2, ndof=ndof, pvalue=pvalue)


def goodness_of_fit(pred: Prediction, data: np.ndarray,
                    data_cov: np.ndarray) -> ChiSquareResult:
    """Chi-square of a prediction against experimental ``data``.

    The prediction's own MC covariance is added to the data covariance so the MC
    statistical uncertainty is not double counted as agreement.
    """
    delta = pred.values - np.asarray(data, dtype=float)
    cov = np.asarray(data_cov, dtype=float) + pred.covariance
    # Only the prediction term is a bootstrap estimate; debias with its replica
    # count (the exact data covariance makes this mildly conservative).
    chi2 = _chi2_quadratic_form(delta, cov, n_samples=pred.n_boot)
    ndof = delta.shape[0]
    pvalue = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else float("nan")
    return ChiSquareResult(chi2=chi2, ndof=ndof, pvalue=pvalue)


def empirical_pvalue(observed_chi2: float, null_chi2: Sequence[float]) -> float:
    """Non-Gaussian fallback: p from a bootstrap null distribution of chi-square.

    Fraction of the null ensemble at least as extreme as ``observed_chi2``.
    """
    null = np.asarray(null_chi2, dtype=float)
    # +1 in numerator and denominator: never report an impossible p == 0.
    return float((np.sum(null >= observed_chi2) + 1) / (null.size + 1))


# ---------------------------------------------------------------------------
# Multiple-comparison combination
# ---------------------------------------------------------------------------

def bonferroni(pvalues: Sequence[float]) -> float:
    """Family-wise overall p-value: ``min(1, N * min_i p_i)``.

    Standard in particle physics for the look-elsewhere / multiple-testing effect.
    Flag overall concern when this is < 0.05.
    """
    p = np.asarray(pvalues, dtype=float)
    if p.size == 0:
        return float("nan")
    return float(min(1.0, p.size * np.min(p)))


def bonferroni_threshold(n: int, alpha: float = 0.05) -> float:
    """Per-measurement Bonferroni threshold ``alpha / N``."""
    return alpha / n if n > 0 else alpha


# ---------------------------------------------------------------------------
# Self-test: verifies bootstrap calibration (the core statistical assumption)
# ---------------------------------------------------------------------------

def _selftest() -> int:
    rng = np.random.default_rng(1234)
    nbins = 12
    edges = np.linspace(0.0, 1.0, nbins + 1)

    def sample_prediction(n_events: int, shift: float,
                          seed_rng: np.random.Generator) -> Prediction:
        """Draw a weighted event sample from a (optionally shifted) distribution."""
        x = np.clip(seed_rng.normal(0.5 + shift, 0.18, size=n_events), 0, 0.999)
        bin_index = np.digitize(x, edges) - 1
        weights = seed_rng.uniform(0.5, 1.5, size=n_events)
        return bootstrap_covariance(bin_index, weights, nbins, n_boot=150, rng=rng)

    # --- Calibration under the null: same distribution, independent samples. -----
    # p_compat should be ~uniform on [0, 1] and reject at ~alpha.
    alpha = 0.05
    trials = 300
    pvals = []
    for _ in range(trials):
        srng = np.random.default_rng(rng.integers(1 << 30))
        a = sample_prediction(4000, 0.0, srng)
        b = sample_prediction(4000, 0.0, srng)
        pvals.append(compatibility(a, b).pvalue)
    pvals = np.asarray(pvals)
    reject_rate = float(np.mean(pvals < alpha))
    mean_p = float(np.mean(pvals))
    ks_p = float(stats.kstest(pvals, "uniform").pvalue)

    print(f"[null] false-flag rate={reject_rate:.3f} (target ~{alpha}), "
          f"mean p={mean_p:.3f} (target ~0.5), uniformity KS p={ks_p:.3f}")
    ok_null = (0.01 <= reject_rate <= 0.12) and (0.40 <= mean_p <= 0.60) and (ks_p > 0.01)

    # --- Power under a real shift: distribution moved, expect small p_compat. -----
    srng = np.random.default_rng(7)
    a = sample_prediction(4000, 0.0, srng)
    b = sample_prediction(4000, 0.06, srng)
    shifted = compatibility(a, b)
    print(f"[shift] chi2/ndof={shifted.chi2_per_ndof:.2f} p_compat={shifted.pvalue:.2e} "
          "(expect << 0.05)")
    ok_power = shifted.pvalue < 0.05

    # --- Bonferroni sanity. -------------------------------------------------------
    p_overall = bonferroni([0.6, 0.4, 0.013, 0.2])
    print(f"[bonferroni] p_overall={p_overall:.3f} (== min(1, 4*0.013)=0.052)")
    ok_bonf = abs(p_overall - 0.052) < 1e-9

    ok = ok_null and ok_power and ok_bonf
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    print(__doc__)
