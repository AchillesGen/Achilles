#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Statistics for the Achilles physics validation (no NUISANCE/Achilles dependency).

* ``trial_covariance``   MC covariance of a prediction from the generator's trial counts
* ``compatibility``      main-vs-PR correlated chi-square (``p_compat``)
* ``shape_compatibility`` the same with the normalisation fitted out
* ``goodness_of_fit``    prediction vs published data
* ``benjamini_hochberg`` / ``bonferroni``  multiple-comparison corrections
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy import stats


@dataclass
class Prediction:
    """A binned prediction and its MC covariance."""

    values: np.ndarray       # shape (nbins,)
    covariance: np.ndarray   # shape (nbins, nbins)

    @property
    def nbins(self) -> int:
        return self.values.shape[0]

    def to_dict(self) -> dict:
        return {"values": self.values.tolist(), "covariance": self.covariance.tolist()}

    @classmethod
    def from_dict(cls, d: dict) -> "Prediction":
        return cls(values=np.asarray(d["values"], dtype=float),
                   covariance=np.asarray(d["covariance"], dtype=float))


def weighted_histogram(bin_index: np.ndarray, weights: np.ndarray,
                       nbins: int) -> np.ndarray:
    return np.bincount(bin_index, weights=weights, minlength=nbins)[:nbins]


def trial_covariance(bin_index: np.ndarray, weights: np.ndarray, nbins: int,
                     n_nonzero_trials: int, rel_xsec_err: float = 0.0,
                     response: Optional[np.ndarray] = None) -> Prediction:
    """Analytic MC covariance from the generator's own trial counts.

    The file holds exactly N non-zero trials, so binning them is multinomial over all
    N, including those the selection rejects. With ``S_k`` the weight sum of bin k:

        C_jk = delta_jk * sum_{i in k} w_i^2 - S_j S_k / N  +  rel_xsec_err^2 S_j S_k

    Resampling only the selected events instead fixes their count and drops the
    acceptance fluctuation, which under-covers the normalisation by ~3x at 19%
    acceptance. ``response`` (a Wiener-SVD A_C) is applied as R S and R C R^T.
    """
    if not n_nonzero_trials or n_nonzero_trials <= 0:
        raise ValueError(f"n_nonzero_trials must be positive, got {n_nonzero_trials!r}")
    S = weighted_histogram(bin_index, np.asarray(weights, dtype=float), nbins)
    Q = weighted_histogram(bin_index, np.asarray(weights, dtype=float) ** 2, nbins)
    cov = np.diag(Q) - np.outer(S, S) / float(n_nonzero_trials)
    if rel_xsec_err:
        cov = cov + float(rel_xsec_err) ** 2 * np.outer(S, S)
    if response is not None:
        response = np.asarray(response, dtype=float)
        if response.shape != (nbins, nbins):
            raise ValueError(f"response matrix is {response.shape}, expected "
                             f"({nbins}, {nbins})")
        S = response @ S
        cov = response @ cov @ response.T
    return Prediction(values=S, covariance=np.atleast_2d(cov))


@dataclass
class ChiSquareResult:
    chi2: float
    ndof: int
    pvalue: float

    @property
    def chi2_per_ndof(self) -> float:
        return self.chi2 / self.ndof if self.ndof > 0 else float("nan")


def _chi2(delta: np.ndarray, cov: np.ndarray, ndof: int) -> ChiSquareResult:
    # pinv: an empty bin in both predictions has zero variance.
    chi2 = float(delta @ np.linalg.pinv(cov) @ delta)
    pvalue = float(stats.chi2.sf(chi2, ndof)) if ndof > 0 else float("nan")
    return ChiSquareResult(chi2=chi2, ndof=ndof, pvalue=pvalue)


def compatibility(main: Prediction, feature: Prediction) -> ChiSquareResult:
    """(h_f - h_m)^T (C_m + C_f)^-1 (h_f - h_m); small p means a real change."""
    delta = feature.values - main.values
    return _chi2(delta, main.covariance + feature.covariance, delta.shape[0])


def shape_compatibility(main: Prediction, feature: Prediction) -> ChiSquareResult:
    """``compatibility`` with the feature rescaled onto main, one dof spent on the scale.

    The scale is fitted in the covariance metric: with unequal bin widths a plain-sum
    ratio leaves a residual along the direction the chi-square is most sensitive to.
    """
    cov = main.covariance + feature.covariance
    inv = np.linalg.pinv(cov)
    denom = float(feature.values @ inv @ feature.values)
    if denom <= 0.0:
        return compatibility(main, feature)
    scale = float(feature.values @ inv @ main.values) / denom
    delta = feature.values * scale - main.values
    return _chi2(delta, main.covariance + feature.covariance * scale ** 2,
                 max(delta.shape[0] - 1, 1))


def goodness_of_fit(pred: Prediction, data: np.ndarray,
                    data_cov: np.ndarray) -> ChiSquareResult:
    """Prediction vs data, with the MC covariance added to the data's."""
    delta = pred.values - np.asarray(data, dtype=float)
    return _chi2(delta, np.asarray(data_cov, dtype=float) + pred.covariance,
                 delta.shape[0])


def bonferroni(pvalues: Sequence[float]) -> float:
    """Family-wise p-value ``min(1, N * min p)``."""
    p = np.asarray(pvalues, dtype=float)
    return float(min(1.0, p.size * np.min(p))) if p.size else float("nan")


def benjamini_hochberg(pvalues: Sequence[float]) -> np.ndarray:
    """BH q-values: flagging q < alpha keeps the expected false-flag fraction at alpha."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if not n:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def _runs(rng: np.random.Generator, n_runs: int, shift: float = 0.0, *,
          nbins: int = 25, n_nonzero: int = 50_000, acceptance: float = 0.2,
          rel_xsec_err: float = 0.0015):
    """Runs of one configuration as a generator delivers them, as Predictions.

    N non-zero trials fixed, Binomial(N, acceptance) selected, a per-run xsec error,
    and non-uniform bin widths as published binnings have.
    """
    shape = np.linspace(2.0, 1.0, nbins) + shift * np.linspace(-1.0, 1.0, nbins)
    probs = np.append(shape / shape.sum() * acceptance, 1.0 - acceptance)
    widths = np.where(np.arange(nbins) < nbins // 2, 0.1, 1.0)
    for _ in range(n_runs):
        counts = rng.multinomial(n_nonzero, probs)[:nbins]
        xsec = 1.0 + rel_xsec_err * rng.normal()
        bin_index = np.repeat(np.arange(nbins), counts)
        weights = (xsec / n_nonzero) / widths[bin_index]
        yield trial_covariance(bin_index, weights, nbins, n_nonzero,
                               rel_xsec_err=rel_xsec_err)


def _selftest() -> int:
    rng = np.random.default_rng(98765)
    alpha, pairs = 0.05, 400

    # Calibration: two runs of one configuration must give uniform p_compat.
    runs = list(_runs(rng, 2 * pairs))
    pvals = np.array([compatibility(runs[2 * i], runs[2 * i + 1]).pvalue
                      for i in range(pairs)])
    rate = float(np.mean(pvals < alpha))
    ks_p = float(stats.kstest(pvals, "uniform").pvalue)
    tol = 3 * np.sqrt(alpha * (1 - alpha) / pairs)
    print(f"[null] flag rate={rate:.3f} (target {alpha}±{tol:.3f}), KS p={ks_p:.3f}")

    # Power: a real shape change must be seen.
    a, b = next(_runs(rng, 1)), next(_runs(rng, 1, shift=0.5))
    p_shift = compatibility(a, b).pvalue
    print(f"[shift] p_compat={p_shift:.2e} (expect << 0.05)")

    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
    print(f"[bh] q={np.round(q, 4).tolist()} (expect [0.04, 0.0533, 0.0533, 0.5])")

    ok = (abs(rate - alpha) <= tol and ks_p > 0.01 and p_shift < 0.05
          and np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.5])
          and abs(bonferroni([0.6, 0.013]) - 0.026) < 1e-12)
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    print(__doc__)
