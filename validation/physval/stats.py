#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Statistics for the Achilles physics validation (no NUISANCE/Achilles dependency).

* ``trial_covariance``   MC covariance of a prediction from the generator's trial counts
* ``compatibility``      main-vs-PR correlated chi-square (``p_compat``)
* ``goodness_of_fit``    prediction vs published data
* ``benjamini_hochberg``  multiple-comparison correction
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
    acceptance. ``response`` (the sample's conversion to published units) is applied as
    R S and R C R^T.
    """
    if not n_nonzero_trials or n_nonzero_trials <= 0:
        raise ValueError(f"n_nonzero_trials must be positive, got {n_nonzero_trials!r}")
    S = weighted_histogram(bin_index, np.asarray(weights, dtype=float), nbins)
    Q = weighted_histogram(bin_index, np.asarray(weights, dtype=float) ** 2, nbins)
    cov = np.diag(Q) - np.outer(S, S) / float(n_nonzero_trials)
    if rel_xsec_err:
        cov = cov + float(rel_xsec_err) ** 2 * np.outer(S, S)
    if response is not None:
        S = response @ S
        cov = response @ cov @ response.T
    return Prediction(values=S, covariance=cov)


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


def goodness_of_fit(pred: Prediction, data: np.ndarray,
                    data_cov: np.ndarray) -> ChiSquareResult:
    """Prediction vs data, with the MC covariance added to the data's."""
    delta = pred.values - np.asarray(data, dtype=float)
    return _chi2(delta, np.asarray(data_cov, dtype=float) + pred.covariance,
                 delta.shape[0])


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
