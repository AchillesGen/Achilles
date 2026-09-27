#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Statistics for the Achilles distribution-based physics validation.

This module is deliberately free of any NUISANCE / Achilles dependency so that it
can be unit-tested in isolation (see ``--selftest``).  It implements the three
quantities the physval summary table is built from:

* ``trial_covariance`` - the bin-by-bin MC covariance of a prediction, built
  analytically from the generator's own trial counts and its cross-section
  uncertainty.  The stored ``main`` baseline keeps this covariance so a pull request
  never has to regenerate ``main``.
* ``bootstrap_covariance`` - the same thing by resampling the events with
  replacement.  Kept as a cross-check and as the fallback for event samples that
  carry no trial counts; see ``trial_covariance`` for why it cannot be used on a
  selected subsample.
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
                         rng: Optional[np.random.Generator] = None,
                         response: Optional[np.ndarray] = None) -> Prediction:
    """Central histogram and bootstrap covariance for a weighted event sample.

    The events (``bin_index``/``weights`` pairs) are resampled with replacement
    ``n_boot`` times; the covariance is estimated across the resulting ensemble of
    histograms.  This captures the prediction's MC statistical uncertainty without
    any re-generation, and works for weighted (including negative-weight) events.

    ``response`` is an (nbins, nbins) matrix applied to every histogram, for a
    measurement whose unfolded data is only comparable to A_C * prediction. Applying
    it per replica rather than to the central values alone carries it into the MC
    covariance as A cov A^T, which is what the chi-square then uses.
    """
    if rng is None:
        rng = np.random.default_rng()
    bin_index = np.asarray(bin_index)
    weights = np.asarray(weights, dtype=float)
    n_events = bin_index.shape[0]

    if response is not None:
        response = np.asarray(response, dtype=float)
        if response.shape != (nbins, nbins):
            raise ValueError(f"response matrix is {response.shape}, expected "
                             f"({nbins}, {nbins})")

    def histogram(idx, w):
        h = weighted_histogram(idx, w, nbins)
        return h if response is None else response @ h

    central = histogram(bin_index, weights)

    ensemble = np.empty((n_boot, nbins), dtype=float)
    for b in range(n_boot):
        pick = rng.integers(0, n_events, size=n_events)
        ensemble[b] = histogram(bin_index[pick], weights[pick])

    # rowvar=False: each column is a bin, each row a bootstrap replica.
    cov = np.cov(ensemble, rowvar=False)
    cov = np.atleast_2d(cov)
    return Prediction(values=central, covariance=cov, n_boot=n_boot)


def trial_covariance(bin_index: np.ndarray, weights: np.ndarray, nbins: int,
                     n_nonzero_trials: int, rel_xsec_err: float = 0.0,
                     response: Optional[np.ndarray] = None) -> Prediction:
    """Analytic MC covariance from the generator's own trial counts.

    A generator asked for N events runs trials until N of them have a non-zero
    weight, so the events in the file *are* those N non-zero trials and N is fixed by
    construction.  Binning them is therefore multinomial over all N -- including the
    ones a selection throws away -- which gives, with ``S_k`` the weight sum of bin
    ``k``,

        C_jk = delta_jk * sum_{i in k} w_i^2  -  S_j S_k / N            (a)

    The ``- S_j S_k / N`` term is the constraint that exactly N events were
    generated.  Applying that constraint to the *selected* events instead -- which is
    what resampling only those events does -- asserts that the number passing the
    selection is exact and so drops the acceptance fluctuation entirely.  At a 19%
    acceptance that understates the normalisation uncertainty by a factor ~3, which
    is enough to give two runs of one configuration ``p_compat`` ~ 1e-150 while every
    individual bin still agrees within its error.

    ``rel_xsec_err`` is the generator's own relative uncertainty on the total cross
    section (NuHepMC's GenCrossSection attribute carries it next to the trial
    counters).  It scales every bin together, so it enters fully correlated:

        C_jk += rel_xsec_err^2 * S_j S_k                                (b)

    Both terms are needed: (a) is the larger one and (b) is the only carrier of the
    integrator's own error.  (b) double counts slightly, since the cross-section
    estimate is built from the same trials, which costs ~3% in chi-square -- cheap
    insurance against the mode that (a) cannot see.

    ``response`` is applied as ``R S`` and ``R C R^T``, exactly rather than per
    replica.  The result carries ``n_boot=None``: this is not a finite-sample
    estimate, so no Hartlap correction applies to its inverse and there is no
    ``n_boot >> n_bins`` requirement to satisfy.
    """
    if not n_nonzero_trials or n_nonzero_trials <= 0:
        raise ValueError("n_nonzero_trials must be the positive number of non-zero "
                         f"trials the generator produced, got {n_nonzero_trials!r}")
    bin_index = np.asarray(bin_index)
    weights = np.asarray(weights, dtype=float)

    S = weighted_histogram(bin_index, weights, nbins)
    Q = weighted_histogram(bin_index, weights * weights, nbins)
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

    return Prediction(values=S, covariance=np.atleast_2d(cov), n_boot=None)


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


def _metric_scale(feature: np.ndarray, main: np.ndarray,
                  cov: np.ndarray) -> float:
    """Scale ``a`` minimising ``(a*feature - main)^T cov^-1 (a*feature - main)``.

    Falls back to the ratio of plain sums when the covariance carries no usable
    information in the feature's direction.
    """
    inv = np.linalg.pinv(cov)
    denom = float(feature @ inv @ feature)
    if denom > 0.0:
        return float(feature @ inv @ main) / denom
    total_feat = float(np.sum(feature))
    return float(np.sum(main)) / total_feat if total_feat else 1.0


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


def shape_compatibility(main: Prediction, feature: Prediction) -> ChiSquareResult:
    """Compatibility of two predictions' *shapes*, with the normalisation divided out.

    With ``trial_covariance`` the overall scale is covered by the test itself, so
    ``compatibility`` is the statistic to read and this is a diagnostic: it separates
    "the distribution moved" from "the two runs disagree on the total cross section".
    The feature is rescaled to the reference (its covariance with it, so the MC
    uncertainty is rescaled consistently) and one degree of freedom is given up for
    the scale that was fitted.

    The scale is fitted *in the covariance metric* rather than as a ratio of plain
    bin sums.  The direction a histogram's covariance constrains most weakly is the
    bin-content direction, which is not the plain sum once the bins have unequal
    widths (the prediction is divided by width), so a plain-sum rescale leaves a
    residual along exactly the direction the chi-square is most sensitive to -- it
    reports a broken shape for measurements whose only problem is normalisation, and
    the effect grows with how non-uniform the binning is.
    """
    total_feat = float(np.sum(feature.values))
    if not total_feat:
        return compatibility(main, feature)
    scale = _metric_scale(feature.values, main.values,
                          main.covariance + feature.covariance)
    scaled = Prediction(values=feature.values * scale,
                        covariance=feature.covariance * scale * scale,
                        n_boot=feature.n_boot)
    delta = scaled.values - main.values
    cov = main.covariance + scaled.covariance
    n_eff = _combine_n_boot(main.n_boot, feature.n_boot)
    chi2 = _chi2_quadratic_form(delta, cov, n_samples=n_eff)
    ndof = max(delta.shape[0] - 1, 1)  # one dof spent on the fitted scale
    pvalue = float(stats.chi2.sf(chi2, ndof))
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

def _null_control_runs(rng: np.random.Generator, n_runs: int, *,
                       nbins: int = 25, n_nonzero: int = 50_000,
                       acceptance: float = 0.2, rel_xsec_err: float = 0.0015):
    """Simulate ``n_runs`` runs of ONE configuration, the way a generator delivers them.

    The generator is asked for ``n_nonzero`` events and runs trials until that many
    have a non-zero weight, so that count is fixed and the number passing a selection
    is ``Binomial(n_nonzero, acceptance)``.  The prediction is scaled by the run's own
    total cross section, which carries ``rel_xsec_err``.  Bin widths are deliberately
    non-uniform, as published binnings are: the per-event weight then varies from bin
    to bin, which is what gives a fixed-count resample any spread in the total at all.

    Yields ``(bin_index, weights)`` per run -- everything an estimator is given.
    """
    shape = np.linspace(2.0, 1.0, nbins)
    probs = np.append(shape / shape.sum() * acceptance, 1.0 - acceptance)
    widths = np.where(np.arange(nbins) < nbins // 2, 0.1, 1.0)
    for _ in range(n_runs):
        counts = rng.multinomial(n_nonzero, probs)[:nbins]
        xsec = 1.0 + rel_xsec_err * rng.normal()
        bin_index = np.repeat(np.arange(nbins), counts)
        weights = (xsec / n_nonzero) / widths[bin_index]
        yield bin_index, weights


def _selftest_null_control(rng: np.random.Generator, alpha: float = 0.05):
    """Regression test for the estimator the summary table is built from.

    Two runs of one configuration differ only by seed, so their ``p_compat`` must be
    uniform and flag at ``alpha``.  This also pins the failure it was written for:
    resampling the *selected* events at fixed count holds their number exact, which
    denies the acceptance fluctuation and under-covers the normalisation -- the
    ``covers`` ratio below is ~1 for the trial-count covariance and far under 1 for
    the resample, and that is what turns identical configurations into p ~ 1e-150.
    """
    nbins, n_nonzero, acceptance, rel_err = 25, 50_000, 0.2, 0.0015
    kwargs = dict(nbins=nbins, n_nonzero=n_nonzero, acceptance=acceptance,
                  rel_xsec_err=rel_err)

    def trial_pred(bin_index, weights):
        return trial_covariance(bin_index, weights, nbins, n_nonzero,
                                rel_xsec_err=rel_err)

    def boot_pred(bin_index, weights):
        return bootstrap_covariance(bin_index, weights, nbins, n_boot=200, rng=rng)

    def coverage(make, n_runs):
        """Normalisation uncertainty the estimator reports / the runs' actual scatter."""
        preds = [make(*run) for run in _null_control_runs(rng, n_runs, **kwargs)]
        totals = np.array([p.values.sum() for p in preds])
        one = np.ones(nbins)
        claimed = np.mean([np.sqrt(one @ p.covariance @ one) / p.values.sum()
                           for p in preds])
        return claimed / (totals.std(ddof=1) / totals.mean())

    # --- calibration: independent pairs of runs, nothing different but the seed ---
    pairs = 400
    runs = list(_null_control_runs(rng, 2 * pairs, **kwargs))
    pvals = np.array([compatibility(trial_pred(*runs[2 * i]),
                                    trial_pred(*runs[2 * i + 1])).pvalue
                      for i in range(pairs)])
    rate = float(np.mean(pvals < alpha))
    ks_p = float(stats.kstest(pvals, "uniform").pvalue)
    # +-3 sigma of Binomial(400, alpha), i.e. the rate really does sit at alpha
    lo, hi = alpha - 3 * np.sqrt(alpha * (1 - alpha) / pairs), \
        alpha + 3 * np.sqrt(alpha * (1 - alpha) / pairs)
    ok_rate = lo <= rate <= hi
    print(f"[trial-null] flag rate={rate:.3f} (target {alpha}, allowed "
          f"{lo:.3f}-{hi:.3f}), mean p={pvals.mean():.3f}, uniformity KS p={ks_p:.3f}")

    cover_trial = coverage(trial_pred, 200)
    cover_boot = coverage(boot_pred, 100)
    print(f"[trial-cover] normalisation uncertainty / actual scatter: "
          f"trial counts {cover_trial:.2f} (target ~1), "
          f"resampled selected events {cover_boot:.2f} (the bug: << 1)")

    # cover_boot is ~0.6 here and ~0.3 on real MINERvA events (the synthetic weight
    # spread is wider than a real binning's, which flatters the resample); 0.75 keeps
    # a clear statement without sitting on top of the estimate's own noise.
    ok = (ok_rate and ks_p > 0.01 and 0.8 <= cover_trial <= 1.3
          and cover_boot < 0.75)
    return ok


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

    # --- Null control for the trial-count covariance (see _selftest_null_control).
    ok_trial = _selftest_null_control(np.random.default_rng(98765))

    ok = ok_null and ok_power and ok_bonf and ok_trial
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    print(__doc__)
