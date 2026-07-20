#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles + NUISANCE3 boundary for physval.

Runs inside the achilles-physval container, so ``achilles`` and NUISANCE3 are on
PATH. Events are generated once per experimental setup (``generate``) and reused
across that setup's measurements (``histogram``).

``Nuisance3Adapter`` is the real path (skeleton). ``SyntheticAdapter`` backs
``--dry-run`` and the self-tests; both expose the same three methods, so the driver
just picks one — no base class needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class GeneratedEvents:
    """Events for one experimental setup, reused across its measurements."""

    x: Optional[np.ndarray] = None        # synthetic: in-memory event positions
    weights: Optional[np.ndarray] = None  # synthetic: per-event weights
    path: Optional[str] = None            # nuisance3: NuHepMC-v1.0 event file


@dataclass
class EventSample:
    """One measurement's events, binned and ready for bootstrapping."""

    bin_index: np.ndarray
    weights: np.ndarray
    nbins: int


@dataclass
class DataTable:
    """A published measurement: central values and their covariance."""

    values: np.ndarray
    covariance: np.ndarray


class Nuisance3Adapter:
    """Achilles generation + NUISANCE3 histogramming. NOT YET WIRED UP.

    generate:   achilles <experiment['achilles_run']> -> NuHepMC-v1.0 file (once/setup).
    histogram:  NUISANCE3 <measurement['name']> over that file -> (bin, weight).
    data_table: the sample's published data values + covariance from NUISANCE3.

    ``achilles_run`` is repo-root-relative, so run achilles with cwd=repo root (or put
    it on $ACHILLES_PATH). The flux/data paths inside a card resolve on their own:
    Filesystem::FindFile/FindFlux search cwd, $ACHILLES_PATH, $ACHILLES_DATA_DIR and
    the install tree's share/Achilles, which is where CMake ships flux/ and data/.
    """

    def __init__(self, workdir: str = "physval-work"):
        self.workdir = workdir

    def generate(self, experiment: dict, branch: str, seed: int,
                 n_events: int) -> GeneratedEvents:  # pragma: no cover
        raise NotImplementedError(
            "Nuisance3Adapter.generate: run Achilles for "
            f"'{experiment.get('achilles_run', experiment.get('name'))}'.")

    def histogram(self, generated: GeneratedEvents,
                  measurement: dict) -> EventSample:  # pragma: no cover
        raise NotImplementedError(
            "Nuisance3Adapter.histogram: run NUISANCE3 for "
            f"'{measurement['name']}'.")

    def data_table(self, measurement: dict) -> DataTable:  # pragma: no cover
        raise NotImplementedError(
            "Nuisance3Adapter.data_table: pull NUISANCE3 data for "
            f"'{measurement['name']}'.")


class SyntheticAdapter:
    """Dry-run/self-test stand-in: draws weighted events from a tunable Gaussian.

    ``experiment['dryrun']['feature_shift']`` shifts only the feature branch (a fake
    physics change); ``measurement['dryrun']['nbins']`` sets the histogram binning.
    """

    def __init__(self, base_seed: int = 0):
        self.base_seed = base_seed

    def _nbins(self, measurement: dict) -> int:
        return int(measurement.get("dryrun", {}).get("nbins", 12))

    def _draw(self, shift: float, n_events: int,
              rng: np.random.Generator) -> GeneratedEvents:
        x = np.clip(rng.normal(0.5 + shift, 0.18, size=n_events), 0.0, 0.999)
        weights = rng.uniform(0.5, 1.5, size=n_events)
        return GeneratedEvents(x=x, weights=weights)

    def generate(self, experiment: dict, branch: str, seed: int,
                 n_events: int) -> GeneratedEvents:
        rng = np.random.default_rng(
            self.base_seed + seed + (0 if branch == "main" else 1))
        shift = float(experiment.get("dryrun", {}).get("feature_shift", 0.0)) \
            if branch == "feature" else 0.0
        return self._draw(shift, n_events, rng)

    def histogram(self, generated: GeneratedEvents,
                  measurement: dict) -> EventSample:
        nbins = self._nbins(measurement)
        bin_index = np.digitize(generated.x, np.linspace(0.0, 1.0, nbins + 1)) - 1
        return EventSample(bin_index=bin_index, weights=generated.weights,
                           nbins=nbins)

    def data_table(self, measurement: dict) -> DataTable:
        nbins = self._nbins(measurement)
        rng = np.random.default_rng(self.base_seed + 999)
        gen = self._draw(0.0, 200000, rng)  # data = the unshifted truth
        bin_index = np.digitize(gen.x, np.linspace(0.0, 1.0, nbins + 1)) - 1
        values = np.bincount(bin_index, weights=gen.weights, minlength=nbins)[:nbins]
        err = np.sqrt(np.maximum(values, 1.0)) * 0.05 + 0.02 * values
        return DataTable(values=values, covariance=np.diag(err ** 2))
