#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Adapter boundary between the physval statistics and Achilles / NUISANCE3.

Everything that actually shells out to Achilles (event generation) or NUISANCE3
(histogramming an analysis and pulling its data + covariance) lives behind the
``Adapter`` interface so the rest of the pipeline stays testable.

Two implementations:

* ``SyntheticAdapter`` - no external tools; draws weighted events from a tunable
  distribution.  Used by ``--dry-run`` and the self-tests.
* ``Nuisance3Adapter`` - the real path.  It is intentionally a documented skeleton:
  the exact NUISANCE3 CLI, the container image, and the per-measurement sample
  names are the external hand-offs called out in the plan and must be filled in
  once the NUISANCE3 image and the Achilles NuHepMC-v1.0 output are wired up.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass
class EventSample:
    """Weighted events assigned to analysis bins, ready for bootstrapping."""

    bin_index: np.ndarray   # shape (n_events,), int bin assignment
    weights: np.ndarray     # shape (n_events,), float event weights
    nbins: int


@dataclass
class DataTable:
    """The published measurement: central data values and their covariance."""

    values: np.ndarray      # shape (nbins,)
    covariance: np.ndarray  # shape (nbins, nbins)


class Adapter(abc.ABC):
    """Produces, for a measurement, the feature/main event samples and the data."""

    @abc.abstractmethod
    def data_table(self, measurement: dict) -> DataTable:
        ...

    @abc.abstractmethod
    def generate(self, measurement: dict, branch: str, seed: int,
                 n_events: int) -> EventSample:
        """Generate ``n_events`` for ``measurement`` on ``branch`` ('main'|'feature')."""
        ...


# ---------------------------------------------------------------------------
# Synthetic adapter (no external tools) - used for --dry-run and tests
# ---------------------------------------------------------------------------

class SyntheticAdapter(Adapter):
    """A stand-in that fakes generation with a Gaussian-ish distribution.

    ``measurement`` may carry a ``dryrun`` block::

        dryrun:
          nbins: 12
          feature_shift: 0.0   # bin-space shift applied only to the feature branch

    A non-zero ``feature_shift`` simulates a physics change so the injected-regression
    path can be exercised end to end.
    """

    def __init__(self, base_seed: int = 0):
        self.base_seed = base_seed

    def _nbins(self, measurement: dict) -> int:
        return int(measurement.get("dryrun", {}).get("nbins", 12))

    def _edges(self, nbins: int) -> np.ndarray:
        return np.linspace(0.0, 1.0, nbins + 1)

    def data_table(self, measurement: dict) -> DataTable:
        nbins = self._nbins(measurement)
        rng = np.random.default_rng(self.base_seed + 999)
        # "Data" = the unshifted truth with a modest diagonal uncertainty.
        sample = self._draw(measurement, shift=0.0, n_events=200000, rng=rng)
        values = np.bincount(sample.bin_index, weights=sample.weights,
                             minlength=nbins)[:nbins]
        err = np.sqrt(np.maximum(values, 1.0)) * 0.05 + 0.02 * values
        cov = np.diag(err ** 2)
        return DataTable(values=values, covariance=cov)

    def _draw(self, measurement: dict, shift: float, n_events: int,
              rng: np.random.Generator) -> EventSample:
        nbins = self._nbins(measurement)
        edges = self._edges(nbins)
        x = np.clip(rng.normal(0.5 + shift, 0.18, size=n_events), 0.0, 0.999)
        bin_index = np.digitize(x, edges) - 1
        weights = rng.uniform(0.5, 1.5, size=n_events)
        return EventSample(bin_index=bin_index, weights=weights, nbins=nbins)

    def generate(self, measurement: dict, branch: str, seed: int,
                 n_events: int) -> EventSample:
        rng = np.random.default_rng(self.base_seed + seed + (0 if branch == "main" else 1))
        shift = 0.0
        if branch == "feature":
            shift = float(measurement.get("dryrun", {}).get("feature_shift", 0.0))
        return self._draw(measurement, shift=shift, n_events=n_events, rng=rng)


# ---------------------------------------------------------------------------
# Real NUISANCE3 adapter - documented skeleton (external hand-off)
# ---------------------------------------------------------------------------

class Nuisance3Adapter(Adapter):
    """Real generation + NUISANCE3 analysis.  NOT YET WIRED UP.

    Intended flow per measurement (see validation/physval/README.md):

      1. Render an Achilles run card for ``measurement['achilles_run']`` with the
         pinned ``seed`` and ``n_events``; run ``achilles <card>`` to produce a
         NuHepMC-v1.0 event file (stream / gzip; delete after step 2).
      2. Run the NUISANCE3 comparison app (with the NUISANCE legacy interface
         enabled) for ``measurement['nuisance_sample']`` over that event file to
         obtain the per-bin prediction, and pull the sample's published data values
         and covariance.
      3. Return the per-event ``(bin_index, weight)`` arrays so the caller can form
         the bootstrap covariance, plus the ``DataTable``.

    TODOs (the external hand-offs from the plan):
      * NUISANCE3 container image + version pin.
      * Exact NUISANCE3 CLI / Python API for "run sample, dump per-event bin+weight".
      * Mapping measurement -> Achilles run card + NUISANCE sample name.
    """

    def __init__(self, workdir: str, nuisance_image: str):
        self.workdir = workdir
        self.nuisance_image = nuisance_image

    def data_table(self, measurement: dict) -> DataTable:  # pragma: no cover
        raise NotImplementedError(
            "Nuisance3Adapter.data_table: wire up NUISANCE3 sample data extraction "
            f"for '{measurement.get('nuisance_sample', measurement.get('name'))}'.")

    def generate(self, measurement: dict, branch: str, seed: int,
                 n_events: int) -> EventSample:  # pragma: no cover
        raise NotImplementedError(
            "Nuisance3Adapter.generate: wire up Achilles generation + NUISANCE3 "
            "histogramming (container image, sample names). See README.md.")


def make_adapter(kind: str, **kwargs) -> Adapter:
    if kind == "synthetic":
        return SyntheticAdapter(base_seed=int(kwargs.get("base_seed", 0)))
    if kind == "nuisance3":
        return Nuisance3Adapter(workdir=kwargs["workdir"],
                                nuisance_image=kwargs["nuisance_image"])
    raise ValueError(f"unknown adapter kind: {kind!r}")
