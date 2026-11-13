#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2018-2026 Achilles Developers
# SPDX-License-Identifier: GPL-3.0-or-later
"""Achilles + NUISANCE3 boundary for physval.

``Nuisance3Adapter`` runs achilles and bins with NUISANCE3 (inside the physval
container); ``SyntheticAdapter`` backs ``--dry-run`` and the tests. Both expose
``generate``, ``histogram_many``, ``response_matrix`` and ``data_table``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional

import numpy as np
import yaml


# The run cards use `!include <path>` for shared blocks. Round-trip the tag so a
# rendered card keeps its includes instead of tripping up the YAML loader.
class _Include:
    def __init__(self, value: str):
        self.value = value


class _CardLoader(yaml.SafeLoader):
    pass


class _CardDumper(yaml.SafeDumper):
    pass


_CardLoader.add_constructor(
    "!include", lambda loader, node: _Include(loader.construct_scalar(node)))
_CardDumper.add_representer(
    _Include, lambda dumper, data: dumper.represent_scalar("!include", data.value))


def read_gen_cross_section(path: str, tail: int = 1 << 24) -> "tuple[float, float, int]":
    """``(xsec, xsec_err, n_nonzero_trials)`` from the file's end; walking it costs minutes."""
    with open(path, "rb") as fh:
        fh.seek(max(0, os.path.getsize(path) - tail))
        lines = [l for l in fh.read().decode("utf-8", errors="replace").splitlines()
                 if l.startswith("A 0 GenCrossSection")]
    if not lines:
        raise RuntimeError(
            f"no GenCrossSection attribute at the end of {path}; the MC uncertainty cannot "
            "be built from trial counts. A compressed event file looks like this -- "
            "generate with Main.Output.Zipped: false.")
    xsec, err, nonzero = lines[-1].split()[3:6]
    return float(xsec), float(err), int(nonzero)


@dataclass
class GeneratedEvents:
    """Events for one experimental setup, reused across its measurements."""

    x: Optional[np.ndarray] = None        # synthetic: in-memory event positions
    weights: Optional[np.ndarray] = None  # synthetic: per-event weights
    path: Optional[str] = None            # nuisance3: NuHepMC-v1.0 event file
    xsec_divisor: Optional[float] = None  # per-atom xsec / this = published convention
    # What stats.trial_covariance needs.
    n_nonzero_trials: Optional[int] = None
    rel_xsec_err: Optional[float] = None

    def cleanup(self) -> None:
        """Drop the on-disk event file once every measurement has been histogrammed."""
        if self.path and os.path.exists(self.path):
            os.remove(self.path)


@dataclass
class EventSample:
    """One measurement's selected events: bin index and scaled weight per event."""

    bin_index: np.ndarray
    weights: np.ndarray
    nbins: int


@dataclass
class DataTable:
    """A published measurement: central values and their covariance."""

    values: np.ndarray
    covariance: np.ndarray
    edges: Optional[np.ndarray] = None  # bin edges, for plotting


class Nuisance3Adapter:
    """Achilles generation + NUISANCE3 (legacy NUISANCE2 record) binning.

    Cards are repo-root-relative and achilles runs from the repo root, so the flux and
    data paths inside a card resolve on their own.
    """

    name = "nuisance3"

    def __init__(self, workdir: str = "physval-work", repo_root: Optional[str] = None,
                 achilles: str = "achilles"):
        self.workdir = os.path.abspath(workdir)
        # Default: this file lives at <repo>/validation/physval/adapters.py.
        self.repo_root = os.path.abspath(
            repo_root or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      os.pardir, os.pardir))
        self.achilles = achilles
        self._record = None

    # -- NUISANCE handles ----------------------------------------------------

    def _pn(self):
        import pyNUISANCE as pn  # imported lazily: only present inside the container
        return pn

    def _analysis(self, name: str):
        pn = self._pn()
        if self._record is None:
            self._record = pn.RecordFactory().make_record({"type": "nuisance2"})
        return self._record.analysis(name)

    # -- generation ----------------------------------------------------------

    def _render_card(self, experiment: dict, branch: str, seed: int, n_events: int,
                     out_path: str) -> str:
        """Copy the run card with our seed, event count and output path pinned."""
        os.makedirs(self.workdir, exist_ok=True)
        card_path = os.path.join(self.repo_root, experiment["achilles_run"])
        if not os.path.isfile(card_path):  # achilles segfaults on a missing card
            raise FileNotFoundError(f"run card not found: {card_path}")
        with open(card_path) as fh:
            card = yaml.load(fh, Loader=_CardLoader)

        card.setdefault("Main", {})["NEvents"] = int(n_events)
        card["Main"]["Output"] = {"Format": "NuHepMC", "Name": out_path,
                                  "Zipped": False}

        # The seed lives in the `!include`d Options block; expand it to pin the seed.
        options = card.get("Options")
        if isinstance(options, _Include):
            with open(os.path.join(self.repo_root, options.value)) as fh:
                options = yaml.load(fh, Loader=_CardLoader)
        card["Options"] = options or {}
        card["Options"].setdefault("Initialize", {})["Seed"] = int(seed)
        # A loaded integrator cache moves the RNG stream, so the same seed would give
        # different events; keep it off whatever the card says.
        card["Cache"] = {"Save": False, "Load": False}

        rendered = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.card.yml")
        with open(rendered, "w") as fh:
            yaml.dump(card, fh, Dumper=_CardDumper, sort_keys=False)
        return rendered


    def generate(self, experiment: dict, branch: str, seed: int,
                 n_events: int) -> GeneratedEvents:
        os.makedirs(self.workdir, exist_ok=True)
        # main and PR draw independent streams.
        seed = int(seed) + (0 if branch == "main" else 1)
        out_path = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.hepmc")
        card = self._render_card(experiment, branch, seed, n_events, out_path)

        exe = shutil.which(self.achilles) or self.achilles
        log_path = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.achilles.log")
        with open(log_path, "w") as log:
            proc = subprocess.run([exe, card], cwd=self.repo_root,
                                  stdout=log, stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0:
            with open(log_path) as fh:
                tail = "\n".join(fh.read().splitlines()[-25:])
            how = (f"signal {-proc.returncode}" if proc.returncode < 0
                   else f"exit {proc.returncode}")
            raise RuntimeError(
                f"achilles failed for {experiment['name']} ({how}); card={card} "
                f"log={log_path}\n{tail or '<no output captured>'}")
        if not os.path.exists(out_path):
            raise RuntimeError(
                f"achilles produced no event file at {out_path} for {experiment['name']}")
        xsec, xsec_err, n_nonzero = read_gen_cross_section(out_path)
        return GeneratedEvents(path=out_path,
                               xsec_divisor=float(experiment["xsec_divisor"]),
                               n_nonzero_trials=n_nonzero,
                               rel_xsec_err=(xsec_err / xsec) if xsec else None)


    # -- analysis ------------------------------------------------------------

    # 1 pb in cm^2: the frame reports fatx/sumw in pb, the published data is in cm^2.
    _PB_TO_CM2 = 1e-36

    def response_matrix(self, measurement: dict, nbins: int) -> np.ndarray:
        """Binned rates to the published quantity, as R S and R C R^T.

        Call after binning, since a sample's conversion can depend on the events it saw.
        """
        analysis = self._analysis(measurement["name"])
        solid_angle = measurement.get("solid_angle")
        if solid_angle:  # the Durham samples' own conversion has no normalisation
            widths = np.asarray(list(analysis.get_data()[0].binning.bin_sizes()),
                                dtype=float).reshape(-1)[:nbins]
            return np.diag(np.divide(1.0 / float(solid_angle), widths,
                                     out=np.zeros(nbins), where=widths != 0))
        return np.column_stack([analysis.convert_event_rates(unit, 1.0)
                                for unit in np.eye(nbins)])

    def histogram_many(self, generated: GeneratedEvents,
                       measurements) -> "dict[str, EventSample]":
        """Bin every measurement of a setup in one blocked pass over the event file."""
        pn = self._pn()
        measurements = list(measurements)
        if not measurements:
            return {}

        evs = pn.EventSource(generated.path)
        if not evs:
            raise RuntimeError(f"could not open event file {generated.path}")

        fg = pn.EventFrameGen(evs, 250_000)
        specs = []
        for i, measurement in enumerate(measurements):
            analysis = self._analysis(measurement["name"])
            # The legacy record does not register its columns; add them under our own
            # names, since a sample's fname can collide with another's.
            projections = analysis.get_projections()
            sel_col = f"sel{i}"
            proj_cols = [f"proj{i}_{j}" for j in range(len(projections))]
            fg.add_int_column(sel_col, analysis.get_selection().op)
            for name, projection in zip(proj_cols, projections):
                fg.add_double_column(name, projection.op)

            binned = analysis.get_data()[0]
            specs.append({
                "measurement": measurement,
                "binning": binned.binning,
                "nbins": int(np.asarray(binned.values).reshape(-1).shape[0]),
                "sel": sel_col,
                "projs": proj_cols,
                "bins": [],
                "weights": [],
            })

        fatx_per_sumw = None
        block = fg.first()
        while block is not None:
            table = np.asarray(block.table, dtype=float)
            if table.size == 0:
                break
            cols = {n: i for i, n in enumerate(block.column_names)}
            # A running estimate: the newest row wins.
            fatx_per_sumw = float(table[-1, cols["fatx_per_sumw.pb_per_target.estimate"]])
            weights = table[:, cols["weight.cv"]]
            for spec in specs:
                bins = np.asarray(spec["binning"].find_bins(block, spec["projs"]),
                                  dtype=np.int64).reshape(-1)
                # Under/overflow is outside the published binning.
                keep = (table[:, cols[spec["sel"]]] != 0) & (bins >= 0) & (bins < spec["nbins"])
                spec["bins"].append(bins[keep])
                spec["weights"].append(weights[keep])
            block = fg.next()

        # Per-atom fatx divided by our own nucleon count: NUISANCE's PerNucleon uses
        # the struck nucleus' A (12 for CH, not 13). The frame's own estimate matches
        # weight.cv; EventSource.norm_info does not.
        scale = (fatx_per_sumw or 0.0) * self._PB_TO_CM2 / generated.xsec_divisor
        return {spec["measurement"]["name"]: EventSample(
                    bin_index=np.concatenate(spec["bins"] or [np.empty(0, dtype=np.int64)]),
                    weights=np.concatenate(spec["weights"] or [np.empty(0)]) * scale,
                    nbins=spec["nbins"])
                for spec in specs}

    def data_table(self, measurement: dict) -> DataTable:
        pn = self._pn()
        analysis = self._analysis(measurement["name"])
        binned = analysis.get_data()[0]

        values = np.asarray(binned.values, dtype=float).reshape(-1)
        errors = np.asarray(binned.errors, dtype=float).reshape(-1)
        covariance = np.asarray(analysis.get_covariance_matrix(), dtype=float)

        data_scale = float(measurement.get("data_scale", 1.0))
        if data_scale != 1.0:
            values = values * data_scale
            errors = errors * data_scale
            covariance = covariance * data_scale ** 2

        if covariance.shape != (values.size, values.size):
            covariance = np.diag(errors ** 2)
        else:
            # The covariance comes in the sample's published units (e.g. 1e-38 cm^2
            # squared); rescale it onto the errors, which match the values.
            sd = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
            usable = (sd > 0) & (errors > 0)
            if usable.any():
                covariance = covariance * float(np.median(errors[usable] / sd[usable])) ** 2

        try:
            edges = np.asarray(pn.Binning.get_bin_edges1D(binned.binning.bins),
                               dtype=float)
        except Exception:
            edges = None  # multi-dimensional binning: plot against bin number

        return DataTable(values=values, covariance=covariance, edges=edges)


class SyntheticAdapter:
    """Weighted Gaussian events. A setup's ``feature_shift`` fakes a physics change in
    the PR."""

    name = "synthetic"
    _NBINS = 12

    def __init__(self, base_seed: int = 0):
        self.base_seed = base_seed

    # Every prediction is normalised to this total, as the real cross section would be.
    _TOTAL = 2.0e5

    def _draw(self, shift: float, n_events: int,
              rng: np.random.Generator) -> GeneratedEvents:
        x = np.clip(rng.normal(0.5 + shift, 0.18, size=n_events), 0.0, 0.999)
        return GeneratedEvents(x=x, weights=rng.uniform(0.5, 1.5, size=n_events))

    def generate(self, experiment: dict, branch: str, seed: int,
                 n_events: int) -> GeneratedEvents:
        rng = np.random.default_rng(self.base_seed + seed + (0 if branch == "main" else 1))
        shift = experiment.get("feature_shift", 0.0) if branch == "feature" else 0.0
        gen = self._draw(shift, n_events, rng)
        gen.n_nonzero_trials = int(gen.weights.size)
        gen.rel_xsec_err = 0.0
        return gen


    def histogram_many(self, generated: GeneratedEvents,
                       measurements) -> "dict[str, EventSample]":
        total = float(np.sum(np.abs(generated.weights)))
        scale = self._TOTAL / total if total > 0.0 else 1.0
        bin_index = np.digitize(generated.x, np.linspace(0.0, 1.0, self._NBINS + 1)) - 1
        return {m["name"]: EventSample(bin_index, generated.weights * scale, self._NBINS)
                for m in measurements}

    @staticmethod
    def response_matrix(measurement: dict, nbins: int) -> Optional[np.ndarray]:
        return None  # the synthetic weights are already in published units

    def data_table(self, measurement: dict) -> DataTable:
        nbins = self._NBINS
        rng = np.random.default_rng(self.base_seed + 999)
        gen = self._draw(0.0, 200000, rng)  # data = the unshifted truth
        bin_index = np.digitize(gen.x, np.linspace(0.0, 1.0, nbins + 1)) - 1
        values = np.bincount(bin_index, weights=gen.weights, minlength=nbins)[:nbins]
        values = values * (self._TOTAL / float(np.sum(values)))
        err = np.sqrt(np.maximum(values, 1.0)) * 0.05 + 0.02 * values
        return DataTable(values=values, covariance=np.diag(err ** 2),
                         edges=np.linspace(0.0, 1.0, nbins + 1))
