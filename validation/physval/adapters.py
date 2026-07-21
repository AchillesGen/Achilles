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


@dataclass
class GeneratedEvents:
    """Events for one experimental setup, reused across its measurements."""

    x: Optional[np.ndarray] = None        # synthetic: in-memory event positions
    weights: Optional[np.ndarray] = None  # synthetic: per-event weights
    path: Optional[str] = None            # nuisance3: NuHepMC-v1.0 event file
    target_nucleons: Optional[float] = None  # nucleons per target atom, for scaling

    def cleanup(self) -> None:
        """Drop the on-disk event file once every measurement has been histogrammed."""
        if self.path and os.path.exists(self.path):
            os.remove(self.path)


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
    edges: Optional[np.ndarray] = None  # bin edges, for plotting
    xlabel: str = "bin"
    ylabel: str = "d$\\sigma$/dx"


class Nuisance3Adapter:
    """Achilles generation + NUISANCE3 (legacy NUISANCE2 record) histogramming.

    generate:   achilles <experiment['achilles_run']> -> NuHepMC file (once/setup).
    histogram:  NUISANCE3 <measurement['name']> over that file -> (bin, weight).
    data_table: the sample's published values + covariance from NUISANCE3.

    ``achilles_run`` is repo-root-relative and achilles runs with cwd=repo root, so the
    flux/data paths inside a card resolve on their own (Filesystem::FindFile/FindFlux
    search the cwd, $ACHILLES_PATH, $ACHILLES_DATA_DIR and share/Achilles).

    Normalisation is delegated to NUISANCE: ``IAnalysis.process`` produces the
    cross-section-scaled, bin-width-divided prediction, and the per-event weights are
    rescaled so they sum to it. The bootstrap then measures MC uncertainty directly in
    the data's units without this code reimplementing any of the scaling.
    """

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
        if not os.path.isfile(card_path):
            # achilles segfaults rather than erroring on a missing card, so check here.
            raise FileNotFoundError(f"run card not found: {card_path}")
        with open(card_path) as fh:
            card = yaml.load(fh, Loader=_CardLoader)

        card.setdefault("Main", {})["NEvents"] = int(n_events)
        card["Main"]["Output"] = {"Format": "NuHepMC", "Name": out_path,
                                  "Zipped": False}

        # Seed lives in Options, which the cards pull in via `!include`; expand that
        # include so the seed can be pinned without editing the shared defaults.
        options = card.get("Options")
        if isinstance(options, _Include):
            with open(os.path.join(self.repo_root, options.value)) as fh:
                options = yaml.load(fh, Loader=_CardLoader)
        card["Options"] = options or {}
        card["Options"].setdefault("Initialize", {})["Seed"] = int(seed)

        rendered = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.card.yml")
        with open(rendered, "w") as fh:
            yaml.dump(card, fh, Dumper=_CardDumper, sort_keys=False)
        return rendered

    @staticmethod
    def _runtime_env(exe: str) -> dict:
        """Environment for achilles with its own libraries ahead of the image's.

        The physval image puts /opt/nuisance2/lib on LD_LIBRARY_PATH, which outranks
        the binary's RUNPATH. NUISANCE2 ships spdlog built against fmt v10 while
        Achilles bundles fmt v11, so leaving that ordering alone loads both and
        Achilles segfaults inside InitializeLogging with no output at all. Putting
        Achilles' own lib dir first keeps the pair consistent.
        """
        env = dict(os.environ)
        libdir = os.path.join(os.path.dirname(os.path.dirname(
            os.path.realpath(exe))), "lib")
        parts = [libdir, libdir + "64"]
        if env.get("LD_LIBRARY_PATH"):
            parts.append(env["LD_LIBRARY_PATH"])
        env["LD_LIBRARY_PATH"] = ":".join(parts)
        return env

    def generate(self, experiment: dict, branch: str, seed: int,
                 n_events: int) -> GeneratedEvents:
        os.makedirs(self.workdir, exist_ok=True)
        # Offset the seed by branch so an inline 'main' is not the identical stream.
        seed = int(seed) + (0 if branch == "main" else 1)
        out_path = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.hepmc")
        card = self._render_card(experiment, branch, seed, n_events, out_path)

        exe = shutil.which(self.achilles) or self.achilles
        log_path = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.achilles.log")
        with open(log_path, "w") as log:
            proc = subprocess.run([exe, card], cwd=self.repo_root,
                                  env=self._runtime_env(exe),
                                  stdout=log, stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0:
            with open(log_path) as fh:
                tail = "\n".join(fh.read().splitlines()[-25:])
            # A negative code is a signal (e.g. -11 = SIGSEGV), which usually dies
            # without flushing anything useful, so say so rather than show a blank.
            how = (f"signal {-proc.returncode}" if proc.returncode < 0
                   else f"exit {proc.returncode}")
            raise RuntimeError(
                f"achilles failed for {experiment['name']} ({how}); card={card} "
                f"log={log_path}\n{tail or '<no output captured>'}")
        if not os.path.exists(out_path):
            raise RuntimeError(
                f"achilles produced no event file at {out_path} for {experiment['name']}")
        if "target_nucleons" not in experiment:
            raise KeyError(
                f"experiment {experiment['name']!r} needs 'target_nucleons' (nucleons "
                "per target atom, e.g. 13 for CH) to normalise the cross section")
        return GeneratedEvents(path=out_path,
                               target_nucleons=float(experiment["target_nucleons"]))

    # -- analysis ------------------------------------------------------------

    # 1 pb in cm^2: the frame reports fatx/sumw in pb, the published data is in cm^2.
    _PB_TO_CM2 = 1e-36

    def histogram(self, generated: GeneratedEvents,
                  measurement: dict) -> EventSample:
        pn = self._pn()
        analysis = self._analysis(measurement["name"])

        evs = pn.EventSource(generated.path)
        if not evs:
            raise RuntimeError(f"could not open event file {generated.path}")

        # The legacy NUISANCE2 record's add_to_framegen does not register the sample's
        # columns, so add the selection and projection operators explicitly.
        selection = analysis.get_selection()
        projections = analysis.get_projections()
        fg = pn.EventFrameGen(evs)
        fg.add_int_column(selection.fname, selection.op)
        for projection in projections:
            fg.add_double_column(projection.fname, projection.op)
        frame = fg.all()

        cols = {n: i for i, n in enumerate(frame.column_names)}
        missing = [n for n in [selection.fname, "weight.cv"] +
                   [p.fname for p in projections] if n not in cols]
        if missing:
            raise RuntimeError(f"event frame is missing columns {missing}; "
                               f"got {list(cols)}")

        binned = analysis.get_data()[0]
        binning = binned.binning
        nbins = int(np.asarray(binned.values).reshape(-1).shape[0])

        table = np.asarray(frame.table, dtype=float)
        if table.size == 0:
            return EventSample(bin_index=np.empty(0, int),
                               weights=np.empty(0, float), nbins=nbins)

        proj_cols = [cols[p.fname] for p in projections]
        selected = table[table[:, cols[selection.fname]] != 0]

        bin_index, weights = [], []
        for row in selected:
            values = [float(row[c]) for c in proj_cols]
            b = (binning.find_bin(values[0]) if len(values) == 1
                 else binning.find_bin(pn.Vector_double(values)))
            if b is None or b < 0 or b >= nbins:
                continue  # under/overflow: outside the published binning
            bin_index.append(int(b))
            weights.append(float(row[cols["weight.cv"]]))

        bin_index = np.asarray(bin_index, dtype=int)
        weights = np.asarray(weights, dtype=float)

        # Cross-section normalisation, per the recipe that reproduces the published
        # results: take the flux-averaged total xsec *per atom* and divide by the
        # target's nucleon count ourselves.
        #
        # NUISANCE's own PerNucleon conversion must NOT be used here: for a composite
        # target it divides by the struck nucleus' A, which for CH measures out at
        # exactly 12 (carbon only), whereas the published data uses 13 (CH). That is a
        # silent 13/12 = 8% normalisation error.
        if generated.target_nucleons is None:
            raise ValueError("GeneratedEvents.target_nucleons is required to "
                             "normalise; set 'target_nucleons' on the experiment")
        fatx_col = cols.get("fatx_per_sumw.pb_per_target.estimate")
        if fatx_col is None:
            raise RuntimeError("event frame carries no per-target fatx estimate; "
                               f"got {list(cols)}")
        # Use the frame's own per-TARGET estimate: it is normalised against the same
        # weight.cv column used above. (EventSource.norm_info reports a differently
        # normalised sumweights -- 0.45 where the frame's weights sum to ~92000 --
        # so mixing the two overstates the prediction by ~4 orders of magnitude.)
        fatx_per_sumw = float(table[-1, fatx_col])

        widths = np.asarray(list(binning.bin_sizes()), dtype=float).reshape(-1)[:nbins]
        scale = np.divide(
            fatx_per_sumw * self._PB_TO_CM2 / generated.target_nucleons, widths,
            out=np.zeros(nbins), where=widths != 0)
        if bin_index.size:
            weights = weights * scale[bin_index]

        return EventSample(bin_index=bin_index, weights=weights, nbins=nbins)

    def data_table(self, measurement: dict) -> DataTable:
        pn = self._pn()
        analysis = self._analysis(measurement["name"])
        binned = analysis.get_data()[0]

        values = np.asarray(binned.values, dtype=float).reshape(-1)
        errors = np.asarray(binned.errors, dtype=float).reshape(-1)
        covariance = np.asarray(analysis.get_covariance_matrix(), dtype=float)

        if covariance.shape != (values.size, values.size):
            # No published covariance: fall back to the per-bin errors.
            covariance = np.diag(errors ** 2)
        else:
            # NUISANCE returns the covariance in the sample's published units (e.g.
            # 1e-38 cm^2 squared) while values/errors are absolute, so the two are not
            # directly comparable -- for these samples the diagonal is off by 1e76.
            # Recover the factor from the errors, which are in the values' units, so
            # the correlation structure is kept and no unit convention is assumed.
            sd = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
            usable = (sd > 0) & (errors > 0)
            if usable.any():
                covariance = covariance * float(np.median(errors[usable] / sd[usable])) ** 2

        try:
            edges = np.asarray(pn.Binning.get_bin_edges1D(binned.binning.bins),
                               dtype=float)
        except Exception:
            edges = None  # multi-dimensional binning: plot against bin number

        projections = analysis.get_projections()
        xlabel = "bin"
        if projections:
            p = projections[0]
            xlabel = p.prettyname or p.fname
            if p.units:
                xlabel = f"{xlabel} [{p.units}]"

        return DataTable(values=values, covariance=covariance, edges=edges,
                         xlabel=xlabel,
                         ylabel="d$\\sigma$/d" + (projections[0].fname
                                                  if projections else "x"))


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
        return DataTable(values=values, covariance=np.diag(err ** 2),
                         edges=np.linspace(0.0, 1.0, nbins + 1), xlabel="x")
