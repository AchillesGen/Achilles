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
import re
import shutil
import subprocess
import types
from dataclasses import dataclass
from typing import Optional

import numpy as np
import yaml


def _as_float(text: str) -> Optional[float]:
    """Parse a scraped number, discarding the -nan/0 that empty groups print."""
    try:
        value = float(text)
    except ValueError:
        return None
    return value if np.isfinite(value) and value > 0.0 else None


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
    xsec_divisor: Optional[float] = None  # per-atom xsec / this = published convention
    run: Optional["RunStats"] = None      # generator-side cost of this sample

    def cleanup(self) -> None:
        """Drop the on-disk event file once every measurement has been histogrammed."""
        if self.path and os.path.exists(self.path):
            os.remove(self.path)


@dataclass
class RunStats:
    """What the generation itself cost, independent of any one measurement.

    ``unweight_eff`` is Achilles' own estimate (<xsec>/max_weight per process group,
    the acceptance an ideal unweighter would reach for the cap it settled on); the
    smallest finite group value is kept, since the slowest group paces the run.
    """

    seconds: Optional[float] = None
    unweight_eff: Optional[float] = None


@dataclass
class EventSample:
    """One measurement's events, binned and ready for bootstrapping."""

    bin_index: np.ndarray
    weights: np.ndarray
    nbins: int
    raw_weights: Optional[np.ndarray] = None  # pre-normalisation weight.cv

    def ess_fraction(self) -> float:
        """Kish effective sample size as a fraction of the selected events.

        ``(Σ|w|)² / (N Σw²)`` — 1 for perfectly unit-weight events, and the factor by
        which a scheme's statistical power falls short of its raw event count. Scale
        free, so it is computed on the raw generator weights rather than the
        bin-width-divided ones (whose spread is binning, not unweighting).
        """
        w = self.raw_weights if self.raw_weights is not None else self.weights
        w = np.abs(np.asarray(w, dtype=float))
        denom = float(w.size) * float(np.sum(w ** 2))
        if denom <= 0.0:
            return float("nan")
        return float(np.sum(w) ** 2 / denom)

    def max_over_mean(self) -> float:
        """Heaviest selected weight in units of the mean — the overweight tail."""
        w = np.abs(np.asarray(
            self.raw_weights if self.raw_weights is not None else self.weights,
            dtype=float))
        mean = float(np.mean(w)) if w.size else 0.0
        return float(np.max(w) / mean) if mean > 0.0 else float("nan")


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
                     out_path: str, unweighting: Optional[dict] = None) -> str:
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
        # The unweighting scan swaps this whole block per variant. Replace rather than
        # merge: the schemes take different keys (percentile vs epsilon), so a leftover
        # key from the card's default would silently apply to the wrong scheme.
        if unweighting is not None:
            card["Options"]["Unweighting"] = dict(unweighting)

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

    def generate(self, experiment: dict, branch: str, seed: int, n_events: int, *,
                 unweighting: Optional[dict] = None,
                 seed_offset: Optional[int] = None) -> GeneratedEvents:
        os.makedirs(self.workdir, exist_ok=True)
        # Offset the seed by branch so an inline 'main' is not the identical stream.
        # The unweighting scan pins seed_offset=0 instead: every variant then starts
        # from the same stream, so what differs between them is only the scheme.
        seed = int(seed) + (seed_offset if seed_offset is not None
                            else (0 if branch == "main" else 1))
        out_path = os.path.join(self.workdir,
                                f"{experiment['name']}_{branch}.hepmc")
        card = self._render_card(experiment, branch, seed, n_events, out_path,
                                 unweighting=unweighting)

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
        return GeneratedEvents(path=out_path,
                               xsec_divisor=self._xsec_divisor(experiment),
                               run=self._run_stats(log_path))

    # Achilles' own end-of-run numbers. Both are printed rather than written to a
    # machine-readable file, so they are scraped; a miss just leaves the field blank.
    _RE_DURATION = re.compile(r"Run Duration:\s*(?:(\d+)h\s*)?(?:(\d+)m\s*)?(\d+)s")
    _RE_EFF = re.compile(r"Estimated unweighting eff for this group:\s*(\S+)")

    @classmethod
    def _run_stats(cls, log_path: str) -> RunStats:
        try:
            with open(log_path, errors="replace") as fh:
                text = fh.read()
        except OSError:
            return RunStats()

        seconds = None
        m = cls._RE_DURATION.search(text)
        if m:
            h, mi, s = (int(g or 0) for g in m.groups())
            seconds = float(3600 * h + 60 * mi + s)

        # Process groups with no allowed states print -nan; keep the least efficient
        # real group, which is the one that paces the run.
        effs = [v for v in (_as_float(x) for x in cls._RE_EFF.findall(text)) if v]
        return RunStats(seconds=seconds,
                        unweight_eff=min(effs) if effs else None)

    @staticmethod
    def _xsec_divisor(experiment: dict) -> float:
        """How much to divide the per-ATOM cross section by for this experiment.

        Achilles reports per atom, but published data does not use one convention:
        T2K/MINERvA quote CH per nucleon (divide by 13), while MicroBooNE quotes per
        argon atom (divide by 1 -- *not* by 40). Getting this wrong is silent and
        costs a factor of the target's mass number, so it is declared per experiment
        rather than guessed.
        """
        per = experiment.get("data_per")
        if per not in ("nucleon", "atom"):
            raise KeyError(
                f"experiment {experiment['name']!r} needs data_per: nucleon|atom "
                f"(the published cross-section convention), got {per!r}")
        if per == "atom":
            return 1.0
        if "target_nucleons" not in experiment:
            raise KeyError(
                f"experiment {experiment['name']!r} has data_per: nucleon and so "
                "needs 'target_nucleons' (e.g. 13 for CH)")
        return float(experiment["target_nucleons"])

    # -- analysis ------------------------------------------------------------

    # 1 pb in cm^2: the frame reports fatx/sumw in pb, the published data is in cm^2.
    _PB_TO_CM2 = 1e-36

    @staticmethod
    def _extra_bin_scale(measurement: dict, nbins: int) -> np.ndarray:
        """Per-bin factors the NUISANCE binning does not already account for.

        Dividing by the widths of the NUISANCE binning is all a sample with a physical
        axis needs. Two cases are not: a histogram indexed by bin *number*, whose real
        widths come from ``bin_edges`` (the sample divides by them in its own
        ConvertEventRates, which the legacy record never calls), and a cross section
        published per steradian, whose ``solid_angle`` the selection integrates over.
        """
        scale = np.ones(nbins)
        edges = measurement.get("bin_edges")
        if edges is not None:
            widths = np.diff(np.asarray(edges, dtype=float))
            if widths.size != nbins:
                raise ValueError(
                    f"{measurement['name']}: bin_edges gives {widths.size} bins, "
                    f"NUISANCE reports {nbins}")
            scale = scale / widths
        omega = measurement.get("solid_angle")
        if omega:
            scale = scale / float(omega)
        return scale

    @staticmethod
    def response_matrix(measurement: dict, nbins: int) -> Optional[np.ndarray]:
        """The regularisation matrix A_C for ``measurement``, or None.

        A Wiener-SVD unfolded measurement is only comparable to A_C * prediction, so
        the matrix is applied to the prediction (in stats.bootstrap_covariance, which
        carries it into the MC covariance too), never to the data.
        """
        path = measurement.get("smearing")
        if not path:
            return None
        matrix = np.loadtxt(path, delimiter=",", ndmin=2)
        if matrix.shape != (nbins, nbins):
            raise ValueError(f"{measurement['name']}: smearing matrix is "
                             f"{matrix.shape}, NUISANCE reports {nbins} bins")
        return matrix

    def histogram(self, generated: GeneratedEvents,
                  measurement: dict) -> EventSample:
        """One measurement; a thin wrapper over the batched pass."""
        return self.histogram_many(generated, [measurement])[measurement["name"]]

    def histogram_many(self, generated: GeneratedEvents, measurements,
                       block_size: int = 250_000) -> "dict[str, EventSample]":
        """Bin every measurement of a setup in a single pass over the event file.

        The legacy record evaluates one sample's selection and projections per column,
        so a frame can carry all of a setup's samples at once (this is the notebook's
        "lots of projections" pattern): the file is opened, parsed and walked once
        instead of once per measurement, which is where the time goes for a setup with
        eighteen of them. Events are pulled in blocks rather than with ``all()`` so
        peak memory stays independent of the run length.
        """
        pn = self._pn()
        measurements = list(measurements)
        if not measurements:
            return {}

        evs = pn.EventSource(generated.path)
        if not evs:
            raise RuntimeError(f"could not open event file {generated.path}")

        fg = pn.EventFrameGen(evs, block_size)
        specs = []
        for i, measurement in enumerate(measurements):
            analysis = self._analysis(measurement["name"])
            # The legacy NUISANCE2 record's add_to_framegen does not register the
            # sample's columns, so add the selection and projection operators
            # explicitly. The names are ours: a sample's own fname would collide with
            # another's in a shared frame.
            selection = analysis.get_selection()
            projections = analysis.get_projections()
            sel_col = f"sel{i}"
            proj_cols = [f"proj{i}_{j}" for j in range(len(projections))]
            fg.add_int_column(sel_col, selection.op)
            for name, projection in zip(proj_cols, projections):
                fg.add_double_column(name, projection.op)

            binned = analysis.get_data()[0]
            nbins = int(np.asarray(binned.values).reshape(-1).shape[0])
            specs.append({
                "measurement": measurement,
                "binning": binned.binning,
                "nbins": nbins,
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
            missing = [c for c in ["weight.cv"] +
                       [c for s in specs for c in [s["sel"], *s["projs"]]]
                       if c not in cols]
            if missing:
                raise RuntimeError(f"event frame is missing columns {missing}; "
                                   f"got {list(cols)}")
            fatx_col = cols.get("fatx_per_sumw.pb_per_target.estimate")
            if fatx_col is None:
                raise RuntimeError("event frame carries no per-target fatx estimate; "
                                   f"got {list(cols)}")
            # A running estimate over the file, so the newest block's last row is the
            # one to keep.
            fatx_per_sumw = float(table[-1, fatx_col])

            weight_col = cols["weight.cv"]
            for spec in specs:
                selected = table[table[:, cols[spec["sel"]]] != 0]
                if not selected.size:
                    continue
                proj_cols = [cols[c] for c in spec["projs"]]
                binning, nbins = spec["binning"], spec["nbins"]
                for row in selected:
                    values = [float(row[c]) for c in proj_cols]
                    b = (binning.find_bin(values[0]) if len(values) == 1
                         else binning.find_bin(pn.Vector_double(values)))
                    if b is None or b < 0 or b >= nbins:
                        continue  # under/overflow: outside the published binning
                    spec["bins"].append(int(b))
                    spec["weights"].append(float(row[weight_col]))
            block = fg.next()

        return {spec["measurement"]["name"]: self._finish(spec, generated,
                                                          fatx_per_sumw)
                for spec in specs}

    def _finish(self, spec: dict, generated: GeneratedEvents,
                fatx_per_sumw: Optional[float]) -> EventSample:
        """Scale one measurement's binned events into the published cross section."""
        nbins = spec["nbins"]
        bin_index = np.asarray(spec["bins"], dtype=int)
        weights = np.asarray(spec["weights"], dtype=float)
        # Keep the generator's own weights: the scaling below folds in the bin width,
        # which would otherwise show up as weight spread in the unweighting metrics.
        raw_weights = weights.copy()

        if not bin_index.size:
            return EventSample(bin_index=bin_index, weights=weights, nbins=nbins,
                               raw_weights=raw_weights)

        # Cross-section normalisation, per the recipe that reproduces the published
        # results: take the flux-averaged total xsec *per atom* and divide by the
        # target's nucleon count ourselves.
        #
        # NUISANCE's own PerNucleon conversion must NOT be used here: for a composite
        # target it divides by the struck nucleus' A, which for CH measures out at
        # exactly 12 (carbon only), whereas the published data uses 13 (CH). That is a
        # silent 13/12 = 8% normalisation error.
        if generated.xsec_divisor is None:
            raise ValueError("GeneratedEvents.xsec_divisor is required to normalise; "
                             "set data_per (and target_nucleons) on the experiment")
        if fatx_per_sumw is None:
            raise RuntimeError("no per-target fatx estimate was seen in any block")
        # Use the frame's own per-TARGET estimate: it is normalised against the same
        # weight.cv column used above. (EventSource.norm_info reports a differently
        # normalised sumweights -- 0.45 where the frame's weights sum to ~92000 --
        # so mixing the two overstates the prediction by ~4 orders of magnitude.)
        widths = np.asarray(list(spec["binning"].bin_sizes()),
                            dtype=float).reshape(-1)[:nbins]
        scale = np.divide(
            fatx_per_sumw * self._PB_TO_CM2 / generated.xsec_divisor, widths,
            out=np.zeros(nbins), where=widths != 0)
        scale = scale * self._extra_bin_scale(spec["measurement"], nbins)

        return EventSample(bin_index=bin_index, weights=weights * scale[bin_index],
                           nbins=nbins, raw_weights=raw_weights)

    def data_table(self, measurement: dict) -> DataTable:
        pn = self._pn()
        analysis = self._analysis(measurement["name"])
        binned = analysis.get_data()[0]

        values = np.asarray(binned.values, dtype=float).reshape(-1)
        errors = np.asarray(binned.errors, dtype=float).reshape(-1)
        covariance = np.asarray(analysis.get_covariance_matrix(), dtype=float)

        # A shipped data table in the wrong units, put back on the prediction's
        # footing before anything reads it (see data_scale in the config).
        data_scale = float(measurement.get("data_scale", 1.0))
        if data_scale != 1.0:
            values = values * data_scale
            errors = errors * data_scale
            covariance = covariance * data_scale ** 2

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

        configured = measurement.get("bin_edges")
        if configured is not None:
            # A bin-number histogram carries no usable edges of its own; these are the
            # real ones, and what the prediction has been made differential in.
            edges = np.asarray(configured, dtype=float)
        else:
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


def unweighting_cap(weights: np.ndarray, options: dict) -> Optional[float]:
    """The weight cap an Achilles unweighting scheme settles on for ``weights``.

    Mirrors ``SortedWeightUnweighter::ComputeCap`` for each registered scheme, so the
    synthetic adapter's ``--dry-run`` scan behaves like the real one and the rules
    have a reference implementation outside C++. ``None`` means "no cap" (``None``
    unweighter: events keep their weights).
    """
    name = options.get("Name", "None")
    if name == "None":
        return None

    w = np.sort(np.abs(np.asarray(weights, dtype=float)))
    if w.size == 0:
        return None
    total = float(np.sum(w))

    if name == "Percentile":
        idx = min(int(w.size * float(options["percentile"]) / 100.0), w.size - 1)
        return float(w[idx])

    eps = float(options["epsilon"])
    target = eps * total
    # suffix[i] = sum of the i largest-or-equal weights, i.e. sum_{j>=i} w[j].
    suffix = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])

    if name == "TailFraction":
        # Smallest cap C with sum_{|w|>C} |w| <= eps*total.
        i = int(np.argmax(suffix[:w.size] <= target)) if np.any(
            suffix[:w.size] <= target) else w.size
        return float(w[i - 1]) if i > 0 else float(w[0])

    if name == "Excess":
        # Smallest cap C with sum_j max(|w_j|-C, 0) <= eps*total. Between w[i-1] and
        # w[i] the excess is suffix[i] - (N-i)*C, so solve that for the first i whose
        # own weight already satisfies the bound.
        n = w.size
        excess_at = suffix[:n] - (n - np.arange(n)) * w
        i = int(np.argmax(excess_at <= target)) if np.any(excess_at <= target) else n - 1
        cap = (suffix[i] - target) / float(n - i)
        return float(np.clip(cap, w[0], w[-1]))

    raise KeyError(f"unknown unweighting scheme {name!r}")


class SyntheticAdapter:
    """Dry-run/self-test stand-in: draws weighted events from a tunable Gaussian.

    ``experiment['dryrun']['feature_shift']`` shifts only the feature branch (a fake
    physics change); ``measurement['dryrun']['nbins']`` sets the histogram binning.
    An ``unweighting`` block is applied for real (see ``unweighting_cap``), so the
    scan's plumbing and its metrics can be checked without Achilles.
    """

    def __init__(self, base_seed: int = 0):
        self.base_seed = base_seed

    def _nbins(self, measurement: dict) -> int:
        return int(measurement.get("dryrun", {}).get("nbins", 12))

    # Every prediction is normalised to this total, mirroring the real adapter's
    # flux-averaged cross-section scaling: a scheme that throws away events must not
    # come out smaller, only noisier.
    _TOTAL = 2.0e5

    def _draw(self, shift: float, n_events: int, rng: np.random.Generator, *,
              tail: bool = False) -> GeneratedEvents:
        x = np.clip(rng.normal(0.5 + shift, 0.18, size=n_events), 0.0, 0.999)
        # A long overweight tail is what the schemes differ on, so the scan draws
        # lognormal weights; the branch-comparison path keeps the mild uniform ones.
        weights = (rng.lognormal(0.0, 0.9, size=n_events) if tail
                   else rng.uniform(0.5, 1.5, size=n_events))
        return GeneratedEvents(x=x, weights=weights)

    @staticmethod
    def _unweight(gen: GeneratedEvents, cap: float,
                  rng: np.random.Generator) -> GeneratedEvents:
        """Accept with probability |w|/cap; overweights keep their excess."""
        prob = np.abs(gen.weights) / cap
        keep = prob >= rng.uniform(0.0, 1.0, size=prob.size)
        return GeneratedEvents(x=gen.x[keep], weights=np.maximum(prob[keep], 1.0))

    def generate(self, experiment: dict, branch: str, seed: int, n_events: int, *,
                 unweighting: Optional[dict] = None,
                 seed_offset: Optional[int] = None) -> GeneratedEvents:
        offset = (seed_offset if seed_offset is not None
                  else (0 if branch == "main" else 1))
        rng = np.random.default_rng(self.base_seed + seed + offset)
        shift = float(experiment.get("dryrun", {}).get("feature_shift", 0.0)) \
            if branch == "feature" else 0.0

        # One pilot draw stands in for Achilles' optimisation pass: it fixes the cap
        # and, with it, the acceptance rate.
        pilot = self._draw(shift, min(n_events, 100_000), rng, tail=True)
        cap = (unweighting_cap(pilot.weights, unweighting)
               if unweighting is not None else None)

        if cap is None or not cap > 0.0:
            gen = self._draw(shift, n_events, rng, tail=unweighting is not None)
            eff = 1.0
        else:
            # Achilles generates until it has n_events *accepted*, so oversample by
            # the acceptance rate rather than letting a harsher cap yield fewer
            # events -- otherwise the schemes are compared at different statistics.
            eff = float(np.mean(np.minimum(np.abs(pilot.weights) / cap, 1.0)))
            trials = int(n_events / max(eff, 1e-3) * 1.15) + 1000
            gen = self._unweight(self._draw(shift, trials, rng, tail=True), cap, rng)
            gen = GeneratedEvents(x=gen.x[:n_events], weights=gen.weights[:n_events])

        gen.run = RunStats(seconds=float(n_events) / 5e4 / max(eff, 1e-3),
                           unweight_eff=eff)
        return gen

    def histogram(self, generated: GeneratedEvents,
                  measurement: dict) -> EventSample:
        nbins = self._nbins(measurement)
        bin_index = np.digitize(generated.x, np.linspace(0.0, 1.0, nbins + 1)) - 1
        total = float(np.sum(np.abs(generated.weights)))
        scale = self._TOTAL / total if total > 0.0 else 1.0
        return EventSample(bin_index=bin_index, weights=generated.weights * scale,
                           nbins=nbins, raw_weights=generated.weights)

    def histogram_many(self, generated: GeneratedEvents, measurements,
                       block_size: int = 250_000) -> "dict[str, EventSample]":
        # In-memory events: nothing to save by sharing a pass, but the driver calls
        # only this.
        return {m["name"]: self.histogram(generated, m) for m in measurements}

    @staticmethod
    def response_matrix(measurement: dict, nbins: int) -> Optional[np.ndarray]:
        return None  # the synthetic path has no smeared measurements

    def data_table(self, measurement: dict) -> DataTable:
        nbins = self._nbins(measurement)
        rng = np.random.default_rng(self.base_seed + 999)
        gen = self._draw(0.0, 200000, rng)  # data = the unshifted truth
        bin_index = np.digitize(gen.x, np.linspace(0.0, 1.0, nbins + 1)) - 1
        values = np.bincount(bin_index, weights=gen.weights, minlength=nbins)[:nbins]
        values = values * (self._TOTAL / float(np.sum(values)))
        err = np.sqrt(np.maximum(values, 1.0)) * 0.05 + 0.02 * values
        return DataTable(values=values, covariance=np.diag(err ** 2),
                         edges=np.linspace(0.0, 1.0, nbins + 1), xlabel="x")


# ---------------------------------------------------------------------------
# Self-test: the batched frame walk, against a stand-in for pyNUISANCE
# ---------------------------------------------------------------------------

def _selftest() -> int:
    """Check histogram_many over a fake pyNUISANCE, since the real one needs the image.

    What is worth pinning down here is the frame bookkeeping: that each sample reads
    its own selection and projection columns out of a shared frame, that blocks are
    stitched together, and that the fatx estimate comes from the last block.
    """
    class Binning:
        def __init__(self, edges):
            self.edges = np.asarray(edges, dtype=float)
        def find_bin(self, value):
            b = int(np.digitize([value], self.edges)[0]) - 1
            return b if 0 <= b < self.edges.size - 1 else -1
        def bin_sizes(self):
            return list(np.diff(self.edges))

    class Binned:
        def __init__(self, edges):
            self.binning = Binning(edges)
            self.values = np.zeros(len(edges) - 1)

    def column(tag):
        """Stand-in for a sample's selection/projection handle (.fname, .op)."""
        return types.SimpleNamespace(fname=tag, op=object(), prettyname=tag, units="")

    class Analysis:
        def __init__(self, edges):
            self._binned = Binned(edges)
        def get_selection(self): return column("sel")
        def get_projections(self): return [column("proj")]
        def get_data(self): return [self._binned]

    class Frame:
        def __init__(self, names, table):
            self.column_names = names
            self.table = table

    class FrameGen:
        """Serves pre-baked blocks; column order follows the add_* call order."""
        def __init__(self, blocks): self.blocks, self.i, self.names = blocks, 0, []
        def add_int_column(self, name, op): self.names.append(name)
        def add_double_column(self, name, op): self.names.append(name)
        def first(self): self.i = 0; return self.next()
        def next(self):
            if self.i >= len(self.blocks):
                return Frame(self._names(), np.empty((0, len(self._names()))))
            block = self.blocks[self.i]; self.i += 1
            return Frame(self._names(), block)
        def _names(self):
            return ["weight.cv", "fatx_per_sumw.pb_per_target.estimate"] + self.names

    # Two samples on one frame: [weight, fatx, sel0, proj0, sel1, proj1].
    # Bin edges are 0,1,2,3 for both, so a projection value is its own bin.
    blocks = [np.array([
        # w    fatx  sel0 proj0  sel1 proj1
        [2.0,  10.0,  1,   0.5,   1,   2.5],   # both select, different bins
        [3.0,  10.0,  0,   9.9,   1,   1.5],   # only the second selects
        [1.0,  10.0,  1,   7.0,   0,   0.5],   # first selects, out of range
    ]), np.array([
        [4.0,  20.0,  1,   2.5,   0,   0.5],   # second block; fatx has moved
    ])]

    adapter = Nuisance3Adapter.__new__(Nuisance3Adapter)
    fake_pn = types.SimpleNamespace(
        EventSource=lambda path: object(),
        EventFrameGen=lambda evs, block: FrameGen(blocks),
        Vector_double=list)
    adapter._pn = lambda: fake_pn
    adapter._analysis = lambda name: Analysis([0.0, 1.0, 2.0, 3.0])

    gen = GeneratedEvents(path="fake.hepmc", xsec_divisor=2.0)
    got = adapter.histogram_many(gen, [{"name": "A"}, {"name": "B"}])

    # fatx from the LAST block (20), per-atom divisor 2, unit bin widths, pb->cm^2.
    scale = 20.0 * Nuisance3Adapter._PB_TO_CM2 / 2.0
    checks = {
        "both samples returned": sorted(got) == ["A", "B"],
        "A took its own rows": list(got["A"].bin_index) == [0, 2],
        "A skipped the out-of-range row": got["A"].bin_index.size == 2,
        "A weights scaled": np.allclose(got["A"].weights, np.array([2.0, 4.0]) * scale),
        "A raw weights untouched": np.allclose(got["A"].raw_weights, [2.0, 4.0]),
        "B took its own rows": list(got["B"].bin_index) == [2, 1],
        "B weights scaled": np.allclose(got["B"].weights, np.array([2.0, 3.0]) * scale),
        "nbins from the binning": got["A"].nbins == 3 and got["B"].nbins == 3,
    }

    # solid_angle/bin_edges ride on top of the same scaling.
    wide = adapter.histogram_many(gen, [{"name": "A", "solid_angle": 4.0}])["A"]
    checks["solid_angle divides"] = np.allclose(
        wide.weights, np.array([2.0, 4.0]) * scale / 4.0)
    single = adapter.histogram(gen, {"name": "A"})
    checks["histogram() matches the batch"] = np.allclose(single.weights,
                                                          got["A"].weights)
    for name, ok in checks.items():
        print(f"[{'ok' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        raise SystemExit(_selftest())
    print(__doc__)
