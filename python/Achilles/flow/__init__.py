"""Normalizing-flow sampling backends for the Achilles ``Flow`` integrator.

The C++ ``FlowIntegrator`` (via ``PyFlowSampler``) drives a flow object that replaces the
per-channel Vegas grids and channel-selection weights of the default ``MultiChannel``
integrator, while reusing the existing analytic channel ``Mapper``s. The flow is selected
from the run card and may be a built-in class, or any class you supply in your own file::

    Options:
      Integrator:
        Type: Flow
        Flow:
          File: /path/to/my_flow.py   # optional: load the class from this file
          Class: MyFlow               # class name (default NormFlow)
          Epochs: 2000                # max training iterations
          NCalls: 20000               # samples per iteration
          # ... any other keys are forwarded to the class constructor as kwargs
          hidden: 128
          lr: 1.0e-3

Protocol
--------
A flow is **type-agnostic** -- coupling-layer, continuous (CNF), flow-matching, or a
pre-trained model loaded from disk all work, as long as the class accepts
``(ndims, nchannels, **params)`` and implements these **batched** methods (arrays are
:class:`numpy.ndarray`; ``n`` is the batch size and ``d == ndims``):

``sample(n) -> (channels[n] int, rans[n, d] float, density[n] float)``   (required)
    Draw ``n`` samples: selected channel, unit-hypercube coordinates for that channel's
    ``Mapper``, and the sampling pdf ``g(channel, rans) = P(channel) * p(rans | channel)``.

``density(channels[n], rans[n, d]) -> density[n]``                       (required)
    The same pdf evaluated on externally supplied points.

``train(channels[n], rans[n, d], values[n]) -> loss``                   (optional)
    One optimization step; ``values`` are the integrand weights ``f/g``. Defaults to a
    no-op for inference-only models.

``save(path)`` / ``load(path) -> bool``                                 (optional)
    Persist / restore parameters. Default to no-op / ``False``.

``requires_training``                                                   (optional attr)
    ``False`` for a pre-trained / inference-only model so the integrator only estimates
    the cross section instead of running the training loop. Defaults to ``True``.

Subclassing :class:`Flow` gives you the optional methods for free, but it is not
required -- ``load_flow`` validates the protocol by duck typing.
"""

import abc
import importlib.util
import sys

import numpy as np


class Flow(abc.ABC):
    """Convenience base class implementing the optional half of the flow protocol.

    Subclassing is optional; only :meth:`sample` and :meth:`density` are required.
    """

    requires_training = True

    @abc.abstractmethod
    def sample(self, n):
        ...

    @abc.abstractmethod
    def density(self, channels, rans):
        ...

    def train(self, channels, rans, values):
        return 0.0

    def save(self, path):
        return None

    def load(self, path):
        return False


class UniformFlow(Flow):
    """Trivial, dependency-free reference backend.

    Samples channels uniformly and ``rans`` uniformly on the unit hypercube, so the
    ``FlowIntegrator`` reduces to plain multichannel Monte Carlo with no learned mapping.
    Useful for smoke-testing the C++/Python plumbing without PyTorch, and the behavioural
    analogue of the C++ ``MockFlowSampler`` used in the unit tests.
    """

    requires_training = False

    def __init__(self, ndims, nchannels, **_unused):
        self.ndims = int(ndims)
        self.nchannels = int(nchannels)
        self._density = 1.0 / self.nchannels  # uniform over channels and the unit cube

    def sample(self, n):
        n = int(n)
        channels = np.random.randint(0, self.nchannels, size=n).astype(np.int64)
        rans = np.random.random_sample((n, self.ndims)).astype(np.float64)
        density = np.full(n, self._density, dtype=np.float64)
        return channels, rans, density

    def density(self, channels, rans):
        channels = np.asarray(channels)
        return np.full(channels.shape[0], self._density, dtype=np.float64)


_REQUIRED = ("sample", "density")
_OPTIONAL_DEFAULTS = {
    "train": lambda *a, **k: 0.0,
    "save": lambda *a, **k: None,
    "load": lambda *a, **k: False,
}


def load_flow(file, class_name, ndims, nchannels, **params):
    """Load, instantiate and validate a flow class.

    If ``file`` is given the class is loaded from that Python file; otherwise it is taken
    from the built-in classes in this module. Raises ``TypeError`` if the instance does
    not implement the required protocol; fills in no-op defaults for missing optional
    methods so the C++ side can always call them.
    """
    module = _load_user_module(file) if file else sys.modules[__name__]
    try:
        cls = getattr(module, class_name)
    except AttributeError as exc:
        where = f"file '{file}'" if file else "Achilles.flow"
        raise AttributeError(f"Flow class '{class_name}' not found in {where}.") from exc

    instance = cls(int(ndims), int(nchannels), **params)

    missing = [name for name in _REQUIRED if not callable(getattr(instance, name, None))]
    if missing:
        raise TypeError(
            f"Flow class '{class_name}' does not implement required method(s): "
            f"{', '.join(missing)}. A flow must provide sample(n) and "
            "density(channels, rans)."
        )
    for name, default in _OPTIONAL_DEFAULTS.items():
        if not callable(getattr(instance, name, None)):
            setattr(instance, name, default)
    return instance


def _load_user_module(path):
    spec = importlib.util.spec_from_file_location("achilles_user_flow", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load a Python module from '{path}'.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def __getattr__(name):
    # Lazily expose the heavier backends so importing this package does not require
    # torch/normflows, zuko, flow_matching or onnxruntime unless that backend is used.
    if name == "NormFlow":
        from .normflow import NormFlow

        return NormFlow
    if name == "SplineFlow":
        from .normflow import SplineFlow

        return SplineFlow
    if name == "ZukoFlow":
        from .zukoflow import ZukoFlow

        return ZukoFlow
    if name == "FlowMatchingFlow":
        from .flowmatching import FlowMatchingFlow

        return FlowMatchingFlow
    if name == "OnnxFlow":
        from .onnxflow import OnnxFlow

        return OnnxFlow
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "Flow",
    "UniformFlow",
    "NormFlow",
    "SplineFlow",
    "ZukoFlow",
    "FlowMatchingFlow",
    "OnnxFlow",
    "load_flow",
]
