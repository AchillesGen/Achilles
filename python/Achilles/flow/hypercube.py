"""Unit-hypercube spline flows (i-flow style) and a VEGAS-then-flow warm start.

:class:`HypercubeSplineFlow` puts rational-quadratic spline coupling layers directly on
``[0, 1]^d`` with a uniform base, so there are no linear tails and no ``Φ`` squash: the
whole cube, including its edges, is covered by spline bins. The couplings use the binary
mask set (dimension ``i`` is transformed when bit ``b`` of ``i`` is set, plus the
complements), so every dimension is conditioned on every other one at least once in
``2 * ceil(log2 d)`` layers. The last layer of each conditioner is zero-initialized, which
makes the untrained flow the identity, i.e. pure analytic-``Mapper`` sampling.

:class:`VegasSplineFlow` first adapts a per-channel VEGAS grid (the same algorithm as the
C++ ``AdaptiveMap``) for ``vegas_iterations`` training calls, then freezes it and trains
the spline flow on top of it: ``rans = V(u)`` with ``u`` drawn from the flow. The flow is
still the identity at the switch, so training starts from the VEGAS-adapted density.
"""

import math

import numpy as np

from ._torch_base import ChannelFlow, _EPS
from .vegas import VegasWarmStart


def _require_splines():
    try:
        from normflows.utils.splines import rational_quadratic_spline, DEFAULT_MIN_DERIVATIVE
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "HypercubeSplineFlow requires normflows. Install the flow extra with "
            "`pip install achilles[flow]`."
        ) from exc
    return rational_quadratic_spline, DEFAULT_MIN_DERIVATIVE


def binary_masks(ndims):
    """Masks (1 = transformed) such that every ordered pair of dimensions is covered."""
    if ndims == 1:
        return [[1]]
    masks = []
    for bit in range(math.ceil(math.log2(ndims))):
        mask = [(i >> bit) & 1 for i in range(ndims)]
        masks += [mask, [1 - m for m in mask]]
    return masks


def _module_classes(torch):
    rqs, min_derivative = _require_splines()
    # Raw derivative of zero maps to a knot slope of exactly one (identity spline).
    deriv_offset = math.log(math.expm1(1.0 - min_derivative))

    class ResidualBlock(torch.nn.Module):
        def __init__(self, hidden):
            super().__init__()
            self.net = torch.nn.Sequential(
                torch.nn.SiLU(), torch.nn.Linear(hidden, hidden),
                torch.nn.SiLU(), torch.nn.Linear(hidden, hidden),
            )

        def forward(self, x):
            return x + self.net(x)

    class RQSCoupling(torch.nn.Module):
        def __init__(self, mask, hidden, bins, blocks):
            super().__init__()
            mask = torch.as_tensor(mask, dtype=torch.bool)
            self.register_buffer("trans", mask.nonzero(as_tuple=True)[0])
            self.register_buffer("cond", (~mask).nonzero(as_tuple=True)[0])
            self.bins = int(bins)
            n_out = self.trans.numel() * (3 * self.bins + 1)
            self.net = torch.nn.Sequential(
                torch.nn.Linear(max(1, self.cond.numel()), hidden),
                *(ResidualBlock(hidden) for _ in range(blocks)),
                torch.nn.SiLU(),
                torch.nn.Linear(hidden, n_out),
            )
            torch.nn.init.zeros_(self.net[-1].weight)
            torch.nn.init.zeros_(self.net[-1].bias)

        def forward(self, x, inverse=False):
            n = x.shape[0]
            if self.cond.numel() == 0:
                cond = torch.zeros(n, 1, dtype=x.dtype)
            else:
                cond = 2.0 * x[:, self.cond] - 1.0
            params = self.net(cond).view(n, self.trans.numel(), 3 * self.bins + 1)
            widths = params[..., : self.bins]
            heights = params[..., self.bins : 2 * self.bins]
            derivs = params[..., 2 * self.bins :] + deriv_offset
            out, logdet = rqs(
                x[:, self.trans].clamp(0.0, 1.0), widths, heights, derivs, inverse=inverse
            )
            y = x.clone()
            y[:, self.trans] = out.clamp(0.0, 1.0)
            return y, logdet.sum(dim=1)

    class HypercubeCore(torch.nn.Module):
        """Density direction ``x -> u`` through the layers; sampling runs them inverted."""

        def __init__(self, ndims, hidden, bins, blocks, repeats):
            super().__init__()
            self.ndims = ndims
            masks = binary_masks(ndims) * max(1, int(repeats))
            self.layers = torch.nn.ModuleList(RQSCoupling(m, hidden, bins, blocks) for m in masks)

        def log_prob(self, x):
            logdet = torch.zeros(x.shape[0], dtype=x.dtype)
            for layer in self.layers:
                x, ld = layer(x)
                logdet = logdet + ld
            return logdet  # uniform base: log q(x) = log|du/dx|

        def sample(self, n):
            u = torch.rand(n, self.ndims)
            log_q = torch.zeros(n, dtype=u.dtype)
            for layer in reversed(self.layers):
                u, ld = layer(u, inverse=True)  # ld = log|dx/du|
                log_q = log_q - ld
            return u, log_q

    return HypercubeCore


class HypercubeSplineFlow(ChannelFlow):
    """Rational-quadratic spline coupling flow living directly on the unit hypercube."""

    def _build_cores(self, hidden=64, bins=16, res_blocks=2, repeats=1, conditional=False,
                     **_unused):
        if conditional:
            raise ValueError(
                f"{self.__class__.__name__} does not support a channel-conditional flow; "
                "use one flow per channel (conditional: false)."
            )
        torch = self.torch
        Core = _module_classes(torch)
        return torch.nn.ModuleList(
            Core(self.ndims, int(hidden), int(bins), int(res_blocks), int(repeats))
            for _ in range(self.nchannels)
        )

    def _core_sample(self, channel, n):
        return self.cores[channel].sample(n)[0]

    def _core_sample_and_log_prob(self, channel, n):
        return self.cores[channel].sample(n)

    def _core_log_prob(self, channel, z):
        return self.cores[channel].log_prob(z)

    # The core already lives on the unit cube: no squash.
    def latent(self, channel, rans_t):
        return rans_t.clamp(_EPS, 1.0 - _EPS)

    def unit_cube(self, channel, z):
        return z.clamp(_EPS, 1.0 - _EPS)

    def _squash_log_det(self, channel, z):
        return self.torch.zeros(z.shape[0], dtype=z.dtype)


class VegasSplineFlow(VegasWarmStart, HypercubeSplineFlow):
    """VEGAS grid adaptation for the first ``vegas_iterations`` calls, then a spline flow."""
