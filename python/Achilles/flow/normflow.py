"""Coupling-layer flow (``normflows``) backend.

A per-channel RealNVP-style coupling flow in logit space, squashed to the unit cube by
the shared :class:`ChannelFlow` base. Coupling flows give exact, cheap ``log_prob``, which
is ideal for the importance-sampling weight -- the natural baseline to compare continuous
/ flow-matching backends against.
"""

from ._torch_base import ChannelFlow, checkerboard_mask


def _require_normflows():
    try:
        import normflows as nf  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "NormFlow requires normflows. Install the flow extra with "
            "`pip install achilles[flow]`."
        ) from exc
    return nf


class NormFlow(ChannelFlow):
    def _build_cores(self, hidden=64, blocks=4, conditional=False, **_unused):
        if conditional:
            raise ValueError(
                "NormFlow does not support a channel-conditional flow; use one flow per "
                "channel (conditional: false), or use ZukoFlow/FlowMatchingFlow to "
                "benchmark a conditional flow."
            )
        torch = self.torch
        nf = _require_normflows()

        def make_flow():
            layers = []
            for i in range(int(blocks)):
                param_map = nf.nets.MLP(
                    [self.ndims, int(hidden), int(hidden), self.ndims], init_zeros=True
                )
                layers.append(
                    nf.flows.MaskedAffineFlow(checkerboard_mask(torch, self.ndims, i), param_map)
                )
                layers.append(nf.flows.ActNorm(self.ndims))
            base = nf.distributions.DiagGaussian(self.ndims)
            return nf.NormalizingFlow(base, layers)

        return torch.nn.ModuleList([make_flow() for _ in range(self.nchannels)])

    def _core_sample(self, channel, n):
        z, _ = self.cores[channel].sample(n)
        return z

    def _core_log_prob(self, channel, z):
        return self.cores[channel].log_prob(z)


class SplineFlow(ChannelFlow):
    """Rational-quadratic neural spline flow (RQS) via ``normflows``.

    Uses RQS coupling layers (``CoupledRationalQuadraticSpline``) instead of affine
    coupling -- the i-flow-style transform, more expressive for sharply peaked integrands.
    It sits in the same Gaussian-latent / Φ-squash framework as the other backends, so it
    has the same effective uniform base on the unit cube and starts from the analytic
    Mapper. (``zuko``'s ``NSF`` is the zuko-side RQS; this is the normflows counterpart.)
    """

    def _build_cores(self, hidden=64, blocks=4, bins=8, conditional=False, **_unused):
        if conditional:
            raise ValueError(
                "SplineFlow does not support a channel-conditional flow; use one flow per "
                "channel (conditional: false), or use ZukoFlow/FlowMatchingFlow."
            )
        torch = self.torch
        nf = _require_normflows()

        def make_flow():
            layers = []
            for i in range(int(blocks)):
                layers.append(
                    nf.flows.CoupledRationalQuadraticSpline(
                        self.ndims, 2, int(hidden), num_bins=int(bins), reverse_mask=bool(i % 2)
                    )
                )
            base = nf.distributions.DiagGaussian(self.ndims)
            return nf.NormalizingFlow(base, layers)

        return torch.nn.ModuleList([make_flow() for _ in range(self.nchannels)])

    def _core_sample(self, channel, n):
        z, _ = self.cores[channel].sample(n)
        return z

    def _core_log_prob(self, channel, z):
        return self.cores[channel].log_prob(z)
