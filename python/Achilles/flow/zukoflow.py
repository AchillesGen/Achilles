"""``zuko`` flow backend.

``zuko`` exposes a wide range of flows -- coupling/spline (``NSF``, ``MAF``) and
continuous (``CNF``) -- behind one API: ``flow()`` builds a distribution with ``sample``
and exact ``log_prob``. The flow class is configurable so the same backend can benchmark
several architectures::

    Flow:
      Class: ZukoFlow
      flow: CNF          # any class in zuko.flows (CNF, NSF, MAF, ...)
      hidden: 64
      transforms: 3      # ignored by CNF
      conditional: false # false: one flow per channel; true: one flow conditioned on
                         #        the channel (one-hot context) -- a benchmark axis

Note: continuous flows (``CNF``) evaluate ``log_prob`` by integrating an ODE with the
divergence trace, which is far more expensive than coupling flows -- the cost to watch in
a generation-time benchmark, since density is queried for every point and channel.
"""

from ._torch_base import ChannelFlow


def _require_zuko():
    try:
        import zuko  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "ZukoFlow requires zuko. Install the flow extra with "
            "`pip install achilles[flow]`."
        ) from exc
    return zuko


class ZukoFlow(ChannelFlow):
    def _build_cores(self, flow="CNF", hidden=64, transforms=3, conditional=False, **_unused):
        torch = self.torch
        zuko = _require_zuko()
        flow_cls = getattr(zuko.flows, flow)
        hidden_features = (int(hidden), int(hidden))
        self._conditional = bool(conditional)
        # A conditional flow takes the channel one-hot as zuko `context`; per-channel uses
        # independent unconditional flows.
        context = self.nchannels if self._conditional else 0

        def make_flow():
            kwargs = {"features": self.ndims, "context": context, "hidden_features": hidden_features}
            try:
                return flow_cls(transforms=int(transforms), **kwargs)  # coupling/spline
            except TypeError:
                return flow_cls(**kwargs)  # CNF and others without `transforms`

        ncores = 1 if self._conditional else self.nchannels
        return torch.nn.ModuleList([make_flow() for _ in range(ncores)])

    def _core_sample(self, channel, n):
        if self._conditional:
            return self.cores[0](self._channel_onehot(channel, n)).sample()
        return self.cores[channel]().sample((n,))

    def _core_log_prob(self, channel, z):
        if self._conditional:
            return self.cores[0](self._channel_onehot(channel, z.shape[0])).log_prob(z)
        return self.cores[channel]().log_prob(z)
