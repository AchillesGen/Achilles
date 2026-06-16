"""Flow-matching backend built on ``facebookresearch/flow_matching``.

Per channel we train a velocity field and integrate it as a continuous flow:

* ``sample``  -- integrate the ODE forward from a Gaussian base (``ODESolver.sample``);
* ``density`` -- exact log-likelihood via the instantaneous change of variables
  (``ODESolver.compute_likelihood(..., exact_divergence=True)``). Exact divergence is used
  on purpose: the Hutchinson estimator is stochastic, and since the importance weight is
  ``f/g`` (not ``log g``) a noisy density biases unweighting. This exact ODE likelihood is
  the cost to watch in a generation-time benchmark.

Training (overrides the base importance-weighted MLE, which is impractical through the
ODE) uses **weighted conditional flow matching**: the velocity is regressed toward the
straight-line (CondOT) path from the Gaussian base to the importance-weighted samples, so
the learned distribution moves toward ``p ∝ f``. This reweighted objective is the
experimental piece to validate/tune for the benchmark.

NOTE: targets the current ``flow_matching`` API (``flow_matching.solver.ODESolver`` with
``sample`` / ``compute_likelihood`` and ``flow_matching.utils.ModelWrapper``); adjust if
your installed version differs.
"""

import math

import numpy as np

from ._torch_base import ChannelFlow, _EPS


def _require_flow_matching():
    try:
        import flow_matching  # noqa: F401
        from flow_matching.solver import ODESolver  # noqa: F401
        from flow_matching.utils import ModelWrapper  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "FlowMatchingFlow requires the flow_matching package. Install the flow extra "
            "with `pip install achilles[flow]`."
        ) from exc
    return ODESolver, ModelWrapper


class FlowMatchingFlow(ChannelFlow):
    def _build_cores(self, hidden=64, step_size=0.05, method="midpoint", conditional=False,
                     **_unused):
        torch = self.torch
        ODESolver, ModelWrapper = _require_flow_matching()
        self._ODESolver = ODESolver
        self._step_size = float(step_size)
        self._method = str(method)
        # conditional: one shared velocity net conditioned on the channel one-hot, vs one
        # independent net per channel -- a benchmark axis.
        self._conditional = bool(conditional)

        ndims = self.ndims
        cond_dim = self.nchannels if self._conditional else 0

        class VelocityNet(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.net = torch.nn.Sequential(
                    torch.nn.Linear(ndims + 1 + cond_dim, int(hidden)),
                    torch.nn.SiLU(),
                    torch.nn.Linear(int(hidden), int(hidden)),
                    torch.nn.SiLU(),
                    torch.nn.Linear(int(hidden), ndims),
                )

            def forward(self, x, t, cond=None):
                if not torch.is_tensor(t):
                    t = torch.as_tensor(t, dtype=x.dtype)
                t = t.reshape(-1, 1).expand(x.shape[0], 1).to(x.dtype)
                parts = [x, t] if cond is None else [x, t, cond]
                return self.net(torch.cat(parts, dim=-1))

        # The solver calls model(x, t); a per-channel wrapper injects that channel's
        # one-hot (conditional) or simply forwards to that channel's net (per-channel).
        flow = self

        class Wrapper(ModelWrapper):
            def __init__(self, model, channel):
                super().__init__(model)
                self._channel = channel

            def forward(self, x, t, **extras):
                return flow._velocity(self._channel, x, t)

        ncores = 1 if self._conditional else self.nchannels
        cores = torch.nn.ModuleList([VelocityNet() for _ in range(ncores)])
        self._wrappers = [Wrapper(cores[0], c) for c in range(self.nchannels)]
        return cores

    def _velocity(self, channel, x, t):
        # Dispatch to the right net, injecting the channel one-hot when conditional.
        if self._conditional:
            cond = self._channel_onehot(channel, x.shape[0])
            return self.cores[0](x, t, cond)
        return self.cores[channel](x, t)

    # Standard-normal base log-density, used by compute_likelihood.
    def _log_p0(self, x):
        torch = self.torch
        return -0.5 * (x**2).sum(dim=-1) - 0.5 * self.ndims * math.log(2.0 * math.pi)

    def _solver(self, channel):
        return self._ODESolver(velocity_model=self._wrappers[channel])

    def _core_sample(self, channel, n):
        torch = self.torch
        x0 = torch.randn(n, self.ndims)
        return self._solver(channel).sample(
            x_init=x0, step_size=self._step_size, method=self._method
        )

    def _core_log_prob(self, channel, z):
        result = self._solver(channel).compute_likelihood(
            x_1=z,
            log_p0=self._log_p0,
            step_size=self._step_size,
            method=self._method,
            exact_divergence=True,
        )
        # compute_likelihood returns (x_0, log_likelihood); accept a bare tensor too.
        return result[1] if isinstance(result, (tuple, list)) else result

    def train(self, channels, rans, values):
        torch = self.torch
        channels = np.asarray(channels).astype(np.int64)
        rans_t = torch.as_tensor(np.asarray(rans, dtype=np.float64), dtype=torch.float32)
        rans_t = rans_t.clamp(_EPS, 1.0 - _EPS)
        z1_all = self.latent(rans_t)  # targets in latent space (probit, matches the base squash)
        weights = self._normalized_weights(values)
        if weights is None:
            return 0.0

        self.optimizer.zero_grad()

        # Discrete part: importance-weighted categorical NLL on the channel logits.
        ch_t = torch.as_tensor(channels, dtype=torch.long)
        loss = -(weights * torch.log_softmax(self.logits, dim=0)[ch_t]).sum()

        # Continuous part: weighted conditional flow matching per channel.
        for c in range(self.nchannels):
            idx = np.nonzero(channels == c)[0]
            if idx.size == 0:
                continue
            z1 = z1_all[idx]
            w = weights[idx]
            t = torch.rand(z1.shape[0], 1)
            x0 = torch.randn_like(z1)
            xt = (1.0 - t) * x0 + t * z1          # CondOT straight-line path
            target = z1 - x0                       # path velocity
            pred = self._velocity(c, xt, t.squeeze(-1))
            loss = loss + (w * ((pred - target) ** 2).mean(dim=1)).sum()

        loss.backward()
        self.optimizer.step()
        return float(loss.detach().cpu().item())
