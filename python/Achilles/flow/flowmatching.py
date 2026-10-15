"""Flow-matching backend built on ``facebookresearch/flow_matching``.

Per channel we train a velocity field and integrate it from a Gaussian base:

* fixed-step methods (``euler``, ``midpoint``): the sampler is a discrete map, and both
  ``sample`` and ``density`` use its exact density (per-step Jacobian log-dets; ``density``
  inverts the steps by Newton's method). The continuous ODE likelihood would differ by
  O(h^p) and bias ``f/g``;
* other (adaptive) methods: ``ODESolver.sample`` and the exact-divergence
  ``ODESolver.compute_likelihood``.

Training (overrides the base importance-weighted MLE, which is impractical through the
ODE) uses **weighted conditional flow matching**: the velocity is regressed toward the
straight-line (CondOT) path from the Gaussian base to the importance-weighted samples, so
the learned distribution moves toward ``p ∝ f``. Note ``loss_power`` changes the target to
``f^p q^(1-p)``: the fixed point is still ``q = f``, but it converges only for ``p < 2``.
:class:`VegasFlowMatchingFlow` (VEGAS warm start + the tuned settings) is the recommended
variant.

NOTE: targets the current ``flow_matching`` API (``flow_matching.solver.ODESolver`` with
``sample`` / ``compute_likelihood`` and ``flow_matching.utils.ModelWrapper``); adjust if
your installed version differs.
"""

import math
import warnings

import numpy as np

from ._torch_base import ChannelFlow, _EPS
from .vegas import VegasWarmStart


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
    # Objective is weighted conditional flow matching + categorical NLL, not the base MLE.
    loss_ylabel = "weighted CFM + categorical NLL loss"

    def _build_cores(self, hidden=64, step_size=0.05, method="midpoint", conditional=False,
                     train_steps=1, zero_init=False, time_freqs=0, layers=2, ot_batch=0,
                     buffer_size=0, **_unused):
        torch = self.torch
        ODESolver, ModelWrapper = _require_flow_matching()
        self._ODESolver = ODESolver
        self._step_size = float(step_size)
        self._method = str(method)
        self._train_steps = max(1, int(train_steps))
        self._ot_batch = int(ot_batch)
        # Replay buffer of past (channel, rans, f/g): f/g stays a valid weight for the CFM
        # regression after the flow moves, so old batches remain usable training data.
        self._buffer_size = int(buffer_size)
        self._buffer = None
        # conditional: one shared velocity net conditioned on the channel one-hot, vs one
        # independent net per channel -- a benchmark axis.
        self._conditional = bool(conditional)

        ndims = self.ndims
        cond_dim = self.nchannels if self._conditional else 0
        # Time enters as t plus sin/cos(pi k t), k = 1..time_freqs (Fourier features).
        freqs = math.pi * torch.arange(1, int(time_freqs) + 1, dtype=torch.float32)
        time_dim = 1 + 2 * freqs.numel()
        hidden = int(hidden)

        class VelocityNet(torch.nn.Module):
            def __init__(self):
                super().__init__()
                body = [torch.nn.Linear(ndims + time_dim + cond_dim, hidden), torch.nn.SiLU()]
                for _ in range(int(layers) - 1):
                    body += [torch.nn.Linear(hidden, hidden), torch.nn.SiLU()]
                self.net = torch.nn.Sequential(*body, torch.nn.Linear(hidden, ndims))

            def forward(self, x, t, cond=None):
                if not torch.is_tensor(t):
                    t = torch.as_tensor(t, dtype=x.dtype)
                t = t.reshape(-1, 1).expand(x.shape[0], 1).to(x.dtype)
                parts = [x, t, torch.sin(t * freqs), torch.cos(t * freqs)]
                if cond is not None:
                    parts.append(cond)
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
        if zero_init:
            # Zero velocity is the identity flow, so training starts from the analytic Mapper.
            for core in cores:
                torch.nn.init.zeros_(core.net[-1].weight)
                torch.nn.init.zeros_(core.net[-1].bias)
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

    # -- exact density of the discrete map ---------------------------------------------
    # With a fixed-step method the sampler is a discrete map, so its density is exactly
    # p0(x0) / prod_k |det d step_k / dx|. The continuous ODE likelihood differs from it by
    # O(h^p), which biases f/g (about -0.15% on the cross section at h = 0.05 midpoint).
    _DISCRETE = ("euler", "midpoint")

    def _grid(self):
        n = max(1, round(1.0 / self._step_size))
        return [(k / n, 1.0 / n) for k in range(n)]

    def _step(self, channel, x, t, h):
        tt = self.torch.full((x.shape[0],), t, dtype=x.dtype)
        v = self._velocity(channel, x, tt)
        if self._method == "euler":
            return x + h * v
        return x + h * self._velocity(channel, x + 0.5 * h * v, tt + 0.5 * h)

    def _step_jacobian(self, channel, x, t, h):
        func = self.torch.func
        single = lambda xi: self._step(channel, xi[None], t, h)[0]  # noqa: E731
        return func.vmap(func.jacrev(single))(x)

    def _core_sample_and_log_prob(self, channel, n):
        if self._method not in self._DISCRETE:
            return super()._core_sample_and_log_prob(channel, n)
        torch = self.torch
        x = torch.randn(n, self.ndims)
        log_q = self._log_p0(x)
        for t, h in self._grid():
            log_q = log_q - torch.linalg.slogdet(self._step_jacobian(channel, x, t, h))[1]
            x = self._step(channel, x, t, h)
        return x.detach(), log_q.detach()

    def _discrete_log_prob(self, channel, z, newton_iters=20, tol=1e-6):
        # Invert each step by Newton's method (reusing the step Jacobian), last step first.
        torch = self.torch
        log_det = torch.zeros(z.shape[0], dtype=z.dtype)
        y = z.detach()
        for t, h in reversed(self._grid()):
            x = y - h * self._velocity(channel, y, torch.full((y.shape[0],), t, dtype=y.dtype))
            for _ in range(newton_iters):
                jac = self._step_jacobian(channel, x, t, h)
                resid = self._step(channel, x, t, h) - y
                if resid.abs().max() < tol:
                    break
                x = (x - torch.linalg.solve(jac, resid.unsqueeze(-1)).squeeze(-1)).detach()
            else:
                # A step that folds space is not invertible and its density is wrong.
                warnings.warn(f"FlowMatchingFlow: step inversion did not converge (max "
                              f"residual {float(resid.abs().max()):.1e}); reduce step_size")
            log_det = log_det + torch.linalg.slogdet(self._step_jacobian(channel, x, t, h))[1]
            y = x
        return (self._log_p0(y) - log_det).detach()

    def _core_sample(self, channel, n):
        torch = self.torch
        x0 = torch.randn(n, self.ndims)
        return self._solver(channel).sample(
            x_init=x0, step_size=self._step_size, method=self._method
        )

    def _core_log_prob(self, channel, z):
        if self._method in self._DISCRETE:
            return self._discrete_log_prob(channel, z)
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
        n_draw = len(channels)
        if self._buffer_size:
            channels, rans, values = self._buffered(channels, rans, values)
        channels = np.asarray(channels).astype(np.int64)
        rans_t = torch.as_tensor(np.asarray(rans, dtype=np.float64), dtype=torch.float32)
        rans_t = rans_t.clamp(_EPS, 1.0 - _EPS)
        weights = self._normalized_weights(values)
        if weights is None:
            return 0.0

        ch_t = torch.as_tensor(channels, dtype=torch.long)
        groups = [(c, idx) for c in range(self.nchannels)
                  if (idx := np.nonzero(channels == c)[0]).size > 0]
        z1s = {c: self.latent(c, rans_t[idx]) for c, idx in groups}  # latent targets
        w_cs = {c: weights[idx] for c, idx in groups}
        resample = self._ot_batch or self._buffer_size
        pairs = ({c: self._pairs(z1s[c], w_cs[c], n_draw * idx.size // len(channels))
                  for c, idx in groups} if resample else {})

        # The CFM loss is simulation-free, so the (expensive) batch is reused for several
        # steps, each with a fresh t (and a fresh x0 unless the pairs are OT-coupled).
        for _ in range(self._train_steps):
            self.optimizer.zero_grad()
            # Discrete part: importance-weighted categorical NLL on the channel logits.
            loss = -(weights * torch.log_softmax(self.logits, dim=0)[ch_t]).sum()
            # Continuous part: weighted conditional flow matching per channel.
            for c, idx in groups:
                if resample:
                    x0, z1 = pairs[c]
                    w = torch.full((z1.shape[0],), float(w_cs[c].sum()) / z1.shape[0])
                else:
                    z1, w = z1s[c], w_cs[c]
                    x0 = torch.randn_like(z1)
                t = torch.rand(z1.shape[0], 1)
                xt = (1.0 - t) * x0 + t * z1          # CondOT straight-line path
                target = z1 - x0                       # path velocity
                pred = self._velocity(c, xt, t.squeeze(-1))
                loss = loss + (w * ((pred - target) ** 2).mean(dim=1)).sum()
            loss.backward()
            self._optimizer_step()

        # Record for plot_loss (the base train() does this; this override must too).
        self.loss.append(loss.detach().cpu().item())
        return float(loss.detach().cpu().item())

    def _buffered(self, channels, rans, values):
        batch = (np.asarray(channels), np.asarray(rans), np.asarray(values))
        if self._buffer is not None:
            batch = tuple(np.concatenate((old, new))[-self._buffer_size:]
                          for old, new in zip(self._buffer, batch))
        self._buffer = batch
        return batch

    def _pairs(self, z1, w, n):
        # Resample n targets by weight; with ot_batch, pair each chunk with base draws by
        # exact assignment (minibatch OT-CFM): straighter paths, lower-variance targets.
        torch = self.torch
        z1 = z1[torch.multinomial(w / w.sum(), max(1, n), replacement=True)]
        x0 = torch.randn_like(z1)
        if not self._ot_batch:
            return x0, z1
        from scipy.optimize import linear_sum_assignment

        for lo in range(0, z1.shape[0], self._ot_batch):
            sl = slice(lo, lo + self._ot_batch)
            cost = torch.cdist(x0[sl], z1[sl]).pow(2).numpy()
            _, col = linear_sum_assignment(cost)
            z1[sl] = z1[sl][torch.as_tensor(col)]
        return x0, z1


class VegasFlowMatchingFlow(VegasWarmStart, FlowMatchingFlow):
    """VEGAS grid adaptation for the first ``vegas_iterations`` calls, then flow matching."""

    # Tuned on the QE benchmark (validation/sweep_fm); zero_init makes the flow the
    # identity at the VEGAS -> flow switch.
    # step_size 0.25: four steps are as good as twenty and 2x cheaper; a single step can fold
    # space (non-invertible, wrong density).
    _DEFAULTS = {"zero_init": True, "train_steps": 10, "ot_batch": 256, "buffer_size": 50000,
                 "loss_power": 1.5, "lr_decay": 0.999, "step_size": 0.25}

    def __init__(self, ndims, nchannels, **params):
        super().__init__(ndims, nchannels, **{**self._DEFAULTS, **params})

    def _vegas_phase_loss(self, channels, rans_t, weights):
        # The core is frozen in this phase, so skip the ODE likelihood: only the logits step.
        ch_t = self.torch.as_tensor(channels, dtype=self.torch.long)
        return -(weights * self.torch.log_softmax(self.logits, dim=0)[ch_t]).sum()
