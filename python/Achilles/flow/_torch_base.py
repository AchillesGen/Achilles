"""Shared machinery for PyTorch-based flow backends.

All trainable backends model the same space -- a discrete channel plus continuous unit-
hypercube coordinates ``rans`` -- in the same way, so that a Vegas-vs-flow benchmark only
varies the *continuous core*:

* channel choice: a learned categorical (``logits``);
* continuous part: an independent per-channel flow in an unbounded latent space, squashed
  to the unit cube with the standard-normal CDF ``Φ`` (its Jacobian is folded into the
  reported density);
* training: importance-weighted maximum likelihood toward ``p(x) ∝ f(x)`` using the
  integrand weights handed back from C++.

The squash is ``Φ`` (not a sigmoid) on purpose: the cores have a standard-normal base, so
an identity-initialized flow gives ``rans = Φ(z)`` with ``z ~ N(0, 1)``, i.e. *uniform*
``rans``. That reproduces the pure analytic-``Mapper`` sampling, so the flow starts from
(and refines) the analytic mapping rather than wasting training relearning it. A sigmoid
would instead peak ``rans`` at 0.5 and force the flow to unlearn that distortion first.

A subclass only implements the continuous core via :meth:`_build_cores`,
:meth:`_core_sample` and :meth:`_core_log_prob`. Backends whose training objective is not
likelihood-based (e.g. flow matching) override :meth:`train`.
"""

import math

import numpy as np

from . import Flow

_EPS = 1.0e-6


def require_torch():
    try:
        import torch  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "This flow backend requires PyTorch. Install the flow extra with "
            "`pip install achilles[flow]`."
        ) from exc
    return torch


def checkerboard_mask(torch, ndims, block):
    mask = torch.zeros(ndims)
    mask[block % 2 :: 2] = 1
    return mask


class ChannelFlow(Flow):
    """Base for flows over ``(channel, unit-cube rans)`` trained by importance weighting."""

    requires_training = True

    def __init__(self, ndims, nchannels, lr=1.0e-3, channel_select="learned", **params):
        torch = require_torch()
        self.torch = torch
        self.ndims = int(ndims)
        self.nchannels = int(nchannels)
        self.cores = self._build_cores(**params)  # one core per channel, or one conditional core

        # Channel-selection strategy (a benchmark axis):
        #   "learned"  -- a small learned categorical (logits) picks the channel;
        #   "uniform"  -- channels are drawn uniformly and the analytic multichannel
        #                 weighting (done in C++ over the Mappers) handles channel
        #                 importance, i.e. reuse the existing multichannel selection.
        self._learn_channels = str(channel_select).lower() == "learned"
        self.logits = torch.nn.Parameter(
            torch.zeros(self.nchannels), requires_grad=self._learn_channels
        )
        opt_params = list(self.cores.parameters())
        if self._learn_channels:
            opt_params.append(self.logits)
        self.optimizer = torch.optim.Adam(opt_params, lr=float(lr))

    # -- subclass hooks --------------------------------------------------------------
    def _build_cores(self, **params):
        """Return an ``nn.ModuleList`` with one continuous core per channel."""
        raise NotImplementedError

    def _core_sample(self, channel, n):
        """Return ``n`` latent samples ``z`` (tensor ``[n, ndims]``) from the core."""
        raise NotImplementedError

    def _core_log_prob(self, channel, z):
        """Return ``log q(z)`` (tensor ``[n]``) under the core."""
        raise NotImplementedError

    # One-hot channel context, for backends that condition a single shared core on the
    # channel rather than keeping one core per channel.
    def _channel_onehot(self, channel, n):
        torch = self.torch
        ctx = torch.zeros(n, self.nchannels)
        ctx[:, channel] = 1.0
        return ctx

    # -- shared density in unit-cube space -------------------------------------------
    def _log_density(self, channel, rans_t):
        torch = self.torch
        rans_t = rans_t.clamp(_EPS, 1.0 - _EPS)
        z = self.latent(rans_t)  # probit: z = Phi^{-1}(rans), matches the Gaussian base
        log_q = self._core_log_prob(channel, z)
        # Change of variables for rans = Phi(z): log|dz/d rans| = -log phi(z).
        log_phi = torch.sum(-0.5 * z * z - 0.5 * math.log(2.0 * math.pi), dim=1)
        log_p_channel = torch.log_softmax(self.logits, dim=0)[channel]
        return log_p_channel + log_q - log_phi

    # Latent <-> unit-cube transforms (standard-normal CDF and its inverse).
    def latent(self, rans_t):
        return self.torch.special.ndtri(rans_t.clamp(_EPS, 1.0 - _EPS))

    def unit_cube(self, z):
        return self.torch.special.ndtr(z).clamp(_EPS, 1.0 - _EPS)

    # -- protocol --------------------------------------------------------------------
    def sample(self, n):
        torch = self.torch
        n = int(n)
        with torch.no_grad():
            probs = torch.softmax(self.logits, dim=0)
            channels = torch.multinomial(probs, n, replacement=True)
            rans = torch.empty((n, self.ndims))
            density = torch.empty(n)
            for c in range(self.nchannels):
                idx = (channels == c).nonzero(as_tuple=True)[0]
                if idx.numel() == 0:
                    continue
                z = self._core_sample(c, int(idx.numel()))
                r = self.unit_cube(z)
                rans[idx] = r
                density[idx] = torch.exp(self._log_density(c, r))
        return (
            channels.cpu().numpy().astype(np.int64),
            rans.cpu().numpy().astype(np.float64),
            density.cpu().numpy().astype(np.float64),
        )

    def density(self, channels, rans):
        torch = self.torch
        channels = np.asarray(channels).astype(np.int64)
        rans_t = torch.as_tensor(np.asarray(rans, dtype=np.float64), dtype=torch.float32)
        out = torch.empty(channels.shape[0])
        with torch.no_grad():
            for c in range(self.nchannels):
                idx = np.nonzero(channels == c)[0]
                if idx.size == 0:
                    continue
                out[idx] = torch.exp(self._log_density(c, rans_t[idx]))
        return out.cpu().numpy().astype(np.float64)

    def train(self, channels, rans, values):
        # Importance-weighted MLE toward p ∝ f: maximize sum_i w_i log q(x_i).
        torch = self.torch
        channels = np.asarray(channels).astype(np.int64)
        rans_t = torch.as_tensor(np.asarray(rans, dtype=np.float64), dtype=torch.float32)
        weights = self._normalized_weights(values)
        if weights is None:
            return 0.0

        self.optimizer.zero_grad()
        log_q = torch.zeros(channels.shape[0])
        for c in range(self.nchannels):
            idx = np.nonzero(channels == c)[0]
            if idx.size == 0:
                continue
            log_q[idx] = self._log_density(c, rans_t[idx].clamp(_EPS, 1.0 - _EPS))
        loss = -(weights * log_q).sum()
        loss.backward()
        self.optimizer.step()
        return float(loss.detach().cpu().item())

    # -- helpers ---------------------------------------------------------------------
    def _normalized_weights(self, values):
        torch = self.torch
        weights = torch.as_tensor(np.asarray(values, dtype=np.float64), dtype=torch.float32).abs()
        norm = weights.sum()
        if norm <= 0:
            return None
        return (weights / norm).detach()

    def save(self, path):
        self.torch.save({"cores": self.cores.state_dict(), "logits": self.logits.detach()}, path)

    def load(self, path):
        torch = self.torch
        try:
            state = torch.load(path, map_location="cpu")
        except (FileNotFoundError, RuntimeError):
            return False
        self.cores.load_state_dict(state["cores"])
        with torch.no_grad():
            self.logits.copy_(state["logits"])
        return True
