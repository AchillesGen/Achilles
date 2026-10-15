"""VEGAS warm start shared by the flow backends.

:class:`VegasWarmStart` adapts a per-channel VEGAS grid (the same algorithm as the C++
``AdaptiveMap``) for ``vegas_iterations`` training calls, then freezes it and trains the
backend's flow on top of it: ``rans = V(c(z))``, where ``c`` is the backend's own
latent-to-cube map. The flow should be the identity at the switch, so training starts
from the VEGAS-adapted density.
"""

import numpy as np

from ._torch_base import _EPS


class VegasWarmStart:
    """Mixin: put in front of a :class:`ChannelFlow` subclass in the bases."""

    def __init__(self, ndims, nchannels, vegas_iterations=10, vegas_bins=50, alpha=1.5,
                 **params):
        super().__init__(ndims, nchannels, **params)
        self.vegas_iterations = int(vegas_iterations)
        self.alpha = float(alpha)
        grid = self.torch.linspace(0.0, 1.0, int(vegas_bins) + 1, dtype=self.torch.float64)
        self.edges = grid.repeat(self.nchannels, self.ndims, 1)  # [channel, dim, bin edge]
        self.iteration = 0

    @property
    def in_vegas_phase(self):
        return self.iteration < self.vegas_iterations

    # -- VEGAS map composed with the backend squash: rans = V(c(z)) ------------------
    def _bin(self, channel, x):
        # Bin index of x[:, dim] in this channel's grid, shape [n, d].
        edges = self.edges[channel]
        idx = self.torch.searchsorted(edges, x.T.contiguous(), right=True).T - 1
        return idx.clamp(0, edges.shape[1] - 2)

    def _left_width(self, channel, idx):
        edges = self.edges[channel].expand(idx.shape[0], -1, -1)
        left = edges.gather(2, idx.unsqueeze(-1)).squeeze(-1)
        right = edges.gather(2, (idx + 1).unsqueeze(-1)).squeeze(-1)
        return left, right - left

    def latent(self, channel, rans_t):
        x = rans_t.to(self.torch.float64).clamp(_EPS, 1.0 - _EPS)
        idx = self._bin(channel, x)
        left, width = self._left_width(channel, idx)
        nbins = self.edges.shape[2] - 1
        u = (idx + (x - left) / width) / nbins
        return super().latent(channel, u.to(rans_t.dtype).clamp(_EPS, 1.0 - _EPS))

    def unit_cube(self, channel, z):
        u = super().unit_cube(channel, z)
        nbins = self.edges.shape[2] - 1
        pos = u.to(self.torch.float64) * nbins
        idx = pos.floor().long().clamp(0, nbins - 1)
        left, width = self._left_width(channel, idx)
        return (left + (pos - idx) * width).to(z.dtype).clamp(_EPS, 1.0 - _EPS)

    def _squash_log_det(self, channel, z):
        u = super().unit_cube(channel, z)
        nbins = self.edges.shape[2] - 1
        idx = (u.to(self.torch.float64) * nbins).floor().long().clamp(0, nbins - 1)
        _, width = self._left_width(channel, idx)
        vegas = self.torch.log(width * nbins).sum(dim=1).to(z.dtype)
        return super()._squash_log_det(channel, z) + vegas

    # -- training --------------------------------------------------------------------
    def _vegas_phase_loss(self, channels, rans_t, weights):
        # Same NLL as the flow phase, so the loss curve is continuous.
        log_q = self.torch.zeros(channels.shape[0])
        for c in range(self.nchannels):
            idx = np.nonzero(channels == c)[0]
            if idx.size:
                log_q[idx] = self._log_density(c, rans_t[idx])
        return -(weights * log_q).sum()

    def train(self, channels, rans, values):
        if not self.in_vegas_phase:
            return super().train(channels, rans, values)

        torch = self.torch
        self.iteration += 1
        channels = np.asarray(channels).astype(np.int64)
        rans_t = torch.as_tensor(np.asarray(rans, dtype=np.float64), dtype=torch.float32)
        weights = self._normalized_weights(values)
        if weights is None:
            return 0.0

        # Only the channel logits step; the grid is adapted by VEGAS, not by gradient.
        self.optimizer.zero_grad()
        loss = self._vegas_phase_loss(channels, rans_t, weights)
        loss.backward()
        for p in self.cores.parameters():
            p.grad = None
        self.optimizer.step()
        self.loss.append(loss.detach().cpu().item())

        self._adapt_grids(channels, rans, np.abs(np.asarray(values, dtype=np.float64)))
        return float(loss.detach().cpu().item())

    def _adapt_grids(self, channels, rans, values):
        rans = np.asarray(rans, dtype=np.float64)
        nbins = self.edges.shape[2] - 1
        for c in range(self.nchannels):
            sel = channels == c
            if not sel.any():
                continue
            idx = self._bin(c, self.torch.as_tensor(rans[sel])).numpy()
            val2 = values[sel] ** 2
            for d in range(self.ndims):
                data = np.bincount(idx[:, d], weights=val2, minlength=nbins)
                edges = self.edges[c, d].numpy()
                self.edges[c, d] = self.torch.as_tensor(_adapt(edges, data, self.alpha))

    # -- persistence -----------------------------------------------------------------
    def _state(self):
        state = super()._state()
        state["vegas_edges"] = self.edges
        state["vegas_iteration"] = self.iteration
        return state

    def _set_state(self, state):
        super()._set_state(state)
        self.edges = state["vegas_edges"].to(self.torch.float64)
        self.iteration = int(state["vegas_iteration"])


def _adapt(edges, data, alpha):
    """One VEGAS refinement of a 1D grid; mirrors ``AdaptiveMap::Adapt``."""
    smooth = np.empty_like(data)
    smooth[0] = (data[0] + data[1]) / 2
    smooth[-1] = (data[-2] + data[-1]) / 2
    smooth[1:-1] = (data[:-2] + data[1:-1] + data[2:]) / 3
    norm = smooth.sum()
    if norm == 0:
        return edges

    r = smooth / norm
    fac = np.zeros_like(r)
    nz = r > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        fac[nz] = np.where(r[nz] < 1, ((r[nz] - 1) / np.log(r[nz])) ** alpha, 1.0)

    # New edges split the cumulative importance into equal parts, linear within old bins.
    cumulative = np.concatenate(([0.0], np.cumsum(fac)))
    targets = cumulative[-1] * np.arange(1, len(data)) / len(data)
    new_edges = np.empty_like(edges)
    new_edges[0], new_edges[-1] = 0.0, 1.0
    new_edges[1:-1] = np.interp(targets, cumulative, edges)
    return new_edges
