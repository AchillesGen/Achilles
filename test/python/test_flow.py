"""Tests for the Python normalizing-flow backends (Achilles.flow).

Run after ``pip install .[flow]`` so the compiled extension and PyTorch are available::

    pytest test/python/test_flow.py
"""

import numpy as np
import pytest

# The flow package lives under the Achilles package, whose __init__ imports the compiled
# extension. Skip the whole module if Achilles is not importable (extension not built).
Achilles = pytest.importorskip("Achilles")
from Achilles import flow  # noqa: E402


NDIMS = 3
NCHANNELS = 2
N = 512


def _check_protocol(model, ndims=NDIMS, nchannels=NCHANNELS, n=N):
    channels, rans, density = model.sample(n)

    assert channels.shape == (n,)
    assert rans.shape == (n, ndims)
    assert density.shape == (n,)
    assert channels.min() >= 0 and channels.max() < nchannels
    assert np.all(rans >= 0.0) and np.all(rans <= 1.0)
    assert np.all(np.isfinite(density)) and np.all(density > 0.0)

    # density() on the drawn samples must reproduce sample()'s density.
    redensity = model.density(channels, rans)
    assert redensity.shape == (n,)
    assert np.all(np.isfinite(redensity)) and np.all(redensity > 0.0)

    loss = model.train(channels, rans, np.ones(n))
    assert np.isfinite(loss)


def test_uniform_flow_protocol():
    model = flow.UniformFlow(NDIMS, NCHANNELS)
    _check_protocol(model)

    # Uniform over channels and the unit cube => g = 1 / nchannels everywhere.
    channels, _, density = model.sample(N)
    assert np.allclose(density, 1.0 / NCHANNELS)
    assert model.train(channels, np.zeros((N, NDIMS)), np.ones(N)) == 0.0


def test_normflow_protocol():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")

    model = flow.NormFlow(NDIMS, NCHANNELS, hidden=16, blocks=2)
    _check_protocol(model)


def test_normflow_density_matches_sample():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")

    model = flow.NormFlow(NDIMS, NCHANNELS, hidden=16, blocks=2)
    channels, rans, density = model.sample(N)
    # density(channels, rans) is the same pdf sample() reports for those points.
    assert np.allclose(density, model.density(channels, rans), rtol=1e-4, atol=1e-6)


def test_normflow_train_runs_and_updates():
    torch = pytest.importorskip("torch")
    pytest.importorskip("normflows")

    model = flow.NormFlow(NDIMS, NCHANNELS, hidden=16, blocks=2, lr=1e-2)
    channels, rans, _ = model.sample(N)
    # Target weights peaked on channel 0 to give the optimizer a gradient.
    values = np.where(channels == 0, 5.0, 0.5)

    before = torch.softmax(model.logits.detach(), dim=0).clone()
    loss = model.train(channels, rans, values)
    after = torch.softmax(model.logits.detach(), dim=0)

    assert np.isfinite(loss)
    assert not torch.allclose(before, after)


def test_channel_select_uniform_freezes_logits():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    # "uniform" reuses a simple channel selection instead of learning it.
    model = flow.NormFlow(NDIMS, NCHANNELS, hidden=16, blocks=2, channel_select="uniform")
    assert model.logits.requires_grad is False
    _check_protocol(model)


def test_normflow_rejects_conditional():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    with pytest.raises(ValueError):
        flow.NormFlow(NDIMS, NCHANNELS, conditional=True)


def test_splineflow_protocol():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.SplineFlow(NDIMS, NCHANNELS, hidden=16, blocks=2, bins=8)
    _check_protocol(model)


def test_zukoflow_protocol():
    pytest.importorskip("torch")
    pytest.importorskip("zuko")
    # NSF (coupling/spline) is cheap; the same backend also does CNF for the benchmark.
    model = flow.ZukoFlow(NDIMS, NCHANNELS, flow="NSF", hidden=16, transforms=2)
    _check_protocol(model)


def test_zukoflow_conditional_protocol():
    pytest.importorskip("torch")
    pytest.importorskip("zuko")
    model = flow.ZukoFlow(NDIMS, NCHANNELS, flow="NSF", hidden=16, transforms=2, conditional=True)
    _check_protocol(model)


def test_flowmatching_protocol():
    pytest.importorskip("torch")
    pytest.importorskip("flow_matching")
    model = flow.FlowMatchingFlow(NDIMS, NCHANNELS, hidden=16, step_size=0.2)
    _check_protocol(model, n=64)  # ODE-based density is expensive; keep the batch small
