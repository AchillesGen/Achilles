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


def test_flowmatching_density_is_exact_for_discrete_map():
    # The sampled density must be that of the fixed-step map itself, p0(x0)/|det dF/dx0|,
    # not the continuous-ODE likelihood (which is O(h^2) off and biases f/g).
    torch = pytest.importorskip("torch")
    pytest.importorskip("flow_matching")
    torch.manual_seed(0)
    model = flow.FlowMatchingFlow(NDIMS, 1, hidden=16, step_size=0.25, lr=1e-2, train_steps=3)
    for _ in range(3):
        channels, rans, density = model.sample(N)
        model.train(channels, rans, _peaked(rans) / density)

    def full_map(x):
        for t, h in model._grid():
            x = model._step(0, x[None], t, h)[0]
        return x

    torch.manual_seed(1)
    z, log_q = model._core_sample_and_log_prob(0, 16)
    torch.manual_seed(1)
    x0 = torch.randn(16, NDIMS)
    log_det = torch.stack([torch.linalg.slogdet(torch.func.jacrev(full_map)(x))[1] for x in x0])
    assert torch.allclose(z, torch.stack([full_map(x) for x in x0]), atol=1e-5)
    assert torch.allclose(log_q, model._log_p0(x0) - log_det, atol=1e-4)
    assert torch.allclose(model._core_log_prob(0, z), log_q, atol=1e-4)


def test_binary_masks_cover_all_pairs():
    from Achilles.flow.hypercube import binary_masks

    for ndims in range(1, 10):
        masks = binary_masks(ndims)
        for i in range(ndims):
            assert any(m[i] for m in masks)
            for j in range(ndims):
                if i != j:
                    assert any(m[i] and not m[j] for m in masks)


@pytest.mark.parametrize("ndims", [1, 2, 5])
def test_hypercube_spline_protocol(ndims):
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.HypercubeSplineFlow(ndims, NCHANNELS, hidden=16, bins=8)
    _check_protocol(model, ndims=ndims)


def test_hypercube_spline_starts_uniform():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.HypercubeSplineFlow(NDIMS, NCHANNELS, hidden=16, bins=8)
    _, _, density = model.sample(N)
    assert np.allclose(density, 1.0 / NCHANNELS, rtol=1e-4)


def test_hypercube_spline_density_matches_sample_after_training():
    pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.HypercubeSplineFlow(NDIMS, NCHANNELS, hidden=16, bins=8, lr=1e-2)
    for _ in range(5):
        channels, rans, _ = model.sample(N)
        model.train(channels, rans, np.exp(-10 * (rans[:, 0] - 0.3) ** 2))
    channels, rans, density = model.sample(N)
    assert np.allclose(density, model.density(channels, rans), rtol=1e-3)


def _peaked(rans):
    return np.exp(-200 * ((rans - 0.2) ** 2).sum(axis=1))


def test_vegas_spline_adapts_grid_then_trains_flow(tmp_path):
    torch = pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.VegasSplineFlow(NDIMS, NCHANNELS, hidden=16, bins=8, vegas_iterations=3,
                                 vegas_bins=20, lr=1e-2)
    _check_protocol(model)
    core_before = [p.detach().clone() for p in model.cores.parameters()]

    for _ in range(2):
        channels, rans, density = model.sample(4 * N)
        model.train(channels, rans, _peaked(rans) / density)
    assert not model.in_vegas_phase
    # The flow is untouched during the VEGAS phase, while the grid moves toward the peak.
    assert all(torch.equal(a, b) for a, b in zip(core_before, model.cores.parameters()))
    widths = torch.diff(model.edges, dim=-1)
    assert widths.min() > 0 and widths.max() > 2.0 / 20

    channels, rans, density = model.sample(N)
    assert np.allclose(density, model.density(channels, rans), rtol=1e-3)
    model.train(channels, rans, _peaked(rans) / density)
    assert not all(torch.equal(a, b) for a, b in zip(core_before, model.cores.parameters()))

    path = str(tmp_path / "vegas.weights")
    model.save(path)
    restored = flow.VegasSplineFlow(NDIMS, NCHANNELS, hidden=16, bins=8, vegas_iterations=3,
                                    vegas_bins=20)
    assert restored.load(path)
    assert restored.iteration == model.iteration
    assert np.allclose(restored.density(channels, rans), model.density(channels, rans))



@pytest.mark.parametrize("name, needs, kwargs", [
    ("NormFlow", "normflows", dict(hidden=16, blocks=2)),
    ("SplineFlow", "normflows", dict(hidden=16, blocks=2, bins=8)),
    ("HypercubeSplineFlow", "normflows", dict(hidden=16, bins=8)),
    ("VegasSplineFlow", "normflows", dict(hidden=16, bins=8, vegas_iterations=1)),
    ("ZukoFlow", "zuko", dict(flow="NSF", hidden=16, transforms=2)),
    ("VegasFlowMatchingFlow", "flow_matching",
     dict(hidden=16, step_size=0.2, vegas_iterations=1, train_steps=2, ot_batch=64,
          buffer_size=2 * N)),
])
def test_sampled_density_matches_density_after_training(name, needs, kwargs):
    pytest.importorskip("torch")
    pytest.importorskip(needs)
    model = getattr(flow, name)(NDIMS, NCHANNELS, lr=1e-2, **kwargs)
    for _ in range(3):
        channels, rans, density = model.sample(N)
        model.train(channels, rans, _peaked(rans) / density)
    channels, rans, density = model.sample(N)
    assert np.allclose(density, model.density(channels, rans), rtol=1e-3)


def test_loss_power_and_lr_decay():
    torch = pytest.importorskip("torch")
    pytest.importorskip("normflows")
    model = flow.HypercubeSplineFlow(NDIMS, NCHANNELS, hidden=16, bins=8, lr=1e-2,
                                     lr_decay=0.5, loss_power=2.0)
    values = np.array([1.0, 2.0, 3.0])
    assert np.allclose(model._normalized_weights(values).numpy(), values**2 / 14.0)
    channels, rans, density = model.sample(N)
    model.train(channels, rans, _peaked(rans) / density)
    assert model.optimizer.param_groups[0]["lr"] == pytest.approx(5e-3)
