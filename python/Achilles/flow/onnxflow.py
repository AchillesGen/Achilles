"""Inference-only flow backed by ONNX models (via onnxruntime).

This is a reference for using a **pre-trained** flow exported to ONNX: it does not train
(``requires_training = False``), it only samples and evaluates densities. Because ONNX
graphs are model-specific, the input/output tensor names are configurable and the
expected contract is documented below; adapt the names to match your export.

Two ONNX models are used (commonly how a flow is exported):

* **sample model** (``sample_model``): given Gaussian latents ``z`` of shape ``[n, d]``,
  returns ``channels`` (int64 ``[n]``), ``rans`` (float ``[n, d]`` in the unit cube) and
  ``density`` (float ``[n]``), the sampling pdf ``g(channel, rans)``.
* **density model** (``density_model``): given ``channels`` (int64 ``[n]``) and ``rans``
  (float ``[n, d]``), returns ``density`` (float ``[n]``).

Run card example::

    Flow:
      Class: OnnxFlow
      File: ""                       # built-in, so no file needed
      sample_model: /path/sample.onnx
      density_model: /path/density.onnx
"""

import numpy as np


def _require_onnxruntime():
    try:
        import onnxruntime as ort  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "The OnnxFlow backend requires onnxruntime. Install the flow extra with "
            "`pip install achilles[flow]`."
        ) from exc
    return ort


class OnnxFlow:
    requires_training = False

    def __init__(
        self,
        ndims,
        nchannels,
        sample_model,
        density_model,
        latent_input="z",
        sample_outputs=("channels", "rans", "density"),
        density_inputs=("channels", "rans"),
        density_output="density",
        **_unused,
    ):
        ort = _require_onnxruntime()
        self.ndims = int(ndims)
        self.nchannels = int(nchannels)
        self._latent_input = latent_input
        self._sample_outputs = list(sample_outputs)
        self._density_inputs = list(density_inputs)
        self._density_output = density_output
        self._sample = ort.InferenceSession(sample_model)
        self._density = ort.InferenceSession(density_model)

    def sample(self, n):
        n = int(n)
        z = np.random.standard_normal((n, self.ndims)).astype(np.float32)
        channels, rans, density = self._sample.run(
            self._sample_outputs, {self._latent_input: z}
        )
        return (
            np.asarray(channels).reshape(n).astype(np.int64),
            np.asarray(rans).reshape(n, self.ndims).astype(np.float64),
            np.asarray(density).reshape(n).astype(np.float64),
        )

    def density(self, channels, rans):
        channels = np.asarray(channels).astype(np.int64)
        rans = np.asarray(rans, dtype=np.float32)
        feeds = {
            self._density_inputs[0]: channels,
            self._density_inputs[1]: rans,
        }
        (density,) = self._density.run([self._density_output], feeds)
        return np.asarray(density).reshape(channels.shape[0]).astype(np.float64)

    # Inference-only: training is a no-op; nothing to persist.
    def train(self, channels, rans, values):
        return 0.0

    def save(self, path):
        return None

    def load(self, path):
        return True
