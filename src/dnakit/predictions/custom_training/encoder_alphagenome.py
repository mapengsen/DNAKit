"""Differentiate AlphaGenome's JAX backbone through the PyTorch trainer.

The PyTorch parameters are the optimizer's source of truth. DLPack transfers
values between runtimes, and an explicit JAX VJP supplies their gradients;
there is no NumPy conversion on the differentiable path. Haiku state and the
PRNG key are saved with the weights so evaluation can be restored independently
of the original checkpoint. The outer trainer owns progress reporting.
"""

from __future__ import annotations

import json
import math
import os
import warnings
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.representations import DNAEmbeddingModel

from .encoder import _require

_MANIFEST = "alphagenome_encoder.json"
_WEIGHTS = "alphagenome_weights"


def _build_training_forward(*, force_float32: bool) -> tuple[Any, int]:
    """Build a stateful trunk using the official Haiku module and module names."""
    hk = _require("haiku")
    jmp = _require("jmp")
    model = _require("alphagenome_research.model.model")
    dna_model = _require("alphagenome_research.model.dna_model")
    api = _require("alphagenome.models.dna_model")
    modern_api = hasattr(model.AlphaGenome, "forward_trunk")
    if not modern_api:
        warnings.warn(
            "This AlphaGenome source only exposes inference normalization statistics. "
            "Full fine-tuning updates all backbone weights, including normalization "
            "scale/offset, while keeping running statistics at checkpoint values. "
            "A newer alphagenome_research with forward_trunk also updates these statistics.",
            UserWarning,
            stacklevel=2,
        )
    compute = "float32" if force_float32 else "bfloat16"
    policy = jmp.get_policy(f"params=float32,compute={compute},output={compute}")

    @hk.transform_with_state
    def forward(dna: Any, organism: Any, is_training: bool) -> Any:
        with hk.mixed_precision.push_policy(model.AlphaGenome, policy):
            # Empty output metadata skips unrelated pretrained biological heads;
            # all trunk module names stay identical to the official checkpoint.
            backbone = model.AlphaGenome({})
            if modern_api:
                return backbone.forward_trunk(dna, organism, is_training=is_training).embeddings_1bp
            _, embeddings = backbone(dna, organism)
            return embeddings.embeddings_1bp

    human = api.Organism.HOMO_SAPIENS
    converter = getattr(dna_model, "convert_to_organism_index", None)
    organism_index = int(converter(human)) if converter is not None else int(human.value)
    forward.apply._dnakit_state_policy = (
        "updated_during_training" if modern_api else "fixed_checkpoint_statistics"
    )
    return forward.apply, organism_index


def _split_head_parameters(params: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """Keep unused original task heads out of AdamW and its weight decay."""
    backbone, heads = {}, {}
    for name, values in params.items():
        target = heads if name.startswith("alphagenome/head/") else backbone
        target[name] = values
    return backbone, heads


def _autograd_bridge(torch: Any) -> Any:
    """Create the optional Torch subclass only after torch has been imported."""

    class _JaxBackbone(torch.autograd.Function):
        @staticmethod
        def forward(ctx: Any, encoder: Any, dna: Any, *parameters: Any) -> Any:
            # Keeping the tensors also makes Torch detect optimizer updates that
            # incorrectly happen before their matching backward pass.
            ctx.save_for_backward(*parameters)
            ctx.encoder = encoder
            key = encoder._next_key()
            values = tuple(encoder._to_jax(parameter) for parameter in parameters)
            state = encoder.state

            def apply(*leaves: Any) -> Any:
                return encoder._forward(leaves, state, dna, key, encoder.training)

            output, ctx.pullback, updated_state = encoder.jax.vjp(apply, *values, has_aux=True)
            if encoder.training:
                encoder.state = updated_state
            return encoder._to_torch(output)

        @staticmethod
        def backward(ctx: Any, gradient: Any) -> tuple[Any, ...]:
            encoder = ctx.encoder
            parameters = ctx.saved_tensors
            gradients = ctx.pullback(encoder._to_jax(gradient))
            result = tuple(
                encoder._to_torch(value).to(dtype=parameter.dtype)
                for value, parameter in zip(gradients, parameters, strict=True)
            )
            return (None, None, *result)

    return _JaxBackbone


class AlphaGenomeSequenceEncoder:
    """Train all floating-point backbone weights with JAX reverse-mode autodiff.

    The pretrained biological output heads are retained in the checkpoint, but
    the user-defined Torch task head consumes the pooled 1bp trunk embeddings.
    CPU uses float32; GPU ``auto`` follows the official bfloat16 compute policy
    while keeping master parameters in float32.
    """

    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: str | Path,
        *,
        device: str,
        dtype: str,
        pooling: str,
        max_length: int | None,
        allow_remote_code: bool,
        model_source_path: str | Path | None = None,
        seed: int = 42,
    ) -> None:
        del allow_remote_code  # The official installed package supplies the implementation.
        from dnakit.representations.backends import (
            _cast_jax_float32,
            _prepend_source_path,
            _resolve_torch_device,
        )

        if spec.loader != "alphagenome":
            raise ConfigurationError(
                "The AlphaGenome encoder requires an AlphaGenome checkpoint.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        if dtype not in {"auto", "float32"}:
            raise ConfigurationError(
                "AlphaGenome supports dtype='auto' or dtype='float32'.",
                code="INVALID_REPRESENTATION_DTYPE",
            )
        if pooling not in {"mean", "max", "cls", "last"}:
            raise ConfigurationError(
                "Unknown AlphaGenome pooling method.",
                code="INVALID_REPRESENTATION_POOLING",
            )
        _prepend_source_path(None if model_source_path is None else str(model_source_path))
        self.spec = spec
        self.torch = _require("torch")
        # JAX and Torch must share GPU memory. Respect an explicit user setting,
        # but avoid JAX reserving most of the GPU before AdamW allocates its state.
        os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
        self.jax = _require("jax")
        self.jnp = _require("jax.numpy")
        self.ocp = _require("orbax.checkpoint")
        self.device = _resolve_torch_device(self.torch, device)
        if self.device.type not in {"cpu", "cuda"}:
            raise ConfigurationError(
                "The AlphaGenome JAX/Torch bridge requires CPU or CUDA.",
                code="INVALID_REPRESENTATION_DEVICE",
            )
        platform = "gpu" if self.device.type == "cuda" else "cpu"
        try:
            devices = self.jax.devices(platform)
            index = self.device.index
            if index is None:
                index = self.torch.cuda.current_device() if platform == "gpu" else 0
            self.jax_device = devices[index]
        except (RuntimeError, IndexError) as exc:
            raise ConfigurationError(
                "The requested device must be available in both JAX and PyTorch.",
                code="INVALID_REPRESENTATION_DEVICE",
                context={"device": str(self.device)},
            ) from exc

        self.checkpoint = Path(checkpoint).expanduser().resolve()
        manifest_path = self.checkpoint / _MANIFEST
        saved: dict[str, Any] = {}
        if manifest_path.is_file():
            saved = json.loads(manifest_path.read_text(encoding="utf-8"))
            if saved.get("format_version") != 1 or saved.get("model") != spec.name:
                raise BackendExecutionError(
                    "Invalid AlphaGenome training checkpoint metadata.", code="MODEL_LOAD_FAILED"
                )
        self.pooling = str(saved.get("pooling", pooling))
        self.chunk_length = int(saved.get("chunk_length", max_length or spec.chunk_length))
        if not 1 <= self.chunk_length <= spec.chunk_length:
            raise ConfigurationError(
                "AlphaGenome max_length must be within its supported sequence length.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        self.seed = int(saved.get("seed", seed))
        with self.jax.default_device(self.jax_device):
            self._rng = (
                self.jnp.asarray(saved["rng"], dtype=self.jnp.uint32)
                if "rng" in saved
                else self.jax.random.PRNGKey(self.seed)
            )
        weights = self.checkpoint / _WEIGHTS if saved else self.checkpoint
        with self.ocp.StandardCheckpointer() as checkpointer:
            params, state = checkpointer.restore(str(weights), target=None, strict=False)
        params, heads = _split_head_parameters(params)
        self._head_params = self.jax.device_get(heads)
        params = _cast_jax_float32(params, jax=self.jax, jnp=self.jnp)
        if platform == "cpu":
            state = _cast_jax_float32(state, jax=self.jax, jnp=self.jnp)
        self.state = self.jax.device_put(state, self.jax_device)
        leaves, self._tree = self.jax.tree_util.tree_flatten(params)
        if not leaves or any(
            not self.jnp.issubdtype(leaf.dtype, self.jnp.floating) for leaf in leaves
        ):
            raise BackendExecutionError(
                "AlphaGenome parameters must be a nonempty floating-point tree.",
                code="MODEL_LOAD_FAILED",
            )
        self.model = self.torch.nn.ParameterList(
            [
                self.torch.nn.Parameter(
                    self._to_torch(self.jax.device_put(leaf, self.jax_device)).clone()
                )
                for leaf in leaves
            ]
        )
        forward_apply, organism = _build_training_forward(
            force_float32=platform == "cpu" or dtype == "float32"
        )
        self.normalization_state_policy = getattr(
            forward_apply, "_dnakit_state_policy", "updated_during_training"
        )
        with self.jax.default_device(self.jax_device):
            self.organism = self.jnp.asarray([organism], dtype=self.jnp.int32)

        def forward(flat_params: Any, state: Any, dna: Any, rng: Any, is_training: bool) -> Any:
            current_params = self.jax.tree_util.tree_unflatten(self._tree, flat_params)
            embeddings, updated = forward_apply(
                current_params, state, rng, dna, self.organism, is_training
            )
            values = embeddings.astype(self.jnp.float32)
            if self.pooling == "mean":
                pooled = values.mean(axis=1)
            elif self.pooling == "max":
                pooled = values.max(axis=1)
            elif self.pooling == "cls":
                pooled = values[:, 0, :]
            else:
                pooled = values[:, -1, :]
            return pooled, updated

        self._forward = self.jax.jit(forward, static_argnums=(4,))
        self._bridge = _autograd_bridge(self.torch)
        self.training = False

    def _to_jax(self, tensor: Any) -> Any:
        return self.jax.dlpack.from_dlpack(tensor.detach().contiguous(), device=self.jax_device)

    def _to_torch(self, array: Any) -> Any:
        return self.torch.utils.dlpack.from_dlpack(array).to(self.device)

    def _next_key(self) -> Any:
        if self.training:
            self._rng, key = self.jax.random.split(self._rng)
            return key
        return self._rng

    def parameters(self) -> Any:
        return self.model.parameters()

    def train(self) -> None:
        self.training = True
        self.model.train()

    def eval(self) -> None:
        self.training = False
        self.model.eval()

    def encode(self, sequences: list[str]) -> Any:
        from dnakit.representations.backends import _AlphaGenomeBackend

        if not sequences or any(not sequence for sequence in sequences):
            raise ConfigurationError(
                "AlphaGenome requires nonempty DNA sequences.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        rows = []
        parameters = tuple(self.parameters())
        for sequence in sequences:
            # Dataset ingestion permits IUPAC ambiguity codes. Match frozen
            # representation extraction by masking them before one-hot encoding.
            sequence = "".join(
                base if base in "ACGTN" else "N" for base in sequence.upper().replace("U", "T")
            )
            chunks = []
            for start in range(0, len(sequence), self.chunk_length):
                chunk = sequence[start : start + self.chunk_length]
                size = max(2048, math.ceil(len(chunk) / 2048) * 2048)
                padded = chunk + "N" * (size - len(chunk))
                dna = self.jax.device_put(_AlphaGenomeBackend._one_hot(padded), self.jax_device)
                if self.torch.is_grad_enabled():
                    vector = self._bridge.apply(self, dna, *parameters)
                else:
                    vector_jax, updated = self._forward(
                        tuple(self._to_jax(parameter) for parameter in parameters),
                        self.state,
                        dna,
                        self._next_key(),
                        self.training,
                    )
                    if self.training:
                        self.state = updated
                    vector = self._to_torch(vector_jax)
                chunks.append(vector.squeeze(0))
            rows.append(self.torch.stack(chunks).mean(dim=0))
        return self.torch.stack(rows)

    def save_pretrained(self, directory: Path) -> None:
        directory = Path(directory).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        params = self.jax.tree_util.tree_unflatten(
            self._tree, [self._to_jax(parameter) for parameter in self.parameters()]
        )
        params = {**self._head_params, **params}
        with self.ocp.StandardCheckpointer() as checkpointer:
            checkpointer.save(str(directory / _WEIGHTS), (params, self.state), force=True)
            checkpointer.wait_until_finished()
        manifest = {
            "format_version": 1,
            "model": self.spec.name,
            "pooling": self.pooling,
            "chunk_length": self.chunk_length,
            "seed": self.seed,
            "rng": self.jax.device_get(self._rng).tolist(),
            "normalization_state_policy": self.normalization_state_policy,
        }
        (directory / _MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    def load_pretrained(self, directory: Path) -> None:
        """Restore the best checkpoint without allocating a second backbone."""
        directory = Path(directory).expanduser().resolve()
        saved = json.loads((directory / _MANIFEST).read_text(encoding="utf-8"))
        if (
            saved.get("format_version") != 1
            or saved.get("model") != self.spec.name
            or saved.get("pooling") != self.pooling
            or saved.get("chunk_length") != self.chunk_length
        ):
            raise BackendExecutionError(
                "The AlphaGenome checkpoint configuration does not match this encoder.",
                code="MODEL_LOAD_FAILED",
            )
        with self.ocp.StandardCheckpointer() as checkpointer:
            params, state = checkpointer.restore(
                str(directory / _WEIGHTS), target=None, strict=False
            )
        params, heads = _split_head_parameters(params)
        leaves, tree = self.jax.tree_util.tree_flatten(params)
        if tree != self._tree or any(
            tuple(leaf.shape) != tuple(parameter.shape)
            for leaf, parameter in zip(leaves, self.parameters(), strict=True)
        ):
            raise BackendExecutionError(
                "The AlphaGenome checkpoint parameter tree does not match this encoder.",
                code="MODEL_LOAD_FAILED",
            )
        with self.torch.no_grad():
            for parameter, leaf in zip(self.parameters(), leaves, strict=True):
                parameter.copy_(self._to_torch(self.jax.device_put(leaf, self.jax_device)))
        self.state = self.jax.device_put(state, self.jax_device)
        self._head_params = self.jax.device_get(heads)
        self.seed = int(saved["seed"])
        with self.jax.default_device(self.jax_device):
            self._rng = self.jnp.asarray(saved["rng"], dtype=self.jnp.uint32)


__all__ = ["AlphaGenomeSequenceEncoder"]
