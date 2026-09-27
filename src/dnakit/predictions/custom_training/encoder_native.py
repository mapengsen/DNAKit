"""Differentiable JanusDNA and Evo 2 backbones for supervised fine-tuning.

Native encoders truncate to ``max_length`` bases. Variable-length sequences are
run separately because bidirectional state-space layers do not consistently
honour padding masks. The outer trainer supplies batch/epoch progress bars.
"""

from __future__ import annotations

import json
import pkgutil
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.representations import DNAEmbeddingModel
from dnakit.representations.backends import (
    _clean_janusdna_config,
    _load_janusdna_checkpoint,
    _pool_hidden_states,
    _prepend_source_path,
    _resolve_hidden_states,
    _resolve_torch_device,
    _torch_dtype,
)

from .encoder import _require


def _load_weights(torch: Any, path: Path) -> Any:
    """Load weights without permitting arbitrary checkpoint-provided globals."""
    return torch.load(path, map_location="cpu", weights_only=True)


def _restore_native_state(model: Any, state: dict[str, Any]) -> None:
    """Preserve saved floating precision before load_state_dict copies tensors."""
    for name, tensor in (*model.named_parameters(), *model.named_buffers()):
        saved = state.get(name)
        if saved is not None and hasattr(saved, "dtype") and tensor.dtype != saved.dtype:
            tensor.data = tensor.data.to(dtype=saved.dtype)
    defaults = model.state_dict()
    for name in defaults:
        if name.endswith("._extra_state"):
            state[name] = defaults[name]
    model.load_state_dict(state)


def _require_evo2_cuda(device: Any) -> None:
    if device.type != "cuda":
        raise ConfigurationError(
            "The official Evo 2 / Vortex training implementation requires CUDA.",
            code="INVALID_CUSTOM_PROPERTY_CONFIG",
            hint="Use device='cuda' and CUDA_VISIBLE_DEVICES to select the training GPU.",
        )


class _NativeSequenceEncoder:
    def _configure(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: str | Path,
        *,
        device: str,
        dtype: str,
        pooling: str,
        max_length: int | None,
        model_source_path: str | Path | None,
    ) -> None:
        self.spec = spec
        self.checkpoint = Path(checkpoint).expanduser().resolve()
        self.pooling = pooling
        self.max_length = int(max_length or min(spec.chunk_length, 4_096))
        self.torch = _require("torch")
        self.device = _resolve_torch_device(self.torch, device)
        if self.device.type == "cpu" and dtype in {"float16", "bfloat16"}:
            raise ConfigurationError(
                "Half-precision full fine-tuning is not supported on the CPU path.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
                hint="Use dtype='float32' or dtype='auto' on CPU.",
            )
        self.requested_dtype = _torch_dtype(self.torch, dtype)
        _prepend_source_path(None if model_source_path is None else str(model_source_path))

    def _prepare_model(self) -> None:
        # Vortex's official loader creates inference tensors while permuting
        # projections. Re-enable autograd by materializing ordinary tensors,
        # retaining shared parameter identities when the backbone ties weights.
        replacements: dict[int, Any] = {}
        for module in self.model.modules():
            for name, parameter in tuple(module.named_parameters(recurse=False)):
                if self.torch.is_inference(parameter):
                    replacement = replacements.get(id(parameter))
                    if replacement is None:
                        replacement = self.torch.nn.Parameter(parameter.detach().clone())
                        replacements[id(parameter)] = replacement
                    setattr(module, name, replacement)
            for name, buffer in tuple(module.named_buffers(recurse=False)):
                if self.torch.is_inference(buffer):
                    setattr(module, name, buffer.clone())
        if self.requested_dtype is not None:
            self.model = self.model.to(dtype=self.requested_dtype)
        self.model = self.model.to(self.device)
        for parameter in self.model.parameters():
            parameter.requires_grad_(True)

    def parameters(self) -> Any:
        return self.model.parameters()

    def train(self) -> None:
        self.model.train()

    def eval(self) -> None:
        self.model.eval()

    def _sequences(self, sequences: list[str]) -> list[str]:
        if not sequences or any(not sequence for sequence in sequences):
            raise ConfigurationError(
                "Native training encoders require non-empty DNA sequences.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
            )
        return [
            "".join(base if base in "ACGTN" else "N" for base in sequence[: self.max_length])
            for sequence in sequences
        ]


class JanusDNASequenceEncoder(_NativeSequenceEncoder):
    """Official JanusDNA backbone with autograd and portable trained weights."""

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
    ) -> None:
        del allow_remote_code  # The installed/source-path package supplies the implementation.
        self._configure(
            spec,
            checkpoint,
            device=device,
            dtype=dtype,
            pooling=pooling,
            max_length=max_length,
            model_source_path=model_source_path,
        )
        configuration = _require("janusdna.configuration_janusdna")
        modeling = _require("janusdna.modeling_janusdna")
        payload = json.loads((self.checkpoint / "model_config.json").read_text(encoding="utf-8"))
        self.model_config = configuration.JanusDNAConfig(
            **_clean_janusdna_config(dict(payload["config"]))
        )
        self.model_config._attn_implementation = "eager"
        self.model_config.final_attention_class = "eager"
        if payload.get("dnakit_encoder") == "janusdna":
            self.model = modeling.JanusDNAModel(self.model_config)
            _restore_native_state(
                self.model, _load_weights(self.torch, self.checkpoint / "weights.pt")
            )
        else:
            causal_model = modeling.JanusDNAForCausalLM(self.model_config)
            checkpoint_state = _load_janusdna_checkpoint(
                self.torch, self.checkpoint / "checkpoints" / "last.ckpt", spec
            )
            state = dict(checkpoint_state["state_dict"])
            self.torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(state, "model.")
            for key in tuple(state):
                if ".self_attn.o_projs.0." in key:
                    state[key.replace(".self_attn.o_projs.0.", ".self_attn.o_proj.")] = state.pop(
                        key
                    )
                elif "torchmetrics" in key:
                    state.pop(key)
            missing, unexpected = causal_model.load_state_dict(state, strict=False)
            missing = [key for key in missing if not key.startswith("lm_head.")]
            if missing or unexpected:
                raise BackendExecutionError(
                    "JanusDNA checkpoint did not match the official model architecture.",
                    code="MODEL_LOAD_FAILED",
                    context={"missing": tuple(missing[:10]), "unexpected": tuple(unexpected[:10])},
                )
            self.model = causal_model.model
        self._prepare_model()
        self.base_ids = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}

    def encode(self, sequences: list[str]) -> Any:
        vectors = []
        for sequence in self._sequences(sequences):
            ids = [0, *(self.base_ids[base] for base in sequence), 1]
            input_ids = self.torch.tensor([ids], dtype=self.torch.long, device=self.device)
            # No padding: reversed Mamba streams otherwise ingest other samples' padding.
            outputs = self.model(input_ids=input_ids, use_cache=False, return_dict=True)
            hidden = _resolve_hidden_states(outputs, self.torch)
            special = self.torch.zeros_like(input_ids)
            special[:, 0] = special[:, -1] = 1
            vectors.append(
                _pool_hidden_states(
                    hidden,
                    torch=self.torch,
                    pooling=self.pooling,
                    special_tokens_mask=special,
                )
            )
        return self.torch.cat(vectors, dim=0)

    def save_pretrained(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "model_config.json").write_text(
            json.dumps(
                {"dnakit_encoder": "janusdna", "config": self.model_config.to_dict()}, indent=2
            ),
            encoding="utf-8",
        )
        self.torch.save(self.model.state_dict(), directory / "weights.pt")

    def load_pretrained(self, directory: Path) -> None:
        _restore_native_state(self.model, _load_weights(self.torch, directory / "weights.pt"))


class Evo2SequenceEncoder(_NativeSequenceEncoder):
    """Train all Evo 2 backbone blocks using their final normalized states.

    Frozen representation extraction uses the registry's intermediate layer.
    Full fine-tuning deliberately uses the final norm so every backbone block
    participates in the loss. The inference wrapper's no_grad/detach path and
    language-model unembedding are never called.
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
    ) -> None:
        del allow_remote_code
        self._configure(
            spec,
            checkpoint,
            device=device,
            dtype=dtype,
            pooling=pooling,
            max_length=max_length,
            model_source_path=model_source_path,
        )
        _require_evo2_cuda(self.device)
        modeling = _require("vortex.model.model")
        utilities = _require("vortex.model.utils")
        tokenizer_module = _require("vortex.model.tokenizer")
        saved_config = self.checkpoint / "dnakit_evo2_config.json"
        if saved_config.is_file():
            config = json.loads(saved_config.read_text(encoding="utf-8"))
        else:
            evo_utilities = _require("evo2.utils")
            yaml = _require("yaml")
            raw = pkgutil.get_data("evo2", evo_utilities.CONFIG_MAP["evo2_7b"])
            config = yaml.safe_load(raw)
        # Disable inference-only kernels/FP8 and keep the ordinary autograd path.
        config.update(
            {
                "inference_mode": False,
                "use_fp8_input_projections": False,
                "use_flashfft": False,
                "use_flash_depthwise": False,
                "use_flash_rmsnorm": False,
                "use_hcs_kernel": False,
                "use_hcm_kernel": False,
                "use_hcl_kernel": False,
            }
        )
        if self.device.type != "cuda" or dtype == "float32":
            config["use_flash_attn"] = False
        self.model_config = config
        self.model = modeling.StripedHyena(utilities.dotdict(config))
        if saved_config.is_file():
            # Native state is already permuted. The official loader would permute
            # attention projections a second time and silently corrupt predictions.
            state = _load_weights(self.torch, self.checkpoint / "weights.pt")
            _restore_native_state(self.model, state)
        else:
            utilities.load_checkpoint(self.model, str(self.checkpoint / "evo2_7b.pt"))
        self._prepare_model()
        # Vortex stores a device map separately from nn.Module parameters.
        self.model.block_idx_to_device = {
            index: str(self.device) for index in range(len(self.model.blocks))
        }
        self.tokenizer = tokenizer_module.CharLevelTokenizer(int(config.get("vocab_size", 512)))
        self._clear_inference_state()

    def _clear_inference_state(self) -> None:
        for module in self.model.modules():
            if hasattr(module, "inference_mode"):
                module.inference_mode = False
            # Long Hyena filters must be recomputed from trainable poles/residues.
            # Short/medium filters store h as a Parameter and must remain intact.
            if hasattr(module, "log_poles") and hasattr(module, "residues"):
                if not isinstance(getattr(module, "h", None), self.torch.nn.Parameter):
                    module.h = None
                if hasattr(module, "t"):
                    module.t = None

    def train(self) -> None:
        self._clear_inference_state()
        super().train()

    def encode(self, sequences: list[str]) -> Any:
        vectors = []
        for sequence in self._sequences(sequences):
            input_ids = self.torch.tensor(
                [self.tokenizer.tokenize(sequence)],
                dtype=self.torch.long,
                device=self.device,
            )
            hidden = self.model.embedding_layer(input_ids)
            hidden, _ = self.model.stateless_forward(hidden)
            hidden = hidden.to(self.device)
            if self.model.norm is not None:
                hidden = self.model.norm(hidden)
            vectors.append(_pool_hidden_states(hidden, torch=self.torch, pooling=self.pooling))
        return self.torch.cat(vectors, dim=0)

    def save_pretrained(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "dnakit_evo2_config.json").write_text(
            json.dumps(self.model_config, indent=2),
            encoding="utf-8",
        )
        # FP8 is disabled; its opaque engine metadata is unnecessary and can make
        # torch.load(weights_only=True) depend on the original GPU/runtime.
        state = {
            name: value
            for name, value in self.model.state_dict().items()
            if not name.endswith("._extra_state")
        }
        self.torch.save(state, directory / "weights.pt")

    def load_pretrained(self, directory: Path) -> None:
        _restore_native_state(self.model, _load_weights(self.torch, directory / "weights.pt"))
        self._clear_inference_state()


__all__ = ["Evo2SequenceEncoder", "JanusDNASequenceEncoder"]
