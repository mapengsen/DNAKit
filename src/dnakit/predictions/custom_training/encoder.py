"""Trainable PyTorch encoders for full fine-tuning of Transformers DNA models."""

from __future__ import annotations

import importlib
from contextlib import suppress
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.representations import DNAEmbeddingModel


@runtime_checkable
class TrainableSequenceEncoder(Protocol):
    """Minimal encoder interface consumed by the custom-property trainer."""

    device: Any

    def encode(self, sequences: list[str]) -> Any:
        """Return a differentiable ``(batch, dimensions)`` tensor."""

    def parameters(self) -> Any:
        """Return trainable parameters."""

    def train(self) -> None:
        """Enable training behavior."""

    def eval(self) -> None:
        """Enable evaluation behavior."""

    def save_pretrained(self, directory: Path) -> None:
        """Persist a reloadable encoder and tokenizer."""


def _require(name: str) -> Any:
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        raise BackendUnavailableError(
            "Full custom-property fine-tuning dependencies are not installed.",
            code="MISSING_NEURAL_DEPENDENCY",
            context={"dependency": name},
            hint='Install the training extra with: python -m pip install "dnakit[training]"',
        ) from exc


class TransformersSequenceEncoder:
    """A pooled, differentiable wrapper around a local Transformers checkpoint."""

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
    ) -> None:
        if spec.loader != "transformers":
            raise ConfigurationError(
                "Full fine-tuning currently requires a Transformers DNA model.",
                code="CUSTOM_FULL_FINETUNE_UNSUPPORTED",
                context={"model": spec.name, "loader": spec.loader},
            )
        if spec.trust_remote_code and not allow_remote_code:
            raise ConfigurationError(
                "The selected DNA model requires checkpoint-provided Python code.",
                code="MODEL_REMOTE_CODE_NOT_ALLOWED",
                context={"model": spec.name},
                hint=(
                    "Review the official model repository and checkpoint, then set "
                    "allow_remote_code=True."
                ),
            )
        self.spec = spec
        self.checkpoint = Path(checkpoint).expanduser().resolve()
        self.pooling = pooling
        self.torch = _require("torch")
        transformers = _require("transformers")
        from dnakit.representations.backends import (
            _disable_fast_paths,
            _module_from_file,
            _pool_hidden_states,
            _resolve_hidden_states,
            _resolve_torch_device,
            _torch_dtype,
            _unwrap_base_model,
        )

        self._pool_hidden_states = _pool_hidden_states
        self._resolve_hidden_states = _resolve_hidden_states
        self.device = _resolve_torch_device(self.torch, device)
        if self.device.type == "cpu" and dtype in {"float16", "bfloat16"}:
            raise ConfigurationError(
                "Half-precision full fine-tuning is not supported on the CPU path.",
                code="INVALID_CUSTOM_PROPERTY_CONFIG",
                hint="Use dtype='float32' or dtype='auto' on CPU.",
            )
        if spec.name == "caduceus":
            tokenizer_module = _module_from_file(
                f"_dnakit_train_caduceus_tokenizer_{abs(hash(self.checkpoint))}",
                self.checkpoint / "tokenization_caduceus.py",
            )
            self.tokenizer = tokenizer_module.CaduceusTokenizer(
                model_max_length=int(max_length or min(spec.chunk_length, 4_096)),
                padding_side="right",
            )
        else:
            self.tokenizer = transformers.AutoTokenizer.from_pretrained(
                str(self.checkpoint),
                trust_remote_code=spec.trust_remote_code,
                local_files_only=True,
            )
            self.tokenizer.padding_side = "right"
        if getattr(self.tokenizer, "pad_token_id", None) is None:
            fallback = (
                getattr(self.tokenizer, "eos_token", None)
                or getattr(self.tokenizer, "sep_token", None)
                or getattr(self.tokenizer, "cls_token", None)
                or getattr(self.tokenizer, "unk_token", None)
            )
            if fallback is None:
                raise BackendExecutionError(
                    "The selected tokenizer has no padding token.",
                    code="MODEL_LOAD_FAILED",
                    context={"model": spec.name},
                )
            self.tokenizer.pad_token = fallback
        configured_length = int(max_length or min(spec.chunk_length, 4_096))
        tokenizer_limit = getattr(self.tokenizer, "model_max_length", configured_length)
        with suppress(TypeError, ValueError, OverflowError):
            parsed_limit = int(tokenizer_limit)
            if 0 < parsed_limit <= 10_000_000:
                configured_length = min(configured_length, parsed_limit)
        self.max_length = configured_length

        model_config = None
        errors: list[str] = []
        if spec.name == "caduceus":
            try:
                model_config = transformers.AutoConfig.from_pretrained(
                    str(self.checkpoint),
                    trust_remote_code=True,
                    local_files_only=True,
                )
                model_config.fused_add_norm = False
            except Exception as exc:
                errors.append(f"AutoConfig: {type(exc).__name__}: {exc}")
        kwargs: dict[str, object] = {
            "trust_remote_code": spec.trust_remote_code,
            "local_files_only": True,
        }
        if model_config is not None:
            kwargs["config"] = model_config
        model = None
        for loader_name in ("AutoModel", "AutoModelForMaskedLM"):
            loader = getattr(transformers, loader_name)
            try:
                model = loader.from_pretrained(str(self.checkpoint), **kwargs)
                if loader_name == "AutoModelForMaskedLM":
                    model = _unwrap_base_model(model, self.torch)
                break
            except Exception as exc:
                errors.append(f"{loader_name}: {type(exc).__name__}: {exc}")
        if model is None:
            raise BackendExecutionError(
                "Could not load the Transformers DNA checkpoint for fine-tuning.",
                code="MODEL_LOAD_FAILED",
                context={"model": spec.name, "errors": tuple(errors)},
            )
        _disable_fast_paths(model, spec)
        requested_dtype = _torch_dtype(self.torch, dtype)
        if requested_dtype is not None:
            model = model.to(dtype=requested_dtype)
        self.model = model.to(self.device)
        for parameter in self.model.parameters():
            parameter.requires_grad_(True)

    def parameters(self) -> Any:
        return self.model.parameters()

    def train(self) -> None:
        self.model.train()

    def eval(self) -> None:
        self.model.eval()

    def encode(self, sequences: list[str]) -> Any:
        try:
            tokens = self.tokenizer(
                sequences,
                add_special_tokens=True,
                padding=True,
                truncation=True,
                max_length=self.max_length,
                return_attention_mask=True,
                return_special_tokens_mask=True,
                return_tensors="pt",
            )
        except (TypeError, NotImplementedError):
            tokens = self.tokenizer(
                sequences,
                add_special_tokens=True,
                padding=True,
                truncation=True,
                max_length=self.max_length,
                return_attention_mask=True,
                return_tensors="pt",
            )
        special_tokens_mask = tokens.pop("special_tokens_mask", None)
        tokens = {key: value.to(self.device) for key, value in tokens.items()}
        try:
            outputs = self.model(**tokens, output_hidden_states=True, return_dict=True)
        except TypeError:
            outputs = self.model(**tokens)
        hidden = self._resolve_hidden_states(outputs, self.torch)
        return self._pool_hidden_states(
            hidden,
            torch=self.torch,
            pooling=self.pooling,
            attention_mask=tokens.get("attention_mask"),
            special_tokens_mask=(
                None if special_tokens_mask is None else special_tokens_mask.to(self.device)
            ),
        )

    def save_pretrained(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        try:
            self.model.save_pretrained(directory, safe_serialization=True)
        except TypeError:
            self.model.save_pretrained(directory)
        self.tokenizer.save_pretrained(directory)


__all__ = ["TrainableSequenceEncoder", "TransformersSequenceEncoder"]
