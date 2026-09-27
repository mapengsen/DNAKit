"""Lazy optional backends for extracting sequence-level DNA representations."""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import json
import math
import pickle
import sys
from collections import defaultdict
from collections.abc import Iterable, Sequence
from contextlib import suppress
from pathlib import Path
from typing import Any, Protocol

from dnakit._dependency_hints import model_dependency_hint
from dnakit._verbosity import progress_enabled
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from .models import DNAEmbeddingModel, RepresentationConfig


class RepresentationBackend(Protocol):
    """Backend contract used by :func:`extract_representations`."""

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        """Return one finite numeric vector per input sequence."""


def _missing_dependency(spec: DNAEmbeddingModel, package: str) -> BackendUnavailableError:
    requirement = spec.required_package or package
    return BackendUnavailableError(
        f"The {spec.display_name} representation backend is not installed.",
        code="MISSING_NEURAL_DEPENDENCY",
        context={"model": spec.name, "dependency": requirement},
        hint=model_dependency_hint(spec.name, requirement),
    )


def _require_module(name: str, spec: DNAEmbeddingModel) -> Any:
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        raise _missing_dependency(spec, name) from exc


def _progress(values: Sequence[str], *, enabled: bool, description: str) -> Iterable[str]:
    if not progress_enabled(enabled):
        return values
    from dnakit._progress import track

    return track(values, description=description)


def _prepend_source_path(value: str | None) -> None:
    if value is None:
        return
    root = Path(value).expanduser().resolve()
    for candidate in (root, root / "src"):
        if not candidate.is_dir():
            continue
        candidate_text = str(candidate)
        if candidate_text not in sys.path:
            sys.path.insert(0, candidate_text)


def _resolve_torch_device(torch: Any, requested: str) -> Any:
    normalized = requested.strip().lower()
    if normalized == "auto":
        if bool(torch.cuda.is_available()):
            normalized = "cuda"
        elif hasattr(torch.backends, "mps") and bool(torch.backends.mps.is_available()):
            normalized = "mps"
        else:
            normalized = "cpu"
    try:
        device = torch.device(normalized)
    except (RuntimeError, TypeError) as exc:
        raise ConfigurationError(
            "Invalid PyTorch device.",
            code="INVALID_REPRESENTATION_DEVICE",
            context={"device": requested},
        ) from exc
    if device.type == "cuda" and not bool(torch.cuda.is_available()):
        raise ConfigurationError(
            "Requested CUDA device is unavailable.",
            code="INVALID_REPRESENTATION_DEVICE",
            context={"device": requested},
        )
    if (
        device.type == "cuda"
        and device.index is not None
        and device.index >= torch.cuda.device_count()
    ):
        raise ConfigurationError(
            "Requested CUDA device index is unavailable.",
            code="INVALID_REPRESENTATION_DEVICE",
            context={"device": requested, "device_count": torch.cuda.device_count()},
        )
    if device.type == "mps" and (
        not hasattr(torch.backends, "mps") or not bool(torch.backends.mps.is_available())
    ):
        raise ConfigurationError(
            "Requested MPS device is unavailable.",
            code="INVALID_REPRESENTATION_DEVICE",
            context={"device": requested},
        )
    return device


def _torch_dtype(torch: Any, value: str) -> Any | None:
    if value == "auto":
        return None
    return {
        "float32": torch.float32,
        "float16": torch.float16,
        "bfloat16": torch.bfloat16,
    }[value]


def _resolve_hidden_states(outputs: Any, torch: Any) -> Any:
    candidate = getattr(outputs, "last_hidden_state", None)
    if candidate is not None:
        return candidate
    hidden_states = getattr(outputs, "hidden_states", None)
    if hidden_states:
        return hidden_states[-1]
    if isinstance(outputs, tuple) and outputs and isinstance(outputs[0], torch.Tensor):
        return outputs[0]
    raise BackendExecutionError(
        "Model output does not contain hidden representations.",
        code="MISSING_MODEL_HIDDEN_STATES",
    )


def _pool_hidden_states(
    hidden_states: Any,
    *,
    torch: Any,
    pooling: str,
    attention_mask: Any | None = None,
    special_tokens_mask: Any | None = None,
) -> Any:
    if hidden_states.ndim != 3:
        raise BackendExecutionError(
            "Expected hidden states with shape (batch, tokens, dimensions).",
            code="INVALID_MODEL_HIDDEN_STATES",
            context={"shape": tuple(int(value) for value in hidden_states.shape)},
        )
    if attention_mask is None:
        mask = torch.ones(
            hidden_states.shape[:2],
            device=hidden_states.device,
            dtype=torch.bool,
        )
    else:
        mask = attention_mask.to(hidden_states.device).bool()
    if special_tokens_mask is not None:
        without_special = mask & (~special_tokens_mask.to(hidden_states.device).bool())
        if without_special.any(dim=1).all():
            mask = without_special
    if pooling == "cls":
        return hidden_states[:, 0, :]
    if pooling == "last":
        positions = torch.arange(hidden_states.shape[1], device=hidden_states.device)[None, :]
        indices = positions.masked_fill(~mask, -1).max(dim=1).values.clamp_min(0)
        gather_index = indices[:, None, None].expand(-1, 1, hidden_states.shape[-1])
        return hidden_states.gather(1, gather_index).squeeze(1)
    weights = mask.to(hidden_states.dtype).unsqueeze(-1)
    if pooling == "mean":
        return (hidden_states * weights).sum(dim=1) / weights.sum(dim=1).clamp_min(1.0)
    if pooling == "max":
        minimum = torch.finfo(hidden_states.dtype).min
        return hidden_states.masked_fill(~mask.unsqueeze(-1), minimum).max(dim=1).values
    raise AssertionError(f"Unsupported pooling method: {pooling}")


def _module_from_file(module_name: str, file_path: Path) -> Any:
    module_spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    if module_spec is None or module_spec.loader is None:
        raise ImportError(f"Could not construct module {module_name} from {file_path}.")
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[module_name] = module
    module_spec.loader.exec_module(module)
    return module


def _unwrap_base_model(model: Any, torch: Any) -> Any:
    for name in ("esm", "bert", "roberta", "transformer", "base_model", "model"):
        candidate = getattr(model, name, None)
        if isinstance(candidate, torch.nn.Module):
            return candidate
    return model


def _disable_fast_paths(model: Any, spec: DNAEmbeddingModel) -> None:
    if spec.name == "dnabert2":
        model_modules = {type(layer).__module__ for layer in model.modules()}
        for module in tuple(sys.modules.values()):
            module_file = getattr(module, "__file__", "") or ""
            if (
                spec.cache_name not in module_file
                and "DNABERT-2-117M" not in module_file
                and getattr(module, "__name__", "") not in model_modules
            ):
                continue
            if hasattr(module, "flash_attn_qkvpacked_func"):
                module.__dict__["flash_attn_qkvpacked_func"] = None
    if spec.name == "caduceus":
        for module in model.modules():
            if hasattr(module, "use_fast_path"):
                with suppress(AttributeError, RuntimeError):
                    module.use_fast_path = False


def _transformer_special_tokens_mask(tokenizer: Any, tokens: Any, torch: Any) -> Any:
    """Repair missing/malformed special masks without dropping unknown DNA tokens."""

    provided = tokens.get("special_tokens_mask")
    if provided is not None and provided.shape == tokens["input_ids"].shape:
        return provided
    # NT-v2's EsmTokenizer can include a phantom EOS position in this mask
    # despite eos_token=None. Rebuild from the actual IDs in that case.
    special_mask = torch.zeros_like(tokens["input_ids"], dtype=torch.long)
    unknown = getattr(tokenizer, "unk_token_id", None)
    for token_id in tuple(getattr(tokenizer, "all_special_ids", ()) or ()):
        if token_id != unknown:
            special_mask |= tokens["input_ids"] == int(token_id)
    return special_mask


def _forward_transformers_batch(
    model: Any,
    inputs: dict[str, Any],
    *,
    torch: Any,
    pooling: str,
    special_tokens_mask: Any | None,
) -> Any:
    """Pool differentiable outputs, respecting models without attention masks.

    Bidirectional state-space encoders cannot ignore padded positions during
    their forward pass. Run their variable-length inputs separately so a
    sequence has the same representation regardless of its batch neighbours.
    """

    signature = inspect.signature(model.forward)
    names = signature.parameters
    accepts_kwargs = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in names.values())
    accepts_mask = "attention_mask" in names or accepts_kwargs
    attention_mask = inputs.get("attention_mask")

    def forward(batch: dict[str, Any], special_mask: Any | None) -> Any:
        arguments = {key: value for key, value in batch.items() if key in names or accepts_kwargs}
        # Token type IDs are not a generic **kwargs option (e.g. Llama routes
        # extra keyword arguments into its attention implementation).
        if "token_type_ids" not in names:
            arguments.pop("token_type_ids", None)
        if "return_dict" in names:
            arguments["return_dict"] = True
        if "use_cache" in names:
            arguments["use_cache"] = False
        outputs = model(**arguments)
        hidden = (
            outputs if isinstance(outputs, torch.Tensor) else _resolve_hidden_states(outputs, torch)
        )
        return _pool_hidden_states(
            hidden,
            torch=torch,
            pooling=pooling,
            attention_mask=batch.get("attention_mask"),
            special_tokens_mask=special_mask,
        )

    if accepts_mask or attention_mask is None or bool(attention_mask.bool().all()):
        return forward(inputs, special_tokens_mask)
    vectors = []
    for row in range(inputs["input_ids"].shape[0]):
        valid = attention_mask[row].bool()
        batch = {key: value[row : row + 1, valid] for key, value in inputs.items()}
        special = None if special_tokens_mask is None else special_tokens_mask[row : row + 1, valid]
        vectors.append(forward(batch, special))
    return torch.cat(vectors, dim=0)


class _TransformersBackend:
    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: Path,
        config: RepresentationConfig,
    ) -> None:
        self.spec = spec
        self.checkpoint = checkpoint
        self.config = config
        self.torch = _require_module("torch", spec)
        transformers = _require_module("transformers", spec)
        self.device = _resolve_torch_device(self.torch, config.device)
        self.tokenizer = self._load_tokenizer(transformers)
        self.model = self._load_model(transformers)
        requested_dtype = _torch_dtype(self.torch, config.dtype)
        if requested_dtype is not None:
            self.model = self.model.to(dtype=requested_dtype)
        self.model = self.model.to(self.device)
        self.model.eval()
        self.window_length = self._window_length()

    def _load_tokenizer(self, transformers: Any) -> Any:
        if self.spec.name == "caduceus":
            module = _module_from_file(
                f"_dnakit_caduceus_tokenizer_{abs(hash(self.checkpoint))}",
                self.checkpoint / "tokenization_caduceus.py",
            )
            tokenizer_class = module.CaduceusTokenizer
            return tokenizer_class(
                model_max_length=int(self.config.max_length or self.spec.chunk_length),
                padding_side="right",
            )
        tokenizer = transformers.AutoTokenizer.from_pretrained(
            str(self.checkpoint),
            trust_remote_code=self.spec.trust_remote_code,
            local_files_only=True,
        )
        tokenizer.padding_side = "right"
        if getattr(tokenizer, "pad_token_id", None) is None:
            fallback = (
                getattr(tokenizer, "eos_token", None)
                or getattr(tokenizer, "sep_token", None)
                or getattr(tokenizer, "cls_token", None)
                or getattr(tokenizer, "unk_token", None)
            )
            if fallback is None:
                raise BackendExecutionError(
                    "The selected tokenizer has no usable padding token.",
                    code="MODEL_LOAD_FAILED",
                    context={"model": self.spec.name},
                )
            tokenizer.pad_token = fallback
        return tokenizer

    def _load_model(self, transformers: Any) -> Any:
        errors: list[str] = []
        model_config = None
        if self.spec.name == "caduceus":
            try:
                model_config = transformers.AutoConfig.from_pretrained(
                    str(self.checkpoint),
                    trust_remote_code=True,
                    local_files_only=True,
                )
                model_config.fused_add_norm = False
            except Exception as exc:  # third-party configuration errors are reported together.
                errors.append(f"AutoConfig: {type(exc).__name__}: {exc}")
        kwargs = {
            "trust_remote_code": self.spec.trust_remote_code,
            "local_files_only": True,
        }
        if model_config is not None:
            kwargs["config"] = model_config
        for loader_name in ("AutoModel", "AutoModelForMaskedLM"):
            loader = getattr(transformers, loader_name)
            try:
                model = loader.from_pretrained(str(self.checkpoint), **kwargs)
                if loader_name == "AutoModelForMaskedLM":
                    model = _unwrap_base_model(model, self.torch)
                _disable_fast_paths(model, self.spec)
                return model
            except Exception as exc:  # preserve both official loader diagnostics.
                errors.append(f"{loader_name}: {type(exc).__name__}: {exc}")
        raise BackendExecutionError(
            "Could not load the Transformers checkpoint.",
            code="MODEL_LOAD_FAILED",
            context={"model": self.spec.name, "errors": tuple(errors)},
        )

    def _window_length(self) -> int:
        configured = int(self.config.max_length or self.spec.chunk_length)
        model_max = getattr(self.tokenizer, "model_max_length", configured)
        try:
            model_max_int = int(model_max)
        except (TypeError, ValueError, OverflowError):
            model_max_int = configured
        if not 0 < model_max_int <= 10_000_000:
            model_max_int = configured
        special_count = 0
        with suppress(AttributeError, TypeError, ValueError):
            special_count = int(self.tokenizer.num_special_tokens_to_add(pair=False))
        return max(1, min(configured, model_max_int) - special_count)

    def _special_token_count(self) -> int:
        method = getattr(self.tokenizer, "num_special_tokens_to_add", None)
        if method is None:
            return 0
        try:
            return int(method(pair=False))
        except (TypeError, ValueError):
            return 0

    def _chunks(self, sequence: str) -> tuple[list[Any], str]:
        if self.spec.chunk_unit == "bp":
            bp_chunks: list[Any] = [
                sequence[start : start + self.window_length]
                for start in range(0, len(sequence), self.window_length)
            ]
            return bp_chunks, "bp"
        encoded = self.tokenizer(
            sequence,
            add_special_tokens=False,
            return_attention_mask=False,
            return_token_type_ids=False,
        )
        token_ids = list(encoded["input_ids"])
        if not token_ids:
            raise BackendExecutionError(
                "The selected tokenizer produced no tokens for a non-empty sequence.",
                code="EMPTY_MODEL_TOKENIZATION",
                context={"model": self.spec.name},
            )
        if getattr(self.tokenizer, "is_fast", False):
            # prepare_for_model() does not apply a fast tokenizer's Rust
            # post-processor. Keep its Encoding objects so each window receives
            # the official CLS/SEP tokens and their exact special-token mask.
            encoding = encoded.encodings[0]
            encoding.truncate(self.window_length, stride=0, direction="right")
            return [encoding, *encoding.overflowing], "encoding"
        token_chunks: list[Any] = [
            token_ids[start : start + self.window_length]
            for start in range(0, len(token_ids), self.window_length)
        ]
        return token_chunks, "token"

    def _tokenize(self, chunks: list[Any], kind: str) -> dict[str, Any]:
        if kind == "bp":
            kwargs = {
                "add_special_tokens": True,
                "padding": True,
                "truncation": True,
                "max_length": self.window_length + self._special_token_count(),
                "return_attention_mask": True,
                "return_special_tokens_mask": True,
                "return_tensors": "pt",
            }
            try:
                batch = self.tokenizer(chunks, **kwargs)
            except (TypeError, NotImplementedError):
                kwargs.pop("return_special_tokens_mask")
                batch = self.tokenizer(chunks, **kwargs)
        elif kind == "encoding":
            prepared = []
            for chunk in chunks:
                encoding = self.tokenizer.backend_tokenizer.post_process(
                    chunk, add_special_tokens=True
                )
                prepared.append(
                    {
                        "input_ids": encoding.ids,
                        "attention_mask": encoding.attention_mask,
                        "special_tokens_mask": encoding.special_tokens_mask,
                    }
                )
            batch = self.tokenizer.pad(prepared, padding=True, return_tensors="pt")
        else:
            prepared = [
                self.tokenizer.prepare_for_model(
                    token_ids,
                    add_special_tokens=True,
                    return_attention_mask=True,
                    return_token_type_ids=False,
                    return_special_tokens_mask=not getattr(self.tokenizer, "is_fast", False),
                    truncation=True,
                    max_length=self.window_length + self._special_token_count(),
                )
                for token_ids in chunks
            ]
            batch = self.tokenizer.pad(prepared, padding=True, return_tensors="pt")
        values = dict(batch)
        values.pop("token_type_ids", None)
        values["special_tokens_mask"] = _transformer_special_tokens_mask(
            self.tokenizer, values, self.torch
        )
        return values

    def _encode_chunk_batch(self, chunks: list[Any], kind: str) -> Any:
        batch = self._tokenize(chunks, kind)
        special_tokens_mask = batch.pop("special_tokens_mask", None)
        inputs = {
            key: value.to(self.device)
            for key, value in batch.items()
            if isinstance(value, self.torch.Tensor)
        }
        with self.torch.inference_mode():
            return _forward_transformers_batch(
                self.model,
                inputs,
                torch=self.torch,
                pooling=self.config.pooling,
                special_tokens_mask=(
                    None if special_tokens_mask is None else special_tokens_mask.to(self.device)
                ),
            )

    def _encode_sequence(self, sequence: str) -> Any:
        chunks, kind = self._chunks(sequence)
        vectors: list[Any] = []
        for start in range(0, len(chunks), self.config.batch_size):
            vectors.append(
                self._encode_chunk_batch(chunks[start : start + self.config.batch_size], kind)
            )
        return self.torch.cat(vectors, dim=0).mean(dim=0)

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        vectors = [
            self._encode_sequence(sequence)
            for sequence in _progress(
                sequences,
                enabled=show_progress,
                description=f"Extracting {self.spec.display_name} rep",
            )
        ]
        return self.torch.stack(vectors, dim=0).detach().float().cpu().numpy()


class _EnformerBackend:
    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: Path,
        config: RepresentationConfig,
    ) -> None:
        self.spec = spec
        self.config = config
        self.torch = _require_module("torch", spec)
        enformer_pytorch = _require_module("enformer_pytorch", spec)
        self.device = _resolve_torch_device(self.torch, config.device)
        self.chunk_length = int(config.max_length or spec.chunk_length)
        try:
            self.model = enformer_pytorch.from_pretrained(
                str(checkpoint),
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Enformer checkpoint.",
                code="MODEL_LOAD_FAILED",
                context={"model": spec.name},
            ) from exc
        # ``enformer_pytorch.from_pretrained`` applies this correction only when
        # it receives the remote repository ID.  DNAKit intentionally loads from
        # the resolved local directory, so enable the same official setting here.
        for module in self.model.modules():
            if hasattr(module, "use_tf_gamma"):
                module.use_tf_gamma = self.chunk_length == 196_608
        self.downsample = 2 ** int(self.model.config.num_downsamples)
        self.short_context = self.chunk_length != spec.chunk_length
        if self.short_context:
            self.model.crop_final.target_length = -1
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            self.model = self.model.to(dtype=dtype)
        self.model = self.model.to(self.device).eval()
        lookup = self.torch.zeros(256, 4, dtype=self.torch.float32, device=self.device)
        for base, vector in {
            "A": (1.0, 0.0, 0.0, 0.0),
            "C": (0.0, 1.0, 0.0, 0.0),
            "G": (0.0, 0.0, 1.0, 0.0),
            "T": (0.0, 0.0, 0.0, 1.0),
            "N": (0.0, 0.0, 0.0, 0.0),
        }.items():
            lookup[ord(base)] = self.torch.tensor(vector, device=self.device)
        self.lookup = lookup

    def _encode_sequence(self, sequence: str) -> Any:
        raw_chunks = [
            sequence[start : start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        padded_length = max(
            2 * self.downsample,
            math.ceil(self.chunk_length / self.downsample) * self.downsample,
        )
        chunks = [chunk + "N" * (padded_length - len(chunk)) for chunk in raw_chunks]
        vectors: list[Any] = []
        for start in range(0, len(chunks), self.config.batch_size):
            batch_chunks = chunks[start : start + self.config.batch_size]
            encoded = self.torch.tensor(
                [list(chunk.encode("ascii")) for chunk in batch_chunks],
                dtype=self.torch.long,
                device=self.device,
            )
            one_hot = self.lookup[encoded]
            with self.torch.inference_mode():
                hidden = self.model(
                    one_hot.to(dtype=next(self.model.parameters()).dtype),
                    return_only_embeddings=True,
                )
            mask = None
            if self.short_context:
                lengths = self.torch.tensor(
                    [
                        math.ceil(len(chunk) / self.downsample)
                        for chunk in raw_chunks[start : start + self.config.batch_size]
                    ],
                    device=self.device,
                )
                mask = (
                    self.torch.arange(hidden.shape[1], device=self.device)[None, :]
                    < (lengths[:, None])
                )
            vectors.append(
                _pool_hidden_states(
                    hidden,
                    torch=self.torch,
                    pooling=self.config.pooling,
                    attention_mask=mask,
                )
            )
        return self.torch.cat(vectors, dim=0).mean(dim=0)

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        vectors = [
            self._encode_sequence(sequence)
            for sequence in _progress(
                sequences,
                enabled=show_progress,
                description="Extracting Enformer rep",
            )
        ]
        return self.torch.stack(vectors).detach().float().cpu().numpy()


def _clean_janusdna_config(config: dict[str, Any]) -> dict[str, Any]:
    cleaned = dict(config)
    cleaned.pop("_target_", None)
    for key in (
        "bidirectional",
        "bidirectional_weight_tie",
        "layer_fusion",
        "final_attention",
        "mid_single_direction_attention",
        "bidirectional_attn_tie",
        "gradient_checkpointing",
        "output_router_logits",
    ):
        value = cleaned.get(key)
        if isinstance(value, str):
            cleaned[key] = value.lower().strip(", ") == "true"
    cleaned.update(
        {
            "gradient_checkpointing": False,
            "output_router_logits": False,
            "use_cache": False,
            "use_mamba_kernels": False,
            "attn_implementation": "eager",
            "final_attention_class": "eager",
        }
    )
    return cleaned


def _load_janusdna_checkpoint(torch: Any, checkpoint: Path, spec: DNAEmbeddingModel) -> Any:
    """Read official Lightning weights with a fixed, scoped metadata allowlist.

    OmegaConf's container/metadata classes restore plain data and annotations;
    they do not invoke configuration resolvers during deserialization. Never
    extend this list from checkpoint-provided names or enable unrestricted pickle.
    DNAKit's tensor-only ``weights.pt`` artifacts do not use this loader.
    """
    try:
        return torch.load(checkpoint, map_location="cpu", weights_only=True)
    except pickle.UnpicklingError:
        configuration = _require_module("omegaconf", spec)
        base = _require_module("omegaconf.base", spec)
        nodes = _require_module("omegaconf.nodes", spec)
        metadata_types = [
            base.Metadata,
            base.ContainerMetadata,
            nodes.AnyNode,
            configuration.DictConfig,
            configuration.ListConfig,
            list,
            dict,
            int,
            defaultdict,
            Any,
        ]
        existing = torch.serialization.get_safe_globals()
        additions = [item for item in metadata_types if item not in existing]
        with torch.serialization.safe_globals(additions):
            return torch.load(checkpoint, map_location="cpu", weights_only=True)


class _JanusDNABackend:
    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: Path,
        config: RepresentationConfig,
    ) -> None:
        self.spec = spec
        self.config = config
        source_path = None if config.model_source_path is None else str(config.model_source_path)
        _prepend_source_path(source_path)
        self.torch = _require_module("torch", spec)
        try:
            configuration_module = importlib.import_module("janusdna.configuration_janusdna")
            modeling_module = importlib.import_module("janusdna.modeling_janusdna")
        except ImportError as exc:
            raise BackendUnavailableError(
                "JanusDNA source code is required to load its official checkpoint.",
                code="MISSING_NEURAL_DEPENDENCY",
                hint=model_dependency_hint(spec.name, spec.required_package),
            ) from exc
        payload = json.loads((checkpoint / "model_config.json").read_text(encoding="utf-8"))
        model_config = configuration_module.JanusDNAConfig(
            **_clean_janusdna_config(dict(payload["config"]))
        )
        model_config._attn_implementation = "eager"
        model_config.final_attention_class = "eager"
        causal_model = modeling_module.JanusDNAForCausalLM(model_config)
        state = _load_janusdna_checkpoint(
            self.torch, checkpoint / "checkpoints" / "last.ckpt", spec
        )
        state_dict = dict(state["state_dict"])
        self.torch.nn.modules.utils.consume_prefix_in_state_dict_if_present(state_dict, "model.")
        for key in tuple(state_dict):
            if ".self_attn.o_projs.0." in key:
                state_dict[key.replace(".self_attn.o_projs.0.", ".self_attn.o_proj.")] = (
                    state_dict.pop(key)
                )
            elif "torchmetrics" in key:
                state_dict.pop(key)
        missing, unexpected = causal_model.load_state_dict(state_dict, strict=False)
        missing = [key for key in missing if not key.startswith("lm_head.")]
        if missing or unexpected:
            raise BackendExecutionError(
                "JanusDNA checkpoint did not match the official model architecture.",
                code="MODEL_LOAD_FAILED",
                context={"missing": tuple(missing[:10]), "unexpected": tuple(unexpected[:10])},
            )
        self.device = _resolve_torch_device(self.torch, config.device)
        self.model = causal_model.model
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            self.model = self.model.to(dtype=dtype)
        self.model = self.model.to(self.device).eval()
        self.chunk_length = int(config.max_length or spec.chunk_length)
        self.base_ids = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}

    def _encode_sequence(self, sequence: str) -> Any:
        chunks = [
            sequence[start : start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        vectors: list[Any] = []
        # JanusDNA reverses state-space streams internally. Padding a short tail
        # alongside a full chunk alters its representation even with a mask.
        chunks_by_length: dict[int, list[str]] = {}
        for chunk in chunks:
            chunks_by_length.setdefault(len(chunk), []).append(chunk)
        batches = [
            group[start : start + self.config.batch_size]
            for group in chunks_by_length.values()
            for start in range(0, len(group), self.config.batch_size)
        ]
        for batch_chunks in batches:
            max_size = max(len(chunk) for chunk in batch_chunks) + 2
            input_ids = self.torch.full(
                (len(batch_chunks), max_size),
                4,
                dtype=self.torch.long,
                device=self.device,
            )
            attention_mask = self.torch.zeros_like(input_ids)
            special_mask = self.torch.ones_like(input_ids)
            for row, chunk in enumerate(batch_chunks):
                ids = [0, *(self.base_ids[base] for base in chunk), 1]
                input_ids[row, : len(ids)] = self.torch.tensor(ids, device=self.device)
                attention_mask[row, : len(ids)] = 1
                special_mask[row, 1 : len(ids) - 1] = 0
            with self.torch.inference_mode():
                outputs = self.model(input_ids=input_ids, attention_mask=attention_mask)
            hidden = _resolve_hidden_states(outputs, self.torch)
            vectors.append(
                _pool_hidden_states(
                    hidden,
                    torch=self.torch,
                    pooling=self.config.pooling,
                    attention_mask=attention_mask,
                    special_tokens_mask=special_mask,
                )
            )
        return self.torch.cat(vectors).mean(dim=0)

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        vectors = [
            self._encode_sequence(sequence)
            for sequence in _progress(
                sequences,
                enabled=show_progress,
                description="Extracting JanusDNA rep",
            )
        ]
        return self.torch.stack(vectors).detach().float().cpu().numpy()


class _Evo2Backend:
    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: Path,
        config: RepresentationConfig,
    ) -> None:
        self.spec = spec
        self.config = config
        _prepend_source_path(
            None if config.model_source_path is None else str(config.model_source_path)
        )
        self.torch = _require_module("torch", spec)
        evo2_module = _require_module("evo2", spec)
        checkpoint_file = checkpoint / "evo2_7b.pt"
        try:
            self.wrapper = evo2_module.Evo2(
                model_name="evo2_7b",
                local_path=str(checkpoint_file),
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Evo 2 checkpoint.",
                code="MODEL_LOAD_FAILED",
            ) from exc
        self.model = self.wrapper.model.eval()
        self.device = _resolve_torch_device(self.torch, config.device)
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            self.model = self.model.to(dtype=dtype)
        self.model = self.model.to(self.device)
        self.chunk_length = int(config.max_length or spec.chunk_length)
        if spec.evo2_layer is None:
            raise AssertionError("Evo 2 registry entry is missing an embedding layer.")
        self.layer = spec.evo2_layer

    def _encode_sequence(self, sequence: str) -> Any:
        chunks = [
            sequence[start : start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        vectors: list[Any] = []
        for chunk in chunks:
            input_ids = self.torch.tensor(
                self.wrapper.tokenizer.tokenize(chunk),
                dtype=self.torch.int,
                device=self.device,
            ).unsqueeze(0)
            with self.torch.inference_mode():
                _, embeddings = self.wrapper(
                    input_ids,
                    return_embeddings=True,
                    layer_names=[self.layer],
                )
            vectors.append(
                _pool_hidden_states(
                    embeddings[self.layer],
                    torch=self.torch,
                    pooling=self.config.pooling,
                ).squeeze(0)
            )
        return self.torch.stack(vectors).mean(dim=0)

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        vectors = [
            self._encode_sequence(sequence)
            for sequence in _progress(
                sequences,
                enabled=show_progress,
                description="Extracting Evo 2 rep",
            )
        ]
        return self.torch.stack(vectors).detach().float().cpu().numpy()


def _build_alphagenome_apply_fn(module: Any, metadata: Any, *, force_float32: bool) -> Any:
    if not force_float32:
        # The official API added trunk/head apply functions after the original
        # three-item return value. The prediction apply function remains second.
        return module.create_model(metadata)[1]
    hk = module.hk
    jmp = module.jmp
    research_model = module.model
    policy = jmp.get_policy("params=float32,compute=float32,output=float32")

    @hk.transform_with_state  # type: ignore[untyped-decorator]
    def forward(dna_sequence: Any, organism_index: Any) -> Any:
        with hk.mixed_precision.push_policy(research_model.AlphaGenome, policy):
            return research_model.AlphaGenome(
                metadata,
                num_splice_sites=research_model.DEFAULT_NUM_SPLICE_SITES,
                splice_site_threshold=research_model.DEFAULT_SPLICE_SITE_THRESHOLD,
            )(dna_sequence, organism_index)

    def apply(params: Any, state: Any, dna_sequence: Any, organism_index: Any) -> Any:
        (predictions, _), _ = forward.apply(
            params,
            state,
            None,
            dna_sequence,
            organism_index,
        )
        return predictions

    return apply


def _cast_jax_float32(tree: Any, *, jax: Any, jnp: Any) -> Any:
    def cast(value: Any) -> Any:
        dtype = getattr(value, "dtype", None)
        if dtype is None:
            return value
        try:
            if jnp.issubdtype(dtype, jnp.floating):
                return value.astype(jnp.float32)
        except TypeError:
            pass
        return value

    return jax.tree_util.tree_map(cast, tree)


class _AlphaGenomeBackend:
    def __init__(
        self,
        spec: DNAEmbeddingModel,
        checkpoint: Path,
        config: RepresentationConfig,
    ) -> None:
        self.spec = spec
        self.config = config
        _prepend_source_path(
            None if config.model_source_path is None else str(config.model_source_path)
        )
        try:
            jax = importlib.import_module("jax")
            jnp = importlib.import_module("jax.numpy")
            ocp = importlib.import_module("orbax.checkpoint")
            api_dna_model = importlib.import_module("alphagenome.models.dna_model")
            dna_model = importlib.import_module("alphagenome_research.model.dna_model")
            metadata_lib = importlib.import_module("alphagenome_research.model.metadata.metadata")
        except ImportError as exc:
            raise _missing_dependency(spec, "alphagenome_research") from exc
        self.jax = jax
        self.jnp = jnp
        self.api_dna_model = api_dna_model
        convert_organism = getattr(dna_model, "convert_to_organism_index", None)
        human = api_dna_model.Organism.HOMO_SAPIENS
        self.organism_index = (
            int(convert_organism(human)) if convert_organism is not None else int(human.value)
        )
        requested = config.device.strip().lower()
        if requested == "auto":
            devices = list(jax.devices())
            self.device = next(
                (device for device in devices if device.platform == "gpu"),
                devices[0],
            )
        else:
            platform, separator, index_text = requested.partition(":")
            if platform == "cuda":
                platform = "gpu"
            try:
                matching = list(jax.devices(platform))
            except RuntimeError as exc:
                raise ConfigurationError(
                    "Requested JAX device is unavailable.",
                    code="INVALID_REPRESENTATION_DEVICE",
                    context={"device": config.device},
                ) from exc
            try:
                device_index = 0 if not separator else int(index_text)
            except ValueError as exc:
                raise ConfigurationError(
                    "Invalid JAX device index.",
                    code="INVALID_REPRESENTATION_DEVICE",
                    context={"device": config.device},
                ) from exc
            if not matching or not 0 <= device_index < len(matching):
                raise ConfigurationError(
                    "Requested JAX device is unavailable.",
                    code="INVALID_REPRESENTATION_DEVICE",
                    context={"device": config.device},
                )
            self.device = matching[device_index]
        if config.dtype not in {"auto", "float32"}:
            raise ConfigurationError(
                "AlphaGenome controls mixed precision in its official JAX implementation.",
                code="INVALID_REPRESENTATION_DTYPE",
                hint='Use dtype="auto" or dtype="float32" for AlphaGenome.',
            )
        metadata = {organism: metadata_lib.load(organism) for organism in api_dna_model.Organism}
        params, state = ocp.StandardCheckpointer().restore(
            str(checkpoint),
            target=None,
            strict=False,
        )
        if self.device.platform != "gpu":
            params = _cast_jax_float32(params, jax=jax, jnp=jnp)
            state = _cast_jax_float32(state, jax=jax, jnp=jnp)
        self.params = jax.device_put(params, self.device)
        self.state = jax.device_put(state, self.device)
        apply_fn = _build_alphagenome_apply_fn(
            dna_model,
            metadata,
            force_float32=self.device.platform != "gpu" or config.dtype == "float32",
        )
        self.apply_fn = jax.jit(apply_fn, device=self.device)
        self.chunk_length = int(config.max_length or spec.chunk_length)
        self.minimum_length = 2_048
        self.alignment = 2_048

    def _pad(self, sequence: str) -> str:
        size = max(len(sequence), self.minimum_length)
        size = int(math.ceil(size / self.alignment) * self.alignment)
        return sequence[:size] + "N" * max(0, size - len(sequence))

    @staticmethod
    def _one_hot(sequence: str) -> Any:
        import numpy as np

        lookup = {
            "A": (1.0, 0.0, 0.0, 0.0),
            "C": (0.0, 1.0, 0.0, 0.0),
            "G": (0.0, 0.0, 1.0, 0.0),
            "T": (0.0, 0.0, 0.0, 1.0),
            "N": (0.0, 0.0, 0.0, 0.0),
        }
        return np.asarray([lookup[base] for base in sequence], dtype=np.float32)[None, ...]

    def _pool(self, embeddings: Any) -> Any:
        import numpy as np

        # Reduce on the model device, matching the differentiable encoder and
        # avoiding a host copy of every base's embedding for long sequences.
        values = self.jnp.asarray(embeddings, dtype=self.jnp.float32)
        if self.config.pooling == "mean":
            pooled = values.mean(axis=1).squeeze(0)
        elif self.config.pooling == "max":
            pooled = values.max(axis=1).squeeze(0)
        elif self.config.pooling == "cls":
            pooled = values[:, 0, :].squeeze(0)
        else:
            pooled = values[:, -1, :].squeeze(0)
        return np.asarray(pooled, dtype=np.float32)

    def _encode_sequence(self, sequence: str) -> Any:
        import numpy as np

        chunks = [
            sequence[start : start + self.chunk_length]
            for start in range(0, len(sequence), self.chunk_length)
        ]
        vectors = []
        for chunk in chunks:
            padded = self._pad(chunk)
            dna = self.jax.device_put(self._one_hot(padded), self.device)
            organism = self.jax.device_put(
                self.jnp.asarray(
                    [self.organism_index],
                    dtype=self.jnp.int32,
                ),
                self.device,
            )
            predictions = self.apply_fn(self.params, self.state, dna, organism)
            vectors.append(self._pool(predictions["embeddings_1bp"]))
        return np.stack(vectors).mean(axis=0)

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        import numpy as np

        return np.stack(
            [
                self._encode_sequence(sequence)
                for sequence in _progress(
                    sequences,
                    enabled=show_progress,
                    description="Extracting AlphaGenome rep",
                )
            ]
        ).astype(np.float32, copy=False)


def create_representation_backend(
    spec: DNAEmbeddingModel,
    checkpoint: str | Path,
    config: RepresentationConfig,
    *,
    verbose: bool = True,
) -> RepresentationBackend:
    """Create one lazy backend for a locally resolved official checkpoint."""

    path = Path(checkpoint).expanduser().resolve()
    if spec.loader == "transformers":
        return _TransformersBackend(spec, path, config)
    if spec.loader == "enformer":
        return _EnformerBackend(spec, path, config)
    if spec.loader == "janusdna":
        return _JanusDNABackend(spec, path, config)
    if spec.loader == "evo2":
        return _Evo2Backend(spec, path, config)
    if spec.loader == "alphagenome":
        return _AlphaGenomeBackend(spec, path, config)
    raise AssertionError(f"Unsupported model loader: {spec.loader}")


__all__ = ["RepresentationBackend", "create_representation_backend"]

_register_public_verbose(globals())
