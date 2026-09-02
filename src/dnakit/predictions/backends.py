"""Lazy adapters for official pretrained property-prediction implementations."""

from __future__ import annotations

import ast
import csv
import importlib
import math
import os
import re
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, TypeVar, cast

from dnakit.backends import execute_bounded_command
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    SequenceError,
)

from .checkpoints import PredictionCheckpointInfo
from .enformer_benchmarks import get_enformer_benchmark_task, is_enformer_benchmark_task
from .models import (
    BiologicalSequence,
    BiologicalSequencePair,
    PredictionInput,
    PropertyPredictionConfig,
    VariantContext,
)
from .results import PredictionOutput

_T = TypeVar("_T")
_CANONICAL_BASES = frozenset("ACGT")
_ALPHAGENOME_LENGTHS = (16_384, 131_072, 524_288, 1_048_576)
_LUCA_SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,254}\Z")
_ENFORMER_TARGET_URLS = {
    "human": (
        "https://raw.githubusercontent.com/calico/basenji/0.5/"
        "manuscripts/cross2020/targets_human.txt"
    ),
    "mouse": (
        "https://raw.githubusercontent.com/calico/basenji/0.5/"
        "manuscripts/cross2020/targets_mouse.txt"
    ),
}
_ENFORMER_BENCHMARK_FORMAT = "enformer_full_finetune_best_valid_v2_epoch_validation"
_ENFORMER_BENCHMARK_MAX_CHECKPOINT_BYTES = 2_000_000_000
_ENFORMER_BENCHMARK_MAX_SEQUENCE_LENGTH = 196_608
_ENFORMER_EMBEDDING_DIM = 3_072
_ENFORMER_DOWNSAMPLE = 128


class PredictionBackend(Protocol):
    """Backend contract used by :func:`dnakit.predictions.predict_properties`."""

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        """Return one numeric output per input."""


def _require_module(name: str, model: str, requirement: str | None = None) -> Any:
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        dependency = requirement or name
        raise BackendUnavailableError(
            f"The {model} prediction backend is not installed.",
            code="MISSING_NEURAL_DEPENDENCY",
            context={"model": model, "dependency": dependency},
            hint=f"Install the official prediction dependency: {dependency}",
        ) from exc


def _prepend_source_path(value: str | os.PathLike[str] | None) -> None:
    if value is None:
        return
    root = Path(value).expanduser().resolve()
    for candidate in (root, root / "src"):
        if not candidate.is_dir():
            continue
        text = str(candidate)
        if text not in sys.path:
            sys.path.insert(0, text)


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
            "Invalid prediction device.",
            code="INVALID_PREDICTION_DEVICE",
            context={"device": requested},
        ) from exc
    if device.type == "cuda" and not bool(torch.cuda.is_available()):
        raise ConfigurationError(
            "Requested CUDA device is unavailable.",
            code="INVALID_PREDICTION_DEVICE",
            context={"device": requested},
        )
    if device.type == "mps" and (
        not hasattr(torch.backends, "mps") or not bool(torch.backends.mps.is_available())
    ):
        raise ConfigurationError(
            "Requested MPS device is unavailable.",
            code="INVALID_PREDICTION_DEVICE",
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


def _batches(
    values: Sequence[_T],
    *,
    batch_size: int,
    enabled: bool,
    description: str,
) -> Iterable[Sequence[_T]]:
    starts: Iterable[int] = range(0, len(values), batch_size)
    if enabled:
        from rich.progress import track

        starts = track(starts, description=description)
    for start in starts:
        yield values[start : start + batch_size]


def _dna(sequence: str, config: PropertyPredictionConfig, record_id: str) -> str:
    value = sequence.replace("U", "T")
    ambiguous = sorted(set(value) - _CANONICAL_BASES)
    if ambiguous and config.ambiguity_policy == "error":
        raise SequenceError(
            "The selected ambiguity policy rejects non-ACGT symbols.",
            code="MODEL_AMBIGUOUS_INPUT",
            context={"record_id": record_id, "symbols": tuple(ambiguous)},
            hint='Use ambiguity_policy="replace_with_n" to map IUPAC symbols to N.',
        )
    if ambiguous:
        value = "".join(base if base in _CANONICAL_BASES else "N" for base in value)
    if config.max_length is not None and len(value) > config.max_length:
        raise SequenceError(
            "Sequence length exceeds max_length; prediction inputs are never silently truncated.",
            code="PREDICTION_SEQUENCE_TOO_LONG",
            context={
                "record_id": record_id,
                "length": len(value),
                "max_length": config.max_length,
            },
        )
    return value


def _path(checkpoints: PredictionCheckpointInfo, index: int = 0) -> Path:
    try:
        return Path(checkpoints.paths[index]).expanduser().resolve()
    except IndexError as exc:
        raise BackendExecutionError(
            "The prediction backend did not receive its required checkpoint path.",
            code="MODEL_CHECKPOINT_NOT_FOUND",
        ) from exc


class _SegmentNTBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        self.torch = _require_module("torch", "SegmentNT", "torch")
        transformers = _require_module("transformers", "SegmentNT", "transformers")
        self.device = _resolve_torch_device(self.torch, config.device)
        checkpoint = _path(checkpoints)
        try:
            self.tokenizer = transformers.AutoTokenizer.from_pretrained(
                str(checkpoint),
                trust_remote_code=True,
                local_files_only=True,
            )
            self.model = transformers.AutoModel.from_pretrained(
                str(checkpoint),
                trust_remote_code=True,
                local_files_only=True,
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the official SegmentNT checkpoint.",
                code="MODEL_LOAD_FAILED",
                context={"model": "segmentnt"},
            ) from exc
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            self.model = self.model.to(dtype=dtype)
        self.model = self.model.to(self.device).eval()
        raw_features = getattr(self.model.config, "features", None)
        if not isinstance(raw_features, (list, tuple)) or len(raw_features) != 14:
            raise BackendExecutionError(
                "SegmentNT checkpoint does not declare its 14 genomic features.",
                code="MODEL_LOAD_FAILED",
            )
        self.features = tuple(str(value) for value in raw_features)
        self.maximum_bp = int(config.max_length or 30_000)
        if self.maximum_bp > 30_000:
            raise ConfigurationError(
                (
                    "The integrated human SegmentNT checkpoint is bounded to its "
                    "30-kb training length."
                ),
                code="PREDICTION_SEQUENCE_TOO_LONG",
                hint="Use 30000 or less for max_length.",
            )

    def _predict_batch(self, batch: Sequence[PredictionInput]) -> list[PredictionOutput]:
        sequences: list[str] = []
        lengths: list[int] = []
        for item in batch:
            if not isinstance(item, BiologicalSequence):
                raise AssertionError("SegmentNT received a non-sequence input.")
            sequence = _dna(item.sequence, self.config, item.id)
            if len(sequence) > self.maximum_bp:
                raise SequenceError(
                    "SegmentNT input exceeds 30,000 bp.",
                    code="PREDICTION_SEQUENCE_TOO_LONG",
                    context={"record_id": item.id, "length": len(sequence)},
                )
            sequences.append(sequence)
            lengths.append(len(sequence))
        try:
            unpadded = self.tokenizer.batch_encode_plus(
                sequences,
                add_special_tokens=True,
                padding=False,
            )["input_ids"]
            dna_tokens = max(len(row) - 1 for row in unpadded)
            padded_dna_tokens = int(math.ceil(dna_tokens / 4.0) * 4)
            token_length = padded_dna_tokens + 1
            encoded = self.tokenizer.batch_encode_plus(
                sequences,
                add_special_tokens=True,
                return_tensors="pt",
                padding="max_length",
                max_length=token_length,
                truncation=False,
            )
            input_ids = encoded["input_ids"].to(self.device)
            attention_mask = encoded.get("attention_mask")
            if attention_mask is None:
                attention_mask = input_ids != self.tokenizer.pad_token_id
            else:
                attention_mask = attention_mask.to(self.device)
            with self.torch.inference_mode():
                outputs = self.model(input_ids, attention_mask=attention_mask)
                logits = outputs.logits
                probabilities = self.torch.softmax(logits, dim=-1)[..., 1]
        except (ConfigurationError, SequenceError):
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "SegmentNT inference failed.",
                code="PROPERTY_PREDICTION_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc
        array = probabilities.detach().float().cpu().numpy()
        return [
            PredictionOutput(
                array[index, :length, :],
                self.features,
                {
                    "axes": ("nucleotide", "feature"),
                    "probability": "feature_present",
                    "sequence_length": length,
                    "checkpoint": "InstaDeepAI/segment_nt",
                },
            )
            for index, length in enumerate(lengths)
        ]

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        results: list[PredictionOutput] = []
        for batch in _batches(
            inputs,
            batch_size=self.config.batch_size,
            enabled=show_progress,
            description="Predicting SegmentNT annotations",
        ):
            results.extend(self._predict_batch(batch))
        return results


def _read_enformer_targets(checkpoint: Path, organism: str) -> tuple[str, ...]:
    path = checkpoint / f"targets_{organism}.txt"
    if not path.is_file():
        return ()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            return tuple(
                f"{row.get('identifier', '').strip()} {row.get('description', '').strip()}".strip()
                for row in rows
            )
    except (OSError, UnicodeDecodeError, csv.Error):
        return ()


class _EnformerBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        self.torch = _require_module("torch", "Enformer", "torch")
        enformer_pytorch = _require_module(
            "enformer_pytorch", "Enformer", "enformer-pytorch>=0.8.11"
        )
        self.device = _resolve_torch_device(self.torch, config.device)
        self.checkpoint = _path(checkpoints)
        try:
            self.model = enformer_pytorch.from_pretrained(str(self.checkpoint))
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Enformer checkpoint.",
                code="MODEL_LOAD_FAILED",
            ) from exc
        for module in self.model.modules():
            if hasattr(module, "use_tf_gamma"):
                module.use_tf_gamma = True
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            self.model = self.model.to(dtype=dtype)
        self.model = self.model.to(self.device).eval()
        self.input_length = int(config.max_length or 196_608)
        if not 114_688 <= self.input_length <= 393_216 or self.input_length % 128 != 0:
            raise ConfigurationError(
                "Enformer max_length must be a multiple of 128 in [114688, 393216].",
                code="INVALID_PREDICTION_CONFIG",
            )
        self.organism = "human" if config.task == "human_tracks" else "mouse"
        self.target_names = _read_enformer_targets(self.checkpoint, self.organism)
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

    def _pad(self, sequence: str, record_id: str) -> tuple[str, int, int]:
        if len(sequence) > self.input_length:
            raise SequenceError(
                "Enformer input exceeds the configured context length.",
                code="PREDICTION_SEQUENCE_TOO_LONG",
                context={"record_id": record_id, "length": len(sequence)},
            )
        missing = self.input_length - len(sequence)
        left = missing // 2
        right = missing - left
        return "N" * left + sequence + "N" * right, left, right

    def _predict_sequences(self, sequences: Sequence[str]) -> Any:
        encoded = self.torch.tensor(
            [list(sequence.encode("ascii")) for sequence in sequences],
            dtype=self.torch.long,
            device=self.device,
        )
        with self.torch.inference_mode():
            predictions = self.model(self.lookup[encoded], head=self.organism)
        return predictions.detach().float().cpu().numpy()

    def _names(self, count: int) -> tuple[str, ...]:
        return self.target_names if len(self.target_names) == count else ()

    def _metadata(
        self,
        *,
        original_length: int,
        left: int,
        right: int,
        shape: Sequence[int],
        variant: bool,
    ) -> dict[str, object]:
        return {
            "axes": (("allele", "bin", "track") if variant else ("bin", "track")),
            "allele_axis": ("reference", "alternate", "alternate_minus_reference")
            if variant
            else (),
            "organism": self.organism,
            "bin_size_bp": 128,
            "input_length": self.input_length,
            "original_length": original_length,
            "padding_left": left,
            "padding_right": right,
            "output_bins": int(shape[-2]),
            "track_metadata_url": _ENFORMER_TARGET_URLS[self.organism],
        }

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        import numpy as np

        results: list[PredictionOutput] = []
        for batch in _batches(
            inputs,
            batch_size=self.config.batch_size,
            enabled=show_progress,
            description="Predicting Enformer tracks",
        ):
            padded: list[str] = []
            padding: list[tuple[int, int, int]] = []
            variant_batch = isinstance(batch[0], VariantContext)
            for item in batch:
                if isinstance(item, BiologicalSequence):
                    sequence = _dna(item.sequence, self.config, item.id)
                    value, left, right = self._pad(sequence, item.id)
                    padded.append(value)
                    padding.append((len(sequence), left, right))
                elif isinstance(item, VariantContext):
                    reference = _dna(item.reference_sequence, self.config, item.id)
                    alternate = _dna(item.alternate_sequence, self.config, item.id)
                    ref_padded, left, right = self._pad(reference, item.id)
                    alt_padded, _, _ = self._pad(alternate, item.id)
                    padded.extend((ref_padded, alt_padded))
                    padding.append((len(reference), left, right))
                else:
                    raise AssertionError("Enformer received an unsupported input.")
            try:
                raw = self._predict_sequences(padded)
            except (ConfigurationError, SequenceError):
                raise
            except Exception as exc:
                raise BackendExecutionError(
                    "Enformer inference failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={"error_type": type(exc).__name__},
                ) from exc
            if variant_batch:
                for index, (length, left, right) in enumerate(padding):
                    reference = raw[index * 2]
                    alternate = raw[index * 2 + 1]
                    values = np.stack((reference, alternate, alternate - reference), axis=0)
                    item = batch[index]
                    assert isinstance(item, VariantContext)
                    assert item.variant_index is not None
                    metadata = self._metadata(
                        original_length=length,
                        left=left,
                        right=right,
                        shape=values.shape,
                        variant=True,
                    )
                    metadata.update(
                        {
                            "variant_index": item.variant_index,
                            "padded_variant_index": left + item.variant_index,
                            "reference_base": item.reference_base,
                            "alternate_base": item.alternate_base,
                        }
                    )
                    results.append(
                        PredictionOutput(
                            values,
                            self._names(int(values.shape[-1])),
                            metadata,
                        )
                    )
            else:
                for values, (length, left, right) in zip(raw, padding, strict=True):
                    results.append(
                        PredictionOutput(
                            values,
                            self._names(int(values.shape[-1])),
                            self._metadata(
                                original_length=length,
                                left=left,
                                right=right,
                                shape=values.shape,
                                variant=False,
                            ),
                        )
                    )
        return results


def _checkpoint_labels(payload: Mapping[str, object], num_classes: int) -> tuple[str, ...]:
    mapping = payload.get("label_mapping")
    if not isinstance(mapping, Mapping) or len(mapping) != num_classes:
        raise BackendExecutionError(
            "The Enformer task checkpoint has an invalid label mapping.",
            code="INVALID_MODEL_CHECKPOINT",
        )
    labels = [""] * num_classes
    for label, index in mapping.items():
        if (
            not isinstance(label, str)
            or not label
            or isinstance(index, bool)
            or not isinstance(index, int)
            or not 0 <= index < num_classes
            or labels[index]
        ):
            raise BackendExecutionError(
                "The Enformer task checkpoint has an invalid label mapping.",
                code="INVALID_MODEL_CHECKPOINT",
            )
        labels[index] = label
    if any(not label for label in labels):
        raise BackendExecutionError(
            "The Enformer task checkpoint label mapping is incomplete.",
            code="INVALID_MODEL_CHECKPOINT",
        )
    return tuple(labels)


class _EnformerBenchmarkBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        self.spec = get_enformer_benchmark_task(config.task)
        self.torch = _require_module("torch", "Enformer benchmark", "torch>=2.5")
        enformer_pytorch = _require_module(
            "enformer_pytorch",
            "Enformer benchmark",
            "enformer-pytorch>=0.8.11",
        )
        self.device = _resolve_torch_device(self.torch, config.device)
        self.checkpoint = _path(checkpoints)
        self.payload = self._load_checkpoint()
        self.num_classes = self._integer_metadata("num_classes")
        if self.num_classes != self.spec.num_classes:
            raise BackendExecutionError(
                "The checkpoint class count does not match the selected Enformer task.",
                code="INVALID_MODEL_CHECKPOINT",
                context={
                    "task": self.spec.name,
                    "expected": self.spec.num_classes,
                    "actual": self.num_classes,
                },
            )
        self.labels = _checkpoint_labels(self.payload, self.num_classes)
        self._validate_metadata()
        self.checkpoint_metadata = {
            "format": self.payload["format"],
            "saved_at_utc": self.payload.get("saved_at_utc", ""),
            "fold_id_one_based": self.payload.get("fold_id_one_based"),
        }
        self.model = self._create_model(enformer_pytorch)
        # On CUDA, the assigned memory-mapped CPU state is no longer needed after
        # ``to(device)``. Keep only compact provenance instead of retaining a second
        # roughly 0.9-GB checkpoint mapping for the lifetime of the backend.
        del self.payload

    def _load_checkpoint(self) -> Mapping[str, object]:
        try:
            size = self.checkpoint.stat().st_size
        except OSError as exc:
            raise BackendExecutionError(
                "Could not inspect the Enformer task checkpoint.",
                code="MODEL_CHECKPOINT_NOT_FOUND",
                context={"path": str(self.checkpoint)},
            ) from exc
        if not 1 <= size <= _ENFORMER_BENCHMARK_MAX_CHECKPOINT_BYTES:
            raise BackendExecutionError(
                "The Enformer task checkpoint exceeds the allowed size or is empty.",
                code="INVALID_MODEL_CHECKPOINT",
                context={"path": str(self.checkpoint), "bytes": size},
            )
        try:
            payload = self.torch.load(
                self.checkpoint,
                map_location="cpu",
                weights_only=True,
                mmap=True,
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not safely load the Enformer task checkpoint.",
                code="MODEL_LOAD_FAILED",
                context={"path": str(self.checkpoint), "error_type": type(exc).__name__},
            ) from exc
        if not isinstance(payload, Mapping):
            raise BackendExecutionError(
                "The Enformer task checkpoint root must be a mapping.",
                code="INVALID_MODEL_CHECKPOINT",
            )
        return cast(Mapping[str, object], payload)

    def _integer_metadata(self, name: str) -> int:
        value = self.payload.get(name)
        if isinstance(value, bool) or not isinstance(value, int):
            raise BackendExecutionError(
                "The Enformer task checkpoint metadata is incomplete.",
                code="INVALID_MODEL_CHECKPOINT",
                context={"field": name},
            )
        return value

    def _validate_metadata(self) -> None:
        checkpoint_format = self.payload.get("format")
        family = self.payload.get("family")
        task = self.payload.get("task")
        if checkpoint_format != _ENFORMER_BENCHMARK_FORMAT:
            raise BackendExecutionError(
                "The file is not a supported Enformer full-finetune checkpoint.",
                code="INVALID_MODEL_CHECKPOINT",
                context={"format": checkpoint_format},
            )
        if (
            family != self.spec.family
            or not isinstance(task, str)
            or (task.lower() != self.spec.checkpoint_task.lower())
        ):
            raise BackendExecutionError(
                "The checkpoint does not belong to the selected Enformer task.",
                code="MODEL_CHECKPOINT_TASK_MISMATCH",
                context={
                    "selected_task": self.spec.name,
                    "checkpoint_family": family,
                    "checkpoint_task": task,
                },
            )

    def _create_model(self, enformer_pytorch: Any) -> Any:
        torch = self.torch
        num_classes = self.num_classes

        class EnformerSequenceClassifier(torch.nn.Module):  # type: ignore[misc,name-defined]
            def __init__(self) -> None:
                super().__init__()
                self.backbone = enformer_pytorch.Enformer.from_hparams(
                    output_heads={},
                    target_length=-1,
                    use_tf_gamma=False,
                    use_checkpointing=False,
                )
                self.classifier = torch.nn.Linear(_ENFORMER_EMBEDDING_DIM, num_classes)

            def forward(self, one_hot: Any, lengths: Any) -> Any:
                embeddings = self.backbone(one_hot, return_only_embeddings=True)
                valid_bins = (lengths + _ENFORMER_DOWNSAMPLE - 1) // _ENFORMER_DOWNSAMPLE
                positions = torch.arange(embeddings.shape[1], device=embeddings.device)
                mask = positions.unsqueeze(0) < valid_bins.unsqueeze(1)
                pooled = (embeddings * mask.unsqueeze(-1)).sum(dim=1)
                pooled = pooled / valid_bins.unsqueeze(-1)
                return self.classifier(pooled)

        state = self.payload.get("model_state_dict")
        if not isinstance(state, Mapping):
            raise BackendExecutionError(
                "The Enformer task checkpoint has no model_state_dict.",
                code="INVALID_MODEL_CHECKPOINT",
            )
        classifier_weight = state.get("classifier.weight")
        classifier_bias = state.get("classifier.bias")
        if (
            not isinstance(classifier_weight, torch.Tensor)
            or tuple(classifier_weight.shape) != (num_classes, _ENFORMER_EMBEDDING_DIM)
            or not isinstance(classifier_bias, torch.Tensor)
            or tuple(classifier_bias.shape) != (num_classes,)
        ):
            raise BackendExecutionError(
                "The Enformer task checkpoint classification head is incompatible.",
                code="INVALID_MODEL_CHECKPOINT",
            )
        try:
            with torch.device("meta"):
                model = EnformerSequenceClassifier()
            model.load_state_dict(state, strict=True, assign=True)
            model = model.to(self.device).eval()
        except Exception as exc:
            raise BackendExecutionError(
                "The Enformer task checkpoint parameters are incompatible.",
                code="MODEL_LOAD_FAILED",
                context={"task": self.spec.name, "error_type": type(exc).__name__},
            ) from exc
        return model

    def _autocast(self) -> Any:
        requested = self.config.dtype
        if self.device.type == "cuda":
            if requested == "float32":
                return nullcontext()
            if requested == "auto":
                dtype = (
                    self.torch.bfloat16
                    if bool(self.torch.cuda.is_bf16_supported())
                    else self.torch.float16
                )
            else:
                dtype = _torch_dtype(self.torch, requested)
            return self.torch.autocast("cuda", dtype=dtype)
        if self.device.type == "cpu" and requested == "bfloat16":
            return self.torch.autocast("cpu", dtype=self.torch.bfloat16)
        return nullcontext()

    def _encode(self, sequences: Sequence[str]) -> tuple[Any, Any]:
        maximum = max(len(sequence) for sequence in sequences)
        encoded = self.torch.full(
            (len(sequences), maximum),
            ord("N"),
            dtype=self.torch.long,
            device=self.device,
        )
        for row, sequence in enumerate(sequences):
            encoded[row, : len(sequence)] = self.torch.tensor(
                list(sequence.encode("ascii")),
                dtype=self.torch.long,
                device=self.device,
            )
        lookup = self.torch.zeros(256, 4, dtype=self.torch.float32, device=self.device)
        for base, index in zip("ACGT", range(4), strict=True):
            lookup[ord(base), index] = 1.0
        lengths = self.torch.tensor(
            [len(sequence) for sequence in sequences],
            dtype=self.torch.long,
            device=self.device,
        )
        return lookup[encoded], lengths

    def _output(self, probabilities: Any, logits: Any, length: int) -> PredictionOutput:
        predicted_index = int(probabilities.argmax())
        metadata = {
            "axes": ("class",),
            "predicted_index": predicted_index,
            "predicted_label": self.labels[predicted_index],
            "confidence": float(probabilities[predicted_index]),
            "logits": tuple(float(value) for value in logits),
            "dataset_family": self.spec.family,
            "dataset_name": self.spec.dataset_name,
            "checkpoint_filename": self.spec.checkpoint_filename,
            "checkpoint_format": self.checkpoint_metadata["format"],
            "checkpoint_saved_at_utc": self.checkpoint_metadata["saved_at_utc"],
            "fold_id_one_based": self.checkpoint_metadata["fold_id_one_based"],
            "sequence_length": length,
            "valid_bins": (length + _ENFORMER_DOWNSAMPLE - 1) // _ENFORMER_DOWNSAMPLE,
            "pooling": "valid-bin mean",
            "checkpoint_is_fine_tuned": True,
        }
        return PredictionOutput(probabilities, self.labels, metadata)

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        prepared: list[tuple[int, str]] = []
        limit = int(self.config.max_length or _ENFORMER_BENCHMARK_MAX_SEQUENCE_LENGTH)
        for index, item in enumerate(inputs):
            if not isinstance(item, BiologicalSequence):
                raise AssertionError("Enformer benchmark received a non-sequence input.")
            sequence = _dna(item.sequence, self.config, item.id)
            if len(sequence) > limit:
                raise SequenceError(
                    "Enformer benchmark input exceeds the configured sequence limit.",
                    code="PREDICTION_SEQUENCE_TOO_LONG",
                    context={"record_id": item.id, "length": len(sequence), "max_length": limit},
                )
            prepared.append((index, sequence))

        ordered = sorted(prepared, key=lambda item: len(item[1]))
        outputs: list[PredictionOutput | None] = [None] * len(inputs)
        for batch in _batches(
            ordered,
            batch_size=self.config.batch_size,
            enabled=show_progress,
            description=f"Predicting {self.spec.display_name}",
        ):
            indices = [item[0] for item in batch]
            sequences = [item[1] for item in batch]
            one_hot, lengths = self._encode(sequences)
            try:
                with self.torch.inference_mode(), self._autocast():
                    logits = self.model(one_hot, lengths).float()
                    probabilities = logits.softmax(dim=-1)
            except (ConfigurationError, SequenceError):
                raise
            except Exception as exc:
                raise BackendExecutionError(
                    "Enformer benchmark inference failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={"task": self.spec.name, "error_type": type(exc).__name__},
                ) from exc
            logits_cpu = logits.detach().cpu()
            probabilities_cpu = probabilities.detach().cpu()
            for row, index in enumerate(indices):
                outputs[index] = self._output(
                    probabilities_cpu[row].numpy(),
                    logits_cpu[row].numpy(),
                    len(sequences[row]),
                )
        if any(output is None for output in outputs):
            raise BackendExecutionError(
                "Enformer benchmark returned incomplete predictions.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        return cast(Sequence[PredictionOutput], outputs)


def _alphagenome_device(jax: Any, requested: str) -> Any:
    devices = list(jax.devices())
    normalized = requested.strip().lower()
    if normalized == "auto":
        return next((device for device in devices if device.platform == "gpu"), devices[0])
    platform, separator, index_text = normalized.partition(":")
    if platform == "cuda":
        platform = "gpu"
    matching = [device for device in devices if device.platform == platform]
    try:
        index = 0 if not separator else int(index_text)
    except ValueError as exc:
        raise ConfigurationError(
            "Invalid JAX device index.", code="INVALID_PREDICTION_DEVICE"
        ) from exc
    if not matching or not 0 <= index < len(matching):
        raise ConfigurationError(
            "Requested JAX device is unavailable.",
            code="INVALID_PREDICTION_DEVICE",
            context={"device": requested},
        )
    return matching[index]


class _AlphaGenomeBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        _prepend_source_path(config.model_source_path)
        jax = _require_module(
            "jax",
            "AlphaGenome",
            "git+https://github.com/google-deepmind/alphagenome_research.git",
        )
        try:
            api_model = importlib.import_module("alphagenome.models.dna_model")
            output_module = importlib.import_module("alphagenome.models.dna_output")
            research_model = importlib.import_module("alphagenome_research.model.dna_model")
        except ImportError as exc:
            raise BackendUnavailableError(
                "AlphaGenome research code is required for local property prediction.",
                code="MISSING_NEURAL_DEPENDENCY",
                hint=(
                    "Install https://github.com/google-deepmind/alphagenome_research "
                    "or pass model_source_path to its checkout."
                ),
            ) from exc
        if config.dtype not in {"auto", "float32"}:
            raise ConfigurationError(
                "AlphaGenome controls mixed precision in its official JAX implementation.",
                code="INVALID_PREDICTION_DTYPE",
                hint='Use dtype="auto" or dtype="float32".',
            )
        self.device = _alphagenome_device(jax, config.device)
        settings = {
            api_model.Organism.HOMO_SAPIENS: research_model.OrganismSettings(),
            api_model.Organism.MUS_MUSCULUS: research_model.OrganismSettings(),
        }
        try:
            self.model = research_model.create(
                str(_path(checkpoints)),
                organism_settings=settings,
                device=self.device,
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the local AlphaGenome all-folds checkpoint.",
                code="MODEL_LOAD_FAILED",
                context={"error_type": type(exc).__name__},
            ) from exc
        self.organism = (
            api_model.Organism.HOMO_SAPIENS
            if config.organism == "human"
            else api_model.Organism.MUS_MUSCULUS
        )
        self.output_type = output_module.OutputType[config.task.upper()]
        if config.max_length is not None and config.max_length not in _ALPHAGENOME_LENGTHS:
            raise ConfigurationError(
                "AlphaGenome max_length must be 16384, 131072, 524288, or 1048576.",
                code="INVALID_PREDICTION_CONFIG",
            )

    def _pad(self, sequence: str, record_id: str) -> tuple[str, int, int]:
        if self.config.max_length is not None:
            target = self.config.max_length
        else:
            target = next((size for size in _ALPHAGENOME_LENGTHS if len(sequence) <= size), 0)
        if target == 0 or len(sequence) > target:
            raise SequenceError(
                "AlphaGenome input exceeds its supported context length.",
                code="PREDICTION_SEQUENCE_TOO_LONG",
                context={"record_id": record_id, "length": len(sequence)},
            )
        missing = target - len(sequence)
        left = missing // 2
        right = missing - left
        return "N" * left + sequence + "N" * right, left, right

    @staticmethod
    def _track_names(track: Any) -> tuple[str, ...]:
        try:
            names = tuple(str(value) for value in track.metadata["name"].tolist())
        except (AttributeError, KeyError, TypeError):
            return ()
        try:
            strands = tuple(str(value) for value in track.metadata["strand"].tolist())
        except (AttributeError, KeyError, TypeError):
            return names
        if len(names) != len(strands):
            return names
        return tuple(f"{name} [{strand}]" for name, strand in zip(names, strands, strict=True))

    def _axes(self, *, variant: bool) -> tuple[str, ...]:
        axes: tuple[str, ...]
        if self.config.task == "contact_maps":
            axes = ("position_1", "position_2", "track")
        elif self.config.task == "splice_junctions":
            axes = ("junction", "track")
        else:
            axes = ("position", "track")
        return ("allele", *axes) if variant else axes

    def _predict_track(self, sequence: str) -> tuple[Any, tuple[str, ...], dict[str, object]]:
        output = self.model.predict_sequence(
            sequence=sequence,
            organism=self.organism,
            requested_outputs=[self.output_type],
            ontology_terms=list(self.config.ontology_terms) or None,
        )
        track = output.get(self.output_type)
        if track is None or not hasattr(track, "values"):
            raise BackendExecutionError(
                "AlphaGenome did not return the requested output type.",
                code="INVALID_PREDICTION_OUTPUT",
                context={"task": self.config.task},
            )
        metadata: dict[str, object] = {
            "resolution_bp": int(getattr(track, "resolution", 1)),
        }
        junctions = getattr(track, "junctions", None)
        if junctions is not None:
            metadata["junctions"] = tuple(str(value) for value in junctions)
        return track.values, self._track_names(track), metadata

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        import numpy as np

        iterable: Iterable[PredictionInput] = inputs
        if show_progress:
            from rich.progress import track

            iterable = track(inputs, description="Predicting AlphaGenome tracks")
        results: list[PredictionOutput] = []
        for item in iterable:
            if isinstance(item, BiologicalSequence):
                sequence = _dna(item.sequence, self.config, item.id)
                padded, left, right = self._pad(sequence, item.id)
                try:
                    values, names, metadata = self._predict_track(padded)
                except (BackendExecutionError, ConfigurationError, SequenceError):
                    raise
                except Exception as exc:
                    raise BackendExecutionError(
                        "AlphaGenome inference failed.",
                        code="PROPERTY_PREDICTION_FAILED",
                        context={"error_type": type(exc).__name__},
                    ) from exc
                metadata.update(
                    {
                        "axes": self._axes(variant=False),
                        "organism": self.config.organism,
                        "ontology_terms": self.config.ontology_terms,
                        "input_length": len(padded),
                        "original_length": len(sequence),
                        "padding_left": left,
                        "padding_right": right,
                    }
                )
                results.append(PredictionOutput(values, names, metadata))
            elif isinstance(item, VariantContext):
                reference = _dna(item.reference_sequence, self.config, item.id)
                alternate = _dna(item.alternate_sequence, self.config, item.id)
                ref_padded, left, right = self._pad(reference, item.id)
                alt_padded, _, _ = self._pad(alternate, item.id)
                try:
                    ref_values, ref_names, metadata = self._predict_track(ref_padded)
                    alt_values, alt_names, _ = self._predict_track(alt_padded)
                except (BackendExecutionError, ConfigurationError, SequenceError):
                    raise
                except Exception as exc:
                    raise BackendExecutionError(
                        "AlphaGenome variant-track inference failed.",
                        code="PROPERTY_PREDICTION_FAILED",
                        context={"error_type": type(exc).__name__},
                    ) from exc
                ref_array = np.asarray(ref_values)
                alt_array = np.asarray(alt_values)
                if ref_array.shape != alt_array.shape or ref_names != alt_names:
                    raise BackendExecutionError(
                        "AlphaGenome reference and alternate track outputs do not align.",
                        code="INVALID_PREDICTION_OUTPUT",
                    )
                values = np.stack(
                    (ref_array, alt_array, alt_array - ref_array),
                    axis=0,
                )
                assert item.variant_index is not None
                metadata.update(
                    {
                        "axes": self._axes(variant=True),
                        "allele_axis": (
                            "reference",
                            "alternate",
                            "alternate_minus_reference",
                        ),
                        "organism": self.config.organism,
                        "ontology_terms": self.config.ontology_terms,
                        "input_length": len(ref_padded),
                        "original_length": len(reference),
                        "padding_left": left,
                        "padding_right": right,
                        "variant_index": item.variant_index,
                        "padded_variant_index": left + item.variant_index,
                    }
                )
                results.append(PredictionOutput(values, ref_names, metadata))
            else:
                raise AssertionError("AlphaGenome received an unsupported input.")
        return results


class _Evo2VariantBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        _prepend_source_path(config.model_source_path)
        evo2 = _require_module("evo2", "Evo 2", "evo2")
        checkpoint = _path(checkpoints) / "evo2_7b.pt"
        try:
            self.wrapper = evo2.Evo2(model_name="evo2_7b", local_path=str(checkpoint))
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Evo 2 7B checkpoint.",
                code="MODEL_LOAD_FAILED",
            ) from exc

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        results: list[PredictionOutput] = []
        for batch in _batches(
            inputs,
            batch_size=self.config.batch_size,
            enabled=show_progress,
            description="Scoring Evo 2 variants",
        ):
            variants = [item for item in batch if isinstance(item, VariantContext)]
            if len(variants) != len(batch):
                raise AssertionError("Evo 2 VEP received a non-variant input.")
            references = [_dna(item.reference_sequence, self.config, item.id) for item in variants]
            alternates = [_dna(item.alternate_sequence, self.config, item.id) for item in variants]
            try:
                ref_scores = self.wrapper.score_sequences(
                    references,
                    batch_size=len(references),
                )
                alt_scores = self.wrapper.score_sequences(
                    alternates,
                    batch_size=len(alternates),
                )
            except Exception as exc:
                raise BackendExecutionError(
                    "Evo 2 variant scoring failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={"error_type": type(exc).__name__},
                ) from exc
            for item, ref_score, alt_score in zip(variants, ref_scores, alt_scores, strict=True):
                reference = float(ref_score)
                alternate = float(alt_score)
                results.append(
                    PredictionOutput(
                        (reference, alternate, alternate - reference),
                        (
                            "reference_likelihood",
                            "alternate_likelihood",
                            "alternate_minus_reference",
                        ),
                        {
                            "score_type": "zero_shot_sequence_likelihood",
                            "reduction": "mean_log_likelihood",
                            "calibrated_probability": False,
                            "sequence_length": len(item.reference_sequence),
                            "variant_index": item.variant_index,
                            "reference_base": item.reference_base,
                            "alternate_base": item.alternate_base,
                        },
                    )
                )
        return results


class _Evo2ExonBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        _prepend_source_path(config.model_source_path)
        self.torch = _require_module("torch", "Evo 2 exon classifier", "torch")
        evo2 = _require_module("evo2", "Evo 2 exon classifier", "evo2")
        transformers = _require_module("transformers", "Evo 2 exon classifier", "transformers")
        self.device = _resolve_torch_device(self.torch, config.device)
        base = _path(checkpoints, 0) / "evo2_7b_base.pt"
        head = _path(checkpoints, 1)
        try:
            self.wrapper = evo2.Evo2(
                model_name="evo2_7b_base",
                local_path=str(base),
            )
            self.wrapper.model = self.wrapper.model.to(self.device).eval()
            self.classifier = (
                transformers.AutoModel.from_pretrained(
                    str(head),
                    trust_remote_code=True,
                    local_files_only=True,
                )
                .to(self.device)
                .eval()
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the Evo 2 exon-classifier checkpoints.",
                code="MODEL_LOAD_FAILED",
            ) from exc
        self.layer = "blocks.26"

    def _embedding(self, sequence: str) -> Any:
        input_ids = self.torch.tensor(
            self.wrapper.tokenizer.tokenize(sequence),
            dtype=self.torch.int,
            device=self.device,
        ).unsqueeze(0)
        with self.torch.inference_mode():
            _, embeddings = self.wrapper(
                input_ids,
                return_embeddings=True,
                layer_names=[self.layer],
            )
        return embeddings[self.layer][0, -1, :].float()

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        iterable: Iterable[PredictionInput] = inputs
        if show_progress:
            from rich.progress import track

            iterable = track(inputs, description="Predicting Evo 2 exon probability")
        results: list[PredictionOutput] = []
        for item in iterable:
            if not isinstance(item, BiologicalSequencePair):
                raise AssertionError("Evo 2 exon classifier received a non-pair input.")
            forward = _dna(item.first.sequence, self.config, item.first.id)
            reverse = _dna(item.second.sequence, self.config, item.second.id)
            try:
                combined = self.torch.cat(
                    (self._embedding(forward), self._embedding(reverse)), dim=-1
                ).reshape(1, 1, -1)
                with self.torch.inference_mode():
                    probability = float(self.classifier(combined)["logits"].item())
            except Exception as exc:
                raise BackendExecutionError(
                    "Evo 2 exon classification failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={"error_type": type(exc).__name__},
                ) from exc
            results.append(
                PredictionOutput(
                    (probability,),
                    ("exon_probability",),
                    {
                        "threshold": self.config.threshold,
                        "predicted_label": "exon"
                        if probability >= self.config.threshold
                        else "non_exon",
                        "context_order": ("forward", "reverse"),
                    },
                )
            )
        return results


class _GeneratorVariantBackend:
    def __init__(
        self,
        checkpoints: PredictionCheckpointInfo,
        config: PropertyPredictionConfig,
    ) -> None:
        self.config = config
        self.torch = _require_module("torch", "GENERator", "torch")
        transformers = _require_module("transformers", "GENERator", "transformers")
        self.device = _resolve_torch_device(self.torch, config.device)
        checkpoint = _path(checkpoints)
        load_kwargs: dict[str, object] = {
            "trust_remote_code": True,
            "local_files_only": True,
        }
        dtype = _torch_dtype(self.torch, config.dtype)
        if dtype is not None:
            load_kwargs["torch_dtype"] = dtype
        try:
            self.tokenizer = transformers.AutoTokenizer.from_pretrained(
                str(checkpoint),
                trust_remote_code=True,
                local_files_only=True,
            )
            self.tokenizer.padding_side = "left"
            self.model = (
                transformers.AutoModelForCausalLM.from_pretrained(str(checkpoint), **load_kwargs)
                .to(self.device)
                .eval()
            )
        except Exception as exc:
            raise BackendExecutionError(
                "Could not load the GENERator variant-effect checkpoint.",
                code="MODEL_LOAD_FAILED",
            ) from exc
        sorted_vocabulary = sorted(
            ((token_id, token) for token, token_id in self.tokenizer.get_vocab().items()),
            key=lambda item: item[0],
        )
        indices: dict[str, list[int]] = {}
        for token_id, token in sorted_vocabulary:
            if isinstance(token, str) and token:
                indices.setdefault(token[0], []).append(int(token_id))
        self.allele_indices = indices
        self.context_length = int(config.max_length or 8_192)

    def _prefix(self, item: VariantContext) -> str:
        assert item.variant_index is not None
        start = max(0, item.variant_index - self.context_length)
        sequence = _dna(item.reference_sequence[start : item.variant_index], self.config, item.id)
        sequence = sequence.lstrip("N")
        remainder = len(sequence) % 6
        return sequence[remainder:] if remainder else sequence

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        results: list[PredictionOutput] = []
        for batch in _batches(
            inputs,
            batch_size=self.config.batch_size,
            enabled=show_progress,
            description="Scoring GENERator variants",
        ):
            variants = [item for item in batch if isinstance(item, VariantContext)]
            if len(variants) != len(batch):
                raise AssertionError("GENERator received a non-variant input.")
            prefixes = ["<s>" + self._prefix(item) for item in variants]
            try:
                encoded = self.tokenizer(
                    prefixes,
                    add_special_tokens=False,
                    return_tensors="pt",
                    padding=True,
                )
                encoded = {key: value.to(self.device) for key, value in encoded.items()}
                with self.torch.inference_mode():
                    logits = self.model(**encoded).logits[:, -1, :]
                    probabilities = self.torch.softmax(logits, dim=-1).float().cpu()
            except Exception as exc:
                raise BackendExecutionError(
                    "GENERator variant scoring failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={"error_type": type(exc).__name__},
                ) from exc
            for row, item, prefix in zip(probabilities, variants, prefixes, strict=True):
                p_ref = float(row[self.allele_indices.get(item.reference_base, [])].sum().item())
                p_alt = float(row[self.allele_indices.get(item.alternate_base, [])].sum().item())
                score = math.log(max(p_ref, 1e-30) / (p_alt + 1e-10))
                results.append(
                    PredictionOutput(
                        (p_ref, p_alt, score),
                        ("reference_probability", "alternate_probability", "log_ratio_score"),
                        {
                            "score_formula": "log(p_ref / (p_alt + 1e-10))",
                            "calibrated_probability": False,
                            "context_bases": len(prefix) - 3,
                            "variant_index": item.variant_index,
                            "reference_base": item.reference_base,
                            "alternate_base": item.alternate_base,
                        },
                    )
                )
        return results


@dataclass(frozen=True, slots=True)
class _LucaPreset:
    dataset_name: str
    dataset_type: str
    task_type: str
    model_type: str
    input_mode: str
    time_str: str
    step: str
    sequence_types: tuple[str, ...]


_LUCA_PRESETS = {
    "central_dogma": _LucaPreset(
        "CentralDogma",
        "gene_protein",
        "binary_class",
        "lucappi2",
        "pair",
        "20240406173806",
        "64000",
        ("gene", "protein"),
    ),
    "genustax": _LucaPreset(
        "GenusTax",
        "gene",
        "multi_class",
        "luca_base",
        "single",
        "20240412100337",
        "24500",
        ("gene",),
    ),
    "influenza_antigenicity": _LucaPreset(
        "InfA",
        "gene_gene",
        "binary_class",
        "lucappi",
        "pair",
        "20240214105653",
        "9603",
        ("gene", "gene"),
    ),
    "ncrna_family": _LucaPreset(
        "ncRNAFam",
        "gene",
        "multi_class",
        "luca_base",
        "single",
        "20240414155526",
        "1958484",
        ("gene",),
    ),
    "ncrna_protein_interaction": _LucaPreset(
        "ncRPI",
        "gene_protein",
        "binary_class",
        "lucappi2",
        "pair",
        "20240404105148",
        "716380",
        ("gene", "protein"),
    ),
    "protein_interaction": _LucaPreset(
        "PPI",
        "protein",
        "binary_class",
        "lucappi",
        "pair",
        "20240216205421",
        "52304",
        ("protein", "protein"),
    ),
    "protein_location": _LucaPreset(
        "ProtLoc",
        "protein",
        "multi_class",
        "luca_base",
        "single",
        "20240412140824",
        "466005",
        ("protein",),
    ),
    "protein_stability": _LucaPreset(
        "ProtStab",
        "protein",
        "regression",
        "luca_base",
        "single",
        "20240404104215",
        "70371",
        ("protein",),
    ),
    "speciestax": _LucaPreset(
        "SpeciesTax",
        "gene",
        "multi_class",
        "luca_base",
        "single",
        "20240411144916",
        "24000",
        ("gene",),
    ),
    "supktax": _LucaPreset(
        "SupKTax",
        "gene",
        "multi_class",
        "luca_base",
        "single",
        "20240212202328",
        "37000",
        ("gene",),
    ),
}

_LUCA_SEQUENCE_TYPES = {"gene": "gene", "protein": "prot"}


class _LucaOneTasksBackend:
    def __init__(self, config: PropertyPredictionConfig) -> None:
        self.config = config
        if config.model_source_path is None:
            raise BackendUnavailableError(
                "LucaOneTasks prediction requires an explicit official source checkout.",
                code="MISSING_NEURAL_DEPENDENCY",
                hint=(
                    "Clone https://github.com/LucaOne/LucaOneTasks, install its official "
                    "environment, and pass model_source_path."
                ),
            )
        self.root = Path(config.model_source_path).expanduser().resolve()
        self.source = self.root / "src"
        self.script = self.source / "predict_v1.py"
        if not self.script.is_file():
            raise BackendUnavailableError(
                "model_source_path is not a LucaOneTasks checkout.",
                code="MISSING_NEURAL_DEPENDENCY",
                context={"path": str(self.root)},
            )
        self.preset = _LUCA_PRESETS[config.task]

    @staticmethod
    def _safe_id(value: str) -> str:
        if _LUCA_SAFE_ID.fullmatch(value) is None or ".." in value:
            raise SequenceError(
                "LucaOneTasks sequence IDs must use only letters, digits, '.', '_', or '-'.",
                code="INVALID_PREDICTION_INPUT",
                context={"record_id": value},
            )
        return value

    def _gpu_id(self) -> int:
        device = self.config.device.strip().lower()
        if device == "cpu":
            return -1
        if device == "auto":
            torch = _require_module("torch", "LucaOneTasks", "its official environment")
            return 0 if bool(torch.cuda.is_available()) else -1
        if device == "cuda":
            return 0
        if device.startswith("cuda:"):
            try:
                return int(device.split(":", 1)[1])
            except ValueError as exc:
                raise ConfigurationError(
                    "Invalid LucaOneTasks CUDA device.",
                    code="INVALID_PREDICTION_DEVICE",
                ) from exc
        raise ConfigurationError(
            "LucaOneTasks supports device='auto', 'cpu', 'cuda', or 'cuda:N'.",
            code="INVALID_PREDICTION_DEVICE",
        )

    def _write_input(self, path: Path, inputs: Sequence[PredictionInput]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            if self.preset.input_mode == "single":
                writer.writerow(("seq_id", "seq_type", "seq"))
                for item in inputs:
                    if not isinstance(item, BiologicalSequence):
                        raise AssertionError("LucaOne single task received another input kind.")
                    if item.sequence_type != self.preset.sequence_types[0]:
                        raise SequenceError(
                            "Sequence type does not match the selected LucaOneTasks head.",
                            code="PREDICTION_INPUT_KIND_MISMATCH",
                            context={
                                "task": self.config.task,
                                "expected": self.preset.sequence_types,
                                "actual": item.sequence_type,
                            },
                        )
                    writer.writerow(
                        (
                            self._safe_id(item.id),
                            _LUCA_SEQUENCE_TYPES[item.sequence_type],
                            item.sequence,
                        )
                    )
            else:
                writer.writerow(
                    ("seq_id_a", "seq_id_b", "seq_type_a", "seq_type_b", "seq_a", "seq_b")
                )
                for item in inputs:
                    if not isinstance(item, BiologicalSequencePair):
                        raise AssertionError("LucaOne pair task received another input kind.")
                    actual = (item.first.sequence_type, item.second.sequence_type)
                    if actual != self.preset.sequence_types:
                        raise SequenceError(
                            "Sequence-pair types do not match the selected LucaOneTasks head.",
                            code="PREDICTION_INPUT_KIND_MISMATCH",
                            context={
                                "task": self.config.task,
                                "expected": self.preset.sequence_types,
                                "actual": actual,
                            },
                        )
                    writer.writerow(
                        (
                            self._safe_id(item.first.id),
                            self._safe_id(item.second.id),
                            _LUCA_SEQUENCE_TYPES[item.first.sequence_type],
                            _LUCA_SEQUENCE_TYPES[item.second.sequence_type],
                            item.first.sequence,
                            item.second.sequence,
                        )
                    )

    @staticmethod
    def _literal_list(value: str, field_name: str) -> tuple[object, ...]:
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError) as exc:
            raise BackendExecutionError(
                f"LucaOneTasks returned invalid {field_name}.",
                code="INVALID_PREDICTION_OUTPUT",
            ) from exc
        if not isinstance(parsed, (list, tuple)) or len(parsed) > 100:
            raise BackendExecutionError(
                f"LucaOneTasks returned invalid {field_name}.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        return tuple(parsed)

    def _read_output(self, path: Path) -> list[PredictionOutput]:
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
        except (OSError, UnicodeDecodeError, csv.Error) as exc:
            raise BackendExecutionError(
                "Could not read LucaOneTasks prediction output.",
                code="INVALID_PREDICTION_OUTPUT",
            ) from exc
        results: list[PredictionOutput] = []
        for row in rows:
            if self.preset.task_type == "multi_class" and self.config.top_k > 1:
                probabilities = tuple(
                    float(cast(Any, value))
                    for value in self._literal_list(
                        row[f"top{self.config.top_k}_probs"], "top-k probabilities"
                    )
                )
                labels = tuple(
                    str(value)
                    for value in self._literal_list(
                        row[f"top{self.config.top_k}_labels"], "top-k labels"
                    )
                )
                results.append(
                    PredictionOutput(
                        probabilities,
                        labels,
                        {
                            "predicted_label": row["top1_label"],
                            "predicted_probability": float(row["top1_prob"]),
                            "top_k": self.config.top_k,
                        },
                    )
                )
            else:
                name = "prediction" if self.preset.task_type == "regression" else "probability"
                results.append(
                    PredictionOutput(
                        (float(row["prob"]),),
                        (name,),
                        {
                            "predicted_label": row["label"],
                            "threshold": self.config.threshold
                            if self.preset.task_type == "binary_class"
                            else None,
                        },
                    )
                )
        return results

    def predict(
        self,
        inputs: Sequence[PredictionInput],
        *,
        show_progress: bool,
    ) -> Sequence[PredictionOutput]:
        with tempfile.TemporaryDirectory(prefix=".dnakit-lucaone-", dir=self.source) as raw:
            workspace = Path(raw)
            input_path = workspace / "input.csv"
            output_path = workspace / "output.csv"
            self._write_input(input_path, inputs)
            arguments = [
                str(self.script),
                "--input_file",
                str(input_path),
                "--llm_truncation_seq_length",
                str(self.config.max_length or 4096),
                "--model_path",
                str(self.root),
                "--save_path",
                str(output_path),
                "--dataset_name",
                self.preset.dataset_name,
                "--dataset_type",
                self.preset.dataset_type,
                "--task_type",
                self.preset.task_type,
                "--task_level_type",
                "seq_level",
                "--model_type",
                self.preset.model_type,
                "--input_type",
                "matrix",
                "--input_mode",
                self.preset.input_mode,
                "--time_str",
                self.preset.time_str,
                "--step",
                self.preset.step,
                "--threshold",
                str(self.config.threshold),
                "--print_per_num",
                str(self.config.batch_size),
                "--gpu_id",
                str(self._gpu_id()),
            ]
            if self.preset.task_type == "multi_class":
                arguments.extend(("--topk", str(self.config.top_k)))
            progress_ui: Any | None = None
            if show_progress:
                from rich.progress import Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

                progress_ui = Progress(
                    SpinnerColumn(),
                    TextColumn("Running official LucaOneTasks predictor"),
                    TimeElapsedColumn(),
                )
                progress_ui.start()
                progress_ui.add_task("lucaone", total=None)
            try:
                result = execute_bounded_command(
                    sys.executable,
                    arguments,
                    backend_id="lucaone-tasks",
                    cwd=self.source,
                    timeout_seconds=self.config.timeout_seconds,
                    max_output_bytes=self.config.max_backend_output_bytes,
                    monitored_output_paths=(output_path,),
                    max_monitored_output_bytes=100_000_000,
                )
            finally:
                if progress_ui is not None:
                    progress_ui.stop()
            if result.return_code != 0:
                raise BackendExecutionError(
                    "The official LucaOneTasks predictor failed.",
                    code="PROPERTY_PREDICTION_FAILED",
                    context={
                        "return_code": result.return_code,
                        "output_excerpt": result.output[-2_000:],
                    },
                )
            outputs = self._read_output(output_path)
        if len(outputs) != len(inputs):
            raise BackendExecutionError(
                "LucaOneTasks returned the wrong number of predictions.",
                code="INVALID_PREDICTION_OUTPUT",
                context={"expected": len(inputs), "actual": len(outputs)},
            )
        return outputs


def create_prediction_backend(
    checkpoints: PredictionCheckpointInfo,
    config: PropertyPredictionConfig,
) -> PredictionBackend:
    """Create the lazy official adapter selected by ``config``."""

    if config.model == "segmentnt":
        return _SegmentNTBackend(checkpoints, config)
    if config.model == "enformer" and is_enformer_benchmark_task(config.task):
        return _EnformerBenchmarkBackend(checkpoints, config)
    if config.model == "enformer":
        return _EnformerBackend(checkpoints, config)
    if config.model == "alphagenome":
        return _AlphaGenomeBackend(checkpoints, config)
    if config.model == "generator":
        return _GeneratorVariantBackend(checkpoints, config)
    if config.model == "evo2" and config.task == "variant_effect":
        return _Evo2VariantBackend(checkpoints, config)
    if config.model == "evo2" and config.task == "exon_probability":
        return _Evo2ExonBackend(checkpoints, config)
    if config.model == "lucaone":
        return _LucaOneTasksBackend(config)
    raise AssertionError(f"Unsupported direct-prediction backend: {config.model}:{config.task}")


__all__ = ["PredictionBackend", "create_prediction_backend"]
