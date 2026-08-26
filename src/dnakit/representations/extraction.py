"""Public sequence-level representation extraction API."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from dnakit.core import DNARecord
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    SequenceError,
    UnsupportedGapOperationError,
)

from .backends import RepresentationBackend, create_representation_backend
from .checkpoints import ensure_model_checkpoint
from .models import RepresentationConfig
from .results import RepresentationResult

_CANONICAL_BASES = frozenset("ACGT")


def _materialize_records(
    records: Iterable[DNARecord],
    *,
    max_records: int,
) -> tuple[DNARecord, ...]:
    iterator = iter(records)
    materialized: list[DNARecord] = []
    for index, record in enumerate(iterator):
        if index >= max_records:
            raise ConfigurationError(
                "Representation extraction exceeded max_records.",
                code="REPRESENTATION_RECORD_LIMIT",
                context={"max_records": max_records},
            )
        if not isinstance(record, DNARecord):
            raise SequenceError(
                "Representation extraction accepts only DNARecord objects.",
                code="INVALID_RECORD_SEQUENCE",
                context={"index": index, "type": type(record).__name__},
            )
        materialized.append(record)
    if not materialized:
        raise ConfigurationError(
            "At least one DNARecord is required for representation extraction.",
            code="EMPTY_REPRESENTATION_INPUT",
        )
    return tuple(materialized)


def _model_sequence(record: DNARecord, *, ambiguity_policy: str) -> str:
    sequence = record.sequence
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Neural representation extraction does not accept explicit Gap objects.",
            code="MODEL_GAPPED_INPUT",
            context={"record_id": record.id},
            hint="Resolve/remove gaps before extracting model representations.",
        )
    symbols = sequence.symbols
    if not symbols:
        raise SequenceError(
            "Neural representation extraction does not accept empty sequences.",
            code="MODEL_EMPTY_INPUT",
            context={"record_id": record.id},
        )
    ambiguous = sorted(set(symbols) - _CANONICAL_BASES)
    if ambiguous and ambiguity_policy == "error":
        raise SequenceError(
            "The selected ambiguity policy rejects non-ACGT symbols.",
            code="MODEL_AMBIGUOUS_INPUT",
            context={"record_id": record.id, "symbols": tuple(ambiguous)},
            hint='Use ambiguity_policy="replace_with_n" to map IUPAC symbols to N.',
        )
    if not ambiguous:
        return symbols
    return "".join(base if base in _CANONICAL_BASES else "N" for base in symbols)


def _require_numpy() -> object:
    try:
        import numpy as np
    except ImportError as exc:
        raise BackendUnavailableError(
            "Representation extraction requires NumPy.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
        ) from exc
    return np


def extract_representations(
    records: Iterable[DNARecord],
    *,
    config: RepresentationConfig | None = None,
    backend: RepresentationBackend | None = None,
) -> RepresentationResult:
    """Extract one fixed-size representation per DNA record.

    With the standard backend, a missing checkpoint is downloaded once to
    ``./ckpt/<model>``.  Passing ``backend`` supports controlled custom models and
    unit tests without downloading or importing a foundation-model runtime.
    """

    resolved = RepresentationConfig() if config is None else config
    if not isinstance(resolved, RepresentationConfig):
        raise ConfigurationError(
            "config must be RepresentationConfig or None.",
            code="INVALID_REPRESENTATION_CONFIG",
        )
    materialized = _materialize_records(records, max_records=resolved.max_records)
    sequences = tuple(
        _model_sequence(record, ambiguity_policy=resolved.ambiguity_policy)
        for record in materialized
    )

    checkpoint_path: str | None = None
    selected_backend = backend
    if selected_backend is None:
        if resolved.model_spec.trust_remote_code and not resolved.allow_remote_code:
            raise ConfigurationError(
                "The selected model requires loading checkpoint-provided Python code.",
                code="MODEL_REMOTE_CODE_NOT_ALLOWED",
                context={"model": resolved.model},
                hint=(
                    "Review the official repository and checkpoint, then set "
                    "allow_remote_code=True to opt in."
                ),
            )
        if resolved.checkpoint_path is None:
            checkpoint = ensure_model_checkpoint(
                resolved.model,
                checkpoint_dir=resolved.checkpoint_dir,
                hf_token=resolved.hf_token,
                show_progress=resolved.show_progress,
            )
            checkpoint_path = checkpoint.path
        else:
            path = Path(resolved.checkpoint_path).expanduser().resolve()
            if not path.is_dir():
                raise ConfigurationError(
                    "checkpoint_path must be an existing checkpoint directory.",
                    code="MODEL_CHECKPOINT_NOT_FOUND",
                    context={"path": str(path)},
                )
            checkpoint_path = str(path)
        try:
            selected_backend = create_representation_backend(
                resolved.model_spec,
                checkpoint_path,
                resolved,
            )
        except (
            BackendExecutionError,
            BackendUnavailableError,
            ConfigurationError,
            SequenceError,
        ):
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "Could not initialize the representation backend.",
                code="MODEL_LOAD_FAILED",
                context={"model": resolved.model, "error_type": type(exc).__name__},
            ) from exc
    elif not hasattr(selected_backend, "extract"):
        raise ConfigurationError(
            "backend must provide extract(sequences, show_progress=...).",
            code="INVALID_REPRESENTATION_BACKEND",
        )

    try:
        raw_matrix = selected_backend.extract(sequences, show_progress=resolved.show_progress)
    except (BackendExecutionError, BackendUnavailableError, ConfigurationError, SequenceError):
        raise
    except Exception as exc:
        raise BackendExecutionError(
            "The representation backend failed while encoding DNA sequences.",
            code="REPRESENTATION_EXTRACTION_FAILED",
            context={"model": resolved.model, "error_type": type(exc).__name__},
        ) from exc

    np = _require_numpy()
    try:
        matrix = np.asarray(raw_matrix, dtype=np.float32)  # type: ignore[attr-defined]
    except (TypeError, ValueError) as exc:
        raise BackendExecutionError(
            "The representation backend returned a non-numeric matrix.",
            code="INVALID_REPRESENTATION_MATRIX",
        ) from exc
    if matrix.ndim != 2 or matrix.shape[0] != len(materialized) or matrix.shape[1] < 1:
        raise BackendExecutionError(
            "The representation backend returned an invalid matrix shape.",
            code="INVALID_REPRESENTATION_MATRIX",
            context={
                "shape": tuple(int(value) for value in matrix.shape),
                "expected_rows": len(materialized),
            },
        )
    if not bool(np.isfinite(matrix).all()):  # type: ignore[attr-defined]
        raise BackendExecutionError(
            "The representation backend returned non-finite values.",
            code="INVALID_REPRESENTATION_MATRIX",
        )
    return RepresentationResult(
        tuple(record.id for record in materialized),
        matrix,
        resolved.model,
        checkpoint_path,
        resolved.pooling,
        int(matrix.shape[1]),
        len(materialized),
        resolved.ambiguity_policy,
    )


__all__ = ["extract_representations"]
