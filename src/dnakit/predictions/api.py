"""Public API for official pretrained property-prediction tasks."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from dnakit.core import DNARecord
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    SequenceError,
    UnsupportedGapOperationError,
)

from .backends import PredictionBackend, create_prediction_backend
from .checkpoints import PredictionCheckpointInfo, ensure_prediction_checkpoint
from .models import (
    BiologicalSequence,
    BiologicalSequencePair,
    PredictionInput,
    PredictionInputKind,
    PropertyPredictionConfig,
    VariantContext,
)
from .results import PredictionOutput, PredictionRecord, PropertyPredictionResult


def _record_input(record: DNARecord) -> BiologicalSequence:
    if record.sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Deep-learning property prediction does not accept explicit Gap objects.",
            code="MODEL_GAPPED_INPUT",
            context={"record_id": record.id},
            hint="Resolve or remove gaps before property prediction.",
        )
    return BiologicalSequence(record.id, record.sequence.symbols, "gene")


def _materialize_inputs(
    inputs: Iterable[PredictionInput | DNARecord],
    *,
    max_records: int,
) -> tuple[PredictionInput, ...]:
    raw_inputs: object = inputs
    if isinstance(raw_inputs, (str, bytes)):
        raise SequenceError(
            "Prediction inputs must be an iterable of typed inputs, not raw text.",
            code="INVALID_PREDICTION_INPUT",
        )
    try:
        iterator = iter(inputs)
    except TypeError as exc:
        raise SequenceError(
            "Prediction inputs must be iterable.", code="INVALID_PREDICTION_INPUT"
        ) from exc
    materialized: list[PredictionInput] = []
    for index, item in enumerate(iterator):
        if index >= max_records:
            raise ConfigurationError(
                "Property prediction exceeded max_records.",
                code="PREDICTION_RECORD_LIMIT",
                context={"max_records": max_records},
            )
        if isinstance(item, DNARecord):
            materialized.append(_record_input(item))
        elif isinstance(item, (BiologicalSequence, BiologicalSequencePair, VariantContext)):
            materialized.append(item)
        else:
            raise SequenceError(
                "Unsupported property-prediction input object.",
                code="INVALID_PREDICTION_INPUT",
                context={"index": index, "type": type(item).__name__},
            )
    if not materialized:
        raise ConfigurationError(
            "At least one property-prediction input is required.",
            code="EMPTY_PREDICTION_INPUT",
        )
    ids = tuple(_input_id(item) for item in materialized)
    if len(set(ids)) != len(ids):
        raise ConfigurationError(
            "Property-prediction input IDs must be unique.",
            code="DUPLICATE_PREDICTION_ID",
        )
    return tuple(materialized)


def _input_id(item: PredictionInput) -> str:
    return item.id


def _input_kind(item: PredictionInput) -> PredictionInputKind:
    if isinstance(item, BiologicalSequence):
        return "sequence"
    if isinstance(item, BiologicalSequencePair):
        return "pair"
    return "variant"


def _validate_input_contract(
    inputs: Sequence[PredictionInput],
    config: PropertyPredictionConfig,
) -> PredictionInputKind:
    kinds = {_input_kind(item) for item in inputs}
    if len(kinds) != 1:
        raise ConfigurationError(
            "One prediction call cannot mix sequence, pair, and variant inputs.",
            code="PREDICTION_INPUT_KIND_MISMATCH",
            context={"input_kinds": tuple(sorted(kinds))},
        )
    kind = next(iter(kinds))
    if kind not in config.task_spec.input_kinds:
        raise ConfigurationError(
            "Input kind does not match the selected pretrained task.",
            code="PREDICTION_INPUT_KIND_MISMATCH",
            context={
                "model": config.model,
                "task": config.task,
                "expected": config.task_spec.input_kinds,
                "actual": kind,
            },
        )
    if config.model != "lucaone":
        for item in inputs:
            if isinstance(item, BiologicalSequence) and item.sequence_type != "gene":
                raise ConfigurationError(
                    "The selected DNA model accepts gene sequences, not proteins.",
                    code="PREDICTION_INPUT_KIND_MISMATCH",
                )
            if isinstance(item, BiologicalSequencePair) and (
                item.first.sequence_type != "gene" or item.second.sequence_type != "gene"
            ):
                raise ConfigurationError(
                    "The selected DNA model requires two gene-sequence contexts.",
                    code="PREDICTION_INPUT_KIND_MISMATCH",
                )
    return kind


def predict_properties(
    inputs: Iterable[PredictionInput | DNARecord],
    *,
    config: PropertyPredictionConfig,
    backend: PredictionBackend | None = None,
) -> PropertyPredictionResult:
    """Run one official pretrained property task without fitting or fine-tuning.

    ``DNARecord`` is accepted for single-gene tasks. Pairwise and variant tasks
    use :class:`BiologicalSequencePair` and :class:`VariantContext`, respectively.
    Passing ``backend`` provides a controlled custom adapter and bypasses model
    download and checkpoint-provided remote code.
    """

    if not isinstance(config, PropertyPredictionConfig):
        raise ConfigurationError(
            "config must be PropertyPredictionConfig.",
            code="INVALID_PREDICTION_CONFIG",
        )
    materialized = _materialize_inputs(inputs, max_records=config.max_records)
    input_kind = _validate_input_contract(materialized, config)

    checkpoints = PredictionCheckpointInfo((), False, ())
    selected_backend = backend
    if selected_backend is None:
        if config.task_spec.requires_remote_code and not config.allow_remote_code:
            raise ConfigurationError(
                "The selected official task loads checkpoint-provided Python code.",
                code="MODEL_REMOTE_CODE_NOT_ALLOWED",
                context={"model": config.model, "task": config.task},
                hint=(
                    "Review the official repository and checkpoint, then set "
                    "allow_remote_code=True to opt in."
                ),
            )
        checkpoints = ensure_prediction_checkpoint(config)
        try:
            selected_backend = create_prediction_backend(checkpoints, config)
        except (
            BackendExecutionError,
            BackendUnavailableError,
            ConfigurationError,
            SequenceError,
        ):
            raise
        except Exception as exc:
            raise BackendExecutionError(
                "Could not initialize the property-prediction backend.",
                code="MODEL_LOAD_FAILED",
                context={
                    "model": config.model,
                    "task": config.task,
                    "error_type": type(exc).__name__,
                },
            ) from exc
    elif not hasattr(selected_backend, "predict"):
        raise ConfigurationError(
            "backend must provide predict(inputs, show_progress=...).",
            code="INVALID_PREDICTION_BACKEND",
        )

    try:
        raw_outputs = tuple(
            selected_backend.predict(materialized, show_progress=config.show_progress)
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
            "The property-prediction backend failed.",
            code="PROPERTY_PREDICTION_FAILED",
            context={
                "model": config.model,
                "task": config.task,
                "error_type": type(exc).__name__,
            },
        ) from exc
    if len(raw_outputs) != len(materialized) or any(
        not isinstance(item, PredictionOutput) for item in raw_outputs
    ):
        raise BackendExecutionError(
            "The backend must return one PredictionOutput per input.",
            code="INVALID_PREDICTION_OUTPUT",
            context={"expected": len(materialized), "actual": len(raw_outputs)},
        )
    records = tuple(
        PredictionRecord(_input_id(item), output)
        for item, output in zip(materialized, raw_outputs, strict=True)
    )
    return PropertyPredictionResult(
        config.model,
        config.task,
        input_kind,
        records,
        checkpoints.paths,
        {
            "model_display_name": config.model_spec.display_name,
            "task_display_name": config.task_spec.display_name,
            "task_description": config.task_spec.description,
            "output_kind": config.task_spec.output_kind,
            "source_repository": config.model_spec.source_repository,
            "checkpoint_sources": checkpoints.sources,
            "checkpoint_downloaded": checkpoints.downloaded,
            "fine_tuning_performed": False,
            "ambiguity_policy": config.ambiguity_policy,
        },
    )


def predict_sequence_properties(
    inputs: Iterable[BiologicalSequence | DNARecord],
    *,
    config: PropertyPredictionConfig,
    backend: PredictionBackend | None = None,
) -> PropertyPredictionResult:
    """Convenience wrapper for a single-sequence pretrained task."""

    return predict_properties(inputs, config=config, backend=backend)


def predict_pair_properties(
    inputs: Iterable[BiologicalSequencePair],
    *,
    config: PropertyPredictionConfig,
    backend: PredictionBackend | None = None,
) -> PropertyPredictionResult:
    """Convenience wrapper for a pairwise pretrained task."""

    return predict_properties(inputs, config=config, backend=backend)


def predict_variant_effects(
    inputs: Iterable[VariantContext],
    *,
    config: PropertyPredictionConfig,
    backend: PredictionBackend | None = None,
) -> PropertyPredictionResult:
    """Convenience wrapper for an official sequence-context variant task."""

    return predict_properties(inputs, config=config, backend=backend)


__all__ = [
    "predict_pair_properties",
    "predict_properties",
    "predict_sequence_properties",
    "predict_variant_effects",
]
