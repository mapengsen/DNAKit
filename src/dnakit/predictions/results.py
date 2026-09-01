"""Immutable outputs returned by direct property prediction."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - prediction extras provide NumPy.
        raise BackendUnavailableError(
            "Property prediction requires NumPy.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the neural extra with: python -m pip install "dnakit[neural]"',
        ) from exc
    return np


@dataclass(frozen=True, slots=True, eq=False)
class PredictionOutput:
    """One finite read-only numeric output plus names for its final axis."""

    values: Any
    output_names: tuple[str, ...] = ()
    metadata: FrozenDict | Mapping[str, object] = field(default_factory=FrozenDict)

    def __post_init__(self) -> None:
        np = _numpy()
        try:
            values = np.asarray(self.values, dtype=np.float32)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Prediction values must be numeric.",
                code="INVALID_PREDICTION_OUTPUT",
            ) from exc
        if values.ndim == 0:
            values = values.reshape(1)
        if values.size < 1:
            raise ConfigurationError(
                "Prediction values must not be empty.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        if not bool(np.isfinite(values).all()):
            raise ConfigurationError(
                "Prediction values must contain only finite numbers.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        names = tuple(self.output_names)
        if any(not isinstance(name, str) or not name.strip() for name in names):
            raise ConfigurationError(
                "output_names must contain non-empty strings.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        if names and len(names) != int(values.shape[-1]):
            raise ConfigurationError(
                "output_names must align with the final values axis.",
                code="INVALID_PREDICTION_OUTPUT",
                context={"name_count": len(names), "final_axis": int(values.shape[-1])},
            )
        metadata = self.metadata
        if not isinstance(metadata, Mapping):
            raise ConfigurationError(
                "Prediction metadata must be a mapping.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        frozen_values = np.array(values, dtype=np.float32, copy=True, order="C")
        frozen_values.setflags(write=False)
        object.__setattr__(self, "values", frozen_values)
        object.__setattr__(self, "output_names", names)
        object.__setattr__(
            self,
            "metadata",
            metadata if isinstance(metadata, FrozenDict) else FrozenDict(metadata),
        )

    def to_dict(
        self,
        *,
        include_values: bool = True,
        include_metadata: bool = True,
    ) -> dict[str, object]:
        """Return a JSON-compatible output envelope."""

        payload: dict[str, object] = {
            "shape": tuple(int(value) for value in self.values.shape),
            "output_names": self.output_names,
        }
        if include_values:
            payload["values"] = self.values.tolist()
        if include_metadata:
            payload["metadata"] = to_json_compatible(self.metadata)
        return payload


@dataclass(frozen=True, slots=True, eq=False)
class PredictionRecord:
    """A model output aligned with one caller-provided input ID."""

    record_id: str
    output: PredictionOutput

    def __post_init__(self) -> None:
        if not isinstance(self.record_id, str) or not self.record_id.strip():
            raise ConfigurationError(
                "record_id must be non-empty text.", code="INVALID_PREDICTION_OUTPUT"
            )
        if not isinstance(self.output, PredictionOutput):
            raise ConfigurationError(
                "output must be PredictionOutput.", code="INVALID_PREDICTION_OUTPUT"
            )

    def to_dict(
        self,
        *,
        include_values: bool = True,
        include_metadata: bool = True,
    ) -> dict[str, object]:
        """Return this record as JSON-compatible data."""

        return {
            "record_id": self.record_id,
            "output": self.output.to_dict(
                include_values=include_values,
                include_metadata=include_metadata,
            ),
        }


@dataclass(frozen=True, slots=True, eq=False)
class PropertyPredictionResult:
    """Auditable outputs for one model/task invocation."""

    model_name: str
    task_name: str
    input_kind: str
    records: tuple[PredictionRecord, ...]
    checkpoint_paths: tuple[str, ...] = ()
    metadata: FrozenDict | Mapping[str, object] = field(default_factory=FrozenDict)

    def __post_init__(self) -> None:
        for field_name in ("model_name", "task_name", "input_kind"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(
                    f"{field_name} must be non-empty text.",
                    code="INVALID_PREDICTION_OUTPUT",
                )
        records = tuple(self.records)
        if not records or any(not isinstance(item, PredictionRecord) for item in records):
            raise ConfigurationError(
                "records must contain at least one PredictionRecord.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        record_ids = tuple(item.record_id for item in records)
        if len(set(record_ids)) != len(record_ids):
            raise ConfigurationError(
                "Prediction record IDs must be unique.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        paths = tuple(self.checkpoint_paths)
        if any(not isinstance(path, str) or not path.strip() for path in paths):
            raise ConfigurationError(
                "checkpoint_paths must contain non-empty paths.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        metadata = self.metadata
        if not isinstance(metadata, Mapping):
            raise ConfigurationError(
                "Prediction result metadata must be a mapping.",
                code="INVALID_PREDICTION_OUTPUT",
            )
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "checkpoint_paths", paths)
        object.__setattr__(
            self,
            "metadata",
            metadata if isinstance(metadata, FrozenDict) else FrozenDict(metadata),
        )

    @property
    def record_ids(self) -> tuple[str, ...]:
        """Return input IDs in prediction order."""

        return tuple(item.record_id for item in self.records)

    def to_dict(
        self,
        *,
        include_values: bool = True,
        include_output_metadata: bool = True,
    ) -> dict[str, object]:
        """Return JSON-compatible metadata and optionally the numeric arrays."""

        return {
            "model_name": self.model_name,
            "task_name": self.task_name,
            "input_kind": self.input_kind,
            "record_count": len(self.records),
            "checkpoint_paths": self.checkpoint_paths,
            "metadata": to_json_compatible(self.metadata),
            "records": [
                record.to_dict(
                    include_values=include_values,
                    include_metadata=include_output_metadata,
                )
                for record in self.records
            ],
        }


__all__ = ["PredictionOutput", "PredictionRecord", "PropertyPredictionResult"]
