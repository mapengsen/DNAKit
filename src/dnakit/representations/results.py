"""Results returned by DNA representation extraction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True, eq=False)
class RepresentationResult:
    """A read-only float32 matrix aligned with the input record IDs."""

    record_ids: tuple[str, ...]
    representations: Any
    model_name: str
    checkpoint_path: str | None
    pooling: str
    embedding_dimension: int
    input_count: int
    ambiguity_policy: str

    def __post_init__(self) -> None:
        try:
            import numpy as np
        except ImportError as exc:  # pragma: no cover - reached through the neural extra.
            raise ConfigurationError(
                "RepresentationResult requires NumPy.",
                code="MISSING_NEURAL_DEPENDENCY",
            ) from exc
        matrix = np.asarray(self.representations, dtype=np.float32)
        if matrix.ndim != 2:
            raise ConfigurationError(
                "representations must be a two-dimensional matrix.",
                code="INVALID_REPRESENTATION_MATRIX",
            )
        if matrix.shape != (self.input_count, self.embedding_dimension):
            raise ConfigurationError(
                "Representation matrix shape does not match result metadata.",
                code="INVALID_REPRESENTATION_MATRIX",
                context={
                    "shape": tuple(int(value) for value in matrix.shape),
                    "input_count": self.input_count,
                    "embedding_dimension": self.embedding_dimension,
                },
            )
        if len(self.record_ids) != self.input_count:
            raise ConfigurationError(
                "record_ids must align with representation rows.",
                code="INVALID_REPRESENTATION_MATRIX",
            )
        if not np.isfinite(matrix).all():
            raise ConfigurationError(
                "Representations must contain only finite numbers.",
                code="INVALID_REPRESENTATION_MATRIX",
            )
        matrix = np.array(matrix, dtype=np.float32, copy=True, order="C")
        matrix.setflags(write=False)
        object.__setattr__(self, "representations", matrix)

    def to_dict(self, *, include_representations: bool = True) -> dict[str, object]:
        """Return JSON-compatible metadata and, by default, the extracted vectors."""

        payload: dict[str, object] = {
            "record_ids": self.record_ids,
            "model_name": self.model_name,
            "checkpoint_path": self.checkpoint_path,
            "pooling": self.pooling,
            "embedding_dimension": self.embedding_dimension,
            "input_count": self.input_count,
            "ambiguity_policy": self.ambiguity_policy,
        }
        if include_representations:
            payload["representations"] = self.representations.tolist()
        return payload


__all__ = ["RepresentationResult"]
