"""Deterministic train/validation/test splitting for custom property data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from dnakit.exceptions import BackendUnavailableError, InputFormatError

from .data import PropertyData

_SPLIT_ALIASES = {
    "train": "train",
    "training": "train",
    "valid": "validation",
    "validation": "validation",
    "val": "validation",
    "dev": "validation",
    "test": "test",
    "testing": "test",
}


def _numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - provided by the training extra.
        raise BackendUnavailableError(
            "Custom property splitting requires NumPy.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the training extra with: python -m pip install "dnakit[training]"',
        ) from exc
    return np


@dataclass(frozen=True, slots=True)
class SplitAssignment:
    """Row indices for the three non-overlapping supervised splits."""

    train: tuple[int, ...]
    validation: tuple[int, ...]
    test: tuple[int, ...]
    source: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        groups = (self.train, self.validation, self.test)
        if any(not group for group in groups):
            raise InputFormatError(
                "Training, validation, and test splits must all be non-empty.",
                code="INVALID_CUSTOM_SPLIT",
            )
        flattened = (*self.train, *self.validation, *self.test)
        if len(set(flattened)) != len(flattened):
            raise InputFormatError("Custom-property splits overlap.", code="INVALID_CUSTOM_SPLIT")

    @property
    def counts(self) -> dict[str, int]:
        return {
            "train": len(self.train),
            "validation": len(self.validation),
            "test": len(self.test),
        }


def _counts(total: int, ratios: tuple[float, float, float]) -> tuple[int, int, int]:
    if total < 3:
        raise InputFormatError(
            "At least 3 records are required for an 8:1:1-style split.",
            code="CUSTOM_DATA_TOO_SMALL",
            context={"record_count": total},
        )
    raw = [total * ratio for ratio in ratios]
    counts = [max(1, int(value)) for value in raw]
    while sum(counts) > total:
        candidates = [index for index, value in enumerate(counts) if value > 1]
        if not candidates:
            break
        index = max(candidates, key=lambda item: (counts[item] - raw[item], counts[item]))
        counts[index] -= 1
    while sum(counts) < total:
        index = max(range(3), key=lambda item: (raw[item] - counts[item], ratios[item]))
        counts[index] += 1
    return counts[0], counts[1], counts[2]


def _supplied_assignment(data: PropertyData) -> SplitAssignment | None:
    if all(value is None for value in data.supplied_splits):
        return None
    groups: dict[str, list[int]] = {"train": [], "validation": [], "test": []}
    for index, raw in enumerate(data.supplied_splits):
        key = "" if raw is None else str(raw).strip().lower()
        canonical = _SPLIT_ALIASES.get(key)
        if canonical is None:
            raise InputFormatError(
                "The mapped split field contains an unsupported value.",
                code="INVALID_CUSTOM_SPLIT",
                context={"record_index": index, "value": raw},
                hint="Use train, validation/valid/val, or test.",
            )
        groups[canonical].append(index)
    return SplitAssignment(
        tuple(groups["train"]),
        tuple(groups["validation"]),
        tuple(groups["test"]),
        "source",
    )


def _stratification_labels(data: PropertyData) -> list[str] | None:
    if len(data.label_names) == 1:
        return [str(row[0]) for row in data.labels]
    binary: list[list[int]] = []
    for row in data.labels:
        values: list[int] = []
        for raw in row:
            text = str(raw).strip().lower()
            if text in {"0", "false", "no", "negative", "neg"}:
                values.append(0)
            elif text in {"1", "true", "yes", "positive", "pos"}:
                values.append(1)
            else:
                return None
        binary.append(values)
    if all(sum(row) == 1 for row in binary):
        return [str(row.index(1)) for row in binary]
    return None


def _stratified_indices(
    total: int,
    counts: tuple[int, int, int],
    labels: list[str],
    *,
    seed: int,
) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    try:
        from sklearn.model_selection import (  # type: ignore[import-not-found]
            StratifiedShuffleSplit,
        )
    except ImportError as exc:
        raise BackendUnavailableError(
            "Stratified custom-property splitting requires scikit-learn.",
            code="MISSING_NEURAL_DEPENDENCY",
            hint='Install the training extra with: python -m pip install "dnakit[training]"',
        ) from exc
    np = _numpy()
    indices = np.arange(total)
    train_count, validation_count, test_count = counts
    first = StratifiedShuffleSplit(
        n_splits=1,
        train_size=train_count,
        test_size=validation_count + test_count,
        random_state=seed,
    )
    train_positions, remaining_positions = next(first.split(indices, labels))
    remaining = indices[remaining_positions]
    remaining_labels = np.asarray(labels, dtype=object)[remaining_positions]
    second = StratifiedShuffleSplit(
        n_splits=1,
        train_size=validation_count,
        test_size=test_count,
        random_state=seed + 1,
    )
    validation_positions, test_positions = next(second.split(remaining, remaining_labels))
    return (
        tuple(int(value) for value in indices[train_positions]),
        tuple(int(value) for value in remaining[validation_positions]),
        tuple(int(value) for value in remaining[test_positions]),
    )


def split_property_data(
    data: PropertyData,
    *,
    train_ratio: float = 0.8,
    validation_ratio: float = 0.1,
    test_ratio: float = 0.1,
    seed: int = 42,
    task_type: str = "classification",
    stratify: bool = True,
) -> SplitAssignment:
    """Use an existing split field or produce a seeded random 8:1:1 split."""

    supplied = _supplied_assignment(data)
    if supplied is not None:
        return supplied
    total = data.record_count
    counts = _counts(total, (train_ratio, validation_ratio, test_ratio))
    warnings: list[str] = []
    if task_type == "classification" and stratify:
        labels = _stratification_labels(data)
        if labels is not None:
            try:
                train, validation, test = _stratified_indices(
                    total,
                    counts,
                    labels,
                    seed=seed,
                )
                return SplitAssignment(train, validation, test, "random_stratified")
            except ValueError:
                warnings.append(
                    "Stratification was infeasible for the class counts and split sizes; "
                    "a seeded random split was used."
                )
        else:
            warnings.append(
                "Automatic stratification supports one categorical target or one-hot "
                "multiclass targets; a seeded random split was used."
            )
    np = _numpy()
    permutation = np.random.default_rng(seed).permutation(total)
    train_count, validation_count, _ = counts
    train = tuple(int(value) for value in permutation[:train_count])
    validation = tuple(
        int(value) for value in permutation[train_count : train_count + validation_count]
    )
    test = tuple(int(value) for value in permutation[train_count + validation_count :])
    return SplitAssignment(train, validation, test, "random", tuple(warnings))


__all__ = ["SplitAssignment", "split_property_data"]
