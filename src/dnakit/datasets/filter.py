"""Explicit metadata-based exclusion for immutable DNA record sets."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Literal, TypeAlias, overload

from dnakit.core import DNA, DNARecord, DNASet
from dnakit.exceptions import ConfigurationError, InputFormatError

MetadataExclusionField: TypeAlias = Literal["chromosome", "species"]
MetadataValues: TypeAlias = object | Iterable[object]
MetadataDataset: TypeAlias = DNA | DNASet

_EXCLUSION_FIELDS = frozenset({"chromosome", "species"})


def _copy_set(dataset: DNASet, records: Iterable[DNARecord]) -> DNASet:
    return DNASet(
        records,
        name=dataset.name,
        source=dataset.source,
        version=dataset.version,
        metadata=dataset.metadata,
    )


def _materialize_values(values: MetadataValues) -> tuple[object, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Iterable):
        materialized = (values,)
    else:
        materialized = tuple(values)
    if not materialized:
        raise ConfigurationError(
            "At least one metadata value must be excluded.",
            code="EMPTY_METADATA_EXCLUSION_VALUES",
        )
    return materialized


@overload
def exclude_by_metadata(
    dataset: DNA,
    field: MetadataExclusionField,
    values: MetadataValues,
) -> DNA: ...


@overload
def exclude_by_metadata(
    dataset: DNASet,
    field: MetadataExclusionField,
    values: MetadataValues,
) -> DNASet: ...


def exclude_by_metadata(
    dataset: MetadataDataset,
    field: MetadataExclusionField,
    values: MetadataValues,
) -> DNA | DNASet:
    """Return a copy without records whose chromosome or species matches.

    Matching is exact and compares the stored JSON-compatible metadata value.
    Every record must contain ``field``; missing fields raise before any record
    is removed so that an incomplete dataset cannot be silently filtered.
    """

    if not isinstance(dataset, (DNA, DNASet)):
        raise TypeError("dataset must be DNA or DNASet.")
    source = dataset.dataset if isinstance(dataset, DNA) else dataset
    if not isinstance(field, str) or field not in _EXCLUSION_FIELDS:
        raise ConfigurationError(
            "Metadata exclusion field must be 'chromosome' or 'species'.",
            code="INVALID_METADATA_EXCLUSION_FIELD",
            context={"field": field},
        )
    excluded_values = _materialize_values(values)
    missing_ids = tuple(record.id for record in source if field not in record.metadata)
    if missing_ids:
        raise InputFormatError(
            "Every DNASet record must contain the metadata field used for exclusion.",
            code="MISSING_METADATA_FIELD_FOR_EXCLUSION",
            context={
                "field": field,
                "missing_record_ids": missing_ids[:10],
                "missing_count": len(missing_ids),
            },
            hint=f"Populate metadata[{field!r}] for every record before exclusion.",
        )
    selected = tuple(
        (index, record)
        for index, record in enumerate(source)
        if not any(record.metadata[field] == value for value in excluded_values)
    )
    if isinstance(dataset, DNA):
        return dataset._derive(selected)
    return _copy_set(dataset, (record for _, record in selected))


@overload
def exclude_chromosomes(dataset: DNA, chromosomes: MetadataValues) -> DNA: ...


@overload
def exclude_chromosomes(dataset: DNASet, chromosomes: MetadataValues) -> DNASet: ...


def exclude_chromosomes(dataset: MetadataDataset, chromosomes: MetadataValues) -> DNA | DNASet:
    """Return a copy without records from the selected chromosome value(s)."""

    return exclude_by_metadata(dataset, "chromosome", chromosomes)


@overload
def exclude_species(dataset: DNA, species: MetadataValues) -> DNA: ...


@overload
def exclude_species(dataset: DNASet, species: MetadataValues) -> DNASet: ...


def exclude_species(dataset: MetadataDataset, species: MetadataValues) -> DNA | DNASet:
    """Return a copy without records from the selected species value(s)."""

    return exclude_by_metadata(dataset, "species", species)


__all__ = [
    "MetadataExclusionField",
    "MetadataValues",
    "exclude_by_metadata",
    "exclude_chromosomes",
    "exclude_species",
]
