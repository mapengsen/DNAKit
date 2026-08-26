from __future__ import annotations

import pytest

from dnakit import DNARecord, DNASequence, DNASet
from dnakit.datasets import exclude_by_metadata, exclude_chromosomes, exclude_species
from dnakit.exceptions import ConfigurationError, InputFormatError


def _dataset() -> DNASet:
    return DNASet(
        (
            DNARecord(
                DNASequence("A"),
                "human-chr1",
                metadata={"species": "human", "chromosome": "chr1"},
            ),
            DNARecord(
                DNASequence("C"),
                "human-chr2",
                metadata={"species": "human", "chromosome": "chr2"},
            ),
            DNARecord(
                DNASequence("G"),
                "mouse-chr1",
                metadata={"species": "mouse", "chromosome": "chr1"},
            ),
        ),
        name="genomes",
        source="test",
        version="v1",
        metadata={"owner": "test"},
    )


def test_exclude_chromosomes_accepts_one_or_many_values_and_preserves_dataset_metadata() -> None:
    dataset = _dataset()

    one = exclude_chromosomes(dataset, "chr2")
    many = exclude_chromosomes(dataset, ["chr1"])

    assert one.ids == ("human-chr1", "mouse-chr1")
    assert many.ids == ("human-chr2",)
    assert dataset.ids == ("human-chr1", "human-chr2", "mouse-chr1")
    assert (one.name, one.source, one.version, one.metadata) == (
        dataset.name,
        dataset.source,
        dataset.version,
        dataset.metadata,
    )


def test_exclude_species_and_generic_function_use_exact_metadata_values() -> None:
    dataset = _dataset()

    assert exclude_species(dataset, ["human"]).ids == ("mouse-chr1",)
    assert exclude_by_metadata(dataset, "species", "mouse").ids == (
        "human-chr1",
        "human-chr2",
    )


@pytest.mark.parametrize(
    "field, function, values",
    [
        ("chromosome", exclude_chromosomes, "chr1"),
        ("species", exclude_species, "human"),
    ],
)
def test_exclusion_errors_when_any_record_lacks_the_target_field(
    field: str,
    function: object,
    values: object,
) -> None:
    dataset = DNASet(
        (
            DNARecord(DNASequence("A"), "complete", metadata={field: values}),
            DNARecord(DNASequence("C"), "missing"),
        )
    )

    with pytest.raises(InputFormatError) as exc_info:
        function(dataset, values)  # type: ignore[operator]

    assert exc_info.value.code == "MISSING_METADATA_FIELD_FOR_EXCLUSION"
    assert exc_info.value.context["field"] == field
    assert exc_info.value.context["missing_record_ids"] == ("missing",)


def test_exclusion_rejects_unsupported_fields_and_empty_values() -> None:
    with pytest.raises(ConfigurationError) as field_error:
        exclude_by_metadata(_dataset(), "batch", "1")  # type: ignore[call-overload]
    assert field_error.value.code == "INVALID_METADATA_EXCLUSION_FIELD"

    with pytest.raises(ConfigurationError) as values_error:
        exclude_species(_dataset(), [])
    assert values_error.value.code == "EMPTY_METADATA_EXCLUSION_VALUES"
