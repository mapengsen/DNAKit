"""User-facing tests for the unified DNA entry point."""

from io import StringIO

import pytest

import dnakit
from dnakit.alignment import align_pairwise
from dnakit.core import DNA, DNAFeature, DNARecord, DNASequence, Gap, Interval
from dnakit.datasets import exclude_species
from dnakit.descriptors import gc_at_content
from dnakit.evaluation import create_reference_library
from dnakit.exceptions import ConfigurationError, SequenceError
from dnakit.fingerprints import integer_encode
from dnakit.io import filter_by_metadata, with_metadata
from dnakit.ops import RecordOperationResult, delete, reverse_complement
from dnakit.patterns import scan_motif
from dnakit.similarity import exact_similarity
from dnakit.standardize import DatasetValidationReport, ValidationReport, validate
from dnakit.visualization import plot_sequence


def test_dna_builds_one_normalized_record_with_optional_information() -> None:
    data = DNA(
        " acgu ",
        id="seq-1",
        description="example",
        topology="circular",
        metadata={"species": "synthetic"},
        features=[
            {
                "type": "motif",
                "start": 0,
                "end": 2,
                "strand": "forward",
                "label": "start",
            }
        ],
    )

    assert data.is_single
    assert data.record_count == 1
    assert data.symbols == "ACG"
    assert data.id == "seq-1"
    assert data.topology.value == "circular"
    assert data.metadata["species"] == "synthetic"
    assert data.features == (DNAFeature("motif", Interval(0, 2), strand="forward", label="start"),)
    assert data.normalization is not None
    assert data.normalization.was_modified


def test_dna_builds_multiple_records_with_generated_or_explicit_ids() -> None:
    generated = DNA(["ACGT", "TTAA"])
    described = DNA(
        [
            {"sequence": "ACGT", "id": "first", "metadata": {"group": "a"}},
            {
                "sequence": "ttaa",
                "id": "second",
                "topology": "circular",
                "features": [{"type": "site", "start": 1, "end": 3}],
            },
        ],
        name="example",
    )

    assert isinstance(generated, DNA)
    assert generated.ids == ("sequence_1", "sequence_2")
    assert described.ids == ("first", "second")
    assert described.name == "example"
    assert isinstance(described[0], DNA)
    assert isinstance(described[:1], DNA)
    assert described[:1].ids == ("first",)
    assert described[0].metadata["group"] == "a"
    assert described[1].sequence.symbols == "TTAA"
    assert described[1].features[0].location == Interval(1, 3)
    with pytest.raises(SequenceError) as error:
        _ = described.sequence
    assert error.value.code == "SINGLE_DNA_REQUIRED"


def test_dna_accepts_legacy_values_and_explicit_gapped_parts() -> None:
    record = DNARecord(DNASequence("ACGT"), "legacy")
    wrapped_record = DNA(record)
    wrapped_set = DNA(dnakit.DNASet([record]))
    gapped = DNA(["AC", Gap(3), "GT"], id="gapped")

    assert wrapped_record.record == record
    assert wrapped_set.record == record
    assert gapped.sequence.parts == ("AC", Gap(3), "GT")


def test_dna_rejects_ambiguous_multi_record_options_and_unbounded_input() -> None:
    with pytest.raises(ConfigurationError) as shared_id:
        DNA(["AC", "GT"], id="same")
    assert shared_id.value.code == "MULTIPLE_DNA_SHARED_ID"

    with pytest.raises(ConfigurationError) as unknown_field:
        DNA({"sequence": "AC", "unknown": True})
    assert unknown_field.value.code == "UNKNOWN_DNA_RECORD_FIELD"

    with pytest.raises(ConfigurationError) as limit:
        DNA(("A" for _ in range(3)), max_records=2)
    assert limit.value.code == "DNA_RECORD_LIMIT_EXCEEDED"


def test_dna_works_with_common_analysis_entry_points() -> None:
    left = DNA("ACGT", id="left")
    right = DNA("ACGA", id="right")

    report = validate(left)
    assert isinstance(report, ValidationReport)
    assert report.record_id == "left"
    assert gc_at_content(left).gc_fraction == 0.5
    assert integer_encode(left).values == (0, 1, 2, 3)
    assert len(scan_motif(left, "CG").hits) == 2
    assert exact_similarity(left, right).value == 0.0
    assert align_pairwise(left, right).score == 2.0
    assert plot_sequence(left).kind == "sequence"

    dataset_report = validate(DNA(["ACGT", "TTAA"]))
    assert isinstance(dataset_report, DatasetValidationReport)
    assert dataset_report.record_count == 2


def test_unified_operations_dispatch_without_record_suffixes() -> None:
    data = DNA(
        "AACCGG",
        id="record-1",
        features=[{"type": "site", "start": 1, "end": 5}],
    )

    deleted = delete(data, 2, 4)
    reversed_data = reverse_complement(data)

    assert isinstance(deleted, DNA)
    assert deleted.symbols == "AAGG"
    assert isinstance(reversed_data, DNA)
    assert reversed_data.symbols == "CCGGTT"

    legacy = delete(DNASequence("AACCGG"), 2, 4)
    assert legacy.sequence.symbols == "AAGG"
    audited = delete(data.record, 2, 4)
    assert isinstance(audited, RecordOperationResult)
    assert audited.record.sequence.symbols == "AAGG"


def test_metadata_filters_and_reference_library_accept_the_facade() -> None:
    data = DNA(
        [
            {"sequence": "ACGT", "id": "a", "metadata": {"species": "human"}},
            {"sequence": "TTAA", "id": "b", "metadata": {"species": "mouse"}},
        ],
        name="source",
    )

    annotated = with_metadata(DNA("ACGT", id="one"), {"species": "human"})
    human = filter_by_metadata(data, {"species": "human"})
    without_mouse = exclude_species(data, "mouse")
    reference = create_reference_library(
        data,
        name="reference",
        version="1",
        source="local",
    )

    assert isinstance(annotated, DNA)
    assert annotated.metadata["species"] == "human"
    assert annotated.normalization is not None
    assert isinstance(human, DNA)
    assert human.ids == ("a",)
    assert human.normalization == data.normalizations[0]
    assert isinstance(without_mouse, DNA)
    assert without_mouse.ids == ("a",)
    assert reference.records.ids == ("a", "b")


def test_read_mode_dna_uses_the_same_facade_for_one_or_many_records() -> None:
    loaded = dnakit.read(
        StringIO(">a\nACGT\n>b\nTTAA\n"),
        format="fasta",
        mode="dna",
    )

    assert isinstance(loaded, DNA)
    assert loaded.ids == ("a", "b")
    assert [record.sequence.symbols for record in loaded] == ["ACGT", "TTAA"]
