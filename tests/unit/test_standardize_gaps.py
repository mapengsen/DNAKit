"""Tests for strict, auditable STD-007 gap standardization."""

from __future__ import annotations

import itertools
import json
from typing import cast

import pytest

from dnakit.core.coordinates import Interval
from dnakit.core.enums import DNAAlphabet, GapKind, Topology
from dnakit.core.feature import DNAFeature
from dnakit.core.gap import Gap
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, SequenceError, UnsupportedGapOperationError
from dnakit.io.annotations import AGPComponent, AGPDocument, AGPGap
from dnakit.standardize import (
    AGPAssemblyConfig,
    GapNormalizationConfig,
    normalize_gaps,
    sequence_from_agp,
)


def _component(
    start: int,
    end: int,
    part: int,
    component_id: str,
    *,
    component_start: int = 0,
    orientation: str = "+",
    object_id: str = "chr1",
) -> AGPComponent:
    return AGPComponent(
        object_id,
        Interval(start, end),
        part,
        "W",
        component_id,
        Interval(component_start, component_start + end - start),
        orientation,
    )


def _known_gap(start: int, end: int, part: int, *, object_id: str = "chr1") -> AGPGap:
    gap = Gap(
        end - start,
        GapKind.SCAFFOLD,
        True,
        ("paired-ends",),
        {"agp_gap_type": "scaffold"},
    )
    return AGPGap(
        object_id,
        Interval(start, end),
        part,
        "N",
        gap,
        "scaffold",
        True,
        ("paired-ends",),
    )


def _unknown_gap(start: int, part: int, *, object_id: str = "chr1") -> AGPGap:
    gap = Gap(
        None,
        GapKind.CONTIG,
        False,
        (),
        {"agp_gap_type": "contig"},
    )
    return AGPGap(
        object_id,
        Interval(start, start + 100),
        part,
        "U",
        gap,
        "contig",
        False,
        (),
    )


def test_n_runs_are_converted_with_config_and_zero_based_audit() -> None:
    explicit = Gap(3, GapKind.CONTIG, False, ("manual",), {"source": "curator"})
    source = DNASequence(
        ["AANNNNNCCNN", explicit, "NNNNGG"],
        alphabet=DNAAlphabet.IUPAC,
    )
    result = normalize_gaps(
        source,
        config=GapNormalizationConfig(
            min_run_length=4,
            kind=GapKind.SCAFFOLD,
            crossable=True,
            evidence=("alignment",),
            metadata={"confidence": "high"},
        ),
    )

    assert result.was_modified
    assert len(result.source_sha256) == 64
    assert result.preserved_gap_count == 1
    assert len(result.changes) == 2
    first, second = result.changes
    assert first.part_interval == Interval(2, 7)
    assert first.symbol_interval == Interval(2, 7)
    assert first.coordinate_interval == Interval(2, 7)
    assert second.symbol_interval == Interval(11, 15)
    assert second.coordinate_interval == Interval(14, 18)
    assert first.original_symbol == "N"
    assert first.original_length == 5
    assert len(first.original_sha256) == 64
    assert first.replacement.kind is GapKind.SCAFFOLD
    assert first.replacement.crossable is True
    assert first.replacement.evidence == ("alignment",)
    assert first.replacement.metadata["confidence"] == "high"
    audit = first.replacement.metadata["dnakit_gap_normalization"]
    assert audit == {
        "source": "N-run",
        "original_symbol": "N",
        "original_length": 5,
        "source_part_index": 0,
        "part_interval": {"start": 2, "end": 7},
        "symbol_interval": {"start": 2, "end": 7},
    }
    assert result.sequence.parts[3] is explicit
    assert result.sequence.parts == (
        "AA",
        first.replacement,
        "CCNN",
        explicit,
        second.replacement,
        "GG",
    )
    assert result.sequence.coordinate_span == source.coordinate_span
    assert result.provenance.implementation.label.value == "native"
    json.dumps(result.to_dict(), ensure_ascii=False, sort_keys=True)


def test_short_n_runs_and_existing_gaps_are_preserved_by_identity() -> None:
    gap = Gap(None, GapKind.UNKNOWN, None, ("assembly",))
    source = DNASequence(["ANNN", gap, "NNGT"], alphabet=DNAAlphabet.IUPAC)

    result = normalize_gaps(source, config=GapNormalizationConfig(min_run_length=4))

    assert result.output is source
    assert result.sequence is source
    assert not result.was_modified
    assert result.sequence.parts[1] is gap
    assert result.preserved_gap_count == 1


def test_change_after_unknown_gap_has_unresolved_coordinate_but_known_symbol_interval() -> None:
    source = DNASequence(
        ["AC", Gap(None), "NNNN"],
        alphabet=DNAAlphabet.IUPAC,
    )

    result = normalize_gaps(source, config=GapNormalizationConfig(min_run_length=4))

    assert result.changes[0].symbol_interval == Interval(2, 6)
    assert result.changes[0].coordinate_interval is None
    assert result.sequence.coordinate_span is None


def test_record_features_are_retained_because_known_gap_span_is_unchanged() -> None:
    feature = DNAFeature("region", Interval(6, 8), id="after-gap")
    record = DNARecord(
        DNASequence("AANNNNCC", alphabet=DNAAlphabet.IUPAC),
        "record-1",
        description="audited",
        features=(feature,),
        metadata={"sample": "S1"},
    )

    result = normalize_gaps(record, config=GapNormalizationConfig(min_run_length=4))

    assert result.record is not None
    assert result.output is result.record
    assert result.record.id == "record-1"
    assert result.record.description == "audited"
    assert result.record.features == (feature,)
    assert result.record.metadata == {"sample": "S1"}
    assert result.record.sequence.coordinate_span == 8


def test_letter_annotations_are_rejected_only_when_symbols_would_be_removed() -> None:
    changed = DNARecord(
        DNASequence("AANNNNCC", alphabet=DNAAlphabet.IUPAC),
        "changed",
        letter_annotations={"quality": range(8)},
    )
    unchanged = DNARecord(
        DNASequence("AANNCC", alphabet=DNAAlphabet.IUPAC),
        "unchanged",
        letter_annotations={"quality": range(6)},
    )

    with pytest.raises(UnsupportedGapOperationError) as error:
        normalize_gaps(changed, config=GapNormalizationConfig(min_run_length=4))
    assert error.value.code == "GAP_NORMALIZATION_LETTER_ANNOTATIONS_UNSUPPORTED"
    unchanged_result = normalize_gaps(
        unchanged,
        config=GapNormalizationConfig(min_run_length=4),
    )
    assert unchanged_result.output is unchanged


def test_letter_annotation_rejection_precedes_partial_iteration_side_effects() -> None:
    record = DNARecord(
        DNASequence("NNNN", alphabet=DNAAlphabet.IUPAC),
        "annotated",
        letter_annotations={"quality": range(4)},
    )

    with pytest.raises(UnsupportedGapOperationError) as error:
        normalize_gaps(record, config=GapNormalizationConfig(min_run_length=4))

    assert error.value.code == "GAP_NORMALIZATION_LETTER_ANNOTATIONS_UNSUPPORTED"
    assert record.sequence.symbols == "NNNN"


def test_circular_origin_run_is_never_silently_merged() -> None:
    source = DNASequence(
        "NNACNNN",
        alphabet=DNAAlphabet.IUPAC,
        topology=Topology.CIRCULAR,
    )

    separate = normalize_gaps(source, config=GapNormalizationConfig(min_run_length=4))
    assert separate.sequence is source
    assert [issue.code for issue in separate.issues] == ["STD_CIRCULAR_N_RUN_SPLIT_AT_ORIGIN"]

    with pytest.raises(UnsupportedGapOperationError) as error:
        normalize_gaps(
            source,
            config=GapNormalizationConfig(
                min_run_length=4,
                circular_boundary_policy="error",
            ),
        )
    assert error.value.code == "GAP_NORMALIZATION_CIRCULAR_BOUNDARY_UNRESOLVED"


def test_circular_all_n_output_and_resource_limits_are_explicit_errors() -> None:
    all_n = DNASequence("NNNN", alphabet=DNAAlphabet.IUPAC, topology=Topology.CIRCULAR)
    with pytest.raises(UnsupportedGapOperationError) as circular_error:
        normalize_gaps(all_n, config=GapNormalizationConfig(min_run_length=1))
    assert circular_error.value.code == "GAP_NORMALIZATION_EMPTY_CIRCULAR_OUTPUT"

    with pytest.raises(SequenceError) as input_error:
        normalize_gaps(
            DNASequence("ACGT"),
            config=GapNormalizationConfig(max_input_symbols=3),
        )
    assert input_error.value.code == "GAP_NORMALIZATION_INPUT_LIMIT_EXCEEDED"

    with pytest.raises(SequenceError) as run_error:
        normalize_gaps(
            DNASequence("NNNN", alphabet=DNAAlphabet.IUPAC),
            config=GapNormalizationConfig(min_run_length=2, max_n_run_length=3),
        )
    assert run_error.value.code == "GAP_NORMALIZATION_RUN_LIMIT_EXCEEDED"

    explicit_only = DNASequence((Gap(1), Gap(1), Gap(1)))
    with pytest.raises(SequenceError) as part_error:
        normalize_gaps(
            explicit_only,
            config=GapNormalizationConfig(max_output_parts=2),
        )
    assert part_error.value.code == "GAP_NORMALIZATION_OUTPUT_LIMIT_EXCEEDED"


def test_gap_normalization_config_is_strict_and_metadata_is_immutable() -> None:
    with pytest.raises(ConfigurationError):
        GapNormalizationConfig(min_run_length=True)
    with pytest.raises(ConfigurationError):
        GapNormalizationConfig(crossable=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        GapNormalizationConfig(evidence=("x", "x"))
    with pytest.raises(ConfigurationError):
        GapNormalizationConfig(metadata={"dnakit_gap_normalization": {}})
    with pytest.raises(ConfigurationError) as invalid_metadata:
        GapNormalizationConfig(metadata={"bad": object()})
    assert invalid_metadata.value.code == "INVALID_GAP_NORMALIZATION_METADATA"
    with pytest.raises(ConfigurationError):
        GapNormalizationConfig(circular_boundary_policy="merge")  # type: ignore[arg-type]

    original = {"nested": {"values": [1, 2]}}
    config = GapNormalizationConfig(metadata=original)
    original["nested"] = {"values": [99]}
    assert config.metadata["nested"] == {"values": (1, 2)}


def test_agp_known_gap_and_reverse_component_build_a_fragmented_record() -> None:
    gap_entry = _known_gap(4, 6, 2)
    document = AGPDocument(
        (
            _component(0, 4, 1, "left"),
            gap_entry,
            _component(6, 10, 3, "right", orientation="-"),
        )
    )

    result = sequence_from_agp(
        document,
        {"left": DNASequence("AACC"), "right": DNASequence("AGTC")},
    )

    assert result.record.id == "chr1"
    assert result.sequence.parts == ("AACC", gap_entry.gap, "GACT")
    assert result.sequence.coordinate_span == 10
    assert result.used_components == ("left", "right")
    assert result.segments[0].object_interval == Interval(0, 4)
    assert result.segments[0].agp_component_type == "W"
    assert result.segments[2].component_interval == Interval(0, 4)
    assert result.segments[1].gap_type == "scaffold"
    assert result.segments[1].linkage is True
    assert result.segments[1].linkage_evidence == ("paired-ends",)
    assert result.parameters["coordinate_system"] == "0-based-half-open"
    assert result.provenance.implementation.label.value == "reimplementation"
    assert not result.coordinate_span_unresolved
    json.dumps(result.to_dict(), ensure_ascii=False, sort_keys=True)


def test_agp_gap_type_is_authoritative_when_optional_metadata_is_absent() -> None:
    gap = Gap(2, GapKind.CONTIG, False)
    entry = AGPGap(
        "chr1",
        Interval(0, 2),
        1,
        "N",
        gap,
        "contig",
        False,
        (),
    )

    result = sequence_from_agp((entry,), {})

    assert result.sequence.parts == (gap,)
    assert result.segments[0].segment_type == "known_gap"


def test_agp_unknown_gap_retains_u_semantics_and_later_components() -> None:
    unknown = _unknown_gap(2, 2)
    document = AGPDocument(
        (
            _component(0, 2, 1, "left"),
            unknown,
            _component(102, 104, 3, "right"),
        )
    )

    result = sequence_from_agp(
        document,
        {"left": DNASequence("AA"), "right": DNASequence("TT")},
    )

    assert result.sequence.parts == ("AA", unknown.gap, "TT")
    assert result.sequence.parts[1] is unknown.gap
    assert result.sequence.coordinate_span is None
    assert result.coordinate_span_unresolved
    assert result.segments[1].segment_type == "unknown_gap"
    assert result.segments[1].object_interval == Interval(2, 102)
    assert result.segments[2].object_interval == Interval(102, 104)
    assert [issue.code for issue in result.issues] == ["STD_AGP_COORDINATE_SPAN_UNRESOLVED"]
    assert result.issues[0].details["agp_placeholder_intervals"] == ({"start": 2, "end": 102},)


def test_agp_iupac_and_record_components_preserve_declared_domain() -> None:
    document = AGPDocument((_component(0, 3, 1, "iupac"),))
    component = DNARecord(
        DNASequence("ANR", alphabet=DNAAlphabet.IUPAC),
        "component-record",
        metadata={"ignored": True},
    )

    result = sequence_from_agp(document, {"iupac": component}, record_id="assembled")

    assert result.record.id == "assembled"
    assert result.sequence.symbols == "ANR"
    assert result.sequence.alphabet is DNAAlphabet.IUPAC
    assert result.record.metadata["dnakit_agp_assembly"] == {
        "object_id": "chr1",
        "coordinate_system": "0-based-half-open",
        "entry_count": 1,
        "component_count": 1,
        "known_gap_count": 0,
        "unknown_gap_count": 0,
    }


def test_agp_object_selection_is_explicit_and_blocks_must_not_recur() -> None:
    entries = (
        _component(0, 2, 1, "a", object_id="chrA"),
        _component(0, 2, 1, "b", object_id="chrB"),
    )
    components = {"a": DNASequence("AA"), "b": DNASequence("CC")}

    with pytest.raises(ConfigurationError) as required:
        sequence_from_agp(entries, components)
    assert required.value.code == "AGP_OBJECT_ID_REQUIRED"
    assert sequence_from_agp(entries, components, object_id="chrB").sequence.symbols == "CC"

    recurring = (*entries, _component(2, 4, 2, "a", object_id="chrA"))
    with pytest.raises(SequenceError) as discontiguous:
        sequence_from_agp(recurring, components, object_id="chrA")
    assert discontiguous.value.code == "AGP_OBJECT_BLOCK_DISCONTIGUOUS"


def test_agp_continuity_missing_component_and_orientation_fail_without_guessing() -> None:
    bad_start = (_component(1, 3, 1, "a"),)
    with pytest.raises(SequenceError) as continuity:
        sequence_from_agp(bad_start, {"a": DNASequence("AAA")})
    assert continuity.value.code == "INVALID_AGP_ASSEMBLY_CONTINUITY"

    entry = _component(0, 2, 1, "missing")
    with pytest.raises(SequenceError) as missing:
        sequence_from_agp((entry,), {})
    assert missing.value.code == "AGP_COMPONENT_NOT_FOUND"

    unresolved = _component(0, 2, 1, "a", orientation="?")
    with pytest.raises(UnsupportedGapOperationError) as orientation:
        sequence_from_agp((unresolved,), {"a": DNASequence("AA")})
    assert orientation.value.code == "AGP_COMPONENT_ORIENTATION_UNRESOLVED"


def test_agp_component_bounds_gaps_and_circular_topology_are_rejected() -> None:
    interval = _component(0, 3, 1, "a", component_start=1)
    with pytest.raises(SequenceError) as bounds:
        sequence_from_agp((interval,), {"a": DNASequence("AAA")})
    assert bounds.value.code == "AGP_COMPONENT_INTERVAL_OUT_OF_BOUNDS"

    with pytest.raises(UnsupportedGapOperationError) as gapped:
        sequence_from_agp(
            (_component(0, 2, 1, "a"),),
            {"a": DNASequence(["A", Gap(1), "T"])},
        )
    assert gapped.value.code == "AGP_GAPPED_COMPONENT_UNSUPPORTED"

    with pytest.raises(UnsupportedGapOperationError) as circular:
        sequence_from_agp(
            (_component(0, 2, 1, "a"),),
            {"a": DNASequence("AT", topology=Topology.CIRCULAR)},
        )
    assert circular.value.code == "AGP_CIRCULAR_COMPONENT_UNSUPPORTED"


def test_agp_resources_and_output_alphabet_are_bounded() -> None:
    entry = _component(0, 2, 1, "a")
    with pytest.raises(SequenceError) as entry_limit:
        sequence_from_agp(
            itertools.repeat(entry),
            {"a": DNASequence("AA")},
            config=AGPAssemblyConfig(max_entries=2),
        )
    assert entry_limit.value.code == "AGP_ENTRY_LIMIT_EXCEEDED"

    with pytest.raises(SequenceError) as source_limit:
        sequence_from_agp(
            (entry,),
            {"a": DNASequence("AAA")},
            config=AGPAssemblyConfig(max_component_symbols=2),
        )
    assert source_limit.value.code == "AGP_COMPONENT_SYMBOL_LIMIT_EXCEEDED"

    with pytest.raises(SequenceError) as output_limit:
        sequence_from_agp(
            (entry,),
            {"a": DNASequence("AA")},
            config=AGPAssemblyConfig(max_output_symbols=1),
        )
    assert output_limit.value.code == "AGP_OUTPUT_SYMBOL_LIMIT_EXCEEDED"
    assert output_limit.value.context["requested_fragment_symbols"] == 2
    assert output_limit.value.context["remaining_output_symbols"] == 1

    with pytest.raises(SequenceError) as span_limit:
        sequence_from_agp(
            (entry,),
            {"a": DNASequence("AA")},
            config=AGPAssemblyConfig(max_output_span=1),
        )
    assert span_limit.value.code == "AGP_OUTPUT_SPAN_LIMIT_EXCEEDED"

    with pytest.raises(SequenceError) as alphabet:
        sequence_from_agp(
            (entry,),
            {"a": DNASequence("AN", alphabet=DNAAlphabet.IUPAC)},
            config=AGPAssemblyConfig(output_alphabet=DNAAlphabet.STRICT),
        )
    assert alphabet.value.code == "AGP_OUTPUT_ALPHABET_MISMATCH"

    header_document = AGPDocument((entry,), ("#one", "#two"))
    with pytest.raises(SequenceError) as header_count:
        sequence_from_agp(
            header_document,
            {"a": DNASequence("AA")},
            config=AGPAssemblyConfig(max_header_lines=1),
        )
    assert header_count.value.code == "AGP_HEADER_LIMIT_EXCEEDED"

    long_header_document = AGPDocument((entry,), ("#long",))
    with pytest.raises(SequenceError) as header_length:
        sequence_from_agp(
            long_header_document,
            {"a": DNASequence("AA")},
            config=AGPAssemblyConfig(max_header_length=4),
        )
    assert header_length.value.code == "AGP_HEADER_LENGTH_LIMIT_EXCEEDED"


def test_agp_all_gap_object_and_config_metadata_remain_json_safe() -> None:
    known = _known_gap(0, 5, 1)
    metadata = {"project": {"id": "P1"}}
    config = AGPAssemblyConfig(
        record_description="gap-only",
        record_metadata=metadata,
    )
    metadata["project"] = {"id": "changed"}

    result = sequence_from_agp((known,), {}, config=config)

    assert result.sequence.symbol_length == 0
    assert result.sequence.coordinate_span == 5
    assert result.record.description == "gap-only"
    assert result.record.metadata["project"] == {"id": "P1"}
    mutable_view = cast(dict[str, object], result.record.metadata)
    with pytest.raises(TypeError):
        mutable_view["project"] = "mutable"
    json.dumps(result.to_dict())


def test_agp_config_and_empty_input_have_stable_errors() -> None:
    with pytest.raises(ConfigurationError):
        AGPAssemblyConfig(max_entries=False)
    with pytest.raises(ConfigurationError):
        AGPAssemblyConfig(record_metadata={"dnakit_agp_assembly": {}})
    with pytest.raises(ConfigurationError) as invalid_metadata:
        AGPAssemblyConfig(record_metadata={"bad": object()})
    assert invalid_metadata.value.code == "INVALID_AGP_ASSEMBLY_METADATA"
    with pytest.raises(ConfigurationError):
        sequence_from_agp((), {}, object_id="")
    with pytest.raises(SequenceError) as bad_mapping:
        sequence_from_agp(
            (_component(0, 1, 1, "a"),),
            {1: DNASequence("A")},  # type: ignore[dict-item]
        )
    assert bad_mapping.value.code == "AGP_COMPONENT_NOT_FOUND"
    with pytest.raises(SequenceError) as empty:
        sequence_from_agp((), {})
    assert empty.value.code == "EMPTY_AGP_ASSEMBLY"
