"""Public immutable operations with unified sequence/record method names."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TypeAlias, overload

from dnakit.core import DNA, DNARecord, DNASequence
from dnakit.ops.circular import (
    CircularOperationResult,
)
from dnakit.ops.circular import (
    canonical_origin as _canonical_origin_sequence,
)
from dnakit.ops.circular import (
    rotate as _rotate_sequence,
)
from dnakit.ops.concat import concat, concat_overlap
from dnakit.ops.direction import (
    complement,
    reverse,
)
from dnakit.ops.direction import (
    reverse_complement as _reverse_complement_sequence,
)
from dnakit.ops.edit import (
    Edit,
    EditResult,
    circular_subsequence,
    subsequence,
)
from dnakit.ops.edit import (
    delete as _delete_sequence,
)
from dnakit.ops.edit import (
    insert as _insert_sequence,
)
from dnakit.ops.edit import (
    mask as _mask_sequence,
)
from dnakit.ops.edit import (
    substitute as _substitute_sequence,
)
from dnakit.ops.edit import (
    trim as _trim_sequence,
)
from dnakit.ops.evolution import (
    DEFAULT_EVOLUTION_AUGMENTATIONS,
    CrossoverResult,
    EvolutionAugmentation,
    EvolutionGenerationResult,
    EvolutionStep,
    IndelOperation,
    KmerShuffleResult,
    RearrangementOperation,
    RearrangementResult,
    crossover,
    evolution_generate,
    indel_generate,
    kmer_shuffle,
    rearrange_generate,
)
from dnakit.ops.mutation import MutationResult, mutate
from dnakit.ops.records import (
    FeatureChange,
    FeatureChangeAction,
    FeatureOverlapPolicy,
    LetterAnnotationAction,
    LetterAnnotationPolicy,
    RecordOperationResult,
    canonical_origin_record,
    delete_record,
    insert_record,
    mask_record,
    reverse_complement_record,
    rotate_record,
    substitute_record,
    trim_record,
)
from dnakit.ops.translation import transcribe, translate

OperationInput: TypeAlias = DNA | DNARecord | DNASequence
FragmentInput: TypeAlias = DNA | DNARecord | DNASequence | str
OperationResult: TypeAlias = DNA | EditResult | CircularOperationResult | RecordOperationResult


def _record_input(value: OperationInput) -> DNARecord | None:
    if isinstance(value, DNA):
        return value.record
    if isinstance(value, DNARecord):
        return value
    return None


def _fragment_input(value: FragmentInput) -> DNASequence | str:
    if isinstance(value, DNA):
        return value.sequence
    if isinstance(value, DNARecord):
        return value.sequence
    return value


def _sequence_input(value: OperationInput) -> DNASequence:
    if isinstance(value, DNASequence):
        return value
    raise TypeError("value must be DNA, DNARecord, or DNASequence.")


def _public_record_result(
    value: OperationInput, result: RecordOperationResult
) -> DNA | RecordOperationResult:
    if not isinstance(value, DNA):
        return result
    return DNA(
        result.record,
        name=value.name,
        source=value.source,
        version=value.version,
        collection_metadata=value.collection_metadata,
    )


@overload
def insert(
    value: DNA,
    position: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA: ...


@overload
def insert(
    value: DNARecord,
    position: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult: ...


@overload
def insert(
    value: DNASequence,
    position: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> EditResult: ...


def insert(
    value: OperationInput,
    position: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA | EditResult | RecordOperationResult:
    """Insert into a sequence or annotated record using one public name."""

    record = _record_input(value)
    resolved_fragment = _fragment_input(fragment)
    if record is not None:
        result = insert_record(
            record,
            position,
            resolved_fragment,
            feature_policy=feature_policy,
            replacement_annotations=replacement_annotations,
            letter_annotation_policy=letter_annotation_policy,
        )
        return _public_record_result(value, result)
    return _insert_sequence(_sequence_input(value), position, resolved_fragment)


@overload
def delete(
    value: DNA,
    start: int,
    end: int,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA: ...


@overload
def delete(
    value: DNARecord,
    start: int,
    end: int,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult: ...


@overload
def delete(
    value: DNASequence,
    start: int,
    end: int,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> EditResult: ...


def delete(
    value: OperationInput,
    start: int,
    end: int,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA | EditResult | RecordOperationResult:
    """Delete from a sequence or annotated record using one public name."""

    record = _record_input(value)
    if record is not None:
        result = delete_record(
            record,
            start,
            end,
            feature_policy=feature_policy,
            letter_annotation_policy=letter_annotation_policy,
        )
        return _public_record_result(value, result)
    return _delete_sequence(_sequence_input(value), start, end)


@overload
def substitute(
    value: DNA,
    start: int,
    end: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA: ...


@overload
def substitute(
    value: DNARecord,
    start: int,
    end: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult: ...


@overload
def substitute(
    value: DNASequence,
    start: int,
    end: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> EditResult: ...


def substitute(
    value: OperationInput,
    start: int,
    end: int,
    fragment: FragmentInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
    replacement_annotations: Mapping[str, Iterable[int | float]] | None = None,
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA | EditResult | RecordOperationResult:
    """Substitute sequence or record content using one public name."""

    record = _record_input(value)
    resolved_fragment = _fragment_input(fragment)
    if record is not None:
        result = substitute_record(
            record,
            start,
            end,
            resolved_fragment,
            feature_policy=feature_policy,
            replacement_annotations=replacement_annotations,
            letter_annotation_policy=letter_annotation_policy,
        )
        return _public_record_result(value, result)
    return _substitute_sequence(_sequence_input(value), start, end, resolved_fragment)


@overload
def mask(
    value: DNA,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA: ...


@overload
def mask(
    value: DNARecord,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult: ...


@overload
def mask(
    value: DNASequence,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> EditResult: ...


def mask(
    value: OperationInput,
    intervals: Iterable[tuple[int, int]],
    *,
    symbol: str = "N",
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA | EditResult | RecordOperationResult:
    """Mask a sequence or annotated record using one public name."""

    record = _record_input(value)
    if record is not None:
        result = mask_record(
            record,
            intervals,
            symbol=symbol,
            feature_policy=feature_policy,
            letter_annotation_policy=letter_annotation_policy,
        )
        return _public_record_result(value, result)
    return _mask_sequence(_sequence_input(value), intervals, symbol=symbol)


@overload
def trim(
    value: DNA,
    *,
    left: int = 0,
    right: int = 0,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA: ...


@overload
def trim(
    value: DNARecord,
    *,
    left: int = 0,
    right: int = 0,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> RecordOperationResult: ...


@overload
def trim(
    value: DNASequence,
    *,
    left: int = 0,
    right: int = 0,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> EditResult: ...


def trim(
    value: OperationInput,
    *,
    left: int = 0,
    right: int = 0,
    feature_policy: FeatureOverlapPolicy = "preserve",
    letter_annotation_policy: LetterAnnotationPolicy = "error",
) -> DNA | EditResult | RecordOperationResult:
    """Trim a sequence or annotated record using one public name."""

    record = _record_input(value)
    if record is not None:
        result = trim_record(
            record,
            left=left,
            right=right,
            feature_policy=feature_policy,
            letter_annotation_policy=letter_annotation_policy,
        )
        return _public_record_result(value, result)
    return _trim_sequence(_sequence_input(value), left=left, right=right)


@overload
def reverse_complement(
    value: DNASequence,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
) -> DNASequence: ...


@overload
def reverse_complement(
    value: DNA,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
) -> DNA: ...


@overload
def reverse_complement(
    value: DNARecord,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
) -> RecordOperationResult: ...


def reverse_complement(
    value: OperationInput,
    *,
    feature_policy: FeatureOverlapPolicy = "preserve",
) -> DNA | DNASequence | RecordOperationResult:
    """Reverse-complement sequence or record content using one public name."""

    record = _record_input(value)
    if record is not None:
        result = reverse_complement_record(record, feature_policy=feature_policy)
        return _public_record_result(value, result)
    return _reverse_complement_sequence(_sequence_input(value))


@overload
def rotate(
    value: DNASequence,
    offset: int,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> CircularOperationResult: ...


@overload
def rotate(
    value: DNA,
    offset: int,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> DNA: ...


@overload
def rotate(
    value: DNARecord,
    offset: int,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> RecordOperationResult: ...


def rotate(
    value: OperationInput,
    offset: int,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> DNA | CircularOperationResult | RecordOperationResult:
    """Rotate a circular sequence or record using one public name."""

    record = _record_input(value)
    if record is not None:
        result = rotate_record(record, offset, feature_policy=feature_policy)
        return _public_record_result(value, result)
    return _rotate_sequence(_sequence_input(value), offset)


@overload
def canonical_origin(
    value: DNASequence,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> CircularOperationResult: ...


@overload
def canonical_origin(
    value: DNA,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> DNA: ...


@overload
def canonical_origin(
    value: DNARecord,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> RecordOperationResult: ...


def canonical_origin(
    value: OperationInput,
    *,
    feature_policy: FeatureOverlapPolicy = "split",
) -> DNA | CircularOperationResult | RecordOperationResult:
    """Canonicalize a circular sequence or record using one public name."""

    record = _record_input(value)
    if record is not None:
        result = canonical_origin_record(record, feature_policy=feature_policy)
        return _public_record_result(value, result)
    return _canonical_origin_sequence(_sequence_input(value))


__all__ = [
    "DEFAULT_EVOLUTION_AUGMENTATIONS",
    "CircularOperationResult",
    "CrossoverResult",
    "Edit",
    "EditResult",
    "EvolutionAugmentation",
    "EvolutionGenerationResult",
    "EvolutionStep",
    "FeatureChange",
    "FeatureChangeAction",
    "FeatureOverlapPolicy",
    "IndelOperation",
    "KmerShuffleResult",
    "LetterAnnotationAction",
    "LetterAnnotationPolicy",
    "MutationResult",
    "RearrangementOperation",
    "RearrangementResult",
    "RecordOperationResult",
    "canonical_origin",
    "canonical_origin_record",
    "circular_subsequence",
    "complement",
    "concat",
    "concat_overlap",
    "crossover",
    "delete",
    "delete_record",
    "evolution_generate",
    "indel_generate",
    "insert",
    "insert_record",
    "kmer_shuffle",
    "mask",
    "mask_record",
    "mutate",
    "rearrange_generate",
    "reverse",
    "reverse_complement",
    "reverse_complement_record",
    "rotate",
    "rotate_record",
    "subsequence",
    "substitute",
    "substitute_record",
    "transcribe",
    "translate",
    "trim",
    "trim_record",
]
