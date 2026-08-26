"""Dependency-free position encodings for DNA sequences."""

from __future__ import annotations

from collections.abc import Sequence

from dnakit.core._json import FrozenDict
from dnakit.core.gap import Gap
from dnakit.exceptions import (
    ConfigurationError,
    InvalidAlphabetError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)
from dnakit.fingerprints._shared import (
    GapEncodingPolicy,
    IntegerAmbiguityPolicy,
    OneHotAmbiguityPolicy,
    SequenceInput,
    coerce_enum,
    sequence_and_id,
)
from dnakit.fingerprints.results import IntegerEncodingResult, OneHotEncodingResult

CANONICAL_BASE_ORDER = ("A", "C", "G", "T")
IUPAC_INTEGER_ORDER = (
    "A",
    "C",
    "G",
    "T",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "B",
    "D",
    "H",
    "V",
    "N",
)
INTEGER_GAP_CODE = -1
INTEGER_AMBIGUITY_CODE = -2
DEFAULT_MAX_OUTPUT_LENGTH = 10_000_000

_IUPAC_BASES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def _gap_counts(value: SequenceInput) -> tuple[int, int]:
    sequence, _ = sequence_and_id(value)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return len(gaps), sum(gap.length is None for gap in gaps)


def _raise_ambiguity(symbol: str, symbol_index: int, policy_name: str) -> None:
    raise InvalidAlphabetError(
        f"Ambiguous IUPAC symbol {symbol!r} is not allowed by this encoding policy.",
        code="ENCODING_AMBIGUITY_NOT_ALLOWED",
        context={"symbol": symbol, "symbol_index": symbol_index, "policy": policy_name},
        hint="Choose an ambiguity policy that explicitly represents or masks IUPAC symbols.",
    )


def _expand_gap(gap: Gap, gap_index: int) -> int:
    if gap.length is None:
        raise UnknownLengthError(
            "An unknown-length Gap cannot be expanded into positional encoding rows.",
            code="UNKNOWN_GAP_ENCODING_LENGTH",
            context={"gap_index": gap_index},
            hint="Use gap_policy='omit' or resolve the Gap length first.",
        )
    return gap.length


def _validate_max_output_length(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            "max_output_length must be a positive integer.",
            code="INVALID_ENCODING_SIZE_LIMIT",
            context={"max_output_length": value},
        )


def _planned_output_length(
    value: SequenceInput,
    *,
    gap_policy: GapEncodingPolicy,
) -> int:
    sequence, _ = sequence_and_id(value)
    if gap_policy is not GapEncodingPolicy.EXPAND:
        return sequence.symbol_length
    length = sequence.symbol_length
    gap_index = 0
    for part in sequence.parts:
        if isinstance(part, Gap):
            length += _expand_gap(part, gap_index)
            gap_index += 1
    return length


def integer_encode(
    value: SequenceInput,
    *,
    ambiguity_policy: IntegerAmbiguityPolicy | str = IntegerAmbiguityPolicy.IUPAC,
    gap_policy: GapEncodingPolicy | str = GapEncodingPolicy.ERROR,
    max_output_length: int = DEFAULT_MAX_OUTPUT_LENGTH,
) -> IntegerEncodingResult:
    """Encode DNA symbols with a versioned integer codebook.

    Canonical bases always map to ``A=0, C=1, G=2, T=3``. Under ``iupac``,
    the remaining IUPAC symbols map to 4 through 14 in
    :data:`IUPAC_INTEGER_ORDER`; under ``sentinel`` they all map to ``-2``.
    ``gap_policy='expand'`` emits ``-1`` once per known Gap coordinate.
    """

    sequence, sequence_id = sequence_and_id(value)
    resolved_ambiguity = coerce_enum(
        ambiguity_policy,
        IntegerAmbiguityPolicy,
        "integer ambiguity policy",
    )
    resolved_gap = coerce_enum(gap_policy, GapEncodingPolicy, "Gap encoding policy")
    _validate_max_output_length(max_output_length)
    planned_length = _planned_output_length(value, gap_policy=resolved_gap)
    if planned_length > max_output_length:
        raise ConfigurationError(
            "Integer encoding output exceeds max_output_length.",
            code="ENCODING_SIZE_LIMIT",
            context={
                "planned_output_length": planned_length,
                "max_output_length": max_output_length,
            },
            hint=(
                "Use gap_policy='omit', process a smaller sequence, or explicitly raise the limit."
            ),
        )
    symbols = (
        IUPAC_INTEGER_ORDER
        if resolved_ambiguity is IntegerAmbiguityPolicy.IUPAC
        else CANONICAL_BASE_ORDER
    )
    codebook_values: dict[str, int] = {symbol: index for index, symbol in enumerate(symbols)}
    if resolved_ambiguity is IntegerAmbiguityPolicy.SENTINEL:
        codebook_values["<IUPAC>"] = INTEGER_AMBIGUITY_CODE
    if resolved_gap is GapEncodingPolicy.EXPAND:
        codebook_values["<GAP>"] = INTEGER_GAP_CODE

    encoded: list[int] = []
    symbol_index = 0
    gap_index = 0
    ambiguity_count = 0
    expanded_gap_length = 0
    omitted_gap_count = 0
    for part in sequence.parts:
        if isinstance(part, str):
            for symbol in part:
                if symbol in "ACGT":
                    encoded.append(codebook_values[symbol])
                elif resolved_ambiguity is IntegerAmbiguityPolicy.ERROR:
                    _raise_ambiguity(symbol, symbol_index, resolved_ambiguity.value)
                elif resolved_ambiguity is IntegerAmbiguityPolicy.IUPAC:
                    encoded.append(codebook_values[symbol])
                    ambiguity_count += 1
                else:
                    encoded.append(INTEGER_AMBIGUITY_CODE)
                    ambiguity_count += 1
                symbol_index += 1
            continue

        if resolved_gap is GapEncodingPolicy.ERROR:
            raise UnsupportedGapOperationError(
                "Integer encoding does not accept explicit Gaps under the error policy.",
                code="ENCODING_GAP_NOT_ALLOWED",
                context={"gap_index": gap_index},
                hint="Choose gap_policy='omit' or 'expand'.",
            )
        if resolved_gap is GapEncodingPolicy.OMIT:
            omitted_gap_count += 1
        else:
            length = _expand_gap(part, gap_index)
            encoded.extend((INTEGER_GAP_CODE,) * length)
            expanded_gap_length += length
        gap_index += 1

    gap_count, unknown_gap_count = _gap_counts(value)
    return IntegerEncodingResult(
        name="integer",
        method="fixed_iupac_integer_encoding",
        schema_version="dnakit.integer.v1",
        sequence_id=sequence_id,
        symbol_length=sequence.symbol_length,
        gap_count=gap_count,
        unknown_gap_count=unknown_gap_count,
        values=tuple(encoded),
        codebook=FrozenDict(codebook_values),
        ambiguity_policy=resolved_ambiguity,
        gap_policy=resolved_gap,
        output_length=len(encoded),
        encoded_ambiguity_count=ambiguity_count,
        expanded_gap_length=expanded_gap_length,
        omitted_gap_count=omitted_gap_count,
        max_output_length=max_output_length,
    )


def _resolve_base_order(value: Sequence[str]) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)):
        raise ConfigurationError("base_order must be a sequence containing A, C, G, and T.")
    resolved = tuple(value)
    if (
        len(resolved) != 4
        or set(resolved) != set("ACGT")
        or any(not isinstance(base, str) for base in resolved)
    ):
        raise ConfigurationError(
            "base_order must contain each of A, C, G, and T exactly once.",
            context={"base_order": resolved},
        )
    return resolved


def one_hot_encode(
    value: SequenceInput,
    *,
    ambiguity_policy: OneHotAmbiguityPolicy | str = OneHotAmbiguityPolicy.ERROR,
    gap_policy: GapEncodingPolicy | str = GapEncodingPolicy.ERROR,
    base_order: Sequence[str] = CANONICAL_BASE_ORDER,
    max_output_length: int = DEFAULT_MAX_OUTPUT_LENGTH,
) -> OneHotEncodingResult:
    """Encode DNA as immutable A/C/G/T rows without requiring NumPy.

    ``fractional`` distributes mass uniformly over the canonical bases allowed
    by an IUPAC symbol, while ``zero`` emits an all-zero row. Known-length Gaps
    can likewise be expanded to zero rows; unknown-length Gaps cannot be
    expanded because their output dimension is undefined.
    """

    sequence, sequence_id = sequence_and_id(value)
    resolved_ambiguity = coerce_enum(
        ambiguity_policy,
        OneHotAmbiguityPolicy,
        "one-hot ambiguity policy",
    )
    resolved_gap = coerce_enum(gap_policy, GapEncodingPolicy, "Gap encoding policy")
    _validate_max_output_length(max_output_length)
    planned_length = _planned_output_length(value, gap_policy=resolved_gap)
    if planned_length > max_output_length:
        raise ConfigurationError(
            "One-hot encoding output exceeds max_output_length.",
            code="ENCODING_SIZE_LIMIT",
            context={
                "planned_output_length": planned_length,
                "max_output_length": max_output_length,
            },
            hint=(
                "Use gap_policy='omit', process a smaller sequence, or explicitly raise the limit."
            ),
        )
    resolved_order = _resolve_base_order(base_order)
    column_by_base = {base: index for index, base in enumerate(resolved_order)}
    zero_row = (0.0, 0.0, 0.0, 0.0)

    rows: list[tuple[float, ...]] = []
    symbol_index = 0
    gap_index = 0
    ambiguity_count = 0
    expanded_gap_length = 0
    omitted_gap_count = 0
    for part in sequence.parts:
        if isinstance(part, str):
            for symbol in part:
                if symbol not in "ACGT" and resolved_ambiguity is OneHotAmbiguityPolicy.ERROR:
                    _raise_ambiguity(symbol, symbol_index, resolved_ambiguity.value)
                if symbol not in "ACGT" and resolved_ambiguity is OneHotAmbiguityPolicy.ZERO:
                    rows.append(zero_row)
                    ambiguity_count += 1
                    symbol_index += 1
                    continue
                allowed_bases = _IUPAC_BASES[symbol]
                weight = 1.0 / len(allowed_bases)
                row = [0.0, 0.0, 0.0, 0.0]
                for base in allowed_bases:
                    row[column_by_base[base]] = weight
                rows.append(tuple(row))
                if symbol not in "ACGT":
                    ambiguity_count += 1
                symbol_index += 1
            continue

        if resolved_gap is GapEncodingPolicy.ERROR:
            raise UnsupportedGapOperationError(
                "One-hot encoding does not accept explicit Gaps under the error policy.",
                code="ENCODING_GAP_NOT_ALLOWED",
                context={"gap_index": gap_index},
                hint="Choose gap_policy='omit' or 'expand'.",
            )
        if resolved_gap is GapEncodingPolicy.OMIT:
            omitted_gap_count += 1
        else:
            length = _expand_gap(part, gap_index)
            rows.extend((zero_row,) * length)
            expanded_gap_length += length
        gap_index += 1

    gap_count, unknown_gap_count = _gap_counts(value)
    return OneHotEncodingResult(
        name="one_hot",
        method="iupac_acgt_one_hot_encoding",
        schema_version="dnakit.one_hot.acgt.v1",
        sequence_id=sequence_id,
        symbol_length=sequence.symbol_length,
        gap_count=gap_count,
        unknown_gap_count=unknown_gap_count,
        values=tuple(rows),
        feature_names=resolved_order,
        ambiguity_policy=resolved_ambiguity,
        gap_policy=resolved_gap,
        output_length=len(rows),
        encoded_ambiguity_count=ambiguity_count,
        expanded_gap_length=expanded_gap_length,
        omitted_gap_count=omitted_gap_count,
        max_output_length=max_output_length,
    )


__all__ = [
    "CANONICAL_BASE_ORDER",
    "DEFAULT_MAX_OUTPUT_LENGTH",
    "INTEGER_AMBIGUITY_CODE",
    "INTEGER_GAP_CODE",
    "IUPAC_INTEGER_ORDER",
    "integer_encode",
    "one_hot_encode",
]
