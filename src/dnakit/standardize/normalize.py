"""Deterministic DNA input normalization with a complete in-memory audit."""

from __future__ import annotations

import hashlib
import json
import unicodedata
from collections import defaultdict
from collections.abc import Iterable
from typing import TypeAlias

from dnakit.core._json import to_json_compatible
from dnakit.core.enums import DNAAlphabet, IssueSeverity
from dnakit.core.gap import Gap
from dnakit.core.issues import Issue
from dnakit.core.provenance import Provenance
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError

from ._shared import CANONICAL_BASES, IUPAC_AMBIGUITY, IUPAC_BASES, ambiguity_report
from .config import AmbiguityPolicy, NormalizationConfig, UPolicy
from .results import (
    InputPosition,
    InvalidSymbol,
    NormalizationChange,
    NormalizationResult,
    NormalizationStep,
    RawInputSnapshot,
)

RawSequenceInput: TypeAlias = str | bytes | Iterable[str | Gap] | DNASequence


def normalize(
    raw: RawSequenceInput,
    *,
    keep_ambiguous: bool = True,
    keep_u: bool = False,
    keep_other: bool = False,
    config: NormalizationConfig | None = None,
) -> NormalizationResult:
    """Normalize raw DNA while retaining an audit entry for every edit.

    By default IUPAC ambiguity symbols are retained, while U and all other
    non-DNA characters are deleted.  Set the three ``keep_*`` arguments to
    change those simple policies.  Advanced callers can pass
    :class:`NormalizationConfig`; non-default keep arguments must not be mixed
    with an explicit config.  Configuration errors and undecodable bytes raise
    immediately.
    """
    simple_flags = {
        "keep_ambiguous": keep_ambiguous,
        "keep_u": keep_u,
        "keep_other": keep_other,
    }
    for field_name, value in simple_flags.items():
        if not isinstance(value, bool):
            raise ConfigurationError(
                f"{field_name} must be a boolean.",
                code="INVALID_NORMALIZATION_FLAG",
                context={"field": field_name},
            )
    if config is None:
        resolved = NormalizationConfig(
            keep_ambiguous=keep_ambiguous,
            keep_u=keep_u,
            keep_other=keep_other,
        )
    elif not isinstance(config, NormalizationConfig):
        raise ConfigurationError(
            "config must be NormalizationConfig or None.",
            code="INVALID_NORMALIZATION_CONFIG",
        )
    else:
        if (keep_ambiguous, keep_u, keep_other) != (True, False, False):
            raise ConfigurationError(
                "Non-default keep arguments cannot be combined with config.",
                code="NORMALIZATION_ARGUMENT_CONFLICT",
                hint="Set keep_ambiguous, keep_u and keep_other on NormalizationConfig.",
            )
        resolved = config

    snapshot, source_parts = _snapshot_and_parts(raw)
    source_sequence = raw if isinstance(raw, DNASequence) else None
    u_policy = resolved.effective_u_policy
    ambiguity_policy = resolved.effective_ambiguity_policy

    output_parts: list[str | Gap] = []
    output_positions: list[tuple[InputPosition, ...] | None] = []
    changes: list[NormalizationChange] = []
    issues: list[Issue] = []
    u_positions: list[InputPosition] = []
    gap_count = 0
    normalized_symbol_offset = 0
    absolute_offset = 0

    for part_index, part in enumerate(source_parts):
        if isinstance(part, Gap):
            gap_count += 1
            output_parts.append(part)
            output_positions.append(None)
            if not resolved.allow_gaps:
                issues.append(
                    _issue(
                        "STD_GAP_NOT_ALLOWED",
                        IssueSeverity.ERROR,
                        "The input contains an explicit Gap but gaps are disabled.",
                        {"part_index": part_index},
                    )
                )
            continue

        normalized_characters: list[str] = []
        normalized_positions: list[InputPosition] = []
        for offset, original in enumerate(part):
            position = InputPosition(part_index, offset, absolute_offset + offset)
            symbol = original

            if symbol.isspace():
                if resolved.remove_whitespace:
                    changes.append(
                        _change(
                            "remove_whitespace",
                            position,
                            symbol,
                            "",
                            normalized_symbol_offset + len(normalized_characters),
                            "Configured whitespace removal.",
                        )
                    )
                    continue
            elif _is_invisible(symbol) and resolved.remove_invisible:
                changes.append(
                    _change(
                        "remove_invisible",
                        position,
                        symbol,
                        "",
                        normalized_symbol_offset + len(normalized_characters),
                        f"Removed Unicode category {unicodedata.category(symbol)}.",
                    )
                )
                continue
            elif symbol in resolved.removable_separators:
                changes.append(
                    _change(
                        "remove_separator",
                        position,
                        symbol,
                        "",
                        normalized_symbol_offset + len(normalized_characters),
                        "Configured separator removal.",
                    )
                )
                continue

            if resolved.uppercase and "a" <= symbol <= "z":
                upper = symbol.upper()
                changes.append(
                    _change(
                        "uppercase",
                        position,
                        symbol,
                        upper,
                        normalized_symbol_offset + len(normalized_characters),
                        "ASCII DNA characters use uppercase canonical form.",
                    )
                )
                symbol = upper

            if symbol == "U":
                u_positions.append(position)
                if u_policy is UPolicy.DELETE:
                    changes.append(
                        _change(
                            "delete_u",
                            position,
                            symbol,
                            "",
                            normalized_symbol_offset + len(normalized_characters),
                            "Configured DNA uracil deletion.",
                        )
                    )
                    continue
                if u_policy is UPolicy.REPLACE:
                    changes.append(
                        _change(
                            "replace_u",
                            position,
                            symbol,
                            "T",
                            normalized_symbol_offset + len(normalized_characters),
                            "Configured DNA uracil-to-thymine conversion.",
                        )
                    )
                    symbol = "T"

            if symbol in IUPAC_AMBIGUITY:
                if ambiguity_policy is AmbiguityPolicy.DELETE:
                    changes.append(
                        _change(
                            "delete_ambiguity",
                            position,
                            symbol,
                            "",
                            normalized_symbol_offset + len(normalized_characters),
                            "Configured IUPAC ambiguity deletion.",
                        )
                    )
                    continue
                if ambiguity_policy is AmbiguityPolicy.MASK and symbol != resolved.ambiguity_mask:
                    changes.append(
                        _change(
                            "mask_ambiguity",
                            position,
                            symbol,
                            resolved.ambiguity_mask,
                            normalized_symbol_offset + len(normalized_characters),
                            "Configured ambiguity masking.",
                        )
                    )
                    symbol = resolved.ambiguity_mask

            if symbol not in IUPAC_BASES and symbol != "U" and not resolved.keep_other:
                changes.append(
                    _change(
                        "delete_other",
                        position,
                        symbol,
                        "",
                        normalized_symbol_offset + len(normalized_characters),
                        "Configured deletion of a non-DNA character.",
                    )
                )
                continue

            normalized_characters.append(symbol)
            normalized_positions.append(position)

        output_parts.append("".join(normalized_characters))
        output_positions.append(tuple(normalized_positions))
        normalized_symbol_offset += len(normalized_characters)
        absolute_offset += len(part)

    normalized_parts = _coalesce_parts(tuple(output_parts))
    alphabet = _resolve_output_alphabet(resolved.alphabet, source_sequence, output_parts)
    invalid = _invalid_symbols(output_parts, output_positions, alphabet)

    if u_positions and u_policy in {UPolicy.ERROR, UPolicy.WARN}:
        severity = IssueSeverity.ERROR if u_policy is UPolicy.ERROR else IssueSeverity.WARNING
        issues.append(
            _issue(
                "STD_U_PRESENT",
                severity,
                f"Found {len(u_positions)} uracil symbol(s); policy is {u_policy.value}.",
                {"positions": [position.absolute_offset for position in u_positions]},
            )
        )

    ambiguity = ambiguity_report(
        normalized_parts,
        base_priors=(
            resolved.base_priors if ambiguity_policy is AmbiguityPolicy.PROBABILITY else None
        ),
    )
    if ambiguity.total_count and ambiguity_policy is AmbiguityPolicy.ERROR:
        issues.append(
            _issue(
                "STD_AMBIGUITY_NOT_ALLOWED",
                IssueSeverity.ERROR,
                f"Found {ambiguity.total_count} IUPAC ambiguity symbol(s).",
                {"symbols": {item.symbol: item.count for item in ambiguity.by_symbol}},
            )
        )

    for item in invalid:
        issues.append(
            _issue(
                "STD_INVALID_SYMBOL",
                IssueSeverity.ERROR,
                f"Character {item.symbol!r} ({item.codepoint}) is not valid for "
                f"the {alphabet.value} DNA alphabet.",
                {
                    "symbol": item.symbol,
                    "codepoint": item.codepoint,
                    "positions": [position.absolute_offset for position in item.positions],
                },
            )
        )

    has_error = any(issue.severity is IssueSeverity.ERROR for issue in issues)
    sequence: DNASequence | None = None
    if not has_error:
        if (
            source_sequence is not None
            and normalized_parts == source_sequence.parts
            and alphabet is source_sequence.alphabet
        ):
            sequence = source_sequence
        else:
            if source_sequence is None:
                sequence = DNASequence(normalized_parts, alphabet=alphabet)
            else:
                sequence = DNASequence(
                    normalized_parts,
                    alphabet=alphabet,
                    topology=source_sequence.topology,
                    strandedness=source_sequence.strandedness,
                )

    steps = _steps(
        resolved,
        alphabet,
        u_policy,
        ambiguity_policy,
        changes,
        invalid,
        gap_count,
        has_error,
    )
    result = NormalizationResult(
        raw_input=snapshot,
        sequence=sequence,
        config=resolved,
        algorithm_version="std-normalize-v2",
        provenance=Provenance(),
        normalized_parts=normalized_parts,
        steps=steps,
        changes=tuple(changes),
        issues=tuple(issues),
        ambiguity=ambiguity,
        invalid_symbols=invalid,
        u_positions=tuple(u_positions),
    )
    if has_error and resolved.raise_on_error:
        raise InvalidAlphabetError(
            "Sequence normalization failed alphabet validation.",
            context={"issue_codes": [issue.code for issue in issues]},
            hint="Inspect NormalizationConfig and the reported invalid input characters.",
        )
    return result


def _snapshot_and_parts(
    raw: RawSequenceInput,
) -> tuple[RawInputSnapshot, tuple[str | Gap, ...]]:
    if isinstance(raw, DNASequence):
        parts = raw.parts
        payload = _parts_payload(parts, raw)
        return (
            RawInputSnapshot(
                input_type="DNASequence",
                content=raw,
                sha256=hashlib.sha256(payload).hexdigest(),
                character_count=sum(len(part) for part in parts if isinstance(part, str)),
            ),
            parts,
        )
    if isinstance(raw, str):
        payload = b"str\0" + raw.encode("utf-8")
        return (
            RawInputSnapshot("str", raw, hashlib.sha256(payload).hexdigest(), len(raw)),
            (raw,),
        )
    if isinstance(raw, bytes):
        text = raw.decode("utf-8", errors="strict")
        payload = b"bytes\0" + raw
        return (
            RawInputSnapshot("bytes", raw, hashlib.sha256(payload).hexdigest(), len(text)),
            (text,),
        )
    if not isinstance(raw, Iterable):
        raise TypeError("raw must be str, bytes, DNASequence or an iterable of str/Gap parts.")

    materialized = tuple(raw)
    for index, part in enumerate(materialized):
        if not isinstance(part, (str, Gap)):
            raise TypeError(
                f"raw part {index} has type {type(part).__name__}; expected str or Gap."
            )
    parts = materialized
    payload = _parts_payload(parts)
    return (
        RawInputSnapshot(
            "parts",
            parts,
            hashlib.sha256(payload).hexdigest(),
            sum(len(part) for part in parts if isinstance(part, str)),
        ),
        parts,
    )


def _parts_payload(parts: tuple[str | Gap, ...], sequence: DNASequence | None = None) -> bytes:
    encoded_parts: list[dict[str, object]] = []
    for part in parts:
        if isinstance(part, str):
            encoded_parts.append({"sequence": part})
        else:
            encoded_parts.append(
                {
                    "gap": {
                        "length": part.length,
                        "kind": getattr(part.kind, "value", part.kind),
                        "crossable": part.crossable,
                        "evidence": list(part.evidence),
                        "metadata": to_json_compatible(part.metadata),
                    }
                }
            )
    payload: dict[str, object] = {"parts": encoded_parts}
    if sequence is not None:
        payload.update(
            alphabet=sequence.alphabet.value,
            topology=sequence.topology.value,
            strandedness=sequence.strandedness.value,
        )
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _resolve_output_alphabet(
    configured: DNAAlphabet | None,
    source_sequence: DNASequence | None,
    parts: list[str | Gap],
) -> DNAAlphabet:
    if configured is not None:
        return configured
    if source_sequence is not None:
        return source_sequence.alphabet
    if any(symbol in IUPAC_AMBIGUITY for part in parts if isinstance(part, str) for symbol in part):
        return DNAAlphabet.IUPAC
    return DNAAlphabet.STRICT


def _invalid_symbols(
    parts: list[str | Gap],
    positions: list[tuple[InputPosition, ...] | None],
    alphabet: DNAAlphabet,
) -> tuple[InvalidSymbol, ...]:
    allowed = CANONICAL_BASES if alphabet is DNAAlphabet.STRICT else IUPAC_BASES
    grouped: dict[str, list[InputPosition]] = defaultdict(list)
    for part, part_positions in zip(parts, positions, strict=True):
        if isinstance(part, Gap):
            continue
        assert part_positions is not None
        for symbol, position in zip(part, part_positions, strict=True):
            if symbol not in allowed:
                grouped[symbol].append(position)
    return tuple(
        InvalidSymbol(
            symbol=symbol,
            codepoint=f"U+{ord(symbol):04X}",
            positions=tuple(items),
        )
        for symbol, items in sorted(grouped.items(), key=lambda item: ord(item[0]))
    )


def _coalesce_parts(parts: tuple[str | Gap, ...]) -> tuple[str | Gap, ...]:
    result: list[str | Gap] = []
    pending: list[str] = []
    for part in parts:
        if isinstance(part, str):
            pending.append(part)
            continue
        if pending:
            result.append("".join(pending))
            pending.clear()
        result.append(part)
    if pending:
        result.append("".join(pending))
    return tuple(result)


def _steps(
    config: NormalizationConfig,
    alphabet: DNAAlphabet,
    u_policy: UPolicy,
    ambiguity_policy: AmbiguityPolicy,
    changes: list[NormalizationChange],
    invalid: tuple[InvalidSymbol, ...],
    gap_count: int,
    failed: bool,
) -> tuple[NormalizationStep, ...]:
    counts: dict[str, int] = defaultdict(int)
    for change in changes:
        counts[change.operation] += 1
    return (
        _step("whitespace", counts["remove_whitespace"], ("enabled", config.remove_whitespace)),
        _step("invisible", counts["remove_invisible"], ("enabled", config.remove_invisible)),
        _step(
            "separators",
            counts["remove_separator"],
            ("count", len(config.removable_separators)),
        ),
        _step("case", counts["uppercase"], ("uppercase", config.uppercase)),
        _step(
            "uracil",
            counts["replace_u"] + counts["delete_u"],
            ("policy", u_policy.value),
        ),
        _step(
            "ambiguity",
            counts["mask_ambiguity"] + counts["delete_ambiguity"],
            ("policy", ambiguity_policy.value),
        ),
        _step("other_characters", counts["delete_other"], ("keep", config.keep_other)),
        NormalizationStep(
            name="explicit_gap",
            status="applied" if gap_count else "no-op",
            change_count=0,
            parameters=(("count", str(gap_count)), ("allowed", str(config.allow_gaps).lower())),
        ),
        NormalizationStep(
            name="alphabet",
            status="failed" if invalid else "no-op",
            change_count=0,
            parameters=(("mode", alphabet.value),),
        ),
        NormalizationStep(
            name="construct",
            status="failed" if failed else "applied",
            change_count=0,
        ),
    )


def _step(name: str, count: int, parameter: tuple[str, object]) -> NormalizationStep:
    return NormalizationStep(
        name=name,
        status="applied" if count else "no-op",
        change_count=count,
        parameters=((parameter[0], str(parameter[1]).lower()),),
    )


def _is_invisible(symbol: str) -> bool:
    return unicodedata.category(symbol) in {"Cc", "Cf"}


def _change(
    operation: str,
    position: InputPosition,
    before: str,
    after: str,
    normalized_offset: int,
    reason: str,
) -> NormalizationChange:
    return NormalizationChange(
        operation=operation,
        position=position,
        before=before,
        after=after,
        normalized_offset=normalized_offset,
        reason=reason,
    )


def _issue(
    code: str,
    severity: IssueSeverity,
    message: str,
    details: dict[str, object],
) -> Issue:
    return Issue(code=code, severity=severity, message=message, details=details)
