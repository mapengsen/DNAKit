"""Golden Gate design and reaction simulation through an optional backend."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Literal, TypeAlias

from dnakit.backends.scientific import _run_scientific_function
from dnakit.core import DNA, DNARecord, DNASequence, ProviderResult
from dnakit.exceptions import ConfigurationError
from dnakit.molbio._shared import materialize_bounded, require_sequence

GoldenGateInput: TypeAlias = DNA | DNARecord | DNASequence
GoldenGateDesignEnzyme: TypeAlias = Literal["BsaI", "BbsI"]
GoldenGateAssemblyEnzyme: TypeAlias = Literal["BsaI", "BbsI", "Esp3I", "BsmBI", "SapI"]
_DESIGN_ENZYMES = frozenset({"BsaI", "BbsI"})
_ASSEMBLY_ENZYMES = frozenset({"BsaI", "BbsI", "Esp3I", "BsmBI", "SapI"})


def _enzyme(value: str, *, assembly: bool) -> str:
    choices = _ASSEMBLY_ENZYMES if assembly else _DESIGN_ENZYMES
    if not isinstance(value, str) or value not in choices:
        raise ConfigurationError(
            "Unsupported Golden Gate enzyme.",
            code="INVALID_GOLDEN_GATE_ENZYME",
            context={"enzyme": value, "supported": tuple(sorted(choices))},
        )
    return value


def _sequences(values: Iterable[GoldenGateInput], name: str) -> tuple[str, ...]:
    items = materialize_bounded(values, max_items=64, name=name)
    if len(items) < 2:
        raise ConfigurationError(
            f"{name} must contain at least two sequences.",
            code="GOLDEN_GATE_INPUT_REQUIRED",
        )
    resolved: list[str] = []
    for index, item in enumerate(items):
        sequence = item.sequence if isinstance(item, DNARecord) else item
        resolved.append(
            require_sequence(
                sequence,
                operation=f"Golden Gate {name}[{index}]",
                max_length=1_000_000,
                canonical=True,
            )
        )
    if sum(map(len, resolved)) > 5_000_000:
        raise ConfigurationError(
            "Golden Gate inputs exceed the 5,000,000 bp total limit.",
            code="GOLDEN_GATE_TOTAL_SEQUENCE_LIMIT",
        )
    return tuple(resolved)


def _labels(value: Sequence[str] | None, count: int) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes, bytearray)):
        raise ConfigurationError(
            "labels must be a sequence of text values or None.",
            code="INVALID_GOLDEN_GATE_LABELS",
        )
    resolved = tuple(value)
    if len(resolved) != count or any(
        not isinstance(label, str) or not label.strip() or "\x00" in label for label in resolved
    ):
        raise ConfigurationError(
            "labels must contain one non-empty label per fragment.",
            code="INVALID_GOLDEN_GATE_LABELS",
        )
    return tuple(label.strip() for label in resolved)


def design_golden_gate(
    parts: Iterable[GoldenGateInput],
    *,
    enzyme: GoldenGateDesignEnzyme = "BsaI",
) -> ProviderResult:
    """Assign junction overhangs and construct Golden Gate part sequences."""

    sequences = _sequences(parts, "parts")
    resolved_enzyme = _enzyme(enzyme, assembly=False)
    return _run_scientific_function(
        "design_golden_gate",
        {"parts": list(sequences), "enzyme": resolved_enzyme},
        parameters={
            "part_count": len(sequences),
            "part_lengths": tuple(map(len, sequences)),
            "enzyme": resolved_enzyme,
        },
    )


def assemble_golden_gate(
    fragments: Iterable[GoldenGateInput],
    *,
    enzyme: GoldenGateAssemblyEnzyme = "BsaI",
    circular: bool = True,
    labels: Sequence[str] | None = None,
) -> ProviderResult:
    """Digest and assemble supplied Golden Gate reaction fragments."""

    sequences = _sequences(fragments, "fragments")
    resolved_enzyme = _enzyme(enzyme, assembly=True)
    if not isinstance(circular, bool):
        raise ConfigurationError(
            "circular must be boolean.",
            code="INVALID_GOLDEN_GATE_CIRCULAR",
        )
    resolved_labels = _labels(labels, len(sequences))
    arguments: dict[str, object] = {
        "fragments": list(sequences),
        "enzyme": resolved_enzyme,
        "circular": circular,
    }
    if resolved_labels is not None:
        arguments["labels"] = list(resolved_labels)
    return _run_scientific_function(
        "assemble_golden_gate",
        arguments,
        parameters={
            "fragment_count": len(sequences),
            "fragment_lengths": tuple(map(len, sequences)),
            "enzyme": resolved_enzyme,
            "circular": circular,
            "labels": resolved_labels,
        },
    )


__all__ = [
    "GoldenGateAssemblyEnzyme",
    "GoldenGateDesignEnzyme",
    "GoldenGateInput",
    "assemble_golden_gate",
    "design_golden_gate",
]
