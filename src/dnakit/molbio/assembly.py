"""Deterministic ordered assembly simulations with explicit method boundaries."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from itertools import pairwise

from dnakit.core import DNAAlphabet, DNASequence, Issue, IssueSeverity, Topology
from dnakit.exceptions import ConfigurationError

from ._shared import (
    freeze_parameters,
    materialize_bounded,
    native_provenance,
    require_sequence,
    validate_positive_int,
    validate_text,
)
from .restriction import LigationFragment, as_ligation_fragment, ligate_fragments
from .results import AssemblyResult, AssemblyStep, DigestFragment


@dataclass(frozen=True)
class AssemblyFragment:
    """Named linear fragment accepted by overlap-based assembly."""

    id: str
    sequence: DNASequence

    def __post_init__(self) -> None:
        validate_text(self.id, "assembly fragment id")
        require_sequence(
            self.sequence,
            operation="assembly fragment construction",
            max_length=100_000_000,
            allow_circular=False,
        )


def _overlap_fragment(
    value: DNASequence | AssemblyFragment | DigestFragment | LigationFragment,
    index: int,
) -> AssemblyFragment:
    if isinstance(value, AssemblyFragment):
        return value
    if isinstance(value, DNASequence):
        return AssemblyFragment(f"fragment_{index}", value)
    if isinstance(value, (DigestFragment, LigationFragment)):
        return AssemblyFragment(value.id, value.sequence)
    raise ConfigurationError("Unsupported assembly fragment type.")


def _longest_overlap(left: str, right: str, minimum: int, maximum: int) -> str:
    upper = min(maximum, len(left), len(right))
    for length in range(upper, minimum - 1, -1):
        if left[-length:] == right[:length]:
            return right[:length]
    return ""


def _resolve_overlap(
    left: str,
    right: str,
    *,
    explicit: str | None,
    minimum: int,
    maximum: int,
) -> str:
    if explicit is not None:
        if not isinstance(explicit, str) or set(explicit) - set("ACGTRYSWKMBDHVN"):
            raise ConfigurationError("Explicit overlaps must be uppercase DNA IUPAC strings.")
        if not minimum <= len(explicit) <= maximum:
            raise ConfigurationError("Explicit overlap length lies outside configured bounds.")
        if not left.endswith(explicit) or not right.startswith(explicit):
            raise ConfigurationError("Explicit overlap is not an exact suffix/prefix junction.")
        return explicit
    overlap = _longest_overlap(left, right, minimum, maximum)
    if not overlap:
        raise ConfigurationError(
            "No exact overlap satisfies the configured assembly boundary.",
            code="ASSEMBLY_OVERLAP_NOT_FOUND",
        )
    return overlap


def simulate_assembly(
    fragments: Iterable[DNASequence | AssemblyFragment | DigestFragment | LigationFragment],
    *,
    method: str,
    overlaps: Iterable[str] | None = None,
    min_overlap: int = 20,
    max_overlap: int = 200,
    circularize: bool = False,
    allow_blunt: bool = False,
    cycles: int = 1,
    max_fragments: int = 1_000,
    max_product_length: int = 100_000_000,
) -> AssemblyResult:
    """Simulate ordered Gibson/LCR overlaps or pre-digested GG/BioBrick ligation.

    Golden Gate and BioBrick inputs must already carry validated end descriptors;
    enzyme cutting, reaction kinetics, scar optimization, and cycle yield are not
    inferred by this function.
    """

    normalized_method = method.lower() if isinstance(method, str) else ""
    if normalized_method not in ("gibson", "lcr", "golden_gate", "biobrick"):
        raise ConfigurationError("method must be gibson, lcr, golden_gate, or biobrick.")
    validate_positive_int(max_fragments, "max_fragments", maximum=100_000)
    validate_positive_int(max_product_length, "max_product_length", maximum=100_000_000)
    validate_positive_int(cycles, "cycles", maximum=10_000)
    if not isinstance(circularize, bool) or not isinstance(allow_blunt, bool):
        raise ConfigurationError("Assembly boolean controls must be booleans.")
    materialized: list[DNASequence | AssemblyFragment | DigestFragment | LigationFragment] = []
    for item in fragments:
        if len(materialized) >= max_fragments:
            raise ConfigurationError("Assembly input exceeds max_fragments.")
        materialized.append(item)
    if len(materialized) < 2:
        raise ConfigurationError("Assembly requires at least two fragments.")
    issues = (
        Issue(
            "ASSEMBLY_KINETICS_NOT_MODELED",
            IssueSeverity.INFO,
            "The result is a deterministic sequence construction, not a reaction-yield model.",
        ),
    )
    steps: list[AssemblyStep] = []
    if normalized_method in ("golden_gate", "biobrick"):
        ligation_inputs: list[LigationFragment] = []
        for item in materialized:
            if not isinstance(item, (DigestFragment, LigationFragment)):
                raise ConfigurationError(
                    "Golden Gate/BioBrick simulation requires pre-digested fragments with ends.",
                    code="PREDIGESTED_ASSEMBLY_INPUT_REQUIRED",
                )
            ligation_inputs.append(as_ligation_fragment(item))
        ligation = ligate_fragments(
            ligation_inputs,
            circularize=circularize,
            allow_blunt=allow_blunt,
            max_fragments=max_fragments,
            max_product_length=max_product_length,
        )
        pairs = list(pairwise(ligation_inputs))
        if circularize:
            pairs.append((ligation_inputs[-1], ligation_inputs[0]))
        for index, (left, right) in enumerate(pairs, 1):
            steps.append(
                AssemblyStep(
                    index=index,
                    left_fragment_id=left.id,
                    right_fragment_id=right.id,
                    operation="compatible_end_ligation",
                    junction_sequence=left.right_end.overhang_sequence_5to3,
                )
            )
        product = ligation.product
    else:
        validate_positive_int(min_overlap, "min_overlap", maximum=100_000)
        validate_positive_int(max_overlap, "max_overlap", maximum=100_000)
        if min_overlap > max_overlap:
            raise ConfigurationError("min_overlap cannot exceed max_overlap.")
        overlap_fragments = tuple(
            _overlap_fragment(item, index) for index, item in enumerate(materialized, 1)
        )
        expected = len(overlap_fragments) if circularize else len(overlap_fragments) - 1
        overlap_values = (
            None
            if overlaps is None
            else materialize_bounded(
                overlaps,
                max_items=expected + 1,
                name="assembly overlaps",
            )
        )
        if overlap_values is not None and len(overlap_values) != expected:
            raise ConfigurationError(
                "overlaps must provide exactly one sequence per assembly junction."
            )
        if any(
            fragment.sequence.strandedness is not overlap_fragments[0].sequence.strandedness
            for fragment in overlap_fragments[1:]
        ):
            raise ConfigurationError("All assembly fragments must use the same strandedness.")
        product_symbols = overlap_fragments[0].sequence.symbols
        if len(product_symbols) > max_product_length:
            raise ConfigurationError("Assembly product exceeds max_product_length.")
        for index, (overlap_left, overlap_right) in enumerate(pairwise(overlap_fragments), 1):
            overlap = _resolve_overlap(
                overlap_left.sequence.symbols,
                overlap_right.sequence.symbols,
                explicit=None if overlap_values is None else overlap_values[index - 1],
                minimum=min_overlap,
                maximum=max_overlap,
            )
            addition = overlap_right.sequence.symbols[len(overlap) :]
            if len(product_symbols) + len(addition) > max_product_length:
                raise ConfigurationError("Assembly product exceeds max_product_length.")
            product_symbols += addition
            steps.append(
                AssemblyStep(
                    index=index,
                    left_fragment_id=overlap_left.id,
                    right_fragment_id=overlap_right.id,
                    operation="exact_overlap_merge",
                    junction_sequence=overlap,
                )
            )
        if circularize:
            closing = _resolve_overlap(
                overlap_fragments[-1].sequence.symbols,
                overlap_fragments[0].sequence.symbols,
                explicit=None if overlap_values is None else overlap_values[-1],
                minimum=min_overlap,
                maximum=max_overlap,
            )
            if not product_symbols.endswith(closing):
                raise ConfigurationError(
                    "Closing overlap was consumed or altered by previous junctions."
                )
            product_symbols = product_symbols[: -len(closing)]
            steps.append(
                AssemblyStep(
                    index=len(steps) + 1,
                    left_fragment_id=overlap_fragments[-1].id,
                    right_fragment_id=overlap_fragments[0].id,
                    operation="circular_exact_overlap_merge",
                    junction_sequence=closing,
                )
            )
        if len(product_symbols) > max_product_length:
            raise ConfigurationError("Assembly product exceeds max_product_length.")
        alphabet = (
            DNAAlphabet.IUPAC
            if any(item.sequence.alphabet is DNAAlphabet.IUPAC for item in overlap_fragments)
            else DNAAlphabet.STRICT
        )
        product = DNASequence(
            product_symbols,
            alphabet=alphabet,
            topology=Topology.CIRCULAR if circularize else Topology.LINEAR,
            strandedness=overlap_fragments[0].sequence.strandedness,
        )
    return AssemblyResult(
        assembly_method=normalized_method,
        product=product,
        steps=tuple(steps),
        circularized=circularize,
        complete=True,
        method="ordered_bounded_assembly_simulation",
        algorithm_version="dnakit-assembly-v1",
        parameters=freeze_parameters(
            {
                "assembly_method": normalized_method,
                "fragment_count": len(materialized),
                "min_overlap": min_overlap if normalized_method in ("gibson", "lcr") else None,
                "max_overlap": max_overlap if normalized_method in ("gibson", "lcr") else None,
                "explicit_overlaps": overlaps is not None,
                "circularize": circularize,
                "allow_blunt": allow_blunt,
                "cycles": cycles,
                "cycles_affect_sequence_result": False,
                "kinetics_modeled": False,
                "max_fragments": max_fragments,
                "max_product_length": max_product_length,
            }
        ),
        provenance=native_provenance(
            reimplementation=True,
            reference_name=f"{normalized_method} ordered assembly abstraction",
            reference_version="boundary-v1",
        ),
        issues=issues,
    )


__all__ = ["AssemblyFragment", "simulate_assembly"]
