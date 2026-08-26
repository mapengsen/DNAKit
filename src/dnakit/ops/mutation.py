"""Specified and reproducible single-site mutation operations."""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Literal, TypeAlias, cast

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.ops._common import require_sequence
from dnakit.ops.edit import Edit, substitute

MutationMode: TypeAlias = Literal["specified", "random"]
RandomSource: TypeAlias = Literal["none", "seed", "rng"]
RandomState: TypeAlias = tuple[int, tuple[int, ...], float | None]


@dataclass(frozen=True)
class MutationResult:
    """A new sequence plus a reproducible, single-site mutation audit."""

    sequence: DNASequence
    mode: MutationMode
    edit: Edit
    seed: int | None
    random_source: RandomSource = "none"
    rng_algorithm_version: int | None = None
    rng_state_before: RandomState | None = None
    rng_state_after: RandomState | None = None
    allowed_bases: str | None = None


def mutate(
    sequence: DNASequence,
    *,
    position: int | None = None,
    replacement: str | None = None,
    seed: int | None = None,
    rng: random.Random | None = None,
    allowed_bases: str = "ACGT",
) -> MutationResult:
    """Apply one specified substitution or one seeded random SNV.

    Supplying ``position`` and ``replacement`` selects specified mode.  Omitting
    both selects random mode and requires exactly one of ``seed`` or ``rng``.
    Random mode chooses a nucleotide position (gaps are never selected), then a
    different base from ``allowed_bases``.  Combinatorial mutation enumeration
    is intentionally outside the MVP.
    """

    source = require_sequence(sequence)
    specified = position is not None or replacement is not None
    if specified:
        if position is None or replacement is None:
            raise ConfigurationError(
                "Specified mutation requires both position and replacement.",
                code="INCOMPLETE_SPECIFIED_MUTATION",
            )
        if seed is not None or rng is not None:
            raise ConfigurationError(
                "seed and rng are only valid for random mutation.",
                code="UNEXPECTED_RANDOM_SOURCE",
            )
        if not isinstance(replacement, str) or len(replacement) != 1:
            raise ConfigurationError(
                "Specified mutation replacement must be one normalized nucleotide symbol.",
                code="INVALID_MUTATION_REPLACEMENT",
                context={"replacement": replacement},
            )
        result = substitute(source, position, position + 1, replacement)
        if result.edits[0].removed_symbols == replacement:
            raise ConfigurationError(
                "Specified mutation replacement must differ from the original symbol.",
                code="NO_OP_MUTATION",
                context={"position": position, "symbol": replacement},
            )
        return MutationResult(result.sequence, "specified", result.edits[0], None)

    if (seed is None) == (rng is None):
        raise ConfigurationError(
            "Random mutation requires exactly one of seed or rng.",
            code="RANDOM_SOURCE_REQUIRED",
        )
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ConfigurationError(
            "seed must be an integer.",
            code="INVALID_RANDOM_SEED",
            context={"seed": seed},
        )
    if rng is not None and not isinstance(rng, random.Random):
        raise ConfigurationError(
            "rng must be random.Random.",
            code="INVALID_RANDOM_GENERATOR",
            context={"type": type(rng).__name__},
        )
    if (
        not isinstance(allowed_bases, str)
        or not allowed_bases
        or len(set(allowed_bases)) != len(allowed_bases)
        or any(base not in "ACGT" for base in allowed_bases)
    ):
        raise ConfigurationError(
            "allowed_bases must contain unique uppercase characters from ACGT.",
            code="INVALID_ALLOWED_BASES",
            context={"allowed_bases": allowed_bases},
        )

    candidates: list[tuple[int, str]] = []
    coordinate = 0
    for part in source.parts:
        if isinstance(part, str):
            candidates.extend((coordinate + index, symbol) for index, symbol in enumerate(part))
            coordinate += len(part)
        else:
            if part.length is None:
                raise ConfigurationError(
                    "Random mutation requires known Gap lengths.",
                    code="UNKNOWN_RANDOM_MUTATION_COORDINATES",
                )
            coordinate += part.length
    viable = [
        (candidate_position, base)
        for candidate_position, base in candidates
        if any(replacement_base != base for replacement_base in allowed_bases)
    ]
    if not viable:
        raise ConfigurationError(
            "No nucleotide position has an allowed alternative base.",
            code="NO_RANDOM_MUTATION_CANDIDATE",
        )
    generator = random.Random(seed) if rng is None else rng
    state_before = cast(RandomState, generator.getstate())
    selected_position, original = generator.choice(viable)
    choices = tuple(base for base in allowed_bases if base != original)
    selected_replacement = generator.choice(choices)
    result = substitute(source, selected_position, selected_position + 1, selected_replacement)
    return MutationResult(
        result.sequence,
        "random",
        result.edits[0],
        seed,
        "seed" if rng is None else "rng",
        random.Random.VERSION,
        state_before,
        cast(RandomState, generator.getstate()),
        allowed_bases,
    )


__all__ = ["MutationResult", "RandomSource", "RandomState", "mutate"]
