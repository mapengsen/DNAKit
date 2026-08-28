"""Reproducible DNA sequence generation and structural variants."""

from __future__ import annotations

import math
import random
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import Literal, TypeAlias, cast

from dnakit.core import DNASequence, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.ops._common import require_sequence
from dnakit.ops.mutation import RandomSource, RandomState

EvolutionAugmentation: TypeAlias = Literal[
    "mutation",
    "deletion",
    "translocation",
    "insertion",
    "inversion",
    "reverse_complement",
]
IndelOperation: TypeAlias = Literal["insertion", "deletion"]
RearrangementOperation: TypeAlias = Literal["exchange", "inversion", "duplication"]

DEFAULT_EVOLUTION_AUGMENTATIONS: tuple[EvolutionAugmentation, ...] = (
    "deletion",
    "insertion",
    "translocation",
    "inversion",
    "reverse_complement",
    "mutation",
)

_EVOLUTION_PRIORITY: tuple[EvolutionAugmentation, ...] = (
    "inversion",
    "deletion",
    "translocation",
    "insertion",
    "reverse_complement",
    "mutation",
)
_EVOLUTION_PRIORITY_INDEX = {name: index for index, name in enumerate(_EVOLUTION_PRIORITY)}
_CANONICAL_BASES = "ACGT"
_MUTATION_ALTERNATIVES = {
    base: tuple(candidate for candidate in _CANONICAL_BASES if candidate != base)
    for base in _CANONICAL_BASES
}


@dataclass(frozen=True, slots=True)
class EvolutionStep:
    """Audit details for one selected EvoAug operation."""

    augmentation: EvolutionAugmentation
    applied: bool
    start: int | None = None
    end: int | None = None
    length: int | None = None
    shift: int | None = None
    mutation_attempts: int | None = None


@dataclass(frozen=True, slots=True)
class EvolutionGenerationResult:
    """Generated sequence, selected operations, and reproducibility metadata."""

    sequence: DNASequence
    steps: tuple[EvolutionStep, ...]
    seed: int | None
    random_source: RandomSource
    rng_algorithm_version: int
    rng_state_before: RandomState
    rng_state_after: RandomState
    algorithm_version: str = "dnakit-evoaug-v3"

    @property
    def augmentations(self) -> tuple[EvolutionAugmentation, ...]:
        """Return the selected operations in their applied priority order."""

        return tuple(step.augmentation for step in self.steps)


@dataclass(frozen=True, slots=True)
class RearrangementResult:
    """A rearranged sequence plus the generated segment layout audit."""

    sequence: DNASequence
    operation: RearrangementOperation
    breakpoints: tuple[int, ...]
    permutation: tuple[int, ...] | None
    selected_segment: int | None
    seed: int | None
    random_source: RandomSource
    rng_algorithm_version: int
    rng_state_before: RandomState
    rng_state_after: RandomState
    algorithm_version: str = "dnakit-rearrangement-v1"


@dataclass(frozen=True, slots=True)
class KmerShuffleResult:
    """A k-mer-count-preserving shuffle plus its reproducibility audit."""

    sequence: DNASequence
    k: int
    kmer_counts: tuple[tuple[str, int], ...]
    attempts: int
    seed: int | None
    random_source: RandomSource
    rng_algorithm_version: int
    rng_state_before: RandomState
    rng_state_after: RandomState
    algorithm_version: str = "dnakit-kmer-shuffle-v1"


@dataclass(frozen=True, slots=True)
class CrossoverResult:
    """A one-point crossover child plus its parent-coordinate audit."""

    sequence: DNASequence
    position: int
    first_parent_length: int
    second_parent_length: int
    seed: int | None
    random_source: RandomSource = "none"
    rng_algorithm_version: int | None = None
    rng_state_before: RandomState | None = None
    rng_state_after: RandomState | None = None
    algorithm_version: str = "dnakit-crossover-v1"


@dataclass(frozen=True, slots=True)
class _EvolutionParameters:
    mut_frac: float
    insert_frac: float
    delete_frac: float
    insert_min: int
    insert_max: int
    shift_min: int
    shift_max: int
    invert_min: int
    invert_max: int
    rc_prob: float


def _validate_nonnegative_int(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(
            f"{name} must be a non-negative integer.",
            code="INVALID_EVOLUTION_PARAMETER",
            context={name: value},
        )
    return value


def _validate_range(name: str, minimum: int, maximum: int) -> tuple[int, int]:
    resolved_minimum = _validate_nonnegative_int(f"{name}_min", minimum)
    resolved_maximum = _validate_nonnegative_int(f"{name}_max", maximum)
    if resolved_minimum > resolved_maximum:
        raise ConfigurationError(
            f"{name}_min cannot exceed {name}_max.",
            code="INVALID_EVOLUTION_PARAMETER_RANGE",
            context={f"{name}_min": minimum, f"{name}_max": maximum},
        )
    return resolved_minimum, resolved_maximum


def _validate_probability(name: str, value: float) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or not 0 <= value <= 1
    ):
        raise ConfigurationError(
            f"{name} must be a finite probability between 0 and 1.",
            code="INVALID_EVOLUTION_PROBABILITY",
            context={name: value},
        )
    return float(value)


def _resolve_parameters(
    *,
    mut_frac: float,
    insert_frac: float,
    delete_frac: float,
    insert_min: int,
    insert_max: int,
    shift_min: int,
    shift_max: int,
    invert_min: int,
    invert_max: int,
    rc_prob: float,
) -> _EvolutionParameters:
    insert_bounds = _validate_range("insert", insert_min, insert_max)
    if insert_bounds[0] == 0:
        raise ConfigurationError(
            "insert_min must be at least 1.",
            code="INVALID_EVOLUTION_PARAMETER",
            context={"insert_min": insert_min},
        )
    shift_bounds = _validate_range("shift", shift_min, shift_max)
    invert_bounds = _validate_range("invert", invert_min, invert_max)
    return _EvolutionParameters(
        mut_frac=_validate_probability("mut_frac", mut_frac),
        insert_frac=_validate_probability("insert_frac", insert_frac),
        delete_frac=_validate_probability("delete_frac", delete_frac),
        insert_min=insert_bounds[0],
        insert_max=insert_bounds[1],
        shift_min=shift_bounds[0],
        shift_max=shift_bounds[1],
        invert_min=invert_bounds[0],
        invert_max=invert_bounds[1],
        rc_prob=_validate_probability("rc_prob", rc_prob),
    )


def _normalize_augmentations(
    augmentations: Sequence[EvolutionAugmentation] | str,
) -> tuple[EvolutionAugmentation, ...]:
    raw: tuple[str, ...]
    if isinstance(augmentations, str):
        raw = (augmentations,)
    else:
        try:
            raw = tuple(augmentations)
        except TypeError as exc:
            raise ConfigurationError(
                "augmentations must be an iterable of supported operation names.",
                code="INVALID_EVOLUTION_AUGMENTATIONS",
            ) from exc
    resolved: list[EvolutionAugmentation] = []
    for name in raw:
        if not isinstance(name, str) or name not in _EVOLUTION_PRIORITY_INDEX:
            raise ConfigurationError(
                "Unsupported EvoAug operation.",
                code="INVALID_EVOLUTION_AUGMENTATION",
                context={"augmentation": name},
                hint="Use mutation, deletion, insertion, translocation, inversion, "
                "or reverse_complement.",
            )
        resolved.append(name)
    if len(set(resolved)) != len(resolved):
        raise ConfigurationError(
            "augmentations cannot contain duplicate operation names.",
            code="DUPLICATE_EVOLUTION_AUGMENTATION",
        )
    return tuple(resolved)


def _resolve_random_source(
    seed: int | None,
    rng: random.Random | None,
) -> tuple[random.Random, int | None, RandomSource, RandomState]:
    if (seed is None) == (rng is None):
        raise ConfigurationError(
            "Evolution generation requires exactly one of seed or rng.",
            code="EVOLUTION_RANDOM_SOURCE_REQUIRED",
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
    generator = random.Random(seed) if rng is None else rng
    return (
        generator,
        seed,
        "seed" if rng is None else "rng",
        cast(RandomState, generator.getstate()),
    )


def _validate_input(sequence: DNASequence) -> tuple[DNASequence, str]:
    source = require_sequence(sequence)
    if source.topology is not Topology.LINEAR:
        raise ConfigurationError(
            "Evolution generation requires a linear DNASequence.",
            code="EVOLUTION_CIRCULAR_NOT_SUPPORTED",
        )
    if source.is_gapped:
        raise UnsupportedGapOperationError(
            "Evolution generation does not operate across explicit Gaps.",
            code="EVOLUTION_GAPPED_SEQUENCE_NOT_SUPPORTED",
            hint="Resolve or split Gap parts before applying EvoAug operations.",
        )
    symbols = source.symbols
    if not symbols:
        raise ConfigurationError(
            "Evolution generation requires at least one nucleotide.",
            code="EVOLUTION_EMPTY_SEQUENCE",
        )
    if any(symbol not in _CANONICAL_BASES for symbol in symbols):
        raise ConfigurationError(
            "Evolution generation requires canonical A, C, G, and T symbols.",
            code="EVOLUTION_CANONICAL_DNA_REQUIRED",
            hint="Resolve IUPAC ambiguity before using the one-hot EvoAug operations.",
        )
    return source, symbols


def _random_dna(generator: random.Random, length: int) -> str:
    return "".join(generator.choice(_CANONICAL_BASES) for _ in range(length))


def _reverse_complement_text(text: str, source: DNASequence) -> str:
    return DNASequence(text, alphabet=source.alphabet).reverse_complement().symbols


def _require_segment_maximum(name: str, maximum: int, sequence_length: int) -> None:
    if maximum > sequence_length:
        raise ConfigurationError(
            f"{name}_max cannot exceed the current sequence length.",
            code="EVOLUTION_SEGMENT_TOO_LONG",
            context={
                "operation": name,
                f"{name}_max": maximum,
                "sequence_length": sequence_length,
            },
            hint=f"Set {name}_max to at most {sequence_length} for this sequence.",
        )


def _apply_mutation(
    symbols: str,
    generator: random.Random,
    mut_frac: float,
) -> tuple[str, EvolutionStep]:
    positions = tuple(position for position in range(len(symbols)) if generator.random() < mut_frac)
    chars = list(symbols)
    for position in positions:
        chars[position] = generator.choice(_MUTATION_ALTERNATIVES[chars[position]])
    mutation_count = len(positions)
    return (
        "".join(chars),
        EvolutionStep(
            "mutation",
            mutation_count > 0,
            length=mutation_count,
            mutation_attempts=mutation_count,
        ),
    )


def _apply_deletion(
    symbols: str,
    generator: random.Random,
    delete_frac: float,
) -> tuple[str, EvolutionStep]:
    retained = tuple(symbol for symbol in symbols if generator.random() >= delete_frac)
    deletion_count = len(symbols) - len(retained)
    return (
        "".join(retained),
        EvolutionStep("deletion", deletion_count > 0, length=deletion_count),
    )


def _apply_insertion(
    symbols: str,
    generator: random.Random,
    parameters: _EvolutionParameters,
) -> tuple[str, EvolutionStep]:
    parts: list[str] = []
    inserted_length = 0
    for symbol in symbols:
        parts.append(symbol)
        if generator.random() < parameters.insert_frac:
            length = generator.randint(parameters.insert_min, parameters.insert_max)
            parts.append(_random_dna(generator, length))
            inserted_length += length
    return (
        "".join(parts),
        EvolutionStep("insertion", inserted_length > 0, length=inserted_length),
    )


def _apply_translocation(
    symbols: str,
    generator: random.Random,
    parameters: _EvolutionParameters,
) -> tuple[str, EvolutionStep]:
    if not symbols:
        return symbols, EvolutionStep("translocation", False, shift=0)
    magnitude = generator.randint(parameters.shift_min, parameters.shift_max)
    shift = -magnitude if generator.random() < 0.5 else magnitude
    normalized = shift % len(symbols)
    generated = symbols if normalized == 0 else symbols[-normalized:] + symbols[:-normalized]
    return generated, EvolutionStep("translocation", normalized != 0, shift=shift)


def _apply_inversion(
    symbols: str,
    source: DNASequence,
    generator: random.Random,
    parameters: _EvolutionParameters,
) -> tuple[str, EvolutionStep]:
    _require_segment_maximum("invert", parameters.invert_max, len(symbols))
    length = generator.randint(parameters.invert_min, parameters.invert_max)
    start = generator.randint(0, len(symbols) - parameters.invert_max)
    segment = _reverse_complement_text(symbols[start : start + length], source)
    generated = symbols[:start] + segment + symbols[start + length :]
    return (
        generated,
        EvolutionStep("inversion", length > 0, start=start, end=start + length, length=length),
    )


def _apply_reverse_complement(
    symbols: str,
    source: DNASequence,
    generator: random.Random,
    probability: float,
) -> tuple[str, EvolutionStep]:
    applied = generator.random() < probability
    generated = _reverse_complement_text(symbols, source) if applied else symbols
    return generated, EvolutionStep("reverse_complement", applied)


def _apply_augmentation(
    symbols: str,
    source: DNASequence,
    generator: random.Random,
    augmentation: EvolutionAugmentation,
    parameters: _EvolutionParameters,
) -> tuple[str, EvolutionStep]:
    if augmentation == "mutation":
        return _apply_mutation(symbols, generator, parameters.mut_frac)
    if augmentation == "deletion":
        return _apply_deletion(symbols, generator, parameters.delete_frac)
    if augmentation == "insertion":
        return _apply_insertion(symbols, generator, parameters)
    if augmentation == "translocation":
        return _apply_translocation(symbols, generator, parameters)
    if augmentation == "inversion":
        return _apply_inversion(symbols, source, generator, parameters)
    return _apply_reverse_complement(symbols, source, generator, parameters.rc_prob)


def evolution_generate(
    sequence: DNASequence,
    *,
    augmentations: Sequence[EvolutionAugmentation] = DEFAULT_EVOLUTION_AUGMENTATIONS,
    max_augmentations: int = 1,
    hard_aug: bool = True,
    seed: int | None = None,
    rng: random.Random | None = None,
    mut_frac: float = 0.05,
    insert_frac: float = 0.05,
    delete_frac: float = 0.05,
    insert_min: int = 1,
    insert_max: int = 1,
    shift_min: int = 0,
    shift_max: int = 20,
    invert_min: int = 0,
    invert_max: int = 20,
    rc_prob: float = 0.5,
) -> EvolutionGenerationResult:
    """Generate one EvoAug-inspired DNA sequence variant.

    The candidate operations are sampled without replacement and applied in
    EvoAug priority order: inversion, deletion, translocation, insertion,
    reverse-complement, then mutation.  ``hard_aug=True`` applies exactly
    ``max_augmentations`` selected operations; ``False`` samples a count from
    one through that maximum.  Exactly one of ``seed`` and ``rng`` is required.
    ``mut_frac``, ``insert_frac``, and ``delete_frac`` are independent
    per-nucleotide probabilities for their respective selected operations.
    A mutation always changes the selected base.  An insertion is placed after
    the selected base and has a random inclusive length from ``insert_min`` to
    ``insert_max``.  A deletion removes the selected base.

    The official EvoAug package operates on fixed-shape one-hot tensors and
    includes Gaussian noise.  This sequence-level API keeps only operations
    that produce valid DNA symbols.  Use ``indel_generate()`` when exactly one
    random contiguous insertion or deletion is required.
    """

    source, symbols = _validate_input(sequence)
    names = _normalize_augmentations(augmentations)
    if isinstance(max_augmentations, bool) or not isinstance(max_augmentations, int):
        raise ConfigurationError(
            "max_augmentations must be a non-negative integer.",
            code="INVALID_MAX_AUGMENTATIONS",
            context={"max_augmentations": max_augmentations},
        )
    if max_augmentations < 0 or max_augmentations > len(names):
        raise ConfigurationError(
            "max_augmentations must be between zero and the number of augmentations.",
            code="INVALID_MAX_AUGMENTATIONS",
            context={"max_augmentations": max_augmentations, "augmentation_count": len(names)},
        )
    if not isinstance(hard_aug, bool):
        raise ConfigurationError(
            "hard_aug must be boolean.",
            code="INVALID_HARD_AUGMENTATION_MODE",
            context={"hard_aug": hard_aug},
        )
    parameters = _resolve_parameters(
        mut_frac=mut_frac,
        insert_frac=insert_frac,
        delete_frac=delete_frac,
        insert_min=insert_min,
        insert_max=insert_max,
        shift_min=shift_min,
        shift_max=shift_max,
        invert_min=invert_min,
        invert_max=invert_max,
        rc_prob=rc_prob,
    )
    generator, resolved_seed, random_source, state_before = _resolve_random_source(seed, rng)
    if max_augmentations == 0:
        selected: tuple[EvolutionAugmentation, ...] = ()
    else:
        selected_count = max_augmentations if hard_aug else generator.randint(1, max_augmentations)
        selected = tuple(
            sorted(
                generator.sample(names, selected_count),
                key=lambda name: _EVOLUTION_PRIORITY_INDEX[name],
            )
        )

    generated = symbols
    steps: list[EvolutionStep] = []
    for augmentation in selected:
        generated, step = _apply_augmentation(
            generated,
            source,
            generator,
            augmentation,
            parameters,
        )
        steps.append(step)
    result = DNASequence(
        generated,
        alphabet=source.alphabet,
        topology=Topology.LINEAR,
        strandedness=source.strandedness,
    )
    return EvolutionGenerationResult(
        sequence=result,
        steps=tuple(steps),
        seed=resolved_seed,
        random_source=random_source,
        rng_algorithm_version=random.Random.VERSION,
        rng_state_before=state_before,
        rng_state_after=cast(RandomState, generator.getstate()),
    )


def _validate_positive_int(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code="INVALID_SEQUENCE_GENERATION_PARAMETER",
            context={name: value},
        )
    return value


def _build_linear_sequence(source: DNASequence, symbols: str) -> DNASequence:
    return DNASequence(
        symbols,
        alphabet=source.alphabet,
        topology=Topology.LINEAR,
        strandedness=source.strandedness,
    )


def _apply_contiguous_deletion(
    symbols: str,
    generator: random.Random,
    minimum: int,
    maximum: int,
    *,
    pad_indels: bool,
) -> tuple[str, EvolutionStep]:
    _require_segment_maximum("delete", maximum, len(symbols))
    length = generator.randint(minimum, maximum)
    start = generator.randint(0, len(symbols) - maximum)
    generated = symbols[:start] + symbols[start + length :]
    if pad_indels:
        padding = _random_dna(generator, maximum)
        pad_begin = length // 2
        pad_end = length - pad_begin
        generated = padding[:pad_begin] + generated + padding[maximum - pad_end :]
    return (
        generated,
        EvolutionStep("deletion", length > 0, start=start, end=start + length, length=length),
    )


def _apply_contiguous_insertion(
    symbols: str,
    generator: random.Random,
    minimum: int,
    maximum: int,
    *,
    pad_indels: bool,
) -> tuple[str, EvolutionStep]:
    length = generator.randint(minimum, maximum)
    start = generator.randrange(len(symbols) + 1)
    if pad_indels:
        padding = _random_dna(generator, maximum)
        pad_begin = (maximum - length) // 2
        generated = (
            padding[:pad_begin]
            + symbols[:start]
            + padding[pad_begin : pad_begin + length]
            + symbols[start:]
            + padding[pad_begin + length :]
        )
    else:
        generated = symbols[:start] + _random_dna(generator, length) + symbols[start:]
    return (
        generated,
        EvolutionStep("insertion", length > 0, start=start, end=start, length=length),
    )


def indel_generate(
    sequence: DNASequence,
    *,
    operation: IndelOperation,
    min_length: int = 1,
    max_length: int = 20,
    seed: int | None = None,
    rng: random.Random | None = None,
    pad_indels: bool = False,
) -> EvolutionGenerationResult:
    """Generate one explicit random insertion or deletion.

    ``pad_indels=False`` returns a naturally variable-length sequence.  Set it
    to ``True`` to use the fixed-shape random-DNA padding used by EvoAug.
    """

    if not isinstance(operation, str) or operation not in {"insertion", "deletion"}:
        raise ConfigurationError(
            "operation must be 'insertion' or 'deletion'.",
            code="INVALID_INDEL_OPERATION",
            context={"operation": operation},
        )
    if not isinstance(pad_indels, bool):
        raise ConfigurationError(
            "pad_indels must be boolean.",
            code="INVALID_EVOLUTION_PARAMETER",
            context={"pad_indels": pad_indels},
        )
    source, symbols = _validate_input(sequence)
    minimum, maximum = _validate_range("indel", min_length, max_length)
    generator, resolved_seed, random_source, state_before = _resolve_random_source(seed, rng)
    if operation == "deletion":
        generated, step = _apply_contiguous_deletion(
            symbols,
            generator,
            minimum,
            maximum,
            pad_indels=pad_indels,
        )
    else:
        generated, step = _apply_contiguous_insertion(
            symbols,
            generator,
            minimum,
            maximum,
            pad_indels=pad_indels,
        )
    return EvolutionGenerationResult(
        sequence=_build_linear_sequence(source, generated),
        steps=(step,),
        seed=resolved_seed,
        random_source=random_source,
        rng_algorithm_version=random.Random.VERSION,
        rng_state_before=state_before,
        rng_state_after=cast(RandomState, generator.getstate()),
    )


def _random_segment_boundaries(
    sequence_length: int,
    segment_count: int,
    generator: random.Random,
) -> tuple[int, ...]:
    if segment_count == 1:
        return (0, sequence_length)
    internal = sorted(generator.sample(range(1, sequence_length), segment_count - 1))
    return (0, *internal, sequence_length)


def rearrange_generate(
    sequence: DNASequence,
    *,
    operation: RearrangementOperation = "exchange",
    segment_count: int = 3,
    seed: int | None = None,
    rng: random.Random | None = None,
) -> RearrangementResult:
    """Rearrange random contiguous segments by exchange, inversion, or duplication.

    ``exchange`` randomly permutes the segments.  ``inversion`` replaces one
    segment with its reverse complement, and ``duplication`` inserts a copy of
    one segment immediately after itself.
    """

    source, symbols = _validate_input(sequence)
    if not isinstance(operation, str) or operation not in {
        "exchange",
        "inversion",
        "duplication",
    }:
        raise ConfigurationError(
            "operation must be 'exchange', 'inversion', or 'duplication'.",
            code="INVALID_REARRANGEMENT_OPERATION",
            context={"operation": operation},
        )
    resolved_segment_count = _validate_positive_int("segment_count", segment_count)
    if resolved_segment_count > len(symbols):
        raise ConfigurationError(
            "segment_count cannot exceed the sequence length.",
            code="INVALID_REARRANGEMENT_SEGMENT_COUNT",
            context={"segment_count": segment_count, "sequence_length": len(symbols)},
        )
    generator, resolved_seed, random_source, state_before = _resolve_random_source(seed, rng)
    boundaries = _random_segment_boundaries(len(symbols), resolved_segment_count, generator)
    segments = [symbols[start:end] for start, end in pairwise(boundaries)]
    permutation: tuple[int, ...] | None = None
    selected_segment: int | None = None

    if operation == "exchange":
        order = list(range(resolved_segment_count))
        generator.shuffle(order)
        if resolved_segment_count > 1 and order == list(range(resolved_segment_count)):
            order[0], order[1] = order[1], order[0]
        permutation = tuple(order)
        generated = "".join(segments[index] for index in order)
    else:
        selected_index = generator.randrange(resolved_segment_count)
        selected_segment = selected_index
        if operation == "inversion":
            selected = _reverse_complement_text(segments[selected_index], source)
            generated_segments = segments.copy()
            generated_segments[selected_index] = selected
            generated = "".join(generated_segments)
        else:
            generated = "".join(
                (
                    *segments[: selected_index + 1],
                    segments[selected_index],
                    *segments[selected_index + 1 :],
                )
            )

    return RearrangementResult(
        sequence=_build_linear_sequence(source, generated),
        operation=operation,
        breakpoints=boundaries[1:-1],
        permutation=permutation,
        selected_segment=selected_segment,
        seed=resolved_seed,
        random_source=random_source,
        rng_algorithm_version=random.Random.VERSION,
        rng_state_before=state_before,
        rng_state_after=cast(RandomState, generator.getstate()),
    )


def _kmer_counts(symbols: str, k: int) -> tuple[tuple[str, int], ...]:
    return tuple(
        sorted(Counter(symbols[index : index + k] for index in range(len(symbols) - k + 1)).items())
    )


def _random_eulerian_shuffle(
    symbols: str,
    k: int,
    generator: random.Random,
) -> str:
    adjacency: dict[str, list[str]] = {}
    for index in range(len(symbols) - k + 1):
        kmer = symbols[index : index + k]
        adjacency.setdefault(kmer[:-1], []).append(kmer[-1])
    for edges in adjacency.values():
        generator.shuffle(edges)

    stack_nodes = [symbols[: k - 1]]
    stack_edges: list[str] = []
    trail_edges: list[str] = []
    while stack_nodes:
        node = stack_nodes[-1]
        outgoing = adjacency.get(node)
        if outgoing:
            next_base = outgoing.pop()
            stack_nodes.append(node[1:] + next_base)
            stack_edges.append(next_base)
        else:
            stack_nodes.pop()
            if stack_edges:
                trail_edges.append(stack_edges.pop())

    expected_edge_count = len(symbols) - k + 1
    if len(trail_edges) != expected_edge_count:
        raise ConfigurationError(
            "The source sequence did not yield a complete Eulerian path.",
            code="KMER_SHUFFLE_GRAPH_FAILURE",
            context={"k": k, "edge_count": expected_edge_count},
        )
    return symbols[: k - 1] + "".join(reversed(trail_edges))


def kmer_shuffle(
    sequence: DNASequence,
    *,
    k: int = 2,
    seed: int | None = None,
    rng: random.Random | None = None,
    ensure_different: bool = True,
    max_attempts: int = 100,
) -> KmerShuffleResult:
    """Shuffle a sequence while exactly preserving overlapping k-mer counts.

    For ``k >= 2`` this samples randomized Eulerian paths in the sequence's
    de Bruijn multigraph.  ``ensure_different=True`` retries until the output
    differs from the source, then raises if no alternative was found.
    """

    source, symbols = _validate_input(sequence)
    resolved_k = _validate_positive_int("k", k)
    if resolved_k > len(symbols):
        raise ConfigurationError(
            "k cannot exceed the sequence length.",
            code="INVALID_KMER_LENGTH",
            context={"k": k, "sequence_length": len(symbols)},
        )
    if not isinstance(ensure_different, bool):
        raise ConfigurationError(
            "ensure_different must be boolean.",
            code="INVALID_KMER_SHUFFLE_PARAMETER",
            context={"ensure_different": ensure_different},
        )
    resolved_attempts = _validate_positive_int("max_attempts", max_attempts)
    generator, resolved_seed, random_source, state_before = _resolve_random_source(seed, rng)
    expected_counts = _kmer_counts(symbols, resolved_k)
    generated = symbols
    attempts = 0

    for attempt in range(1, resolved_attempts + 1):
        if resolved_k == 1:
            chars = list(symbols)
            generator.shuffle(chars)
            candidate = "".join(chars)
        else:
            candidate = _random_eulerian_shuffle(symbols, resolved_k, generator)
        if _kmer_counts(candidate, resolved_k) != expected_counts:
            raise ConfigurationError(
                "The generated sequence changed the requested k-mer counts.",
                code="KMER_SHUFFLE_COUNT_FAILURE",
                context={"k": resolved_k},
            )
        generated = candidate
        attempts = attempt
        if not ensure_different or candidate != symbols:
            break
    else:
        raise ConfigurationError(
            "No different sequence was found with the same k-mer counts.",
            code="KMER_SHUFFLE_NO_ALTERNATIVE",
            context={"k": resolved_k, "max_attempts": resolved_attempts},
            hint="Set ensure_different=False when the source has a unique reconstruction.",
        )

    return KmerShuffleResult(
        sequence=_build_linear_sequence(source, generated),
        k=resolved_k,
        kmer_counts=expected_counts,
        attempts=attempts,
        seed=resolved_seed,
        random_source=random_source,
        rng_algorithm_version=random.Random.VERSION,
        rng_state_before=state_before,
        rng_state_after=cast(RandomState, generator.getstate()),
    )


def crossover(
    first: DNASequence,
    second: DNASequence,
    *,
    position: int | None = None,
    seed: int | None = None,
    rng: random.Random | None = None,
) -> CrossoverResult:
    """Create a child by one-point crossover of two equal-length DNA sequences.

    ``position`` is a zero-based boundary and must leave at least one base from
    each parent.  If it is omitted, exactly one of ``seed`` or ``rng`` chooses
    the boundary reproducibly.
    """

    first_source, first_symbols = _validate_input(first)
    second_source, second_symbols = _validate_input(second)
    if first_source.strandedness is not second_source.strandedness:
        raise ConfigurationError(
            "Crossover parents must have matching strandedness.",
            code="CROSSOVER_STRANDEDNESS_MISMATCH",
        )
    if len(first_symbols) != len(second_symbols):
        raise ConfigurationError(
            "Crossover parents must have equal sequence lengths.",
            code="CROSSOVER_LENGTH_MISMATCH",
            context={
                "first_length": len(first_symbols),
                "second_length": len(second_symbols),
            },
        )
    if len(first_symbols) < 2:
        raise ConfigurationError(
            "Crossover requires parents with at least two bases.",
            code="CROSSOVER_SEQUENCE_TOO_SHORT",
        )

    if position is None:
        generator, resolved_seed, random_source, state_before = _resolve_random_source(seed, rng)
        resolved_position = generator.randint(1, len(first_symbols) - 1)
        state_after: RandomState | None = cast(RandomState, generator.getstate())
        rng_version: int | None = random.Random.VERSION
    else:
        if seed is not None or rng is not None:
            raise ConfigurationError(
                "seed and rng are only valid when position is omitted.",
                code="UNEXPECTED_RANDOM_SOURCE",
            )
        if isinstance(position, bool) or not isinstance(position, int):
            raise ConfigurationError(
                "position must be an integer crossover boundary.",
                code="INVALID_CROSSOVER_POSITION",
                context={"position": position},
            )
        if not 1 <= position < len(first_symbols):
            raise ConfigurationError(
                "position must leave at least one base from each parent.",
                code="INVALID_CROSSOVER_POSITION",
                context={"position": position, "sequence_length": len(first_symbols)},
            )
        generator = None
        resolved_seed = None
        resolved_position = position
        random_source = "none"
        state_before = None
        state_after = None
        rng_version = None

    child = first_symbols[:resolved_position] + second_symbols[resolved_position:]
    return CrossoverResult(
        sequence=_build_linear_sequence(first_source, child),
        position=resolved_position,
        first_parent_length=len(first_symbols),
        second_parent_length=len(second_symbols),
        seed=resolved_seed,
        random_source=random_source,
        rng_algorithm_version=rng_version,
        rng_state_before=state_before,
        rng_state_after=state_after,
    )


__all__ = [
    "DEFAULT_EVOLUTION_AUGMENTATIONS",
    "CrossoverResult",
    "EvolutionAugmentation",
    "EvolutionGenerationResult",
    "EvolutionStep",
    "IndelOperation",
    "KmerShuffleResult",
    "RearrangementOperation",
    "RearrangementResult",
    "crossover",
    "evolution_generate",
    "indel_generate",
    "kmer_shuffle",
    "rearrange_generate",
]
