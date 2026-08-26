"""Bounded rule, codon, and mutation-library design algorithms."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace

from dnakit.core import DNASequence, Issue, IssueSeverity
from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.exceptions import ConfigurationError
from dnakit.ops import translate

from ._shared import (
    finite_fraction,
    freeze_parameters,
    iupac_compatible,
    materialize_bounded,
    native_provenance,
    require_sequence,
    reverse_complement_text,
    validate_iupac_text,
    validate_positive_int,
    validate_text,
)
from .results import (
    CodonOptimizationResult,
    Mutation,
    MutationLibraryResult,
    RuleOptimizationResult,
    SequenceChange,
    SequenceVariant,
)

_BASES = "ACGT"
_CODONS = tuple("".join(parts) for parts in itertools.product(_BASES, repeat=3))
_CODON_TO_AA = {codon: translate(codon) for codon in _CODONS}
_AA_TO_CODONS: dict[str, tuple[str, ...]] = {
    amino_acid: tuple(codon for codon in _CODONS if _CODON_TO_AA[codon] == amino_acid)
    for amino_acid in sorted(set(_CODON_TO_AA.values()))
}


@dataclass(frozen=True, init=False)
class RuleOptimizationConfig:
    """Deterministic local-search constraints and hard resource limits."""

    target_gc_range: tuple[float, float] | None
    forbidden_motifs: tuple[str, ...]
    forbidden_motifs_both_strands: bool
    max_homopolymer: int | None
    allowed_positions: tuple[int, ...] | None
    preserve_translation: bool
    max_iterations: int
    max_candidate_evaluations: int
    max_total_scoring_cells: int
    seed: int
    gc_weight: float
    motif_weight: float
    homopolymer_weight: float

    def __init__(
        self,
        *,
        target_gc_range: tuple[float, float] | None = None,
        forbidden_motifs: Iterable[str] = (),
        forbidden_motifs_both_strands: bool = True,
        max_homopolymer: int | None = None,
        allowed_positions: Iterable[int] | None = None,
        preserve_translation: bool = False,
        max_iterations: int = 1_000,
        max_candidate_evaluations: int = 100_000,
        max_total_scoring_cells: int = 500_000_000,
        seed: int = 0,
        gc_weight: float = 1.0,
        motif_weight: float = 100.0,
        homopolymer_weight: float = 10.0,
    ) -> None:
        resolved_gc: tuple[float, float] | None = None
        if target_gc_range is not None:
            if not isinstance(target_gc_range, tuple) or len(target_gc_range) != 2:
                raise ConfigurationError("target_gc_range must contain two bounds.")
            lower = finite_fraction(target_gc_range[0], "minimum GC")
            upper = finite_fraction(target_gc_range[1], "maximum GC")
            if lower > upper:
                raise ConfigurationError("Minimum GC cannot exceed maximum GC.")
            resolved_gc = (lower, upper)
        if isinstance(forbidden_motifs, (str, bytes)):
            raise ConfigurationError("forbidden_motifs must be an iterable of motifs.")
        motif_inputs = materialize_bounded(
            forbidden_motifs,
            max_items=10_000,
            name="forbidden motifs",
        )
        motifs = tuple(validate_iupac_text(item, "forbidden motif") for item in motif_inputs)
        if not isinstance(forbidden_motifs_both_strands, bool) or not isinstance(
            preserve_translation, bool
        ):
            raise ConfigurationError("Rule optimization switches must be booleans.")
        if max_homopolymer is not None:
            validate_positive_int(max_homopolymer, "max_homopolymer", maximum=1_000_000)
        positions: tuple[int, ...] | None = None
        if allowed_positions is not None:
            if isinstance(allowed_positions, (str, bytes)):
                raise ConfigurationError("allowed_positions must be integer positions.")
            copied = materialize_bounded(
                allowed_positions,
                max_items=1_000_000,
                name="allowed positions",
            )
            if any(
                isinstance(item, bool) or not isinstance(item, int) or item < 0 for item in copied
            ):
                raise ConfigurationError("allowed_positions must be non-negative integers.")
            if len(copied) != len(set(copied)):
                raise ConfigurationError("allowed_positions cannot contain duplicates.")
            positions = tuple(sorted(copied))
        validate_positive_int(max_iterations, "max_iterations", maximum=1_000_000)
        validate_positive_int(
            max_candidate_evaluations,
            "max_candidate_evaluations",
            maximum=1_000_000_000,
        )
        validate_positive_int(
            max_total_scoring_cells,
            "max_total_scoring_cells",
            maximum=10_000_000_000,
        )
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise ConfigurationError("seed must be an integer.")
        weights = []
        for value, name in (
            (gc_weight, "gc_weight"),
            (motif_weight, "motif_weight"),
            (homopolymer_weight, "homopolymer_weight"),
        ):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value <= 0
            ):
                raise ConfigurationError(f"{name} must be finite and positive.")
            weights.append(float(value))
        object.__setattr__(self, "target_gc_range", resolved_gc)
        object.__setattr__(self, "forbidden_motifs", motifs)
        object.__setattr__(self, "forbidden_motifs_both_strands", forbidden_motifs_both_strands)
        object.__setattr__(self, "max_homopolymer", max_homopolymer)
        object.__setattr__(self, "allowed_positions", positions)
        object.__setattr__(self, "preserve_translation", preserve_translation)
        object.__setattr__(self, "max_iterations", max_iterations)
        object.__setattr__(self, "max_candidate_evaluations", max_candidate_evaluations)
        object.__setattr__(self, "max_total_scoring_cells", max_total_scoring_cells)
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "gc_weight", weights[0])
        object.__setattr__(self, "motif_weight", weights[1])
        object.__setattr__(self, "homopolymer_weight", weights[2])


def _motif_occurrences(symbols: str, motifs: tuple[str, ...], both_strands: bool) -> int:
    patterns: list[str] = list(motifs)
    if both_strands:
        patterns.extend(reverse_complement_text(motif) for motif in motifs)
    unique = tuple(sorted(set(patterns)))
    count = 0
    for motif in unique:
        for start in range(len(symbols) - len(motif) + 1):
            if all(
                iupac_compatible(symbols[start + offset], pattern)
                for offset, pattern in enumerate(motif)
            ):
                count += 1
    return count


def _homopolymer_excess(symbols: str, maximum: int | None) -> int:
    if maximum is None or not symbols:
        return 0
    excess = 0
    run = 1
    for previous, current in itertools.pairwise(symbols):
        if previous == current:
            run += 1
        else:
            excess += max(0, run - maximum)
            run = 1
    return excess + max(0, run - maximum)


def _rule_components(symbols: str, config: RuleOptimizationConfig) -> tuple[float, int, int]:
    gc_violation = 0.0
    if config.target_gc_range is not None:
        gc = sum(base in "GC" for base in symbols) / len(symbols)
        lower, upper = config.target_gc_range
        gc_violation = max(lower - gc, gc - upper, 0.0) * len(symbols)
    return (
        gc_violation,
        _motif_occurrences(
            symbols,
            config.forbidden_motifs,
            config.forbidden_motifs_both_strands,
        ),
        _homopolymer_excess(symbols, config.max_homopolymer),
    )


def _rule_score(components: tuple[float, int, int], config: RuleOptimizationConfig) -> float:
    return (
        config.gc_weight * components[0]
        + config.motif_weight * components[1]
        + config.homopolymer_weight * components[2]
    )


def optimize_sequence_rules(
    sequence: DNASequence,
    config: RuleOptimizationConfig,
    *,
    max_sequence_length: int = 1_000_000,
) -> RuleOptimizationResult:
    """Greedily lower a documented rule score with deterministic seeded ties."""

    symbols = require_sequence(
        sequence,
        operation="rule-based sequence optimization",
        max_length=max_sequence_length,
        canonical=True,
        allow_circular=False,
    )
    if not isinstance(config, RuleOptimizationConfig):
        raise ConfigurationError("config must be RuleOptimizationConfig.")
    allowed = (
        list(range(len(symbols)))
        if config.allowed_positions is None
        else list(config.allowed_positions)
    )
    if any(position >= len(symbols) for position in allowed):
        raise ConfigurationError("allowed_positions contains a position outside the sequence.")
    original_translation: str | None = None
    if config.preserve_translation:
        if len(symbols) % 3:
            raise ConfigurationError(
                "preserve_translation requires a complete frame-0 coding sequence."
            )
        original_translation = translate(symbols, incomplete_policy="error")
    order = list(allowed)
    random.Random(config.seed).shuffle(order)
    current = symbols
    initial_components = _rule_components(current, config)
    current_score = _rule_score(initial_components, config)
    selected_patterns = set(config.forbidden_motifs)
    if config.forbidden_motifs_both_strands:
        selected_patterns.update(
            reverse_complement_text(motif) for motif in tuple(selected_patterns)
        )
    cells_per_score = (
        (len(symbols) if config.target_gc_range is not None else 0)
        + (len(symbols) if config.max_homopolymer is not None else 0)
        + sum(max(0, len(symbols) - len(motif) + 1) * len(motif) for motif in selected_patterns)
        + (len(symbols) if config.preserve_translation else 0)
    )
    evaluation_bound = min(
        config.max_candidate_evaluations,
        config.max_iterations * len(order) * 3,
    )
    estimated_scoring_cells = cells_per_score * (evaluation_bound + 1)
    if estimated_scoring_cells > config.max_total_scoring_cells:
        raise ConfigurationError(
            "Rule optimization exceeds max_total_scoring_cells.",
            code="RULE_OPTIMIZATION_CELL_LIMIT_EXCEEDED",
            context={
                "estimated_scoring_cells": estimated_scoring_cells,
                "max_total_scoring_cells": config.max_total_scoring_cells,
            },
            hint="Reduce sequence/constraint size, iterations, or candidate evaluations.",
        )
    changes: list[SequenceChange] = []
    evaluations = 0
    limit_reached = False
    for _ in range(config.max_iterations):
        best: tuple[float, int, str, str] | None = None
        before_components = _rule_components(current, config)
        for position in order:
            for alternate in _BASES:
                if alternate == current[position]:
                    continue
                if evaluations >= config.max_candidate_evaluations:
                    limit_reached = True
                    break
                evaluations += 1
                candidate = current[:position] + alternate + current[position + 1 :]
                if (
                    original_translation is not None
                    and translate(candidate, incomplete_policy="error") != original_translation
                ):
                    continue
                score = _rule_score(_rule_components(candidate, config), config)
                proposal = (score, position, alternate, candidate)
                if score < current_score - 1e-12 and (best is None or score < best[0] - 1e-12):
                    best = proposal
            if limit_reached:
                break
        if best is None:
            break
        score, position, alternate, candidate = best
        after_components = _rule_components(candidate, config)
        labels = ("gc_range", "forbidden_motif", "homopolymer")
        reason = next(
            (
                label
                for label, before, after in zip(
                    labels, before_components, after_components, strict=True
                )
                if after < before
            ),
            "weighted_rule_score",
        )
        changes.append(SequenceChange(position, current[position], alternate, reason))
        current = candidate
        current_score = score
        if all(value <= 1e-12 for value in after_components):
            break
        if limit_reached:
            break
    final_components = _rule_components(current, config)
    satisfied = all(value <= 1e-12 for value in final_components)
    issues: list[Issue] = []
    if limit_reached:
        issues.append(
            Issue(
                "RULE_OPTIMIZATION_EVALUATION_LIMIT",
                IssueSeverity.WARNING,
                "Optimization stopped at max_candidate_evaluations.",
            )
        )
    if not satisfied:
        issues.append(
            Issue(
                "RULE_OPTIMIZATION_INCOMPLETE",
                IssueSeverity.WARNING,
                "The bounded greedy search ended with unsatisfied constraints.",
                details={
                    "gc_violation_bases": final_components[0],
                    "forbidden_motif_occurrences": final_components[1],
                    "homopolymer_excess_bases": final_components[2],
                },
            )
        )
    return RuleOptimizationResult(
        original=sequence,
        optimized=DNASequence(current, strandedness=sequence.strandedness),
        changes=tuple(changes),
        initial_score=_rule_score(initial_components, config),
        final_score=current_score,
        iterations=len(changes),
        constraints_satisfied=satisfied,
        method="seeded_deterministic_greedy_rule_search",
        algorithm_version="dnakit-rule-optimization-v1",
        parameters=freeze_parameters(
            {
                "target_gc_range": config.target_gc_range,
                "forbidden_motifs": config.forbidden_motifs,
                "forbidden_motifs_both_strands": config.forbidden_motifs_both_strands,
                "max_homopolymer": config.max_homopolymer,
                "allowed_positions": config.allowed_positions,
                "preserve_translation": config.preserve_translation,
                "max_iterations": config.max_iterations,
                "max_candidate_evaluations": config.max_candidate_evaluations,
                "max_total_scoring_cells": config.max_total_scoring_cells,
                "candidate_evaluations": evaluations,
                "estimated_scoring_cells": estimated_scoring_cells,
                "seed": config.seed,
                "weights": {
                    "gc": config.gc_weight,
                    "motif": config.motif_weight,
                    "homopolymer": config.homopolymer_weight,
                },
                "search_guarantees_global_optimum": False,
            }
        ),
        provenance=native_provenance(),
        issues=tuple(issues),
    )


@dataclass(frozen=True, init=False)
class CodonUsageTable:
    """User-supplied, versioned table of non-negative table-1 codon weights."""

    name: str
    version: str
    frequencies: FrozenDict
    checksum: str

    def __init__(self, name: str, version: str, frequencies: Mapping[str, float]) -> None:
        validate_text(name, "codon usage table name")
        validate_text(version, "codon usage table version")
        if not isinstance(frequencies, Mapping) or not frequencies:
            raise ConfigurationError("frequencies must be a non-empty codon mapping.")
        normalized: dict[str, float] = {}
        for raw_codon, raw_value in frequencies.items():
            codon = validate_iupac_text(raw_codon, "codon")
            if len(codon) != 3 or set(codon) - set(_BASES):
                raise ConfigurationError("Codon usage keys must be canonical DNA triplets.")
            if (
                isinstance(raw_value, bool)
                or not isinstance(raw_value, (int, float))
                or not math.isfinite(raw_value)
                or raw_value < 0
            ):
                raise ConfigurationError("Codon frequencies must be finite and non-negative.")
            normalized[codon] = float(raw_value)
        payload = json.dumps(
            {"name": name, "version": version, "frequencies": normalized},
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "frequencies", freeze_mapping(normalized))
        object.__setattr__(self, "checksum", hashlib.sha256(payload).hexdigest())


def _literal_locked_positions(symbols: str, motifs: tuple[str, ...]) -> frozenset[int]:
    locked: set[int] = set()
    for motif in motifs:
        start = symbols.find(motif)
        while start >= 0:
            locked.update(range(start, start + len(motif)))
            start = symbols.find(motif, start + 1)
    return frozenset(locked)


def _cai(codons: tuple[str, ...], table: CodonUsageTable) -> float:
    logs: list[float] = []
    for codon in codons:
        amino_acid = _CODON_TO_AA[codon]
        synonymous = _AA_TO_CODONS[amino_acid]
        maximum = max(_codon_frequency(table, item) for item in synonymous)
        weight = _codon_frequency(table, codon)
        if maximum <= 0.0 or weight <= 0.0:
            return 0.0
        logs.append(math.log(weight / maximum))
    return math.exp(math.fsum(logs) / len(logs)) if logs else 0.0


def _codon_frequency(table: CodonUsageTable, codon: str) -> float:
    value = table.frequencies.get(codon, 0.0)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AssertionError("Validated codon frequency changed type.")
    return float(value)


def optimize_codons(
    sequence: DNASequence,
    usage_table: CodonUsageTable,
    *,
    gc_range: tuple[float, float] | None = None,
    forbidden_codons: Iterable[str] = (),
    preserve_motifs: Iterable[str] = (),
    max_codons: int = 2_000,
    max_dp_cells: int = 2_000_000,
) -> CodonOptimizationResult:
    """Maximize relative codon weights under an optional global GC interval."""

    validate_positive_int(max_codons, "max_codons", maximum=100_000)
    validate_positive_int(max_dp_cells, "max_dp_cells", maximum=100_000_000)
    symbols = require_sequence(
        sequence,
        operation="codon optimization",
        max_length=max_codons * 3,
        canonical=True,
        allow_circular=False,
    )
    if len(symbols) % 3:
        raise ConfigurationError("Codon optimization requires a complete frame-0 CDS.")
    if not isinstance(usage_table, CodonUsageTable):
        raise ConfigurationError("usage_table must be CodonUsageTable.")
    if isinstance(forbidden_codons, (str, bytes)):
        raise ConfigurationError("forbidden_codons must be an iterable of codons.")
    forbidden_inputs = materialize_bounded(
        forbidden_codons,
        max_items=10_000,
        name="forbidden codons",
    )
    forbidden = frozenset(validate_iupac_text(item, "forbidden codon") for item in forbidden_inputs)
    if any(len(item) != 3 or set(item) - set(_BASES) for item in forbidden):
        raise ConfigurationError("Forbidden codons must be canonical triplets.")
    if isinstance(preserve_motifs, (str, bytes)):
        raise ConfigurationError("preserve_motifs must be an iterable of motifs.")
    motif_inputs = materialize_bounded(
        preserve_motifs,
        max_items=10_000,
        name="preserved motifs",
    )
    motifs = tuple(validate_iupac_text(item, "preserved motif") for item in motif_inputs)
    if any(set(item) - set(_BASES) for item in motifs):
        raise ConfigurationError("Preserved motifs must be exact canonical DNA strings.")
    locked = _literal_locked_positions(symbols, motifs)
    lower_count, upper_count = 0, len(symbols)
    resolved_gc: tuple[float, float] | None = None
    if gc_range is not None:
        if not isinstance(gc_range, tuple) or len(gc_range) != 2:
            raise ConfigurationError("gc_range must contain two bounds.")
        lower = finite_fraction(gc_range[0], "minimum GC")
        upper = finite_fraction(gc_range[1], "maximum GC")
        if lower > upper:
            raise ConfigurationError("Minimum GC cannot exceed maximum GC.")
        lower_count = math.ceil(lower * len(symbols) - 1e-12)
        upper_count = math.floor(upper * len(symbols) + 1e-12)
        if lower_count > upper_count:
            raise ConfigurationError("GC range contains no attainable integer GC count.")
        resolved_gc = (lower, upper)
    original_codons = tuple(symbols[index : index + 3] for index in range(0, len(symbols), 3))
    original_translation = translate(symbols, incomplete_policy="error")

    # value = (log relative-weight score, previous GC count, chosen codon)
    layers: list[dict[int, tuple[float, int | None, str]]] = [{0: (0.0, None, "")}]
    total_cells = 1
    for codon_index, original in enumerate(original_codons):
        amino_acid = _CODON_TO_AA[original]
        codon_positions = range(codon_index * 3, codon_index * 3 + 3)
        is_locked = any(position in locked for position in codon_positions)
        candidates = (original,) if is_locked else _AA_TO_CODONS[amino_acid]
        candidates = tuple(
            codon
            for codon in candidates
            if codon not in forbidden and _codon_frequency(usage_table, codon) > 0.0
        )
        if not candidates:
            raise ConfigurationError(
                "No positive-weight allowed synonymous codon is available.",
                code="CODON_OPTIMIZATION_NO_OPTION",
                context={"codon_index": codon_index, "amino_acid": amino_acid},
            )
        maximum = max(_codon_frequency(usage_table, item) for item in _AA_TO_CODONS[amino_acid])
        current = layers[-1]
        next_layer: dict[int, tuple[float, int | None, str]] = {}
        for previous_gc in sorted(current):
            previous_score = current[previous_gc][0]
            for candidate in sorted(candidates):
                gc_count = previous_gc + sum(base in "GC" for base in candidate)
                score = previous_score + math.log(
                    _codon_frequency(usage_table, candidate) / maximum
                )
                existing = next_layer.get(gc_count)
                if existing is None or score > existing[0] + 1e-15:
                    next_layer[gc_count] = (score, previous_gc, candidate)
        total_cells += len(next_layer)
        if total_cells > max_dp_cells:
            raise ConfigurationError(
                "Codon optimization exceeds max_dp_cells.",
                code="CODON_DP_LIMIT_EXCEEDED",
                context={"dp_cells": total_cells, "max_dp_cells": max_dp_cells},
            )
        layers.append(next_layer)
    feasible = [gc_count for gc_count in layers[-1] if lower_count <= gc_count <= upper_count]
    if not feasible:
        raise ConfigurationError(
            "No synonymous sequence satisfies the requested GC interval.",
            code="CODON_GC_CONSTRAINT_INFEASIBLE",
        )
    final_gc = max(feasible, key=lambda count: (layers[-1][count][0], -count))
    optimized_reversed: list[str] = []
    cursor = final_gc
    for layer_index in range(len(original_codons), 0, -1):
        _, previous, chosen = layers[layer_index][cursor]
        optimized_reversed.append(chosen)
        if previous is None:
            raise AssertionError("Codon DP backpointer is missing.")
        cursor = previous
    optimized_codons = tuple(reversed(optimized_reversed))
    optimized_symbols = "".join(optimized_codons)
    optimized_translation = translate(optimized_symbols, incomplete_policy="error")
    if optimized_translation != original_translation:
        raise AssertionError("Synonymous DP changed the translated protein.")
    changes = tuple(
        SequenceChange(index * 3, before, after, "synonymous_codon_weight")
        for index, (before, after) in enumerate(zip(original_codons, optimized_codons, strict=True))
        if before != after
    )
    return CodonOptimizationResult(
        original=sequence,
        optimized=DNASequence(optimized_symbols, strandedness=sequence.strandedness),
        original_translation=original_translation,
        optimized_translation=optimized_translation,
        changes=changes,
        original_cai=_cai(original_codons, usage_table),
        optimized_cai=_cai(optimized_codons, usage_table),
        original_gc_fraction=sum(base in "GC" for base in symbols) / len(symbols),
        optimized_gc_fraction=sum(base in "GC" for base in optimized_symbols) / len(symbols),
        usage_table_name=usage_table.name,
        usage_table_version=usage_table.version,
        usage_table_checksum=usage_table.checksum,
        method="synonymous_codon_dynamic_programming",
        algorithm_version="dnakit-codon-optimization-v1",
        parameters=freeze_parameters(
            {
                "genetic_code": 1,
                "gc_range": resolved_gc,
                "forbidden_codons": sorted(forbidden),
                "preserve_motifs": motifs,
                "locked_position_count": len(locked),
                "objective": "maximum_sum_log_relative_codon_weight",
                "tie_break": "lowest_gc_then_first_sorted_dp_traversal_path",
                "dp_cells": total_cells,
                "max_dp_cells": max_dp_cells,
                "max_codons": max_codons,
            }
        ),
        provenance=native_provenance(
            reference_name=usage_table.name,
            reference_version=usage_table.version,
            reference_checksum=usage_table.checksum,
        ),
        issues=(),
    )


def optimize_protein_codons(
    protein: str,
    usage_table: CodonUsageTable,
    *,
    gc_range: tuple[float, float] | None = None,
    forbidden_codons: Iterable[str] = (),
    max_amino_acids: int = 2_000,
    max_dp_cells: int = 2_000_000,
) -> CodonOptimizationResult:
    """Reverse-translate a protein and optimize it with the same bounded DP."""

    if not isinstance(protein, str) or not protein:
        raise ConfigurationError("protein must be a non-empty uppercase amino-acid string.")
    validate_positive_int(max_amino_acids, "max_amino_acids", maximum=100_000)
    if len(protein) > max_amino_acids:
        raise ConfigurationError("Protein exceeds max_amino_acids.")
    invalid = sorted(set(protein) - set(_AA_TO_CODONS))
    if invalid:
        raise ConfigurationError(
            "Protein contains residues unsupported by genetic code table 1.",
            context={"invalid_residues": invalid},
        )
    if not isinstance(usage_table, CodonUsageTable):
        raise ConfigurationError("usage_table must be CodonUsageTable.")
    forbidden_inputs = materialize_bounded(
        forbidden_codons,
        max_items=10_000,
        name="forbidden codons",
    )
    forbidden = frozenset(validate_iupac_text(item, "forbidden codon") for item in forbidden_inputs)
    if any(len(item) != 3 or set(item) - set(_BASES) for item in forbidden):
        raise ConfigurationError("Forbidden codons must be canonical triplets.")
    baseline: list[str] = []
    for residue_index, amino_acid in enumerate(protein):
        candidates = tuple(
            codon
            for codon in _AA_TO_CODONS[amino_acid]
            if codon not in forbidden and _codon_frequency(usage_table, codon) > 0.0
        )
        if not candidates:
            raise ConfigurationError(
                "No positive-weight allowed codon can encode a protein residue.",
                code="PROTEIN_CODON_OPTIMIZATION_NO_OPTION",
                context={"residue_index": residue_index, "amino_acid": amino_acid},
            )
        baseline.append(sorted(candidates)[0])
    result = optimize_codons(
        DNASequence("".join(baseline)),
        usage_table,
        gc_range=gc_range,
        forbidden_codons=forbidden,
        max_codons=max_amino_acids,
        max_dp_cells=max_dp_cells,
    )
    parameters = dict(result.parameters)
    parameters.update(
        {
            "input_type": "protein",
            "source_protein": protein,
            "baseline_reverse_translation": result.original.symbols,
        }
    )
    return replace(
        result,
        method="protein_reverse_translation_and_synonymous_codon_dynamic_programming",
        algorithm_version="dnakit-protein-codon-optimization-v1",
        parameters=freeze_parameters(parameters),
        issues=(
            Issue(
                "PROTEIN_REVERSE_TRANSLATION_BASELINE",
                IssueSeverity.INFO,
                "The result.original field is a deterministic synthetic DNA baseline.",
            ),
        ),
    )


def _replacement_map(
    symbols: str,
    positions: Mapping[int, Iterable[str]] | Iterable[int],
    allowed_bases: tuple[str, ...],
    max_positions: int,
) -> dict[int, tuple[str, ...]]:
    items: Iterable[tuple[int, Iterable[str]]]
    if isinstance(positions, Mapping):
        items = positions.items()
    else:
        if isinstance(positions, (str, bytes)):
            raise ConfigurationError("positions must be integer positions or a mapping.")
        items = ((position, allowed_bases) for position in positions)
    result: dict[int, tuple[str, ...]] = {}
    for position, replacements in items:
        if len(result) >= max_positions:
            raise ConfigurationError("Mutation target count exceeds max_positions.")
        if (
            isinstance(position, bool)
            or not isinstance(position, int)
            or not 0 <= position < len(symbols)
        ):
            raise ConfigurationError("Mutation position lies outside the source sequence.")
        copied = materialize_bounded(
            replacements,
            max_items=16,
            name="replacement bases",
            reject_text=False,
        )
        normalized = {validate_iupac_text(base, "replacement base") for base in copied}
        normalized.discard(symbols[position])
        resolved = tuple(sorted(normalized))
        if any(len(base) != 1 or base not in allowed_bases for base in resolved):
            raise ConfigurationError("Replacement bases must be members of allowed_bases.")
        if not resolved:
            raise ConfigurationError("Every mutation position needs at least one alternate base.")
        if position in result:
            raise ConfigurationError("Mutation positions cannot be repeated.")
        result[position] = resolved
    if not result:
        raise ConfigurationError("At least one mutation position is required.")
    return dict(sorted(result.items()))


def _variant_count(replacements: Mapping[int, tuple[str, ...]], max_order: int) -> int:
    counts = [len(values) for values in replacements.values()]
    totals = [0] * (max_order + 1)
    totals[0] = 1
    for count in counts:
        for order in range(max_order, 0, -1):
            totals[order] += totals[order - 1] * count
    return sum(totals[1:])


def generate_mutation_library(
    sequence: DNASequence,
    positions: Mapping[int, Iterable[str]] | Iterable[int],
    *,
    mode: str = "single",
    allowed_bases: Iterable[str] = _BASES,
    max_order: int = 1,
    sample_size: int | None = None,
    seed: int = 0,
    max_positions: int = 1_000,
    max_variants: int = 100_000,
    max_enumerated_variants: int = 1_000_000,
    max_output_bases: int = 100_000_000,
    max_sequence_length: int = 1_000_000,
) -> MutationLibraryResult:
    """Generate bounded single, saturation, or combinatorial substitutions."""

    symbols = require_sequence(
        sequence,
        operation="mutation library generation",
        max_length=max_sequence_length,
        canonical=True,
        allow_circular=False,
    )
    if mode not in ("single", "saturation", "combinatorial"):
        raise ConfigurationError("mode must be single, saturation, or combinatorial.")
    base_values = materialize_bounded(
        allowed_bases,
        max_items=16,
        name="allowed bases",
        reject_text=False,
    )
    allowed = tuple(sorted({validate_iupac_text(base, "allowed base") for base in base_values}))
    if any(len(base) != 1 or base not in _BASES for base in allowed):
        raise ConfigurationError("allowed_bases must contain canonical single bases.")
    validate_positive_int(max_positions, "max_positions", maximum=100_000)
    validate_positive_int(max_variants, "max_variants", maximum=10_000_000)
    validate_positive_int(max_enumerated_variants, "max_enumerated_variants", maximum=100_000_000)
    validate_positive_int(max_output_bases, "max_output_bases", maximum=10_000_000_000)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ConfigurationError("seed must be an integer.")
    replacements = _replacement_map(symbols, positions, allowed, max_positions)
    resolved_order = 1 if mode in ("single", "saturation") else max_order
    validate_positive_int(resolved_order, "max_order", maximum=20)
    resolved_order = min(resolved_order, len(replacements))
    total = _variant_count(replacements, resolved_order)
    if total > max_enumerated_variants:
        raise ConfigurationError(
            "Mutation search space exceeds max_enumerated_variants.",
            code="MUTATION_ENUMERATION_LIMIT_EXCEEDED",
            context={
                "total_possible_variants": total,
                "max_enumerated_variants": max_enumerated_variants,
            },
        )
    if sample_size is not None:
        validate_positive_int(sample_size, "sample_size", maximum=max_variants)
        if sample_size > total:
            raise ConfigurationError("sample_size cannot exceed total possible variants.")
        output_count = sample_size
    else:
        output_count = total
        if output_count > max_variants:
            raise ConfigurationError(
                "Mutation library exceeds max_variants; provide a bounded sample_size.",
                code="MUTATION_OUTPUT_LIMIT_EXCEEDED",
            )
    if output_count * len(symbols) > max_output_bases:
        raise ConfigurationError(
            "Mutation library exceeds max_output_bases.",
            code="MUTATION_OUTPUT_BASE_LIMIT_EXCEEDED",
            context={
                "output_bases": output_count * len(symbols),
                "max_output_bases": max_output_bases,
            },
        )
    selected_indices = (
        None if sample_size is None else set(random.Random(seed).sample(range(total), sample_size))
    )
    generated: list[SequenceVariant] = []
    positions_tuple = tuple(replacements)
    enumeration_index = 0
    for order in range(1, resolved_order + 1):
        for selected_positions in itertools.combinations(positions_tuple, order):
            options = tuple(replacements[position] for position in selected_positions)
            for alternates in itertools.product(*options):
                if selected_indices is not None and enumeration_index not in selected_indices:
                    enumeration_index += 1
                    continue
                mutation_tuple = tuple(
                    Mutation(position, symbols[position], alternate)
                    for position, alternate in zip(selected_positions, alternates, strict=True)
                )
                mutated = list(symbols)
                for mutation in mutation_tuple:
                    mutated[mutation.position] = mutation.alternate
                generated.append(
                    SequenceVariant(
                        id=f"variant_{enumeration_index + 1}",
                        sequence=DNASequence("".join(mutated), strandedness=sequence.strandedness),
                        mutations=mutation_tuple,
                    )
                )
                enumeration_index += 1
    if len(generated) != output_count:
        raise AssertionError("Mutation enumeration did not produce the expected bounded output.")
    return MutationLibraryResult(
        source=sequence,
        variants=tuple(generated),
        total_possible_variants=total,
        sampled=sample_size is not None,
        seed=seed,
        method="bounded_substitution_product_enumeration",
        algorithm_version="dnakit-mutation-library-v1",
        parameters=freeze_parameters(
            {
                "mode": mode,
                "positions": {str(key): value for key, value in replacements.items()},
                "allowed_bases": allowed,
                "max_order": resolved_order,
                "sample_size": sample_size,
                "seed": seed,
                "total_possible_variants": total,
                "max_positions": max_positions,
                "max_variants": max_variants,
                "max_enumerated_variants": max_enumerated_variants,
                "max_output_bases": max_output_bases,
            }
        ),
        provenance=native_provenance(),
        issues=(),
    )


__all__ = [
    "CodonUsageTable",
    "RuleOptimizationConfig",
    "generate_mutation_library",
    "optimize_codons",
    "optimize_protein_codons",
    "optimize_sequence_rules",
]
