from __future__ import annotations

import json

import pytest

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.molbio import (
    CodonUsageTable,
    RuleOptimizationConfig,
    generate_mutation_library,
    optimize_codons,
    optimize_protein_codons,
    optimize_sequence_rules,
)


def test_rule_optimizer_improves_gc_repeat_constraints_deterministically() -> None:
    config = RuleOptimizationConfig(
        target_gc_range=(0.5, 0.5),
        max_homopolymer=2,
        seed=7,
    )

    first = optimize_sequence_rules(DNASequence("AAAAAA"), config)
    second = optimize_sequence_rules(DNASequence("AAAAAA"), config)

    assert first.optimized == second.optimized
    assert first.final_score < first.initial_score
    assert first.constraints_satisfied
    assert first.parameters["seed"] == 7
    assert json.loads(json.dumps(first.to_dict()))["optimized"]["parts"]


def test_rule_optimizer_can_preserve_translation_and_reports_local_incompleteness() -> None:
    preserved = optimize_sequence_rules(
        DNASequence("GCTGCT"),
        RuleOptimizationConfig(
            forbidden_motifs=("GCT",),
            preserve_translation=True,
        ),
    )
    blocked = optimize_sequence_rules(
        DNASequence("AAAA"),
        RuleOptimizationConfig(
            target_gc_range=(1.0, 1.0),
            allowed_positions=(),
        ),
    )

    assert preserved.optimized.symbols != "GCTGCT"
    assert blocked.constraints_satisfied is False
    assert blocked.issues[-1].code == "RULE_OPTIMIZATION_INCOMPLETE"
    with pytest.raises(ConfigurationError, match="complete frame-0"):
        optimize_sequence_rules(
            DNASequence("GCTA"),
            RuleOptimizationConfig(preserve_translation=True),
        )
    with pytest.raises(ConfigurationError, match="max_total_scoring_cells"):
        optimize_sequence_rules(
            DNASequence("A" * 100),
            RuleOptimizationConfig(
                target_gc_range=(0.5, 0.5),
                max_total_scoring_cells=1,
            ),
        )


def test_codon_optimization_uses_versioned_table_and_preserves_protein() -> None:
    table = CodonUsageTable(
        "test-host",
        "2026.1",
        {"GCT": 1.0, "GCC": 10.0, "GCA": 2.0, "GCG": 3.0},
    )

    result = optimize_codons(DNASequence("GCTGCT"), table)

    assert result.optimized.symbols == "GCCGCC"
    assert result.original_translation == result.optimized_translation == "AA"
    assert result.optimized_cai > result.original_cai
    assert result.usage_table_checksum == table.checksum
    assert result.provenance.reference is not None
    assert result.provenance.reference.checksum == table.checksum


def test_codon_gc_dynamic_program_and_preserved_motif_constraints() -> None:
    table = CodonUsageTable(
        "test-host",
        "1",
        {"GCT": 1.0, "GCC": 10.0, "GCA": 2.0, "GCG": 3.0},
    )
    constrained = optimize_codons(
        DNASequence("GCTGCT"),
        table,
        gc_range=(4 / 6, 4 / 6),
    )
    locked = optimize_codons(
        DNASequence("GCTGCT"),
        table,
        preserve_motifs=("GCT",),
    )

    assert constrained.optimized.symbols == "GCAGCA"
    assert locked.optimized.symbols == "GCTGCT"
    with pytest.raises(ConfigurationError, match="No synonymous sequence"):
        optimize_codons(
            DNASequence("GCTGCT"),
            table,
            gc_range=(0.0, 0.0),
        )
    with pytest.raises(ConfigurationError, match="max_dp_cells"):
        optimize_codons(DNASequence("GCTGCT"), table, max_dp_cells=1)


def test_protein_codon_optimization_builds_an_audited_reverse_translation() -> None:
    table = CodonUsageTable(
        "test-host",
        "1",
        {"GCT": 1.0, "GCC": 10.0, "GCA": 2.0, "GCG": 3.0},
    )

    result = optimize_protein_codons("AA", table)

    assert result.optimized.symbols == "GCCGCC"
    assert result.optimized_translation == "AA"
    assert result.parameters["input_type"] == "protein"
    assert result.issues[0].code == "PROTEIN_REVERSE_TRANSLATION_BASELINE"
    with pytest.raises(ConfigurationError, match="unsupported"):
        optimize_protein_codons("AX", table)


def test_mutation_library_single_combinatorial_and_seeded_sample() -> None:
    single = generate_mutation_library(
        DNASequence("AC"),
        {0: ("C", "G"), 1: ("A",)},
    )
    combinatorial = generate_mutation_library(
        DNASequence("AC"),
        {0: ("C", "G"), 1: ("A",)},
        mode="combinatorial",
        max_order=2,
    )
    sampled_one = generate_mutation_library(
        DNASequence("ACGT"),
        range(4),
        mode="combinatorial",
        max_order=2,
        sample_size=5,
        seed=11,
    )
    sampled_two = generate_mutation_library(
        DNASequence("ACGT"),
        range(4),
        mode="combinatorial",
        max_order=2,
        sample_size=5,
        seed=11,
    )

    assert len(single.variants) == 3
    assert combinatorial.total_possible_variants == 5
    assert len(combinatorial.variants) == 5
    assert sampled_one.variants == sampled_two.variants
    assert all(1 <= len(variant.mutations) <= 2 for variant in sampled_one.variants)


def test_mutation_library_limits_are_checked_before_large_materialization() -> None:
    with pytest.raises(ConfigurationError, match="max_enumerated_variants"):
        generate_mutation_library(
            DNASequence("A" * 10),
            range(10),
            mode="combinatorial",
            max_order=5,
            max_enumerated_variants=10,
        )
    with pytest.raises(ConfigurationError, match="max_variants"):
        generate_mutation_library(
            DNASequence("AAAA"),
            range(4),
            max_variants=2,
        )
    with pytest.raises(ConfigurationError, match="outside"):
        generate_mutation_library(DNASequence("AC"), [2])
    with pytest.raises(ConfigurationError, match="max_output_bases"):
        generate_mutation_library(DNASequence("AAAA"), [0], max_output_bases=4)
    with pytest.raises(ConfigurationError, match="alternate"):
        generate_mutation_library(DNASequence("AC"), {0: ("a",)})
