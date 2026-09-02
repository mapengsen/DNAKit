from __future__ import annotations

import itertools
import json

import pytest

from dnakit.core import DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError, UnsupportedGapOperationError
from dnakit.molbio import (
    AssemblyFragment,
    check_end_compatibility,
    classify_restriction_end,
    digest_restriction,
    ligate_fragments,
    simulate_assembly,
)


def test_restriction_digest_classifies_ends_and_religation_reconstructs_template() -> None:
    sequence = DNASequence("AGAATTCT")

    result = digest_restriction(sequence, ["EcoRI"])

    assert [(cut.top_cut, cut.bottom_cut, cut.polarity) for cut in result.cuts] == [
        (2, 6, "5prime")
    ]
    assert [fragment.sequence.symbols for fragment in result.fragments] == ["AG", "AATTCT"]
    assert result.fragments[0].right_end.overhang_sequence_5to3 == "AATT"
    compatibility = check_end_compatibility(
        result.fragments[0].right_end,
        result.fragments[1].left_end,
    )
    assert compatibility.compatible
    ligated = ligate_fragments(result.fragments)
    assert ligated.product == sequence
    assert json.loads(json.dumps(result.to_dict()))["algorithm_version"] == (
        "dnakit-restriction-digest-v1"
    )


def test_blunt_end_requires_explicit_policy_and_phosphorylation() -> None:
    result = digest_restriction(DNASequence("AGGCCA"), ["HaeIII"])
    left, right = result.fragments

    assert left.right_end.polarity == "blunt"
    assert not check_end_compatibility(left.right_end, right.left_end).compatible
    assert check_end_compatibility(
        left.right_end,
        right.left_end,
        allow_blunt=True,
    ).compatible
    explicit = classify_restriction_end(DNASequence("ACGT"), 2, 2)
    assert explicit.end.polarity == "blunt"


def test_circular_single_cut_linearizes_one_full_length_fragment() -> None:
    result = digest_restriction(
        DNASequence("AATTCCCG", topology=Topology.CIRCULAR),
        ["EcoRI"],
    )

    assert len(result.fragments) == 1
    assert result.fragments[0].sequence.symbols == "AATTCCCG"
    assert result.fragments[0].sequence.topology is Topology.LINEAR
    assert result.fragments[0].wraps_origin
    assert result.fragments[0].left_end.cohesive_key == "AATT"


def test_digest_rejects_unmodeled_states_and_resource_truncation() -> None:
    with pytest.raises(InvalidAlphabetError):
        digest_restriction(DNASequence("NGAATTC", alphabet="iupac"), ["EcoRI"])
    ambiguous = digest_restriction(
        DNASequence("NGAATTC", alphabet="iupac"),
        ["EcoRI"],
        allow_ambiguous_template=True,
    )
    assert ambiguous.issues[0].code == "AMBIGUOUS_DIGEST_MATCHES"
    with pytest.raises(UnsupportedGapOperationError):
        digest_restriction(DNASequence(["AAA", Gap(1), "GAATTC"]), ["EcoRI"])
    with pytest.raises(ConfigurationError, match="unmethylated"):
        digest_restriction(DNASequence("GAATTC"), ["EcoRI"], methylation_state="dam")
    with pytest.raises(ConfigurationError):
        digest_restriction(DNASequence("GAATTCGAATTC"), ["EcoRI"], max_fragments=1)
    with pytest.raises(ConfigurationError, match="item limit"):
        digest_restriction(
            DNASequence("GAATTC"),
            itertools.repeat("EcoRI"),
            max_enzymes=2,
        )


def test_gibson_and_predigested_golden_gate_sequence_abstractions() -> None:
    gibson = simulate_assembly(
        [
            AssemblyFragment("left", DNASequence("AAAACCCC")),
            AssemblyFragment("right", DNASequence("CCCCGGGG")),
        ],
        method="gibson",
        min_overlap=4,
        max_overlap=4,
    )
    digest = digest_restriction(DNASequence("AGAATTCT"), ["EcoRI"])
    golden_gate = simulate_assembly(digest.fragments, method="golden_gate")

    assert gibson.product.symbols == "AAAACCCCGGGG"
    assert gibson.steps[0].junction_sequence == "CCCC"
    assert golden_gate.product.symbols == "AGAATTCT"
    assert golden_gate.parameters["kinetics_modeled"] is False


def test_assembly_rejects_missing_overlap_and_raw_golden_gate_fragments() -> None:
    with pytest.raises(ConfigurationError, match="No exact overlap"):
        simulate_assembly(
            [DNASequence("AAAA"), DNASequence("CCCC")],
            method="lcr",
            min_overlap=2,
            max_overlap=2,
        )
    with pytest.raises(ConfigurationError, match="pre-digested"):
        simulate_assembly(
            [DNASequence("AAAA"), DNASequence("AAAA")],
            method="golden_gate",
        )
    with pytest.raises(ConfigurationError, match="item limit"):
        simulate_assembly(
            [DNASequence("AAAA"), DNASequence("AAAA")],
            method="gibson",
            overlaps=itertools.repeat("AA"),
            min_overlap=2,
            max_overlap=2,
        )
