"""Tests for the DNA adaptations of MOSES Frag and SNN metrics."""

from __future__ import annotations

import json
from importlib import import_module
from typing import Any, cast

import pytest

from dnakit.core import DNARecord, DNASequence, DNASet, Gap
from dnakit.evaluation import (
    EvaluationLimits,
    FragmentSimilarityConfig,
    SNNConfig,
    evaluate_fragment_similarity,
    evaluate_snn,
)
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.fingerprints import hashed_kmer_fingerprint
from dnakit.similarity import fingerprint_similarity


def _set(*pairs: tuple[str, str]) -> DNASet:
    return DNASet(DNARecord(DNASequence(symbols), record_id) for record_id, symbols in pairs)


def test_fragment_similarity_matches_known_count_cosine_and_is_symmetric() -> None:
    generated = _set(("a", "AAAA"), ("c", "CCCC"))
    reference = _set(("a2", "AAAA"), ("g", "GGGG"))
    config = FragmentSimilarityConfig(k=2, canonical=False, show_progress=False)

    forward = evaluate_fragment_similarity(generated, reference, config=config)
    reverse = evaluate_fragment_similarity(reference, generated, config=config)

    assert forward.metrics["frag"] == pytest.approx(0.5)
    assert reverse.metrics["frag"] == pytest.approx(0.5)
    assert forward.metrics["generated_fragment_observations"] == 6
    assert forward.metrics["reference_fragment_observations"] == 6
    assert forward.metrics["shared_unique_fragments"] == 1
    assert forward.parameters["higher_is_better"] is True


def test_fragment_similarity_identical_is_one_and_disjoint_is_zero() -> None:
    first = _set(("a", "AAAA"), ("c", "CCCC"))
    second = _set(("g", "GGGG"), ("t", "TTTT"))
    config = FragmentSimilarityConfig(k=2, canonical=False, show_progress=False)

    identical = evaluate_fragment_similarity(first, first, config=config)
    disjoint = evaluate_fragment_similarity(first, second, config=config)

    assert identical.metrics["fragment_similarity"] == pytest.approx(1.0)
    assert disjoint.metrics["fragment_similarity"] == pytest.approx(0.0)
    json.dumps(identical.to_dict())


def test_fragment_similarity_respects_gaps_ambiguity_and_limits() -> None:
    gapped = DNASet([DNARecord(DNASequence(["AA", Gap(3), "AA"]), "gapped")])
    ambiguous = DNASet([DNARecord(DNASequence("AANAAA", alphabet="iupac"), "ambiguous")])
    reference = _set(("reference", "AAAAAA"))
    config = FragmentSimilarityConfig(k=3, show_progress=False)

    report = evaluate_fragment_similarity(ambiguous, reference, config=config)
    assert report.metrics["generated_fragment_observations"] == 1
    assert report.metrics["generated_ignored_ambiguity_count"] == 1

    with pytest.raises(ConfigurationError) as empty:
        evaluate_fragment_similarity(gapped, reference, config=config)
    assert empty.value.code == "EMPTY_FRAGMENT_VECTOR"

    with pytest.raises(InvalidAlphabetError):
        evaluate_fragment_similarity(
            ambiguous,
            reference,
            config=FragmentSimilarityConfig(
                k=3,
                ambiguity_policy="error",
                show_progress=False,
            ),
        )

    with pytest.raises(ConfigurationError) as bounded:
        evaluate_fragment_similarity(
            reference,
            reference,
            config=FragmentSimilarityConfig(
                k=3,
                max_kmer_observations=2,
                show_progress=False,
            ),
        )
    assert bounded.value.code == "FRAGMENT_KMER_LIMIT"


def test_snn_matches_known_binary_tanimoto_and_stable_tie_break() -> None:
    generated = _set(("copy", "AAAA"), ("far", "ATAT"))
    reference = _set(("first", "AAAA"), ("second", "CCCC"))
    report = evaluate_snn(
        generated,
        reference,
        config=SNNConfig(
            k=2,
            n_bits=65_536,
            canonical=False,
            max_fingerprint_elements=300_000,
            show_progress=False,
        ),
    )

    assert report.metrics["snn"] == pytest.approx(0.5)
    assert report.metrics["minimum_nearest_similarity"] == pytest.approx(0.0)
    assert report.entries[0].metrics["nearest_reference_id"] == "first"
    assert report.entries[1].metrics["nearest_reference_id"] == "first"
    assert report.parameters["dna_fingerprint"] == "SHA-256 hashed binary k-mer fingerprint"
    json.dumps(report.to_dict())


def test_snn_canonical_reverse_complements_and_zero_vector_policy() -> None:
    config = SNNConfig(k=2, n_bits=1_024, show_progress=False)
    canonical = evaluate_snn(
        _set(("generated", "GGGG")),
        _set(("reference", "CCCC")),
        config=config,
    )
    empty = evaluate_snn(
        _set(("short", "A")),
        _set(("short-reference", "T")),
        config=config,
    )

    assert canonical.metrics["snn"] == pytest.approx(1.0)
    assert empty.metrics["snn"] == pytest.approx(1.0)
    assert empty.metrics["generated_zero_fingerprint_count"] == 1
    assert empty.parameters["zero_vector_policy"].startswith("two empty")


def test_snn_pair_values_match_public_binary_fingerprint_tanimoto() -> None:
    generated = _set(("generated", "AACCGGTT"))
    reference = _set(("first", "AACCGGTA"), ("second", "TTTTGGGG"))
    config = SNNConfig(k=3, n_bits=2_048, canonical=False, show_progress=False)

    report = evaluate_snn(generated, reference, config=config)
    query_fingerprint = hashed_kmer_fingerprint(
        generated[0],
        k=config.k,
        n_bits=config.n_bits,
        canonical=config.canonical,
        seed=config.seed,
        representation="sparse",
        ambiguity_policy=config.ambiguity_policy,
    )
    expected = max(
        fingerprint_similarity(
            query_fingerprint,
            hashed_kmer_fingerprint(
                record,
                k=config.k,
                n_bits=config.n_bits,
                canonical=config.canonical,
                seed=config.seed,
                representation="sparse",
                ambiguity_policy=config.ambiguity_policy,
            ),
            metric="tanimoto",
        ).value
        for record in reference
    )

    assert report.metrics["snn"] == pytest.approx(expected)
    entry = report.entries[0]
    assert entry.metrics["nearest_union_bit_count"] >= entry.metrics["nearest_shared_bit_count"]


def test_snn_limits_precede_fingerprint_work(monkeypatch: pytest.MonkeyPatch) -> None:
    module = import_module("dnakit.evaluation.generative")
    fingerprint_called = False

    def unexpected_fingerprint(*args: object, **kwargs: object) -> object:
        nonlocal fingerprint_called
        fingerprint_called = True
        raise AssertionError("fingerprints must not be built for rejected workloads")

    monkeypatch.setattr(module, "hashed_kmer_fingerprint", unexpected_fingerprint)
    records = _set(("a", "AAAA"), ("b", "CCCC"))
    with pytest.raises(ConfigurationError) as pair_limit:
        evaluate_snn(
            records,
            records,
            config=SNNConfig(
                show_progress=False,
                limits=EvaluationLimits(max_pairwise_comparisons=3),
            ),
        )
    assert pair_limit.value.code == "EVALUATION_PAIRWISE_LIMIT"
    assert fingerprint_called is False

    with pytest.raises(ConfigurationError) as element_limit:
        evaluate_snn(
            records,
            records,
            config=SNNConfig(
                n_bits=1_024,
                max_fingerprint_elements=4_095,
                show_progress=False,
            ),
        )
    assert element_limit.value.code == "SNN_FINGERPRINT_LIMIT"
    assert fingerprint_called is False


def test_generative_configs_reject_invalid_values() -> None:
    defaults = FragmentSimilarityConfig()
    snn_defaults = SNNConfig()
    assert (defaults.k, defaults.canonical, defaults.show_progress) == (3, True, True)
    assert (snn_defaults.k, snn_defaults.n_bits, snn_defaults.show_progress) == (7, 1_024, True)

    with pytest.raises(ConfigurationError):
        FragmentSimilarityConfig(k=0)
    with pytest.raises(ConfigurationError) as fragment_progress:
        FragmentSimilarityConfig(show_progress=cast(Any, 1))
    assert fragment_progress.value.code == "INVALID_FRAGMENT_SIMILARITY_CONFIG"
    with pytest.raises(ConfigurationError):
        SNNConfig(n_bits=0)
    with pytest.raises(ConfigurationError) as seed:
        SNNConfig(seed=-1)
    assert seed.value.code == "INVALID_SNN_CONFIG"
