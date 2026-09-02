"""Tests for reproducible random, stable-hash, stratified, group, and similarity splits."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.datasets import SplitConfig, split
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError


def _records(count: int) -> list[DNARecord]:
    return [DNARecord(DNASequence("ACGT"), f"r{index}") for index in range(count)]


def test_random_split_has_exact_counts_seed_and_original_subset_order() -> None:
    records = _records(10)
    config = SplitConfig(ratios={"train": 0.6, "valid": 0.2, "test": 0.2}, seed=17)

    first = split(records, config=config)
    second = split(records, config=config)

    assert first.assignments == second.assignments
    assert first.counts == {"train": 6, "valid": 2, "test": 2}
    for subset in first.subsets:
        numeric_ids = [int(record.id[1:]) for record in subset.records]
        assert numeric_ids == sorted(numeric_ids)
    assert json.loads(json.dumps(first.to_dict()))["seed"] == 17


def test_hash_split_is_independent_of_input_order_and_can_preserve_hash_order() -> None:
    records = _records(10)
    config = SplitConfig(
        method="hash",
        ratios={"train": 0.6, "valid": 0.2, "test": 0.2},
        seed=17,
        preserve_order=False,
    )

    first = split(records, config=config)
    second = split(list(reversed(records)), config=config)

    first_by_id = {item.record_id: item.split for item in first.assignments}
    second_by_id = {item.record_id: item.split for item in second.assignments}
    assert first_by_id == second_by_id
    assert first.counts == {"train": 6, "valid": 2, "test": 2}
    split_names = ("train", "valid", "test")
    assert {name: first.get(name).ids for name in split_names} == {
        name: second.get(name).ids for name in split_names
    }
    assert first.assignment_strategy == "sha256_record_id_rank_v1_largest_remainder"


def test_hash_split_rejects_duplicate_record_ids() -> None:
    records = [
        DNARecord(DNASequence("AC"), "duplicate"),
        DNARecord(DNASequence("GT"), "duplicate"),
    ]

    with pytest.raises(ConfigurationError) as exc_info:
        split(
            records,
            config=SplitConfig(method="hash", ratios={"train": 0.5, "test": 0.5}),
        )

    assert exc_info.value.code == "HASH_SPLIT_DUPLICATE_RECORD_ID"


def test_stratified_split_allocates_each_label_independently() -> None:
    records = [
        DNARecord(DNASequence("AC"), f"p{i}", metadata={"label": "positive"}) for i in range(4)
    ] + [DNARecord(DNASequence("GT"), f"n{i}", metadata={"label": "negative"}) for i in range(4)]
    result = split(
        records,
        config=SplitConfig(
            method="stratified",
            ratios={"train": 0.5, "test": 0.5},
            metadata_key="label",
            seed=2,
        ),
    )

    for name in ("train", "test"):
        labels = [record.metadata["label"] for record in result.get(name)]
        assert labels.count("positive") == 2
        assert labels.count("negative") == 2


def test_stratified_split_uses_global_quotas_for_singleton_strata() -> None:
    records = [
        DNARecord(DNASequence("AC"), f"r{i}", metadata={"label": f"label-{i}"}) for i in range(10)
    ]

    result = split(
        records,
        config=SplitConfig(
            method="stratified",
            ratios={"train": 0.8, "test": 0.2},
            metadata_key="label",
            seed=7,
        ),
    )

    assert result.counts == {"train": 8, "test": 2}
    assert result.assignment_strategy == "global_quota_round_robin_stratified"


def test_stratified_split_spreads_minority_assignments_across_small_strata() -> None:
    records = [
        DNARecord(DNASequence("AC"), f"r{label}-{member}", metadata={"label": label})
        for label in range(10)
        for member in range(2)
    ]

    result = split(
        records,
        config=SplitConfig(
            method="stratified",
            ratios={"train": 0.8, "test": 0.2},
            metadata_key="label",
            shuffle=False,
        ),
    )
    test_labels = {record.metadata["label"] for record in result.get("test")}

    assert result.counts == {"train": 16, "test": 4}
    assert len(test_labels) == 4


@pytest.mark.parametrize("metadata_key", ["cluster", "species", "chromosome", "donor", "locus"])
def test_group_split_never_separates_a_metadata_group(metadata_key: str) -> None:
    records = [
        DNARecord(DNASequence("AC"), f"a{i}", metadata={metadata_key: "A"}) for i in range(3)
    ] + [DNARecord(DNASequence("GT"), f"b{i}", metadata={metadata_key: "B"}) for i in range(3)]
    result = split(
        records,
        config=SplitConfig(
            method="group",
            ratios={"train": 0.5, "test": 0.5},
            metadata_key=metadata_key,
            seed=11,
        ),
    )
    split_by_id = {assignment.record_id: assignment.split for assignment in result.assignments}

    assert len({split_by_id[f"a{i}"] for i in range(3)}) == 1
    assert len({split_by_id[f"b{i}"] for i in range(3)}) == 1
    assert split_by_id["a0"] != split_by_id["b0"]


def test_missing_group_metadata_errors_or_becomes_a_separate_unit() -> None:
    records = _records(2)
    with pytest.raises(ConfigurationError) as exc_info:
        split(
            records,
            config=SplitConfig(
                method="group",
                ratios={"train": 0.5, "test": 0.5},
                metadata_key="donor",
            ),
        )
    assert exc_info.value.code == "SPLIT_METADATA_MISSING"

    result = split(
        records,
        config=SplitConfig(
            method="group",
            ratios={"train": 0.5, "test": 0.5},
            metadata_key="donor",
            missing_metadata_policy="separate",
        ),
    )
    assert tuple(result.counts.values()) == (1, 1)


def test_missing_group_key_cannot_collide_with_json_metadata() -> None:
    records = [
        DNARecord(DNASequence("AC"), "missing"),
        DNARecord(
            DNASequence("GT"),
            "json-lookalike",
            metadata={"donor": ["<missing>", 0]},
        ),
    ]

    result = split(
        records,
        config=SplitConfig(
            method="group",
            ratios={"train": 0.5, "test": 0.5},
            metadata_key="donor",
            missing_metadata_policy="separate",
            shuffle=False,
        ),
    )
    assignments = {item.record_id: item.split for item in result.assignments}

    assert assignments["missing"] != assignments["json-lookalike"]


def test_missing_stratum_key_cannot_collide_with_json_metadata() -> None:
    records = [
        DNARecord(DNASequence("AC"), "missing"),
        DNARecord(
            DNASequence("GT"),
            "json-lookalike",
            metadata={"label": ["<missing>", 0]},
        ),
        DNARecord(DNASequence("AA"), "other", metadata={"label": "other"}),
    ]

    result = split(
        records,
        config=SplitConfig(
            method="stratified",
            ratios={"train": 0.8, "test": 0.2},
            metadata_key="label",
            missing_metadata_policy="separate",
            shuffle=False,
        ),
    )
    assignments = {item.record_id: item.split for item in result.assignments}

    assert assignments == {
        "missing": "train",
        "json-lookalike": "train",
        "other": "test",
    }


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (True, 1),
        (True, 1.0),
        (1, 1.0),
        (
            {"outer": [True, {"value": 1.0}]},
            {"outer": [1, {"value": 1}]},
        ),
    ],
)
def test_group_split_metadata_keys_preserve_json_types_recursively(
    left: object,
    right: object,
) -> None:
    records = [
        DNARecord(DNASequence("AC"), "left", metadata={"group": left}),
        DNARecord(DNASequence("GT"), "right", metadata={"group": right}),
    ]

    result = split(
        records,
        config=SplitConfig(
            method="group",
            ratios={"train": 0.5, "test": 0.5},
            metadata_key="group",
            shuffle=False,
        ),
    )
    assignments = {item.record_id: item.split for item in result.assignments}

    assert assignments["left"] != assignments["right"]


@pytest.mark.parametrize(
    ("left", "right"),
    [
        (True, 1),
        (True, 1.0),
        (1, 1.0),
        (
            {"outer": [True, {"value": 1.0}]},
            {"outer": [1, {"value": 1}]},
        ),
    ],
)
def test_stratified_split_metadata_keys_preserve_json_types_recursively(
    left: object,
    right: object,
) -> None:
    records = [
        DNARecord(DNASequence("AC"), "left", metadata={"label": left}),
        DNARecord(DNASequence("GT"), "right-1", metadata={"label": right}),
        DNARecord(DNASequence("AA"), "right-2", metadata={"label": right}),
    ]

    result = split(
        records,
        config=SplitConfig(
            method="stratified",
            ratios={"train": 0.6, "test": 0.4},
            metadata_key="label",
            shuffle=False,
        ),
    )
    assignments = {item.record_id: item.split for item in result.assignments}

    assert assignments == {
        "left": "train",
        "right-1": "train",
        "right-2": "test",
    }


def test_similarity_split_uses_connected_components_to_prevent_cross_split_edges() -> None:
    records = [
        DNARecord(DNASequence("AAAA"), "a"),
        DNARecord(DNASequence("AAAT"), "b"),
        DNARecord(DNASequence("CCCC"), "c"),
        DNARecord(DNASequence("CCCG"), "d"),
    ]
    result = split(
        records,
        config=SplitConfig(
            method="similarity",
            ratios={"train": 0.5, "test": 0.5},
            similarity_k=2,
            similarity_threshold=0.5,
            seed=5,
        ),
    )
    assigned = {item.record_id: item.split for item in result.assignments}

    assert assigned["a"] == assigned["b"]
    assert assigned["c"] == assigned["d"]
    assert assigned["a"] != assigned["c"]
    assert result.component_count == 2
    assert result.pairwise_comparison_count == 6
    assert result.similarity_method == "kmer_jaccard"


def test_similarity_split_ambiguity_gap_size_and_empty_kmer_policies_are_explicit() -> None:
    ambiguous = DNARecord(DNASequence("AN", alphabet=DNAAlphabet.IUPAC), "ambiguous")
    gapped = DNARecord(DNASequence(["A", Gap(None), "T"]), "gapped")

    with pytest.raises(UnsupportedGapOperationError):
        split(
            [gapped],
            config=SplitConfig(method="similarity", ratios={"train": 0.5, "test": 0.5}),
        )
    split(
        [ambiguous],
        config=SplitConfig(
            method="similarity",
            ratios={"train": 0.5, "test": 0.5},
            similarity_ambiguity_policy="ignore",
        ),
    )
    with pytest.raises(ConfigurationError) as size_error:
        split(
            _records(3),
            config=SplitConfig(
                method="similarity",
                ratios={"train": 0.5, "test": 0.5},
                max_pairwise_records=2,
            ),
        )
    assert size_error.value.code == "SIMILARITY_SPLIT_SIZE_LIMIT"


def test_similarity_split_consumes_at_most_limit_plus_one_records() -> None:
    consumed: list[int] = []

    def records() -> object:
        for index in range(100):
            consumed.append(index)
            yield DNARecord(DNASequence("ACGT"), f"r{index}")

    with pytest.raises(ConfigurationError) as exc_info:
        split(
            records(),  # type: ignore[arg-type]
            config=SplitConfig(
                method="similarity",
                ratios={"train": 0.5, "test": 0.5},
                max_pairwise_records=2,
            ),
        )

    assert exc_info.value.code == "SIMILARITY_SPLIT_SIZE_LIMIT"
    assert consumed == [0, 1, 2]


def test_split_result_records_resource_and_assignment_parameters() -> None:
    result = split(
        _records(4),
        config=SplitConfig(
            method="similarity",
            ratios={"train": 0.5, "test": 0.5},
            shuffle=False,
            similarity_k=2,
            max_pairwise_records=9,
        ),
    )
    payload = result.to_dict()

    assert payload["shuffle"] is False
    assert payload["similarity_k"] == 2
    assert payload["max_pairwise_records"] == 9
    assert payload["assignment_strategy"] == ("connected_component_atomic_unit_greedy_target_error")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"ratios": {"train": 0.8, "test": 0.3}},
        {"ratios": {"train": 1.0}},
        {"ratios": {"train": True, "test": 0.0}},
        {"ratios": {"train": 0.5, "test": 0.5}, "seed": True},
        {"method": "group", "ratios": {"train": 0.5, "test": 0.5}},
        {
            "method": "similarity",
            "ratios": {"train": 0.5, "test": 0.5},
            "similarity_threshold": 1.1,
        },
    ],
)
def test_split_configuration_rejects_ambiguous_or_unsafe_values(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ConfigurationError):
        SplitConfig(**kwargs)  # type: ignore[arg-type]
