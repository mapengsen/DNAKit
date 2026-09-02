"""Unit tests for the STD-001..007 and STD-009 normalization contract."""

from __future__ import annotations

import json
import math

import pytest

from dnakit.core.enums import DNAAlphabet, GapKind
from dnakit.core.gap import Gap
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError, InvalidAlphabetError
from dnakit.standardize import (
    AmbiguityPolicy,
    NormalizationConfig,
    UPolicy,
    normalize,
)


def issue_codes(result: object) -> set[str]:
    return {issue.code for issue in result.issues}  # type: ignore[attr-defined]


def test_normalize_records_case_whitespace_u_and_original_coordinates() -> None:
    result = normalize(" a\ncu\t")

    assert result.is_valid
    assert result.sequence is not None
    assert result.sequence.symbols == "AC"
    assert result.algorithm_version == "std-normalize-v2"
    assert result.raw_input.content == " a\ncu\t"
    assert result.raw_input.character_count == 6
    assert [change.operation for change in result.changes] == [
        "remove_whitespace",
        "uppercase",
        "remove_whitespace",
        "uppercase",
        "uppercase",
        "delete_u",
        "remove_whitespace",
    ]
    assert result.u_positions[0].absolute_offset == 4
    assert result.raw_input.to_dict() == {
        "input_type": "str",
        "sha256": result.raw_input.sha256,
        "character_count": 6,
    }
    assert result.raw_input.to_dict(include_content=True)["content"] == " a\ncu\t"
    serialized = result.to_dict()
    assert "content" not in serialized["raw_input"]
    json.dumps(serialized)


def test_default_deletes_other_characters_but_audits_original_positions() -> None:
    fullwidth_a = "\uff21"
    result = normalize(f"AC\u200bG{fullwidth_a}T")

    assert result.sequence is not None
    assert result.sequence.symbols == "ACGT"
    assert result.invalid_symbols == ()
    invisible = next(change for change in result.changes if change.operation == "remove_invisible")
    assert invisible.position.absolute_offset == 2
    deleted = next(change for change in result.changes if change.operation == "delete_other")
    assert deleted.before == fullwidth_a
    assert deleted.position.absolute_offset == 4


def test_keep_other_retains_invalid_character_without_fabricating_dna() -> None:
    result = normalize("ACXG", keep_other=True)

    assert result.sequence is None
    assert result.normalized_parts == ("ACXG",)
    assert result.invalid_symbols[0].symbol == "X"
    assert result.invalid_symbols[0].positions[0].absolute_offset == 2
    assert "STD_INVALID_SYMBOL" in issue_codes(result)


def test_simple_character_flags_have_requested_defaults_and_deletion_modes() -> None:
    default = normalize("ANURX-")
    deleted_ambiguity = normalize("ANURX-", keep_ambiguous=False)
    kept_u = normalize("AU", keep_u=True)

    assert default.sequence is not None
    assert default.sequence.symbols == "ANR"
    assert default.sequence.alphabet is DNAAlphabet.IUPAC
    assert default.config.keep_ambiguous is True
    assert default.config.keep_u is False
    assert default.config.keep_other is False
    assert [change.operation for change in default.changes] == [
        "delete_u",
        "delete_other",
        "delete_other",
    ]
    assert default.u_positions[0].absolute_offset == 2
    uracil_step = next(step for step in default.steps if step.name == "uracil")
    other_step = next(step for step in default.steps if step.name == "other_characters")
    assert uracil_step.parameters == (("policy", "delete"),)
    assert other_step.change_count == 2
    assert deleted_ambiguity.sequence is not None
    assert deleted_ambiguity.sequence.symbols == "A"
    assert [change.operation for change in deleted_ambiguity.changes] == [
        "delete_ambiguity",
        "delete_u",
        "delete_ambiguity",
        "delete_other",
        "delete_other",
    ]
    assert kept_u.sequence is None
    assert kept_u.normalized_parts == ("AU",)


def test_advanced_u_replace_policy_remains_available() -> None:
    result = normalize("AU", config=NormalizationConfig(u_policy=UPolicy.REPLACE))

    assert result.sequence is not None
    assert result.sequence.symbols == "AT"
    assert result.changes[0].operation == "replace_u"


@pytest.mark.parametrize("policy", [UPolicy.WARN, UPolicy.KEEP])
def test_retained_u_never_constructs_a_dna_sequence(policy: UPolicy) -> None:
    result = normalize("AU", config=NormalizationConfig(u_policy=policy))

    assert result.sequence is None
    assert result.normalized_parts == ("AU",)
    assert "STD_INVALID_SYMBOL" in issue_codes(result)
    assert ("STD_U_PRESENT" in issue_codes(result)) is (policy is UPolicy.WARN)


def test_u_error_can_be_returned_or_raised() -> None:
    result = normalize("AU", config=NormalizationConfig(u_policy=UPolicy.ERROR))
    assert result.sequence is None
    assert {"STD_U_PRESENT", "STD_INVALID_SYMBOL"} <= issue_codes(result)

    with pytest.raises(InvalidAlphabetError) as error:
        normalize(
            "AU",
            config=NormalizationConfig(u_policy=UPolicy.ERROR, raise_on_error=True),
        )
    assert error.value.code == "INVALID_ALPHABET"


def test_strict_and_iupac_alphabets_are_explicit() -> None:
    strict = normalize("ACGN", config=NormalizationConfig(alphabet=DNAAlphabet.STRICT))
    iupac = normalize("ACGN")

    assert strict.sequence is None
    assert strict.invalid_symbols[0].positions[0].absolute_offset == 3
    assert iupac.sequence is not None
    assert iupac.sequence.alphabet is DNAAlphabet.IUPAC
    assert iupac.ambiguity.count("N") == 1


def test_all_iupac_symbols_have_counts_positions_and_probability_resolution() -> None:
    symbols = "RYSWKMBDHVN"
    result = normalize(
        symbols,
        config=NormalizationConfig(
            alphabet=DNAAlphabet.IUPAC,
            ambiguity_policy=AmbiguityPolicy.PROBABILITY,
        ),
    )

    assert result.is_valid
    assert result.ambiguity.total_count == len(symbols)
    assert result.ambiguity.fraction == 1.0
    assert tuple(item.symbol_offset for item in result.ambiguity.occurrences) == tuple(
        range(len(symbols))
    )
    assert len(result.ambiguity.probability_resolutions) == len(symbols)
    assert all(
        math.isclose(sum(value for _, value in item.probabilities), 1.0)
        for item in result.ambiguity.probability_resolutions
    )


def test_ambiguity_error_and_mask_are_distinct_policies() -> None:
    errored = normalize(
        "ARY",
        config=NormalizationConfig(
            alphabet=DNAAlphabet.IUPAC,
            ambiguity_policy=AmbiguityPolicy.ERROR,
        ),
    )
    masked = normalize(
        "ARY",
        config=NormalizationConfig(
            alphabet=DNAAlphabet.IUPAC,
            ambiguity_policy=AmbiguityPolicy.MASK,
        ),
    )

    assert errored.sequence is None
    assert "STD_AMBIGUITY_NOT_ALLOWED" in issue_codes(errored)
    assert masked.sequence is not None
    assert masked.sequence.symbols == "ANN"
    assert [item.operation for item in masked.changes] == [
        "mask_ambiguity",
        "mask_ambiguity",
    ]


def test_explicit_gap_is_preserved_and_raw_generator_is_materialized_once() -> None:
    gap = Gap(None, kind=GapKind.SCAFFOLD, crossable=False, evidence=("manual",))
    source_parts: list[str | Gap] = [" ac ", gap, "gt"]
    source = (part for part in source_parts)
    result = normalize(source)

    assert result.is_valid
    assert result.sequence is not None
    assert result.sequence.parts == ("AC", gap, "GT")
    assert result.raw_input.content == (" ac ", gap, "gt")
    assert result.sequence.coordinate_span is None
    gap_step = next(step for step in result.steps if step.name == "explicit_gap")
    assert gap_step.status == "applied"


def test_nested_gap_metadata_has_stable_digest_and_json_snapshot() -> None:
    source: list[str | Gap] = [
        "A",
        Gap(
            2,
            metadata={"nested": {"items": [1, {"accepted": True}]}},
        ),
        "T",
    ]

    first = normalize(source)
    second = normalize(source)
    payload = first.to_dict(include_raw_content=True)

    assert first.raw_input.sha256 == second.raw_input.sha256
    assert payload["raw_input"]["content"][1]["gap"]["metadata"] == {
        "nested": {"items": [1, {"accepted": True}]}
    }
    json.dumps(payload, ensure_ascii=False, sort_keys=True)


def test_gap_can_be_rejected_without_losing_the_normalized_parts() -> None:
    gap = Gap(10)
    result = normalize(["AC", gap, "GT"], config=NormalizationConfig(allow_gaps=False))

    assert result.sequence is None
    assert result.normalized_parts == ("AC", gap, "GT")
    assert "STD_GAP_NOT_ALLOWED" in issue_codes(result)


def test_dna_sequence_normalization_is_audited_no_op_and_preserves_identity() -> None:
    sequence = DNASequence("ACGN", alphabet=DNAAlphabet.IUPAC)
    result = normalize(sequence)

    assert result.sequence is sequence
    assert not result.was_modified
    assert result.raw_input.content is sequence
    alphabet_step = next(step for step in result.steps if step.name == "alphabet")
    assert alphabet_step.parameters == (("mode", "iupac"),)


def test_empty_dna_sequence_no_op_preserves_identity_and_core_parts() -> None:
    sequence = DNASequence("")
    result = normalize(sequence)

    assert result.sequence is sequence
    assert result.normalized_parts == ()
    assert not result.was_modified


def test_empty_raw_sequence_is_normalized_but_left_for_qc_policy() -> None:
    result = normalize("")

    assert result.is_valid
    assert result.sequence is not None
    assert result.sequence.parts == ()


@pytest.mark.parametrize("bad_prior", [float("nan"), float("inf"), -1.0])
def test_nonfinite_or_negative_base_prior_is_a_configuration_error(bad_prior: float) -> None:
    with pytest.raises(ConfigurationError) as error:
        NormalizationConfig(base_priors={"A": bad_prior, "C": 1, "G": 1, "T": 1})
    assert error.value.code == "INVALID_BASE_PRIORS"


def test_normalization_config_rejects_truthy_flags_and_invalid_audit_fields() -> None:
    with pytest.raises(ConfigurationError):
        NormalizationConfig(uppercase=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(remove_whitespace="false")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(keep_ambiguous=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(keep_u="false")  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(keep_other=[])  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(allow_gaps=0)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(raise_on_error=[])  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(operator=7)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        NormalizationConfig(run_id="")
    with pytest.raises(ConfigurationError):
        NormalizationConfig(removable_separators="-")  # type: ignore[arg-type]


def test_boolean_base_prior_is_not_silently_converted_to_one() -> None:
    with pytest.raises(ConfigurationError) as error:
        NormalizationConfig(base_priors={"A": True, "C": 1, "G": 1, "T": 1})

    assert error.value.code == "INVALID_BASE_PRIORS"


def test_normalize_rejects_falsy_non_config_and_case_colliding_priors() -> None:
    invalid_configs: tuple[object, ...] = ({}, [], 0, False)
    for value in invalid_configs:
        with pytest.raises(ConfigurationError):
            normalize("AC", config=value)  # type: ignore[arg-type]

    with pytest.raises(ConfigurationError) as error:
        NormalizationConfig(base_priors={"a": 0.1, "A": 0.2, "C": 0.3, "G": 0.3, "T": 0.1})
    assert error.value.code == "INVALID_BASE_PRIORS"


def test_normalize_rejects_invalid_or_conflicting_direct_keep_arguments() -> None:
    with pytest.raises(ConfigurationError) as flag_error:
        normalize("AC", keep_u=1)  # type: ignore[arg-type]
    assert flag_error.value.code == "INVALID_NORMALIZATION_FLAG"

    with pytest.raises(ConfigurationError) as conflict_error:
        normalize("AC", keep_other=True, config=NormalizationConfig())
    assert conflict_error.value.code == "NORMALIZATION_ARGUMENT_CONFLICT"


def test_undecodable_bytes_raise_instead_of_replacing_source_data() -> None:
    with pytest.raises(UnicodeDecodeError):
        normalize(b"AC\xffGT")
