"""Human-verifiable tests for integer and one-hot DNA encodings."""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap
from dnakit.exceptions import (
    ConfigurationError,
    InvalidAlphabetError,
    UnknownLengthError,
    UnsupportedGapOperationError,
)
from dnakit.fingerprints import integer_encode, one_hot_encode


def test_integer_encoding_has_a_fixed_a_c_g_t_codebook_and_record_id() -> None:
    result = integer_encode(DNARecord(DNASequence("ACGT"), "record-1"))

    assert result.values == (0, 1, 2, 3)
    assert result.codebook == {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
        "R": 4,
        "Y": 5,
        "S": 6,
        "W": 7,
        "K": 8,
        "M": 9,
        "B": 10,
        "D": 11,
        "H": 12,
        "V": 13,
        "N": 14,
    }
    assert result.sequence_id == "record-1"
    assert result.schema_version == "dnakit.integer.v1"
    assert result.dimension == 15


def test_integer_encoding_makes_each_ambiguity_policy_explicit() -> None:
    sequence = DNASequence("ARN", alphabet=DNAAlphabet.IUPAC)

    iupac = integer_encode(sequence)
    sentinel = integer_encode(sequence, ambiguity_policy="sentinel")

    assert iupac.values == (0, 4, 14)
    assert sentinel.values == (0, -2, -2)
    assert sentinel.codebook == {"A": 0, "C": 1, "G": 2, "T": 3, "<IUPAC>": -2}
    assert sentinel.encoded_ambiguity_count == 2
    with pytest.raises(InvalidAlphabetError, match="Ambiguous IUPAC"):
        integer_encode(sequence, ambiguity_policy="error")


def test_integer_gap_policies_cover_known_and_unknown_spans() -> None:
    known = DNASequence(["AC", Gap(2), "GT"])
    unknown = DNASequence(["A", Gap(None), "T"])

    expanded = integer_encode(known, gap_policy="expand")
    omitted = integer_encode(known, gap_policy="omit")

    assert expanded.values == (0, 1, -1, -1, 2, 3)
    assert expanded.expanded_gap_length == 2
    assert expanded.codebook["<GAP>"] == -1
    assert omitted.values == (0, 1, 2, 3)
    assert omitted.omitted_gap_count == 1
    assert integer_encode(unknown, gap_policy="omit").unknown_gap_count == 1
    with pytest.raises(UnknownLengthError, match="unknown-length Gap"):
        integer_encode(unknown, gap_policy="expand")
    with pytest.raises(UnsupportedGapOperationError, match="explicit Gaps"):
        integer_encode(known)


def test_one_hot_canonical_and_custom_column_order_are_deterministic() -> None:
    standard = one_hot_encode(DNASequence("ACGT"))
    reversed_order = one_hot_encode(DNASequence("AT"), base_order=("T", "G", "C", "A"))

    assert standard.feature_names == ("A", "C", "G", "T")
    assert standard.values == (
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    )
    assert standard.dimension == 4
    assert reversed_order.values == ((0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 0.0))


def test_one_hot_fractional_zero_and_error_ambiguity_policies() -> None:
    sequence = DNASequence("RN", alphabet=DNAAlphabet.IUPAC)

    fractional = one_hot_encode(sequence, ambiguity_policy="fractional")
    zero = one_hot_encode(sequence, ambiguity_policy="zero")

    assert fractional.values == ((0.5, 0.0, 0.5, 0.0), (0.25, 0.25, 0.25, 0.25))
    assert fractional.encoded_ambiguity_count == 2
    assert zero.values == ((0.0, 0.0, 0.0, 0.0),) * 2
    with pytest.raises(InvalidAlphabetError):
        one_hot_encode(sequence)


def test_one_hot_gap_expansion_uses_zero_rows_and_empty_input_is_supported() -> None:
    known = DNASequence(["A", Gap(2), "T"])
    unknown = DNASequence(["A", Gap(None), "T"])

    expanded = one_hot_encode(known, gap_policy="expand")
    omitted = one_hot_encode(unknown, gap_policy="omit")
    empty = one_hot_encode(DNASequence(""))

    assert expanded.values == (
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    )
    assert expanded.expanded_gap_length == 2
    assert omitted.output_length == 2
    assert empty.values == ()
    assert empty.output_length == 0
    with pytest.raises(UnknownLengthError):
        one_hot_encode(unknown, gap_policy="expand")


def test_positional_encodings_enforce_an_explicit_output_size_limit() -> None:
    sequence = DNASequence(["A", Gap(10), "T"])

    with pytest.raises(ConfigurationError) as integer_error:
        integer_encode(sequence, gap_policy="expand", max_output_length=5)
    assert integer_error.value.code == "ENCODING_SIZE_LIMIT"
    with pytest.raises(ConfigurationError) as one_hot_error:
        one_hot_encode(sequence, gap_policy="expand", max_output_length=5)
    assert one_hot_error.value.code == "ENCODING_SIZE_LIMIT"
    assert integer_encode(sequence, gap_policy="omit", max_output_length=2).values == (0, 3)
    with pytest.raises(ConfigurationError):
        integer_encode(DNASequence("A"), max_output_length=True)


@pytest.mark.parametrize(
    ("function_name", "kwargs"),
    [
        ("integer", {"ambiguity_policy": "guess"}),
        ("integer", {"gap_policy": "mask"}),
        ("one_hot", {"base_order": ("A", "A", "G", "T")}),
        ("one_hot", {"base_order": "ACGT"}),
    ],
)
def test_encodings_reject_invalid_configuration(
    function_name: str,
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ConfigurationError):
        if function_name == "integer":
            integer_encode(DNASequence("A"), **kwargs)  # type: ignore[arg-type]
        else:
            one_hot_encode(DNASequence("A"), **kwargs)  # type: ignore[arg-type]


def test_encoding_results_are_immutable_and_json_serializable() -> None:
    result = one_hot_encode(DNASequence("AN", alphabet=DNAAlphabet.IUPAC), ambiguity_policy="zero")

    with pytest.raises(FrozenInstanceError):
        result.output_length = 100  # type: ignore[misc]
    assert json.loads(json.dumps(result.to_dict()))["ambiguity_policy"] == "zero"
