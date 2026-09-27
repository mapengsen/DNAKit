"""Real Goldman codec round trips and damaged-packet rejection.

These tests check digital encoding and integrity, not physical DNA storage or
sequencing-error correction. Optional Chamaeleo is never replaced by a fake codec.
"""

from __future__ import annotations

import copy
import hashlib
import json
from itertools import pairwise

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.genomics.storage import decode_dna_storage, encode_dna_storage


@pytest.fixture
def real_codec() -> None:
    pytest.importorskip("Chamaeleo.methods.fixed")


@pytest.fixture
def packet(real_codec: None) -> dict:
    return encode_dna_storage(bytes(range(64)).hex(), segment_bytes=16, verbose=False)


@pytest.mark.usefixtures("real_codec")
@pytest.mark.parametrize("segment_bytes", [1, 31, 256, 4096])
def test_all_byte_values_round_trip_across_segment_boundaries(segment_bytes: int) -> None:
    data = bytes(range(256)) * 17
    packet = encode_dna_storage(data.hex().upper(), segment_bytes=segment_bytes, verbose=False)
    assert packet["byte_count"] == len(data)
    assert packet["sha256"] == hashlib.sha256(data).hexdigest()
    assert packet["backend_version"]
    assert packet["error_correction"] == "none"
    assert "JSON" in packet["index_storage"]
    for segment in packet["segments"]:
        dna = segment["sequence"]
        assert set(dna) <= set("ACGT")
        # Goldman rotating ternary code forbids adjacent repeated nucleotides.
        assert all(left != right for left, right in pairwise(dna))
    decoded = decode_dna_storage(json.loads(json.dumps(packet)), verbose=False)
    assert decoded["data_hex"] == data.hex()
    assert decoded["integrity"] == "passed"


@pytest.mark.usefixtures("real_codec")
def test_empty_packet_has_verified_empty_payload() -> None:
    packet = encode_dna_storage("", verbose=False)
    assert packet["segments"] == []
    decoded = decode_dna_storage(packet, verbose=False)
    assert decoded["data_hex"] == ""
    assert decoded["byte_count"] == 0
    assert decoded["sha256"] == hashlib.sha256(b"").hexdigest()


def test_unordered_segments_are_restored_without_mutating_packet(packet: dict) -> None:
    packet["segments"].reverse()
    before = copy.deepcopy(packet)
    decoded = decode_dna_storage(packet, verbose=False)
    assert decoded["data_hex"] == bytes(range(64)).hex()
    assert packet == before


@pytest.mark.parametrize(
    "damage", ["missing", "duplicate", "wrong_length", "reindex", "swap_payload"]
)
def test_damaged_segment_collection_fails(packet: dict, damage: str) -> None:
    segments = packet["segments"]
    if damage == "missing":
        segments.pop()
    elif damage == "duplicate":
        segments[1]["index"] = segments[0]["index"]
    elif damage == "wrong_length":
        segments[0]["byte_count"] += 1
    elif damage == "reindex":
        segments[-1]["index"] += 1
    else:
        segments[0]["sequence"], segments[1]["sequence"] = (
            segments[1]["sequence"],
            segments[0]["sequence"],
        )
    with pytest.raises(InputFormatError):
        decode_dna_storage(packet, verbose=False)


@pytest.mark.parametrize("damage", ["homopolymer", "incomplete_codeword", "truncate", "checksum"])
def test_backend_silent_skips_and_partial_codewords_are_not_success(
    packet: dict, damage: str
) -> None:
    segment = packet["segments"][0]
    if damage == "homopolymer":
        segment["sequence"] = "A" + segment["sequence"][1:]
    elif damage == "incomplete_codeword":
        # One rotating base is a partial ternary codeword; Chamaeleo silently
        # ignores its leftover trit, so the adapter must verify canonical re-encoding.
        segment["sequence"] += next(base for base in "ACGT" if base != segment["sequence"][-1])
    elif damage == "truncate":
        segment["sequence"] = segment["sequence"][:-1]
    else:
        packet["sha256"] = "0" * 64
    with pytest.raises(InputFormatError):
        decode_dna_storage(packet, verbose=False)


@pytest.mark.parametrize("value", [None, 42, "0", "0g", "00 01", "0x00", "A" * 2_000_002])
def test_invalid_hex_is_rejected_before_optional_import(
    value: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unexpected_import(*args: object) -> None:
        pytest.fail("Invalid input must not load the optional backend")

    monkeypatch.setattr("dnakit.genomics.storage.require", unexpected_import)
    with pytest.raises(ConfigurationError):
        encode_dna_storage(value, verbose=False)  # type: ignore[arg-type]


@pytest.mark.parametrize("value", [True, 0, -1, 4097, 1.5, "32"])
def test_invalid_segment_size(value: object) -> None:
    with pytest.raises(ConfigurationError):
        encode_dna_storage("00", segment_bytes=value, verbose=False)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "field,value",
    [
        ("byte_count", True),
        ("byte_count", -1),
        ("sha256", "A" * 64),
        ("segments", None),
        ("scheme", "unknown"),
    ],
)
def test_packet_metadata_contract(packet: dict, field: str, value: object) -> None:
    packet[field] = value
    with pytest.raises(InputFormatError):
        decode_dna_storage(packet, verbose=False)
