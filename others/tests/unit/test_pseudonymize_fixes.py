"""Regression tests for the shared metadata pseudonymisation (it replaced two misleadingly named functions that
hashed six top-level keys with an unsalted SHA-256)."""

from __future__ import annotations

import copy
import hashlib
import hmac
import json

import pytest

import dnakit._candidate_runtime as runtime
from dnakit import knowledge_graph as kg
from dnakit import workflow_security as ws
from dnakit.exceptions import ConfigurationError, InputFormatError

FUNCTIONS = [ws.pseudonymize_metadata_fields, kg.pseudonymize_genomic_metadata_fields]
NAMES = ["pseudonymize_metadata_fields", "pseudonymize_genomic_metadata_fields"]
SALT = "unit-test-salt"


def run(function, data, **options):
    envelope = function(data, options=options or None, verbose=False)
    assert envelope["status"] == "ok"
    return envelope


def hmac_token(salt: str, message: str, chars: int = 16) -> str:
    return hmac.new(salt.encode(), message.encode(), hashlib.sha256).hexdigest()[:chars]


def test_old_names_are_gone():
    old_ws, old_kg = "remove_phi" + "_from_metadata", "deidentify" + "_genomic_metadata"  # spelled apart on purpose
    for module, old in ((ws, old_ws), (kg, old_kg)):
        assert not hasattr(module, old) and old not in module.__all__
        assert old not in runtime._REAL_OPERATIONS
    for module, new in ((ws, NAMES[0]), (kg, NAMES[1])):
        assert new in module.__all__ and new in module.__candidate_adapter_operations__
    assert runtime.run_candidate_operation("workflow_security", old_ws, {"name": "x"})["status"] == "not_implemented"


@pytest.mark.parametrize("function", FUNCTIONS)
def test_redact_nested_lists_and_case_variants(function):
    data = {
        "Patient_ID": "P-0042",
        "FULL NAME": "Ada Lovelace",
        "sample": "S1",
        "visits": [{"Date Of Birth": "1815-12-10", "site": "lab"}, {"contact": {"E-mail": "a@b.org", "n": 3}}],
        "dx": ["flu", {"mrn": 99}],
    }
    result = run(function, data, scan_text=False)["result"]
    assert result["metadata"] == {
        "Patient_ID": "[REDACTED]",
        "FULL NAME": "[REDACTED]",
        "sample": "S1",
        "visits": [{"Date Of Birth": "[REDACTED]", "site": "lab"}, {"contact": {"E-mail": "[REDACTED]", "n": 3}}],
        "dx": ["flu", {"mrn": "[REDACTED]"}],
    }
    assert result["mode"] == "redact" and result["n_redacted"] == 5
    assert result["redacted_paths"] == [
        "$.Patient_ID", "$['FULL NAME']", "$.visits[0]['Date Of Birth']", "$.visits[1].contact['E-mail']", "$.dx[1].mrn",
    ]


def test_sensitive_key_with_container_value_is_replaced_whole():
    result = run(ws.pseudonymize_metadata_fields, {"address": {"street": "1 Main St", "city": "X"}, "id": 1})["result"]
    assert result["metadata"] == {"address": "[REDACTED]", "id": 1}
    assert result["redacted_paths"] == ["$.address"]


def test_extra_keys_and_keep_exceptions():
    data = {"clinic": "North", "name": "Ada", "Sample_Owner": "bob"}
    result = run(ws.pseudonymize_metadata_fields, data, keys=["clinic", "sample owner"], keep=["NAME"])["result"]
    assert result["metadata"] == {"clinic": "[REDACTED]", "name": "Ada", "Sample_Owner": "[REDACTED]"}
    assert run(ws.pseudonymize_metadata_fields, data, keys="clinic")["result"]["metadata"]["clinic"] == "[REDACTED]"


def test_scan_text_replaces_patterns_anywhere():
    data = {"note": "mail ada@example.org or call +44 20 7946 0958 / (415) 555-0100; ssn 123-45-6789 ok", "rows": ["x@y.io"]}
    result = run(ws.pseudonymize_metadata_fields, data)["result"]
    assert result["metadata"]["note"] == "mail [REDACTED:email] or call [REDACTED:phone] / [REDACTED:phone]; ssn [REDACTED:ssn] ok"
    assert result["metadata"]["rows"] == ["[REDACTED:email]"]
    assert result["n_text_matches"] == 5 and result["n_key_matches"] == 0 and result["n_redacted"] == 5
    assert result["text_matches_by_pattern"] == {"email": 2, "phone": 2, "ssn": 1}
    assert result["redacted_paths"] == ["$.note", "$.rows[0]"]
    # harmless text and plain numbers are left alone, and scan_text=False disables the scan
    plain = {"note": "GC 0.41, reads 12345, 2024-01-15", "n": 5551234567}
    assert run(ws.pseudonymize_metadata_fields, plain)["result"]["metadata"] == plain
    off = run(ws.pseudonymize_metadata_fields, data, scan_text=False)
    assert off["result"]["metadata"] == data and off["result"]["n_redacted"] == 0
    assert any("scan_text is disabled" in w for w in off["warnings"])


def test_hash_mode_matches_independent_hmac_and_depends_on_salt():
    data = {"name": "Ada", "n": 7, "meta": {"phone": 5551234, "address": {"city": "X"}}, "Zip": "02139"}
    result = run(ws.pseudonymize_metadata_fields, data, mode="hash", salt=SALT, scan_text=False)["result"]
    metadata = result["metadata"]
    assert metadata["name"] == hmac_token(SALT, "Ada")
    assert metadata["Zip"] == hmac_token(SALT, "02139")
    assert metadata["meta"]["phone"] == hmac_token(SALT, "5551234")  # non-strings: canonical JSON text
    assert metadata["meta"]["address"] == hmac_token(SALT, '{"city":"X"}')
    assert metadata["n"] == 7 and len(metadata["name"]) == 16
    assert metadata["name"] != hashlib.sha256(b"Ada").hexdigest()[:16]  # not the old unsalted token
    other = run(ws.pseudonymize_metadata_fields, data, mode="hash", salt="another-salt!", scan_text=False)["result"]["metadata"]
    assert other["name"] != metadata["name"] and other["name"] == hmac_token("another-salt!", "Ada")
    again = run(kg.pseudonymize_genomic_metadata_fields, data, mode="hash", salt=SALT, scan_text=False)["result"]["metadata"]
    assert again == metadata  # deterministic, and both functions share the implementation
    long = run(ws.pseudonymize_metadata_fields, {"name": "Ada"}, mode="hash", salt=SALT, digest_chars=64)["result"]
    assert long["metadata"]["name"] == hmac_token(SALT, "Ada", 64)
    assert any("linkable" in w for w in long["warnings"])


def test_drop_mode_removes_keys():
    data = {"name": "Ada", "keep": 1, "visits": [{"email": "a@b.org", "x": 2}], "note": "mail a@b.org"}
    result = run(ws.pseudonymize_metadata_fields, data, mode="drop")["result"]
    assert result["metadata"] == {"keep": 1, "visits": [{"x": 2}], "note": "mail [REDACTED:email]"}
    assert result["redacted_paths"] == ["$.name", "$.visits[0].email", "$.note"] and result["n_redacted"] == 3


@pytest.mark.parametrize("function", FUNCTIONS)
def test_input_forms(function):
    as_list = [{"name": "A"}, {"name": "B", "k": 1}]
    assert run(function, as_list)["result"]["metadata"] == [{"name": "[REDACTED]"}, {"name": "[REDACTED]", "k": 1}]
    assert run(function, as_list)["input_summary"] == {"record_count": 2, "n_keys": 3}
    text = json.dumps({"Name": "A", "x": [1]})
    assert run(function, text)["result"]["metadata"] == {"Name": "[REDACTED]", "x": [1]}
    assert run(function, json.dumps(as_list))["result"]["metadata"][1] == {"name": "[REDACTED]", "k": 1}
    assert run(function, [])["result"]["metadata"] == []


@pytest.mark.parametrize("function", FUNCTIONS)
@pytest.mark.parametrize("bad", ["not json", "42", '"just text"', "[1, 2]", 42, None, [1, 2], [{"a": 1}, "x"], b"\xff\xfe"])
def test_invalid_input_raises_input_format_error(function, bad):
    with pytest.raises(InputFormatError) as caught:
        function(bad, verbose=False)
    assert caught.value.code in {"INVALID_METADATA_JSON", "INVALID_METADATA_INPUT"} and caught.value.hint


def test_unsupported_value_and_cycle_raise_without_leaking():
    with pytest.raises(InputFormatError) as caught:
        ws.pseudonymize_metadata_fields({"a": {1, 2}}, verbose=False)
    assert caught.value.code == "UNSUPPORTED_METADATA_VALUE"
    cyclic: dict = {}
    cyclic["self"] = cyclic
    with pytest.raises(InputFormatError) as caught:
        ws.pseudonymize_metadata_fields(cyclic, verbose=False)
    assert caught.value.code == "METADATA_TOO_DEEP"


@pytest.mark.parametrize(
    "options, code",
    [
        ({"mode": "hash"}, "PSEUDONYMIZE_SALT_REQUIRED"),
        ({"mode": "hash", "salt": ""}, "PSEUDONYMIZE_SALT_REQUIRED"),
        ({"mode": "hash", "salt": "short"}, "PSEUDONYMIZE_SALT_REQUIRED"),
        ({"mode": "hash", "salt": "        "}, "PSEUDONYMIZE_SALT_REQUIRED"),
        ({"mode": "hash", "salt": 12345678}, "PSEUDONYMIZE_SALT_REQUIRED"),
        ({"mode": "hash", "salt": SALT, "digest_chars": 65}, "INVALID_DIGEST_CHARS"),
        ({"mode": "hash", "salt": SALT, "digest_chars": 2}, "INVALID_DIGEST_CHARS"),
        ({"mode": "hash", "salt": SALT, "digest_chars": True}, "INVALID_DIGEST_CHARS"),
        ({"mode": "sha"}, "INVALID_PSEUDONYMIZE_MODE"),
        ({"keys": [1]}, "INVALID_PSEUDONYMIZE_KEYS"),
        ({"keep": ["--"]}, "INVALID_PSEUDONYMIZE_KEYS"),
        ({"scan_text": "yes"}, "INVALID_PSEUDONYMIZE_OPTIONS"),
    ],
)
def test_option_errors(options, code):
    with pytest.raises(ConfigurationError) as caught:
        ws.pseudonymize_metadata_fields({"name": "Ada"}, options=options, verbose=False)
    assert caught.value.code == code and caught.value.hint


def test_options_must_be_mapping():
    with pytest.raises(ConfigurationError):
        ws.pseudonymize_metadata_fields({"name": "Ada"}, options=["mode"], verbose=False)


@pytest.mark.parametrize("function", FUNCTIONS)
def test_input_is_not_mutated(function):
    data = {"name": "Ada", "visits": [{"email": "a@b.org", "note": "x@y.io"}], "tags": ["a"]}
    snapshot = copy.deepcopy(data)
    for options in ({}, {"mode": "drop"}, {"mode": "hash", "salt": SALT}):
        run(function, data, **options)
    assert data == snapshot


@pytest.mark.parametrize("function", FUNCTIONS)
@pytest.mark.parametrize("options", [{}, {"mode": "drop"}, {"mode": "hash", "salt": SALT}])
def test_no_plaintext_in_output_warnings_summary_or_paths(function, options):
    secrets = ["Ada Lovelace", "P-0042", "ada@example.org", "555-0100", "123-45-6789", "Graceland", "1815-12-10"]
    data = {
        "name": "Ada Lovelace", "patient_id": "P-0042", "dob": "1815-12-10",
        "contact": {"note": "mail ada@example.org tel (415) 555-0100 ssn 123-45-6789"},
        "rows": [{"address": "Graceland"}],
    }
    envelope = run(function, data, **options)
    dumped = json.dumps({k: v for k, v in envelope.items() if k != "result"}) + json.dumps(
        {k: v for k, v in envelope["result"].items() if k != "metadata"}
    )
    dumped_all = json.dumps(envelope)
    for secret in secrets:
        assert secret not in dumped_all, secret
    assert SALT not in dumped_all
    assert envelope["input_summary"] == {"record_count": 1, "n_keys": 7}
    assert "sha256" not in envelope["input_summary"] and "text_length" not in envelope["input_summary"]
    assert dumped  # bookkeeping is present


@pytest.mark.parametrize("function", FUNCTIONS)
def test_warnings_state_the_limits(function):
    envelope = run(function, {"name": "Ada"})
    joined = " ".join(envelope["result"]["warnings"])
    assert envelope["warnings"] == envelope["result"]["warnings"]
    for phrase in ("key-name and pattern based", "NOT a guarantee", "not shifted or generalised", "ages over 89",
                   "small geographic units", "Free text may still identify", "Review the output"):
        assert phrase in joined
    assert envelope["operation"] in NAMES and envelope["module"] in {"workflow_security", "knowledge_graph"}


def test_envelope_override_default_is_unchanged():
    old = runtime.candidate_ok("m", "op", {"a": 1}, 1)
    assert set(old["input_summary"]) == {"record_count", "text_length", "sequence_bases", "lengths", "sha256"}
    new = runtime.candidate_ok("m", "op", {"a": 1}, 1, input_summary={"record_count": 1, "n_keys": 1})
    assert new["input_summary"] == {"record_count": 1, "n_keys": 1}
