"""The envelope's ``input_summary`` counts records and bases by the real format of the input."""

from __future__ import annotations

import hashlib

import pytest

from dnakit._candidate_runtime import _input_summary, candidate_not_implemented, candidate_ok


@pytest.mark.parametrize(
    ("value", "records", "bases", "lengths"),
    [
        # two reads of 4 and 2 bases; the quality letters must not be counted as bases
        ("@r1\nACGT\n+\nIIII\n@r2\nGG\n+\nII\n", 2, 6, [4, 2]),
        # a multi-line FASTA record is the sum of its lines: 4 + 2 = 6, and 4
        (">a\nACGT\nAC\n>b\nGGAA\n", 2, 10, [6, 4]),
        # a Newick tree has one record and no sequence; its labels are not bases
        ("((A:0.1,B:0.2):0.3,C:0.4);", 1, 0, []),
        (None, 0, 0, []),  # the word "None" used to count as 4 bases
        ("run1.pod5", 1, 0, []),  # a file name is not a sequence
        ("hello world", 1, 0, []),  # prose is not a sequence
        ("ACGT\nGGAA", 2, 8, [4, 4]),
        ([100, 60, 40], 3, 200, [100, 60, 40]),  # contig lengths
        ([{"id": "a", "sequence": "ACGT"}, {"id": "b", "sequence": "GG"}], 2, 6, [4, 2]),
        ({"a": 1}, 1, 0, []),
        (b"\xff\xfe\x00binary", 1, 0, []),
        (7, 1, 0, []),
    ],
)
def test_input_summary_counts_by_format(value, records, bases, lengths):
    summary = _input_summary(value)
    assert summary["record_count"] == records
    assert summary["sequence_bases"] == bases
    assert summary["lengths"] == lengths


def test_input_summary_hash_and_length_are_of_the_input():
    text = "@r1\nACGT\n+\nIIII\n"
    summary = _input_summary(text)
    assert summary["text_length"] == len(text)
    assert summary["sha256"] == hashlib.sha256(text.encode()).hexdigest()
    assert set(summary) == {"record_count", "text_length", "sequence_bases", "lengths", "sha256"}


def test_envelopes_use_the_new_summary():
    ok = candidate_ok("m", "op", "@r1\nACGT\n+\nIIII\n", {"x": 1})
    assert ok["input_summary"]["record_count"] == 1 and ok["input_summary"]["sequence_bases"] == 4
    ni = candidate_not_implemented("m", "op", "run1.pod5")
    assert ni["input_summary"]["sequence_bases"] == 0 and ni["input_summary"]["record_count"] == 1
