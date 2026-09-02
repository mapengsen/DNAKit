from __future__ import annotations

import io
from collections.abc import Mapping
from pathlib import Path

import pytest

from dnakit.core import (
    CompoundLocation,
    DNAAlphabet,
    DNAFeature,
    DNARecord,
    DNASequence,
    Interval,
    Strand,
    Strandedness,
    Topology,
)
from dnakit.exceptions import InputFormatError
from dnakit.io import ReadConfig, WriteConfig, read_one, write


def _record() -> DNARecord:
    return DNARecord(
        DNASequence(
            "ACGTNACG",
            alphabet=DNAAlphabet.IUPAC,
            topology=Topology.CIRCULAR,
            strandedness=Strandedness.DOUBLE,
        ),
        "TEST.1",
        description="auditable example",
        features=(
            DNAFeature(
                "misc_feature",
                CompoundLocation((Interval(0, 2), Interval(5, 8))),
                id="feature-1",
                label="joined",
                strand=Strand.REVERSE,
                score=1.25,
                phase=1,
                qualifiers={"note": "example", "tag": ["a", "b"]},
                source="manual",
            ),
        ),
        metadata={
            "genbank": {
                "locus": "TEST",
                "accession": "TEST",
                "version": "TEST.1",
                "division": "SYN",
                "date": "13-AUG-2026",
            }
        },
    )


def test_genbank_subset_semantic_roundtrip() -> None:
    stream = io.StringIO()

    result = write(_record(), stream, format="genbank")
    stream.seek(0)
    restored = read_one(stream, format="genbank")

    assert result.record_count == 1
    assert restored.id == "TEST.1"
    assert restored.sequence == _record().sequence
    assert restored.description == "auditable example"
    assert restored.features == _record().features
    metadata = restored.metadata["genbank"]
    assert isinstance(metadata, Mapping)
    assert metadata["codec"] == "dnakit.genbank.subset.v1"
    assert stream.getvalue().endswith("//\n")


def test_genbank_path_inference_and_atomic_write(tmp_path: Path) -> None:
    path = tmp_path / "record.gbk"

    write(_record(), path)

    assert read_one(path).sequence.to_string() == "ACGTNACG"
    with pytest.raises(FileExistsError):
        write(_record(), path)
    write(_record(), path, config=WriteConfig(overwrite=True))


@pytest.mark.parametrize(
    ("text", "code"),
    [
        (
            "LOCUS       X 4 bp DNA linear UNK 01-JAN-1980\n"
            "FEATURES             Location/Qualifiers\n"
            "     gene            <1..4\nORIGIN\n        1 acgt\n//\n",
            "UNSUPPORTED_GENBANK_LOCATION",
        ),
        (
            "LOCUS       X 5 bp DNA linear UNK 01-JAN-1980\nORIGIN\n        1 acgt\n//\n",
            "GENBANK_LENGTH_MISMATCH",
        ),
        (
            "LOCUS       X 4 bp DNA linear UNK 01-JAN-1980\nORIGIN\n        1 acgt\n",
            "TRUNCATED_GENBANK_RECORD",
        ),
    ],
)
def test_genbank_subset_rejects_unsupported_or_malformed_input(text: str, code: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_one(io.StringIO(text), format="genbank")

    assert exc_info.value.code == code


def test_genbank_rejects_gap_and_nested_qualifier_loss() -> None:
    feature = DNAFeature("gene", Interval(0, 1), qualifiers={"nested": {"x": 1}})
    record = DNARecord(DNASequence("A"), "x", features=(feature,))

    with pytest.raises(InputFormatError) as exc_info:
        write(record, io.StringIO(), format="genbank")

    assert exc_info.value.code == "UNSUPPORTED_GENBANK_QUALIFIER"


def test_genbank_rejects_missing_origin_and_record_line_limit() -> None:
    missing = "LOCUS       X 0 bp DNA linear UNK 01-JAN-1980\n//\n"
    with pytest.raises(InputFormatError) as missing_error:
        read_one(io.StringIO(missing), format="genbank")
    assert missing_error.value.code == "GENBANK_MISSING_ORIGIN"

    text = "LOCUS       X 1 bp DNA linear UNK 01-JAN-1980\nORIGIN\n        1 a\n//\n"
    with pytest.raises(InputFormatError) as limit_error:
        read_one(
            io.StringIO(text),
            format="genbank",
            config=ReadConfig(max_record_lines=1),
        )
    assert limit_error.value.code == "GENBANK_RECORD_LINE_LIMIT_EXCEEDED"

    bounded = io.StringIO()
    write(_record(), bounded, format="genbank")
    bounded.seek(0)
    with pytest.raises(InputFormatError) as symbol_error:
        read_one(
            bounded,
            format="genbank",
            config=ReadConfig(max_sequence_symbols=7),
        )
    assert symbol_error.value.code == "SEQUENCE_SYMBOL_LIMIT_EXCEEDED"


def test_genbank_semantics_match_biopython_for_simple_record() -> None:
    seqio = pytest.importorskip("Bio.SeqIO")
    stream = io.StringIO()
    write(DNARecord(DNASequence("ACGT"), "SIMPLE"), stream, format="genbank")

    stream.seek(0)
    parsed = seqio.read(stream, "genbank")

    assert parsed.id == "SIMPLE"
    assert str(parsed.seq) == "ACGT"
    assert parsed.annotations["topology"] == "linear"
