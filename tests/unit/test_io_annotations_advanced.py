from __future__ import annotations

import io

import pytest

from dnakit.core import GapKind, Interval, Strand
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.io import (
    AGPComponent,
    AGPGap,
    read_agp,
    read_bed,
    read_gff3,
    write_agp,
    write_bed,
    write_gff3,
)


def test_gff3_roundtrip_and_coordinate_conversion() -> None:
    text = (
        "##gff-version 3\n"
        "##sequence-region chr1 1 20\n"
        "chr1\ttest\tCDS\t2\t8\t0.5\t-\t1\tID=g1;Name=alpha;Alias=a,b\n"
    )

    document = read_gff3(io.StringIO(text))
    feature = document.entries[0].feature

    assert feature.location == Interval(1, 8)
    assert feature.strand is Strand.REVERSE
    assert feature.id == "g1"
    assert feature.qualifiers["Alias"] == ("a", "b")
    output = io.StringIO()
    assert write_gff3(document, output) == 1
    assert read_gff3(io.StringIO(output.getvalue())).entries == document.entries


@pytest.mark.parametrize(
    ("text", "code"),
    [
        ("chr1\tx\tgene\t1\t2\t.\t+\t.\t.\n", "GFF3_VERSION_REQUIRED"),
        ("##gff-version 3\nchr1\tx\tgene\t1\t2\t.\t+\t.\n", "INVALID_GFF3_COLUMNS"),
        ("##gff-version 3\n##FASTA\n", "UNSUPPORTED_GFF3_FASTA"),
    ],
)
def test_gff3_rejects_missing_version_columns_and_embedded_fasta(text: str, code: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_gff3(io.StringIO(text))

    assert exc_info.value.code == code


def test_gff3_accepts_version_revision_and_requires_cds_phase() -> None:
    document = read_gff3(
        io.StringIO("##gff-version 3.1.26\nchr1\tx\tCDS\t1\t3\t.\t+\t0\tID=cds1\n")
    )
    assert document.entries[0].feature.phase == 0

    with pytest.raises(InputFormatError) as phase_error:
        read_gff3(io.StringIO("##gff-version 3\nchr1\tx\tCDS\t1\t3\t.\t+\t.\t.\n"))
    assert phase_error.value.code == "INVALID_GFF3_PHASE"


def test_annotation_resource_limits_are_enforced() -> None:
    text = "##gff-version 3\nchr1\tx\tgene\t1\t2\t.\t+\t.\t.\n"

    with pytest.raises(ConfigurationError) as count_error:
        read_gff3(io.StringIO(text), max_records=1 - 1)
    assert count_error.value.code == "INVALID_ANNOTATION_LIMIT"

    with pytest.raises(InputFormatError) as line_error:
        read_gff3(io.StringIO(text), max_line_length=10)
    assert line_error.value.code == "LINE_TOO_LONG"


def test_bed6_roundtrip_and_bed12_explicit_rejection() -> None:
    document = read_bed(io.StringIO("track name=demo\nchr1\t2\t9\tregion-1\t500\t+\n"))

    assert document.entries[0].feature.location == Interval(2, 9)
    assert document.entries[0].feature.score == 500.0
    output = io.StringIO()
    assert write_bed(document, output) == 1
    assert read_bed(io.StringIO(output.getvalue())).entries == document.entries

    with pytest.raises(InputFormatError) as exc_info:
        read_bed(io.StringIO("chr1\t0\t10\tx\t0\t+\t0\t10\t0\t1\t10\t0\n"))
    assert exc_info.value.code == "UNSUPPORTED_BED_COLUMNS"


def test_agp_components_and_gaps_roundtrip() -> None:
    text = (
        "##agp-version 2.1\n"
        "scaf1\t1\t4\t1\tW\tcontig1\t1\t4\t+\n"
        "scaf1\t5\t7\t2\tN\t3\tscaffold\tyes\tpaired-ends\n"
        "scaf1\t8\t107\t3\tU\t100\tcontig\tno\tna\n"
    )

    document = read_agp(io.StringIO(text))

    assert isinstance(document.entries[0], AGPComponent)
    assert isinstance(document.entries[1], AGPGap)
    assert document.entries[1].gap.kind is GapKind.SCAFFOLD
    assert document.entries[1].gap.length == 3
    assert isinstance(document.entries[2], AGPGap)
    assert document.entries[2].gap.length is None
    output = io.StringIO()
    assert write_agp(document, output) == 3
    assert read_agp(io.StringIO(output.getvalue())) == document


@pytest.mark.parametrize(
    "text",
    [
        "obj\t2\t4\t1\tW\tc\t1\t3\t+\n",
        "obj\t1\t4\t1\tW\tc\t1\t3\t+\n",
        "obj\t1\t4\t1\tN\t4\tnot-a-gap-type\tno\tna\n",
    ],
)
def test_agp_rejects_invalid_continuity_span_and_gap_type(text: str) -> None:
    with pytest.raises(InputFormatError) as exc_info:
        read_agp(io.StringIO(text))

    assert getattr(exc_info.value, "code", "") != ""


def test_agp_rejects_nonstandard_unknown_gap_and_body_comment() -> None:
    with pytest.raises(InputFormatError) as gap_error:
        read_agp(io.StringIO("obj\t1\t99\t1\tU\t99\tcontig\tno\tna\n"))
    assert gap_error.value.code == "INVALID_AGP_GAP_SPAN"

    with pytest.raises(InputFormatError) as header_error:
        read_agp(io.StringIO("obj\t1\t4\t1\tW\tcontig\t1\t4\t+\n# comment inside body\n"))
    assert header_error.value.code == "INVALID_AGP_HEADER_POSITION"
