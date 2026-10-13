"""Contract tests for G1b wrappers: pyfaidx, dnaio, screed, and the pyfastx/py2bit layers."""

from __future__ import annotations

import inspect
import io
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import (
    dnaio_functions as dnaio_wrappers,
)
from dnakit.sequence import (
    py2bit_functions,
    pyfastx_functions,
    reverse_complement_dna,
)
from dnakit.sequence import (
    pyfaidx_functions as pyfaidx_wrappers,
)
from dnakit.sequence import (
    screed_functions as screed_wrappers,
)
from dnakit.sequence.backends import build_screed_database

FASTA_TEXT = ">r1 first read\nACGTNACGT\n>r2\nggcc\n"
FASTQ_TEXT = "@q1 desc\nACGTAC\n+\nIIII#!\n@q2\nGGCCAA\n+\nABCDEF\n"
MATE1_TEXT = "@p/1 x\nACGT\n+\nIIII\n@q/1\nGG\n+\nII\n"
MATE2_TEXT = "@p/2 y\nTTTT\n+\nIIII\n@q/2\nCC\n+\nII\n"


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="ascii")
    return path


def test_wrapped_functions_take_keyword_only_verbose_defaulting_to_true() -> None:
    for module in (pyfaidx_wrappers, dnaio_wrappers, screed_wrappers):
        assert module.__all__, module.__name__
        for name in module.__all__:
            parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
            assert parameter.default is True, name


# pyfaidx -------------------------------------------------------------------------------------


def test_translate_chromosome_names_maps_explicit_pairs_and_reports_unmapped() -> None:
    pyfaidx = pytest.importorskip("pyfaidx")
    result = pyfaidx_wrappers.translate_chromosome_names(
        ["chr1", "chrM", "chr2", "chrX"],
        from_names=["chr1", "chr2", "chrM"],
        to_names=["1", "2", "MT"],
        verbose=False,
    )
    assert result["translations"] == [
        {"name": "chr1", "translated": "1"},
        {"name": "chrM", "translated": "MT"},
        {"name": "chr2", "translated": "2"},
    ]
    assert result["unmapped"] == ["chrX"]
    assert (result["requested_count"], result["translated_count"]) == (4, 3)
    upstream = pyfaidx.translate_chr_name(["chr1", "chr2", "chrM"], ["1", "2", "MT"])
    assert upstream("chr2") == "2"


def test_translate_chromosome_names_rejects_repeated_or_mismatched_pairs() -> None:
    pytest.importorskip("pyfaidx")
    with pytest.raises(ConfigurationError, match="repeat"):
        pyfaidx_wrappers.translate_chromosome_names(
            ["chr1"], from_names=["chr1", "chr1"], to_names=["1", "2"], verbose=False
        )
    with pytest.raises(ConfigurationError, match="same length"):
        pyfaidx_wrappers.translate_chromosome_names(
            ["chr1"], from_names=["chr1"], to_names=["1", "2"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        pyfaidx_wrappers.translate_chromosome_names(
            "chr1", from_names=["chr1"], to_names=["1"], verbose=False
        )


def test_split_bed_entry_matches_pyfaidx_and_hand_length() -> None:
    pyfaidx = pytest.importorskip("pyfaidx")
    result = pyfaidx_wrappers.split_bed_entry("chr1\t10\t20\tgeneA\n", verbose=False)
    assert result["chromosome"] == "chr1"
    assert (result["start"], result["end"], result["length"]) == (10, 20, 10)
    assert result["is_comment"] is False
    assert pyfaidx.bed_split("chr1\t10\t20\tgeneA\n") == ("chr1", 10, 20)

    comment = pyfaidx_wrappers.split_bed_entry("# header", verbose=False)
    assert comment["is_comment"] is True
    assert comment["chromosome"] is None


def test_split_bed_entry_rejects_malformed_lines() -> None:
    pytest.importorskip("pyfaidx")
    for bad in ["chr1 10", "chr1 20 10", "chr1 -5 10", "   ", "track name=x"]:
        with pytest.raises(ConfigurationError):
            pyfaidx_wrappers.split_bed_entry(bad, verbose=False)


def test_split_ucsc_region_converts_one_based_inclusive_to_half_open() -> None:
    pyfaidx = pytest.importorskip("pyfaidx")
    result = pyfaidx_wrappers.split_ucsc_region("chr1:101-200", verbose=False)
    assert result["chromosome"] == "chr1"
    assert (result["start"], result["end"], result["length"]) == (100, 200, 100)
    assert result["is_whole_chromosome"] is False
    assert pyfaidx.ucsc_split("chr1:101-200") == ("chr1", 100, 200)

    whole = pyfaidx_wrappers.split_ucsc_region("chr7", verbose=False)
    assert whole["is_whole_chromosome"] is True
    assert whole["start"] is None and whole["end"] is None


def test_split_ucsc_region_rejects_bad_coordinates() -> None:
    pytest.importorskip("pyfaidx")
    for bad in ["chr1:0-5", "chr1:1,000-2000", "chr1:10-5", "chr1:x", "chr1:5"]:
        with pytest.raises(ConfigurationError):
            pyfaidx_wrappers.split_ucsc_region(bad, verbose=False)


# dnaio ---------------------------------------------------------------------------------------


def test_read_fastx_records_matches_direct_dnaio_for_fastq(tmp_path: Path) -> None:
    dnaio = pytest.importorskip("dnaio")
    source = _write(tmp_path, "reads.fq", FASTQ_TEXT)
    result = dnaio_wrappers.read_fastx_records(source, verbose=False)
    assert result["format"] == "fastq"
    assert result["entry_unit"] == "record"
    assert result["entry_count"] == 2
    assert result["records"][0] == {
        "name": "q1 desc",
        "id": "q1",
        "comment": "desc",
        "sequence": "ACGTAC",
        "length": 6,
        "qualities": "IIII#!",
    }
    with dnaio.open(str(source)) as reader:
        direct = [(read.name, read.sequence, read.qualities) for read in reader]
    assert direct == [(e["name"], e["sequence"], e["qualities"]) for e in result["records"]]


def test_read_fastx_records_reads_fasta_without_qualities(tmp_path: Path) -> None:
    pytest.importorskip("dnaio")
    source = _write(tmp_path, "reads.fa", FASTA_TEXT)
    result = dnaio_wrappers.read_fastx_records(source, verbose=False)
    assert result["format"] == "fasta"
    assert result["records"][1] == {
        "name": "r2",
        "id": "r2",
        "comment": "",
        "sequence": "ggcc",
        "length": 4,
    }


def test_read_fastx_records_pairs_mates_and_interleaved_input(tmp_path: Path) -> None:
    pytest.importorskip("dnaio")
    first = _write(tmp_path, "r1.fq", MATE1_TEXT)
    second = _write(tmp_path, "r2.fq", MATE2_TEXT)
    paired = dnaio_wrappers.read_fastx_records(first, mate_source=second, verbose=False)
    assert paired["layout"] == "paired-end"
    assert paired["entry_unit"] == "pair"
    assert paired["entry_count"] == 2
    assert paired["records"][0]["read1"]["name"] == "p/1 x"
    assert paired["records"][0]["read2"]["name"] == "p/2 y"
    assert paired["records"][1]["read2"]["sequence"] == "CC"

    interleaved_path = _write(
        tmp_path, "interleaved.fq", "@p/1\nACGT\n+\nIIII\n@p/2\nTTTT\n+\nIIII\n"
    )
    interleaved = dnaio_wrappers.read_fastx_records(
        interleaved_path, interleaved=True, verbose=False
    )
    assert interleaved["layout"] == "interleaved"
    assert interleaved["entry_count"] == 1
    assert interleaved["records"][0]["read1"]["sequence"] == "ACGT"
    assert interleaved["records"][0]["read2"]["sequence"] == "TTTT"


def test_read_fastx_records_rejects_mismatched_mates_and_bad_limits(tmp_path: Path) -> None:
    pytest.importorskip("dnaio")
    first = _write(tmp_path, "r1.fq", MATE1_TEXT)
    wrong = _write(tmp_path, "wrong.fq", "@z/2\nTTTT\n+\nIIII\n@q/2\nCC\n+\nII\n")
    with pytest.raises(ConfigurationError):
        dnaio_wrappers.read_fastx_records(first, mate_source=wrong, verbose=False)
    with pytest.raises(ConfigurationError, match="cannot be combined"):
        dnaio_wrappers.read_fastx_records(first, mate_source=first, interleaved=True, verbose=False)
    source = _write(tmp_path, "reads.fq", FASTQ_TEXT)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        dnaio_wrappers.read_fastx_records(source, max_records=1, verbose=False)


def test_check_read_mates_matches_dnaio_records_are_mates() -> None:
    dnaio = pytest.importorskip("dnaio")
    assert dnaio_wrappers.check_read_mates(["p/1", "p/2"], verbose=False)["mates"] is True
    assert (
        dnaio_wrappers.check_read_mates(["read 1:N:0", "read 2:N:0"], verbose=False)["mates"]
        is True
    )
    unrelated = dnaio_wrappers.check_read_mates(["p/1", "q/2"], verbose=False)
    assert unrelated["mates"] is False
    direct = dnaio.records_are_mates(
        dnaio.SequenceRecord("p/1", "", None), dnaio.SequenceRecord("q/2", "", None)
    )
    assert direct is False


def test_check_read_mates_rejects_single_name_and_bad_text() -> None:
    pytest.importorskip("dnaio")
    with pytest.raises(ConfigurationError, match="at least two"):
        dnaio_wrappers.check_read_mates(["p/1"], verbose=False)
    with pytest.raises(ConfigurationError):
        dnaio_wrappers.check_read_mates(["p/1", ""], verbose=False)
    with pytest.raises(ConfigurationError, match="line breaks"):
        dnaio_wrappers.check_read_mates(["p/1", "p\n/2"], verbose=False)


# screed --------------------------------------------------------------------------------------


def test_export_screed_database_fasta_matches_upstream_tofasta(tmp_path: Path) -> None:
    screed = pytest.importorskip("screed")
    source = _write(tmp_path, "reads.fa", FASTA_TEXT)
    built = build_screed_database(source, database=tmp_path / "reads", verbose=False)
    output = tmp_path / "exported.fa"
    result = screed_wrappers.export_screed_database(built["database"], output, verbose=False)
    assert result["record_count"] == 2
    assert result["output_format"] == "fasta"
    # screed writes a trailing space after a header that has no description.
    assert output.read_text(encoding="utf-8") == ">r1 first read\nACGTNACGT\n>r2 \nggcc\n"

    upstream = tmp_path / "upstream.fa"
    screed.ToFasta(built["database"], str(upstream))
    assert output.read_text(encoding="utf-8") == upstream.read_text(encoding="utf-8")


def test_export_screed_database_fastq_keeps_quality_values(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write(tmp_path, "reads.fq", FASTQ_TEXT)
    built = build_screed_database(source, database=tmp_path / "reads_fq", verbose=False)
    output = tmp_path / "exported.fq"
    result = screed_wrappers.export_screed_database(
        built["database"], output, output_format="fastq", verbose=False
    )
    assert result["record_count"] == 2
    assert output.read_text(encoding="utf-8") == (
        "@q1 desc\nACGTAC\n+\nIIII#!\n@q2 \nGGCCAA\n+\nABCDEF\n"
    )


def test_export_screed_database_refuses_existing_output_and_extra_records(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write(tmp_path, "reads.fa", FASTA_TEXT)
    built = build_screed_database(source, database=tmp_path / "reads", verbose=False)
    existing = _write(tmp_path, "taken.fa", "keep\n")
    with pytest.raises(ConfigurationError, match="already exists"):
        screed_wrappers.export_screed_database(built["database"], existing, verbose=False)
    assert existing.read_text(encoding="ascii") == "keep\n"
    with pytest.raises(ConfigurationError, match="max_records=1"):
        screed_wrappers.export_screed_database(
            built["database"], tmp_path / "new.fa", max_records=1, verbose=False
        )
    with pytest.raises(ConfigurationError, match="does not exist"):
        screed_wrappers.export_screed_database(
            tmp_path / "missing", tmp_path / "new2.fa", verbose=False
        )


def test_read_screed_fastx_records_matches_upstream_fasta_iterator(tmp_path: Path) -> None:
    screed = pytest.importorskip("screed")
    source = _write(tmp_path, "reads.fa", FASTA_TEXT)
    result = screed_wrappers.read_screed_fastx_records(source, verbose=False)
    assert result["format"] == "fasta"
    assert result["records"] == [
        {"name": "r1", "description": "first read", "length": 9, "sequence": "ACGTNACGT"},
        {"name": "r2", "description": "", "length": 4, "sequence": "ggcc"},
    ]
    with source.open("rb") as handle:
        upstream = [
            (record.name, record.sequence)
            for record in screed.fasta.fasta_iter(handle, parse_description=True)
        ]
    assert upstream == [(e["name"], e["sequence"]) for e in result["records"]]


def test_read_screed_fastx_records_reports_fastq_quality_and_annotations(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    source = _write(tmp_path, "reads.fq", FASTQ_TEXT)
    result = screed_wrappers.read_screed_fastx_records(source, file_format="fastq", verbose=False)
    assert result["records"] == [
        {
            "name": "q1",
            "description": "desc",
            "length": 6,
            "sequence": "ACGTAC",
            "quality": "IIII#!",
        },
        {"name": "q2", "description": "", "length": 6, "sequence": "GGCCAA", "quality": "ABCDEF"},
    ]


def test_read_screed_fastx_records_rejects_wrong_format_and_limits(tmp_path: Path) -> None:
    pytest.importorskip("screed")
    fasta = _write(tmp_path, "reads.fa", FASTA_TEXT)
    fastq = _write(tmp_path, "reads.fq", FASTQ_TEXT)
    with pytest.raises(ConfigurationError, match="as fastq"):
        screed_wrappers.read_screed_fastx_records(fasta, file_format="fastq", verbose=False)
    with pytest.raises(ConfigurationError, match="as fasta"):
        screed_wrappers.read_screed_fastx_records(fastq, verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        screed_wrappers.read_screed_fastx_records(fasta, max_records=1, verbose=False)


def test_format_fastx_records_writes_hand_checked_text() -> None:
    pytest.importorskip("screed")
    fasta = screed_wrappers.format_fastx_records(
        [
            {"name": "r1", "description": "first read", "sequence": "ACGTN"},
            {"name": "r2", "sequence": "ggcc"},
        ],
        verbose=False,
    )
    assert fasta["text"] == ">r1 first read\nACGTN\n>r2\nggcc\n"
    assert fasta["record_count"] == 2

    fastq = screed_wrappers.format_fastx_records(
        [{"name": "q1", "description": "desc", "sequence": "ACGTAC", "quality": "IIII#!"}],
        output_format="fastq",
        verbose=False,
    )
    assert fastq["text"] == "@q1 desc\nACGTAC\n+\nIIII#!\n"


def test_format_fastx_records_matches_direct_write_fastx() -> None:
    screed = pytest.importorskip("screed")
    buffer = io.BytesIO()
    screed.screedRecord.write_fastx(
        screed.Record(name="r1", description="first read", sequence="ACGTN"), buffer
    )
    result = screed_wrappers.format_fastx_records(
        [{"name": "r1", "description": "first read", "sequence": "ACGTN"}], verbose=False
    )
    assert result["text"] == buffer.getvalue().decode("utf-8")


def test_format_fastx_records_rejects_invalid_records() -> None:
    pytest.importorskip("screed")
    cases = [
        ([], "fasta", "non-empty"),
        ([{"name": "q", "sequence": "AC", "quality": "I"}], "fastq", "one symbol per base"),
        ([{"name": "q", "sequence": "AC", "quality": "II"}], "fasta", "only when"),
        ([{"name": "r 1", "sequence": "AC"}], "fasta", "whitespace"),
        ([{"name": "r1", "sequence": "AC1"}], "fasta", "ASCII letters"),
        (["AC"], "fasta", "mapping"),
    ]
    for records, output_format, message in cases:
        with pytest.raises(ConfigurationError, match=message):
            screed_wrappers.format_fastx_records(
                records, output_format=output_format, verbose=False
            )


# pyfastx and py2bit -------------------------------------------------------------------------


def test_pyfastx_reverse_complement_agrees_with_dnakit_duplicate() -> None:
    pyfastx = pytest.importorskip("pyfastx")
    assert pyfastx.reverse_complement("ACGTN") == reverse_complement_dna("ACGTN") == "NACGT"


def test_pyfastx_and_py2bit_layers_export_no_new_wrappers() -> None:
    pytest.importorskip("py2bit")
    assert pyfastx_functions.__all__ == ()
    assert py2bit_functions.__all__ == ()
