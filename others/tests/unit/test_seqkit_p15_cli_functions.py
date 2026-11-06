from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import seqkit_p15_cli_functions as sk

needs_seqkit = pytest.mark.skipif(shutil.which("seqkit") is None, reason="seqkit is not installed")

FASTA = ">a first\nACGTACGTAC\n>b second\nGGGGCCCCAA\n>c third\nTTTTTTTTTT\n"
OTHER = ">a other\nACGTACGTAC\n>x only\nCCCCCCCCCC\n>c third\nAAAAAAAAAA\n"


def write(path: Path, text: str) -> str:
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_missing_program_has_install_hint(tmp_path: Path) -> None:
    fasta = write(tmp_path / "x.fa", FASTA)
    with pytest.raises(BackendUnavailableError):
        sk.seqkit_faidx(fasta, seqkit_path=str(tmp_path / "no_such_seqkit"), verbose=False)


def test_argument_validation(tmp_path: Path) -> None:
    fasta = write(tmp_path / "x.fa", FASTA)
    with pytest.raises(ConfigurationError):
        sk.seqkit_faidx(fasta, regions=["a"], verbose=False)  # output_path missing
    with pytest.raises(ConfigurationError):
        sk.seqkit_faidx(fasta, output_path=str(tmp_path / "o.fa"), verbose=False)  # regions missing
    with pytest.raises(ConfigurationError):
        sk.seqkit_common_records([fasta], str(tmp_path / "o.fa"), verbose=False)  # needs two files
    with pytest.raises(ConfigurationError):
        sk.seqkit_common_records(
            [fasta, fasta],
            str(tmp_path / "o.fa"),
            by="id",
            only_positive_strand=True,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        sk.seqkit_stream_directory(
            str(tmp_path), str(tmp_path / "o.fq"), verbose=False
        )  # no time_limit and not find_only
    with pytest.raises(ConfigurationError):
        sk.seqkit_stream_directory(
            str(tmp_path), str(tmp_path / "o.fq"), time_limit="soon", verbose=False
        )
    with pytest.raises(ConfigurationError):
        sk.seqkit_bam_dump_records(
            fasta, str(tmp_path / "d.tsv"), fields=["NotAField"], verbose=False
        )
    existing = write(tmp_path / "exists.fa", FASTA)
    with pytest.raises(ConfigurationError):
        sk.seqkit_sanitize_fastx(fasta, existing, verbose=False)  # outputs are never overwritten


@needs_seqkit
def test_faidx_index_and_regions(tmp_path: Path) -> None:
    fasta = write(tmp_path / "x.fa", FASTA)
    result = sk.seqkit_faidx(fasta, verbose=False)
    assert result["status"] == "completed"
    assert result["sequence_count"] == 3
    assert [row["name"] for row in result["index"]] == ["a", "b", "c"]
    assert result["index"][0]["length"] == 10
    extracted = sk.seqkit_faidx(
        fasta,
        regions=["a:2:4".replace(":2:4", ":2-4"), "b"],
        output_path=str(tmp_path / "r.fa"),
        verbose=False,
    )
    assert extracted["extracted_records"] == 2
    assert (tmp_path / "r.fa").read_text().splitlines()[1] == "CGT"


@needs_seqkit
def test_sanitize_drops_malformed_records(tmp_path: Path) -> None:
    broken = write(
        tmp_path / "b.fq", "@r1\nACGT\n+\nIIII\n@bad\nAC\n+\nIIII\n@r3\nACGTN\n+\nIIIII\nnonsense\n"
    )
    result = sk.seqkit_sanitize_fastx(broken, str(tmp_path / "clean.fq"), verbose=False)
    assert result["passed_records"] == 2
    assert result["discarded_lines"] == 5
    assert result["output_record_count"] == 2


@needs_seqkit
def test_common_records_by_id_name_and_sequence(tmp_path: Path) -> None:
    first, second = write(tmp_path / "a.fa", FASTA), write(tmp_path / "b.fa", OTHER)
    by_id = sk.seqkit_common_records([first, second], str(tmp_path / "id.fa"), verbose=False)
    assert by_id["shared_record_count"] == 2  # a and c
    by_name = sk.seqkit_common_records(
        [first, second], str(tmp_path / "name.fa"), by="name", verbose=False
    )
    assert by_name["shared_record_count"] == 1  # only ">c third"
    by_seq = sk.seqkit_common_records(
        [first, second], str(tmp_path / "seq.fa"), by="sequence", verbose=False
    )
    assert by_seq["shared_record_count"] == 2  # ACGTACGTAC and T10 (the reverse complement of A10)
    assert by_seq["key"] == "sequence"
    positive = sk.seqkit_common_records(
        [first, second],
        str(tmp_path / "pos.fa"),
        by="sequence",
        only_positive_strand=True,
        verbose=False,
    )
    assert positive["shared_record_count"] == 1


@needs_seqkit
def test_replace_header_with_numbering_and_kv(tmp_path: Path) -> None:
    fasta = write(tmp_path / "x.fa", FASTA)
    result = sk.seqkit_replace_header(
        fasta, str(tmp_path / "r.fa"), "^(\\w)", "seq{nr}_$1", number_width=2, verbose=False
    )
    headers = [
        line for line in (tmp_path / "r.fa").read_text().splitlines() if line.startswith(">")
    ]
    assert headers == [">seq01_a first", ">seq02_b second", ">seq03_c third"]
    assert result["output_record_count"] == 3
    kv = write(tmp_path / "kv.tsv", "a\talpha\nb\tbeta\n")
    mapped = sk.seqkit_replace_header(
        fasta,
        str(tmp_path / "m.fa"),
        "^(\\w)",
        "${1}_{kv}",
        key_value_file=kv,
        keep_untouched=True,
        verbose=False,
    )
    names = [line for line in (tmp_path / "m.fa").read_text().splitlines() if line.startswith(">")]
    assert names[0].startswith(">a_alpha") and names[1].startswith(">b_beta")
    assert mapped["status"] == "completed"


@needs_seqkit
def test_stream_directory_find_only(tmp_path: Path) -> None:
    folder = tmp_path / "run"
    folder.mkdir()
    write(folder / "one.fq", "@r1\nACGT\n+\nIIII\n")
    write(folder / "two.fastq", "@r2\nTTTT\n+\nIIII\n")
    result = sk.seqkit_stream_directory(
        str(folder), str(tmp_path / "all.fq"), find_only=True, verbose=False
    )
    assert result["output_record_count"] == 2
    assert result["passed_records"] == 2


def _bam_with_nm(path: Path) -> Path:
    pysam = pytest.importorskip("pysam")
    reference = path.parent / "ref.fa"
    sequence = "ACGTTGCAAGGCTTAACCGGATATCGCGATTTAAACCCGGGTTTACGTACGTGGCCAATT"
    reference.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
    pysam.faidx(str(reference))
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": len(sequence)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for index, (start, edit) in enumerate(((5, None), (10, 3))):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            query = list(sequence[start : start + 20])
            if edit is not None:
                query[edit] = "A" if query[edit] != "A" else "C"
            read.query_sequence = "".join(query)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigarstring = "20M"
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read.set_tag("NM", 0 if edit is None else 1)
            out.write(read)
    return reference


@needs_seqkit
def test_bam_toolbox_dump_accuracy_and_context(tmp_path: Path) -> None:
    bam = tmp_path / "x.bam"
    reference = _bam_with_nm(bam)
    dump = sk.seqkit_bam_dump_records(
        str(bam),
        str(tmp_path / "dump.tsv"),
        fields=["Read", "Ref", "Pos", "Acc", "MapQual"],
        verbose=False,
    )
    assert dump["row_count"] == 2
    assert dump["rows"][0]["Read"] == "read0" and dump["rows"][0]["Acc"] == 100.0
    assert dump["rows"][1]["Acc"] == 95.0
    stats = sk.seqkit_bam_accuracy_stats(
        str(bam), tsv_path=str(tmp_path / "acc.tsv"), verbose=False
    )
    assert stats["accuracy_mean"] == pytest.approx(97.5)
    assert stats["accuracy_weighted_mean"] == pytest.approx(97.5)
    context = sk.seqkit_bam_filter_by_alignment_context(
        str(bam),
        str(reference),
        str(tmp_path / "kept.bam"),
        regex_start="^AAC",
        left_shift=0,
        right_shift=3,
        tsv_path=str(tmp_path / "ctx.tsv"),
        verbose=False,
    )
    assert len(context["context_rows"]) == 2
    assert context["context_rows"][0]["StartSeq"] == "GCA"
    import pysam

    kept = pysam.AlignmentFile(str(tmp_path / "kept.bam"), check_sq=False).count(until_eof=True)
    assert kept == 0  # neither window starts with AAC
