"""Paired chunks retain forward coordinates, including a discarded tail."""

from pathlib import Path

import pytest

from dnakit.chunking import ChunkingConfig, export_paired_chunks, iter_paired_chunks
from dnakit.exceptions import InputFormatError
from dnakit.io import iter_rows


def sources(tmp_path: Path) -> tuple[Path, Path]:
    forward = tmp_path / "forward.fasta"
    reverse = tmp_path / "reverse.fasta"
    forward.write_text(">chr1 description\nAACCGGTTAN\n>short\nAC\n")
    reverse.write_text(">chr1 description\nNTAACCGGTT\n>short\nGT\n")
    return forward, reverse


def test_given_reverse_matches_generated_reverse_and_coordinates(tmp_path: Path) -> None:
    forward, reverse = sources(tmp_path)
    config = ChunkingConfig(length=4)
    supplied = list(iter_paired_chunks(forward, reverse_source=reverse, config=config))
    generated = list(iter_paired_chunks(forward, config=config))
    assert supplied == generated
    assert [(r["sequence"], r["source_start"], r["source_end"], r["strand"]) for r in supplied] == [
        ("AACC", 0, 4, "+"),
        ("GGTT", 0, 4, "-"),
        ("GGTT", 4, 8, "+"),
        ("AACC", 4, 8, "-"),
    ]
    assert [r["id"] for r in supplied] == ["chr1:0-4:+", "chr1:0-4:-", "chr1:4-8:+", "chr1:4-8:-"]
    result = export_paired_chunks(
        forward,
        [tmp_path / "out.csv", tmp_path / "out.fasta"],
        reverse_source=reverse,
        config=config,
    )
    assert result == {
        "pair_count": 2,
        "record_count": 4,
        "length": 4,
        "coordinates": "0-based-half-open",
    }
    assert list(iter_rows(tmp_path / "out.csv"))[3]["pair_id"] == "chr1:4-8"


def test_partial_iupac_and_filter(tmp_path: Path) -> None:
    forward, reverse = sources(tmp_path)
    rows = list(
        iter_paired_chunks(
            forward,
            reverse_source=reverse,
            config=ChunkingConfig(length=4, include_partial=True),
            source_ids=["chr1"],
        )
    )
    assert [(r["sequence"], r["source_start"], r["source_end"]) for r in rows[-2:]] == [
        ("AN", 8, 10),
        ("NT", 8, 10),
    ]
    result = export_paired_chunks(forward, tmp_path / "empty.csv", source_ids=["absent"])
    assert result["record_count"] == 0
    assert list(iter_rows(tmp_path / "empty.csv")) == []


@pytest.mark.parametrize("content", [">wrong\nNTAACCGGTT\n", ">chr1\nACGT\n"])
def test_mismatched_sources_fail_without_replacing_output(tmp_path: Path, content: str) -> None:
    forward, reverse = sources(tmp_path)
    reverse.write_text(content)
    output = tmp_path / "chunks.csv"
    output.write_text("original")
    with pytest.raises(InputFormatError, match="matching IDs and lengths"):
        export_paired_chunks(
            forward, output, reverse_source=reverse, config=ChunkingConfig(length=4), overwrite=True
        )
    assert output.read_text() == "original"


def test_record_count_mismatch_is_not_truncated(tmp_path: Path) -> None:
    forward, reverse = sources(tmp_path)
    reverse.write_text(">chr1\nNTAACCGGTT\n")
    with pytest.raises(ValueError, match="shorter"):
        export_paired_chunks(
            forward,
            tmp_path / "chunks.csv",
            reverse_source=reverse,
            config=ChunkingConfig(length=4),
        )
    assert not (tmp_path / "chunks.csv").exists()
