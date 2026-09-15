"""File-region reverse complements agree byte-for-byte with the general codec."""

import random
from pathlib import Path

import pytest

from dnakit.exceptions import InputFormatError
from dnakit.io import ReadConfig, WriteConfig
from dnakit.ops import _fasta_direction as fast
from dnakit.ops import reverse_complement


@pytest.fixture(autouse=True)
def small_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fast, "_MIN_BYTES", 0)
    monkeypatch.setattr(fast, "_BLOCK_BYTES", 67)


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("newline", [b"\n", b"\r\n"])
@pytest.mark.parametrize("final_newline", [True, False])
def test_regions_match_codec(
    tmp_path: Path, workers: int, newline: bytes, final_newline: bool
) -> None:
    source = tmp_path / "in.fna"
    records = [
        ("chr1 description", b"aCGTryswkmbdhvn" * 23, 19),
        ("empty", b"", 9),
        ("chr2", b"TGCAN" * 31, 11),
    ]
    payload = b""
    for header, sequence, width in records:
        payload += b">" + header.encode() + newline
        payload += b"".join(
            sequence[i : i + width] + newline for i in range(0, len(sequence), width)
        )
    if not final_newline:
        payload = payload[: -len(newline)]
    source.write_bytes(payload)
    writer = WriteConfig(line_width=17)
    reference = tmp_path / "reference.fna"
    reverse_complement(
        source, reference, read_config=ReadConfig(), write_config=writer, workers=1, verbose=False
    )
    result = fast.try_reverse_fasta(
        source,
        tmp_path / "actual.fna",
        read_config=None,
        write_config=writer,
        sequence_column=None,
        workers=workers,
        show_progress=False,
    )
    assert result is not None and result.record_count == 3
    assert (tmp_path / "actual.fna").read_bytes() == reference.read_bytes()
    assert result.byte_count == reference.stat().st_size
    assert source.read_bytes() == payload
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize(
    "body",
    [b"AAGC", b"AAAA\nAA\nACGT\nAAAA", b"AA\n\nAAGCT\nA", b"AAAA\r\nCCCC\n", b"AAAA\nCCCC\r\n"],
)
def test_unwrapped_and_irregular_layouts(tmp_path: Path, body: bytes) -> None:
    source = tmp_path / "input.fna"
    source.write_bytes(b">chr\n" + body)
    reverse_complement(source, tmp_path / "actual.fna", workers=1, verbose=False)
    reverse_complement(
        source, tmp_path / "reference.fna", read_config=ReadConfig(), workers=1, verbose=False
    )
    assert (tmp_path / "actual.fna").read_bytes() == (tmp_path / "reference.fna").read_bytes()


def test_failure_preserves_existing_output(tmp_path: Path) -> None:
    source = tmp_path / "input.fna"
    source.write_bytes(b">valid\n" + b"ACGT\n" * 100 + b">invalid\nACGX\n")
    target = tmp_path / "out.fna"
    target.write_bytes(b"previous output")
    with pytest.raises(InputFormatError):
        reverse_complement(
            source, target, write_config=WriteConfig(overwrite=True), workers=2, verbose=False
        )
    assert target.read_bytes() == b"previous output"
    assert not list(tmp_path.glob("*.tmp"))


def test_size_limit_before_commit(tmp_path: Path) -> None:
    source = tmp_path / "input.fna"
    source.write_bytes(b">chr\nACGT\nACGT\n")
    with pytest.raises(InputFormatError, match="max_output_bytes"):
        reverse_complement(
            source,
            tmp_path / "out.fna",
            write_config=WriteConfig(max_output_bytes=5),
            workers=1,
            verbose=False,
        )
    assert not (tmp_path / "out.fna").exists()


@pytest.mark.parametrize("seed", range(20))
def test_mixed_wrapping_never_loses_bases(tmp_path: Path, seed: int) -> None:
    rng = random.Random(seed)
    lines = []
    for index in range(3):
        lines.append(f">chr{index} description\n".encode())
        sequence = bytes(rng.choices(b"ACGTryswkmbdhvn", k=rng.randrange(1, 800)))
        start = 0
        while start < len(sequence):
            width = rng.randrange(1, 60)
            lines.append(sequence[start : start + width] + rng.choice([b"\n", b"\r\n", b"\r"]))
            start += width
    source = tmp_path / "mixed.fna"
    source.write_bytes(b"".join(lines))
    reverse_complement(source, tmp_path / "actual.fna", workers=1, verbose=False)
    reverse_complement(
        source, tmp_path / "reference.fna", read_config=ReadConfig(), workers=1, verbose=False
    )
    assert (tmp_path / "actual.fna").read_bytes() == (tmp_path / "reference.fna").read_bytes()
