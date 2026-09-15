"""Split progress reports actual file work without changing outputs or quiet mode."""

import gzip
from io import StringIO
from pathlib import Path

import pytest
from rich.console import Console

from dnakit._io_progress import read_observer
from dnakit.datasets import SplitConfig, split
from dnakit.io import read_inputs


@pytest.mark.parametrize("suffix,workers", [("fna", 1), ("fna", 2), ("txt", 1), ("txt.gz", 1)])
def test_split_shows_measured_totals_and_eta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, suffix: str, workers: int
) -> None:
    from dnakit._processing import progress

    captured = StringIO()
    monkeypatch.setattr(
        progress,
        "Console",
        lambda **kwargs: Console(
            file=captured,
            width=180,
            color_system=None,
        ),
    )
    source = tmp_path / f"input.{suffix}"
    text = "".join(
        (f">s{i}\n" if suffix == "fna" else "") + "ACGT" * (i + 1) + "\n" for i in range(10)
    )
    if suffix.endswith(".gz"):
        with gzip.open(source, "wt") as handle:
            handle.write(text)
    else:
        source.write_text(text)
    with read_inputs(source, verbose=False) as records:
        original = list(records)
    config = SplitConfig(ratios={"train": 0.8, "test": 0.2}, seed=42)
    native = split(original, config=config, workers=1, verbose=False)
    result = split(source, config=config, output_path=tmp_path / "split", workers=workers)
    assert result.to_dict() == native.to_dict()
    for name in result.counts:
        assert tuple(result.get(name)) == tuple(native.get(name))
    output = captured.getvalue()
    assert "ETA" in output and "elapsed" in output and "100%" in output
    assert "Saved train.jsonl" in output and "Saved test.jsonl" in output
    assert read_observer.get() is None
    assert not list(tmp_path.rglob("*.progress"))


def test_quiet_split_does_not_render_progress(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    source = tmp_path / "input.txt"
    source.write_text("ACGT\n" * 10)
    split(source, config=SplitConfig(ratios={"train": 0.8, "test": 0.2}), workers=1, verbose=False)
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("boundary", [-1, 0])
def test_index_reports_byte_progress_across_header_boundary(tmp_path: Path, boundary: int) -> None:
    from dnakit.ops._fasta_direction import _index

    source = tmp_path / "input.fna"
    source.write_bytes(b">a\n" + b"A" * (8 * 1024**2 - 3 + boundary) + b"\n>b\nACGT\n")
    events = []
    token = read_observer.set(lambda path, completed: events.append(completed))
    try:
        records, _ = _index(source, 80)
    finally:
        read_observer.reset(token)
    assert [record.identifier for record in records] == ["a", "b"]
    assert events == sorted(events)
    assert events[0] < source.stat().st_size
    assert events[-1] == source.stat().st_size
