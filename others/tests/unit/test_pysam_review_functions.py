"""Tests for dnakit.alignment.pysam_review_functions (samtools tview text view).

The fake-samtools tests run a small POSIX shell script in place of samtools, so they run without
samtools installed. The real-samtools test is skipped unless samtools is on PATH.
"""

from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from dnakit.alignment.pysam_review_functions import pysam_tview_text
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

posix_only = pytest.mark.skipif(os.name != "posix", reason="fake samtools is a POSIX shell script")

_FAKE = (
    "#!/bin/sh\n"
    'if [ "$1" = "--version" ]; then echo "samtools 9.9-stub"; exit 0; fi\n'
    'echo "ARGS: $*"\n'
    'echo "stub-log" >&2\n'
    "exit 0\n"
)


def _fake_samtools(directory: Path) -> Path:
    script = directory / "fake_samtools"
    script.write_text(_FAKE, encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def _alignment(directory: Path) -> Path:
    path = directory / "reads.bam"
    path.write_bytes(b"not a real bam; only the path is used by the stub")
    return path


def test_missing_samtools_raises_backend_unavailable(tmp_path: Path) -> None:
    bam = _alignment(tmp_path)
    with pytest.raises(BackendUnavailableError) as info:
        pysam_tview_text(
            bam,
            "chr1:10",
            output=tmp_path / "view.txt",
            samtools=tmp_path / "no_such_samtools",
            verbose=False,
        )
    assert "samtools" in str(info.value)
    assert not (tmp_path / "view.txt").exists()


@pytest.mark.parametrize("position", ["chr1", "chr1:0", "chr1:-5", "chr1:1.5", "", ":10"])
def test_invalid_position_is_rejected(tmp_path: Path, position: str) -> None:
    with pytest.raises(ConfigurationError):
        pysam_tview_text(
            _alignment(tmp_path),
            position,
            output=tmp_path / "view.txt",
            samtools=tmp_path / "unused",
            verbose=False,
        )


@pytest.mark.parametrize("width", [0, 19, 1001, True, "80"])
def test_invalid_width_is_rejected(tmp_path: Path, width: object) -> None:
    with pytest.raises(ConfigurationError):
        pysam_tview_text(
            _alignment(tmp_path),
            "chr1:10",
            output=tmp_path / "view.txt",
            width=width,  # type: ignore[arg-type]
            samtools=tmp_path / "unused",
            verbose=False,
        )


def test_existing_output_is_not_overwritten(tmp_path: Path) -> None:
    existing = tmp_path / "view.txt"
    existing.write_text("keep me", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        pysam_tview_text(
            _alignment(tmp_path),
            "chr1:10",
            output=existing,
            samtools=tmp_path / "unused",
            verbose=False,
        )
    assert existing.read_text(encoding="utf-8") == "keep me"


@posix_only
def test_fixed_tview_arguments_and_outputs(tmp_path: Path) -> None:
    bam = _alignment(tmp_path)
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    result = pysam_tview_text(
        bam,
        "chr1:10",
        output=tmp_path / "view.txt",
        width=60,
        reference=reference,
        samtools=_fake_samtools(tmp_path),
        verbose=False,
    )
    assert result["status"] == "ok"
    assert result["backend"] == "samtools"
    assert result["version"] == "samtools 9.9-stub"
    assert result["position"] == "chr1:10"
    assert result["width"] == 60
    assert result["reference"] == str(reference.resolve())
    assert result["lines"] == [
        f"ARGS: tview -d T -p chr1:10 -w 60 {bam.resolve()} {reference.resolve()}"
    ]
    assert result["line_count"] == 1
    assert "stub-log" in result["log"]
    assert (tmp_path / "view.txt").read_text(encoding="utf-8").splitlines() == result["lines"]


@posix_only
def test_reference_is_optional(tmp_path: Path) -> None:
    bam = _alignment(tmp_path)
    result = pysam_tview_text(
        bam,
        "chr2:5",
        output=tmp_path / "view.txt",
        samtools=_fake_samtools(tmp_path),
        verbose=False,
    )
    assert result["reference"] is None
    assert result["lines"] == [f"ARGS: tview -d T -p chr2:5 -w 80 {bam.resolve()}"]


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools is not on PATH")
def test_real_samtools_matches_direct_tview(tmp_path: Path) -> None:
    sam = tmp_path / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n"
        "r1\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"
        "r2\t0\tchr1\t15\t60\t20M\t*\t0\t0\tTACGTACGTACGTACGTACG\tIIIIIIIIIIIIIIIIIIII\n",
        encoding="utf-8",
    )
    bam = tmp_path / "reads.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    direct = subprocess.run(
        ["samtools", "tview", "-d", "T", "-p", "chr1:20", "-w", "40", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = pysam_tview_text(bam, "chr1:20", output=tmp_path / "view.txt", width=40, verbose=False)
    assert result["lines"] == direct.stdout.splitlines()
    assert result["version"] and result["version"].startswith("samtools")
