from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import trimmomatic_p15_cli_functions as tm

JAR = os.environ.get("TRIMMOMATIC_JAR")
RUNNABLE = shutil.which("java") is not None and (
    (JAR is not None and Path(JAR).is_file()) or shutil.which("trimmomatic") is not None
)
needs_trimmomatic = pytest.mark.skipif(
    not RUNNABLE,
    reason="Trimmomatic (jar via TRIMMOMATIC_JAR or launcher) and java are not available",
)
JAR_ARGUMENT = JAR if JAR is not None and Path(JAR).is_file() else None


def fastq(names: list[str], sequences: list[str], quality: str = "I") -> str:
    return "".join(
        f"@{n}\n{s}\n+\n{quality * len(s)}\n" for n, s in zip(names, sequences, strict=True)
    )


def test_missing_program_has_install_hint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = tmp_path / "r.fq"
    reads.write_text(fastq(["a"], ["ACGTACGTAC"]))
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.setattr("sys.executable", str(tmp_path / "python"))
    with pytest.raises(BackendUnavailableError):
        tm.trimmomatic_process_reads(str(reads), str(tmp_path / "out"), min_length=3, verbose=False)
    with pytest.raises(BackendUnavailableError):
        tm.trimmomatic_process_reads(
            str(reads),
            str(tmp_path / "out2"),
            min_length=3,
            trimmomatic_jar=str(reads),
            verbose=False,
        )


def test_argument_validation(tmp_path: Path) -> None:
    reads = tmp_path / "r.fq"
    reads.write_text(fastq(["a"], ["ACGTACGTAC"]))
    with pytest.raises(ConfigurationError):
        tm.trimmomatic_process_reads(str(reads), str(tmp_path / "o1"), verbose=False)  # no step
    with pytest.raises(ConfigurationError):
        tm.trimmomatic_process_reads(
            str(reads), str(tmp_path / "o2"), sliding_window=[4], verbose=False
        )
    with pytest.raises(ConfigurationError):
        tm.trimmomatic_process_reads(
            str(reads), str(tmp_path / "o3"), base_count=["N!", 0], verbose=False
        )
    with pytest.raises(ConfigurationError):
        tm.trimmomatic_process_reads(
            str(reads), str(tmp_path / "o4"), max_info=[40, 2.0], verbose=False
        )
    with pytest.raises(ConfigurationError):
        tm.trimmomatic_pair_reads(
            str(reads), str(reads), str(tmp_path / "o5"), trimmomatic_jar="x", verbose=False
        )


@needs_trimmomatic
def test_steps_and_survival_counts(tmp_path: Path) -> None:
    reads = tmp_path / "r.fq"
    long_read = "ACGT" * 10
    reads.write_text(
        fastq(
            ["good", "short", "withN"],
            [long_read, "ACGTACGT", "ACGTNNNNACGTACGTACGTACGTACGTACGTACGTACGT"],
        )
        + f"@lowq\n{long_read}\n+\n"
        + "#" * 40
        + "\n"
    )
    result = tm.trimmomatic_process_reads(
        str(reads),
        str(tmp_path / "out"),
        min_length=20,
        max_length=100,
        average_quality=20,
        base_count=["N", 0, 2],
        trimmomatic_jar=JAR_ARGUMENT,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["steps"] == ["AVGQUAL:20", "BASECOUNT:N:0:2", "MINLEN:20", "MAXLEN:100"]
    assert result["survival"]["input_reads"] == 4
    assert (
        result["survival"]["surviving_reads"] == 1
    )  # short, N-rich and low-quality reads are dropped


@needs_trimmomatic
def test_sliding_window_convert_quality_and_pairs(tmp_path: Path) -> None:
    first, second = tmp_path / "r1.fq", tmp_path / "r2.fq"
    first.write_text(fastq(["p1/1", "p2/1"], ["ACGT" * 8, "TTGCA" * 6]))
    second.write_text(fastq(["p1/2", "p2/2"], ["GGCC" * 8, "ACACA" * 6]))
    result = tm.trimmomatic_process_reads(
        str(first),
        str(tmp_path / "pe"),
        read2=str(second),
        sliding_window=[4, 20],
        convert_quality=64,
        input_phred=33,
        trimmomatic_jar=JAR_ARGUMENT,
        verbose=False,
    )
    assert result["survival"]["input_read_pairs"] == 2
    assert result["survival"]["both_surviving_reads"] == 2
    assert set(result["read_outputs"]) == {
        "R1.paired.fastq.gz",
        "R1.unpaired.fastq.gz",
        "R2.paired.fastq.gz",
        "R2.unpaired.fastq.gz",
    }
    import gzip

    with gzip.open(result["read_outputs"]["R1.paired.fastq.gz"], "rt") as handle:
        quality = handle.read().splitlines()[3]
    assert quality == "h" * 32  # 'I' (Phred 40) re-encoded as Phred+64


@pytest.mark.skipif(
    JAR_ARGUMENT is None or shutil.which("java") is None, reason="needs TRIMMOMATIC_JAR and java"
)
def test_pairomatic_splits_unmatched_mates(tmp_path: Path) -> None:
    first, second = tmp_path / "r1.fq", tmp_path / "r2.fq"
    first.write_text(fastq(["m1/1", "m2/1", "m3/1"], ["ACGTAC", "ACGTAC", "ACGTAC"]))
    second.write_text(fastq(["m1/2", "m3/2"], ["GGGGGG", "TTTTTT"]))
    result = tm.trimmomatic_pair_reads(
        str(first),
        str(second),
        str(tmp_path / "paired"),
        delimiter="/",
        trimmomatic_jar=JAR_ARGUMENT,
        verbose=False,
    )
    counts = result["read_counts"]
    assert counts["R1.paired.fastq"] == 2 and counts["R2.paired.fastq"] == 2
    assert counts["R1.unpaired.fastq"] == 1 and counts["R2.unpaired.fastq"] == 0
