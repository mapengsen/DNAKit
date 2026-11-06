from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import afterqc_p15_cli_functions as aq

needs_afterqc = pytest.mark.skipif(
    shutil.which("after.py") is None, reason="AfterQC (after.py) is not installed"
)

INSERT = "ACGTTGCAAGGCTTAACCGGATATCGCGATTTAAACCCGGGTTTACGTACGTGGCCAATTGGAACCT"


def revcomp(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def fastq(prefix: str, sequences: list[str]) -> str:
    return "".join(
        f"@{prefix}{i} 1:N:0:1\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(sequences)
    )


def test_missing_program_has_install_hint(tmp_path: Path) -> None:
    reads = tmp_path / "r_R1.fq"
    reads.write_text(fastq("a", ["ACGT"]))
    other = tmp_path / "r_R2.fq"
    other.write_text(fastq("a", ["ACGT"]))
    with pytest.raises(BackendUnavailableError):
        aq.afterqc_overlap_analysis(
            str(reads),
            str(other),
            str(tmp_path / "out"),
            afterqc_path=str(tmp_path / "missing"),
            verbose=False,
        )


def test_argument_validation(tmp_path: Path) -> None:
    reads = tmp_path / "r_R1.fq"
    reads.write_text(fastq("a", ["ACGT"]))
    with pytest.raises(ConfigurationError):
        aq.afterqc_overlap_analysis(str(reads), str(reads), str(tmp_path / "o1"), verbose=False)
    with pytest.raises(ConfigurationError):
        aq.afterqc_move_barcodes_to_read_names(
            str(reads), str(tmp_path / "o2"), verbose=False
        )  # name lacks 'barcode'
    with pytest.raises(ConfigurationError):
        aq.afterqc_overlap_analysis(
            str(reads), str(tmp_path / "missing.fq"), str(tmp_path / "o3"), verbose=False
        )


@needs_afterqc
def test_overlap_analysis_trims_adapter_and_stores_overlap(tmp_path: Path) -> None:
    adapter = "AGATCGGAAGAGC"
    first = [INSERT + adapter + "TTTT"] * 4
    second = [revcomp(INSERT) + adapter + "TTTT"] * 4
    r1, r2 = tmp_path / "s_R1.fq", tmp_path / "s_R2.fq"
    r1.write_text(fastq("p", first))
    r2.write_text(fastq("p", second))
    result = aq.afterqc_overlap_analysis(
        str(r1),
        str(r2),
        str(tmp_path / "out"),
        store_overlap=True,
        trim_front=0,
        trim_tail=0,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["good_files"] and result["qc_files"] and result["overlap_files"]
    trimmed = next(f for f in result["good_files"] if f.endswith("R1.fq") or "R1" in Path(f).name)
    lines = Path(trimmed).read_text().splitlines()
    assert len(lines[1]) <= len(INSERT) + 2  # the adapter that sticks out of the overlap is removed
