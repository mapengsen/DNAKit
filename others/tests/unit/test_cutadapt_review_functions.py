"""Tests for dnakit.fasta_fastq.cutadapt_review_functions (cutadapt read-modification facades)."""

from __future__ import annotations

import inspect
import shutil
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.fasta_fastq import cutadapt_review_functions as facade

_LOCAL_CUTADAPT = Path(sys.executable).parent / "cutadapt"
_CUTADAPT = shutil.which("cutadapt") or (
    str(_LOCAL_CUTADAPT) if _LOCAL_CUTADAPT.is_file() else None
)
needs_cutadapt = pytest.mark.skipif(
    _CUTADAPT is None, reason="cutadapt executable is not installed"
)

_READS = [
    ("r1", "ACGTACGTAG", "IIIIIIIIII"),
    ("r2", "NNACGTNN", "IIIIIIII"),
    ("r3", "ACGTGGGGG", "IIIIIIIII"),
]


def _write_fastq(path: Path) -> Path:
    path.write_text(
        "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in _READS), encoding="utf-8"
    )
    return path


def _sequences(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 4)}


@needs_cutadapt
def test_unconditional_cut_removes_bases_from_both_ends(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    target = tmp_path / "cut.fastq"
    result = facade.cutadapt_unconditional_cut(
        str(source), str(target), cut_5p=2, cut_3p=3, verbose=False
    )
    assert result["status"] == "completed"
    assert result["backend"] == "cutadapt"
    assert result["command"][1:3] == ["--report", "minimal"]
    assert "-u" in result["command"] and "-3" in result["command"] and "2" in result["command"]
    # r1 ACGTACGTAG -> GTACGTAG -> GTACG ; r2 NNACGTNN -> ACGTNN -> ACG ; r3 ACGTGGGGG -> GTGGGGG -> GTGG
    assert _sequences(target) == {"r1": "GTACG", "r2": "ACG", "r3": "GTGG"}


@needs_cutadapt
def test_unconditional_cut_five_prime_only(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    target = tmp_path / "cut5.fastq"
    facade.cutadapt_unconditional_cut(str(source), str(target), cut_5p=4, verbose=False)
    assert _sequences(target) == {"r1": "ACGTAG", "r2": "GTNN", "r3": "GGGGG"}


@needs_cutadapt
def test_shorten_keeps_first_or_last_bases(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    first = tmp_path / "first4.fastq"
    last = tmp_path / "last4.fastq"
    facade.cutadapt_shorten(str(source), str(first), length=4, verbose=False)
    facade.cutadapt_shorten(str(source), str(last), length=-4, verbose=False)
    assert _sequences(first) == {"r1": "ACGT", "r2": "NNAC", "r3": "ACGT"}
    assert _sequences(last) == {"r1": "GTAG", "r2": "GTNN", "r3": "GGGG"}


@needs_cutadapt
@pytest.mark.parametrize(
    ("action", "expected_r3"),
    [
        ("trim", "ACGT"),
        ("retain", "ACGTGGGGG"),
        ("mask", "ACGTNNNNN"),
        ("lowercase", "ACGTggggg"),
        ("crop", "GGGGG"),
        ("none", "ACGTGGGGG"),
    ],
)
def test_adapter_action_applies_the_requested_action(
    tmp_path: Path, action: str, expected_r3: str
) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    target = tmp_path / f"{action}.fastq"
    result = facade.cutadapt_adapter_action(
        str(source), str(target), adapter_3p="ggggg", action=action, verbose=False
    )
    assert result["status"] == "completed"
    sequences = _sequences(target)
    assert sequences["r3"] == expected_r3
    assert sequences["r1"] == "ACGTACGTAG"


@needs_cutadapt
def test_adapter_action_reverse_complement_matches_the_other_strand(tmp_path: Path) -> None:
    source = tmp_path / "strand.fastq"
    # The 3' adapter GGGGG is absent from the read but present at the 3' end of its reverse complement (ACGTGGGGG).
    source.write_text("@s1\nCCCCCACGT\n+\nIIIIIIIII\n", encoding="utf-8")
    plain = tmp_path / "plain.fastq"
    both = tmp_path / "both.fastq"
    facade.cutadapt_adapter_action(str(source), str(plain), adapter_3p="GGGGG", verbose=False)
    facade.cutadapt_adapter_action(
        str(source), str(both), adapter_3p="GGGGG", reverse_complement=True, verbose=False
    )
    assert _sequences(plain) == {"s1": "CCCCCACGT"}
    # cutadapt appends " rc" to the name of a read that is written as its reverse complement.
    assert _sequences(both) == {"s1 rc": "ACGT"}


@needs_cutadapt
def test_trim_n_ends_removes_leading_and_trailing_n(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    target = tmp_path / "trimn.fastq"
    result = facade.cutadapt_trim_n_ends(str(source), str(target), verbose=False)
    assert result["status"] == "completed"
    assert "--trim-n" in result["command"]
    assert _sequences(target) == {"r1": "ACGTACGTAG", "r2": "ACGT", "r3": "ACGTGGGGG"}


def test_missing_program_raises_backend_unavailable(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    with pytest.raises(BackendUnavailableError):
        facade.cutadapt_shorten(
            str(source),
            str(tmp_path / "out.fastq"),
            length=4,
            cutadapt_path=str(tmp_path / "missing-cutadapt"),
            verbose=False,
        )


def test_invalid_arguments_are_rejected_before_running(tmp_path: Path) -> None:
    source = _write_fastq(tmp_path / "reads.fastq")
    target = tmp_path / "out.fastq"
    with pytest.raises(ConfigurationError):
        facade.cutadapt_unconditional_cut(str(source), str(target), verbose=False)
    with pytest.raises(ConfigurationError):
        facade.cutadapt_unconditional_cut(str(source), str(target), cut_5p=-1, verbose=False)
    with pytest.raises(ConfigurationError):
        facade.cutadapt_shorten(str(source), str(target), length=0, verbose=False)
    with pytest.raises(ConfigurationError):
        facade.cutadapt_shorten(str(source), str(target), length=True, verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        facade.cutadapt_adapter_action(str(source), str(target), verbose=False)
    with pytest.raises(ConfigurationError):
        facade.cutadapt_adapter_action(
            str(source), str(target), adapter_3p="GGG", action="delete", verbose=False
        )
    with pytest.raises(ConfigurationError):
        facade.cutadapt_adapter_action(
            str(source), str(target), adapter_3p="GGG", times=0, verbose=False
        )
    with pytest.raises(ConfigurationError):
        facade.cutadapt_adapter_action(str(source), str(target), adapter_3p="GXZ", verbose=False)
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim_n_ends(str(tmp_path / "absent.fastq"), str(target), verbose=False)
    existing = tmp_path / "existing.fastq"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        facade.cutadapt_trim_n_ends(str(source), str(existing), verbose=False)
    assert not target.exists()


def test_public_functions_take_keyword_only_verbose() -> None:
    names = facade.__all__
    assert names == [
        "cutadapt_adapter_action",
        "cutadapt_shorten",
        "cutadapt_trim_n_ends",
        "cutadapt_unconditional_cut",
    ]
    for name in names:
        parameter = inspect.signature(getattr(facade, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True
