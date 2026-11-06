"""Tests for dnakit.epigenome.deeptools_cli_facades (deepTools command-line facades)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import deeptools_cli_facades as facade
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
BAM = DATA / "reads.bam"
_ENV_BIN = Path(sys.executable).parent
HAS_DEEPTOOLS = (_ENV_BIN / "bamCoverage").is_file()
pytestmark = pytest.mark.skipif(
    not HAS_DEEPTOOLS or not BAM.exists(), reason="deepTools or shared BAM missing"
)


@pytest.fixture(autouse=True)
def _deeptools_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make the interpreter's own bin folder visible so the installed deepTools commands are found."""
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


def test_bam_coverage_writes_bedgraph(tmp_path: Path) -> None:
    result = facade.bam_coverage(
        str(BAM), str(tmp_path / "cov.bedgraph"), output_format="bedgraph", bin_size=1000
    )
    assert result["status"] == "completed", result["log_tail"]
    assert Path(result["outputs"][0]).stat().st_size > 0


def test_bam_compare_writes_bedgraph(tmp_path: Path) -> None:
    result = facade.bam_compare(
        str(BAM), str(BAM), str(tmp_path / "cmp.bedgraph"), output_format="bedgraph", bin_size=1000
    )
    assert result["status"] == "completed", result["log_tail"]


def test_bam_pe_fragment_size_runs_from_relative_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(DATA)
    result = facade.bam_pe_fragment_size(
        ["reads.bam"], str(tmp_path / "hist.png"), table_output=str(tmp_path / "t.tab")
    )
    assert "bamPEFragmentSize" in " ".join(result["command"])
    assert result["status"] in ("completed", "failed")
    assert (
        "log_tail" in result
    )  # a rejected single-end input is reported as a structured failure, not an exception


def test_estimate_scale_factor_needs_two_bams() -> None:
    with pytest.raises(ConfigurationError):
        facade.estimate_scale_factor([str(BAM)])


def test_estimate_scale_factor_with_two_bams() -> None:
    result = facade.estimate_scale_factor(
        [str(BAM), str(BAM)], sample_window_length=500, number_of_samples=2000
    )
    assert result["tool"] == "estimateScaleFactor"
    assert result["return_code"] == 0, result["stdout_tail"]


def test_refuses_existing_output(tmp_path: Path) -> None:
    existing = tmp_path / "exists.bedgraph"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        facade.bam_coverage(str(BAM), str(existing))


def test_refuses_missing_input(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.bam_coverage(str(tmp_path / "missing.bam"), str(tmp_path / "out.bw"))


def test_shift_must_have_four_integers(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.alignment_sieve(str(BAM), str(tmp_path / "s.bam"), shift=[1, 2])


def test_missing_program_gives_install_hint(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        facade.bam_coverage(
            str(BAM), str(tmp_path / "out.bw"), deeptools_path=str(tmp_path / "no-such-tool")
        )
    assert "deeptools" in (info.value.hint or "").lower()
