"""Tests for dnakit.epigenome.deeptools_cli_facades3."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import deeptools_cli_facades3 as facade
from dnakit.exceptions import ConfigurationError

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
BAM = DATA / "reads.bam"
BED = DATA / "regions_a.bed"
_ENV_BIN = Path(sys.executable).parent
pytestmark = pytest.mark.skipif(not (_ENV_BIN / "plotEnrichment").is_file() or not BAM.exists(),
                                reason="deepTools or shared BAM missing")


@pytest.fixture(autouse=True)
def _deeptools_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


def test_plot_enrichment_runs(tmp_path: Path) -> None:
    result = facade.plot_enrichment([str(BAM)], [str(BED)], str(tmp_path / "enr.png"), labels=["s1"])
    assert result["status"] == "completed", result["log_tail"]


def test_estimate_read_filtering_runs(tmp_path: Path) -> None:
    result = facade.estimate_read_filtering([str(BAM)], str(tmp_path / "filter.tab"), bin_size=1000)
    assert result["status"] == "completed", result["log_tail"]


def test_labels_must_match(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.plot_enrichment([str(BAM)], [str(BED)], str(tmp_path / "x.png"), labels=["a", "b"])
