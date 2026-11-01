"""Tests for dnakit.epigenome.deeptools_gc_bias (computeGCBias / correctGCBias) with a 2bit genome built from the shared data."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import deeptools_gc_bias as facade
from dnakit.exceptions import ConfigurationError

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
BAM = DATA / "reads.bam"
GENOME_FA = DATA / "genome.fa"
_ENV_BIN = Path(sys.executable).parent
pytestmark = pytest.mark.skipif(not (_ENV_BIN / "computeGCBias").is_file() or not BAM.exists(),
                                reason="deepTools or shared BAM missing")


@pytest.fixture(autouse=True)
def _deeptools_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


@pytest.fixture()
def genome_2bit(tmp_path: Path) -> Path:
    from _twobit import write_twobit  # test helper beside this file

    sequences: dict[str, str] = {}
    name, buf = None, []
    for line in GENOME_FA.read_text().splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    sequences[name] = "".join(buf)
    return Path(write_twobit(str(tmp_path / "genome.2bit"), sequences))


def test_compute_then_correct_gc_bias(tmp_path: Path, genome_2bit: Path) -> None:
    freq = tmp_path / "freq.txt"
    computed = facade.compute_gc_bias(str(BAM), str(genome_2bit), 120000, str(freq), fragment_length=150)
    assert computed["status"] == "completed", computed["log_tail"]
    corrected = facade.correct_gc_bias(str(BAM), str(genome_2bit), 120000, str(freq), str(tmp_path / "corr.bam"))
    assert corrected["status"] == "completed", corrected["log_tail"]


def test_effective_genome_size_must_be_positive(tmp_path: Path, genome_2bit: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.compute_gc_bias(str(BAM), str(genome_2bit), 0, str(tmp_path / "f.txt"))
