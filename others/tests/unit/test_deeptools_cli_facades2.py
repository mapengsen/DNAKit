"""Tests for dnakit.epigenome.deeptools_cli_facades2 (second batch of deepTools facades)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import deeptools_cli_facades as base
from dnakit.epigenome import deeptools_cli_facades2 as facade
from dnakit.exceptions import ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
BAM = DATA / "reads.bam"
BED = DATA / "regions_a.bed"
_ENV_BIN = Path(sys.executable).parent
pytestmark = pytest.mark.skipif(
    not (_ENV_BIN / "computeMatrix").is_file() or not BAM.exists(),
    reason="deepTools or shared BAM missing",
)


@pytest.fixture(autouse=True)
def _deeptools_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


@pytest.fixture()
def bigwig(tmp_path: Path) -> Path:
    out = tmp_path / "cov.bw"
    result = base.bam_coverage(str(BAM), str(out), output_format="bigwig", bin_size=500)
    assert result["status"] == "completed", result["log_tail"]
    return out


def test_compute_matrix_scale_regions(tmp_path: Path, bigwig: Path) -> None:
    result = facade.compute_matrix_scale_regions(
        [str(bigwig)],
        [str(BED)],
        str(tmp_path / "m.gz"),
        region_body_length=500,
        before_region_length=200,
        after_region_length=200,
        bin_size=100,
    )
    assert result["status"] == "completed", result["log_tail"]


def test_compute_matrix_reference_point_then_heatmap(tmp_path: Path, bigwig: Path) -> None:
    matrix = tmp_path / "ref.gz"
    result = facade.compute_matrix_reference_point(
        [str(bigwig)],
        [str(BED)],
        str(matrix),
        reference_point="center",
        before_region_length=200,
        after_region_length=200,
        bin_size=100,
    )
    assert result["status"] == "completed", result["log_tail"]
    heat = facade.plot_heatmap(str(matrix), str(tmp_path / "heat.png"))
    assert heat["status"] == "completed", heat["log_tail"]


def test_compute_matrix_profile_plot(tmp_path: Path, bigwig: Path) -> None:
    matrix = tmp_path / "prof.gz"
    facade.compute_matrix_scale_regions(
        [str(bigwig)],
        [str(BED)],
        str(matrix),
        region_body_length=500,
        before_region_length=200,
        after_region_length=200,
        bin_size=100,
    )
    result = facade.plot_profile(str(matrix), str(tmp_path / "prof.png"))
    assert result["status"] == "completed", result["log_tail"]


def test_multi_bam_summary_and_correlation(tmp_path: Path) -> None:
    npz = tmp_path / "summary.npz"
    result = facade.multi_bam_summary(
        [str(BAM), str(BAM)], str(npz), bin_size=1000, labels=["a", "b"]
    )
    assert result["status"] == "completed", result["log_tail"]
    corr = facade.plot_correlation(
        str(npz), str(tmp_path / "corr.png"), method="spearman", plot_type="heatmap"
    )
    assert corr["status"] == "completed", corr["log_tail"]


def test_multi_bigwig_summary_and_bigwig_average(tmp_path: Path, bigwig: Path) -> None:
    npz = tmp_path / "bw.npz"
    result = facade.multi_bigwig_summary([str(bigwig), str(bigwig)], str(npz), bin_size=1000)
    assert result["status"] == "completed", result["log_tail"]
    avg = facade.bigwig_average([str(bigwig), str(bigwig)], str(tmp_path / "avg.bw"), bin_size=500)
    assert avg["status"] == "completed", avg["log_tail"]


def test_pca_plot(tmp_path: Path) -> None:
    npz = tmp_path / "pca.npz"
    facade.multi_bam_summary([str(BAM), str(BAM)], str(npz), bin_size=1000)
    result = facade.plot_pca(str(npz), str(tmp_path / "pca.png"))
    assert result["status"] in ("completed", "failed")
    assert result["command"][0] == "plotPCA"


def test_labels_must_match_number_of_files(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.multi_bam_summary([str(BAM), str(BAM)], str(tmp_path / "x.npz"), labels=["only-one"])


def test_kmeans_must_be_positive(tmp_path: Path, bigwig: Path) -> None:
    matrix = tmp_path / "k.gz"
    built = facade.compute_matrix_scale_regions(
        [str(bigwig)],
        [str(BED)],
        str(matrix),
        region_body_length=500,
        before_region_length=200,
        after_region_length=200,
        bin_size=100,
    )
    assert built["status"] == "completed", built["log_tail"]
    with pytest.raises(ConfigurationError, match="kmeans"):
        facade.plot_heatmap(str(matrix), str(tmp_path / "h.png"), kmeans=0)
