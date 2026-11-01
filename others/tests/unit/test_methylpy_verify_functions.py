"""Tests for dnakit.epigenome.methylpy_verify_functions (methylpy filter-allc wrapper)."""

from __future__ import annotations

import gzip
import inspect

import pytest

pytest.importorskip("methylpy.utilities")

from dnakit.epigenome import methylpy_verify_functions as module  # noqa: E402
from dnakit.exceptions import ConfigurationError  # noqa: E402

ROWS = (
    "chr1\t100\t+\tCGG\t3\t5\t1\t4,2\t1,0\n"
    "chr1\t200\t+\tCGA\t0\t4\t0\t3,4\t2,0\n"
    "chr2\t300\t-\tCHH\t1\t6\t1\t5,1\t0,3\n"
)


@pytest.fixture()
def allc(tmp_path):
    path = tmp_path / "in.allc"
    path.write_text(ROWS)
    return path


def _sites(path):
    opener = gzip.open if str(path).endswith("gz") else open
    with opener(path, "rt") as handle:
        return [line.split("\t")[0] + ":" + line.split("\t")[1] for line in handle if line.strip()]


def test_exports_and_keyword_only_verbose():
    assert module.__all__ == ["methylpy_filter_allc_files"]
    params = inspect.signature(module.methylpy_filter_allc_files).parameters
    assert params["verbose"].kind is inspect.Parameter.KEYWORD_ONLY
    assert params["verbose"].default is True
    assert "max_mismatch" not in params


def test_min_cov_keeps_rows_at_or_above_threshold(tmp_path, allc):
    result = module.methylpy_filter_allc_files([allc], [tmp_path / "out.allc"], min_cov=5)
    assert _sites(tmp_path / "out.allc") == ["chr1:100", "chr2:300"]
    assert result["outputs"][0]["records_in"] == 3
    assert result["outputs"][0]["records_out"] == 2


def test_max_cov_is_inclusive(tmp_path, allc):
    module.methylpy_filter_allc_files([allc], [tmp_path / "out.allc"], max_cov=4)
    assert _sites(tmp_path / "out.allc") == ["chr1:200"]


def test_chromosome_and_context_filters(tmp_path, allc):
    module.methylpy_filter_allc_files([allc], [tmp_path / "chr2.allc"], chroms=["chr2"])
    assert _sites(tmp_path / "chr2.allc") == ["chr2:300"]
    module.methylpy_filter_allc_files([allc], [tmp_path / "cg.allc"], mc_type=["CGN"])
    assert _sites(tmp_path / "cg.allc") == ["chr1:100", "chr1:200"]
    module.methylpy_filter_allc_files([allc], [tmp_path / "chh.allc"], mc_type=["CHH"])
    assert _sites(tmp_path / "chh.allc") == ["chr2:300"]
    module.methylpy_filter_allc_files([allc], [tmp_path / "both.allc"], mc_type=["CGN", "CHH"])
    assert _sites(tmp_path / "both.allc") == ["chr1:100", "chr1:200", "chr2:300"]


@pytest.mark.parametrize("bad", ["CGN", ["CG"], ["cgn"], [], ["CGN"] * 17])
def test_refuses_malformed_contexts(tmp_path, allc, bad):
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc"], mc_type=bad)


def test_empty_result_still_writes_the_output(tmp_path, allc):
    result = module.methylpy_filter_allc_files([allc], [tmp_path / "none.allc"], min_cov=100)
    assert (tmp_path / "none.allc").exists()
    assert _sites(tmp_path / "none.allc") == []
    assert result["outputs"][0]["records_out"] == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"min_cov": 5},
        {"mc_type": ["CGN"], "min_cov": 4},
        {"chroms": ["chr1"], "max_cov": 5},
    ],
)
def test_matches_direct_upstream_call(tmp_path, allc, kwargs):
    utilities = pytest.importorskip("methylpy.utilities")
    wrapped = tmp_path / "wrapped.allc"
    direct = tmp_path / "direct.allc"
    module.methylpy_filter_allc_files([allc], [wrapped], **kwargs)
    utilities.filter_allc_files(
        [str(allc)],
        [str(direct)],
        num_procs=1,
        mc_type=kwargs.get("mc_type"),  # list form, as the methylpy command line passes it
        chroms=kwargs.get("chroms"),
        compress_output=False,
        min_cov=kwargs.get("min_cov", 0),
        max_cov=kwargs.get("max_cov"),
    )
    assert wrapped.read_text() == direct.read_text()


def test_compress_output_appends_gz(tmp_path, allc):
    result = module.methylpy_filter_allc_files([allc], [tmp_path / "out.allc"], compress_output=True, min_cov=5)
    assert result["outputs"][0]["output"].endswith("out.allc.gz")
    assert _sites(tmp_path / "out.allc.gz") == ["chr1:100", "chr2:300"]


def test_refuses_existing_output_and_bad_arguments(tmp_path, allc):
    existing = tmp_path / "exists.allc"
    existing.write_text("")
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [existing])
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([tmp_path / "missing.allc"], [tmp_path / "x.allc"])
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [])
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc"], min_cov=True)
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc"], min_cov=5, max_cov=4)
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc"], mc_type="CG;rm")
    with pytest.raises(ConfigurationError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc.gz"])
    with pytest.raises(TypeError):
        module.methylpy_filter_allc_files([allc], [tmp_path / "x.allc"], False)
