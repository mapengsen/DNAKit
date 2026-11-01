"""Tests for dnakit.epigenome.cooler_verify_functions (BINS argument resolution of cooler load/cload)."""

from __future__ import annotations

import inspect
import warnings
from pathlib import Path

import pytest

pytest.importorskip("cooler")
pd = pytest.importorskip("pandas")

from dnakit.epigenome import cooler_verify_functions as cv  # noqa: E402
from dnakit.exceptions import ConfigurationError  # noqa: E402

DATA = Path(__file__).resolve().parents[3] / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
CHROMSIZES = DATA / "genome.fa.sizes"


def _upstream_parse_bins(arg: str) -> pd.DataFrame:
    from cooler.cli._util import parse_bins

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, bins = parse_bins(arg)
    return bins


def _read_bed(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end"], dtype={"chrom": str})


def test_public_functions_take_keyword_only_verbose() -> None:
    parameter = inspect.signature(cv.cooler_resolve_bins).parameters["verbose"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is True


@pytest.mark.skipif(not CHROMSIZES.exists(), reason="shared consistency chromsizes not present")
def test_chromsizes_and_binsize_match_upstream_bins(tmp_path: Path) -> None:
    spec = f"{CHROMSIZES}:10000"
    out = tmp_path / "bins.bed"
    result = cv.cooler_resolve_bins(spec, str(out), verbose=False)
    expected = _upstream_parse_bins(spec)
    written = _read_bed(out)
    assert result["source"] == "chromsizes"
    assert result["binsize"] == 10000
    assert result["bin_count"] == len(expected) == len(written)
    expected = expected.reset_index(drop=True).astype({"chrom": str, "start": "int64", "end": "int64"})
    assert written.equals(expected)
    assert result["chromosomes"] == ["chr1", "chr2"]
    assert result["package"].startswith("cooler ")


@pytest.mark.skipif(not CHROMSIZES.exists(), reason="shared consistency chromsizes not present")
def test_bin_size_that_does_not_divide_chromosome_keeps_partial_last_bin(tmp_path: Path) -> None:
    out = tmp_path / "bins.bed"
    cv.cooler_resolve_bins(f"{CHROMSIZES}:30000", str(out), verbose=False)
    written = _read_bed(out)
    chr2 = written[written["chrom"] == "chr2"]
    assert chr2["end"].max() == 80000
    assert (chr2["end"] - chr2["start"]).max() == 30000
    assert int(chr2.iloc[-1]["end"] - chr2.iloc[-1]["start"]) == 20000


def test_bed_bins_match_upstream_and_are_written(tmp_path: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t100\nchr1\t100\t250\nchr2\t0\t500\n", encoding="utf-8")
    out = tmp_path / "resolved.bed"
    result = cv.cooler_resolve_bins(str(bed), str(out), verbose=False)
    expected = _upstream_parse_bins(str(bed))
    assert result["source"] == "bed" and result["binsize"] is None
    assert result["bin_count"] == 3
    assert result["chromosomes"] == ["chr1", "chr2"]
    assert _read_bed(out).equals(expected.reset_index(drop=True))


def test_missing_chromsizes_file_is_configuration_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="not found"):
        cv.cooler_resolve_bins(f"{tmp_path / 'absent.sizes'}:1000", str(tmp_path / "o.bed"), verbose=False)


@pytest.mark.skipif(not CHROMSIZES.exists(), reason="shared consistency chromsizes not present")
def test_non_integer_binsize_is_configuration_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="integer binsize"):
        cv.cooler_resolve_bins(f"{CHROMSIZES}:1k", str(tmp_path / "o.bed"), verbose=False)


def test_unknown_bins_argument_is_configuration_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="Expected bins"):
        cv.cooler_resolve_bins(str(tmp_path / "missing.bed"), str(tmp_path / "o.bed"), verbose=False)


def test_malformed_bed_is_configuration_error_where_upstream_raises_attribute_error(tmp_path: Path) -> None:
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\t0\n\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        cv.cooler_resolve_bins(str(bad), str(tmp_path / "o.bed"), verbose=False)
    # Upstream defect, reproduced: its error handler names pd.parser, which pandas 2.x does not have.
    from cooler.cli._util import parse_bins

    with pytest.raises(AttributeError, match="parser"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            parse_bins(str(bad))


def test_existing_output_is_refused(tmp_path: Path) -> None:
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t100\n", encoding="utf-8")
    existing = tmp_path / "exists.bed"
    existing.write_text("", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="must not exist"):
        cv.cooler_resolve_bins(str(bed), str(existing), verbose=False)
