"""Tests for dnakit.variants.bed_reader_gap_functions (bed-reader cloud reads, cloud checks, sample
downloads).

Cloud reads and checks use file:// URLs of local files, so bed-reader's native cloud readers run
without any
network access. The pooch-backed sample wrappers are tested with fakes; nothing is downloaded.
"""

from __future__ import annotations

import inspect
import socket
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import bed_reader_gap_functions as bg

# bed-reader's documented sample: -127 is the int8 missing value
INT8_VALUES = [[1, 0, -127, 0], [2, 0, -127, 2], [0, 1, 2, 0]]
FLOAT_VALUES = [[1.0, 0.0, None, 0.0], [2.0, 0.0, None, 2.0], [0.0, 1.0, 2.0, 0.0]]


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    monkeypatch.setattr(socket, "socket", refuse)


@pytest.fixture()
def bed_path(tmp_path: Path) -> Path:
    bed_reader = pytest.importorskip("bed_reader")
    np = pytest.importorskip("numpy")
    path = tmp_path / "small.bed"
    bed_reader.to_bed(str(path), np.array(INT8_VALUES, dtype=np.int8))
    return path


@pytest.fixture()
def bed_uri(bed_path: Path) -> str:
    return bed_path.as_uri()


# ---- cloud genotype reads (file:// URLs run through the native cloud readers) -----------------


@pytest.mark.parametrize("dtype", ["int8", "float32", "float64"])
def test_read_cloud_genotypes_each_dtype(bed_uri: str, dtype: str) -> None:
    pytest.importorskip("bed_reader")
    result = bg.read_cloud_genotypes(bed_uri, dtype=dtype, verbose=False)
    assert result["status"] == "completed"
    assert result["operation"] == "read_cloud_genotypes"
    assert result["dtype"] == dtype
    assert result["shape"] == [3, 4]
    assert result["missing_calls"] == 2
    expected = [[1, 0, None, 0], [2, 0, None, 2], [0, 1, 2, 0]] if dtype == "int8" else FLOAT_VALUES
    assert result["values"] == expected


def test_read_cloud_genotypes_matches_local_open_bed(bed_path: Path) -> None:
    bed_reader = pytest.importorskip("bed_reader")
    np = pytest.importorskip("numpy")
    local = bed_reader.open_bed(str(bed_path)).read(dtype="float32")
    cloud = bg.read_cloud_genotypes(bed_path.as_uri(), dtype="float32", verbose=False)
    expected = [[None if np.isnan(value) else float(value) for value in row] for row in local]
    assert cloud["values"] == expected


def test_read_cloud_genotypes_index_subset(bed_uri: str) -> None:
    pytest.importorskip("bed_reader")
    one_sample = bg.read_cloud_genotypes(bed_uri, iid_index=[2], sid_index=[1, 3], verbose=False)
    assert one_sample["values"] == [[1.0, 0.0]]
    assert one_sample["shape"] == [1, 2]
    assert one_sample["missing_calls"] == 0
    missing_column = bg.read_cloud_genotypes(bed_uri, sid_index=[2], verbose=False)
    assert missing_column["values"] == [[None], [None], [2.0]]
    assert missing_column["missing_calls"] == 2


@pytest.mark.parametrize(
    ("location", "kwargs", "error"),
    [
        ("plain/path/small.bed", {}, "scheme"),
        ("https://", {}, "scheme"),
        ("", {}, "URL string"),
        ("file:///a\nb.bed", {}, "control characters"),
        ("file:///x/small.bed", {"dtype": "float16"}, "dtype"),
        ("file:///x/small.bed", {"cloud_options": {"timeout": 10}}, "cloud_options"),
        ("file:///x/small.bed", {"max_cells": 0}, "max_cells"),
    ],
)
def test_read_cloud_genotypes_rejects_bad_arguments(
    location: str, kwargs: dict[str, Any], error: str
) -> None:
    with pytest.raises(ConfigurationError, match=error):
        bg.read_cloud_genotypes(location, verbose=False, **kwargs)


def test_read_cloud_genotypes_rejects_bad_selection_and_size(bed_uri: str) -> None:
    pytest.importorskip("bed_reader")
    with pytest.raises(ConfigurationError, match="iid_index"):
        bg.read_cloud_genotypes(bed_uri, iid_index=[3], verbose=False)
    with pytest.raises(ConfigurationError, match="sid_index"):
        bg.read_cloud_genotypes(bed_uri, sid_index=[], verbose=False)
    with pytest.raises(ConfigurationError, match="sid_index"):
        bg.read_cloud_genotypes(bed_uri, sid_index=[True], verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells"):
        bg.read_cloud_genotypes(bed_uri, max_cells=3, verbose=False)


def test_read_cloud_genotypes_missing_file_is_reported(tmp_path: Path) -> None:
    pytest.importorskip("bed_reader")
    missing = (tmp_path / "absent.bed").as_uri()
    with pytest.raises(ConfigurationError, match="could not be opened"):
        bg.read_cloud_genotypes(missing, verbose=False)


def test_read_cloud_genotypes_missing_dependency_is_reported(
    monkeypatch: pytest.MonkeyPatch, bed_uri: str
) -> None:
    monkeypatch.setitem(sys.modules, "bed_reader", None)
    with pytest.raises(BackendUnavailableError, match="bio-bed-reader"):
        bg.read_cloud_genotypes(bed_uri, verbose=False)


# ---- cloud location check -----------------------------------------------------------------------


def test_check_cloud_location_true_for_readable_file_and_false_for_missing(
    bed_uri: str, tmp_path: Path
) -> None:
    pytest.importorskip("bed_reader")
    assert bg.check_cloud_location(bed_uri, verbose=False) is True
    assert bg.check_cloud_location((tmp_path / "absent.bed").as_uri(), verbose=False) is False


def test_check_cloud_location_passes_options_and_maps_value_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    native = pytest.importorskip("bed_reader.bed_reader")
    calls: list[tuple[str, dict[str, str]]] = []

    def accept(location: str, options: dict[str, str]) -> None:
        calls.append((location, options))

    monkeypatch.setattr(native, "check_file_cloud", accept)
    assert (
        bg.check_cloud_location(
            "https://example.org/x.bed", cloud_options={"timeout": "10s"}, verbose=False
        )
        is True
    )
    assert calls == [("https://example.org/x.bed", {"timeout": "10s"})]

    def refuse(location: str, options: dict[str, str]) -> None:
        raise ValueError("Object store error: not found")

    monkeypatch.setattr(native, "check_file_cloud", refuse)
    assert bg.check_cloud_location("https://example.org/x.bed", verbose=False) is False


def test_check_cloud_location_rejects_plain_path_and_reports_missing_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(ConfigurationError, match="scheme"):
        bg.check_cloud_location("/local/small.bed", verbose=False)
    monkeypatch.setitem(sys.modules, "bed_reader.bed_reader", None)
    with pytest.raises(BackendUnavailableError, match="bio-bed-reader"):
        bg.check_cloud_location("https://example.org/x.bed", verbose=False)


# ---- pooch-backed sample wrappers (fakes only) ---------------------------------------------------


def test_sample_bed_file_returns_local_bed_fam_and_bim(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    bed_reader = pytest.importorskip("bed_reader")
    bed = tmp_path / "small.bed"
    bed.write_bytes(b"")
    (tmp_path / "small.fam").write_text("fam\n")
    seen: list[str] = []

    def fake_sample_file(name: str) -> str:
        seen.append(name)
        return str(bed)

    monkeypatch.setattr(bed_reader, "sample_file", fake_sample_file)
    result = bg.sample_bed_file("small.bed", verbose=False)
    assert seen == ["small.bed"]
    assert result["service"] == "bed_reader.sample_file"
    assert result["source"] == "small.bed"
    assert result["result"] == {
        "bed": str(bed.resolve()),
        "fam": str((tmp_path / "small.fam").resolve()),
        "bim": None,
    }


def test_sample_bed_file_reports_missing_pooch(monkeypatch: pytest.MonkeyPatch) -> None:
    bed_reader = pytest.importorskip("bed_reader")

    def without_pooch(name: str) -> str:
        raise ImportError("The function sample_file() requires pooch.")

    monkeypatch.setattr(bed_reader, "sample_file", without_pooch)
    with pytest.raises(BackendUnavailableError, match="bed-reader\\[samples\\]"):
        bg.sample_bed_file("small.bed", verbose=False)


def test_sample_bed_url_returns_the_url(monkeypatch: pytest.MonkeyPatch) -> None:
    bed_reader = pytest.importorskip("bed_reader")
    seen: list[str] = []

    def fake_sample_url(name: str) -> str:
        seen.append(name)
        return "file:///cache/bed_reader/small.bed"

    monkeypatch.setattr(bed_reader, "sample_url", fake_sample_url)
    result = bg.sample_bed_url("small.bed", verbose=False)
    assert seen == ["small.bed"]
    assert result == {
        "service": "bed_reader.sample_url",
        "source": "small.bed",
        "result": {"url": "file:///cache/bed_reader/small.bed"},
    }


@pytest.mark.parametrize("name", ["", "   ", "small\n.bed", "x" * 300])
def test_sample_wrappers_reject_bad_names(name: str) -> None:
    with pytest.raises(ConfigurationError):
        bg.sample_bed_file(name, verbose=False)
    with pytest.raises(ConfigurationError):
        bg.sample_bed_url(name, verbose=False)


# ---- module contract --------------------------------------------------------------------------


def test_public_functions_end_with_keyword_only_verbose() -> None:
    assert sorted(bg.__all__) == [
        "check_cloud_location",
        "read_cloud_genotypes",
        "sample_bed_file",
        "sample_bed_url",
    ]
    for name in bg.__all__:
        parameters = list(inspect.signature(getattr(bg, name)).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name
