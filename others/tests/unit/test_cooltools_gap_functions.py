"""Tests for dnakit.epigenome.cooltools_gap_functions (cooltools example data, file downloads and saddle plots).

Download functions replace the one upstream file fetch (download_file) with a recording fake that writes known bytes,
so no request is sent. Dataset listings are read from a small local table, which cooltools also accepts, so those
tests need no fake at all. The saddle plot runs for real on a synthetic track built in the test. Tests skip cleanly
when cooltools or matplotlib is missing, and the socket guard turns any network attempt into a test failure.
"""

from __future__ import annotations

import hashlib
import inspect
import random
import socket
import sys
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from dnakit.epigenome import cooltools_gap_functions as cg
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

@pytest.fixture(autouse=True)
def _close_figures_left_by_other_tests():
    """Other tests in the same process may leave figures open; these tests check that the function closes its own."""
    plt = pytest.importorskip("matplotlib.pyplot")
    plt.close("all")
    yield


PAYLOAD = b"cooltools test payload\n"
OTHER = b"second dataset\n"
MD5 = hashlib.md5(PAYLOAD).hexdigest()
DEMO_URL = "https://example.invalid/demo.txt"
OTHER_URL = "https://example.invalid/other.txt"


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


class Recorder:
    """A fake upstream call that records its arguments and returns a canned value."""

    def __init__(self, value: Any = None, *, effect: Callable[..., Any] | None = None) -> None:
        self.value = value
        self.effect = effect
        self.calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append((args, kwargs))
        if self.effect is not None:
            return self.effect(*args, **kwargs)
        return self.value


@pytest.fixture
def cooltools_io() -> Any:
    pytest.importorskip("cooltools", reason="cooltools is not installed")
    return pytest.importorskip("cooltools.lib.io", reason="cooltools is not installed")


def _write_payload(url: str, local_filename: str) -> str:
    """Stands in for cooltools' download_file: writes the known bytes for the requested URL."""
    Path(local_filename).write_bytes(PAYLOAD if "demo" in url else OTHER)
    return local_filename


def _listing(path: Path, *, demo_checksum: str = MD5) -> Path:
    """A dataset table in cooltools' format: a '#' header line, then tab-separated rows."""
    path.write_text(
        "# key\tcomment\tlink\tfilename\tchecksum\n"
        f"demo\tdemo file\t{DEMO_URL}\tdemo.txt\t{demo_checksum}\n"
        f"other\tsecond file\t{OTHER_URL}\tother.txt\t{hashlib.md5(OTHER).hexdigest()}\n",
        encoding="utf-8",
    )
    return path


# --- module contract ---------------------------------------------------------------------------


def test_public_functions_have_keyword_only_verbose_and_docstrings() -> None:
    assert len(cg.__all__) == len(set(cg.__all__))
    for name in cg.__all__:
        function = getattr(cg, name)
        assert callable(function), name
        signature = inspect.signature(function)
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert signature.parameters["verbose"].default is True, name
        assert function.__doc__ and function.__doc__.strip(), name
        for parameter in signature.parameters.values():
            assert parameter.annotation is not inspect.Parameter.empty, f"{name}.{parameter.name}"


def test_data_dir_is_required_so_no_file_lands_in_the_package_tree() -> None:
    parameter = inspect.signature(cg.download_cooltools_example_data).parameters["data_dir"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is inspect.Parameter.empty


def test_aliases_are_the_same_upstream_objects(cooltools_io: Any) -> None:
    import cooltools

    assert cooltools.download_data is cooltools_io.download_data
    assert cooltools.lib.download_data is cooltools_io.download_data
    assert cooltools.lib.io.download_data is cooltools_io.download_data
    assert cooltools.print_available_datasets is cooltools_io.print_available_datasets
    assert cooltools.lib.print_available_datasets is cooltools_io.print_available_datasets
    assert cooltools.lib.io.print_available_datasets is cooltools_io.print_available_datasets
    assert cooltools.download_file is cooltools_io.download_file
    assert cooltools.lib.download_file is cooltools_io.download_file
    assert cooltools.lib.io.download_file is cooltools_io.download_file


def test_cli_helpers_are_shared_across_command_modules() -> None:
    """The validate_csv and sniff_for_header rows are one helper each, re-imported by several CLI modules."""
    cooltools_cli_util = pytest.importorskip(
        "cooltools.cli.util", reason="cooltools is not installed"
    )
    from cooltools.cli import dots, eigs_cis, pileup

    assert dots.validate_csv is cooltools_cli_util.validate_csv
    assert pileup.validate_csv is cooltools_cli_util.validate_csv
    assert eigs_cis.sniff_for_header is cooltools_cli_util.sniff_for_header


def test_missing_cooltools_reports_the_optional_extra(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # A None entry makes the import fail the way a missing package does.
    monkeypatch.setitem(sys.modules, "cooltools.lib.io", None)
    with pytest.raises(BackendUnavailableError) as raised:
        cg.download_cooltools_file("https://example.invalid/x.txt", str(tmp_path / "x.txt"))
    assert raised.value.hint is not None and "bio-cooltools" in raised.value.hint
    assert not (tmp_path / "x.txt").exists()


# --- example datasets: download (fake fetch) ----------------------------------------------------


def test_named_dataset_is_fetched_into_the_folder_and_checked(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    table = _listing(tmp_path / "table.tsv")
    folder = tmp_path / "data"

    result = cg.download_cooltools_example_data(
        "demo", data_dir=str(folder), info_url=str(table), verbose=False
    )

    assert len(fake.calls) == 1
    (url, local), _ = fake.calls[0]
    assert url == DEMO_URL
    assert Path(local) == folder / "demo.txt"
    assert (folder / "demo.txt").read_bytes() == PAYLOAD
    assert result["status"] == "complete"
    assert result["backend"] == "cooltools"
    assert result["data"]["files"] == [str(folder / "demo.txt")]
    assert result["data"]["file_count"] == 1
    assert result["request"]["info_url"] == str(table)
    assert result["request"]["cache"] is True


def test_cached_dataset_is_reused_without_a_second_fetch(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    table = _listing(tmp_path / "table.tsv")
    folder = tmp_path / "data"
    first = cg.download_cooltools_example_data("demo", data_dir=str(folder), info_url=str(table))
    second = cg.download_cooltools_example_data("demo", data_dir=str(folder), info_url=str(table))

    assert len(fake.calls) == 1
    assert first["data"]["files"] == second["data"]["files"]


def test_checksum_mismatch_is_reported_as_backend_error(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cooltools_io, "download_file", Recorder(effect=_write_payload))
    table = _listing(tmp_path / "table.tsv", demo_checksum="0" * 32)
    with pytest.raises(BackendExecutionError, match="differs from original test|failed"):
        cg.download_cooltools_example_data(
            "demo", data_dir=str(tmp_path / "data"), info_url=str(table)
        )


def test_all_datasets_are_fetched_and_listed(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    table = _listing(tmp_path / "table.tsv")

    result = cg.download_cooltools_example_data(
        "all", data_dir=str(tmp_path / "data"), info_url=str(table)
    )

    assert len(fake.calls) == 2
    assert result["data"]["file_count"] == 2
    assert {Path(path).name for path in result["data"]["files"]} == {"demo.txt", "other.txt"}


def test_unknown_dataset_name_is_reported_by_cooltools(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    table = _listing(tmp_path / "table.tsv")
    with pytest.raises(BackendExecutionError, match="is not available"):
        cg.download_cooltools_example_data(
            "missing", data_dir=str(tmp_path / "data"), info_url=str(table)
        )
    assert fake.calls == []


def test_download_arguments_are_refused_before_any_fetch(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    table = _listing(tmp_path / "table.tsv")
    not_a_folder = tmp_path / "occupied"
    not_a_folder.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        cg.download_cooltools_example_data("", data_dir=str(tmp_path / "data"), info_url=str(table))
    with pytest.raises(ConfigurationError):
        cg.download_cooltools_example_data("demo", data_dir=str(not_a_folder), info_url=str(table))
    with pytest.raises(ConfigurationError):
        cg.download_cooltools_example_data(
            "demo", data_dir=str(tmp_path / "data"), info_url="ftp://example.invalid/t"
        )
    with pytest.raises(ConfigurationError):
        cg.download_cooltools_example_data(
            "demo",
            data_dir=str(tmp_path / "data"),
            info_url=str(table),
            cache="yes",  # type: ignore[arg-type]
        )
    assert fake.calls == []


# --- example datasets: listing ------------------------------------------------------------------


def test_listing_from_a_local_table_needs_no_network(cooltools_io: Any, tmp_path: Path) -> None:
    table = _listing(tmp_path / "table.tsv")

    result = cg.list_cooltools_example_datasets(url=str(table), verbose=False)

    datasets = result["data"]["datasets"]
    assert result["data"]["count"] == 2
    assert [row["key"] for row in datasets] == ["demo", "other"]
    assert datasets[0] == {
        "index": datasets[0]["index"],
        "key": "demo",
        "comment": "demo file",
        "link": DEMO_URL,
        "filename": "demo.txt",
        "checksum": MD5,
    }
    assert isinstance(datasets[0]["index"], int)
    assert result["request"] == {"url": str(table)}


def test_default_listing_reads_cooltools_published_table(
    cooltools_io: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(
        [
            {
                "index": 1,
                "key": "remote",
                "comment": "c",
                "link": DEMO_URL,
                "filename": "r.txt",
                "checksum": MD5,
            }
        ]
    )
    monkeypatch.setattr(cooltools_io, "_get_datasets_info", fake)

    result = cg.list_cooltools_example_datasets()

    assert fake.calls == [((cooltools_io.URL_DATA,), {})]
    assert result["data"]["datasets"][0]["key"] == "remote"


def test_listing_refuses_missing_local_files_and_corrupt_tables(
    cooltools_io: Any, tmp_path: Path
) -> None:
    with pytest.raises(ConfigurationError):
        cg.list_cooltools_example_datasets(url=str(tmp_path / "absent.tsv"))
    corrupt = tmp_path / "corrupt.tsv"
    corrupt.write_text("demo\tno header line\n", encoding="utf-8")
    with pytest.raises(BackendExecutionError, match="corrupted"):
        cg.list_cooltools_example_datasets(url=str(corrupt))


# --- single-file download -----------------------------------------------------------------------


def test_single_file_download_records_size_and_md5(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    target = tmp_path / "nested" / "demo.txt"

    result = cg.download_cooltools_file(DEMO_URL, str(target), verbose=False)

    assert fake.calls == [((DEMO_URL, str(target)), {})]
    assert target.read_bytes() == PAYLOAD
    assert result["data"] == {"output_path": str(target), "bytes": len(PAYLOAD), "md5": MD5}
    assert result["request"] == {"url": DEMO_URL, "output_path": str(target)}


def test_single_file_download_refuses_existing_output_and_non_http_urls(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = Recorder(effect=_write_payload)
    monkeypatch.setattr(cooltools_io, "download_file", fake)
    existing = tmp_path / "existing.txt"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="already exists"):
        cg.download_cooltools_file(DEMO_URL, str(existing))
    for bad in ("ftp://example.invalid/x", "file:///etc/hosts", str(tmp_path / "local.txt"), ""):
        with pytest.raises(ConfigurationError):
            cg.download_cooltools_file(bad, str(tmp_path / "new.txt"))
    assert fake.calls == []
    assert existing.read_text(encoding="utf-8") == "keep"


def test_single_file_download_failures_leave_no_partial_output(
    cooltools_io: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def _http_error(url: str, local_filename: str) -> str:
        raise OSError(f"404 Client Error for {url}")

    monkeypatch.setattr(cooltools_io, "download_file", Recorder(effect=_http_error))
    target = tmp_path / "missing.txt"
    with pytest.raises(BackendExecutionError, match="404"):
        cg.download_cooltools_file(DEMO_URL, str(target))
    assert not target.exists()

    monkeypatch.setattr(cooltools_io, "download_file", Recorder(effect=lambda url, local: local))
    with pytest.raises(BackendExecutionError, match="wrote no data"):
        cg.download_cooltools_file(DEMO_URL, str(tmp_path / "empty.txt"))


# --- saddle plot --------------------------------------------------------------------------------


@pytest.fixture
def pyplot() -> Any:
    pytest.importorskip("cooltools", reason="cooltools is not installed")
    return pytest.importorskip("matplotlib.pyplot", reason="matplotlib is not installed")


def _track(rows: int = 200, *, missing: int = 5) -> list[list[object]]:
    """A sorted, non-overlapping synthetic track with some missing values."""
    rng = random.Random(0)
    values: list[object] = [round(rng.uniform(0.2, 3.0), 3) for _ in range(rows)]
    for index in range(missing):
        values[index * 7] = None
    return [["chr1", index * 100, index * 100 + 100, value] for index, value in enumerate(values)]


def _saddle(side: int) -> list[list[float]]:
    return [[1.0 + 0.05 * (row + column) for column in range(side)] for row in range(side)]


def test_saddle_plot_writes_a_png_and_the_digitised_counts(pyplot: Any, tmp_path: Path) -> None:
    track = _track()
    output = tmp_path / "saddle.png"

    result = cg.plot_cooltools_saddle(track, _saddle(6), 4, str(output), verbose=False)

    assert result["status"] == "completed"
    assert result["image_path"] == str(output)
    assert output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    data = result["data"]
    assert data["track_rows"] == 200
    assert len(data["bin_counts"]) == 4
    total = (
        sum(data["bin_counts"])
        + data["left_outlier_count"]
        + data["right_outlier_count"]
        + data["missing_count"]
    )
    assert total == 200
    assert data["missing_count"] == 5
    assert data["heatmap_shape"] == [4, 4]
    assert data["heatmap_values"][0] == pytest.approx(
        _saddle(6)[1][1:5]
    )  # outlier flanks are trimmed
    assert any("deprecat" in message.lower() for message in result["upstream_warnings"])
    assert pyplot.get_fignums() == []  # the figure is closed after saving


def test_upstream_saddleplot_fails_without_outlier_flanks(pyplot: Any) -> None:
    """Reproduces why the wrapper requires the flanked matrix: upstream fails on an n x n matrix."""
    saddle_api = pytest.importorskip("cooltools.api.saddle", reason="cooltools is not installed")
    import numpy as np
    import pandas as pd

    track = pd.DataFrame(_track(), columns=["chrom", "start", "end", "value"])
    with warnings.catch_warnings(), pytest.raises(ValueError, match="shape mismatch"):
        warnings.simplefilter(
            "ignore", DeprecationWarning
        )  # upstream's notice that saddleplots will be deprecated
        saddle_api.saddleplot(track, np.array(_saddle(4)), 4, qrange=(0.0, 1.0))
    pyplot.close(
        "all"
    )  # upstream leaves its half-drawn figure open; the wrapper closes its own figures


def test_unflanked_matrix_is_refused_before_drawing(pyplot: Any, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="outlier flanks"):
        cg.plot_cooltools_saddle(_track(), _saddle(4), 4, str(tmp_path / "plain.png"))
    assert not (tmp_path / "plain.png").exists()


def test_value_limits_replace_quantile_limits(pyplot: Any, tmp_path: Path) -> None:
    import numpy as np

    result = cg.plot_cooltools_saddle(
        _track(),
        _saddle(6),
        4,
        str(tmp_path / "vrange.png"),
        vrange=[0.2, 3.0],
        qrange=None,
        verbose=False,
    )
    assert result["data"]["vrange"] == [0.2, 3.0]
    assert result["data"]["qrange"] is None
    assert result["data"]["bin_edges"] == pytest.approx(list(np.linspace(0.2, 3.0, 5)))


def test_linear_scale_draws_too(pyplot: Any, tmp_path: Path) -> None:
    result = cg.plot_cooltools_saddle(
        _track(),
        _saddle(6),
        4,
        str(tmp_path / "linear.png"),
        scale="linear",
        vmin=0.5,
        vmax=2.0,
        verbose=False,
    )
    assert result["data"]["scale"] == "linear"
    assert (tmp_path / "linear.png").is_file()


def test_saddle_plot_argument_checks_run_before_drawing(pyplot: Any, tmp_path: Path) -> None:
    track, matrix = _track(), _saddle(6)
    with pytest.raises(ConfigurationError, match="exactly one"):
        cg.plot_cooltools_saddle(track, matrix, 4, str(tmp_path / "a.png"), qrange=None)
    with pytest.raises(ConfigurationError, match="exactly one"):
        cg.plot_cooltools_saddle(track, matrix, 4, str(tmp_path / "a.png"), vrange=[0.2, 3.0])
    with pytest.raises(ConfigurationError, match="outlier flanks"):
        cg.plot_cooltools_saddle(track, _saddle(5), 4, str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="square"):
        cg.plot_cooltools_saddle(
            track, [[1.0] * 4, [1.0] * 3, [1.0] * 4, [1.0] * 4], 2, str(tmp_path / "a.png")
        )
    with pytest.raises(ConfigurationError, match="scale"):
        cg.plot_cooltools_saddle(track, matrix, 4, str(tmp_path / "a.png"), scale="sqrt")
    with pytest.raises(ConfigurationError, match="above zero"):
        cg.plot_cooltools_saddle(track, matrix, 4, str(tmp_path / "a.png"), vmin=0.0)
    with pytest.raises(ConfigurationError, match="\\.png"):
        cg.plot_cooltools_saddle(track, matrix, 4, str(tmp_path / "a.pdf"))
    with pytest.raises(ConfigurationError, match="end greater than start"):
        cg.plot_cooltools_saddle([["chr1", 100, 100, 1.0]], matrix, 4, str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="four values"):
        cg.plot_cooltools_saddle([["chr1", 0, 100]], matrix, 4, str(tmp_path / "a.png"))
    with pytest.raises(ConfigurationError, match="positive integer"):
        cg.plot_cooltools_saddle(track, matrix, 0, str(tmp_path / "a.png"))
    assert not any(tmp_path.glob("a.*"))


def test_existing_saddle_output_is_refused(pyplot: Any, tmp_path: Path) -> None:
    output = tmp_path / "exists.png"
    output.write_bytes(b"keep")
    with pytest.raises(ConfigurationError, match="already exists"):
        cg.plot_cooltools_saddle(_track(), _saddle(6), 4, str(output))
    assert output.read_bytes() == b"keep"


def test_missing_matplotlib_reports_the_plot_extra(
    cooltools_io: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", None)
    with pytest.raises(BackendUnavailableError) as raised:
        cg.plot_cooltools_saddle(_track(), _saddle(6), 4, str(tmp_path / "x.png"))
    assert raised.value.hint is not None and "bio-plots" in raised.value.hint
    assert not (tmp_path / "x.png").exists()
