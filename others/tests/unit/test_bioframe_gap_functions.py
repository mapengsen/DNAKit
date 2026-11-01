"""Tests for bioframe_gap_functions: UCSC metadata wrappers, pairix reads, UCSC bigBed/bigWig
facades and interval plots.

Network wrappers are tested with fakes that replace bioframe's fetch functions; the socket guard
below blocks any real connection. Local-provider assembly tables ship with bioframe and are read
offline. The UCSC programs and pypairix are not required: tests that need them are skipped with a
reason, while the missing-program paths run everywhere. Facade argument tests use a small stand-in
executable written under tmp_path.
"""

from __future__ import annotations

import shutil
import socket
import sys
import types
from pathlib import Path

import pytest

from dnakit.bioinformatics import bioframe_gap_functions as module
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


def _bioframe():
    return pytest.importorskip("bioframe")


def _is_png(path: Path) -> bool:
    return path.is_file() and path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def _write_stand_in(directory: Path, name: str, body: str) -> Path:
    """Write a POSIX shell script that stands in for a UCSC program; it records its arguments next
    to itself.
    """
    script = directory / name
    script.write_text("#!/bin/sh\n" + body)
    script.chmod(0o755)
    return script


def _empty_path(monkeypatch, tmp_path: Path) -> None:
    empty = tmp_path / "empty_bin"
    empty.mkdir(exist_ok=True)
    monkeypatch.setenv("PATH", str(empty))


def test_public_names_are_exported_and_callable():
    for name in module.__all__:
        assert callable(getattr(module, name)), name
    assert set(module.__all__) == {
        "fetch_centromeres",
        "fetch_chromsizes",
        "plot_intervals",
        "plot_intervals_arr",
        "read_pairix",
        "to_bigbed",
        "to_bigwig",
    }


def test_fetch_centromeres_local_provider_reads_bundled_tables(tmp_path):
    _bioframe()
    out = tmp_path / "centromeres.tsv"
    result = module.fetch_centromeres("hg38", provider="local", output_path=str(out), verbose=False)
    assert result["status"] == "completed"
    assert result["request"] == {"db": "hg38", "provider": "local"}
    assert result["retrieved_at"].endswith("+00:00")
    chr1 = next(row for row in result["data"] if row["chrom"] == "chr1")
    assert chr1["start"] < chr1["mid"] < chr1["end"]
    assert result["row_count"] == len(result["data"])
    assert out.read_text().splitlines()[0].split("\t") == ["chrom", "start", "end", "mid"]


def test_fetch_centromeres_ucsc_uses_bioframe_with_fake(monkeypatch, tmp_path):
    bioframe = _bioframe()
    import pandas as pd

    calls = []

    def fake_fetch(db, provider="local"):
        calls.append((db, provider))
        return pd.DataFrame({"chrom": ["chrX"], "start": [1], "end": [3], "mid": [2]})

    monkeypatch.setattr(bioframe, "fetch_centromeres", fake_fetch)
    result = module.fetch_centromeres("hg19", output_path=str(tmp_path / "cen.tsv"), verbose=False)
    assert calls == [("hg19", "ucsc")]
    assert result["data"] == [{"chrom": "chrX", "start": 1, "end": 3, "mid": 2}]
    assert result["service"].startswith("UCSC")
    assert result["output_path"] == str(tmp_path / "cen.tsv")
    assert (tmp_path / "cen.tsv").is_file()


def test_fetch_centromeres_maps_download_failure_to_backend_error(monkeypatch):
    bioframe = _bioframe()

    def failing(db, provider="local"):
        raise OSError("connection refused")

    monkeypatch.setattr(bioframe, "fetch_centromeres", failing)
    with pytest.raises(BackendExecutionError, match="connection refused"):
        module.fetch_centromeres("hg38", verbose=False)


def test_fetch_centromeres_rejects_bad_input():
    with pytest.raises(ConfigurationError):
        module.fetch_centromeres("../hg38", verbose=False)
    with pytest.raises(ConfigurationError):
        module.fetch_centromeres("hg38", provider="ftp", verbose=False)


def test_fetch_chromsizes_local_provider_as_mapping_and_bed(tmp_path):
    _bioframe()
    mapping = module.fetch_chromsizes("hg38", provider="local", verbose=False)
    assert mapping["data"]["chr1"] == 248956422
    assert mapping["row_count"] == len(mapping["data"])
    table = module.fetch_chromsizes("hg38", provider="local", as_bed=True, verbose=False)
    first = table["data"][0]
    assert first == {"chrom": "chr1", "start": 0, "end": 248956422}

    out = tmp_path / "hg38.chrom.sizes"
    module.fetch_chromsizes("hg38", provider="local", output_path=str(out), verbose=False)
    assert out.read_text().splitlines()[0] == "chr1\t248956422"


def test_fetch_chromsizes_ucsc_uses_bioframe_with_fake(monkeypatch):
    bioframe = _bioframe()
    import pandas as pd

    seen = {}

    def fake_fetch(db, provider="local", as_bed=False, filter_chroms=True):
        seen.update(db=db, provider=provider, as_bed=as_bed, filter_chroms=filter_chroms)
        return pd.Series({"chr1": 10, "chr2": 20}, dtype="int64")

    monkeypatch.setattr(bioframe, "fetch_chromsizes", fake_fetch)
    result = module.fetch_chromsizes("mm10", filter_chroms=False, verbose=False)
    assert seen == {"db": "mm10", "provider": "ucsc", "as_bed": False, "filter_chroms": False}
    assert result["data"] == {"chr1": 10, "chr2": 20}
    assert result["request"]["provider"] == "ucsc"


def test_read_pairix_reports_missing_bindings(monkeypatch, tmp_path):
    source = tmp_path / "pairs.px2.gz"
    source.write_bytes(b"")
    monkeypatch.setitem(sys.modules, "pypairix", None)
    with pytest.raises(BackendUnavailableError) as caught:
        module.read_pairix(str(source), "chr1:0-100", verbose=False)
    assert "dnakit[bio-pairix]" in (caught.value.hint or "")


def test_read_pairix_wrapper_passes_arguments_through(monkeypatch, tmp_path):
    bioframe = _bioframe()
    pd = pytest.importorskip("pandas")
    pytest.importorskip("cytoolz")
    fake_pypairix = types.ModuleType("pypairix")
    monkeypatch.setitem(sys.modules, "pypairix", fake_pypairix)
    source = tmp_path / "pairs.px2.gz"
    source.write_bytes(b"")
    seen = {}

    def fake_read(path, region1, region2=None, chromsizes=None, usecols=None, **kwargs):
        seen.update(
            path=path, region1=region1, region2=region2, chromsizes=chromsizes, usecols=usecols
        )
        return pd.DataFrame({"chr1": ["chr1"], "pos1": [150]})

    monkeypatch.setattr(bioframe, "read_pairix", fake_read)
    result = module.read_pairix(
        str(source),
        "chr1:0-200",
        region2="chr1:0-200",
        chromsizes={"chr1": 1000},
        usecols=["chr1", "pos1"],
        verbose=False,
    )
    assert seen["path"] == str(source)
    assert seen["region2"] == "chr1:0-200"
    assert list(seen["chromsizes"].to_dict().items()) == [("chr1", 1000)]
    assert seen["usecols"] == ["chr1", "pos1"]
    assert result["data"] == [{"chr1": "chr1", "pos1": 150}]
    assert result["row_count"] == 1


def test_read_pairix_real_file_when_pypairix_present(tmp_path):
    pypairix = pytest.importorskip("pypairix")
    pysam = pytest.importorskip("pysam")
    _bioframe()
    text = tmp_path / "pairs.txt"
    text.write_text("r1\tchr1\t100\tchr1\t500\nr2\tchr1\t300\tchr2\t900\nr3\tchr2\t50\tchr2\t80\n")
    gz = tmp_path / "pairs.txt.gz"
    pysam.tabix_compress(str(text), str(gz))
    pypairix.build_index(str(gz), force=1)
    result = module.read_pairix(str(gz), "chr1:0-1000", verbose=False)
    assert result["row_count"] >= 1
    assert result["backend"] == "pypairix"


def test_to_bigbed_reports_missing_program(monkeypatch, tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    _empty_path(monkeypatch, tmp_path)
    with pytest.raises(BackendUnavailableError) as caught:
        module.to_bigbed(str(bed), str(sizes), str(tmp_path / "out.bb"), verbose=False)
    assert "ucsc-bedtobigbed" in (caught.value.hint or "")


def test_to_bigwig_reports_missing_program(monkeypatch, tmp_path):
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1.5\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    _empty_path(monkeypatch, tmp_path)
    with pytest.raises(BackendUnavailableError) as caught:
        module.to_bigwig(str(graph), str(sizes), str(tmp_path / "out.bw"), verbose=False)
    assert "ucsc-bedgraphtobigwig" in (caught.value.hint or "")


def test_to_bigbed_builds_fixed_arguments_with_stand_in(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\tn\t5\t+\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    out = tmp_path / "out.bb"
    tool = _write_stand_in(
        tmp_path,
        "bedToBigBed",
        'for last in "$@"; do :; done\n'
        'printf "bigbed" > "$last"\n'
        'printf "%s\\n" "$@" > "$(dirname "$0")/argv.txt"\n',
    )
    result = module.to_bigbed(
        str(bed), str(sizes), str(out), bed_format="bed6", bedtobigbed_path=str(tool), verbose=False
    )
    assert result["status"] == "completed"
    assert result["return_code"] == 0
    assert result["outputs"] == [str(out)]
    assert result["command"][0] == "bedToBigBed"
    assert (tmp_path / "argv.txt").read_text().splitlines() == [
        "-type=bed6",
        str(bed.resolve()),
        str(sizes.resolve()),
        str(out.resolve()),
    ]


def test_to_bigwig_reports_failed_run_with_stand_in(tmp_path):
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1.5\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    tool = _write_stand_in(tmp_path, "bedGraphToBigWig", 'echo "bad input" >&2\nexit 3\n')
    result = module.to_bigwig(
        str(graph),
        str(sizes),
        str(tmp_path / "out.bw"),
        bedgraphtobigwig_path=str(tool),
        verbose=False,
    )
    assert result["status"] == "failed"
    assert result["return_code"] == 3
    assert result["outputs"] == []
    assert "bad input" in result["log_tail"]


def test_bigwig_bigbed_refuse_bad_inputs(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    existing = tmp_path / "exists.bb"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        module.to_bigbed(str(bed), str(sizes), str(existing), verbose=False)
    with pytest.raises(ConfigurationError):
        module.to_bigbed(
            str(bed), str(sizes), str(tmp_path / "new.bb"), bed_format="bed6; rm", verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.to_bigwig(
            str(tmp_path / "missing.bedGraph"), str(sizes), str(tmp_path / "x.bw"), verbose=False
        )


@pytest.mark.skipif(shutil.which("bedToBigBed") is None, reason="UCSC bedToBigBed is not installed")
def test_to_bigbed_real_program_writes_bigbed(tmp_path):
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\tfeat\t0\t+\nchr1\t20\t30\tfeat2\t0\t-\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    out = tmp_path / "out.bb"
    result = module.to_bigbed(str(bed), str(sizes), str(out), bed_format="bed6", verbose=False)
    assert result["status"] == "completed"
    assert out.read_bytes()[:4] == bytes([0xEB, 0xF2, 0x89, 0x87])


@pytest.mark.skipif(
    shutil.which("bedGraphToBigWig") is None, reason="UCSC bedGraphToBigWig is not installed"
)
def test_to_bigwig_real_program_writes_bigwig(tmp_path):
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1.5\nchr1\t10\t20\t2.5\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    out = tmp_path / "out.bw"
    result = module.to_bigwig(str(graph), str(sizes), str(out), verbose=False)
    assert result["status"] == "completed"
    assert out.read_bytes()[:4] == bytes([0x26, 0xFC, 0x8F, 0x88])


def _interval_file(path: Path) -> Path:
    path.write_text("chr1\t100\t500\tA\nchr1\t300\t900\tB\nchr2\t0\t50\tC\n")
    return path


def test_plot_intervals_writes_png_and_data(tmp_path):
    _bioframe()
    source = _interval_file(tmp_path / "iv.bed")
    out = tmp_path / "tracks.png"
    result = module.plot_intervals(str(source), str(out), label_column=3, verbose=False)
    assert result["status"] == "completed"
    assert result["image_path"] == str(out)
    assert _is_png(out)
    data = result["data"]
    assert data["chromosomes"] == ["chr1", "chr2"]
    assert data["interval_count"] == 3
    chr1_levels = sorted(row["level"] for row in data["intervals"] if row["chrom"] == "chr1")
    assert chr1_levels == [0, 1]
    assert {row["label"] for row in data["intervals"]} == {"A", "B", "C"}


def test_plot_intervals_refuses_existing_output_and_bad_columns(tmp_path):
    _bioframe()
    source = _interval_file(tmp_path / "iv.bed")
    existing = tmp_path / "exists.png"
    existing.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        module.plot_intervals(str(source), str(existing), verbose=False)
    with pytest.raises(ConfigurationError):
        module.plot_intervals(str(source), str(tmp_path / "new.png"), label_column=9, verbose=False)
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\t500\t100\n")
    with pytest.raises(ConfigurationError):
        module.plot_intervals(str(bad), str(tmp_path / "bad.png"), verbose=False)


def test_plot_intervals_arr_auto_levels_and_options(tmp_path):
    _bioframe()
    out = tmp_path / "arr.png"
    result = module.plot_intervals_arr(
        [0, 2, 20],
        [10, 5, 30],
        str(out),
        labels=["x", None, "z"],
        colors=["red", "blue", "green"],
        xlim=[-1, 35],
        verbose=False,
    )
    assert _is_png(out)
    data = result["data"]
    assert data["levels"] == [0, 1, 0]
    assert data["xlim"] == [-1, 35]
    assert data["colors"] == ["red", "blue", "green"]
    assert data["ylim"] == [-0.5, 1.5]


def test_plot_intervals_arr_explicit_levels_and_default_colour(tmp_path):
    _bioframe()
    out = tmp_path / "levels.png"
    result = module.plot_intervals_arr([1, 4], [3, 9], str(out), levels=[2, 0], verbose=False)
    assert _is_png(out)
    assert result["data"]["levels"] == [2, 0]
    assert result["data"]["colors"] == ["skyblue", "skyblue"]
    assert result["data"]["xlim"] == [-0.5, 9.5]


def test_plot_intervals_arr_rejects_bad_input(tmp_path):
    _bioframe()
    with pytest.raises(ConfigurationError):
        module.plot_intervals_arr([0, 1], [5], str(tmp_path / "a.png"), verbose=False)
    with pytest.raises(ConfigurationError):
        module.plot_intervals_arr([9], [3], str(tmp_path / "b.png"), verbose=False)
    with pytest.raises(ConfigurationError):
        module.plot_intervals_arr([0], [3], str(tmp_path / "c.png"), levels=[0, 1], verbose=False)
    with pytest.raises(ConfigurationError):
        module.plot_intervals_arr([0], [3], str(tmp_path / "d.png"), xlim=[5, 1], verbose=False)


def test_plotting_leaves_pyplot_figures_alone(tmp_path):
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    _bioframe()
    before = plt.get_fignums()
    module.plot_intervals_arr([0], [5], str(tmp_path / "pure.png"), verbose=False)
    module.plot_intervals(
        str(_interval_file(tmp_path / "iv.bed")), str(tmp_path / "pure2.png"), verbose=False
    )
    assert plt.get_fignums() == before
