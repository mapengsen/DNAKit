"""Tests for dnakit.epigenome.deeptools_plot_functions (deepTools plots, matrix heatmap, template
lengths, blacklist).

Every input is written by the test itself under tmp_path; nothing is downloaded and no socket is
opened.
"""

from __future__ import annotations

import gzip
import inspect
import json
import math
import socket
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.epigenome import deeptools_plot_functions as plots
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(autouse=True)
def _no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("network access is not allowed in this test")

    monkeypatch.setattr(socket, "socket", refuse)


def _assert_png(path: str) -> None:
    target = Path(path)
    assert target.is_absolute()
    assert target.is_file()
    assert target.read_bytes()[:8] == PNG_SIGNATURE


def _open_figures() -> list[int]:
    plt = pytest.importorskip("matplotlib.pyplot")
    return list(plt.get_fignums())


def _write_matrix(path: Path, *, groups: tuple[int, int] = (2, 2)) -> Path:
    """A gzipped computeMatrix-style file: a JSON header line, then BED-like rows with values."""

    rows_per_group = list(groups)
    total = sum(rows_per_group)
    params = {
        "group_boundaries": [0, rows_per_group[0], total],
        "sample_boundaries": [0, 3],
        "group_labels": ["g1", "g2"],
        "sample_labels": ["s1"],
        "ref point": ["TSS"],
        "bin size": [10],
        "upstream": [30],
        "downstream": [30],
        "body": [0],
        "unscaled 5 prime": [0],
        "unscaled 3 prime": [0],
        "sort regions": "no",
        "sort using": "mean",
    }
    lines = ["@" + json.dumps(params, separators=(",", ":"))]
    for index in range(total):
        values = [float(index + column) for column in range(3)]
        fields = ["chr1", str(index * 100), str(index * 100 + 60), f"r{index}", "0", "+"]
        lines.append("\t".join(fields + [f"{value:f}" for value in values]))
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _aligned_read(
    pysam: Any, header: Any, name: str, start: int, cigar: list[tuple[int, int]], tlen: int
) -> Any:
    segment = pysam.AlignedSegment(header)
    query_length = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    segment.query_name = name
    segment.query_sequence = "A" * query_length
    segment.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    segment.flag = 0
    segment.reference_id = 0
    segment.reference_start = start
    segment.cigartuples = cigar
    segment.mapping_quality = 30
    segment.template_length = tlen
    return segment


def _write_bam(path: Path) -> Path:
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        handle.write(_aligned_read(pysam, handle.header, "with_tlen", 10, [(0, 50)], -200))
        handle.write(
            _aligned_read(pysam, handle.header, "gapped", 10, [(0, 20), (3, 100), (0, 30)], 0)
        )
        handle.write(
            _aligned_read(pysam, handle.header, "inserted", 10, [(0, 10), (1, 5), (0, 10)], 0)
        )
    return path


# ---- fragment-length density -------------------------------------------------------------


def test_fragment_length_density_known_values() -> None:
    pytest.importorskip("deeptools.bamPEFragmentSize")
    result = plots.fragment_length_density([150, 150, 150, 250], 0, 300, verbose=False)
    assert result["status"] == "completed"
    assert result["fragments"] == 4
    assert result["bins"] == 100
    assert result["range"] == [0.0, 300.0]
    edges = result["bin_edges"]
    assert len(edges) == 101
    assert edges[0] == 0.0 and edges[-1] == 300.0
    density = result["density"]
    # bins are 3 bp wide: 150 falls in bin 50 (three fragments), 250 in bin 83 (one fragment)
    assert density[50] == pytest.approx(0.25)
    assert density[83] == pytest.approx(1 / 12)
    assert sum(density) * 3.0 == pytest.approx(1.0)


def test_fragment_length_density_does_not_clear_callers_figure() -> None:
    plt = pytest.importorskip("matplotlib.pyplot")
    pytest.importorskip("deeptools.bamPEFragmentSize")
    figure = plt.figure()
    try:
        axes = figure.add_subplot()
        axes.plot([0, 1], [0, 1])
        before = _open_figures()
        plots.fragment_length_density([100, 200], 0, 300, verbose=False)
        assert _open_figures() == before
        assert len(figure.axes) == 1
        assert len(figure.axes[0].lines) == 1
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("lengths", "low", "high"),
    [
        ([], 0, 300),
        ([100, -5], 0, 300),
        ([100], 300, 300),
        ([100], 0, True),
        ("100,200", 0, 300),
    ],
)
def test_fragment_length_density_rejects_bad_input(
    lengths: object, low: object, high: object
) -> None:
    with pytest.raises(ConfigurationError):
        plots.fragment_length_density(lengths, low, high, verbose=False)  # type: ignore[arg-type]


def test_fragment_length_density_missing_dependency_is_reported(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(sys.modules, "deeptools.bamPEFragmentSize", None)
    with pytest.raises(BackendUnavailableError, match="bio-deeptools"):
        plots.fragment_length_density([100], 0, 300, verbose=False)


# ---- GC-bias plot ------------------------------------------------------------------------


def _gc_inputs() -> tuple[list[list[float]], list[list[float]]]:
    frequencies = [[100.0, 50.0, 0.5 + index * 0.01] for index in range(101)]
    # one region per GC step from 0.2 to 0.7, skipping the step that lands on 0.5
    reads = [[float(k + 1), 0.2 + k * 0.01] for k in range(51) if k != 30]
    # three regions at exactly GC 0.5 with reads 4, 6 and 8: that bin has median 6
    reads += [[4.0, 0.5], [6.0, 0.5], [8.0, 0.5]]
    return frequencies, reads


def test_plot_gc_bias_writes_png_and_summary(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.computeGCBias")
    frequencies, reads = _gc_inputs()
    result = plots.plot_gc_bias(
        frequencies, reads, 200, str(tmp_path / "gc" / "bias.png"), verbose=False
    )
    _assert_png(result["path"])
    data = result["data"]
    assert data["region_size"] == 200
    assert len(data["log2_observed_expected"]) == 101
    assert data["log2_observed_expected"][0] == pytest.approx(-1.0)
    assert data["log2_observed_expected"][-1] == pytest.approx(math.log2(1.5))
    assert data["gc_fraction"][0] == 0.0 and data["gc_fraction"][-1] == 1.0
    assert data["window"] == [0.2, 0.7]
    gc_values = [entry["gc"] for entry in data["gc_bins"]]
    assert all(0.2 <= gc <= 0.7 for gc in gc_values)
    at_half = [entry for entry in data["gc_bins"] if entry["gc"] == pytest.approx(0.5)]
    assert len(at_half) == 1
    assert at_half[0]["regions"] == 3
    assert at_half[0]["median_reads"] == pytest.approx(6.0)


def test_plot_gc_bias_refuses_existing_output_and_non_png(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.computeGCBias")
    frequencies, reads = _gc_inputs()
    existing = tmp_path / "bias.png"
    existing.write_bytes(b"old")
    with pytest.raises(ConfigurationError, match="already exists"):
        plots.plot_gc_bias(frequencies, reads, 200, str(existing), verbose=False)
    with pytest.raises(ConfigurationError, match=r"\.png"):
        plots.plot_gc_bias(frequencies, reads, 200, str(tmp_path / "bias.pdf"), verbose=False)
    assert existing.read_bytes() == b"old"


@pytest.mark.parametrize(
    ("frequencies", "reads", "size"),
    [
        ([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], [[1.0, 0.5]], 200),
        ([[1.0, 1.0, 1.0]] * 3, [[1.0, 0.5]], 0),
        ([[1.0, 1.0, 0.0]] * 3, [[1.0, 0.5]], 200),
        ([[1.0, 1.0]] * 3, [[1.0, 0.5]], 200),
        ([[1.0, 1.0, 1.0]] * 3, [[1.0, 1.5]], 200),
        ([[1.0, 1.0, 1.0]] * 3, [], 200),
    ],
)
def test_plot_gc_bias_rejects_bad_input(
    tmp_path: Path, frequencies: list[list[float]], reads: list[list[float]], size: int
) -> None:
    with pytest.raises(ConfigurationError):
        plots.plot_gc_bias(frequencies, reads, size, str(tmp_path / "bias.png"), verbose=False)
    assert not (tmp_path / "bias.png").exists()


# ---- matrix heatmap ----------------------------------------------------------------------


def test_plot_matrix_heatmap_reads_computematrix_output(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotHeatmap")
    matrix = _write_matrix(tmp_path / "matrix.gz")
    result = plots.plot_matrix_heatmap(
        str(matrix), str(tmp_path / "out" / "heat.png"), verbose=False
    )
    _assert_png(result["path"])
    data = result["data"]
    assert data["samples"] == ["s1"]
    assert data["groups"] == ["g1", "g2"]
    assert data["regions"] == 4
    assert data["group_boundaries"] == [0, 2, 4]
    assert data["group_sizes"] == [2, 2]
    assert data["columns_per_sample"] == [3]
    assert data["parameters"]["ref point"] == ["TSS"]
    assert data["average_type"] == "median"


def test_plot_matrix_heatmap_accepts_uneven_groups(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotHeatmap")
    matrix = _write_matrix(tmp_path / "matrix.gz", groups=(3, 1))
    result = plots.plot_matrix_heatmap(
        str(matrix), str(tmp_path / "heat.png"), z_min=0.5, z_max=4.0, per_group=True, verbose=False
    )
    _assert_png(result["path"])
    assert result["data"]["group_sizes"] == [3, 1]
    assert result["data"]["z_min"] == 0.5 and result["data"]["z_max"] == 4.0


def test_plot_matrix_heatmap_rejects_non_matrix_and_bad_options(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotHeatmap")
    plain = tmp_path / "plain.gz"
    with gzip.open(plain, "wt") as handle:
        handle.write("chr1\t0\t10\tr\t0\t+\t1.0\t2.0\n")
    with pytest.raises(ConfigurationError, match="not a computeMatrix output"):
        plots.plot_matrix_heatmap(str(plain), str(tmp_path / "a.png"), verbose=False)
    matrix = _write_matrix(tmp_path / "matrix.gz")
    with pytest.raises(ConfigurationError, match="z_max"):
        plots.plot_matrix_heatmap(
            str(matrix), str(tmp_path / "b.png"), z_min=2, z_max=1, verbose=False
        )
    with pytest.raises(ConfigurationError, match="average_type"):
        plots.plot_matrix_heatmap(
            str(matrix), str(tmp_path / "c.png"), average_type="mode", verbose=False
        )
    with pytest.raises(ConfigurationError, match="what_to_show"):
        plots.plot_matrix_heatmap(
            str(matrix), str(tmp_path / "d.png"), what_to_show="all", verbose=False
        )
    with pytest.raises(ConfigurationError, match="does not exist"):
        plots.plot_matrix_heatmap(
            str(tmp_path / "missing.gz"), str(tmp_path / "e.png"), verbose=False
        )


def test_plot_matrix_heatmap_unknown_colormap_closes_figures(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotHeatmap")
    matrix = _write_matrix(tmp_path / "matrix.gz")
    before = _open_figures()
    with pytest.raises(ConfigurationError, match="could not draw"):
        plots.plot_matrix_heatmap(
            str(matrix), str(tmp_path / "heat.png"), colormap="no-such-colormap", verbose=False
        )
    assert _open_figures() == before
    assert not (tmp_path / "heat.png").exists()


def test_importing_plotting_module_keeps_numpy_error_state(monkeypatch: pytest.MonkeyPatch) -> None:
    np = pytest.importorskip("numpy")
    plt = pytest.importorskip("matplotlib.pyplot")
    pytest.importorskip("deeptools.plotHeatmap")
    # deepTools' plotHeatmap silences all numpy floating-point errors and turns interactive mode off
    # when
    # it is imported; _plotting_module must undo both, so a fresh import is needed to see the
    # effect.
    monkeypatch.delitem(sys.modules, "deeptools.plotHeatmap", raising=False)
    saved_errors = np.geterr()
    np.seterr(all="warn")
    plt.ion()
    try:
        plots._plotting_module("deeptools.plotHeatmap")
        assert np.geterr()["divide"] == "warn"
        assert np.geterr()["invalid"] == "warn"
        assert plt.isinteractive() is True
    finally:
        plt.ioff()
        np.seterr(**saved_errors)


# ---- enrichment bars ---------------------------------------------------------------------


def test_plot_enrichment_counts_percentages_and_png(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotEnrichment")
    result = plots.plot_enrichment_counts(
        ["peak", "promoter"],
        ["a", "b"],
        [{"peak": 30, "promoter": 10}, {"peak": 5, "promoter": 45}],
        [100, 100],
        str(tmp_path / "enrich.png"),
        verbose=False,
    )
    _assert_png(result["path"])
    samples = result["data"]["samples"]
    assert samples[0]["features"]["peak"] == {"count": 30, "percent": 30.0}
    assert samples[1]["features"]["promoter"]["percent"] == 45.0
    assert result["data"]["bars_per_panel"] == 2


def test_plot_enrichment_counts_per_sample_layout_and_colors(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.plotEnrichment")
    result = plots.plot_enrichment_counts(
        ["peak", "promoter", "exon"],
        ["a"],
        [{"peak": 1, "promoter": 2, "exon": 3}],
        [10],
        str(tmp_path / "enrich.png"),
        per_sample=True,
        colors=["red", "#00ff00", "blue"],
        plot_title="title",
        verbose=False,
    )
    _assert_png(result["path"])
    assert result["data"]["bars_per_panel"] == 3
    assert result["data"]["samples"][0]["features"]["exon"]["percent"] == 30.0
    with pytest.raises(ConfigurationError, match="colors"):
        plots.plot_enrichment_counts(
            ["peak", "promoter"],
            ["a"],
            [{"peak": 1, "promoter": 2}],
            [10],
            str(tmp_path / "other.png"),
            per_sample=True,
            colors=["red"],
            verbose=False,
        )


@pytest.mark.parametrize(
    ("features", "labels", "counts", "totals"),
    [
        (["peak"], ["a"], [{"promoter": 1}], [10]),
        (["peak"], ["a", "b"], [{"peak": 1}], [10, 10]),
        (["peak"], ["a"], [{"peak": 1}], [0]),
        (["peak"], ["a"], [{"peak": -1}], [10]),
        (["peak"], ["a"], [[1]], [10]),
    ],
)
def test_plot_enrichment_counts_rejects_bad_input(
    tmp_path: Path, features: list[str], labels: list[str], counts: list, totals: list[int]
) -> None:
    with pytest.raises(ConfigurationError):
        plots.plot_enrichment_counts(
            features, labels, counts, totals, str(tmp_path / "e.png"), verbose=False
        )
    assert not (tmp_path / "e.png").exists()


# ---- template lengths --------------------------------------------------------------------


def test_template_lengths_from_bam_known_values(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.utilities")
    bam = _write_bam(tmp_path / "reads.bam")
    result = plots.template_lengths_from_bam(str(bam), verbose=False)
    # with TLEN -200 -> 200; no TLEN and 20M 100N 30M -> 50 (N is skipped); no TLEN 10M 5I 10M -> 20
    assert result["template_lengths"] == [200, 50, 20]
    assert result["alignments"] == 3
    assert result["truncated"] is False
    assert result["summary"] == {"min": 20, "median": 50, "max": 200}


def test_template_lengths_from_bam_signed_and_truncated(tmp_path: Path) -> None:
    pytest.importorskip("deeptools.utilities")
    bam = _write_bam(tmp_path / "reads.bam")
    signed = plots.template_lengths_from_bam(str(bam), not_abs=True, verbose=False)
    assert signed["template_lengths"] == [-200, 50, 20]
    capped = plots.template_lengths_from_bam(str(bam), max_reads=2, verbose=False)
    assert capped["alignments"] == 2
    assert capped["truncated"] is True
    assert capped["template_lengths"] == [200, 50]


def test_template_lengths_from_bam_rejects_non_alignment_file(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    text = tmp_path / "plain.txt"
    text.write_text("this is not an alignment file\n")
    with pytest.raises(ConfigurationError, match="could not be read as alignments"):
        plots.template_lengths_from_bam(str(text), verbose=False)
    with pytest.raises(ConfigurationError):
        plots.template_lengths_from_bam(str(tmp_path / "missing.bam"), verbose=False)


def test_template_lengths_missing_dependency_is_reported(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    bam = tmp_path / "reads.bam"
    bam.write_bytes(b"")
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(BackendUnavailableError, match="bio-deeptools"):
        plots.template_lengths_from_bam(str(bam), verbose=False)


# ---- blacklist subtraction ---------------------------------------------------------------


@pytest.fixture()
def blacklist(tmp_path: Path) -> Path:
    path = tmp_path / "blacklist.bed"
    path.write_text("chr1\t100\t200\nchr1\t300\t350\n")
    return path


def test_subtract_blacklist_intervals_known_values(blacklist: Path) -> None:
    pytest.importorskip("deeptoolsintervals")
    result = plots.subtract_blacklist_intervals(str(blacklist), "chr1", 0, 400, verbose=False)
    assert result["blacklisted_overlaps"] == [[100, 200], [300, 350]]
    assert result["remaining"] == [[0, 100], [200, 300], [350, 400]]
    assert result["remaining_bases"] == 250


def test_subtract_blacklist_intervals_edge_cases(blacklist: Path) -> None:
    pytest.importorskip("deeptoolsintervals")
    covered = plots.subtract_blacklist_intervals(str(blacklist), "chr1", 120, 180, verbose=False)
    assert covered["remaining"] == []
    assert covered["remaining_bases"] == 0
    other = plots.subtract_blacklist_intervals(str(blacklist), "chr2", 0, 400, verbose=False)
    assert other["blacklisted_overlaps"] == []
    assert other["remaining"] == [[0, 400]]


@pytest.mark.parametrize(("start", "end"), [(200, 200), (300, 100), (-1, 10)])
def test_subtract_blacklist_intervals_rejects_bad_coordinates(
    blacklist: Path, start: int, end: int
) -> None:
    with pytest.raises(ConfigurationError):
        plots.subtract_blacklist_intervals(str(blacklist), "chr1", start, end, verbose=False)


def test_subtract_blacklist_intervals_missing_dependency_is_reported(
    monkeypatch: pytest.MonkeyPatch, blacklist: Path
) -> None:
    monkeypatch.setitem(sys.modules, "deeptoolsintervals", None)
    with pytest.raises(BackendUnavailableError, match="bio-deeptools"):
        plots.subtract_blacklist_intervals(str(blacklist), "chr1", 0, 400, verbose=False)


# ---- module contract ---------------------------------------------------------------------


def test_public_functions_end_with_keyword_only_verbose() -> None:
    assert sorted(plots.__all__) == [
        "fragment_length_density",
        "plot_enrichment_counts",
        "plot_gc_bias",
        "plot_matrix_heatmap",
        "subtract_blacklist_intervals",
        "template_lengths_from_bam",
    ]
    for name in plots.__all__:
        parameters = list(inspect.signature(getattr(plots, name)).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name
