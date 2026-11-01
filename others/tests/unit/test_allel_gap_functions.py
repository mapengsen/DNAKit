"""Tests for dnakit.population.allel_gap_functions (scikit-allel plots and the ROH wrapper).

Each plot is drawn from a small known input into tmp_path and the returned numbers are checked
against values worked out for that input. Optional packages are stubbed in sys.modules to test the
guard paths; nothing touches the network.
"""

from __future__ import annotations

import importlib.util
import inspect
import json
import socket
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.population import allel_gap_functions as gap

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

HAPLOTYPES = [
    [0, 0, 1, 1, 0],
    [0, 0, 1, 1, 1],
    [1, 0, 0, 1, 1],
    [1, 0, 0, 0, 1],
    [1, 1, 0, 0, 0],
    [0, 1, 0, 0, 0],
]
POSITIONS = [1, 5, 9, 13, 20, 28]


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network access is disabled in unit tests")

    monkeypatch.setattr(socket, "socket", _blocked)


def _allel():
    return pytest.importorskip("allel")


def _assert_png(result: dict, expected_path: Path) -> None:
    path = Path(result["path"])
    assert path == expected_path.resolve()
    assert path.is_file()
    assert path.read_bytes()[:8] == PNG_SIGNATURE
    json.dumps(result)


def test_every_public_function_takes_keyword_only_verbose() -> None:
    assert gap.__all__, "the module must export its functions"
    for name in gap.__all__:
        signature = inspect.signature(getattr(gap, name))
        parameters = list(signature.parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name


def test_plot_sfs_draws_the_clipped_spectrum(tmp_path: Path) -> None:
    _allel()
    target = tmp_path / "sfs.png"
    result = gap.plot_sfs([5, 3, 2, 1, 4], str(target), verbose=False)
    _assert_png(result, target)
    data = result["data"]
    assert data["spectrum"] == [5, 3, 2, 1, 4]
    # scikit-allel drops the first and last entries by default.
    assert data["x"] == [1, 2, 3]
    assert data["y"] == [3, 2, 1]
    assert data["xlabel"] == "derived allele count"
    assert data["ylabel"] == "site frequency"
    assert data["yscale"] == "log"


def test_plot_sfs_frequency_axis_and_bins(tmp_path: Path) -> None:
    _allel()
    target = tmp_path / "sfs_binned.png"
    result = gap.plot_sfs([5, 3, 2, 1, 4], str(target), n=10, bins=2)
    _assert_png(result, target)
    data = result["data"]
    # Bins over allele counts 1..3 hold [3] and [2, 1]; x is the bin midpoint divided by n.
    assert data["y"] == pytest.approx([3.0, 3.0])
    assert data["x"] == pytest.approx([0.15, 0.25])
    assert data["xlabel"] == "derived allele frequency"


def test_plot_sfs_folded_and_scaled_labels(tmp_path: Path) -> None:
    _allel()
    spectrum = [5, 3, 2, 1, 4]
    folded = gap.plot_sfs_folded(spectrum, str(tmp_path / "folded.png"))
    assert folded["data"]["xlabel"] == "minor allele count"
    assert folded["data"]["y"] == [3, 2, 1]

    scaled = gap.plot_sfs_scaled(spectrum, str(tmp_path / "scaled.png"))
    assert scaled["data"]["ylabel"] == "scaled site frequency"
    # scikit-allel draws the values it is given; the wrapper reports them unchanged.
    assert scaled["data"]["y"] == [3, 2, 1]

    both = gap.plot_sfs_folded_scaled(spectrum, str(tmp_path / "both.png"))
    assert both["data"]["xlabel"] == "minor allele count"
    assert both["data"]["ylabel"] == "scaled site frequency"


def test_plot_joint_sfs_variants_draw_the_transposed_matrix(tmp_path: Path) -> None:
    _allel()
    spectrum = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    transposed = [[1, 4, 7], [2, 5, 8], [3, 6, 9]]
    plain = gap.plot_joint_sfs(spectrum, str(tmp_path / "joint.png"))
    _assert_png(plain, tmp_path / "joint.png")
    assert plain["data"]["spectrum"] == spectrum
    assert plain["data"]["shape"] == [3, 3]
    assert plain["data"]["drawn_matrix"] == transposed
    assert plain["data"]["xlabel"] == "derived allele count (population 1)"

    folded = gap.plot_joint_sfs_folded(spectrum, str(tmp_path / "joint_folded.png"))
    assert folded["data"]["drawn_matrix"] == transposed
    assert folded["data"]["xlabel"] == "minor allele count (population 1)"

    scaled = gap.plot_joint_sfs_scaled(spectrum, str(tmp_path / "joint_scaled.png"))
    assert scaled["data"]["drawn_matrix"] == transposed

    both = gap.plot_joint_sfs_folded_scaled(spectrum, str(tmp_path / "joint_both.png"))
    assert both["data"]["drawn_matrix"] == transposed
    assert both["data"]["ylabel"] == "minor allele count (population 2)"


def test_fig_voight_painting_cuts_flanks_at_the_middle_variant(tmp_path: Path) -> None:
    allel = _allel()
    pytest.importorskip("seaborn")
    target = tmp_path / "voight.png"
    result = gap.fig_voight_painting(HAPLOTYPES, str(target))
    _assert_png(result, target)
    data = result["data"]
    assert data["index"] == 3
    assert data["n_variants"] == 6
    assert data["n_haplotypes"] == 5
    left, _ = allel.voight_painting(HAPLOTYPES[:4][::-1])
    right, _ = allel.voight_painting(HAPLOTYPES[3:])
    assert data["left_painting"] == left.tolist()
    assert data["right_painting"] == right.tolist()


def test_fig_voight_painting_rejects_bad_index_and_multiallelic_input(tmp_path: Path) -> None:
    _allel()
    pytest.importorskip("seaborn")
    with pytest.raises(ConfigurationError):
        gap.fig_voight_painting(HAPLOTYPES, str(tmp_path / "bad_index.png"), index=6)
    multiallelic = [[0, 2], [1, 0], [0, 1]]
    with pytest.raises(ConfigurationError, match="could not draw"):
        gap.fig_voight_painting(multiallelic, str(tmp_path / "multi.png"))
    assert not (tmp_path / "bad_index.png").exists()


def test_plot_voight_painting_reports_colours_and_flank(tmp_path: Path) -> None:
    _allel()
    pytest.importorskip("seaborn")
    painting = [[1, 1, 1, 2, 2], [1, 1, 0, 2, 2], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    target = tmp_path / "painting.png"
    result = gap.plot_voight_painting(painting, str(target), flank="left")
    _assert_png(result, target)
    assert result["data"]["painting"] == painting
    assert result["data"]["flank"] == "left"
    # Group ids 1 and 2 are the shared prefixes, so two colours after the white singletons.
    assert result["data"]["n_colors"] == 2
    with pytest.raises(ConfigurationError):
        gap.plot_voight_painting(painting, str(tmp_path / "flank.png"), flank="middle")


def test_plot_haplotype_frequencies_counts_distinct_haplotypes(tmp_path: Path) -> None:
    _allel()
    pytest.importorskip("seaborn")
    # Columns 0, 1 and 3 are the same haplotype; columns 2 and 4 are singletons.
    haplotypes = [
        [0, 0, 1, 0, 0],
        [1, 1, 0, 1, 0],
        [0, 0, 1, 0, 1],
        [1, 1, 0, 1, 1],
    ]
    target = tmp_path / "haplotypes.png"
    result = gap.plot_haplotype_frequencies(haplotypes, str(target))
    _assert_png(result, target)
    data = result["data"]
    assert sorted(data["haplotype_counts"], reverse=True) == [3, 1, 1]
    assert sum(data["haplotype_counts"]) == 5
    assert data["singleton_count"] == 2
    assert data["n_variants"] == 4
    assert data["n_haplotypes"] == 5


def test_plot_moving_haplotype_frequencies_window_bounds_and_ranks(tmp_path: Path) -> None:
    _allel()
    pytest.importorskip("seaborn")
    target = tmp_path / "moving.png"
    result = gap.plot_moving_haplotype_frequencies(POSITIONS, HAPLOTYPES, 3, str(target))
    _assert_png(result, target)
    data = result["data"]
    assert data["window_size"] == 3
    assert data["window_bounds"] == [[1, 9], [13, 28]]
    # Window 0 has no repeated haplotype. In window 1 the repeated pair is rank 1 and fills
    # the first two slots of the frequency-sorted order.
    assert data["rank_matrix"] == [[0, 0, 0, 0, 0], [1, 1, 0, 0, 0]]


def test_plot_moving_haplotype_frequencies_checks_positions(tmp_path: Path) -> None:
    _allel()
    pytest.importorskip("seaborn")
    with pytest.raises(ConfigurationError, match="one entry per variant"):
        gap.plot_moving_haplotype_frequencies(POSITIONS[:4], HAPLOTYPES, 3, str(tmp_path / "x.png"))


def test_plot_variant_locator_default_step_for_large_and_small_contigs(tmp_path: Path) -> None:
    _allel()
    large = list(range(1, 451, 3))
    result = gap.plot_variant_locator(large, str(tmp_path / "large.png"))
    _assert_png(result, tmp_path / "large.png")
    assert result["data"]["n_variants"] == 150
    assert result["data"]["step"] == 1
    assert result["data"]["lines_drawn"] == 150

    # scikit-allel raises "slice step cannot be zero" here; the wrapper never uses a zero step.
    small = gap.plot_variant_locator([10, 20, 30, 40, 50], str(tmp_path / "small.png"))
    _assert_png(small, tmp_path / "small.png")
    assert small["data"]["step"] == 1
    assert small["data"]["lines_drawn"] == 5
    assert small["data"]["start"] == pytest.approx(10.0)
    assert small["data"]["stop"] == pytest.approx(50.0)


def test_plot_variant_locator_explicit_step_and_order_check(tmp_path: Path) -> None:
    _allel()
    result = gap.plot_variant_locator([10, 20, 30, 40, 50], str(tmp_path / "step.png"), step=2)
    assert result["data"]["step"] == 2
    assert result["data"]["lines_drawn"] == 3
    with pytest.raises(ConfigurationError, match="sorted"):
        gap.plot_variant_locator([30, 10, 20], str(tmp_path / "unsorted.png"))


def test_plot_pairwise_distance_expands_condensed_form(tmp_path: Path) -> None:
    _allel()
    target = tmp_path / "distance.png"
    result = gap.plot_pairwise_distance(
        [1.0, 2.0, 3.0], str(target), labels=["a", "b", "c"], verbose=False
    )
    _assert_png(result, target)
    data = result["data"]
    assert data["n_samples"] == 3
    assert data["labels"] == ["a", "b", "c"]
    assert data["square_matrix"] == [[0, 1, 2], [1, 0, 3], [2, 3, 0]]

    square = gap.plot_pairwise_distance(
        [[0, 1, 2], [1, 0, 3], [2, 3, 0]], str(tmp_path / "square.png")
    )
    assert square["data"]["square_matrix"] == [[0, 1, 2], [1, 0, 3], [2, 3, 0]]


def test_plot_pairwise_distance_rejects_bad_shapes(tmp_path: Path) -> None:
    _allel()
    with pytest.raises(ConfigurationError, match="n\\*\\(n-1\\)/2"):
        gap.plot_pairwise_distance([1.0, 2.0], str(tmp_path / "bad.png"))
    with pytest.raises(ConfigurationError, match="names"):
        gap.plot_pairwise_distance([1.0, 2.0, 3.0], str(tmp_path / "labels.png"), labels=["a"])


def test_plot_pairwise_ld_draws_lower_triangle_flipped(tmp_path: Path) -> None:
    _allel()
    target = tmp_path / "ld.png"
    result = gap.plot_pairwise_ld([0.5, 0.2, 0.8], str(target))
    _assert_png(result, target)
    data = result["data"]
    assert data["n_variants"] == 3
    assert data["ld_values"] == [0.5, 0.2, 0.8]
    expected = [[0.2, 0.8, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.0]]
    drawn = data["drawn_matrix"]
    assert len(drawn) == 3
    for row, expected_row in zip(drawn, expected, strict=True):
        assert row == pytest.approx(expected_row)


def test_output_path_rules(tmp_path: Path) -> None:
    _allel()
    with pytest.raises(ConfigurationError, match="\\.png"):
        gap.plot_sfs([1, 2, 3, 4], str(tmp_path / "plot.pdf"))
    existing = tmp_path / "exists.png"
    existing.write_bytes(b"keep")
    with pytest.raises(ConfigurationError, match="not overwritten"):
        gap.plot_sfs([1, 2, 3, 4], str(existing))
    assert existing.read_bytes() == b"keep"
    nested = tmp_path / "new" / "dir" / "sfs.png"
    result = gap.plot_sfs([1, 2, 3, 4], str(nested))
    assert Path(result["path"]).is_file()


def test_input_validation(tmp_path: Path) -> None:
    _allel()
    with pytest.raises(ConfigurationError):
        gap.plot_sfs([5, -3, 2, 1], str(tmp_path / "negative.png"))
    with pytest.raises(ConfigurationError):
        gap.plot_joint_sfs([[1, 2], [3]], str(tmp_path / "ragged.png"))
    with pytest.raises(ConfigurationError):
        gap.plot_sfs([1, 2, 3, 4], str(tmp_path / "yscale.png"), yscale="")
    with pytest.raises(ConfigurationError):
        gap.plot_haplotype_frequencies([[0, 1.5], [1, 0]], str(tmp_path / "fraction.png"))
    with pytest.raises(ConfigurationError):
        gap.plot_moving_haplotype_frequencies(POSITIONS, HAPLOTYPES, 0, str(tmp_path / "w.png"))


def test_missing_allel_gives_backend_unavailable(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setitem(sys.modules, "allel", None)
    with pytest.raises(BackendUnavailableError) as raised:
        gap.plot_sfs([1, 2, 3, 4], str(tmp_path / "no_allel.png"))
    assert "bio-scikit-allel" in (raised.value.hint or "")
    assert not (tmp_path / "no_allel.png").exists()


def test_missing_seaborn_gives_backend_unavailable(monkeypatch, tmp_path: Path) -> None:
    _allel()
    monkeypatch.setitem(sys.modules, "seaborn", None)
    with pytest.raises(BackendUnavailableError) as raised:
        gap.plot_haplotype_frequencies(HAPLOTYPES, str(tmp_path / "no_seaborn.png"))
    assert "seaborn" in raised.value.message
    assert "bio-scikit-allel" in (raised.value.hint or "")
    with pytest.raises(BackendUnavailableError):
        gap.plot_voight_painting([[1, 0]], str(tmp_path / "no_seaborn2.png"))


def test_roh_guard_without_protopunica(monkeypatch) -> None:
    _allel()
    monkeypatch.setitem(sys.modules, "protopunica", None)
    with pytest.raises(BackendUnavailableError) as raised:
        gap.roh_poissonhmm(
            [[0, 0], [0, 0], [0, 1]], [100, 200, 300], contig_size=1000, verbose=False
        )
    assert "protopunica" in raised.value.message
    assert "bio-scikit-allel" in (raised.value.hint or "")


def test_roh_input_rules(monkeypatch) -> None:
    _allel()
    monkeypatch.setitem(sys.modules, "protopunica", None)
    with pytest.raises(ConfigurationError, match="contig_size"):
        gap.roh_poissonhmm([[0, 0], [0, 1]], [100, 200])
    with pytest.raises(ConfigurationError, match="one entry per genotype row"):
        gap.roh_poissonhmm([[0, 0], [0, 1]], [100], contig_size=1000)


@pytest.mark.skipif(
    importlib.util.find_spec("protopunica") is None,
    reason="optional package protopunica is not installed; the Poisson HMM ROH caller cannot run",
)
def test_roh_poissonhmm_homozygous_contig() -> None:
    _allel()
    genotypes = [[0, 0] for _ in range(100)]
    positions = [10 * i + 1 for i in range(100)]
    result = gap.roh_poissonhmm(genotypes, positions, contig_size=1000)
    assert result["froh"] > 0.5
    assert result["block_count"] >= 1
