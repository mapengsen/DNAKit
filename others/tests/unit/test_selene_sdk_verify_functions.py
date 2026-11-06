"""Tests for selene_sdk_verify_functions: the sequence logo and the variant max-difference scatter.

The logo is checked against selene's own ``sequence_logo``: upstream returns its axes without rendering
(its failure is at draw time), so the bar heights and the glyph outlines of the two figures can be
compared directly. The scatter is checked by hand on a small matrix; the upstream scatter is asserted to
fail on plotly 6 (``titlefont``). Plots are written into tmp_path. Socket connections are blocked.
"""

from __future__ import annotations

import importlib.util
import socket

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.predictions import selene_sdk_verify_functions as module

requires_selene = pytest.mark.skipif(
    importlib.util.find_spec("selene_sdk") is None, reason="selene-sdk is not installed"
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

MATRICES = {
    "positive": [[0.5, 0.1, 0.2, 0.9], [0.3, 0.8, 0.1, 0.2], [0.6, 0.2, 0.7, 0.05]],
    "mixed": [
        [0.5, -0.2, 0.1, 0.0],
        [-0.3, 0.8, -0.1, 0.2],
        [0.1, 0.1, 0.9, -0.4],
        [-0.6, -0.5, 0.3, 0.3],
    ],
    "all_negative": [[-0.5, -0.2, -0.1, -0.9], [-0.3, -0.8, -0.1, -0.2]],
}


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)


@requires_selene
@pytest.mark.parametrize(
    "name,order",
    [("positive", "value"), ("mixed", "value"), ("mixed", "alpha"), ("all_negative", "value")],
)
def test_logo_bars_and_glyphs_match_selene(tmp_path, name, order):
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from selene_sdk.interpret import vis

    matrix = MATRICES[name]
    upstream_ax = vis.sequence_logo(np.array(matrix, dtype=float), order=order)
    try:
        ours = module.plot_sequence_logo_scores(
            matrix, str(tmp_path / "logo.png"), order=order, verbose=False
        )
        assert (tmp_path / "logo.png").read_bytes()[:8] == PNG_SIGNATURE
        length = len(matrix)
        bars = [patch for patch in upstream_ax.patches if isinstance(patch, Rectangle)][
            : 4 * length
        ]
        heights = np.array([bar.get_height() for bar in bars]).reshape(4, length)
        assert np.allclose(heights, np.array(ours["scores"]))
        glyphs = [patch for patch in upstream_ax.patches if not isinstance(patch, Rectangle)]
        assert len(glyphs) == 4 * length
        for index, patch in enumerate(glyphs):
            base = ours["bases"][index // length][index % length]
            expected = vis._SVG_PATHS[base][0]
            assert np.array_equal(np.asarray(patch.get_path().vertices), np.asarray(expected)), (
                index,
                base,
            )
    finally:
        plt.close("all")


@requires_selene
def test_logo_rejects_wrong_column_count(tmp_path):
    with pytest.raises(ConfigurationError, match="columns"):
        module.plot_sequence_logo_scores(
            [[0.1, 0.2, 0.3]], str(tmp_path / "bad.png"), verbose=False
        )


@requires_selene
def test_logo_refuses_an_existing_file(tmp_path):
    target = tmp_path / "exists.png"
    target.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        module.plot_sequence_logo_scores(MATRICES["positive"], str(target), verbose=False)


LABELS = [["chr2", 500, "v1", "A", "C"], ["chr1", 20, "v2", "G", "T"], ["chr1", 5, "v3", "C", "A"]]
SCORES = [[0.1, 0.5, 0.2], [0.9, 0.3, 0.1], [0.2, 0.05, 0.6]]
FEATURES = ["f1", "f2", "f3"]


@requires_selene
def test_scatter_orders_by_chromosome_and_takes_feature_maximum(tmp_path):
    target = tmp_path / "scatter.html"
    result = module.plot_variant_diffs_scatter(SCORES, LABELS, FEATURES, str(target), verbose=False)
    assert target.is_file()
    assert [(v["chrom"], v["pos"]) for v in result["variants"]] == [
        ("chr1", 5),
        ("chr1", 20),
        ("chr2", 500),
    ]
    assert [v["max_abs_diff"] for v in result["variants"]] == [0.6, 0.9, 0.5]


@requires_selene
def test_scatter_feature_filter_keeps_each_variant_with_its_own_scores(tmp_path):
    result = module.plot_variant_diffs_scatter(
        SCORES,
        LABELS,
        FEATURES,
        str(tmp_path / "filtered.html"),
        feature_indices=[1, 2],
        verbose=False,
    )
    assert result["features"] == ["f2", "f3"]
    expected = {("chr1", 5): 0.6, ("chr1", 20): 0.3, ("chr2", 500): 0.5}
    assert {(v["chrom"], v["pos"]): v["max_abs_diff"] for v in result["variants"]} == expected


@requires_selene
def test_scatter_percentile_and_threshold(tmp_path):
    result = module.plot_variant_diffs_scatter(
        SCORES,
        LABELS,
        FEATURES,
        str(tmp_path / "pct.html"),
        nth_percentile=40,
        threshold_line=0.4,
        verbose=False,
    )
    # maxima in label order are 0.6, 0.9, 0.5; the 40th percentile is 0.58, so variants at 5 and 20 remain
    assert [v["pos"] for v in result["variants"]] == [5, 20]
    text = (tmp_path / "pct.html").read_text()
    assert "Max probability difference scores" in text


@requires_selene
def test_upstream_scatter_fails_on_plotly_titlefont(tmp_path):
    import numpy as np
    from selene_sdk.interpret import vis

    with pytest.raises(ValueError, match="titlefont"):
        vis.variant_diffs_scatter_plot(
            np.array(SCORES), LABELS, FEATURES, str(tmp_path / "up.html")
        )
