"""Tests for dnakit.predictions.alphagenome_gap_functions.

Network calls are replaced with fakes on the alphagenome module attributes, and a socket guard makes
any accidental connection fail. Plot tests run alphagenome's own drawing code offline and check the
PNG files they write. Guard tests (missing alphagenome module) run without alphagenome installed.
"""

from __future__ import annotations

import importlib
import inspect
import json
import os
import socket
import struct
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.predictions import alphagenome_gap_functions as gap

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
LIVE = os.environ.get("DNAKIT_LIVE_TESTS") == "1"


@pytest.fixture(autouse=True)
def _block_network(monkeypatch, request):
    if LIVE and request.node.name.startswith("test_live_"):
        return

    original = socket.socket

    class BlockedSocket(original):
        """Subclassable so that third-party imports still work; creating a socket fails."""

        def __init__(self, *args, **kwargs):
            raise RuntimeError("network access is blocked in unit tests")

    monkeypatch.setattr(socket, "socket", BlockedSocket)


def _alphagenome():
    pytest.importorskip("alphagenome")
    pytest.importorskip("matplotlib")
    return importlib.import_module("alphagenome")


def _module(name: str):
    _alphagenome()
    return importlib.import_module(name)


def _png_size(path: Path) -> tuple[int, int]:
    data = path.read_bytes()
    assert data[:8] == PNG_SIGNATURE
    return struct.unpack(">II", data[16:24])


def _assert_png(result: dict, tmp_path: Path, name: str) -> dict:
    path = Path(result["path"])
    assert path.is_absolute()
    assert path == (tmp_path / name).resolve()
    assert path.exists()
    assert path.read_bytes()[:8] == PNG_SIGNATURE
    return result["data"]


# ---------------------------------------------------------------- guards (no alphagenome needed)


def test_missing_atlas_module_raises_backend_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, "alphagenome.atlas.atlas", None)
    with pytest.raises(BackendUnavailableError, match="bio-alphagenome"):
        gap.create_alphagenome_atlas_client("key", verbose=False)


def test_missing_dna_client_module_raises_backend_unavailable(monkeypatch):
    monkeypatch.setitem(sys.modules, "alphagenome.models.dna_client", None)
    with pytest.raises(BackendUnavailableError, match="bio-alphagenome"):
        gap.create_alphagenome_dna_client("key", verbose=False)


@pytest.mark.parametrize(
    "call",
    [
        lambda p: gap.get_alphagenome_junctions_to_plot([], [], ["T"], "T", "+", verbose=False),
        lambda p: gap.pad_alphagenome_track([[1.0]], 3, verbose=False),
        lambda p: gap.plot_alphagenome_seqlogo(
            [[1.0, 0.0, 0.0, 0.0]], str(p / "a.png"), verbose=False
        ),
        lambda p: gap.plot_alphagenome_contact_map(
            [[1.0, 2.0], [2.0, 1.0]], str(p / "a.png"), verbose=False
        ),
        lambda p: gap.plot_alphagenome_track([1.0, 2.0], str(p / "a.png"), verbose=False),
        lambda p: gap.plot_alphagenome_tracks({"a": [1.0, 2.0]}, str(p / "a.png"), verbose=False),
        lambda p: gap.plot_alphagenome_sashimi(
            [{"chromosome": "chr1", "start": 1, "end": 5, "strand": "+", "k": 2}],
            str(p / "a.png"),
            verbose=False,
        ),
        lambda p: gap.plot_alphagenome_components(
            [{"type": "empty"}], str(p / "a.png"), chromosome="chr1", start=0, end=1, verbose=False
        ),
        lambda p: gap.plot_alphagenome_transcripts(
            [{"transcript_id": "T", "strand": "+", "exons": [[1, 5]]}],
            str(p / "a.png"),
            chromosome="chr1",
            start=0,
            end=10,
            verbose=False,
        ),
    ],
)
def test_missing_alphagenome_package_guards_every_data_and_plot_wrapper(
    monkeypatch, tmp_path, call
):
    for module_name in (
        "alphagenome.visualization.plot",
        "alphagenome.visualization.plot_components",
        "alphagenome.visualization.plot_transcripts",
        "alphagenome.data.junction_data",
        "alphagenome.data.genome",
    ):
        monkeypatch.setitem(sys.modules, module_name, None)
    with pytest.raises(BackendUnavailableError) as info:
        call(tmp_path)
    assert "bio-alphagenome" in str(info.value)


def test_verbose_is_last_keyword_only_parameter():
    for name in gap.__all__:
        parameters = list(inspect.signature(getattr(gap, name)).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters[-1].default is True, name


# ---------------------------------------------------------------- network wrappers (fakes)


class _FakeClient:
    pass


def test_atlas_client_calls_upstream_create_and_never_returns_the_key(monkeypatch):
    atlas = _module("alphagenome.atlas.atlas")
    calls = []

    def fake_create(api_key, *, timeout=None, address=None):
        calls.append((api_key, timeout, address))
        return _FakeClient()

    monkeypatch.setattr(atlas, "create", fake_create)
    result = gap.create_alphagenome_atlas_client(
        "SECRET-KEY-123", timeout=5, address="localhost:5005", verbose=False
    )
    assert calls == [("SECRET-KEY-123", 5.0, "localhost:5005")]
    assert result == {
        "service": "alphagenome.atlas",
        "source": "localhost:5005",
        "result": {
            "client_type": "_FakeClient",
            "address": "localhost:5005",
            "timeout_seconds": 5.0,
        },
    }
    assert "SECRET-KEY" not in json.dumps(result)


def test_atlas_client_default_address_is_the_hosted_service(monkeypatch):
    atlas = _module("alphagenome.atlas.atlas")
    monkeypatch.setattr(atlas, "create", lambda api_key, **kwargs: _FakeClient())
    result = gap.create_alphagenome_atlas_client("key", verbose=False)
    assert result["source"] == "dns:///gdmscience.googleapis.com:443"
    assert result["result"]["timeout_seconds"] is None


def test_atlas_client_rejects_empty_key_without_calling_upstream(monkeypatch):
    atlas = _module("alphagenome.atlas.atlas")

    def fail(*args, **kwargs):
        raise AssertionError("upstream create must not run for an empty key")

    monkeypatch.setattr(atlas, "create", fail)
    with pytest.raises(ConfigurationError):
        gap.create_alphagenome_atlas_client("  ", verbose=False)


def test_dna_client_passes_requested_model_version(monkeypatch):
    _alphagenome()
    dna_client = importlib.import_module("alphagenome.models.dna_client")
    dna_model = importlib.import_module("alphagenome.models.v1.dna_model")
    seen = {}

    def fake_create(api_key, *, model_version=None, timeout=None, address=None):
        seen.update(api_key=api_key, model_version=model_version, timeout=timeout, address=address)
        return _FakeClient()

    monkeypatch.setattr(dna_client, "create", fake_create)
    result = gap.create_alphagenome_dna_client(
        "KEY", model_version="fold_2", timeout=3.5, verbose=False
    )
    assert seen == {
        "api_key": "KEY",
        "model_version": dna_model.ModelVersion.FOLD_2,
        "timeout": 3.5,
        "address": None,
    }
    assert result["service"] == "alphagenome.models.dna_client"
    assert result["result"]["model_version"] == "fold_2"
    assert result["result"]["model_version_source"] == "requested"
    assert "KEY" not in json.dumps(result)


def test_dna_client_default_model_version_is_all_folds(monkeypatch):
    _alphagenome()
    dna_client = importlib.import_module("alphagenome.models.dna_client")
    seen = {}

    def fake_create(api_key, *, model_version=None, timeout=None, address=None):
        seen["model_version"] = model_version
        return _FakeClient()

    monkeypatch.setattr(dna_client, "create", fake_create)
    result = gap.create_alphagenome_dna_client("KEY", verbose=False)
    assert seen["model_version"] is None
    assert result["result"]["model_version"] == "all_folds"
    assert result["result"]["model_version_source"] == "upstream default"


def test_dna_client_rejects_unknown_model_version_before_connecting(monkeypatch):
    _alphagenome()
    dna_client = importlib.import_module("alphagenome.models.dna_client")

    def fail(*args, **kwargs):
        raise AssertionError("upstream create must not run for an invalid model version")

    monkeypatch.setattr(dna_client, "create", fail)
    with pytest.raises(ConfigurationError):
        gap.create_alphagenome_dna_client("KEY", model_version="fold_9", verbose=False)


@pytest.mark.skipif(
    not LIVE, reason="live test: set DNAKIT_LIVE_TESTS=1 to contact the AlphaGenome service"
)
@pytest.mark.skipif(
    not os.environ.get("ALPHAGENOME_API_KEY"), reason="live test needs ALPHAGENOME_API_KEY"
)
def test_live_dna_client_connects_with_api_key():
    _alphagenome()
    key = os.environ["ALPHAGENOME_API_KEY"]
    result = gap.create_alphagenome_dna_client(key, timeout=60, verbose=False)
    assert result["service"] == "alphagenome.models.dna_client"
    assert result["result"]["model_version"] == "all_folds"
    assert key not in json.dumps(result)


# ---------------------------------------------------------------- junction data


JUNCTIONS = [
    {"chromosome": "chr1", "start": 100, "end": 200, "strand": "+"},
    {"chromosome": "chr1", "start": 150, "end": 300, "strand": "+"},
    {"chromosome": "chr1", "start": 150, "end": 300, "strand": "-"},
]


def test_junction_selection_applies_strand_and_threshold():
    _alphagenome()
    result = gap.get_alphagenome_junctions_to_plot(
        JUNCTIONS, [[5.0], [0.5], [3.0]], ["T"], "T", "+", k_threshold=1.0, verbose=False
    )
    assert result["count"] == 1
    assert result["junctions"] == [
        {"chromosome": "chr1", "start": 100, "end": 200, "strand": "+", "k": 5.0}
    ]


def test_junction_selection_with_zero_threshold_keeps_all_strand_junctions():
    _alphagenome()
    result = gap.get_alphagenome_junctions_to_plot(
        JUNCTIONS, [[5.0], [0.5], [3.0]], ["T"], "T", "-", k_threshold=0.0, verbose=False
    )
    assert result["junctions"] == [
        {"chromosome": "chr1", "start": 150, "end": 300, "strand": "-", "k": 3.0}
    ]


def test_junction_selection_none_threshold_uses_five_percent_of_maximum():
    _alphagenome()
    # Maximum over all junctions of the track is 5.0, so the cut-off is 0.25.
    result = gap.get_alphagenome_junctions_to_plot(
        JUNCTIONS, [[5.0], [0.5], [3.0]], ["T"], "T", "+", k_threshold=None, verbose=False
    )
    assert result["k_threshold"] is None
    assert [item["k"] for item in result["junctions"]] == [5.0, 0.5]


def test_junction_selection_rounds_counts_to_two_decimals_and_reports_unknown_track():
    _alphagenome()
    rounded = gap.get_alphagenome_junctions_to_plot(
        JUNCTIONS[:1], [[2.456]], ["T"], "T", "+", verbose=False
    )
    assert rounded["junctions"][0]["k"] == 2.46
    missing = gap.get_alphagenome_junctions_to_plot(
        JUNCTIONS, [[5.0], [0.5], [3.0]], ["T"], "other", "+", verbose=False
    )
    assert missing["count"] == 0


def test_junction_selection_rejects_duplicate_track_names_and_bad_shapes():
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.get_alphagenome_junctions_to_plot(
            JUNCTIONS[:1], [[1.0, 2.0]], ["T", "T"], "T", "+", verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.get_alphagenome_junctions_to_plot(JUNCTIONS, [[1.0]], ["T"], "T", "+", verbose=False)
    with pytest.raises(ConfigurationError):
        gap.get_alphagenome_junctions_to_plot(
            [{"chromosome": "chr1", "start": 9, "end": 3, "strand": "+"}],
            [[1.0]],
            ["T"],
            "T",
            "+",
            verbose=False,
        )


# ---------------------------------------------------------------- plots


def test_seqlogo_writes_png_of_expected_size_and_reports_heights(tmp_path):
    _alphagenome()
    heights = [[1.0, 0.0, 0.0, 0.0], [0.0, 2.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.5]]
    result = gap.plot_alphagenome_seqlogo(
        heights, str(tmp_path / "nested" / "logo.png"), verbose=False
    )
    data = _assert_png(result, tmp_path / "nested", "logo.png")
    assert _png_size(Path(result["path"])) == (
        600,
        375,
    )  # max(4.0, 0.25 * 3) in x 2.5 in at 150 dpi
    assert data == {
        "letter_heights": heights,
        "alphabet": "ACGT",
        "one_based": True,
        "start": 0,
        "letter_colors": "default",
        "positions": 3,
    }


def test_seqlogo_rejects_alphabet_mismatch_and_non_png_path(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_seqlogo(
            [[1.0, 0.0]], str(tmp_path / "a.png"), alphabet="ACGT", verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_seqlogo([[1.0, 0.0, 0.0, 0.0]], str(tmp_path / "a.pdf"), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_seqlogo(
            [[1.0, 0.0, 0.0, 0.0]], str(tmp_path / "a.png"), letter_colors="rainbow", verbose=False
        )


def test_contact_map_reports_scale_taken_from_positive_values(tmp_path):
    _alphagenome()
    matrix = [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.5]]
    result = gap.plot_alphagenome_contact_map(matrix, str(tmp_path / "cm.png"), verbose=False)
    data = _assert_png(result, tmp_path, "cm.png")
    assert data["vmin"] == 0.5
    assert data["vmax"] == 3.0
    assert data["contact_map"] == matrix


def test_contact_map_rejects_non_square_and_bad_scales(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_contact_map([[1.0, 2.0, 3.0]], str(tmp_path / "a.png"), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_contact_map(
            [[0.0, 0.0], [0.0, 0.0]], str(tmp_path / "a.png"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_contact_map(
            [[1.0, 2.0], [2.0, 1.0]], str(tmp_path / "a.png"), vmin=3.0, vmax=2.0, verbose=False
        )


def test_track_rendering_follows_alphagenome_shape_rules(tmp_path):
    _alphagenome()
    line = gap.plot_alphagenome_track(
        [1.0, 2.0, 3.0, 2.0], str(tmp_path / "line.png"), verbose=False
    )
    assert line["data"]["rendering"] == "line"
    assert line["data"]["x"] is None
    logo = gap.plot_alphagenome_track(
        [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]], str(tmp_path / "logo.png"), verbose=False
    )
    assert logo["data"]["rendering"] == "sequence_logo"
    strands = gap.plot_alphagenome_track(
        [[1.0, -1.0], [2.0, -2.0], [3.0, -3.0]],
        str(tmp_path / "strands.png"),
        legend=True,
        ylim=[-4.0, 4.0],
        color=["red", "blue"],
        verbose=False,
    )
    assert strands["data"]["rendering"] == "two_strand_lines"
    assert strands["data"]["ylim"] == [-4.0, 4.0]
    for name in ("line.png", "logo.png", "strands.png"):
        assert (tmp_path / name).read_bytes()[:8] == PNG_SIGNATURE


def test_track_rejects_unsupported_channels_and_length_mismatch(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_track([[1.0, 2.0, 3.0]] * 2, str(tmp_path / "a.png"), verbose=False)
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_track([1.0, 2.0], str(tmp_path / "a.png"), x=[1.0], verbose=False)
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_track([1.0, float("nan")], str(tmp_path / "a.png"), verbose=False)


def test_tracks_plot_uses_shared_scale_and_per_track_heights(tmp_path):
    _alphagenome()
    result = gap.plot_alphagenome_tracks(
        {"a": [1.0, 2.0, 3.0, 4.0], "b": [4.0, 3.0, 2.0, 1.0]},
        str(tmp_path / "tracks.png"),
        fig_track_height={"a": 1.0, "b": 2.0},
        ylim="same",
        color={"a": "red", "b": "blue"},
        filled_tracks=["a"],
        verbose=False,
    )
    data = _assert_png(result, tmp_path, "tracks.png")
    assert data["track_names"] == ["a", "b"]
    assert data["tracks"] == {"a": [1.0, 2.0, 3.0, 4.0], "b": [4.0, 3.0, 2.0, 1.0]}
    assert data["filled_tracks"] == ["a"]
    assert _png_size(Path(result["path"])) == (3000, 450)  # 20 in x (1 + 2) in at 150 dpi


def test_tracks_plot_rejects_unequal_lengths_and_bad_heights(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_tracks(
            {"a": [1.0, 2.0], "b": [1.0]}, str(tmp_path / "t.png"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_tracks(
            {"a": [1.0, 2.0]}, str(tmp_path / "t.png"), fig_track_height={"b": 1.0}, verbose=False
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_tracks(
            {"a": [1.0, 2.0]}, str(tmp_path / "t.png"), ylim="fast", verbose=False
        )


SASHIMI_JUNCTIONS = [
    {"chromosome": "chr1", "start": 100, "end": 200, "strand": "+", "k": 5},
    {"chromosome": "chr1", "start": 150, "end": 300, "strand": "+", "k": 0.5},
    {"chromosome": "chr1", "start": 150, "end": 300, "strand": "-", "k": 3},
]


def test_sashimi_reports_filtered_junctions_and_is_reproducible_with_a_seed(tmp_path):
    _alphagenome()
    first = gap.plot_alphagenome_sashimi(
        SASHIMI_JUNCTIONS,
        str(tmp_path / "one.png"),
        interval={"chromosome": "chr1", "start": 0, "end": 1000},
        verbose=False,
    )
    second = gap.plot_alphagenome_sashimi(
        SASHIMI_JUNCTIONS,
        str(tmp_path / "two.png"),
        interval={"chromosome": "chr1", "start": 0, "end": 1000},
        verbose=False,
    )
    data = _assert_png(first, tmp_path, "one.png")
    assert data["total_k"] == 8.5
    assert data["filter_threshold"] == 0.01
    assert data["seed"] == 0
    assert data["drawn_count"] == 3
    assert [item["k"] for item in data["drawn"]] == [5.0, 0.5, 3.0]
    assert (tmp_path / "one.png").read_bytes() == (tmp_path / "two.png").read_bytes()
    assert second["data"]["drawn_count"] == 3


def test_sashimi_threshold_removes_small_junctions_from_the_drawn_list(tmp_path):
    _alphagenome()
    result = gap.plot_alphagenome_sashimi(
        SASHIMI_JUNCTIONS, str(tmp_path / "s.png"), filter_threshold=1.0, verbose=False
    )
    assert [item["k"] for item in result["data"]["drawn"]] == [5.0, 3.0]
    assert len(result["data"]["junctions"]) == 3


def test_sashimi_rejects_negative_counts_and_all_zero_counts(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_sashimi(
            [{"chromosome": "chr1", "start": 1, "end": 5, "strand": "+", "k": -1}],
            str(tmp_path / "s.png"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_sashimi(
            [{"chromosome": "chr1", "start": 1, "end": 5, "strand": "+", "k": 0}],
            str(tmp_path / "s.png"),
            verbose=False,
        )


def test_pad_track_centres_even_and_odd_padding():
    _alphagenome()
    track = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
    even = gap.pad_alphagenome_track(track, 7, verbose=False)
    assert even["offset"] == 2
    assert even["padded"] == [
        [0.0, 0.0],
        [0.0, 0.0],
        [1.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
    ]
    odd = gap.pad_alphagenome_track(track, 6, value=-1.0, verbose=False)
    assert odd["offset"] == 1
    assert odd["padded"] == [
        [-1.0, -1.0],
        [1.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0],
        [-1.0, -1.0],
        [-1.0, -1.0],
    ]
    assert odd["length"] == 3 and odd["new_length"] == 6


def test_pad_track_rejects_shrinking_and_one_dimensional_input():
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.pad_alphagenome_track([[1.0], [2.0], [3.0]], 2, verbose=False)
    with pytest.raises(ConfigurationError):
        gap.pad_alphagenome_track([1.0, 2.0], 4, verbose=False)


def test_components_plot_tracks_logo_and_spacer(tmp_path):
    _alphagenome()
    components = [
        {
            "type": "tracks",
            "names": ["ctcf", "atac"],
            "values": [[float(i), i + 0.5] for i in range(10)],
            "strands": ["+", "-"],
        },
        {"type": "seq_logo", "scores": [[1.0, 0.0, 0.0, 0.0]] * 10},
        {"type": "empty", "fig_height": 0.5},
    ]
    result = gap.plot_alphagenome_components(
        components,
        str(tmp_path / "comp.png"),
        chromosome="chr1",
        start=1000,
        end=1010,
        verbose=False,
    )
    data = _assert_png(result, tmp_path, "comp.png")
    assert data["interval"] == {"chromosome": "chr1", "start": 1000, "end": 1010}
    assert data["components"][0]["names"] == ["ctcf", "atac"]
    assert data["components"][0]["strands"] == ["+", "-"]
    assert data["components"][1]["alphabet"] == "ACGT"
    assert data["components"][2] == {"type": "empty", "fig_height": 0.5}


def test_components_reject_wrong_width_bad_keys_and_duplicate_names(tmp_path):
    _alphagenome()
    base = {"chromosome": "chr1", "start": 0, "end": 4}
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_components(
            [{"type": "tracks", "names": ["a"], "values": [[1.0]] * 3}],
            str(tmp_path / "c.png"),
            verbose=False,
            **base,
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_components(
            [{"type": "empty", "fig_height": 1.0, "colour": "red"}],
            str(tmp_path / "c.png"),
            verbose=False,
            **base,
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_components(
            [{"type": "transcripts"}], str(tmp_path / "c.png"), verbose=False, **base
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_components(
            [{"type": "tracks", "names": ["a", "a"], "values": [[1.0, 2.0]] * 4}],
            str(tmp_path / "c.png"),
            verbose=False,
            **base,
        )


TRANSCRIPTS = [
    {
        "transcript_id": "T1",
        "strand": "+",
        "exons": [[1000, 1040], [1060, 1100]],
        "cds": [[1010, 1040]],
        "info": {"gene_name": "GENE1"},
    },
    {"transcript_id": "T2", "strand": "-", "exons": [[1020, 1090]], "info": {"gene_name": "GENE2"}},
]


def test_transcripts_plot_draws_both_coding_and_non_coding_models(tmp_path):
    _alphagenome()
    result = gap.plot_alphagenome_transcripts(
        TRANSCRIPTS,
        str(tmp_path / "tx.png"),
        chromosome="chr1",
        start=1000,
        end=1100,
        label_name="gene_name",
        verbose=False,
    )
    data = _assert_png(result, tmp_path, "tx.png")
    assert data["transcripts"][0]["exons"] == [[1000, 1040], [1060, 1100]]
    assert data["transcripts"][0]["cds"] == [[1010, 1040]]
    assert data["transcripts"][1]["cds"] is None
    assert data["label_name"] == "gene_name"


def test_transcripts_plot_rejects_missing_label_bad_strand_and_empty_list(tmp_path):
    _alphagenome()
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_transcripts(
            [{"transcript_id": "T", "strand": "+", "exons": [[1, 5]]}],
            str(tmp_path / "tx.png"),
            chromosome="chr1",
            start=0,
            end=10,
            label_name="gene_name",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_transcripts(
            [{"transcript_id": "T", "strand": ".", "exons": [[1, 5]]}],
            str(tmp_path / "tx.png"),
            chromosome="chr1",
            start=0,
            end=10,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        gap.plot_alphagenome_transcripts(
            [], str(tmp_path / "tx.png"), chromosome="chr1", start=0, end=10, verbose=False
        )
