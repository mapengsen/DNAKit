"""Tests for the P19 Bio.Graphics figures (PDF and SVG; PNG needs rlPyCairo and is skipped without it)."""

from __future__ import annotations

import pytest

pytest.importorskip("Bio")
pytest.importorskip("reportlab")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.visualization import biopython_p19_graphics_functions as m  # noqa: E402

TRACKS = [
    {"name": "genomeA", "length": 1000, "features": [
        {"start": 100, "end": 300, "strand": 1, "label": "geneA", "color": "steelblue"},
        {"start": 500, "end": 700, "strand": -1, "label": "geneB", "color": "#cc3333"}]},
    {"name": "genomeB", "length": 800, "features": [{"start": 50, "end": 250, "strand": -1, "label": "geneA2"}]},
    {"name": "genomeC", "length": 900},
]
LINKS = [
    {"track_a": "genomeA", "start_a": 100, "end_a": 300, "track_b": "genomeB", "start_b": 50, "end_b": 250, "flip": True, "color": "red"},
    {"track_a": 0, "start_a": 500, "end_a": 700, "track_b": 2, "start_b": 100, "end_b": 300},
]
CHROMOSOMES = [
    {"name": "chr1", "length": 1000, "features": [
        {"start": 100, "end": 200, "strand": 1, "label": "a", "color": "red"},
        {"start": 400, "end": 420, "strand": -1, "label": "b"}]},
    {"name": "chr2", "length": 700},
]


def _has_cairo() -> bool:
    try:
        import rlPyCairo  # noqa: F401
    except ImportError:
        return False
    return True


@pytest.mark.parametrize("fmt,marker", [("pdf", b"%PDF"), ("svg", b"<?xml"), ("eps", b"%!PS")])
def test_genome_diagram_formats(tmp_path, fmt, marker) -> None:
    target = tmp_path / f"g.{fmt}"
    result = m.genome_diagram_crosslinks(TRACKS, LINKS, target, image_format=fmt, verbose=False)
    assert target.read_bytes().lstrip()[:5].startswith(marker[:4])
    assert result["bytes"] == target.stat().st_size > 500
    assert [t["name"] for t in result["tracks"]] == ["genomeA", "genomeB", "genomeC"]
    assert result["links"][0]["flip"] is True and result["links"][1]["track_b"] == "genomeC"


def test_genome_diagram_svg_contains_labels_and_links(tmp_path) -> None:
    target = tmp_path / "g.svg"
    m.genome_diagram_crosslinks(TRACKS, LINKS, target, image_format="svg", title="demo", verbose=False)
    text = target.read_text()
    assert "geneA" in text and "genomeB" in text and "<path" in text


def test_genome_diagram_errors(tmp_path) -> None:
    out = tmp_path / "e.pdf"
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS[:1], LINKS, out, verbose=False)
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, [], out, verbose=False)
    same = [{"track_a": 0, "start_a": 1, "end_a": 5, "track_b": 0, "start_b": 1, "end_b": 5}]
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, same, out, verbose=False)
    beyond = [{"track_a": 0, "start_a": 1, "end_a": 5000, "track_b": 1, "start_b": 1, "end_b": 5}]
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, beyond, out, verbose=False)
    unknown_track = [{"track_a": "zzz", "start_a": 1, "end_a": 5, "track_b": 1, "start_b": 1, "end_b": 5}]
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, unknown_track, out, verbose=False)
    bad_color = [{**TRACKS[1], "features": [{"start": 1, "end": 5, "color": "not-a-colour"}]}]
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks([TRACKS[0], *bad_color], LINKS[:1], out, verbose=False)
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, LINKS, out, image_format="gif", verbose=False)
    out.write_bytes(b"x")
    with pytest.raises(ConfigurationError):
        m.genome_diagram_crosslinks(TRACKS, LINKS, out, verbose=False)  # exists


@pytest.mark.parametrize("fmt", ["pdf", "svg"])
def test_chromosome_ideogram(tmp_path, fmt) -> None:
    result = m.chromosome_ideogram(CHROMOSOMES, tmp_path / f"c.{fmt}", title="Demo", image_format=fmt, verbose=False)
    assert result["bytes"] > 1000 and result["telomere_length"] == 20
    assert [c["name"] for c in result["chromosomes"]] == ["chr1", "chr2"]
    if fmt == "svg":
        text = (tmp_path / "c.svg").read_text()
        assert "chr1" in text and "Demo" in text


def test_chromosome_ideogram_errors(tmp_path) -> None:
    with pytest.raises(ConfigurationError):
        m.chromosome_ideogram([], tmp_path / "a.pdf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.chromosome_ideogram([{"name": "c", "length": 10, "features": [{"start": 5, "end": 50}]}], tmp_path / "b.pdf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.chromosome_ideogram([{"name": "c", "length": 10, "extra": 1}], tmp_path / "c.pdf", verbose=False)


@pytest.mark.skipif(_has_cairo(), reason="rlPyCairo is installed")
def test_png_without_cairo_gives_install_hint(tmp_path) -> None:
    with pytest.raises(BackendUnavailableError):
        m.genome_diagram_crosslinks(TRACKS, LINKS, tmp_path / "g.png", image_format="png", verbose=False)
    assert not (tmp_path / "g.png").exists()
    with pytest.raises(BackendUnavailableError):
        m.chromosome_ideogram(CHROMOSOMES, tmp_path / "c.png", image_format="png", verbose=False)
