"""Tests for the dependency-free sequence SVG renderer."""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap, Topology
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import Highlight, LimitPolicy, SequencePlotConfig, plot_sequence

_SVG = "{http://www.w3.org/2000/svg}"


def _nodes(artifact_svg: str, tag: str, class_name: str) -> list[ET.Element]:
    root = ET.fromstring(artifact_svg)
    return [
        node
        for node in root.iter(f"{_SVG}{tag}")
        if class_name in node.attrib.get("class", "").split()
    ]


def test_sequence_plot_is_valid_deterministic_svg_for_iupac() -> None:
    sequence = DNASequence("ARYN", alphabet=DNAAlphabet.IUPAC)

    first = plot_sequence(sequence)
    second = plot_sequence(sequence)

    assert first.kind == "sequence"
    assert first.svg == second.svg
    assert first.sha256 == second.sha256
    assert first.width == first.height
    assert "matplotlib" not in first.svg
    assert not _nodes(first.svg, "text", "plot-title")
    bases = _nodes(first.svg, "text", "base")
    assert [node.text for node in bases] == list("ARYN")


def test_sequence_plot_escapes_record_id_title_and_highlight_label() -> None:
    record = DNARecord(DNASequence("ACGT"), "sample<&")
    artifact = plot_sequence(
        record,
        highlights=[Highlight(0, 2, label="tag<&", color="#abc")],
    )

    root = ET.fromstring(artifact.svg)
    title = root.find(f"{_SVG}title")
    assert title is not None
    assert title.text == "DNA sequence: sample<&"
    assert _nodes(artifact.svg, "text", "highlight-label")[0].text == "tag<&"
    assert "sample&lt;&amp;" in artifact.svg


def test_sequence_plot_shows_known_and_unknown_gaps_without_expanding_them() -> None:
    sequence = DNASequence(
        ["AC", Gap(500), "T", Gap(None), "G"],
        alphabet=DNAAlphabet.STRICT,
    )

    artifact = plot_sequence(sequence)

    assert [node.text for node in _nodes(artifact.svg, "text", "gap-label")] == [
        "[500 bp]",
        "[… bp]",
    ]
    assert len(_nodes(artifact.svg, "text", "base")) == 4
    final_base = _nodes(artifact.svg, "text", "base")[-1]
    assert final_base.text == "G"
    assert final_base.attrib["data-coordinate"] == "?"
    assert artifact.metadata["unknown_gap_count"] == 1


def test_sequence_plot_can_show_iupac_complement() -> None:
    sequence = DNASequence("ARYK", alphabet=DNAAlphabet.IUPAC)

    artifact = plot_sequence(sequence, config=SequencePlotConfig(show_complement=True))

    complements = _nodes(artifact.svg, "text", "complement")
    assert [node.text for node in complements] == list("TYRM")


def test_sequence_plot_can_render_custom_symbols_without_changing_semantics() -> None:
    source_map = {"A": "*", "T": "-", "C": "+", "G": "]"}
    config = SequencePlotConfig(
        bases_per_line=4,
        show_complement=True,
        symbol_map=source_map,
    )
    source_map["A"] = "!"

    artifact = plot_sequence(DNASequence("ATCG"), config=config)
    bases = _nodes(artifact.svg, "text", "base")
    complements = _nodes(artifact.svg, "text", "complement")

    assert [node.text for node in bases] == ["*", "-", "+", "]"]
    assert [node.attrib["data-source-symbol"] for node in bases] == list("ATCG")
    assert [node.attrib["fill"] for node in bases] == [
        "#16a34a",
        "#dc2626",
        "#2563eb",
        "#d97706",
    ]
    assert [node.attrib["data-coordinate"] for node in bases] == ["0", "1", "2", "3"]
    assert [node.text for node in complements] == ["-", "*", "]", "+"]
    assert [node.attrib["data-source-symbol"] for node in complements] == list("TAGC")
    assert artifact.metadata["symbol_map"] == {"A": "*", "C": "+", "G": "]", "T": "-"}
    assert config.symbol_map["A"] == "*"
    assert hash(config)


def test_sequence_plot_custom_symbols_fall_back_and_escape_xml() -> None:
    artifact = plot_sequence(
        DNASequence("AN", alphabet=DNAAlphabet.IUPAC),
        config=SequencePlotConfig(symbol_map={"A": "<"}),
    )

    assert [node.text for node in _nodes(artifact.svg, "text", "base")] == ["<", "N"]
    assert "&lt;" in artifact.svg


@pytest.mark.parametrize(
    "symbol_map",
    [
        [("A", "*")],
        {"a": "*"},
        {"U": "*"},
        {"AA": "*"},
        {"A": ""},
        {"A": "++"},
        {"A": " "},
        {"A": "\x00"},
        {"A": 1},
    ],
)
def test_sequence_plot_rejects_invalid_symbol_map(symbol_map: object) -> None:
    with pytest.raises(ConfigurationError) as error:
        SequencePlotConfig(symbol_map=symbol_map)  # type: ignore[arg-type]
    assert error.value.code == "INVALID_VISUALIZATION_SYMBOL_MAP"


def test_sequence_plot_controls_line_and_column_spacing() -> None:
    sequence = DNASequence("ACGTACGT")
    default = plot_sequence(sequence, config=SequencePlotConfig(bases_per_line=4))
    spaced = plot_sequence(
        sequence,
        config=SequencePlotConfig(
            bases_per_line=4,
            column_spacing=7,
            line_spacing=19,
        ),
    )

    bases = _nodes(spaced.svg, "text", "base")
    assert float(bases[1].attrib["x"]) - float(bases[0].attrib["x"]) == 19
    assert float(bases[4].attrib["y"]) - float(bases[0].attrib["y"]) == 81
    assert spaced.width == spaced.height
    assert default.width == default.height
    assert spaced.metadata["column_spacing"] == 7
    assert spaced.metadata["line_spacing"] == 19

    wide = plot_sequence(
        DNASequence("A" * 40),
        config=SequencePlotConfig(
            bases_per_line=40,
            column_spacing=7,
            show_coordinates=False,
        ),
    )
    assert wide.width == 40 * 12 + 39 * 7 + 2 * 24
    assert wide.height == wide.width


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("column_spacing", -1),
        ("column_spacing", True),
        ("line_spacing", -1),
        ("line_spacing", 1.5),
    ],
)
def test_sequence_plot_rejects_invalid_spacing(name: str, value: object) -> None:
    with pytest.raises(ConfigurationError) as error:
        SequencePlotConfig(**{name: value})  # type: ignore[arg-type]
    assert error.value.code == "INVALID_VISUALIZATION_CONFIG"


def test_sequence_highlight_priority_then_input_order_is_deterministic() -> None:
    artifact = plot_sequence(
        DNASequence("AAAA"),
        highlights=(
            Highlight(0, 4, color="#111111", priority=1),
            Highlight(1, 3, color="#222222", priority=1),
            Highlight(2, 4, color="#333333", priority=2),
        ),
    )

    rectangles = _nodes(artifact.svg, "rect", "highlight")
    assert [node.attrib["data-highlight-index"] for node in rectangles] == ["0", "1", "2", "2"]


def test_sequence_plot_handles_empty_and_circular_inputs() -> None:
    empty = plot_sequence(DNASequence(""))
    circular = plot_sequence(DNASequence("AC", topology=Topology.CIRCULAR))

    assert _nodes(empty.svg, "text", "empty-state")[0].text == "Empty sequence"
    assert any(issue.code == "VIZ_CIRCULAR_LINEARIZED" for issue in circular.issues)


def test_sequence_size_limit_errors_unless_truncation_is_explicit() -> None:
    sequence = DNASequence("ACGTAC")

    with pytest.raises(ConfigurationError) as error:
        plot_sequence(sequence, config=SequencePlotConfig(max_symbols=4))
    assert error.value.code == "VISUALIZATION_SIZE_LIMIT"

    artifact = plot_sequence(
        sequence,
        highlights=[Highlight(3, 6)],
        config=SequencePlotConfig(max_symbols=4, limit_policy=LimitPolicy.TRUNCATE),
    )
    assert len(_nodes(artifact.svg, "text", "base")) == 4
    assert {issue.code for issue in artifact.issues} == {
        "VIZ_HIGHLIGHT_CLIPPED",
        "VIZ_SEQUENCE_TRUNCATED",
    }


def test_sequence_plot_rejects_invalid_highlights_and_too_many_gaps() -> None:
    with pytest.raises(ConfigurationError) as range_error:
        plot_sequence(DNASequence("AC"), highlights=[Highlight(1, 3)])
    assert range_error.value.code == "HIGHLIGHT_OUT_OF_RANGE"

    gapped = DNASequence(["A", Gap(1), "C", Gap(1), "G"])
    with pytest.raises(ConfigurationError) as gap_error:
        plot_sequence(gapped, config=SequencePlotConfig(max_gaps=1))
    assert gap_error.value.code == "VISUALIZATION_SIZE_LIMIT"


def test_sequence_plot_rejects_xml_control_characters_with_stable_error() -> None:
    record = DNARecord(DNASequence("AC"), "bad\x00id")

    with pytest.raises(ConfigurationError) as error:
        plot_sequence(record)

    assert error.value.code == "INVALID_VISUALIZATION_TEXT"


def test_sequence_plot_does_not_mutate_input() -> None:
    sequence = DNASequence(["AC", Gap(None), "GT"])
    parts_before = sequence.parts

    plot_sequence(sequence)

    assert sequence.parts == parts_before
