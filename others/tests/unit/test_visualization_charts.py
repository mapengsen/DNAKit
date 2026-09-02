"""Tests for similarity-matrix SVG rendering."""

from __future__ import annotations

import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.similarity import similarity_matrix
from dnakit.similarity.results import SimilarityMatrixResult
from dnakit.visualization import (
    HeatmapConfig,
    LimitPolicy,
    plot_similarity_matrix,
)

_SVG = "{http://www.w3.org/2000/svg}"


def _nodes(artifact_svg: str, tag: str, class_name: str) -> list[ET.Element]:
    root = ET.fromstring(artifact_svg)
    return [
        node
        for node in root.iter(f"{_SVG}{tag}")
        if class_name in node.attrib.get("class", "").split()
    ]


def test_similarity_heatmap_preserves_directional_matrix_values_and_order() -> None:
    result = similarity_matrix(
        [DNASequence("A"), DNASequence("AC")],
        method="kmer_containment",
        k=1,
    )

    artifact = plot_similarity_matrix(result, config=HeatmapConfig(show_values=True))

    cells = _nodes(artifact.svg, "rect", "heatmap-cell")
    triples = [
        (node.attrib["data-row"], node.attrib["data-column"], node.attrib["data-value"])
        for node in cells
    ]
    assert triples == [("0", "0", "1"), ("0", "1", "1"), ("1", "0", "0.5"), ("1", "1", "1")]
    assert artifact.metadata["symmetric"] is False
    assert artifact.width == artifact.height
    assert not _nodes(artifact.svg, "text", "plot-title")


def test_similarity_heatmap_escapes_labels_and_keeps_full_label_in_data() -> None:
    result = similarity_matrix(
        [
            DNARecord(DNASequence("A"), "left<&"),
            DNARecord(DNASequence("C"), "right"),
        ],
        method="exact",
    )

    artifact = plot_similarity_matrix(result)

    labels = _nodes(artifact.svg, "text", "row-label")
    assert labels[0].text == "left<&"
    assert labels[0].attrib["data-label"] == "left<&"
    assert "left&lt;&amp;" in artifact.svg


def test_similarity_heatmap_handles_empty_matrix() -> None:
    result = similarity_matrix([], method="exact")

    artifact = plot_similarity_matrix(result)

    assert _nodes(artifact.svg, "text", "empty-state")[0].text == "Empty matrix"
    assert any(issue.code == "VIZ_EMPTY_MATRIX" for issue in artifact.issues)


def test_similarity_heatmap_limit_errors_or_uses_same_stride_on_both_axes() -> None:
    records = [DNARecord(DNASequence("A"), f"r{index}") for index in range(5)]
    result = similarity_matrix(records, method="exact")

    with pytest.raises(ConfigurationError) as error:
        plot_similarity_matrix(result, config=HeatmapConfig(max_rows=2, max_columns=2))
    assert error.value.code == "VISUALIZATION_SIZE_LIMIT"

    artifact = plot_similarity_matrix(
        result,
        config=HeatmapConfig(
            max_rows=2,
            max_columns=2,
            limit_policy=LimitPolicy.STRIDE,
        ),
    )
    assert artifact.metadata["display_indices"] == (0, 3)
    cells = _nodes(artifact.svg, "rect", "heatmap-cell")
    assert {(node.attrib["data-row"], node.attrib["data-column"]) for node in cells} == {
        ("0", "0"),
        ("0", "3"),
        ("3", "0"),
        ("3", "3"),
    }


def test_similarity_heatmap_rejects_xml_controls_and_oversized_labels() -> None:
    control = SimilarityMatrixResult(
        name="matrix",
        method="manual",
        value_kind="similarity",
        labels=("bad\x00",),
        values=((1.0,),),
        symmetric=True,
        max_items=1,
        parameters=FrozenDict(),
    )
    oversized = SimilarityMatrixResult(
        name="matrix",
        method="manual",
        value_kind="similarity",
        labels=("x" * 1_025,),
        values=((1.0,),),
        symmetric=True,
        max_items=1,
        parameters=FrozenDict(),
    )

    with pytest.raises(ConfigurationError) as control_error:
        plot_similarity_matrix(control)
    assert control_error.value.code == "INVALID_VISUALIZATION_LABEL"
    with pytest.raises(ConfigurationError) as size_error:
        plot_similarity_matrix(oversized)
    assert size_error.value.code == "VISUALIZATION_SIZE_LIMIT"


def test_visualization_import_does_not_require_matplotlib() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; sys.modules['matplotlib'] = None; import dnakit.visualization",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
