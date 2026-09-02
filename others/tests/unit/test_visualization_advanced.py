"""Tests for advanced dependency-free SVG renderers."""

import xml.etree.ElementTree as ET

import pytest

from dnakit.alignment import align_pairwise
from dnakit.core import DNAFeature, DNARecord, DNASequence, Interval, Topology
from dnakit.exceptions import ConfigurationError
from dnakit.visualization import (
    plot_alignment,
    plot_circular_map,
    plot_linear_map,
)


def test_linear_and_circular_maps_are_valid_deterministic_svg() -> None:
    linear = DNARecord(
        DNASequence("ACGTACGT"), "linear", features=[DNAFeature("gene", Interval(1, 6))]
    )
    circular = DNARecord(
        DNASequence("ACGTACGT", topology=Topology.CIRCULAR),
        "circular",
        features=[DNAFeature("gene", Interval(6, 8))],
    )

    linear_artifact = plot_linear_map(linear)
    circular_artifact = plot_circular_map(circular)

    ET.fromstring(linear_artifact.svg)
    ET.fromstring(circular_artifact.svg)
    assert linear_artifact.kind == "linear-map"
    assert circular_artifact.kind == "circular-map"
    assert linear_artifact.width == linear_artifact.height
    assert circular_artifact.width == circular_artifact.height
    assert plot_circular_map(circular).svg == circular_artifact.svg


def test_alignment_plot_uses_precomputed_result() -> None:
    result = align_pairwise(DNASequence("ACGT"), DNASequence("AGT"))
    artifact = plot_alignment(result, columns_per_line=3)

    root = ET.fromstring(artifact.svg)
    assert root.attrib["data-kind"] == "alignment"
    assert artifact.metadata["columns"] == result.alignment_length
    assert artifact.width == artifact.height


def test_full_circle_feature_is_rendered_and_plot_parameters_are_validated() -> None:
    circular = DNARecord(
        DNASequence("ACGT", topology=Topology.CIRCULAR),
        "circle",
        features=[DNAFeature("source", Interval(0, 4))],
    )
    root = ET.fromstring(plot_circular_map(circular).svg)
    feature_circles = [
        item
        for item in root
        if item.tag.endswith("circle") and item.attrib.get("class") == "feature"
    ]
    assert len(feature_circles) == 1

    with pytest.raises(ConfigurationError):
        plot_linear_map(DNASequence("AC"), width=True)
    with pytest.raises(ConfigurationError):
        plot_alignment(align_pairwise(DNASequence("A"), DNASequence("A")), max_columns=0)
