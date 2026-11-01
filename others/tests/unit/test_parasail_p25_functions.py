"""parasail traceback rendering (text view and CIGAR)."""

from __future__ import annotations

import pytest

pytest.importorskip("parasail")
import parasail

from dnakit.alignment.parasail_family_functions import _ALGORITHMS
from dnakit.alignment.parasail_p25_functions import parasail_render_alignment
from dnakit.exceptions import ConfigurationError

QUERY = "ACGTACGGT"
TARGET = "TTACGTCCGTAA"
MATRIX = parasail.matrix_create("ACGTN", 5, -4)


@pytest.mark.parametrize("algorithm", _ALGORITHMS)
def test_matches_direct_parasail_call(algorithm: str) -> None:
    result = parasail_render_alignment(QUERY, TARGET, algorithm=algorithm, verbose=False)
    direct = getattr(parasail, f"{algorithm}_trace_striped_sat")(QUERY, TARGET, 10, 1, MATRIX)
    trace = direct.get_traceback()
    assert result["score"] == direct.score
    assert result["query_line"] == trace.query
    assert result["comparison_line"] == trace.comp
    assert result["target_line"] == trace.ref
    assert result["cigar"] == direct.get_cigar().decode.decode()
    assert result["text"].split("\n") == [trace.query, trace.comp, trace.ref]


def test_views() -> None:
    text = parasail_render_alignment(QUERY, TARGET, view="text", verbose=False)
    cigar = parasail_render_alignment(QUERY, TARGET, view="cigar", verbose=False)
    assert "cigar" not in text and "text" in text
    assert "text" not in cigar and "cigar" in cigar


def test_identical_sequences_edge_case() -> None:
    result = parasail_render_alignment("ACGT", "ACGT", algorithm="nw", verbose=False)
    assert result["cigar"] == "4=" and result["comparison_line"] == "||||"


def test_invalid_arguments() -> None:
    with pytest.raises(ConfigurationError):
        parasail_render_alignment(QUERY, TARGET, algorithm="bogus", verbose=False)
    with pytest.raises(ConfigurationError):
        parasail_render_alignment(QUERY, TARGET, view="x", verbose=False)
    with pytest.raises(ConfigurationError):
        parasail_render_alignment(QUERY, "ACGU!", verbose=False)
