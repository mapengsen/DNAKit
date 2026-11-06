"""Tests for dnakit.variants.hgvs_gap_functions (sequence-context text, full house, alignment context).

Text helpers are compared with the hgvs output on parsed HGVS strings. The mapper-based functions use a stub
mapper because UTA and a sequence store are not available offline; the alignment context is only checked for
its dependency error here.
"""

from __future__ import annotations

import inspect
import sys

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pytest.importorskip("hgvs")
context = pytest.importorskip("hgvs.utils.context")

from dnakit.variants import hgvs_gap_functions as facades

INSERTION = "NM_000399.3:c.2689_2690insA"
DUPLICATION = "NC_000010.10:g.64572045dupT"


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_sequence_line_matches_hgvs_layout():
    from hgvs.parser import Parser

    line = facades.hgvs_sequence_line(
        INSERTION, (2680, 2700), "ACGTACGTAC", direction=">", verbose=False
    )
    expected = context.seq_line_fmt(
        var=Parser().parse_hgvs_variant(INSERTION), span=(2680, 2700), content="ACGTACGTAC", dir=">"
    )
    assert line == expected
    assert "NM_000399.3" in line and "2680 >" in line and line.endswith(INSERTION)


def test_pointer_line_marks_the_edit():
    line = facades.hgvs_pointer_line(DUPLICATION, (64572025, 64572065), verbose=False)
    assert "*" in line
    assert line.endswith("NC_000010.10:g.64572045dup")


def test_pointer_line_marks_insertions_with_arrow():
    line = facades.hgvs_pointer_line(INSERTION, (2680, 2700), verbose=False)
    assert "><" in line


def test_text_helpers_reject_bad_input():
    with pytest.raises(ConfigurationError):
        facades.hgvs_sequence_line("not an hgvs variant", (1, 2), "AC", verbose=False)
    with pytest.raises(ConfigurationError):
        facades.hgvs_sequence_line(INSERTION, (10, 5), "AC", verbose=False)
    with pytest.raises(ConfigurationError):
        facades.hgvs_sequence_line(INSERTION, (1, 2), "AC", direction="x", verbose=False)
    with pytest.raises(ConfigurationError):
        facades.hgvs_pointer_line("NP_000390.2:p.Ser10Phe", (1, 20), verbose=False)


def test_format_sequence_blocks_matches_hgvs():
    sequence = "ACGT" * 30
    text = facades.format_sequence_blocks(sequence, start=0, end=120, verbose=False)
    assert text == "".join(context.format_sequence(sequence, 0, 120))
    first_row = text.splitlines()[1]
    assert first_row.strip().startswith("1 :")
    assert "ACG TAC GTA" in first_row


def test_format_sequence_blocks_validates_arguments():
    with pytest.raises(ConfigurationError):
        facades.format_sequence_blocks("", verbose=False)
    with pytest.raises(ConfigurationError):
        facades.format_sequence_blocks("ACGT", start=2, end=2, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.format_sequence_blocks("ACGT", group_size=0, verbose=False)
    with pytest.raises(ConfigurationError):
        facades.format_sequence_blocks("ACGT", end=5, verbose=False)


class StubMapper:
    """Minimal stand-in for an AlignmentMapper: returns fixed strings for the conversions."""

    def __init__(self, transcripts=("NM_000399.3",)):
        self.transcripts = list(transcripts)

    def relevant_transcripts(self, var_g):
        return list(self.transcripts)

    def c_to_g(self, var_c):
        return "NC_000010.10:g.64572045dupT"

    def c_to_n(self, var_c):
        return "NM_000399.3:n.2696dupA"

    def c_to_p(self, var_c):
        return "NP_000390.2:p.Ser10Phe"

    def g_to_n(self, var_g, tx_ac):
        return "NM_000399.3:n.2696dupA"

    def n_to_c(self, var_n):
        return "NM_000399.3:c.2689_2690insA"

    def n_to_g(self, var_n):
        return "NC_000010.10:g.64572045dupT"


def test_full_house_returns_four_forms_from_a_coding_variant():
    result = facades.hgvs_full_house(StubMapper(), INSERTION, verbose=False)
    assert set(result) == {"g", "c", "n", "p"}
    assert result["c"] == INSERTION
    assert result["g"] == "NC_000010.10:g.64572045dupT"


def test_full_house_needs_a_transcript_when_several_apply():
    mapper = StubMapper(transcripts=("NM_1", "NM_2"))
    with pytest.raises(ConfigurationError):
        facades.hgvs_full_house(mapper, DUPLICATION, verbose=False)
    picked = facades.hgvs_full_house(mapper, DUPLICATION, transcript="NM_1", verbose=False)
    assert picked["n"] == "NM_000399.3:n.2696dupA"


def test_full_house_requires_a_mapper():
    with pytest.raises(ConfigurationError):
        facades.hgvs_full_house(None, INSERTION, verbose=False)


def test_alignment_context_reports_missing_uta_align(monkeypatch):
    monkeypatch.setitem(sys.modules, "uta_align", None)
    monkeypatch.setitem(sys.modules, "uta_align.align", None)
    monkeypatch.setitem(sys.modules, "uta_align.align.algorithms", None)
    with pytest.raises(BackendUnavailableError) as info:
        facades.hgvs_variant_context_with_alignment(StubMapper(), INSERTION, verbose=False)
    assert "uta-align" in (info.value.hint or "")
