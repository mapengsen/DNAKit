"""Real Pyrodigal coordinate and translation checks, plus input contracts."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNASequence, Topology
from dnakit.exceptions import ConfigurationError
from dnakit.genomics.genes import predict_prokaryotic_genes

# A synthetic, non-organism-derived CDS with explicit stop codons and an RBS.
# The coding region's translation is known by construction, independently of
# the adapter. No external genome or model download is needed for these tests.
_CDS = "ATG" + "GCTGACGAACTGAAAGGT" * 20 + "TAA"
_SEQUENCE = "TAA" * 10 + "AGGAGG" + "A" * 7 + _CDS + "TAA" * 10


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@pytest.fixture
def real_pyrodigal() -> None:
    pytest.importorskip("pyrodigal", minversion="3.4")


@pytest.mark.usefixtures("real_pyrodigal")
def test_known_cds_and_json_provenance() -> None:
    result = predict_prokaryotic_genes(
        _SEQUENCE, sequence_id="synthetic", closed=True, verbose=False
    )
    assert len(result.genes) == 1
    gene = result.genes[0]
    assert (gene.start, gene.end, gene.strand) == (43, 409, 1)
    assert gene.cds == _CDS
    assert gene.protein == "M" + "ADELKG" * 20 + "*"
    assert gene.translation_table == 11
    assert not gene.partial_begin and not gene.partial_end
    assert 0 <= gene.confidence_percent <= 100
    payload = json.loads(json.dumps(result.to_dict(), allow_nan=False))
    assert payload["genes"][0]["original_begin"] == gene.start + 1
    assert payload["genes"][0]["original_end"] == gene.end
    assert payload["backend_version"]
    assert payload["prodigal_version"]
    assert payload["coordinate_system"] == "0-based half-open"
    assert len(payload["citations"]) == 2


@pytest.mark.usefixtures("real_pyrodigal")
def test_reverse_strand_coordinates_and_oriented_cds() -> None:
    forward = predict_prokaryotic_genes(_SEQUENCE, closed=True, verbose=False)
    reverse_sequence = _reverse_complement(_SEQUENCE)
    reverse = predict_prokaryotic_genes(reverse_sequence, closed=True, verbose=False)
    assert len(forward.genes) == len(reverse.genes) == 1
    gene, reverse_gene = forward.genes[0], reverse.genes[0]
    assert reverse_gene.strand == -1
    assert reverse_gene.start == len(_SEQUENCE) - gene.end
    assert reverse_gene.end == len(_SEQUENCE) - gene.start
    assert reverse_gene.cds == _reverse_complement(
        reverse_sequence[reverse_gene.start : reverse_gene.end]
    )
    assert reverse_gene.cds == gene.cds
    assert reverse_gene.protein == gene.protein


@pytest.mark.usefixtures("real_pyrodigal")
def test_terminal_stop_is_optional_without_coordinate_changes() -> None:
    stopped = predict_prokaryotic_genes(_SEQUENCE, closed=True, verbose=False)
    unstopped = predict_prokaryotic_genes(_SEQUENCE, closed=True, include_stop=False, verbose=False)
    assert unstopped.genes[0].protein == stopped.genes[0].protein.removesuffix("*")
    assert unstopped.genes[0].cds == stopped.genes[0].cds
    assert unstopped.genes[0].end == stopped.genes[0].end


@pytest.mark.usefixtures("real_pyrodigal")
def test_ambiguous_input_preserves_positions_and_empty_predictions() -> None:
    result = predict_prokaryotic_genes("RYSWKMBDHVN" * 30, verbose=False)
    assert result.ambiguous_bases_to_n == 330
    assert result.sequence_length == 330
    assert result.genes == ()


@pytest.mark.usefixtures("real_pyrodigal")
def test_short_training_is_rejected_without_falling_back_to_meta() -> None:
    with pytest.raises(ConfigurationError, match="Pyrodigal rejected"):
        predict_prokaryotic_genes(_SEQUENCE, meta=False, training_sequence=_SEQUENCE, verbose=False)


@pytest.mark.usefixtures("real_pyrodigal")
def test_explicit_training_mode_uses_supplied_sequence_and_table() -> None:
    result = predict_prokaryotic_genes(
        _SEQUENCE,
        meta=False,
        training_sequence=_SEQUENCE * 250,
        translation_table=11,
        closed=True,
        verbose=False,
    )
    assert result.meta is False
    assert result.training_translation_table == 11
    assert result.training_sha256 is not None
    assert result.training_sha256 != result.input_sha256
    assert len(result.genes) == 1
    assert result.genes[0].cds == _CDS
    assert result.genes[0].protein == "M" + "ADELKG" * 20 + "*"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"meta": False},
        {"training_sequence": _SEQUENCE},
        {"translation_table": 11},
        {"meta": 1},
        {"closed": "yes"},
        {"mask": None},
        {"sequence_id": " "},
        {"meta": False, "training_sequence": _SEQUENCE, "translation_table": True},
        {"meta": False, "training_sequence": _SEQUENCE, "translation_table": 0},
    ],
)
def test_invalid_mode_contract(kwargs: dict[str, object]) -> None:
    with pytest.raises(ConfigurationError):
        predict_prokaryotic_genes(_SEQUENCE, verbose=False, **kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize("sequence", ["", "AC-GT", "AC GT", "ACUGT", ">seq\nACGT"])
def test_invalid_sequence_is_not_silently_accepted(sequence: str) -> None:
    with pytest.raises(ConfigurationError):
        predict_prokaryotic_genes(sequence, verbose=False)


def test_circular_object_requires_explicit_linearization() -> None:
    sequence = DNASequence(_SEQUENCE, topology=Topology.CIRCULAR)
    with pytest.raises(ConfigurationError, match="linearize"):
        predict_prokaryotic_genes(sequence, verbose=False)
