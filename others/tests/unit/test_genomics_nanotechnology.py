"""Real scadnano complement, geometry, and oxDNA topology checks.

No molecular dynamics simulation or experimental folding validation is claimed.
"""

from __future__ import annotations

import copy
import json
import math

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.genomics.nanotechnology import design_dna_duplex, export_oxdna_design


@pytest.fixture
def real_scadnano() -> None:
    pytest.importorskip("scadnano")


@pytest.fixture
def duplex(real_scadnano: None) -> dict:
    return design_dna_duplex("acgttgca", verbose=False)


def test_duplex_assigns_reverse_complement_and_antiparallel_domains(duplex: dict) -> None:
    forward, reverse = duplex["strands"]
    assert forward["sequence"] == "ACGTTGCA"
    assert reverse["sequence"] == "TGCAACGT"
    assert forward["domains"][0] == {"helix": 0, "forward": True, "start": 0, "end": 8}
    assert reverse["domains"][0] == {"helix": 0, "forward": False, "start": 0, "end": 8}
    json.dumps(duplex, allow_nan=False)


def test_export_has_valid_sequences_bonds_and_finite_geometry_without_writes(
    duplex: dict, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    before = copy.deepcopy(duplex)
    result = export_oxdna_design(duplex, verbose=False)
    assert duplex == before
    assert list(tmp_path.iterdir()) == []
    assert result["simulation_executed"] is False
    assert result["nucleotide_count"] == 16
    assert result["strand_count"] == 2
    assert result["backend_version"]
    topology = result["topology"].splitlines()
    assert topology[0] == "16 2"
    rows = [line.split() for line in topology[1:]]
    assert "".join(row[1] for row in rows[:8]) == "ACGTTGCA"
    assert "".join(row[1] for row in rows[8:]) == "TGCAACGT"
    for index, row in enumerate(rows):
        strand, _, n3, n5 = row
        assert int(strand) == index // 8 + 1
        for neighbor, reciprocal in [(int(n3), 3), (int(n5), 2)]:
            if neighbor != -1:
                assert 0 <= neighbor < 16
                assert rows[neighbor][0] == strand
                assert int(rows[neighbor][reciprocal]) == index
    assert rows[0][3] == rows[8][3] == "-1"
    assert rows[7][2] == rows[15][2] == "-1"
    positions = []
    configuration = result["configuration"].splitlines()
    assert len(configuration) == 19
    for line in configuration[3:]:
        values = [float(value) for value in line.split()]
        assert len(values) == 15 and all(math.isfinite(value) for value in values)
        positions.append(tuple(values[:3]))
        backbone, normal = values[3:6], values[6:9]
        assert sum(x * y for x, y in zip(backbone, normal, strict=True)) == pytest.approx(
            0, abs=1e-10
        )
        assert sum(x * x for x in backbone) == pytest.approx(1)
        assert sum(x * x for x in normal) == pytest.approx(1)
    assert len(set(positions)) == 16


@pytest.mark.parametrize("sequence", [None, 42, "", "ACNT", "ACUT", "AC GT", "A" * 100_001])
def test_invalid_duplex_sequence_fails_before_optional_import(
    sequence: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    def unexpected_import(*args: object) -> None:
        pytest.fail("Invalid sequence must not load optional scadnano")

    monkeypatch.setattr("dnakit.genomics.nanotechnology.require", unexpected_import)
    with pytest.raises(ConfigurationError):
        design_dna_duplex(sequence, verbose=False)  # type: ignore[arg-type]


@pytest.mark.parametrize("sequence", [None, "ACGT", "ACGNTGCA", "ACG?TGCA"])
def test_unassigned_incomplete_or_wildcard_strands_cannot_export(
    duplex: dict, sequence: object
) -> None:
    duplex["strands"][0]["sequence"] = sequence
    with pytest.raises(ConfigurationError):
        export_oxdna_design(duplex, verbose=False)


def test_circular_strand_cannot_silently_export_as_linear(duplex: dict) -> None:
    duplex["strands"][0]["circular"] = True
    with pytest.raises(ConfigurationError, match="Circular"):
        export_oxdna_design(duplex, verbose=False)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_nonfinite_input_geometry_is_rejected(duplex: dict, value: float) -> None:
    duplex["helices"][0]["roll"] = value
    with pytest.raises(ConfigurationError, match="finite"):
        export_oxdna_design(duplex, verbose=False)


def test_nonfinite_derived_geometry_is_rejected(duplex: dict) -> None:
    duplex["geometry"] = {"rise_per_base_pair": 1e308}
    with pytest.raises(ConfigurationError, match="finite"):
        export_oxdna_design(duplex, verbose=False)


@pytest.mark.parametrize("damage", ["zero_turn", "invalid_helix", "no_strands"])
def test_unexportable_geometry_has_an_actionable_domain_error(duplex: dict, damage: str) -> None:
    if damage == "zero_turn":
        duplex["geometry"] = {"bases_per_turn": 0}
    elif damage == "invalid_helix":
        duplex["strands"][0]["domains"][0]["helix"] = 99
    else:
        duplex["strands"] = []
    with pytest.raises(ConfigurationError):
        export_oxdna_design(duplex, verbose=False)


def test_nonzero_min_offset_never_leaks_upstream_index_error(duplex: dict) -> None:
    # scadnano 0.21.1 has an offset-indexing bug for nonzero min_offset. Future
    # corrected backends may export this valid design; do not require that bug.
    duplex["helices"][0].update(min_offset=8, max_offset=16)
    for strand in duplex["strands"]:
        strand["domains"][0].update(start=8, end=16)
    before = copy.deepcopy(duplex)
    try:
        result = export_oxdna_design(duplex, verbose=False)
    except ConfigurationError as exc:
        assert "scadnano cannot export" in str(exc)
    else:
        assert result["nucleotide_count"] == 16
        assert result["strand_count"] == 2
    assert duplex == before
