"""Tests for bounded PDB geometry and standard 3DNA/DSSR parsers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dnakit.exceptions import InputFormatError
from dnakit.structure3d import (
    analyze_ensemble_flexibility,
    analyze_structure,
    load_pdb,
    load_pdb_ensemble,
    read_3dna_bp_step,
    read_dssr_json,
)


def _pdb_atom(
    serial: int,
    name: str,
    residue: str,
    chain: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain:1s}{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}  "
    )


def _minimal_duplex(model: int, shift: float = 0.0) -> list[str]:
    lines = [f"MODEL     {model:4d}"]
    residues = (
        ("DA", "A", 1, 0.0, 0.0, 0.0),
        ("DT", "A", 2, 0.0, 0.0, 3.4),
        ("DA", "B", 1, 10.0, 0.0, 3.4),
        ("DT", "B", 2, 10.0, 0.0, 0.0),
    )
    serial = 1
    for residue, chain, number, x, y, z in residues:
        for name, dx, dy, element in (
            ("P", -1.0, 0.0, "P"),
            ("C4'", 0.0, 0.0, "C"),
            ("C1'", 1.0, 0.0, "C"),
            ("N1", 2.0, 0.0, "N"),
        ):
            lines.append(
                _pdb_atom(
                    serial,
                    name,
                    residue,
                    chain,
                    number,
                    x + dx + shift,
                    y + dy,
                    z,
                    element,
                )
            )
            serial += 1
    lines.append("ENDMDL")
    return lines


def test_pdb_parser_geometry_progress_and_ensemble_flexibility(tmp_path: Path) -> None:
    pdb = tmp_path / "fixture.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HEADER    DNA                                     01-JAN-00   TST1",
                "TITLE     SYNTHETIC DNA TEST FIXTURE",
                *_minimal_duplex(1),
                *_minimal_duplex(2, shift=0.2),
                "END",
            ]
        )
        + "\n",
        encoding="ascii",
    )
    models = load_pdb_ensemble(pdb)
    progress: list[tuple[str, int, int]] = []
    result = analyze_structure(
        models[0],
        sasa_points_per_atom=24,
        volume_grid_spacing_angstrom=2.0,
        progress=lambda stage, completed, total: progress.append((stage, completed, total)),
    )
    flexibility = analyze_ensemble_flexibility(models)

    assert len(models) == 2
    assert models[0].pdb_id == "TST1"
    assert models[0].sequence_by_chain == ("AT", "AT")
    assert result.atom_count == 16
    assert result.residue_count == 4
    assert result.chain_count == 2
    assert result.shape.radius_of_gyration_angstrom > 0.0
    assert result.solvent_accessible_surface_area_angstrom2 > 0.0
    assert result.molecular_volume_angstrom3 > 0.0
    assert result.helix is not None
    assert {stage for stage, _, _ in progress} == {"sasa", "volume"}
    assert flexibility.model_count == 2
    assert flexibility.mean_atomic_rmsf_angstrom == pytest.approx(0.0, abs=1e-12)


def test_pdb_parser_rejects_non_dna_coordinate_files(tmp_path: Path) -> None:
    pdb = tmp_path / "protein.pdb"
    pdb.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C") + "\n",
        encoding="ascii",
    )

    with pytest.raises(InputFormatError) as error:
        load_pdb(pdb)
    assert error.value.code == "PDB_DNA_ATOMS_NOT_FOUND"


def test_3dna_parameter_parser_preserves_all_twelve_standard_values(tmp_path: Path) -> None:
    parameters = tmp_path / "bp_step.par"
    parameters.write_text(
        "2 # base-pairs\n"
        "# Shear Stretch Stagger Buckle Prop-T Opening Shift Slide Rise Tilt Roll Twist\n"
        "A-T 0.1 0.2 0.3 1 2 3 0.4 0.5 3.4 4 5 36\n"
        "T-A 0.2 0.3 0.4 2 3 4 0.5 0.6 3.5 5 6 35\n",
        encoding="utf-8",
    )
    result = read_3dna_bp_step(parameters)

    assert result.base_pair_count == 2
    assert result.parameters[0].shear_angstrom == 0.1
    assert result.parameters[1].twist_degree == 35.0
    assert result.mean_rise_angstrom == 3.5
    assert result.mean_twist_degree == 35.0
    assert result.base_pairs_per_turn == pytest.approx(360.0 / 35.0)


def test_dssr_json_summary_supports_current_and_legacy_torsion_keys(tmp_path: Path) -> None:
    payload = {
        "dssr_ver": "DSSR test",
        "num_nts": 2,
        "nts": [{}, {}],
        "num_pairs": 1,
        "pairs": [{}],
        "num_helices": 1,
        "helices": [{}],
        "num_stems": 1,
        "stems": [{}],
        "num_hairpins": 0,
        "hairpins": [],
        "num_hbonds": 2,
        "hbonds": [{}, {}],
        "ntPars": [{}, {}],
    }
    source = tmp_path / "dssr.json"
    source.write_text(json.dumps(payload), encoding="utf-8")
    result = read_dssr_json(source)

    assert result.nucleotide_count == 2
    assert result.pair_count == 1
    assert result.hydrogen_bond_count == 2
    assert result.backbone_torsion_record_count == 2
