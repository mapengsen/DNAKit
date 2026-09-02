"""Integration checks for the checksummed public DNA structure fixtures in temp/."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from dnakit.structure3d import (
    analyze_ensemble_flexibility,
    analyze_structure,
    load_pdb,
    load_pdb_ensemble,
    read_dssr_json,
)

_DATA = Path(__file__).resolve().parents[3] / "temp" / "dna_structures"


def _sample(name: str) -> Path:
    path = _DATA / name
    if not path.is_file():
        pytest.skip(f"downloaded validation fixture is absent: {path}")
    return path


def test_downloaded_b_dna_dodecamer_geometry() -> None:
    structure = load_pdb(_sample("1BNA.pdb"))
    result = analyze_structure(
        structure,
        sasa_points_per_atom=24,
        volume_grid_spacing_angstrom=1.5,
    )

    assert structure.source_sha256 == (
        "df42f1506792f191b957227b061360652adcf6f813eb69d9ec553067ea584670"
    )
    assert structure.sequence_by_chain == ("CGCGAATTCGCG", "CGCGAATTCGCG")
    assert result.residue_count == 24
    assert result.chain_count == 2
    assert result.helix is not None
    assert result.helix.mean_rise_angstrom == pytest.approx(3.3514, abs=0.01)
    assert result.helix.mean_twist_degree == pytest.approx(35.797, abs=0.1)
    assert result.shape.radius_of_gyration_angstrom == pytest.approx(13.228, abs=0.02)


@pytest.mark.parametrize(
    ("filename", "model_count", "residue_count", "expected_sequence"),
    (
        ("1AC7.pdb", 10, 16, "ATCCTAGTTATAGGAT"),
        ("139D.pdb", 4, 28, "TTGGGGT"),
    ),
)
def test_downloaded_nmr_ensembles_have_nonzero_flexibility(
    filename: str,
    model_count: int,
    residue_count: int,
    expected_sequence: str,
) -> None:
    models = load_pdb_ensemble(_sample(filename))
    flexibility = analyze_ensemble_flexibility(models)

    assert len(models) == model_count
    assert len(models[0].residues) == residue_count
    assert models[0].sequence_by_chain[0] == expected_sequence
    assert flexibility.mean_atomic_rmsf_angstrom > 0.0
    assert flexibility.max_atomic_rmsf_angstrom > flexibility.mean_atomic_rmsf_angstrom


def test_downloaded_dssr_schema_example_is_parsed() -> None:
    result = read_dssr_json(_sample("1ehz-dssr-example.json"))

    assert result.nucleotide_count == 76
    assert result.pair_count == 34
    assert result.helix_count == 2
    assert result.stem_count == 4
    assert result.hairpin_count == 3
    assert result.hydrogen_bond_count == 116
    assert result.backbone_torsion_record_count == 76


def test_downloaded_manifest_hashes_and_secondary_annotations() -> None:
    manifest_path = _sample("manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checked = 0
    for item in manifest["files"]:
        expected = item.get("sha256")
        if expected is None:
            continue
        path = _sample(item["filename"])
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
        checked += 1

    dbn = _sample("1AC7_secondary_structure.dbn").read_text(encoding="utf-8").splitlines()
    assert checked == 4
    assert dbn[1] == "ATCCTAGTTATAGGAT"
    assert dbn[2] == "((((((....))))))"
    assert "annotation_not_MFE" in dbn[0]
