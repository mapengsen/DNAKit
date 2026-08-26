"""Immutable coordinate, geometry, flexibility, and 3DNA/DSSR result objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import Provenance
from dnakit.core._json import to_json_compatible


class _SerializableResult:
    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class Atom3D(_SerializableResult):
    serial: int
    name: str
    element: str
    residue_name: str
    chain_id: str
    residue_number: int
    insertion_code: str
    x_angstrom: float
    y_angstrom: float
    z_angstrom: float
    occupancy: float
    model_index: int

    @property
    def coordinates(self) -> tuple[float, float, float]:
        return (self.x_angstrom, self.y_angstrom, self.z_angstrom)


@dataclass(frozen=True)
class Residue3D(_SerializableResult):
    name: str
    base: str
    chain_id: str
    number: int
    insertion_code: str
    atoms: tuple[Atom3D, ...]


@dataclass(frozen=True)
class DNA3DStructure(_SerializableResult):
    pdb_id: str | None
    title: str | None
    model_index: int
    source_path: str
    source_sha256: str
    atoms: tuple[Atom3D, ...]
    residues: tuple[Residue3D, ...]
    chain_ids: tuple[str, ...]
    sequence_by_chain: tuple[str, ...]
    parser: str


@dataclass(frozen=True)
class ChainGeometry(_SerializableResult):
    chain_id: str
    residue_count: int
    contour_length_angstrom: float
    end_to_end_distance_angstrom: float
    anchor_atom: str


@dataclass(frozen=True)
class ShapeGeometry(_SerializableResult):
    center_of_mass_angstrom: tuple[float, float, float]
    radius_of_gyration_angstrom: float
    principal_moments_dalton_angstrom2: tuple[float, float, float]
    gyration_eigenvalues_angstrom2: tuple[float, float, float]
    sphericity: float
    eccentricity: float
    asphericity_angstrom2: float
    relative_shape_anisotropy: float


@dataclass(frozen=True)
class HydrogenBondGeometry(_SerializableResult):
    explicit_hydrogen_bond_count: int | None
    geometric_base_n_o_contact_count: int
    heavy_atom_distance_cutoff_angstrom: float
    method: str
    limitation: str


@dataclass(frozen=True)
class BackboneTorsion(_SerializableResult):
    chain_id: str
    residue_number: int
    insertion_code: str
    alpha_degree: float | None
    beta_degree: float | None
    gamma_degree: float | None
    delta_degree: float | None
    epsilon_degree: float | None
    zeta_degree: float | None


@dataclass(frozen=True)
class HelixGeometry(_SerializableResult):
    chain_ids: tuple[str, str]
    base_pair_count: int
    helical_axis_unit_vector: tuple[float, float, float]
    helical_contour_length_angstrom: float
    mean_rise_angstrom: float
    mean_twist_degree: float
    helical_pitch_angstrom: float
    base_pairs_per_turn: float
    bending_angle_degree: float
    mean_curvature_inverse_angstrom: float
    method: str
    limitation: str


@dataclass(frozen=True)
class Conditional3DFeature(_SerializableResult):
    name: str
    status: str
    required_input_or_backend: str
    reason: str


@dataclass(frozen=True)
class Structure3DAnalysisResult(_SerializableResult):
    pdb_id: str | None
    model_index: int
    atom_count: int
    residue_count: int
    chain_count: int
    chain_geometry: tuple[ChainGeometry, ...]
    shape: ShapeGeometry
    solvent_accessible_surface_area_angstrom2: float
    molecular_volume_angstrom3: float
    sasa_probe_radius_angstrom: float
    sasa_points_per_atom: int
    volume_grid_spacing_angstrom: float
    hydrogen_bonds: HydrogenBondGeometry
    backbone_torsions: tuple[BackboneTorsion, ...]
    helix: HelixGeometry | None
    conditional_features: tuple[Conditional3DFeature, ...]
    method: str
    provenance: Provenance


@dataclass(frozen=True)
class ResidueFlexibility(_SerializableResult):
    chain_id: str
    residue_number: int
    insertion_code: str
    rmsf_angstrom: float
    common_atom_count: int


@dataclass(frozen=True)
class EnsembleFlexibilityResult(_SerializableResult):
    model_count: int
    common_atom_count: int
    mean_atomic_rmsf_angstrom: float
    max_atomic_rmsf_angstrom: float
    residue_flexibility: tuple[ResidueFlexibility, ...]
    method: str
    limitation: str
    provenance: Provenance


@dataclass(frozen=True)
class BasePairStepParameters(_SerializableResult):
    index: int
    base_pair: str
    shear_angstrom: float
    stretch_angstrom: float
    stagger_angstrom: float
    buckle_degree: float
    propeller_degree: float
    opening_degree: float
    shift_angstrom: float
    slide_angstrom: float
    rise_angstrom: float
    tilt_degree: float
    roll_degree: float
    twist_degree: float


@dataclass(frozen=True)
class ThreeDNAParameterResult(_SerializableResult):
    source_path: str
    base_pair_count: int
    parameters: tuple[BasePairStepParameters, ...]
    mean_rise_angstrom: float | None
    mean_twist_degree: float | None
    helical_pitch_angstrom: float | None
    base_pairs_per_turn: float | None
    method: str
    provenance: Provenance


@dataclass(frozen=True)
class DSSRSummaryResult(_SerializableResult):
    source_path: str
    program_version: str
    nucleotide_count: int
    pair_count: int
    helix_count: int
    stem_count: int
    hairpin_count: int
    hydrogen_bond_count: int
    backbone_torsion_record_count: int
    method: str
    provenance: Provenance


__all__ = [
    "Atom3D",
    "BackboneTorsion",
    "BasePairStepParameters",
    "ChainGeometry",
    "Conditional3DFeature",
    "DNA3DStructure",
    "DSSRSummaryResult",
    "EnsembleFlexibilityResult",
    "HelixGeometry",
    "HydrogenBondGeometry",
    "Residue3D",
    "ResidueFlexibility",
    "ShapeGeometry",
    "Structure3DAnalysisResult",
    "ThreeDNAParameterResult",
]
