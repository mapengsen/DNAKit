"""Deterministic geometry descriptors derived from explicit DNA coordinates."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Callable, Iterable
from contextlib import suppress
from itertools import combinations, pairwise
from typing import TypeAlias

from dnakit.exceptions import ConfigurationError

from ._shared import native_structure_provenance
from .results import (
    Atom3D,
    BackboneTorsion,
    ChainGeometry,
    Conditional3DFeature,
    DNA3DStructure,
    EnsembleFlexibilityResult,
    HelixGeometry,
    HydrogenBondGeometry,
    Residue3D,
    ResidueFlexibility,
    ShapeGeometry,
    Structure3DAnalysisResult,
)

Vector3: TypeAlias = tuple[float, float, float]
StructureProgress: TypeAlias = Callable[[str, int, int], None]

_ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
}
_VDW_RADIUS = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "S": 1.80}
_BASE_HETERO_ATOMS = frozenset({"N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"})
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _add(first: Vector3, second: Vector3) -> Vector3:
    return (first[0] + second[0], first[1] + second[1], first[2] + second[2])


def _subtract(first: Vector3, second: Vector3) -> Vector3:
    return (first[0] - second[0], first[1] - second[1], first[2] - second[2])


def _scale(value: Vector3, factor: float) -> Vector3:
    return (value[0] * factor, value[1] * factor, value[2] * factor)


def _dot(first: Vector3, second: Vector3) -> float:
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


def _cross(first: Vector3, second: Vector3) -> Vector3:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _norm(value: Vector3) -> float:
    return math.sqrt(_dot(value, value))


def _unit(value: Vector3) -> Vector3:
    length = _norm(value)
    if length <= 1e-12:
        raise ConfigurationError(
            "A required geometric vector has zero length.", code="DEGENERATE_3D_GEOMETRY"
        )
    return _scale(value, 1.0 / length)


def _distance(first: Vector3, second: Vector3) -> float:
    return _norm(_subtract(first, second))


def _angle(first: Vector3, second: Vector3) -> float:
    denominator = _norm(first) * _norm(second)
    if denominator <= 1e-12:
        return 0.0
    cosine = min(1.0, max(-1.0, _dot(first, second) / denominator))
    return math.degrees(math.acos(cosine))


def _atom(residue: Residue3D, *names: str) -> Atom3D | None:
    mapping = {item.name.replace("*", "'"): item for item in residue.atoms}
    return next(
        (mapping[name.replace("*", "'")] for name in names if name.replace("*", "'") in mapping),
        None,
    )


def _dihedral(
    first: Atom3D | None, second: Atom3D | None, third: Atom3D | None, fourth: Atom3D | None
) -> float | None:
    if first is None or second is None or third is None or fourth is None:
        return None
    b0 = _scale(_subtract(second.coordinates, first.coordinates), -1.0)
    b1 = _unit(_subtract(third.coordinates, second.coordinates))
    b2 = _subtract(fourth.coordinates, third.coordinates)
    v = _subtract(b0, _scale(b1, _dot(b0, b1)))
    w = _subtract(b2, _scale(b1, _dot(b2, b1)))
    if _norm(v) <= 1e-12 or _norm(w) <= 1e-12:
        return None
    return math.degrees(math.atan2(_dot(_cross(b1, v), w), _dot(v, w)))


def _eigenvalues_symmetric(matrix: tuple[Vector3, Vector3, Vector3]) -> tuple[float, float, float]:
    values = [list(row) for row in matrix]
    for _ in range(64):
        first, second = max(
            ((0, 1), (0, 2), (1, 2)), key=lambda pair: abs(values[pair[0]][pair[1]])
        )
        off_diagonal = values[first][second]
        if abs(off_diagonal) < 1e-12:
            break
        angle = 0.5 * math.atan2(2.0 * off_diagonal, values[second][second] - values[first][first])
        cosine, sine = math.cos(angle), math.sin(angle)
        old_first = values[first][first]
        old_second = values[second][second]
        values[first][first] = (
            cosine * cosine * old_first
            - 2.0 * sine * cosine * off_diagonal
            + sine * sine * old_second
        )
        values[second][second] = (
            sine * sine * old_first
            + 2.0 * sine * cosine * off_diagonal
            + cosine * cosine * old_second
        )
        values[first][second] = values[second][first] = 0.0
        for index in range(3):
            if index in {first, second}:
                continue
            old_index_first = values[index][first]
            old_index_second = values[index][second]
            values[index][first] = values[first][index] = (
                cosine * old_index_first - sine * old_index_second
            )
            values[index][second] = values[second][index] = (
                sine * old_index_first + cosine * old_index_second
            )
    result = sorted(
        (max(0.0, values[index][index]) for index in range(3)),
        reverse=True,
    )
    return (result[0], result[1], result[2])


def _shape(atoms: tuple[Atom3D, ...]) -> ShapeGeometry:
    masses = tuple(_ATOMIC_MASS.get(atom.element, 12.011) for atom in atoms)
    total_mass = math.fsum(masses)
    center: Vector3 = (
        math.fsum(mass * atom.coordinates[0] for mass, atom in zip(masses, atoms, strict=True))
        / total_mass,
        math.fsum(mass * atom.coordinates[1] for mass, atom in zip(masses, atoms, strict=True))
        / total_mass,
        math.fsum(mass * atom.coordinates[2] for mass, atom in zip(masses, atoms, strict=True))
        / total_mass,
    )
    centered = tuple(_subtract(atom.coordinates, center) for atom in atoms)
    radius_squared = (
        math.fsum(
            mass * _dot(coordinate, coordinate)
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    xx = (
        math.fsum(
            mass * coordinate[0] * coordinate[0]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    yy = (
        math.fsum(
            mass * coordinate[1] * coordinate[1]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    zz = (
        math.fsum(
            mass * coordinate[2] * coordinate[2]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    xy = (
        math.fsum(
            mass * coordinate[0] * coordinate[1]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    xz = (
        math.fsum(
            mass * coordinate[0] * coordinate[2]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    yz = (
        math.fsum(
            mass * coordinate[1] * coordinate[2]
            for mass, coordinate in zip(masses, centered, strict=True)
        )
        / total_mass
    )
    gyration = _eigenvalues_symmetric(((xx, xy, xz), (xy, yy, yz), (xz, yz, zz)))
    inertia = _eigenvalues_symmetric(
        (
            ((yy + zz) * total_mass, -xy * total_mass, -xz * total_mass),
            (-xy * total_mass, (xx + zz) * total_mass, -yz * total_mass),
            (-xz * total_mass, -yz * total_mass, (xx + yy) * total_mass),
        )
    )
    eigen_sum = math.fsum(gyration)
    sphericity = (
        0.0
        if eigen_sum == 0.0
        else 3.0 * (gyration[0] * gyration[1] * gyration[2]) ** (1.0 / 3.0) / eigen_sum
    )
    eccentricity = (
        0.0 if gyration[0] == 0.0 else math.sqrt(max(0.0, 1.0 - gyration[2] / gyration[0]))
    )
    asphericity = gyration[0] - 0.5 * (gyration[1] + gyration[2])
    pair_product = gyration[0] * gyration[1] + gyration[1] * gyration[2] + gyration[2] * gyration[0]
    anisotropy = 0.0 if eigen_sum == 0.0 else 1.0 - 3.0 * pair_product / (eigen_sum**2)
    return ShapeGeometry(
        center_of_mass_angstrom=center,
        radius_of_gyration_angstrom=math.sqrt(max(0.0, radius_squared)),
        principal_moments_dalton_angstrom2=inertia,
        gyration_eigenvalues_angstrom2=gyration,
        sphericity=min(1.0, max(0.0, sphericity)),
        eccentricity=min(1.0, max(0.0, eccentricity)),
        asphericity_angstrom2=asphericity,
        relative_shape_anisotropy=min(1.0, max(0.0, anisotropy)),
    )


def _chain_geometry(structure: DNA3DStructure) -> tuple[ChainGeometry, ...]:
    results: list[ChainGeometry] = []
    for chain_id in structure.chain_ids:
        residues = tuple(item for item in structure.residues if item.chain_id == chain_id)
        anchor_name = next(
            (
                name
                for name in ("C4'", "C1'", "P")
                if all(_atom(residue, name) is not None for residue in residues)
            ),
            None,
        )
        if anchor_name is None:
            continue
        anchors = tuple(
            item.coordinates
            for residue in residues
            if (item := _atom(residue, anchor_name)) is not None
        )
        contour = math.fsum(_distance(first, second) for first, second in pairwise(anchors))
        results.append(
            ChainGeometry(
                chain_id=chain_id,
                residue_count=len(residues),
                contour_length_angstrom=contour,
                end_to_end_distance_angstrom=(
                    0.0 if len(anchors) < 2 else _distance(anchors[0], anchors[-1])
                ),
                anchor_atom=anchor_name,
            )
        )
    return tuple(results)


def _sphere_points(count: int) -> tuple[Vector3, ...]:
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    return tuple(
        (
            math.sqrt(max(0.0, 1.0 - (1.0 - 2.0 * (index + 0.5) / count) ** 2))
            * math.cos(index * golden_angle),
            math.sqrt(max(0.0, 1.0 - (1.0 - 2.0 * (index + 0.5) / count) ** 2))
            * math.sin(index * golden_angle),
            1.0 - 2.0 * (index + 0.5) / count,
        )
        for index in range(count)
    )


def _emit(progress: StructureProgress | None, stage: str, completed: int, total: int) -> None:
    if progress is not None:
        with suppress(Exception):
            progress(stage, completed, total)


def _sasa(
    atoms: tuple[Atom3D, ...],
    *,
    probe_radius: float,
    points_per_atom: int,
    progress: StructureProgress | None,
) -> float:
    expanded = tuple(_VDW_RADIUS.get(atom.element, 1.70) + probe_radius for atom in atoms)
    cell_size = 2.0 * max(expanded)
    grid: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for index, atom in enumerate(atoms):
        cell = (
            math.floor(atom.coordinates[0] / cell_size),
            math.floor(atom.coordinates[1] / cell_size),
            math.floor(atom.coordinates[2] / cell_size),
        )
        grid[cell].append(index)
    unit_points = _sphere_points(points_per_atom)
    area = 0.0
    for atom_index, (atom, radius) in enumerate(zip(atoms, expanded, strict=True), start=1):
        exposed = 0
        for unit_point in unit_points:
            point = _add(atom.coordinates, _scale(unit_point, radius))
            cell = (
                math.floor(point[0] / cell_size),
                math.floor(point[1] / cell_size),
                math.floor(point[2] / cell_size),
            )
            occluded = False
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for other_index in grid.get((cell[0] + dx, cell[1] + dy, cell[2] + dz), ()):
                            if other_index == atom_index - 1:
                                continue
                            if (
                                _distance(point, atoms[other_index].coordinates)
                                < expanded[other_index]
                            ):
                                occluded = True
                                break
                        if occluded:
                            break
                    if occluded:
                        break
                if occluded:
                    break
            if not occluded:
                exposed += 1
        area += 4.0 * math.pi * radius * radius * exposed / points_per_atom
        _emit(progress, "sasa", atom_index, len(atoms))
    return area


def _volume(
    atoms: tuple[Atom3D, ...],
    *,
    spacing: float,
    progress: StructureProgress | None,
) -> float:
    occupied: set[tuple[int, int, int]] = set()
    for atom_index, atom in enumerate(atoms, start=1):
        radius = _VDW_RADIUS.get(atom.element, 1.70)
        lower = tuple(math.floor((value - radius) / spacing) for value in atom.coordinates)
        upper = tuple(math.ceil((value + radius) / spacing) for value in atom.coordinates)
        radius_squared = radius * radius
        for x_index in range(lower[0], upper[0] + 1):
            for y_index in range(lower[1], upper[1] + 1):
                for z_index in range(lower[2], upper[2] + 1):
                    center = (
                        (x_index + 0.5) * spacing,
                        (y_index + 0.5) * spacing,
                        (z_index + 0.5) * spacing,
                    )
                    if (
                        _dot(
                            _subtract(center, atom.coordinates), _subtract(center, atom.coordinates)
                        )
                        <= radius_squared
                    ):
                        occupied.add((x_index, y_index, z_index))
        _emit(progress, "volume", atom_index, len(atoms))
    return len(occupied) * spacing**3


def _hydrogen_bonds(structure: DNA3DStructure) -> HydrogenBondGeometry:
    def residue_key(atom: Atom3D) -> tuple[str, int, str]:
        return (atom.chain_id, atom.residue_number, atom.insertion_code)

    base_hetero = tuple(
        atom
        for atom in structure.atoms
        if atom.name in _BASE_HETERO_ATOMS and atom.element in {"N", "O"}
    )
    contacts = 0
    for first, second in combinations(base_hetero, 2):
        if residue_key(first) == residue_key(second):
            continue
        if (
            first.chain_id == second.chain_id
            and abs(first.residue_number - second.residue_number) <= 1
        ):
            continue
        distance = _distance(first.coordinates, second.coordinates)
        if 2.2 <= distance <= 3.5:
            contacts += 1
    hydrogens = tuple(atom for atom in structure.atoms if atom.element == "H")
    explicit_count: int | None = None
    if hydrogens:
        explicit = 0
        for hydrogen in hydrogens:
            donors = tuple(
                atom
                for atom in base_hetero
                if residue_key(atom) == residue_key(hydrogen)
                and _distance(atom.coordinates, hydrogen.coordinates) <= 1.3
            )
            for donor in donors:
                for acceptor in base_hetero:
                    if residue_key(acceptor) == residue_key(donor):
                        continue
                    if (
                        _distance(hydrogen.coordinates, acceptor.coordinates) <= 2.5
                        and _angle(
                            _subtract(donor.coordinates, hydrogen.coordinates),
                            _subtract(acceptor.coordinates, hydrogen.coordinates),
                        )
                        >= 120.0
                    ):
                        explicit += 1
        explicit_count = explicit
    return HydrogenBondGeometry(
        explicit_hydrogen_bond_count=explicit_count,
        geometric_base_n_o_contact_count=contacts,
        heavy_atom_distance_cutoff_angstrom=3.5,
        method="explicit-D-H-A-geometry-or-base-N-O-contact-screen-v1",
        limitation=(
            "Actual hydrogen bonds require explicit hydrogen atoms and protonation; when absent, "
            "only geometric base N/O contacts are reported and must not be called actual bonds."
        ),
    )


def _torsions(structure: DNA3DStructure) -> tuple[BackboneTorsion, ...]:
    results: list[BackboneTorsion] = []
    for chain_id in structure.chain_ids:
        residues = tuple(item for item in structure.residues if item.chain_id == chain_id)
        for index, residue in enumerate(residues):
            previous = None if index == 0 else residues[index - 1]
            following = None if index + 1 == len(residues) else residues[index + 1]
            results.append(
                BackboneTorsion(
                    chain_id=chain_id,
                    residue_number=residue.number,
                    insertion_code=residue.insertion_code,
                    alpha_degree=_dihedral(
                        None if previous is None else _atom(previous, "O3'"),
                        _atom(residue, "P"),
                        _atom(residue, "O5'"),
                        _atom(residue, "C5'"),
                    ),
                    beta_degree=_dihedral(
                        _atom(residue, "P"),
                        _atom(residue, "O5'"),
                        _atom(residue, "C5'"),
                        _atom(residue, "C4'"),
                    ),
                    gamma_degree=_dihedral(
                        _atom(residue, "O5'"),
                        _atom(residue, "C5'"),
                        _atom(residue, "C4'"),
                        _atom(residue, "C3'"),
                    ),
                    delta_degree=_dihedral(
                        _atom(residue, "C5'"),
                        _atom(residue, "C4'"),
                        _atom(residue, "C3'"),
                        _atom(residue, "O3'"),
                    ),
                    epsilon_degree=_dihedral(
                        _atom(residue, "C4'"),
                        _atom(residue, "C3'"),
                        _atom(residue, "O3'"),
                        None if following is None else _atom(following, "P"),
                    ),
                    zeta_degree=_dihedral(
                        _atom(residue, "C3'"),
                        _atom(residue, "O3'"),
                        None if following is None else _atom(following, "P"),
                        None if following is None else _atom(following, "O5'"),
                    ),
                )
            )
    return tuple(results)


def _project_perpendicular(value: Vector3, axis: Vector3) -> Vector3:
    return _subtract(value, _scale(axis, _dot(value, axis)))


def _helix(structure: DNA3DStructure) -> HelixGeometry | None:
    chains = {
        chain_id: tuple(item for item in structure.residues if item.chain_id == chain_id)
        for chain_id in structure.chain_ids
    }
    selected: tuple[str, str] | None = None
    for first, second in combinations(structure.chain_ids, 2):
        first_sequence = "".join(item.base for item in chains[first])
        second_sequence = "".join(item.base for item in chains[second])
        if (
            len(first_sequence) >= 2
            and second_sequence == first_sequence.translate(_COMPLEMENT)[::-1]
        ):
            selected = (first, second)
            break
    if selected is None:
        return None
    first_residues = chains[selected[0]]
    second_residues = chains[selected[1]][::-1]
    anchors: list[tuple[Vector3, Vector3]] = []
    for first_residue, second_residue in zip(first_residues, second_residues, strict=True):
        first_atom = _atom(first_residue, "C1'")
        second_atom = _atom(second_residue, "C1'")
        if first_atom is None or second_atom is None:
            return None
        anchors.append((first_atom.coordinates, second_atom.coordinates))
    centers = tuple(_scale(_add(first, second), 0.5) for first, second in anchors)
    axis = _unit(_subtract(centers[-1], centers[0]))
    orientations = tuple(
        _project_perpendicular(_subtract(second, first), axis) for first, second in anchors
    )
    if any(_norm(value) <= 1e-12 for value in orientations):
        return None
    twist_values = tuple(abs(_angle(first, second)) for first, second in pairwise(orientations))
    rise_values = tuple(
        abs(_dot(_subtract(second, first), axis)) for first, second in pairwise(centers)
    )
    mean_twist = math.fsum(twist_values) / len(twist_values)
    mean_rise = math.fsum(rise_values) / len(rise_values)
    contour = math.fsum(_distance(first, second) for first, second in pairwise(centers))
    first_direction = _subtract(centers[min(2, len(centers) - 1)], centers[0])
    last_direction = _subtract(centers[-1], centers[max(0, len(centers) - 3)])
    bending = _angle(first_direction, last_direction)
    pairs_per_turn = 0.0 if mean_twist == 0.0 else 360.0 / mean_twist
    return HelixGeometry(
        chain_ids=selected,
        base_pair_count=len(centers),
        helical_axis_unit_vector=axis,
        helical_contour_length_angstrom=contour,
        mean_rise_angstrom=mean_rise,
        mean_twist_degree=mean_twist,
        helical_pitch_angstrom=mean_rise * pairs_per_turn,
        base_pairs_per_turn=pairs_per_turn,
        bending_angle_degree=bending,
        mean_curvature_inverse_angstrom=(
            0.0 if contour == 0.0 else math.radians(bending) / contour
        ),
        method="paired-C1-prime-centers-global-axis-v1",
        limitation=(
            "Approximate global-axis descriptor for two equal reverse-complement chains; "
            "use 3DNA/DSSR or Curves+ for standard local base-pair reference frames."
        ),
    )


def analyze_structure(
    structure: DNA3DStructure,
    *,
    sasa_probe_radius_angstrom: float = 1.4,
    sasa_points_per_atom: int = 96,
    volume_grid_spacing_angstrom: float = 0.75,
    progress: StructureProgress | None = None,
) -> Structure3DAnalysisResult:
    """Calculate native geometry from one explicit PDB model."""

    if not isinstance(structure, DNA3DStructure):
        raise ConfigurationError("structure must be DNA3DStructure.", code="INVALID_3D_STRUCTURE")
    if not structure.atoms or not structure.residues:
        raise ConfigurationError("structure is empty.", code="EMPTY_3D_STRUCTURE")
    if (
        isinstance(sasa_probe_radius_angstrom, bool)
        or not isinstance(sasa_probe_radius_angstrom, (int, float))
        or not math.isfinite(sasa_probe_radius_angstrom)
        or not 0.0 <= sasa_probe_radius_angstrom <= 5.0
    ):
        raise ConfigurationError(
            "sasa_probe_radius_angstrom must be in [0, 5].", code="INVALID_SASA_CONFIG"
        )
    if (
        isinstance(sasa_points_per_atom, bool)
        or not isinstance(sasa_points_per_atom, int)
        or not 24 <= sasa_points_per_atom <= 4_096
    ):
        raise ConfigurationError(
            "sasa_points_per_atom must be in [24, 4096].", code="INVALID_SASA_CONFIG"
        )
    if (
        isinstance(volume_grid_spacing_angstrom, bool)
        or not isinstance(volume_grid_spacing_angstrom, (int, float))
        or not math.isfinite(volume_grid_spacing_angstrom)
        or not 0.25 <= volume_grid_spacing_angstrom <= 5.0
    ):
        raise ConfigurationError(
            "volume_grid_spacing_angstrom must be in [0.25, 5].",
            code="INVALID_VOLUME_CONFIG",
        )
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.", code="INVALID_3D_PROGRESS")
    conditional = (
        Conditional3DFeature(
            "major_minor_groove_width_depth",
            "conditional",
            "3DNA/DSSR or Curves+ output",
            "Standard groove geometry requires local base-pair frames and phosphate tracing.",
        ),
        Conditional3DFeature(
            "base_stacking_area",
            "conditional",
            "3DNA/DSSR stacking analysis",
            "Projected ring overlap depends on standard base frames and atom classification.",
        ),
        Conditional3DFeature(
            "electrostatic_potential_charge_distribution",
            "conditional",
            "PQR charges plus an electrostatics backend such as APBS",
            "A PDB coordinate file does not contain force-field charges or solvent electrostatics.",
        ),
        Conditional3DFeature(
            "persistence_length_and_mechanical_stiffness",
            "conditional",
            "trajectory or structural ensemble with a declared mechanical model",
            "A single conformation cannot identify persistence, bend, twist, or stretch moduli.",
        ),
        Conditional3DFeature(
            "standard_base_pair_and_step_parameters",
            "conditional",
            "3DNA/DSSR bp_step.par output",
            "Use read_3dna_bp_step() to preserve the standard rigid-body definitions.",
        ),
    )
    return Structure3DAnalysisResult(
        pdb_id=structure.pdb_id,
        model_index=structure.model_index,
        atom_count=len(structure.atoms),
        residue_count=len(structure.residues),
        chain_count=len(structure.chain_ids),
        chain_geometry=_chain_geometry(structure),
        shape=_shape(structure.atoms),
        solvent_accessible_surface_area_angstrom2=_sasa(
            structure.atoms,
            probe_radius=float(sasa_probe_radius_angstrom),
            points_per_atom=sasa_points_per_atom,
            progress=progress,
        ),
        molecular_volume_angstrom3=_volume(
            structure.atoms,
            spacing=float(volume_grid_spacing_angstrom),
            progress=progress,
        ),
        sasa_probe_radius_angstrom=float(sasa_probe_radius_angstrom),
        sasa_points_per_atom=sasa_points_per_atom,
        volume_grid_spacing_angstrom=float(volume_grid_spacing_angstrom),
        hydrogen_bonds=_hydrogen_bonds(structure),
        backbone_torsions=_torsions(structure),
        helix=_helix(structure),
        conditional_features=conditional,
        method="dnakit-explicit-coordinate-geometry-v1",
        provenance=native_structure_provenance(),
    )


def analyze_ensemble_flexibility(
    structures: Iterable[DNA3DStructure],
) -> EnsembleFlexibilityResult:
    """Calculate translation-centered RMSF for a pre-oriented PDB model ensemble."""

    if not isinstance(structures, Iterable):
        raise ConfigurationError("structures must be an iterable.", code="INVALID_3D_ENSEMBLE")
    models = tuple(structures)
    if not 2 <= len(models) <= 10_000 or any(
        not isinstance(item, DNA3DStructure) for item in models
    ):
        raise ConfigurationError(
            "A 3D ensemble requires 2-10000 DNA3DStructure models.",
            code="INVALID_3D_ENSEMBLE",
        )

    def atom_map(model: DNA3DStructure) -> dict[tuple[str, int, str, str, str], Atom3D]:
        return {
            (
                atom.chain_id,
                atom.residue_number,
                atom.insertion_code,
                atom.name,
                atom.element,
            ): atom
            for atom in model.atoms
        }

    maps = tuple(atom_map(model) for model in models)
    common = set(maps[0])
    for mapping in maps[1:]:
        common.intersection_update(mapping)
    ordered_keys = tuple(sorted(common))
    if len(ordered_keys) < 3:
        raise ConfigurationError(
            "Ensemble models have fewer than three common atoms.",
            code="INSUFFICIENT_COMMON_ENSEMBLE_ATOMS",
        )
    centered_models: list[dict[tuple[str, int, str, str, str], Vector3]] = []
    for mapping in maps:
        centroid: Vector3 = (
            math.fsum(mapping[key].coordinates[0] for key in ordered_keys) / len(ordered_keys),
            math.fsum(mapping[key].coordinates[1] for key in ordered_keys) / len(ordered_keys),
            math.fsum(mapping[key].coordinates[2] for key in ordered_keys) / len(ordered_keys),
        )
        centered_models.append(
            {key: _subtract(mapping[key].coordinates, centroid) for key in ordered_keys}
        )
    rmsf_by_key: dict[tuple[str, int, str, str, str], float] = {}
    for key in ordered_keys:
        mean: Vector3 = (
            math.fsum(model[key][0] for model in centered_models) / len(centered_models),
            math.fsum(model[key][1] for model in centered_models) / len(centered_models),
            math.fsum(model[key][2] for model in centered_models) / len(centered_models),
        )
        rmsf_by_key[key] = math.sqrt(
            math.fsum(
                _dot(_subtract(model[key], mean), _subtract(model[key], mean))
                for model in centered_models
            )
            / len(centered_models)
        )
    grouped: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    for key, value in rmsf_by_key.items():
        grouped[key[:3]].append(value)
    residue_flexibility = tuple(
        ResidueFlexibility(
            chain_id=key[0],
            residue_number=key[1],
            insertion_code=key[2],
            rmsf_angstrom=math.sqrt(math.fsum(item * item for item in values) / len(values)),
            common_atom_count=len(values),
        )
        for key, values in sorted(grouped.items())
    )
    atomic_values = tuple(rmsf_by_key.values())
    return EnsembleFlexibilityResult(
        model_count=len(models),
        common_atom_count=len(ordered_keys),
        mean_atomic_rmsf_angstrom=math.fsum(atomic_values) / len(atomic_values),
        max_atomic_rmsf_angstrom=max(atomic_values),
        residue_flexibility=residue_flexibility,
        method="common-atom-translation-centered-rmsf-v1",
        limitation=(
            "Models must already share a common rotational frame; DNAKit removes translation "
            "but does not perform Kabsch rotational superposition or infer force constants."
        ),
        provenance=native_structure_provenance(),
    )


__all__ = [
    "StructureProgress",
    "analyze_ensemble_flexibility",
    "analyze_structure",
]
