"""Bounded parser for canonical DNA atoms in legacy PDB coordinate files."""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

from .results import Atom3D, DNA3DStructure, Residue3D

_DNA_RESIDUES = {
    "DA": "A",
    "DC": "C",
    "DG": "G",
    "DT": "T",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "ADE": "A",
    "CYT": "C",
    "GUA": "G",
    "THY": "T",
}
_MAX_PDB_BYTES = 100_000_000
_MAX_ATOMS = 1_000_000
_MAX_MODELS = 10_000


def _integer(text: str, field: str, line_number: int) -> int:
    try:
        return int(text.strip())
    except ValueError as exc:
        raise InputFormatError(
            f"PDB {field} is not an integer.",
            code="INVALID_PDB_FIELD",
            context={"field": field, "line_number": line_number, "value": text},
        ) from exc


def _number(text: str, field: str, line_number: int, *, default: float | None = None) -> float:
    if not text.strip() and default is not None:
        return default
    try:
        value = float(text.strip())
    except ValueError as exc:
        raise InputFormatError(
            f"PDB {field} is not numeric.",
            code="INVALID_PDB_FIELD",
            context={"field": field, "line_number": line_number, "value": text},
        ) from exc
    if not math.isfinite(value):
        raise InputFormatError(
            f"PDB {field} must be finite.",
            code="INVALID_PDB_FIELD",
            context={"field": field, "line_number": line_number},
        )
    return value


def _element(atom_name: str, explicit: str) -> str:
    if explicit.strip():
        return explicit.strip().upper()
    letters = "".join(character for character in atom_name if character.isalpha())
    return letters[:1].upper() or "X"


def _structure(
    *,
    path: Path,
    digest: str,
    pdb_id: str | None,
    title: str | None,
    model_index: int,
    atoms: tuple[Atom3D, ...],
) -> DNA3DStructure:
    residue_atoms: dict[tuple[str, int, str, str], list[Atom3D]] = {}
    residue_order: list[tuple[str, int, str, str]] = []
    for atom in atoms:
        key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)
        if key not in residue_atoms:
            residue_atoms[key] = []
            residue_order.append(key)
        residue_atoms[key].append(atom)
    residues = tuple(
        Residue3D(
            name=residue_name,
            base=_DNA_RESIDUES[residue_name],
            chain_id=chain_id,
            number=number,
            insertion_code=insertion_code,
            atoms=tuple(residue_atoms[(chain_id, number, insertion_code, residue_name)]),
        )
        for chain_id, number, insertion_code, residue_name in residue_order
    )
    chain_ids = tuple(dict.fromkeys(residue.chain_id for residue in residues))
    sequence_by_chain = tuple(
        "".join(residue.base for residue in residues if residue.chain_id == chain_id)
        for chain_id in chain_ids
    )
    return DNA3DStructure(
        pdb_id=pdb_id,
        title=title,
        model_index=model_index,
        source_path=str(path),
        source_sha256=digest,
        atoms=atoms,
        residues=residues,
        chain_ids=chain_ids,
        sequence_by_chain=sequence_by_chain,
        parser="dnakit-bounded-legacy-pdb-parser-v1",
    )


def load_pdb_ensemble(path: str | Path) -> tuple[DNA3DStructure, ...]:
    """Load all PDB MODEL records, retaining canonical DNA ATOM records only."""

    if not isinstance(path, (str, Path)):
        raise ConfigurationError("path must be text or Path.", code="INVALID_PDB_PATH")
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise InputFormatError(
            "PDB path is not a readable file.",
            code="PDB_FILE_NOT_FOUND",
            context={"path": str(source)},
        )
    size = source.stat().st_size
    if not 1 <= size <= _MAX_PDB_BYTES:
        raise InputFormatError(
            f"PDB file size must be in [1, {_MAX_PDB_BYTES}] bytes.",
            code="PDB_FILE_SIZE_LIMIT",
            context={"size_bytes": size},
        )
    payload = source.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "Legacy PDB files must contain ASCII text.", code="INVALID_PDB_ENCODING"
        ) from exc
    pdb_id: str | None = None
    title_parts: list[str] = []
    current_model = 1
    explicit_models = False
    selected_atoms: dict[int, dict[tuple[str, int, str, str, str], tuple[int, float, Atom3D]]] = (
        defaultdict(dict)
    )
    atom_records = 0
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.ljust(80)
        record = line[:6]
        if record == "HEADER":
            candidate = line[62:66].strip().upper()
            if candidate:
                pdb_id = candidate
        elif record == "TITLE ":
            part = line[10:80].strip()
            if part:
                title_parts.append(part)
        elif record == "MODEL ":
            explicit_models = True
            current_model = _integer(line[10:14], "model", line_number)
            if not 1 <= current_model <= _MAX_MODELS:
                raise InputFormatError(
                    "PDB model index exceeds the safety limit.",
                    code="PDB_MODEL_LIMIT",
                    context={"model_index": current_model},
                )
        elif record == "ATOM  ":
            residue_name = line[17:20].strip().upper()
            if residue_name not in _DNA_RESIDUES:
                continue
            atom_records += 1
            if atom_records > _MAX_ATOMS:
                raise InputFormatError(
                    "PDB DNA atom count exceeds the safety limit.",
                    code="PDB_ATOM_LIMIT",
                )
            altloc = line[16].strip()
            if altloc not in {"", "A"}:
                continue
            serial = _integer(line[6:11], "serial", line_number)
            atom_name = line[12:16].strip()
            chain_id = line[21].strip() or "_"
            residue_number = _integer(line[22:26], "residue_number", line_number)
            insertion_code = line[26].strip()
            occupancy = _number(line[54:60], "occupancy", line_number, default=1.0)
            atom = Atom3D(
                serial=serial,
                name=atom_name,
                element=_element(atom_name, line[76:78]),
                residue_name=residue_name,
                chain_id=chain_id,
                residue_number=residue_number,
                insertion_code=insertion_code,
                x_angstrom=_number(line[30:38], "x", line_number),
                y_angstrom=_number(line[38:46], "y", line_number),
                z_angstrom=_number(line[46:54], "z", line_number),
                occupancy=occupancy,
                model_index=current_model,
            )
            key = (chain_id, residue_number, insertion_code, residue_name, atom_name)
            priority = 2 if not altloc else 1
            previous = selected_atoms[current_model].get(key)
            if previous is None or (priority, occupancy) > (previous[0], previous[1]):
                selected_atoms[current_model][key] = (priority, occupancy, atom)
    if not selected_atoms:
        raise InputFormatError(
            "PDB file contains no canonical DNA ATOM records.",
            code="PDB_DNA_ATOMS_NOT_FOUND",
        )
    if explicit_models and len(selected_atoms) > _MAX_MODELS:
        raise InputFormatError("PDB model count exceeds the safety limit.", code="PDB_MODEL_LIMIT")
    title = " ".join(title_parts) or None
    return tuple(
        _structure(
            path=source,
            digest=digest,
            pdb_id=pdb_id,
            title=title,
            model_index=model_index,
            atoms=tuple(value[2] for value in model_atoms.values()),
        )
        for model_index, model_atoms in sorted(selected_atoms.items())
    )


def load_pdb(path: str | Path, *, model_index: int = 1) -> DNA3DStructure:
    """Load one selected PDB model."""

    if isinstance(model_index, bool) or not isinstance(model_index, int) or model_index <= 0:
        raise ConfigurationError(
            "model_index must be a positive integer.", code="INVALID_PDB_MODEL_INDEX"
        )
    structures = load_pdb_ensemble(path)
    for structure in structures:
        if structure.model_index == model_index:
            return structure
    raise InputFormatError(
        "Requested PDB model is absent.",
        code="PDB_MODEL_NOT_FOUND",
        context={
            "model_index": model_index,
            "available_models": tuple(item.model_index for item in structures),
        },
    )


__all__ = ["load_pdb", "load_pdb_ensemble"]
