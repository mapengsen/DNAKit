"""Biotite nucleic-acid structure functions (integration pack P30).

Base pairs, base stacking, dot-bracket notation, pseudoknot order and secondary-structure
drawing for nucleic-acid 3D structures, wrapping ``biotite.structure``. Structures are read from
PDB or PDBx/mmCIF files (``biotite.structure.io.load_structure``). Atom indices in Biotite's
base-pair and stacking results are the index of the first atom of each residue; these functions
also report the chain, residue number and residue name. Biotite is imported only when a function
runs; the tested version is biotite 2.0.0. Drawing needs the ViennaRNA program ``RNAplot``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-biotite"
_PACKAGE = "biotite"
_TESTED_VERSION = "2.0.0"
_STRUCTURE_SUFFIXES = (".pdb", ".cif", ".pdbx", ".bcif", ".mmtf", ".ent")
_EDGES = {0: "invalid", 1: "watson_crick", 2: "hoogsteen", 3: "sugar"}
_ORIENTATIONS = {0: "invalid", 1: "cis", 2: "trans"}
_LAYOUTS = ("radial", "naview", "circular", "rnaturtle", "rnapuzzler")
_MAX_LENGTH = 20_000
_MAX_PAIRS = 100_000


def _provenance(**notes: Any) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {
        "package": _PACKAGE,
        "package_version": installed,
        "tested_version": _TESTED_VERSION,
        **notes,
    }


def _int(value: object, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _load(structure_file: str | Path, model: int, chain_id: str | None) -> Any:
    """Load one model of a structure file, keeping only nucleotides (optionally one chain)."""

    path = Path(structure_file).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"structure_file does not exist: {path}")
    if path.suffix.lower() not in _STRUCTURE_SUFFIXES:
        raise ConfigurationError(f"structure_file must end in one of {list(_STRUCTURE_SUFFIXES)}.")
    model_number = _int(model, "model", minimum=1)
    structure = require("biotite.structure", _EXTRA)
    structure_io = require("biotite.structure.io", _EXTRA)
    try:
        atoms = structure_io.load_structure(str(path), model=model_number)
    except Exception as exc:
        raise ConfigurationError(f"{path.name} could not be read as model {model_number}: {exc}") from exc
    if chain_id is not None:
        if not isinstance(chain_id, str) or not chain_id:
            raise ConfigurationError("chain_id must be a non-empty string.")
        atoms = atoms[atoms.chain_id == chain_id]
    atoms = atoms[structure.filter_nucleotides(atoms)]
    if atoms.array_length() == 0:
        raise ConfigurationError("The structure contains no nucleotides for the chosen model/chain.")
    return atoms


def _residue(atoms: Any, atom_index: int) -> dict[str, Any]:
    return {
        "atom_index": int(atom_index),
        "chain_id": str(atoms.chain_id[atom_index]),
        "res_id": int(atoms.res_id[atom_index]),
        "res_name": str(atoms.res_name[atom_index]),
    }


def _pairs(values: object, length: int | None = None) -> list[tuple[int, int]]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError("base_pairs must be a list of [i, j] index pairs.")
    if len(values) > _MAX_PAIRS:
        raise ConfigurationError(f"base_pairs may hold at most {_MAX_PAIRS} pairs.")
    pairs = []
    for row in values:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 2:
            raise ConfigurationError("each base pair must be [i, j].")
        first, second = _int(row[0], "base pair index", minimum=0), _int(row[1], "base pair index", minimum=0)
        if first == second or (length is not None and max(first, second) >= length):
            raise ConfigurationError(f"invalid base pair [{first}, {second}] for length {length}.")
        pairs.append((first, second))
    flat = [index for pair in pairs for index in pair]
    if len(set(flat)) != len(flat):
        raise ConfigurationError("a base can take part in only one base pair.")
    return pairs


def nucleic_acid_base_pairs(
    structure_file: str,
    *,
    model: int = 1,
    chain_id: str | None = None,
    min_atoms_per_base: int = 3,
    unique: bool = True,
    annotate_edges: bool = False,
    annotate_glycosidic_bond: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the base pairs of a nucleic-acid 3D structure (``biotite.structure.base_pairs``).

    Pairs follow the DSSR criteria. ``annotate_edges=True`` adds the Leontis-Westhof edge of
    each base (``base_pairs_edge``: watson_crick, hoogsteen, sugar); ``annotate_glycosidic_bond``
    adds the cis/trans orientation (``base_pairs_glycosidic_bond``). The edge annotation needs
    hydrogen atoms in the structure (for example NMR models); without them Biotite raises an
    error that is reported as a configuration error.
    """

    atoms = _load(structure_file, model, chain_id)
    edges_wanted = _flag(annotate_edges, "annotate_edges")
    glycosidic_wanted = _flag(annotate_glycosidic_bond, "annotate_glycosidic_bond")
    structure = require("biotite.structure", _EXTRA)
    pairs = structure.base_pairs(
        atoms,
        min_atoms_per_base=_int(min_atoms_per_base, "min_atoms_per_base", minimum=1),
        unique=_flag(unique, "unique"),
    )
    edges = glycosidic = None
    try:
        if edges_wanted and len(pairs):
            edges = structure.base_pairs_edge(atoms, pairs)
        if glycosidic_wanted and len(pairs):
            glycosidic = structure.base_pairs_glycosidic_bond(atoms, pairs)
    except Exception as exc:
        raise ConfigurationError(
            f"Edge or orientation annotation failed (hydrogens are required): {exc}"
        ) from exc
    rows = []
    for number, (first, second) in enumerate(pairs.tolist()):
        row: dict[str, Any] = {"base_1": _residue(atoms, first), "base_2": _residue(atoms, second)}
        if edges is not None:
            row["edge_1"] = _EDGES.get(int(edges[number][0]), "invalid")
            row["edge_2"] = _EDGES.get(int(edges[number][1]), "invalid")
        if glycosidic is not None:
            row["glycosidic_bond"] = _ORIENTATIONS.get(int(glycosidic[number]), "invalid")
        rows.append(row)
    return {
        "status": "completed",
        "operation": "nucleic_acid_base_pairs",
        "pair_count": len(rows),
        "base_pairs": rows,
        "provenance": _provenance(method="base_pairs"),
    }


def nucleic_acid_base_stacking(
    structure_file: str,
    *,
    model: int = 1,
    chain_id: str | None = None,
    min_atoms_per_base: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find stacked neighbouring bases of a nucleic-acid 3D structure (``base_stacking``)."""

    atoms = _load(structure_file, model, chain_id)
    structure = require("biotite.structure", _EXTRA)
    stacking = structure.base_stacking(
        atoms, min_atoms_per_base=_int(min_atoms_per_base, "min_atoms_per_base", minimum=1)
    )
    rows = [
        {"base_1": _residue(atoms, first), "base_2": _residue(atoms, second)}
        for first, second in stacking.tolist()
    ]
    return {
        "status": "completed",
        "operation": "nucleic_acid_base_stacking",
        "stack_count": len(rows),
        "stacked_bases": rows,
        "provenance": _provenance(method="base_stacking"),
    }


def _scores(values: Sequence[int] | None, count: int) -> Any:
    if values is None:
        return None
    numpy = require("numpy", _EXTRA)
    if isinstance(values, (str, bytes)) or len(values) != count:
        raise ConfigurationError("scores must list one integer per base pair.")
    return numpy.asarray([_int(item, "score") for item in values])


def base_pairs_to_dot_bracket(
    base_pairs: Sequence[Sequence[int]],
    length: int,
    *,
    scores: Sequence[int] | None = None,
    max_pseudoknot_order: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write base pairs as dot-bracket notation (``biotite.structure.dot_bracket``).

    ``base_pairs`` are 0-based position pairs ``[i, j]`` of a strand of ``length`` bases.
    Pseudoknots use bracket types ``()``, ``[]``, ``{}``, ``<>`` and letters. ``scores`` rank
    pairs when the pseudoknot assignment is ambiguous; ``max_pseudoknot_order`` caps the
    pseudoknot order. More than one equivalent notation can be returned.
    """

    size = _int(length, "length", minimum=1)
    if size > _MAX_LENGTH:
        raise ConfigurationError(f"length may be at most {_MAX_LENGTH}.")
    pairs = _pairs(base_pairs, size)
    numpy = require("numpy", _EXTRA)
    structure = require("biotite.structure", _EXTRA)
    array = numpy.asarray(pairs, dtype=int).reshape(-1, 2)
    order = None if max_pseudoknot_order is None else _int(max_pseudoknot_order, "max_pseudoknot_order", minimum=0)
    notation = structure.dot_bracket(
        array, size, scores=_scores(scores, len(pairs)), max_pseudoknot_order=order
    )
    return {
        "status": "completed",
        "operation": "base_pairs_to_dot_bracket",
        "length": size,
        "dot_bracket": [str(item) for item in notation],
        "provenance": _provenance(method="dot_bracket"),
    }


def nucleic_acid_dot_bracket(
    structure_file: str,
    *,
    model: int = 1,
    chain_id: str | None = None,
    scores: Sequence[int] | None = None,
    max_pseudoknot_order: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dot-bracket secondary structure of a nucleic-acid 3D structure.

    Finds the base pairs with ``base_pairs`` and writes them with
    ``dot_bracket_from_structure``. Chains of a duplex are concatenated in file order;
    ``chain_id`` restricts the strand. Returns the residue name of every nucleotide used.
    """

    atoms = _load(structure_file, model, chain_id)
    structure = require("biotite.structure", _EXTRA)
    order = None if max_pseudoknot_order is None else _int(max_pseudoknot_order, "max_pseudoknot_order", minimum=0)
    try:
        notation = structure.dot_bracket_from_structure(
            atoms,
            scores=_scores(scores, len(structure.base_pairs(atoms))),
            max_pseudoknot_order=order,
        )
    except Exception as exc:
        raise ConfigurationError(f"The dot-bracket notation could not be built: {exc}") from exc
    starts = structure.get_residue_starts(atoms)
    names = [str(atoms.res_name[index]) for index in starts]
    return {
        "status": "completed",
        "operation": "nucleic_acid_dot_bracket",
        "nucleotide_count": len(names),
        "residue_names": names,
        "dot_bracket": [str(item) for item in notation],
        "provenance": _provenance(method="dot_bracket_from_structure"),
    }


def pseudoknot_order(
    base_pairs: Sequence[Sequence[int]],
    *,
    scores: Sequence[int] | None = None,
    max_pseudoknot_order: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign each base pair its pseudoknot order (``biotite.structure.pseudoknots``).

    Order 0 pairs are nested; order 1 pairs cross them, and so on. Returns one assignment per
    optimal solution (``orders[solution][pair]``); -1 marks pairs above ``max_pseudoknot_order``.
    ``scores`` give an integer weight per pair for choosing among equal solutions.
    """

    pairs = _pairs(base_pairs)
    numpy = require("numpy", _EXTRA)
    structure = require("biotite.structure", _EXTRA)
    if not pairs:
        raise ConfigurationError("base_pairs must contain at least one pair.")
    order = None if max_pseudoknot_order is None else _int(max_pseudoknot_order, "max_pseudoknot_order", minimum=0)
    result = structure.pseudoknots(
        numpy.asarray(pairs, dtype=int),
        scores=_scores(scores, len(pairs)),
        max_pseudoknot_order=order,
    )
    return {
        "status": "completed",
        "operation": "pseudoknot_order",
        "base_pairs": [list(pair) for pair in pairs],
        "solution_count": int(result.shape[0]),
        "orders": [[int(value) for value in row] for row in result.tolist()],
        "provenance": _provenance(method="pseudoknots"),
    }


def plot_nucleotide_secondary_structure(
    bases: str,
    dot_bracket: str,
    output_path: str,
    *,
    layout: str = "naview",
    draw_pseudoknots: bool = True,
    angle: float = 0.0,
    rnaplot_path: str = "RNAplot",
    dpi: int = 150,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a nucleic-acid secondary structure to an image file.

    Uses ``plot_nucleotide_secondary_structure``: the 2D layout comes from ViennaRNA's
    ``RNAplot`` (``layout`` is ``radial``, ``naview``, ``circular``, ``rnaturtle`` or
    ``rnapuzzler``). ``bases`` is the sequence and ``dot_bracket`` the structure of the same
    length. Returns the image path, the base pairs and the layout.
    """

    if not isinstance(bases, str) or not bases.strip():
        raise ConfigurationError("bases must be a non-empty string.")
    text = bases.strip().upper()
    notation = dot_bracket.strip() if isinstance(dot_bracket, str) else ""
    if len(text) != len(notation) or not notation or len(text) > _MAX_LENGTH:
        raise ConfigurationError("bases and dot_bracket must have the same non-zero length.")
    way = layout.lower() if isinstance(layout, str) else ""
    if way not in _LAYOUTS:
        raise ConfigurationError(f"layout must be one of {list(_LAYOUTS)}.")
    out = Path(output_path).expanduser().resolve()
    if out.exists() or not out.parent.is_dir():
        raise ConfigurationError(f"output_path must be a new file in an existing folder: {out}")
    if out.suffix.lower() not in (".png", ".svg", ".pdf"):
        raise ConfigurationError("output_path must end in .png, .svg or .pdf.")
    from dnakit.bioinformatics import _gap_cli

    executable = _gap_cli.program(
        "RNAplot",
        explicit=None if rnaplot_path == "RNAplot" else rnaplot_path,
        hint=(
            "Install ViennaRNA (for example `conda install -c bioconda viennarna`) and put "
            "RNAplot on PATH, or pass rnaplot_path."
        ),
    )
    structure = require("biotite.structure", _EXTRA)
    viennarna = require("biotite.application.viennarna", _EXTRA)
    graphics = require("biotite.structure.graphics", _EXTRA)
    matplotlib = require("matplotlib", _EXTRA)
    matplotlib.use("Agg", force=True)
    pyplot = require("matplotlib.pyplot", _EXTRA)
    try:
        pairs = structure.base_pairs_from_dot_bracket(notation)
    except Exception as exc:
        raise ConfigurationError(f"dot_bracket is not valid: {exc}") from exc
    layout_enum = viennarna.PlotApp.Layout[way.upper()]
    figure, axes = pyplot.subplots(figsize=(8.0, 8.0))
    try:
        graphics.plot_nucleotide_secondary_structure(
            axes,
            list(text),
            pairs,
            len(text),
            layout_type=layout_enum,
            draw_pseudoknots=_flag(draw_pseudoknots, "draw_pseudoknots"),
            angle=float(angle),
            bin_path=executable,
        )
        figure.savefig(out, dpi=_int(dpi, "dpi", minimum=30), bbox_inches="tight")
    except Exception as exc:
        raise ConfigurationError(f"The structure could not be drawn: {exc}") from exc
    finally:
        pyplot.close(figure)
    return {
        "status": "completed",
        "operation": "plot_nucleotide_secondary_structure",
        "image": str(out),
        "layout": way,
        "bases": text,
        "dot_bracket": notation,
        "base_pairs": [[int(a), int(b)] for a, b in pairs.tolist()],
        "provenance": _provenance(method="plot_nucleotide_secondary_structure", rnaplot=executable),
    }


__all__ = [
    "base_pairs_to_dot_bracket",
    "nucleic_acid_base_pairs",
    "nucleic_acid_base_stacking",
    "nucleic_acid_dot_bracket",
    "plot_nucleotide_secondary_structure",
    "pseudoknot_order",
]
_register_public_verbose(globals(), tuple(__all__))
