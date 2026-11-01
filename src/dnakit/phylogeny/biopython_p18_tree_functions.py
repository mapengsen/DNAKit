"""Biopython tree-file and distance-matrix writers of integration pack P18 (phylogeny area).

``write_tree_file`` writes phylogenetic trees given as Newick strings to Newick, NEXUS,
phyloXML, NeXML or CDAO files (``Bio.Phylo.write``). ``format_phylip_distance_matrix`` prints a
distance matrix in the PHYLIP square format (``DistanceMatrix.format_phylip``).

Biopython is imported only when a function runs. Writers never overwrite an existing file.
"""

from __future__ import annotations

import io
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _biopython_p18_common as _c
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = _c.EXTRA
_TREE_FORMATS = ("newick", "nexus", "phyloxml", "nexml", "cdao")
_MAX_TREES = 1_000
_MAX_TAXA = 2_000


def write_tree_file(
    newick_trees: Sequence[str],
    output_path: str | Path,
    *,
    file_format: str = "nexus",
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one or more trees, given as Newick strings, to a new tree file.

    ``file_format`` is ``newick`` (one tree per line), ``nexus``, ``phyloxml``, ``nexml`` or
    ``cdao`` (RDF). Internal node labels that are numbers become branch confidences, as
    Biopython's Newick reader does, and are written back as such. The output file must not
    exist.
    """

    fmt = _c.choice(file_format, _TREE_FORMATS, "file_format")
    if isinstance(newick_trees, str) or not isinstance(newick_trees, Sequence):
        raise ConfigurationError("newick_trees must be a list of Newick strings.")
    if not 1 <= len(newick_trees) <= _MAX_TREES:
        raise ConfigurationError(f"newick_trees must hold 1 to {_MAX_TREES} trees.")
    destination = _c.new_output(output_path)
    phylo = require("Bio.Phylo", _EXTRA)
    trees = []
    for index, text in enumerate(newick_trees):
        value = _c.text(text, f"newick_trees[{index}]")
        try:
            trees.append(phylo.read(io.StringIO(value), "newick"))
        except Exception as exc:
            raise ConfigurationError(f"newick_trees[{index}] is not valid Newick: {exc}") from exc
    try:
        with destination.open("x", encoding="utf-8") as handle:
            phylo.write(trees, handle, fmt)
    except Exception as exc:
        destination.unlink(missing_ok=True)
        raise ConfigurationError(f"Biopython could not write {fmt}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_tree_file",
        "output_path": str(destination),
        "file_format": fmt,
        "tree_count": len(trees),
        "bytes_written": destination.stat().st_size,
        "provenance": _c.provenance(reader="Bio.Phylo.read (newick)", writer="Bio.Phylo.write"),
    }


def format_phylip_distance_matrix(
    names: Sequence[str],
    distances: Sequence[Sequence[float]],
    *,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format a symmetric distance matrix as PHYLIP text (``DistanceMatrix.format_phylip``).

    ``names`` lists the taxa and ``distances`` is the full square matrix (zero diagonal,
    symmetric, non-negative). The text has the number of taxa on the first line, then one line
    per taxon with the name and its distances to all taxa, to four decimals, as Biopython
    prints them. When ``output_path`` is given the text is also written to that new file.
    """

    if isinstance(names, str) or not isinstance(names, Sequence):
        raise ConfigurationError("names must be a list of taxon names.")
    if not 2 <= len(names) <= _MAX_TAXA:
        raise ConfigurationError(f"names must hold 2 to {_MAX_TAXA} taxa.")
    labels = [_c.identifier(name, "taxon name") for name in names]
    if len(set(labels)) != len(labels):
        raise ConfigurationError("taxon names must be unique.")
    size = len(labels)
    table: object = distances
    if (
        isinstance(table, str)
        or not isinstance(table, Sequence)
        or len(table) != size
        or any(
            isinstance(row, str) or not isinstance(row, Sequence) or len(row) != size
            for row in table
        )
    ):
        raise ConfigurationError(f"distances must be a {size} x {size} matrix.")
    matrix = [[_c.number(value, "distance") for value in row] for row in distances]
    for i in range(size):
        if matrix[i][i] != 0.0:
            raise ConfigurationError("the diagonal of distances must be zero.")
        for j in range(i):
            if matrix[i][j] < 0 or abs(matrix[i][j] - matrix[j][i]) > 1e-12:
                raise ConfigurationError("distances must be symmetric and non-negative.")
    destination = None if output_path is None else _c.new_output(output_path)
    tree_construction = require("Bio.Phylo.TreeConstruction", _EXTRA)
    lower = [row[: i + 1] for i, row in enumerate(matrix)]
    buffer = io.StringIO()
    tree_construction.DistanceMatrix(labels, lower).format_phylip(buffer)
    content = buffer.getvalue()
    if destination is not None:
        with destination.open("x", encoding="utf-8") as handle:
            handle.write(content)
    return {
        "status": "completed",
        "operation": "format_phylip_distance_matrix",
        "taxon_count": size,
        "text": content,
        "output_path": None if destination is None else str(destination),
        "provenance": _c.provenance(
            method="Bio.Phylo.TreeConstruction.DistanceMatrix.format_phylip"
        ),
    }


__all__ = ["format_phylip_distance_matrix", "write_tree_file"]

_register_public_verbose(globals(), tuple(__all__))
