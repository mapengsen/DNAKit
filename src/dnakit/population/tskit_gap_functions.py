"""tskit drawing and random-sequence functions for the coverage-gap rows (wave 1d).

- ``draw_tree_sequence_svg()`` writes an SVG drawing of a whole tree sequence, or of one tree, with
  tskit's ``TreeSequence.draw_svg()`` and ``Tree.draw_svg()``.
- ``draw_tree_sequence_html()``, ``draw_tree_html()`` and ``draw_variant_html()`` write tskit's HTML
  views of a tree sequence (``tskit.tree_sequence_html``), of one tree (``tskit.tree_html``) and of
  one variant (``tskit.variant_html``) to a new file.
- ``tskit_random_nucleotides()`` draws a seeded random nucleotide string with
  ``tskit.random_nucleotides()``. Its values differ from ``generate_random_dna()`` for the same seed
  because the generators differ.

Drawing inputs are the explicit table rows of ``build_tree_sequence()``. Output files are never
overwritten, no network is used, and tskit is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import islice
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.population.backends import (
    _DEFAULT_MAX_ROWS,
    _TSKIT_EXTRA,
    _integer,
    _number,
    _provenance,
    _tree_sequence,
)
from dnakit.population.tskit_functions import _new_file, _tree_index

_MAX_RANDOM_LENGTH = 10_000_000
_MAX_SEED = 2**32 - 1
_RANDOM_CONVENTION = (
    "tskit.random_nucleotides: numpy RandomState(seed) draws uniformly from the alphabet ACTG; "
    "the same seed gives the same string as tskit, not as generate_random_dna"
)
_DRAW_CONVENTION = (
    "SVG from tskit draw_svg(); the whole tree sequence when tree_index is null, otherwise the "
    "tree at that index; intervals are half-open [left, right)"
)
_HTML_CONVENTION = (
    "HTML from tskit's tree_sequence_html, tree_html or variant_html (the notebook view of the "
    "object), written as UTF-8; the file must not exist yet"
)


def tskit_random_nucleotides(
    *,
    length: int,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a seeded random nucleotide string of ``length`` bases, using tskit's generator.

    Characters are drawn uniformly from ``ACTG`` by ``tskit.random_nucleotides()``. The same
    ``seed`` reproduces tskit's string exactly; ``seed=None`` gives a fresh random string.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    count = _integer(length, "length", minimum=0, maximum=_MAX_RANDOM_LENGTH)
    chosen = None if seed is None else _integer(seed, "seed", minimum=0, maximum=_MAX_SEED)
    sequence = tskit.random_nucleotides(count, seed=chosen)
    return {
        "status": "completed",
        "operation": "tskit_random_nucleotides",
        "length": len(sequence),
        "seed": chosen,
        "alphabet": "ACTG",
        "sequence": sequence,
        "provenance": _provenance("tskit", _RANDOM_CONVENTION),
    }


def _x_lim(value: Sequence[float], length: float) -> list[float]:
    """Two numbers ``[left, right]`` with 0 <= left < right <= sequence length."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError("x_lim must be two numbers [left, right].")
    left = _number(value[0], "x_lim left")
    right = _number(value[1], "x_lim right", positive=True)
    if not left < right <= length:
        raise ConfigurationError(f"x_lim must satisfy 0 <= left < right <= {length:g}.")
    return [left, right]


def _write_text(text: str, target: Path) -> int:
    try:
        target.write_text(text, encoding="utf-8")
    except OSError as exc:
        raise BackendExecutionError(f"could not write {target.name}: {exc}") from exc
    if not target.is_file() or target.stat().st_size == 0:
        raise BackendExecutionError(f"tskit produced an empty file for {target}.")
    return target.stat().st_size


def draw_tree_sequence_svg(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    tree_index: int | None = None,
    x_lim: Sequence[float] | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a tree sequence, or one of its trees, as an SVG file, using tskit's ``draw_svg()``.

    ``tables`` takes the same explicit rows as ``build_tree_sequence()``. With ``tree_index`` null
    the whole tree sequence is drawn (``TreeSequence.draw_svg()``); with an index only that tree is
    drawn (``Tree.draw_svg()``). ``x_lim`` (two numbers ``[left, right]``, whole sequence only) sets
    the drawn interval, which must lie inside the sequence length. The SVG file is written to
    ``destination``, which must not exist yet. The result returns the image path and the tree
    intervals.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    target = _new_file(destination, "destination")
    applied_x_lim: list[float] | None = None
    if tree_index is None:
        kwargs: dict[str, Any] = {}
        if x_lim is not None:
            applied_x_lim = _x_lim(x_lim, float(ts.sequence_length))
            kwargs["x_lim"] = applied_x_lim
        ts.draw_svg(path=str(target), **kwargs)
        scope = "tree_sequence"
    else:
        if x_lim is not None:
            raise ConfigurationError("x_lim applies to the whole tree sequence, not to one tree.")
        index = _tree_index(tree_index, int(ts.num_trees))
        ts.at_index(index).draw_svg(path=str(target))
        scope = "tree"
    intervals = [[float(tree.interval.left), float(tree.interval.right)] for tree in ts.trees()]
    return {
        "status": "completed",
        "operation": "draw_tree_sequence_svg",
        "image_path": str(target),
        "format": "svg",
        "scope": scope,
        "tree_index": tree_index,
        "x_lim": applied_x_lim,
        "bytes": target.stat().st_size,
        "sequence_length": float(ts.sequence_length),
        "num_trees": int(ts.num_trees),
        "intervals": intervals,
        "num_samples": int(ts.num_samples),
        "num_sites": int(ts.num_sites),
        "num_mutations": int(ts.num_mutations),
        "provenance": _provenance("tskit", _DRAW_CONVENTION),
    }


def draw_tree_sequence_html(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write tskit's HTML summary of a tree sequence (``tskit.tree_sequence_html``) to a new file.

    ``tables`` takes the explicit rows of ``build_tree_sequence()``. The summary lists the table
    sizes, the provenance records and the metadata of the tree sequence. The result returns the HTML
    path and the counts of the tree sequence.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    target = _new_file(destination, "destination")
    size = _write_text(tskit.tree_sequence_html(ts), target)
    return {
        "status": "completed",
        "operation": "draw_tree_sequence_html",
        "image_path": str(target),
        "format": "html",
        "bytes": size,
        "sequence_length": float(ts.sequence_length),
        "num_trees": int(ts.num_trees),
        "num_samples": int(ts.num_samples),
        "num_sites": int(ts.num_sites),
        "num_mutations": int(ts.num_mutations),
        "num_provenances": int(ts.num_provenances),
        "provenance": _provenance("tskit", _HTML_CONVENTION),
    }


def draw_tree_html(
    tables: Mapping[str, Any],
    tree_index: int,
    destination: str | Path,
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write tskit's HTML view of one tree of a tree sequence (``tskit.tree_html``) to a new file.

    ``tree_index`` selects the tree in genomic order (0 to num_trees - 1). The result returns the
    HTML path, the tree interval ``[left, right)``, the number of roots and sites of that tree.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    index = _tree_index(tree_index, int(ts.num_trees))
    tree = ts.at_index(index)
    target = _new_file(destination, "destination")
    size = _write_text(tskit.tree_html(tree), target)
    return {
        "status": "completed",
        "operation": "draw_tree_html",
        "image_path": str(target),
        "format": "html",
        "bytes": size,
        "tree_index": index,
        "interval": [float(tree.interval.left), float(tree.interval.right)],
        "num_roots": int(tree.num_roots),
        "num_sites": int(tree.num_sites),
        "provenance": _provenance("tskit", _HTML_CONVENTION),
    }


def draw_variant_html(
    tables: Mapping[str, Any],
    site_index: int,
    destination: str | Path,
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write tskit's HTML view of one variant (``tskit.variant_html``) to a new file.

    ``site_index`` selects the site in position order (0 to num_sites - 1), which is the order of
    ``TreeSequence.variants()``. The result returns the HTML path, the site position, its alleles
    (null for a missing allele) and the number of samples in the genotype row.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    index = _integer(site_index, "site_index", minimum=0, maximum=int(ts.num_sites) - 1)
    variant = next(islice(ts.variants(), index, None))
    target = _new_file(destination, "destination")
    size = _write_text(tskit.variant_html(variant), target)
    return {
        "status": "completed",
        "operation": "draw_variant_html",
        "image_path": str(target),
        "format": "html",
        "bytes": size,
        "site_index": index,
        "position": float(variant.site.position),
        "alleles": [None if allele is None else str(allele) for allele in variant.alleles],
        "num_samples": int(ts.num_samples),
        "provenance": _provenance("tskit", _HTML_CONVENTION),
    }


__all__ = [
    "draw_tree_html",
    "draw_tree_sequence_html",
    "draw_tree_sequence_svg",
    "draw_variant_html",
    "tskit_random_nucleotides",
]

_register_public_verbose(globals())
