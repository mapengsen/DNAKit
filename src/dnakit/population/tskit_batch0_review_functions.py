"""tskit text-tree drawing for the batch-0 coverage-gap review.

``draw_tree_text()`` writes one tree of a tree sequence as tskit's text drawing to a new UTF-8
file. The text drawing is ``Tree.draw(format="unicode")`` or ``Tree.draw(format="ascii")``; the
same ``format`` argument of tskit's ``draw_tree`` dispatcher also selects SVG and HTML, which
DNAKit already covers with ``draw_tree_sequence_svg()`` and ``draw_tree_html()``.

Inputs are the explicit table rows accepted by the other tskit drawing functions. Output files are
never overwritten, no network is used, and tskit is imported only when the function runs.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DEFAULT_MAX_ROWS, _TSKIT_EXTRA, _provenance, _tree_sequence
from dnakit.population.tskit_functions import _new_file, _tree_index
from dnakit.population.tskit_gap_functions import _write_text

_TEXT_FORMATS = ("unicode", "ascii")
_TEXT_CONVENTION = (
    "text drawing from tskit Tree.draw(format=...): one line per level of the tree, drawn with "
    "unicode box characters or with ascii characters; written as UTF-8; the file must not exist yet"
)


def draw_tree_text(
    tables: Mapping[str, Any],
    tree_index: int,
    destination: str | Path,
    *,
    text_format: str = "unicode",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one tree of a tree sequence as tskit's text drawing (``Tree.draw``) to a new file.

    ``tree_index`` selects the tree in genomic order (0 to num_trees - 1). ``text_format`` is
    ``"unicode"`` (box characters, tskit's default text drawing) or ``"ascii"``. The result returns
    the file path, the drawing text, the tree's Newick string and its interval ``[left, right)``.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if text_format not in _TEXT_FORMATS:
        raise ConfigurationError(f"text_format must be one of {list(_TEXT_FORMATS)}.")
    ts = _tree_sequence(tskit, tables, max_rows)
    index = _tree_index(tree_index, int(ts.num_trees))
    tree = ts.at_index(index)
    drawing = str(tree.draw(format=text_format))
    target = _new_file(destination, "destination")
    size = _write_text(drawing, target)
    return {
        "status": "completed",
        "operation": "draw_tree_text",
        "image_path": str(target),
        "format": text_format,
        "bytes": size,
        "tree_index": index,
        "interval": [float(tree.interval.left), float(tree.interval.right)],
        "newick": tree.newick(),
        "num_roots": int(tree.num_roots),
        "num_sites": int(tree.num_sites),
        "text": drawing,
        "provenance": _provenance("tskit", _TEXT_CONVENTION),
    }


__all__ = ["draw_tree_text"]

_register_public_verbose(globals())
