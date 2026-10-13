"""Row-level overlap selection through bioframe (family-variant check, group H4).

``bioframe.select`` returns the overlapping intervals; the existing ``select_bioframe_region``
wrapper returns those records. The sibling functions ``select_mask``, ``select_indices`` and
``select_labels`` return positions in the input table instead, which the records do not carry,
so this module exposes them. bioframe is imported only when a function runs; a missing package
raises ``BackendUnavailableError`` with the install hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import _DEFAULT_ROWS, _package, _rows_limit, _scan_bed, _text
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_REGION_MAX = 200


def select_bioframe_overlap_rows(
    source: str | Path,
    region: str,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return which rows of a BED file overlap a region, as a boolean mask and row positions.

    ``mask`` has one boolean per data row of the file (comment lines skipped, file order);
    ``indices`` are the positions of the True entries. ``labels`` are bioframe's index labels
    for the selected rows, which equal the positions because the file is read with a default
    integer index.
    """

    path = input_path(source)
    text = _text(region, "region", maximum=_REGION_MAX)
    limit = _rows_limit(max_rows)
    bioframe = require("bioframe", "bio-bioframe")
    schema = f"bed{_scan_bed(path, limit)}"
    try:
        frame = bioframe.read_table(str(path), schema=schema, comment="#")
        mask = bioframe.select_mask(frame, text)
        indices = bioframe.select_indices(frame, text)
        labels = bioframe.select_labels(frame, text)
    except (IndexError, KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"bioframe rejected the input: {exc}") from exc
    return {
        "status": "completed",
        "operation": "select_bioframe_overlap_rows",
        "region": text,
        "row_count": len(frame),
        "selected_count": len(indices),
        "mask": [bool(value) for value in mask],
        "indices": [int(value) for value in indices],
        "labels": [int(value) for value in labels],
        "provenance": {
            "package": _package("bioframe"),
            "convention": "0-based half-open intervals; an interval is selected when it overlaps "
            "the region; a bare chromosome name selects every interval on it",
        },
    }


__all__ = ["select_bioframe_overlap_rows"]

_register_public_verbose(globals(), tuple(__all__))
