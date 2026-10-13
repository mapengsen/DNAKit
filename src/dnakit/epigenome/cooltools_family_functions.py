"""Cooltools wrappers for sandbox variants that the existing cooltools wrappers do not reach.

``pairwise_region_sums`` adds the inter-chromosomal block row to the diagonal sums, and
``compute_region_contact_scaling`` computes a distance-decay curve with the possible-pair areas of
the regions from explicit contact positions. Packages are imported only when a function runs. A
missing package raises ``BackendUnavailableError`` with the ``pip install "dnakit[bio-cooltools]"``
hint.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.cooltools_review_functions import (
    _COOLTOOLS_EXTRA,
    _MAX_POSITION,
    _NUMPY_EXTRA,
    _existing,
    _integer,
    _open,
    _plain,
    _provenance,
    _records,
    _weight,
)
from dnakit.exceptions import ConfigurationError

_MAX_CONTACTS = 5_000_000
_MAX_DISTANCE_BINS = 1_000


def pairwise_region_sums(
    cool_path: str | Path,
    *,
    weight: str = "weight",
    ignore_diags: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Diagonal sums within each region pair plus one total row for each pair of chromosomes.

    Calls ``cooltools.sandbox.obs_over_exp_cooler.sum_pairwise`` on the chromosome view of the
    cooler. Rows with ``region1 == region2`` and a distance of 0 or more equal the rows of
    ``pairwise_diagonal_sums`` (``count.sum`` is empty for the ignored diagonals). In addition,
    every pair of different regions has one row with ``dist = -1`` that holds the total count of
    the block (``count.sum``) and the number of valid pixels in it (``n_valid``); ``n_total`` is
    empty for those rows.
    """

    source = _existing(cool_path, "cool_path")
    diags = _integer(ignore_diags, "ignore_diags", 0, 10)
    _, clr = _open(source)
    name = _weight(clr, weight)
    common = require("cooltools.lib.common", _COOLTOOLS_EXTRA)
    sandbox = require("cooltools.sandbox.obs_over_exp_cooler", _COOLTOOLS_EXTRA)
    view = common.make_cooler_view(clr)
    try:
        frame = sandbox.sum_pairwise(clr, view, clr_weight_name=name, ignore_diags=diags)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not sum the regions: {exc}") from exc
    table = frame.sort_values(["region1", "region2", "dist"]).reset_index(drop=True)
    table = table[["region1", "region2", "dist", "n_total", "n_valid", "count.sum"]]
    return {
        "status": "completed",
        "operation": "pairwise_region_sums",
        "cool_path": str(source),
        "weight": name,
        "ignore_diags": diags,
        "row_count": len(table),
        "block_row_count": int((table["dist"] == -1).sum()),
        "rows": _records(table),
        "provenance": _provenance(
            "cooltools",
            "diagonal sums of intra-chromosomal blocks and total sums of inter-chromosomal "
            "blocks (sandbox.obs_over_exp_cooler.sum_pairwise)",
            order="sorted by region1, region2 and distance; dist -1 marks a whole block",
        ),
    }


def _region(value: object, name: str) -> tuple[int, int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 2:
        raise ConfigurationError(f"{name} must be a pair [start, end].")
    start = _integer(value[0], f"{name} start", 0, _MAX_POSITION)
    end = _integer(value[1], f"{name} end", 1, _MAX_POSITION)
    if end <= start:
        raise ConfigurationError(f"{name} must have end greater than start.")
    return start, end


def compute_region_contact_scaling(
    positions1: Sequence[int],
    positions2: Sequence[int],
    *,
    region1: Sequence[int],
    region2: Sequence[int] | None = None,
    min_distance: int = 10,
    max_distance: int = 10_000_000,
    n_bins: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Raw contact counts and possible-pair areas per log-spaced distance bin for one region pair.

    Calls ``cooltools.sandbox.pairs_scaling_functions.compute_scaling``. Contact ``i`` joins
    ``positions1[i]`` and ``positions2[i]`` on one chromosome. Only contacts with
    ``positions1`` inside ``region1`` and ``positions2`` inside ``region2`` (both given as
    ``[start, end)``; ``region2`` defaults to ``region1``) are counted, in the distance bin of
    ``positions2 - positions1``, for distances from ``min_distance`` up to but not including
    ``max_distance``. The area of a bin is the integral of the number of possible position pairs
    over the bin, so ``frequency = observed / area`` is a contact density. The counts are raw:
    unlike ``compute_contact_scaling`` there are no chromosome sizes and no pair-file input.
    """

    first = _region(region1, "region1")
    second = first if region2 is None else _region(region2, "region2")
    low = _integer(min_distance, "min_distance", 1, _MAX_POSITION)
    high = _integer(max_distance, "max_distance", low + 1, _MAX_POSITION)
    count = _integer(n_bins, "n_bins", 2, _MAX_DISTANCE_BINS)
    for name, values in (("positions1", positions1), ("positions2", positions2)):
        if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
            raise ConfigurationError(f"{name} must be a list of integer positions.")
    if len(positions1) != len(positions2) or not positions1:
        raise ConfigurationError("positions1 and positions2 must have equal, non-zero length.")
    if len(positions1) > _MAX_CONTACTS:
        raise ConfigurationError(f"At most {_MAX_CONTACTS} contacts can be given.")
    first_positions = [_integer(v, "positions1 entries", 0, _MAX_POSITION) for v in positions1]
    second_positions = [_integer(v, "positions2 entries", 0, _MAX_POSITION) for v in positions2]
    pandas = require("pandas", _NUMPY_EXTRA)
    scaling = require("cooltools.sandbox.pairs_scaling_functions", _COOLTOOLS_EXTRA)
    frame = pandas.DataFrame({"pos1": first_positions, "pos2": second_positions})
    try:
        edges, observed, areas = scaling.compute_scaling(
            frame,
            first,
            None if region2 is None else second,
            dmin=low,
            dmax=high,
            n_bins=count,
        )
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"cooltools could not compute the scaling: {exc}") from exc
    bins = []
    for index in range(len(observed)):
        area = float(areas[index])
        bins.append(
            {
                "min_dist": int(edges[index]),
                "max_dist": int(edges[index + 1]),
                "observed": int(observed[index]),
                "area": _plain(area),
                "frequency": None if area == 0 else _plain(float(observed[index]) / area),
            }
        )
    return {
        "status": "completed",
        "operation": "compute_region_contact_scaling",
        "region1": list(first),
        "region2": list(second),
        "dist_range": [low, high],
        "n_bins": count,
        "contacts_counted": int(sum(item["observed"] for item in bins)),
        "bins": bins,
        "provenance": _provenance(
            "cooltools",
            "histogram of contact distances over log-spaced bins divided by the possible-pair "
            "area of the regions (sandbox.pairs_scaling_functions.compute_scaling)",
            convention="bin edges are rounded to integers; the area is a continuous integral",
        ),
    }


__all__ = ["compute_region_contact_scaling", "pairwise_region_sums"]

_register_public_verbose(globals(), tuple(__all__))
