"""deepTools wrapper for a region-splitting variant that the existing wrappers do not reach.

``chop_regions_from_middle`` wraps ``deeptools.heatmapper.chopRegionsFromMiddle``. Packages are
imported only when the function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-deeptools]"`` hint.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.deeptools_functions import (
    _DEEPTOOLS_EXTRA,
    _MAX_LENGTH,
    _integer,
    _interval_list,
    _provenance,
    _upstream_call,
)


def chop_regions_from_middle(
    exons: Sequence[Sequence[int]],
    *,
    left: int = 0,
    right: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cut a fixed number of bases on each side of the centre of a set of exons.

    Calls ``deeptools.heatmapper.chopRegionsFromMiddle``. The centre is the point after
    ``total_exon_length // 2`` exonic bases; an exon that spans it is split. The last ``left``
    bases before the centre and the first ``right`` bases after it are returned as lists of
    ``[start, end]`` intervals. When a side has fewer bases than requested, ``pad_left`` or
    ``pad_right`` holds the missing count. Unlike ``chop_regions_for_profile`` the intervals are
    cut from the centre and there is no body.
    """

    pairs = _interval_list(exons, "exons")
    left_size = _integer(left, "left", 0, _MAX_LENGTH)
    right_size = _integer(right, "right", 0, _MAX_LENGTH)
    heatmapper = require("deeptools.heatmapper", _DEEPTOOLS_EXTRA)
    left_bins, right_bins, pad_left, pad_right = _upstream_call(
        heatmapper.chopRegionsFromMiddle, pairs, left_size, right_size
    )
    return {
        "status": "completed",
        "operation": "chop_regions_from_middle",
        "left_bins": [[int(s), int(e)] for s, e in left_bins],
        "right_bins": [[int(s), int(e)] for s, e in right_bins],
        "pad_left": int(pad_left),
        "pad_right": int(pad_right),
        "provenance": _provenance(
            "deeptools",
            "flanks of the given sizes are cut on each side of the centre of the exons",
            note="pad_left and pad_right count bases missing from a side that is too short",
        ),
    }


__all__ = ["chop_regions_from_middle"]

_register_public_verbose(globals(), tuple(__all__))
