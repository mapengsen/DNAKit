"""Gimmemotifs wrappers for variants that the existing gimmemotifs wrappers do not reach.

``normalize_peak_widths`` wraps ``gimmemotifs.preprocessing.peakfile2bedfile``. Unlike
``equal_size_intervals`` it is strand aware, drops intervals that would start below zero and always
writes six columns. ``count_motif_matches_at_fpr`` wraps ``gimmemotifs.rocmetrics.matches_at_fpr``,
which returns two counts instead of one number and so cannot go through ``motif_score_metric``.
"""

from __future__ import annotations

import math
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import _fraction, _gimme_provenance, _integer

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_MAX_ROWS = 1_000_000


def normalize_peak_widths(
    path: str | Path,
    *,
    window: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resize peaks to a fixed window around their midpoint with ``peakfile2bedfile``.

    The midpoint of each peak is ``(start + end) // 2``. Rows on the ``+`` strand (the default when
    the strand column is missing) become ``midpoint - window // 2`` to ``midpoint + window // 2``.
    For other strands gimmemotifs writes the two coordinates the other way round, so ``start`` is
    larger than ``end``; such rows carry ``start_after_end=True``. Rows whose start would be below
    zero are dropped, not clamped. Missing name, score and strand columns become ``""``, ``0`` and
    ``+``. Lines starting with ``#`` or ``track`` are skipped.
    """

    peaks = require("gimmemotifs.preprocessing", _GIMMEMOTIFS_EXTRA)
    source = input_path(path)
    width = _integer(window, "window", 1, 1_000_000)
    cleaned: list[str] = []
    with source.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ConfigurationError(f"{source.name} line {number} needs at least 3 columns.")
            try:
                int(fields[1])
                int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{source.name} line {number} has non-integer coordinates."
                ) from exc
            cleaned.append(line.strip())
            if len(cleaned) > _MAX_ROWS:
                raise ConfigurationError(f"{source.name} has more than {_MAX_ROWS} rows.")
    rows = len(cleaned)
    if rows == 0:
        raise ConfigurationError(f"{source.name} contains no peaks.")
    with tempfile.TemporaryDirectory(prefix="dnakit-gimmemotifs-") as scratch:
        staged = Path(scratch) / "input.bed"
        staged.write_text("\n".join(cleaned) + "\n", encoding="utf-8")
        output = Path(scratch) / "peaks.bed"
        try:
            peaks.peakfile2bedfile(str(staged), str(output), window=width)
        except (KeyboardInterrupt, SystemExit):
            raise
        except BaseException as exc:  # gimmemotifs raises several exception types
            raise BackendExecutionError(f"gimmemotifs failed: {exc}") from exc
        written = output.read_text(encoding="utf-8").splitlines()
    intervals = []
    for text in written:
        chrom, start, end, name, score, strand = text.split("\t")
        intervals.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "name": name,
                "score": score,
                "strand": strand,
                "start_after_end": int(start) > int(end),
            }
        )
    return {
        "status": "completed",
        "operation": "normalize_peak_widths",
        "window": width,
        "input_rows": rows,
        "interval_count": len(intervals),
        "dropped_rows": rows - len(intervals),
        "intervals": intervals,
        "provenance": _gimme_provenance(
            function="gimmemotifs.preprocessing.peakfile2bedfile",
            coordinate_system="0-based half-open BED coordinates as written by gimmemotifs",
            convention="strand aware; start > end on non-plus strands; start < 0 rows dropped",
        ),
    }


def count_motif_matches_at_fpr(
    foreground: Sequence[float],
    background: Sequence[float],
    *,
    fpr: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count foreground and background sequences scoring above the score cutoff at a given FPR.

    Calls ``gimmemotifs.rocmetrics.matches_at_fpr``: the cutoff is the ``100 - fpr * 100``
    percentile of the background scores, and the result is the number of foreground scores and the
    number of background scores at or above it. ``motif_score_metric`` cannot return these two
    counts, because it expects one number per metric.
    """

    rocmetrics = require("gimmemotifs.rocmetrics", _GIMMEMOTIFS_EXTRA)
    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    rate = _fraction(fpr, "fpr")
    if rate <= 0.0:
        raise ConfigurationError("fpr must be greater than 0.")
    arrays = [
        np.asarray(_finite_scores(values, name), dtype=float)
        for name, values in (("foreground", foreground), ("background", background))
    ]
    try:
        counts = rocmetrics.matches_at_fpr(arrays[0], arrays[1], fpr=rate)
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:  # gimmemotifs raises several exception types
        raise BackendExecutionError(f"gimmemotifs failed: {exc}") from exc
    foreground_matches, background_matches = (int(value) for value in counts)
    return {
        "status": "completed",
        "operation": "count_motif_matches_at_fpr",
        "fpr": rate,
        "foreground_count": int(arrays[0].size),
        "background_count": int(arrays[1].size),
        "foreground_matches": foreground_matches,
        "background_matches": background_matches,
        "foreground_fraction": foreground_matches / int(arrays[0].size),
        "background_fraction": background_matches / int(arrays[1].size),
        "provenance": _gimme_provenance(
            function="gimmemotifs.rocmetrics.matches_at_fpr",
            convention="cutoff is the (100 - fpr * 100) percentile of the background scores; "
            "scores equal to the cutoff count as matches",
        ),
    }


def calculate_partial_motif_roc_auc(
    foreground: Sequence[float],
    background: Sequence[float],
    *,
    fpr_limit: float = 0.1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Area under the motif-score ROC curve up to a false-positive-rate limit.

    Calls ``gimmemotifs.rocmetrics.roc_auc_xlim`` with the foreground and background score lists
    (the upstream parameters are named ``x_bla`` and ``y_bla``, but they hold the scores, not
    curve points). The area is the integral of the true-positive rate over false-positive rates
    from 0 to ``fpr_limit``, so it is at most ``fpr_limit``. ``motif_score_metric`` cannot return
    it: that wrapper feeds the function the ROC curve points, which makes gimmemotifs fail.
    """

    rocmetrics = require("gimmemotifs.rocmetrics", _GIMMEMOTIFS_EXTRA)
    limit = _fraction(fpr_limit, "fpr_limit")
    if limit <= 0.0:
        raise ConfigurationError("fpr_limit must be greater than 0.")
    lists = []
    for name, values in (("foreground", foreground), ("background", background)):
        lists.append(_finite_scores(values, name))
    try:
        area = float(rocmetrics.roc_auc_xlim(lists[0][:], lists[1][:], xlim=limit))
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:  # gimmemotifs raises several exception types
        raise BackendExecutionError(f"gimmemotifs failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "calculate_partial_motif_roc_auc",
        "fpr_limit": limit,
        "foreground_count": len(lists[0]),
        "background_count": len(lists[1]),
        "area": area,
        "area_fraction_of_limit": area / limit,
        "provenance": _gimme_provenance(
            function="gimmemotifs.rocmetrics.roc_auc_xlim",
            convention="integral of the true-positive rate over false-positive rates 0 to "
            "fpr_limit; divide by fpr_limit for a value between 0 and 1",
        ),
    }


def _finite_scores(values: object, name: str) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of numbers.")
    if len(values) > _MAX_ROWS:
        raise ConfigurationError(f"{name} has more than {_MAX_ROWS} scores.")
    out = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigurationError(f"{name} must contain only numbers.")
        if not math.isfinite(value):
            raise ConfigurationError(f"{name} must contain only finite numbers.")
        out.append(float(value))
    return out


__all__ = ["calculate_partial_motif_roc_auc", "count_motif_matches_at_fpr", "normalize_peak_widths"]

_register_public_verbose(globals(), tuple(__all__))
