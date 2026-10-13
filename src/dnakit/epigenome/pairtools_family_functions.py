"""Pairtools wrappers for variants that the existing pairtools wrappers do not reach.

``deduplicate_pairs_online`` wraps the Cython online deduplication, which marks duplicates against
kept pairs instead of merging chains of near-identical pairs, so its counts can differ from
``deduplicate_pairs``. ``normalize_scaling_curve_unity`` wraps the area normalisation of a
scaling curve. Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-pairtools]"`` hint.
"""

from __future__ import annotations

import io
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.backends import (
    _MAX_ROWS,
    _PAIRTOOLS_EXTRA,
    _integer,
    _new_path,
    _path,
    _plain,
    _provenance,
    _read_pairs_header,
    _writer,
)
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_MAX_LISTED_IDS = 10_000
_MAX_CURVE = 10_000
_MAX_DISTANCE = 4_000_000_000
_REQUIRED = ("readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2")


def deduplicate_pairs_online(
    pairs_path: str | Path,
    *,
    output_path: str | Path | None = None,
    duplicates_path: str | Path | None = None,
    unmapped_path: str | Path | None = None,
    max_mismatch: int = 3,
    method: str = "max",
    mark_duplicates: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mark PCR duplicate pairs with pairtools' Cython online algorithm.

    Calls ``pairtools.lib.dedup.streaming_dedup_cython`` on a sorted ``.pairs`` file. A pair is a
    duplicate when both of its sides lie within ``max_mismatch`` bp of a pair that was kept
    earlier (``method="max"`` uses the larger and ``"sum"`` the summed difference of the two
    positions). Duplicates are not chained: ``deduplicate_pairs`` joins pairs that are close to a
    duplicate into one cluster, so on a run of pairs spaced 2 bp apart it keeps fewer pairs than
    this function does. ``mark_duplicates=True`` rewrites the pair type of duplicates written to
    ``duplicates_path`` as ``DD``. Optional output files must not exist yet; the duplicate and
    unique read IDs are also listed (up to 10000 each).
    """

    source = _path(pairs_path, "pairs_path")
    outputs = {
        "output_path": output_path,
        "duplicates_path": duplicates_path,
        "unmapped_path": unmapped_path,
    }
    targets = {
        name: _new_path(value, name) if value is not None else None
        for name, value in outputs.items()
    }
    named = [path for path in targets.values() if path is not None]
    if len(set(named)) != len(named):
        raise ConfigurationError("output_path, duplicates_path and unmapped_path must differ.")
    mismatch = _integer(max_mismatch, "max_mismatch", 0, 1_000)
    if method not in ("max", "sum"):
        raise ConfigurationError("method must be 'max' or 'sum'.")
    if not isinstance(mark_duplicates, bool):
        raise ConfigurationError("mark_duplicates must be a bool.")
    pairsam_format = require("pairtools.lib.pairsam_format", _PAIRTOOLS_EXTRA)
    dedup = require("pairtools.lib.dedup", _PAIRTOOLS_EXTRA)
    counter_class = require("pairtools.lib.stats", _PAIRTOOLS_EXTRA).PairCounter
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    header, columns = _read_pairs_header(source)
    missing = [name for name in _REQUIRED if name not in columns]
    if missing:
        raise ConfigurationError(f"the pairs file lacks the columns {missing}.")
    if "pair_type" not in columns:
        raise ConfigurationError("the pairs file needs a pair_type column.")
    index = {name: columns.index(name) for name in _REQUIRED}
    counter = counter_class(
        bytile_dups=False, filters=None, startup_code="", type_cast=(), engine="pandas"
    )
    with source.open("r", encoding="utf-8") as handle:
        _, body = headerops.get_header(handle)
        lines = [line if line.endswith("\n") else line + "\n" for line in body]
    if len(lines) > _MAX_ROWS:
        raise ConfigurationError(f"the pairs file has more than {_MAX_ROWS} pairs.")
    with (
        _writer(targets["output_path"]) as unique_file,
        _writer(targets["duplicates_path"]) as duplicate_file,
        _writer(targets["unmapped_path"]) as unmapped_file,
    ):
        for name, file in (
            ("output_path", unique_file),
            ("duplicates_path", duplicate_file),
            ("unmapped_path", unmapped_file),
        ):
            if targets[name] is not None:
                file.writelines(line + "\n" for line in header)
        streams = {
            "unique": io.StringIO(),
            "duplicates": io.StringIO(),
            "unmapped": io.StringIO(),
        }
        tees = {
            "unique": _Tee(streams["unique"], unique_file),
            "duplicates": _Tee(streams["duplicates"], duplicate_file),
            "unmapped": _Tee(streams["unmapped"], unmapped_file),
        }
        try:
            dedup.streaming_dedup_cython(
                method,
                mismatch,
                "\t",
                index["chrom1"],
                index["chrom2"],
                index["pos1"],
                index["pos2"],
                index["strand1"],
                index["strand2"],
                [],
                [],
                pairsam_format.UNMAPPED_CHROM,
                iter(lines),
                tees["unique"],
                tees["duplicates"],
                tees["unmapped"],
                counter,
                mark_duplicates,
                False,
                index["readID"],
            )
        except (ValueError, IndexError, KeyError) as exc:
            raise BackendExecutionError(f"pairtools could not deduplicate: {exc}") from exc

    def read_ids(name: str) -> list[str]:
        text = streams[name].getvalue()
        return [line.split("\t", 1)[0] for line in text.splitlines() if line]

    unique_ids = read_ids("unique")
    duplicate_ids = read_ids("duplicates")
    mapped = int(counter["total_mapped"])
    duplicates = int(counter["total_dups"])
    return {
        "status": "completed",
        "operation": "deduplicate_pairs_online",
        "pairs": str(source),
        "output": None if targets["output_path"] is None else str(targets["output_path"]),
        "duplicates_output": (
            None if targets["duplicates_path"] is None else str(targets["duplicates_path"])
        ),
        "unmapped_output": (
            None if targets["unmapped_path"] is None else str(targets["unmapped_path"])
        ),
        "max_mismatch": mismatch,
        "method": method,
        "mark_duplicates": mark_duplicates,
        "total": int(counter["total"]),
        "total_mapped": mapped,
        "total_nodups": int(counter["total_nodups"]),
        "total_dups": duplicates,
        "duplicate_fraction": _plain(duplicates / mapped) if mapped else None,
        "unique_read_ids": unique_ids[:_MAX_LISTED_IDS],
        "duplicate_read_ids": duplicate_ids[:_MAX_LISTED_IDS],
        "ids_truncated": len(unique_ids) > _MAX_LISTED_IDS or len(duplicate_ids) > _MAX_LISTED_IDS,
        "provenance": _provenance(
            "pairtools",
            "online duplicate detection against kept pairs (pairtools.lib.dedup."
            "streaming_dedup_cython); duplicates are not chained",
            note="counts can differ from deduplicate_pairs, which joins chains of close pairs",
        ),
    }


class _Tee:
    """Write each chunk to an in-memory buffer and to an optional file handle."""

    def __init__(self, buffer: io.StringIO, handle: Any) -> None:
        self._buffer = buffer
        self._handle = handle

    def write(self, text: str) -> int:
        self._handle.write(text)
        return self._buffer.write(text)

    def __bool__(self) -> bool:
        return True


def normalize_scaling_curve_unity(
    bin_edges: Sequence[int],
    frequencies: Sequence[float],
    *,
    norm_range: Sequence[int] | None = (10_000, 1_000_000_000),
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a contact-frequency curve so that its area under the curve is 1.

    Calls ``pairtools.lib.scaling.unity_norm_scaling``. ``bin_edges`` are the ``n + 1`` edges of
    the distance bins and ``frequencies`` the ``n`` contact frequencies. Every frequency is divided
    by ``sum(bin_length * frequency)`` over the finite bins whose geometric midpoint
    ``sqrt(low * high)`` lies strictly inside ``norm_range``; ``norm_range=None`` uses all finite
    bins. This differs from ``normalize_scaling_curve``, which divides by the mean frequency in
    a window.
    """

    edges_in = _sequence(bin_edges, "bin_edges")
    freq_in = _sequence(frequencies, "frequencies")
    if len(edges_in) != len(freq_in) + 1 or len(freq_in) < 1:
        raise ConfigurationError("bin_edges must have exactly one more value than frequencies.")
    edges = [_integer(value, "bin edge", 1, _MAX_DISTANCE) for value in edges_in]
    if edges != sorted(set(edges)):
        raise ConfigurationError("bin_edges must be unique and ascending.")
    values = []
    for value in freq_in:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigurationError("frequencies must be numbers.")
        number = float(value)
        if math.isnan(number) or number < 0:
            raise ConfigurationError("frequencies must be non-negative numbers.")
        values.append(number)
    window: tuple[int, int] | None = None
    if norm_range is not None:
        if isinstance(norm_range, (str, bytes)) or not isinstance(norm_range, Sequence):
            raise ConfigurationError("norm_range must be None or a pair (low, high).")
        if len(norm_range) != 2:
            raise ConfigurationError("norm_range must be None or a pair (low, high).")
        low = _integer(norm_range[0], "norm_range low", 1, _MAX_DISTANCE)
        high = _integer(norm_range[1], "norm_range high", low, _MAX_DISTANCE)
        window = (low, high)
    scaling = require("pairtools.lib.scaling", _PAIRTOOLS_EXTRA)
    numpy = require("numpy", _PAIRTOOLS_EXTRA)
    with numpy.errstate(divide="ignore", invalid="ignore"):
        normalized = scaling.unity_norm_scaling(
            numpy.asarray(edges, dtype=float), numpy.asarray(values, dtype=float), window
        )
    mids = numpy.sqrt(
        numpy.asarray(edges[1:], dtype=float) * numpy.asarray(edges[:-1], dtype=float)
    )
    return {
        "status": "completed",
        "operation": "normalize_scaling_curve_unity",
        "norm_range": None if window is None else list(window),
        "bins": [
            {
                "min_dist": edges[index],
                "max_dist": edges[index + 1],
                "midpoint": _plain(float(mids[index])),
                "frequency": _plain(float(normalized[index])),
            }
            for index in range(len(values))
        ],
        "provenance": _provenance(
            "pairtools",
            "frequencies divided by the area sum(bin_length * frequency) inside norm_range "
            "(pairtools.lib.scaling.unity_norm_scaling)",
            note="a window without finite bins gives non-finite values, reported as null",
        ),
    }


def _sequence(value: object, name: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if len(value) > _MAX_CURVE + 1:
        raise ConfigurationError(f"{name} may contain at most {_MAX_CURVE + 1} values.")
    return value


__all__ = ["deduplicate_pairs_online", "normalize_scaling_curve_unity"]

_register_public_verbose(globals(), tuple(__all__))
