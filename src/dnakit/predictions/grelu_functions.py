"""Offline gReLU utilities for DNA encodings, sequence statistics, intervals, motifs and variants.

Each public function imports ``grelu`` (and ``torch`` where tensors are needed) only when
it runs. Inputs are local: sequences, nested lists of encoded sequences, interval and
variant mappings, PWM matrices and FASTA files. Model training, model zoo downloads,
genome-based lookups through genomepy, and plotting are not wrapped here.
"""

from __future__ import annotations

import contextlib
import io
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _clean,
    _dna,
    _dna_batch,
    _integer,
    _provenance,
    _rejected,
    _text,
)

_GRELU = "grelu"
_EXTRAS = "bio-grelu"
_ENCODINGS = ("strings", "indices", "one_hot")
_CONVERSIONS = {
    ("indices", "strings"),
    ("one_hot", "strings"),
    ("one_hot", "indices"),
    ("indices", "one_hot"),
}
_BLACKLIST_GENOMES = ("hg19", "hg38", "mm10")
_GC_MODES = ("distribution", "fraction")
_OVERLAP_METHODS = ("any", "all")
_MAX_ROWS = 100_000
_MAX_SEQUENCE = 1_000_000
_MAX_MOTIFS = 200
_MAX_MOTIF_LENGTH = 1_000
_MAX_SEQUENCES = 10_000
_MAX_COORDINATE = 10**12
_MAX_CHROM = 128
_MAX_FASTA_SEQUENCES = 10_000
_ALLELE_CHARS = frozenset("ACGTN-")


def _quiet(function: Callable[..., Any], *args: Any, verbose: bool, **kwargs: Any) -> Any:
    """Run an upstream call that prints progress; keep its messages when verbose is True."""

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        result = function(*args, **kwargs)
    messages = [line.strip() for line in buffer.getvalue().splitlines() if line.strip()]
    return result, (messages if verbose else [])


def _rows(values: object, name: str, *, limit: int = _MAX_ROWS) -> list[Any]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} may contain at most {limit} entries.")
    return list(values)


def _interval_frame(values: object, name: str, *, extra: Sequence[str] = ()) -> Any:
    """Build the chrom/start/end table that gReLU's interval helpers expect."""

    pandas = require("pandas", _EXTRAS)
    records = []
    for index, row in enumerate(_rows(values, name)):
        if not isinstance(row, Mapping):
            raise ConfigurationError(f"{name}[{index}] must be a mapping with chrom, start, end.")
        label = f"{name}[{index}]"
        chrom = _text(row.get("chrom"), f"{label}.chrom", maximum=_MAX_CHROM)
        start = _integer(row.get("start"), f"{label}.start", minimum=0, maximum=_MAX_COORDINATE)
        end = _integer(row.get("end"), f"{label}.end", minimum=1, maximum=_MAX_COORDINATE)
        if end <= start:
            raise ConfigurationError(f"{label}.end must be greater than start.")
        record: dict[str, Any] = {"chrom": chrom, "start": start, "end": end}
        for key in extra:
            if key not in row:
                raise ConfigurationError(f"{label} is missing the column {key!r}.")
            record[key] = _clean(row[key])
        records.append(record)
    return pandas.DataFrame(records, columns=["chrom", "start", "end", *extra])


def _records(frame: Any, extra: Sequence[str] = ()) -> list[dict[str, Any]]:
    """Sorted, JSON-safe interval records with the extra columns kept as given."""

    output = []
    for row in frame.to_dict(orient="records"):
        record: dict[str, Any] = {
            "chrom": str(row["chrom"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
        }
        for key in extra:
            record[key] = _clean(row.get(key))
        output.append(record)
    output.sort(
        key=lambda item: (
            item["chrom"],
            item["start"],
            item["end"],
            repr([item.get(key) for key in extra]),
        )
    )
    return output


def _matrix(values: object, name: str, *, rows: int = 4) -> Any:
    """Return a finite float array of shape (rows, L) from nested lists."""

    np = require("numpy", _EXTRAS)
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or len(values) != rows:
        raise ConfigurationError(f"{name} must be a list of {rows} rows.")
    lengths = set()
    for row in values:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise ConfigurationError(f"{name} must be a list of {rows} numeric rows.")
        lengths.add(len(row))
    if len(lengths) != 1 or 0 in lengths or max(lengths) > _MAX_MOTIF_LENGTH:
        raise ConfigurationError(f"{name} rows must have one length between 1 and 1000.")
    for row in values:
        for item in row:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise ConfigurationError(f"{name} entries must be finite numbers.")
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ConfigurationError(f"{name} entries must be finite numbers.")
    return array


def _motif_set(values: object, name: str) -> dict[str, Any]:
    if not isinstance(values, Mapping) or not values:
        raise ConfigurationError(f"{name} must be a non-empty mapping of names to 4xL matrices.")
    if len(values) > _MAX_MOTIFS:
        raise ConfigurationError(f"{name} may contain at most {_MAX_MOTIFS} motifs.")
    return {
        _text(key, f"{name} key", maximum=_MAX_CHROM): _matrix(value, f"{name}[{key!r}]")
        for key, value in sorted(values.items(), key=lambda item: str(item[0]))
    }


def _allele(value: object, name: str) -> str:
    text = _text(value, name, maximum=_MAX_SEQUENCE).upper()
    if set(text) - _ALLELE_CHARS:
        raise ConfigurationError(f"{name} may contain only A, C, G, T, N and '-'.")
    return text


def _json_matrix(array: Any) -> list[list[float]]:
    return [[float(item) + 0.0 for item in row] for row in array.tolist()]


def convert_dna_encoding(
    values: Sequence[Any],
    *,
    source: str,
    target: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a batch of encoded DNA between gReLU's indices, one-hot and string forms.

    Indices are integers 0 to 3 for A, C, G, T and 4 for N. One-hot input and output
    are lists of 4 x L matrices with channels A, C, G, T; an N column is all zeros.
    All sequences in a batch need the same length. String input is handled by
    ``one_hot_encode_dna`` and ``integer_encode`` in dnakit.
    """

    src = _choice(source, "source", ("indices", "one_hot"))
    dst = _choice(target, "target", ("strings", "indices", "one_hot"))
    if (src, dst) not in _CONVERSIONS:
        raise ConfigurationError(f"Unsupported conversion {src} -> {dst}.")
    grelu_format = require("grelu.sequence.format", _EXTRAS)
    torch = require("torch", _EXTRAS)
    np = require("numpy", _EXTRAS)
    batch = _rows(values, "values", limit=_MAX_SEQUENCES)
    if src == "indices":
        indices = _index_batch(batch)
        try:
            if dst == "strings":
                output: Any = grelu_format.indices_to_strings(indices)
            else:
                output = grelu_format.indices_to_one_hot(indices).numpy()
        except (ValueError, TypeError) as exc:
            raise _rejected(_GRELU, exc) from exc
    else:
        one_hot = _one_hot_batch(batch, torch, np)
        try:
            if dst == "strings":
                output = grelu_format.one_hot_to_strings(one_hot)
            else:
                output = grelu_format.one_hot_to_indices(one_hot)
        except (ValueError, TypeError) as exc:
            raise _rejected(_GRELU, exc) from exc
    if dst == "strings":
        result_values: Any = [str(item) for item in output]
    elif dst == "indices":
        result_values = [[int(item) for item in row] for row in np.asarray(output).tolist()]
    else:
        result_values = [_json_matrix(np.asarray(item)) for item in output]
    return {
        "status": "completed",
        "operation": "convert_dna_encoding",
        "source": src,
        "target": dst,
        "count": len(result_values),
        "values": result_values,
        "provenance": _provenance(
            _GRELU,
            index_convention="A=0, C=1, G=2, T=3, N=4",
            one_hot_convention="channels A, C, G, T; N is all zeros",
        ),
    }


def _index_rows(batch: list[Any]) -> list[list[int]]:
    """Validate container shape and integer types; range checks are left to callers."""

    rows = []
    for index, row in enumerate(batch):
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row:
            raise ConfigurationError(f"values[{index}] must be a non-empty list of indices.")
        for item in row:
            if isinstance(item, bool) or not isinstance(item, int):
                raise ConfigurationError(f"values[{index}] must contain integers only.")
        rows.append(list(row))
    if len({len(row) for row in rows}) != 1:
        raise ConfigurationError("All indices in a batch must have the same length.")
    return rows


def _index_batch(batch: list[Any]) -> Any:
    """Indices for conversion: gReLU expects int8 arrays with values 0 (A) to 4 (N)."""

    np = require("numpy", _EXTRAS)
    rows = _index_rows(batch)
    for index, row in enumerate(rows):
        if any(item < 0 or item > 4 for item in row):
            raise ConfigurationError(f"values[{index}] indices must be integers from 0 to 4.")
    return np.asarray(rows, dtype=np.int8)


def _one_hot_columns_ok(matrices: Any, np: Any) -> bool:
    """gReLU's one-hot check only tests dtype and shape; each column must also be 0/1
    with at most one 1 (all zeros is N)."""

    binary = bool(np.isin(matrices, (0.0, 1.0)).all())
    return binary and bool((matrices.sum(axis=1) <= 1).all())


def _one_hot_batch(batch: list[Any], torch: Any, np: Any) -> Any:
    grelu_format = require("grelu.sequence.format", _EXTRAS)
    matrices = [_matrix(item, f"values[{index}]") for index, item in enumerate(batch)]
    if len({matrix.shape for matrix in matrices}) != 1 or matrices[0].shape[0] != 4:
        raise ConfigurationError("All one-hot matrices in a batch must be 4 x L with equal L.")
    stacked = np.stack(matrices)
    if not _one_hot_columns_ok(stacked, np):
        raise ConfigurationError("One-hot input must have one 1 per column, or all zeros for N.")
    tensor = torch.tensor(stacked, dtype=torch.float32)
    if not grelu_format.check_one_hot(tensor):
        raise ConfigurationError("One-hot input must be a float tensor shaped (B, 4, L).")
    return tensor


def check_dna_encoding(
    values: Sequence[Any],
    *,
    encoding: str = "strings",
    verbose: bool = True,
) -> dict[str, Any]:
    """Check whether a batch is valid DNA strings, indices or one-hot matrices.

    ``valid`` combines gReLU's dtype, shape and value checks with the per-column rule for
    one-hot input. Strings are uppercased first. Ragged lists and non-integer or
    non-numeric entries raise ``ConfigurationError``; out-of-range values return False.
    """

    kind = _choice(encoding, "encoding", _ENCODINGS)
    grelu_format = require("grelu.sequence.format", _EXTRAS)
    batch = _rows(values, "values", limit=_MAX_SEQUENCES)
    if kind == "strings":
        if any(not isinstance(item, str) for item in batch):
            raise ConfigurationError("String encoding needs a list of DNA strings.")
        valid = bool(grelu_format.check_string_dna([item.upper() for item in batch]))
    elif kind == "indices":
        np = require("numpy", _EXTRAS)
        rows = _index_rows(batch)
        # int8 is gReLU's index dtype; values outside it cannot be valid indices.
        if any(not -128 <= item <= 127 for row in rows for item in row):
            valid = False
        else:
            valid = bool(grelu_format.check_indices(np.asarray(rows, dtype=np.int8)))
    else:
        torch = require("torch", _EXTRAS)
        np = require("numpy", _EXTRAS)
        matrices = [_matrix(item, f"values[{index}]") for index, item in enumerate(batch)]
        if len({matrix.shape for matrix in matrices}) != 1:
            raise ConfigurationError("All one-hot matrices must have the same shape.")
        stacked = np.stack(matrices)
        valid = bool(
            grelu_format.check_one_hot(torch.tensor(stacked, dtype=torch.float32))
        ) and _one_hot_columns_ok(stacked, np)
    return {
        "status": "completed",
        "operation": "check_dna_encoding",
        "encoding": kind,
        "count": len(batch),
        "valid": valid,
        "provenance": _provenance(_GRELU, check="grelu.sequence.format"),
    }


def gc_content_distribution(
    sequences: str | Sequence[str],
    *,
    mode: str = "distribution",
    binwidth: float = 0.1,
    normalize: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """GC fraction per sequence, or a histogram of GC fraction across sequences.

    ``mode="fraction"`` returns one GC value per sequence (N counts as non-GC and
    non-AT in gReLU's denominator convention). ``mode="distribution"`` bins the GC
    values into bins of width ``binwidth`` starting at 0, as gReLU does, and optionally
    normalizes the counts to sum to 1.
    """

    kind = _choice(mode, "mode", _GC_MODES)
    values = _dna_batch(sequences, "sequences")
    metrics = require("grelu.sequence.metrics", _EXTRAS)
    np = require("numpy", _EXTRAS)
    if isinstance(binwidth, bool) or not isinstance(binwidth, (int, float)):
        raise ConfigurationError("binwidth must be a number in (0, 1].")
    if not 0 < float(binwidth) <= 1:
        raise ConfigurationError("binwidth must be in (0, 1].")
    if not isinstance(normalize, bool):
        raise ConfigurationError("normalize must be boolean.")
    if kind == "fraction":
        try:
            fractions = metrics.gc(values)
        except (ValueError, TypeError) as exc:
            raise _rejected(_GRELU, exc) from exc
        return {
            "status": "completed",
            "operation": "gc_content_distribution",
            "mode": kind,
            "count": len(values),
            "gc": [_clean(float(item)) for item in np.atleast_1d(fractions)],
            "provenance": _provenance(_GRELU, denominator="all bases including N"),
        }
    try:
        counts = metrics.gc_distribution(values, binwidth=float(binwidth), normalize=normalize)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    edges = [float(item) for item in np.arange(0, 1, float(binwidth))]
    return {
        "status": "completed",
        "operation": "gc_content_distribution",
        "mode": kind,
        "count": len(values),
        "binwidth": float(binwidth),
        "bin_edges": edges,
        "counts": [_clean(float(item)) for item in counts],
        "provenance": _provenance(_GRELU, binning="np.arange(0, 1, binwidth); lower-edge bins"),
    }


def summarize_sequence_lengths(
    sequences: str | Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lengths of each sequence, whether all are equal, and the shared length if so."""

    values = _dna_batch(sequences, "sequences")
    utils = require("grelu.sequence.utils", _EXTRAS)
    try:
        lengths = [int(item) for item in utils.get_lengths(values)]
        equal = bool(utils.check_equal_lengths(values))
        unique = int(utils.get_unique_length(values)) if equal else None
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "summarize_sequence_lengths",
        "count": len(lengths),
        "lengths": lengths,
        "all_equal": equal,
        "unique_length": unique,
        "provenance": _provenance(_GRELU, length_unit="bases"),
    }


def dna_sequence_differences(
    sequence_a: str,
    sequence_b: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """List positions where two equal-length DNA sequences differ, with both bases."""

    left = _dna(sequence_a, "sequence_a")
    right = _dna(sequence_b, "sequence_b")
    if len(left) != len(right):
        raise ConfigurationError("sequence_a and sequence_b must have the same length.")
    mutate = require("grelu.sequence.mutate", _EXTRAS)
    try:
        positions = [int(item) for item in mutate.seq_differences(left, right, verbose=False)]
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "dna_sequence_differences",
        "length": len(left),
        "count": len(positions),
        "differences": [
            {"position": item, "reference": left[item], "alternate": right[item]}
            for item in positions
        ],
        "provenance": _provenance(_GRELU, position_convention="0-based"),
    }


def extend_intervals_from_centers(
    intervals: Sequence[Mapping[str, Any]],
    *,
    seq_len: int,
    center_col: str = "summit",
    verbose: bool = True,
) -> dict[str, Any]:
    """Make fixed-length intervals around a summit given as an offset from each start.

    This follows gReLU's convention: ``start + summit - seq_len // 2`` (the narrowPeak
    summit column is an offset from the peak start). Extended starts below 0 are rejected.
    """

    length = _integer(seq_len, "seq_len", minimum=1, maximum=_MAX_SEQUENCE)
    key = _text(center_col, "center_col", maximum=_MAX_CHROM)
    frame = _interval_frame(intervals, "intervals", extra=(key,))
    for index, value in enumerate(frame[key].tolist()):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            raise ConfigurationError(f"intervals[{index}].{key} must be a non-negative integer.")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        extended = preprocess.extend_from_coord(frame, seq_len=length, center_col=key)
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    if (extended["start"] < 0).any():
        raise ConfigurationError("An extended interval starts before coordinate 0.")
    return {
        "status": "completed",
        "operation": "extend_intervals_from_centers",
        "seq_len": length,
        "count": len(extended),
        "intervals": _records(extended),
        "provenance": _provenance(_GRELU, center_convention=f"{key} is an offset from start"),
    }


def clip_intervals_to_bounds(
    intervals: Sequence[Mapping[str, Any]],
    *,
    start: int | None = None,
    end: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip each interval's start up to ``start`` and its end down to ``end``.

    Clipping can leave an interval with start >= end; those rows are kept and counted
    in ``invalid_count`` so the caller can drop them.
    """

    lower = None if start is None else _integer(start, "start", minimum=0, maximum=_MAX_COORDINATE)
    upper = None if end is None else _integer(end, "end", minimum=0, maximum=_MAX_COORDINATE)
    frame = _interval_frame(intervals, "intervals")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        clipped = preprocess.clip_intervals(frame, start=lower, end=upper)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(clipped)
    invalid = sum(1 for item in records if item["start"] >= item["end"])
    return {
        "status": "completed",
        "operation": "clip_intervals_to_bounds",
        "count": len(records),
        "invalid_count": invalid,
        "intervals": records,
        "provenance": _provenance(_GRELU, bounds="start raised to start, end lowered to end"),
    }


def merge_intervals_by_group(
    intervals: Sequence[Mapping[str, Any]],
    *,
    group_col: str,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge intervals that share a group value into one span per value.

    Each group must map to a single chromosome. The span runs from the smallest start to
    the largest end in the group.
    """

    key = _text(group_col, "group_col", maximum=_MAX_CHROM)
    frame = _interval_frame(intervals, "intervals", extra=(key,))
    chrom_per_group = frame.groupby(key)["chrom"].nunique()
    if (chrom_per_group > 1).any():
        raise ConfigurationError("Each group value must map to a single chromosome.")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        merged = preprocess.merge_intervals_by_column(frame, group_col=key)
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "merge_intervals_by_group",
        "group_col": key,
        "input_count": len(frame),
        "count": len(merged),
        "intervals": _records(merged, extra=(key,)),
        "provenance": _provenance(_GRELU, merge="min start and max end per group"),
    }


def filter_intervals_overlapping_reference(
    intervals: Sequence[Mapping[str, Any]],
    reference: Sequence[Mapping[str, Any]],
    *,
    window: int = 0,
    invert: bool = False,
    method: str = "any",
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep intervals that overlap reference intervals (or drop them with ``invert``).

    ``method="any"`` keeps any overlap; ``method="all"`` requires the interval to be
    covered by the reference. ``window`` widens each reference interval by that many
    bases on both sides.
    """

    extension = _integer(window, "window", minimum=0, maximum=_MAX_COORDINATE)
    choice = _choice(method, "method", _OVERLAP_METHODS)
    if not isinstance(invert, bool):
        raise ConfigurationError("invert must be boolean.")
    frame = _interval_frame(intervals, "intervals")
    ref = _interval_frame(reference, "reference")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        kept, messages = _quiet(
            preprocess.filter_overlapping,
            frame,
            ref,
            window=extension,
            invert=invert,
            inplace=False,
            method=choice,
            verbose=verbose,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(kept)
    return {
        "status": "completed",
        "operation": "filter_intervals_overlapping_reference",
        "method": choice,
        "invert": invert,
        "window": extension,
        "input_count": len(frame),
        "count": len(records),
        "intervals": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, overlap="half-open 0-based coordinates"),
    }


def split_intervals_by_chromosome(
    intervals: Sequence[Mapping[str, Any]],
    *,
    train_chroms: Sequence[str] | None = None,
    val_chroms: Sequence[str] | None = ("chr10",),
    test_chroms: Sequence[str] | None = ("chr11",),
    verbose: bool = True,
) -> dict[str, Any]:
    """Split intervals into train, validation and test sets by chromosome.

    Defaults follow gReLU: validation is chr10 and test is chr11. ``train_chroms=None``
    gives training every chromosome not used for validation or test. The three lists must
    not share a chromosome.
    """

    frame = _interval_frame(intervals, "intervals")
    lists = {
        "train_chroms": _chrom_list(train_chroms, "train_chroms", allow_none=True),
        "val_chroms": _chrom_list(val_chroms, "val_chroms", allow_none=True),
        "test_chroms": _chrom_list(test_chroms, "test_chroms", allow_none=True),
    }
    named = [set(value) for value in lists.values() if value is not None]
    for index, first in enumerate(named):
        for second in named[index + 1 :]:
            if first & second:
                raise ConfigurationError("train, validation and test chromosomes must be disjoint.")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        (train, val, test), messages = _quiet(
            preprocess.split,
            frame,
            train_chroms=lists["train_chroms"],
            val_chroms=lists["val_chroms"],
            test_chroms=lists["test_chroms"],
            sample=[],
            seed=None,
            verbose=verbose,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "split_intervals_by_chromosome",
        "counts": {"train": len(train), "validation": len(val), "test": len(test)},
        "train": _records(train),
        "validation": _records(val),
        "test": _records(test),
        "messages": messages,
        "provenance": _provenance(_GRELU, split="by chromosome"),
    }


def _chrom_list(values: object, name: str, *, allow_none: bool) -> list[str] | None:
    if values is None and allow_none:
        return None
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of chromosome names.")
    return [
        _text(item, f"{name}[{index}]", maximum=_MAX_CHROM) for index, item in enumerate(values)
    ]


def filter_intervals_blacklist(
    intervals: Sequence[Mapping[str, Any]],
    *,
    genome: str = "hg38",
    window: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove intervals that overlap the ENCODE blacklist bundled with gReLU.

    Bundled blacklists exist for ``hg19``, ``hg38`` and ``mm10``; no download is needed.
    ``window`` pads each blacklist region by that many bases.
    """

    name = _choice(genome, "genome", _BLACKLIST_GENOMES)
    extension = _integer(window, "window", minimum=0, maximum=_MAX_COORDINATE)
    frame = _interval_frame(intervals, "intervals")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        kept, messages = _quiet(
            preprocess.filter_blacklist,
            frame,
            genome=name,
            window=extension,
            inplace=False,
            verbose=verbose,
        )
    except (ValueError, TypeError, KeyError, OSError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(kept)
    return {
        "status": "completed",
        "operation": "filter_intervals_blacklist",
        "genome": name,
        "window": extension,
        "input_count": len(frame),
        "count": len(records),
        "removed_count": len(frame) - len(records),
        "intervals": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, blacklist=f"ENCODE v2 {name} (bundled)"),
    }


def filter_intervals_by_chromosome(
    intervals: Sequence[Mapping[str, Any]],
    *,
    include: Sequence[str] | None = None,
    exclude: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep intervals on the ``include`` chromosomes and drop those on ``exclude``."""

    if include is None and exclude is None:
        raise ConfigurationError("Pass include, exclude, or both.")
    included = _chrom_list(include, "include", allow_none=True)
    excluded = _chrom_list(exclude, "exclude", allow_none=True)
    frame = _interval_frame(intervals, "intervals")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        kept, messages = _quiet(
            preprocess.filter_chromosomes,
            frame,
            include=included,
            exclude=excluded,
            inplace=False,
            verbose=verbose,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(kept)
    return {
        "status": "completed",
        "operation": "filter_intervals_by_chromosome",
        "input_count": len(frame),
        "count": len(records),
        "intervals": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, filter="chromosome names as given"),
    }


def pwm_to_consensus(
    pwm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Consensus DNA string of a 4 x L position weight matrix (rows A, C, G, T).

    Each position takes its highest-weight base; ties go to the earlier base in A, C, G, T.
    """

    matrix = _matrix(pwm, "pwm")
    motifs = require("grelu.interpret.motifs", _EXTRAS)
    try:
        consensus = motifs.motifs_to_strings(matrix)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "pwm_to_consensus",
        "length": int(matrix.shape[1]),
        "consensus": str(consensus),
        "provenance": _provenance(_GRELU, rule="argmax per position"),
    }


def trim_motif_pwm(
    pwm: Sequence[Sequence[float]],
    *,
    trim_threshold: float = 0.3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim the flanks of a 4 x L PWM where per-position weight is low.

    The position score is the sum of absolute values per column. A column is kept when its
    score is at least ``trim_threshold`` times the largest score. Use log-odds or centred
    PWMs: for probability columns every score is 1 and nothing is trimmed.
    """

    matrix = _matrix(pwm, "pwm")
    threshold = _real(trim_threshold, "trim_threshold")
    motifs = require("grelu.interpret.motifs", _EXTRAS)
    try:
        start, end = motifs.trim_pwm(matrix, trim_threshold=threshold, return_indices=True)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    start_index, end_index = int(start), int(end)
    trimmed = matrix[:, start_index:end_index]
    return {
        "status": "completed",
        "operation": "trim_motif_pwm",
        "input_length": int(matrix.shape[1]),
        "start": start_index,
        "end": end_index,
        "trim_threshold": threshold,
        "pwm": _json_matrix(trimmed),
        "provenance": _provenance(_GRELU, positions="0-based half-open"),
    }


def _real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [0, 1].")
    if not 0 <= float(value) <= 1:
        raise ConfigurationError(f"{name} must be in [0, 1].")
    return float(value)


def scan_motifs_with_pwms(
    sequences: Mapping[str, str] | Sequence[str],
    motifs: Mapping[str, Sequence[Sequence[float]]],
    *,
    pthresh: float = 0.001,
    rc: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan DNA sequences with several PWMs, FIMO-style, on both strands by default.

    ``sequences`` is an ID-to-sequence mapping or a list (IDs ``seq0``, ``seq1``, ...).
    Each motif is a 4 x L PWM of probabilities or log-odds with rows A, C, G, T. Hits are
    sorted by sequence, start, strand and motif. Only hits with p-value at or below
    ``pthresh`` are returned.
    """

    ids, seqs = _sequence_items(sequences)
    pwms = _motif_set(motifs, "motifs")
    threshold = _real_positive(pthresh, "pthresh")
    if not isinstance(rc, bool):
        raise ConfigurationError("rc must be boolean.")
    motif_module = require("grelu.interpret.motifs", _EXTRAS)
    try:
        frame = motif_module.scan_sequences(
            seqs,
            pwms,
            names=list(pwms),
            seq_ids=ids,
            pthresh=threshold,
            rc=rc,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    hits = []
    for row in frame.to_dict(orient="records"):
        hits.append(
            {
                "motif": str(row["motif"]),
                "sequence": str(row["sequence"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": str(row["strand"]),
                "score": _clean(float(row["fimo_score"])),
                "p_value": _clean(float(row["fimo_p-value"])),
                "matched_sequence": str(row["matched_seq"]),
            }
        )
    hits.sort(key=lambda item: (item["sequence"], item["start"], item["strand"], item["motif"]))
    return {
        "status": "completed",
        "operation": "scan_motifs_with_pwms",
        "sequence_count": len(ids),
        "motif_count": len(pwms),
        "count": len(hits),
        "hits": hits,
        "provenance": _provenance(_GRELU, method="FIMO-style scan", pthresh=str(threshold)),
    }


def _sequence_items(values: object) -> tuple[list[str], list[str]]:
    if isinstance(values, Mapping):
        if not values:
            raise ConfigurationError("sequences must not be empty.")
        ids = [_text(key, "sequence id", maximum=_MAX_CHROM) for key in values]
        seqs = [_dna(value, f"sequences[{key!r}]") for key, value in values.items()]
    else:
        items = _rows(values, "sequences", limit=_MAX_SEQUENCES)
        ids = [f"seq{index}" for index in range(len(items))]
        seqs = [_dna(item, f"sequences[{index}]") for index, item in enumerate(items)]
    if len(set(ids)) != len(ids):
        raise ConfigurationError("sequence IDs must be unique.")
    return ids, seqs


def _real_positive(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in (0, 1].")
    if not 0 < float(value) <= 1:
        raise ConfigurationError(f"{name} must be in (0, 1].")
    return float(value)


def compare_motif_pwm_alleles(
    reference: str,
    alternate: str,
    motifs: Mapping[str, Sequence[Sequence[float]]],
    *,
    pthresh: float = 0.001,
    rc: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare PWM hits between a reference and an alternate allele sequence.

    Sites are reported once per motif, start and strand, with the reference and alternate
    p-values and scores and their difference (alternate minus reference). Missing values
    are ``None``.
    """

    ref_seq = _dna(reference, "reference")
    alt_seq = _dna(alternate, "alternate")
    pwms = _motif_set(motifs, "motifs")
    threshold = _real_positive(pthresh, "pthresh")
    if not isinstance(rc, bool):
        raise ConfigurationError("rc must be boolean.")
    motif_module = require("grelu.interpret.motifs", _EXTRAS)
    try:
        frame = motif_module.compare_motifs(
            ref_seq,
            alt_seq,
            pwms,
            names=list(pwms),
            pthresh=threshold,
            rc=rc,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    sites = []
    for row in frame.to_dict(orient="records"):
        sites.append(
            {
                "motif": str(row["motif"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": str(row["strand"]),
                "p_value_reference": _clean(float(row["fimo_p-value_ref"])),
                "p_value_alternate": _clean(float(row["fimo_p-value_alt"])),
                "score_reference": _clean(float(row["fimo_score_ref"])),
                "score_alternate": _clean(float(row["fimo_score_alt"])),
                "score_difference": _clean(float(row["fimo_score_diff"])),
            }
        )
    sites.sort(key=lambda item: (item["motif"], item["start"], item["strand"]))
    return {
        "status": "completed",
        "operation": "compare_motif_pwm_alleles",
        "count": len(sites),
        "sites": sites,
        "provenance": _provenance(_GRELU, comparison="alternate minus reference score"),
    }


def filter_variants_by_length(
    variants: Sequence[Mapping[str, Any]],
    *,
    standard_bases: bool = True,
    max_insert_len: int = 0,
    max_del_len: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter variants by allele type and indel length, as gReLU does for VCF-like rows.

    Each variant needs ``chrom``, ``pos`` (1-based), ``ref`` and ``alt``. With
    ``standard_bases`` only A, C, G, T, and N alleles are kept. Insertions longer than
    ``max_insert_len`` and deletions longer than ``max_del_len`` are removed.
    """

    if not isinstance(standard_bases, bool):
        raise ConfigurationError("standard_bases must be boolean.")
    insert_limit = _integer(max_insert_len, "max_insert_len", minimum=0, maximum=_MAX_SEQUENCE)
    delete_limit = _integer(max_del_len, "max_del_len", minimum=0, maximum=_MAX_SEQUENCE)
    pandas = require("pandas", _EXTRAS)
    rows = []
    for index, item in enumerate(_rows(variants, "variants")):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"variants[{index}] must be a mapping.")
        rows.append(
            {
                "chrom": _text(item.get("chrom"), f"variants[{index}].chrom", maximum=_MAX_CHROM),
                "pos": _integer(
                    item.get("pos"), f"variants[{index}].pos", minimum=1, maximum=_MAX_COORDINATE
                ),
                "ref": _allele(item.get("ref"), f"variants[{index}].ref"),
                "alt": _allele(item.get("alt"), f"variants[{index}].alt"),
            }
        )
    frame = pandas.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"])
    filtering = require("grelu.variant", _EXTRAS)
    try:
        kept, messages = _quiet(
            filtering.filter_variants,
            frame,
            standard_bases=standard_bases,
            max_insert_len=insert_limit,
            max_del_len=delete_limit,
            inplace=False,
            verbose=verbose,
        )
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = [
        {
            "chrom": str(row["chrom"]),
            "pos": int(row["pos"]),
            "ref": str(row["ref"]),
            "alt": str(row["alt"]),
        }
        for row in kept.to_dict(orient="records")
    ]
    records.sort(key=lambda item: (item["chrom"], item["pos"], item["ref"], item["alt"]))
    return {
        "status": "completed",
        "operation": "filter_variants_by_length",
        "input_count": len(rows),
        "count": len(records),
        "removed_count": len(rows) - len(records),
        "variants": records,
        "messages": messages,
        "provenance": _provenance(_GRELU, position_convention="1-based, as in VCF"),
    }


def variants_to_sequence_intervals(
    variants: Sequence[Mapping[str, Any]],
    *,
    seq_len: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fixed-length intervals centred on each variant position, as gReLU defines them.

    The rule is ``start = pos - ceil(seq_len / 2)`` with ``pos`` used as given, and
    ``end = start + seq_len``. With ``seq_len=1`` the interval is ``[pos - 1, pos)``, the
    0-based base of a 1-based VCF position.
    """

    length = _integer(seq_len, "seq_len", minimum=1, maximum=_MAX_SEQUENCE)
    pandas = require("pandas", _EXTRAS)
    rows = []
    for index, item in enumerate(_rows(variants, "variants")):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"variants[{index}] must be a mapping.")
        rows.append(
            {
                "chrom": _text(item.get("chrom"), f"variants[{index}].chrom", maximum=_MAX_CHROM),
                "pos": _integer(
                    item.get("pos"), f"variants[{index}].pos", minimum=1, maximum=_MAX_COORDINATE
                ),
                "ref": _allele(item.get("ref"), f"variants[{index}].ref"),
                "alt": _allele(item.get("alt"), f"variants[{index}].alt"),
            }
        )
    frame = pandas.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"])
    filtering = require("grelu.variant", _EXTRAS)
    try:
        intervals = filtering.variants_to_intervals(frame, seq_len=length, inplace=False)
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(intervals)
    if any(item["start"] < 0 for item in records):
        raise ConfigurationError("A variant interval starts before coordinate 0.")
    return {
        "status": "completed",
        "operation": "variants_to_sequence_intervals",
        "seq_len": length,
        "count": len(records),
        "intervals": records,
        "provenance": _provenance(_GRELU, rule="start = pos - ceil(seq_len / 2)"),
    }


def read_fasta_sequences(
    fasta_path: str | Path,
    *,
    max_sequences: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read sequences from a local FASTA (plain or gzipped) with gReLU's reader.

    gReLU returns the sequences without their headers, so the result has the sequences
    and their lengths in file order. Files with more than ``max_sequences`` records are
    rejected.
    """

    path = input_path(fasta_path)
    limit = _integer(max_sequences, "max_sequences", minimum=1, maximum=_MAX_FASTA_SEQUENCES)
    fasta = require("grelu.io.fasta", _EXTRAS)
    try:
        sequences = fasta.read_fasta(str(path))
    except (ValueError, OSError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    if len(sequences) > limit:
        raise ConfigurationError(f"FASTA has {len(sequences)} records; max_sequences is {limit}.")
    values = [str(item) for item in sequences]
    return {
        "status": "completed",
        "operation": "read_fasta_sequences",
        "source": str(path),
        "count": len(values),
        "lengths": [len(item) for item in values],
        "sequences": values,
        "provenance": _provenance(_GRELU, headers="not returned by gReLU"),
    }


def sample_intervals(
    intervals: Sequence[Mapping[str, Any]],
    *,
    n: int,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Random subset of ``n`` intervals without replacement, reproducible with ``seed``.

    Output is sorted by coordinates. gReLU's ``filter_random`` returns nothing unless it
    works in place, so the table is filtered in place on a copy made here.
    """

    count = _integer(n, "n", minimum=1, maximum=_MAX_ROWS)
    if seed is not None:
        seed = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    frame = _interval_frame(intervals, "intervals").copy()
    if count > len(frame):
        raise ConfigurationError(f"n={count} exceeds the {len(frame)} intervals provided.")
    preprocess = require("grelu.data.preprocess", _EXTRAS)
    try:
        _quiet(preprocess.filter_random, frame, count, seed=seed, inplace=True, verbose=verbose)
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    records = _records(frame)
    return {
        "status": "completed",
        "operation": "sample_intervals",
        "input_count": len(_rows(intervals, "intervals")),
        "seed": seed,
        "count": len(records),
        "intervals": records,
        "provenance": _provenance(_GRELU, sampling="without replacement"),
    }


__all__ = [
    "check_dna_encoding",
    "clip_intervals_to_bounds",
    "compare_motif_pwm_alleles",
    "convert_dna_encoding",
    "dna_sequence_differences",
    "extend_intervals_from_centers",
    "filter_intervals_blacklist",
    "filter_intervals_by_chromosome",
    "filter_intervals_overlapping_reference",
    "filter_variants_by_length",
    "gc_content_distribution",
    "merge_intervals_by_group",
    "pwm_to_consensus",
    "read_fasta_sequences",
    "sample_intervals",
    "scan_motifs_with_pwms",
    "split_intervals_by_chromosome",
    "summarize_sequence_lengths",
    "trim_motif_pwm",
    "variants_to_sequence_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
