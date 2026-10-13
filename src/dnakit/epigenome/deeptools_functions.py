"""Optional deepTools wrappers for GC content, region parsing, scale factors and bigWig tracks.

The functions call deepTools' Python helpers directly and avoid its command-line entry points.
deepTools calls ``sys.exit`` on some bad inputs; those exits are converted into
``ConfigurationError``. Packages are imported only when a function runs.
"""

from __future__ import annotations

import contextlib
import io
import itertools
import math
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_DEEPTOOLS_EXTRA = "bio-deeptools"
_TWOBIT_EXTRA = "bio-py2bit"
_BIGWIG_EXTRA = "bio-pybigwig"
_MAX_TEXT = 1_000_000
_MAX_CHROMOSOMES = 10_000
_MAX_LENGTH = 2_000_000_000
_MAX_FILES = 32
_MAX_WINDOW = 10_000_000
_MAX_SAMPLES = 1_000_000
_MAX_ROWS = 1_000_000
_MAX_COLUMNS = 1_000
_MAX_INTERVALS = 10_000_000
_AVERAGE_METHODS = ("mean", "median")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str, limit: int = _MAX_TEXT) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > limit
        or any(char in value for char in "\x00\r\n\t")
    ):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _sequence(value: object, name: str, maximum: int) -> Sequence[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if not 1 <= len(value) <= maximum:
        raise ConfigurationError(f"{name} must contain 1 to {maximum} entries.")
    return value


def _path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new_path(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _chromsizes(value: object, name: str = "chrom_sizes") -> dict[str, int]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CHROMOSOMES:
        raise ConfigurationError(
            f"{name} must map 1 to {_MAX_CHROMOSOMES} chromosome names to lengths."
        )
    sizes: dict[str, int] = {}
    for chrom, length in value.items():
        sizes[_text(chrom, "chromosome name", limit=255)] = _integer(
            length, f"length of {chrom}", 1, _MAX_LENGTH
        )
    return sizes


def _plain(value: Any) -> Any:
    """Convert numpy values to Python values and non-finite floats to None."""

    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [_plain(item) for item in sorted(value, key=repr)]
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _provenance(distribution: str, statistic: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": statistic, **notes}


def _upstream_call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run a deepTools function that prints progress and may call sys.exit on bad input."""

    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            return function(*args, **kwargs)
    except SystemExit as exc:
        lines = captured.getvalue().strip().splitlines()
        detail = lines[-1] if lines else str(exc.code)
        raise ConfigurationError(f"deepTools rejected the input: {detail}") from exc


def gc_fraction_from_twobit(
    twobit_path: str | Path,
    chrom: str,
    start: int,
    end: int,
    *,
    fraction: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """GC fraction (or G+C base count) of one interval of a 2bit genome."""

    source = _path(twobit_path, "twobit_path")
    name = _text(chrom, "chrom", limit=255)
    left = _integer(start, "start", 0, _MAX_LENGTH)
    right = _integer(end, "end", 1, _MAX_LENGTH)
    if right <= left or right - left > _MAX_WINDOW:
        raise ConfigurationError(f"end must follow start by at most {_MAX_WINDOW} bases.")
    as_fraction = _flag(fraction, "fraction")
    py2bit = require("py2bit", _TWOBIT_EXTRA)
    utilities = require("deeptools.utilities", _DEEPTOOLS_EXTRA)
    handle = py2bit.open(str(source))
    try:
        lengths = handle.chroms()
        if name not in lengths:
            raise ConfigurationError(f"{name} is not a sequence in {source.name}.")
        if right > lengths[name]:
            raise ConfigurationError(f"end {right} is beyond the {lengths[name]} bp of {name}.")
        try:
            value = utilities.getGC_content(handle, name, left, right, as_fraction)
        except Exception as exc:
            raise ConfigurationError(str(exc)) from exc
    finally:
        handle.close()
    return {
        "status": "completed",
        "operation": "gc_fraction_from_twobit",
        "twobit": str(source),
        "chrom": name,
        "start": left,
        "end": right,
        "unit": "fraction" if as_fraction else "bases",
        "gc": _plain(value),
        "provenance": _provenance(
            "deeptools",
            "(G + C) over the interval; intervals with more than 5% N bases are rejected",
            coordinate_system="0-based half-open",
        ),
    }


def parse_deeptools_region(
    chrom_sizes: Mapping[str, int],
    region: str,
    *,
    max_chunk_size: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a ``chrom:start:end[:tile]`` region against chromosome lengths, as deepTools does."""

    sizes = _chromsizes(chrom_sizes)
    text = _text(region, "region", limit=512)
    chunk_limit = _integer(max_chunk_size, "max_chunk_size", 1, 1_000_000_000)
    chrom = text.split(":", 1)[0]
    if chrom not in sizes:
        raise ConfigurationError(f"Unknown chromosome in region: {chrom}")
    mapreduce = require("deeptools.mapReduce", _DEEPTOOLS_EXTRA)
    chunks, start, stop, chunk = _upstream_call(
        mapreduce.getUserRegion, sizes, text, float(chunk_limit)
    )
    return {
        "status": "completed",
        "operation": "parse_deeptools_region",
        "region": text,
        "regions": [{"chrom": name, "length": int(length)} for name, length in chunks],
        "start": int(start),
        "end": int(stop),
        "chunk_size": int(chunk),
        "provenance": _provenance(
            "deeptools",
            "start is moved to a multiple of the tile size, end is clipped to the chromosome",
            coordinate_system="0-based, end exclusive",
        ),
    }


def estimate_bam_scale_factors(
    bam_paths: Sequence[str | Path],
    *,
    bin_length: int = 500,
    number_of_samples: int = 1000,
    normalization_length: int = 1000,
    avg_method: str = "median",
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale factors for two indexed BAM files from samples across the genome (SES method)."""

    if isinstance(bam_paths, (str, Path)) or not isinstance(bam_paths, Sequence):
        raise ConfigurationError("bam_paths must be a list of two BAM files.")
    if len(bam_paths) != 2:
        raise ConfigurationError("estimate_bam_scale_factors compares exactly two BAM files.")
    bams = [_path(item, "bam_paths") for item in bam_paths]
    if len(set(bams)) != 2:
        raise ConfigurationError("bam_paths must name two different files.")
    for bam in bams:
        if not (bam.with_name(bam.name + ".bai").is_file() or bam.with_suffix(".bai").is_file()):
            raise ConfigurationError(f"{bam.name} needs a .bai index next to it.")
    window = _integer(bin_length, "bin_length", 1, _MAX_WINDOW)
    samples = _integer(number_of_samples, "number_of_samples", 1, _MAX_SAMPLES)
    normalization = _integer(normalization_length, "normalization_length", 1, _MAX_WINDOW)
    if avg_method not in _AVERAGE_METHODS:
        raise ConfigurationError("avg_method must be 'mean' or 'median'.")
    estimator = require("deeptools.estimateScaleFactor", _DEEPTOOLS_EXTRA)
    result = _upstream_call(
        estimator.estimateScaleFactor,
        [str(bam) for bam in bams],
        window,
        samples,
        normalization,
        avg_method=avg_method,
        numberOfProcessors=1,
        verbose=False,
    )
    values = _plain(result)
    return {
        "status": "completed",
        "operation": "estimate_bam_scale_factors",
        "bams": [str(bam) for bam in bams],
        "bin_length": window,
        "number_of_samples": samples,
        "normalization_length": normalization,
        "avg_method": avg_method,
        "size_factors": values["size_factors"],
        "sites_sampled": values["sites_sampled"],
        "estimates": {
            key: values[key]
            for key in sorted(values)
            if key not in ("size_factors", "sites_sampled")
        },
        "provenance": _provenance(
            "deeptools",
            "genome sampled at evenly spaced sites; size factors from the sampled read counts",
            note="requires exactly two BAM files",
        ),
    }


def bedgraph_to_bigwig(
    chrom_sizes: Mapping[str, int],
    bedgraph_paths: Sequence[str | Path],
    bigwig_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write sorted bedGraph files, in chromosome order, to one bigWig file."""

    sizes = _chromsizes(chrom_sizes)
    files = _sequence(bedgraph_paths, "bedgraph_paths", _MAX_FILES)
    sources = [_path(item, "bedgraph_paths") for item in files]
    target = _new_path(bigwig_path, "bigwig_path")
    order = {name: index for index, name in enumerate(sizes)}
    last_chrom = -1
    last_end = 0
    intervals = 0
    for source in sources:
        with source.open("r", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                fields = line.split()
                if not fields:
                    continue
                if len(fields) < 4 or fields[0] not in order:
                    raise ConfigurationError(
                        f"{source.name} line {number} is not a bedGraph record of a known "
                        "chromosome."
                    )
                chrom_index = order[fields[0]]
                start, end = int(fields[1]), int(fields[2])
                if chrom_index < last_chrom:
                    raise ConfigurationError("bedGraph files must follow the order of chrom_sizes.")
                if chrom_index != last_chrom:
                    last_chrom = chrom_index
                    last_end = 0
                if start < last_end or end <= start or end > sizes[fields[0]]:
                    raise ConfigurationError(
                        f"{source.name} line {number} is unsorted, overlapping or out of range."
                    )
                last_end = end
                intervals += 1
                if intervals > _MAX_INTERVALS:
                    raise ConfigurationError(
                        f"More than {_MAX_INTERVALS} intervals are not allowed."
                    )
    writer = require("deeptools.writeBedGraph", _DEEPTOOLS_EXTRA)
    require("pyBigWig", _BIGWIG_EXTRA)
    _upstream_call(
        writer.bedGraphToBigWig,
        list(sizes.items()),
        [str(source) for source in sources],
        str(target),
    )
    return {
        "status": "completed",
        "operation": "bedgraph_to_bigwig",
        "bigwig": str(target),
        "chromosomes": list(sizes),
        "files": len(sources),
        "intervals": intervals,
        "provenance": _provenance(
            "deeptools",
            "bedGraph intervals written with pyBigWig; input must be sorted and non-overlapping",
        ),
    }


def _interval_list(value: object, name: str) -> list[tuple[int, int]]:
    items = _sequence(value, name, 100_000)
    pairs: list[tuple[int, int]] = []
    for item in items:
        if isinstance(item, (str, bytes, Mapping)) or not isinstance(item, Sequence):
            raise ConfigurationError(f"{name} must contain (start, end) pairs.")
        if len(item) != 2:
            raise ConfigurationError(f"{name} must contain (start, end) pairs.")
        start = _integer(item[0], f"{name} start", 0, _MAX_LENGTH)
        end = _integer(item[1], f"{name} end", 1, _MAX_LENGTH)
        if end <= start:
            raise ConfigurationError(f"{name} must have end greater than start.")
        pairs.append((start, end))
    if any(later[0] < earlier[1] for earlier, later in itertools.pairwise(pairs)):
        raise ConfigurationError(f"{name} must be sorted and non-overlapping.")
    return pairs


def chop_regions_for_profile(
    exons: Sequence[Sequence[int]],
    *,
    left: int = 0,
    right: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split exon intervals into unscaled 5' and 3' flanks and a scaled body, as for profiles."""

    pairs = _interval_list(exons, "exons")
    left_size = _integer(left, "left", 0, _MAX_LENGTH)
    right_size = _integer(right, "right", 0, _MAX_LENGTH)
    heatmapper = require("deeptools.heatmapper", _DEEPTOOLS_EXTRA)
    left_bins, body_bins, right_bins, pad_left, pad_right = _upstream_call(
        heatmapper.chopRegions, pairs, left_size, right_size
    )
    return {
        "status": "completed",
        "operation": "chop_regions_for_profile",
        "left_bins": [[int(s), int(e)] for s, e in left_bins],
        "body_bins": [[int(s), int(e)] for s, e in body_bins],
        "right_bins": [[int(s), int(e)] for s, e in right_bins],
        "pad_left": int(pad_left),
        "pad_right": int(pad_right),
        "provenance": _provenance(
            "deeptools",
            "unscaled flanks of the given sizes are cut from the exon ends; the rest is the body",
            note="pad_left and pad_right count bases missing from a flank that is too short",
        ),
    }


def match_chromosome_name(
    chrom: str, chromosome_names: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Find the chromosome name that matches a request with or without the chr prefix."""

    name = _text(chrom, "chrom", limit=255)
    names = [
        _text(item, "chromosome_names entry", limit=255)
        for item in _sequence(chromosome_names, "chromosome_names", _MAX_CHROMOSOMES)
    ]
    utilities = require("deeptools.utilities", _DEEPTOOLS_EXTRA)
    matched = _upstream_call(utilities.mungeChromosome, name, names)
    return {
        "status": "completed",
        "operation": "match_chromosome_name",
        "requested": name,
        "matched": matched,
        "provenance": _provenance(
            "deeptools",
            "adds or removes the chr prefix until the name is in the list",
            note="matched is None when no variant is listed",
        ),
    }


def map_twobit_to_bam_chromosomes(
    twobit_names: Sequence[str], bam_names: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Map 2bit chromosome names to BAM reference names, adding or removing the chr prefix."""

    twobit = [
        _text(item, "twobit_names entry", limit=255)
        for item in _sequence(twobit_names, "twobit_names", _MAX_CHROMOSOMES)
    ]
    bam = [
        _text(item, "bam_names entry", limit=255)
        for item in _sequence(bam_names, "bam_names", _MAX_CHROMOSOMES)
    ]
    utilities = require("deeptools.utilities", _DEEPTOOLS_EXTRA)
    mapping = _upstream_call(utilities.tbitToBamChrName, twobit, bam)
    if mapping is None:
        raise ConfigurationError("The 2bit and BAM chromosome names cannot be matched.")
    pairs = {str(key): str(value) for key, value in mapping.items()}
    return {
        "status": "completed",
        "operation": "map_twobit_to_bam_chromosomes",
        "mapping": {key: pairs[key] for key in sorted(pairs)},
        "complete": set(pairs) == set(twobit),
        "provenance": _provenance(
            "deeptools",
            "names are equal, or differ only by a chr prefix",
        ),
    }


def fingerprint_jensen_shannon(
    chip_counts: Sequence[float], input_counts: Sequence[float], *, verbose: bool = True
) -> dict[str, Any]:
    """Jensen-Shannon distance between ChIP and input cumulative coverage curves."""

    chip = _sequence(chip_counts, "chip_counts", _MAX_ROWS)
    control = _sequence(input_counts, "input_counts", _MAX_ROWS)
    if len(chip) != len(control) or len(chip) < 2:
        raise ConfigurationError("chip_counts and input_counts need equal length of at least 2.")
    chip_values = [_number(item, "chip_counts entry", 0.0, 1e12) for item in chip]
    input_values = [_number(item, "input_counts entry", 0.0, 1e12) for item in control]
    if sum(chip_values) <= 0 or sum(input_values) <= 0:
        raise ConfigurationError("Both count vectors must contain at least one count.")
    plot = require("deeptools.plotFingerprint", _DEEPTOOLS_EXTRA)
    value = _upstream_call(plot.getJSDcommon, chip_values, input_values)
    return {
        "status": "completed",
        "operation": "fingerprint_jensen_shannon",
        "bins": len(chip_values),
        "jsd": _plain(value),
        "provenance": _provenance(
            "deeptools",
            "Jensen-Shannon measure between the cumulative coverage curves of ChIP and input",
            note="the same value plotFingerprint reports; bins are used in the given order",
        ),
    }


def estimate_count_size_factors(
    counts: Sequence[Sequence[float]], *, verbose: bool = True
) -> dict[str, Any]:
    """DESeq2-style size factors for a bin-by-sample count matrix, inverted for bamCoverage."""

    rows = _sequence(counts, "counts", _MAX_ROWS)
    matrix: list[list[float]] = []
    width: int | None = None
    for row in rows:
        if isinstance(row, (str, bytes, Mapping)) or not isinstance(row, Sequence):
            raise ConfigurationError("counts must be a list of equal-length rows.")
        if width is None:
            width = len(row)
            if not 2 <= width <= _MAX_COLUMNS:
                raise ConfigurationError(f"counts must have 2 to {_MAX_COLUMNS} samples.")
        if len(row) != width:
            raise ConfigurationError("counts must be a list of equal-length rows.")
        matrix.append([_number(item, "counts entry", 0.0, 1e12) for item in row])
    numpy = require("numpy", _DEEPTOOLS_EXTRA)
    counter = require("deeptools.countReadsPerBin", _DEEPTOOLS_EXTRA)
    factors = _upstream_call(counter.estimateSizeFactors, numpy.array(matrix, dtype=float))
    values = _plain(factors)
    if any(item is None for item in values):
        raise ConfigurationError("No row without zero counts; size factors cannot be estimated.")
    return {
        "status": "completed",
        "operation": "estimate_count_size_factors",
        "rows": len(matrix),
        "samples": width,
        "size_factors": values,
        "provenance": _provenance(
            "deeptools",
            "DESeq2 size factors on rows without zeros, returned as their inverse",
            compatibility="values can be used directly as bamCoverage scale factors",
        ),
    }


def bigwig_chromosome_sizes(
    bigwig_paths: Sequence[str | Path], *, verbose: bool = True
) -> dict[str, Any]:
    """Chromosome names and lengths shared by bigWig files, and any that differ between them."""

    files = _sequence(bigwig_paths, "bigwig_paths", _MAX_FILES)
    sources = [_path(item, "bigwig_paths") for item in files]
    bigwig = require("deeptools.getScorePerBigWigBin", _DEEPTOOLS_EXTRA)
    common, non_common = _upstream_call(bigwig.getChromSizes, [str(item) for item in sources])
    return {
        "status": "completed",
        "operation": "bigwig_chromosome_sizes",
        "files": [str(item) for item in sources],
        "chromosomes": [{"chrom": name, "length": int(length)} for name, length in sorted(common)],
        "non_matching": [
            {"chrom": name, "length": int(length)} for name, length in sorted(non_common)
        ],
        "provenance": _provenance(
            "deeptools",
            "chromosome names and lengths read with pyBigWig from each file",
        ),
    }


__all__ = [
    "bedgraph_to_bigwig",
    "bigwig_chromosome_sizes",
    "chop_regions_for_profile",
    "estimate_bam_scale_factors",
    "estimate_count_size_factors",
    "fingerprint_jensen_shannon",
    "gc_fraction_from_twobit",
    "map_twobit_to_bam_chromosomes",
    "match_chromosome_name",
    "parse_deeptools_region",
]

_register_public_verbose(globals(), tuple(__all__))
