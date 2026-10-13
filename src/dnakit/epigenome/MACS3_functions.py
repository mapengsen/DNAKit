"""Optional MACS3 wrappers for fragment pileups, peak summits, tags and variant likelihoods.

The functions call MACS3's Cython signal, peak and statistics functions directly. They do not
run the MACS3 command line, so no peak-calling files are written. Packages are imported only
when a function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-macs3]"`` hint.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MACS3_EXTRA = "bio-macs3"
_MAX_POSITIONS = 10_000_000
_MAX_SIGNAL = 10_000_000
_MAX_COORDINATE = 2_147_483_647
_MAX_SHIFT = 100_000
_MAX_TEXT = 1_000_000
_MAX_SAMPLES = 1_000_000
_HMM_CLASS_COUNT = 4
_HMM_CLASSES = ("short", "mono", "di", "tri")
_REDUCERS = ("max", "min", "mean")
_MODELS = ("homozygous", "heterozygous")


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


def _text(value: object, name: str, limit: int = _MAX_TEXT) -> str:
    if (
        not isinstance(value, str)
        or not value.strip()
        or len(value) > limit
        or any(char in value for char in "\x00\r\n\t")
    ):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _sequence(value: object, name: str, maximum: int) -> Sequence[object]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} may contain at most {maximum} values.")
    return value


def _positions(value: object, name: str, maximum: int = _MAX_POSITIONS) -> list[int]:
    items = _sequence(value, name, maximum)
    return [_integer(item, f"{name} entries", 0, _MAX_COORDINATE) for item in items]


def _signal(value: object, name: str) -> list[float]:
    items = _sequence(value, name, _MAX_SIGNAL)
    if len(items) < 3:
        raise ConfigurationError(f"{name} must contain at least 3 values.")
    return [_number(item, f"{name} entries", -1e12, 1e12) for item in items]


def _require_sorted(values: list[int], name: str) -> None:
    if any(later < earlier for earlier, later in itertools.pairwise(values)):
        raise ConfigurationError(f"{name} must be sorted in ascending order.")


def _plain(value: Any) -> Any:
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


def _pileup_result(pileup: Any) -> dict[str, list[Any]]:
    """Convert a MACS3 legacy ``[positions, values]`` pair into plain lists."""

    positions, values = pileup
    return {
        "positions": [int(_plain(item)) for item in positions],
        "values": [float(_plain(item)) for item in values],
    }


def build_paired_fragment_pileup(
    left_positions: Sequence[int],
    right_positions: Sequence[int],
    *,
    scale_factor: float = 1.0,
    baseline_value: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up paired-end fragments given their left and right coordinates."""

    left = _positions(left_positions, "left_positions")
    right = _positions(right_positions, "right_positions")
    if len(left) != len(right) or not left:
        raise ConfigurationError(
            "left_positions and right_positions must have equal, non-zero length."
        )
    if any(end <= start for start, end in zip(left, right, strict=True)):
        raise ConfigurationError(
            "Every fragment must have its right coordinate after its left one."
        )
    scale = _number(scale_factor, "scale_factor", 0.0, 1_000_000.0)
    baseline = _number(baseline_value, "baseline_value", -1_000_000.0, 1_000_000.0)
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.PileupV2", _MACS3_EXTRA)
    fragments = numpy.array(list(zip(left, right, strict=True)), dtype=[("l", "<i4"), ("r", "<i4")])
    pileup = pileup_module.pileup_from_LR_as_list(fragments, scale, baseline)
    return {
        "status": "completed",
        "operation": "build_paired_fragment_pileup",
        "fragments": len(left),
        **_pileup_result(pileup),
        "provenance": _provenance(
            "MACS3",
            "count pileup of fragments; positions are segment ends, values cover the segment "
            "that ends at each position",
        ),
    }


def build_shifted_read_pileup(
    plus_positions: Sequence[int],
    minus_positions: Sequence[int],
    *,
    five_shift: int,
    three_shift: int,
    chrom_length: int,
    scale_factor: float = 1.0,
    baseline_value: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up single-end reads as shifted fragments on each strand."""

    plus = _positions(plus_positions, "plus_positions")
    minus = _positions(minus_positions, "minus_positions")
    if not plus and not minus:
        raise ConfigurationError("plus_positions or minus_positions must not both be empty.")
    shift_five = _integer(five_shift, "five_shift", -_MAX_SHIFT, _MAX_SHIFT)
    shift_three = _integer(three_shift, "three_shift", -_MAX_SHIFT, _MAX_SHIFT)
    length = _integer(chrom_length, "chrom_length", 1, _MAX_COORDINATE)
    scale = _number(scale_factor, "scale_factor", 0.0, 1_000_000.0)
    baseline = _number(baseline_value, "baseline_value", -1_000_000.0, 1_000_000.0)
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.PileupV2", _MACS3_EXTRA)
    plus_array = numpy.array(plus, dtype=numpy.int32)
    minus_array = numpy.array(minus, dtype=numpy.int32)
    pileup = pileup_module.pileup_from_PN_shifted(
        plus_array, minus_array, shift_five, shift_three, length, scale, baseline
    )
    return {
        "status": "completed",
        "operation": "build_shifted_read_pileup",
        "plus_reads": len(plus),
        "minus_reads": len(minus),
        "five_shift": shift_five,
        "three_shift": shift_three,
        "chrom_length": length,
        **_pileup_result(pileup),
        "provenance": _provenance(
            "MACS3",
            "plus reads cover [p - five_shift, p + three_shift), minus reads cover "
            "[p - three_shift, p + five_shift); coordinates are clipped to the chromosome",
        ),
    }


def find_signal_maxima(
    signal: Sequence[float], *, window_size: int = 51, verbose: bool = True
) -> dict[str, Any]:
    """Find local maxima of a signal after a second-order Savitzky-Golay smoothing."""

    values = _signal(signal, "signal")
    window = _integer(window_size, "window_size", 3, _MAX_SIGNAL)
    if window % 2 == 0 or window > len(values):
        raise ConfigurationError("window_size must be odd and no longer than the signal.")
    numpy = require("numpy", _MACS3_EXTRA)
    processing = require("MACS3.Signal.SignalProcessing", _MACS3_EXTRA)
    maxima = processing.maxima(numpy.asarray(values, dtype=numpy.float32), window)
    indices = [int(item) for item in maxima]
    return {
        "status": "completed",
        "operation": "find_signal_maxima",
        "signal_length": len(values),
        "window_size": window,
        "maxima": indices,
        "provenance": _provenance(
            "MACS3", "local maxima after second-order Savitzky-Golay smoothing", note="0-based"
        ),
    }


def _summit_indices(values: object, signal_length: int, name: str) -> list[int]:
    summits = _positions(values, name, maximum=_MAX_SIGNAL)
    if any(item >= signal_length for item in summits):
        raise ConfigurationError(f"{name} must index positions inside the signal.")
    _require_sorted(summits, name)
    if len(set(summits)) != len(summits):
        raise ConfigurationError(f"{name} must not repeat an index.")
    return summits


def filter_summits_by_peakiness(
    signal: Sequence[float], summits: Sequence[int], *, verbose: bool = True
) -> dict[str, Any]:
    """Keep summits whose surrounding peak meets MACS3's width and shape criteria."""

    values = _signal(signal, "signal")
    indices = _summit_indices(summits, len(values), "summits")
    numpy = require("numpy", _MACS3_EXTRA)
    processing = require("MACS3.Signal.SignalProcessing", _MACS3_EXTRA)
    kept = processing.enforce_peakyness(
        numpy.asarray(values, dtype=numpy.float32), numpy.asarray(indices, dtype=numpy.int32)
    )
    return {
        "status": "completed",
        "operation": "filter_summits_by_peakiness",
        "input_summits": len(indices),
        "summits": [int(item) for item in kept],
        "provenance": _provenance(
            "MACS3",
            "keeps summits whose peak is at least 50 bp wide and has at least 6 distinct values",
            note="0-based indices into the signal",
        ),
    }


def filter_summits_by_valleys(
    signal: Sequence[float],
    summits: Sequence[int],
    *,
    min_valley: float = 0.8,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep summits whose neighbouring valley is low enough compared with the lower summit."""

    values = _signal(signal, "signal")
    indices = _summit_indices(summits, len(values), "summits")
    ratio = _number(min_valley, "min_valley", 0.0, 1.0)
    numpy = require("numpy", _MACS3_EXTRA)
    processing = require("MACS3.Signal.SignalProcessing", _MACS3_EXTRA)
    kept = processing.enforce_valleys(
        numpy.asarray(values, dtype=numpy.float32),
        numpy.asarray(indices, dtype=numpy.int32),
        ratio,
    )
    return {
        "status": "completed",
        "operation": "filter_summits_by_valleys",
        "input_summits": len(indices),
        "min_valley": ratio,
        "summits": [int(item) for item in kept],
        "provenance": _provenance(
            "MACS3",
            "a valley between two summits must be at most min_valley times the lower summit",
            note="0-based indices into the signal",
        ),
    }


def refine_peak_summit(
    chrom: str,
    plus_positions: Sequence[int],
    minus_positions: Sequence[int],
    peak_start: int,
    peak_end: int,
    *,
    name: str = "peak",
    window_size: int = 100,
    cutoff: float = 5.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate the summit of one called peak from the 5' read positions on each strand."""

    chromosome = _text(chrom, "chrom", limit=255)
    plus = _positions(plus_positions, "plus_positions")
    minus = _positions(minus_positions, "minus_positions")
    start = _integer(peak_start, "peak_start", 0, _MAX_COORDINATE)
    end = _integer(peak_end, "peak_end", 1, _MAX_COORDINATE)
    if end <= start:
        raise ConfigurationError("peak_end must be greater than peak_start.")
    label = _text(name, "name", limit=128)
    window = _integer(window_size, "window_size", 1, 10_000)
    threshold = _number(cutoff, "cutoff", 0.0, 1_000.0)
    numpy = require("numpy", _MACS3_EXTRA)
    refine = require("MACS3.Commands.refinepeak_cmd", _MACS3_EXTRA)
    result = refine.find_summit(
        chromosome.encode("ascii"),
        numpy.asarray(plus, dtype=numpy.int32),
        numpy.asarray(minus, dtype=numpy.int32),
        start,
        end,
        label.encode("ascii"),
        window,
        threshold,
    )
    if result is None:
        raise ConfigurationError("MACS3 found no summit inside the peak.")
    chrom_out, summit_start, summit_end, name_out, score = result
    return {
        "status": "completed",
        "operation": "refine_peak_summit",
        "chrom": chrom_out.decode("ascii") if isinstance(chrom_out, bytes) else chrom_out,
        "summit_start": int(summit_start),
        "summit_end": int(summit_end),
        "name": name_out.decode("ascii") if isinstance(name_out, bytes) else name_out,
        "score": _plain(score),
        "provenance": _provenance(
            "MACS3",
            "the summit is the position with the most 5' read ends after window smoothing",
            note="coordinates are as returned by MACS3 (refinepeak)",
        ),
    }


def max_duplicate_tags(
    genome_size: int, tag_count: int, *, p_value: float = 1e-5, verbose: bool = True
) -> dict[str, Any]:
    """Largest number of identical tags expected by chance at a given P-value."""

    size = _integer(genome_size, "genome_size", 1, 10_000_000_000)
    tags = _integer(tag_count, "tag_count", 1, 10_000_000_000)
    p = _number(p_value, "p_value", 1e-300, 0.999999)
    callpeak = require("MACS3.Commands.callpeak_cmd", _MACS3_EXTRA)
    limit = callpeak.cal_max_dup_tags(size, tags, p)
    return {
        "status": "completed",
        "operation": "max_duplicate_tags",
        "genome_size": size,
        "tag_count": tags,
        "p_value": p,
        "max_duplicate_tags": int(limit),
        "provenance": _provenance(
            "MACS3",
            "reverse binomial CDF of the tags per position, computed by brute force",
        ),
    }


def guess_tag_file_format(
    path: str | Path, *, buffer_size: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """Identify the tag file format MACS3 would read, from the start of the file."""

    source = input_path(path)
    buffer = _integer(buffer_size, "buffer_size", 1_000, 10_000_000)
    parser_module = require("MACS3.IO.Parser", _MACS3_EXTRA)
    try:
        parser = parser_module.guess_parser(str(source), buffer)
    except Exception as exc:
        raise ConfigurationError(
            f"MACS3 could not recognise the tag file format of {source.name}: {exc}"
        ) from exc
    class_name = type(parser).__name__
    return {
        "status": "completed",
        "operation": "guess_tag_file_format",
        "path": str(source),
        "parser": class_name,
        "format": class_name[: -len("Parser")] if class_name.endswith("Parser") else class_name,
        "provenance": _provenance("MACS3", "first parser that accepts the leading records"),
    }


def calculate_genotype_quality(
    log_likelihood_1: float,
    log_likelihood_2: float,
    log_likelihood_3: float,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phred-scaled genotype quality from the natural-log likelihoods of three genotypes."""

    values = [
        _number(item, name, -1e9, 1e9)
        for item, name in (
            (log_likelihood_1, "log_likelihood_1"),
            (log_likelihood_2, "log_likelihood_2"),
            (log_likelihood_3, "log_likelihood_3"),
        )
    ]
    variant = require("MACS3.Signal.VariantStat", _MACS3_EXTRA)
    quality = variant.calculate_GQ(*values)
    return {
        "status": "completed",
        "operation": "calculate_genotype_quality",
        "genotype_quality": int(quality),
        "provenance": _provenance(
            "MACS3",
            "GQ1 = -10 log10((L2 + L3) / (L1 + L2 + L3)) with L = exp(log-likelihood)",
            note="MACS3 returns the value as an integer",
        ),
    }


def _count_array(value: object, name: str, numpy: Any) -> Any:
    items = _sequence(value, name, _MAX_SAMPLES)
    if len(items) < 1:
        raise ConfigurationError(f"{name} must contain at least one read.")
    counts = [_integer(item, f"{name} entries", 0, 1_000_000) for item in items]
    return numpy.array(counts, dtype=numpy.int32)


def calculate_allele_model_likelihood(
    top1_bq_t: Sequence[int],
    top1_bq_c: Sequence[int],
    top2_bq_t: Sequence[int],
    top2_bq_c: Sequence[int],
    *,
    model: str = "homozygous",
    verbose: bool = True,
) -> dict[str, Any]:
    """Log-likelihood and BIC of a homozygous or heterozygous allele model at one site."""

    if model not in _MODELS:
        raise ConfigurationError("model must be 'homozygous' or 'heterozygous'.")
    numpy = require("numpy", _MACS3_EXTRA)
    arrays = [
        _count_array(top1_bq_t, "top1_bq_t", numpy),
        _count_array(top1_bq_c, "top1_bq_c", numpy),
        _count_array(top2_bq_t, "top2_bq_t", numpy),
        _count_array(top2_bq_c, "top2_bq_c", numpy),
    ]
    if any(len(array) != len(arrays[0]) for array in arrays):
        raise ConfigurationError("All four base-quality arrays must have the same length.")
    variant = require("MACS3.Signal.VariantStat", _MACS3_EXTRA)
    if model == "homozygous":
        log_likelihood, bic = variant.CalModel_Homo(*arrays)
    else:
        log_likelihood, bic = variant.CalModel_Heter_noAS(*arrays)
    return {
        "status": "completed",
        "operation": "calculate_allele_model_likelihood",
        "model": model,
        "reads": len(arrays[0]),
        "log_likelihood": _plain(log_likelihood),
        "bic": _plain(bic),
        "provenance": _provenance(
            "MACS3",
            "allele-count likelihood of a single-site genotype model, with BIC",
        ),
    }


def hmmratac_fragment_weights(
    fragment_lengths: Sequence[int],
    means: Sequence[float],
    stddevs: Sequence[float],
    *,
    min_fragment_probability: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weight of each fragment length for the HMMRATAC short, mono, di and tri signal tracks."""

    lengths = _positions(fragment_lengths, "fragment_lengths", maximum=10_000)
    if not lengths or any(item < 1 for item in lengths):
        raise ConfigurationError("fragment_lengths must be positive integers.")
    mean_values = _sequence(means, "means", _HMM_CLASS_COUNT)
    sd_values = _sequence(stddevs, "stddevs", _HMM_CLASS_COUNT)
    if len(mean_values) != _HMM_CLASS_COUNT or len(sd_values) != _HMM_CLASS_COUNT:
        raise ConfigurationError("means and stddevs must have four values: short, mono, di, tri.")
    mu = [_number(item, "means entries", 0.0, 1_000_000.0) for item in mean_values]
    sigma = [_number(item, "stddevs entries", 1e-9, 1_000_000.0) for item in sd_values]
    floor = _number(min_fragment_probability, "min_fragment_probability", 0.0, 1.0)
    hmm = require("MACS3.Signal.HMMR_Signal_Processing", _MACS3_EXTRA)
    mappings = hmm.generate_weight_mapping(sorted(set(lengths)), mu, sigma, floor)
    rows = []
    for length in sorted(set(lengths)):
        row: dict[str, Any] = {"fragment_length": length}
        for label, mapping in zip(_HMM_CLASSES, mappings, strict=True):
            row[label] = _plain(mapping.get(length, 0.0))
        rows.append(row)
    return {
        "status": "completed",
        "operation": "hmmratac_fragment_weights",
        "classes": list(_HMM_CLASSES),
        "weights": rows,
        "provenance": _provenance(
            "MACS3",
            "Gaussian fragment-length mixture weight per HMMRATAC signal class",
            note="a length absent from the upstream mapping is reported as 0.0",
        ),
    }


def merge_pileup_arrays(
    first: Mapping[str, Sequence[float]],
    second: Mapping[str, Sequence[float]],
    *,
    reducer: str = "max",
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge two pileup tracks (positions and values) with a pointwise reducer."""

    if reducer not in _REDUCERS:
        raise ConfigurationError("reducer must be 'max', 'min' or 'mean'.")
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.Pileup", _MACS3_EXTRA)
    pairs = []
    for label, track in (("first", first), ("second", second)):
        if not isinstance(track, Mapping) or set(track) != {"positions", "values"}:
            raise ConfigurationError(f"{label} must be a mapping with positions and values.")
        positions = _positions(track["positions"], f"{label} positions")
        values = _sequence(track["values"], f"{label} values", _MAX_POSITIONS)
        if len(values) != len(positions) or not positions:
            raise ConfigurationError(f"{label} positions and values must match and be non-empty.")
        numbers = [_number(item, f"{label} values entries", -1e12, 1e12) for item in values]
        pairs.append(
            [
                numpy.asarray(positions, dtype=numpy.int32),
                numpy.asarray(numbers, dtype=numpy.float32),
            ]
        )
    merged = pileup_module.over_two_pv_array(pairs[0], pairs[1], reducer)
    return {
        "status": "completed",
        "operation": "merge_pileup_arrays",
        "reducer": reducer,
        **_pileup_result(merged),
        "provenance": _provenance(
            "MACS3",
            "pointwise reducer over the union of two [positions, values] tracks",
        ),
    }


__all__ = [
    "build_paired_fragment_pileup",
    "build_shifted_read_pileup",
    "calculate_allele_model_likelihood",
    "calculate_genotype_quality",
    "filter_summits_by_peakiness",
    "filter_summits_by_valleys",
    "find_signal_maxima",
    "guess_tag_file_format",
    "hmmratac_fragment_weights",
    "max_duplicate_tags",
    "merge_pileup_arrays",
    "refine_peak_summit",
]

_register_public_verbose(globals(), tuple(__all__))
