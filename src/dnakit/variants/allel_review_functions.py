"""Offline scikit-allel wrappers for population-genetic statistics, windows and archives.

Inputs are plain Python lists (or NumPy arrays) of allele counts, phased haplotypes, genotype
allele counts or site frequency spectra, all from local data. Positions are 1-based and must be
sorted. Window bounds are 1-based and inclusive; moving windows count variants, not bases.
Statistics that cannot be computed for a window are returned as ``None``.

scikit-allel is imported only when a function runs; a missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-scikit-allel]"`` hint. HDF5
archives also need h5py, which allel installs as a dependency.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_ALLEL_EXTRA = "bio-scikit-allel"
_NUMPY_EXTRA = "bio-analysis"
_MAX_VARIANTS = 1_000_000
_MAX_COUNT = 2_000_000
_MAX_POSITION = 2_000_000_000
_MAX_WINDOW = 1_000_000_000
_MAX_SAMPLE_SIZE = 2_000
_MAX_SAMPLES = 500
_MAX_POPULATIONS = 50
_MAX_FIELDS = 100
_CONTIG_LINE = re.compile(r"^##contig=<(.*)>\s*$")


def _numpy() -> Any:
    return require("numpy", _NUMPY_EXTRA)


def _allel() -> Any:
    return require("allel", _ALLEL_EXTRA)


def _provenance() -> dict[str, str]:
    return {
        "package": f"scikit-allel {metadata.version('scikit-allel')}",
        "reference": "none: offline arrays or local files",
    }


def _clean(value: float) -> float | None:
    if value != value or value in (float("inf"), float("-inf")):
        return None
    return round(float(value), 10) + 0.0


def _values(value: Any) -> list[Any]:
    """Flatten a NumPy array into a JSON list with None for undefined statistics."""
    np = _numpy()
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return _clean(float(array))
    return [_values(item) if array.ndim > 1 else _clean(float(item)) for item in array]


def _int_array(
    value: object,
    name: str,
    *,
    ndim: int,
    minimum: int = 0,
    maximum: int = _MAX_COUNT,
    allow_missing: bool = False,
) -> Any:
    np = _numpy()
    try:
        array = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular array.") from exc
    if array.dtype.kind not in "iu" or array.ndim != ndim or array.size == 0:
        raise ConfigurationError(f"{name} must be a non-empty {ndim}-dimensional integer array.")
    low = -1 if allow_missing else minimum
    if (array < low).any() or (array > maximum).any():
        raise ConfigurationError(f"{name} must hold integers in [{low}, {maximum}].")
    return array.astype(np.int64)


def _float_array(value: object, name: str, *, ndim: int, allow_nan: bool = False) -> Any:
    np = _numpy()
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a numeric array.") from exc
    if array.ndim != ndim or array.size == 0:
        raise ConfigurationError(f"{name} must be a non-empty {ndim}-dimensional array.")
    finite = np.isfinite(array) | (np.isnan(array) if allow_nan else False)
    if not finite.all():
        raise ConfigurationError(f"{name} must contain finite numbers.")
    return array


def _bool_array(value: object, name: str) -> Any:
    np = _numpy()
    array = np.asarray(value)
    if array.dtype != bool or array.ndim != 1 or array.size == 0:
        raise ConfigurationError(f"{name} must be a non-empty one-dimensional boolean array.")
    return array


def _positions(value: object, name: str = "positions") -> Any:
    np = _numpy()
    array = _int_array(value, name, ndim=1, minimum=1, maximum=_MAX_POSITION)
    if array.shape[0] > _MAX_VARIANTS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_VARIANTS} variants.")
    if (np.diff(array) < 0).any():
        raise ConfigurationError(f"{name} must be sorted in ascending order.")
    return array


def _counts(value: object, name: str) -> Any:
    array = _int_array(value, name, ndim=2, maximum=_MAX_COUNT)
    if array.shape[1] < 2 or array.shape[0] > _MAX_VARIANTS:
        raise ConfigurationError(
            f"{name} must have shape (variants, alleles) with at least two alleles."
        )
    return array


def _same_rows(arrays: Sequence[Any], names: Sequence[str]) -> None:
    rows = arrays[0].shape[0]
    for array, label in zip(arrays[1:], names[1:], strict=False):
        if array.shape[0] != rows:
            raise ConfigurationError(f"{label} must describe the same number of variants.")


def _haplotypes(value: object, name: str) -> Any:
    allel = _allel()
    array = _int_array(value, name, ndim=2, minimum=0, maximum=1)
    if array.shape[1] < 2:
        raise ConfigurationError(f"{name} must have at least two haplotypes.")
    return allel.HaplotypeArray(array.astype("int8"))


def _window_size(value: object, name: str = "window_size") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_WINDOW:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _optional_window(value: object, name: str) -> int | None:
    return None if value is None else _window_size(value, name)


def _optional_position(value: object, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_POSITION:
        raise ConfigurationError(f"{name} must be a positive integer position.")
    return value


def _index_bound(value: object, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _MAX_VARIANTS:
        raise ConfigurationError(f"{name} must be a non-negative integer index.")
    return value


def _bounds(left: int | None, right: int | None) -> None:
    if left is not None and right is not None and right < left:
        raise ConfigurationError("stop must be at least start.")


def _window_rows(windows: Any, label: str, values: Any, extra: Any = None) -> list[dict[str, Any]]:
    np = _numpy()
    rows = []
    listed_values = np.asarray(values).tolist()
    listed_extra = None if extra is None else np.asarray(extra).tolist()
    for index, (start, stop) in enumerate(np.asarray(windows).tolist()):
        row: dict[str, Any] = {"start": int(start), "stop": int(stop)}
        if listed_extra is not None:
            row["sites"] = int(listed_extra[index])
        row[label] = _clean(float(listed_values[index]))
        rows.append(row)
    return rows


def _square(condensed: Any, samples: int) -> list[list[float | None]]:
    np = _numpy()
    square = np.zeros((samples, samples))
    upper = np.triu_indices(samples, k=1)
    square[upper] = np.asarray(condensed, dtype=float)
    square = square + square.T
    return _values(square)


def _fresh_file(value: object, name: str, suffix: str | None = None) -> Path:
    if not isinstance(value, (str, Path)) or not str(value):
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if suffix is not None and target.suffix != suffix:
        raise ConfigurationError(f"{name} must end with {suffix}.")
    if target.exists():
        raise ConfigurationError(f"{name} already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"{name} directory does not exist: {target.parent}")
    return target


def _fields(value: object) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("fields must be a non-empty list of allel field names.")
    if len(value) > _MAX_FIELDS or not all(isinstance(item, str) and item for item in value):
        raise ConfigurationError(f"fields must hold at most {_MAX_FIELDS} non-empty names.")
    return list(value)


def locate_population_private_alleles(
    allele_counts: Sequence[Sequence[Sequence[int]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Alleles present in exactly one of several populations, from per-population allele counts."""

    allel = _allel()
    if isinstance(allele_counts, (str, bytes)) or not 2 <= len(allele_counts) <= _MAX_POPULATIONS:
        raise ConfigurationError(f"allele_counts must hold 2 to {_MAX_POPULATIONS} populations.")
    arrays = [_counts(item, f"allele_counts[{index}]") for index, item in enumerate(allele_counts)]
    _same_rows(arrays, [f"allele_counts[{index}]" for index in range(len(arrays))])
    with _quiet():
        flags = allel.locate_private_alleles(*[allel.AlleleCountsArray(a) for a in arrays])
    private = [
        {"variant_index": int(variant), "allele_index": int(allele)}
        for variant, allele in zip(*_nonzero(flags), strict=False)
    ]
    return {
        "status": "completed",
        "operation": "locate_population_private_alleles",
        "populations": len(arrays),
        "variants": int(arrays[0].shape[0]),
        "private_alleles": private,
        "private_variant_indices": sorted({item["variant_index"] for item in private}),
        "provenance": _provenance(),
    }


def _nonzero(flags: Any) -> tuple[list[int], list[int]]:
    np = _numpy()
    rows, columns = np.nonzero(np.asarray(flags))
    return rows.tolist(), columns.tolist()


def _quiet() -> Any:
    np = _numpy()
    return np.errstate(divide="ignore", invalid="ignore")


def locate_population_fixed_differences(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variants with no allele shared between two populations (fixed differences)."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    flags = allel.locate_fixed_differences(
        allel.AlleleCountsArray(first), allel.AlleleCountsArray(second)
    )
    indices = [int(index) for index, flag in enumerate(flags.tolist()) if flag]
    return {
        "status": "completed",
        "operation": "locate_population_fixed_differences",
        "variants": int(first.shape[0]),
        "fixed_indices": indices,
        "fixed_count": len(indices),
        "provenance": _provenance(),
    }


def count_variants_in_windows(
    positions: Sequence[int],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of variants in each window of ``window_size`` base pairs over one contig."""

    allel = _allel()
    pos = _positions(positions)
    size = _window_size(window_size)
    left = _optional_position(start, "start")
    right = _optional_position(stop, "stop")
    stride = _optional_window(step, "step")
    _bounds(left, right)
    with _quiet():
        counts, windows = allel.windowed_count(
            allel.SortedIndex(pos), size=size, start=left, stop=right, step=stride
        )
    rows = [
        {"start": int(start_pos), "stop": int(stop_pos), "variants": int(count)}
        for (start_pos, stop_pos), count in zip(windows.tolist(), counts.tolist(), strict=False)
    ]
    return {
        "status": "completed",
        "operation": "count_variants_in_windows",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def windowed_values_per_base(
    values: Sequence[float],
    windows: Sequence[Sequence[int]],
    *,
    is_accessible: Sequence[bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Divide one windowed statistic by the number of bases in each window (or accessible bases)."""

    allel = _allel()
    np = _numpy()
    vals = _float_array(values, "values", ndim=1)
    bounds = _int_array(windows, "windows", ndim=2, minimum=1, maximum=_MAX_POSITION)
    if bounds.shape[1] != 2 or bounds.shape[0] != vals.shape[0]:
        raise ConfigurationError("windows must have one (start, stop) row per value.")
    if (bounds[:, 1] < bounds[:, 0]).any():
        raise ConfigurationError("each window stop must be at least its start.")
    access = None
    if is_accessible is not None:
        access = _bool_array(is_accessible, "is_accessible")
        if access.shape[0] < int(bounds[:, 1].max()):
            raise ConfigurationError("is_accessible must cover every window stop.")
    with _quiet():
        per_base, bases = allel.per_base(vals, bounds, is_accessible=access, fill=np.nan)
    rows = [
        {
            "start": int(start),
            "stop": int(stop),
            "callable_bases": int(count),
            "value_per_base": _clean(float(value)),
        }
        for (start, stop), count, value in zip(
            bounds.tolist(), np.asarray(bases).tolist(), np.asarray(per_base).tolist(), strict=False
        )
    ]
    return {
        "status": "completed",
        "operation": "windowed_values_per_base",
        "windows": rows,
        "provenance": _provenance(),
    }


def build_equally_accessible_windows(
    is_accessible: Sequence[bool],
    *,
    window_size: int,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Windows that each hold the same number of accessible bases, from a base-1 mask."""

    allel = _allel()
    np = _numpy()
    access = _bool_array(is_accessible, "is_accessible")
    size = _window_size(window_size)
    stride = _optional_window(step, "step")
    windows = allel.equally_accessible_windows(access, size, step=stride)
    rows = [
        {
            "start": int(start),
            "stop": int(stop),
            "accessible_bases": int(np.count_nonzero(access[start - 1 : stop])),
        }
        for start, stop in np.asarray(windows).tolist()
    ]
    return {
        "status": "completed",
        "operation": "build_equally_accessible_windows",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def calculate_between_population_pi(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean number of pairwise differences between two populations, for each variant."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    with _quiet():
        values = allel.mean_pairwise_difference_between(
            allel.AlleleCountsArray(first), allel.AlleleCountsArray(second)
        )
    return {
        "status": "completed",
        "operation": "calculate_between_population_pi",
        "variants": int(first.shape[0]),
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_windowed_fixed_difference_density(
    positions: Sequence[int],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    is_accessible: Sequence[bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Density of fixed differences between two populations in windows over one contig."""

    allel = _allel()
    pos = _positions(positions)
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    if first.shape[0] != pos.shape[0]:
        raise ConfigurationError("first_counts must have one row per position.")
    size = _window_size(window_size)
    left = _optional_position(start, "start")
    right = _optional_position(stop, "stop")
    stride = _optional_window(step, "step")
    _bounds(left, right)
    access = None
    if is_accessible is not None:
        access = _bool_array(is_accessible, "is_accessible")
        if access.shape[0] < int(pos[-1]):
            raise ConfigurationError("is_accessible must cover every position.")
    with _quiet():
        density, windows, bases, counts = allel.windowed_df(
            allel.SortedIndex(pos),
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            size=size,
            start=left,
            stop=right,
            step=stride,
            is_accessible=access,
        )
    rows = _window_rows(windows, "density", density, counts)
    for row, callable_bases in zip(rows, _numpy().asarray(bases).tolist(), strict=False):
        row["callable_bases"] = int(callable_bases)
    return {
        "status": "completed",
        "operation": "calculate_windowed_fixed_difference_density",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def calculate_patterson_fst(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's Fst between two populations: per-variant parts and their ratio."""

    allel = _allel()
    np = _numpy()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    with _quiet():
        numerator, denominator = allel.patterson_fst(
            allel.AlleleCountsArray(first), allel.AlleleCountsArray(second)
        )
    total = float(np.nansum(denominator))
    return {
        "status": "completed",
        "operation": "calculate_patterson_fst",
        "variants": int(first.shape[0]),
        "numerator": _values(numerator),
        "denominator": _values(denominator),
        "fst": None if total == 0 else _clean(float(np.nansum(numerator)) / total),
        "provenance": _provenance(),
    }


def calculate_windowed_patterson_fst(
    positions: Sequence[int],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's Fst averaged in windows of ``window_size`` base pairs over one contig."""

    allel = _allel()
    pos = _positions(positions)
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    if first.shape[0] != pos.shape[0]:
        raise ConfigurationError("first_counts must have one row per position.")
    size = _window_size(window_size)
    left = _optional_position(start, "start")
    right = _optional_position(stop, "stop")
    stride = _optional_window(step, "step")
    _bounds(left, right)
    with _quiet():
        fst, windows, counts = allel.windowed_patterson_fst(
            allel.SortedIndex(pos),
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            size=size,
            start=left,
            stop=right,
            step=stride,
        )
    return {
        "status": "completed",
        "operation": "calculate_windowed_patterson_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": _window_rows(windows, "fst", fst, counts),
        "provenance": _provenance(),
    }


def calculate_blockwise_patterson_fst(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    block_length: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's Fst with a block-jackknife standard error over blocks of variants."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    blen = _window_size(block_length, "block_length")
    with _quiet():
        fst, error, block_values, jackknife = allel.blockwise_patterson_fst(
            allel.AlleleCountsArray(first), allel.AlleleCountsArray(second), blen
        )
    return {
        "status": "completed",
        "operation": "calculate_blockwise_patterson_fst",
        "variants": int(first.shape[0]),
        "block_length": blen,
        "fst": _clean(float(fst)),
        "standard_error": _clean(float(error)),
        "block_values": _values(block_values),
        "jackknife_values": _values(jackknife),
        "provenance": _provenance(),
    }


def calculate_moving_patterson_fst(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's Fst in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    size = _window_size(window_size)
    left = _index_bound(start, "start") or 0
    right = _index_bound(stop, "stop")
    stride = _optional_window(step, "step")
    with _quiet():
        values = allel.moving_patterson_fst(
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            size=size,
            start=left,
            stop=right,
            step=stride,
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_patterson_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_pairwise_dxy(
    positions: Sequence[int],
    genotype_allele_counts: Sequence[Sequence[Sequence[int]]],
    *,
    start: int | None = None,
    stop: int | None = None,
    is_accessible: Sequence[bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise nucleotide divergence (dxy) between samples, as a square matrix over one region."""

    allel = _allel()
    pos = _positions(positions)
    gac = _int_array(genotype_allele_counts, "genotype_allele_counts", ndim=3, maximum=_MAX_COUNT)
    if gac.shape[0] != pos.shape[0] or gac.shape[2] < 2:
        raise ConfigurationError(
            "genotype_allele_counts must have shape (variants, samples, alleles)."
        )
    if not 2 <= gac.shape[1] <= _MAX_SAMPLES:
        raise ConfigurationError(f"genotype_allele_counts must hold 2 to {_MAX_SAMPLES} samples.")
    left = _optional_position(start, "start")
    right = _optional_position(stop, "stop")
    _bounds(left, right)
    access = None
    if is_accessible is not None:
        access = _bool_array(is_accessible, "is_accessible")
        if access.shape[0] < int(pos[-1]):
            raise ConfigurationError("is_accessible must cover every position.")
    with _quiet():
        condensed = allel.pairwise_dxy(
            allel.SortedIndex(pos), gac, start=left, stop=right, is_accessible=access
        )
    samples = int(gac.shape[1])
    return {
        "status": "completed",
        "operation": "calculate_pairwise_dxy",
        "samples": samples,
        "matrix": _square(condensed, samples),
        "provenance": _provenance(),
    }


def calculate_patterson_f2(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's f2 branch-length estimate between two populations: per variant and its mean."""

    allel = _allel()
    np = _numpy()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    with _quiet():
        values = allel.patterson_f2(allel.AlleleCountsArray(first), allel.AlleleCountsArray(second))
    defined = np.asarray(values, dtype=float)
    defined = defined[~np.isnan(defined)]
    return {
        "status": "completed",
        "operation": "calculate_patterson_f2",
        "variants": int(first.shape[0]),
        "values": _values(values),
        "mean": _clean(float(defined.mean())) if defined.size else None,
        "provenance": _provenance(),
    }


def calculate_patterson_f3(
    target_counts: Sequence[Sequence[int]],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's f3(C; A, B) admixture test: summed numerator and heterozygosity, and f3*."""

    allel = _allel()
    np = _numpy()
    target = _counts(target_counts, "target_counts")
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([target, first, second], ["target_counts", "first_counts", "second_counts"])
    with _quiet():
        numerator, heterozygosity = allel.patterson_f3(
            allel.AlleleCountsArray(target),
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
        )
    t_sum = float(np.nansum(numerator))
    b_sum = float(np.nansum(heterozygosity))
    return {
        "status": "completed",
        "operation": "calculate_patterson_f3",
        "variants": int(target.shape[0]),
        "numerator": _values(numerator),
        "heterozygosity": _values(heterozygosity),
        "numerator_sum": _clean(t_sum),
        "heterozygosity_sum": _clean(b_sum),
        "f3_star": None if b_sum == 0 else _clean(t_sum / b_sum),
        "provenance": _provenance(),
    }


def calculate_blockwise_patterson_f3(
    target_counts: Sequence[Sequence[int]],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    block_length: int,
    normed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's f3 with a block-jackknife standard error and z-score."""

    allel = _allel()
    target = _counts(target_counts, "target_counts")
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([target, first, second], ["target_counts", "first_counts", "second_counts"])
    blen = _window_size(block_length, "block_length")
    if not isinstance(normed, bool):
        raise ConfigurationError("normed must be True or False.")
    with _quiet():
        f3, error, z_score, block_values, jackknife = allel.blockwise_patterson_f3(
            allel.AlleleCountsArray(target),
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            blen,
            normed=normed,
        )
    return {
        "status": "completed",
        "operation": "calculate_blockwise_patterson_f3",
        "variants": int(target.shape[0]),
        "block_length": blen,
        "normed": normed,
        "f3": _clean(float(f3)),
        "standard_error": _clean(float(error)),
        "z_score": _clean(float(z_score)),
        "block_values": _values(block_values),
        "jackknife_values": _values(jackknife),
        "provenance": _provenance(),
    }


def calculate_moving_patterson_f3(
    target_counts: Sequence[Sequence[int]],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    normed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's f3 in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    target = _counts(target_counts, "target_counts")
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([target, first, second], ["target_counts", "first_counts", "second_counts"])
    size = _window_size(window_size)
    left = _index_bound(start, "start") or 0
    right = _index_bound(stop, "stop")
    stride = _optional_window(step, "step")
    if not isinstance(normed, bool):
        raise ConfigurationError("normed must be True or False.")
    with _quiet():
        values = allel.moving_patterson_f3(
            allel.AlleleCountsArray(target),
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            size,
            start=left,
            stop=right,
            step=stride,
            normed=normed,
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_patterson_f3",
        "window_size": size,
        "step": size if stride is None else stride,
        "normed": normed,
        "values": _values(values),
        "provenance": _provenance(),
    }


def paint_shared_haplotype_prefixes(
    haplotypes: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Paint haplotypes so that each shared prefix gets one integer, after sorting by prefix."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    painting, order = allel.voight_painting(haps)
    return {
        "status": "completed",
        "operation": "paint_shared_haplotype_prefixes",
        "variants": int(haps.shape[0]),
        "haplotypes": int(haps.shape[1]),
        "painting": [[int(value) for value in row] for row in painting.tolist()],
        "haplotype_order": [int(index) for index in order.tolist()],
        "provenance": _provenance(),
    }


def calculate_haplotype_diversity(
    haplotypes: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Haplotype diversity: one minus the sum of squared haplotype frequencies, bias-corrected."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    value = allel.haplotype_diversity(haps)
    return {
        "status": "completed",
        "operation": "calculate_haplotype_diversity",
        "variants": int(haps.shape[0]),
        "haplotypes": int(haps.shape[1]),
        "diversity": _clean(float(value)),
        "provenance": _provenance(),
    }


def calculate_moving_haplotype_diversity(
    haplotypes: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Haplotype diversity in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    size = _window_size(window_size)
    left = _index_bound(start, "start") or 0
    right = _index_bound(stop, "stop")
    stride = _optional_window(step, "step")
    values = allel.moving_haplotype_diversity(haps, size, start=left, stop=right, step=stride)
    return {
        "status": "completed",
        "operation": "calculate_moving_haplotype_diversity",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_xpnsl(
    haplotypes_a: Sequence[Sequence[int]],
    haplotypes_b: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unstandardized cross-population nSL score per variant for two phased haplotype sets."""

    allel = _allel()
    first = _haplotypes(haplotypes_a, "haplotypes_a")
    second = _haplotypes(haplotypes_b, "haplotypes_b")
    if first.shape[0] != second.shape[0]:
        raise ConfigurationError("haplotypes_b must describe the same number of variants.")
    with _quiet():
        values = allel.xpnsl(first, second, use_threads=False)
    return {
        "status": "completed",
        "operation": "calculate_xpnsl",
        "variants": int(first.shape[0]),
        "xpnsl": _values(values),
        "provenance": _provenance(),
    }


def standardize_score_by_allele_count(
    score: Sequence[float],
    allele_counts: Sequence[int],
    *,
    n_bins: int | None = None,
    bins: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Standardize a selection score (such as iHS) within allele-count bins."""

    allel = _allel()
    np = _numpy()
    values = _float_array(score, "score", ndim=1, allow_nan=True)
    counts = _int_array(allele_counts, "allele_counts", ndim=1, maximum=_MAX_VARIANTS)
    if counts.shape[0] != values.shape[0]:
        raise ConfigurationError("allele_counts must have one entry per score.")
    if n_bins is not None and (
        isinstance(n_bins, bool) or not isinstance(n_bins, int) or not 2 <= n_bins <= 1000
    ):
        raise ConfigurationError("n_bins must be an integer in [2, 1000].")
    bin_array = None
    if bins is not None:
        bin_array = _int_array(bins, "bins", ndim=1, maximum=_MAX_VARIANTS)
        if (np.diff(bin_array) <= 0).any():
            raise ConfigurationError("bins must be strictly increasing.")
    with _quiet():
        standardized, used = allel.standardize_by_allele_count(
            values, counts, bins=bin_array, n_bins=n_bins, diagnostics=False
        )
    return {
        "status": "completed",
        "operation": "standardize_score_by_allele_count",
        "variants": int(values.shape[0]),
        "values": _values(standardized),
        "bins": [int(item) for item in np.asarray(used).tolist()],
        "provenance": _provenance(),
    }


def calculate_moving_delta_tajimas_d(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Difference in Tajima's D between two populations in moving variant windows."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    size = _window_size(window_size)
    left = _index_bound(start, "start") or 0
    right = _index_bound(stop, "stop")
    stride = _optional_window(step, "step")
    with _quiet():
        values = allel.moving_delta_tajima_d(
            allel.AlleleCountsArray(first),
            allel.AlleleCountsArray(second),
            size,
            start=left,
            stop=right,
            step=stride,
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_delta_tajimas_d",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def fold_joint_site_frequency_spectrum(
    spectrum: Sequence[Sequence[int]],
    *,
    first_sample_size: int,
    second_sample_size: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fold a joint site frequency spectrum of two populations onto the minor-allele corner."""

    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=2, maximum=_MAX_COUNT)
    n1 = _sample_size(first_sample_size, "first_sample_size")
    n2 = _sample_size(second_sample_size, "second_sample_size")
    _spectrum_fits(counts, n1, n2)
    folded = allel.fold_joint_sfs(counts, n1, n2)
    return {
        "status": "completed",
        "operation": "fold_joint_site_frequency_spectrum",
        "first_sample_size": n1,
        "second_sample_size": n2,
        "folded": [[int(item) for item in row] for row in folded.tolist()],
        "provenance": _provenance(),
    }


def scale_site_frequency_spectrum(
    spectrum: Sequence[int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale each class of a site frequency spectrum by its derived allele count (i times xi_i)."""

    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=1, maximum=_MAX_COUNT)
    scaled = allel.scale_sfs(counts)
    return {
        "status": "completed",
        "operation": "scale_site_frequency_spectrum",
        "classes": int(counts.shape[0]),
        "scaled": [int(item) for item in scaled.tolist()],
        "provenance": _provenance(),
    }


def scale_folded_site_frequency_spectrum(
    spectrum: Sequence[int],
    *,
    sample_size: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a folded site frequency spectrum of ``sample_size`` chromosomes by i(n-i)/n."""

    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=1, maximum=_MAX_COUNT)
    n = _sample_size(sample_size, "sample_size")
    if counts.shape[0] > n + 1:
        raise ConfigurationError("spectrum has more classes than sample_size allows.")
    scaled = allel.scale_sfs_folded(counts, n)
    return {
        "status": "completed",
        "operation": "scale_folded_site_frequency_spectrum",
        "sample_size": n,
        "scaled": _values(scaled),
        "provenance": _provenance(),
    }


def scale_joint_site_frequency_spectrum(
    spectrum: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale each cell of a joint site frequency spectrum by i times j."""

    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=2, maximum=_MAX_COUNT)
    scaled = allel.scale_joint_sfs(counts)
    return {
        "status": "completed",
        "operation": "scale_joint_site_frequency_spectrum",
        "shape": [int(item) for item in counts.shape],
        "scaled": [[int(item) for item in row] for row in scaled.tolist()],
        "provenance": _provenance(),
    }


def scale_folded_joint_site_frequency_spectrum(
    spectrum: Sequence[Sequence[int]],
    *,
    first_sample_size: int,
    second_sample_size: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a folded joint site frequency spectrum by i j (n1 - i)(n2 - j) per cell."""

    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=2, maximum=_MAX_COUNT)
    n1 = _sample_size(first_sample_size, "first_sample_size")
    n2 = _sample_size(second_sample_size, "second_sample_size")
    _spectrum_fits(counts, n1, n2)
    scaled = allel.scale_joint_sfs_folded(counts, n1, n2)
    return {
        "status": "completed",
        "operation": "scale_folded_joint_site_frequency_spectrum",
        "first_sample_size": n1,
        "second_sample_size": n2,
        "scaled": [[int(item) for item in row] for row in scaled.tolist()],
        "provenance": _provenance(),
    }


def _sample_size(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_SAMPLE_SIZE:
        raise ConfigurationError(f"{name} must be an integer in [1, {_MAX_SAMPLE_SIZE}].")
    return value


def _spectrum_fits(counts: Any, n1: int, n2: int) -> None:
    if counts.shape[0] > n1 + 1 or counts.shape[1] > n2 + 1:
        raise ConfigurationError("spectrum has more classes than the sample sizes allow.")


def paint_parental_transmission(
    parent_haplotypes: Sequence[Sequence[int]],
    progeny_haplotypes: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Paint progeny haplotypes by which parental haplotype each allele was inherited from."""

    allel = _allel()
    parent = _int_array(
        parent_haplotypes, "parent_haplotypes", ndim=2, maximum=1, allow_missing=True
    )
    progeny = _int_array(
        progeny_haplotypes, "progeny_haplotypes", ndim=2, maximum=1, allow_missing=True
    )
    if parent.shape[1] != 2:
        raise ConfigurationError("parent_haplotypes must hold the two haplotypes of one parent.")
    if progeny.shape[0] != parent.shape[0]:
        raise ConfigurationError("progeny_haplotypes must describe the same variants as parent.")
    painting = allel.paint_transmission(parent.astype("int8"), progeny.astype("int8"))
    return {
        "status": "completed",
        "operation": "paint_parental_transmission",
        "variants": int(parent.shape[0]),
        "progeny": int(progeny.shape[1]),
        "painting": [[int(code) for code in row] for row in painting.tolist()],
        "codes": {
            "0": "undetermined",
            "1": "first parental haplotype",
            "2": "second parental haplotype",
            "3": "reference allele on both parental haplotypes",
            "4": "non-reference allele on both parental haplotypes",
            "5": "non-parental allele",
            "6": "a parental allele is missing",
            "7": "progeny allele is missing",
        },
        "provenance": _provenance(),
    }


def read_vcf_header_summary(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Samples, INFO and FORMAT fields, filters and contig lengths in a local VCF header."""

    allel = _allel()
    source = input_path(path)
    headers = allel.read_vcf_headers(str(source))
    contigs = []
    for line in headers.headers:
        match = _CONTIG_LINE.match(line.strip())
        if match is None:
            continue
        fields = dict(part.split("=", 1) for part in match.group(1).split(",") if "=" in part)
        length = fields.get("length")
        contigs.append(
            {
                "id": fields.get("ID", ""),
                "length": int(length) if length and length.isdigit() else None,
            }
        )
    contigs.sort(key=lambda item: item["id"])
    return {
        "status": "completed",
        "operation": "read_vcf_header_summary",
        "path": str(source),
        "samples": [str(item) for item in headers.samples],
        "info_fields": sorted(str(key) for key in headers.infos),
        "format_fields": sorted(str(key) for key in headers.formats),
        "filters": sorted(str(key) for key in headers.filters),
        "contigs": contigs,
        "provenance": _provenance(),
    }


def vcf_to_npz_archive(
    source: str,
    output: str,
    *,
    fields: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local VCF into NumPy arrays and save them as a new compressed ``.npz`` archive."""

    np = _numpy()
    allel = _allel()
    src = input_path(source)
    target = _fresh_file(output, "output", ".npz")
    selected = _fields(fields)
    allel.vcf_to_npz(str(src), str(target), compressed=True, overwrite=False, fields=selected)
    with np.load(target) as archive:
        arrays = {
            name: [int(item) for item in archive[name].shape] for name in sorted(archive.files)
        }
    return {
        "status": "completed",
        "operation": "vcf_to_npz_archive",
        "source": str(src),
        "output": str(target),
        "arrays": arrays,
        "provenance": _provenance(),
    }


def vcf_to_hdf5_archive(
    source: str,
    output: str,
    *,
    fields: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a local VCF into NumPy arrays and load them into a new HDF5 file."""

    allel = _allel()
    h5py = require("h5py", _ALLEL_EXTRA)
    src = input_path(source)
    target = _fresh_file(output, "output")
    selected = _fields(fields)
    allel.vcf_to_hdf5(str(src), str(target), overwrite=False, fields=selected)
    datasets: dict[str, list[int]] = {}
    with h5py.File(target, "r") as handle:

        def collect(name: str, item: Any) -> None:
            if isinstance(item, h5py.Dataset):
                datasets[name] = [int(size) for size in item.shape]

        handle.visititems(collect)
    return {
        "status": "completed",
        "operation": "vcf_to_hdf5_archive",
        "source": str(src),
        "output": str(target),
        "datasets": {name: datasets[name] for name in sorted(datasets)},
        "provenance": _provenance(),
    }


__all__ = [
    "build_equally_accessible_windows",
    "calculate_between_population_pi",
    "calculate_blockwise_patterson_f3",
    "calculate_blockwise_patterson_fst",
    "calculate_haplotype_diversity",
    "calculate_moving_delta_tajimas_d",
    "calculate_moving_haplotype_diversity",
    "calculate_moving_patterson_f3",
    "calculate_moving_patterson_fst",
    "calculate_pairwise_dxy",
    "calculate_patterson_f2",
    "calculate_patterson_f3",
    "calculate_patterson_fst",
    "calculate_windowed_fixed_difference_density",
    "calculate_windowed_patterson_fst",
    "calculate_xpnsl",
    "count_variants_in_windows",
    "fold_joint_site_frequency_spectrum",
    "locate_population_fixed_differences",
    "locate_population_private_alleles",
    "paint_parental_transmission",
    "paint_shared_haplotype_prefixes",
    "read_vcf_header_summary",
    "scale_folded_joint_site_frequency_spectrum",
    "scale_folded_site_frequency_spectrum",
    "scale_joint_site_frequency_spectrum",
    "scale_site_frequency_spectrum",
    "standardize_score_by_allele_count",
    "vcf_to_hdf5_archive",
    "vcf_to_npz_archive",
    "windowed_values_per_base",
]

_register_public_verbose(globals(), tuple(__all__))
