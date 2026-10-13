"""Genotype-matrix operations on PLINK 1 data through bed-reader.

- ``standardize_genotypes`` wraps bed-reader's ``standardize_f64`` and ``standardize_f32`` kernels.
- ``calculate_genotype_relationship_matrix`` reads a ``.bed`` file with ``open_bed``, standardizes
  it with ``standardize_f64`` and returns ``Z Z' / m`` computed with NumPy.
- ``subset_genotype_matrix`` wraps bed-reader's ``subset_f64_f64`` and ``subset_f32_f32`` kernels.

Genotype matrices have one row per sample and one column per variant. Missing calls are ``None``
in inputs and results. Standardization uses population standard deviations (ddof 0). Packages are
imported only when a function runs; a missing package raises ``BackendUnavailableError``.
"""

from __future__ import annotations

import math
import numbers
from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.variants.genotype_formats import _check_cells, _limit, _plink_paths, _provenance

_BED_EXTRA = "bio-bed-reader"
_NUMPY_EXTRA = "bio-analysis"
_DEFAULT_MAX_CELLS = 5_000_000
_MAX_CELLS_LIMIT = 100_000_000
_PRECISION = {"float32": "f32", "float64": "f64"}


def _precision(value: object) -> tuple[str, str]:
    if not isinstance(value, str) or value not in _PRECISION:
        raise ConfigurationError("dtype must be 'float32' or 'float64'.")
    return value, _PRECISION[value]


def _number(value: object, name: str) -> float:
    if value is None:
        return math.nan
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ConfigurationError(f"{name} values must be finite numbers or None.")
    number = float(value)
    if not math.isfinite(number):
        raise ConfigurationError(f"{name} values must be finite numbers or None.")
    return number


def _rows(values: object, name: str) -> list[list[float]]:
    if not isinstance(values, Sequence) or isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of sample rows.")
    rows: list[list[float]] = []
    width = -1
    for row in values:
        if not isinstance(row, Sequence) or isinstance(row, str) or not row:
            raise ConfigurationError(f"each row of {name} must be a non-empty list.")
        if width < 0:
            width = len(row)
        if len(row) != width:
            raise ConfigurationError(f"all rows of {name} must have {width} variant columns.")
        rows.append([_number(item, name) for item in row])
    return rows


def _indices(values: object, limit: int, name: str) -> list[int]:
    if not isinstance(values, Sequence) or isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of integers.")
    selected: list[int] = []
    for item in values:
        if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < limit:
            raise ConfigurationError(f"{name} must contain integers in [0, {limit}).")
        selected.append(item)
    return selected


def _finite_or_none(value: object) -> float | None:
    number = float(value)  # type: ignore[arg-type]
    return number if math.isfinite(number) else None


def _matrix_rows(values: Any) -> list[list[float | None]]:
    return [[_finite_or_none(item) for item in row] for row in values.tolist()]


def standardize_genotypes(
    genotypes: Sequence[Sequence[float | None]],
    *,
    dtype: str = "float64",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Standardize each variant column to mean 0 and population SD 1 with bed-reader."""

    name, suffix = _precision(dtype)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    rows = _rows(genotypes, "genotypes")
    _check_cells(len(rows), len(rows[0]), cells)
    np = require("numpy", _NUMPY_EXTRA)
    module = require("bed_reader", _BED_EXTRA)
    values = np.asfortranarray(np.array(rows, dtype=name))
    observed = np.isfinite(values).any(axis=0)
    if not observed.all():
        first = int(np.flatnonzero(~observed)[0])
        raise ConfigurationError(f"variant column {first} has no observed genotypes.")
    stats = np.zeros((values.shape[1], 2), dtype=name, order="F")
    kernel = getattr(module, f"standardize_{suffix}")
    # bed-reader writes the standardized values in place and the (mean, SD) pairs into stats.
    kernel(values, False, 1.0, 25.0, True, False, stats, module.get_num_threads())
    return {
        "status": "completed",
        "operation": "standardize_genotypes",
        "dtype": name,
        "sample_count": int(values.shape[0]),
        "variant_count": int(values.shape[1]),
        "standardized": _matrix_rows(values),
        "variant_statistics": [
            {
                "index": index,
                "mean": _finite_or_none(stats[index, 0]),
                "std": _finite_or_none(stats[index, 1]),
            }
            for index in range(values.shape[1])
        ],
        "provenance": _provenance(
            "bed-reader",
            kernel=f"standardize_{suffix}",
            convention="population SD (ddof 0); missing and zero-variance values become 0",
        ),
    }


def calculate_genotype_relationship_matrix(
    path: object,
    *,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a standardized genotype relationship matrix Z Z' / m from a PLINK 1 .bed file."""

    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_path, _bim, _fam = _plink_paths(path)
    np = require("numpy", _NUMPY_EXTRA)
    module = require("bed_reader", _BED_EXTRA)
    with module.open_bed(str(bed_path)) as handle:
        sample_count = int(handle.iid_count)
        variant_count = int(handle.sid_count)
        _check_cells(sample_count, variant_count, cells)
        raw = np.asarray(handle.read(dtype="float64"), dtype=np.float64)
        samples = [
            {"family_id": str(family), "sample_id": str(individual)}
            for family, individual in zip(handle.fid, handle.iid, strict=True)
        ]
    observed = np.isfinite(raw).any(axis=0)
    if not observed.any():
        raise ConfigurationError("no variant has observed genotypes.")
    values = np.asfortranarray(raw[:, observed])
    stats = np.zeros((values.shape[1], 2), dtype=np.float64, order="F")
    module.standardize_f64(values, False, 1.0, 25.0, True, False, stats, module.get_num_threads())
    spread = stats[:, 1]
    polymorphic = int(np.count_nonzero(np.isfinite(spread) & (spread > 0)))
    if polymorphic == 0:
        raise ConfigurationError("no polymorphic variant is available for the matrix.")
    matrix = np.dot(values, values.T) / polymorphic
    return {
        "status": "completed",
        "operation": "calculate_genotype_relationship_matrix",
        "samples": samples,
        "matrix": matrix.tolist(),
        "variant_count": variant_count,
        "observed_variants": int(values.shape[1]),
        "polymorphic_variants": polymorphic,
        "provenance": _provenance(
            "bed-reader",
            kernel="standardize_f64",
            formula="Z Z' divided by the number of polymorphic variants",
            convention="unobserved variants are skipped; missing genotypes contribute 0",
        ),
    }


def subset_genotype_matrix(
    genotypes: Sequence[Sequence[float | None]],
    *,
    sample_indices: Sequence[int],
    variant_indices: Sequence[int],
    dtype: str = "float64",
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select samples (rows) and variants (columns) from an in-memory genotype matrix."""

    name, suffix = _precision(dtype)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    rows = _rows(genotypes, "genotypes")
    samples = _indices(sample_indices, len(rows), "sample_indices")
    variants = _indices(variant_indices, len(rows[0]), "variant_indices")
    _check_cells(len(samples), len(variants), cells)
    np = require("numpy", _NUMPY_EXTRA)
    module = require("bed_reader", _BED_EXTRA)
    # bed-reader's subset kernels take three-dimensional arrays; the last axis has length 1.
    source = np.asfortranarray(np.array(rows, dtype=name)[:, :, None])
    target = np.zeros((len(samples), len(variants), 1), dtype=name, order="F")
    kernel = getattr(module, f"subset_{suffix}_{suffix}")
    kernel(
        source,
        np.asarray(samples, dtype=np.uintp),
        np.asarray(variants, dtype=np.uintp),
        target,
        module.get_num_threads(),
    )
    return {
        "status": "completed",
        "operation": "subset_genotype_matrix",
        "dtype": name,
        "sample_indices": samples,
        "variant_indices": variants,
        "genotypes": _matrix_rows(target[:, :, 0]),
        "provenance": _provenance(
            "bed-reader",
            kernel=f"subset_{suffix}_{suffix}",
            convention="indices select rows and columns in the given order and may repeat",
        ),
    }


__all__ = [
    "calculate_genotype_relationship_matrix",
    "standardize_genotypes",
    "subset_genotype_matrix",
]

_register_public_verbose(globals(), tuple(__all__))
