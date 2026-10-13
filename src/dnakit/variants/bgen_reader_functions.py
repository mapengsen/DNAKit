"""Allele expectations, dosages and metafile creation for BGEN files through bgen-reader.

- ``read_bgen_allele_expectation`` wraps ``bgen_reader.allele_expectation`` for unphased variants.
- ``read_bgen_dosage`` wraps ``bgen_reader.compute_dosage`` on those expectations.
- ``create_bgen_metafile`` wraps ``bgen_reader.create_metafile`` and checks the result by reading
  the variant table back with ``bgen_reader.read_bgen``.

Expectations are the expected number of copies of each allele per sample. Missing samples are
``None``. Phased variants are rejected because bgen-reader defines these values for unphased
genotypes only. Packages are imported only when a function runs.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.variants.genotype_formats import (
    _DEFAULT_MAX_CELLS,
    _DEFAULT_MAX_RECORDS,
    _MAX_CELLS_LIMIT,
    _bgen_indices,
    _bgen_variant_record,
    _check_cells,
    _limit,
    _open_bgen,
    _path,
    _provenance,
)

_BGEN_EXTRA = "bio-bgen-reader"
_NUMPY_EXTRA = "bio-analysis"


def _unphased_indices(handle: Any, variant_indices: object, operation: str) -> list[int]:
    selected = _bgen_indices(variant_indices, int(handle.nvariants))
    if any(bool(handle.phased[index]) for index in selected):
        raise ConfigurationError(f"{operation} supports unphased variants only.")
    return selected


def _allele_total(handle: Any, selected: Sequence[int]) -> int:
    return sum(len(str(handle.allele_ids[index]).split(",")) for index in selected)


def _alt_index(value: object) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError("alt_allele_index must be None or a non-negative integer.")
    return value


def _finite_row(values: Sequence[float]) -> list[float] | None:
    row = [float(item) for item in values]
    return row if all(math.isfinite(item) for item in row) else None


def _finite_scalar(value: object) -> float | None:
    number = float(value)  # type: ignore[arg-type]
    return number if math.isfinite(number) else None


def _metafile_path(value: object, bgen_path: Path, overwrite: object) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError("metafile_path must be a non-empty file path.")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be True or False.")
    target = Path(value).expanduser().resolve()
    if target == bgen_path:
        raise ConfigurationError("metafile_path must differ from the BGEN file.")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    if target.exists() and not overwrite:
        raise ConfigurationError(f"Refusing to overwrite existing file: {target}")
    return target


def read_bgen_allele_expectation(
    path: object,
    *,
    variant_indices: Sequence[int] | None = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_variants: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the expected copies of each allele per sample for unphased BGEN variants."""

    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    handle = _open_bgen(path, max_variants)
    np = require("numpy", _NUMPY_EXTRA)
    with handle:
        selected = _unphased_indices(handle, variant_indices, "read_bgen_allele_expectation")
        sample_count = int(handle.nsamples)
        _check_cells(sample_count, _allele_total(handle, selected), cells)
        # bgen-reader returns samples x variants x alleles for the selected variant indices.
        expectation = np.asarray(handle.allele_expectation(selected), dtype=float)
        samples = [str(item) for item in handle.samples]
        variants = [_bgen_variant_record(handle, index) for index in selected]
    expectations = [
        [_finite_row(expectation[sample, column, :]) for sample in range(sample_count)]
        for column in range(len(selected))
    ]
    return {
        "status": "completed",
        "operation": "read_bgen_allele_expectation",
        "samples": samples,
        "variants": variants,
        "expectations": expectations,
        "provenance": _provenance(
            "bgen-reader",
            layout="expectations are ordered as variant, sample, allele",
            convention="unphased genotypes only; a missing sample is None",
        ),
    }


def read_bgen_dosage(
    path: object,
    *,
    variant_indices: Sequence[int] | None = None,
    alt_allele_index: int | None = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_variants: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the dosage of one alternative allele per sample for unphased BGEN variants."""

    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    alt = _alt_index(alt_allele_index)
    handle = _open_bgen(path, max_variants)
    np = require("numpy", _NUMPY_EXTRA)
    bgen = require("bgen_reader", _BGEN_EXTRA)
    with handle:
        selected = _unphased_indices(handle, variant_indices, "read_bgen_dosage")
        sample_count = int(handle.nsamples)
        _check_cells(sample_count, _allele_total(handle, selected), cells)
        expectation = np.asarray(handle.allele_expectation(selected), dtype=float)
        samples = [str(item) for item in handle.samples]
        variants: list[dict[str, Any]] = []
        dosages: list[list[float | None]] = []
        for column, index in enumerate(selected):
            record = _bgen_variant_record(handle, index)
            alleles = record["alleles"]
            if alt is not None and alt >= len(alleles):
                raise ConfigurationError(
                    f"alt_allele_index {alt} is out of range for variant {record['id']}, "
                    f"which has {len(alleles)} alleles."
                )
            # bgen-reader's compute_dosage with alt=None uses the last allele column.
            used = len(alleles) - 1 if alt is None else alt
            values = bgen.compute_dosage(expectation[:, column, :], alt=alt)
            record["alt_allele_index"] = used
            record["alt_allele"] = alleles[used]
            variants.append(record)
            dosages.append([_finite_scalar(value) for value in np.asarray(values, dtype=float)])
    return {
        "status": "completed",
        "operation": "read_bgen_dosage",
        "samples": samples,
        "variants": variants,
        "dosages": dosages,
        "provenance": _provenance(
            "bgen-reader",
            statistic="expected copies of the alternative allele",
            convention="alt_allele_index None uses the last allele, bgen-reader's default",
        ),
    }


def create_bgen_metafile(
    path: object,
    metafile_path: object,
    *,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a variant metadata file that speeds up later read_bgen calls on the same BGEN file."""

    bgen_path = input_path(_path(path, ".bgen"))
    target = _metafile_path(metafile_path, bgen_path, overwrite)
    module = require("bgen_reader", _BGEN_EXTRA)
    if target.exists():
        target.unlink()  # bgen-reader refuses to replace an existing metafile
    module.create_metafile(str(bgen_path), str(target), verbose=False)
    frame = module.read_bgen(str(bgen_path), metafile_filepath=str(target), verbose=False)
    return {
        "status": "completed",
        "operation": "create_bgen_metafile",
        "path": str(bgen_path),
        "metafile": str(target),
        "variants": len(frame["variants"]),
        "bytes": target.stat().st_size,
        "provenance": _provenance(
            "bgen-reader",
            note="metafile is checked by reading the variant table with read_bgen",
        ),
    }


__all__ = [
    "create_bgen_metafile",
    "read_bgen_allele_expectation",
    "read_bgen_dosage",
]

_register_public_verbose(globals(), tuple(__all__))
