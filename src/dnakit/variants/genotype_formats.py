"""Optional genotype-file adapters for PLINK 1 and BGEN files.

The functions in this module are explicit adapters over optional packages:

- ``bed-reader`` reads, writes and summarizes PLINK 1 ``.bed`` genotype matrices.
- ``pandas-plink`` reads PLINK 1 metadata and genotype matrices, and writes ``.bed`` files.
- ``bgen-reader`` reads BGEN variant metadata, probabilities and allele frequencies.

Genotype matrices have one row per sample and one column per variant. A PLINK dosage counts
one allele: ``count_allele="allele_1"`` counts the first (A1) allele and ``"allele_2"`` counts
the second. Missing calls are returned as ``None``. Every function bounds the number of records
and matrix cells it reads, and raises ``ConfigurationError`` instead of loading more.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_BED_EXTRA = "bio-bed-reader"
_PANDAS_PLINK_EXTRA = "bio-pandas-plink"
_BGEN_EXTRA = "bio-bgen-reader"
_COUNT_ALLELES = ("allele_1", "allele_2")
_PLINK_SUFFIXES = (".bed", ".bim", ".fam")
_DEFAULT_MAX_RECORDS = 100_000
_MAX_RECORDS_LIMIT = 1_000_000
_DEFAULT_MAX_CELLS = 5_000_000
_MAX_CELLS_LIMIT = 100_000_000


def _limit(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
    return value


def _count_allele(value: object) -> str:
    if value not in _COUNT_ALLELES:
        raise ConfigurationError("count_allele must be 'allele_1' or 'allele_2'.")
    return str(value)


def _text(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or "\x00" in value
        or any(character.isspace() for character in value)
    ):
        raise ConfigurationError(f"{name} must contain non-empty strings without whitespace.")
    return value


def _strings(values: object, name: str) -> list[str]:
    if not isinstance(values, Sequence) or isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of strings.")
    return [_text(item, name) for item in values]


def _unique(values: list[str], name: str) -> list[str]:
    if len(set(values)) != len(values):
        raise ConfigurationError(f"{name} must be unique.")
    return values


def _same_length(values: list[str], expected: int, name: str) -> None:
    if len(values) != expected:
        raise ConfigurationError(f"{name} must have {expected} entries.")


def _positions(values: object, expected: int) -> list[int]:
    if not isinstance(values, Sequence) or isinstance(values, str) or len(values) != expected:
        raise ConfigurationError(f"positions must be a list of {expected} integers.")
    result = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ConfigurationError("positions must contain non-negative integers.")
        result.append(value)
    return result


def _path(value: object, suffix: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip() or "\x00" in str(value):
        raise ConfigurationError("path must be a non-empty file path.")
    path = Path(value).expanduser().resolve()
    if path.suffix != suffix:
        raise ConfigurationError(f"path must end with {suffix}.")
    return path


def _new_plink_path(value: object, overwrite: bool) -> Path:
    """Return the .bed path for a PLINK file set; existing files need overwrite=True."""

    bed_path = _path(value, ".bed")
    if not bed_path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {bed_path.parent}")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be True or False.")
    for suffix in _PLINK_SUFFIXES:
        target = bed_path.with_suffix(suffix)
        if target.exists() and not overwrite:
            raise ConfigurationError(f"Refusing to overwrite existing file: {target}")
    return bed_path


def _plink_paths(value: object) -> tuple[Path, Path, Path]:
    bed = input_path(_path(value, ".bed"))
    return bed, input_path(bed.with_suffix(".bim")), input_path(bed.with_suffix(".fam"))


def _bounded_lines(path: Path, limit: int) -> int:
    """Count non-blank lines while streaming; stop as soon as the limit is exceeded."""

    count = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                count += 1
                if count > limit:
                    raise ConfigurationError(f"{path.name} has more than {limit} records.")
    return count


def _check_cells(rows: int, columns: int, limit: int) -> None:
    if rows * columns > limit:
        raise ConfigurationError(
            f"Matrix has {rows * columns} cells, above max_cells={limit}; "
            "select fewer samples or variants."
        )


def _select(available: Iterable[object], requested: object, name: str) -> list[int]:
    names = [str(item) for item in available]
    if requested is None:
        return list(range(len(names)))
    if not isinstance(requested, Sequence) or isinstance(requested, str) or not requested:
        raise ConfigurationError(f"{name} must be a non-empty list of IDs or None.")
    first: dict[str, int] = {}
    for index, item in enumerate(names):
        first.setdefault(item, index)
    selected = []
    for item in requested:
        identifier = _text(item, name)
        if identifier not in first:
            raise ConfigurationError(f"{name} contains an ID that is not in the file: {item!r}.")
        selected.append(first[identifier])
    return selected


def _dosage(value: object) -> int | None:
    number = float(value)
    if not math.isfinite(number):
        return None
    if number not in (0.0, 1.0, 2.0):
        raise ConfigurationError(f"Genotype dosage must be 0, 1, 2 or missing; got {number}.")
    return int(number)


def _finite_float(value: object) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None


def _matrix_values(
    genotypes: object, rows: int, columns: int, limit: int
) -> list[list[int | None]]:
    if not isinstance(genotypes, Sequence) or isinstance(genotypes, str):
        raise ConfigurationError("genotypes must be a list of sample rows.")
    if len(genotypes) != rows:
        raise ConfigurationError(f"genotypes must have {rows} sample rows.")
    _check_cells(rows, columns, limit)
    matrix: list[list[int | None]] = []
    for row in genotypes:
        if not isinstance(row, Sequence) or isinstance(row, str) or len(row) != columns:
            raise ConfigurationError(f"each genotype row must have {columns} variant columns.")
        values: list[int | None] = []
        for value in row:
            if value is None:
                values.append(None)
            elif isinstance(value, bool) or value not in (0, 1, 2):
                raise ConfigurationError("genotype values must be 0, 1, 2 or None.")
            else:
                values.append(int(value))
        matrix.append(values)
    return matrix


def _provenance(package: str, **notes: str) -> dict[str, str]:
    return {"package": f"{package} {metadata.version(package)}", **notes}


def _read_bed_matrix(
    path: object,
    *,
    samples: object,
    variants: object,
    count_allele: str,
    max_records: int,
    max_cells: int,
) -> tuple[Any, list[int], list[int], Any]:
    records = _limit(max_records, "max_records", _MAX_RECORDS_LIMIT)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_path, bim_path, fam_path = _plink_paths(path)
    sample_count = _bounded_lines(fam_path, records)
    variant_count = _bounded_lines(bim_path, records)
    module = require("bed_reader", _BED_EXTRA)
    handle = module.open_bed(str(bed_path), count_A1=count_allele == "allele_1")
    if len(handle.iid) != sample_count or len(handle.sid) != variant_count:
        raise ConfigurationError("PLINK .fam and .bim record counts do not match the .bed file.")
    sample_index = _select(handle.iid, samples, "samples")
    variant_index = _select(handle.sid, variants, "variants")
    _check_cells(len(sample_index), len(variant_index), cells)
    values = handle.read(index=(sample_index, variant_index), dtype="float32")
    return handle, sample_index, variant_index, values


def read_plink_bed(
    path: object,
    *,
    samples: Sequence[str] | None = None,
    variants: Sequence[str] | None = None,
    count_allele: str = "allele_1",
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a PLINK 1 .bed genotype matrix, optionally selecting samples and variants by ID."""

    count = _count_allele(count_allele)
    handle, sample_index, variant_index, values = _read_bed_matrix(
        path,
        samples=samples,
        variants=variants,
        count_allele=count,
        max_records=max_records,
        max_cells=max_cells,
    )
    return {
        "status": "completed",
        "operation": "read_plink_bed",
        "samples": [str(handle.iid[index]) for index in sample_index],
        "variants": [
            {
                "id": str(handle.sid[index]),
                "chromosome": str(handle.chromosome[index]),
                "position": int(handle.bp_position[index]),
                "allele_1": str(handle.allele_1[index]),
                "allele_2": str(handle.allele_2[index]),
            }
            for index in variant_index
        ],
        "genotypes": [[_dosage(value) for value in row] for row in values.tolist()],
        "count_allele": count,
        "provenance": _provenance("bed-reader", convention="dosage counts the selected allele"),
    }


def write_plink_bed(
    path: object,
    genotypes: Sequence[Sequence[int | None]],
    *,
    sample_ids: Sequence[str],
    variant_ids: Sequence[str],
    chromosomes: Sequence[str],
    positions: Sequence[int],
    allele_1: Sequence[str],
    allele_2: Sequence[str],
    count_allele: str = "allele_1",
    overwrite: bool = False,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a genotype matrix as PLINK 1 .bed, .bim and .fam files with bed-reader."""

    count = _count_allele(count_allele)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_path = _new_plink_path(path, overwrite)
    samples = _unique(_strings(sample_ids, "sample_ids"), "sample_ids")
    variants = _unique(_strings(variant_ids, "variant_ids"), "variant_ids")
    _same_length(_strings(chromosomes, "chromosomes"), len(variants), "chromosomes")
    pos = _positions(positions, len(variants))
    first = _strings(allele_1, "allele_1")
    second = _strings(allele_2, "allele_2")
    _same_length(first, len(variants), "allele_1")
    _same_length(second, len(variants), "allele_2")
    matrix = _matrix_values(genotypes, len(samples), len(variants), cells)

    np = require("numpy", "bio-analysis")
    values = np.array(
        [[math.nan if value is None else value for value in row] for row in matrix],
        dtype=np.float32,
    )
    module = require("bed_reader", _BED_EXTRA)
    module.to_bed(
        str(bed_path),
        values,
        properties={
            "iid": samples,
            "sid": variants,
            "chromosome": _strings(chromosomes, "chromosomes"),
            "bp_position": pos,
            "allele_1": first,
            "allele_2": second,
        },
        count_A1=count == "allele_1",
    )
    return {
        "status": "completed",
        "operation": "write_plink_bed",
        "path": str(bed_path),
        "files": [str(bed_path.with_suffix(suffix)) for suffix in _PLINK_SUFFIXES],
        "samples": len(samples),
        "variants": len(variants),
        "count_allele": count,
        "provenance": _provenance("bed-reader", convention="dosage counts the selected allele"),
    }


def summarize_plink_bed(
    path: object,
    *,
    count_allele: str = "allele_1",
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count called and missing genotypes and the allele frequency of each PLINK variant."""

    count = _count_allele(count_allele)
    handle, sample_index, variant_index, values = _read_bed_matrix(
        path,
        samples=None,
        variants=None,
        count_allele=count,
        max_records=max_records,
        max_cells=max_cells,
    )
    rows = [[_dosage(value) for value in row] for row in values.tolist()]
    summaries = []
    missing_total = 0
    for column, index in enumerate(variant_index):
        called = [row[column] for row in rows if row[column] is not None]
        missing = len(rows) - len(called)
        missing_total += missing
        frequency = sum(called) / (2 * len(called)) if called else None
        summaries.append(
            {
                "id": str(handle.sid[index]),
                "called": len(called),
                "missing": missing,
                "allele_frequency": frequency,
            }
        )
    return {
        "status": "completed",
        "operation": "summarize_plink_bed",
        "samples": len(sample_index),
        "variants": summaries,
        "missing_genotypes": missing_total,
        "count_allele": count,
        "provenance": _provenance(
            "bed-reader",
            statistic="allele frequency = counted allele copies / (2 x called genotypes)",
        ),
    }


def read_plink_metadata(
    path: object,
    *,
    max_records: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read PLINK .fam sample rows and .bim variant rows without reading genotypes."""

    records = _limit(max_records, "max_records", _MAX_RECORDS_LIMIT)
    bed_path, bim_path, fam_path = _plink_paths(path)
    _bounded_lines(fam_path, records)
    _bounded_lines(bim_path, records)
    module = require("pandas_plink", _PANDAS_PLINK_EXTRA)
    bim, fam, _bed = module.read_plink(str(bed_path.with_suffix("")), verbose=False)
    return {
        "status": "completed",
        "operation": "read_plink_metadata",
        "samples": [
            {
                "family_id": str(row["fid"]),
                "sample_id": str(row["iid"]),
                "father_id": str(row["father"]),
                "mother_id": str(row["mother"]),
                "sex": str(row["gender"]),
                "phenotype": str(row["trait"]),
            }
            for row in fam.to_dict("records")
        ],
        "variants": [
            {
                "chromosome": str(row["chrom"]),
                "id": str(row["snp"]),
                "genetic_position": _finite_float(row["cm"]),
                "position": int(row["pos"]),
                "allele_1": str(row["a0"]),
                "allele_2": str(row["a1"]),
            }
            for row in bim.to_dict("records")
        ],
        "provenance": _provenance("pandas-plink", source="read_plink .fam and .bim tables"),
    }


def read_plink_genotypes(
    path: object,
    *,
    count_allele: str = "allele_1",
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a PLINK 1 .bed genotype matrix through pandas-plink's read_plink1_bin."""

    count = _count_allele(count_allele)
    records = _limit(max_records, "max_records", _MAX_RECORDS_LIMIT)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_path, bim_path, fam_path = _plink_paths(path)
    sample_count = _bounded_lines(fam_path, records)
    variant_count = _bounded_lines(bim_path, records)
    _check_cells(sample_count, variant_count, cells)
    module = require("pandas_plink", _PANDAS_PLINK_EXTRA)
    # pandas-plink's ref names the allele whose dosage is returned.
    ref = "a0" if count == "allele_1" else "a1"
    array = module.read_plink1_bin(str(bed_path), verbose=False, ref=ref)
    coords = array.coords
    variant_rows = zip(
        coords["snp"].values,
        coords["chrom"].values,
        coords["pos"].values,
        coords["a0"].values,
        coords["a1"].values,
        strict=True,
    )
    return {
        "status": "completed",
        "operation": "read_plink_genotypes",
        "samples": [str(item) for item in coords["sample"].values],
        "variants": [
            {
                "id": str(snp),
                "chromosome": str(chrom),
                "position": int(pos),
                "allele_1": str(first),
                "allele_2": str(second),
            }
            for snp, chrom, pos, first, second in variant_rows
        ],
        "genotypes": [[_dosage(value) for value in row] for row in array.values.tolist()],
        "count_allele": count,
        "provenance": _provenance("pandas-plink", ref=ref),
    }


def write_plink_genotypes(
    path: object,
    genotypes: Sequence[Sequence[int | None]],
    *,
    sample_ids: Sequence[str],
    variant_ids: Sequence[str],
    chromosomes: Sequence[str],
    positions: Sequence[int],
    allele_1: Sequence[str],
    allele_2: Sequence[str],
    count_allele: str = "allele_1",
    overwrite: bool = False,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a genotype matrix as PLINK 1 files through pandas-plink's write_plink1_bin."""

    count = _count_allele(count_allele)
    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    bed_path = _new_plink_path(path, overwrite)
    samples = _unique(_strings(sample_ids, "sample_ids"), "sample_ids")
    variants = _unique(_strings(variant_ids, "variant_ids"), "variant_ids")
    chrom = _strings(chromosomes, "chromosomes")
    _same_length(chrom, len(variants), "chromosomes")
    pos = _positions(positions, len(variants))
    first = _strings(allele_1, "allele_1")
    second = _strings(allele_2, "allele_2")
    _same_length(first, len(variants), "allele_1")
    _same_length(second, len(variants), "allele_2")
    matrix = _matrix_values(genotypes, len(samples), len(variants), cells)

    np = require("numpy", "bio-analysis")
    xarray = require("xarray", _PANDAS_PLINK_EXTRA)
    module = require("pandas_plink", _PANDAS_PLINK_EXTRA)
    # pandas-plink's writer stores counts of the second allele, so A1 counts are flipped.
    values = np.array(
        [
            [
                math.nan if value is None else (2 - value if count == "allele_1" else value)
                for value in row
            ]
            for row in matrix
        ],
        dtype=np.float32,
    )
    array = xarray.DataArray(
        values,
        dims=["sample", "variant"],
        coords={
            "sample": samples,
            "fid": ("sample", samples),
            "variant": variants,
            "snp": ("variant", variants),
            "chrom": ("variant", chrom),
            "pos": ("variant", pos),
            "a0": ("variant", first),
            "a1": ("variant", second),
        },
    )
    module.write_plink1_bin(array, str(bed_path), verbose=False)
    return {
        "status": "completed",
        "operation": "write_plink_genotypes",
        "path": str(bed_path),
        "files": [str(bed_path.with_suffix(suffix)) for suffix in _PLINK_SUFFIXES],
        "samples": len(samples),
        "variants": len(variants),
        "count_allele": count,
        "provenance": _provenance("pandas-plink", convention="writer stores second-allele counts"),
    }


def _open_bgen(path: object, max_variants: int) -> Any:
    """Open a BGEN file after checking its variant count against the limit."""

    limit = _limit(max_variants, "max_variants", _MAX_RECORDS_LIMIT)
    bgen_path = input_path(_path(path, ".bgen"))
    module = require("bgen_reader", _BGEN_EXTRA)
    handle = module.open_bgen(str(bgen_path), verbose=False)
    if int(handle.nvariants) > limit:
        handle.close()
        raise ConfigurationError(
            f"BGEN file has {handle.nvariants} variants, above max_variants={limit}."
        )
    return handle


def _bgen_indices(variant_indices: object, total: int) -> list[int]:
    if variant_indices is None:
        return list(range(total))
    if not isinstance(variant_indices, Sequence) or isinstance(variant_indices, str):
        raise ConfigurationError("variant_indices must be a list of integers or None.")
    if not variant_indices:
        raise ConfigurationError("variant_indices must not be empty.")
    selected = []
    for item in variant_indices:
        if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < total:
            raise ConfigurationError(f"variant_indices must be integers in [0, {total}).")
        selected.append(item)
    return selected


def _bgen_variant_record(handle: Any, index: int) -> dict[str, Any]:
    return {
        "id": str(handle.ids[index]),
        "rsid": str(handle.rsids[index]),
        "chromosome": str(handle.chromosomes[index]),
        "position": int(handle.positions[index]),
        "alleles": str(handle.allele_ids[index]).split(","),
    }


def read_bgen_variants(
    path: object, *, max_variants: int = _DEFAULT_MAX_RECORDS, verbose: bool = True
) -> dict[str, Any]:
    """Read BGEN variant metadata with bgen-reader, without reading probabilities."""

    handle = _open_bgen(path, max_variants)
    with handle:
        variants = []
        for index in range(int(handle.nvariants)):
            record = _bgen_variant_record(handle, index)
            record["phased"] = bool(handle.phased[index])
            record["genotype_combinations"] = int(handle.ncombinations[index])
            variants.append(record)
        sample_count = int(handle.nsamples)
    return {
        "status": "completed",
        "operation": "read_bgen_variants",
        "sample_count": sample_count,
        "variants": variants,
        "provenance": _provenance("bgen-reader", note="may create a .metadata2.mmm cache file"),
    }


def read_bgen_probabilities(
    path: object,
    *,
    variant_indices: Sequence[int] | None = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_variants: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read BGEN genotype probabilities as variant, sample and genotype-combination lists."""

    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    handle = _open_bgen(path, max_variants)
    with handle:
        combinations = [int(item) for item in handle.ncombinations]
        selected = _bgen_indices(variant_indices, len(combinations))
        samples = [str(item) for item in handle.samples]
        total_values = len(samples) * sum(combinations[index] for index in selected)
        if total_values > cells:
            raise ConfigurationError(
                f"Probability block has {total_values} values, above max_cells={cells}."
            )
        probabilities, missing = handle.read(selected, dtype="float64", return_missings=True)
        variants = [_bgen_variant_record(handle, index) for index in selected]
    values = []
    for column in range(len(selected)):
        block = []
        for sample in range(len(samples)):
            if missing[sample, column]:
                block.append(None)
            else:
                block.append([float(item) for item in probabilities[sample, column]])
        values.append(block)
    return {
        "status": "completed",
        "operation": "read_bgen_probabilities",
        "samples": samples,
        "variants": variants,
        "probabilities": values,
        "provenance": _provenance(
            "bgen-reader",
            layout="probabilities are ordered as variant, sample, genotype combination",
        ),
    }


def calculate_bgen_allele_frequency(
    path: object,
    *,
    variant_indices: Sequence[int] | None = None,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_variants: int = _DEFAULT_MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute allele frequencies of unphased BGEN variants from expected allele copies."""

    cells = _limit(max_cells, "max_cells", _MAX_CELLS_LIMIT)
    handle = _open_bgen(path, max_variants)
    np = require("numpy", "bio-analysis")
    with handle:
        selected = _bgen_indices(variant_indices, int(handle.nvariants))
        if any(bool(handle.phased[index]) for index in selected):
            raise ConfigurationError(
                "calculate_bgen_allele_frequency supports unphased variants only."
            )
        nsamples = int(handle.nsamples)
        allele_total = sum(len(str(handle.allele_ids[index]).split(",")) for index in selected)
        if nsamples * allele_total > cells:
            raise ConfigurationError(f"Allele block is above max_cells={cells}.")
        expectation = np.asarray(handle.allele_expectation(selected), dtype=float)
        missing = handle.read(selected, return_probabilities=False, return_missings=True)
        variants = []
        for column, index in enumerate(selected):
            present = ~np.asarray(missing[:, column], dtype=bool)
            copies = np.nansum(expectation[present, column, :], axis=0)
            total = float(copies.sum())
            record = _bgen_variant_record(handle, index)
            record["called_samples"] = int(present.sum())
            record["allele_copies"] = [float(item) for item in copies]
            record["allele_frequencies"] = [
                float(item / total) if total > 0 else None for item in copies
            ]
            variants.append(record)
    return {
        "status": "completed",
        "operation": "calculate_bgen_allele_frequency",
        "variants": variants,
        "provenance": _provenance(
            "bgen-reader",
            statistic="expected allele copies from allele_expectation / total copies",
            note="bgen_reader.allele_frequency returns copy sums divided by allele count",
        ),
    }


__all__ = [
    "calculate_bgen_allele_frequency",
    "read_bgen_probabilities",
    "read_bgen_variants",
    "read_plink_bed",
    "read_plink_genotypes",
    "read_plink_metadata",
    "summarize_plink_bed",
    "write_plink_bed",
    "write_plink_genotypes",
]

_register_public_verbose(globals(), tuple(__all__))
