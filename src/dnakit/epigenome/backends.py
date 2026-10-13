"""Optional backend adapters for contact maps, pairs statistics, coverage tracks and peaks.

The functions in this module call upstream Python APIs of optional packages:

- ``build_contact_matrix()``, ``fetch_contact_matrix()`` and
  ``balance_contact_matrix()`` use cooler to bin a .pairs contact list into a
  .cool matrix, read blocks of it and store iterative-correction weights.
- ``calculate_coverage()``, ``calculate_expected_contacts()`` and
  ``calculate_insulation_score()`` use cooltools on a cooler file.
- ``read_pairs_header()``, ``summarize_pairs()`` and ``deduplicate_pairs()`` use
  pairtools' header, statistics and deduplication functions.
- ``build_fragment_pileup()``, ``build_extended_tag_pileup()`` and
  ``call_fragment_peaks()`` use MACS3's pileup and naive peak-calling functions.
- ``filter_allc_sites()``, ``calculate_non_conversion_rate()`` and
  ``calculate_site_binomial_pvalues()`` use methylpy's allc functions for
  bisulfite cytosine counts.
- ``summarize_bigwig_bins()``, ``count_bam_reads_per_bin()`` and
  ``estimate_fragment_lengths()`` use deeptools.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import contextlib
import io
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_COOLTOOLS_EXTRA = "bio-cooltools"
_PAIRTOOLS_EXTRA = "bio-pairtools"
_MACS3_EXTRA = "bio-macs3"
_METHYLPY_EXTRA = "bio-methylpy"
_SHELL_SAFE_PATH = re.compile(r"^[A-Za-z0-9_./@+=,:-]+$")
_SORT_BUFFER = re.compile(r"^[0-9]+[KMGTkmgt%]?$")
_CONTROL_REGION = re.compile(r"^[^:\s]+(?::[0-9]+(?:-[0-9]+)?)?$")
_DEEPTOOLS_EXTRA = "bio-deeptools"
_MAX_FILES = 32
_MAX_CHROMOSOMES = 10_000
_MAX_CHROM_LENGTH = 2_000_000_000
_MAX_BIN_SIZE = 10_000_000
_DEFAULT_MAX_BINS = 100_000
_MAX_BINS_LIMIT = 1_000_000
_MAX_PROCESSORS = 16
_MAX_MAPPING_QUALITY = 255
_MAX_FRAGMENTS = 10_000_000
_DEFAULT_MAX_FRAGMENTS = 1_000_000
_MAX_INT32 = 2_147_483_647
_MAX_ROWS = 10_000_000


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


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
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


def _paths(value: object, name: str) -> list[Path]:
    if (
        isinstance(value, (str, Path))
        or not isinstance(value, Sequence)
        or not 1 <= len(value) <= _MAX_FILES
    ):
        raise ConfigurationError(f"{name} must be a list of 1 to {_MAX_FILES} file paths.")
    paths = [_path(item, name) for item in value]
    if len(set(paths)) != len(paths):
        raise ConfigurationError(f"{name} must not list the same file twice.")
    return paths


def _chromsizes(value: object) -> dict[str, int]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CHROMOSOMES:
        raise ConfigurationError(
            f"chromsizes must map 1 to {_MAX_CHROMOSOMES} chromosome names to lengths."
        )
    sizes: dict[str, int] = {}
    for name, length in value.items():
        _text(name, "chromosome name")
        sizes[name] = _integer(length, f"length of {name}", 1, _MAX_CHROM_LENGTH)
    return sizes


def _check_bin_count(chromosomes: Sequence[tuple[str, int]], bin_size: int, max_bins: int) -> int:
    count = sum(-(-length // bin_size) for _, length in chromosomes)
    if count > max_bins:
        raise ConfigurationError(
            f"The bin layout has {count} bins, above max_bins={max_bins}; "
            "increase bin_size or max_bins."
        )
    return count


def _plain(value: Any) -> Any:
    """Convert numpy scalars to Python values and non-finite floats to None."""

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
    return {"package": f"{distribution} {installed}", "statistic": statistic, **notes}


def _shared_bigwig_chromosomes(paths: Sequence[Path]) -> list[tuple[str, int]]:
    """Chromosomes with the same name and length in every bigWig file."""

    pybigwig = require("pyBigWig", "bio-pybigwig")
    shared: dict[str, int] | None = None
    for path in paths:
        try:
            handle = pybigwig.open(str(path))
        except RuntimeError as exc:
            raise ConfigurationError(f"pyBigWig could not open {path}: {exc}") from exc
        if handle is None or not handle.isBigWig():
            raise ConfigurationError(f"Input is not a bigWig file: {path}")
        chromosomes = {name: int(length) for name, length in handle.chroms().items()}
        handle.close()
        shared = chromosomes if shared is None else _intersect(shared, chromosomes)
    if not shared:
        raise ConfigurationError(
            "The bigWig files share no chromosome with the same name and length."
        )
    return sorted(shared.items())


def _shared_bam_chromosomes(paths: Sequence[Path]) -> list[tuple[str, int]]:
    """Chromosomes with the same name and length in every indexed BAM file."""

    pysam = require("pysam", "bio-reads")
    shared: dict[str, int] | None = None
    for path in paths:
        try:
            handle = pysam.AlignmentFile(str(path), "rb")
        except (OSError, ValueError) as exc:
            raise ConfigurationError(f"Input is not a readable BAM file: {path}") from exc
        with handle:
            try:
                indexed = handle.check_index()
            except (OSError, ValueError):
                indexed = False
            if not indexed:
                raise ConfigurationError(f"BAM file needs a .bai index: {path}")
            chromosomes = {
                name: int(length)
                for name, length in zip(handle.references, handle.lengths, strict=True)
            }
        shared = chromosomes if shared is None else _intersect(shared, chromosomes)
    if not shared:
        raise ConfigurationError("The BAM files share no chromosome with the same name and length.")
    return sorted(shared.items())


def _intersect(left: Mapping[str, int], right: Mapping[str, int]) -> dict[str, int]:
    return {name: length for name, length in left.items() if right.get(name) == length}


def _chromosome_rows(chromosomes: Sequence[tuple[str, int]]) -> list[dict[str, Any]]:
    return [{"name": name, "length": length} for name, length in chromosomes]


def _read_bins(raw: Path, columns: int, convert: Any) -> list[dict[str, Any]]:
    """Parse deeptools' raw per-bin table: chrom, start, end, then one value per file."""

    rows: list[dict[str, Any]] = []
    for line in raw.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) != 3 + columns:
            raise BackendExecutionError(
                "deeptools wrote a bin row with an unexpected column count."
            )
        rows.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "values": [convert(value) for value in fields[3:]],
            }
        )
    return rows


def _open_cooler(cooler: Any, path: Path) -> Any:
    try:
        return cooler.Cooler(str(path))
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {path} ({exc})") from exc


def _weighted(clr: Any, balanced: bool) -> None:
    if balanced and "weight" not in clr.bins()[:1].columns:
        raise ConfigurationError(
            "This cooler has no balancing weights; run balance_contact_matrix() first."
        )


def _matrix_rows(values: Any) -> list[list[Any]]:
    return [[_plain(float(value)) for value in row] for row in values]


def build_contact_matrix(
    pairs_path: str | Path,
    cool_path: str | Path,
    *,
    chromsizes: Mapping[str, int],
    binsize: int,
    max_bins: int = _DEFAULT_MAX_BINS,
    chunksize: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bin a .pairs contact list into a new cooler (.cool) matrix, using cooler."""

    source = _path(pairs_path, "pairs_path")
    output = _new_path(cool_path, "cool_path")
    sizes = _chromsizes(chromsizes)
    size = _integer(binsize, "binsize", 1, _MAX_BIN_SIZE)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    rows_per_chunk = _integer(chunksize, "chunksize", 1, _MAX_ROWS)
    chromosomes = list(sizes.items())
    expected = _check_bin_count(chromosomes, size, limit)

    cooler = require("cooler", _COOLER_EXTRA)
    pandas = require("pandas", _COOLER_EXTRA)
    bins = cooler.binnify(pandas.Series(sizes, dtype="int64"), size)
    sanitize = cooler.create.sanitize_records(
        bins,
        schema="pairs",
        decode_chroms=True,
        is_one_based=True,
        tril_action="reflect",
        sort=True,
        validate=True,
    )
    aggregate = cooler.create.aggregate_records(agg={}, count=True, sort=False)

    def pipeline(frame: Any) -> Any:
        return aggregate(sanitize(frame))

    try:
        with pandas.read_csv(
            source,
            sep="\t",
            comment="#",
            header=None,
            usecols=[1, 2, 3, 4],
            names=["chrom1", "pos1", "chrom2", "pos2"],
            dtype={"chrom1": str, "pos1": "int64", "chrom2": str, "pos2": "int64"},
            iterator=True,
            chunksize=rows_per_chunk,
        ) as reader:
            cooler.create_cooler(
                str(output),
                bins,
                map(pipeline, reader),
                columns=["count"],
                dtypes={},
                ordered=False,
                symmetric_upper=True,
                boundscheck=False,
                triucheck=False,
                dupcheck=False,
                ensure_sorted=False,
                temp_dir=str(output.parent),
                delete_temp=True,
            )
        info = cooler.Cooler(str(output)).info
    except (ValueError, IndexError, KeyError, OSError) as exc:
        output.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler rejected the pairs input: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_contact_matrix",
        "pairs": str(source),
        "cool_path": str(output),
        "binsize": size,
        "chromosomes": _chromosome_rows(chromosomes),
        "bins": expected,
        "nonzero_pixels": int(info["nnz"]),
        "total_contacts": int(info["sum"]),
        "provenance": _provenance(
            "cooler",
            "contact counts per bin pair; each pair is listed once, lower triangle reflected",
            coordinate_system="pair positions are 1-based; bins are 0-based half-open",
            input_columns="pairs columns 2 to 5: chrom1, pos1, chrom2, pos2",
        ),
    }


def fetch_contact_matrix(
    cool_path: str | Path,
    region: str,
    *,
    balance: bool = False,
    max_region_bins: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the dense contact matrix for one genomic region from a cooler file."""

    path = _path(cool_path, "cool_path")
    text = _text(region, "region")
    use_weights = _flag(balance, "balance")
    limit = _integer(max_region_bins, "max_region_bins", 1, 1_000)
    cooler = require("cooler", _COOLER_EXTRA)
    clr = _open_cooler(cooler, path)
    _weighted(clr, use_weights)
    try:
        lo, hi = clr.extent(text)
    except (KeyError, ValueError) as exc:
        raise ConfigurationError(f"Region is not valid for this cooler: {text} ({exc})") from exc
    if hi - lo > limit:
        raise ConfigurationError(
            f"Region spans {hi - lo} bins, above max_region_bins={limit}; choose a smaller region."
        )
    matrix = clr.matrix(balance=use_weights).fetch(text)
    bins = clr.bins()[lo:hi]
    return {
        "status": "completed",
        "operation": "fetch_contact_matrix",
        "cool_path": str(path),
        "region": text,
        "binsize": int(clr.binsize),
        "balanced": use_weights,
        "bins": [
            {"chrom": str(row.chrom), "start": int(row.start), "end": int(row.end)}
            for row in bins.itertuples(index=False)
        ],
        "matrix": _matrix_rows(matrix),
        "provenance": _provenance(
            "cooler",
            "balanced matrix values use the stored weights; null marks bins without a weight"
            if use_weights
            else "raw contact counts",
            coordinate_system="0-based half-open bins",
        ),
    }


def balance_contact_matrix(
    cool_path: str | Path,
    *,
    cis_only: bool = False,
    ignore_diags: int = 2,
    min_nnz: int = 10,
    min_count: int = 0,
    mad_max: int = 5,
    tol: float = 1e-5,
    max_iters: int = 200,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute iterative-correction weights with cooler and store them as a weight column."""

    path = _path(cool_path, "cool_path")
    options = {
        "cis_only": _flag(cis_only, "cis_only"),
        "ignore_diags": _integer(ignore_diags, "ignore_diags", 0, 100),
        "min_nnz": _integer(min_nnz, "min_nnz", 0, 1_000_000),
        "min_count": _integer(min_count, "min_count", 0, 1_000_000_000),
        "mad_max": _integer(mad_max, "mad_max", 0, 100),
        "tol": _number(tol, "tol", 0.0, 1.0),
        "max_iters": _integer(max_iters, "max_iters", 1, 100_000),
    }
    if options["tol"] <= 0.0:
        raise ConfigurationError("tol must be greater than 0.")
    cooler = require("cooler", _COOLER_EXTRA)
    np = require("numpy", "bio-analysis")
    clr = _open_cooler(cooler, path)
    weights, _ = cooler.balance_cooler(clr, store=False, **options)
    finite = int(np.isfinite(weights).sum())
    if finite == 0:
        raise ConfigurationError(
            "No bin passed the balancing filters; lower min_nnz or min_count. "
            "The file was not changed."
        )
    try:
        _, stats = cooler.balance_cooler(clr, store=True, **options)
    except OSError as exc:
        raise ConfigurationError(f"The cooler file could not be updated: {path} ({exc})") from exc
    return {
        "status": "completed",
        "operation": "balance_contact_matrix",
        "cool_path": str(path),
        "weight_column": "weight",
        "bins": len(weights),
        "finite_weights": finite,
        "scale": _plain(stats["scale"]),
        "converged": bool(stats["converged"]),
        "variance": _plain(stats["var"]),
        "parameters": {
            "cis_only": options["cis_only"],
            "ignore_diags": options["ignore_diags"],
            "min_nnz": options["min_nnz"],
            "min_count": options["min_count"],
            "mad_max": options["mad_max"],
            "tol": options["tol"],
            "max_iters": options["max_iters"],
        },
        "provenance": _provenance(
            "cooler",
            "iterative correction (Knight-Ruiz style balancing) with marginals rescaled by cooler",
            note="an existing weight column is replaced",
        ),
    }


def calculate_coverage(
    cool_path: str | Path,
    *,
    balanced: bool = False,
    ignore_diags: int = 0,
    max_bins: int = _DEFAULT_MAX_BINS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cis and genome-wide contact sums for each bin of a cooler file, using cooltools."""

    path = _path(cool_path, "cool_path")
    use_weights = _flag(balanced, "balanced")
    diagonals = _integer(ignore_diags, "ignore_diags", 0, 100)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    cooler = require("cooler", _COOLER_EXTRA)
    cooltools = require("cooltools", _COOLTOOLS_EXTRA)
    clr = _open_cooler(cooler, path)
    _weighted(clr, use_weights)
    if int(clr.info["nbins"]) > limit:
        raise ConfigurationError(
            f"The cooler has {clr.info['nbins']} bins, above max_bins={limit}."
        )
    cis, total = cooltools.coverage(
        clr,
        ignore_diags=diagonals,
        clr_weight_name="weight" if use_weights else None,
    )
    bins = clr.bins()[:]
    rows = [
        {
            "chrom": str(row.chrom),
            "start": int(row.start),
            "end": int(row.end),
            "cis": _plain(float(cis_value)),
            "total": _plain(float(total_value)),
        }
        for row, cis_value, total_value in zip(
            bins.itertuples(index=False), cis, total, strict=True
        )
    ]
    return {
        "status": "completed",
        "operation": "calculate_coverage",
        "cool_path": str(path),
        "balanced": use_weights,
        "ignore_diags": diagonals,
        "bins": rows,
        "provenance": _provenance(
            "cooltools",
            "sum of contacts in each bin; cis counts same-chromosome partners only",
            note=(
                "raw counts; each contact is counted on both of its bins"
                if not use_weights
                else "sum of balanced values using the stored weights"
            ),
        ),
    }


def calculate_expected_contacts(
    cool_path: str | Path,
    *,
    balanced: bool = False,
    ignore_diags: int = 2,
    max_bins: int = _DEFAULT_MAX_BINS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average contact frequency at each genomic distance on each chromosome, using cooltools."""

    path = _path(cool_path, "cool_path")
    use_weights = _flag(balanced, "balanced")
    diagonals = _integer(ignore_diags, "ignore_diags", 0, 100)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    cooler = require("cooler", _COOLER_EXTRA)
    cooltools = require("cooltools", _COOLTOOLS_EXTRA)
    clr = _open_cooler(cooler, path)
    _weighted(clr, use_weights)
    if int(clr.info["nbins"]) > limit:
        raise ConfigurationError(
            f"The cooler has {clr.info['nbins']} bins, above max_bins={limit}."
        )
    table = cooltools.expected_cis(
        clr,
        intra_only=True,
        smooth=False,
        aggregate_smoothed=False,
        clr_weight_name="weight" if use_weights else None,
        ignore_diags=diagonals,
    )
    prefix = "balanced" if use_weights else "count"
    rows: list[dict[str, Any]] = []
    sum_column = f"{prefix}.sum"
    avg_column = f"{prefix}.avg"
    for record in table.to_dict("records"):
        rows.append(
            {
                "chrom": str(record["region1"]),
                "distance_bins": int(record["dist"]),
                "distance_bp": int(record["dist_bp"]),
                "valid_pixels": int(record["n_valid"]),
                "sum": _plain(float(record[sum_column])),
                "mean": _plain(float(record[avg_column])),
            }
        )
    return {
        "status": "completed",
        "operation": "calculate_expected_contacts",
        "cool_path": str(path),
        "balanced": use_weights,
        "ignore_diags": diagonals,
        "rows": rows,
        "provenance": _provenance(
            "cooltools",
            "mean contacts on each diagonal: sum divided by the number of valid pixels on it",
            coordinate_system="distance in bins and in bp; intra-chromosomal pixels only",
        ),
    }


def calculate_insulation_score(
    cool_path: str | Path,
    window_bp: int,
    *,
    ignore_diags: int | None = None,
    max_bins: int = _DEFAULT_MAX_BINS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Diamond insulation score and boundary calls for a balanced cooler, using cooltools."""

    path = _path(cool_path, "cool_path")
    window = _integer(window_bp, "window_bp", 1, 100_000_000)
    if ignore_diags is not None:
        ignore_diags = _integer(ignore_diags, "ignore_diags", 0, 100)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    cooler = require("cooler", _COOLER_EXTRA)
    cooltools = require("cooltools", _COOLTOOLS_EXTRA)
    clr = _open_cooler(cooler, path)
    _weighted(clr, True)
    if int(clr.info["nbins"]) > limit:
        raise ConfigurationError(
            f"The cooler has {clr.info['nbins']} bins, above max_bins={limit}."
        )
    try:
        table = cooltools.insulation(clr, [window], ignore_diags=ignore_diags)
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"cooltools rejected the insulation request: {exc}") from exc
    suffix = f"_{window}"
    rows = [
        {
            "chrom": str(record["chrom"]),
            "start": int(record["start"]),
            "end": int(record["end"]),
            "is_bad_bin": bool(record["is_bad_bin"]),
            "log2_insulation_score": _plain(float(record[f"log2_insulation_score{suffix}"])),
            "valid_pixels": _plain(float(record[f"n_valid_pixels{suffix}"])),
            "boundary_strength": _plain(float(record[f"boundary_strength{suffix}"])),
            "is_boundary": bool(record[f"is_boundary{suffix}"]),
        }
        for record in table.to_dict("records")
    ]
    return {
        "status": "completed",
        "operation": "calculate_insulation_score",
        "cool_path": str(path),
        "window_bp": window,
        "boundary_count": sum(1 for row in rows if row["is_boundary"]),
        "bins": rows,
        "provenance": _provenance(
            "cooltools",
            "diamond insulation score on balanced pixels; boundaries thresholded by the Li method",
            coordinate_system="0-based half-open bins",
        ),
    }


def _read_pairs_header(source: Path) -> tuple[list[str], list[str]]:
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    with source.open("r", encoding="utf-8") as handle:
        header, _ = headerops.get_header(handle)
    return header, headerops.extract_column_names(header)


def read_pairs_header(pairs_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read columns, chromosome sizes and sort order from a .pairs header, using pairtools."""

    source = _path(pairs_path, "pairs_path")
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    header, columns = _read_pairs_header(source)
    sorted_values = [
        line.split(":", 1)[1].strip() for line in header if line.startswith("#sorted:")
    ]
    chromsizes = (
        headerops.extract_chromsizes(header)
        if any(line.startswith("#chromsize:") for line in header)
        else {}
    )
    return {
        "status": "completed",
        "operation": "read_pairs_header",
        "pairs": str(source),
        "header_lines": len(header),
        "columns": columns,
        "sorted": sorted_values[0] if sorted_values else None,
        "chromsizes": {str(name): int(length) for name, length in chromsizes.items()},
        "provenance": _provenance(
            "pairtools", "header lines starting with '#' parsed by pairtools"
        ),
    }


def summarize_pairs(
    pairs_path: str | Path, *, chunksize: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """Count mapped, cis, trans and pair-type totals in a .pairs file, using pairtools."""

    source = _path(pairs_path, "pairs_path")
    rows = _integer(chunksize, "chunksize", 1, _MAX_ROWS)
    pandas = require("pandas", _PAIRTOOLS_EXTRA)
    counter_class = require("pairtools.lib.stats", _PAIRTOOLS_EXTRA).PairCounter
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    with source.open("r", encoding="utf-8") as handle:
        header, body = headerops.get_header(handle)
        columns = headerops.extract_column_names(header)
        if "pair_type" not in columns:
            raise ConfigurationError(
                "summarize_pairs() needs a pair_type column, as written by pairtools parse."
            )
        counter = counter_class(
            n_dist_bins_decade=8,
            bytile_dups=False,
            filters=None,
            startup_code="",
            type_cast=(),
            engine="pandas",
        )
        for chunk in pandas.read_table(body, names=columns, chunksize=rows):
            counter.add_pairs_from_dataframe(chunk)
    return {
        "status": "completed",
        "operation": "summarize_pairs",
        "pairs": str(source),
        "total": int(counter["total"]),
        "total_mapped": int(counter["total_mapped"]),
        "total_single_sided_mapped": int(counter["total_single_sided_mapped"]),
        "total_unmapped": int(counter["total_unmapped"]),
        "cis": int(counter["cis"]),
        "trans": int(counter["trans"]),
        "pair_types": {str(key): int(value) for key, value in counter["pair_types"].items()},
        "provenance": _provenance(
            "pairtools",
            "pair counts from pairtools PairCounter; cis pairs share one chromosome on both sides",
        ),
    }


def deduplicate_pairs(
    pairs_path: str | Path,
    *,
    output_path: str | Path | None = None,
    unmapped_path: str | Path | None = None,
    max_mismatch: int = 3,
    method: str = "max",
    chunksize: int = 10_000,
    carryover: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mark PCR duplicate pairs and count them, using pairtools' streaming deduplication."""

    source = _path(pairs_path, "pairs_path")
    target = _new_path(output_path, "output_path") if output_path is not None else None
    unmapped = _new_path(unmapped_path, "unmapped_path") if unmapped_path is not None else None
    mismatch = _integer(max_mismatch, "max_mismatch", 0, 1_000)
    if method not in ("max", "sum"):
        raise ConfigurationError("method must be 'max' or 'sum'.")
    rows = _integer(chunksize, "chunksize", 1, _MAX_ROWS)
    carry = _integer(carryover, "carryover", 0, 1_000_000)
    pairsam_format = require("pairtools.lib.pairsam_format", _PAIRTOOLS_EXTRA)
    dedup = require("pairtools.lib.dedup", _PAIRTOOLS_EXTRA)
    counter_class = require("pairtools.lib.stats", _PAIRTOOLS_EXTRA).PairCounter
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    header, columns = _read_pairs_header(source)
    needed = ("chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2")
    if [name for name in needed if name not in columns]:
        raise ConfigurationError("deduplicate_pairs() needs the chrom, pos and strand columns.")
    counter = counter_class(
        bytile_dups=False, filters=None, startup_code="", type_cast=(), engine="pandas"
    )
    with source.open("r", encoding="utf-8") as handle:
        _, body = headerops.get_header(handle)
        with _writer(target) as unique_out, _writer(unmapped) as unmapped_out:
            if target is not None:
                unique_out.writelines(line + "\n" for line in header)
            if unmapped is not None:
                unmapped_out.writelines(line + "\n" for line in header)
            dedup.streaming_dedup(
                in_stream=body,
                colnames=columns,
                chunksize=rows,
                carryover=carry,
                method=method,
                mark_dups=False,
                max_mismatch=mismatch,
                extra_col_pairs=[],
                unmapped_chrom=pairsam_format.UNMAPPED_CHROM,
                outstream=unique_out,
                outstream_dups=None,
                outstream_unmapped=unmapped_out if unmapped is not None else None,
                keep_parent_id=False,
                out_stat=counter,
                backend="scipy",
                n_proc=1,
            )
    mapped = int(counter["total_mapped"])
    duplicates = int(counter["total_dups"])
    return {
        "status": "completed",
        "operation": "deduplicate_pairs",
        "pairs": str(source),
        "output": str(target) if target is not None else None,
        "unmapped_output": str(unmapped) if unmapped is not None else None,
        "max_mismatch": mismatch,
        "method": method,
        "total": int(counter["total"]),
        "total_mapped": mapped,
        "total_single_sided_mapped": int(counter["total_single_sided_mapped"]),
        "total_nodups": int(counter["total_nodups"]),
        "total_dups": duplicates,
        "duplicate_fraction": _plain(duplicates / mapped) if mapped else None,
        "provenance": _provenance(
            "pairtools",
            "pairs within max_mismatch bp on both sides form one cluster; one pair is kept",
            note="single-sided pairs are counted but written only when unmapped_path is given",
        ),
    }


class _writer:
    """Context manager for an optional output file; None writes nowhere."""

    def __init__(self, path: Path | None) -> None:
        self._path = path
        self._handle: Any = None

    def __enter__(self) -> Any:
        if self._path is None:
            self._handle = open(os.devnull, "w", encoding="utf-8")
        else:
            self._handle = self._path.open("x", encoding="utf-8")
        return self._handle

    def __exit__(self, *exc: object) -> None:
        self._handle.close()


def _fragment_arrays(starts: object, ends: object, max_fragments: int) -> tuple[Any, Any]:
    np = require("numpy", "bio-analysis")
    if isinstance(starts, (str, bytes)) or isinstance(ends, (str, bytes)):
        raise ConfigurationError("starts and ends must be lists of integers.")
    if not isinstance(starts, Sequence) or not isinstance(ends, Sequence):
        raise ConfigurationError("starts and ends must be lists of integers.")
    if len(starts) != len(ends) or not 1 <= len(starts) <= max_fragments:
        raise ConfigurationError(
            f"starts and ends must have the same length, between 1 and {max_fragments}."
        )
    pairs = []
    for start, end in zip(starts, ends, strict=True):
        left = _integer(start, "fragment start", 0, _MAX_INT32 - 1)
        right = _integer(end, "fragment end", 1, _MAX_INT32)
        if right <= left:
            raise ConfigurationError("every fragment end must be greater than its start.")
        pairs.append((left, right))
    start_array = np.array(sorted(left for left, _ in pairs), dtype="i4")
    end_array = np.array(sorted(right for _, right in pairs), dtype="i4")
    return start_array, end_array


def _pileup_segments(pileup: tuple[Any, Any]) -> list[dict[str, Any]]:
    ends, values = pileup
    segments: list[dict[str, Any]] = []
    previous = 0
    for end, value in zip(ends.tolist(), values.tolist(), strict=True):
        segments.append({"start": int(previous), "end": int(end), "value": _plain(float(value))})
        previous = int(end)
    return segments


def build_fragment_pileup(
    starts: Sequence[int],
    ends: Sequence[int],
    *,
    scale_factor: float = 1.0,
    baseline_value: float = 0.0,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Piecewise-constant coverage of fragments on one chromosome, using MACS3's quick_pileup."""

    limit = _integer(max_fragments, "max_fragments", 1, _MAX_FRAGMENTS)
    scale = _number(scale_factor, "scale_factor", 0.0, 1e9)
    baseline = _number(baseline_value, "baseline_value", 0.0, 1e9)
    start_array, end_array = _fragment_arrays(starts, ends, limit)
    pileup = require("MACS3.Signal.Pileup", _MACS3_EXTRA).quick_pileup(
        start_array, end_array, scale, baseline
    )
    segments = _pileup_segments(pileup)
    return {
        "status": "completed",
        "operation": "build_fragment_pileup",
        "fragments": len(start_array),
        "segments": segments,
        "provenance": _provenance(
            "macs3",
            "fragment coverage as contiguous segments from position 0 to the largest fragment end",
            coordinate_system="0-based half-open segments",
            baseline="baseline_value is the minimum value of every segment",
        ),
    }


def build_extended_tag_pileup(
    positions: Sequence[int],
    extension: int,
    *,
    max_tags: int = _DEFAULT_MAX_FRAGMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coverage of single-end tags extended symmetrically by a fixed length, using MACS3."""

    limit = _integer(max_tags, "max_tags", 1, _MAX_FRAGMENTS)
    width = _integer(extension, "extension", 0, 1_000_000)
    if isinstance(positions, (str, bytes)) or not isinstance(positions, Sequence):
        raise ConfigurationError("positions must be a list of integers.")
    if not 1 <= len(positions) <= limit:
        raise ConfigurationError(f"positions must contain between 1 and {limit} tags.")
    values = sorted(_integer(position, "tag position", 0, _MAX_INT32 - 1) for position in positions)
    np = require("numpy", "bio-analysis")
    array = np.array(values, dtype="i4")
    pileup = require("MACS3.Signal.Pileup", _MACS3_EXTRA).naive_quick_pileup(array, width)
    return {
        "status": "completed",
        "operation": "build_extended_tag_pileup",
        "tags": len(values),
        "extension": width,
        "segments": _pileup_segments(pileup),
        "provenance": _provenance(
            "macs3",
            "each tag covers [position - extension, position + extension); coverage starts at 0",
            coordinate_system="0-based half-open segments",
        ),
    }


def call_fragment_peaks(
    starts: Sequence[int],
    ends: Sequence[int],
    *,
    min_value: float,
    max_gap: int = 50,
    min_length: int = 200,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call coverage summits from fragments on one chromosome, using MACS3's naive peak caller."""

    limit = _integer(max_fragments, "max_fragments", 1, _MAX_FRAGMENTS)
    threshold = _number(min_value, "min_value", 0.0, 1e9)
    gap = _integer(max_gap, "max_gap", 0, 1_000_000)
    length = _integer(min_length, "min_length", 1, 1_000_000_000)
    start_array, end_array = _fragment_arrays(starts, ends, limit)
    macs = require("MACS3.Signal.Pileup", _MACS3_EXTRA)
    pileup = macs.quick_pileup(start_array, end_array, 1.0, 0.0)
    summits = macs.naive_call_peaks([pileup[0], pileup[1]], threshold, 1e30, gap, length)
    return {
        "status": "completed",
        "operation": "call_fragment_peaks",
        "fragments": len(start_array),
        "min_value": threshold,
        "max_gap": gap,
        "min_length": length,
        "peaks": [
            {"summit": int(summit), "height": _plain(float(height))} for summit, height in summits
        ],
        "provenance": _provenance(
            "macs3",
            "segments above min_value merge across gaps up to max_gap; kept if long enough",
            note="a naive caller on the fragment pileup; MACS3 statistical peak calling is not run",
        ),
    }


def summarize_bigwig_bins(
    bigwig_files: Sequence[str | Path],
    *,
    bin_size: int = 1000,
    max_bins: int = _DEFAULT_MAX_BINS,
    number_of_processors: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean signal of one or more bigWig files in fixed-width genome bins, using deeptools."""

    paths = _paths(bigwig_files, "bigwig_files")
    size = _integer(bin_size, "bin_size", 1, _MAX_BIN_SIZE)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    processors = _integer(number_of_processors, "number_of_processors", 1, _MAX_PROCESSORS)
    chromosomes = _shared_bigwig_chromosomes(paths)
    expected = _check_bin_count(chromosomes, size, limit)
    get_scores = require("deeptools.getScorePerBigWigBin", _DEEPTOOLS_EXTRA).getScorePerBin
    with tempfile.TemporaryDirectory(prefix="dnakit-deeptools-") as scratch:
        raw = Path(scratch) / "bins.bed"
        matrix = get_scores(
            [str(path) for path in paths],
            size,
            numberOfProcessors=processors,
            out_file_for_raw_data=str(raw),
        )
        bins = _read_bins(raw, len(paths), lambda text: _plain(float(text)))
    if len(bins) != expected or len(matrix) != expected:
        raise BackendExecutionError(
            f"deeptools returned {len(bins)} bins; expected {expected} from the chromosome sizes."
        )
    return {
        "status": "completed",
        "operation": "summarize_bigwig_bins",
        "bin_size": size,
        "files": [str(path) for path in paths],
        "chromosomes": _chromosome_rows(chromosomes),
        "bins": bins,
        "provenance": _provenance(
            "deeptools",
            "mean signal over bases that have data; null when a bin has no data",
            coordinate_system="0-based half-open; the last bin of a chromosome may be shorter",
            scope="chromosomes with the same name and length in every bigWig file",
        ),
    }


def count_bam_reads_per_bin(
    bam_files: Sequence[str | Path],
    *,
    bin_size: int = 1000,
    max_bins: int = _DEFAULT_MAX_BINS,
    min_mapping_quality: int = 0,
    number_of_processors: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count indexed BAM alignments that overlap each fixed-width genome bin, using deeptools."""

    paths = _paths(bam_files, "bam_files")
    size = _integer(bin_size, "bin_size", 1, _MAX_BIN_SIZE)
    limit = _integer(max_bins, "max_bins", 1, _MAX_BINS_LIMIT)
    quality = _integer(min_mapping_quality, "min_mapping_quality", 0, _MAX_MAPPING_QUALITY)
    processors = _integer(number_of_processors, "number_of_processors", 1, _MAX_PROCESSORS)
    chromosomes = _shared_bam_chromosomes(paths)
    expected = _check_bin_count(chromosomes, size, limit)
    counter = require("deeptools.countReadsPerBin", _DEEPTOOLS_EXTRA).CountReadsPerBin
    with tempfile.TemporaryDirectory(prefix="dnakit-deeptools-") as scratch:
        raw = Path(scratch) / "bins.bed"
        try:
            matrix = counter(
                [str(path) for path in paths],
                binLength=size,
                stepSize=size,
                numberOfProcessors=processors,
                minMappingQuality=quality,
                out_file_for_raw_data=str(raw),
            ).run()
        except ValueError as exc:
            raise ConfigurationError(f"deeptools rejected the bin layout: {exc}") from exc
        bins = _read_bins(raw, len(paths), lambda text: int(float(text)))
    if len(bins) != expected or len(matrix) != expected:
        raise BackendExecutionError(
            f"deeptools returned {len(bins)} bins; expected {expected} from the chromosome sizes."
        )
    for row in bins:
        row["counts"] = row.pop("values")
    return {
        "status": "completed",
        "operation": "count_bam_reads_per_bin",
        "bin_size": size,
        "files": [str(path) for path in paths],
        "chromosomes": _chromosome_rows(chromosomes),
        "bins": bins,
        "provenance": _provenance(
            "deeptools",
            "alignments overlapping each bin; each mate counts on its own, without read extension",
            coordinate_system="0-based half-open; the last bin of a chromosome may be shorter",
            scope="chromosomes with the same name and length in every BAM file",
        ),
    }


def estimate_fragment_lengths(
    bam_file: str | Path,
    *,
    number_of_processors: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate fragment and read length distributions from an indexed BAM, using deeptools."""

    path = _path(bam_file, "bam_file")
    processors = _integer(number_of_processors, "number_of_processors", 1, _MAX_PROCESSORS)
    _shared_bam_chromosomes([path])
    estimate = require(
        "deeptools.getFragmentAndReadSize", _DEEPTOOLS_EXTRA
    ).get_read_and_fragment_length
    fragment, read = estimate(
        str(path), return_lengths=False, numberOfProcessors=processors, verbose=False
    )
    return {
        "status": "completed",
        "operation": "estimate_fragment_lengths",
        "bam": str(path),
        "fragment_length": _summary(fragment),
        "read_length": _summary(read),
        "provenance": _provenance(
            "deeptools",
            "template lengths of sampled proper-pair read1 alignments; "
            "falls back to all mapped alignments when the file has no proper pairs",
            sampling="estimates come from a sample of genome regions, not from every alignment",
        ),
    }


def _shell_safe_output(path: Path, name: str) -> str:
    """methylpy runs sort through the shell, so output paths use a restricted character set."""

    text = str(path)
    if not _SHELL_SAFE_PATH.fullmatch(text):
        raise ConfigurationError(
            f"{name} may only contain letters, digits and ._/@+=,:- characters: {text}"
        )
    return text


def _upstream_call(function: Any, *args: Any, **kwargs: Any) -> tuple[Any, str]:
    """Run a methylpy function that prints progress and may call sys.exit on bad input."""

    captured = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured), contextlib.redirect_stderr(captured):
            result = function(*args, **kwargs)
    except SystemExit as exc:
        lines = captured.getvalue().strip().splitlines()
        detail = lines[-1] if lines else str(exc.code)
        raise ConfigurationError(f"methylpy rejected the input: {detail}") from exc
    return result, captured.getvalue()


def _line_count(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def filter_allc_sites(
    allc_path: str | Path,
    output_path: str | Path,
    *,
    chromosomes: Sequence[str] | None = None,
    mc_type: str | None = None,
    min_coverage: int = 0,
    max_coverage: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep allc cytosine records that match a chromosome list, context and coverage range."""

    source = _path(allc_path, "allc_path")
    target = _new_path(output_path, "output_path")
    names: list[str] | None = None
    if chromosomes is not None:
        if isinstance(chromosomes, (str, Path)) or not isinstance(chromosomes, Sequence):
            raise ConfigurationError("chromosomes must be a list of chromosome names.")
        names = [_text(name, "chromosome name") for name in chromosomes]
        if not names:
            raise ConfigurationError("chromosomes must not be empty when given.")
    context = None if mc_type is None else _text(mc_type, "mc_type")
    minimum = _integer(min_coverage, "min_coverage", 0, 1_000_000_000)
    maximum = None
    if max_coverage is not None:
        maximum = _integer(max_coverage, "max_coverage", minimum, 1_000_000_000)
    utilities = require("methylpy.utilities", _METHYLPY_EXTRA)
    _upstream_call(
        utilities.filter_allc_files,
        [str(source)],
        [str(target)],
        num_procs=1,
        mc_type=[context] if context is not None else None,
        chroms=names,
        compress_output=False,
        min_cov=minimum,
        max_cov=maximum,
    )
    return {
        "status": "completed",
        "operation": "filter_allc_sites",
        "allc": str(source),
        "output": str(target),
        "input_sites": _line_count(source),
        "kept_sites": _line_count(target),
        "chromosomes": names,
        "mc_type": context,
        "min_coverage": minimum,
        "max_coverage": maximum,
        "provenance": _provenance(
            "methylpy",
            "allc records kept when chromosome, cytosine context and total coverage all match",
            coordinate_system="allc positions as written in the input (1-based)",
        ),
    }


def calculate_non_conversion_rate(
    allc_path: str | Path, control_region: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Estimate the bisulfite non-conversion rate from unmethylated control cytosines."""

    source = _path(allc_path, "allc_path")
    region = _text(control_region, "control_region")
    if not _CONTROL_REGION.fullmatch(region):
        raise ConfigurationError("control_region must look like 'chrom' or 'chrom:start-end'.")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    with tempfile.TemporaryDirectory(prefix="dnakit-methylpy-") as scratch:
        # methylpy writes its index next to the file it reads, so use a linked copy in scratch.
        link = Path(scratch) / ("allc" + "".join(source.suffixes))
        try:
            os.symlink(source, link)
        except OSError:
            shutil.copyfile(source, link)
        rate, _ = _upstream_call(
            call_mc.calculate_non_conversion_rate, region, str(link), None, False
        )
    if rate is None:
        raise ConfigurationError("methylpy found no cytosine reads in the control region.")
    return {
        "status": "completed",
        "operation": "calculate_non_conversion_rate",
        "allc": str(source),
        "control_region": region,
        "non_conversion_rate": _plain(float(rate)),
        "provenance": _provenance(
            "methylpy",
            "methylated reads divided by all cytosine reads in the control region, all contexts",
            note="chromosome names are matched exactly as they are written in the allc file",
        ),
    }


def calculate_site_binomial_pvalues(
    allc_path: str | Path,
    output_path: str | Path,
    *,
    non_conversion: float,
    min_coverage: int = 1,
    sort_buffer: str = "64M",
    max_report: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """One-sided binomial P-value per cytosine, testing methylation above background."""

    source = _path(allc_path, "allc_path")
    target = _new_path(output_path, "output_path")
    rate = _number(non_conversion, "non_conversion", 0.0, 1.0)
    minimum = _integer(min_coverage, "min_coverage", 1, 1_000_000)
    report = _integer(max_report, "max_report", 0, 1_000)
    buffer = _text(sort_buffer, "sort_buffer")
    if not _SORT_BUFFER.fullmatch(buffer):
        raise ConfigurationError("sort_buffer must look like '64M' or '1G'.")
    output_text = _shell_safe_output(target, "output_path")
    call_mc = require("methylpy.call_mc_se", _METHYLPY_EXTRA)
    counts, _ = _upstream_call(
        call_mc.allc_run_binom_tests,
        str(source),
        output_text,
        rate,
        min_cov=minimum,
        sort_mem=buffer,
    )
    lowest: list[dict[str, Any]] = []
    tested = 0
    written = 0
    with target.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            written += 1
            p_value = float(fields[6])
            if p_value < 2.0:
                tested += 1
            if len(lowest) < report:
                lowest.append(
                    {
                        "chrom": fields[0],
                        "pos": int(fields[1]),
                        "strand": fields[2],
                        "mc_class": fields[3],
                        "methylated": int(fields[4]),
                        "coverage": int(fields[5]),
                        "p_value": _plain(p_value),
                    }
                )
    return {
        "status": "completed",
        "operation": "calculate_site_binomial_pvalues",
        "allc": str(source),
        "output": str(target),
        "non_conversion_rate": rate,
        "min_coverage": minimum,
        "sites_written": written,
        "tested_sites": tested,
        "mc_class_counts": {str(key): int(value) for key, value in counts.items()},
        "lowest_p_values": lowest,
        "provenance": _provenance(
            "methylpy",
            "P(X >= methylated reads) for a binomial with p = non-conversion rate",
            note=(
                "cytosines with no methylated reads or coverage below min_coverage get p-value 2.0 "
                "and are not tested; the output is sorted by p-value"
            ),
        ),
    }


def _summary(values: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if values is None:
        return None
    return {key: _plain(value) for key, value in values.items()}


__all__ = [
    "balance_contact_matrix",
    "build_contact_matrix",
    "build_extended_tag_pileup",
    "build_fragment_pileup",
    "calculate_coverage",
    "calculate_expected_contacts",
    "calculate_insulation_score",
    "calculate_non_conversion_rate",
    "calculate_site_binomial_pvalues",
    "call_fragment_peaks",
    "count_bam_reads_per_bin",
    "deduplicate_pairs",
    "estimate_fragment_lengths",
    "fetch_contact_matrix",
    "filter_allc_sites",
    "read_pairs_header",
    "summarize_bigwig_bins",
    "summarize_pairs",
]

_register_public_verbose(globals(), tuple(__all__))
