"""Chromosome-size, interval, FASTA-length and contact-matrix file functions from cooler.

Each public function imports cooler (or pyfaidx) only when it is called. Missing
dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Inputs are checked before upstream calls; upstream ``ValueError``, ``KeyError``
and ``OSError`` are reported as ``ConfigurationError``.

Output files are created as new paths and are never overwritten.
"""

from __future__ import annotations

import re
import shutil
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_COOLER_EXTRA = "bio-cooler"
_PYFAIDX_EXTRA = "bio-pyfaidx"
_MAX_CHROMOSOMES = 10_000
_MAX_CHROM_LENGTH = 2_000_000_000
_MAX_PATTERNS = 16
_MAX_FACTOR = 1_000
_MAX_CHUNK = 100_000_000
_MAX_MERGE_INPUTS = 64
_MAX_RESOLUTIONS = 32
_MAX_RESOLUTION = 1_000_000_000
_MAX_PIXEL_ROWS = 10_000
_MAX_FASTA_NAMES = 10_000
_DATABASE = re.compile(r"^[A-Za-z0-9_]{2,64}$")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\x00\r\n\t"):
        raise ConfigurationError(f"{name} must be a non-empty string without tabs or newlines.")
    return value


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _new(value: object, name: str) -> Path:
    """Output path that must not exist yet; its parent directory must exist."""

    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The directory for {name} does not exist: {path.parent}")
    return path


def _names(value: object, name: str, limit: int) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of chromosome names.")
    if not 1 <= len(value) <= limit:
        raise ConfigurationError(f"{name} must contain 1 to {limit} names.")
    names = [_text(item, name) for item in value]
    if len(set(names)) != len(names):
        raise ConfigurationError(f"{name} must not repeat a name.")
    return names


def _lengths(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CHROMOSOMES:
        raise ConfigurationError(
            f"{name} must map 1 to {_MAX_CHROMOSOMES} chromosome names to lengths."
        )
    return {
        _text(key, f"{name} chromosome"): _integer(length, f"length of {key}", 1, _MAX_CHROM_LENGTH)
        for key, length in value.items()
    }


def _patterns(value: object) -> tuple[str, ...] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("name_patterns must be a list of regular expressions.")
    if not 1 <= len(value) <= _MAX_PATTERNS:
        raise ConfigurationError(f"name_patterns must contain 1 to {_MAX_PATTERNS} patterns.")
    patterns: list[str] = []
    for item in value:
        text = _text(item, "name_patterns")
        try:
            re.compile(text)
        except re.error as exc:
            raise ConfigurationError(f"Invalid name pattern {text!r}: {exc}") from exc
        patterns.append(text)
    return tuple(patterns)


def _plain(value: Any) -> Any:
    """Convert numpy scalars to Python values and non-finite floats to None."""

    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not (value == value and abs(value) != float("inf")):
        return None
    return value


def _records(frame: Any) -> list[dict[str, Any]]:
    return [
        {str(key): _plain(value) for key, value in row.items()}
        for row in frame.to_dict(orient="records")
    ]


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _sizes_rows(sizes: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [{"name": str(name), "length": int(length)} for name, length in sizes.items()]


def _cooler_info(cool_path: Path) -> dict[str, Any]:
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        clr = cooler.Cooler(str(cool_path))
        info = clr.info
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(
            f"Input is not a readable cooler file: {cool_path} ({exc})"
        ) from exc
    return {
        "binsize": int(clr.binsize),
        "bins": int(info["nbins"]),
        "nonzero_pixels": int(info["nnz"]),
        "total_contacts": int(info["sum"]),
        "chromosomes": list(clr.chromnames),
    }


def read_chrom_sizes(
    path: str | Path,
    *,
    all_names: bool = False,
    name_patterns: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read chromosome lengths from a UCSC chrom.sizes or chromInfo file, using cooler."""

    source = _existing(path, "path")
    keep_all = _flag(all_names, "all_names")
    patterns = _patterns(name_patterns)
    cooler_util = require("cooler.util", _COOLER_EXTRA)
    kwargs: dict[str, Any] = {"all_names": keep_all}
    if patterns is not None:
        kwargs["name_patterns"] = patterns
    try:
        sizes = cooler_util.read_chromsizes(str(source), **kwargs)
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"Could not read chromosome sizes: {source} ({exc})") from exc
    rows = _sizes_rows(sizes)
    if not rows:
        raise ConfigurationError("No chromosome passed the name filter.")
    return {
        "status": "completed",
        "operation": "read_chrom_sizes",
        "path": str(source),
        "all_names": keep_all,
        "chromosome_count": len(rows),
        "chromosomes": rows,
        "provenance": _provenance(
            "cooler",
            "two-column or UCSC chromInfo table read by cooler.util.read_chromsizes",
            order="file order is kept",
        ),
    }


def fetch_chrom_sizes(
    database: str,
    *,
    all_names: bool = False,
    name_patterns: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download chromosome sizes for a UCSC genome assembly, using cooler (needs network)."""

    db = _text(database, "database")
    if not _DATABASE.match(db):
        raise ConfigurationError("database must be a UCSC assembly name such as hg38.")
    keep_all = _flag(all_names, "all_names")
    patterns = _patterns(name_patterns)
    cooler_util = require("cooler.util", _COOLER_EXTRA)
    kwargs: dict[str, Any] = {"all_names": keep_all}
    if patterns is not None:
        kwargs["name_patterns"] = patterns
    try:
        sizes = cooler_util.fetch_chromsizes(db, **kwargs)
    except OSError as exc:
        raise BackendUnavailableError(
            f"UCSC chromosome sizes for {db} could not be downloaded: {exc}",
            hint="Check the network connection, or use read_chrom_sizes() on a local file.",
        ) from exc
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"UCSC chromosome sizes for {db} were not usable: {exc}") from exc
    rows = _sizes_rows(sizes)
    if not rows:
        raise ConfigurationError("No chromosome passed the name filter.")
    return {
        "status": "completed",
        "operation": "fetch_chrom_sizes",
        "database": db,
        "source": f"https://hgdownload.soe.ucsc.edu/goldenPath/{db}/database/chromInfo.txt.gz",
        "all_names": keep_all,
        "chromosome_count": len(rows),
        "chromosomes": rows,
        "provenance": _provenance(
            "cooler",
            "UCSC chromInfo table downloaded by cooler.util.fetch_chromsizes",
            network="one HTTP GET per call; no cache is written",
        ),
    }


def parse_genomic_region(
    region: str,
    *,
    chromsizes: Mapping[str, int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a UCSC-style region such as chr1:100-200 into half-open 0-based coordinates."""

    text = _text(region, "region")
    sizes = None if chromsizes is None else _lengths(chromsizes, "chromsizes")
    cooler_util = require("cooler.util", _COOLER_EXTRA)
    try:
        chrom, start, end = cooler_util.parse_region(text, sizes)
    except (ValueError, KeyError, TypeError) as exc:
        raise ConfigurationError(f"Region is not valid: {text} ({exc})") from exc
    start, end = int(start), int(end)
    return {
        "status": "completed",
        "operation": "parse_genomic_region",
        "region": text,
        "chrom": str(chrom),
        "start": start,
        "end": end,
        "length": end - start,
        "provenance": _provenance(
            "cooler",
            "half-open intervals with 0-based starts and 1-based ends (cooler.util.parse_region)",
        ),
    }


def query_bed_intervals(
    bed_path: str | Path,
    region: str,
    *,
    chromsizes: Mapping[str, int],
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the BED intervals of one region from a sorted, non-overlapping BED file."""

    source = _existing(bed_path, "bed_path")
    text = _text(region, "region")
    sizes = _lengths(chromsizes, "chromsizes")
    pandas = require("pandas", _COOLER_EXTRA)
    cooler_util = require("cooler.util", _COOLER_EXTRA)
    try:
        table = pandas.read_csv(
            source,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 1, 2],
            names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": "int64", "end": "int64"},
        )
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"BED file needs chrom, start and end columns: {exc}") from exc
    if table.empty:
        raise ConfigurationError("The BED file has no intervals.")
    if (table["end"] <= table["start"]).any() or (table["start"] < 0).any():
        raise ConfigurationError("Every BED interval must have 0 <= start < end.")
    if (table.groupby("chrom", sort=False)["start"].diff().fillna(0) < 0).any():
        raise ConfigurationError("BED intervals must be sorted by start within each chromosome.")
    previous_end = table.groupby("chrom", sort=False)["end"].shift(1)
    if (table["start"] < previous_end.fillna(0)).any():
        raise ConfigurationError("BED intervals must not overlap within a chromosome.")
    chrom = text.split(":", 1)[0]
    if chrom not in sizes:
        raise ConfigurationError(f"Chromosome is not in chromsizes: {chrom}")
    if chrom not in set(table["chrom"]):
        rows: list[dict[str, Any]] = []
    else:
        try:
            grouped = table.groupby("chrom", sort=False)
            selected = cooler_util.bedslice(grouped, sizes, text)
        except (ValueError, KeyError, TypeError) as exc:
            raise ConfigurationError(
                f"Region is not valid for this BED file: {text} ({exc})"
            ) from exc
        rows = _records(selected)
    return {
        "status": "completed",
        "operation": "query_bed_intervals",
        "bed_path": str(source),
        "region": text,
        "interval_count": len(rows),
        "intervals": rows,
        "provenance": _provenance(
            "cooler",
            "range query on non-overlapping intervals (cooler.util.bedslice)",
            coordinates="half-open, 0-based starts",
        ),
    }


def rename_cooler_chromosomes(
    cool_path: str | Path,
    output_path: str | Path,
    rename_map: Mapping[str, str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy a cooler file and rename its chromosomes in the copy, using cooler."""

    source = _existing(cool_path, "cool_path")
    target = _new(output_path, "output_path")
    if not isinstance(rename_map, Mapping) or not 1 <= len(rename_map) <= _MAX_CHROMOSOMES:
        raise ConfigurationError("rename_map must map 1 to 10000 old names to new names.")
    mapping = {_text(old, "old name"): _text(new, "new name") for old, new in rename_map.items()}
    cooler = require("cooler", _COOLER_EXTRA)
    before = list(_cooler_info(source)["chromosomes"])
    unknown = sorted(name for name in mapping if name not in before)
    if unknown:
        raise ConfigurationError(f"rename_map names are not chromosomes of this cooler: {unknown}")
    after = [mapping.get(name, name) for name in before]
    if len(set(after)) != len(after):
        raise ConfigurationError("The renamed chromosome names must be unique.")
    shutil.copyfile(source, target)
    try:
        clr = cooler.Cooler(str(target))
        cooler.rename_chroms(clr, mapping)
    except (ValueError, KeyError, OSError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler could not rename the chromosomes: {exc}") from exc
    return {
        "status": "completed",
        "operation": "rename_cooler_chromosomes",
        "input": str(source),
        "output": str(target),
        "renamed": {old: mapping[old] for old in sorted(mapping)},
        "chromosomes_before": before,
        "chromosomes_after": after,
        "provenance": _provenance(
            "cooler",
            "chromosome names rewritten in a copy by cooler.rename_chroms",
            note="the input file is not modified",
        ),
    }


def coarsen_contact_matrix(
    cool_path: str | Path,
    output_path: str | Path,
    *,
    factor: int,
    chunksize: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Coarsen a cooler to a lower resolution by an integer factor, using cooler."""

    source = _existing(cool_path, "cool_path")
    target = _new(output_path, "output_path")
    k = _integer(factor, "factor", 2, _MAX_FACTOR)
    rows_per_chunk = _integer(chunksize, "chunksize", 1, _MAX_CHUNK)
    cooler = require("cooler", _COOLER_EXTRA)
    before = _cooler_info(source)
    try:
        cooler.coarsen_cooler(str(source), str(target), factor=k, chunksize=rows_per_chunk, nproc=1)
        after = _cooler_info(target)
    except (ValueError, KeyError, OSError, ConfigurationError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler could not coarsen this file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "coarsen_contact_matrix",
        "input": str(source),
        "output": str(target),
        "factor": k,
        "input_binsize": before["binsize"],
        "output_binsize": after["binsize"],
        "output_bins": after["bins"],
        "output_nonzero_pixels": after["nonzero_pixels"],
        "output_total_contacts": after["total_contacts"],
        "provenance": _provenance(
            "cooler",
            "k-by-k pooling of pixels, aggregated by sum (cooler.coarsen_cooler)",
            note="each chromosome block is coarsened separately",
        ),
    }


def merge_contact_matrices(
    cool_paths: Sequence[str | Path],
    output_path: str | Path,
    *,
    mergebuf: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sum several coolers that share identical bins into one file, using cooler."""

    if isinstance(cool_paths, (str, bytes, Path)) or not isinstance(cool_paths, Sequence):
        raise ConfigurationError("cool_paths must be a list of cooler file paths.")
    if not 2 <= len(cool_paths) <= _MAX_MERGE_INPUTS:
        raise ConfigurationError(f"cool_paths must contain 2 to {_MAX_MERGE_INPUTS} files.")
    sources = [_existing(item, "cool_paths") for item in cool_paths]
    if len(set(sources)) != len(sources):
        raise ConfigurationError("cool_paths must not list the same file twice.")
    target = _new(output_path, "output_path")
    buffer_size = _integer(mergebuf, "mergebuf", 1, 1_000_000_000)
    cooler = require("cooler", _COOLER_EXTRA)
    infos = [_cooler_info(path) for path in sources]
    try:
        cooler.merge_coolers(
            str(target),
            [str(path) for path in sources],
            mergebuf=buffer_size,
        )
        after = _cooler_info(target)
    except (ValueError, KeyError, OSError, ConfigurationError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler could not merge these files: {exc}") from exc
    return {
        "status": "completed",
        "operation": "merge_contact_matrices",
        "inputs": [str(path) for path in sources],
        "output": str(target),
        "input_total_contacts": sum(info["total_contacts"] for info in infos),
        "output_total_contacts": after["total_contacts"],
        "output_bins": after["bins"],
        "output_nonzero_pixels": after["nonzero_pixels"],
        "provenance": _provenance(
            "cooler",
            "pixel counts summed over inputs with identical bins (cooler.merge_coolers)",
        ),
    }


def zoomify_contact_matrix(
    cool_paths: str | Path | Sequence[str | Path],
    output_path: str | Path,
    *,
    resolutions: Sequence[int],
    chunksize: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a multi-resolution (.mcool) file from one base cooler, using cooler."""

    if isinstance(cool_paths, (str, Path)):
        bases = [_existing(cool_paths, "cool_paths")]
    elif isinstance(cool_paths, Sequence) and cool_paths:
        bases = [_existing(item, "cool_paths") for item in cool_paths]
    else:
        raise ConfigurationError("cool_paths must be a cooler path or a non-empty list of paths.")
    if len(bases) > _MAX_MERGE_INPUTS:
        raise ConfigurationError(f"cool_paths must contain at most {_MAX_MERGE_INPUTS} files.")
    if isinstance(resolutions, (str, bytes)) or not isinstance(resolutions, Sequence):
        raise ConfigurationError("resolutions must be a list of bin sizes in base pairs.")
    if not 1 <= len(resolutions) <= _MAX_RESOLUTIONS:
        raise ConfigurationError(f"resolutions must contain 1 to {_MAX_RESOLUTIONS} sizes.")
    sizes = [_integer(r, "resolution", 1, _MAX_RESOLUTION) for r in resolutions]
    if sizes != sorted(set(sizes)):
        raise ConfigurationError("resolutions must be unique and in ascending order.")
    target = _new(output_path, "output_path")
    rows_per_chunk = _integer(chunksize, "chunksize", 1, _MAX_CHUNK)
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        cooler.zoomify_cooler(
            [str(path) for path in bases],
            str(target),
            sizes,
            rows_per_chunk,
            nproc=1,
        )
    except (ValueError, KeyError, OSError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"cooler could not zoomify these files: {exc}") from exc
    levels = []
    for size in sizes:
        info = _cooler_info_uri(target, size)
        levels.append({"resolution": size, **info})
    return {
        "status": "completed",
        "operation": "zoomify_contact_matrix",
        "inputs": [str(path) for path in bases],
        "output": str(target),
        "resolutions": levels,
        "provenance": _provenance(
            "cooler",
            "recursive coarsening into an MCOOL v2 file (cooler.zoomify_cooler)",
        ),
    }


def _cooler_info_uri(target: Path, resolution: int) -> dict[str, Any]:
    cooler = require("cooler", _COOLER_EXTRA)
    clr = cooler.Cooler(f"{target}::/resolutions/{resolution}")
    info = clr.info
    return {
        "binsize": int(clr.binsize),
        "bins": int(info["nbins"]),
        "nonzero_pixels": int(info["nnz"]),
        "total_contacts": int(info["sum"]),
    }


def load_fasta_lengths(
    fasta_path: str | Path,
    names: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the lengths of named FASTA records without reading sequence into memory."""

    source = _existing(fasta_path, "fasta_path")
    wanted = _names(names, "names", _MAX_FASTA_NAMES)
    require("pyfaidx", _PYFAIDX_EXTRA)
    cooler_util = require("cooler.util", _COOLER_EXTRA)
    try:
        records = cooler_util.load_fasta(wanted, str(source))
    except KeyError as exc:
        raise ConfigurationError(f"A requested name is not in the FASTA file: {exc}") from exc
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Could not index the FASTA file: {exc}") from exc
    rows = [{"name": name, "length": len(records[name])} for name in wanted]
    return {
        "status": "completed",
        "operation": "load_fasta_lengths",
        "fasta_path": str(source),
        "sequence_count": len(rows),
        "sequences": rows,
        "provenance": _provenance(
            "cooler",
            "lazy FASTA records opened by cooler.util.load_fasta (pyfaidx index)",
            note="sequence bases are not returned",
        ),
    }


def read_cooler_pixels(
    cool_path: str | Path,
    *,
    limit: int = 100,
    join: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the first nonzero pixels of a cooler, optionally joined to bin coordinates."""

    source = _existing(cool_path, "cool_path")
    rows_limit = _integer(limit, "limit", 1, _MAX_PIXEL_ROWS)
    use_join = _flag(join, "join")
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        clr = cooler.Cooler(str(source))
        total = int(clr.info["nnz"])
        frame = clr.pixels(join=use_join)[0 : min(rows_limit, total)]
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {source} ({exc})") from exc
    records = _records(frame)
    return {
        "status": "completed",
        "operation": "read_cooler_pixels",
        "cool_path": str(source),
        "join": use_join,
        "nonzero_pixels": total,
        "returned": len(records),
        "columns": [str(column) for column in frame.columns],
        "pixels": records,
        "provenance": _provenance(
            "cooler",
            "upper-triangle pixel table in bin order (cooler.Cooler.pixels, api.pixels)",
            order="rows follow the stored bin1_id, bin2_id order",
        ),
    }


def read_cooler_chroms(cool_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """List the chromosomes stored in a cooler file, in matrix order, using cooler."""

    source = _existing(cool_path, "cool_path")
    cooler = require("cooler", _COOLER_EXTRA)
    try:
        frame = cooler.Cooler(str(source)).chroms()[:]
    except (OSError, KeyError, ValueError) as exc:
        raise ConfigurationError(f"Input is not a readable cooler file: {source} ({exc})") from exc
    rows = [
        {"name": str(name), "length": int(length)}
        for name, length in zip(frame["name"], frame["length"], strict=True)
    ]
    return {
        "status": "completed",
        "operation": "read_cooler_chroms",
        "cool_path": str(source),
        "chromosome_count": len(rows),
        "chromosomes": rows,
        "provenance": _provenance(
            "cooler",
            "chromosome table of the cooler (cooler.Cooler.chroms, api.chroms)",
            order="file order is the heatmap order",
        ),
    }


__all__ = [
    "coarsen_contact_matrix",
    "fetch_chrom_sizes",
    "load_fasta_lengths",
    "merge_contact_matrices",
    "parse_genomic_region",
    "query_bed_intervals",
    "read_chrom_sizes",
    "read_cooler_chroms",
    "read_cooler_pixels",
    "rename_cooler_chromosomes",
    "zoomify_contact_matrix",
]

_register_public_verbose(globals(), tuple(__all__))
