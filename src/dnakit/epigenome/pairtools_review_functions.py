"""SAM-to-pairs conversion, coverage filtering, duplicate flags and SAM header sizes from pairtools.

Each public function imports pairtools (and pysam, through pairtools) only when it is called.
Missing dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Inputs are checked before upstream calls; upstream ``ValueError``, ``KeyError`` and
``OSError`` are reported as ``ConfigurationError``.

Functions that write a file create a new path and never overwrite one.
"""

from __future__ import annotations

from collections.abc import Iterator
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_PAIRTOOLS_EXTRA = "bio-pairtools"
_MAX_CHROMOSOMES = 10_000
_MAX_POSITION = 2_000_000_000
_MAX_PAIRS_LINES = 5_000_000
_MAX_SAM_LINE = 1_000_000
_WALKS_POLICIES = ("mask", "5any", "5unique", "3any", "3unique", "all")
_COVERAGE_METHODS = ("sum", "max")
_PAIR_COLUMNS = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
]


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _choice(value: object, name: str, options: tuple[str, ...]) -> str:
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of: {', '.join(options)}.")
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


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def _sam_reader(source: Path) -> Any:
    parse_pysam = require("pairtools.lib.parse_pysam", _PAIRTOOLS_EXTRA)
    try:
        reader = parse_pysam.AlignmentFilePairtoolized(str(source), "r")
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Input is not a readable SAM or BAM file: {exc}") from exc
    if not reader.header or "SQ" not in reader.header.to_dict():
        reader.close()
        raise ConfigurationError("The SAM input needs an @SQ header; use samtools view -h.")
    return reader


def _count_pairs(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def sam_to_pairs(
    sam_path: str | Path,
    output_path: str | Path,
    *,
    min_mapq: int = 40,
    max_molecule_size: int = 500,
    max_inter_align_gap: int = 20,
    walks_policy: str = "mask",
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify SAM alignments of each read into Hi-C pairs and write a .pairs file (pairtools)."""

    source = _existing(sam_path, "sam_path")
    target = _new(output_path, "output_path")
    quality = _integer(min_mapq, "min_mapq", 0, 255)
    molecule = _integer(max_molecule_size, "max_molecule_size", 1, _MAX_POSITION)
    gap = _integer(max_inter_align_gap, "max_inter_align_gap", 0, _MAX_POSITION)
    policy = _choice(walks_policy, "walks_policy", _WALKS_POLICIES)
    parse = require("pairtools.lib.parse", _PAIRTOOLS_EXTRA)
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    reader = _sam_reader(source)
    try:
        sizes = headerops.get_chromsizes_from_pysam_header(reader.header)
        chromosomes = list(sizes.keys())
        if len(chromosomes) > _MAX_CHROMOSOMES:
            raise ConfigurationError(
                f"The SAM header lists more than {_MAX_CHROMOSOMES} sequences."
            )
        header = headerops.make_standard_pairsheader(
            assembly=None,
            chromsizes=[(chrom, sizes[chrom]) for chrom in chromosomes],
            columns=_PAIR_COLUMNS,
            shape="upper triangle",
        )
        kwargs = {
            "min_mapq": quality,
            "max_molecule_size": molecule,
            "max_inter_align_gap": gap,
            "walks_policy": policy,
            "drop_readid": False,
            "drop_seq": True,
            "drop_sam": True,
            "add_pair_index": False,
            "add_columns": "",
            "parse2": False,
            "flip": False,
            "readid_transform": None,
        }
        try:
            with target.open("w", encoding="utf-8") as handle:
                handle.writelines(line + "\n" for line in header)
                parse.streaming_classify(iter(reader), handle, chromosomes, None, None, **kwargs)
        except (ValueError, KeyError, OSError) as exc:
            target.unlink(missing_ok=True)
            raise ConfigurationError(f"pairtools could not classify the SAM file: {exc}") from exc
    finally:
        reader.close()
    pairs = _count_pairs(target)
    return {
        "status": "completed",
        "operation": "sam_to_pairs",
        "sam_path": str(source),
        "output_path": str(target),
        "chromosome_count": len(chromosomes),
        "pair_count": pairs,
        "columns": list(_PAIR_COLUMNS),
        "min_mapq": quality,
        "walks_policy": policy,
        "provenance": _provenance(
            "pairtools",
            "regular Hi-C read parser over SAM records (pairtools.lib.parse.streaming_classify)",
            output="unique alignments only; SAM text is dropped; header lists chromosome sizes",
        ),
    }


def _body_lines(path: Path) -> tuple[list[str], Iterator[str]]:
    handle = path.open(encoding="utf-8")
    header: list[str] = []
    for line in handle:
        if line.startswith("#"):
            header.append(line.rstrip("\n"))
        else:
            first = line
            break
    else:
        handle.close()
        raise ConfigurationError(f"The pairs file has no pair lines: {path}")

    def body() -> Iterator[str]:
        try:
            yield first
            yield from handle
        finally:
            handle.close()

    return header, body()


def filter_pairs_by_coverage(
    pairs_path: str | Path,
    output_path: str | Path,
    *,
    method: str = "max",
    max_dist: int = 500,
    max_cov: int = 8,
    removed_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove pairs whose fragment ends are covered too often by nearby pairs, using pairtools."""

    source = _existing(pairs_path, "pairs_path")
    target = _new(output_path, "output_path")
    kind = _choice(method, "method", _COVERAGE_METHODS)
    distance = _integer(max_dist, "max_dist", 0, _MAX_POSITION)
    coverage = _integer(max_cov, "max_cov", 0, 1_000_000)
    removed = None if removed_path is None else _new(removed_path, "removed_path")
    filterbycov = require("pairtools.lib.filterbycov", _PAIRTOOLS_EXTRA)
    header, body = _body_lines(source)
    columns = None
    counted = 0
    for line in header:
        if line.startswith("#columns:"):
            columns = line.split()[1:]
    if columns is not None and len(columns) < 7:
        body.close()
        raise ConfigurationError("The pairs file needs the eight standard pair columns.")
    try:
        with (
            target.open("w", encoding="utf-8") as kept,
            removed.open("w", encoding="utf-8") if removed is not None else _NullFile() as dropped,
        ):
            kept.writelines(line + "\n" for line in header)
            filterbycov.streaming_filterbycov(
                kind,
                distance,
                coverage,
                "\t",
                1,
                3,
                2,
                4,
                5,
                6,
                "!",
                _lines_without_newline(body),
                kept,
                dropped if removed is not None else None,
                None,
                None,
                False,
            )
    except (ValueError, KeyError, OSError, IndexError) as exc:
        target.unlink(missing_ok=True)
        if removed is not None:
            removed.unlink(missing_ok=True)
        raise ConfigurationError(f"pairtools could not filter the pairs: {exc}") from exc
    finally:
        body.close()
    kept_count = _count_pairs(target)
    counted = kept_count
    return {
        "status": "completed",
        "operation": "filter_pairs_by_coverage",
        "pairs_path": str(source),
        "output_path": str(target),
        "removed_path": None if removed is None else str(removed),
        "method": kind,
        "max_dist": distance,
        "max_cov": coverage,
        "kept_count": counted,
        "removed_count": None if removed is None else _count_pairs(removed),
        "provenance": _provenance(
            "pairtools",
            "pairs whose fragment ends have more than max_cov partners within max_dist are "
            "removed (pairtools.lib.filterbycov)",
            coverage_rule="sum adds partner counts; max uses the largest partner count per end",
        ),
    }


class _NullFile:
    def __enter__(self) -> _NullFile:
        return self

    def __exit__(self, *exc_info: object) -> None:
        return None


def _lines_without_newline(lines: Iterator[str]) -> Iterator[str]:
    for line in lines:
        yield line.rstrip("\n")


def mark_sam_duplicate(sam_line: str, *, verbose: bool = True) -> dict[str, Any]:
    """Set the duplicate flag (1024) of one SAM record and mark its pair type DD (pairtools)."""

    dedup = require("pairtools.lib.dedup", _PAIRTOOLS_EXTRA)
    pairsam_format = require("pairtools.lib.pairsam_format", _PAIRTOOLS_EXTRA)
    if not isinstance(sam_line, str) or not sam_line.strip() or len(sam_line) > _MAX_SAM_LINE:
        raise ConfigurationError("sam_line must be one SAM record as text.")
    if "\n" in sam_line.rstrip("\n"):
        raise ConfigurationError("sam_line must contain a single SAM record.")
    record = sam_line.rstrip("\n")
    columns = record.split("\t")
    if len(columns) < 11:
        raise ConfigurationError("A SAM record needs at least 11 tab-separated columns.")
    try:
        before = int(columns[1])
    except ValueError as exc:
        raise ConfigurationError("The SAM FLAG column must be an integer.") from exc
    separator = pairsam_format.SAM_SEP
    if separator in record:
        raise ConfigurationError("sam_line must not contain the pairsam field separator.")
    marked = dedup.mark_sam_as_dup(record.replace("\t", separator)).replace(separator, "\t")
    after = int(marked.split("\t")[1])
    return {
        "status": "completed",
        "operation": "mark_sam_duplicate",
        "sam": marked,
        "flag_before": before,
        "flag_after": after,
        "duplicate_flag_set": bool(after & 1024),
        "provenance": _provenance(
            "pairtools",
            "FLAG bit 1024 is set and optional Yt:Z tags are changed to DD "
            "(pairtools.lib.dedup.mark_sam_as_dup)",
        ),
    }


def chromsizes_from_sam_header(sam_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Read sequence names and lengths from the @SQ lines of a SAM or BAM header (pairtools)."""

    source = _existing(sam_path, "sam_path")
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    reader = _sam_reader(source)
    try:
        sizes = headerops.get_chromsizes_from_pysam_header(reader.header)
    finally:
        reader.close()
    rows = [{"name": str(name), "length": int(length)} for name, length in sizes.items()]
    if len(rows) > _MAX_CHROMOSOMES:
        raise ConfigurationError(f"The header lists more than {_MAX_CHROMOSOMES} sequences.")
    return {
        "status": "completed",
        "operation": "chromsizes_from_sam_header",
        "sam_path": str(source),
        "chromosome_count": len(rows),
        "chromosomes": rows,
        "provenance": _provenance(
            "pairtools",
            "@SQ names and lengths in header order "
            "(pairtools.lib.headerops.get_chromsizes_from_pysam_header)",
            order="header order is kept",
        ),
    }


__all__ = [
    "chromsizes_from_sam_header",
    "filter_pairs_by_coverage",
    "mark_sam_duplicate",
    "sam_to_pairs",
]

_register_public_verbose(globals(), tuple(__all__))
