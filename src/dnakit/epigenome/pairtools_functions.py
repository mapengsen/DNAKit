"""Restriction-fragment, alignment, region, scaling, header and phasing functions from pairtools.

Each public function imports pairtools (and pysam where needed) only when it is called.
Missing dependencies raise ``BackendUnavailableError`` with the ``pip install "dnakit[...]"``
hint. Inputs are checked before upstream calls, and upstream ``ValueError``, ``KeyError``
and ``OSError`` are reported as ``ConfigurationError``.

Results are plain dictionaries with ``status``, ``operation``, the result fields and
``provenance``. Grouped results are sorted explicitly.
"""

from __future__ import annotations

import math
from array import array
from collections.abc import Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_PAIRTOOLS_EXTRA = "bio-pairtools"
_NUMPY_EXTRA = "bio-analysis"
_PANDAS_EXTRA = "bio-pairtools"
_MAX_CHROMOSOMES = 10_000
_MAX_CHROM_LENGTH = 2_000_000_000
_MAX_SITES = 1_000_000
_MAX_POSITIONS = 100_000
_MAX_REGIONS = 1_000_000
_MAX_RECORDS = 100_000
_MAX_PAIRS_ROWS = 1_000
_MAX_HEADER_COLUMNS = 64
_MAX_SEQUENCE = 100_000
_MAX_CONTACT_DIST = 1_000_000_000
_MAX_BINS_PER_DECADE = 50
_MAX_CURVE = 100_000
_VALID_BASES = frozenset("ACGTNacgtn")
_PHASE_VARIANTS = ("XA", "XB")
_SHAPES = ("upper triangle", "whole matrix")


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


def _existing(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a file path.")
    return input_path(value)


def _sequence(value: object, name: str, limit: int) -> list[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    if not 1 <= len(value) <= limit:
        raise ConfigurationError(f"{name} must contain 1 to {limit} items.")
    return list(value)


def _chromsizes(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= _MAX_CHROMOSOMES:
        raise ConfigurationError(
            f"{name} must map 1 to {_MAX_CHROMOSOMES} chromosome names to lengths."
        )
    return {
        _text(key, f"{name} chromosome"): _integer(length, f"length of {key}", 1, _MAX_CHROM_LENGTH)
        for key, length in value.items()
    }


def _plain(value: Any) -> Any:
    """Convert numpy scalars to Python values and non-finite floats to None."""

    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _provenance(distribution: str, method: str, **notes: str) -> dict[str, Any]:
    try:
        installed = version(distribution)
    except PackageNotFoundError:
        installed = "unknown"
    return {"package": f"{distribution} {installed}", "method": method, **notes}


def find_restriction_fragment(
    sites: Mapping[str, Sequence[int]],
    chrom: str,
    position: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the restriction fragment that contains a position, using pairtools."""

    if not isinstance(sites, Mapping) or not 1 <= len(sites) <= _MAX_CHROMOSOMES:
        raise ConfigurationError("sites must map chromosome names to restriction site positions.")
    name = _text(chrom, "chrom")
    if name not in sites:
        raise ConfigurationError(f"No restriction sites are listed for {name}.")
    values = _sequence(sites[name], f"sites[{name!r}]", _MAX_SITES)
    ints = [_integer(site, "site position", 0, _MAX_CHROM_LENGTH) for site in values]
    if len(ints) < 2 or ints != sorted(set(ints)):
        raise ConfigurationError(
            "Restriction sites must be at least two unique, ascending positions."
        )
    pos = _integer(position, "position", ints[0], ints[-1])
    numpy = require("numpy", _NUMPY_EXTRA)
    find_rfrag = require("pairtools.lib.restrict", _PAIRTOOLS_EXTRA).find_rfrag
    table = {name: numpy.asarray(ints, dtype=numpy.int64)}
    index, start, end = find_rfrag(table, name, pos)
    return {
        "status": "completed",
        "operation": "find_restriction_fragment",
        "chrom": name,
        "position": pos,
        "fragment_index": int(index),
        "start": int(start),
        "end": int(end),
        "length": int(end) - int(start),
        "provenance": _provenance(
            "pairtools",
            "binary search of sorted restriction sites (pairtools.lib.restrict.find_rfrag)",
            coordinates="fragments are half-open between consecutive sites",
        ),
    }


def compute_contact_scaling(
    pairs_path: str | Path,
    *,
    chromsizes: Mapping[str, int],
    min_distance: int = 1,
    max_distance: int = 1_000_000_000,
    bins_per_decade: int = 8,
    ignore_trans: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Contact frequency against genomic distance (scaling curve) of a .pairs file."""

    source = _existing(pairs_path, "pairs_path")
    sizes = _chromsizes(chromsizes, "chromsizes")
    low = _integer(min_distance, "min_distance", 1, _MAX_CONTACT_DIST)
    high = _integer(max_distance, "max_distance", low + 1, _MAX_CONTACT_DIST)
    per_decade = _integer(bins_per_decade, "bins_per_decade", 1, _MAX_BINS_PER_DECADE)
    skip_trans = _flag(ignore_trans, "ignore_trans")
    scaling = require("pairtools.lib.scaling", _PAIRTOOLS_EXTRA)
    try:
        cis, trans = scaling.compute_scaling(
            str(source),
            chromsizes=sizes,
            dist_range=(low, high),
            n_dist_bins_decade=per_decade,
            ignore_trans=skip_trans,
            nproc_in=1,
        )
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"pairtools could not compute the scaling curve: {exc}") from exc
    grouped = (
        cis.groupby(["min_dist", "max_dist"], as_index=False)[["n_pairs", "n_bp2"]]
        .sum()
        .sort_values(["min_dist", "max_dist"])
    )
    bins = []
    for row in grouped.to_dict(orient="records"):
        pairs_count = int(row["n_pairs"])
        possible = float(row["n_bp2"])
        bins.append(
            {
                "min_dist": int(row["min_dist"]),
                "max_dist": int(row["max_dist"]),
                "n_pairs": pairs_count,
                "n_bp2": _plain(possible),
                "contact_frequency": None if possible == 0 else _plain(pairs_count / possible),
            }
        )
    return {
        "status": "completed",
        "operation": "compute_contact_scaling",
        "pairs_path": str(source),
        "dist_range": [low, high],
        "bins_per_decade": per_decade,
        "ignore_trans": skip_trans,
        "cis_pairs": int(cis["n_pairs"].sum()),
        "trans_pairs": 0 if skip_trans else int(trans["n_pairs"].sum()),
        "bins": bins,
        "provenance": _provenance(
            "pairtools",
            "cis contacts per distance bin divided by possible pairs (pairtools.lib.scaling)",
            note="sums are taken over all regions and strand pairs",
        ),
    }


def normalize_scaling_curve(
    distances: Sequence[int],
    frequencies: Sequence[float],
    *,
    norm_window: Sequence[int] = (10_000, 1_000_000_000),
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a contact-frequency curve so its mean inside a distance window is 1."""

    dist = _sequence(distances, "distances", _MAX_CURVE)
    freq = _sequence(frequencies, "frequencies", _MAX_CURVE)
    if len(dist) != len(freq):
        raise ConfigurationError("distances and frequencies must have the same length.")
    bins = [_integer(value, "distance", 1, _MAX_CONTACT_DIST) for value in dist]
    if bins != sorted(set(bins)):
        raise ConfigurationError("distances must be unique and ascending.")
    values = [_number(value, "frequency", 0.0, 1e12) for value in freq]
    if isinstance(norm_window, (str, bytes)) or not isinstance(norm_window, Sequence):
        raise ConfigurationError("norm_window must be a pair of distances (low, high).")
    if len(norm_window) != 2:
        raise ConfigurationError("norm_window must be a pair of distances (low, high).")
    low = _integer(norm_window[0], "norm_window low", 1, _MAX_CONTACT_DIST)
    high = _integer(norm_window[1], "norm_window high", low, _MAX_CONTACT_DIST)
    scaling = require("pairtools.lib.scaling", _PAIRTOOLS_EXTRA)
    numpy = require("numpy", _NUMPY_EXTRA)
    try:
        factor = scaling.norm_scaling_factor(
            numpy.asarray(bins, dtype=float),
            numpy.asarray(values, dtype=float),
            (low, high),
        )
        normalized = scaling.norm_scaling(
            numpy.asarray(bins, dtype=float),
            numpy.asarray(values, dtype=float),
            (low, high),
        )
    except (ValueError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"The curve could not be normalised: {exc}") from exc
    return {
        "status": "completed",
        "operation": "normalize_scaling_curve",
        "norm_window": [low, high],
        "factor": _plain(float(factor)),
        "points": [
            {"distance": dist_value, "frequency": _plain(float(normalized_value))}
            for dist_value, normalized_value in zip(bins, normalized.tolist(), strict=True)
        ],
        "provenance": _provenance(
            "pairtools",
            "mean frequency in norm_window set to 1 (pairtools.lib.scaling.norm_scaling)",
        ),
    }


def assign_genomic_regions(
    positions: Sequence[Mapping[str, Any]],
    regions: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign genomic positions to the region that contains them, using pairtools."""

    points = _sequence(positions, "positions", _MAX_POSITIONS)
    spans = _sequence(regions, "regions", _MAX_REGIONS)
    chroms: list[str] = []
    coords: list[int] = []
    for item in points:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Each position must be a mapping with chrom and pos.")
        chroms.append(_text(item.get("chrom"), "chrom"))
        coords.append(_integer(item.get("pos"), "pos", 0, _MAX_CHROM_LENGTH))
    table = []
    for item in spans:
        if not isinstance(item, Mapping):
            raise ConfigurationError("Each region must be a mapping with chrom, start and end.")
        start = _integer(item.get("start"), "start", 0, _MAX_CHROM_LENGTH)
        end = _integer(item.get("end"), "end", start + 1, _MAX_CHROM_LENGTH)
        table.append({"chrom": _text(item.get("chrom"), "chrom"), "start": start, "end": end})
    pandas = require("pandas", _PANDAS_EXTRA)
    numpy = require("numpy", _NUMPY_EXTRA)
    scaling = require("pairtools.lib.scaling", _PAIRTOOLS_EXTRA)
    frame = pandas.DataFrame(table, columns=["chrom", "start", "end"])
    try:
        assigned = scaling.assign_regs(
            numpy.asarray(chroms, dtype=object),
            numpy.asarray(coords, dtype=numpy.int64),
            frame,
        )
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"pairtools could not assign the regions: {exc}") from exc
    rows = []
    for chrom, pos, found in zip(chroms, coords, assigned.tolist(), strict=True):
        index, start, end = (int(value) for value in found)
        rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "region_index": None if index < 0 else index,
                "region_start": None if index < 0 else start,
                "region_end": None if index < 0 else end,
            }
        )
    return {
        "status": "completed",
        "operation": "assign_genomic_regions",
        "position_count": len(rows),
        "assigned_count": sum(1 for row in rows if row["region_index"] is not None),
        "positions": rows,
        "provenance": _provenance(
            "pairtools",
            "region lookup by chromosome and start (pairtools.lib.scaling.assign_regs)",
            index="region_index counts regions of the same chromosome after sorting by start",
        ),
    }


def find_read_mismatches(
    sequence: str,
    qualities: Sequence[int],
    aligned_pairs: Sequence[Sequence[Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """List read bases that differ from the reference, given aligned pairs, using pairtools."""

    seq = _text(sequence, "sequence")
    if len(seq) > _MAX_SEQUENCE or set(seq) - _VALID_BASES:
        raise ConfigurationError("sequence must contain only A, C, G, T or N.")
    quals_in = _sequence(qualities, "qualities", _MAX_SEQUENCE)
    if len(quals_in) != len(seq):
        raise ConfigurationError("qualities must have one value per base.")
    quals = array("B", [_integer(value, "quality", 0, 255) for value in quals_in])
    pairs = _sequence(aligned_pairs, "aligned_pairs", _MAX_SEQUENCE)
    clean: list[tuple[int, int, str]] = []
    for item in pairs:
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 3:
            raise ConfigurationError("Each aligned pair must be (read_pos, ref_pos, ref_base).")
        read_pos = _integer(item[0], "read_pos", 0, len(seq) - 1)
        ref_pos = _integer(item[1], "ref_pos", 0, _MAX_CHROM_LENGTH)
        ref_base = _text(item[2], "ref_base")
        if len(ref_base) != 1 or ref_base.upper() not in _VALID_BASES:
            raise ConfigurationError("ref_base must be a single A, C, G, T or N.")
        clean.append((read_pos, ref_pos, ref_base))
    mismatch = require("pairtools.lib.parse_pysam", _PAIRTOOLS_EXTRA).get_mismatches_c
    try:
        found = mismatch(seq, quals, clean)
    except (ValueError, IndexError, TypeError) as exc:
        raise ConfigurationError(f"pairtools could not compare the read: {exc}") from exc
    rows = [
        {
            "ref_base": str(ref),
            "read_base": str(read),
            "quality": int(quality),
            "ref_pos": int(ref_pos),
            "read_pos": int(read_pos),
        }
        for ref, read, quality, ref_pos, read_pos in found
    ]
    rows.sort(key=lambda row: (row["read_pos"], row["ref_pos"]))
    return {
        "status": "completed",
        "operation": "find_read_mismatches",
        "read_length": len(seq),
        "mismatch_count": len(rows),
        "mismatches": rows,
        "provenance": _provenance(
            "pairtools",
            "base-by-base comparison of the read with the reference (pairtools.lib.parse_pysam)",
            order="sorted by read position, then reference position",
        ),
    }


def parse_sam_alignments(
    sam_path: str | Path,
    *,
    min_mapq: int = 30,
    max_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the first SAM or BAM records as pairtools alignments (strand, 5' and 3' ends)."""

    source = _existing(sam_path, "sam_path")
    quality = _integer(min_mapq, "min_mapq", 0, 255)
    limit = _integer(max_records, "max_records", 1, _MAX_RECORDS)
    parse_pysam = require("pairtools.lib.parse_pysam", _PAIRTOOLS_EXTRA)
    parse = require("pairtools.lib.parse", _PAIRTOOLS_EXTRA)
    try:
        reader = parse_pysam.AlignmentFilePairtoolized(str(source), "r")
    except (ValueError, OSError) as exc:
        raise ConfigurationError(f"Input is not a readable SAM or BAM file: {exc}") from exc
    try:
        if not reader.header or "SQ" not in reader.header.to_dict():
            raise ConfigurationError("The SAM input needs an @SQ header; use samtools view -h.")
        alignments: list[dict[str, Any]] = []
        for index, record in enumerate(reader):
            if index >= limit:
                break
            parsed = parse.parse_pysam_entry(record, quality, sam_tags=None, store_seq=False)
            alignments.append(
                {
                    "read_id": str(record.query_name),
                    "flag": int(record.flag),
                    **{str(key): _plain(value) for key, value in parsed.items()},
                }
            )
    except (ValueError, KeyError, OSError, AttributeError) as exc:
        raise ConfigurationError(f"pairtools could not parse the SAM records: {exc}") from exc
    finally:
        reader.close()
    return {
        "status": "completed",
        "operation": "parse_sam_alignments",
        "sam_path": str(source),
        "min_mapq": quality,
        "record_count": len(alignments),
        "alignments": alignments,
        "provenance": _provenance(
            "pairtools",
            "strand, 5' and 3' ends, CIGAR spans and uniqueness (pairtools.lib.parse)",
            order="records are returned in file order",
        ),
    }


def read_pairs_table(
    pairs_path: str | Path,
    *,
    limit: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a .pairs file: header, chromosome sizes, row count and first rows."""

    source = _existing(pairs_path, "pairs_path")
    rows_limit = _integer(limit, "limit", 1, _MAX_PAIRS_ROWS)
    pairsio = require("pairtools.lib.pairsio", _PAIRTOOLS_EXTRA)
    try:
        frame, header, chromsizes = pairsio.read_pairs(str(source), nproc=1)
    except (ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"Input is not a readable pairs file: {source} ({exc})") from exc
    head = frame.head(rows_limit)
    rows = [
        {str(key): _plain(value) for key, value in row.items()}
        for row in head.to_dict(orient="records")
    ]
    return {
        "status": "completed",
        "operation": "read_pairs_table",
        "pairs_path": str(source),
        "row_count": len(frame),
        "columns": [str(column) for column in frame.columns],
        "chromsizes": {str(name): int(size) for name, size in sorted(chromsizes.items())},
        "header": [str(line) for line in header],
        "rows": rows,
        "provenance": _provenance(
            "pairtools",
            "header and body of a .pairs file read by pairtools.lib.pairsio.read_pairs",
            note=f"at most {rows_limit} rows are returned, in file order",
        ),
    }


def make_pairs_header(
    *,
    chromsizes: Mapping[str, int],
    assembly: str = "unknown",
    columns: Sequence[str] | None = None,
    shape: str = "upper triangle",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the standard .pairs header lines for a genome, using pairtools."""

    sizes = _chromsizes(chromsizes, "chromsizes")
    label = _text(assembly, "assembly")
    if shape not in _SHAPES:
        raise ConfigurationError(f"shape must be one of {list(_SHAPES)}.")
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    kwargs: dict[str, Any] = {}
    if columns is not None:
        names = _sequence(columns, "columns", _MAX_HEADER_COLUMNS)
        kwargs["columns"] = [_text(name, "column") for name in names]
    try:
        lines = headerops.make_standard_pairsheader(
            assembly=label,
            chromsizes=dict(sizes),
            shape=shape,
            **kwargs,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"pairtools could not build the header: {exc}") from exc
    return {
        "status": "completed",
        "operation": "make_pairs_header",
        "assembly": label,
        "shape": shape,
        "chromosome_count": len(sizes),
        "lines": [str(line) for line in lines],
        "provenance": _provenance(
            "pairtools",
            "pairs format header with assembly, chromsize and columns lines "
            "(pairtools.lib.headerops.make_standard_pairsheader)",
        ),
    }


def subset_pairs_header(
    header_lines: Sequence[str],
    chroms: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep only the chromosome entries of a .pairs header that are listed, using pairtools."""

    lines = _sequence(header_lines, "header_lines", 100_000)
    text_lines = [_text(line, "header line") for line in lines]
    names = _sequence(chroms, "chroms", _MAX_CHROMOSOMES)
    subset = [_text(name, "chroms") for name in names]
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    try:
        kept = headerops.subset_chroms_in_pairsheader(text_lines, subset)
    except (ValueError, IndexError) as exc:
        raise ConfigurationError(f"pairtools could not subset the header: {exc}") from exc
    return {
        "status": "completed",
        "operation": "subset_pairs_header",
        "chroms": sorted(set(subset)),
        "input_line_count": len(text_lines),
        "output_line_count": len(kept),
        "lines": [str(line) for line in kept],
        "provenance": _provenance(
            "pairtools",
            "chromsize and chromosomes lines filtered to the subset "
            "(pairtools.lib.headerops.subset_chroms_in_pairsheader)",
        ),
    }


def read_chrom_order(
    chroms_path: str | Path,
    *,
    sam_chroms: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number the chromosomes of a two-column chromsizes file in file order, using pairtools."""

    source = _existing(chroms_path, "chroms_path")
    extra: list[str] | None = None
    if sam_chroms is not None:
        extra = [_text(name, "sam_chroms") for name in _sequence(sam_chroms, "sam_chroms", 100_000)]
    headerops = require("pairtools.lib.headerops", _PAIRTOOLS_EXTRA)
    try:
        order = headerops.get_chrom_order(str(source), extra)
    except (ValueError, IndexError, OSError) as exc:
        raise ConfigurationError(f"Chromosome file is not readable: {source} ({exc})") from exc
    rows = [{"chrom": name, "order": int(index)} for name, index in order.items()]
    rows.sort(key=lambda row: (row["order"], row["chrom"]))
    return {
        "status": "completed",
        "operation": "read_chrom_order",
        "chroms_path": str(source),
        "chromosome_count": len(rows),
        "chromosomes": rows,
        "provenance": _provenance(
            "pairtools",
            "1-based chromosome enumeration from the chromsizes file "
            "(pairtools.lib.headerops.get_chrom_order)",
            note="chromosomes present only in sam_chroms are appended in sorted order",
        ),
    }


def assign_haplotype_phase(
    chrom: str,
    *,
    tag_value: str,
    primary_score: int,
    secondary_score: int,
    variant: str = "XB",
    mismatches: int | None = None,
    phase_suffixes: Sequence[str] = ("_1", "_2"),
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign a read side to a haplotype from its XA or XB alternative-hit tag, using pairtools."""

    name = _text(chrom, "chrom")
    if variant not in _PHASE_VARIANTS:
        raise ConfigurationError("variant must be 'XA' or 'XB'.")
    tag = tag_value if isinstance(tag_value, str) else None
    if tag is None:
        raise ConfigurationError("tag_value must be the XA or XB tag text.")
    primary = _integer(primary_score, "primary_score", -1_000_000, 1_000_000)
    secondary = _integer(secondary_score, "secondary_score", -1_000_000, 1_000_000)
    if (
        isinstance(phase_suffixes, (str, bytes))
        or not isinstance(phase_suffixes, Sequence)
        or len(phase_suffixes) != 2
    ):
        raise ConfigurationError("phase_suffixes must be two suffixes, for example ('_1', '_2').")
    suffixes = (
        _text(phase_suffixes[0], "phase suffix"),
        _text(phase_suffixes[1], "phase suffix"),
    )
    if suffixes[0] == suffixes[1]:
        raise ConfigurationError("The two phase suffixes must differ.")
    phase_module = require("pairtools.lib.phase", _PAIRTOOLS_EXTRA)
    try:
        if variant == "XB":
            result = phase_module.phase_side_XB(name, tag, primary, secondary, suffixes)
        else:
            if mismatches is None:
                raise ConfigurationError("mismatches (NM) is required for the XA variant.")
            nm = _integer(mismatches, "mismatches", 0, 1_000_000)
            result = phase_module.phase_side_XA(name, tag, primary, secondary, nm, suffixes)
    except (ValueError, IndexError) as exc:
        raise ConfigurationError(f"The {variant} tag could not be parsed: {exc}") from exc
    phase, base, score_one, score_two, score_three = result
    return {
        "status": "completed",
        "operation": "assign_haplotype_phase",
        "chrom": name,
        "variant": variant,
        "phase": str(phase),
        "chromosome_base": str(base),
        "scores": [int(score_one), int(score_two), int(score_three)],
        "phase_suffixes": list(suffixes),
        "provenance": _provenance(
            "pairtools",
            "haplotype from the alternative-hit tag and the primary and secondary scores "
            "(pairtools.lib.phase)",
            codes="phase '0' or '1' from a suffix, '.' for homologous hits, '!' when ambiguous",
        ),
    }


__all__ = [
    "assign_genomic_regions",
    "assign_haplotype_phase",
    "compute_contact_scaling",
    "find_read_mismatches",
    "find_restriction_fragment",
    "make_pairs_header",
    "normalize_scaling_curve",
    "parse_sam_alignments",
    "read_chrom_order",
    "read_pairs_table",
    "subset_pairs_header",
]

_register_public_verbose(globals(), tuple(__all__))
