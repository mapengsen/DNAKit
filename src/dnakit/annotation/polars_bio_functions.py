"""Polars-bio wrappers for interval, gene, variant, sequence, alignment, and contact-pair files.

Each public function reads local files with ``polars_bio``, runs one upstream operation, and
returns a JSON-serializable dictionary with ``status``, ``operation``, result fields, and
``provenance``. Coordinates follow the BED convention (0-based, half-open) unless a docstring
says otherwise. Rows are sorted by their values so repeated runs return identical lists.

Functions are importable as ``dnakit.annotation.polars_bio_functions.<name>``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, named_sequence_items, require, sequence_text
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_FASTQC_INTEGER_METRICS = frozenset({"n_seq", "total_bases", "min_len", "max_len"})


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int = 512) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _names(value: Sequence[str] | None, name: str) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise ConfigurationError(f"{name} must be a list of field names.")
    items = [_text(item, name, maximum=128) for item in value]
    if not items:
        raise ConfigurationError(f"{name} must not be empty when given.")
    return items


def _lengths(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value or len(value) > 100_000:
        raise ConfigurationError(
            f"{name} must be a non-empty mapping from chromosome name to length."
        )
    return {
        _text(chrom, f"{name} chromosome"): _integer(length, f"{name} length", 1, 10**12)
        for chrom, length in value.items()
    }


def _new_output(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output file already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _plain(value: Any) -> Any:
    """Convert numpy, polars, and missing values to JSON-compatible Python objects."""

    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
        return _plain(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _record_order(record: Mapping[str, Any]) -> tuple[tuple[int, Any], ...]:
    """Type-safe sort key so rows come back in one order regardless of hash seed."""

    keys: list[tuple[int, Any]] = []
    for value in record.values():
        if value is None:
            keys.append((3, ""))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            keys.append((0, float(value)))
        else:
            keys.append((2, str(value)))
    return tuple(keys)


def _ordered(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(records, key=_record_order)


def _table(frame: Any, limit: int, what: str) -> list[dict[str, Any]]:
    """Materialize a polars frame as sorted, JSON-compatible records within the row limit."""

    if hasattr(frame, "collect"):
        frame = frame.collect()
    if frame.height > limit:
        raise ConfigurationError(f"{what} has {frame.height} rows; max_rows={limit} is the limit.")
    return _ordered([{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()])


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("polars-bio"), "convention": convention}


_BED_CONVENTION = "0-based half-open coordinates (BED)"


# Interval arithmetic on BED files.


def overlap_intervals(
    left: str | Path,
    right: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return every pair of overlapping intervals between two BED files, using polars-bio."""

    left_path = input_path(left)
    right_path = input_path(right)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.overlap(
        pb.read_bed(str(left_path), use_zero_based=True),
        pb.read_bed(str(right_path), use_zero_based=True),
        suffixes=("_1", "_2"),
        output_type="polars.DataFrame",
    )
    pairs = _table(frame, limit, "Overlap result")
    return {
        "status": "completed",
        "operation": "overlap_intervals",
        "pair_count": len(pairs),
        "pairs": pairs,
        "provenance": _provenance(_BED_CONVENTION + "; columns suffixed _1 (left) and _2 (right)"),
    }


def count_interval_overlaps(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count how many intervals of ``other`` overlap each interval of ``source``."""

    source_path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.count_overlaps(
        pb.read_bed(str(source_path), use_zero_based=True),
        pb.read_bed(str(other_path), use_zero_based=True),
        output_type="polars.DataFrame",
    )
    intervals = _table(frame, limit, "Count result")
    return {
        "status": "completed",
        "operation": "count_interval_overlaps",
        "interval_count": len(intervals),
        "total_overlaps": sum(int(row["count"]) for row in intervals),
        "intervals": intervals,
        "provenance": _provenance(_BED_CONVENTION + "; one row per interval of source"),
    }


def cluster_intervals(
    source: str | Path,
    *,
    min_distance: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group overlapping or nearby intervals into clusters and number the clusters in order."""

    path = input_path(source)
    distance = _integer(min_distance, "min_distance", 0, 10**9)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.cluster(
        pb.read_bed(str(path), use_zero_based=True),
        min_dist=distance,
        output_type="polars.DataFrame",
    )
    rows = _table(frame, limit, "Cluster result")
    # Upstream cluster numbers depend on processing order; renumber by coordinate order.
    groups = sorted(
        {(row["cluster"], row["cluster_start"], row["cluster_end"], row["chrom"]) for row in rows},
        key=lambda item: (item[3], item[1], item[2], item[0]),
    )
    renumber = {group[0]: index for index, group in enumerate(groups)}
    for row in rows:
        row["cluster"] = renumber[row["cluster"]]
    rows = _ordered(rows)
    return {
        "status": "completed",
        "operation": "cluster_intervals",
        "min_distance": distance,
        "cluster_count": len(groups),
        "intervals": rows,
        "provenance": _provenance(_BED_CONVENTION + "; cluster ids numbered by coordinate order"),
    }


def complement_intervals(
    source: str | Path,
    *,
    chrom_sizes: Mapping[str, int],
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the gaps between intervals inside each chromosome of ``chrom_sizes``."""

    path = input_path(source)
    sizes = _lengths(chrom_sizes, "chrom_sizes")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    pl = require("polars", "bio-polarsbio")
    view = pl.DataFrame(
        {
            "chrom": list(sizes),
            "start": [0] * len(sizes),
            "end": [sizes[chrom] for chrom in sizes],
        }
    )
    frame = pb.complement(
        pb.read_bed(str(path), use_zero_based=True),
        view_df=view,
        output_type="polars.DataFrame",
    )
    gaps = _table(frame, limit, "Complement result")
    return {
        "status": "completed",
        "operation": "complement_intervals",
        "gap_count": len(gaps),
        "gaps": gaps,
        "provenance": _provenance(_BED_CONVENTION + "; gaps are inside the given chromosomes"),
    }


def subtract_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the bases covered by ``other`` from the intervals of ``source``."""

    source_path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.subtract(
        pb.read_bed(str(source_path), use_zero_based=True),
        pb.read_bed(str(other_path), use_zero_based=True),
        output_type="polars.DataFrame",
    )
    pieces = _table(frame, limit, "Subtract result")
    return {
        "status": "completed",
        "operation": "subtract_intervals",
        "piece_count": len(pieces),
        "pieces": pieces,
        "provenance": _provenance(_BED_CONVENTION + "; pieces keep the source columns"),
    }


# Gene and variant files.


def read_gene_annotation(
    path: str | Path,
    *,
    file_format: str = "gtf",
    attributes: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GTF or GFF3 gene annotation file and return one record per feature."""

    source = input_path(path)
    if file_format not in {"gtf", "gff"}:
        raise ConfigurationError("file_format must be 'gtf' or 'gff'.")
    fields = _names(attributes, "attributes")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    reader = pb.read_gtf if file_format == "gtf" else pb.read_gff
    frame = reader(str(source), attr_fields=fields, use_zero_based=True)
    features = _table(frame, limit, "Annotation")
    return {
        "status": "completed",
        "operation": "read_gene_annotation",
        "file_format": file_format,
        "feature_count": len(features),
        "features": features,
        "provenance": _provenance(_BED_CONVENTION + "; start is GTF/GFF start-1; end is exclusive"),
    }


def read_variant_calls(
    path: str | Path,
    *,
    file_format: str = "vcf",
    info_fields: Sequence[str] | None = None,
    format_fields: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF or BCF file into one record per variant with INFO and FORMAT fields."""

    source = input_path(path)
    if file_format not in {"vcf", "bcf"}:
        raise ConfigurationError("file_format must be 'vcf' or 'bcf'.")
    info = _names(info_fields, "info_fields")
    formats = _names(format_fields, "format_fields")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    reader = pb.read_vcf if file_format == "vcf" else pb.read_bcf
    frame = reader(
        str(source),
        info_fields=info,
        format_fields=formats,
        use_zero_based=True,
    )
    variants = _table(frame, limit, "Variant file")
    return {
        "status": "completed",
        "operation": "read_variant_calls",
        "file_format": file_format,
        "variant_count": len(variants),
        "variants": variants,
        "provenance": _provenance(_BED_CONVENTION + "; start is POS-1 and end is start+1 for SNVs"),
    }


def write_vcf_variants(
    source: str | Path,
    destination: str | Path,
    *,
    info_fields: Sequence[str] | None = None,
    format_fields: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a VCF file with polars-bio and write its records to a new VCF file."""

    source_path = input_path(source)
    target = _new_output(destination)
    info = _names(info_fields, "info_fields")
    formats = _names(format_fields, "format_fields")
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_vcf(
        str(source_path),
        info_fields=info,
        format_fields=formats,
        use_zero_based=True,
    )
    written = int(pb.write_vcf(frame, str(target)))
    return {
        "status": "completed",
        "operation": "write_vcf_variants",
        "records_written": written,
        "destination": str(target),
        "provenance": _provenance(_BED_CONVENTION + "; POS is written as start+1"),
    }


# Sequence files.


def read_fasta_sequences(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a FASTA file with polars-bio and return each sequence with its length."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_fasta(str(source))
    records = _table(frame, limit, "FASTA file")
    for record in records:
        record["length"] = len(record["sequence"])
    return {
        "status": "completed",
        "operation": "read_fasta_sequences",
        "sequence_count": len(records),
        "sequences": _ordered(records),
        "provenance": _provenance("sequence text is returned as stored (case preserved)"),
    }


def read_fastq_reads(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a FASTQ file with polars-bio and return each read with its quality string."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_fastq(str(source))
    reads = _table(frame, limit, "FASTQ file")
    for read in reads:
        read["length"] = len(read["sequence"])
    return {
        "status": "completed",
        "operation": "read_fastq_reads",
        "read_count": len(reads),
        "reads": _ordered(reads),
        "provenance": _provenance("Phred+33 quality strings are returned as stored"),
    }


def write_fasta_sequences(
    sequences: Any,
    destination: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write DNA sequences, given as ``DNA`` records or an ID-to-sequence mapping, to FASTA."""

    target = _new_output(destination)
    items = named_sequence_items(sequences)
    if not items:
        raise ConfigurationError("At least one sequence is required.")
    names = [_text(name, "sequence id", maximum=256) for name, _ in items]
    values = [sequence_text(sequence) for _, sequence in items]
    pb = require("polars_bio", "bio-polarsbio")
    pl = require("polars", "bio-polarsbio")
    frame = pl.DataFrame(
        {
            "name": names,
            "description": pl.Series("description", [None] * len(names), dtype=pl.String),
            "sequence": values,
        }
    )
    written = int(pb.write_fasta(frame, str(target)))
    return {
        "status": "completed",
        "operation": "write_fasta_sequences",
        "records_written": written,
        "destination": str(target),
        "provenance": _provenance("sequences are validated as DNA and written uppercase"),
    }


# Alignments, pileups, and read quality.


def read_alignment_records(
    path: str | Path,
    *,
    file_format: str = "sam",
    tag_fields: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read SAM or BAM alignments into one record per read with CIGAR and sequence."""

    source = input_path(path)
    if file_format not in {"sam", "bam"}:
        raise ConfigurationError("file_format must be 'sam' or 'bam'.")
    tags = _names(tag_fields, "tag_fields")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    reader = pb.read_sam if file_format == "sam" else pb.read_bam
    frame = reader(str(source), tag_fields=tags, use_zero_based=True)
    alignments = _table(frame, limit, "Alignment file")
    return {
        "status": "completed",
        "operation": "read_alignment_records",
        "file_format": file_format,
        "record_count": len(alignments),
        "alignments": alignments,
        "provenance": _provenance(_BED_CONVENTION + "; end is the exclusive reference end"),
    }


def alignment_depth(
    path: str | Path,
    *,
    min_mapping_quality: int = 0,
    filter_flag: int = 1796,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return runs of constant read depth from a SAM or BAM file (samtools-style filters)."""

    source = input_path(path)
    quality = _integer(min_mapping_quality, "min_mapping_quality", 0, 255)
    flag = _integer(filter_flag, "filter_flag", 0, 65535)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.depth(
        str(source),
        filter_flag=flag,
        min_mapping_quality=quality,
        use_zero_based=True,
        output_type="polars.DataFrame",
    )
    runs = _table(frame, limit, "Depth result")
    renamed = [
        {
            "chrom": row["contig"],
            "start": row["pos_start"],
            "end": row["pos_end"],
            "depth": row["coverage"],
        }
        for row in runs
    ]
    return {
        "status": "completed",
        "operation": "alignment_depth",
        "run_count": len(renamed),
        "runs": _ordered(renamed),
        "provenance": _provenance(
            _BED_CONVENTION + f"; filter_flag={flag}; min_mapping_quality={quality}"
        ),
    }


def fastqc_summary(
    path: str | Path,
    *,
    modules: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run FastQC-style read statistics over a FASTQ file in one streaming pass."""

    source = input_path(path)
    selected = _names(modules, "modules")
    pb = require("polars_bio", "bio-polarsbio")
    result = pb.fastqc(str(source), modules=selected, group=True)
    basic: dict[str, Any] = {}
    for row in result.basic_stats.collect().to_dicts():
        value = float(row["value"])
        metric = str(row["metric"])
        basic[metric] = int(value) if metric in _FASTQC_INTEGER_METRICS else round(value, 6)
    statuses = {
        str(row["module"]): str(row["status"]) for row in result.summary().collect().to_dicts()
    }
    per_gc = [
        {"gc_pct": int(row["gc_pct"]), "count": int(row["count"])}
        for row in result.per_seq_gc.collect().to_dicts()
        if int(row["count"]) > 0
    ]
    return {
        "status": "completed",
        "operation": "fastqc_summary",
        "basic_statistics": dict(sorted(basic.items())),
        "module_status": dict(sorted(statuses.items())),
        "per_sequence_gc": per_gc,
        "provenance": _provenance("GC percent excludes N bases; FastQC-style module statuses"),
    }


# Contact pairs, BED files, cooler matrices, and alignments of sequences.


def read_contact_pairs(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a 4DN Pairs file of Hi-C read pairs (chromosome and 0-based position per side)."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_pairs(str(source), use_zero_based=True)
    pairs = _table(frame, limit, "Pairs file")
    return {
        "status": "completed",
        "operation": "read_contact_pairs",
        "pair_count": len(pairs),
        "pairs": pairs,
        "provenance": _provenance("Pairs header-declared columns; positions converted to 0-based"),
    }


def read_bed_intervals(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a BED file into one record per interval with its optional name and score columns."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_bed(str(source), use_zero_based=True)
    intervals = _table(frame, limit, "BED file")
    return {
        "status": "completed",
        "operation": "read_bed_intervals",
        "interval_count": len(intervals),
        "intervals": intervals,
        "provenance": _provenance(_BED_CONVENTION),
    }


def read_cooler_matrix(
    path: str | Path,
    *,
    resolution: int | None = None,
    include_weights: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a Cooler (.cool or .mcool) Hi-C contact matrix as one record per non-zero bin pair."""

    source = input_path(path)
    if resolution is not None:
        resolution = _integer(resolution, "resolution", 1, 10**10)
    weights = _flag(include_weights, "include_weights")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_cool(
        str(source),
        resolution=resolution,
        include_weights=weights,
        use_zero_based=True,
    )
    contacts = _table(frame, limit, "Cooler file")
    return {
        "status": "completed",
        "operation": "read_cooler_matrix",
        "contact_count": len(contacts),
        "total_counts": sum(int(row["count"]) for row in contacts),
        "contacts": contacts,
        "provenance": _provenance(_BED_CONVENTION + "; one row per stored bin pair"),
    }


def read_stockholm_alignments(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a Stockholm multiple-sequence alignment, one record per aligned sequence."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_sto(str(source))
    rows = _table(frame, limit, "Stockholm file")
    lengths = {len(row["sequence"]) for row in rows}
    return {
        "status": "completed",
        "operation": "read_stockholm_alignments",
        "sequence_count": len(rows),
        "alignment_columns": lengths.pop() if len(lengths) == 1 else None,
        "alignments": rows,
        "provenance": _provenance("gap characters '-' and '.' are kept in the aligned sequences"),
    }


__all__ = [
    "alignment_depth",
    "cluster_intervals",
    "complement_intervals",
    "count_interval_overlaps",
    "fastqc_summary",
    "overlap_intervals",
    "read_alignment_records",
    "read_bed_intervals",
    "read_contact_pairs",
    "read_cooler_matrix",
    "read_fasta_sequences",
    "read_fastq_reads",
    "read_gene_annotation",
    "read_stockholm_alignments",
    "read_variant_calls",
    "subtract_intervals",
    "write_fasta_sequences",
    "write_vcf_variants",
]

_register_public_verbose(globals(), tuple(__all__))
