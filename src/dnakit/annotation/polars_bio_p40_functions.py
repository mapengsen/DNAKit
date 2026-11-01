"""polars-bio readers for alignment-style MSA files and BigWig tracks (pack P40).

polars-bio's A2M/A3M and BigWig readers as DNAKit functions. The VCF, BCF, GFF, GTF, BAM, SAM readers, the BAM
writer, ``nearest`` and ``coverage`` of the same package are already available as ``read_variant_calls``,
``read_gene_annotation``, ``read_alignment_records``, ``write_alignment_file``,
``nearest_polars_bio_intervals`` and ``coverage_polars_bio_intervals``. Files are explicit local paths;
polars-bio is imported only when a function runs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.polars_bio_functions import _DEFAULT_ROWS, _plain, _provenance, _rows_limit
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MSA_FORMATS = ("a2m", "a3m")


def _rows(frame: Any, limit: int, what: str) -> list[dict[str, Any]]:
    if hasattr(frame, "collect"):
        frame = frame.collect()
    if frame.height > limit:
        raise ConfigurationError(f"{what} has {frame.height} rows; max_rows={limit} is the limit.")
    return [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dicts()]


def read_msa_records(
    path: str | Path,
    *,
    file_format: str = "a2m",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an A2M or A3M alignment file with polars-bio (read_a2m, read_a3m).

    ``file_format`` is ``"a2m"`` or ``"a3m"``. One record per sequence in file order, with name, description
    and the aligned sequence exactly as written (gap and insertion characters are kept).
    """

    source = input_path(path)
    kind = _h.choice(file_format, _MSA_FORMATS, "file_format")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_a2m(str(source)) if kind == "a2m" else pb.read_a3m(str(source))
    records = _rows(frame, limit, f"{kind.upper()} file")
    return {
        "status": "completed",
        "operation": "read_msa_records",
        "file_format": kind,
        "record_count": len(records),
        "records": records,
        "provenance": _provenance("sequences are returned as written in the file"),
    }


def read_bigwig_intervals(
    path: str | Path,
    *,
    zero_based: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the intervals of a BigWig signal track with polars-bio (read_bigwig).

    Returns chrom, start, end and value (float32) for every stored interval. With ``zero_based=True`` (the
    DNAKit default) coordinates are 0-based half-open as in BED; ``False`` gives polars-bio's 1-based
    closed coordinates.
    """

    source = input_path(path)
    if not isinstance(zero_based, bool):
        raise ConfigurationError("zero_based must be True or False.")
    limit = _rows_limit(max_rows)
    pb = require("polars_bio", "bio-polarsbio")
    frame = pb.read_bigwig(str(source), use_zero_based=zero_based)
    records = _rows(frame, limit, "BigWig file")
    return {
        "status": "completed",
        "operation": "read_bigwig_intervals",
        "interval_count": len(records),
        "intervals": records,
        "provenance": _provenance(
            "0-based half-open coordinates (BED)" if zero_based else "1-based closed coordinates"),
    }


__all__ = ["read_bigwig_intervals", "read_msa_records"]

_register_public_verbose(globals(), tuple(__all__))
