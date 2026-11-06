"""BAM alignments as genomic intervals with pysam, with the columns and filters of pyranges1.

The pyranges1.read_bam function is the reference for the columns and filters.

``pyranges1.read_bam`` needs the separate bamread package, which the DNAKit environment does not
install (upstream then exits with status 1). ``read_pyranges1_bam`` reads the same alignments with
pysam (``dnakit[bio-reads]``) and builds a ``pyranges1.PyRanges`` (``dnakit[bio-pyranges1]``) with
the upstream columns: Chromosome, Start, End, Strand and Flag, plus QueryStart, QueryEnd,
QuerySequence, Name, Cigar and Quality when ``sparse`` is false. Coordinates are 0-based, End
exclusive, and the end comes from the CIGAR reference span. Results are JSON-compatible records.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_SPARSE_COLUMNS = ["Chromosome", "Start", "End", "Strand", "Flag"]
_FULL_COLUMNS = ["QueryStart", "QueryEnd", "QuerySequence", "Name", "Cigar", "Quality"]
_REVERSE_FLAG = 16
_DEFAULT_ROWS = 1000
_MAX_ROWS = 1_000_000
_CONVENTION = (
    "Start is the 0-based leftmost reference position of the alignment and End is exclusive, "
    "taken from the "
    "CIGAR reference span (M, D, N, = and X). Strand is '-' for flag 16. Flag is the SAM flag. "
    "Unmapped reads "
    "and reads without a reference span are not reported. QueryStart and QueryEnd (full columns "
    "only) are "
    "pysam's query_alignment_start and query_alignment_end: 0-based, soft clips excluded, end "
    "exclusive."
)


def _existing_bam(path: str | Path) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"path must be an existing BAM file: {resolved}")
    return resolved


def _native(value: Any) -> Any:
    if isinstance(value, float) and value != value:  # NaN is not valid JSON
        return None
    return value.item() if hasattr(value, "item") else value


def _alignment_rows(
    source: Path,
    *,
    mapq: int,
    required_flag: int,
    filter_flag: int,
    sparse: bool,
) -> list[dict[str, Any]]:
    pysam = require("pysam", "bio-reads")
    rows: list[dict[str, Any]] = []
    try:
        bam = pysam.AlignmentFile(str(source), "rb")
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"Could not open {source} as a BAM file: {exc}") from exc
    with bam:
        for read in bam.fetch(until_eof=True):
            flag = int(read.flag)
            if int(read.mapping_quality) < mapq:
                continue
            if flag & required_flag != required_flag:
                continue
            if flag & filter_flag:
                continue
            if read.reference_id < 0 or read.reference_end is None:
                continue
            row: dict[str, Any] = {
                "Chromosome": read.reference_name,
                "Start": int(read.reference_start),
                "End": int(read.reference_end),
                "Strand": "-" if flag & _REVERSE_FLAG else "+",
                "Flag": flag,
            }
            if not sparse:
                row.update(
                    {
                        "QueryStart": int(read.query_alignment_start),
                        "QueryEnd": int(read.query_alignment_end),
                        "QuerySequence": read.query_sequence,
                        "Name": read.query_name,
                        "Cigar": read.cigarstring,
                        "Quality": read.qual,
                    }
                )
            rows.append(row)
    return rows


def read_pyranges1_bam(
    path: str | Path,
    *,
    mapq: int = 0,
    required_flag: int = 0,
    filter_flag: int = 1540,
    sparse: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read BAM alignments as pyranges1 intervals (Chromosome, Start, End, Strand, Flag) with pysam.

    ``mapq`` keeps reads with mapping quality at least this value. ``required_flag`` keeps reads
    that have all of these SAM flag bits; ``filter_flag`` drops reads that have any of them (default
    1540: unmapped, QC fail, duplicate). ``sparse=False`` adds QueryStart, QueryEnd, QuerySequence,
    Name, Cigar and Quality. The whole file is read; the result holds the first ``max_rows``
    intervals with the total count and a ``truncated`` flag.
    """
    source = _existing_bam(path)
    quality = _cli.integer(mapq, "mapq", minimum=0, maximum=255)
    required = _cli.integer(required_flag, "required_flag", minimum=0, maximum=65535)
    filtered = _cli.integer(filter_flag, "filter_flag", minimum=0, maximum=65535)
    if not isinstance(sparse, bool):
        raise ConfigurationError("sparse must be a boolean.")
    limit = _cli.integer(max_rows, "max_rows", minimum=1, maximum=_MAX_ROWS)

    rows = _alignment_rows(
        source,
        mapq=quality,
        required_flag=required,
        filter_flag=filtered,
        sparse=sparse,
    )
    pd = require("pandas", "bio-pyranges1")  # lazy: pandas is not a core dependency
    columns = _SPARSE_COLUMNS if sparse else [*_SPARSE_COLUMNS, *_FULL_COLUMNS]
    frame = pd.DataFrame.from_records(rows, columns=columns)
    frame = frame.astype({"Start": "int64", "End": "int64", "Flag": "uint16"})
    frame["Chromosome"] = frame["Chromosome"].astype("category")
    frame["Strand"] = frame["Strand"].astype("category")

    pyranges1 = require("pyranges1", "bio-pyranges1")
    intervals = pyranges1.PyRanges(frame)
    head = pd.DataFrame(intervals).head(limit).astype(object)
    records = [
        {key: _native(value) for key, value in record.items()} for record in head.to_dict("records")
    ]
    return {
        "status": "completed",
        "operation": "read_pyranges1_bam",
        "provenance": {
            "package": f"pyranges1 {getattr(pyranges1, '__version__', 'unknown')}",
            "backend": f"pysam {require('pysam', 'bio-reads').__version__}",
            "convention": _CONVENTION,
        },
        "columns": list(columns),
        "record_count": len(frame),
        "truncated": bool(len(frame) > limit),
        "records": records,
    }


__all__ = ["read_pyranges1_bam"]
_register_public_verbose(globals())
