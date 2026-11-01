"""Further bioframe operations as DNAKit functions (pack P40): UCSC gap, cytoband and mRNA tables, the UCSC
region-string completeness test and lazy FASTA loading.

Network wrappers call bioframe's UCSC client and record the service, the request, a UTC retrieval time and the
data; ordinary tests never use the network. bioframe is imported only when a function runs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.bioinformatics.bioframe_gap_functions import (
    _assembly,
    _bioframe_version,
    _optional_output,
    _records,
    _utc_now,
    _wrap_fetch_error,
)
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-bioframe"
_TABLES = ("gap", "cytoband", "mrna", "centromeres")
_ENGINES = ("pyfaidx", "pysam")
_URLS = {"gap": "database/gap.txt.gz", "cytoband": "database/cytoBand.txt.gz",
         "mrna": "database/all_mrna.txt.gz", "centromeres": "database/centromeres.txt.gz"}


def bioframe_fetch_ucsc_table(
    db: str,
    *,
    table: str = "gap",
    ideo: bool = False,
    max_rows: int = 10_000,
    output_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch a UCSC assembly table through bioframe's UCSCClient.

    ``table`` is ``"centromeres"`` (UCSCClient.fetch_centromeres), ``"gap"`` (assembly gaps; UCSCClient.fetch_gaps), ``"cytoband"`` (cytogenetic bands;
    fetch_cytoband, with ``ideo=True`` for the ideogram variant cytoBandIdeo) or ``"mrna"`` (all GenBank mRNA
    alignments; fetch_mrna, a very large table). ``max_rows`` bounds the returned records (and, for gap and
    mrna without ``output_path``, the rows read from the server); ``output_path``
    (new file) receives the full table as TSV.
    """

    assembly = _assembly(db)
    kind = _h.choice(table, _TABLES, "table")
    limit = _h.integer(max_rows, "max_rows", minimum=1, maximum=_h.MAX_ROWS * 10)
    target = _optional_output(output_path, "output_path")
    resources = require("bioframe.io.resources", _EXTRA)
    client = resources.UCSCClient(assembly)
    kwargs: dict[str, Any] = {} if target is not None or kind == "cytoband" else {"nrows": limit}
    try:
        if kind == "centromeres":
            frame = client.fetch_centromeres(**kwargs)
        elif kind == "gap":
            frame = client.fetch_gaps(**kwargs)
        elif kind == "cytoband":
            frame = client.fetch_cytoband(ideo=bool(ideo))
        else:
            frame = client.fetch_mrna(**kwargs)
    except (OSError, ValueError, KeyError) as exc:
        raise _wrap_fetch_error(exc, f"fetch the UCSC {kind} table for {assembly}") from exc
    if target is not None:
        frame.to_csv(target, sep="\t", index=False)
    shown = frame.head(limit)
    return {
        "status": "completed",
        "service": "UCSC Genome Browser database tables (hgdownload.soe.ucsc.edu)",
        "request": {"db": assembly, "table": kind, "ideo": bool(ideo), "file": _URLS[kind]},
        "retrieved_at": _utc_now(),
        "bioframe_version": _bioframe_version(),
        "row_count": len(frame),
        "data": _records(shown),
        "truncated": len(frame) > limit,
        "output_path": None if target is None else str(target),
    }


def bioframe_is_complete_ucsc_string(region: str, *, verbose: bool = True) -> dict[str, Any]:
    """Test whether a UCSC-style region string names a chromosome with a known end (bioframe.is_complete_ucsc_string).

    ``chr1:100-200`` and ``chr1:100-`` style strings: True only when chrom, start and end can all be read;
    ``chr1`` or ``chr1:100-`` are not complete.
    """

    text = _h.text(region, "region")
    bioframe = require("bioframe", _EXTRA)
    try:
        complete = bool(bioframe.is_complete_ucsc_string(text))
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"region is not a UCSC-style string: {exc}") from exc
    return {"status": "completed", "backend": "bioframe", "version": _bioframe_version(), "region": text,
            "is_complete": complete}


def bioframe_load_fasta(
    fasta: str | Path,
    *,
    engine: str = "pyfaidx",
    regions: list[str] | None = None,
    max_sequence_chars: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Open an indexed FASTA lazily with bioframe.load_fasta and report its records.

    ``engine`` is ``pyfaidx`` (default; the ``.fai`` index is created if missing, next to the FASTA) or
    ``pysam``. Returns the name and length of every record and, for ``regions`` given as UCSC strings
    (``chr1:101-200``), the sequence of each region (bounded by ``max_sequence_chars``). As in bioframe, the
    parsed numbers are used directly as the 0-based half-open slice ``[start:end]`` of the record.
    """

    source = _h.existing(fasta, "fasta")
    which = _h.choice(engine, _ENGINES, "engine")
    cap = _h.integer(max_sequence_chars, "max_sequence_chars", minimum=1, maximum=10_000_000)
    bioframe = require("bioframe", _EXTRA)
    require("pyfaidx" if which == "pyfaidx" else "pysam", "bio-pyfaidx" if which == "pyfaidx" else "bio-pysam")
    try:
        records = bioframe.load_fasta(str(source), engine=which)
    except (OSError, ValueError, KeyError) as exc:
        raise ConfigurationError(f"bioframe could not load {source}: {exc}") from exc
    names = [{"name": str(name), "length": len(record)} for name, record in records.items()]
    pieces = []
    for region in regions or []:
        chrom, start, end = bioframe.parse_region_string(_h.text(region, "regions"))
        if chrom not in records:
            raise ConfigurationError(f"{chrom!r} is not a record of the FASTA file.")
        record = records[chrom]
        lo = 0 if start is None else int(start)
        hi = len(record) if end is None else int(end)
        sequence = str(record[lo:hi])
        pieces.append({"region": region, "chrom": str(chrom), "start": lo, "end": hi,
                       "sequence": sequence[:cap], "truncated": len(sequence) > cap})
    return {"status": "completed", "backend": "bioframe", "version": _bioframe_version(), "engine": which,
            "record_count": len(names), "records": names, "regions": pieces}


__all__ = ["bioframe_fetch_ucsc_table", "bioframe_is_complete_ucsc_string", "bioframe_load_fasta"]

_register_public_verbose(globals(), tuple(__all__))
