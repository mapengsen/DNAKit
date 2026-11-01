"""genomicranges operations: UCSC gene-annotation URLs and table reading (network, or a local GTF file).

``ucsc_gtf_url`` builds the UCSC download address of a gene model GTF. ``read_ucsc_gtf`` reads that GTF into a
table: it downloads from UCSC (network) unless a local copy is given, writes the whole table as TSV to a new
file, and returns a summary with the source and the retrieval date. Results are plain dictionaries.
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-genomicranges"
_UCSC_BASE = "http://hgdownload.cse.ucsc.edu/goldenPath"
_ANNOTATIONS = ("refGene", "ensGene", "knownGene", "ncbiRefSeq")


def _annotation(value: str) -> str:
    if value not in _ANNOTATIONS:
        raise ConfigurationError(f"annotation must be one of {', '.join(_ANNOTATIONS)}; got {value!r}.")
    return value


def _genome(value: str) -> str:
    if not isinstance(value, str) or not value.strip() or "/" in value or value.startswith("-"):
        raise ConfigurationError("genome must be a UCSC genome code such as hg38 or mm10.")
    return value.strip()


def ucsc_gtf_url(genome: str, annotation: str = "refGene", *, verbose: bool = True) -> str:
    """Return the UCSC download address of a genome's gene-model GTF (genomicranges access_gtf_ucsc).

    Example: ``ucsc_gtf_url("hg38")`` gives the refGene GTF under goldenPath/hg38/bigZips/genes/. The upstream
    address contains a doubled slash before the file name; this function writes a single slash.
    """
    code = _genome(genome)
    kind = _annotation(annotation)
    return f"{_UCSC_BASE}/{code}/bigZips/genes/{code}.{kind}.gtf.gz"


def read_ucsc_gtf(
    genome: str,
    output_table: str,
    *,
    annotation: str = "refGene",
    local_gtf: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a UCSC gene-model GTF into a table and write it as TSV (genomicranges read_ucsc).

    Without ``local_gtf`` the GTF is downloaded from UCSC (network access needed). With ``local_gtf`` the given
    file (plain or ``.gz``) is read and nothing is downloaded; ``genome`` and ``annotation`` then only label the
    result. ``output_table`` must be a new file. The summary has the row count, the columns, the number of rows
    per chromosome, the source, and the retrieval date.
    """
    code = _genome(genome)
    kind = _annotation(annotation)
    target = Path(output_table).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output_table must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for output_table does not exist: {target.parent}")
    parse_gtf = require("genomicranges.io.gtf", _EXTRA).parse_gtf
    if local_gtf is not None:
        source = Path(local_gtf).expanduser().resolve()
        if not source.is_file():
            raise ConfigurationError(f"local_gtf must be an existing file: {source}")
        source_text = str(source)
        compressed = source.suffix == ".gz"
        retrieved = datetime.datetime.fromtimestamp(source.stat().st_mtime, tz=datetime.timezone.utc).date()
        table = parse_gtf(str(source), compressed=compressed)
    else:
        source_text = ucsc_gtf_url(code, kind)
        retrieved = datetime.date.today()
        table = parse_gtf(source_text, compressed=True)
    table.to_csv(target, sep="\t", index=False)
    counts = table["seqnames"].value_counts()
    return {
        "status": "completed",
        "backend": "genomicranges",
        "genome": code,
        "annotation": kind,
        "source": source_text,
        "source_date": retrieved.isoformat(),
        "row_count": int(len(table)),
        "columns": [str(column) for column in table.columns],
        "rows_per_chromosome": {str(name): int(count) for name, count in counts.items()},
        "output_table": str(target),
    }


__all__ = [
    "read_ucsc_gtf",
    "ucsc_gtf_url",
]
_register_public_verbose(globals())
