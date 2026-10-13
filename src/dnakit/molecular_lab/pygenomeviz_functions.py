"""pyGenomeViz GenBank summaries and accession fetch (optional).

pyGenomeViz and Biopython are imported only inside each function. The fetch
function calls NCBI Entrez over the network; its tests replace the Entrez client
with a fake. The example-dataset loaders, the figure workflows and the CLI entry
points are not wrapped.
"""

from __future__ import annotations

import importlib.metadata
import io
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_MAX_SOURCES = 100
_MAX_TEXT = 256
_ACCESSION = re.compile(r"^[A-Za-z]{1,4}_?[0-9]{5,9}(\.[0-9]{1,3})?$")
_EMAIL = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version("pygenomeviz")
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": "pygenomeviz", "version": version, **notes}


def _sources(values: object) -> list[str]:
    if isinstance(values, (str, bytes, Path)) or not isinstance(values, Sequence):
        raise ConfigurationError("sources must be a list of GenBank file paths.")
    if not 1 <= len(values) <= _MAX_SOURCES:
        raise ConfigurationError(f"sources must contain 1 to {_MAX_SOURCES} files.")
    specs: list[str] = []
    for value in values:
        if not isinstance(value, (str, Path)):
            raise ConfigurationError("each source must be a path or a 'path[:min-max[:-1]]' spec.")
        spec = str(value)
        if len(spec) > 4096 or "\x00" in spec:
            raise ConfigurationError("source specification is too long or contains NUL.")
        input_path(spec.split(":", 1)[0])
        specs.append(spec)
    return specs


def _genome_row(genome: Any) -> dict[str, Any]:
    records = sorted(genome.records, key=lambda record: str(record.id))
    if not records:
        raise ConfigurationError(f"{genome.name} contains no GenBank records.")
    return {
        "name": str(genome.name),
        "full_length": int(genome.full_genome_length),
        "range": [int(genome.min_range), int(genome.max_range)],
        "range_size": int(genome.range_size),
        "reverse": bool(genome.reverse),
        "gc_percent": round(float(genome.calc_genome_gc_content()), 6),
        "records": [
            {
                "seqid": str(record.id),
                "length": len(record.seq),
                "feature_count": len(record.features),
            }
            for record in records
        ],
    }


def summarize_genbank_records(
    sources: Sequence[str | Path], *, verbose: bool = True
) -> dict[str, Any]:
    """Summarize GenBank files through pyGenomeViz's ``Genbank`` parser.

    Each source is a local path, optionally followed by ``:min-max`` and
    ``:-1``. ``:-1`` reverse-complements the genome. The output lists the genome
    length, the range, the GC percentage and every record, sorted by sequence ID.
    """

    scripts = require("pygenomeviz.scripts", "bio-pygenomeviz")
    specs = _sources(sources)
    try:
        genomes = scripts.gbk_files2gbk_objects(specs)
    except (OSError, ValueError, IndexError, TypeError) as exc:
        raise ConfigurationError(f"pyGenomeViz could not read GenBank input: {exc}") from exc
    rows = [_genome_row(genome) for genome in genomes]
    return {
        "status": "completed",
        "operation": "summarize_genbank_records",
        "genome_count": len(rows),
        "genomes": rows,
        "provenance": _provenance(parser="pygenomeviz.parser.genbank.Genbank"),
    }


def _accession(value: object) -> str:
    if not isinstance(value, str) or not _ACCESSION.match(value.strip()):
        raise ConfigurationError("accession must look like NC_000913.3, JX128258.1 or CP000001.")
    return value.strip()


def _new_output(value: object) -> Path:
    text = _text_path(value)
    path = Path(text).expanduser().resolve()
    if not path.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"output already exists; choose a new file: {path}")
    return path


def _text_path(value: object) -> str:
    if isinstance(value, Path):
        value = str(value)
    if not isinstance(value, str) or not value.strip() or len(value) > 4096:
        raise ConfigurationError("output must be a file path of at most 4096 characters.")
    return value


def fetch_genbank_by_accession(
    accession: str,
    *,
    output: str | Path | None = None,
    email: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch one nucleotide GenBank record from NCBI by accession (network call).

    ``email`` is passed to NCBI as the contact address. ``output`` must be a new
    file path; the fetched text is written there. The record is parsed with
    Biopython, so a reply that is not exactly one GenBank record raises
    ``ConfigurationError``. A failed request raises ``BackendExecutionError``.
    """

    utils = require("pygenomeviz.utils", "bio-pygenomeviz")
    seqio = require("Bio.SeqIO", "bio-pygenomeviz")
    acc = _accession(accession)
    target = _new_output(output) if output is not None else None
    contact = None
    if email is not None:
        if not isinstance(email, str) or len(email) > _MAX_TEXT or not _EMAIL.match(email):
            raise ConfigurationError("email must be a plain address such as name@example.org.")
        contact = email
    # Upstream may overwrite the global Bio.Entrez.email; restore the caller's value afterwards.
    entrez = require("Bio.Entrez", "bio-pygenomeviz")
    previous_email = getattr(entrez, "email", None)
    try:
        handle = utils.fetch_genbank_by_accid(
            acc, gbk_outfile=str(target) if target is not None else None, email=contact
        )
        text = handle.read()
        handle.close()
    except (OSError, ValueError, RuntimeError) as exc:
        raise BackendExecutionError(f"NCBI fetch failed for {acc}: {exc}") from exc
    finally:
        entrez.email = previous_email
    try:
        record = seqio.read(io.StringIO(text), "genbank")
    except ValueError as exc:
        raise ConfigurationError(
            f"NCBI did not return one GenBank record for {acc}: {exc}"
        ) from exc
    return {
        "status": "completed",
        "operation": "fetch_genbank_by_accession",
        "accession": acc,
        "record_id": str(record.id),
        "name": str(record.name),
        "description": str(record.description),
        "length": len(record.seq),
        "output": str(target) if target is not None else None,
        "bytes": len(text.encode("utf-8")),
        "provenance": _provenance(
            source="NCBI Entrez efetch, nucleotide, gbwithparts", network=True
        ),
    }


__all__ = ["fetch_genbank_by_accession", "summarize_genbank_records"]

_register_public_verbose(globals(), tuple(__all__))
