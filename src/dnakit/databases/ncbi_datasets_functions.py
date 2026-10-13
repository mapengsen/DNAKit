"""NCBI Datasets assembly metadata lookups and offline field selectors for genes and assemblies.

The assembly lookups call the NCBI Datasets REST API through ncbi-datasets-pylib and need the
network. The field selectors build openapi descriptors from a plain mapping and make no request.
Upstream prints warnings and API errors to stdout; they are returned as ``messages`` instead.
The ``print_*`` helpers and the gRPC servicer registrations are not wrapped.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-ncbi-datasets]"`` hint.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import io
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-ncbi-datasets"
_DISTRIBUTION = "ncbi-datasets-pylib"
_MAX_ACCESSIONS = 100
_MAX_RECORDS = 10_000
_MAX_TAXON_NAME = 64
_MAX_TAX_ID = 2**63 - 1
_ASSEMBLY = re.compile(r"GC[AF]_[0-9]{9}(?:\.[0-9]+)?")
_BIOPROJECT = re.compile(r"PRJ[A-Z]{1,2}[0-9]{1,12}")
_TAXON_NAME = re.compile(r"[A-Za-z][A-Za-z0-9 _.-]{0,63}")
_GENE_ID = re.compile(r"[0-9]{1,20}")
_SYMBOL = re.compile(r"[A-Za-z0-9_.@-]{1,64}")
_GENE_FIELDS = ("gene_id", "symbol")


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


def _pattern(value: object, name: str, pattern: re.Pattern[str]) -> str:
    text = _text(value, name, limit=128)
    if not pattern.fullmatch(text):
        raise ConfigurationError(f"{name} has an unsupported format: {text!r}")
    return text


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _sequence_of(values: object, name: str, *, limit: int, minimum: int = 0) -> list[object]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list or tuple.")
    if len(values) < minimum:
        raise ConfigurationError(f"{name} must contain at least {minimum} item.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} has {len(values)} items; the limit is {limit}.")
    return list(values)


def _plain(value: object) -> Any:
    """Convert openapi models, numpy scalars and non-finite floats to JSON-compatible data."""
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        value = to_dict()
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if value == value and abs(value) != float("inf") else None
    return str(value)


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version(_DISTRIBUTION)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"{_DISTRIBUTION} {version}", **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": provenance}


def _upstream(operation: str, call: Callable[[], Any]) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"ncbi-datasets-pylib {operation} failed: {exc}") from exc


def _model(name: str, build: Callable[[], Any]) -> Any:
    try:
        return build()
    except Exception as exc:
        raise ConfigurationError(f"{name} rejected the record: {exc}") from exc


def _collect(
    operation: str,
    source: Callable[[], Iterable[object]],
    *,
    limit: int,
    provenance: dict[str, Any],
    **fields: object,
) -> dict[str, Any]:
    """Drain an upstream iterator with a record limit; capture its stdout messages."""
    records: list[Any] = []
    buffer = io.StringIO()

    def run() -> None:
        with contextlib.redirect_stdout(buffer):
            iterator = iter(source())
            try:
                for item in iterator:
                    if len(records) >= limit:
                        raise ConfigurationError(f"more than max_records={limit} records returned.")
                    records.append(_plain(item))
            finally:
                close = getattr(iterator, "close", None)
                if callable(close):
                    close()

    _upstream(operation, run)
    messages = [line.strip() for line in buffer.getvalue().splitlines() if line.strip()]
    return _document(
        operation,
        provenance,
        record_count=len(records),
        records=records,
        messages=messages,
        **fields,
    )


def fetch_ncbi_assembly_metadata(
    accessions: Sequence[str],
    *,
    max_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch NCBI Datasets assembly metadata for GCF or GCA accessions (network)."""
    items = _sequence_of(accessions, "accessions", limit=_MAX_ACCESSIONS, minimum=1)
    ids = [_pattern(item, "accession", _ASSEMBLY) for item in items]
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    genome = require("ncbi.datasets.metadata.genome", _EXTRA)
    return _collect(
        "fetch_ncbi_assembly_metadata",
        lambda: genome.get_assembly_metadata_by_asm_accessions(ids),
        limit=limit,
        provenance=_provenance(source="NCBI Datasets assembly descriptors by accession"),
        requested=ids,
    )


def fetch_ncbi_assemblies_by_taxon(
    taxon: int | str,
    *,
    max_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """List NCBI Datasets assembly metadata for a taxonomy ID or name (network)."""
    if isinstance(taxon, bool) or not isinstance(taxon, (int, str)):
        raise ConfigurationError("taxon must be an integer tax ID or a taxon name.")
    if isinstance(taxon, int):
        value: int | str = _int(taxon, "taxon", 1, _MAX_TAX_ID)
    else:
        value = _pattern(taxon, "taxon", _TAXON_NAME)
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    genome = require("ncbi.datasets.metadata.genome", _EXTRA)
    return _collect(
        "fetch_ncbi_assemblies_by_taxon",
        lambda: genome.get_assembly_metadata_by_taxon(value),
        limit=limit,
        provenance=_provenance(source="NCBI Datasets assembly descriptors by taxon"),
        requested=str(value),
    )


def fetch_ncbi_assemblies_by_bioproject(
    bioproject_accessions: Sequence[str],
    *,
    max_records: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """List NCBI Datasets assembly metadata for BioProject accessions such as PRJNA... (network)."""
    items = _sequence_of(
        bioproject_accessions, "bioproject_accessions", limit=_MAX_ACCESSIONS, minimum=1
    )
    ids = [_pattern(item, "bioproject accession", _BIOPROJECT) for item in items]
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    genome = require("ncbi.datasets.metadata.genome", _EXTRA)
    return _collect(
        "fetch_ncbi_assemblies_by_bioproject",
        lambda: genome.get_assembly_metadata_by_bioproject_accessions(ids),
        limit=limit,
        provenance=_provenance(source="NCBI Datasets assembly descriptors by BioProject"),
        requested=ids,
    )


def select_ncbi_gene_fields(
    gene: Mapping[str, object],
    fields: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select gene_id and symbol from one NCBI gene record (no network)."""
    if not isinstance(gene, Mapping):
        raise ConfigurationError("gene must be a mapping with gene_id and symbol.")
    names = [
        _text(item, "fields item", limit=32) for item in _sequence_of(fields, "fields", limit=2)
    ]
    unknown = sorted(name for name in names if name not in _GENE_FIELDS)
    if unknown or not names:
        raise ConfigurationError(f"fields must be a non-empty subset of {list(_GENE_FIELDS)}.")
    gene_id = _pattern(gene.get("gene_id"), "gene_id", _GENE_ID)
    symbol = _pattern(gene.get("symbol"), "symbol", _SYMBOL)
    models = require("ncbi.datasets.openapi.model.v1_gene_descriptor", _EXTRA)
    metadata = require("ncbi.datasets.metadata.gene", _EXTRA)
    descriptor = _model(
        "V1GeneDescriptor",
        lambda: models.V1GeneDescriptor(gene_id=gene_id, symbol=symbol),
    )
    values = _upstream(
        "gene_values_by_fields",
        lambda: metadata.gene_values_by_fields(descriptor, names),
    )
    return _document(
        "select_ncbi_gene_fields",
        _provenance(source="openapi V1GeneDescriptor; gene_id and symbol only"),
        fields=names,
        values=_plain(values),
    )


def select_ncbi_assembly_fields(
    assembly: Mapping[str, object],
    fields: Sequence[str] = (),
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select top-level fields from one NCBI assembly record; all fields when none are named."""
    if not isinstance(assembly, Mapping) or not all(
        isinstance(key, str) and key for key in assembly
    ):
        raise ConfigurationError("assembly must be a mapping with string keys.")
    names = [
        _text(item, "fields item", limit=64) for item in _sequence_of(fields, "fields", limit=256)
    ]
    models = require("ncbi.datasets.openapi.model.v1_assembly_dataset_descriptor", _EXTRA)
    metadata = require("ncbi.datasets.metadata.genome", _EXTRA)
    descriptor = _model(
        "V1AssemblyDatasetDescriptor",
        lambda: models.V1AssemblyDatasetDescriptor(**dict(assembly)),
    )
    values = _upstream(
        "assembly_values_by_fields",
        lambda: metadata.assembly_values_by_fields(descriptor, names),
    )
    selected = _plain(values)
    missing = sorted(name for name in names if name not in selected)
    return _document(
        "select_ncbi_assembly_fields",
        _provenance(source="openapi V1AssemblyDatasetDescriptor; top-level fields"),
        fields=names,
        values=selected,
        missing_fields=missing,
    )


__all__ = [
    "fetch_ncbi_assemblies_by_bioproject",
    "fetch_ncbi_assemblies_by_taxon",
    "fetch_ncbi_assembly_metadata",
    "select_ncbi_assembly_fields",
    "select_ncbi_gene_fields",
]

_register_public_verbose(globals(), tuple(__all__))
