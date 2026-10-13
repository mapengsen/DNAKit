"""Offline pyensembl helpers for FASTA records, genomic intervals, identifiers and Ensembl URLs.

The FASTA readers read only the file they are given. Interval, locus, identifier and alias
helpers call pyensembl's own functions on the values passed in. The URL builders format
Ensembl FTP paths and make no network request. Genome objects that download and index
annotations are not wrapped here; ``gene_backends`` covers GTF-indexed gene queries.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-pyensembl]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-pyensembl"
_DISTRIBUTION = "pyensembl"
_DEFAULT_MAX_FASTA_BYTES = 64_000_000
_MAX_FASTA_BYTES = 1_000_000_000
_MAX_INTERVALS = 1_000_000
_MAX_LOCI = 100_000
_MAX_NAMES = 100_000
_MAX_SEQUENCES = 1_000_000
_MAX_SEQUENCE_CHARS = 100_000_000
_MAX_POSITION = 2**53
_MAX_RELEASE = 1000
_RELEASE_NAME = re.compile(r"[0-9]{4}_[0-9]{2}")
_FASTA_TYPES = ("toplevel", "primary_assembly")
_MASKS = ("none", "soft", "hard")
_SEQUENCE_TYPES = ("cdna", "ncrna", "pep")
_STRANDS = ("+", "-")


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _sequence_of(values: object, name: str, *, limit: int) -> list[object]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list or tuple.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} has {len(values)} items; the limit is {limit}.")
    return list(values)


def _release(value: object) -> int | str:
    """Accept a numbered Ensembl release or an annotation date such as 2026_04."""
    if isinstance(value, str) and _RELEASE_NAME.fullmatch(value):
        return value
    return _int(value, "ensembl_release", 1, _MAX_RELEASE)


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version(_DISTRIBUTION)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"{_DISTRIBUTION} {version}", **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": provenance}


def _upstream(operation: str, call: Callable[[], Any]) -> Any:
    """Run one upstream call and report unexpected failures as backend errors."""
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"pyensembl {operation} failed: {exc}") from exc


def _sequence_text(value: object, name: str) -> str:
    if not isinstance(value, str) or "\x00" in value or len(value) > _MAX_SEQUENCE_CHARS:
        raise ConfigurationError(f"{name} must be a sequence string without NUL characters.")
    return value


def _input_error(call: Callable[[], Any]) -> Any:
    """Map pyensembl's ValueError and TypeError on user-supplied values to ConfigurationError."""
    try:
        return call()
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(str(exc)) from exc


def read_fasta_records(
    fasta_path: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_FASTA_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a FASTA or gzipped FASTA file into identifier and sequence records."""
    path = input_path(fasta_path)
    limit = _int(max_bytes, "max_bytes", 1, _MAX_FASTA_BYTES)
    size = path.stat().st_size
    if size > limit:
        raise ConfigurationError(f"FASTA file has {size} bytes; max_bytes is {limit}.")
    fasta = require("pyensembl.fasta", _EXTRA)
    sequences = _upstream("parse_fasta_dictionary", lambda: fasta.parse_fasta_dictionary(str(path)))
    records = [
        {"id": identifier, "length": len(sequences[identifier]), "sequence": sequences[identifier]}
        for identifier in sorted(sequences)
    ]
    return _document(
        "read_fasta_records",
        _provenance(input_format="FASTA or gzipped FASTA", ordering="sorted by identifier"),
        records=records,
        record_count=len(records),
    )


def read_first_fasta_id(fasta_path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Return the identifier of the first record in a FASTA file, or None when it has none."""
    path = input_path(fasta_path)
    fasta = require("pyensembl.fasta", _EXTRA)
    identifier = _upstream("first_fasta_id", lambda: fasta.first_fasta_id(str(path)))
    return _document(
        "read_first_fasta_id",
        _provenance(input_format="FASTA or gzipped FASTA"),
        identifier=identifier,
    )


def _interval(item: object, position: int) -> tuple[int, int]:
    if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 2:
        raise ConfigurationError(f"intervals[{position}] must be a (start, end) pair.")
    start = _int(item[0], f"intervals[{position}].start", 1, _MAX_POSITION)
    end = _int(item[1], f"intervals[{position}].end", 1, _MAX_POSITION)
    if end < start:
        raise ConfigurationError(f"intervals[{position}] has end before start.")
    return start, end


def merge_genomic_intervals(
    intervals: Sequence[Sequence[int]],
    *,
    max_intervals: int = _MAX_INTERVALS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge 1-based inclusive intervals that overlap or touch, in ascending order."""
    limit = _int(max_intervals, "max_intervals", 1, _MAX_INTERVALS)
    items = _sequence_of(intervals, "intervals", limit=limit)
    pairs = [_interval(item, position) for position, item in enumerate(items)]
    common = require("pyensembl.common", _EXTRA)
    merged = _upstream("merge_intervals", lambda: common.merge_intervals(pairs))
    return _document(
        "merge_genomic_intervals",
        _provenance(coordinate_system="1-based inclusive", merge_rule="overlap or adjacent"),
        input_count=len(pairs),
        intervals=[{"start": start, "end": end} for start, end in merged],
    )


def _locus_record(locus: Any) -> dict[str, Any]:
    return {"contig": locus.contig, "start": locus.start, "end": locus.end, "strand": locus.strand}


def _make_locus(module: Any, contig: str, start: int, end: int, strand: str) -> Any:
    return _input_error(lambda: module.Locus(contig, start, end, strand))


def find_nearest_genomic_locus(
    start: int,
    end: int,
    loci: Sequence[Mapping[str, object]],
    *,
    max_loci: int = _MAX_LOCI,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the locus nearest to a 1-based inclusive interval, using pyensembl distances."""
    begin = _int(start, "start", 1, _MAX_POSITION)
    finish = _int(end, "end", 1, _MAX_POSITION)
    if finish < begin:
        raise ConfigurationError("end must not be before start.")
    limit = _int(max_loci, "max_loci", 1, _MAX_LOCI)
    items = _sequence_of(loci, "loci", limit=limit)
    locus_module = require("pyensembl.locus", _EXTRA)
    search = require("pyensembl.search", _EXTRA)
    objects = []
    for position, item in enumerate(items):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"loci[{position}] must be a mapping.")
        contig = _text(item.get("contig"), f"loci[{position}].contig", limit=128)
        strand = item.get("strand")
        if strand not in _STRANDS:
            raise ConfigurationError(f"loci[{position}].strand must be '+' or '-'.")
        locus_start = _int(item.get("start"), f"loci[{position}].start", 1, _MAX_POSITION)
        locus_end = _int(item.get("end"), f"loci[{position}].end", 1, _MAX_POSITION)
        objects.append(_make_locus(locus_module, contig, locus_start, locus_end, strand))
    distance, nearest = _upstream(
        "find_nearest_locus", lambda: search.find_nearest_locus(begin, finish, objects)
    )
    nearest_index = None
    if nearest is not None:
        nearest_index = next(index for index, item in enumerate(objects) if item is nearest)
    finite = distance if isinstance(distance, int) else None
    return _document(
        "find_nearest_genomic_locus",
        _provenance(coordinate_system="1-based inclusive"),
        start=begin,
        end=finish,
        distance=finite,
        nearest_index=nearest_index,
        nearest=None if nearest is None else _locus_record(nearest),
    )


def stable_ensembl_id(identifier: str, *, verbose: bool = True) -> dict[str, Any]:
    """Strip an Ensembl version suffix such as .15 and keep other identifier formats unchanged."""
    text = _text(identifier, "identifier")
    gene_aliases = require("pyensembl.gene_name_aliases", _EXTRA)
    stable = _upstream("stable_ensembl_gene_id", lambda: gene_aliases.stable_ensembl_gene_id(text))
    return _document(
        "stable_ensembl_id",
        _provenance(version_rule="strip a trailing .N suffix"),
        identifier=text,
        stable_id=stable,
    )


def normalize_chromosome_name(name: str | int, *, verbose: bool = True) -> dict[str, Any]:
    """Normalize a contig name, for example chrx to chrX and 7 to 7."""
    if isinstance(name, bool) or not isinstance(name, (str, int)):
        raise ConfigurationError("name must be a string or an integer.")
    if isinstance(name, str):
        _text(name, "name", limit=128)
    normalization = require("pyensembl.normalization", _EXTRA)
    chromosome = _input_error(lambda: normalization.normalize_chromosome(name))
    return _document(
        "normalize_chromosome_name",
        _provenance(rule="chr prefix upper-cases the suffix; alphabetic names upper-case"),
        input=str(name),
        chromosome=chromosome,
    )


def normalize_strand_symbol(strand: str | int, *, verbose: bool = True) -> dict[str, Any]:
    """Map strand notation ('+', '-', '1', '-1', 1 or -1) to '+' or '-'."""
    if isinstance(strand, bool) or not isinstance(strand, (str, int)):
        raise ConfigurationError("strand must be a string or an integer.")
    normalization = require("pyensembl.normalization", _EXTRA)
    symbol = _input_error(lambda: normalization.normalize_strand(strand))
    return _document(
        "normalize_strand_symbol",
        _provenance(),
        input=str(strand),
        strand=symbol,
    )


def normalize_gene_aliases(
    aliases: Mapping[str, str | Sequence[str]],
    *,
    max_names: int = _MAX_NAMES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalize a gene-name to gene-ID mapping; each name keeps a sorted set of IDs."""
    if not isinstance(aliases, Mapping):
        raise ConfigurationError("aliases must be a mapping from gene names to gene IDs.")
    limit = _int(max_names, "max_names", 1, _MAX_NAMES)
    if len(aliases) > limit:
        raise ConfigurationError(f"aliases has {len(aliases)} names; max_names is {limit}.")
    gene_aliases = require("pyensembl.gene_name_aliases", _EXTRA)
    normalized = _input_error(lambda: gene_aliases.normalize_aliases(dict(aliases)))
    return _document(
        "normalize_gene_aliases",
        _provenance(ordering="sorted by name; IDs sorted and de-duplicated"),
        records=[{"name": name, "gene_ids": list(ids)} for name, ids in sorted(normalized.items())],
        name_count=len(normalized),
    )


def match_identifier_version(
    identifier: str,
    installed: Mapping[str, int | None],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pick the stored form of an Ensembl ID from installed IDs and their recorded versions."""
    text = _text(identifier, "identifier", limit=128)
    if not isinstance(installed, Mapping):
        raise ConfigurationError("installed must map stored IDs to versions or None.")
    stored: dict[str, int | None] = {}
    for key, version in installed.items():
        stored_id = _text(key, "installed key", limit=128)
        if version is not None:
            version = _int(version, f"installed[{stored_id}]", 0, _MAX_POSITION)
        stored[stored_id] = version
    gene_aliases = require("pyensembl.gene_name_aliases", _EXTRA)
    versioned_ids = require("pyensembl.versioned_ids", _EXTRA)
    # Upstream returns the only installed form for any bare ID, so compare stable IDs first.
    stable = gene_aliases.stable_ensembl_gene_id(text)
    same_stable = {
        key: version
        for key, version in stored.items()
        if gene_aliases.stable_ensembl_gene_id(key) == stable
    }
    matched = None
    if same_stable:
        matched = _input_error(lambda: versioned_ids.match_version(text, same_stable))
    return _document(
        "match_identifier_version",
        _provenance(rule="compare stable IDs first; exact match, then version check"),
        identifier=text,
        matched=matched,
        installed_count=len(stored),
        compared_count=len(same_stable),
    )


def lookup_sequence_version_fallback(
    sequences: Mapping[str, str],
    identifier: str,
    *,
    version: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up a sequence by ID, tolerating .N version suffixes in either direction."""
    if not isinstance(sequences, Mapping):
        raise ConfigurationError("sequences must be a mapping from IDs to sequences.")
    if len(sequences) > _MAX_SEQUENCES:
        raise ConfigurationError(f"sequences has more than {_MAX_SEQUENCES} entries.")
    stored: dict[str, str] = {}
    for key, value in sequences.items():
        stored[_text(key, "sequence ID", limit=256)] = _sequence_text(value, f"sequence {key}")
    text = _text(identifier, "identifier", limit=256)
    recorded = None if version is None else _int(version, "version", 0, _MAX_POSITION)
    sequence_data = require("pyensembl.sequence_data", _EXTRA)
    found = _upstream(
        "lookup_sequence_with_version_fallback",
        lambda: sequence_data.lookup_sequence_with_version_fallback(stored, text, recorded),
    )
    return _document(
        "lookup_sequence_version_fallback",
        _provenance(
            resolution="literal, then strip .N, then stripped index",
            note="plain mappings have no stripped index, so a bare ID misses a versioned key",
        ),
        identifier=text,
        version=recorded,
        found=found is not None,
        sequence=found,
    )


def normalize_genome_reference_name(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the canonical spelling of a genome reference name, for example grch38 to GRCh38."""
    text = _text(name, "name", limit=64)
    reference = require("pyensembl.reference_name", _EXTRA)
    canonical = _input_error(lambda: reference.normalize_reference_name(text))
    return _document(
        "normalize_genome_reference_name",
        _provenance(match="case-insensitive"),
        input=text,
        reference=canonical,
    )


def ensembl_reference_name(
    species_name: str,
    ensembl_release: int | str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the reference assembly name that an Ensembl release uses for a species."""
    species = _text(species_name, "species_name", limit=128)
    release = _release(ensembl_release)
    reference = require("pyensembl.reference_name", _EXTRA)
    assembly = _input_error(lambda: reference.which_reference(species, release))
    return _document(
        "ensembl_reference_name",
        _provenance(source="pyensembl species registry"),
        species=species,
        ensembl_release=release,
        reference=assembly,
    )


def ensembl_gtf_url(
    ensembl_release: int,
    species: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the Ensembl FTP URL of the gene annotation GTF for a species and release."""
    release = _int(ensembl_release, "ensembl_release", 1, _MAX_RELEASE)
    name = _text(species, "species", limit=128)
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    url = _input_error(lambda: templates.make_gtf_url(release, name))
    return _document(
        "ensembl_gtf_url",
        _provenance(network="none; URL string only"),
        species=name,
        ensembl_release=release,
        url=url,
    )


def ensembl_genome_fasta_url(
    ensembl_release: int,
    species: str,
    *,
    fasta_type: str = "toplevel",
    mask: str = "none",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the Ensembl FTP URL of the combined reference DNA FASTA for a species and release."""
    release = _int(ensembl_release, "ensembl_release", 1, _MAX_RELEASE)
    name = _text(species, "species", limit=128)
    if fasta_type not in _FASTA_TYPES:
        raise ConfigurationError(f"fasta_type must be one of {', '.join(_FASTA_TYPES)}.")
    if mask not in _MASKS:
        raise ConfigurationError(f"mask must be one of {', '.join(_MASKS)}.")
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    url = _input_error(lambda: templates.make_genome_fasta_url(release, name, fasta_type, mask))
    return _document(
        "ensembl_genome_fasta_url",
        _provenance(network="none; URL string only"),
        species=name,
        ensembl_release=release,
        fasta_type=fasta_type,
        mask=mask,
        url=url,
    )


def ensembl_sequence_fasta_url(
    ensembl_release: int,
    species: str,
    sequence_type: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the Ensembl FTP URL of the cDNA, ncRNA or protein FASTA for a species and release."""
    release = _int(ensembl_release, "ensembl_release", 1, _MAX_RELEASE)
    name = _text(species, "species", limit=128)
    if sequence_type not in _SEQUENCE_TYPES:
        raise ConfigurationError(f"sequence_type must be one of {', '.join(_SEQUENCE_TYPES)}.")
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    url = _input_error(lambda: templates.make_fasta_url(release, name, sequence_type))
    return _document(
        "ensembl_sequence_fasta_url",
        _provenance(network="none; URL string only"),
        species=name,
        ensembl_release=release,
        sequence_type=sequence_type,
        url=url,
    )


__all__ = [
    "ensembl_genome_fasta_url",
    "ensembl_gtf_url",
    "ensembl_reference_name",
    "ensembl_sequence_fasta_url",
    "find_nearest_genomic_locus",
    "lookup_sequence_version_fallback",
    "match_identifier_version",
    "merge_genomic_intervals",
    "normalize_chromosome_name",
    "normalize_gene_aliases",
    "normalize_genome_reference_name",
    "normalize_strand_symbol",
    "read_fasta_records",
    "read_first_fasta_id",
    "stable_ensembl_id",
]

_register_public_verbose(globals(), tuple(__all__))
