"""Offline pyensembl helpers for dated Ensembl URLs, reference DNA source options and species.

The URL and source helpers match or format strings and make no network request. The species helper
reads pyensembl's built-in species registry. Release objects that depend on a local annotation
cache (``cached_release``, ``genome_for_reference_name``) are not wrapped here, and neither are the
cache path and maintenance helpers.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-pyensembl]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Callable
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-pyensembl"
_DISTRIBUTION = "pyensembl"
_MAX_TEXT = 2048
_MASKS = ("none", "soft", "hard")
_SPECIES_BY = ("name", "reference")
_ACCESSION = re.compile(r"GC[AF]_[0-9]{9}\.[0-9]+")
_PROVIDER = re.compile(r"[a-z][a-z0-9_]{0,31}")


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


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


def _input_error(call: Callable[[], Any]) -> Any:
    """Map pyensembl's lookup and validation errors on user values to ConfigurationError."""
    try:
        return call()
    except (KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(str(exc)) from exc


def ensembl_dated_release_urls(
    assembly_accession: str,
    provider: str,
    annotation_date: str,
    *,
    include_alt: bool = False,
    genome_fasta_mask: str = "none",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the URLs of one dated Ensembl annotation release: GTF, cDNA, peptide, reference DNA."""
    accession = _text(assembly_accession, "assembly_accession", limit=32)
    if _ACCESSION.fullmatch(accession) is None:
        raise ConfigurationError(
            "assembly_accession must look like GCA_000001405.29 or GCF_000001405.40."
        )
    source = _text(provider, "provider", limit=32)
    if _PROVIDER.fullmatch(source) is None:
        raise ConfigurationError("provider must be a lower-case source name such as 'ensembl'.")
    versions = require("pyensembl.ensembl_versions", _EXTRA)
    if not versions.is_dated_release(annotation_date):
        raise ConfigurationError(
            "annotation_date must be a dated release in YYYY_MM form, for example '2026_04'."
        )
    if not isinstance(include_alt, bool):
        raise ConfigurationError("include_alt must be a boolean.")
    if genome_fasta_mask not in _MASKS:
        raise ConfigurationError(f"genome_fasta_mask must be one of {', '.join(_MASKS)}.")
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    urls = _input_error(
        lambda: templates.make_dated_release_urls(
            accession,
            source,
            annotation_date,
            include_alt=include_alt,
            genome_fasta_mask=genome_fasta_mask,
        )
    )
    return _document(
        "ensembl_dated_release_urls",
        _provenance(network="none; URL string only", server=templates.ENSEMBL_PLATFORM_FTP_SERVER),
        assembly_accession=accession,
        provider=source,
        annotation_date=annotation_date,
        include_alt=include_alt,
        genome_fasta_mask=genome_fasta_mask,
        directory=urls.directory,
        gtf=urls.gtf,
        cdna=urls.cdna,
        pep=urls.pep,
        genome_fasta=urls.genome_fasta,
    )


def ensembl_dated_genome_fasta_mask(url: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the soft or hard mask named by a dated-release reference DNA URL, or None."""
    text = _text(url, "url", limit=_MAX_TEXT)
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    mask = _upstream("dated_genome_fasta_mask", lambda: templates.dated_genome_fasta_mask(text))
    return _document(
        "ensembl_dated_genome_fasta_mask",
        _provenance(network="none; URL suffix match"),
        url=text,
        genome_fasta_mask=mask,
        matched=mask is not None,
    )


def ensembl_genome_fasta_source_options(source: str, *, verbose: bool = True) -> dict[str, Any]:
    """Read the FASTA type and mask of an official Ensembl reference DNA URL."""
    text = _text(source, "source", limit=_MAX_TEXT)
    genome_fasta = require("pyensembl.genome_fasta_cache", _EXTRA)
    options = _upstream(
        "canonical_source_options", lambda: genome_fasta.canonical_source_options(text)
    )
    fasta_type, mask = (None, None) if options is None else options
    return _document(
        "ensembl_genome_fasta_source_options",
        _provenance(network="none; URL pattern match"),
        source=text,
        canonical=options is not None,
        genome_fasta_type=fasta_type,
        genome_fasta_mask=mask,
    )


def ensembl_species_record(
    query: str,
    *,
    by: str = "name",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the species record for a species name or a genome reference name."""
    if by not in _SPECIES_BY:
        raise ConfigurationError("by must be 'name' or 'reference'.")
    text = _text(query, "query", limit=128)
    if by == "name":
        registry = require("pyensembl.species", _EXTRA)
        species = _input_error(lambda: registry.find_species_by_name(text))
    else:
        reference = require("pyensembl.reference_name", _EXTRA)
        species = _input_error(lambda: reference.find_species_by_reference(text))
    dated = species.dated_releases
    return _document(
        "ensembl_species_record",
        _provenance(source="pyensembl species registry", network="none"),
        query=text,
        by=by,
        latin_name=species.latin_name,
        synonyms=sorted(str(name) for name in species.synonyms),
        division=species.division,
        ensembl_genomes=bool(species.ensembl_genomes),
        reference_assemblies=[
            {"reference": name, "first_release": first, "last_release": last}
            for name, (first, last) in sorted(species.reference_assemblies.items())
        ],
        dated_releases=(
            None if dated is None else {"assembly_accession": dated[0], "provider": dated[1]}
        ),
    )


__all__ = [
    "ensembl_dated_genome_fasta_mask",
    "ensembl_dated_release_urls",
    "ensembl_genome_fasta_source_options",
    "ensembl_species_record",
]

_register_public_verbose(globals(), tuple(__all__))
