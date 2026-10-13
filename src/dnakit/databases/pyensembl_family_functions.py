"""Ensembl filename and source-URL checks through pyensembl (family-variant check, group H4).

Three pyensembl functions are reached by no parameter of the existing wrappers:
``make_gtf_filename`` and ``make_fasta_filename`` return only the file name that the URL
builders append to a path, and ``is_canonical_source`` is true for every official Ensembl
reference FASTA directory path, including file names that are not toplevel or primary-assembly
DNA (the ``canonical`` flag of ``ensembl_genome_fasta_source_options`` is false for those).
All three work on strings and make no network request. pyensembl is imported only when a
function runs; a missing package raises ``BackendUnavailableError`` with the install hint.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Callable
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-pyensembl"
_MAX_RELEASE = 1000
_SEQUENCE_TYPES = ("cdna", "ncrna", "pep")


def _text(value: object, name: str, *, limit: int) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


def _release(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_RELEASE:
        raise ConfigurationError(f"ensembl_release must be an integer from 1 to {_MAX_RELEASE}.")
    return value


def _provenance() -> dict[str, Any]:
    try:
        version = importlib.metadata.version("pyensembl")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"pyensembl {version}", "network": "none; file name or pattern only"}


def _input_error(call: Callable[[], Any]) -> Any:
    """Map pyensembl's ValueError and TypeError on user-supplied values to ConfigurationError."""
    try:
        return call()
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(str(exc)) from exc
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"pyensembl call failed: {exc}") from exc


def ensembl_gtf_filename(
    ensembl_release: int,
    species: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the file name of the Ensembl gene annotation GTF for a species and release."""
    release = _release(ensembl_release)
    name = _text(species, "species", limit=128)
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    filename = _input_error(lambda: templates.make_gtf_filename(release, name))
    return {
        "status": "completed",
        "operation": "ensembl_gtf_filename",
        "species": name,
        "ensembl_release": release,
        "filename": filename,
        "provenance": _provenance(),
    }


def ensembl_sequence_fasta_filename(
    ensembl_release: int,
    species: str,
    sequence_type: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the file name of the Ensembl cDNA, ncRNA or protein FASTA of a species and release."""
    release = _release(ensembl_release)
    name = _text(species, "species", limit=128)
    if sequence_type not in _SEQUENCE_TYPES:
        raise ConfigurationError(f"sequence_type must be one of {', '.join(_SEQUENCE_TYPES)}.")
    templates = require("pyensembl.ensembl_url_templates", _EXTRA)
    filename = _input_error(lambda: templates.make_fasta_filename(release, name, sequence_type))
    return {
        "status": "completed",
        "operation": "ensembl_sequence_fasta_filename",
        "species": name,
        "ensembl_release": release,
        "sequence_type": sequence_type,
        "filename": filename,
        "provenance": _provenance(),
    }


def check_ensembl_canonical_source(source: str, *, verbose: bool = True) -> dict[str, Any]:
    """Tell whether a URL is an official Ensembl reference FASTA path (https, known host, layout).

    The check is on the host and the ``.../fasta/<species>/dna/<file>`` path only. A file in
    that directory whose name is not a toplevel or primary-assembly DNA FASTA is also reported
    as official; use ``ensembl_genome_fasta_source_options`` to read the DNA type and mask.
    """
    text = _text(source, "source", limit=4096)
    genome_fasta = require("pyensembl.genome_fasta_cache", _EXTRA)
    official = _input_error(lambda: genome_fasta.is_canonical_source(text))
    return {
        "status": "completed",
        "operation": "check_ensembl_canonical_source",
        "source": text,
        "official_source": bool(official),
        "provenance": _provenance(),
    }


__all__ = [
    "check_ensembl_canonical_source",
    "ensembl_gtf_filename",
    "ensembl_sequence_fasta_filename",
]

_register_public_verbose(globals(), tuple(__all__))
