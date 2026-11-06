"""Candidate entry points of the databases capability family.

Only ``query`` and ``map_ids`` delegate to existing DNAKit database functions; the other operations have no
matching DNAKit function (they would need a particular database, release catalogue or cache layout) and
return ``status: "not_implemented"``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from dnakit._candidate_runtime import (
    candidate_not_implemented,
    candidate_ok,
    run_candidate_operation,
)
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import BackendUnavailableError, ConfigurationError, InputFormatError


def _not_implemented(
    operation: str, input: object, options: Mapping[str, object] | None, requires: str, message: str
) -> object:
    return candidate_not_implemented(
        "databases", operation, input, options=options, requires=requires, message=message
    )


def cache(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: no database-result cache is defined for the database adapters."""

    return _not_implemented(
        "cache",
        input,
        options,
        "a database-result cache with a defined layout and expiry policy",
        "The database adapters keep no cache; there is nothing for this function to populate or inspect.",
    )


def download(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: downloading needs a named database release; use the specific download/fetch functions."""

    return _not_implemented(
        "download",
        input,
        options,
        "a named database release and network access",
        "No generic database download exists; use dnakit.databases.search to find the specific operation.",
    )


def export(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: there is no generic export of database content."""

    return _not_implemented(
        "export",
        input,
        options,
        "a database snapshot and a target format",
        "No generic database export exists; export the result of a query with the io functions.",
    )


def fetch(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: fetching needs a named database; use query or the specific fetch_* functions."""

    return _not_implemented(
        "fetch",
        input,
        options,
        "a named database and network access",
        "fetch does not name a database; call dnakit.databases.query with an operation name instead.",
    )


def list_versions(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Not implemented: DNAKit keeps no catalogue of database release versions."""

    return _not_implemented(
        "list_versions",
        input,
        options,
        "a release catalogue of the database",
        "No DNAKit function lists database releases.",
    )


def map_ids(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Map gene identifiers to MyGene.info records (network) by delegating to ``map_gene_identifiers``.

    ``input``: a list of identifiers or a comma-separated string.  Options: ``scopes`` (identifier type, default
    "entrezgene"), ``species`` (default "human"), ``fields`` (default "symbol,name,taxid").  Needs network access and
    the optional ``mygene`` package; without it the result is ``not_implemented``.  The result is the
    ``map_gene_identifiers`` document (one row per identifier under ``results``).
    """

    from dnakit.databases import gene_backends

    settings = dict(options or {})
    if isinstance(input, str):
        identifiers = [item.strip() for item in input.split(",") if item.strip()]
    elif isinstance(input, Sequence):
        identifiers = list(input)
    else:
        raise InputFormatError(
            "input must be a list of identifiers or a comma-separated string.", code="INVALID_INPUT"
        )
    if not identifiers:
        raise InputFormatError("No identifiers were given.", code="INVALID_INPUT")
    keywords = {key: settings[key] for key in ("scopes", "species", "fields") if key in settings}
    try:
        result = gene_backends.map_gene_identifiers(identifiers, verbose=False, **keywords)
    except BackendUnavailableError as exc:
        return candidate_not_implemented(
            "databases",
            "map_ids",
            input,
            options=options,
            requires="pip install mygene",
            message=str(exc),
        )
    return candidate_ok(
        "databases",
        "map_ids",
        input,
        result,
        options=options,
        implementation="dnakit.databases.map_gene_identifiers",
    )


def query(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    """Call one bundled DNAKit database operation by name (network), via ``dnakit.databases.call``.

    ``input``: the operation name (see ``dnakit.databases.list_operations`` / ``search``) with its arguments in
    ``options["arguments"]``, or a mapping ``{"name": operation, "arguments": {...}}``.  Arguments are validated
    against the operation's input schema; provider errors raise.  Option ``max_result_bytes`` (default 20,000,000)
    bounds the result.  The result is the provider's JSON document.
    """

    from dnakit.databases import _runtime

    settings = dict(options or {})
    if isinstance(input, Mapping):
        if not isinstance(input.get("name"), str) or set(input) - {"name", "arguments"}:
            raise InputFormatError(
                "A query mapping needs 'name' and optionally 'arguments'.", code="INVALID_INPUT"
            )
        name, arguments = input["name"], input.get("arguments")
    elif isinstance(input, str) and input.strip():
        name, arguments = input.strip(), settings.get("arguments")
    else:
        raise InputFormatError(
            "input must be an operation name or a {'name', 'arguments'} mapping.",
            code="INVALID_INPUT",
        )
    if arguments is not None and not isinstance(arguments, Mapping):
        raise ConfigurationError("arguments must be a mapping.", code="INVALID_DATABASE_ARGUMENTS")
    extra = (
        {"max_result_bytes": settings["max_result_bytes"]} if "max_result_bytes" in settings else {}
    )
    result = _runtime.call(name, arguments, verbose=False, **extra)
    return candidate_ok(
        "databases", "query", input, result, options=options, implementation="dnakit.databases.call"
    )


def search(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return run_candidate_operation("databases", "search", input, options=options)


__all__ = ["cache", "download", "export", "fetch", "list_versions", "map_ids", "query"]

_register_public_verbose(globals())
