"""Local candidate adapters for the databases capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def cache(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'cache', input, options=options)

def download(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'download', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'export', input, options=options)

def fetch(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'fetch', input, options=options)

def list_versions(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'list_versions', input, options=options)

def map_ids(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'map_ids', input, options=options)

def query(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'query', input, options=options)

def search(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('databases', 'search', input, options=options)

__all__ = ['cache', 'download', 'export', 'fetch', 'list_versions', 'map_ids', 'query']

from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
