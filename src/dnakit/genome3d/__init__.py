"""Local candidate adapters for the genome3d capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def aggregate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'aggregate', input, options=options)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'annotate', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'calculate', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'compare', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'export', input, options=options)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'normalize', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'plot', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'read', input, options=options)

def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('genome3d', 'write', input, options=options)

__all__ = ['aggregate', 'annotate', 'calculate', 'compare', 'export', 'normalize', 'plot', 'read', 'write']

__candidate_adapter_operations__ = set(['aggregate', 'annotate', 'calculate', 'compare', 'export', 'normalize', 'plot', 'read', 'write'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
