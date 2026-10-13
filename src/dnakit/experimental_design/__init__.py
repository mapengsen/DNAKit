"""Local candidate adapters for the experimental design capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def assign(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'assign', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'calculate', input, options=options)

def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'check', input, options=options)

def design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'design', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'export', input, options=options)

def generate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'generate', input, options=options)

def randomize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'randomize', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'simulate', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('experimental_design', 'validate', input, options=options)

__all__ = ['assign', 'calculate', 'check', 'design', 'export', 'generate', 'randomize', 'simulate', 'validate']

__candidate_adapter_operations__ = set(['assign', 'calculate', 'check', 'design', 'export', 'generate', 'randomize', 'simulate', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
