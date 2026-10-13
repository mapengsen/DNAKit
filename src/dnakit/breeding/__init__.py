"""Local candidate adapters for the breeding capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'calculate', input, options=options)

def fit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'fit', input, options=options)

def optimize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'optimize', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'plot', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'predict', input, options=options)

def select(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'select', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('breeding', 'simulate', input, options=options)

__all__ = ['calculate', 'fit', 'optimize', 'plot', 'predict', 'select', 'simulate']

__candidate_adapter_operations__ = set(['calculate', 'fit', 'optimize', 'plot', 'predict', 'select', 'simulate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
