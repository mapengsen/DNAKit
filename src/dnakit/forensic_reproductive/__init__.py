"""Local candidate adapters for the forensic reproductive capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'calculate', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'compare', input, options=options)

def deconvolute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'deconvolute', input, options=options)

def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'infer', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'plot', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'predict', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'read', input, options=options)

def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'report', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('forensic_reproductive', 'validate', input, options=options)

__all__ = ['calculate', 'compare', 'deconvolute', 'infer', 'plot', 'predict', 'read', 'report', 'validate']

__candidate_adapter_operations__ = set(['calculate', 'compare', 'deconvolute', 'infer', 'plot', 'predict', 'read', 'report', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
