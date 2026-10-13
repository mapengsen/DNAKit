"""Local candidate adapters for the replication repair capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'annotate', input, options=options)

def build_recombination_map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'build_recombination_map', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'calculate', input, options=options)

def classify_repair_pathway(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'classify_repair_pathway', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'compare', input, options=options)

def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'detect', input, options=options)

def map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'map', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'plot', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'predict', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('replication_repair', 'simulate', input, options=options)

__all__ = ['annotate', 'build_recombination_map', 'calculate', 'classify_repair_pathway', 'compare', 'detect', 'map', 'plot', 'predict', 'simulate']

__candidate_adapter_operations__ = set(['annotate', 'build_recombination_map', 'calculate', 'classify_repair_pathway', 'compare', 'detect', 'map', 'plot', 'predict', 'simulate'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
