"""Local candidate adapters for the external adapters capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def build_command(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'build_command', input, options=options)

def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'check', input, options=options)

def convert_input(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'convert_input', input, options=options)

def convert_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'convert_output', input, options=options)

def make_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'make_report', input, options=options)

def parse_log(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'parse_log', input, options=options)

def parse_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'parse_output', input, options=options)

def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'run', input, options=options)

def validate_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('external_adapters', 'validate_output', input, options=options)

__all__ = ['build_command', 'check', 'convert_input', 'convert_output', 'make_report', 'parse_log', 'parse_output', 'run', 'validate_output']

__candidate_adapter_operations__ = set(['build_command', 'check', 'convert_input', 'convert_output', 'make_report', 'parse_log', 'parse_output', 'run', 'validate_output'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
