"""Local candidate adapters for the workflow security capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def benchmark(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'benchmark', input, options=options)

def cache(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'cache', input, options=options)

def cancel(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'cancel', input, options=options)

def check_dependencies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'check_dependencies', input, options=options)

def check_external_tools(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'check_external_tools', input, options=options)

def check_numerical_reproducibility(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'check_numerical_reproducibility', input, options=options)

def checkpoint(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'checkpoint', input, options=options)

def cleanup_temp_workspace(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'cleanup_temp_workspace', input, options=options)

def clear(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'clear', input, options=options)

def create(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'create', input, options=options)

def decrypt_genomic_file(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'decrypt_genomic_file', input, options=options)

def encrypt_genomic_file(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'encrypt_genomic_file', input, options=options)

def estimate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'estimate', input, options=options)

def monitor(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'monitor', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'read', input, options=options)

def record_provenance(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'record_provenance', input, options=options)

def remove_phi_from_metadata(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'remove_phi_from_metadata', input, options=options)

def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'report', input, options=options)

def resume(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'resume', input, options=options)

def retry(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'retry', input, options=options)

def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'run', input, options=options)

def save_config(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'save_config', input, options=options)

def set_random_seed(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'set_random_seed', input, options=options)

def submit_cloud_job(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'submit_cloud_job', input, options=options)

def submit_hpc_job(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'submit_hpc_job', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'validate', input, options=options)

def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('workflow_security', 'write', input, options=options)

__all__ = ['benchmark', 'cache', 'cancel', 'check_dependencies', 'check_external_tools', 'check_numerical_reproducibility', 'checkpoint', 'cleanup_temp_workspace', 'clear', 'create', 'decrypt_genomic_file', 'encrypt_genomic_file', 'estimate', 'monitor', 'read', 'record_provenance', 'remove_phi_from_metadata', 'report', 'resume', 'retry', 'run', 'save_config', 'set_random_seed', 'submit_cloud_job', 'submit_hpc_job', 'validate', 'write']

__candidate_adapter_operations__ = set(['benchmark', 'cache', 'cancel', 'check_dependencies', 'check_external_tools', 'check_numerical_reproducibility', 'checkpoint', 'cleanup_temp_workspace', 'clear', 'create', 'decrypt_genomic_file', 'encrypt_genomic_file', 'estimate', 'monitor', 'read', 'record_provenance', 'remove_phi_from_metadata', 'report', 'resume', 'retry', 'run', 'save_config', 'set_random_seed', 'submit_cloud_job', 'submit_hpc_job', 'validate', 'write'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
