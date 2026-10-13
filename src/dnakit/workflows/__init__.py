"""Strict local workflows, manifests, and auditable execution helpers."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.workflows.manifest import (
    RunManifestBuilder,
    artifact_from_path,
    load_manifest,
    save_manifest,
)
from dnakit.workflows.runner import (
    ProgressCallback,
    WorkflowProgress,
    WorkflowRunResult,
    WorkflowStepResult,
    run_workflow,
)
from dnakit.workflows.schema import (
    LoadedWorkflow,
    WorkflowErrorPolicy,
    WorkflowInput,
    WorkflowLimits,
    WorkflowOperation,
    WorkflowSpec,
    WorkflowStep,
    load_workflow,
)

from .execution import resource_path, run, run_support, scripts
from .inheritance import simulate_inheritance
from .reporting import write_wes_report

__all__ = [
    "LoadedWorkflow",
    "ProgressCallback",
    "RunManifestBuilder",
    "WorkflowErrorPolicy",
    "WorkflowInput",
    "WorkflowLimits",
    "WorkflowOperation",
    "WorkflowProgress",
    "WorkflowRunResult",
    "WorkflowSpec",
    "WorkflowStep",
    "WorkflowStepResult",
    "artifact_from_path",
    "load_manifest",
    "load_workflow",
    "resource_path",
    "run",
    "run_support",
    "run_workflow",
    "save_manifest",
    "scripts",
    "simulate_inheritance",
    "write_wes_report",
]

# Optional-backend exports (see docs/api/features) 
from .backends import (
    check_wdl_inputs,
    list_galaxy_histories,
    list_galaxy_tools,
    list_galaxy_workflows,
    list_nfcore_processes,
    parse_wdl_document,
    plan_snakemake_jobs,
    read_cwl_workflow,
    read_snakefile,
)

__all__ = [*__all__, "check_wdl_inputs", "list_galaxy_histories", "list_galaxy_tools", "list_galaxy_workflows", "list_nfcore_processes", "parse_wdl_document", "plan_snakemake_jobs", "read_cwl_workflow", "read_snakefile"]

_register_public_verbose(globals())
