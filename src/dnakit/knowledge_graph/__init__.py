"""Local candidate adapters for the knowledge graph capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def add_graph_relationship(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'add_graph_relationship', input, options=options)

def aggregate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'aggregate', input, options=options)

def apply_access_policy(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'apply_access_policy', input, options=options)

def build_data_lineage(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'build_data_lineage', input, options=options)

def build_knowledge_graph(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'build_knowledge_graph', input, options=options)

def compare_database_versions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'compare_database_versions', input, options=options)

def create(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'create', input, options=options)

def deidentify_genomic_metadata(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'deidentify_genomic_metadata', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'export', input, options=options)

def link(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'link', input, options=options)

def map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'map', input, options=options)

def query_graph_neighbors(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'query_graph_neighbors', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'read', input, options=options)

def search(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'search', input, options=options)

def snapshot_database_version(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'snapshot_database_version', input, options=options)

def update(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'update', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'validate', input, options=options)

def version(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('knowledge_graph', 'version', input, options=options)

__all__ = ['add_graph_relationship', 'aggregate', 'apply_access_policy', 'build_data_lineage', 'build_knowledge_graph', 'compare_database_versions', 'create', 'deidentify_genomic_metadata', 'export', 'link', 'map', 'query_graph_neighbors', 'read', 'search', 'snapshot_database_version', 'update', 'validate', 'version']

__candidate_adapter_operations__ = set(['add_graph_relationship', 'aggregate', 'apply_access_policy', 'build_data_lineage', 'build_knowledge_graph', 'compare_database_versions', 'create', 'deidentify_genomic_metadata', 'export', 'link', 'map', 'query_graph_neighbors', 'read', 'search', 'snapshot_database_version', 'update', 'validate', 'version'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
