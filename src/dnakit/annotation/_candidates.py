"""Local candidate adapters for the annotation capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def annotate(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'annotate', input, options=options)

def assign_gene_ontology(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'assign_gene_ontology', input, options=options)

def assign_pathways(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'assign_pathways', input, options=options)

def assign_protein_domains(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'assign_protein_domains', input, options=options)

def build_gene_model(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'build_gene_model', input, options=options)

def check_annotation_consistency(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'check_annotation_consistency', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'compare', input, options=options)

def convert_annotation_coordinates(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'convert_annotation_coordinates', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'export', input, options=options)

def extract(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'extract', input, options=options)

def find_homologous_genes(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'find_homologous_genes', input, options=options)

def index(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'index', input, options=options)

def make_intergenic_regions(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'make_intergenic_regions', input, options=options)

def make_promoter_regions(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'make_promoter_regions', input, options=options)

def merge_gene_models(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'merge_gene_models', input, options=options)

def parse(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'parse', input, options=options)

def summarize(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'summarize', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True) -> object:
    return run_candidate_operation('annotation', 'validate', input, options=options)

__all__ = ['annotate', 'assign_gene_ontology', 'assign_pathways', 'assign_protein_domains', 'build_gene_model', 'check_annotation_consistency', 'compare', 'convert_annotation_coordinates', 'export', 'extract', 'find_homologous_genes', 'index', 'make_intergenic_regions', 'make_promoter_regions', 'merge_gene_models', 'parse', 'summarize', 'validate']

from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
