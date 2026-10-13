"""Local candidate adapters for the transcriptomics capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'align', input, options=options)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'annotate', input, options=options)

def calculate_differential_expression(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'calculate_differential_expression', input, options=options)

def calculate_differential_splicing(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'calculate_differential_splicing', input, options=options)

def calculate_translation_efficiency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'calculate_translation_efficiency', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'compare', input, options=options)

def deconvolute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'deconvolute', input, options=options)

def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'detect', input, options=options)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'normalize', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'plot', input, options=options)

def predict_protein_function(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'predict_protein_function', input, options=options)

def predict_protein_localization(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'predict_protein_localization', input, options=options)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'qc', input, options=options)

def quantify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'quantify', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'read', input, options=options)

def run_featurecounts(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'run_featurecounts', input, options=options)

def run_kallisto_quantification(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'run_kallisto_quantification', input, options=options)

def run_salmon_quantification(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('transcriptomics', 'run_salmon_quantification', input, options=options)

__all__ = ['align', 'annotate', 'calculate_differential_expression', 'calculate_differential_splicing', 'calculate_translation_efficiency', 'compare', 'deconvolute', 'detect', 'normalize', 'plot', 'predict_protein_function', 'predict_protein_localization', 'qc', 'quantify', 'read', 'run_featurecounts', 'run_kallisto_quantification', 'run_salmon_quantification']

__candidate_adapter_operations__ = set(['align', 'annotate', 'calculate_differential_expression', 'calculate_differential_splicing', 'calculate_translation_efficiency', 'compare', 'deconvolute', 'detect', 'normalize', 'plot', 'predict_protein_function', 'predict_protein_localization', 'qc', 'quantify', 'read', 'run_featurecounts', 'run_kallisto_quantification', 'run_salmon_quantification'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    aggregate_pseudobulk,
    calculate_cell_qc_metrics,
    calculate_gene_qc_metrics,
    filter_genes_by_min_cells,
    rank_marker_genes,
    read_anndata_summary,
    select_highly_variable_genes,
    summarize_mudata,
    train_scvi_model,
)

__all__ = [*__all__, "aggregate_pseudobulk", "calculate_cell_qc_metrics", "calculate_gene_qc_metrics", "filter_genes_by_min_cells", "rank_marker_genes", "read_anndata_summary", "select_highly_variable_genes", "summarize_mudata", "train_scvi_model"]

_register_public_verbose(globals())
