"""Local candidate adapters for the singlecell spatial capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def analyze_crispr_perturbations(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'analyze_crispr_perturbations', input, options=options)

def annotate_cells(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'annotate_cells', input, options=options)

def annotate_cells_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'annotate_cells_spatial', input, options=options)

def build_cell_gene_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'build_cell_gene_matrix', input, options=options)

def build_cell_peak_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'build_cell_peak_matrix', input, options=options)

def calculate_cell_cycle_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_cell_cycle_score', input, options=options)

def calculate_cell_interactions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_cell_interactions', input, options=options)

def calculate_cell_interactions_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_cell_interactions_spatial', input, options=options)

def calculate_markers(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_markers', input, options=options)

def calculate_markers_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_markers_spatial', input, options=options)

def calculate_perturbation_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'calculate_perturbation_effect', input, options=options)

def cluster(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'cluster', input, options=options)

def cluster_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'cluster_spatial', input, options=options)

def correct_umis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'correct_umis', input, options=options)

def correct_umis_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'correct_umis_spatial', input, options=options)

def deconvolute_spots(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'deconvolute_spots', input, options=options)

def deconvolute_spots_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'deconvolute_spots_spatial', input, options=options)

def demultiplex(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'demultiplex', input, options=options)

def demultiplex_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'demultiplex_spatial', input, options=options)

def extract_barcodes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'extract_barcodes', input, options=options)

def extract_barcodes_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'extract_barcodes_spatial', input, options=options)

def extract_umis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'extract_umis', input, options=options)

def extract_umis_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'extract_umis_spatial', input, options=options)

def find_differential_features(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'find_differential_features', input, options=options)

def find_differential_features_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'find_differential_features_spatial', input, options=options)

def infer_cell_type_from_markers(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_cell_type_from_markers', input, options=options)

def infer_cnv(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_cnv', input, options=options)

def infer_cnv_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_cnv_spatial', input, options=options)

def infer_trajectory(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_trajectory', input, options=options)

def infer_trajectory_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_trajectory_spatial', input, options=options)

def infer_velocity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_velocity', input, options=options)

def infer_velocity_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'infer_velocity_spatial', input, options=options)

def integrate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'integrate', input, options=options)

def integrate_single_cell_batches(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'integrate_single_cell_batches', input, options=options)

def integrate_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'integrate_spatial', input, options=options)

def map_cells(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'map_cells', input, options=options)

def map_cells_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'map_cells_spatial', input, options=options)

def map_single_cell_to_reference(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'map_single_cell_to_reference', input, options=options)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'normalize', input, options=options)

def normalize_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'normalize_spatial', input, options=options)

def plot_embedding(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'plot_embedding', input, options=options)

def plot_embedding_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'plot_embedding_spatial', input, options=options)

def plot_heatmap(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'plot_heatmap', input, options=options)

def plot_heatmap_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'plot_heatmap_spatial', input, options=options)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'qc', input, options=options)

def qc_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'qc_spatial', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'read', input, options=options)

def read_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'read_spatial', input, options=options)

def remove_ambient(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'remove_ambient', input, options=options)

def remove_ambient_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'remove_ambient_spatial', input, options=options)

def remove_doublets(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'remove_doublets', input, options=options)

def remove_doublets_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'remove_doublets_spatial', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'validate', input, options=options)

def validate_spatial(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('singlecell_spatial', 'validate_spatial', input, options=options)

__all__ = ['analyze_crispr_perturbations', 'annotate_cells', 'annotate_cells_spatial', 'build_cell_gene_matrix', 'build_cell_peak_matrix', 'calculate_cell_cycle_score', 'calculate_cell_interactions', 'calculate_cell_interactions_spatial', 'calculate_markers', 'calculate_markers_spatial', 'calculate_perturbation_effect', 'cluster', 'cluster_spatial', 'correct_umis', 'correct_umis_spatial', 'deconvolute_spots', 'deconvolute_spots_spatial', 'demultiplex', 'demultiplex_spatial', 'extract_barcodes', 'extract_barcodes_spatial', 'extract_umis', 'extract_umis_spatial', 'find_differential_features', 'find_differential_features_spatial', 'infer_cell_type_from_markers', 'infer_cnv', 'infer_cnv_spatial', 'infer_trajectory', 'infer_trajectory_spatial', 'infer_velocity', 'infer_velocity_spatial', 'integrate', 'integrate_single_cell_batches', 'integrate_spatial', 'map_cells', 'map_cells_spatial', 'map_single_cell_to_reference', 'normalize', 'normalize_spatial', 'plot_embedding', 'plot_embedding_spatial', 'plot_heatmap', 'plot_heatmap_spatial', 'qc', 'qc_spatial', 'read', 'read_spatial', 'remove_ambient', 'remove_ambient_spatial', 'remove_doublets', 'remove_doublets_spatial', 'validate', 'validate_spatial']

__candidate_adapter_operations__ = set(['analyze_crispr_perturbations', 'annotate_cells', 'annotate_cells_spatial', 'build_cell_gene_matrix', 'build_cell_peak_matrix', 'calculate_cell_cycle_score', 'calculate_cell_interactions', 'calculate_cell_interactions_spatial', 'calculate_markers', 'calculate_markers_spatial', 'calculate_perturbation_effect', 'cluster', 'cluster_spatial', 'correct_umis', 'correct_umis_spatial', 'deconvolute_spots', 'deconvolute_spots_spatial', 'demultiplex', 'demultiplex_spatial', 'extract_barcodes', 'extract_barcodes_spatial', 'extract_umis', 'extract_umis_spatial', 'find_differential_features', 'find_differential_features_spatial', 'infer_cell_type_from_markers', 'infer_cnv', 'infer_cnv_spatial', 'infer_trajectory', 'infer_trajectory_spatial', 'infer_velocity', 'infer_velocity_spatial', 'integrate', 'integrate_single_cell_batches', 'integrate_spatial', 'map_cells', 'map_cells_spatial', 'map_single_cell_to_reference', 'normalize', 'normalize_spatial', 'plot_embedding', 'plot_embedding_spatial', 'plot_heatmap', 'plot_heatmap_spatial', 'qc', 'qc_spatial', 'read', 'read_spatial', 'remove_ambient', 'remove_ambient_spatial', 'remove_doublets', 'remove_doublets_spatial', 'validate', 'validate_spatial'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
