"""Local candidate adapters for the population capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'calculate', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'compare', input, options=options)

def estimate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'estimate', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'export', input, options=options)

def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'infer', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'plot', input, options=options)

def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'simulate', input, options=options)

def test(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('population', 'test', input, options=options)

__all__ = ['calculate', 'compare', 'estimate', 'export', 'infer', 'plot', 'simulate', 'test']

__candidate_adapter_operations__ = set(['calculate', 'compare', 'estimate', 'export', 'infer', 'plot', 'simulate', 'test'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    build_tree_sequence,
    calculate_diversity,
    calculate_fst,
    calculate_sfs_log_likelihood,
    calculate_site_frequency_spectrum,
    calculate_tajimas_d,
    estimate_sfs_theta,
    expected_neutral_sfs,
    fold_sfs,
    read_tree_sequence,
    write_tree_sequence,
)

__all__ = [*__all__, "build_tree_sequence", "calculate_diversity", "calculate_fst", "calculate_sfs_log_likelihood", "calculate_site_frequency_spectrum", "calculate_tajimas_d", "estimate_sfs_theta", "expected_neutral_sfs", "fold_sfs", "read_tree_sequence", "write_tree_sequence"]

# Full-pass exports (see docs/api/features)
from .dadi_functions import (
    annotate_snps_from_annovar,
    apply_ancestral_misidentification,
    bootstrap_snp_spectra,
    expected_sfs_1d,
    expected_sfs_2d,
    expected_sfs_3d,
    expected_sfs_from_demes,
    expected_sfs_selection,
    fit_demographic_model,
    parameter_uncertainty_fim,
    pool_sfs_populations,
    read_sfs_file,
    read_snp_data,
    read_vcf_data,
    sfs_residuals,
    snp_configuration_counts,
    trinucleotide_frequency_array,
    write_sfs_file,
)

from .tskit_functions import (
    complement_intervals,
    count_sample_topologies,
    export_ms_format,
    export_sample_sequences,
    map_allele_indices,
    random_tree_newick,
    read_plink_fam,
    read_text_table,
    read_text_tree_sequence,
    resolve_polytomies,
    tree_to_newick,
    write_text_tree_sequence,
)

__all__ = [*__all__, "annotate_snps_from_annovar", "apply_ancestral_misidentification", "bootstrap_snp_spectra", "expected_sfs_1d", "expected_sfs_2d", "expected_sfs_3d", "expected_sfs_from_demes", "expected_sfs_selection", "fit_demographic_model", "parameter_uncertainty_fim", "pool_sfs_populations", "read_sfs_file", "read_snp_data", "read_vcf_data", "sfs_residuals", "snp_configuration_counts", "trinucleotide_frequency_array", "write_sfs_file", "complement_intervals", "count_sample_topologies", "export_ms_format", "export_sample_sequences", "map_allele_indices", "random_tree_newick", "read_plink_fam", "read_text_table", "read_text_tree_sequence", "resolve_polytomies", "tree_to_newick", "write_text_tree_sequence"]

# Review-pass exports (see docs/api/features)
from .dadi_review_functions import (
    expected_neutral_sfs_2d,
    fit_demographic_model_to_residual,
    godambe_matrices,
    ms_command_line,
    ms_core_command_2d,
    nested_model_statistic,
    parameter_uncertainty_gim,
    slim_output_data_dict,
    total_substitution_rate,
    trinucleotide_misid_table,
    vcf_subsample_bootstrap_spectra,
)

__all__ = [*__all__, "expected_neutral_sfs_2d", "fit_demographic_model_to_residual", "godambe_matrices", "ms_command_line", "ms_core_command_2d", "nested_model_statistic", "parameter_uncertainty_gim", "slim_output_data_dict", "total_substitution_rate", "trinucleotide_misid_table", "vcf_subsample_bootstrap_spectra"]

# Family-pass exports (see docs/api/features)
from .dadi_family_functions import (
    calculate_model_log_likelihood,
    fit_demographic_model_constrained,
    fit_demographic_model_grid,
    fit_demographic_model_natural_scale,
    fit_demographic_model_nlopt,
    scale_model_sfs_to_data,
)

from .tskit_family_functions import (
    count_tree_topologies,
    export_tree_sequence_nexus,
)

__all__ = [*__all__, "calculate_model_log_likelihood", "fit_demographic_model_constrained", "fit_demographic_model_grid", "fit_demographic_model_natural_scale", "fit_demographic_model_nlopt", "scale_model_sfs_to_data", "count_tree_topologies", "export_tree_sequence_nexus"]

_register_public_verbose(globals())
