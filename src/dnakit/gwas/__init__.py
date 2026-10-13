"""Local candidate adapters for the gwas capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def annotate_lead_snps(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'annotate_lead_snps', input, options=options)

def apply_genomic_control(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'apply_genomic_control', input, options=options)

def build_credible_set(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'build_credible_set', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'calculate', input, options=options)

def calibrate_prs_risk(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'calibrate_prs_risk', input, options=options)

def clump_and_threshold_prs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'clump_and_threshold_prs', input, options=options)

def compute_locus_plot_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'compute_locus_plot_data', input, options=options)

def compute_manhattan_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'compute_manhattan_data', input, options=options)

def compute_qq_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'compute_qq_data', input, options=options)

def conditional_gwas(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'conditional_gwas', input, options=options)

def define_association_loci(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'define_association_loci', input, options=options)

def estimate_genetic_correlation(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'estimate_genetic_correlation', input, options=options)

def estimate_snp_heritability(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'estimate_snp_heritability', input, options=options)

def evaluate_prs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'evaluate_prs', input, options=options)

def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'export', input, options=options)

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'filter', input, options=options)

def fixed_effect_meta_analysis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'fixed_effect_meta_analysis', input, options=options)

def harmonize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'harmonize', input, options=options)

def joint_gwas_analysis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'joint_gwas_analysis', input, options=options)

def merge_nearby_significant_hits(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'merge_nearby_significant_hits', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'plot', input, options=options)

def random_effect_meta_analysis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'random_effect_meta_analysis', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'read', input, options=options)

def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'run', input, options=options)

def select_lead_snps(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'select_lead_snps', input, options=options)

def simulate_gwas_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'simulate_gwas_data', input, options=options)

def standardize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'standardize', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'validate', input, options=options)

def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('gwas', 'write', input, options=options)

__all__ = ['annotate_lead_snps', 'apply_genomic_control', 'build_credible_set', 'calculate', 'calibrate_prs_risk', 'clump_and_threshold_prs', 'compute_locus_plot_data', 'compute_manhattan_data', 'compute_qq_data', 'conditional_gwas', 'define_association_loci', 'estimate_genetic_correlation', 'estimate_snp_heritability', 'evaluate_prs', 'export', 'filter', 'fixed_effect_meta_analysis', 'harmonize', 'joint_gwas_analysis', 'merge_nearby_significant_hits', 'plot', 'random_effect_meta_analysis', 'read', 'run', 'select_lead_snps', 'simulate_gwas_data', 'standardize', 'validate', 'write']

__candidate_adapter_operations__ = set(['annotate_lead_snps', 'apply_genomic_control', 'build_credible_set', 'calculate', 'calibrate_prs_risk', 'clump_and_threshold_prs', 'compute_locus_plot_data', 'compute_manhattan_data', 'compute_qq_data', 'conditional_gwas', 'define_association_loci', 'estimate_genetic_correlation', 'estimate_snp_heritability', 'evaluate_prs', 'export', 'filter', 'fixed_effect_meta_analysis', 'harmonize', 'joint_gwas_analysis', 'merge_nearby_significant_hits', 'plot', 'random_effect_meta_analysis', 'read', 'run', 'select_lead_snps', 'simulate_gwas_data', 'standardize', 'validate', 'write'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
