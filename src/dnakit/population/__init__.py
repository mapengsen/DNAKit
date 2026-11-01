"""Population-genetics functions: native statistics on dosage genotype matrices plus
optional-backend wrappers.

Genotype input for the generic functions: ``{'genotypes': samples x variants alternate-allele
dosages (0/1/2, fractional imputed dosages allowed where noted, null = missing), 'sample_ids':
[...], 'variant_ids': [...], 'populations': [...]}`` or just the matrix.  Every function returns
the DNAKit candidate envelope.
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit.gwas._impl_common import finish as _finish
from dnakit.gwas._impl_common import get_options as _get_options

from . import _impl_pop as _pop

_M = "population"


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate per-variant allele frequencies, heterozygosity or pairwise linkage
    disequilibrium from a genotype matrix.

    input: genotype mapping/matrix. options: statistic = 'allele_frequency' (default:
    alternate-allele frequency, MAF, called samples, missing rate; per population when
    'populations' is given), 'heterozygosity' (hard calls: Ho, He = 2pq, unbiased He, F = 1 -
    Ho/He) or 'ld' (composite-LD r^2, D and D' for all variant pairs; max_variants=100).
    """

    return _finish(_M, "calculate", input, options, _pop.calculate)


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare populations by Fst (Weir & Cockerham 1984 or Hudson), globally, per variant and
    pairwise.

    input: genotype mapping with 'populations' (hard calls 0/1/2). options: method =
    'weir_cockerham' (default; theta = sum(a)/sum(a+b+c) over variants, all populations
    jointly plus a pairwise table) or 'hudson' (pairwise ratio of sums, Bhatia et al. 2013).
    Variants with undefined estimates are None and excluded from the global ratio.
    """

    return _finish(_M, "compare", input, options, _pop.compare)


def estimate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate genomic relationships (VanRaden GRM or kinship), per-sample inbreeding F, or
    genotype principal components.

    input: genotype mapping/matrix. options: quantity = 'grm' (default; VanRaden method 1),
    'kinship' (GRM/2), 'inbreeding' (PLINK --het method-of-moments F per sample, hard calls)
    or 'pca' (frequency-standardised PCA; components=2). 'grm', 'kinship' and 'pca' need
    numpy and report not_implemented without it.
    """

    return _finish(_M, "estimate", input, options, _pop.estimate)


def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a genotype matrix as TSV/CSV, PLINK text PED+MAP, or an allele-frequency table.

    input: genotype mapping/matrix. options: format = 'tsv' (default), 'csv', 'ped_map' (hard
    calls; ref/alt lists in the input name the alleles, otherwise A/G; path is the file stem)
    or 'frq' (SNP, A1, A2, MAF, NCHROBS); path (existing files are not overwritten unless
    overwrite=True). The text is always returned.
    """

    return _finish(_M, "export", input, options, _pop.export)


def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer genetic ancestry or population structure (not implemented: needs a reference panel
    or a model-fitting program).
    """

    return candidate_not_implemented(
        _M,
        "infer",
        input,
        options=_get_options(options),
        requires=(
            "ancestry reference panel (e.g. 1000 Genomes) or an ADMIXTURE/STRUCTURE-type program"
        ),
        message=(
            "Ancestry/structure inference needs labelled reference genotypes or a model-fitting "
            "backend, neither of which is bundled; "
            "use estimate(quantity='pca') for a PCA of your own samples."
        ),
    )


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a population-genetics plot (MAF histogram, PCA scatter or per-variant Fst) as an SVG
    string.

    input: genotype mapping/matrix ('populations' colours the PCA). options: kind =
    'maf_histogram' (default; bins=10), 'pca' (numpy needed) or 'fst' (method as compare);
    width=640, height=420; path (write the SVG; no overwrite unless overwrite=True).
    """

    return _finish(_M, "plot", input, options, _pop.plot)


def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate population genotypes (Balding-Nichols structured populations) or neutral
    Wright-Fisher allele-frequency drift.

    input: None or a mapping of the options (overrides options). model='balding_nichols'
    (default): n_populations=2, n_samples_per_population=50, n_variants=100, fst=0.1,
    maf_min=0.1; p_k ~ Beta(p(1-F)/F, (1-p)(1-F)/F), genotypes Binomial(2, p_k).
    model='wright_fisher': population_size=100 (diploid), generations=50,
    initial_frequency=0.5, replicates=100. Both are deterministic given seed (default 0). Not
    a coalescent/forward simulator such as msprime or SLiM.
    """

    return _finish(_M, "simulate", input, options, _pop.simulate)


def test(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test Hardy-Weinberg equilibrium per variant (exact or chi-square) or allele-frequency
    differentiation between populations.

    input: genotype mapping (hard calls; 'populations' for differentiation). options: test =
    'hwe' (default) or 'differentiation' (2 x k allele-count chi-square); method = 'exact'
    (Wigginton et al. 2005, default) or 'chisq' for HWE. Each tested variant gets Bonferroni
    (p*m) and Benjamini-Hochberg q-values.
    """

    return _finish(_M, "test", input, options, _pop.test)


__all__ = ["calculate", "compare", "estimate", "export", "infer", "plot", "simulate", "test"]

__candidate_adapter_operations__ = set(
    ["calculate", "compare", "estimate", "export", "infer", "plot", "simulate", "test"]
)


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


# Coverage-gap exports (see paper/write_paper/figure/fig2/task_coverage/coverage_gap/waves/CANONICAL_MODULES.md)
from .allel_gap_functions import (
    fig_voight_painting,
    plot_haplotype_frequencies,
    plot_joint_sfs,
    plot_joint_sfs_folded,
    plot_joint_sfs_folded_scaled,
    plot_joint_sfs_scaled,
    plot_moving_haplotype_frequencies,
    plot_pairwise_distance,
    plot_pairwise_ld,
    plot_sfs,
    plot_sfs_folded,
    plot_sfs_folded_scaled,
    plot_sfs_scaled,
    plot_variant_locator,
    plot_voight_painting,
    roh_poissonhmm,
)
from .allel_verify_functions import (
    allel_vcf_to_zarr,
)
from .dadi_gap2_functions import (
    plot_3d_spectrum_mayavi,
    vourlaki_visual_validation_png,
)
from .dadi_gap_functions import (
    coverage_depth_distribution,
    dfe_mixture_spectrum,
    fit_demes_graph_parameters,
    plot_dfe_bivariate,
    plot_dfe_bivariate_point_pos,
    plot_sfs_1d,
    plot_sfs_1d_comparison,
    plot_sfs_2d,
    plot_sfs_2d_comparison,
    plot_sfs_2d_meta_residual,
    plot_sfs_2d_residual,
    plot_sfs_3d,
    plot_sfs_3d_comparison,
    plot_sfs_3d_pairwise,
    plot_trispectrum,
    plot_trispectrum_comparison,
    plot_two_locus_comparison,
    plot_two_locus_spectrum,
)
from .dadi_verify_functions import (
    dfe_mixture_point_pos_spectrum,
    optimize_demes_graph_parameters,
    plot_triallelic_spectrum,
    plot_triallelic_spectrum_comparison,
    plot_two_locus_spectrum_3d,
)
from .tskit_batch0_review_functions import (
    draw_tree_text,
)
from .tskit_cli_functions import (
    export_tree_file_edges_table,
    export_tree_file_fasta,
    export_tree_file_individuals_table,
    export_tree_file_migrations_table,
    export_tree_file_mutations_table,
    export_tree_file_nodes_table,
    export_tree_file_populations_table,
    export_tree_file_provenances_table,
    export_tree_file_sites_table,
    export_tree_file_vcf,
    tree_file_info,
    tree_file_trees,
)
from .tskit_gap_functions import (
    draw_tree_html,
    draw_tree_sequence_html,
    draw_tree_sequence_svg,
    draw_variant_html,
    tskit_random_nucleotides,
)

__all__ = [*__all__, "fig_voight_painting", "plot_haplotype_frequencies", "plot_joint_sfs", "plot_joint_sfs_folded", "plot_joint_sfs_folded_scaled", "plot_joint_sfs_scaled", "plot_moving_haplotype_frequencies", "plot_pairwise_distance", "plot_pairwise_ld", "plot_sfs", "plot_sfs_folded", "plot_sfs_folded_scaled", "plot_sfs_scaled", "plot_variant_locator", "plot_voight_painting", "roh_poissonhmm", "allel_vcf_to_zarr", "plot_3d_spectrum_mayavi", "vourlaki_visual_validation_png", "coverage_depth_distribution", "dfe_mixture_spectrum", "fit_demes_graph_parameters", "plot_dfe_bivariate", "plot_dfe_bivariate_point_pos", "plot_sfs_1d", "plot_sfs_1d_comparison", "plot_sfs_2d", "plot_sfs_2d_comparison", "plot_sfs_2d_meta_residual", "plot_sfs_2d_residual", "plot_sfs_3d", "plot_sfs_3d_comparison", "plot_sfs_3d_pairwise", "plot_trispectrum", "plot_trispectrum_comparison", "plot_two_locus_comparison", "plot_two_locus_spectrum", "dfe_mixture_point_pos_spectrum", "optimize_demes_graph_parameters", "plot_triallelic_spectrum", "plot_triallelic_spectrum_comparison", "plot_two_locus_spectrum_3d", "draw_tree_text", "export_tree_file_edges_table", "export_tree_file_fasta", "export_tree_file_individuals_table", "export_tree_file_migrations_table", "export_tree_file_mutations_table", "export_tree_file_nodes_table", "export_tree_file_populations_table", "export_tree_file_provenances_table", "export_tree_file_sites_table", "export_tree_file_vcf", "tree_file_info", "tree_file_trees", "draw_tree_html", "draw_tree_sequence_html", "draw_tree_sequence_svg", "draw_variant_html", "tskit_random_nucleotides"]

# Public integration of the reviewed P00–P45 operation packs.
from .biopython_p19_cluster_functions import (
    cluster_centroids,
    cluster_distance,
    cluster_pca,
    cut_cluster_tree,
    hierarchical_cluster,
    kmeans_cluster,
    kmedoids_cluster,
    profile_distance_matrix,
    remove_genepop_loci,
    remove_genepop_populations,
    som_cluster,
    sort_cluster_tree,
    split_genepop_by_locus,
    split_genepop_by_population,
)
from .dadi_p22_dfe_triallele_functions import (
    alternate_mutation_triallelic_spectrum,
    dfe_integrated_spectrum,
    triallelic_neutral_equilibrium,
    triallelic_segregating_sites,
    triallelic_spectrum_transform,
)
from .dadi_p22_spectrum_functions import (
    calculate_sfs_statistic,
    expected_sfs_extrapolated,
    expected_sfs_low_pass,
    fragment_snp_data,
    gatk_no_call_probability,
    integrate_diffusion_phi,
    log_likelihood_derivatives,
    log_spectrum,
    mask_spectrum_corners,
    probability_enough_individuals_covered,
    project_spectrum,
    reorder_spectrum_populations,
    resample_spectrum,
    simulate_gatk_calling_spectrum,
    simulate_low_coverage_reads,
    spectrum_from_simulation_file,
    subsample_genotypes,
    subset_spectrum_populations,
    unfold_spectrum,
)
from .dadi_p22_twolocus_functions import (
    condition_two_locus_spectrum,
    read_multilocus_spectrum_file,
    sample_two_locus_density,
    scale_two_locus_spectrum,
    two_locus_ld_per_bin,
    two_locus_log_likelihood,
    two_locus_log_likelihood_over_rho_bins,
    two_locus_marginal_spectra,
    two_locus_mean_r2,
    two_locus_spectrum_transform,
    write_multilocus_spectrum_file,
)
from .dadi_p23_phi_functions import (
    phi_admix_into_population,
    phi_create_admixed_population,
    phi_one_population,
    phi_select_populations,
    phi_split_population,
    scramble_spectrum_pop_ids,
    spectrum_from_phi,
    spectrum_from_snp_data,
)
from .tskit_p31_functions import (
    tskit_alignments,
    tskit_breakpoints,
    tskit_compute_mutation_parents,
    tskit_compute_mutation_times,
    tskit_concatenate,
    tskit_deduplicate_sites,
    tskit_delete_sites,
    tskit_extend_haplotypes,
    tskit_genealogical_nearest_neighbours,
    tskit_general_stat,
    tskit_genetic_relatedness,
    tskit_genetic_relatedness_vector,
    tskit_genetic_relatedness_weighted,
    tskit_genotype_matrix,
    tskit_haplotypes,
    tskit_ibd_segments,
    tskit_ld_r2,
    tskit_link_ancestors,
    tskit_list_mutations,
    tskit_list_sites,
    tskit_map_mutations,
    tskit_map_to_vcf_model,
    tskit_mean_descendants,
    tskit_pair_coalescence,
    tskit_read_hapmap,
    tskit_recombination_map_lookup,
    tskit_remove_older_genealogy,
    tskit_sample_nodes_by_ploidy,
    tskit_sample_set_statistic,
    tskit_segregating_sites,
    tskit_select_intervals,
    tskit_shift_coordinates,
    tskit_simplify,
    tskit_split_edges,
    tskit_subset,
    tskit_trait_statistic,
    tskit_trim_genome,
    tskit_union,
    tskit_variants,
)

__all__ = [
    *__all__,
    "alternate_mutation_triallelic_spectrum",
    "calculate_sfs_statistic",
    "cluster_centroids",
    "cluster_distance",
    "cluster_pca",
    "condition_two_locus_spectrum",
    "cut_cluster_tree",
    "dfe_integrated_spectrum",
    "expected_sfs_extrapolated",
    "expected_sfs_low_pass",
    "fragment_snp_data",
    "gatk_no_call_probability",
    "hierarchical_cluster",
    "integrate_diffusion_phi",
    "kmeans_cluster",
    "kmedoids_cluster",
    "log_likelihood_derivatives",
    "log_spectrum",
    "mask_spectrum_corners",
    "phi_admix_into_population",
    "phi_create_admixed_population",
    "phi_one_population",
    "phi_select_populations",
    "phi_split_population",
    "probability_enough_individuals_covered",
    "profile_distance_matrix",
    "project_spectrum",
    "read_multilocus_spectrum_file",
    "remove_genepop_loci",
    "remove_genepop_populations",
    "reorder_spectrum_populations",
    "resample_spectrum",
    "sample_two_locus_density",
    "scale_two_locus_spectrum",
    "scramble_spectrum_pop_ids",
    "simulate_gatk_calling_spectrum",
    "simulate_low_coverage_reads",
    "som_cluster",
    "sort_cluster_tree",
    "spectrum_from_phi",
    "spectrum_from_simulation_file",
    "spectrum_from_snp_data",
    "split_genepop_by_locus",
    "split_genepop_by_population",
    "subsample_genotypes",
    "subset_spectrum_populations",
    "triallelic_neutral_equilibrium",
    "triallelic_segregating_sites",
    "triallelic_spectrum_transform",
    "tskit_alignments",
    "tskit_breakpoints",
    "tskit_compute_mutation_parents",
    "tskit_compute_mutation_times",
    "tskit_concatenate",
    "tskit_deduplicate_sites",
    "tskit_delete_sites",
    "tskit_extend_haplotypes",
    "tskit_genealogical_nearest_neighbours",
    "tskit_general_stat",
    "tskit_genetic_relatedness",
    "tskit_genetic_relatedness_vector",
    "tskit_genetic_relatedness_weighted",
    "tskit_genotype_matrix",
    "tskit_haplotypes",
    "tskit_ibd_segments",
    "tskit_ld_r2",
    "tskit_link_ancestors",
    "tskit_list_mutations",
    "tskit_list_sites",
    "tskit_map_mutations",
    "tskit_map_to_vcf_model",
    "tskit_mean_descendants",
    "tskit_pair_coalescence",
    "tskit_read_hapmap",
    "tskit_recombination_map_lookup",
    "tskit_remove_older_genealogy",
    "tskit_sample_nodes_by_ploidy",
    "tskit_sample_set_statistic",
    "tskit_segregating_sites",
    "tskit_select_intervals",
    "tskit_shift_coordinates",
    "tskit_simplify",
    "tskit_split_edges",
    "tskit_subset",
    "tskit_trait_statistic",
    "tskit_trim_genome",
    "tskit_union",
    "tskit_variants",
    "two_locus_ld_per_bin",
    "two_locus_log_likelihood",
    "two_locus_log_likelihood_over_rho_bins",
    "two_locus_marginal_spectra",
    "two_locus_mean_r2",
    "two_locus_spectrum_transform",
    "unfold_spectrum",
    "write_multilocus_spectrum_file",
]

_register_public_verbose(globals())
