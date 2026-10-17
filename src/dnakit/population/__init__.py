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

_register_public_verbose(globals())
