"""Local candidate adapters for the variants capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation
from dnakit._verbosity import register_public as _register_public_verbose

from .external import (
    annotate_vcf_ensembl_vep,
    annotate_vcf_vcfanno,
    call_cutesv,
    call_deepvariant,
    call_delly,
    call_freebayes,
    call_gatk_haplotypecaller,
    call_sniffles2,
    filter_vcf_vcftools,
    load_gemini_database,
    merge_structural_variants_survivor,
    query_gemini_database,
)
from .statistics import (
    calculate_allele_frequency,
    calculate_heterozygosity,
    calculate_missingness,
    calculate_transition_transversion,
    calculate_variant_density,
    calculate_minor_allele_frequency as _calculate_minor_allele_frequency_local,
    compare as _compare_local,
    count as _count_local,
    filter_variants_by_allele_frequency as _filter_variants_by_allele_frequency_local,
    filter_variants_by_depth as _filter_variants_by_depth_local,
    filter_variants_by_genotype_quality as _filter_variants_by_genotype_quality_local,
    filter_variants_by_missingness as _filter_variants_by_missingness_local,
    filter_variants_by_quality as _filter_variants_by_quality_local,
    parse_vcf_records as _parse_vcf_records_local,
    split_multiallelic as _split_multiallelic_local,
)
from .backends import (
    calculate_fst,
    calculate_nucleotide_diversity,
    calculate_observed_heterozygosity,
    calculate_pca,
    calculate_site_frequency_spectrum,
    calculate_tajimas_d,
    calculate_watterson_theta,
    calculate_windowed_ld,
    fetch_vcf_region,
    read_vcf,
)


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate', input, options=options)

def annotate_coding_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_coding_effect', input, options=options)

def annotate_conservation_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_conservation_score', input, options=options)

def annotate_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_mt', input, options=options)

def annotate_pathogenicity_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_pathogenicity_score', input, options=options)

def annotate_population_frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_population_frequency', input, options=options)

def annotate_regulatory_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_regulatory_effect', input, options=options)

def annotate_splice_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_splice_effect', input, options=options)

def annotate_variant_consequence(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'annotate_variant_consequence', input, options=options)

def calculate_copy_number_profile(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_copy_number_profile', input, options=options)

def calculate_genotype_likelihood(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_genotype_likelihood', input, options=options)

def calculate_haplotype_frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_haplotype_frequency', input, options=options)

def calculate_hwe_pvalue(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_hwe_pvalue', input, options=options)

def calculate_identity_by_descent(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_identity_by_descent', input, options=options)

def calculate_ld_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_ld_matrix', input, options=options)

def calculate_runs_of_homozygosity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'calculate_runs_of_homozygosity', input, options=options)

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'call', input, options=options)

def call_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'call_mt', input, options=options)

def compare_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'compare_mt', input, options=options)

def compare_vcfs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'compare_vcfs', input, options=options)

def construct_haplotypes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'construct_haplotypes', input, options=options)

def count_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'count_mt', input, options=options)

def detect_breakpoints(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'detect_breakpoints', input, options=options)

def detect_low_frequency_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'detect_low_frequency_variants', input, options=options)

def detect_mobile_element_insertions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'detect_mobile_element_insertions', input, options=options)

def detect_mosaic_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'detect_mosaic_variants', input, options=options)

def detect_repeat_expansions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'detect_repeat_expansions', input, options=options)

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'filter', input, options=options)

def filter_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'filter_mt', input, options=options)

def genotype_structural_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'genotype_structural_variants', input, options=options)

def impute_missing_genotypes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'impute_missing_genotypes', input, options=options)

def left_align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'left_align', input, options=options)

def left_align_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'left_align_mt', input, options=options)

def left_align_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'left_align_variants', input, options=options)

def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'merge', input, options=options)

def merge_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'merge_mt', input, options=options)

def merge_vcf_files(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'merge_vcf_files', input, options=options)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'normalize', input, options=options)

def normalize_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'normalize_mt', input, options=options)

def normalize_vcf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'normalize_vcf', input, options=options)

def phase_genotypes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'phase_genotypes', input, options=options)

def plot_density(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'plot_density', input, options=options)

def plot_density_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'plot_density_mt', input, options=options)

def prioritize_candidate_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'prioritize_candidate_variants', input, options=options)

def split_multiallelic_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'split_multiallelic_mt', input, options=options)

def split_multiallelic_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'split_multiallelic_variants', input, options=options)

def summarize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'summarize', input, options=options)

def summarize_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'summarize_mt', input, options=options)

__all__ = ['annotate', 'annotate_coding_effect', 'annotate_conservation_score', 'annotate_mt', 'annotate_pathogenicity_score', 'annotate_population_frequency', 'annotate_regulatory_effect', 'annotate_splice_effect', 'annotate_variant_consequence', 'annotate_vcf_ensembl_vep', 'annotate_vcf_vcfanno', 'calculate_allele_frequency', 'call_cutesv', 'call_deepvariant', 'call_delly', 'calculate_copy_number_profile', 'calculate_genotype_likelihood', 'calculate_haplotype_frequency', 'calculate_heterozygosity', 'calculate_hwe_pvalue', 'calculate_identity_by_descent', 'calculate_ld_matrix', 'calculate_minor_allele_frequency', 'calculate_missingness', 'calculate_runs_of_homozygosity', 'calculate_transition_transversion', 'calculate_variant_density', 'call', 'call_freebayes', 'call_gatk_haplotypecaller', 'call_sniffles2', 'call_mt', 'compare', 'compare_mt', 'compare_vcfs', 'construct_haplotypes', 'count', 'count_mt', 'detect_breakpoints', 'detect_low_frequency_variants', 'detect_mobile_element_insertions', 'detect_mosaic_variants', 'detect_repeat_expansions', 'filter', 'filter_mt', 'filter_variants_by_allele_frequency', 'filter_variants_by_depth', 'filter_variants_by_genotype_quality', 'filter_variants_by_missingness', 'filter_variants_by_quality', 'filter_vcf_vcftools', 'genotype_structural_variants', 'impute_missing_genotypes', 'left_align', 'left_align_mt', 'left_align_variants', 'load_gemini_database', 'merge', 'merge_mt', 'merge_vcf_files', 'merge_structural_variants_survivor', 'normalize', 'normalize_mt', 'normalize_vcf', 'parse_vcf_records', 'phase_genotypes', 'plot_density', 'plot_density_mt', 'prioritize_candidate_variants', 'query_gemini_database', 'split_multiallelic', 'split_multiallelic_mt', 'split_multiallelic_variants', 'summarize', 'summarize_mt', 'read_vcf', 'fetch_vcf_region', 'calculate_observed_heterozygosity', 'calculate_fst', 'calculate_pca', 'calculate_tajimas_d', 'calculate_watterson_theta', 'calculate_nucleotide_diversity', 'calculate_site_frequency_spectrum', 'calculate_windowed_ld']

__candidate_adapter_operations__ = set(['annotate', 'annotate_coding_effect', 'annotate_conservation_score', 'annotate_mt', 'annotate_pathogenicity_score', 'annotate_population_frequency', 'annotate_regulatory_effect', 'annotate_splice_effect', 'annotate_variant_consequence', 'calculate_allele_frequency', 'calculate_copy_number_profile', 'calculate_genotype_likelihood', 'calculate_haplotype_frequency', 'calculate_heterozygosity', 'calculate_hwe_pvalue', 'calculate_identity_by_descent', 'calculate_ld_matrix', 'calculate_minor_allele_frequency', 'calculate_missingness', 'calculate_runs_of_homozygosity', 'calculate_transition_transversion', 'calculate_variant_density', 'call', 'call_mt', 'compare', 'compare_mt', 'compare_vcfs', 'construct_haplotypes', 'count', 'count_mt', 'detect_breakpoints', 'detect_low_frequency_variants', 'detect_mobile_element_insertions', 'detect_mosaic_variants', 'detect_repeat_expansions', 'filter', 'filter_mt', 'filter_variants_by_allele_frequency', 'filter_variants_by_depth', 'filter_variants_by_genotype_quality', 'filter_variants_by_missingness', 'filter_variants_by_quality', 'genotype_structural_variants', 'impute_missing_genotypes', 'left_align', 'left_align_mt', 'left_align_variants', 'merge', 'merge_mt', 'merge_vcf_files', 'normalize', 'normalize_mt', 'normalize_vcf', 'parse_vcf_records', 'phase_genotypes', 'plot_density', 'plot_density_mt', 'prioritize_candidate_variants', 'split_multiallelic', 'split_multiallelic_mt', 'split_multiallelic_variants', 'summarize', 'summarize_mt'])

# These summaries are implemented locally and do not require a caller,
# annotation database, or optional backend.
__candidate_adapter_operations__ -= {
    'calculate_allele_frequency',
    'calculate_heterozygosity',
    'calculate_missingness',
    'calculate_transition_transversion',
    'calculate_variant_density',
}

# These record-level operations are deterministic and dependency-free. The
# broad wrappers above remain available for unimplemented operations, while
# these names resolve to the tested local implementations.
calculate_minor_allele_frequency = _calculate_minor_allele_frequency_local
compare = _compare_local
count = _count_local
filter_variants_by_allele_frequency = _filter_variants_by_allele_frequency_local
filter_variants_by_depth = _filter_variants_by_depth_local
filter_variants_by_genotype_quality = _filter_variants_by_genotype_quality_local
filter_variants_by_missingness = _filter_variants_by_missingness_local
filter_variants_by_quality = _filter_variants_by_quality_local
parse_vcf_records = _parse_vcf_records_local
split_multiallelic = _split_multiallelic_local
__candidate_adapter_operations__ -= {
    'calculate_minor_allele_frequency',
    'compare',
    'count',
    'filter_variants_by_allele_frequency',
    'filter_variants_by_depth',
    'filter_variants_by_genotype_quality',
    'filter_variants_by_missingness',
    'filter_variants_by_quality',
    'parse_vcf_records',
    'split_multiallelic',
}


# Optional-backend exports (see docs/api/features) 
from .nomenclature import (
    check_hgvs_syntax,
    describe_variant,
    format_hgvs_variant,
    load_vcf_variants,
    parse_hgvs_variant,
)

from .genotype_formats import (
    calculate_bgen_allele_frequency,
    read_bgen_probabilities,
    read_bgen_variants,
    read_plink_bed,
    read_plink_genotypes,
    read_plink_metadata,
    summarize_plink_bed,
    write_plink_bed,
    write_plink_genotypes,
)

__all__ = [*__all__, "check_hgvs_syntax", "describe_variant", "format_hgvs_variant", "load_vcf_variants", "parse_hgvs_variant", "calculate_bgen_allele_frequency", "read_bgen_probabilities", "read_bgen_variants", "read_plink_bed", "read_plink_genotypes", "read_plink_metadata", "summarize_plink_bed", "write_plink_bed", "write_plink_genotypes"]

# Optional-backend exports (hail)
from .hail_backends import (
    calculate_association_tests,
    calculate_hardy_weinberg_pvalues,
    calculate_hwe_normalized_pca,
    calculate_sample_qc,
    calculate_variant_qc,
    import_vcf_matrix_table,
    prune_linkage_disequilibrium,
    summarize_vcf_matrix_table,
)

__all__ = [*__all__, "calculate_association_tests", "calculate_hardy_weinberg_pvalues", "calculate_hwe_normalized_pca", "calculate_sample_qc", "calculate_variant_qc", "import_vcf_matrix_table", "prune_linkage_disequilibrium", "summarize_vcf_matrix_table"]

# Full-pass exports (see docs/api/features)
from .allel_functions import (
    calculate_ehh_decay,
    calculate_expected_heterozygosity,
    calculate_garud_h,
    calculate_hudson_fst_per_site,
    calculate_ihs,
    calculate_inbreeding_coefficient,
    calculate_mean_pairwise_difference,
    calculate_nsl,
    calculate_pairwise_distance,
    calculate_patterson_d,
    calculate_population_branching_statistic,
    calculate_rogers_huff_r,
    calculate_sequence_divergence,
    calculate_windowed_nucleotide_diversity,
    calculate_windowed_tajimas_d,
    calculate_xpehh,
    call_roh_hmm,
    find_mendel_errors,
    phase_trio_by_transmission,
    select_unlinked_variants,
)

from .bed_reader_functions import (
    calculate_genotype_relationship_matrix,
    standardize_genotypes,
    subset_genotype_matrix,
)

from .bgen_reader_functions import (
    create_bgen_metafile,
    read_bgen_allele_expectation,
    read_bgen_dosage,
)

from .cyvcf2_functions import (
    pearson_genotype_correlation,
    summarize_pairwise_relatedness,
    summarize_sample_heterozygosity,
)

from .hgvs_functions import (
    assemble_repeat_blocks,
    build_transcript_cigar,
    detect_repeat_blocks,
    normalize_indel_alleles,
    trim_common_prefix_strings,
    validate_hgvs_type_accession,
)

from .pandas_plink_functions import (
    read_gcta_relationship_matrix,
    read_plink_relationship_matrix,
)

from .pysam_functions import (
    annotate_vcf_from_tabix,
    build_vcf_consensus,
    compress_bgzip_file,
    concatenate_vcf_files,
    fetch_padded_sequence,
    index_tabix_file,
    intersect_vcf_sites,
    parse_genomic_region,
    parse_genomic_regions,
    query_vcf_fields,
    select_vcf_records,
    sort_vcf_records,
)

from .varcode_functions import (
    breakend_alt_sides,
    breakend_insert_bases,
    breakend_mate_position,
    check_vcf_sample_identity,
    classify_breakend_pair,
    compare_sv_calls,
    detect_germline_allele_overlap,
    enumerate_cryptic_exon_candidates,
    enumerate_phase_hypotheses,
    four_gamete_phase_test,
    load_maf_variants,
    normalize_nucleotide_allele,
    parse_symbolic_sv_alt,
    parse_vcf_genotype,
    partition_germline_by_phase,
    score_splice_site,
    trim_shared_allele_flanks,
    variant_overlaps_interval,
)

from .vcfpy_functions import (
    classify_vcf_alt,
    classify_vcf_substitution,
    format_vcf_field_value,
    parse_vcf_breakend,
    parse_vcf_field_value,
    remove_vcf_header_lines,
)

__all__ = [*__all__, "calculate_ehh_decay", "calculate_expected_heterozygosity", "calculate_garud_h", "calculate_hudson_fst_per_site", "calculate_ihs", "calculate_inbreeding_coefficient", "calculate_mean_pairwise_difference", "calculate_nsl", "calculate_pairwise_distance", "calculate_patterson_d", "calculate_population_branching_statistic", "calculate_rogers_huff_r", "calculate_sequence_divergence", "calculate_windowed_nucleotide_diversity", "calculate_windowed_tajimas_d", "calculate_xpehh", "call_roh_hmm", "find_mendel_errors", "phase_trio_by_transmission", "select_unlinked_variants", "calculate_genotype_relationship_matrix", "standardize_genotypes", "subset_genotype_matrix", "create_bgen_metafile", "read_bgen_allele_expectation", "read_bgen_dosage", "pearson_genotype_correlation", "summarize_pairwise_relatedness", "summarize_sample_heterozygosity", "assemble_repeat_blocks", "build_transcript_cigar", "detect_repeat_blocks", "normalize_indel_alleles", "trim_common_prefix_strings", "validate_hgvs_type_accession", "read_gcta_relationship_matrix", "read_plink_relationship_matrix", "annotate_vcf_from_tabix", "build_vcf_consensus", "compress_bgzip_file", "concatenate_vcf_files", "fetch_padded_sequence", "index_tabix_file", "intersect_vcf_sites", "parse_genomic_region", "parse_genomic_regions", "query_vcf_fields", "select_vcf_records", "sort_vcf_records", "breakend_alt_sides", "breakend_insert_bases", "breakend_mate_position", "check_vcf_sample_identity", "classify_breakend_pair", "compare_sv_calls", "detect_germline_allele_overlap", "enumerate_cryptic_exon_candidates", "enumerate_phase_hypotheses", "four_gamete_phase_test", "load_maf_variants", "normalize_nucleotide_allele", "parse_symbolic_sv_alt", "parse_vcf_genotype", "partition_germline_by_phase", "score_splice_site", "trim_shared_allele_flanks", "variant_overlaps_interval", "classify_vcf_alt", "classify_vcf_substitution", "format_vcf_field_value", "parse_vcf_breakend", "parse_vcf_field_value", "remove_vcf_header_lines"]

# Review-pass exports (see docs/api/features)
from .allel_review_functions import (
    build_equally_accessible_windows,
    calculate_between_population_pi,
    calculate_blockwise_patterson_f3,
    calculate_blockwise_patterson_fst,
    calculate_haplotype_diversity,
    calculate_moving_delta_tajimas_d,
    calculate_moving_haplotype_diversity,
    calculate_moving_patterson_f3,
    calculate_moving_patterson_fst,
    calculate_pairwise_dxy,
    calculate_patterson_f2,
    calculate_patterson_f3,
    calculate_patterson_fst,
    calculate_windowed_fixed_difference_density,
    calculate_windowed_patterson_fst,
    calculate_xpnsl,
    count_variants_in_windows,
    fold_joint_site_frequency_spectrum,
    locate_population_fixed_differences,
    locate_population_private_alleles,
    paint_parental_transmission,
    paint_shared_haplotype_prefixes,
    read_vcf_header_summary,
    scale_folded_joint_site_frequency_spectrum,
    scale_folded_site_frequency_spectrum,
    scale_joint_site_frequency_spectrum,
    scale_site_frequency_spectrum,
    standardize_score_by_allele_count,
    vcf_to_hdf5_archive,
    vcf_to_npz_archive,
    windowed_values_per_base,
)

from .hgvs_review_functions import (
    hgvs_interbase_span,
    repeat_notation_for_alleles,
    trim_common_suffix_strings,
)

from .pysam_review_functions import (
    add_bam_read_group,
    annotate_vcf_consequence,
    calculate_bam_read_depth,
    check_alignment_file_integrity,
    check_vcf_sample_concordance,
    clip_amplicon_primers,
    collate_bam_reads,
    concatenate_bam_files,
    count_bam_reads_per_interval,
    decode_sam_flags,
    detect_bcftools_roh,
    extract_fasta_from_bam,
    extract_fastq_from_bam,
    filter_vcf_by_expression,
    fix_bam_mate_fields,
    import_fastq_to_bam,
    index_alignment_file,
    index_fastq_file,
    index_vcf_csi,
    list_read_group_samples,
    mark_bam_duplicates,
    merge_bam_files,
    merge_indexed_vcfs,
    phase_bam_haplotypes,
    pileup_alignment_positions,
    pileup_bam_to_vcf,
    read_bam_header_and_records,
    read_vcf_header,
    recompute_nm_md_tags,
    reconstruct_reference_from_alignments,
    rename_vcf_samples,
    replace_bam_header,
    reset_alignments_to_unmapped,
    select_alignment_records,
    sort_alignment_file,
    split_bam_by_read_group,
    summarize_alignment_stats,
    summarize_amplicon_stats,
    summarize_bam_coverage,
    summarize_bam_index_stats,
    write_fasta_sequence_dictionary,
)

from .varcode_review_functions import (
    breakend_junction_span,
    has_only_standard_nucleotides,
    infer_reference_assembly,
    is_mitochondrial_contig,
    is_purine_base,
    is_ucsc_reference_name,
    join_breakend_layouts,
    pair_breakend_mates,
    project_intervals_onto_layout,
    rewrite_vcf_records,
)

__all__ = [*__all__, "build_equally_accessible_windows", "calculate_between_population_pi", "calculate_blockwise_patterson_f3", "calculate_blockwise_patterson_fst", "calculate_haplotype_diversity", "calculate_moving_delta_tajimas_d", "calculate_moving_haplotype_diversity", "calculate_moving_patterson_f3", "calculate_moving_patterson_fst", "calculate_pairwise_dxy", "calculate_patterson_f2", "calculate_patterson_f3", "calculate_patterson_fst", "calculate_windowed_fixed_difference_density", "calculate_windowed_patterson_fst", "calculate_xpnsl", "count_variants_in_windows", "fold_joint_site_frequency_spectrum", "locate_population_fixed_differences", "locate_population_private_alleles", "paint_parental_transmission", "paint_shared_haplotype_prefixes", "read_vcf_header_summary", "scale_folded_joint_site_frequency_spectrum", "scale_folded_site_frequency_spectrum", "scale_joint_site_frequency_spectrum", "scale_site_frequency_spectrum", "standardize_score_by_allele_count", "vcf_to_hdf5_archive", "vcf_to_npz_archive", "windowed_values_per_base", "hgvs_interbase_span", "repeat_notation_for_alleles", "trim_common_suffix_strings", "add_bam_read_group", "annotate_vcf_consequence", "calculate_bam_read_depth", "check_alignment_file_integrity", "check_vcf_sample_concordance", "clip_amplicon_primers", "collate_bam_reads", "concatenate_bam_files", "count_bam_reads_per_interval", "decode_sam_flags", "detect_bcftools_roh", "extract_fasta_from_bam", "extract_fastq_from_bam", "filter_vcf_by_expression", "fix_bam_mate_fields", "import_fastq_to_bam", "index_alignment_file", "index_fastq_file", "index_vcf_csi", "list_read_group_samples", "mark_bam_duplicates", "merge_bam_files", "merge_indexed_vcfs", "phase_bam_haplotypes", "pileup_alignment_positions", "pileup_bam_to_vcf", "read_bam_header_and_records", "read_vcf_header", "recompute_nm_md_tags", "reconstruct_reference_from_alignments", "rename_vcf_samples", "replace_bam_header", "reset_alignments_to_unmapped", "select_alignment_records", "sort_alignment_file", "split_bam_by_read_group", "summarize_alignment_stats", "summarize_amplicon_stats", "summarize_bam_coverage", "summarize_bam_index_stats", "write_fasta_sequence_dictionary", "breakend_junction_span", "has_only_standard_nucleotides", "infer_reference_assembly", "is_mitochondrial_contig", "is_purine_base", "is_ucsc_reference_name", "join_breakend_layouts", "pair_breakend_mates", "project_intervals_onto_layout", "rewrite_vcf_records"]

# Family-pass exports (see docs/api/features)
from .allel_family_functions import (
    calculate_blockwise_patterson_d,
    calculate_joint_site_frequency_spectrum_from_counts,
    calculate_moving_garud_h,
    calculate_moving_hudson_fst,
    calculate_moving_patterson_d,
    calculate_moving_tajimas_d,
    calculate_moving_weir_cockerham_fst,
    calculate_randomized_pca,
    calculate_rogers_huff_r_between,
    calculate_site_frequency_spectrum_from_counts,
    calculate_weir_cockerham_fst_per_variant,
    calculate_windowed_divergence,
    calculate_windowed_hudson_fst,
    calculate_windowed_watterson_theta,
    calculate_windowed_weir_cockerham_fst,
    phase_trio_step_by_transmission,
    read_vcf_field_table,
    vcf_to_csv_table,
)

from .pysam_family_functions import (
    export_vcf_to_phasing_format,
)

from .varcode_family_functions import (
    load_maf_table,
    read_vcf_table,
)

__all__ = [*__all__, "calculate_blockwise_patterson_d", "calculate_joint_site_frequency_spectrum_from_counts", "calculate_moving_garud_h", "calculate_moving_hudson_fst", "calculate_moving_patterson_d", "calculate_moving_tajimas_d", "calculate_moving_weir_cockerham_fst", "calculate_randomized_pca", "calculate_rogers_huff_r_between", "calculate_site_frequency_spectrum_from_counts", "calculate_weir_cockerham_fst_per_variant", "calculate_windowed_divergence", "calculate_windowed_hudson_fst", "calculate_windowed_watterson_theta", "calculate_windowed_weir_cockerham_fst", "phase_trio_step_by_transmission", "read_vcf_field_table", "vcf_to_csv_table", "export_vcf_to_phasing_format", "load_maf_table", "read_vcf_table"]

_register_public_verbose(globals())
