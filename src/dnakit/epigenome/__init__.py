"""Local candidate adapters for the epigenome capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation
from dnakit.epigenome.external import extract_methylation_methyldackel

def align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'align', input, options=options)

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'annotate', input, options=options)

def build_chromatin_state_segments(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'build_chromatin_state_segments', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'calculate', input, options=options)

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'call', input, options=options)

def find(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'find', input, options=options)

def link_enhancers_to_genes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'link_enhancers_to_genes', input, options=options)

def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'merge', input, options=options)

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'normalize', input, options=options)

def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'plot', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'predict', input, options=options)

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'qc', input, options=options)

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('epigenome', 'read', input, options=options)

__all__ = ['align', 'annotate', 'build_chromatin_state_segments', 'calculate', 'call', 'extract_methylation_methyldackel', 'find', 'link_enhancers_to_genes', 'merge', 'normalize', 'plot', 'predict', 'qc', 'read']

__candidate_adapter_operations__ = set(['align', 'annotate', 'build_chromatin_state_segments', 'calculate', 'call', 'find', 'link_enhancers_to_genes', 'merge', 'normalize', 'plot', 'predict', 'qc', 'read'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    balance_contact_matrix,
    build_contact_matrix,
    build_extended_tag_pileup,
    build_fragment_pileup,
    calculate_coverage,
    calculate_expected_contacts,
    calculate_insulation_score,
    calculate_non_conversion_rate,
    calculate_site_binomial_pvalues,
    call_fragment_peaks,
    count_bam_reads_per_bin,
    deduplicate_pairs,
    estimate_fragment_lengths,
    fetch_contact_matrix,
    filter_allc_sites,
    read_pairs_header,
    summarize_bigwig_bins,
    summarize_pairs,
)

__all__ = [*__all__, "balance_contact_matrix", "build_contact_matrix", "build_extended_tag_pileup", "build_fragment_pileup", "calculate_coverage", "calculate_expected_contacts", "calculate_insulation_score", "calculate_non_conversion_rate", "calculate_site_binomial_pvalues", "call_fragment_peaks", "count_bam_reads_per_bin", "deduplicate_pairs", "estimate_fragment_lengths", "fetch_contact_matrix", "filter_allc_sites", "read_pairs_header", "summarize_bigwig_bins", "summarize_pairs"]

# Full-pass exports (see docs/api/features)
from .MACS3_functions import (
    build_paired_fragment_pileup,
    build_shifted_read_pileup,
    calculate_allele_model_likelihood,
    calculate_genotype_quality,
    filter_summits_by_peakiness,
    filter_summits_by_valleys,
    find_signal_maxima,
    guess_tag_file_format,
    hmmratac_fragment_weights,
    max_duplicate_tags,
    merge_pileup_arrays,
    refine_peak_summit,
)

from .cooler_functions import (
    coarsen_contact_matrix,
    fetch_chrom_sizes,
    load_fasta_lengths,
    merge_contact_matrices,
    parse_genomic_region,
    query_bed_intervals,
    read_chrom_sizes,
    read_cooler_chroms,
    read_cooler_pixels,
    rename_cooler_chromosomes,
    zoomify_contact_matrix,
)

from .cooltools_functions import (
    align_windows_to_bins,
    assign_features_to_view,
    calculate_compartments,
    calculate_directionality_index,
    calculate_expected_trans,
    calculate_logbin_expected,
    calculate_saddle,
    calculate_trans_compartments,
    calculate_virtual_4c,
    dense_cis_eigenvectors,
    digitize_track,
    expected_from_dense_matrix,
    iterative_correction_dense,
    make_chromosome_view,
    observed_over_expected_matrix,
    pileup_features,
    read_expected_table,
    read_viewframe_bed,
    rearrange_contact_matrix,
    sample_contacts,
)

from .deeptools_functions import (
    bedgraph_to_bigwig,
    bigwig_chromosome_sizes,
    chop_regions_for_profile,
    estimate_bam_scale_factors,
    estimate_count_size_factors,
    fingerprint_jensen_shannon,
    gc_fraction_from_twobit,
    map_twobit_to_bam_chromosomes,
    match_chromosome_name,
    parse_deeptools_region,
)

from .methylpy_functions import (
    benjamini_hochberg_methylation_cutoffs,
    convert_fastq_to_bisulfite,
    count_bisulfite_basecalls,
    decode_read_c_positions,
    encode_read_c_positions,
    expand_methylation_context,
    filter_bam_by_mch_level,
    index_allc_table,
    read_allc_index,
    read_fasta_chromosome,
)

from .pairtools_functions import (
    assign_genomic_regions,
    assign_haplotype_phase,
    compute_contact_scaling,
    find_read_mismatches,
    find_restriction_fragment,
    make_pairs_header,
    normalize_scaling_curve,
    parse_sam_alignments,
    read_chrom_order,
    read_pairs_table,
    subset_pairs_header,
)

__all__ = [*__all__, "build_paired_fragment_pileup", "build_shifted_read_pileup", "calculate_allele_model_likelihood", "calculate_genotype_quality", "filter_summits_by_peakiness", "filter_summits_by_valleys", "find_signal_maxima", "guess_tag_file_format", "hmmratac_fragment_weights", "max_duplicate_tags", "merge_pileup_arrays", "refine_peak_summit", "coarsen_contact_matrix", "fetch_chrom_sizes", "load_fasta_lengths", "merge_contact_matrices", "parse_genomic_region", "query_bed_intervals", "read_chrom_sizes", "read_cooler_chroms", "read_cooler_pixels", "rename_cooler_chromosomes", "zoomify_contact_matrix", "align_windows_to_bins", "assign_features_to_view", "calculate_compartments", "calculate_directionality_index", "calculate_expected_trans", "calculate_logbin_expected", "calculate_saddle", "calculate_trans_compartments", "calculate_virtual_4c", "dense_cis_eigenvectors", "digitize_track", "expected_from_dense_matrix", "iterative_correction_dense", "make_chromosome_view", "observed_over_expected_matrix", "pileup_features", "read_expected_table", "read_viewframe_bed", "rearrange_contact_matrix", "sample_contacts", "bedgraph_to_bigwig", "bigwig_chromosome_sizes", "chop_regions_for_profile", "estimate_bam_scale_factors", "estimate_count_size_factors", "fingerprint_jensen_shannon", "gc_fraction_from_twobit", "map_twobit_to_bam_chromosomes", "match_chromosome_name", "parse_deeptools_region", "benjamini_hochberg_methylation_cutoffs", "convert_fastq_to_bisulfite", "count_bisulfite_basecalls", "decode_read_c_positions", "encode_read_c_positions", "expand_methylation_context", "filter_bam_by_mch_level", "index_allc_table", "read_allc_index", "read_fasta_chromosome", "assign_genomic_regions", "assign_haplotype_phase", "compute_contact_scaling", "find_read_mismatches", "find_restriction_fragment", "make_pairs_header", "normalize_scaling_curve", "parse_sam_alignments", "read_chrom_order", "read_pairs_table", "subset_pairs_header"]

# Review-pass exports (see docs/api/features)
from .MACS3_review_functions import (
    digest_hmmratac_fragments,
    extract_hmmratac_training_data,
    hmmratac_proba_to_bedgraph,
    hmmratac_states_path,
    load_hmmratac_model,
    predict_hmmratac_states,
    save_hmmratac_model,
    train_hmmratac_model,
    write_hmmratac_accessible_regions,
    write_hmmratac_states_bed,
)

from .cooler_review_functions import (
    create_single_cell_cooler,
    infer_cooler_binsize,
    list_cooler_groups,
    list_scool_cells,
)

from .cooltools_review_functions import (
    adaptive_coarsegrain_matrix,
    align_track_to_cooler_bins,
    assign_bins_to_regions,
    assign_paired_features_to_view,
    call_chromatin_dots,
    dense_trans_eigenvectors,
    expand_features_to_bins,
    expected_full_table,
    make_track_view,
    mask_track_bad_bins,
    observed_over_expected_pixels,
    pairwise_block_sums,
    pairwise_diagonal_sums,
    replace_cis_with_decoys,
    saddle_compartment_strength,
    smooth_contact_decay,
)

from .deeptools_review_functions import (
    bigwig_bedgraph_tiles,
    count_reads_per_gc_bins,
    fingerprint_synthetic_jsd,
    tabulate_gc_fragment_counts,
)

from .methylpy_review_functions import (
    add_cytosine_info_to_dmrs,
    add_methylation_levels_to_dmrs,
    collapse_dmr_windows,
    dmr_residual_cutoff,
    estimate_dmr_fdr_cutoff,
    filter_allc_by_pvalue_cutoffs,
    merge_allc_tables,
    merge_dms_to_dmrs,
    run_allc_binomial_test,
    split_allc_by_chromosome,
    split_allc_into_chunks,
    split_mpileup_chunks,
)

from .pairtools_review_functions import (
    chromsizes_from_sam_header,
    filter_pairs_by_coverage,
    mark_sam_duplicate,
    sam_to_pairs,
)

__all__ = [*__all__, "digest_hmmratac_fragments", "extract_hmmratac_training_data", "hmmratac_proba_to_bedgraph", "hmmratac_states_path", "load_hmmratac_model", "predict_hmmratac_states", "save_hmmratac_model", "train_hmmratac_model", "write_hmmratac_accessible_regions", "write_hmmratac_states_bed", "create_single_cell_cooler", "infer_cooler_binsize", "list_cooler_groups", "list_scool_cells", "adaptive_coarsegrain_matrix", "align_track_to_cooler_bins", "assign_bins_to_regions", "assign_paired_features_to_view", "call_chromatin_dots", "dense_trans_eigenvectors", "expand_features_to_bins", "expected_full_table", "make_track_view", "mask_track_bad_bins", "observed_over_expected_pixels", "pairwise_block_sums", "pairwise_diagonal_sums", "replace_cis_with_decoys", "saddle_compartment_strength", "smooth_contact_decay", "bigwig_bedgraph_tiles", "count_reads_per_gc_bins", "fingerprint_synthetic_jsd", "tabulate_gc_fragment_counts", "add_cytosine_info_to_dmrs", "add_methylation_levels_to_dmrs", "collapse_dmr_windows", "dmr_residual_cutoff", "estimate_dmr_fdr_cutoff", "filter_allc_by_pvalue_cutoffs", "merge_allc_tables", "merge_dms_to_dmrs", "run_allc_binomial_test", "split_allc_by_chromosome", "split_allc_into_chunks", "split_mpileup_chunks", "chromsizes_from_sam_header", "filter_pairs_by_coverage", "mark_sam_duplicate", "sam_to_pairs"]

# Family-pass exports (see docs/api/features)
from .MACS3_family_functions import (
    build_counted_fragment_pileup,
    build_fragment_end_pileup,
    build_hmmratac_fragment_pileup,
    calculate_allele_specific_genotype_quality,
    calculate_allele_specific_model_likelihood,
)

from .cooltools_family_functions import (
    compute_region_contact_scaling,
    pairwise_region_sums,
)

from .deeptools_family_functions import (
    chop_regions_from_middle,
)

from .methylpy_family_functions import (
    convert_paired_fastq_to_bisulfite,
)

from .pairtools_family_functions import (
    deduplicate_pairs_online,
    normalize_scaling_curve_unity,
)

__all__ = [*__all__, "build_counted_fragment_pileup", "build_fragment_end_pileup", "build_hmmratac_fragment_pileup", "calculate_allele_specific_genotype_quality", "calculate_allele_specific_model_likelihood", "compute_region_contact_scaling", "pairwise_region_sums", "chop_regions_from_middle", "convert_paired_fastq_to_bisulfite", "deduplicate_pairs_online", "normalize_scaling_curve_unity"]

_register_public_verbose(globals())
