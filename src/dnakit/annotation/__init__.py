"""Reference-based gene annotation, variant annotation, and population frequencies."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation import _candidates as _candidate_module
from dnakit.annotation import external as _external
from dnakit.annotation import gene_models as _gene_models
from dnakit.annotation import genes as _genes
from dnakit.annotation import variants as _variants
from dnakit.annotation.local import (
    convert_annotation_coordinates as convert_annotation_coordinates_local,
    make_intergenic_regions as make_intergenic_regions_local,
    make_promoter_regions as make_promoter_regions_local,
    merge_gene_models as merge_gene_models_local,
)
from dnakit.core import ProviderResult as _ProviderResult

_candidate_all = tuple(_candidate_module.__all__)
globals().update({name: getattr(_candidate_module, name) for name in _candidate_all})

for _module, _names in (
    (_gene_models, ("GeneAnnotationConfig", "GeneAnnotationResult", "GeneRecord", "GeneSegment")),
    (_genes, ("annotate_genes",)),
    (_variants, (
        "annotate_rsid_vep", "annotate_variant_vep", "get_clinvar_significance",
        "get_clinvar_variant", "get_dbsnp_frequencies", "get_dbsnp_variant",
        "get_gnomad_population_frequencies", "get_gnomad_variant", "recode_variant",
        "search_clinvar_variants", "search_gnomad_variants",
    )),
    (_external, ("annotate_bakta", "annotate_prokka", "predict_checkm2")),
):
    globals().update({name: getattr(_module, name) for name in _names})
ProviderResult = _ProviderResult

__all__ = list(dict.fromkeys([
    "GeneAnnotationConfig",
    "GeneAnnotationResult",
    "GeneRecord",
    "GeneSegment",
    "ProviderResult",
    "annotate_genes",
    "annotate_bakta",
    "annotate_prokka",
    "annotate_rsid_vep",
    "annotate_variant_vep",
    "get_clinvar_significance",
    "get_clinvar_variant",
    "get_dbsnp_frequencies",
    "get_dbsnp_variant",
    "get_gnomad_population_frequencies",
    "get_gnomad_variant",
    "recode_variant",
    "predict_checkm2",
    "search_clinvar_variants",
    "search_gnomad_variants",
    "convert_annotation_coordinates",
    "make_intergenic_regions",
    "make_promoter_regions",
    "merge_gene_models",
    *_candidate_all,
]))

convert_annotation_coordinates = convert_annotation_coordinates_local
make_intergenic_regions = make_intergenic_regions_local
make_promoter_regions = make_promoter_regions_local
merge_gene_models = merge_gene_models_local

__candidate_adapter_operations__ = set(_candidate_all) - {
    "convert_annotation_coordinates",
    "make_intergenic_regions",
    "make_promoter_regions",
    "merge_gene_models",
}

# Optional-backend exports (see docs/api/features) 
from .backends import (
    closest_bioframe_intervals,
    closest_pybedtools_intervals,
    complement_pyranges1_intervals,
    count_gffutils_featuretypes,
    count_ncls_overlaps,
    coverage_bioframe_intervals,
    coverage_polars_bio_intervals,
    coverage_pybedtools_intervals,
    create_gffutils_database,
    extract_pybedtools_sequences,
    find_ncls_containing_intervals,
    flank_genomicranges_intervals,
    gaps_gtars_regions,
    intersect_genomicranges_intervals,
    intersect_pyranges1_intervals,
    jaccard_gtars_regions,
    list_gffutils_children,
    merge_bioframe_intervals,
    merge_polars_bio_intervals,
    merge_pybedtools_intervals,
    nearest_polars_bio_intervals,
    nearest_pyranges1_intervals,
    query_ncls_overlaps,
    read_gtfparse_table,
    read_pybigwig_intervals,
    read_pybigwig_values,
    reduce_genomicranges_intervals,
    select_gtfparse_records,
    summarize_gtfparse_annotation,
    summarize_pybigwig_bed_intervals,
    tokenize_gtars_regions,
)

__all__ = [*__all__, "closest_bioframe_intervals", "closest_pybedtools_intervals", "complement_pyranges1_intervals", "count_gffutils_featuretypes", "count_ncls_overlaps", "coverage_bioframe_intervals", "coverage_polars_bio_intervals", "coverage_pybedtools_intervals", "create_gffutils_database", "extract_pybedtools_sequences", "find_ncls_containing_intervals", "flank_genomicranges_intervals", "gaps_gtars_regions", "intersect_genomicranges_intervals", "intersect_pyranges1_intervals", "jaccard_gtars_regions", "list_gffutils_children", "merge_bioframe_intervals", "merge_polars_bio_intervals", "merge_pybedtools_intervals", "nearest_polars_bio_intervals", "nearest_pyranges1_intervals", "query_ncls_overlaps", "read_gtfparse_table", "read_pybigwig_intervals", "read_pybigwig_values", "reduce_genomicranges_intervals", "select_gtfparse_records", "summarize_gtfparse_annotation", "summarize_pybigwig_bed_intervals", "tokenize_gtars_regions"]

# Full-pass exports (see docs/api/features)
from .bioframe_functions import (
    binnify_bioframe_chromosomes,
    cluster_bioframe_intervals,
    expand_bioframe_intervals,
    overlap_bioframe_intervals,
    parse_bioframe_region,
    subtract_bioframe_intervals,
)

from .genomicranges_functions import (
    read_gtf_ranges,
    slide_interval_windows,
    split_interval_bins,
    upstream_downstream_regions,
)

from .gffutils_functions import (
    inspect_gffutils_source,
    parse_gffutils_feature_line,
    to_bed12_gffutils_transcript,
    tss_gffutils_intervals,
)

from .gtars_functions import (
    consensus_regions,
    digest_fasta_sequences,
    fasta_index_records,
    lola_enrichment,
    normalize_variant_allele,
    parse_hgvs_variant,
    partition_region_counts,
    region_dinucleotide_frequencies,
    region_gc_content,
    sequence_digest,
    signal_matrix_summary,
    tokenize_fragments,
    vrs_allele_id,
)

from .gtfparse_functions import (
    expand_gtfparse_attributes,
    parse_gtfparse_fixed_columns,
    write_gtfparse_table,
)

from .polars_bio_functions import (
    alignment_depth,
    cluster_intervals,
    complement_intervals,
    count_interval_overlaps,
    fastqc_summary,
    overlap_intervals,
    read_alignment_records,
    read_bed_intervals,
    read_contact_pairs,
    read_cooler_matrix,
    read_fasta_sequences,
    read_fastq_reads,
    read_gene_annotation,
    read_stockholm_alignments,
    read_variant_calls,
    subtract_intervals,
    write_fasta_sequences,
    write_vcf_variants,
)

from .pyBigWig_functions import (
    read_bigwig_header,
    summarize_bigwig_bins,
    write_bigwig_intervals,
)

from .pybedtools_functions import (
    five_prime_pybedtools_coordinates,
    midpoint_pybedtools_coordinates,
    tag_pybedtools_bedpe,
    write_pybedtools_chromsizes,
)

from .pyranges1_functions import (
    count_pyranges1_overlaps,
    jaccard_pyranges1_intervals,
    parse_pyranges1_string,
    tile_pyranges1_genome,
)

__all__ = [*__all__, "binnify_bioframe_chromosomes", "cluster_bioframe_intervals", "expand_bioframe_intervals", "overlap_bioframe_intervals", "parse_bioframe_region", "subtract_bioframe_intervals", "read_gtf_ranges", "slide_interval_windows", "split_interval_bins", "upstream_downstream_regions", "inspect_gffutils_source", "parse_gffutils_feature_line", "to_bed12_gffutils_transcript", "tss_gffutils_intervals", "consensus_regions", "digest_fasta_sequences", "fasta_index_records", "lola_enrichment", "normalize_variant_allele", "parse_hgvs_variant", "partition_region_counts", "region_dinucleotide_frequencies", "region_gc_content", "sequence_digest", "signal_matrix_summary", "tokenize_fragments", "vrs_allele_id", "expand_gtfparse_attributes", "parse_gtfparse_fixed_columns", "write_gtfparse_table", "alignment_depth", "cluster_intervals", "complement_intervals", "count_interval_overlaps", "fastqc_summary", "overlap_intervals", "read_alignment_records", "read_bed_intervals", "read_contact_pairs", "read_cooler_matrix", "read_fasta_sequences", "read_fastq_reads", "read_gene_annotation", "read_stockholm_alignments", "read_variant_calls", "subtract_intervals", "write_fasta_sequences", "write_vcf_variants", "read_bigwig_header", "summarize_bigwig_bins", "write_bigwig_intervals", "five_prime_pybedtools_coordinates", "midpoint_pybedtools_coordinates", "tag_pybedtools_bedpe", "write_pybedtools_chromsizes", "count_pyranges1_overlaps", "jaccard_pyranges1_intervals", "parse_pyranges1_string", "tile_pyranges1_genome"]

# Review-pass exports (see docs/api/features)
from .bioframe_review_functions import (
    assign_bioframe_views,
    check_bioframe_cataloged,
    check_bioframe_contained,
    check_bioframe_covering,
    check_bioframe_overlapping,
    check_bioframe_sorted,
    check_bioframe_tiling,
    check_bioframe_viewframe,
    complement_bioframe_intervals,
    count_bioframe_overlaps,
    describe_bioframe_assembly,
    digest_bioframe_fasta,
    format_bioframe_region,
    frac_gc_bioframe_intervals,
    frac_gene_coverage_bioframe,
    frac_mapped_bioframe_intervals,
    list_bioframe_assemblies,
    make_bioframe_chromarms,
    mark_bioframe_runs,
    merge_bioframe_runs,
    pair_bioframe_intervals_by_distance,
    read_bioframe_alignments,
    read_bioframe_bigbed,
    read_bioframe_chromsizes,
    read_bioframe_tabix,
    read_bioframe_table,
    resolve_bioframe_region,
    sanitize_bioframe_intervals,
    select_bioframe_region,
    sequence_gc_fraction,
    setdiff_bioframe_intervals,
    sort_bioframe_intervals,
    summarize_bioframe_fasta,
    trim_bioframe_intervals,
    validate_bioframe_bed,
    write_bioframe_bed,
)

from .gffutils_review_functions import (
    asinterval_gffutils_feature,
    bin_gffutils_interval,
    canonical_gffutils_transcript_sequences,
    features_gffutils_to_bedtool,
    gffutils_line_to_seqfeature,
    infer_gffutils_attribute_dialect,
    merge_gffutils_features,
    read_gffutils_records,
    sanitize_gffutils_file,
    seqfeature_to_gffutils_line,
)

from .gtars_review_functions import (
    lola_check_universe,
    lola_redefine_user_sets,
    lola_restricted_universe,
    read_region_tokens,
    vrs_location_digest,
    write_region_tokens,
)

from .gtfparse_review_functions import (
    create_missing_gtfparse_features,
    expand_gtfparse_attribute_strings,
)

from .polars_bio_review_functions import (
    describe_genomic_file,
    query_genomic_tables,
    read_alignment_sequences,
    read_bigbed_intervals,
    read_cram_records,
    read_file_header,
    write_alignment_file,
    write_fastq_reads,
)

from .pybedtools_review_functions import (
    bed_to_gff_pybedtools,
    center_pybedtools_coordinates,
    cis_trans_pybedtools_interactions,
    filter_pybedtools_intervals_by_length,
    gff_to_bed_pybedtools,
    normalize_pybedtools_scores_to_length,
    overlap_pybedtools_bp,
    random_jaccard_pybedtools_intervals,
    scale_pybedtools_bedgraph_values,
    three_prime_pybedtools_coordinates,
)

from .pyranges1_review_functions import (
    calculate_pyranges1_cds_frame,
    calculate_pyranges1_forbes,
    calculate_pyranges1_mcc,
    calculate_pyranges1_relative_distance,
    concatenate_pyranges1_intervals,
    extend_pyranges1_orfs,
    read_pyranges1_bigwig,
    read_pyranges1_gff,
    read_pyranges1_gff3,
    read_pyranges1_gtf_full,
)

__all__ = [*__all__, "assign_bioframe_views", "check_bioframe_cataloged", "check_bioframe_contained", "check_bioframe_covering", "check_bioframe_overlapping", "check_bioframe_sorted", "check_bioframe_tiling", "check_bioframe_viewframe", "complement_bioframe_intervals", "count_bioframe_overlaps", "describe_bioframe_assembly", "digest_bioframe_fasta", "format_bioframe_region", "frac_gc_bioframe_intervals", "frac_gene_coverage_bioframe", "frac_mapped_bioframe_intervals", "list_bioframe_assemblies", "make_bioframe_chromarms", "mark_bioframe_runs", "merge_bioframe_runs", "pair_bioframe_intervals_by_distance", "read_bioframe_alignments", "read_bioframe_bigbed", "read_bioframe_chromsizes", "read_bioframe_tabix", "read_bioframe_table", "resolve_bioframe_region", "sanitize_bioframe_intervals", "select_bioframe_region", "sequence_gc_fraction", "setdiff_bioframe_intervals", "sort_bioframe_intervals", "summarize_bioframe_fasta", "trim_bioframe_intervals", "validate_bioframe_bed", "write_bioframe_bed", "asinterval_gffutils_feature", "bin_gffutils_interval", "canonical_gffutils_transcript_sequences", "features_gffutils_to_bedtool", "gffutils_line_to_seqfeature", "infer_gffutils_attribute_dialect", "merge_gffutils_features", "read_gffutils_records", "sanitize_gffutils_file", "seqfeature_to_gffutils_line", "lola_check_universe", "lola_redefine_user_sets", "lola_restricted_universe", "read_region_tokens", "vrs_location_digest", "write_region_tokens", "create_missing_gtfparse_features", "expand_gtfparse_attribute_strings", "describe_genomic_file", "query_genomic_tables", "read_alignment_sequences", "read_bigbed_intervals", "read_cram_records", "read_file_header", "write_alignment_file", "write_fastq_reads", "bed_to_gff_pybedtools", "center_pybedtools_coordinates", "cis_trans_pybedtools_interactions", "filter_pybedtools_intervals_by_length", "gff_to_bed_pybedtools", "normalize_pybedtools_scores_to_length", "overlap_pybedtools_bp", "random_jaccard_pybedtools_intervals", "scale_pybedtools_bedgraph_values", "three_prime_pybedtools_coordinates", "calculate_pyranges1_cds_frame", "calculate_pyranges1_forbes", "calculate_pyranges1_mcc", "calculate_pyranges1_relative_distance", "concatenate_pyranges1_intervals", "extend_pyranges1_orfs", "read_pyranges1_bigwig", "read_pyranges1_gff", "read_pyranges1_gff3", "read_pyranges1_gtf_full"]

# Family-pass exports (see docs/api/features)
from .bioframe_family_functions import (
    select_bioframe_overlap_rows,
)

from .gtars_family_functions import (
    load_fasta_sequence_collection,
    sequence_digest_record,
)

__all__ = [*__all__, "select_bioframe_overlap_rows", "load_fasta_sequence_collection", "sequence_digest_record"]

_register_public_verbose(globals())
