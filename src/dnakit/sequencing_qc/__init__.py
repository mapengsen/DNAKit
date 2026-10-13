"""Local candidate adapters for the sequencing qc capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

from ._external import run_cutadapt as _run_cutadapt_backend
from ._external import run_fastp as _run_fastp_backend
from ._external import run_fastqc as _run_fastqc_backend
from ._external import run_multiqc as _run_multiqc_backend
from .long_reads import (
    filter_long_reads,
    filter_nanopore_reads,
    summarize_long_reads,
    trim_porechop,
)
from .short_reads import run_afterqc, trim_trimmomatic


def assign_sample_barcodes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'assign_sample_barcodes', input, options=options)

def basecall_nanopore_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'basecall_nanopore_reads', input, options=options)

def calculate_adapter_content(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_adapter_content', input, options=options)

def calculate_adapter_content_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_adapter_content_tenx', input, options=options)

def calculate_barcode_quality(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_barcode_quality', input, options=options)

def calculate_barcode_quality_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_barcode_quality_tenx', input, options=options)

def calculate_basecalling_qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_basecalling_qc', input, options=options)

def calculate_basecalling_qc_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_basecalling_qc_tenx', input, options=options)

def calculate_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_contamination', input, options=options)

def calculate_contamination_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_contamination_tenx', input, options=options)

def calculate_duplicate_rate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_duplicate_rate', input, options=options)

def calculate_duplication(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_duplication', input, options=options)

def calculate_duplication_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_duplication_tenx', input, options=options)

def calculate_error_profile(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_error_profile', input, options=options)

def calculate_error_profile_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_error_profile_tenx', input, options=options)

def calculate_gc_bias(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_gc_bias', input, options=options)

def calculate_gc_bias_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_gc_bias_tenx', input, options=options)

def calculate_index_hopping(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_index_hopping', input, options=options)

def calculate_index_hopping_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_index_hopping_tenx', input, options=options)

def calculate_insert_size(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_insert_size', input, options=options)

def calculate_insert_size_distribution(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_insert_size_distribution', input, options=options)

def calculate_insert_size_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_insert_size_tenx', input, options=options)

def calculate_phasing_qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_phasing_qc', input, options=options)

def calculate_phasing_qc_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_phasing_qc_tenx', input, options=options)

def calculate_qscore_summary(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_qscore_summary', input, options=options)

def calculate_quality(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_quality', input, options=options)

def calculate_quality_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_quality_tenx', input, options=options)

def calculate_read_length(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_read_length', input, options=options)

def calculate_read_length_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_read_length_tenx', input, options=options)

def calculate_yield(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_yield', input, options=options)

def calculate_yield_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'calculate_yield_tenx', input, options=options)

def collapse_umi_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'collapse_umi_reads', input, options=options)

def correct_long_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'correct_long_reads', input, options=options)

def correct_umi_errors(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'correct_umi_errors', input, options=options)

def demultiplex_barcoded_reads(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'demultiplex_barcoded_reads', input, options=options)

def detect_cross_sample_contamination(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'detect_cross_sample_contamination', input, options=options)

def detect_platform(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'detect_platform', input, options=options)

def detect_platform_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'detect_platform_tenx', input, options=options)

def estimate_library_complexity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'estimate_library_complexity', input, options=options)

def estimate_read_error_model(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'estimate_read_error_model', input, options=options)

def extract_umi_sequences(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'extract_umi_sequences', input, options=options)

def generate_qc_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'generate_qc_report', input, options=options)

def generate_qc_report_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'generate_qc_report_tenx', input, options=options)

def plot_per_base_quality(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'plot_per_base_quality', input, options=options)

def plot_qc_summary(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'plot_qc_summary', input, options=options)

def plot_qc_summary_tenx(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('sequencing_qc', 'plot_qc_summary_tenx', input, options=options)

def run_cutadapt(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_cutadapt_backend(input, options=options, verbose=verbose)

def run_fastp(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_fastp_backend(input, options=options, verbose=verbose)

def run_fastqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_fastqc_backend(input, options=options, verbose=verbose)

def run_multiqc(
    input: object, *, options: Mapping[str, object] | None = None, verbose: bool = True
) -> object:
    return _run_multiqc_backend(input, options=options, verbose=verbose)

__all__ = ['assign_sample_barcodes', 'basecall_nanopore_reads', 'calculate_adapter_content', 'calculate_adapter_content_tenx', 'calculate_barcode_quality', 'calculate_barcode_quality_tenx', 'calculate_basecalling_qc', 'calculate_basecalling_qc_tenx', 'calculate_contamination', 'calculate_contamination_tenx', 'calculate_duplicate_rate', 'calculate_duplication', 'calculate_duplication_tenx', 'calculate_error_profile', 'calculate_error_profile_tenx', 'calculate_gc_bias', 'calculate_gc_bias_tenx', 'calculate_index_hopping', 'calculate_index_hopping_tenx', 'calculate_insert_size', 'calculate_insert_size_distribution', 'calculate_insert_size_tenx', 'calculate_phasing_qc', 'calculate_phasing_qc_tenx', 'calculate_qscore_summary', 'calculate_quality', 'calculate_quality_tenx', 'calculate_read_length', 'calculate_read_length_tenx', 'calculate_yield', 'calculate_yield_tenx', 'collapse_umi_reads', 'correct_long_reads', 'correct_umi_errors', 'demultiplex_barcoded_reads', 'detect_cross_sample_contamination', 'detect_platform', 'detect_platform_tenx', 'estimate_library_complexity', 'estimate_read_error_model', 'extract_umi_sequences', 'filter_long_reads', 'filter_nanopore_reads', 'generate_qc_report', 'generate_qc_report_tenx', 'plot_per_base_quality', 'plot_qc_summary', 'plot_qc_summary_tenx', 'run_afterqc', 'run_cutadapt', 'run_fastp', 'run_fastqc', 'run_multiqc', 'summarize_long_reads', 'trim_porechop', 'trim_trimmomatic']

__candidate_adapter_operations__ = set(['assign_sample_barcodes', 'basecall_nanopore_reads', 'calculate_adapter_content', 'calculate_adapter_content_tenx', 'calculate_barcode_quality', 'calculate_barcode_quality_tenx', 'calculate_basecalling_qc', 'calculate_basecalling_qc_tenx', 'calculate_contamination', 'calculate_contamination_tenx', 'calculate_duplicate_rate', 'calculate_duplication', 'calculate_duplication_tenx', 'calculate_error_profile', 'calculate_error_profile_tenx', 'calculate_gc_bias', 'calculate_gc_bias_tenx', 'calculate_index_hopping', 'calculate_index_hopping_tenx', 'calculate_insert_size', 'calculate_insert_size_distribution', 'calculate_insert_size_tenx', 'calculate_phasing_qc', 'calculate_phasing_qc_tenx', 'calculate_qscore_summary', 'calculate_quality', 'calculate_quality_tenx', 'calculate_read_length', 'calculate_read_length_tenx', 'calculate_yield', 'calculate_yield_tenx', 'collapse_umi_reads', 'correct_long_reads', 'correct_umi_errors', 'demultiplex_barcoded_reads', 'detect_cross_sample_contamination', 'detect_platform', 'detect_platform_tenx', 'estimate_library_complexity', 'estimate_read_error_model', 'extract_umi_sequences', 'generate_qc_report', 'generate_qc_report_tenx', 'plot_per_base_quality', 'plot_qc_summary', 'plot_qc_summary_tenx', 'run_cutadapt', 'run_fastp', 'run_fastqc', 'run_multiqc'])


from dnakit._verbosity import register_public as _register_public_verbose

# Optional-backend exports (see docs/api/features) 
from .backends import (
    calculate_expected_errors,
    cluster_umi_counts,
    extract_read_barcodes,
    find_adapter_match,
    find_quality_trim_bounds,
    map_barcodes_to_whitelist,
    read_slow5_signal,
    summarize_slow5_signals,
    trim_fastq_adapters,
    write_slow5_signals,
)

__all__ = [*__all__, "calculate_expected_errors", "cluster_umi_counts", "extract_read_barcodes", "find_adapter_match", "find_quality_trim_bounds", "map_barcodes_to_whitelist", "read_slow5_signal", "summarize_slow5_signals", "trim_fastq_adapters", "write_slow5_signals"]

# Full-pass exports (see docs/api/features)
from .cutadapt_functions import (
    adapter_kmer_positions,
    adapter_kmer_table,
    edit_neighbours,
    expand_adapter_braces,
    hamming_neighbours,
    nextseq_trim_bounds,
    parse_adapter_spec,
    poly_a_tail_start,
    read_adapter_fasta,
    sequence_kmer_chunks,
)

from .umi_tools_functions import (
    average_umi_distance,
    barcode_indel_matches,
    call_cell_whitelist,
    cigar_splice_offset,
    count_umis_per_gene,
    detect_bam_tags,
    extract_barcode_positions,
    find_umi_neighbours,
    mask_low_quality_umi,
    parse_umi_read_name,
    read_barcode_whitelist,
    tag_read_name_with_umi,
    umi_connected_component,
)

__all__ = [*__all__, "adapter_kmer_positions", "adapter_kmer_table", "edit_neighbours", "expand_adapter_braces", "hamming_neighbours", "nextseq_trim_bounds", "parse_adapter_spec", "poly_a_tail_start", "read_adapter_fasta", "sequence_kmer_chunks", "average_umi_distance", "barcode_indel_matches", "call_cell_whitelist", "cigar_splice_offset", "count_umis_per_gene", "detect_bam_tags", "extract_barcode_positions", "find_umi_neighbours", "mask_low_quality_umi", "parse_umi_read_name", "read_barcode_whitelist", "tag_read_name_with_umi", "umi_connected_component"]

# Review-pass exports (see docs/api/features)
from .cutadapt_review_functions import (
    adapter_untrimmed_interval,
)

from .umi_tools_review_functions import (
    barcode_from_bam_tag,
    metacontig_contig_map,
    parse_prefixed_umi_name,
    resolve_above_knee_barcodes,
    soft_clip_adjusted_position,
    tag_reads_by_metacontig,
)

__all__ = [*__all__, "adapter_untrimmed_interval", "barcode_from_bam_tag", "metacontig_contig_map", "parse_prefixed_umi_name", "resolve_above_knee_barcodes", "soft_clip_adjusted_position", "tag_reads_by_metacontig"]

_register_public_verbose(globals())
