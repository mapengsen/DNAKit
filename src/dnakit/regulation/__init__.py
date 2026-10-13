"""Local candidate adapters for the regulation capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'annotate', input, options=options)

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'calculate', input, options=options)

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'classify', input, options=options)

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'compare', input, options=options)

def enrich(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'enrich', input, options=options)

def link(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'link', input, options=options)

def perform_in_silico_mutagenesis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'perform_in_silico_mutagenesis', input, options=options)

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'predict', input, options=options)

def prioritize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'prioritize', input, options=options)

def score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'score', input, options=options)

def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('regulation', 'validate', input, options=options)

__all__ = ['annotate', 'calculate', 'classify', 'compare', 'enrich', 'link', 'perform_in_silico_mutagenesis', 'predict', 'prioritize', 'score', 'validate']

__candidate_adapter_operations__ = set(['annotate', 'calculate', 'classify', 'compare', 'enrich', 'link', 'perform_in_silico_mutagenesis', 'predict', 'prioritize', 'score', 'validate'])


from dnakit._verbosity import register_public as _register_public_verbose
# Optional-backend exports (see docs/api/features) 
from .backends import (
    count_peak_fragments,
    count_tile_matrix,
    deep_lift_shap_attribution,
    marginalize_motif_effect,
    predict_sequences,
    read_fragment_qc,
    read_motif_file,
    scan_motif_hits,
    score_motif_kmers,
)

__all__ = [*__all__, "count_peak_fragments", "count_tile_matrix", "deep_lift_shap_attribution", "marginalize_motif_effect", "predict_sequences", "read_fragment_qc", "read_motif_file", "scan_motif_hits", "score_motif_kmers"]

# Full-pass exports (see docs/api/features)
from .gimmemotifs_functions import (
    cluster_motif_set,
    compare_motif_positions,
    equal_size_intervals,
    motif_enrichment_stats,
    motif_from_alignment,
    motif_from_consensus,
    motif_score_correlation,
    promoter_regions,
    random_genomic_regions,
    read_peak_summits,
    scan_best_motif_hits,
    summit_intervals_from_narrowpeak,
)

from .snapatac2_functions import (
    chromosome_fragment_counts,
    fragment_in_peaks_fraction,
    fragment_size_distribution,
    gene_activity_counts,
    merge_peak_sets,
    region_motif_enrichment,
    regulatory_regions_to_genes,
    tss_enrichment,
)

from .tangermeme_functions import (
    ablate_region,
    call_seqlets,
    deletion_effect,
    design_greedy_substitution,
    extract_loci_sequences,
    insertion_effect,
    kmer_counts,
    matched_background_loci,
    motif_pair_spacing,
    place_motif,
    pwm_information_content,
    saturation_mutagenesis,
    substitution_effect,
)

__all__ = [*__all__, "cluster_motif_set", "compare_motif_positions", "equal_size_intervals", "motif_enrichment_stats", "motif_from_alignment", "motif_from_consensus", "motif_score_correlation", "promoter_regions", "random_genomic_regions", "read_peak_summits", "scan_best_motif_hits", "summit_intervals_from_narrowpeak", "chromosome_fragment_counts", "fragment_in_peaks_fraction", "fragment_size_distribution", "gene_activity_counts", "merge_peak_sets", "region_motif_enrichment", "regulatory_regions_to_genes", "tss_enrichment", "ablate_region", "call_seqlets", "deletion_effect", "design_greedy_substitution", "extract_loci_sequences", "insertion_effect", "kmer_counts", "matched_background_loci", "motif_pair_spacing", "place_motif", "pwm_information_content", "saturation_mutagenesis", "substitution_effect"]

# Review-pass exports (see docs/api/features)
from .gimmemotifs_review_functions import (
    combine_summit_files,
    equal_length_motif_pair,
    fasta_record_count,
    motif_activity_scores,
    motif_position_ks,
    motif_score_metric,
    peak_coverage_table,
    scan_motifs_to_file,
    scan_region_scores,
)

from .snapatac2_review_functions import (
    add_tf_binding_edges,
    call_peaks_macs3,
    export_coverage_tracks,
    export_fragment_files,
    fragments_from_bam,
    import_base_pair_values,
    import_contact_pairs,
    tenx_atac_metrics,
)

from .tangermeme_review_functions import (
    ablate_annotations,
    annotate_seqlets,
    count_seqlet_annotations,
    delete_region,
    example_to_sequence_coords,
    gapped_kmer_features,
    greedy_marginalize_design,
    local_dinucleotide_shuffle,
    marginalize_annotations,
    one_hot_to_fasta_file,
    pairwise_annotation_counts,
    pisa_attribution,
    place_motif_set,
    pwm_consensus_sequence,
    random_one_hot_sequences,
    randomize_region,
    screen_sequence_design,
    shuffle_region,
    tfmodisco_seqlets,
)

__all__ = [*__all__, "combine_summit_files", "equal_length_motif_pair", "fasta_record_count", "motif_activity_scores", "motif_position_ks", "motif_score_metric", "peak_coverage_table", "scan_motifs_to_file", "scan_region_scores", "add_tf_binding_edges", "call_peaks_macs3", "export_coverage_tracks", "export_fragment_files", "fragments_from_bam", "import_base_pair_values", "import_contact_pairs", "tenx_atac_metrics", "ablate_annotations", "annotate_seqlets", "count_seqlet_annotations", "delete_region", "example_to_sequence_coords", "gapped_kmer_features", "greedy_marginalize_design", "local_dinucleotide_shuffle", "marginalize_annotations", "one_hot_to_fasta_file", "pairwise_annotation_counts", "pisa_attribution", "place_motif_set", "pwm_consensus_sequence", "random_one_hot_sequences", "randomize_region", "screen_sequence_design", "shuffle_region", "tfmodisco_seqlets"]

# Family-pass exports (see docs/api/features)
from .gimmemotifs_family_functions import (
    calculate_partial_motif_roc_auc,
    count_motif_matches_at_fpr,
    normalize_peak_widths,
)

from .tangermeme_family_functions import (
    design_beam_substitution,
    insert_motif,
)

__all__ = [*__all__, "calculate_partial_motif_roc_auc", "count_motif_matches_at_fpr", "normalize_peak_widths", "design_beam_substitution", "insert_motif"]

_register_public_verbose(globals())
