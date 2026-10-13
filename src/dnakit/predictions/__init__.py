"""Custom training and official pretrained property-prediction APIs."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.predictions.api import (
    predict_enformer_benchmark,
    predict_pair_properties,
    predict_properties,
    predict_sequence_properties,
    predict_variant_effects,
)
from dnakit.predictions.backends import PredictionBackend
from dnakit.predictions.checkpoints import (
    PredictionCheckpointInfo,
    default_prediction_checkpoint_root,
    enformer_benchmark_checkpoint_path,
    ensure_prediction_checkpoint,
)
from dnakit.predictions.custom_training import (
    CustomPropertyPredictionResult,
    CustomPropertyTrainingConfig,
    CustomPropertyTrainingResult,
    CustomTrainingModel,
    DataInspection,
    DataMapping,
    DataReader,
    LLMSchemaConfig,
    MachineLearningModel,
    OpenAICompatibleSchemaClient,
    PropertyData,
    SchemaLLMClient,
    SchemaLLMProvider,
    SplitAssignment,
    TargetSchema,
    TrainableSequenceEncoder,
    available_custom_training_models,
    available_machine_learning_models,
    available_schema_llm_providers,
    compute_metrics,
    detect_data_format,
    get_custom_training_model,
    get_machine_learning_model,
    get_schema_llm_provider,
    inspect_property_data,
    load_property_data,
    predict_custom_properties,
    split_property_data,
    train_custom_property_model,
)
from dnakit.predictions.enformer_benchmarks import (
    ENFORMER_BENCHMARK_CHECKPOINTS_URL,
    ENFORMER_BENCHMARK_TASKS,
    EnformerBenchmarkFamily,
    EnformerBenchmarkTask,
    available_enformer_benchmark_tasks,
    get_enformer_benchmark_task,
)
from dnakit.predictions.models import (
    BiologicalSequence,
    BiologicalSequencePair,
    DirectPredictionModel,
    DirectPredictionTask,
    PropertyPredictionConfig,
    VariantContext,
    available_prediction_models,
    available_prediction_tasks,
    get_prediction_model,
    get_prediction_task,
)
from dnakit.predictions.results import (
    PredictionOutput,
    PredictionRecord,
    PropertyPredictionResult,
)

__all__ = [
    "ENFORMER_BENCHMARK_CHECKPOINTS_URL",
    "ENFORMER_BENCHMARK_TASKS",
    "BiologicalSequence",
    "BiologicalSequencePair",
    "CustomPropertyPredictionResult",
    "CustomPropertyTrainingConfig",
    "CustomPropertyTrainingResult",
    "CustomTrainingModel",
    "DataInspection",
    "DataMapping",
    "DataReader",
    "DirectPredictionModel",
    "DirectPredictionTask",
    "EnformerBenchmarkFamily",
    "EnformerBenchmarkTask",
    "LLMSchemaConfig",
    "MachineLearningModel",
    "OpenAICompatibleSchemaClient",
    "PredictionBackend",
    "PredictionCheckpointInfo",
    "PredictionOutput",
    "PredictionRecord",
    "PropertyData",
    "PropertyPredictionConfig",
    "PropertyPredictionResult",
    "SchemaLLMClient",
    "SchemaLLMProvider",
    "SplitAssignment",
    "TargetSchema",
    "TrainableSequenceEncoder",
    "VariantContext",
    "available_custom_training_models",
    "available_enformer_benchmark_tasks",
    "available_machine_learning_models",
    "available_prediction_models",
    "available_prediction_tasks",
    "available_schema_llm_providers",
    "compute_metrics",
    "default_prediction_checkpoint_root",
    "detect_data_format",
    "enformer_benchmark_checkpoint_path",
    "ensure_prediction_checkpoint",
    "get_custom_training_model",
    "get_enformer_benchmark_task",
    "get_machine_learning_model",
    "get_prediction_model",
    "get_prediction_task",
    "get_schema_llm_provider",
    "inspect_property_data",
    "load_property_data",
    "predict_custom_properties",
    "predict_enformer_benchmark",
    "predict_pair_properties",
    "predict_properties",
    "predict_sequence_properties",
    "predict_variant_effects",
    "split_property_data",
    "train_custom_property_model",
]

# Optional-backend exports (see docs/api/features) 
from .dna_model_backends import (
    check_genomic_benchmark_cache,
    describe_genomic_benchmark,
    describe_kipoi_model,
    encode_one_hot_enformer,
    encode_one_hot_selene,
    extract_gencode_transcripts,
    fetch_genome_sequences,
    generate_evo2_sequences,
    list_genomic_benchmark_datasets,
    list_kipoi_models,
    merge_alphagenome_intervals,
    mutate_dna_sequence,
    parse_alphagenome_interval,
    parse_alphagenome_variant,
    read_contig_lengths,
    read_enformer_bed_windows,
    read_gencode_gene_intervals,
    resize_dna_sequences,
    reverse_complement_with_grelu,
    shuffle_dinucleotides,
)

__all__ = [*__all__, "check_genomic_benchmark_cache", "describe_genomic_benchmark", "describe_kipoi_model", "encode_one_hot_enformer", "encode_one_hot_selene", "extract_gencode_transcripts", "fetch_genome_sequences", "generate_evo2_sequences", "list_genomic_benchmark_datasets", "list_kipoi_models", "merge_alphagenome_intervals", "mutate_dna_sequence", "parse_alphagenome_interval", "parse_alphagenome_variant", "read_contig_lengths", "read_enformer_bed_windows", "read_gencode_gene_intervals", "resize_dna_sequences", "reverse_complement_with_grelu", "shuffle_dinucleotides"]

# Full-pass exports (see docs/api/features)
from .alphagenome_functions import (
    build_alphagenome_ism_matrix,
    combine_alphagenome_intervals,
    enumerate_alphagenome_ism_variants,
    extract_alphagenome_tss,
    filter_alphagenome_transcripts,
    get_alphagenome_fold_intervals,
    normalize_alphagenome_variant,
    upgrade_alphagenome_annotation_ids,
    validate_alphagenome_sequence_length,
)

from .enformer_pytorch_functions import (
    encode_enformer_indices,
    indices_to_one_hot_enformer,
    reverse_complement_enformer,
)

from .genomic_benchmarks_functions import (
    locate_sequences_in_reference,
)

from .grelu_functions import (
    check_dna_encoding,
    clip_intervals_to_bounds,
    compare_motif_pwm_alleles,
    convert_dna_encoding,
    dna_sequence_differences,
    extend_intervals_from_centers,
    filter_intervals_blacklist,
    filter_intervals_by_chromosome,
    filter_intervals_overlapping_reference,
    filter_variants_by_length,
    gc_content_distribution,
    merge_intervals_by_group,
    pwm_to_consensus,
    read_fasta_sequences,
    sample_intervals,
    scan_motifs_with_pwms,
    split_intervals_by_chromosome,
    summarize_sequence_lengths,
    trim_motif_pwm,
    variants_to_sequence_intervals,
)

from .selene_sdk_functions import (
    apply_one_hot_mutations,
    decode_one_hot_sequence,
    interval_length_probabilities,
    list_in_silico_mutants,
    load_variant_diff_scores,
    order_variant_labels,
    read_selene_vcf,
    rescale_attribution_scores,
    reverse_complement_allele,
    reverse_complement_one_hot,
    sort_standard_chromosomes,
    unpack_bit_packed_sequence,
)

__all__ = [*__all__, "build_alphagenome_ism_matrix", "combine_alphagenome_intervals", "enumerate_alphagenome_ism_variants", "extract_alphagenome_tss", "filter_alphagenome_transcripts", "get_alphagenome_fold_intervals", "normalize_alphagenome_variant", "upgrade_alphagenome_annotation_ids", "validate_alphagenome_sequence_length", "encode_enformer_indices", "indices_to_one_hot_enformer", "reverse_complement_enformer", "locate_sequences_in_reference", "check_dna_encoding", "clip_intervals_to_bounds", "compare_motif_pwm_alleles", "convert_dna_encoding", "dna_sequence_differences", "extend_intervals_from_centers", "filter_intervals_blacklist", "filter_intervals_by_chromosome", "filter_intervals_overlapping_reference", "filter_variants_by_length", "gc_content_distribution", "merge_intervals_by_group", "pwm_to_consensus", "read_fasta_sequences", "sample_intervals", "scan_motifs_with_pwms", "split_intervals_by_chromosome", "summarize_sequence_lengths", "trim_motif_pwm", "variants_to_sequence_intervals", "apply_one_hot_mutations", "decode_one_hot_sequence", "interval_length_probabilities", "list_in_silico_mutants", "load_variant_diff_scores", "order_variant_labels", "read_selene_vcf", "rescale_attribution_scores", "reverse_complement_allele", "reverse_complement_one_hot", "sort_standard_chromosomes", "unpack_bit_packed_sequence"]

# Review-pass exports (see docs/api/features)
from .alphagenome_review_functions import (
    convert_alphagenome_variant_scores,
    list_alphagenome_recommended_scorers,
    merge_alphagenome_stranded_gene_tracks,
    merge_alphagenome_track_metadata,
    tidy_alphagenome_anndata,
    tidy_alphagenome_scores,
    unmerge_alphagenome_stranded_gene_tracks,
)

from .grelu_review_functions import (
    add_negative_intervals,
    check_genomic_interval_table,
    design_sequence_ledidi,
    design_sequences_evolve,
    expand_chromosome_names,
    filter_interval_coverage,
    get_attention_scores,
    get_sequence_attributions,
    grelu_marginalize_variants,
    grelu_predict_variant_effects,
    in_silico_mutagenesis_scores,
    marginalize_pattern_insertions,
    marginalize_pattern_spacing,
    read_modisco_motifs,
    read_narrowpeak_intervals,
    read_tomtom_matches,
    run_modisco_motifs,
    score_motif_pwm_products,
    score_motif_sites,
    shuffle_tile_effects,
)

__all__ = [*__all__, "convert_alphagenome_variant_scores", "list_alphagenome_recommended_scorers", "merge_alphagenome_stranded_gene_tracks", "merge_alphagenome_track_metadata", "tidy_alphagenome_anndata", "tidy_alphagenome_scores", "unmerge_alphagenome_stranded_gene_tracks", "add_negative_intervals", "check_genomic_interval_table", "design_sequence_ledidi", "design_sequences_evolve", "expand_chromosome_names", "filter_interval_coverage", "get_attention_scores", "get_sequence_attributions", "grelu_marginalize_variants", "grelu_predict_variant_effects", "in_silico_mutagenesis_scores", "marginalize_pattern_insertions", "marginalize_pattern_spacing", "read_modisco_motifs", "read_narrowpeak_intervals", "read_tomtom_matches", "run_modisco_motifs", "score_motif_pwm_products", "score_motif_sites", "shuffle_tile_effects"]

_register_public_verbose(globals())
