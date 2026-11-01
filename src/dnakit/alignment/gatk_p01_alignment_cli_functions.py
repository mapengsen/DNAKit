"""GATK alignment and read-file command-line facades (pack P01): BAM/SAM/CRAM editing, duplicate marking, BQSR, downsampling, PathSeq.

Each public function runs one GATK 4.7.0.0 tool (or one family of tools selected by a mode argument) through
``dnakit.backends.external.execute_bounded_command`` with a fixed option whitelist and no shell. The user installs GATK
(conda ``gatk4`` or the release zip) with Java 17 and either puts the ``gatk`` launcher on PATH, sets ``DNAKIT_GATK`` or
passes ``gatk_path``. Outputs are written to explicit new paths that must not exist; nothing is overwritten. Results are
JSON-serialisable dictionaries with status, backend, tool, version, command, outputs, a bounded log tail and provenance.
35 public functions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gatk_p01_core as _core
from dnakit.bioinformatics._gatk_p01_core import Param as _Param
from dnakit.bioinformatics._gatk_p01_core import Spec as _Spec

_SPECS: dict[str, _Spec] = {
    'gatk_pathseq': _Spec(
        'gatk_pathseq',
        {'filter': 'PathSeqFilterSpark', 'align': 'PathSeqBwaSpark', 'score': 'PathSeqScoreSpark', 'pipeline': 'PathSeqPipelineSpark'},
        (
            _Param('input_path', 'in', 'input', v=('filter', 'pipeline'), req=('filter', 'pipeline')),
            _Param('paired_input', 'in', 'paired-input', v=('align', 'score')),
            _Param('unpaired_input', 'in', 'unpaired-input', v=('align', 'score')),
            _Param('microbe_bwa_image', 'in', 'microbe-bwa-image', v=('align', 'pipeline'), req=('align', 'pipeline')),
            _Param('microbe_dict', 'in', 'microbe-dict', v=('align', 'pipeline'), req=('align', 'pipeline')),
            _Param('taxonomy_file', 'in', 'taxonomy-file', v=('score', 'pipeline'), req=('score', 'pipeline')),
            _Param('paired_output', 'out', 'paired-output', v=('filter', 'align')),
            _Param('unpaired_output', 'out', 'unpaired-output', v=('filter', 'align')),
            _Param('scores_output', 'out', 'scores-output', v=('score', 'pipeline'), req=('score', 'pipeline')),
        ),
        selector='step',
        any_of=((('align', 'score'), ('paired_input', 'unpaired_input')), (('filter', 'align'), ('paired_output', 'unpaired_output'))),
    ),
    'gatk_detect_cram_issue_8768': _Spec(
        'gatk_detect_cram_issue_8768',
        {'default': 'CRAMIssue8768Detector'},
        (
            _Param('cram_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_create_bam_splitting_index': _Spec(
        'gatk_create_bam_splitting_index',
        {'default': 'CreateHadoopBamSplittingIndex'},
        (
            _Param('bam_path', 'in', 'input', req=True),
            _Param('create_bai', 'bool', 'create-bai'),
        ),
        derived=(('sbi', 'bam_path', ('.sbi',), True),),
    ),
    'gatk_print_bgzf_block_information': _Spec(
        'gatk_print_bgzf_block_information',
        {'default': 'PrintBGZFBlockInformation'},
        (
            _Param('bgzf_path', 'in', 'bgzf-file', req=True),
        ),
        hook='bgzf_blocks',
        big_output=True,
    ),
    'gatk_print_file_diagnostics': _Spec(
        'gatk_print_file_diagnostics',
        {'default': 'PrintFileDiagnostics'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_read_anonymizer': _Spec(
        'gatk_read_anonymizer',
        {'default': 'ReadAnonymizer'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('ref_base_quality', 'int', 'ref-base-quality', minimum=0),
            _Param('use_simple_cigar', 'bool', 'use-simple-cigar'),
        ),
    ),
    'gatk_add_comments_to_bam': _Spec(
        'gatk_add_comments_to_bam',
        {'default': 'AddCommentsToBam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('comments', 'strs', 'COMMENT', req=True),
        ),
    ),
    'gatk_add_original_alignment_tags': _Spec(
        'gatk_add_original_alignment_tags',
        {'oa_tag': 'AddOATag', 'original_alignment': 'AddOriginalAlignmentTags'},
        (
            _Param('input_path', 'in', {'oa_tag': 'INPUT', 'original_alignment': 'input'}, req=True),
            _Param('output_path', 'out', {'oa_tag': 'OUTPUT', 'original_alignment': 'output'}, req=True),
            _Param('interval_list', 'in', 'INTERVAL_LIST', v=('oa_tag',)),
        ),
        selector='mode',
    ),
    'gatk_base_recalibrator': _Spec(
        'gatk_base_recalibrator',
        {'default': 'BaseRecalibrator'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('known_sites', 'ins', 'known-sites', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_apply_bqsr': _Spec(
        'gatk_apply_bqsr',
        {'default': 'ApplyBQSR'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('recalibration_table', 'in', 'bqsr-recal-file', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('reference_path', 'in', 'reference'),
            _Param('static_quantized_quals', 'strs', 'static-quantized-quals'),
            _Param('use_original_qualities', 'bool', 'use-original-qualities'),
            _Param('emit_original_quals', 'bool', 'emit-original-quals'),
        ),
    ),
    'gatk_bqsr_pipeline_spark': _Spec(
        'gatk_bqsr_pipeline_spark',
        {'default': 'BQSRPipelineSpark'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('known_sites', 'ins', 'known-sites', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_bam_to_bfq': _Spec(
        'gatk_bam_to_bfq',
        {'default': 'BamToBfq'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('analysis_dir', 'outd', 'ANALYSIS_DIR', req=True),
            _Param('output_file_prefix', 'str', 'OUTPUT_FILE_PREFIX'),
            _Param('flowcell_barcode', 'str', 'FLOWCELL_BARCODE'),
            _Param('paired_run', 'bool', 'PAIRED_RUN', req=True),
            _Param('lane', 'int', 'LANE', minimum=1),
        ),
        any_of=((None, ('output_file_prefix', 'flowcell_barcode')),),
    ),
    'gatk_bwa_spark': _Spec(
        'gatk_bwa_spark',
        {'align': 'BwaSpark', 'align_mark_duplicates': 'BwaAndMarkDuplicatesPipelineSpark'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('bwa_mem_index_image', 'in', 'bwa-mem-index-image'),
            _Param('single_end_alignment', 'bool', 'single-end-alignment'),
        ),
        selector='pipeline',
    ),
    'gatk_clean_sam': _Spec(
        'gatk_clean_sam',
        {'default': 'CleanSam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
        ),
    ),
    'gatk_clip_reads': _Spec(
        'gatk_clip_reads',
        {'default': 'ClipReads'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('clip_sequences', 'strs', 'clip-sequence'),
            _Param('cycles_to_trim', 'str', 'cycles-to-trim'),
            _Param('q_trimming_threshold', 'int', 'q-trimming-threshold', minimum=0),
            _Param('clip_representation', 'choice', 'clip-representation', choices=('WRITE_NS', 'WRITE_Q0S', 'WRITE_NS_Q0S', 'SOFTCLIP_BASES', 'HARDCLIP_BASES', 'REVERT_SOFTCLIPPED_BASES')),
            _Param('statistics_path', 'out', 'output-statistics'),
        ),
    ),
    'gatk_convert_headerless_bam_shard': _Spec(
        'gatk_convert_headerless_bam_shard',
        {'default': 'ConvertHeaderlessHadoopBamShardToBam'},
        (
            _Param('shard_path', 'in', 'bam-shard', req=True),
            _Param('header_bam_path', 'in', 'bam-with-header', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_downsample_by_duplicate_set': _Spec(
        'gatk_downsample_by_duplicate_set',
        {'default': 'DownsampleByDuplicateSet'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('fraction_to_keep', 'float', 'fraction-to-keep', minimum=0),
            _Param('min_reads', 'int', 'min-reads', minimum=0),
            _Param('min_per_strand_reads', 'int', 'min-per-strand-reads', minimum=0),
        ),
    ),
    'gatk_downsample_sam': _Spec(
        'gatk_downsample_sam',
        {'default': 'DownsampleSam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('probability', 'float', 'PROBABILITY', minimum=0),
            _Param('strategy', 'choice', 'STRATEGY', choices=('HighAccuracy', 'ConstantMemory', 'Chained')),
            _Param('accuracy', 'float', 'ACCURACY', minimum=0),
            _Param('random_seed', 'int', 'RANDOM_SEED'),
            _Param('metrics_path', 'out', 'METRICS_FILE'),
        ),
    ),
    'gatk_position_based_downsample_sam': _Spec(
        'gatk_position_based_downsample_sam',
        {'default': 'PositionBasedDownsampleSam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('fraction', 'float', 'FRACTION', req=True, minimum=0),
            _Param('stop_after', 'int', 'STOP_AFTER', minimum=0),
        ),
    ),
    'gatk_extract_original_alignment_records_by_name': _Spec(
        'gatk_extract_original_alignment_records_by_name',
        {'default': 'ExtractOriginalAlignmentRecordsByNameSpark'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('read_name_file', 'in', 'read-name-file', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('invert_match', 'bool', 'invert-match'),
        ),
    ),
    'gatk_filter_sam_reads': _Spec(
        'gatk_filter_sam_reads',
        {'default': 'FilterSamReads'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('filter', 'choice', 'FILTER', req=True, choices=('includeAligned', 'excludeAligned', 'includeReadList', 'excludeReadList', 'includePairedIntervals', 'includeTagValues', 'excludeTagValues')),
            _Param('read_list_file', 'in', 'READ_LIST_FILE'),
            _Param('interval_list', 'in', 'INTERVAL_LIST'),
            _Param('tag', 'str', 'TAG'),
            _Param('tag_values', 'strs', 'TAG_VALUE'),
        ),
    ),
    'gatk_fix_misencoded_base_quality_reads': _Spec(
        'gatk_fix_misencoded_base_quality_reads',
        {'default': 'FixMisencodedBaseQualityReads'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_left_align_indels': _Spec(
        'gatk_left_align_indels',
        {'default': 'LeftAlignIndels'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_mark_duplicates_with_mate_cigar': _Spec(
        'gatk_mark_duplicates_with_mate_cigar',
        {'standard': 'MarkDuplicatesWithMateCigar', 'simple': 'SimpleMarkDuplicatesWithMateCigar', 'umi_aware': 'UmiAwareMarkDuplicatesWithMateCigar'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('metrics_path', 'out', 'METRICS_FILE', req=True),
            _Param('umi_metrics_path', 'out', 'UMI_METRICS_FILE', v=('umi_aware',), req=('umi_aware',)),
            _Param('remove_duplicates', 'bool', 'REMOVE_DUPLICATES'),
            _Param('umi_tag_name', 'str', 'UMI_TAG_NAME', v=('umi_aware',)),
        ),
        selector='variant',
        tables={'metrics_path': 'metrics'},
    ),
    'gatk_merge_bam_alignment': _Spec(
        'gatk_merge_bam_alignment',
        {'default': 'MergeBamAlignment'},
        (
            _Param('unmapped_bam', 'in', 'UNMAPPED_BAM', req=True),
            _Param('reference_path', 'in', 'REFERENCE_SEQUENCE', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('aligned_bams', 'ins', 'ALIGNED_BAM'),
            _Param('paired_run', 'bool', 'PAIRED_RUN'),
            _Param('sort_order', 'choice', 'SORT_ORDER', choices=('coordinate', 'queryname', 'unsorted')),
            _Param('add_mate_cigar', 'bool', 'ADD_MATE_CIGAR'),
            _Param('clip_adapters', 'bool', 'CLIP_ADAPTERS'),
        ),
    ),
    'gatk_print_distant_mates': _Spec(
        'gatk_print_distant_mates',
        {'default': 'PrintDistantMates'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_reorder_sam': _Spec(
        'gatk_reorder_sam',
        {'default': 'ReorderSam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('sequence_dictionary', 'in', 'SEQUENCE_DICTIONARY', req=True),
            _Param('allow_incomplete_dict_concordance', 'bool', 'ALLOW_INCOMPLETE_DICT_CONCORDANCE'),
            _Param('allow_contig_length_discordance', 'bool', 'ALLOW_CONTIG_LENGTH_DISCORDANCE'),
        ),
    ),
    'gatk_revert_base_quality_scores': _Spec(
        'gatk_revert_base_quality_scores',
        {'original_qualities': 'RevertBaseQualityScores', 'with_mate_cigar': 'RevertOriginalBaseQualitiesAndAddMateCigar'},
        (
            _Param('input_path', 'in', {'original_qualities': 'input', 'with_mate_cigar': 'INPUT'}, req=True),
            _Param('output_path', 'out', {'original_qualities': 'output', 'with_mate_cigar': 'OUTPUT'}, req=True),
            _Param('restore_original_qualities', 'bool', 'RESTORE_ORIGINAL_QUALITIES', v=('with_mate_cigar',)),
            _Param('sort_order', 'choice', 'SORT_ORDER', v=('with_mate_cigar',), choices=('coordinate', 'queryname', 'unsorted')),
        ),
        selector='mode',
    ),
    'gatk_revert_sam': _Spec(
        'gatk_revert_sam',
        {'default': 'RevertSam'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('restore_original_qualities', 'bool', 'RESTORE_ORIGINAL_QUALITIES'),
            _Param('remove_duplicate_information', 'bool', 'REMOVE_DUPLICATE_INFORMATION'),
            _Param('remove_alignment_information', 'bool', 'REMOVE_ALIGNMENT_INFORMATION'),
            _Param('sanitize', 'bool', 'SANITIZE'),
            _Param('sample_alias', 'str', 'SAMPLE_ALIAS'),
            _Param('library_name', 'str', 'LIBRARY_NAME'),
            _Param('sort_order', 'choice', 'SORT_ORDER', choices=('coordinate', 'queryname', 'unsorted')),
        ),
    ),
    'gatk_convert_sam_format': _Spec(
        'gatk_convert_sam_format',
        {'default': 'SamFormatConverter'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
        ),
    ),
    'gatk_sam_to_fastq_with_tags': _Spec(
        'gatk_sam_to_fastq_with_tags',
        {'default': 'SamToFastqWithTags'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('fastq_path', 'out', 'FASTQ', req=True),
            _Param('sequence_tag_groups', 'strs', 'SEQUENCE_TAG_GROUP', req=True),
            _Param('quality_tag_groups', 'strs', 'QUALITY_TAG_GROUP'),
            _Param('second_end_fastq', 'out', 'SECOND_END_FASTQ'),
            _Param('unpaired_fastq', 'out', 'UNPAIRED_FASTQ'),
            _Param('interleave', 'bool', 'INTERLEAVE'),
        ),
    ),
    'gatk_split_n_cigar_reads': _Spec(
        'gatk_split_n_cigar_reads',
        {'default': 'SplitNCigarReads'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('max_bases_in_overhang', 'int', 'max-bases-in-overhang', minimum=0),
            _Param('max_mismatches_in_overhang', 'int', 'max-mismatches-in-overhang', minimum=0),
            _Param('skip_mapping_quality_transform', 'bool', 'skip-mapping-quality-transform'),
        ),
    ),
    'gatk_split_sam': _Spec(
        'gatk_split_sam',
        {'library': 'SplitSamByLibrary', 'number_of_reads': 'SplitSamByNumberOfReads'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('output_dir', 'outd', 'OUTPUT', req=True),
            _Param('split_to_n_files', 'int', 'SPLIT_TO_N_FILES', v=('number_of_reads',), minimum=1),
            _Param('split_to_n_reads', 'int', 'SPLIT_TO_N_READS', v=('number_of_reads',), minimum=1),
            _Param('total_reads_in_input', 'int', 'TOTAL_READS_IN_INPUT', v=('number_of_reads',), minimum=1),
            _Param('out_prefix', 'str', 'OUT_PREFIX', v=('number_of_reads',)),
        ),
        selector='by',
    ),
    'gatk_transfer_read_tags': _Spec(
        'gatk_transfer_read_tags',
        {'default': 'TransferReadTags'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('unmapped_sam', 'in', 'unmapped-sam', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('read_tags', 'strs', 'read-tags', req=True),
        ),
    ),
    'gatk_unmark_duplicates': _Spec(
        'gatk_unmark_duplicates',
        {'default': 'UnmarkDuplicates'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
}


def gatk_pathseq(
    *,
    step: str = 'filter',
    input_path: str | None = None,
    paired_input: str | None = None,
    unpaired_input: str | None = None,
    microbe_bwa_image: str | None = None,
    microbe_dict: str | None = None,
    taxonomy_file: str | None = None,
    paired_output: str | None = None,
    unpaired_output: str | None = None,
    scores_output: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the PathSeq pathogen-discovery pipeline: host filtering, microbe alignment, abundance scoring, or all steps at once.

    Run `gatk PathSeqFilterSpark`, `PathSeqBwaSpark`, `PathSeqScoreSpark` or `PathSeqPipelineSpark` (Spark local mode) depending on `step`. Host and microbe resources (k-mer file, host BWA image, hash files) are passed as GATK options such as `filter-bwa-image` or `kmer-file`.

    ``step``: `filter` = host-read filtering (PathSeqFilterSpark), `align` = microbe alignment (PathSeqBwaSpark), `score` = abundance scoring (PathSeqScoreSpark), `pipeline` = all steps (PathSeqPipelineSpark).

    Parameters:
    - ``input_path`` (required for some modes): Query BAM/CRAM/SAM of reads (host-filtering input).
    - ``paired_input`` (optional): Paired reads BAM from the filter step.
    - ``unpaired_input`` (optional): Unpaired reads BAM from the filter step.
    - ``microbe_bwa_image`` (required for some modes): BWA-MEM index image of the microbe reference (see `gatk_bwa_mem_index_image`).
    - ``microbe_dict`` (required for some modes): Sequence dictionary of the microbe reference.
    - ``taxonomy_file`` (required for some modes): Taxonomy database built by PathSeqBuildReferenceTaxonomy.
    - ``paired_output`` (optional): Output BAM for paired reads.
    - ``unpaired_output`` (optional): Output BAM for unpaired reads.
    - ``scores_output`` (required for some modes): Pathogen abundance score table (tab-separated).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields, `mode` and the output paths that were requested.
    """
    return _core.execute(
        _SPECS['gatk_pathseq'],
        {'input_path': input_path, 'paired_input': paired_input, 'unpaired_input': unpaired_input, 'microbe_bwa_image': microbe_bwa_image, 'microbe_dict': microbe_dict, 'taxonomy_file': taxonomy_file, 'paired_output': paired_output, 'unpaired_output': unpaired_output, 'scores_output': scores_output},
        mode=step,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_detect_cram_issue_8768(
    cram_path: str,
    reference_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a CRAM file for base corruption caused by GATK issue 8768.

    Run `gatk CRAMIssue8768Detector`; the mismatch-rate threshold and TSV output are GATK options (`mismatch-rate-threshold`, `output-tsv`).

    Parameters:
    - ``cram_path`` (required): CRAM file to check.
    - ``reference_path`` (required): Reference FASTA used for the CRAM (with .fai and .dict).
    - ``output_path`` (required): New report file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the report.
    """
    return _core.execute(
        _SPECS['gatk_detect_cram_issue_8768'],
        {'cram_path': cram_path, 'reference_path': reference_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_create_bam_splitting_index(
    bam_path: str,
    *,
    create_bai: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a Hadoop-BAM splitting index (.sbi) for a BAM file.

    Run `gatk CreateHadoopBamSplittingIndex`; the index `<bam>.sbi` is written next to the BAM and must not exist yet.

    Parameters:
    - ``bam_path`` (required): Coordinate-sorted BAM file.
    - ``create_bai`` (optional): Also write a .bai index.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.sbi` is the splitting index path.
    """
    return _core.execute(
        _SPECS['gatk_create_bam_splitting_index'],
        {'bam_path': bam_path, 'create_bai': create_bai},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_print_bgzf_block_information(
    bgzf_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report the compressed and uncompressed size and offset of every BGZF block of a file.

    Run `gatk PrintBGZFBlockInformation` and parse its `Block #n at file offset ...` listing (the first 10,000 blocks are returned; `block_count` is the total).

    Parameters:
    - ``bgzf_path`` (required): BGZF-compressed file (BAM, bgzipped VCF ...).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields plus `block_count` and `blocks` (block, offset, compressed_size, uncompressed_size).
    """
    return _core.execute(
        _SPECS['gatk_print_bgzf_block_information'],
        {'bgzf_path': bgzf_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_print_file_diagnostics(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write diagnostic information (record counts, header and format details) about a genomics file.

    Run `gatk PrintFileDiagnostics`; `count-limit` limits the number of records examined.

    Parameters:
    - ``input_path`` (required): Genomics file in a format the tool analyses (CRAM; GATK 4.7.0.0 rejects BAM).
    - ``output_path`` (required): New local diagnostics file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the diagnostics file.
    """
    return _core.execute(
        _SPECS['gatk_print_file_diagnostics'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_read_anonymizer(
    input_path: str,
    reference_path: str,
    output_path: str,
    *,
    ref_base_quality: int | None = None,
    use_simple_cigar: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the bases of reads with the reference bases they align to, so that no personal sequence remains.

    Run `gatk ReadAnonymizer`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New output alignment file.
    - ``ref_base_quality`` (optional): Base quality assigned to the replacing reference bases.
    - ``use_simple_cigar`` (optional): Rewrite the CIGAR as a simple full-length match.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the anonymised file.
    """
    return _core.execute(
        _SPECS['gatk_read_anonymizer'],
        {'input_path': input_path, 'reference_path': reference_path, 'output_path': output_path, 'ref_base_quality': ref_base_quality, 'use_simple_cigar': use_simple_cigar},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_add_comments_to_bam(
    input_path: str,
    output_path: str,
    comments: Sequence[str],
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add comment (@CO) lines to the header of a SAM/BAM/CRAM file.

    Run Picard-based `gatk AddCommentsToBam`; one `--COMMENT` per entry of `comments`.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New output file.
    - ``comments`` (required): Comment texts; each becomes a @CO header line.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_add_comments_to_bam'],
        {'input_path': input_path, 'output_path': output_path, 'comments': comments},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_add_original_alignment_tags(
    input_path: str,
    output_path: str,
    *,
    mode: str = 'oa_tag',
    interval_list: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Record the current alignment of every read in the OA tag (and the original mate contig).

    Run Picard `AddOATag` (`mode=oa_tag`) or GATK `AddOriginalAlignmentTags` (`mode=original_alignment`, which also adds the original mate contig tag).

    ``mode``: `oa_tag` = AddOATag, `original_alignment` = AddOriginalAlignmentTags.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New output file.
    - ``interval_list`` (optional): Only reads overlapping these intervals receive the OA tag.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_add_original_alignment_tags'],
        {'input_path': input_path, 'output_path': output_path, 'interval_list': interval_list},
        mode=mode,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_base_recalibrator(
    input_path: str,
    reference_path: str,
    known_sites: Sequence[str],
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the base quality score recalibration (BQSR) table.

    Run `gatk BaseRecalibrator`; covariates and quality settings are GATK options such as `low-quality-tail` or `maximum-cycle-value`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``known_sites`` (required): Known-variant VCF files (indexed) to mask.
    - ``output_path`` (required): New recalibration table.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the recalibration table.
    """
    return _core.execute(
        _SPECS['gatk_base_recalibrator'],
        {'input_path': input_path, 'reference_path': reference_path, 'known_sites': known_sites, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_apply_bqsr(
    input_path: str,
    recalibration_table: str,
    output_path: str,
    *,
    reference_path: str | None = None,
    static_quantized_quals: Sequence[str] | None = None,
    use_original_qualities: bool | None = None,
    emit_original_quals: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply a base quality score recalibration table to the reads.

    Run `gatk ApplyBQSR`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM.
    - ``recalibration_table`` (required): Table from `gatk_base_recalibrator`.
    - ``output_path`` (required): New recalibrated BAM/CRAM.
    - ``reference_path`` (optional): Reference FASTA (needed for CRAM).
    - ``static_quantized_quals`` (optional): Static quality levels to quantize to.
    - ``use_original_qualities`` (optional): Use the OQ tag as the input qualities.
    - ``emit_original_quals`` (optional): Keep the original qualities in the OQ tag.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the recalibrated file.
    """
    return _core.execute(
        _SPECS['gatk_apply_bqsr'],
        {'input_path': input_path, 'recalibration_table': recalibration_table, 'output_path': output_path, 'reference_path': reference_path, 'static_quantized_quals': static_quantized_quals, 'use_original_qualities': use_original_qualities, 'emit_original_quals': emit_original_quals},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_bqsr_pipeline_spark(
    input_path: str,
    reference_path: str,
    known_sites: Sequence[str],
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run both BQSR steps (BaseRecalibrator and ApplyBQSR) in one Spark-parallel pass.

    Run `gatk BQSRPipelineSpark` (beta); pass `spark-master` (default `local[*]`) through `options`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (a .2bit file also works for Spark).
    - ``known_sites`` (required): Known-variant VCF files (indexed).
    - ``output_path`` (required): New recalibrated BAM.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the recalibrated BAM.
    """
    return _core.execute(
        _SPECS['gatk_bqsr_pipeline_spark'],
        {'input_path': input_path, 'reference_path': reference_path, 'known_sites': known_sites, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_bam_to_bfq(
    input_path: str,
    analysis_dir: str,
    paired_run: bool,
    *,
    output_file_prefix: str | None = None,
    flowcell_barcode: str | None = None,
    lane: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a BAM file to BFQ (binary FASTQ) files for the Maq aligner.

    Run Picard `gatk BamToBfq`; the files written to `analysis_dir` are listed in `files`.

    Parameters:
    - ``input_path`` (required): Input BAM/SAM.
    - ``analysis_dir`` (required): New directory for the .bfq files.
    - ``output_file_prefix`` (optional): Prefix of the .bfq file names (use this or `flowcell_barcode`).
    - ``flowcell_barcode`` (optional): Flowcell barcode from which the file names are built (use this or `output_file_prefix`).
    - ``paired_run`` (required): Whether the run is paired-end.
    - ``lane`` (optional): Lane number to write.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `files` (file names inside `analysis_dir`).
    """
    return _core.execute(
        _SPECS['gatk_bam_to_bfq'],
        {'input_path': input_path, 'analysis_dir': analysis_dir, 'output_file_prefix': output_file_prefix, 'flowcell_barcode': flowcell_barcode, 'paired_run': paired_run, 'lane': lane},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_bwa_spark(
    input_path: str,
    reference_path: str,
    output_path: str,
    *,
    pipeline: str = 'align',
    bwa_mem_index_image: str | None = None,
    single_end_alignment: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align reads to a reference with BWA-MEM on Spark, optionally followed by duplicate marking.

    Run `gatk BwaSpark` (`pipeline=align`) or `BwaAndMarkDuplicatesPipelineSpark` (`pipeline=align_mark_duplicates`, name-sorted input); both are beta tools.

    ``pipeline``: `align` = BwaSpark, `align_mark_duplicates` = BwaAndMarkDuplicatesPipelineSpark.

    Parameters:
    - ``input_path`` (required): Unaligned (or name-sorted) BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New aligned BAM.
    - ``bwa_mem_index_image`` (optional): Prebuilt BWA-MEM index image (see `gatk_bwa_mem_index_image`).
    - ``single_end_alignment`` (optional): Align reads as single-end.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the aligned BAM.
    """
    return _core.execute(
        _SPECS['gatk_bwa_spark'],
        {'input_path': input_path, 'reference_path': reference_path, 'output_path': output_path, 'bwa_mem_index_image': bwa_mem_index_image, 'single_end_alignment': single_end_alignment},
        mode=pipeline,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_clean_sam(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clean a SAM/BAM/CRAM file: soft-clip alignments that hang off the end of the reference and set MAPQ to 0 for unmapped reads.

    Run Picard `gatk CleanSam`.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New cleaned file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the cleaned file.
    """
    return _core.execute(
        _SPECS['gatk_clean_sam'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_clip_reads(
    input_path: str,
    output_path: str,
    *,
    clip_sequences: Sequence[str] | None = None,
    cycles_to_trim: str | None = None,
    q_trimming_threshold: int | None = None,
    clip_representation: str | None = None,
    statistics_path: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clip reads by sequence, cycle or quality and report the clipping statistics.

    Run `gatk ClipReads`; `clip_representation` takes GATK's `WRITE_NS` (default), `WRITE_Q0S`, `WRITE_NS_Q0S`, `SOFTCLIP_BASES`, `HARDCLIP_BASES` or `REVERT_SOFTCLIPPED_BASES`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM.
    - ``output_path`` (required): New clipped file.
    - ``clip_sequences`` (optional): Sequences to clip wherever they occur in a read.
    - ``cycles_to_trim`` (optional): Cycle ranges to trim, for example `1-5,10-12`.
    - ``q_trimming_threshold`` (optional): Trim ends below this base quality.
    - ``clip_representation`` (optional): How clipped bases are written.
    - ``statistics_path`` (optional): New file for the clipping statistics.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the clipped file (and `statistics_path` when requested).
    """
    return _core.execute(
        _SPECS['gatk_clip_reads'],
        {'input_path': input_path, 'output_path': output_path, 'clip_sequences': clip_sequences, 'cycles_to_trim': cycles_to_trim, 'q_trimming_threshold': q_trimming_threshold, 'clip_representation': clip_representation, 'statistics_path': statistics_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_convert_headerless_bam_shard(
    shard_path: str,
    header_bam_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Turn a headerless BAM shard (Hadoop-BAM output) into a readable BAM by adding a header.

    Run `gatk ConvertHeaderlessHadoopBamShardToBam` (beta).

    Parameters:
    - ``shard_path`` (required): Headerless BAM shard written by a Spark tool.
    - ``header_bam_path`` (required): BAM file whose header is used.
    - ``output_path`` (required): New readable BAM.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the BAM.
    """
    return _core.execute(
        _SPECS['gatk_convert_headerless_bam_shard'],
        {'shard_path': shard_path, 'header_bam_path': header_bam_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_downsample_by_duplicate_set(
    input_path: str,
    output_path: str,
    *,
    fraction_to_keep: float | None = None,
    min_reads: int | None = None,
    min_per_strand_reads: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Discard a set fraction of the duplicate sets of a UMI-grouped BAM.

    Run `gatk DownsampleByDuplicateSet` (beta).

    Parameters:
    - ``input_path`` (required): UMI-grouped BAM carrying MI tags (`<id>/A`, `<id>/B`) and sorted by MI.
    - ``output_path`` (required): New downsampled BAM.
    - ``fraction_to_keep`` (optional): Fraction of duplicate sets to keep.
    - ``min_reads`` (optional): Minimum reads in a set to be kept.
    - ``min_per_strand_reads`` (optional): Minimum reads per strand in a set.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the downsampled BAM.
    """
    return _core.execute(
        _SPECS['gatk_downsample_by_duplicate_set'],
        {'input_path': input_path, 'output_path': output_path, 'fraction_to_keep': fraction_to_keep, 'min_reads': min_reads, 'min_per_strand_reads': min_per_strand_reads},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_downsample_sam(
    input_path: str,
    output_path: str,
    *,
    probability: float | None = None,
    strategy: str | None = None,
    accuracy: float | None = None,
    random_seed: int | None = None,
    metrics_path: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly keep a fraction of the reads of a SAM/BAM/CRAM file.

    Run Picard `gatk DownsampleSam`; strategy defaults to GATK's own default (`ConstantMemory`).

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New downsampled file.
    - ``probability`` (optional): Probability of keeping each read pair (0-1).
    - ``strategy`` (optional): Downsampling strategy.
    - ``accuracy`` (optional): Accuracy of the kept fraction (not for ConstantMemory).
    - ``random_seed`` (optional): Random seed (default 1).
    - ``metrics_path`` (optional): New metrics file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the downsampled file.
    """
    return _core.execute(
        _SPECS['gatk_downsample_sam'],
        {'input_path': input_path, 'output_path': output_path, 'probability': probability, 'strategy': strategy, 'accuracy': accuracy, 'random_seed': random_seed, 'metrics_path': metrics_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_position_based_downsample_sam(
    input_path: str,
    output_path: str,
    fraction: float,
    *,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Downsample reads by their physical position on the flowcell tile (keeps whole neighbourhoods of clusters, for duplicate studies).

    Run Picard `gatk PositionBasedDownsampleSam`.

    Parameters:
    - ``input_path`` (required): Coordinate-sorted input with Illumina-style read names.
    - ``output_path`` (required): New downsampled file.
    - ``fraction`` (required): Fraction of the read positions (tile areas) to keep (0-1).
    - ``stop_after`` (optional): Stop after this many reads (0 = all).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the downsampled file.
    """
    return _core.execute(
        _SPECS['gatk_position_based_downsample_sam'],
        {'input_path': input_path, 'output_path': output_path, 'fraction': fraction, 'stop_after': stop_after},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_extract_original_alignment_records_by_name(
    input_path: str,
    read_name_file: str,
    output_path: str,
    *,
    invert_match: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Subset a BAM to the reads whose names are listed in a file (with their original alignments).

    Run `gatk ExtractOriginalAlignmentRecordsByNameSpark` (beta).

    Parameters:
    - ``input_path`` (required): Input BAM.
    - ``read_name_file`` (required): Text file with one read name per line.
    - ``output_path`` (required): New BAM with the selected reads.
    - ``invert_match`` (optional): Keep the reads that are NOT listed.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the subset.
    """
    return _core.execute(
        _SPECS['gatk_extract_original_alignment_records_by_name'],
        {'input_path': input_path, 'read_name_file': read_name_file, 'output_path': output_path, 'invert_match': invert_match},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_filter_sam_reads(
    input_path: str,
    output_path: str,
    filter: str,
    *,
    read_list_file: str | None = None,
    interval_list: str | None = None,
    tag: str | None = None,
    tag_values: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Subset the reads of a SAM/BAM/CRAM file by aligned status, read name, interval or tag value.

    Run Picard `gatk FilterSamReads`. The `includeJavascript` filter is not offered because it would execute a user script.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New filtered file.
    - ``filter`` (required): Which filter to apply.
    - ``read_list_file`` (optional): Read names for the read-list filters.
    - ``interval_list`` (optional): Intervals for `includePairedIntervals`.
    - ``tag`` (optional): Tag name for the tag-value filters.
    - ``tag_values`` (optional): Tag values for the tag-value filters.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the filtered file.
    """
    return _core.execute(
        _SPECS['gatk_filter_sam_reads'],
        {'input_path': input_path, 'output_path': output_path, 'filter': filter, 'read_list_file': read_list_file, 'interval_list': interval_list, 'tag': tag, 'tag_values': tag_values},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_fix_misencoded_base_quality_reads(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fix Illumina base quality scores that were stored with the old Phred+64 encoding.

    Run `gatk FixMisencodedBaseQualityReads`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM with Illumina-encoded (Phred+64) qualities.
    - ``output_path`` (required): New file with Phred+33 qualities.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the fixed file.
    """
    return _core.execute(
        _SPECS['gatk_fix_misencoded_base_quality_reads'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_left_align_indels(
    input_path: str,
    reference_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Left-align the indels in the CIGARs of the reads of a SAM/BAM/CRAM file.

    Run `gatk LeftAlignIndels`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New file with left-aligned indels.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_left_align_indels'],
        {'input_path': input_path, 'reference_path': reference_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_mark_duplicates_with_mate_cigar(
    input_path: str,
    output_path: str,
    metrics_path: str,
    *,
    variant: str = 'standard',
    umi_metrics_path: str | None = None,
    remove_duplicates: bool | None = None,
    umi_tag_name: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mark duplicate reads using the mate CIGAR (standard, simplified or UMI-aware).

    Run Picard `MarkDuplicatesWithMateCigar`, `SimpleMarkDuplicatesWithMateCigar` or `UmiAwareMarkDuplicatesWithMateCigar` depending on `variant`; the metrics file is parsed into `tables.metrics_path`.

    ``variant``: `standard` = MarkDuplicatesWithMateCigar, `simple` = SimpleMarkDuplicatesWithMateCigar (experimental), `umi_aware` = UmiAwareMarkDuplicatesWithMateCigar (experimental).

    Parameters:
    - ``input_path`` (required): Coordinate-sorted BAM with MC (mate CIGAR) tags.
    - ``output_path`` (required): New BAM with duplicates marked.
    - ``metrics_path`` (required): New duplication metrics file.
    - ``umi_metrics_path`` (required for some modes): New UMI metrics file.
    - ``remove_duplicates`` (optional): Drop duplicates instead of flagging them.
    - ``umi_tag_name`` (optional): Read tag that holds the UMI (default RX).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields, `mode` and `tables.metrics_path` (the Picard metrics rows).
    """
    return _core.execute(
        _SPECS['gatk_mark_duplicates_with_mate_cigar'],
        {'input_path': input_path, 'output_path': output_path, 'metrics_path': metrics_path, 'umi_metrics_path': umi_metrics_path, 'remove_duplicates': remove_duplicates, 'umi_tag_name': umi_tag_name},
        mode=variant,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_merge_bam_alignment(
    unmapped_bam: str,
    reference_path: str,
    output_path: str,
    *,
    aligned_bams: Sequence[str] | None = None,
    paired_run: bool | None = None,
    sort_order: str | None = None,
    add_mate_cigar: bool | None = None,
    clip_adapters: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the alignment data of an aligner BAM with the original data of the unmapped BAM.

    Run Picard `gatk MergeBamAlignment`.

    Parameters:
    - ``unmapped_bam`` (required): Original unmapped BAM (query-name sorted).
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New merged BAM.
    - ``aligned_bams`` (optional): Aligner output files (query-name sorted); omit for unmapped-only output.
    - ``paired_run`` (optional): Whether the run is paired-end.
    - ``sort_order`` (optional): Sort order of the output.
    - ``add_mate_cigar`` (optional): Add the MC tag.
    - ``clip_adapters`` (optional): Clip adapters marked by MarkIlluminaAdapters.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the merged BAM.
    """
    return _core.execute(
        _SPECS['gatk_merge_bam_alignment'],
        {'unmapped_bam': unmapped_bam, 'reference_path': reference_path, 'output_path': output_path, 'aligned_bams': aligned_bams, 'paired_run': paired_run, 'sort_order': sort_order, 'add_mate_cigar': add_mate_cigar, 'clip_adapters': clip_adapters},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_print_distant_mates(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unmap reads whose mates are mapped far away and write them out.

    Run `gatk PrintDistantMates`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM.
    - ``output_path`` (required): New file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_print_distant_mates'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_reorder_sam(
    input_path: str,
    output_path: str,
    sequence_dictionary: str,
    *,
    allow_incomplete_dict_concordance: bool | None = None,
    allow_contig_length_discordance: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reorder the reads of a SAM/BAM/CRAM file to match the contig order of another reference.

    Run Picard `gatk ReorderSam`.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New reordered file.
    - ``sequence_dictionary`` (required): SAM/BAM/CRAM/dict/fasta whose contig order is the target.
    - ``allow_incomplete_dict_concordance`` (optional): Allow contigs missing from the target.
    - ``allow_contig_length_discordance`` (optional): Allow contigs of different length.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the reordered file.
    """
    return _core.execute(
        _SPECS['gatk_reorder_sam'],
        {'input_path': input_path, 'output_path': output_path, 'sequence_dictionary': sequence_dictionary, 'allow_incomplete_dict_concordance': allow_incomplete_dict_concordance, 'allow_contig_length_discordance': allow_contig_length_discordance},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_revert_base_quality_scores(
    input_path: str,
    output_path: str,
    *,
    mode: str = 'original_qualities',
    restore_original_qualities: bool | None = None,
    sort_order: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Revert base qualities to the original values stored in the OQ tag (optionally adding the mate CIGAR tag).

    Run `gatk RevertBaseQualityScores` (`mode=original_qualities`) or Picard `RevertOriginalBaseQualitiesAndAddMateCigar` (`mode=with_mate_cigar`).

    ``mode``: `original_qualities` = RevertBaseQualityScores, `with_mate_cigar` = RevertOriginalBaseQualitiesAndAddMateCigar.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM that carries OQ tags.
    - ``output_path`` (required): New file.
    - ``restore_original_qualities`` (optional): Restore the OQ values.
    - ``sort_order`` (optional): Sort order of the output.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_revert_base_quality_scores'],
        {'input_path': input_path, 'output_path': output_path, 'restore_original_qualities': restore_original_qualities, 'sort_order': sort_order},
        mode=mode,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_revert_sam(
    input_path: str,
    output_path: str,
    *,
    restore_original_qualities: bool | None = None,
    remove_duplicate_information: bool | None = None,
    remove_alignment_information: bool | None = None,
    sanitize: bool | None = None,
    sample_alias: str | None = None,
    library_name: str | None = None,
    sort_order: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Revert a SAM/BAM/CRAM file to an earlier state (unmapped, original qualities, no duplicate flags) for re-alignment.

    Run Picard `gatk RevertSam` with a single output file; `OUTPUT_BY_READGROUP` is not offered by this function.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New reverted file.
    - ``restore_original_qualities`` (optional): Restore OQ qualities (Picard default true).
    - ``remove_duplicate_information`` (optional): Clear the duplicate flag.
    - ``remove_alignment_information`` (optional): Remove alignment information.
    - ``sanitize`` (optional): Remove reads that cannot be reverted consistently.
    - ``sample_alias`` (optional): Sample name for the read groups.
    - ``library_name`` (optional): Library name for the read groups.
    - ``sort_order`` (optional): Sort order of the output.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the reverted file.
    """
    return _core.execute(
        _SPECS['gatk_revert_sam'],
        {'input_path': input_path, 'output_path': output_path, 'restore_original_qualities': restore_original_qualities, 'remove_duplicate_information': remove_duplicate_information, 'remove_alignment_information': remove_alignment_information, 'sanitize': sanitize, 'sample_alias': sample_alias, 'library_name': library_name, 'sort_order': sort_order},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_convert_sam_format(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a SAM file to BAM or a BAM file to SAM (or CRAM).

    Run Picard `gatk SamFormatConverter`.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``output_path`` (required): New file; the format follows its extension (.sam, .bam, .cram).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the converted file.
    """
    return _core.execute(
        _SPECS['gatk_convert_sam_format'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sam_to_fastq_with_tags(
    input_path: str,
    fastq_path: str,
    sequence_tag_groups: Sequence[str],
    *,
    quality_tag_groups: Sequence[str] | None = None,
    second_end_fastq: str | None = None,
    unpaired_fastq: str | None = None,
    interleave: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a SAM/BAM/CRAM file to FASTQ and also write FASTQ files built from read tags (for example UMIs).

    Run Picard `gatk SamToFastqWithTags`.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM.
    - ``fastq_path`` (required): New FASTQ for the first read (or all reads when interleaved).
    - ``sequence_tag_groups`` (required): Tags whose values are written as extra FASTQ files; join several tags with commas.
    - ``quality_tag_groups`` (optional): Quality tags that go with the sequence tag groups.
    - ``second_end_fastq`` (optional): New FASTQ for the second read.
    - ``unpaired_fastq`` (optional): New FASTQ for unpaired reads.
    - ``interleave`` (optional): Write both reads of a pair to one FASTQ.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs` lists the FASTQ files that were requested.
    """
    return _core.execute(
        _SPECS['gatk_sam_to_fastq_with_tags'],
        {'input_path': input_path, 'fastq_path': fastq_path, 'sequence_tag_groups': sequence_tag_groups, 'quality_tag_groups': quality_tag_groups, 'second_end_fastq': second_end_fastq, 'unpaired_fastq': unpaired_fastq, 'interleave': interleave},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_split_n_cigar_reads(
    input_path: str,
    reference_path: str,
    output_path: str,
    *,
    max_bases_in_overhang: int | None = None,
    max_mismatches_in_overhang: int | None = None,
    skip_mapping_quality_transform: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split RNA-seq reads at the N operators of their CIGARs into exon segments and clip overhangs.

    Run `gatk SplitNCigarReads`.

    Parameters:
    - ``input_path`` (required): RNA-seq BAM/CRAM/SAM with N operators in CIGARs.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New file.
    - ``max_bases_in_overhang`` (optional): Maximum overhang bases before a split piece is hard-clipped.
    - ``max_mismatches_in_overhang`` (optional): Maximum mismatches allowed in an overhang.
    - ``skip_mapping_quality_transform`` (optional): Do not change MAPQ 255 to 60.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the split file.
    """
    return _core.execute(
        _SPECS['gatk_split_n_cigar_reads'],
        {'input_path': input_path, 'reference_path': reference_path, 'output_path': output_path, 'max_bases_in_overhang': max_bases_in_overhang, 'max_mismatches_in_overhang': max_mismatches_in_overhang, 'skip_mapping_quality_transform': skip_mapping_quality_transform},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_split_sam(
    input_path: str,
    output_dir: str,
    *,
    by: str = 'library',
    split_to_n_files: int | None = None,
    split_to_n_reads: int | None = None,
    total_reads_in_input: int | None = None,
    out_prefix: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a SAM/BAM/CRAM file into several files, by library name or by number of reads.

    Run Picard `SplitSamByLibrary` (`by=library`, reads without a library go to `unknown`) or `SplitSamByNumberOfReads` (`by=number_of_reads`); the files written are listed in `files`.

    ``by``: `library` = SplitSamByLibrary, `number_of_reads` = SplitSamByNumberOfReads.

    Parameters:
    - ``input_path`` (required): Input SAM/BAM/CRAM (query-grouped for `number_of_reads`).
    - ``output_dir`` (required): New directory for the split files.
    - ``split_to_n_files`` (optional): Number of output files.
    - ``split_to_n_reads`` (optional): Reads per output file.
    - ``total_reads_in_input`` (optional): Total reads, which saves a counting pass.
    - ``out_prefix`` (optional): File name prefix (default `shard`).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields, `mode` and `files` (file names inside `output_dir`).
    """
    return _core.execute(
        _SPECS['gatk_split_sam'],
        {'input_path': input_path, 'output_dir': output_dir, 'split_to_n_files': split_to_n_files, 'split_to_n_reads': split_to_n_reads, 'total_reads_in_input': total_reads_in_input, 'out_prefix': out_prefix},
        mode=by,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_transfer_read_tags(
    input_path: str,
    unmapped_sam: str,
    output_path: str,
    read_tags: Sequence[str],
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy read tags (for example UMIs) from an unmapped SAM/BAM onto the matching reads of an aligned file.

    Run `gatk TransferReadTags` (experimental).

    Parameters:
    - ``input_path`` (required): Aligned BAM/SAM (query-name sorted) that receives the tags.
    - ``unmapped_sam`` (required): Query-name sorted unmapped SAM/BAM that holds the tags.
    - ``output_path`` (required): New output file.
    - ``read_tags`` (required): Tag names to transfer.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_transfer_read_tags'],
        {'input_path': input_path, 'unmapped_sam': unmapped_sam, 'output_path': output_path, 'read_tags': read_tags},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_unmark_duplicates(
    input_path: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Clear the duplicate flag (0x400) of all reads.

    Run `gatk UnmarkDuplicates`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM/SAM.
    - ``output_path`` (required): New file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the new file.
    """
    return _core.execute(
        _SPECS['gatk_unmark_duplicates'],
        {'input_path': input_path, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


__all__ = [
    "gatk_pathseq",
    "gatk_detect_cram_issue_8768",
    "gatk_create_bam_splitting_index",
    "gatk_print_bgzf_block_information",
    "gatk_print_file_diagnostics",
    "gatk_read_anonymizer",
    "gatk_add_comments_to_bam",
    "gatk_add_original_alignment_tags",
    "gatk_base_recalibrator",
    "gatk_apply_bqsr",
    "gatk_bqsr_pipeline_spark",
    "gatk_bam_to_bfq",
    "gatk_bwa_spark",
    "gatk_clean_sam",
    "gatk_clip_reads",
    "gatk_convert_headerless_bam_shard",
    "gatk_downsample_by_duplicate_set",
    "gatk_downsample_sam",
    "gatk_position_based_downsample_sam",
    "gatk_extract_original_alignment_records_by_name",
    "gatk_filter_sam_reads",
    "gatk_fix_misencoded_base_quality_reads",
    "gatk_left_align_indels",
    "gatk_mark_duplicates_with_mate_cigar",
    "gatk_merge_bam_alignment",
    "gatk_print_distant_mates",
    "gatk_reorder_sam",
    "gatk_revert_base_quality_scores",
    "gatk_revert_sam",
    "gatk_convert_sam_format",
    "gatk_sam_to_fastq_with_tags",
    "gatk_split_n_cigar_reads",
    "gatk_split_sam",
    "gatk_transfer_read_tags",
    "gatk_unmark_duplicates",
]

_register_public_verbose(globals(), tuple(__all__))
