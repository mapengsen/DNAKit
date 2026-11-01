"""GATK variant-calling, GVCF, structural-variant and evaluation command-line facades (pack P01).

Each public function runs one GATK 4.7.0.0 tool (or one family of tools selected by a mode argument) through
``dnakit.backends.external.execute_bounded_command`` with a fixed option whitelist and no shell. The user installs GATK
(conda ``gatk4`` or the release zip) with Java 17 and either puts the ``gatk`` launcher on PATH, sets ``DNAKIT_GATK`` or
passes ``gatk_path``. Outputs are written to explicit new paths that must not exist; nothing is overwritten. Results are
JSON-serialisable dictionaries with status, backend, tool, version, command, outputs, a bounded log tail and provenance.
30 public functions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gatk_p01_core as _core
from dnakit.bioinformatics._gatk_p01_core import Param as _Param
from dnakit.bioinformatics._gatk_p01_core import Spec as _Spec

_SPECS: dict[str, _Spec] = {
    'gatk_reblock_gvcf': _Spec(
        'gatk_reblock_gvcf',
        {'default': 'ReblockGVCF'},
        (
            _Param('gvcf_path', 'in', 'variant', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('gvcf_gq_bands', 'strs', 'gvcf-gq-bands'),
            _Param('drop_low_quals', 'bool', 'drop-low-quals'),
            _Param('do_qual_score_approximation', 'bool', 'do-qual-score-approximation'),
        ),
    ),
    'gatk_calibrate_dragstr_model': _Spec(
        'gatk_calibrate_dragstr_model',
        {'default': 'CalibrateDragstrModel'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('str_table_path', 'in', 'str-table-path', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('threads', 'int', 'threads', minimum=1),
            _Param('min_loci_count', 'int', 'min-loci-count', minimum=1),
        ),
    ),
    'gatk_combine_gvcfs': _Spec(
        'gatk_combine_gvcfs',
        {'default': 'CombineGVCFs'},
        (
            _Param('gvcf_paths', 'ins', 'variant', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('convert_to_base_pair_resolution', 'bool', 'convert-to-base-pair-resolution'),
            _Param('break_bands_at_multiples_of', 'int', 'break-bands-at-multiples-of', minimum=0),
        ),
    ),
    'gatk_genomicsdb_import': _Spec(
        'gatk_genomicsdb_import',
        {'create': 'GenomicsDBImport', 'update': 'GenomicsDBImport'},
        (
            _Param('workspace_path', 'outd_tool', {'create': 'genomicsdb-workspace-path'}, v=('create',), req=('create',)),
            _Param('existing_workspace_path', 'in_dir', {'update': 'genomicsdb-update-workspace-path'}, v=('update',), req=('update',)),
            _Param('gvcf_paths', 'ins', 'variant'),
            _Param('sample_name_map', 'in', 'sample-name-map'),
            _Param('intervals', 'strs', 'intervals', v=('create',), req=('create',)),
            _Param('batch_size', 'int', 'batch-size', minimum=0),
        ),
        selector='mode',
        any_of=((None, ('gvcf_paths', 'sample_name_map')),),
    ),
    'gatk_genotype_gvcfs': _Spec(
        'gatk_genotype_gvcfs',
        {'standard': 'GenotypeGVCFs', 'gnarly': 'GnarlyGenotyper'},
        (
            _Param('variant_source', 'src', 'variant', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('dbsnp', 'in', 'dbsnp'),
            _Param('include_non_variant_sites', 'bool', 'include-non-variant-sites', v=('standard',)),
            _Param('keep_all_sites', 'bool', 'keep-all-sites', v=('gnarly',)),
        ),
        selector='engine',
    ),
    'gatk_haplotype_based_variant_recaller': _Spec(
        'gatk_haplotype_based_variant_recaller',
        {'default': 'HaplotypeBasedVariantRecaller'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('alleles_vcf', 'in', 'alleles-file-vcf', req=True),
            _Param('haplotypes_bam', 'in', 'haplotypes-file-bam', req=True),
            _Param('matrix_csv', 'out', 'matrix-file-csv', req=True),
        ),
    ),
    'gatk_learn_read_orientation_model': _Spec(
        'gatk_learn_read_orientation_model',
        {'default': 'LearnReadOrientationModel'},
        (
            _Param('f1r2_tar_gz_files', 'ins', 'input', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('num_em_iterations', 'int', 'num-em-iterations', minimum=1),
            _Param('convergence_threshold', 'float', 'convergence-threshold', minimum=0),
        ),
    ),
    'gatk_mutect2': _Spec(
        'gatk_mutect2',
        {'default': 'Mutect2'},
        (
            _Param('input_paths', 'ins', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('tumor_sample', 'str', 'tumor-sample'),
            _Param('normal_samples', 'strs', 'normal-sample'),
            _Param('germline_resource', 'in', 'germline-resource'),
            _Param('panel_of_normals', 'in', 'panel-of-normals'),
            _Param('f1r2_tar_gz', 'out', 'f1r2-tar-gz'),
            _Param('bam_output', 'out', 'bam-output'),
        ),
    ),
    'gatk_reads_pipeline_spark': _Spec(
        'gatk_reads_pipeline_spark',
        {'default': 'ReadsPipelineSpark'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('known_sites', 'ins', 'known-sites', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('output_bam', 'out', 'output-bam'),
            _Param('align', 'bool', 'align'),
            _Param('emit_ref_confidence', 'choice', 'emit-ref-confidence', choices=('NONE', 'BP_RESOLUTION', 'GVCF')),
        ),
    ),
    'gatk_collect_sv_evidence': _Spec(
        'gatk_collect_sv_evidence',
        {'default': 'CollectSVEvidence'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('sample_name', 'str', 'sample-name', req=True),
            _Param('pe_file', 'out', 'pe-file'),
            _Param('sr_file', 'out', 'sr-file'),
            _Param('sd_file', 'out', 'sd-file'),
            _Param('site_depth_locs_vcf', 'in', 'site-depth-locs-vcf'),
            _Param('depth_evidence_file', 'out', 'depth-evidence-file'),
        ),
        any_of=((None, ('pe_file', 'sr_file', 'sd_file', 'depth_evidence_file')),),
    ),
    'gatk_condense_depth_evidence': _Spec(
        'gatk_condense_depth_evidence',
        {'default': 'CondenseDepthEvidence'},
        (
            _Param('depth_evidence', 'in', 'depth-evidence', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('max_interval_size', 'int', 'max-interval-size', minimum=1),
            _Param('min_interval_size', 'int', 'min-interval-size', minimum=1),
        ),
    ),
    'gatk_convert_counts_to_depth_file': _Spec(
        'gatk_convert_counts_to_depth_file',
        {'default': 'ConvertCountsToDepthFile'},
        (
            _Param('counts_path', 'in', 'counts-filename', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('sample_name', 'str', 'sample-name'),
        ),
    ),
    'gatk_sv_cluster': _Spec(
        'gatk_sv_cluster',
        {'cluster': 'SVCluster', 'grouped': 'GroupedSVCluster'},
        (
            _Param('variant_paths', 'ins', 'variant', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('ploidy_table', 'in', 'ploidy-table', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('clustering_config', 'in', 'clustering-config', v=('grouped',), req=('grouped',)),
            _Param('stratify_config', 'in', 'stratify-config', v=('grouped',), req=('grouped',)),
            _Param('algorithm', 'choice', 'algorithm', choices=('SINGLE_LINKAGE', 'MAX_CLIQUE')),
            _Param('variant_prefix', 'str', 'variant-prefix'),
        ),
        selector='mode',
    ),
    'gatk_joint_germline_cnv_segmentation': _Spec(
        'gatk_joint_germline_cnv_segmentation',
        {'default': 'JointGermlineCNVSegmentation'},
        (
            _Param('variant_paths', 'ins', 'variant', req=True),
            _Param('pedigree', 'in', 'pedigree', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('model_call_intervals', 'in', 'model-call-intervals'),
            _Param('min_sample_set_fraction_overlap', 'float', 'min-sample-set-fraction-overlap', minimum=0),
        ),
    ),
    'gatk_print_read_counts': _Spec(
        'gatk_print_read_counts',
        {'default': 'PrintReadCounts'},
        (
            _Param('counts_path', 'in', 'input-counts', req=True),
            _Param('output_prefix', 'str', 'output-prefix'),
            _Param('output_file_list', 'out', 'output-file-list'),
        ),
        any_of=((None, ('output_prefix', 'output_file_list')),),
    ),
    'gatk_print_sv_evidence': _Spec(
        'gatk_print_sv_evidence',
        {'default': 'PrintSVEvidence'},
        (
            _Param('evidence_file', 'in', 'evidence-file', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('sample_names', 'strs', 'sample-names'),
        ),
    ),
    'gatk_sv_annotate': _Spec(
        'gatk_sv_annotate',
        {'default': 'SVAnnotate'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('protein_coding_gtf', 'in', 'protein-coding-gtf'),
            _Param('non_coding_bed', 'in', 'non-coding-bed'),
            _Param('promoter_window_length', 'int', 'promoter-window-length', minimum=0),
        ),
    ),
    'gatk_sv_concordance': _Spec(
        'gatk_sv_concordance',
        {'default': 'SVConcordance'},
        (
            _Param('evaluation', 'in', 'evaluation', req=True),
            _Param('truth', 'in', 'truth', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('clustering_config', 'in', 'clustering-config'),
            _Param('stratify_config', 'in', 'stratify-config'),
        ),
    ),
    'gatk_sv_stratify': _Spec(
        'gatk_sv_stratify',
        {'default': 'SVStratify'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('stratify_config', 'in', 'stratify-config', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('track_name', 'strs', 'track-name'),
            _Param('track_intervals', 'ins', 'track-intervals'),
            _Param('split_output', 'bool', 'split-output'),
        ),
    ),
    'gatk_site_depth_to_baf': _Spec(
        'gatk_site_depth_to_baf',
        {'default': 'SiteDepthtoBAF'},
        (
            _Param('site_depth_paths', 'ins', 'site-depth', req=True),
            _Param('baf_sites_vcf', 'in', 'baf-sites-vcf', req=True),
            _Param('baf_evidence_output', 'out', 'baf-evidence-output', req=True),
            _Param('sample_names', 'strs', 'sample-names'),
            _Param('min_het_probability', 'float', 'min-het-probability', minimum=0),
        ),
    ),
    'gatk_structural_variation_discovery_pipeline_spark': _Spec(
        'gatk_structural_variation_discovery_pipeline_spark',
        {'default': 'StructuralVariationDiscoveryPipelineSpark'},
        (
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('aligner_index_image', 'in', 'aligner-index-image', req=True),
            _Param('kmers_to_ignore', 'in', 'kmers-to-ignore', req=True),
            _Param('contig_sam_file', 'out', 'contig-sam-file', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_calculate_genotype_posteriors': _Spec(
        'gatk_calculate_genotype_posteriors',
        {'default': 'CalculateGenotypePosteriors'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('pedigree', 'in', 'pedigree'),
            _Param('supporting_callsets', 'ins', 'supporting-callsets'),
            _Param('de_novo_prior', 'float', 'de-novo-prior', minimum=0),
        ),
    ),
    'gatk_concordance': _Spec(
        'gatk_concordance',
        {'default': 'Concordance'},
        (
            _Param('evaluation', 'in', 'evaluation', req=True),
            _Param('truth', 'in', 'truth', req=True),
            _Param('summary_path', 'out', 'summary', req=True),
            _Param('true_positives_and_false_negatives', 'out', 'true-positives-and-false-negatives'),
            _Param('true_positives_and_false_positives', 'out', 'true-positives-and-false-positives'),
        ),
        tables={'summary_path': 'tsv'},
    ),
    'gatk_count_false_positives': _Spec(
        'gatk_count_false_positives',
        {'default': 'CountFalsePositives'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('intervals', 'strs', 'intervals', req=True),
            _Param('output_path', 'out', 'output', req=True),
        ),
    ),
    'gatk_evaluate_info_field_concordance': _Spec(
        'gatk_evaluate_info_field_concordance',
        {'default': 'EvaluateInfoFieldConcordance'},
        (
            _Param('evaluation', 'in', 'evaluation', req=True),
            _Param('truth', 'in', 'truth', req=True),
            _Param('eval_info_key', 'str', 'eval-info-key', req=True),
            _Param('truth_info_key', 'str', 'truth-info-key', req=True),
            _Param('summary_path', 'out', 'summary', req=True),
            _Param('epsilon', 'float', 'epsilon', minimum=0),
        ),
    ),
    'gatk_reference_block_concordance': _Spec(
        'gatk_reference_block_concordance',
        {'default': 'ReferenceBlockConcordance'},
        (
            _Param('evaluation', 'in', 'evaluation', req=True),
            _Param('truth', 'in', 'truth', req=True),
            _Param('eval_block_histogram', 'out', 'eval-block-histogram', req=True),
            _Param('truth_block_histogram', 'out', 'truth-block-histogram', req=True),
            _Param('confidence_concordance_histogram', 'out', 'confidence-concordance-histogram', req=True),
        ),
    ),
    'gatk_filter_funcotations': _Spec(
        'gatk_filter_funcotations',
        {'default': 'FilterFuncotations'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('ref_version', 'choice', 'ref-version', req=True, choices=('hg19', 'hg38')),
            _Param('allele_frequency_data_source', 'choice', 'allele-frequency-data-source', req=True, choices=('gnomad', 'exac')),
        ),
    ),
    'gatk_find_mendelian_violations': _Spec(
        'gatk_find_mendelian_violations',
        {'default': 'FindMendelianViolations'},
        (
            _Param('input_path', 'in', 'INPUT', req=True),
            _Param('pedigree_path', 'in', 'TRIOS', req=True),
            _Param('output_path', 'out', 'OUTPUT', req=True),
            _Param('min_gq', 'int', 'MIN_GQ', minimum=0),
            _Param('min_dp', 'int', 'MIN_DP', minimum=0),
            _Param('skip_chroms', 'strs', 'SKIP_CHROMS'),
        ),
        tables={'output_path': 'metrics'},
    ),
    'gatk_funcotate': _Spec(
        'gatk_funcotate',
        {'variants': 'Funcotator', 'segments': 'FuncotateSegments'},
        (
            _Param('input_path', 'in', {'variants': 'variant', 'segments': 'segments'}, req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('data_sources_path', 'in_dir', 'data-sources-path', req=True),
            _Param('ref_version', 'str', 'ref-version', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('output_file_format', 'choice', 'output-file-format', req=True, choices=('VCF', 'MAF', 'SEG')),
            _Param('transcript_selection_mode', 'choice', 'transcript-selection-mode', choices=('BEST_EFFECT', 'CANONICAL', 'ALL')),
        ),
        selector='input_kind',
    ),
    'gatk_validate_basic_somatic_short_mutations': _Spec(
        'gatk_validate_basic_somatic_short_mutations',
        {'default': 'ValidateBasicSomaticShortMutations'},
        (
            _Param('variant_path', 'in', 'variant', req=True),
            _Param('input_path', 'in', 'input', req=True),
            _Param('reference_path', 'in', 'reference', req=True),
            _Param('discovery_sample_name', 'str', 'discovery-sample-name', req=True),
            _Param('val_case_sample_name', 'str', 'val-case-sample-name', req=True),
            _Param('val_control_sample_name', 'str', 'val-control-sample-name', req=True),
            _Param('output_path', 'out', 'output', req=True),
            _Param('summary_path', 'out', 'summary'),
        ),
    ),
}


def gatk_reblock_gvcf(
    gvcf_path: str,
    reference_path: str,
    output_path: str,
    *,
    gvcf_gq_bands: Sequence[str] | None = None,
    drop_low_quals: bool | None = None,
    do_qual_score_approximation: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Condense the homozygous-reference blocks of a single-sample GVCF into fewer, wider blocks.

    Run `gatk ReblockGVCF`.

    Parameters:
    - ``gvcf_path`` (required): Single-sample GVCF (indexed).
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New reblocked GVCF.
    - ``gvcf_gq_bands`` (optional): GQ band boundaries.
    - ``drop_low_quals`` (optional): Drop low-quality variants and non-ref blocks.
    - ``do_qual_score_approximation`` (optional): Approximate QUAL from the allele-specific annotations.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the reblocked GVCF.
    """
    return _core.execute(
        _SPECS['gatk_reblock_gvcf'],
        {'gvcf_path': gvcf_path, 'reference_path': reference_path, 'output_path': output_path, 'gvcf_gq_bands': gvcf_gq_bands, 'drop_low_quals': drop_low_quals, 'do_qual_score_approximation': do_qual_score_approximation},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_calibrate_dragstr_model(
    input_path: str,
    reference_path: str,
    str_table_path: str,
    output_path: str,
    *,
    threads: int | None = None,
    min_loci_count: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the parameters of the DRAGstr model (STR-aware indel error model) from a BAM.

    Run `gatk CalibrateDragstrModel`.

    Parameters:
    - ``input_path`` (required): Input BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``str_table_path`` (required): STR table from `gatk_compose_str_table`.
    - ``output_path`` (required): New DRAGstr parameter file.
    - ``threads`` (optional): Number of threads.
    - ``min_loci_count`` (optional): Minimum loci per period/repeat combination.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the parameter file.
    """
    return _core.execute(
        _SPECS['gatk_calibrate_dragstr_model'],
        {'input_path': input_path, 'reference_path': reference_path, 'str_table_path': str_table_path, 'output_path': output_path, 'threads': threads, 'min_loci_count': min_loci_count},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_combine_gvcfs(
    gvcf_paths: Sequence[str],
    reference_path: str,
    output_path: str,
    *,
    convert_to_base_pair_resolution: bool | None = None,
    break_bands_at_multiples_of: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge one or more HaplotypeCaller GVCF files into a single GVCF.

    Run `gatk CombineGVCFs`.

    Parameters:
    - ``gvcf_paths`` (required): GVCF files (indexed) to merge.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New combined GVCF.
    - ``convert_to_base_pair_resolution`` (optional): Write one record per base.
    - ``break_bands_at_multiples_of`` (optional): Break reference bands at multiples of this.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the merged GVCF.
    """
    return _core.execute(
        _SPECS['gatk_combine_gvcfs'],
        {'gvcf_paths': gvcf_paths, 'reference_path': reference_path, 'output_path': output_path, 'convert_to_base_pair_resolution': convert_to_base_pair_resolution, 'break_bands_at_multiples_of': break_bands_at_multiples_of},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_genomicsdb_import(
    *,
    mode: str = 'create',
    workspace_path: str | None = None,
    existing_workspace_path: str | None = None,
    gvcf_paths: Sequence[str] | None = None,
    sample_name_map: str | None = None,
    intervals: Sequence[str] | None = None,
    batch_size: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Import single-sample GVCFs into a GenomicsDB workspace for joint genotyping.

    Run `gatk GenomicsDBImport` to create a workspace (`mode=create`) or to add samples to an existing one (`mode=update`).

    ``mode``: `create` = new workspace, `update` = add samples to an existing workspace.

    Parameters:
    - ``workspace_path`` (required for some modes): New GenomicsDB workspace directory (created by GATK).
    - ``existing_workspace_path`` (required for some modes): Existing workspace to add samples to.
    - ``gvcf_paths`` (optional): GVCF files (indexed) to import; use this or `sample_name_map`.
    - ``sample_name_map`` (optional): Two-column file `sample<TAB>gvcf path`.
    - ``intervals`` (required for some modes): Intervals to import (required for `create`).
    - ``batch_size`` (optional): Samples read per batch (0 = all).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.workspace_path` for a new workspace, with `files` listing it.
    """
    return _core.execute(
        _SPECS['gatk_genomicsdb_import'],
        {'workspace_path': workspace_path, 'existing_workspace_path': existing_workspace_path, 'gvcf_paths': gvcf_paths, 'sample_name_map': sample_name_map, 'intervals': intervals, 'batch_size': batch_size},
        mode=mode,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_genotype_gvcfs(
    variant_source: str,
    reference_path: str,
    output_path: str,
    *,
    engine: str = 'standard',
    dbsnp: str | None = None,
    include_non_variant_sites: bool | None = None,
    keep_all_sites: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Perform joint genotyping of GVCFs (standard or the fast GnarlyGenotyper model).

    Run `gatk GenotypeGVCFs` (`engine=standard`) or the beta `gatk GnarlyGenotyper` (`engine=gnarly`, expects a GenomicsDB workspace or a combined GVCF with allele-specific annotations); give the intervals through the `intervals` option.

    ``engine``: `standard` = GenotypeGVCFs, `gnarly` = GnarlyGenotyper.

    Parameters:
    - ``variant_source`` (required): Combined GVCF, or `gendb://<workspace>` for a GenomicsDB workspace.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New genotyped VCF.
    - ``dbsnp`` (optional): dbSNP VCF for the ID column.
    - ``include_non_variant_sites`` (optional): Also output non-variant sites.
    - ``keep_all_sites`` (optional): Keep all sites, also those without alternate alleles.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the VCF.
    """
    return _core.execute(
        _SPECS['gatk_genotype_gvcfs'],
        {'variant_source': variant_source, 'reference_path': reference_path, 'output_path': output_path, 'dbsnp': dbsnp, 'include_non_variant_sites': include_non_variant_sites, 'keep_all_sites': keep_all_sites},
        mode=engine,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_haplotype_based_variant_recaller(
    input_path: str,
    reference_path: str,
    alleles_vcf: str,
    haplotypes_bam: str,
    matrix_csv: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate the likelihood of each allele of a VCF against the reads, limited by a set of haplotypes.

    Run `gatk HaplotypeBasedVariantRecaller` (experimental).

    Parameters:
    - ``input_path`` (required): Reads BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``alleles_vcf`` (required): VCF of the alleles to evaluate.
    - ``haplotypes_bam`` (required): BAM of the haplotypes.
    - ``matrix_csv`` (required): New CSV with the likelihood matrix.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.matrix_csv` is the matrix.
    """
    return _core.execute(
        _SPECS['gatk_haplotype_based_variant_recaller'],
        {'input_path': input_path, 'reference_path': reference_path, 'alleles_vcf': alleles_vcf, 'haplotypes_bam': haplotypes_bam, 'matrix_csv': matrix_csv},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_learn_read_orientation_model(
    f1r2_tar_gz_files: Sequence[str],
    output_path: str,
    *,
    num_em_iterations: int | None = None,
    convergence_threshold: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Get maximum-likelihood estimates of the artifact priors of the read-orientation bias mixture model used by FilterMutectCalls.

    Run `gatk LearnReadOrientationModel`.

    Parameters:
    - ``f1r2_tar_gz_files`` (required): F1R2 counts archives written by Mutect2 (`f1r2-tar-gz`).
    - ``output_path`` (required): New orientation-bias model archive (.tar.gz).
    - ``num_em_iterations`` (optional): Maximum EM iterations.
    - ``convergence_threshold`` (optional): EM convergence threshold.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the model archive.
    """
    return _core.execute(
        _SPECS['gatk_learn_read_orientation_model'],
        {'f1r2_tar_gz_files': f1r2_tar_gz_files, 'output_path': output_path, 'num_em_iterations': num_em_iterations, 'convergence_threshold': convergence_threshold},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_mutect2(
    input_paths: Sequence[str],
    reference_path: str,
    output_path: str,
    *,
    tumor_sample: str | None = None,
    normal_samples: Sequence[str] | None = None,
    germline_resource: str | None = None,
    panel_of_normals: str | None = None,
    f1r2_tar_gz: str | None = None,
    bam_output: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call somatic SNVs and indels by local assembly of haplotypes (Mutect2).

    Run `gatk Mutect2`; give the regions through the `intervals` option.

    Parameters:
    - ``input_paths`` (required): Tumor (and normal) BAM/CRAM files.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New somatic VCF.
    - ``tumor_sample`` (optional): Tumor sample name (read-group SM).
    - ``normal_samples`` (optional): Normal sample name(s).
    - ``germline_resource`` (optional): Population allele-frequency VCF (for example gnomAD).
    - ``panel_of_normals`` (optional): Panel-of-normals VCF.
    - ``f1r2_tar_gz`` (optional): New archive of F1R2 counts for `gatk_learn_read_orientation_model`.
    - ``bam_output`` (optional): New BAM of the locally re-assembled reads.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the VCF.
    """
    return _core.execute(
        _SPECS['gatk_mutect2'],
        {'input_paths': input_paths, 'reference_path': reference_path, 'output_path': output_path, 'tumor_sample': tumor_sample, 'normal_samples': normal_samples, 'germline_resource': germline_resource, 'panel_of_normals': panel_of_normals, 'f1r2_tar_gz': f1r2_tar_gz, 'bam_output': bam_output},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_reads_pipeline_spark(
    input_path: str,
    reference_path: str,
    known_sites: Sequence[str],
    output_path: str,
    *,
    output_bam: str | None = None,
    align: bool | None = None,
    emit_ref_confidence: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run BWA (optional), MarkDuplicates, BQSR and HaplotypeCaller on reads in one Spark pipeline to produce a VCF/GVCF.

    Run `gatk ReadsPipelineSpark` (beta, Spark local mode).

    Parameters:
    - ``input_path`` (required): Unaligned or aligned BAM/CRAM.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict). A .2bit file also works for Spark.
    - ``known_sites`` (required): Known-variant VCF files (indexed) for BQSR.
    - ``output_path`` (required): New VCF/GVCF.
    - ``output_bam`` (optional): New BAM of the processed reads.
    - ``align`` (optional): Run BWA-MEM first.
    - ``emit_ref_confidence`` (optional): Reference-confidence mode of the caller.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the variant file.
    """
    return _core.execute(
        _SPECS['gatk_reads_pipeline_spark'],
        {'input_path': input_path, 'reference_path': reference_path, 'known_sites': known_sites, 'output_path': output_path, 'output_bam': output_bam, 'align': align, 'emit_ref_confidence': emit_ref_confidence},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_collect_sv_evidence(
    input_path: str,
    sample_name: str,
    *,
    pe_file: str | None = None,
    sr_file: str | None = None,
    sd_file: str | None = None,
    site_depth_locs_vcf: str | None = None,
    depth_evidence_file: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Gather paired-end, split-read and depth evidence files from a BAM for the GATK-SV pipeline.

    Run `gatk CollectSVEvidence` (beta); the evidence files that are requested are written, and at least one must be requested.

    Parameters:
    - ``input_path`` (required): Coordinate-sorted BAM/CRAM.
    - ``sample_name`` (required): Sample name written to the evidence files.
    - ``pe_file`` (optional): New paired-end (discordant pair) evidence file.
    - ``sr_file`` (optional): New split-read evidence file.
    - ``sd_file`` (optional): New site-depth evidence file.
    - ``site_depth_locs_vcf`` (optional): VCF of the sites for site-depth counting.
    - ``depth_evidence_file`` (optional): New read-depth evidence file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs` lists the evidence files.
    """
    return _core.execute(
        _SPECS['gatk_collect_sv_evidence'],
        {'input_path': input_path, 'sample_name': sample_name, 'pe_file': pe_file, 'sr_file': sr_file, 'sd_file': sd_file, 'site_depth_locs_vcf': site_depth_locs_vcf, 'depth_evidence_file': depth_evidence_file},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_condense_depth_evidence(
    depth_evidence: str,
    output_path: str,
    *,
    max_interval_size: int | None = None,
    min_interval_size: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge adjacent depth-evidence records of a GATK-SV depth file.

    Run `gatk CondenseDepthEvidence` (experimental).

    Parameters:
    - ``depth_evidence`` (required): Depth evidence file (`.rd.txt`).
    - ``output_path`` (required): New condensed depth evidence file.
    - ``max_interval_size`` (optional): Maximum merged interval size.
    - ``min_interval_size`` (optional): Minimum interval size.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the condensed file.
    """
    return _core.execute(
        _SPECS['gatk_condense_depth_evidence'],
        {'depth_evidence': depth_evidence, 'output_path': output_path, 'max_interval_size': max_interval_size, 'min_interval_size': min_interval_size},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_convert_counts_to_depth_file(
    counts_path: str,
    output_path: str,
    *,
    sample_name: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a GATK read-count file (`*.counts.tsv`) to a GATK-SV depth file (`.rd.txt`).

    Run `gatk ConvertCountsToDepthFile` (experimental).

    Parameters:
    - ``counts_path`` (required): Read-count file (`*.counts.tsv`).
    - ``output_path`` (required): New depth evidence file (`.rd.txt`).
    - ``sample_name`` (optional): Sample name written to the file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the depth file.
    """
    return _core.execute(
        _SPECS['gatk_convert_counts_to_depth_file'],
        {'counts_path': counts_path, 'output_path': output_path, 'sample_name': sample_name},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sv_cluster(
    variant_paths: Sequence[str],
    reference_path: str,
    ploidy_table: str,
    output_path: str,
    *,
    mode: str = 'cluster',
    clustering_config: str | None = None,
    stratify_config: str | None = None,
    algorithm: str | None = None,
    variant_prefix: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster structural variants by type, size and breakpoints (optionally grouped by stratum and track overlap).

    Run `gatk SVCluster` (`mode=cluster`) or `GroupedSVCluster` (`mode=grouped`, which also needs the clustering and stratification configurations); both are beta tools.

    ``mode``: `cluster` = SVCluster, `grouped` = GroupedSVCluster.

    Parameters:
    - ``variant_paths`` (required): Structural-variant VCF files.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``ploidy_table`` (required): Tab-separated sample/contig ploidy table.
    - ``output_path`` (required): New clustered VCF.
    - ``clustering_config`` (required for some modes): Clustering configuration table (per variant group).
    - ``stratify_config`` (required for some modes): Stratification configuration table.
    - ``algorithm`` (optional): Clustering algorithm.
    - ``variant_prefix`` (optional): Prefix of the new variant IDs.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the clustered VCF.
    """
    return _core.execute(
        _SPECS['gatk_sv_cluster'],
        {'variant_paths': variant_paths, 'reference_path': reference_path, 'ploidy_table': ploidy_table, 'output_path': output_path, 'clustering_config': clustering_config, 'stratify_config': stratify_config, 'algorithm': algorithm, 'variant_prefix': variant_prefix},
        mode=mode,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_joint_germline_cnv_segmentation(
    variant_paths: Sequence[str],
    pedigree: str,
    reference_path: str,
    output_path: str,
    *,
    model_call_intervals: str | None = None,
    min_sample_set_fraction_overlap: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine segmented germline copy-number (gCNV) VCFs of many samples into one cohort VCF.

    Run `gatk JointGermlineCNVSegmentation` (beta).

    Parameters:
    - ``variant_paths`` (required): Segmented gCNV VCFs of the samples.
    - ``pedigree`` (required): Pedigree (PED) file.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``output_path`` (required): New combined VCF.
    - ``model_call_intervals`` (optional): Interval list of the gCNV model calls.
    - ``min_sample_set_fraction_overlap`` (optional): Minimum overlap fraction to merge segments.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the combined VCF.
    """
    return _core.execute(
        _SPECS['gatk_joint_germline_cnv_segmentation'],
        {'variant_paths': variant_paths, 'pedigree': pedigree, 'reference_path': reference_path, 'output_path': output_path, 'model_call_intervals': model_call_intervals, 'min_sample_set_fraction_overlap': min_sample_set_fraction_overlap},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_print_read_counts(
    counts_path: str,
    *,
    output_prefix: str | None = None,
    output_file_list: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print per-sample count files for CNV determination from a read-count file.

    Run `gatk PrintReadCounts` (experimental). Output files are named by GATK from `output_prefix`.

    Parameters:
    - ``counts_path`` (required): Input read-count file.
    - ``output_prefix`` (optional): Prefix of the output paths `<prefix><sample>.counts.tsv`; include a folder that exists.
    - ``output_file_list`` (optional): New two-column list `sample<TAB>output path`.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_file_list` when requested.
    """
    return _core.execute(
        _SPECS['gatk_print_read_counts'],
        {'counts_path': counts_path, 'output_prefix': output_prefix, 'output_file_list': output_file_list},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_print_sv_evidence(
    evidence_file: str,
    output_path: str,
    *,
    sample_names: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge or subset GATK-SV evidence records into one output file.

    Run `gatk PrintSVEvidence` (experimental).

    Parameters:
    - ``evidence_file`` (required): SV evidence file (`.pe.txt`, `.sr.txt`, `.rd.txt`, `.sd.txt`, `.baf.txt`).
    - ``output_path`` (required): New evidence file.
    - ``sample_names`` (optional): Samples to keep.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the evidence file.
    """
    return _core.execute(
        _SPECS['gatk_print_sv_evidence'],
        {'evidence_file': evidence_file, 'output_path': output_path, 'sample_names': sample_names},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sv_annotate(
    variant_path: str,
    output_path: str,
    *,
    protein_coding_gtf: str | None = None,
    non_coding_bed: str | None = None,
    promoter_window_length: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate an SV VCF with gene overlaps and predicted variant consequences.

    Run `gatk SVAnnotate`.

    Parameters:
    - ``variant_path`` (required): Structural-variant VCF from the GATK-SV pipeline.
    - ``output_path`` (required): New annotated VCF.
    - ``protein_coding_gtf`` (optional): GTF of protein-coding genes.
    - ``non_coding_bed`` (optional): BED of non-coding elements.
    - ``promoter_window_length`` (optional): Promoter window upstream of a gene (bp).
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the annotated VCF.
    """
    return _core.execute(
        _SPECS['gatk_sv_annotate'],
        {'variant_path': variant_path, 'output_path': output_path, 'protein_coding_gtf': protein_coding_gtf, 'non_coding_bed': non_coding_bed, 'promoter_window_length': promoter_window_length},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sv_concordance(
    evaluation: str,
    truth: str,
    output_path: str,
    *,
    clustering_config: str | None = None,
    stratify_config: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate structural variants with their genotype concordance against a truth set.

    Run `gatk SVConcordance` (beta).

    Parameters:
    - ``evaluation`` (required): Evaluation SV VCF.
    - ``truth`` (required): Truth SV VCF.
    - ``output_path`` (required): New VCF with concordance annotations.
    - ``clustering_config`` (optional): Clustering configuration table.
    - ``stratify_config`` (optional): Stratification configuration table.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the annotated VCF.
    """
    return _core.execute(
        _SPECS['gatk_sv_concordance'],
        {'evaluation': evaluation, 'truth': truth, 'output_path': output_path, 'clustering_config': clustering_config, 'stratify_config': stratify_config},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_sv_stratify(
    variant_path: str,
    stratify_config: str,
    output_path: str,
    *,
    track_name: Sequence[str] | None = None,
    track_intervals: Sequence[str] | None = None,
    split_output: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate structural variants by SV type, size and reference-track overlap (stratification).

    Run `gatk SVStratify` (beta).

    Parameters:
    - ``variant_path`` (required): SV VCF.
    - ``stratify_config`` (required): Stratification configuration table.
    - ``output_path`` (required): New stratified VCF.
    - ``track_name`` (optional): Names of the reference tracks.
    - ``track_intervals`` (optional): Interval files of the reference tracks.
    - ``split_output`` (optional): Write one file per stratum.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the stratified VCF.
    """
    return _core.execute(
        _SPECS['gatk_sv_stratify'],
        {'variant_path': variant_path, 'stratify_config': stratify_config, 'output_path': output_path, 'track_name': track_name, 'track_intervals': track_intervals, 'split_output': split_output},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_site_depth_to_baf(
    site_depth_paths: Sequence[str],
    baf_sites_vcf: str,
    baf_evidence_output: str,
    *,
    sample_names: Sequence[str] | None = None,
    min_het_probability: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert site-depth evidence to B-allele-frequency (BAF) evidence for GATK-SV.

    Run `gatk SiteDepthtoBAF` (experimental).

    Parameters:
    - ``site_depth_paths`` (required): Site-depth evidence files (`.sd.txt`).
    - ``baf_sites_vcf`` (required): VCF of the candidate heterozygous sites.
    - ``baf_evidence_output`` (required): New B-allele-frequency evidence file.
    - ``sample_names`` (optional): Sample names of the inputs.
    - ``min_het_probability`` (optional): Minimum probability of a heterozygous site.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.baf_evidence_output` is the BAF file.
    """
    return _core.execute(
        _SPECS['gatk_site_depth_to_baf'],
        {'site_depth_paths': site_depth_paths, 'baf_sites_vcf': baf_sites_vcf, 'baf_evidence_output': baf_evidence_output, 'sample_names': sample_names, 'min_het_probability': min_het_probability},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_structural_variation_discovery_pipeline_spark(
    input_path: str,
    reference_path: str,
    aligner_index_image: str,
    kmers_to_ignore: str,
    contig_sam_file: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the structural-variation discovery workflow (evidence, local assembly, breakpoint calling) on one sample.

    Run `gatk StructuralVariationDiscoveryPipelineSpark` (beta, Spark local mode).

    Parameters:
    - ``input_path`` (required): Query-name-grouped or coordinate-sorted BAM of one sample.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict). A .2bit file also works for Spark.
    - ``aligner_index_image`` (required): BWA-MEM index image of the reference.
    - ``kmers_to_ignore`` (required): K-mer list from `gatk_find_bad_genomic_kmers_spark`.
    - ``contig_sam_file`` (required): New SAM of the assembled contig alignments.
    - ``output_path`` (required): New SV VCF.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the SV VCF.
    """
    return _core.execute(
        _SPECS['gatk_structural_variation_discovery_pipeline_spark'],
        {'input_path': input_path, 'reference_path': reference_path, 'aligner_index_image': aligner_index_image, 'kmers_to_ignore': kmers_to_ignore, 'contig_sam_file': contig_sam_file, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_calculate_genotype_posteriors(
    variant_path: str,
    output_path: str,
    *,
    pedigree: str | None = None,
    supporting_callsets: Sequence[str] | None = None,
    de_novo_prior: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate genotype posterior probabilities given family and/or known population genotypes.

    Run `gatk CalculateGenotypePosteriors`.

    Parameters:
    - ``variant_path`` (required): Input VCF with genotype likelihoods.
    - ``output_path`` (required): New VCF with posterior genotype probabilities.
    - ``pedigree`` (optional): Pedigree (PED) file for family priors.
    - ``supporting_callsets`` (optional): Population VCFs that provide allele-count priors.
    - ``de_novo_prior`` (optional): Prior probability of a de novo mutation.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the VCF.
    """
    return _core.execute(
        _SPECS['gatk_calculate_genotype_posteriors'],
        {'variant_path': variant_path, 'output_path': output_path, 'pedigree': pedigree, 'supporting_callsets': supporting_callsets, 'de_novo_prior': de_novo_prior},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_concordance(
    evaluation: str,
    truth: str,
    summary_path: str,
    *,
    true_positives_and_false_negatives: str | None = None,
    true_positives_and_false_positives: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate the concordance of a VCF against a validated truth VCF (true/false positives, false negatives).

    Run `gatk Concordance`; the summary table is parsed into `tables.summary_path`. Give the regions through the `intervals` option.

    Parameters:
    - ``evaluation`` (required): VCF to evaluate.
    - ``truth`` (required): Validated truth VCF.
    - ``summary_path`` (required): New summary table (sensitivity, precision).
    - ``true_positives_and_false_negatives`` (optional): New VCF of TP and FN records.
    - ``true_positives_and_false_positives`` (optional): New VCF of TP and FP records.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `tables.summary_path` (rows per variant type).
    """
    return _core.execute(
        _SPECS['gatk_concordance'],
        {'evaluation': evaluation, 'truth': truth, 'summary_path': summary_path, 'true_positives_and_false_negatives': true_positives_and_false_negatives, 'true_positives_and_false_positives': true_positives_and_false_positives},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_count_false_positives(
    variant_path: str,
    intervals: Sequence[str],
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the PASS variants of a call set (false positives when the sample has no true variants).

    Run `gatk CountFalsePositives` (beta).

    Parameters:
    - ``variant_path`` (required): VCF of calls on a truth sample (all calls are false positives).
    - ``intervals`` (required): Intervals or an interval file to restrict the count to.
    - ``output_path`` (required): New file for the count.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` holds the count.
    """
    return _core.execute(
        _SPECS['gatk_count_false_positives'],
        {'variant_path': variant_path, 'intervals': intervals, 'output_path': output_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_evaluate_info_field_concordance(
    evaluation: str,
    truth: str,
    eval_info_key: str,
    truth_info_key: str,
    summary_path: str,
    *,
    epsilon: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate the concordance of an INFO field of a VCF against the same or another field of a truth VCF.

    Run `gatk EvaluateInfoFieldConcordance` (beta).

    Parameters:
    - ``evaluation`` (required): VCF to evaluate.
    - ``truth`` (required): Truth VCF.
    - ``eval_info_key`` (required): INFO key in the evaluation VCF.
    - ``truth_info_key`` (required): INFO key in the truth VCF.
    - ``summary_path`` (required): New summary file.
    - ``epsilon`` (optional): Tolerance when comparing values.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.summary_path` is the summary.
    """
    return _core.execute(
        _SPECS['gatk_evaluate_info_field_concordance'],
        {'evaluation': evaluation, 'truth': truth, 'eval_info_key': eval_info_key, 'truth_info_key': truth_info_key, 'summary_path': summary_path, 'epsilon': epsilon},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_reference_block_concordance(
    evaluation: str,
    truth: str,
    eval_block_histogram: str,
    truth_block_histogram: str,
    confidence_concordance_histogram: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate the reference-block concordance of a GVCF against a truth GVCF.

    Run `gatk ReferenceBlockConcordance`.

    Parameters:
    - ``evaluation`` (required): GVCF to evaluate.
    - ``truth`` (required): Truth GVCF.
    - ``eval_block_histogram`` (required): New histogram of the evaluation reference blocks.
    - ``truth_block_histogram`` (required): New histogram of the truth reference blocks.
    - ``confidence_concordance_histogram`` (required): New histogram of the confidence concordance.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs` lists the three histograms.
    """
    return _core.execute(
        _SPECS['gatk_reference_block_concordance'],
        {'evaluation': evaluation, 'truth': truth, 'eval_block_histogram': eval_block_histogram, 'truth_block_histogram': truth_block_histogram, 'confidence_concordance_histogram': confidence_concordance_histogram},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_filter_funcotations(
    variant_path: str,
    output_path: str,
    ref_version: str,
    allele_frequency_data_source: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter variants on clinically significant Funcotator annotations (ClinVar, allele frequency, ...).

    Run `gatk FilterFuncotations` (experimental).

    Parameters:
    - ``variant_path`` (required): Funcotator-annotated VCF.
    - ``output_path`` (required): New filtered VCF.
    - ``ref_version`` (required): Reference build of the annotation.
    - ``allele_frequency_data_source`` (required): Population database that holds the allele frequencies.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the filtered VCF.
    """
    return _core.execute(
        _SPECS['gatk_filter_funcotations'],
        {'variant_path': variant_path, 'output_path': output_path, 'ref_version': ref_version, 'allele_frequency_data_source': allele_frequency_data_source},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_find_mendelian_violations(
    input_path: str,
    pedigree_path: str,
    output_path: str,
    *,
    min_gq: int | None = None,
    min_dp: int | None = None,
    skip_chroms: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find Mendelian violations of all types in a VCF with trios.

    Run Picard `gatk FindMendelianViolations`; the metrics file is parsed into `tables.output_path`.

    Parameters:
    - ``input_path`` (required): VCF with the trio members.
    - ``pedigree_path`` (required): PED file of the trios (no genotype columns).
    - ``output_path`` (required): New metrics file.
    - ``min_gq`` (optional): Minimum genotype quality.
    - ``min_dp`` (optional): Minimum depth.
    - ``skip_chroms`` (optional): Contigs to skip.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `tables.output_path` (Picard metrics rows).
    """
    return _core.execute(
        _SPECS['gatk_find_mendelian_violations'],
        {'input_path': input_path, 'pedigree_path': pedigree_path, 'output_path': output_path, 'min_gq': min_gq, 'min_dp': min_dp, 'skip_chroms': skip_chroms},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_funcotate(
    input_path: str,
    reference_path: str,
    data_sources_path: str,
    ref_version: str,
    output_path: str,
    output_file_format: str,
    *,
    input_kind: str = 'variants',
    transcript_selection_mode: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate variants or copy-number segments functionally (genes, consequences, ClinVar, gnomAD ...) with Funcotator.

    Run `gatk Funcotator` (`input_kind=variants`) or `gatk FuncotateSegments` (`input_kind=segments`, beta; output format not well defined upstream).

    ``input_kind``: `variants` = Funcotator, `segments` = FuncotateSegments.

    Parameters:
    - ``input_path`` (required): VCF (`variants`) or segment file (`segments`) to annotate.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``data_sources_path`` (required): Funcotator data-sources folder (see FuncotatorDataSourceDownloader).
    - ``ref_version`` (required): Reference build of the input, for example hg19 or hg38.
    - ``output_path`` (required): New annotated file.
    - ``output_file_format`` (required): Output format.
    - ``transcript_selection_mode`` (optional): How the transcript of a variant is chosen.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields and `mode`; `outputs.output_path` is the annotated file.
    """
    return _core.execute(
        _SPECS['gatk_funcotate'],
        {'input_path': input_path, 'reference_path': reference_path, 'data_sources_path': data_sources_path, 'ref_version': ref_version, 'output_path': output_path, 'output_file_format': output_file_format, 'transcript_selection_mode': transcript_selection_mode},
        mode=input_kind,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


def gatk_validate_basic_somatic_short_mutations(
    variant_path: str,
    input_path: str,
    reference_path: str,
    discovery_sample_name: str,
    val_case_sample_name: str,
    val_control_sample_name: str,
    output_path: str,
    *,
    summary_path: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check somatic short mutations against tumor-normal BAMs of the same samples sequenced differently (validation power).

    Run `gatk ValidateBasicSomaticShortMutations` (experimental).

    Parameters:
    - ``variant_path`` (required): VCF of the discovery calls to validate.
    - ``input_path`` (required): BAM/CRAM of the validation tumor and normal samples.
    - ``reference_path`` (required): Reference FASTA (with .fai and .dict).
    - ``discovery_sample_name`` (required): Sample name in the discovery VCF.
    - ``val_case_sample_name`` (required): Validation case (tumor) sample name in the BAM.
    - ``val_control_sample_name`` (required): Validation control (normal) sample name in the BAM.
    - ``output_path`` (required): New validation table.
    - ``summary_path`` (optional): New summary file.
    - ``options``: further GATK options of the tool as ``{long option name: value}``; only options listed by `gatk <Tool> --help` are accepted.
    - ``gatk_path`` / ``java_max_memory`` / ``verbose``: launcher path (default: DNAKIT_GATK or PATH), Java heap such as ``'4g'``, console progress.

    Output: The common result fields; `outputs.output_path` is the validation table.
    """
    return _core.execute(
        _SPECS['gatk_validate_basic_somatic_short_mutations'],
        {'variant_path': variant_path, 'input_path': input_path, 'reference_path': reference_path, 'discovery_sample_name': discovery_sample_name, 'val_case_sample_name': val_case_sample_name, 'val_control_sample_name': val_control_sample_name, 'output_path': output_path, 'summary_path': summary_path},
        mode=None,
        options=options,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )


__all__ = [
    "gatk_reblock_gvcf",
    "gatk_calibrate_dragstr_model",
    "gatk_combine_gvcfs",
    "gatk_genomicsdb_import",
    "gatk_genotype_gvcfs",
    "gatk_haplotype_based_variant_recaller",
    "gatk_learn_read_orientation_model",
    "gatk_mutect2",
    "gatk_reads_pipeline_spark",
    "gatk_collect_sv_evidence",
    "gatk_condense_depth_evidence",
    "gatk_convert_counts_to_depth_file",
    "gatk_sv_cluster",
    "gatk_joint_germline_cnv_segmentation",
    "gatk_print_read_counts",
    "gatk_print_sv_evidence",
    "gatk_sv_annotate",
    "gatk_sv_concordance",
    "gatk_sv_stratify",
    "gatk_site_depth_to_baf",
    "gatk_structural_variation_discovery_pipeline_spark",
    "gatk_calculate_genotype_posteriors",
    "gatk_concordance",
    "gatk_count_false_positives",
    "gatk_evaluate_info_field_concordance",
    "gatk_reference_block_concordance",
    "gatk_filter_funcotations",
    "gatk_find_mendelian_violations",
    "gatk_funcotate",
    "gatk_validate_basic_somatic_short_mutations",
]

_register_public_verbose(globals(), tuple(__all__))
