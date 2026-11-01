"""DeepVariant command-line facades: end-to-end runners, pipeline stages, training and the helper tools.

Each public function runs one fixed program of DeepVariant (https://github.com/google/deepvariant; 1.10.0 is the
version whose flags were used) with an explicit, validated argument list: the runners ``run_deepvariant``,
``run_deeptrio``, ``run_pangenome_aware_deepvariant``, ``run_deepsomatic`` and ``run_oracle_inference``; the stages
``make_examples`` (and its trio, multi-sample, somatic and pangenome-aware variants), ``call_variants`` and
``postprocess_variants``; ``train``, ``show_examples``, ``vcf_stats_report``, ``labeled_examples_to_vcf`` and the helper
tools ``compare_labelers``, ``preprocess_truth`` and ``print_f1``. No shell is used.

DeepVariant is distributed as container images (Docker or Apptainer/Singularity), not as a pip package. DNAKit therefore
does not install it: the user either runs DNAKit inside the DeepVariant image, or makes the programs executable on the
host (a wrapper script that calls the container is enough) and puts them on PATH or passes ``deepvariant_dir``, the folder
that holds them (``/opt/deepvariant/bin`` in the image). All paths passed to a function must be visible to the program.
``compare_labelers``, ``preprocess_truth`` and ``print_f1`` are source-only tools of the DeepVariant repository; they
must be built (Bazel) or wrapped by the user. A missing program raises ``BackendUnavailableError`` with an install hint.
Results are JSON-serialisable dictionaries with the status, backend, version, command, outputs and a bounded log tail.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p12_cli as _cli
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_TOOL = _cli.Tool(
    name="DeepVariant",
    executable="run_deepvariant",
    hint=(
        "DeepVariant is distributed as a Docker/Apptainer image (https://github.com/google/deepvariant). Run DNAKit "
        "inside the image or make its programs executable on the host (for example wrapper scripts that call the "
        "container), put them on PATH or pass deepvariant_dir."
    ),
    tested_version="1.10.0",
    citation="https://doi.org/10.1038/nbt.4235",
    version_arguments=("--version",),
    version_pattern=r"version\s+(\d[\w.]*)",
)


def _names(text: str) -> frozenset[str]:
    """The white-space separated words of ``text`` as a set (the flag whitelists below)."""
    return frozenset(text.split())


_MODEL_TYPES = ("WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA", "MASSEQ", "RNASEQ")
_TRIO_MODEL_TYPES = ("WGS", "WES", "PACBIO", "ONT")
_PANGENOME_MODEL_TYPES = ("WGS", "WES")
_SOMATIC_MODEL_TYPES = (
    "WGS",
    "WES",
    "PACBIO",
    "ONT",
    "FFPE_WGS",
    "FFPE_WES",
    "FFPE_WGS_TUMOR_ONLY",
    "FFPE_WES_TUMOR_ONLY",
    "WGS_TUMOR_ONLY",
    "WES_TUMOR_ONLY",
    "PACBIO_TUMOR_ONLY",
    "ONT_TUMOR_ONLY",
)
_REGION_WORD = re.compile(r"^[^\s'\"\\\x00]+$")
_CONFIG_NAME = re.compile(r"^[A-Za-z0-9_]+$")
_ORACLE_MODEL_TYPES = ("WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA", "MASSEQ")

# Flags that make_examples (all variants) accepts, from make_examples_options.py of DeepVariant 1.10.0.
_MAKE_EXAMPLES_FLAGS = _names(
    """ref use_ref_for_cram examples checkpoint checkpoint_json candidates mode regions exclude_regions variant_caller
    gvcf gvcf_gq_binsize p_error include_med_dp confident_regions truth_variants task partition_size
    max_reads_per_partition max_reads_for_dynamic_bases_per_region multi_allelic_mode realign_reads
    max_read_length_to_realign trim_reads_for_pileup write_run_info alt_aligned_pileup types_to_alt_align
    hts_block_size min_base_quality min_mapping_quality vsc_min_count_snps vsc_min_count_indels
    vsc_min_fraction_snps vsc_min_fraction_indels vsc_min_fraction_multiplier
    vsc_max_fraction_indels_for_non_target_sample vsc_max_fraction_snps_for_non_target_sample downsample_classes
    training_random_emit_ref_sites pileup_image_width labeler_algorithm customized_classes_labeler_classes_list
    customized_classes_labeler_info_field_name logging_every_n_candidates keep_duplicates
    keep_supplementary_alignments keep_secondary_alignments parse_sam_aux_fields aux_fields_to_keep
    use_original_quality_scores select_variant_types sequencing_type sort_by_haplotypes reverse_haplotypes
    hp_tag_for_assembly_polishing sort_by_alt_allele_support add_hp_channel channels channel_list
    add_supporting_other_alt_color population_vcfs use_allele_frequency mean_coverage_per_sample
    sample_mean_coverage_on_calling_regions runtime_by_region track_ref_reads normalize_reads
    keep_legacy_allele_counter_behavior enable_strict_insertion_filter enable_methylation_calling
    methylation_calling_threshold enable_methylation_aware_phasing exclude_contigs_for_methylation_phasing
    phase_reads min_alleles_to_phase phase_max_candidates enable_joint_realignment output_local_read_phasing
    output_phasing_error_stats discard_non_dna_regions output_sitelist denovo_regions deterministic_serialization
    write_small_model_examples skip_pileup_image_generation call_small_model_examples
    small_model_snp_gq_threshold small_model_indel_gq_threshold small_model_call_multiallelics
    small_model_emit_all_candidates small_model_inference_batch_size small_model_vaf_context_window_size
    stream_examples shm_prefix shm_buffer_size exclude_variants_vcf_filename exclude_variants_af_threshold
    output_phase_info haploid_contigs par_regions_bed create_complex_alleles use_non_uniform_downsampling
    non_uniform_downsampling_threshold output_debug_info assign_phase_from_normal
    vsc_min_indel_fraction_for_small_indels vsc_min_indel_fraction_for_large_indels vsc_small_indel_threshold
    use_rejected_alleles filter_low_vaf_candidates low_vaf_threshold low_vaf_max_base_quality
    low_vaf_max_mapping_quality"""
)
_MAKE_EXAMPLES_PROGRAMS: dict[str, tuple[str, frozenset[str]]] = {
    "make_examples": (
        "make_examples",
        _names(
            "reads sample_name pileup_image_height downsample_fraction proposed_variants candidate_positions "
            "trained_small_model_path"
        ),
    ),
    "make_examples_somatic": (
        "make_examples_somatic",
        _names(
            "reads_tumor reads_normal sample_name_tumor sample_name_normal downsample_fraction_tumor "
            "downsample_fraction_normal pileup_image_height_tumor pileup_image_height_normal candidate_positions "
            "trained_small_model_path"
        ),
    ),
    "multisample_make_examples": (
        "multisample_make_examples",
        _names("reads sample_names downsample_fractions pileup_image_heights"),
    ),
    "deeptrio_make_examples": (
        "deeptrio/make_examples",
        _names(
            "sample_name_to_train reads reads_parent1 reads_parent2 downsample_fraction_child "
            "downsample_fraction_parents sample_name sample_name_parent1 sample_name_parent2 "
            "pileup_image_height_parent pileup_image_height_child proposed_variants_child proposed_variants_parent1 "
            "proposed_variants_parent2 small_model_path_child small_model_path_parent candidate_positions "
            "skip_parent_calling"
        ),
    ),
    "make_examples_pangenome_aware_dv": (
        "make_examples_pangenome_aware_dv",
        _names(
            "reads pangenome sample_name_reads sample_name_pangenome downsample_fraction_reads "
            "pileup_image_height_reads pileup_image_height_pangenome keep_only_window_spanning_haplotypes "
            "variant_types_to_blank ref_name_pangenome ref_chrom_prefix use_loaded_gbz_shared_memory "
            "gbz_shared_memory_name trained_small_model_path"
        ),
    ),
}
_CALL_VARIANTS_FLAGS = _names(
    "batch_size max_batches num_readers include_debug_info activation_layers debugging_true_label_mode "
    "execution_hardware config_string kmp_blocktime writer_threads limit num_input_shards shm_prefix "
    "stream_examples allow_empty_examples"
)
_POSTPROCESS_FLAGS = _names(
    "small_model_cvo_records qual_filter cnn_homref_call_min_gq multi_allelic_qual_filter "
    "nonvariant_site_tfrecord_path phased_reads_input_path checkpoint_json phased_reads_switches_output_path "
    "phased_reads_corrected_output_path gvcf_outfile group_variants sample_name use_multiallelic_model "
    "multiallelic_mode debug_output_all_candidates only_keep_pass haploid_contigs cpus num_partitions "
    "par_regions_bed regions process_somatic pon_filtering resolve_call_variants_outputs_by_model "
    "small_model_gq_threshold"
)
_TRAIN_CONFIG_KEYS = _names(
    "batch_size num_epochs num_validation_examples warmup_steps limit steps_per_iter shuffle_buffer_elements "
    "init_checkpoint best_checkpoint_metric optimizer momentum use_ema ema_momentum optimizer_weight_decay "
    "weight_decay early_stopping_patience learning_rate learning_rate_num_epochs_per_decay "
    "learning_rate_decay_rate label_smoothing backbone_dropout_rate use_mixed_precision beta_1 beta_2 "
    "adaptive_epsilon include_snp_indel_metrics"
)


# ---------------------------------------------------------------------------
# Validation and argument helpers
# ---------------------------------------------------------------------------


def _explicit(deepvariant_dir: str | None, name: str) -> str | None:
    if deepvariant_dir is None:
        return None
    folder = Path(deepvariant_dir).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError(f"deepvariant_dir must be an existing folder: {folder}")
    return str(folder / name)


def _scalar(value: object, label: str) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, str) and value and not any(char in value for char in "\x00\n\r"):
        return value
    raise ConfigurationError(
        f"{label} must be a boolean, a number or non-empty text without line breaks."
    )


def _flag_arguments(
    values: Mapping[str, object] | None, allowed: frozenset[str], label: str
) -> list[str]:
    """``--name=value`` arguments for a mapping whose keys are all in the whitelist ``allowed``."""
    if values is None:
        return []
    if not isinstance(values, Mapping):
        raise ConfigurationError(f"{label} must be a mapping from flag name to value.")
    arguments: list[str] = []
    for key, value in values.items():
        if key not in allowed:
            raise ConfigurationError(
                f"{label}: flag {key!r} is not allowed. Allowed flags: {sorted(allowed)}."
            )
        arguments.append(f"--{key}={_scalar(value, f'{label}[{key}]')}")
    return arguments


def _extra_string(
    values: Mapping[str, object] | None, allowed: frozenset[str], label: str
) -> str | None:
    """The comma-separated ``name=value`` text of the ``*_extra_args`` flags of the runners."""
    if not values:
        return None
    parts: list[str] = []
    for key, value in values.items():
        if key not in allowed:
            raise ConfigurationError(
                f"{label}: flag {key!r} is not allowed. Allowed flags: {sorted(allowed)}."
            )
        text = _scalar(value, f"{label}[{key}]")
        if "," in text:
            raise ConfigurationError(f"{label}[{key}] must not contain a comma.")
        parts.append(f"{key}={text}")
    return ",".join(parts)


def _regions(regions: Sequence[str] | None) -> str | None:
    if regions is None:
        return None
    if isinstance(regions, str) or not regions:
        raise ConfigurationError(
            "regions must be a non-empty list of region literals or BED/BEDPE paths."
        )
    for item in regions:
        if not isinstance(item, str) or not _REGION_WORD.match(item) or item.startswith("-"):
            raise ConfigurationError(
                "each region must be a text without spaces, quotes or a leading '-'."
            )
    return " ".join(regions)


def _optional(arguments: list[str], flag: str, value: object | None) -> None:
    if value is not None:
        arguments.append(f"--{flag}={_scalar(value, flag)}")


def _boolean(arguments: list[str], flag: str, value: bool | None) -> None:
    if value is not None:
        arguments.append(f"--{flag}={'true' if _cli.flag(value, flag) else 'false'}")


def _folder_to_use(value: str, label: str) -> Path:
    """An existing folder, or a new folder whose parent exists (the program creates it)."""
    path = Path(_cli.free_text(value, label)).expanduser().resolve()
    if path.exists() and not path.is_dir():
        raise ConfigurationError(f"{label} must be a folder: {path}")
    if not path.exists() and not path.parent.is_dir():
        raise ConfigurationError(f"The parent folder of {label} does not exist: {path.parent}")
    return path


def _sharded(pattern: str) -> list[Path]:
    """Files named by a DeepVariant sharded path (``name@4.gz`` means ``name-0000?-of-00004.gz``) or a plain path."""
    match = re.search(r"@(\d+)", pattern)
    if match is None:
        path = Path(pattern)
        return [path] if path.is_file() else []
    count = int(match.group(1))
    path = Path(pattern.replace(match.group(0), f"-?????-of-{count:05d}"))
    return sorted(path.parent.glob(path.name))


def _version(deepvariant_dir: str | None) -> str | None:
    try:
        executable = _cli.program(
            "run_deepvariant",
            explicit=_explicit(deepvariant_dir, "run_deepvariant"),
            hint=_TOOL.hint,
        )
    except Exception:
        return None
    return _cli._version_of(executable, ("--version",), r"version\s+(\d[\w.]*)")


def _run(
    name: str,
    arguments: list[str],
    *,
    cwd: Path,
    expected: Sequence[Path],
    deepvariant_dir: str | None,
    operation: str,
    timeout_seconds: float,
    stdout_to: Path | None = None,
    **extra: Any,
) -> dict[str, Any]:
    record = _cli.invoke(
        _TOOL,
        arguments,
        cwd=cwd,
        expected=list(expected),
        stdout_to=stdout_to,
        explicit=_explicit(deepvariant_dir, name),
        executable_name=name.split("/")[-1] if deepvariant_dir is None else name,
        timeout_seconds=timeout_seconds,
    )
    if record.get("version") is None:
        record["version"] = _version(deepvariant_dir)
    resolved = {key: (value() if callable(value) else value) for key, value in extra.items()}
    return _cli.finish(record, operation, **resolved)


def _stats_report_path(vcf: Path) -> Path:
    name = vcf.name
    for suffix in (".gz", ".vcf"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return vcf.with_name(name + ".visual_report.html")


# ---------------------------------------------------------------------------
# End-to-end runners
# ---------------------------------------------------------------------------


def deepvariant_run(
    model_type: str,
    reference_path: str,
    reads_path: str,
    output_vcf: str,
    *,
    output_gvcf: str | None = None,
    num_shards: int = 1,
    regions: Sequence[str] | None = None,
    sample_name: str | None = None,
    haploid_contigs: str | None = None,
    par_regions_bed: str | None = None,
    intermediate_results_dir: str | None = None,
    logging_dir: str | None = None,
    runtime_report: bool = False,
    customized_model: str | None = None,
    customized_model_json: str | None = None,
    disable_small_model: bool = False,
    customized_small_model: str | None = None,
    phase_vcf: bool | None = None,
    vcf_stats_report: bool = False,
    report_title: str | None = None,
    postprocess_cpus: int | None = None,
    make_examples_extra_args: Mapping[str, str | int | float | bool] | None = None,
    call_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    postprocess_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    dry_run: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call germline variants from aligned reads with DeepVariant's ``run_deepvariant`` (make_examples, call_variants, postprocess_variants).

    ``model_type`` is WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA, MASSEQ or RNASEQ. ``reference_path`` is the
    FASTA with its ``.fai`` index and ``reads_path`` an aligned, sorted, indexed BAM or CRAM. ``output_vcf`` (and the
    optional ``output_gvcf``) must not exist. ``phase_vcf=True`` writes read-backed haplotype phasing into the VCF
    (``--phase_vcf``, long-read PACBIO and ONT_R104 models only; None keeps the program default, which is on for those
    models); ``vcf_stats_report=True`` also writes the HTML statistics report ``<output_vcf>.visual_report.html``
    (``--vcf_stats_report``, titled ``report_title``). Other options: ``num_shards`` (parallel make_examples shards),
    ``regions`` (region literals or BED paths), ``sample_name``, ``haploid_contigs`` and ``par_regions_bed`` for
    sex chromosomes, ``intermediate_results_dir``, ``logging_dir`` with ``runtime_report``, ``customized_model`` (a
    checkpoint) with ``customized_model_json``, ``disable_small_model``, ``customized_small_model``, ``postprocess_cpus``
    and mappings ``make_examples_extra_args``, ``call_variants_extra_args`` and ``postprocess_variants_extra_args``
    holding the flags of the corresponding stage (see ``deepvariant_make_examples`` for the names). ``dry_run`` only
    prints the commands the program would run. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report``
    (None when not produced).
    """
    kind = _cli.choice(model_type, _MODEL_TYPES, "model_type")
    reference = _cli.existing_file(reference_path, "reference_path")
    reads = _cli.existing_file(reads_path, "reads_path")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    arguments = [
        f"--model_type={kind}",
        f"--ref={reference}",
        f"--reads={reads}",
        f"--output_vcf={vcf}",
        f"--num_shards={_cli.integer(num_shards, 'num_shards', minimum=1, maximum=1024)}",
    ]
    expected: list[Path] = []
    gvcf = None
    if output_gvcf is not None:
        gvcf = _cli.new_file(output_gvcf, "output_gvcf")
        arguments.append(f"--output_gvcf={gvcf}")
    region_text = _regions(regions)
    _optional(arguments, "regions", region_text)
    _optional(
        arguments,
        "sample_name",
        _cli.text(sample_name, "sample_name") if sample_name is not None else None,
    )
    _optional(
        arguments,
        "haploid_contigs",
        _cli.text(haploid_contigs, "haploid_contigs") if haploid_contigs else None,
    )
    if par_regions_bed is not None:
        arguments.append(
            f"--par_regions_bed={_cli.existing_file(par_regions_bed, 'par_regions_bed')}"
        )
    if intermediate_results_dir is not None:
        arguments.append(
            f"--intermediate_results_dir={_folder_to_use(intermediate_results_dir, 'intermediate_results_dir')}"
        )
    if logging_dir is not None:
        arguments.append(f"--logging_dir={_folder_to_use(logging_dir, 'logging_dir')}")
        _boolean(arguments, "runtime_report", runtime_report)
    elif runtime_report:
        raise ConfigurationError("runtime_report needs logging_dir.")
    if customized_model is not None:
        arguments.append(
            f"--customized_model={_cli.free_text(customized_model, 'customized_model')}"
        )
    if customized_model_json is not None:
        arguments.append(
            f"--customized_model_json={_cli.existing_file(customized_model_json, 'customized_model_json')}"
        )
    _boolean(arguments, "disable_small_model", disable_small_model)
    if customized_small_model is not None:
        arguments.append(
            f"--customized_small_model={_cli.free_text(customized_small_model, 'customized_small_model')}"
        )
    _boolean(arguments, "phase_vcf", phase_vcf)
    _boolean(arguments, "vcf_stats_report", vcf_stats_report)
    _optional(
        arguments,
        "report_title",
        _cli.text(report_title, "report_title") if report_title is not None else None,
    )
    if postprocess_cpus is not None:
        arguments.append(
            f"--postprocess_cpus={_cli.integer(postprocess_cpus, 'postprocess_cpus', minimum=0, maximum=1024)}"
        )
    for flag_name, mapping, allowed in (
        (
            "make_examples_extra_args",
            make_examples_extra_args,
            _MAKE_EXAMPLES_FLAGS | _MAKE_EXAMPLES_PROGRAMS["make_examples"][1],
        ),
        ("call_variants_extra_args", call_variants_extra_args, _CALL_VARIANTS_FLAGS),
        ("postprocess_variants_extra_args", postprocess_variants_extra_args, _POSTPROCESS_FLAGS),
    ):
        text = _extra_string(mapping, allowed, flag_name)
        _optional(arguments, flag_name, text)
    _boolean(arguments, "dry_run", dry_run)
    report = _stats_report_path(vcf) if vcf_stats_report else None
    if not dry_run:
        expected = [vcf] + ([gvcf] if gvcf else []) + ([report] if report else [])
    return _run(
        "run_deepvariant",
        arguments,
        cwd=vcf.parent,
        expected=expected,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_run",
        timeout_seconds=timeout_seconds,
        output_vcf=lambda: str(vcf) if vcf.exists() else None,
        output_gvcf=lambda: str(gvcf) if gvcf and gvcf.exists() else None,
        stats_report=lambda: str(report) if report and report.exists() else None,
    )


def deepvariant_run_deeptrio(
    model_type: str,
    reference_path: str,
    reads_child: str,
    output_vcf_child: str,
    sample_name_child: str,
    *,
    reads_parent1: str | None = None,
    output_vcf_parent1: str | None = None,
    sample_name_parent1: str | None = None,
    reads_parent2: str | None = None,
    output_vcf_parent2: str | None = None,
    sample_name_parent2: str | None = None,
    output_gvcf_child: str | None = None,
    output_gvcf_parent1: str | None = None,
    output_gvcf_parent2: str | None = None,
    num_shards: int = 1,
    regions: Sequence[str] | None = None,
    intermediate_results_dir: str | None = None,
    logging_dir: str | None = None,
    customized_model_child: str | None = None,
    customized_model_parent: str | None = None,
    disable_small_model: bool = False,
    use_candidate_partition: bool | None = None,
    vcf_stats_report: bool = False,
    make_examples_extra_args: Mapping[str, str | int | float | bool] | None = None,
    call_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    postprocess_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    dry_run: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call variants jointly in a trio or a duo from aligned reads with ``run_deeptrio`` (DeepTrio).

    ``model_type`` is WGS, WES, PACBIO or ONT. The child's reads, VCF and sample name are required; give the reads, VCF
    and sample name of one parent for a duo, of both parents for a trio (each parent needs all three). ``*_gvcf``
    paths are optional; every output must not exist. ``customized_model_child`` and ``customized_model_parent`` replace the
    default checkpoints, ``use_candidate_partition`` toggles the candidate-sweep partitioning (None keeps the program
    default), and the mappings ``make_examples_extra_args``, ``call_variants_extra_args`` and
    ``postprocess_variants_extra_args`` pass further stage flags (names as in ``deepvariant_make_examples`` for the
    DeepTrio variant). Other options as in ``deepvariant_run``. Result keys: ``output_vcfs`` (a mapping from ``child``,
    ``parent1``, ``parent2`` to the written VCF) and ``output_gvcfs``.
    """
    kind = _cli.choice(model_type, _TRIO_MODEL_TYPES, "model_type")
    reference = _cli.existing_file(reference_path, "reference_path")
    members: dict[str, tuple[Path, Path, str, Path | None]] = {}
    for role, reads, vcf, name, gvcf in (
        ("child", reads_child, output_vcf_child, sample_name_child, output_gvcf_child),
        ("parent1", reads_parent1, output_vcf_parent1, sample_name_parent1, output_gvcf_parent1),
        ("parent2", reads_parent2, output_vcf_parent2, sample_name_parent2, output_gvcf_parent2),
    ):
        if reads is None and vcf is None and name is None and gvcf is None:
            continue
        if reads is None or vcf is None or name is None:
            raise ConfigurationError(
                f"{role}: reads, output VCF and sample name must be given together."
            )
        members[role] = (
            _cli.existing_file(reads, f"reads_{role}"),
            _cli.new_file(vcf, f"output_vcf_{role}"),
            _cli.text(name, f"sample_name_{role}"),
            _cli.new_file(gvcf, f"output_gvcf_{role}") if gvcf is not None else None,
        )
    if "parent1" not in members and "parent2" not in members:
        raise ConfigurationError("Give the reads, VCF and sample name of at least one parent.")
    arguments = [
        f"--model_type={kind}",
        f"--ref={reference}",
        f"--num_shards={_cli.integer(num_shards, 'num_shards', minimum=1, maximum=1024)}",
    ]
    expected: list[Path] = []
    for role, (member_reads, member_vcf, member_name, member_gvcf) in members.items():
        arguments += [
            f"--reads_{role}={member_reads}",
            f"--output_vcf_{role}={member_vcf}",
            f"--sample_name_{role}={member_name}",
        ]
        expected.append(member_vcf)
        if member_gvcf is not None:
            arguments.append(f"--output_gvcf_{role}={member_gvcf}")
            expected.append(member_gvcf)
    _optional(arguments, "regions", _regions(regions))
    if intermediate_results_dir is not None:
        arguments.append(
            f"--intermediate_results_dir={_folder_to_use(intermediate_results_dir, 'intermediate_results_dir')}"
        )
    if logging_dir is not None:
        arguments.append(f"--logging_dir={_folder_to_use(logging_dir, 'logging_dir')}")
    for flag_name, value in (
        ("customized_model_child", customized_model_child),
        ("customized_model_parent", customized_model_parent),
    ):
        if value is not None:
            arguments.append(f"--{flag_name}={_cli.free_text(value, flag_name)}")
    _boolean(arguments, "disable_small_model", disable_small_model)
    _boolean(arguments, "use_candidate_partition", use_candidate_partition)
    _boolean(arguments, "vcf_stats_report", vcf_stats_report)
    trio_flags = _MAKE_EXAMPLES_FLAGS | _MAKE_EXAMPLES_PROGRAMS["deeptrio_make_examples"][1]
    for flag_name, mapping, allowed in (
        ("make_examples_extra_args", make_examples_extra_args, trio_flags),
        ("call_variants_extra_args", call_variants_extra_args, _CALL_VARIANTS_FLAGS),
        ("postprocess_variants_extra_args", postprocess_variants_extra_args, _POSTPROCESS_FLAGS),
    ):
        _optional(arguments, flag_name, _extra_string(mapping, allowed, flag_name))
    _boolean(arguments, "dry_run", dry_run)
    cwd = members["child"][1].parent
    return _run(
        "run_deeptrio",
        arguments,
        cwd=cwd,
        expected=[] if dry_run else expected,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_run_deeptrio",
        timeout_seconds=timeout_seconds,
        output_vcfs=lambda: {
            role: str(item[1]) for role, item in members.items() if item[1].exists()
        },
        output_gvcfs=lambda: {
            role: str(item[3])
            for role, item in members.items()
            if item[3] is not None and item[3].exists()
        },
    )


def deepvariant_run_pangenome_aware(
    model_type: str,
    reference_path: str,
    reads_path: str,
    pangenome_path: str,
    output_vcf: str,
    *,
    output_gvcf: str | None = None,
    num_shards: int = 1,
    regions: Sequence[str] | None = None,
    sample_name_reads: str | None = None,
    sample_name_pangenome: str | None = None,
    ref_name_pangenome: str | None = None,
    intermediate_results_dir: str | None = None,
    logging_dir: str | None = None,
    customized_model: str | None = None,
    disable_small_model: bool = False,
    customized_small_model: str | None = None,
    vcf_stats_report: bool = False,
    report_title: str | None = None,
    postprocess_cpus: int | None = None,
    gbz_shared_memory_name: str | None = None,
    gbz_shared_memory_size_gb: int | None = None,
    make_examples_extra_args: Mapping[str, str | int | float | bool] | None = None,
    call_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    postprocess_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    dry_run: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call germline variants with pangenome-aware DeepVariant (``run_pangenome_aware_deepvariant``).

    The pileup images also contain haplotypes of a pangenome graph. ``model_type`` is WGS or WES, ``pangenome_path`` the
    pangenome graph in GBZ format and ``reads_path`` aligned, sorted, indexed short reads. ``sample_name_reads``,
    ``sample_name_pangenome`` and ``ref_name_pangenome`` name the read sample, the pangenome sample and the reference
    path inside the graph; ``gbz_shared_memory_name`` and ``gbz_shared_memory_size_gb`` use a graph already loaded into
    shared memory. Other options as in ``deepvariant_run``. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report``.
    """
    kind = _cli.choice(model_type, _PANGENOME_MODEL_TYPES, "model_type")
    reference = _cli.existing_file(reference_path, "reference_path")
    reads = _cli.existing_file(reads_path, "reads_path")
    graph = _cli.existing_file(pangenome_path, "pangenome_path")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    arguments = [
        f"--model_type={kind}",
        f"--ref={reference}",
        f"--reads={reads}",
        f"--pangenome={graph}",
        f"--output_vcf={vcf}",
        f"--num_shards={_cli.integer(num_shards, 'num_shards', minimum=1, maximum=1024)}",
    ]
    gvcf = None
    if output_gvcf is not None:
        gvcf = _cli.new_file(output_gvcf, "output_gvcf")
        arguments.append(f"--output_gvcf={gvcf}")
    _optional(arguments, "regions", _regions(regions))
    for flag_name, value in (
        ("sample_name_reads", sample_name_reads),
        ("sample_name_pangenome", sample_name_pangenome),
        ("ref_name_pangenome", ref_name_pangenome),
        ("gbz_shared_memory_name", gbz_shared_memory_name),
        ("report_title", report_title),
    ):
        _optional(arguments, flag_name, _cli.text(value, flag_name) if value is not None else None)
    if gbz_shared_memory_size_gb is not None:
        arguments.append(
            f"--gbz_shared_memory_size_gb={_cli.integer(gbz_shared_memory_size_gb, 'gbz_shared_memory_size_gb', minimum=1)}"
        )
    if intermediate_results_dir is not None:
        arguments.append(
            f"--intermediate_results_dir={_folder_to_use(intermediate_results_dir, 'intermediate_results_dir')}"
        )
    if logging_dir is not None:
        arguments.append(f"--logging_dir={_folder_to_use(logging_dir, 'logging_dir')}")
    if customized_model is not None:
        arguments.append(
            f"--customized_model={_cli.free_text(customized_model, 'customized_model')}"
        )
    if customized_small_model is not None:
        arguments.append(
            f"--customized_small_model={_cli.free_text(customized_small_model, 'customized_small_model')}"
        )
    _boolean(arguments, "disable_small_model", disable_small_model)
    _boolean(arguments, "vcf_stats_report", vcf_stats_report)
    if postprocess_cpus is not None:
        arguments.append(
            f"--postprocess_cpus={_cli.integer(postprocess_cpus, 'postprocess_cpus', minimum=0, maximum=1024)}"
        )
    pan_flags = (
        _MAKE_EXAMPLES_FLAGS | _MAKE_EXAMPLES_PROGRAMS["make_examples_pangenome_aware_dv"][1]
    )
    for flag_name, mapping, allowed in (
        ("make_examples_extra_args", make_examples_extra_args, pan_flags),
        ("call_variants_extra_args", call_variants_extra_args, _CALL_VARIANTS_FLAGS),
        ("postprocess_variants_extra_args", postprocess_variants_extra_args, _POSTPROCESS_FLAGS),
    ):
        _optional(arguments, flag_name, _extra_string(mapping, allowed, flag_name))
    _boolean(arguments, "dry_run", dry_run)
    report = _stats_report_path(vcf) if vcf_stats_report else None
    expected = [] if dry_run else [vcf] + ([gvcf] if gvcf else []) + ([report] if report else [])
    return _run(
        "run_pangenome_aware_deepvariant",
        arguments,
        cwd=vcf.parent,
        expected=expected,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_run_pangenome_aware",
        timeout_seconds=timeout_seconds,
        output_vcf=lambda: str(vcf) if vcf.exists() else None,
        output_gvcf=lambda: str(gvcf) if gvcf and gvcf.exists() else None,
        stats_report=lambda: str(report) if report and report.exists() else None,
    )


def deepvariant_run_deepsomatic(
    model_type: str,
    reference_path: str,
    reads_tumor: str,
    output_vcf: str,
    sample_name_tumor: str,
    *,
    reads_normal: str | None = None,
    sample_name_normal: str | None = None,
    output_gvcf: str | None = None,
    num_shards: int = 1,
    regions: Sequence[str] | None = None,
    intermediate_results_dir: str | None = None,
    logging_dir: str | None = None,
    customized_model: str | None = None,
    customized_model_json: str | None = None,
    disable_small_model: bool = False,
    customized_small_model: str | None = None,
    process_somatic: bool | None = None,
    pon_filtering: str | None = None,
    use_default_pon_filtering: bool | None = None,
    use_candidate_partition: bool | None = None,
    vcf_stats_report: bool = False,
    report_title: str | None = None,
    make_examples_extra_args: Mapping[str, str | int | float | bool] | None = None,
    call_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    postprocess_variants_extra_args: Mapping[str, str | int | float | bool] | None = None,
    dry_run: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call somatic variants from tumor and normal reads with DeepSomatic (``run_deepsomatic``).

    ``model_type`` is WGS, WES, PACBIO, ONT, FFPE_WGS, FFPE_WES or, for tumor-only calling, one of the ``*_TUMOR_ONLY``
    types (FFPE_WGS_TUMOR_ONLY, FFPE_WES_TUMOR_ONLY, WGS_TUMOR_ONLY, WES_TUMOR_ONLY, PACBIO_TUMOR_ONLY,
    ONT_TUMOR_ONLY). Paired types need ``reads_normal`` and ``sample_name_normal``; tumor-only types must not have them.
    ``process_somatic`` marks the input as somatic (None keeps the program default), ``pon_filtering`` is a panel-of-normals
    VCF used to filter calls and ``use_default_pon_filtering`` turns on the default panel for tumor-only models. Other options as
    in ``deepvariant_run``. Result keys: ``output_vcf``, ``output_gvcf`` and ``stats_report``.
    """
    kind = _cli.choice(model_type, _SOMATIC_MODEL_TYPES, "model_type")
    reference = _cli.existing_file(reference_path, "reference_path")
    tumor = _cli.existing_file(reads_tumor, "reads_tumor")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    tumor_only = kind.endswith("_TUMOR_ONLY")
    if tumor_only and (reads_normal is not None or sample_name_normal is not None):
        raise ConfigurationError(
            "Tumor-only model types take no normal reads or normal sample name."
        )
    if not tumor_only and (reads_normal is None or sample_name_normal is None):
        raise ConfigurationError("Paired model types need reads_normal and sample_name_normal.")
    arguments = [
        f"--model_type={kind}",
        f"--ref={reference}",
        f"--reads_tumor={tumor}",
        f"--sample_name_tumor={_cli.text(sample_name_tumor, 'sample_name_tumor')}",
        f"--output_vcf={vcf}",
        f"--num_shards={_cli.integer(num_shards, 'num_shards', minimum=1, maximum=1024)}",
    ]
    if reads_normal is not None and sample_name_normal is not None:
        arguments += [
            f"--reads_normal={_cli.existing_file(reads_normal, 'reads_normal')}",
            f"--sample_name_normal={_cli.text(sample_name_normal, 'sample_name_normal')}",
        ]
    gvcf = None
    if output_gvcf is not None:
        gvcf = _cli.new_file(output_gvcf, "output_gvcf")
        arguments.append(f"--output_gvcf={gvcf}")
    _optional(arguments, "regions", _regions(regions))
    if intermediate_results_dir is not None:
        arguments.append(
            f"--intermediate_results_dir={_folder_to_use(intermediate_results_dir, 'intermediate_results_dir')}"
        )
    if logging_dir is not None:
        arguments.append(f"--logging_dir={_folder_to_use(logging_dir, 'logging_dir')}")
    if customized_model is not None:
        arguments.append(
            f"--customized_model={_cli.free_text(customized_model, 'customized_model')}"
        )
    if customized_model_json is not None:
        arguments.append(
            f"--customized_model_json={_cli.existing_file(customized_model_json, 'customized_model_json')}"
        )
    if customized_small_model is not None:
        arguments.append(
            f"--customized_small_model={_cli.free_text(customized_small_model, 'customized_small_model')}"
        )
    _boolean(arguments, "disable_small_model", disable_small_model)
    _boolean(arguments, "process_somatic", process_somatic)
    if pon_filtering is not None:
        arguments.append(f"--pon_filtering={_cli.existing_file(pon_filtering, 'pon_filtering')}")
    _boolean(arguments, "use_default_pon_filtering", use_default_pon_filtering)
    _boolean(arguments, "use_candidate_partition", use_candidate_partition)
    _boolean(arguments, "vcf_stats_report", vcf_stats_report)
    _optional(
        arguments,
        "report_title",
        _cli.text(report_title, "report_title") if report_title is not None else None,
    )
    somatic_flags = _MAKE_EXAMPLES_FLAGS | _MAKE_EXAMPLES_PROGRAMS["make_examples_somatic"][1]
    for flag_name, mapping, allowed in (
        ("make_examples_extra_args", make_examples_extra_args, somatic_flags),
        ("call_variants_extra_args", call_variants_extra_args, _CALL_VARIANTS_FLAGS),
        ("postprocess_variants_extra_args", postprocess_variants_extra_args, _POSTPROCESS_FLAGS),
    ):
        _optional(arguments, flag_name, _extra_string(mapping, allowed, flag_name))
    _boolean(arguments, "dry_run", dry_run)
    report = _stats_report_path(vcf) if vcf_stats_report else None
    expected = [] if dry_run else [vcf] + ([gvcf] if gvcf else []) + ([report] if report else [])
    return _run(
        "run_deepsomatic",
        arguments,
        cwd=vcf.parent,
        expected=expected,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_run_deepsomatic",
        timeout_seconds=timeout_seconds,
        output_vcf=lambda: str(vcf) if vcf.exists() else None,
        output_gvcf=lambda: str(gvcf) if gvcf and gvcf.exists() else None,
        stats_report=lambda: str(report) if report and report.exists() else None,
    )


def deepvariant_run_oracle_inference(
    model_type: str,
    reference_path: str,
    reads_path: str,
    truth_variants: str,
    confident_regions: str,
    output_vcf: str,
    *,
    labeler_algorithm: str | None = None,
    haploid_contigs: str | None = None,
    par_regions_bed: str | None = None,
    num_shards: int = 1,
    regions: Sequence[str] | None = None,
    sample_name: str | None = None,
    intermediate_results_dir: str | None = None,
    logging_dir: str | None = None,
    make_examples_extra_args: Mapping[str, str | int | float | bool] | None = None,
    dry_run: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the VCF that a perfect classifier would produce from labeled examples (``run_oracle_inference``).

    The program runs make_examples in training mode against a truth VCF and then converts the labeled examples to a VCF
    with ``labeled_examples_to_vcf``; it measures the upper bound that candidate generation allows. ``model_type`` is
    WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA or MASSEQ; ``truth_variants`` is the truth VCF and
    ``confident_regions`` the BED file of confident regions; ``labeler_algorithm`` selects the labeler (for example
    ``haplotype_labeler``). Other options as in ``deepvariant_run``. Result key: ``output_vcf``.
    """
    kind = _cli.choice(model_type, _ORACLE_MODEL_TYPES, "model_type")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    arguments = [
        f"--model_type={kind}",
        f"--ref={_cli.existing_file(reference_path, 'reference_path')}",
        f"--reads={_cli.existing_file(reads_path, 'reads_path')}",
        f"--truth_variants={_cli.existing_file(truth_variants, 'truth_variants')}",
        f"--confident_regions={_cli.existing_file(confident_regions, 'confident_regions')}",
        f"--output_vcf={vcf}",
        f"--num_shards={_cli.integer(num_shards, 'num_shards', minimum=1, maximum=1024)}",
    ]
    _optional(
        arguments,
        "labeler_algorithm",
        _cli.text(labeler_algorithm, "labeler_algorithm") if labeler_algorithm else None,
    )
    _optional(
        arguments,
        "haploid_contigs",
        _cli.text(haploid_contigs, "haploid_contigs") if haploid_contigs else None,
    )
    if par_regions_bed is not None:
        arguments.append(
            f"--par_regions_bed={_cli.existing_file(par_regions_bed, 'par_regions_bed')}"
        )
    _optional(arguments, "regions", _regions(regions))
    _optional(
        arguments,
        "sample_name",
        _cli.text(sample_name, "sample_name") if sample_name is not None else None,
    )
    if intermediate_results_dir is not None:
        arguments.append(
            f"--intermediate_results_dir={_folder_to_use(intermediate_results_dir, 'intermediate_results_dir')}"
        )
    if logging_dir is not None:
        arguments.append(f"--logging_dir={_folder_to_use(logging_dir, 'logging_dir')}")
    _optional(
        arguments,
        "make_examples_extra_args",
        _extra_string(
            make_examples_extra_args,
            _MAKE_EXAMPLES_FLAGS | _MAKE_EXAMPLES_PROGRAMS["make_examples"][1],
            "make_examples_extra_args",
        ),
    )
    _boolean(arguments, "dry_run", dry_run)
    return _run(
        "run_oracle_inference",
        arguments,
        cwd=vcf.parent,
        expected=[] if dry_run else [vcf],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_run_oracle_inference",
        timeout_seconds=timeout_seconds,
        output_vcf=lambda: str(vcf) if vcf.exists() else None,
    )


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------


def deepvariant_make_examples(
    program: str,
    reference_path: str,
    reads: Sequence[str],
    *,
    examples_path: str | None = None,
    mode: str = "calling",
    sample_names: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
    confident_regions: str | None = None,
    truth_variants: str | None = None,
    candidate_positions: str | None = None,
    pangenome_path: str | None = None,
    gvcf_path: str | None = None,
    task: int | None = None,
    checkpoint: str | None = None,
    extra_flags: Mapping[str, str | int | float | bool] | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate pileup examples (and candidate positions) with the DeepVariant ``make_examples`` stage programs.

    One function covers six operations of the work list through ``program`` and ``mode``. ``program`` selects the
    executable and how ``reads`` and ``sample_names`` are read: ``make_examples`` (one BAM/CRAM; germline), ``make_examples_somatic``
    (tumor, then optionally normal reads), ``multisample_make_examples`` (several BAMs of one run; joined with commas),
    ``deeptrio_make_examples`` (child, then parent 1 and parent 2; the program lives in the ``deeptrio`` sub-folder of
    ``deepvariant_dir``) and ``make_examples_pangenome_aware_dv`` (one BAM plus ``pangenome_path``, a GBZ graph). ``mode``
    is ``calling`` (examples for inference), ``training`` (examples with labels; needs ``truth_variants`` and
    ``confident_regions``) or ``candidate_sweep`` (writes ``candidate_positions`` used to cut the genome into intervals
    with equally many candidates; run ``calling`` afterwards). ``examples_path`` is the output TFRecord (gzip) path; it may
    be sharded as ``name.tfrecord@4.gz``, in which case ``task`` (0 to 3) selects the shard written by this call. It is
    required unless ``mode='candidate_sweep'``. ``gvcf_path`` also writes non-variant sites for gVCF output and
    ``checkpoint`` points to the model whose example format is used. ``regions`` limits the genome. ``extra_flags`` holds
    further flags of the chosen program (``--flag=value``); only the documented flags of DeepVariant 1.10.0 for that
    program are accepted. Result keys: ``examples_files`` (the files written) and ``candidate_positions_files``.
    """
    name = _cli.choice(program, tuple(_MAKE_EXAMPLES_PROGRAMS), "program")
    executable, program_flags = _MAKE_EXAMPLES_PROGRAMS[name]
    chosen_mode = _cli.choice(mode, ("calling", "training", "candidate_sweep"), "mode")
    reference = _cli.existing_file(reference_path, "reference_path")
    if isinstance(reads, str) or not reads:
        raise ConfigurationError("reads must be a non-empty list of BAM or CRAM paths.")
    bams = [str(_cli.existing_file(item, "reads")) for item in reads]
    names = list(sample_names) if sample_names is not None else None
    if names is not None:
        names = [_cli.text(item, "sample_names") for item in names]
        if len(names) != len(bams):
            raise ConfigurationError("sample_names must have one entry per reads file.")
    arguments = [f"--ref={reference}", f"--mode={chosen_mode}"]
    if name == "make_examples":
        if len(bams) != 1:
            raise ConfigurationError("make_examples takes exactly one reads file.")
        arguments.append(f"--reads={bams[0]}")
        if names:
            arguments.append(f"--sample_name={names[0]}")
    elif name == "make_examples_somatic":
        if len(bams) not in (1, 2):
            raise ConfigurationError(
                "make_examples_somatic takes the tumor reads and optionally the normal reads."
            )
        arguments.append(f"--reads_tumor={bams[0]}")
        if len(bams) == 2:
            arguments.append(f"--reads_normal={bams[1]}")
        if names:
            arguments.append(f"--sample_name_tumor={names[0]}")
            if len(names) == 2:
                arguments.append(f"--sample_name_normal={names[1]}")
    elif name == "multisample_make_examples":
        if len(bams) < 2 or any("," in item for item in bams):
            raise ConfigurationError(
                "multisample_make_examples takes two or more reads files without commas in their paths."
            )
        arguments.append(f"--reads={','.join(bams)}")
        if names:
            arguments.append(f"--sample_names={','.join(names)}")
    elif name == "deeptrio_make_examples":
        if len(bams) not in (2, 3):
            raise ConfigurationError(
                "deeptrio_make_examples takes the child reads and one or two parents' reads."
            )
        arguments.append(f"--reads={bams[0]}")
        arguments.append(f"--reads_parent1={bams[1]}")
        if len(bams) == 3:
            arguments.append(f"--reads_parent2={bams[2]}")
        if names:
            arguments.append(f"--sample_name={names[0]}")
            arguments.append(f"--sample_name_parent1={names[1]}")
            if len(names) == 3:
                arguments.append(f"--sample_name_parent2={names[2]}")
        if deepvariant_dir is None:
            raise ConfigurationError(
                "deeptrio_make_examples needs deepvariant_dir (the program is deeptrio/make_examples)."
            )
    else:  # make_examples_pangenome_aware_dv
        if len(bams) != 1 or pangenome_path is None:
            raise ConfigurationError(
                "make_examples_pangenome_aware_dv takes one reads file and pangenome_path."
            )
        arguments.append(f"--reads={bams[0]}")
        arguments.append(f"--pangenome={_cli.existing_file(pangenome_path, 'pangenome_path')}")
        if names:
            arguments.append(f"--sample_name_reads={names[0]}")
    if pangenome_path is not None and name != "make_examples_pangenome_aware_dv":
        raise ConfigurationError("pangenome_path only applies to make_examples_pangenome_aware_dv.")
    shard_pattern = None
    if examples_path is not None:
        shard_pattern = _cli.free_text(examples_path, "examples_path")
        parent = Path(shard_pattern).expanduser().resolve().parent
        if not parent.is_dir():
            raise ConfigurationError(f"The folder of examples_path does not exist: {parent}")
        shard_pattern = str(parent / Path(shard_pattern).name)
        if _sharded(shard_pattern):
            raise ConfigurationError("examples_path already exists and will not be overwritten.")
        arguments.append(f"--examples={shard_pattern}")
    elif chosen_mode != "candidate_sweep":
        raise ConfigurationError("examples_path is required unless mode='candidate_sweep'.")
    cwd = Path(shard_pattern).parent if shard_pattern else reference.parent
    if chosen_mode == "training" and (truth_variants is None or confident_regions is None):
        raise ConfigurationError("mode='training' needs truth_variants and confident_regions.")
    if truth_variants is not None:
        arguments.append(f"--truth_variants={_cli.existing_file(truth_variants, 'truth_variants')}")
    if confident_regions is not None:
        arguments.append(
            f"--confident_regions={_cli.existing_file(confident_regions, 'confident_regions')}"
        )
    _optional(arguments, "regions", _regions(regions))
    candidates = None
    if candidate_positions is not None:
        candidates = _cli.new_file(candidate_positions, "candidate_positions")
        arguments.append(f"--candidate_positions={candidates}")
        cwd = candidates.parent if shard_pattern is None else cwd
    elif chosen_mode == "candidate_sweep":
        raise ConfigurationError("mode='candidate_sweep' needs candidate_positions.")
    if gvcf_path is not None:
        arguments.append(f"--gvcf={_cli.free_text(gvcf_path, 'gvcf_path')}")
    if task is not None:
        arguments.append(f"--task={_cli.integer(task, 'task', minimum=0, maximum=100_000)}")
    if checkpoint is not None:
        arguments.append(f"--checkpoint={_cli.free_text(checkpoint, 'checkpoint')}")
    arguments += _flag_arguments(extra_flags, _MAKE_EXAMPLES_FLAGS | program_flags, "extra_flags")
    result = _run(
        executable,
        arguments,
        cwd=cwd,
        expected=[candidates] if candidates else [],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_make_examples",
        timeout_seconds=timeout_seconds,
        program=name,
        mode=chosen_mode,
    )
    produced = _sharded(shard_pattern) if shard_pattern else []
    if shard_pattern and not produced:
        raise BackendExecutionError(
            "make_examples did not write the examples file(s).",
            context={"examples_path": shard_pattern},
        )
    result.update(
        examples_files=[str(path) for path in produced],
        candidate_positions_files=[str(candidates)] if candidates and candidates.exists() else [],
    )
    result["outputs"] = [*result["examples_files"], *result["candidate_positions_files"]]
    return result


def deepvariant_call_variants(
    examples_path: str,
    output_path: str,
    checkpoint: str,
    *,
    batch_size: int | None = None,
    num_readers: int | None = None,
    writer_threads: int | None = None,
    execution_hardware: str | None = None,
    include_debug_info: bool | None = None,
    allow_empty_examples: bool | None = None,
    extra_flags: Mapping[str, str | int | float | bool] | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run the CNN genotype inference on pileup examples with ``call_variants``.

    ``examples_path`` names the TFRecord example files written by ``deepvariant_make_examples`` (a sharded path such as
    ``examples.tfrecord@4.gz`` or a glob; at least one file must exist), ``checkpoint`` the model (a saved-model folder
    or checkpoint prefix, as understood by DeepVariant) and ``output_path`` the new TFRecord (gzip) file of
    ``CallVariantsOutput`` protos. ``batch_size``, ``num_readers``, ``writer_threads``, ``execution_hardware`` (for example
    ``cpu`` or ``accelerator``), ``include_debug_info`` and ``allow_empty_examples`` are passed on; ``extra_flags`` holds
    further ``call_variants`` flags of DeepVariant 1.10.0. Result keys: ``output_file`` and ``input_example_files``.
    """
    pattern = _cli.free_text(examples_path, "examples_path")
    inputs = _sharded(pattern) or sorted(Path(pattern).parent.glob(Path(pattern).name))
    if not inputs:
        raise ConfigurationError(f"No example file matches examples_path: {pattern}")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        f"--examples={pattern}",
        f"--outfile={target}",
        f"--checkpoint={_cli.free_text(checkpoint, 'checkpoint')}",
    ]
    _optional(
        arguments,
        "batch_size",
        _cli.integer(batch_size, "batch_size", minimum=1) if batch_size is not None else None,
    )
    _optional(
        arguments,
        "num_readers",
        _cli.integer(num_readers, "num_readers", minimum=1) if num_readers is not None else None,
    )
    _optional(
        arguments,
        "writer_threads",
        _cli.integer(writer_threads, "writer_threads", minimum=1)
        if writer_threads is not None
        else None,
    )
    _optional(
        arguments,
        "execution_hardware",
        _cli.text(execution_hardware, "execution_hardware") if execution_hardware else None,
    )
    _boolean(arguments, "include_debug_info", include_debug_info)
    _boolean(arguments, "allow_empty_examples", allow_empty_examples)
    arguments += _flag_arguments(extra_flags, _CALL_VARIANTS_FLAGS, "extra_flags")
    return _run(
        "call_variants",
        arguments,
        cwd=target.parent,
        expected=[target],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_call_variants",
        timeout_seconds=timeout_seconds,
        output_file=str(target),
        input_example_files=[str(path) for path in inputs],
    )


def deepvariant_postprocess_variants(
    reference_path: str,
    call_variants_output: str,
    output_vcf: str,
    *,
    output_gvcf: str | None = None,
    nonvariant_site_tfrecord_path: str | None = None,
    small_model_cvo_records: str | None = None,
    sample_name: str | None = None,
    cpus: int | None = None,
    regions: Sequence[str] | None = None,
    haploid_contigs: str | None = None,
    par_regions_bed: str | None = None,
    process_somatic: bool | None = None,
    extra_flags: Mapping[str, str | int | float | bool] | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the model predictions of ``call_variants`` to a VCF (and gVCF) with ``postprocess_variants``.

    ``call_variants_output`` is the TFRecord written by ``deepvariant_call_variants`` (sharded paths and globs are
    passed on unchanged), ``output_vcf`` the new VCF path (``.vcf.gz`` for compressed output) and ``output_gvcf`` a new gVCF
    path, which also needs ``nonvariant_site_tfrecord_path`` (the gVCF records written by make_examples with ``gvcf_path``).
    ``small_model_cvo_records`` are the calls of the small model, ``sample_name`` overrides the sample name, ``cpus`` the
    worker processes, ``regions`` restricts the output, and ``haploid_contigs``/``par_regions_bed`` handle sex
    chromosomes. ``extra_flags`` holds further ``postprocess_variants`` flags of DeepVariant 1.10.0. Result keys:
    ``output_vcf`` and ``output_gvcf``.
    """
    reference = _cli.existing_file(reference_path, "reference_path")
    infile = _cli.free_text(call_variants_output, "call_variants_output")
    if not (_sharded(infile) or list(Path(infile).parent.glob(Path(infile).name))):
        raise ConfigurationError(f"No file matches call_variants_output: {infile}")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    arguments = [f"--ref={reference}", f"--infile={infile}", f"--outfile={vcf}"]
    expected = [vcf]
    gvcf = None
    if output_gvcf is not None:
        if nonvariant_site_tfrecord_path is None:
            raise ConfigurationError("output_gvcf needs nonvariant_site_tfrecord_path.")
        gvcf = _cli.new_file(output_gvcf, "output_gvcf")
        arguments.append(f"--gvcf_outfile={gvcf}")
        expected.append(gvcf)
    if nonvariant_site_tfrecord_path is not None:
        arguments.append(
            f"--nonvariant_site_tfrecord_path={_cli.free_text(nonvariant_site_tfrecord_path, 'nonvariant_site_tfrecord_path')}"
        )
    if small_model_cvo_records is not None:
        arguments.append(
            f"--small_model_cvo_records={_cli.free_text(small_model_cvo_records, 'small_model_cvo_records')}"
        )
    _optional(
        arguments,
        "sample_name",
        _cli.text(sample_name, "sample_name") if sample_name is not None else None,
    )
    _optional(
        arguments,
        "cpus",
        _cli.integer(cpus, "cpus", minimum=0, maximum=1024) if cpus is not None else None,
    )
    _optional(arguments, "regions", _regions(regions))
    _optional(
        arguments,
        "haploid_contigs",
        _cli.text(haploid_contigs, "haploid_contigs") if haploid_contigs else None,
    )
    if par_regions_bed is not None:
        arguments.append(
            f"--par_regions_bed={_cli.existing_file(par_regions_bed, 'par_regions_bed')}"
        )
    _boolean(arguments, "process_somatic", process_somatic)
    blocked = frozenset(
        {
            "gvcf_outfile",
            "nonvariant_site_tfrecord_path",
            "small_model_cvo_records",
            "sample_name",
            "cpus",
            "regions",
            "haploid_contigs",
            "par_regions_bed",
            "process_somatic",
        }
    )
    arguments += _flag_arguments(extra_flags, _POSTPROCESS_FLAGS - blocked, "extra_flags")
    return _run(
        "postprocess_variants",
        arguments,
        cwd=vcf.parent,
        expected=expected,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_postprocess_variants",
        timeout_seconds=timeout_seconds,
        output_vcf=str(vcf),
        output_gvcf=str(gvcf) if gvcf else None,
    )


# ---------------------------------------------------------------------------
# Reports, training and helper tools
# ---------------------------------------------------------------------------


def deepvariant_vcf_stats_report(
    input_vcf: str,
    output_base: str,
    *,
    title: str | None = None,
    max_records: int | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the visual HTML statistics report of a single-sample VCF with ``vcf_stats_report`` (as ``run_deepvariant --vcf_stats_report``).

    ``input_vcf`` needs exactly one sample and a GT field. The report is written to ``<output_base>.visual_report.html``
    (which must not exist). ``title`` is the heading (default: the sample name) and ``max_records`` limits the number of
    VCF records read (default: all). Result keys: ``report_file``.
    """
    source = _cli.existing_file(input_vcf, "input_vcf")
    base = Path(_cli.free_text(output_base, "output_base")).expanduser().resolve()
    if not base.parent.is_dir():
        raise ConfigurationError(f"The folder of output_base does not exist: {base.parent}")
    report = base.with_name(base.name + ".visual_report.html")
    if report.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {report}")
    arguments = [f"--input_vcf={source}", f"--outfile_base={base}"]
    _optional(arguments, "title", _cli.text(title, "title") if title is not None else None)
    _optional(
        arguments,
        "num_records",
        _cli.integer(max_records, "max_records", minimum=1) if max_records is not None else None,
    )
    return _run(
        "vcf_stats_report",
        arguments,
        cwd=base.parent,
        expected=[report],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_vcf_stats_report",
        timeout_seconds=timeout_seconds,
        report_file=str(report),
    )


def deepvariant_train(
    train_dataset_pbtxt: str,
    tune_dataset_pbtxt: str,
    experiment_dir: str,
    *,
    config_name: str = "base",
    config_file: str = "dv_config.py",
    strategy: str = "mirrored",
    config_overrides: Mapping[str, str | int | float | bool] | None = None,
    debug: bool = False,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Train a DeepVariant CNN model from labeled examples with ``train``.

    ``train_dataset_pbtxt`` and ``tune_dataset_pbtxt`` are the dataset configuration files (``.pbtxt``) that list the
    labeled example TFRecords of the training and tuning sets; ``experiment_dir`` is the folder for checkpoints and logs
    (it may exist or be new). ``config_file`` and ``config_name`` select the configuration (``dv_config.py:base`` is the
    DeepVariant default; other names of that file are also accepted) and ``strategy`` is ``mirrored`` (GPU/CPU) or ``tpu``.
    ``config_overrides`` sets configuration values (``--config.<name>=<value>``) for the names batch_size, num_epochs,
    num_validation_examples, warmup_steps, limit, steps_per_iter, shuffle_buffer_elements, init_checkpoint,
    best_checkpoint_metric, optimizer, momentum, use_ema, ema_momentum, optimizer_weight_decay, weight_decay,
    early_stopping_patience, learning_rate, learning_rate_num_epochs_per_decay, learning_rate_decay_rate,
    label_smoothing, backbone_dropout_rate, use_mixed_precision, beta_1, beta_2, adaptive_epsilon and
    include_snp_indel_metrics. Result keys: ``experiment_dir`` and ``files`` (the files of the experiment folder,
    relative, at most 200).
    """
    train_set = _cli.existing_file(train_dataset_pbtxt, "train_dataset_pbtxt")
    tune_set = _cli.existing_file(tune_dataset_pbtxt, "tune_dataset_pbtxt")
    folder = _folder_to_use(experiment_dir, "experiment_dir")
    config = (
        f"{_cli.free_text(config_file, 'config_file')}:{_cli.free_text(config_name, 'config_name')}"
    )
    if not _CONFIG_NAME.match(config_name):
        raise ConfigurationError("config_name must contain only letters, digits and '_'.")
    arguments = [
        f"--config={config}",
        f"--config.train_dataset_pbtxt={train_set}",
        f"--config.tune_dataset_pbtxt={tune_set}",
        f"--experiment_dir={folder}",
        f"--strategy={_cli.choice(strategy, ('mirrored', 'tpu'), 'strategy')}",
    ]
    arguments += [
        item.replace("--", "--config.", 1)
        for item in _flag_arguments(config_overrides, _TRAIN_CONFIG_KEYS, "config_overrides")
    ]
    _boolean(arguments, "debug", debug)
    cwd = folder if folder.exists() else folder.parent
    result = _run(
        "train",
        arguments,
        cwd=cwd,
        expected=[],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_train",
        timeout_seconds=timeout_seconds,
        experiment_dir=str(folder),
    )
    files = (
        sorted(str(path.relative_to(folder)) for path in folder.rglob("*") if path.is_file())
        if folder.exists()
        else []
    )
    result.update(files=files[:200], file_count=len(files))
    return result


def deepvariant_show_examples(
    examples_path: str,
    output_prefix: str,
    *,
    example_info_json: str | None = None,
    vcf: str | None = None,
    image_type: str = "channels",
    regions: Sequence[str] | None = None,
    num_records: int | None = None,
    annotation: bool | None = None,
    truth_labels: bool | None = None,
    scale: int | None = None,
    max_examples_to_scan: int | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Render pileup examples as images with ``show_examples``.

    ``examples_path`` names the example TFRecord file(s) (sharded path or glob), ``example_info_json`` the
    ``example_info.json`` that describes their channels (needed when it is not next to the examples) and ``vcf``
    restricts the output to the variants of a VCF. ``output_prefix`` is the path prefix of the image files; ``image_type`` is
    ``channels`` (a row of channel images), ``RGB`` (stacked channels), ``both`` or ``none`` (no images). ``regions``
    filters by 1-based region or BED path, ``num_records`` limits how many examples are drawn, ``annotation`` and
    ``truth_labels`` add labels, ``scale`` enlarges the images and ``max_examples_to_scan`` bounds the search. Result keys:
    ``image_files`` (new files that start with the output prefix, at most 1000) and ``image_count``.
    """
    pattern = _cli.free_text(examples_path, "examples_path")
    if not (_sharded(pattern) or list(Path(pattern).parent.glob(Path(pattern).name))):
        raise ConfigurationError(f"No example file matches examples_path: {pattern}")
    prefix = Path(_cli.free_text(output_prefix, "output_prefix")).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"The folder of output_prefix does not exist: {prefix.parent}")
    before = {path for path in prefix.parent.glob(prefix.name + "*")}
    arguments = [
        f"--examples={pattern}",
        f"--output={prefix}",
        f"--image_type={_cli.choice(image_type, ('channels', 'RGB', 'both', 'none'), 'image_type')}",
    ]
    if example_info_json is not None:
        arguments.append(
            f"--example_info_json={_cli.existing_file(example_info_json, 'example_info_json')}"
        )
    if vcf is not None:
        arguments.append(f"--vcf={_cli.existing_file(vcf, 'vcf')}")
    _optional(arguments, "regions", _regions(regions))
    _optional(
        arguments,
        "num_records",
        _cli.integer(num_records, "num_records", minimum=1) if num_records is not None else None,
    )
    _boolean(arguments, "annotation", annotation)
    _boolean(arguments, "truth_labels", truth_labels)
    _optional(
        arguments,
        "scale",
        _cli.integer(scale, "scale", minimum=1, maximum=100) if scale is not None else None,
    )
    _optional(
        arguments,
        "max_examples_to_scan",
        _cli.integer(max_examples_to_scan, "max_examples_to_scan", minimum=1)
        if max_examples_to_scan is not None
        else None,
    )
    result = _run(
        "show_examples",
        arguments,
        cwd=prefix.parent,
        expected=[],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_show_examples",
        timeout_seconds=timeout_seconds,
    )
    created = sorted(
        str(path) for path in prefix.parent.glob(prefix.name + "*") if path not in before
    )
    result.update(image_files=created[:1000], image_count=len(created), outputs=created[:1000])
    return result


def deepvariant_labeled_examples_to_vcf(
    reference_path: str,
    examples_path: str,
    output_vcf: str,
    *,
    sample_name: str | None = None,
    max_records: int | None = None,
    allow_unlabeled_examples: bool | None = None,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert labeled training examples back to a VCF of their true genotypes with ``labeled_examples_to_vcf``.

    ``examples_path`` names the labeled example TFRecord file(s) written by ``deepvariant_make_examples`` in training mode.
    ``sample_name`` is the sample name written in the VCF, ``max_records`` limits the examples read and
    ``allow_unlabeled_examples`` also accepts examples without a label. Result key: ``output_vcf``.
    """
    reference = _cli.existing_file(reference_path, "reference_path")
    pattern = _cli.free_text(examples_path, "examples_path")
    if not (_sharded(pattern) or list(Path(pattern).parent.glob(Path(pattern).name))):
        raise ConfigurationError(f"No example file matches examples_path: {pattern}")
    vcf = _cli.new_file(output_vcf, "output_vcf")
    arguments = [f"--ref={reference}", f"--examples={pattern}", f"--output_vcf={vcf}"]
    _optional(
        arguments,
        "sample_name",
        _cli.text(sample_name, "sample_name") if sample_name is not None else None,
    )
    _optional(
        arguments,
        "max_records",
        _cli.integer(max_records, "max_records", minimum=1) if max_records is not None else None,
    )
    _boolean(arguments, "allow_unlabeled_examples", allow_unlabeled_examples)
    return _run(
        "labeled_examples_to_vcf",
        arguments,
        cwd=vcf.parent,
        expected=[vcf],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_labeled_examples_to_vcf",
        timeout_seconds=timeout_seconds,
        output_vcf=str(vcf),
    )


def deepvariant_compare_labelers(
    labeled_examples_a: str,
    labeled_examples_b: str,
    output_path: str,
    *,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the variant labels of two labeled example sets (two labelers) with ``compare_labelers``.

    The example sets (TFRecord, one per labeler) must come from the same candidates. The comparison report is written to
    ``output_path``. ``compare_labelers`` is a source-only tool of the DeepVariant repository (``//deepvariant/labeler:compare_labelers``
    in Bazel); build or wrap it first. Result keys: ``output_file`` and ``lines`` (the report lines, at most 200).
    """
    first = _cli.free_text(labeled_examples_a, "labeled_examples_a")
    second = _cli.free_text(labeled_examples_b, "labeled_examples_b")
    for value, label in ((first, "labeled_examples_a"), (second, "labeled_examples_b")):
        if not (_sharded(value) or list(Path(value).parent.glob(Path(value).name))):
            raise ConfigurationError(f"No file matches {label}: {value}")
    target = _cli.new_file(output_path, "output_path")
    arguments = [
        f"--labeled_examples_a={first}",
        f"--labeled_examples_b={second}",
        f"--output={target}",
    ]
    result = _run(
        "compare_labelers",
        arguments,
        cwd=target.parent,
        expected=[target],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_compare_labelers",
        timeout_seconds=timeout_seconds,
        output_file=str(target),
    )
    lines = target.read_text(encoding="utf-8", errors="replace").splitlines()
    result["lines"] = lines[:200]
    return result


def deepvariant_preprocess_truth(
    truth_vcf: str,
    output_vcf: str,
    *,
    timeout_seconds: float = 86_400.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalise the alleles and genotypes of a phased truth VCF (star alleles, overlapping variants) with ``preprocess_truth``.

    Use it only for truth VCFs that are fully phased, such as the T2T Q100 truth sets. ``truth_vcf`` is the input and
    ``output_vcf`` the new normalised VCF. ``preprocess_truth`` is a source-only tool of the DeepVariant repository
    (``tools/preprocess_truth.py``); wrap or build it first. Result key: ``output_vcf``.
    """
    source = _cli.existing_file(truth_vcf, "truth_vcf")
    target = _cli.new_file(output_vcf, "output_vcf")
    return _run(
        "preprocess_truth",
        [f"--truth_vcf={source}", f"--output_vcf={target}"],
        cwd=target.parent,
        expected=[target],
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_preprocess_truth",
        timeout_seconds=timeout_seconds,
        output_vcf=str(target),
    )


def deepvariant_print_f1(
    metrics_dir: str,
    output_path: str,
    *,
    timeout_seconds: float = 3600.0,
    deepvariant_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the F1 score of every checkpoint from the evaluation metrics files of a training run with ``print_f1``.

    ``metrics_dir`` holds the ``*.metrics`` JSON files that ``train`` writes (names ending ``ckpt-<number>.metrics``).
    ``print_f1`` prints one tab-separated line per file; it is saved in ``output_path`` and returned as ``records``
    (``checkpoint``, ``tps_plus_fns`` and ``f1_all``). ``print_f1`` is a source-only tool of the DeepVariant repository
    (``tools/print_f1.py``); wrap or build it first.
    """
    folder = _cli.existing_directory(metrics_dir, "metrics_dir")
    target = _cli.new_file(output_path, "output_path")
    result = _run(
        "print_f1",
        [f"--metrics_dir={folder}"],
        cwd=target.parent,
        expected=[],
        stdout_to=target,
        deepvariant_dir=deepvariant_dir,
        operation="deepvariant_print_f1",
        timeout_seconds=timeout_seconds,
        output_file=str(target),
    )
    records: list[dict[str, Any]] = []
    for line in target.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) == 3:
            records.append(
                {
                    "checkpoint": _cli.convert_number(fields[0]),
                    "tps_plus_fns": _cli.convert_number(fields[1]),
                    "f1_all": _cli.convert_number(fields[2]),
                }
            )
    result["records"] = records
    return result


__all__ = [
    "deepvariant_call_variants",
    "deepvariant_compare_labelers",
    "deepvariant_labeled_examples_to_vcf",
    "deepvariant_make_examples",
    "deepvariant_postprocess_variants",
    "deepvariant_preprocess_truth",
    "deepvariant_print_f1",
    "deepvariant_run",
    "deepvariant_run_deepsomatic",
    "deepvariant_run_deeptrio",
    "deepvariant_run_oracle_inference",
    "deepvariant_run_pangenome_aware",
    "deepvariant_show_examples",
    "deepvariant_train",
    "deepvariant_vcf_stats_report",
]
_register_public_verbose(globals())
