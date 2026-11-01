"""Picard SAM/BAM/CRAM and Illumina-run tools as typed DNAKit facades (integration pack P08, part
2).

The user installs Picard (``conda install -c bioconda picard``, or ``picard.jar`` plus a Java
runtime). Each function runs one Picard tool with an explicit, typed ``KEY=value`` argument list (no
shell). Picard is found from ``picard_path`` (a ``picard`` launcher or a ``picard.jar``), then
``picard`` on PATH, then the ``PICARD_JAR`` environment variable; a missing program raises
``BackendUnavailableError`` with an install hint. Outputs are new explicit paths: an existing output
is never overwritten. Results are JSON-serialisable dictionaries with ``status`` (``completed`` only
when the exit code is 0 and every expected output exists), ``backend``, ``tool``, ``version``,
``command``, ``return_code``, ``outputs`` and ``log_tail``.

Optional arguments left as ``None`` are not passed, so the defaults are exactly Picard's own. The
output format of a SAM/BAM/CRAM file follows its extension, as in Picard. Tested with Picard 3.5.0.

Not exposed on purpose: Picard's ``FilterSamReads FILTER=includeJavascript`` (it would run
caller-supplied scripts).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _picard_cli as _p
from dnakit.exceptions import ConfigurationError

_SORT_ORDERS = ("unsorted", "queryname", "coordinate", "duplicate", "unknown")
_ADAPTER_PAIRS = (
    "PAIRED_END",
    "INDEXED",
    "SINGLE_END",
    "NEXTERA_V1",
    "NEXTERA_V2",
    "DUAL_INDEXED",
    "FLUIDIGM",
    "TRUSEQ_SMALLRNA",
    "ALTERNATIVE_SINGLE_END",
)
_FILTERS = (
    "includeAligned",
    "excludeAligned",
    "includeReadList",
    "excludeReadList",
    "includePairedIntervals",
    "includeTagValues",
    "excludeTagValues",
)
_PRIMARY_STRATEGIES = ("BestMapq", "EarliestFragment", "BestEndMapq", "MostDistant")
_TAGGING_POLICIES = ("DontTag", "OpticalOnly", "All")
_SCORING = ("SUM_OF_BASE_QUALITIES", "TOTAL_MAPPED_REFERENCE_LENGTH", "RANDOM")
_READ_NAME_FORMATS = ("CASAVA_1_8", "ILLUMINA")
_DISTANCE_MODES = ("HAMMING", "LENIENT_HAMMING", "FREE")
_BARCODE_STRATEGIES = ("ORPHANS_ONLY", "INEXACT_MATCH", "ALWAYS")
_PAIR_ORIENTATIONS = ("FR", "RF", "TANDEM")
_UNMAPPED_STRATEGIES = ("COPY_TO_TAG", "DO_NOT_CHANGE", "DO_NOT_CHANGE_INVALID", "MOVE_TO_TAG")
_DOWNSAMPLE_STRATEGIES = ("HighAccuracy", "ConstantMemory", "Chained")
_ALIGNED_FORMATS = ("sam", "bam", "cram", "dynamic")


def _sam_common(
    validation_stringency: str | None, reference_sequence: str | Path | None
) -> list[str]:
    return _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency)) + _p.a_file(
        "REFERENCE_SEQUENCE", reference_sequence, "reference_sequence"
    )


def _metrics(result: dict[str, Any], path: Path | None, key: str = "metrics") -> dict[str, Any]:
    if path is not None and path.is_file():
        parsed = _p.read_metrics(path)
        result[key] = parsed["metrics"]
        if parsed["histogram"]:
            result["histogram"] = parsed["histogram"]
    return result


def _files_like(directory: Path, pattern: str) -> list[Path]:
    return sorted(path for path in directory.glob(pattern) if path.is_file())


# ---------------------------------------------------------------------------
# Subsetting and downsampling
# ---------------------------------------------------------------------------


def picard_downsample_sam(
    input_path: str,
    output_path: str,
    *,
    probability: float | None = None,
    strategy: str | None = None,
    random_seed: int | None = None,
    accuracy: float | None = None,
    metrics_file: str | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep a deterministic random subset of the reads of a SAM/BAM/CRAM file, by query name
    (``DownsampleSam``).

    Reads of one template (pairs, secondary and supplementary alignments) are kept or dropped
    together. ``probability`` is the chance of keeping a template (Picard's default is 1.0, which
    keeps everything); ``strategy`` is ``ConstantMemory`` (default), ``HighAccuracy`` or
    ``Chained``; ``random_seed`` makes the choice reproducible. ``metrics_file`` receives a
    ``QualityYieldMetrics`` file of the downsampled output.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    metrics = _p.optional_new_file(metrics_file, "metrics_file")
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_num("PROBABILITY", probability, "probability", minimum=0, maximum=1)
    args += _p.a_choice("STRATEGY", strategy, _DOWNSAMPLE_STRATEGIES, "strategy")
    args += _p.a_int("RANDOM_SEED", random_seed, "random_seed")
    args += _p.a_num("ACCURACY", accuracy, "accuracy", minimum=0)
    args += _p.kv("METRICS_FILE", metrics)
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"output": target, **({"metrics_file": metrics} if metrics is not None else {})}
    result = _p.run(
        "DownsampleSam",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics(result, metrics)


def picard_position_based_downsample_sam(
    input_path: str,
    output_path: str,
    fraction: float,
    *,
    read_name_regex: str | None = None,
    remove_duplicate_information: bool | None = None,
    allow_multiple_downsampling_despite_warnings: bool | None = None,
    stop_after: int | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Downsample a SAM/BAM/CRAM file by flowcell position, keeping approximately ``fraction`` of
    the templates (``PositionBasedDownsampleSam``).

    The choice is deterministic from each read's tile and x/y cluster coordinates (parsed from the
    read name with ``read_name_regex``), not random, so optical-duplicate structure is kept. Reads
    of one template stay together. ``fraction`` must be between 0 and 1.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    args = [
        f"INPUT={source}",
        f"OUTPUT={target}",
        f"FRACTION={_p.number(fraction, 'fraction', minimum=0, maximum=1)}",
    ]
    args += _p.a_text("READ_NAME_REGEX", read_name_regex, "read_name_regex")
    args += _p.a_bool(
        "REMOVE_DUPLICATE_INFORMATION", remove_duplicate_information, "remove_duplicate_information"
    )
    args += _p.a_bool(
        "ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS",
        allow_multiple_downsampling_despite_warnings,
        "allow_multiple_downsampling_despite_warnings",
    )
    args += _p.a_int("STOP_AFTER", stop_after, "stop_after", minimum=0)
    args += _sam_common(validation_stringency, reference_sequence)
    return _p.run(
        "PositionBasedDownsampleSam",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )


def picard_filter_sam_reads(
    input_path: str,
    output_path: str,
    filter: str,
    *,
    read_list_file: str | None = None,
    interval_list: str | None = None,
    tag: str | None = None,
    tag_values: Sequence[str] | None = None,
    sort_order: str | None = None,
    write_reads_files: bool | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Subset the reads of a SAM/BAM/CRAM file by one filter (``FilterSamReads``).

    ``filter`` is one of ``includeAligned`` / ``excludeAligned`` (the input must be queryname sorted
    and both mates must be aligned to be included), ``includeReadList`` / ``excludeReadList``
    (``read_list_file`` holds one read name per line), ``includePairedIntervals``
    (``interval_list``) and ``includeTagValues`` / ``excludeTagValues`` (``tag`` and string
    ``tag_values``). The JavaScript filter is not offered. ``sort_order`` defaults to the input's.
    """
    source = _p.existing_file(input_path, "input_path")
    target = _p.new_file(output_path, "output_path")
    mode = _p.choice(filter, _FILTERS, "filter")
    if mode in ("includeReadList", "excludeReadList") and read_list_file is None:
        raise ConfigurationError("read_list_file is required for the read-list filters.")
    if mode == "includePairedIntervals" and interval_list is None:
        raise ConfigurationError("interval_list is required for filter='includePairedIntervals'.")
    if mode in ("includeTagValues", "excludeTagValues") and (tag is None or not tag_values):
        raise ConfigurationError("tag and tag_values are required for the tag-value filters.")
    args = [f"INPUT={source}", f"OUTPUT={target}", f"FILTER={mode}"]
    args += _p.a_file("READ_LIST_FILE", read_list_file, "read_list_file")
    args += _p.a_file("INTERVAL_LIST", interval_list, "interval_list")
    args += _p.a_text("TAG", tag, "tag")
    args += _p.a_texts("TAG_VALUE", tag_values, "tag_values")
    args += _p.a_choice("SORT_ORDER", sort_order, _SORT_ORDERS, "sort_order")
    args += _p.a_bool("WRITE_READS_FILES", write_reads_files, "write_reads_files")
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "FilterSamReads",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return result


def picard_split_sam_by_number_of_reads(
    input_path: str,
    output_dir: str,
    *,
    split_to_n_files: int | None = None,
    split_to_n_reads: int | None = None,
    out_prefix: str | None = None,
    total_reads_in_input: int | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a query-grouped SAM/BAM/CRAM file into several files with the same sort order
    (``SplitSamByNumberOfReads``).

    Give exactly one of ``split_to_n_files`` or ``split_to_n_reads`` (about that many reads per
    file). ``output_dir`` is an existing directory; the files ``<out_prefix>_N.<ext>`` (``shard`` by
    default, extension of the input) are written into it and must not exist yet. The new files are
    listed in ``files``. ``total_reads_in_input`` saves a first pass over the input.
    """
    source = _p.existing_file(input_path, "input_path")
    folder = _p.existing_directory(output_dir, "output_dir")
    if (split_to_n_files is None) == (split_to_n_reads is None):
        raise ConfigurationError("Give exactly one of split_to_n_files and split_to_n_reads.")
    prefix = "shard" if out_prefix is None else _p.text(out_prefix, "out_prefix")
    if _files_like(folder, f"{prefix}_*"):
        raise ConfigurationError(f"output_dir already holds files named {prefix}_*: {folder}")
    args = [f"INPUT={source}", f"OUTPUT={folder}"]
    args += _p.a_int("SPLIT_TO_N_FILES", split_to_n_files, "split_to_n_files", minimum=1)
    args += _p.a_int("SPLIT_TO_N_READS", split_to_n_reads, "split_to_n_reads", minimum=1)
    args += _p.a_text("OUT_PREFIX", out_prefix, "out_prefix")
    args += _p.a_int(
        "TOTAL_READS_IN_INPUT", total_reads_in_input, "total_reads_in_input", minimum=1
    )
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "SplitSamByNumberOfReads",
        args,
        outputs={},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        cwd=folder,
    )
    files = _files_like(folder, f"{prefix}_*")
    result["files"] = [str(path) for path in files]
    if result["status"] == "completed" and not files:
        result["status"] = "failed"
    return result


# ---------------------------------------------------------------------------
# Reshaping alignments
# ---------------------------------------------------------------------------


def picard_reorder_sam(
    input_path: str,
    output_path: str,
    sequence_dictionary: str,
    *,
    allow_contig_length_discordance: bool | None = None,
    allow_incomplete_dict_concordance: bool | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reorder the reads of a SAM/BAM/CRAM file to match the contig order of a new reference
    (``ReorderSam``).

    Contigs are matched by exact name; reads on contigs absent from the new reference are unmapped.
    ``sequence_dictionary`` may be a SAM/BAM/CRAM, VCF/BCF, interval list, FASTA or ``.dict`` file.
    Both ``allow_*`` flags are dangerous relaxations and default to Picard's false.
    """
    source = _p.existing_file(input_path, "input_path")
    dictionary = _p.existing_file(sequence_dictionary, "sequence_dictionary")
    target = _p.new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"OUTPUT={target}", f"SEQUENCE_DICTIONARY={dictionary}"]
    args += _p.a_bool(
        "ALLOW_CONTIG_LENGTH_DISCORDANCE",
        allow_contig_length_discordance,
        "allow_contig_length_discordance",
    )
    args += _p.a_bool(
        "ALLOW_INCOMPLETE_DICT_CONCORDANCE",
        allow_incomplete_dict_concordance,
        "allow_incomplete_dict_concordance",
    )
    args += _sam_common(validation_stringency, reference_sequence)
    return _p.run(
        "ReorderSam",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )


def picard_revert_sam(
    input_path: str,
    output_path: str,
    *,
    output_by_readgroup: bool | None = None,
    output_by_readgroup_file_format: str | None = None,
    remove_alignment_information: bool | None = None,
    remove_duplicate_information: bool | None = None,
    restore_original_qualities: bool | None = None,
    restore_hardclips: bool | None = None,
    attribute_to_clear: Sequence[str] | None = None,
    sanitize: bool | None = None,
    keep_first_duplicate: bool | None = None,
    max_discard_fraction: float | None = None,
    sample_alias: str | None = None,
    library_name: str | None = None,
    sort_order: str | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Revert an aligned SAM/BAM/CRAM file to an unaligned file, restoring original qualities
    (``RevertSam``).

    With Picard's defaults the alignment and duplicate flags are removed, original qualities
    (``OQ``) are restored, hard-clipped bases are restored and the output is queryname sorted.
    ``attribute_to_clear`` replaces nothing: its tags are added to Picard's default list (``NM UQ PG
    MD MQ SA MC AS``). With ``output_by_readgroup=True``, ``output_path`` is a new directory that
    receives one file per read group (``output_by_readgroup_file_format`` is ``sam``, ``bam``,
    ``cram`` or ``dynamic``, by default the format of the input); its files are listed in ``files``.
    ``sanitize`` is destructive and discards inconsistent reads.
    """
    source = _p.existing_file(input_path, "input_path")
    by_group = output_by_readgroup is not None and _p.boolean(
        output_by_readgroup, "output_by_readgroup"
    )
    target = (
        _p.existing_directory(output_path, "output_path")
        if by_group
        else _p.new_file(output_path, "output_path")
    )
    existing = set(_files_like(target, "*")) if by_group else set()
    args = [f"INPUT={source}", f"OUTPUT={target}"]
    args += _p.a_bool("OUTPUT_BY_READGROUP", output_by_readgroup, "output_by_readgroup")
    args += _p.a_choice(
        "OUTPUT_BY_READGROUP_FILE_FORMAT",
        output_by_readgroup_file_format,
        _ALIGNED_FORMATS,
        "output_by_readgroup_file_format",
    )
    args += _p.a_bool(
        "REMOVE_ALIGNMENT_INFORMATION", remove_alignment_information, "remove_alignment_information"
    )
    args += _p.a_bool(
        "REMOVE_DUPLICATE_INFORMATION", remove_duplicate_information, "remove_duplicate_information"
    )
    args += _p.a_bool(
        "RESTORE_ORIGINAL_QUALITIES", restore_original_qualities, "restore_original_qualities"
    )
    args += _p.a_bool("RESTORE_HARDCLIPS", restore_hardclips, "restore_hardclips")
    args += _p.a_texts("ATTRIBUTE_TO_CLEAR", attribute_to_clear, "attribute_to_clear")
    args += _p.a_bool("SANITIZE", sanitize, "sanitize")
    args += _p.a_bool("KEEP_FIRST_DUPLICATE", keep_first_duplicate, "keep_first_duplicate")
    args += _p.a_num(
        "MAX_DISCARD_FRACTION", max_discard_fraction, "max_discard_fraction", minimum=0, maximum=1
    )
    args += _p.a_text("SAMPLE_ALIAS", sample_alias, "sample_alias")
    args += _p.a_text("LIBRARY_NAME", library_name, "library_name")
    args += _p.a_choice("SORT_ORDER", sort_order, _SORT_ORDERS, "sort_order")
    args += _sam_common(validation_stringency, reference_sequence)
    result = _p.run(
        "RevertSam",
        args,
        outputs={} if by_group else {"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        cwd=target if by_group else target.parent,
    )
    if by_group:
        created = [path for path in _files_like(target, "*") if path not in existing]
        result["files"] = [str(path) for path in created]
        if result["status"] == "completed" and not created:
            result["status"] = "failed"
    return result


def picard_mark_illumina_adapters(
    input_path: str,
    metrics_path: str,
    *,
    output_path: str | None = None,
    adapters: Sequence[str] | None = None,
    five_prime_adapter: str | None = None,
    three_prime_adapter: str | None = None,
    max_error_rate_pe: float | None = None,
    max_error_rate_se: float | None = None,
    min_match_bases_pe: int | None = None,
    min_match_bases_se: int | None = None,
    adapter_truncation_length: int | None = None,
    num_adapters_to_keep: int | None = None,
    prune_adapter_list_after_this_many_adapters_seen: int | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find Illumina adapter sequences in queryname-sorted reads and write the ``XT`` trimming tag
    (``MarkIlluminaAdapters``).

    ``metrics_path`` receives a histogram of bases clipped per read (parsed in ``histogram``).
    Without ``output_path`` only the metrics are produced. ``adapters`` are Picard adapter pairs
    (``PAIRED_END``, ``INDEXED``, ``SINGLE_END``, ``NEXTERA_V1``, ``NEXTERA_V2``, ``DUAL_INDEXED``,
    ``FLUIDIGM``, ``TRUSEQ_SMALLRNA``, ``ALTERNATIVE_SINGLE_END``); by default ``INDEXED``,
    ``DUAL_INDEXED`` and ``PAIRED_END``. Existing ``XT`` tags are cleared.
    """
    source = _p.existing_file(input_path, "input_path")
    metrics = _p.new_file(metrics_path, "metrics_path")
    target = _p.optional_new_file(output_path, "output_path")
    args = [f"INPUT={source}", f"METRICS={metrics}"]
    args += _p.kv("OUTPUT", target)
    args += _p.a_choices("ADAPTERS", adapters, _ADAPTER_PAIRS, "adapters")
    args += _p.a_text("FIVE_PRIME_ADAPTER", five_prime_adapter, "five_prime_adapter")
    args += _p.a_text("THREE_PRIME_ADAPTER", three_prime_adapter, "three_prime_adapter")
    args += _p.a_num(
        "MAX_ERROR_RATE_PE", max_error_rate_pe, "max_error_rate_pe", minimum=0, maximum=1
    )
    args += _p.a_num(
        "MAX_ERROR_RATE_SE", max_error_rate_se, "max_error_rate_se", minimum=0, maximum=1
    )
    args += _p.a_int("MIN_MATCH_BASES_PE", min_match_bases_pe, "min_match_bases_pe", minimum=1)
    args += _p.a_int("MIN_MATCH_BASES_SE", min_match_bases_se, "min_match_bases_se", minimum=1)
    args += _p.a_int(
        "ADAPTER_TRUNCATION_LENGTH",
        adapter_truncation_length,
        "adapter_truncation_length",
        minimum=1,
    )
    args += _p.a_int(
        "NUM_ADAPTERS_TO_KEEP", num_adapters_to_keep, "num_adapters_to_keep", minimum=1
    )
    args += _p.a_int(
        "PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN",
        prune_adapter_list_after_this_many_adapters_seen,
        "prune_adapter_list_after_this_many_adapters_seen",
        minimum=-1,
    )
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"metrics": metrics, **({"output": target} if target is not None else {})}
    result = _p.run(
        "MarkIlluminaAdapters",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )
    return _metrics(result, metrics)


def picard_merge_bam_alignment(
    unmapped_bam: str,
    reference_sequence: str,
    output_path: str,
    *,
    aligned_bams: Sequence[str] | None = None,
    read1_aligned_bams: Sequence[str] | None = None,
    read2_aligned_bams: Sequence[str] | None = None,
    sort_order: str | None = None,
    add_mate_cigar: bool | None = None,
    clip_adapters: bool | None = None,
    clip_overlapping_reads: bool | None = None,
    hard_clip_overlapping_reads: bool | None = None,
    include_secondary_alignments: bool | None = None,
    aligned_reads_only: bool | None = None,
    aligner_proper_pair_flags: bool | None = None,
    is_bisulfite_sequence: bool | None = None,
    max_insertions_or_deletions: int | None = None,
    primary_alignment_strategy: str | None = None,
    attributes_to_retain: Sequence[str] | None = None,
    attributes_to_remove: Sequence[str] | None = None,
    attributes_to_reverse: Sequence[str] | None = None,
    attributes_to_reverse_complement: Sequence[str] | None = None,
    expected_orientations: Sequence[str] | None = None,
    read1_trim: int | None = None,
    read2_trim: int | None = None,
    unmap_contaminant_reads: bool | None = None,
    unmapped_read_strategy: str | None = None,
    min_unclipped_bases: int | None = None,
    add_pg_tag_to_reads: bool | None = None,
    program_record_id: str | None = None,
    program_group_name: str | None = None,
    program_group_version: str | None = None,
    program_group_command_line: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge aligner output with the original unmapped BAM, restoring the tags the aligner dropped
    (``MergeBamAlignment``).

    ``unmapped_bam`` must be queryname sorted and contain only unmapped reads. Give the alignments
    either as ``aligned_bams`` or as separate ``read1_aligned_bams`` and ``read2_aligned_bams`` (not
    both). The merged file is coordinate sorted unless ``sort_order`` says otherwise.
    ``attributes_to_reverse`` and ``attributes_to_reverse_complement`` replace Picard's defaults
    (``OQ U2`` and ``E2 SQ``) only when given. ``primary_alignment_strategy`` is ``BestMapq``
    (Picard's default), ``EarliestFragment``, ``BestEndMapq`` or ``MostDistant``.
    """
    unmapped = _p.existing_file(unmapped_bam, "unmapped_bam")
    reference = _p.existing_file(reference_sequence, "reference_sequence")
    target = _p.new_file(output_path, "output_path")
    if aligned_bams and (read1_aligned_bams or read2_aligned_bams):
        raise ConfigurationError(
            "Use either aligned_bams or read1_aligned_bams/read2_aligned_bams, not both."
        )
    if not aligned_bams and not (read1_aligned_bams and read2_aligned_bams):
        raise ConfigurationError(
            "Give aligned_bams, or both read1_aligned_bams and read2_aligned_bams."
        )
    args = [f"UNMAPPED_BAM={unmapped}", f"REFERENCE_SEQUENCE={reference}", f"OUTPUT={target}"]
    args += _p.a_files("ALIGNED_BAM", aligned_bams, "aligned_bams")
    args += _p.a_files("READ1_ALIGNED_BAM", read1_aligned_bams, "read1_aligned_bams")
    args += _p.a_files("READ2_ALIGNED_BAM", read2_aligned_bams, "read2_aligned_bams")
    args += _p.a_choice("SORT_ORDER", sort_order, _SORT_ORDERS, "sort_order")
    args += _p.a_bool("ADD_MATE_CIGAR", add_mate_cigar, "add_mate_cigar")
    args += _p.a_bool("CLIP_ADAPTERS", clip_adapters, "clip_adapters")
    args += _p.a_bool("CLIP_OVERLAPPING_READS", clip_overlapping_reads, "clip_overlapping_reads")
    args += _p.a_bool(
        "HARD_CLIP_OVERLAPPING_READS", hard_clip_overlapping_reads, "hard_clip_overlapping_reads"
    )
    args += _p.a_bool(
        "INCLUDE_SECONDARY_ALIGNMENTS", include_secondary_alignments, "include_secondary_alignments"
    )
    args += _p.a_bool("ALIGNED_READS_ONLY", aligned_reads_only, "aligned_reads_only")
    args += _p.a_bool(
        "ALIGNER_PROPER_PAIR_FLAGS", aligner_proper_pair_flags, "aligner_proper_pair_flags"
    )
    args += _p.a_bool("IS_BISULFITE_SEQUENCE", is_bisulfite_sequence, "is_bisulfite_sequence")
    args += _p.a_int(
        "MAX_INSERTIONS_OR_DELETIONS",
        max_insertions_or_deletions,
        "max_insertions_or_deletions",
        minimum=-1,
    )
    args += _p.a_choice(
        "PRIMARY_ALIGNMENT_STRATEGY",
        primary_alignment_strategy,
        _PRIMARY_STRATEGIES,
        "primary_alignment_strategy",
    )
    args += _p.a_texts("ATTRIBUTES_TO_RETAIN", attributes_to_retain, "attributes_to_retain")
    args += _p.a_texts("ATTRIBUTES_TO_REMOVE", attributes_to_remove, "attributes_to_remove")
    args += _p.a_texts("ATTRIBUTES_TO_REVERSE", attributes_to_reverse, "attributes_to_reverse")
    args += _p.a_texts(
        "ATTRIBUTES_TO_REVERSE_COMPLEMENT",
        attributes_to_reverse_complement,
        "attributes_to_reverse_complement",
    )
    args += _p.a_choices(
        "EXPECTED_ORIENTATIONS", expected_orientations, _PAIR_ORIENTATIONS, "expected_orientations"
    )
    args += _p.a_int("READ1_TRIM", read1_trim, "read1_trim", minimum=0)
    args += _p.a_int("READ2_TRIM", read2_trim, "read2_trim", minimum=0)
    args += _p.a_bool("UNMAP_CONTAMINANT_READS", unmap_contaminant_reads, "unmap_contaminant_reads")
    args += _p.a_choice(
        "UNMAPPED_READ_STRATEGY",
        unmapped_read_strategy,
        _UNMAPPED_STRATEGIES,
        "unmapped_read_strategy",
    )
    args += _p.a_int("MIN_UNCLIPPED_BASES", min_unclipped_bases, "min_unclipped_bases", minimum=0)
    args += _p.a_bool("ADD_PG_TAG_TO_READS", add_pg_tag_to_reads, "add_pg_tag_to_reads")
    args += _p.a_text("PROGRAM_RECORD_ID", program_record_id, "program_record_id")
    args += _p.a_text("PROGRAM_GROUP_NAME", program_group_name, "program_group_name")
    args += _p.a_text("PROGRAM_GROUP_VERSION", program_group_version, "program_group_version")
    args += _p.a_text(
        "PROGRAM_GROUP_COMMAND_LINE", program_group_command_line, "program_group_command_line"
    )
    args += _p.kv("VALIDATION_STRINGENCY", _p.stringency(validation_stringency))
    return _p.run(
        "MergeBamAlignment",
        args,
        outputs={"output": target},
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
    )


# ---------------------------------------------------------------------------
# Duplicate marking
# ---------------------------------------------------------------------------


def picard_mark_duplicates_with_mate_cigar(
    input_paths: Sequence[str],
    output_path: str,
    metrics_file: str,
    *,
    umi_aware: bool = False,
    umi_metrics_file: str | None = None,
    umi_tag_name: str | None = None,
    max_edit_distance_to_join: int | None = None,
    allow_missing_umis: bool | None = None,
    remove_duplicates: bool | None = None,
    remove_sequencing_duplicates: bool | None = None,
    tagging_policy: str | None = None,
    tag_duplicate_set_members: bool | None = None,
    duplicate_scoring_strategy: str | None = None,
    optical_duplicate_pixel_distance: int | None = None,
    read_name_regex: str | None = None,
    barcode_tag: str | None = None,
    read_one_barcode_tag: str | None = None,
    read_two_barcode_tag: str | None = None,
    molecular_identifier_tag: str | None = None,
    assume_sort_order: str | None = None,
    add_pg_tag_to_reads: bool | None = None,
    reference_sequence: str | None = None,
    validation_stringency: str | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag duplicate molecules in coordinate-sorted, mate-cigar-tagged reads
    (``SimpleMarkDuplicatesWithMateCigar`` or ``UmiAwareMarkDuplicatesWithMateCigar``).

    Reads need ``MC`` mate-cigar tags. ``umi_aware=False`` runs
    ``SimpleMarkDuplicatesWithMateCigar`` (Picard marks it EXPERIMENTAL); ``umi_aware=True`` runs
    ``UmiAwareMarkDuplicatesWithMateCigar``, which also groups reads by their UMI (tag
    ``umi_tag_name``, ``RX`` by default; reads within ``max_edit_distance_to_join`` edits are
    joined) and needs ``umi_metrics_file``. The duplication metrics are parsed in ``metrics`` (and
    ``umi_metrics``). ``tagging_policy`` is ``DontTag``, ``OpticalOnly`` or ``All``;
    ``duplicate_scoring_strategy`` is ``SUM_OF_BASE_QUALITIES`` (default),
    ``TOTAL_MAPPED_REFERENCE_LENGTH`` or ``RANDOM``.
    """
    if not isinstance(umi_aware, bool):
        raise ConfigurationError("umi_aware must be True or False.")
    sources = _p.existing_files(input_paths, "input_paths")
    target = _p.new_file(output_path, "output_path")
    metrics = _p.new_file(metrics_file, "metrics_file")
    umi_metrics = _p.optional_new_file(umi_metrics_file, "umi_metrics_file")
    if umi_aware and umi_metrics is None:
        raise ConfigurationError("umi_metrics_file is required when umi_aware=True.")
    if not umi_aware and any(
        value is not None
        for value in (umi_metrics, umi_tag_name, max_edit_distance_to_join, allow_missing_umis)
    ):
        raise ConfigurationError("The UMI options only apply with umi_aware=True.")
    args = [f"INPUT={path}" for path in sources] + [f"OUTPUT={target}", f"METRICS_FILE={metrics}"]
    if umi_aware:
        args += _p.kv("UMI_METRICS_FILE", umi_metrics)
        args += _p.a_text("UMI_TAG_NAME", umi_tag_name, "umi_tag_name")
        args += _p.a_int(
            "MAX_EDIT_DISTANCE_TO_JOIN",
            max_edit_distance_to_join,
            "max_edit_distance_to_join",
            minimum=0,
        )
        args += _p.a_bool("ALLOW_MISSING_UMIS", allow_missing_umis, "allow_missing_umis")
    args += _p.a_bool("REMOVE_DUPLICATES", remove_duplicates, "remove_duplicates")
    args += _p.a_bool(
        "REMOVE_SEQUENCING_DUPLICATES", remove_sequencing_duplicates, "remove_sequencing_duplicates"
    )
    args += _p.a_choice("TAGGING_POLICY", tagging_policy, _TAGGING_POLICIES, "tagging_policy")
    args += _p.a_bool(
        "TAG_DUPLICATE_SET_MEMBERS", tag_duplicate_set_members, "tag_duplicate_set_members"
    )
    args += _p.a_choice(
        "DUPLICATE_SCORING_STRATEGY",
        duplicate_scoring_strategy,
        _SCORING,
        "duplicate_scoring_strategy",
    )
    args += _p.a_int(
        "OPTICAL_DUPLICATE_PIXEL_DISTANCE",
        optical_duplicate_pixel_distance,
        "optical_duplicate_pixel_distance",
        minimum=0,
    )
    args += _p.a_text("READ_NAME_REGEX", read_name_regex, "read_name_regex")
    args += _p.a_text("BARCODE_TAG", barcode_tag, "barcode_tag")
    args += _p.a_text("READ_ONE_BARCODE_TAG", read_one_barcode_tag, "read_one_barcode_tag")
    args += _p.a_text("READ_TWO_BARCODE_TAG", read_two_barcode_tag, "read_two_barcode_tag")
    args += _p.a_text(
        "MOLECULAR_IDENTIFIER_TAG", molecular_identifier_tag, "molecular_identifier_tag"
    )
    args += _p.a_choice("ASSUME_SORT_ORDER", assume_sort_order, _SORT_ORDERS, "assume_sort_order")
    args += _p.a_bool("ADD_PG_TAG_TO_READS", add_pg_tag_to_reads, "add_pg_tag_to_reads")
    args += _sam_common(validation_stringency, reference_sequence)
    outputs = {"output": target, "metrics_file": metrics}
    if umi_metrics is not None:
        outputs["umi_metrics_file"] = umi_metrics
    tool = (
        "UmiAwareMarkDuplicatesWithMateCigar" if umi_aware else "SimpleMarkDuplicatesWithMateCigar"
    )
    result = _p.run(
        tool, args, outputs=outputs, picard_path=picard_path, timeout_seconds=timeout_seconds
    )
    _metrics(result, metrics)
    if umi_metrics is not None:
        _metrics(result, umi_metrics, "umi_metrics")
    return result


# ---------------------------------------------------------------------------
# Illumina run folders
# ---------------------------------------------------------------------------


def _basecalls_args(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    run_barcode: str,
) -> list[str]:
    directory = _p.existing_directory(basecalls_dir, "basecalls_dir")
    if isinstance(lanes, (str, bytes)) or not lanes:
        raise ConfigurationError("lanes must be a non-empty list of lane numbers.")
    args = [f"BASECALLS_DIR={directory}"]
    args += _p.kv_many("LANE", [_p.integer(lane, "lanes", minimum=1) for lane in lanes])
    args += _p.kv("READ_STRUCTURE", _p.text(read_structure, "read_structure"))
    args += _p.kv("RUN_BARCODE", _p.text(run_barcode, "run_barcode"))
    return args


def picard_illumina_basecalls_to_fastq(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    run_barcode: str,
    *,
    output_prefix: str | None = None,
    multiplex_params: str | None = None,
    machine_name: str | None = None,
    flowcell_barcode: str | None = None,
    read_name_format: str | None = None,
    compress_outputs: bool | None = None,
    include_non_pf_reads: bool | None = None,
    apply_eamss_filter: bool | None = None,
    adapters_to_check: Sequence[str] | None = None,
    five_prime_adapter: str | None = None,
    three_prime_adapter: str | None = None,
    trimming_quality: int | None = None,
    min_trimmed_length: int | None = None,
    minimum_quality: int | None = None,
    max_mismatches: int | None = None,
    min_mismatch_delta: int | None = None,
    max_no_calls: int | None = None,
    minimum_base_quality: int | None = None,
    distance_mode: str | None = None,
    match_barcodes_inline: bool | None = None,
    ignore_unexpected_barcodes: bool | None = None,
    barcodes_dir: str | None = None,
    metrics_file: str | None = None,
    first_tile: int | None = None,
    tile_limit: int | None = None,
    sort: bool | None = None,
    num_processors: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate FASTQ files from the raw base calls of an Illumina run folder
    (``IlluminaBasecallsToFastq``).

    ``basecalls_dir`` is ``Data/Intensities/BaseCalls`` of the run, ``lanes`` the lane numbers and
    ``read_structure`` the description of the clusters (for example ``76T8B76T``). Give
    ``output_prefix`` for a run without demultiplexing (files ``<prefix>.1.fastq``,
    ``<prefix>.2.fastq``, ``<prefix>.barcode_1.fastq`` ...) or ``multiplex_params``, a tab-separated
    file (``OUTPUT_PREFIX``, ``BARCODE_1`` ...) that demultiplexes by barcode. Read names are in
    Casava 1.8 style unless ``read_name_format='ILLUMINA'``; that style needs ``machine_name`` and
    ``flowcell_barcode``. The FASTQ files that exist afterwards are listed in ``files``.
    """
    args = _basecalls_args(basecalls_dir, lanes, read_structure, run_barcode)
    if (output_prefix is None) == (multiplex_params is None):
        raise ConfigurationError("Give exactly one of output_prefix and multiplex_params.")
    prefix: Path | None = None
    watched: Path
    if output_prefix is not None:
        prefix = Path(output_prefix).expanduser().resolve()
        if not prefix.parent.is_dir():
            raise ConfigurationError(
                f"The folder for output_prefix does not exist: {prefix.parent}"
            )
        if _files_like(prefix.parent, f"{prefix.name}*.fastq*"):
            raise ConfigurationError(
                f"output_prefix would overwrite existing FASTQ files: {prefix}*"
            )
        args += _p.kv("OUTPUT_PREFIX", prefix)
        watched = prefix.parent
    else:
        params = _p.existing_file(str(multiplex_params), "multiplex_params")
        args += _p.kv("MULTIPLEX_PARAMS", params)
        watched = params.parent
    metrics = _p.optional_new_file(metrics_file, "metrics_file")
    args += _p.a_text("MACHINE_NAME", machine_name, "machine_name")
    args += _p.a_text("FLOWCELL_BARCODE", flowcell_barcode, "flowcell_barcode")
    args += _p.a_choice(
        "READ_NAME_FORMAT", read_name_format, _READ_NAME_FORMATS, "read_name_format"
    )
    args += _p.a_bool("COMPRESS_OUTPUTS", compress_outputs, "compress_outputs")
    args += _p.a_bool("INCLUDE_NON_PF_READS", include_non_pf_reads, "include_non_pf_reads")
    args += _p.a_bool("APPLY_EAMSS_FILTER", apply_eamss_filter, "apply_eamss_filter")
    args += _p.a_choices(
        "ADAPTERS_TO_CHECK", adapters_to_check, _ADAPTER_PAIRS, "adapters_to_check"
    )
    args += _p.a_text("FIVE_PRIME_ADAPTER", five_prime_adapter, "five_prime_adapter")
    args += _p.a_text("THREE_PRIME_ADAPTER", three_prime_adapter, "three_prime_adapter")
    args += _p.a_int("TRIMMING_QUALITY", trimming_quality, "trimming_quality", minimum=0)
    args += _p.a_int("MIN_TRIMMED_LENGTH", min_trimmed_length, "min_trimmed_length", minimum=0)
    args += _p.a_int("MINIMUM_QUALITY", minimum_quality, "minimum_quality", minimum=0)
    args += _p.a_int("MAX_MISMATCHES", max_mismatches, "max_mismatches", minimum=0)
    args += _p.a_int("MIN_MISMATCH_DELTA", min_mismatch_delta, "min_mismatch_delta", minimum=0)
    args += _p.a_int("MAX_NO_CALLS", max_no_calls, "max_no_calls", minimum=0)
    args += _p.a_int(
        "MINIMUM_BASE_QUALITY", minimum_base_quality, "minimum_base_quality", minimum=0
    )
    args += _p.a_choice("DISTANCE_MODE", distance_mode, _DISTANCE_MODES, "distance_mode")
    args += _p.a_bool("MATCH_BARCODES_INLINE", match_barcodes_inline, "match_barcodes_inline")
    args += _p.a_bool(
        "IGNORE_UNEXPECTED_BARCODES", ignore_unexpected_barcodes, "ignore_unexpected_barcodes"
    )
    if barcodes_dir is not None:
        args += _p.kv("BARCODES_DIR", _p.existing_directory(barcodes_dir, "barcodes_dir"))
    args += _p.kv("METRICS_FILE", metrics)
    args += _p.a_int("FIRST_TILE", first_tile, "first_tile", minimum=1)
    args += _p.a_int("TILE_LIMIT", tile_limit, "tile_limit", minimum=1)
    args += _p.a_bool("SORT", sort, "sort")
    args += _p.a_int("NUM_PROCESSORS", num_processors, "num_processors")
    outputs = {"metrics_file": metrics} if metrics is not None else {}
    before = {path for path in _files_like(watched, "*.fastq*")}
    result = _p.run(
        "IlluminaBasecallsToFastq",
        args,
        outputs=outputs,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        cwd=watched,
    )
    created = [path for path in _files_like(watched, "*.fastq*") if path not in before]
    result["files"] = [str(path) for path in created]
    if result["status"] == "completed" and not created:
        result["status"] = "failed"
    return result


def picard_illumina_basecalls_to_sam(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    run_barcode: str,
    sequencing_center: str,
    *,
    output_path: str | None = None,
    sample_alias: str | None = None,
    library_name: str | None = None,
    library_params: str | None = None,
    platform: str | None = None,
    read_group_id: str | None = None,
    run_start_date: str | None = None,
    include_non_pf_reads: bool | None = None,
    apply_eamss_filter: bool | None = None,
    adapters_to_check: Sequence[str] | None = None,
    five_prime_adapter: str | None = None,
    three_prime_adapter: str | None = None,
    minimum_quality: int | None = None,
    max_mismatches: int | None = None,
    min_mismatch_delta: int | None = None,
    max_no_calls: int | None = None,
    minimum_base_quality: int | None = None,
    distance_mode: str | None = None,
    barcode_population_strategy: str | None = None,
    include_barcode_quality: bool | None = None,
    include_bc_in_rg_tag: bool | None = None,
    match_barcodes_inline: bool | None = None,
    ignore_unexpected_barcodes: bool | None = None,
    barcodes_dir: str | None = None,
    molecular_index_tag: str | None = None,
    molecular_index_base_quality_tag: str | None = None,
    tag_per_molecular_index: Sequence[str] | None = None,
    metrics_file: str | None = None,
    first_tile: int | None = None,
    process_single_tile: int | None = None,
    tile_limit: int | None = None,
    sort: bool | None = None,
    num_processors: int | None = None,
    picard_path: str | None = None,
    timeout_seconds: float = _p.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transform raw Illumina base calls into unmapped SAM/BAM/CRAM files, demultiplexed by barcode
    (``IlluminaBasecallsToSam``).

    Single-sample run: give ``output_path`` (a new file whose extension picks the format),
    ``sample_alias`` and ``library_name``. Multi-sample run: give ``library_params``, a
    tab-separated file with the columns ``OUTPUT``, ``SAMPLE_ALIAS``, ``LIBRARY_NAME`` and
    ``BARCODE_1`` ...; the output files it names are listed in ``files``. ``read_structure``
    describes the clusters (for example ``76T8B76T``); ``sequencing_center`` fills the ``CN`` tag of
    the read group; ``molecular_index_tag`` is ``RX`` by default; ``run_start_date`` is a date such
    as ``2024-05-01``.
    """
    args = _basecalls_args(basecalls_dir, lanes, read_structure, run_barcode)
    args += _p.kv("SEQUENCING_CENTER", _p.text(sequencing_center, "sequencing_center"))
    target = _p.optional_new_file(output_path, "output_path")
    expected: dict[str, Path] = {}
    if library_params is None:
        if target is None or sample_alias is None:
            raise ConfigurationError("Give output_path and sample_alias, or library_params.")
        args += [f"OUTPUT={target}", f"SAMPLE_ALIAS={_p.text(sample_alias, 'sample_alias')}"]
        args += _p.a_text("LIBRARY_NAME", library_name, "library_name")
        expected["output"] = target
        cwd = target.parent
    else:
        if target is not None or sample_alias is not None or library_name is not None:
            raise ConfigurationError(
                "library_params replaces output_path, sample_alias and library_name."
            )
        params = _p.existing_file(library_params, "library_params")
        args += _p.kv("LIBRARY_PARAMS", params)
        rows = [
            line.split("\t")
            for line in params.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if len(rows) < 2:
            raise ConfigurationError(
                "library_params needs a header line and at least one library row."
            )
        for index, row in enumerate(rows[1:]):
            output = Path(row[0]).expanduser()
            expected[f"output_{index + 1}"] = (
                output if output.is_absolute() else (params.parent / output)
            )
            if expected[f"output_{index + 1}"].exists():
                raise ConfigurationError(
                    "library_params would overwrite an existing file: "
                    f"{expected[f'output_{index + 1}']}"
                )
        cwd = params.parent
    metrics = _p.optional_new_file(metrics_file, "metrics_file")
    args += _p.a_text("PLATFORM", platform, "platform")
    args += _p.a_text("READ_GROUP_ID", read_group_id, "read_group_id")
    args += _p.a_text("RUN_START_DATE", run_start_date, "run_start_date")
    args += _p.a_bool("INCLUDE_NON_PF_READS", include_non_pf_reads, "include_non_pf_reads")
    args += _p.a_bool("APPLY_EAMSS_FILTER", apply_eamss_filter, "apply_eamss_filter")
    args += _p.a_choices(
        "ADAPTERS_TO_CHECK", adapters_to_check, _ADAPTER_PAIRS, "adapters_to_check"
    )
    args += _p.a_text("FIVE_PRIME_ADAPTER", five_prime_adapter, "five_prime_adapter")
    args += _p.a_text("THREE_PRIME_ADAPTER", three_prime_adapter, "three_prime_adapter")
    args += _p.a_int("MINIMUM_QUALITY", minimum_quality, "minimum_quality", minimum=0)
    args += _p.a_int("MAX_MISMATCHES", max_mismatches, "max_mismatches", minimum=0)
    args += _p.a_int("MIN_MISMATCH_DELTA", min_mismatch_delta, "min_mismatch_delta", minimum=0)
    args += _p.a_int("MAX_NO_CALLS", max_no_calls, "max_no_calls", minimum=0)
    args += _p.a_int(
        "MINIMUM_BASE_QUALITY", minimum_base_quality, "minimum_base_quality", minimum=0
    )
    args += _p.a_choice("DISTANCE_MODE", distance_mode, _DISTANCE_MODES, "distance_mode")
    args += _p.a_choice(
        "BARCODE_POPULATION_STRATEGY",
        barcode_population_strategy,
        _BARCODE_STRATEGIES,
        "barcode_population_strategy",
    )
    args += _p.a_bool("INCLUDE_BARCODE_QUALITY", include_barcode_quality, "include_barcode_quality")
    args += _p.a_bool("INCLUDE_BC_IN_RG_TAG", include_bc_in_rg_tag, "include_bc_in_rg_tag")
    args += _p.a_bool("MATCH_BARCODES_INLINE", match_barcodes_inline, "match_barcodes_inline")
    args += _p.a_bool(
        "IGNORE_UNEXPECTED_BARCODES", ignore_unexpected_barcodes, "ignore_unexpected_barcodes"
    )
    if barcodes_dir is not None:
        args += _p.kv("BARCODES_DIR", _p.existing_directory(barcodes_dir, "barcodes_dir"))
    args += _p.a_text("MOLECULAR_INDEX_TAG", molecular_index_tag, "molecular_index_tag")
    args += _p.a_text(
        "MOLECULAR_INDEX_BASE_QUALITY_TAG",
        molecular_index_base_quality_tag,
        "molecular_index_base_quality_tag",
    )
    args += _p.a_texts(
        "TAG_PER_MOLECULAR_INDEX", tag_per_molecular_index, "tag_per_molecular_index"
    )
    args += _p.kv("METRICS_FILE", metrics)
    args += _p.a_int("FIRST_TILE", first_tile, "first_tile", minimum=1)
    args += _p.a_int("PROCESS_SINGLE_TILE", process_single_tile, "process_single_tile", minimum=1)
    args += _p.a_int("TILE_LIMIT", tile_limit, "tile_limit", minimum=1)
    args += _p.a_bool("SORT", sort, "sort")
    args += _p.a_int("NUM_PROCESSORS", num_processors, "num_processors")
    if metrics is not None:
        expected["metrics_file"] = metrics
    result = _p.run(
        "IlluminaBasecallsToSam",
        args,
        outputs=expected,
        picard_path=picard_path,
        timeout_seconds=timeout_seconds,
        cwd=cwd,
    )
    result["files"] = [
        str(path) for key, path in expected.items() if key != "metrics_file" and path.is_file()
    ]
    return result


__all__ = [
    "picard_downsample_sam",
    "picard_filter_sam_reads",
    "picard_illumina_basecalls_to_fastq",
    "picard_illumina_basecalls_to_sam",
    "picard_mark_duplicates_with_mate_cigar",
    "picard_mark_illumina_adapters",
    "picard_merge_bam_alignment",
    "picard_position_based_downsample_sam",
    "picard_reorder_sam",
    "picard_revert_sam",
    "picard_split_sam_by_number_of_reads",
]
_register_public_verbose(globals())
