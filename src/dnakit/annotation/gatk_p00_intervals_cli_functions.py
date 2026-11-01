"""GATK command-line facades, part 3: interval lists, genotyping arrays, flow reads, PathSeq.

Each public function runs one ``gatk <Tool>`` command of GATK 4.6.2.0 (Picard 3.4.0 is bundled in
the same launcher) with an explicit argument list: no shell, outputs only to new paths, bounded
time and output, and the version in every result. The user installs GATK and Java 17
(``conda install -c bioconda gatk4`` or the release zip from
https://github.com/broadinstitute/gatk/releases) and puts the ``gatk`` launcher on PATH, sets
``DNAKIT_GATK`` or passes ``gatk_path``. A function types the inputs and outputs of the tool and
the options that change its result; every other long option of the tool can be given through
``options`` and is checked against the list printed by ``gatk <Tool> --help``. Options that are
not given keep the GATK default. Results are JSON-serialisable dictionaries with ``status``,
``tool``, ``version``, ``command``, ``outputs`` (paths), small parsed results and a log tail.

This part covers Picard interval-list tools, the Illumina genotyping-array tools (GTC, manifests,
ADPC, array metrics), the flow-based read tools for Ultima-style data and the PathSeq reference
builders.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics import _gatk_p00_core as _core
from dnakit.exceptions import BackendExecutionError, ConfigurationError

__all__ = [
    "gatk_add_flow_base_quality",
    "gatk_add_flow_snv_quality",
    "gatk_bpm_to_normalization_manifest_csv",
    "gatk_calculate_average_combined_annotations",
    "gatk_collect_arrays_variant_calling_metrics",
    "gatk_combine_genotyping_array_vcfs",
    "gatk_compare_gtc_files",
    "gatk_compare_interval_lists",
    "gatk_convert_interval_list",
    "gatk_create_extended_illumina_manifest",
    "gatk_flow_feature_mapper",
    "gatk_flow_pair_hmm_align_reads_to_haplotypes",
    "gatk_ground_truth_reads_builder",
    "gatk_ground_truth_scorer",
    "gatk_gtc_to_vcf",
    "gatk_interval_list_tools",
    "gatk_lift_over_interval_list",
    "gatk_merge_ped_into_vcf",
    "gatk_pathseq_build_kmers",
    "gatk_pathseq_build_reference_taxonomy",
    "gatk_split_intervals",
    "gatk_vcf_to_adpc",
]


def gatk_convert_interval_list(
    input_path: str,
    output_path: str,
    output_format: str,
    *,
    sequence_dictionary: str | None = None,
    keep_length_zero_intervals: bool | None = None,
    unique: bool | None = None,
    score: int | None = None,
    sort: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert between BED and Picard interval-list formats: output_format='interval_list' runs
    BedToIntervalList, output_format='bed' runs IntervalListToBed.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): BED file (for 'interval_list') or Picard interval list (for
      'bed').
    - ``output_path`` (``--OUTPUT``): New file in the requested format.
    - ``output_format`` (``--None``): 'interval_list' converts BED to an interval list, 'bed'
      converts an interval list to BED.
    - ``sequence_dictionary`` (``--SEQUENCE_DICTIONARY``): Sequence dictionary (.dict) for the
      header of the interval list ('interval_list' only).
    - ``keep_length_zero_intervals`` (``--KEEP_LENGTH_ZERO_INTERVALS``): Keep zero-length BED
      intervals ('interval_list' only).
    - ``unique`` (``--UNIQUE``): Merge overlapping and adjacent intervals ('interval_list'
      only).
    - ``score`` (``--SCORE``): Value written to the BED score column ('bed' only).
    - ``sort`` (``--SORT``): Sort the output by coordinate.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk <tool> --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``output_format`` and ``interval_count`` (data lines of the
    output).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    _core.choice(output_format, ("interval_list", "bed"), "output_format")
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if sequence_dictionary is not None and output_format not in ("interval_list",):
        raise ConfigurationError(
            "sequence_dictionary only applies to output_format in ['interval_list']."
        )
    if sequence_dictionary is None and output_format in ("interval_list",):
        raise ConfigurationError(
            "sequence_dictionary is required when output_format is in ['interval_list']."
        )
    if sequence_dictionary is not None:
        sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
        arguments += ["--SEQUENCE_DICTIONARY", str(sequence_dictionary_path)]
    if keep_length_zero_intervals is not None and output_format not in ("interval_list",):
        raise ConfigurationError(
            "keep_length_zero_intervals only applies to output_format in ['interval_list']."
        )
    if keep_length_zero_intervals is not None:
        arguments += [
            "--KEEP_LENGTH_ZERO_INTERVALS",
            "true"
            if _core.flag(keep_length_zero_intervals, "keep_length_zero_intervals")
            else "false",
        ]
    if unique is not None and output_format not in ("interval_list",):
        raise ConfigurationError("unique only applies to output_format in ['interval_list'].")
    if unique is not None:
        arguments += ["--UNIQUE", "true" if _core.flag(unique, "unique") else "false"]
    if score is not None and output_format not in ("bed",):
        raise ConfigurationError("score only applies to output_format in ['bed'].")
    if score is not None:
        score_value = _core.integer(score, "score", minimum=0)
        arguments += ["--SCORE", str(score_value)]
    if sort is not None:
        arguments += ["--SORT", "true" if _core.flag(sort, "sort") else "false"]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "BedToIntervalList" if output_format == "interval_list" else "IntervalListToBed"
    reserved = _core.names(
        "INPUT KEEP_LENGTH_ZERO_INTERVALS OUTPUT REFERENCE_SEQUENCE SCORE "
        "SEQUENCE_DICTIONARY SORT UNIQUE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_convert_interval_list",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["interval_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_compare_interval_lists(
    first: str,
    second: str,
    reference: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two interval lists for equality and report the intervals that differ.

    Parameters (the GATK option is in parentheses):

    - ``first`` (``--L``): First interval list (file) or locus.
    - ``second`` (``--L2``): Second interval list (file) or locus.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CompareIntervalLists --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``identical`` (GATK exit code 0), ``has_findings`` and ``report_lines`` (the
    comparison lines of the log).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    first_text = _core.interval(first, "first")
    arguments += ["--L", first_text]
    second_text = _core.interval(second, "second")
    arguments += ["--L2", second_text]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    tool = "CompareIntervalLists"
    reserved = _core.names("L L2 reference")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_compare_interval_lists",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
        success_codes=(0, 1),
    )
    extras["has_findings"] = result["return_code"] != 0
    extras["identical"] = result["return_code"] == 0
    extras["report_lines"] = [
        line
        for line in result["_log"].splitlines()
        if line.startswith(("<", ">", "Only", "Intervals"))
        or "differ" in line.lower()
        or "identical" in line.lower()
    ][:100]
    return _core.finish(result, outputs=outs, **extras)


def gatk_interval_list_tools(
    input_paths: Sequence[str],
    *,
    output_path: str | None = None,
    output_directory: str | None = None,
    action: str | None = None,
    second_input_paths: Sequence[str] | None = None,
    padding: int | None = None,
    invert: bool | None = None,
    unique: bool | None = None,
    sort: bool | None = None,
    dont_merge_abutting: bool | None = None,
    break_bands_at_multiples_of: int | None = None,
    scatter_count: int | None = None,
    scatter_content: int | None = None,
    subdivision_mode: str | None = None,
    count_output: str | None = None,
    output_value: str | None = None,
    comments: Sequence[str] | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Manipulate Picard interval lists: concatenate, union, intersect, subtract, symmetric
    difference, pad, invert, merge, break bands and scatter.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--INPUT``): One or more interval lists, BED/VCF files or loci readable
      by Picard.
    - ``output_path`` (``--OUTPUT``): New interval list (when scatter_count is 1).
    - ``output_directory`` (``--OUTPUT``): New folder for the scattered lists (when
      scatter_count is above 1).
    - ``action`` (``--ACTION``): How the inputs are combined; GATK default CONCAT.
    - ``second_input_paths`` (``--SECOND_INPUT``): Second set of intervals for SUBTRACT, SYMDIFF
      and OVERLAPS.
    - ``padding`` (``--PADDING``): Bases added to both ends of each interval (negative shrinks).
    - ``invert`` (``--INVERT``): Produce the intervals that are not covered.
    - ``unique`` (``--UNIQUE``): Merge overlapping and adjacent intervals (implies sorting).
    - ``sort`` (``--SORT``): Sort the result by coordinate.
    - ``dont_merge_abutting`` (``--DONT_MERGE_ABUTTING``): Keep abutting intervals separate.
    - ``break_bands_at_multiples_of`` (``--BREAK_BANDS_AT_MULTIPLES_OF``): Break intervals at
      multiples of this many bases.
    - ``scatter_count`` (``--SCATTER_COUNT``): Number of lists to scatter the result into.
    - ``scatter_content`` (``--SCATTER_CONTENT``): Target number of bases per scattered list.
    - ``subdivision_mode`` (``--SUBDIVISION_MODE``): How the intervals are scattered.
    - ``count_output`` (``--COUNT_OUTPUT``): New file that receives the count requested with
      output_value.
    - ``output_value`` (``--OUTPUT_VALUE``): What to count into count_output or the log.
    - ``comments`` (``--COMMENT``): Comment lines added to the header of the output.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk IntervalListTools --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``interval_count`` (single output), ``scattered_lists`` (scatter
    output) and ``count`` (count_output).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if output_path is not None and output_directory is not None:
        raise ConfigurationError("Give output_path or output_directory, not both.")
    if output_path is None and output_directory is None and count_output is None:
        raise ConfigurationError("Give output_path, output_directory or count_output.")
    if output_directory is not None and (scatter_count is None or scatter_count < 2):
        raise ConfigurationError(
            "output_directory is for scatter_count of 2 or more; use output_path otherwise."
        )
    if output_path is not None and scatter_count is not None and scatter_count > 1:
        raise ConfigurationError("scatter_count above 1 needs output_directory, not output_path.")
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--INPUT", str(item)]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if output_directory is not None:
        output_directory_path = _core.new_directory(output_directory, "output_directory")
        arguments += ["--OUTPUT", str(output_directory_path)]
        made_directories.append(output_directory_path)
        outs["output_directory"] = output_directory_path
    if action is not None:
        action_value = _core.choice(
            action, ("CONCAT", "UNION", "INTERSECT", "SUBTRACT", "SYMDIFF", "OVERLAPS"), "action"
        )
        arguments += ["--ACTION", action_value]
    if second_input_paths is not None:
        second_input_paths_paths = _core.existing_files(second_input_paths, "second_input_paths")
        for item in second_input_paths_paths:
            arguments += ["--SECOND_INPUT", str(item)]
    if padding is not None:
        padding_value = _core.integer(padding, "padding")
        arguments += ["--PADDING", str(padding_value)]
    if invert is not None:
        arguments += ["--INVERT", "true" if _core.flag(invert, "invert") else "false"]
    if unique is not None:
        arguments += ["--UNIQUE", "true" if _core.flag(unique, "unique") else "false"]
    if sort is not None:
        arguments += ["--SORT", "true" if _core.flag(sort, "sort") else "false"]
    if dont_merge_abutting is not None:
        arguments += [
            "--DONT_MERGE_ABUTTING",
            "true" if _core.flag(dont_merge_abutting, "dont_merge_abutting") else "false",
        ]
    if break_bands_at_multiples_of is not None:
        break_bands_at_multiples_of_value = _core.integer(
            break_bands_at_multiples_of, "break_bands_at_multiples_of", minimum=0
        )
        arguments += ["--BREAK_BANDS_AT_MULTIPLES_OF", str(break_bands_at_multiples_of_value)]
    if scatter_count is not None:
        scatter_count_value = _core.integer(scatter_count, "scatter_count", minimum=1)
        arguments += ["--SCATTER_COUNT", str(scatter_count_value)]
    if scatter_content is not None:
        scatter_content_value = _core.integer(scatter_content, "scatter_content", minimum=1)
        arguments += ["--SCATTER_CONTENT", str(scatter_content_value)]
    if subdivision_mode is not None:
        subdivision_mode_value = _core.choice(
            subdivision_mode,
            (
                "INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
                "INTERVAL_COUNT",
                "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER",
            ),
            "subdivision_mode",
        )
        arguments += ["--SUBDIVISION_MODE", subdivision_mode_value]
    if count_output is not None:
        count_output_path = _core.new_file(count_output, "count_output")
        arguments += ["--COUNT_OUTPUT", str(count_output_path)]
        expected.append(count_output_path)
        outs["count_output"] = count_output_path
    if output_value is not None:
        output_value_value = _core.choice(
            output_value, ("NONE", "BASES", "INTERVALS"), "output_value"
        )
        arguments += ["--OUTPUT_VALUE", output_value_value]
    if comments is not None:
        comments_values = _core.texts(comments, "comments")
        for item in comments_values:
            arguments += ["--COMMENT", item]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "IntervalListTools"
    reserved = _core.names(
        "ACTION BREAK_BANDS_AT_MULTIPLES_OF COMMENT COUNT_OUTPUT "
        "DONT_MERGE_ABUTTING INPUT INVERT OUTPUT OUTPUT_VALUE PADDING "
        "REFERENCE_SEQUENCE SCATTER_CONTENT SCATTER_COUNT SECOND_INPUT SORT "
        "SUBDIVISION_MODE UNIQUE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_interval_list_tools",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs:
        extras["interval_count"] = _core.count_data_lines(outs["output_path"])
    if "output_directory" in outs:
        extras["scattered_lists"] = sorted(
            str(item) for item in outs["output_directory"].rglob("*.interval_list")
        )
    if "count_output" in outs and outs["count_output"].is_file():
        extras["count"] = _core.first_integer(outs["count_output"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_lift_over_interval_list(
    input_path: str,
    chain: str,
    sequence_dictionary: str,
    output_path: str,
    *,
    reject: str | None = None,
    min_liftover_pct: float | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift over an interval list from one reference build to another with a chain file; intervals
    that cannot be lifted are rejected.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Interval list on the old reference.
    - ``chain`` (``--CHAIN``): Chain file from the old to the new reference.
    - ``sequence_dictionary`` (``--SEQUENCE_DICTIONARY``): Sequence dictionary of the new
      reference.
    - ``output_path`` (``--OUTPUT``): New lifted interval list.
    - ``reject`` (``--REJECT``): New interval list of intervals that could not be lifted.
    - ``min_liftover_pct`` (``--MIN_LIFTOVER_PCT``): Minimum fraction of an interval that must
      lift over.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk LiftOverIntervalList --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``interval_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    chain_path = _core.existing_file(chain, "chain")
    arguments += ["--CHAIN", str(chain_path)]
    sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
    arguments += ["--SEQUENCE_DICTIONARY", str(sequence_dictionary_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reject is not None:
        reject_path = _core.new_file(reject, "reject")
        arguments += ["--REJECT", str(reject_path)]
        expected.append(reject_path)
        outs["reject"] = reject_path
    if min_liftover_pct is not None:
        min_liftover_pct_value = _core.number(
            min_liftover_pct, "min_liftover_pct", minimum=0.0, maximum=1.0
        )
        arguments += ["--MIN_LIFTOVER_PCT", str(min_liftover_pct_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "LiftOverIntervalList"
    reserved = _core.names(
        "CHAIN INPUT MIN_LIFTOVER_PCT OUTPUT REFERENCE_SEQUENCE REJECT SEQUENCE_DICTIONARY"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_lift_over_interval_list",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["interval_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_split_intervals(
    output_dir: str,
    *,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    scatter_count: int | None = None,
    subdivision_mode: str | None = None,
    dont_mix_contigs: bool | None = None,
    min_contig_size: int | None = None,
    interval_file_prefix: str | None = None,
    extension: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split an interval list (or a whole reference) into sub-interval files for scatter-gather
    processing.

    Parameters (the GATK option is in parentheses):

    - ``output_dir`` (``--output``): New folder for the scattered interval lists.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``scatter_count`` (``--scatter-count``): Number of files to write.
    - ``subdivision_mode`` (``--subdivision-mode``): How the intervals are split.
    - ``dont_mix_contigs`` (``--dont-mix-contigs``): Keep each contig in one file.
    - ``min_contig_size`` (``--min-contig-size``): Contigs shorter than this are ignored.
    - ``interval_file_prefix`` (``--interval-file-prefix``): Prefix of the interval file names.
    - ``extension`` (``--extension``): Extension of the interval file names.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk SplitIntervals --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``interval_files``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if scatter_count is not None:
        scatter_count_value = _core.integer(scatter_count, "scatter_count", minimum=1)
        arguments += ["--scatter-count", str(scatter_count_value)]
    if subdivision_mode is not None:
        subdivision_mode_value = _core.choice(
            subdivision_mode,
            (
                "INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
                "INTERVAL_COUNT",
                "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER",
            ),
            "subdivision_mode",
        )
        arguments += ["--subdivision-mode", subdivision_mode_value]
    if dont_mix_contigs is not None:
        arguments += [
            "--dont-mix-contigs",
            "true" if _core.flag(dont_mix_contigs, "dont_mix_contigs") else "false",
        ]
    if min_contig_size is not None:
        min_contig_size_value = _core.integer(min_contig_size, "min_contig_size", minimum=0)
        arguments += ["--min-contig-size", str(min_contig_size_value)]
    if interval_file_prefix is not None:
        interval_file_prefix_text = _core.text(interval_file_prefix, "interval_file_prefix")
        arguments += ["--interval-file-prefix", interval_file_prefix_text]
    if extension is not None:
        extension_text = _core.text(extension, "extension")
        arguments += ["--extension", extension_text]
    tool = "SplitIntervals"
    reserved = _core.names(
        "dont-mix-contigs extension interval-file-prefix intervals "
        "min-contig-size output reference scatter-count subdivision-mode"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_split_intervals",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["interval_files"] = sorted(str(item) for item in output_dir_path.glob("*"))
    return _core.finish(result, outputs=outs, **extras)


def gatk_bpm_to_normalization_manifest_csv(
    input_path: str,
    cluster_file: str,
    output_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an Illumina bead pool manifest (.bpm) and cluster file into a bpm.csv normalization
    manifest.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Illumina bead pool manifest (.bpm).
    - ``cluster_file`` (``--CLUSTER_FILE``): Illumina cluster file (.egt).
    - ``output_path`` (``--OUTPUT``): New bpm.csv file.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk BpmToNormalizationManifestCsv --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``row_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    cluster_file_path = _core.existing_file(cluster_file, "cluster_file")
    arguments += ["--CLUSTER_FILE", str(cluster_file_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "BpmToNormalizationManifestCsv"
    reserved = _core.names("CLUSTER_FILE INPUT OUTPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_bpm_to_normalization_manifest_csv",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["row_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_combine_genotyping_array_vcfs(
    input_paths: Sequence[str],
    output_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine several genotyping-array VCF files into one VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--INPUT``): Genotyping-array VCF files.
    - ``output_path`` (``--OUTPUT``): New combined VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CombineGenotypingArrayVcfs --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--INPUT", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CombineGenotypingArrayVcfs"
    reserved = _core.names("INPUT OUTPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_combine_genotyping_array_vcfs",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["record_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_compare_gtc_files(
    input_paths: Sequence[str],
    illumina_bead_pool_manifest_file: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two Illumina GTC files (genotypes and intensities) against a bead pool manifest.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--INPUT``): Exactly two GTC files.
    - ``illumina_bead_pool_manifest_file`` (``--ILLUMINA_BEAD_POOL_MANIFEST_FILE``): Illumina
      bead pool manifest (.bpm).
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CompareGtcFiles --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``files_match`` (GATK exit code 0) and ``has_findings``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if isinstance(input_paths, (str, bytes)) or len(input_paths) != 2:
        raise ConfigurationError("input_paths must hold exactly two GTC files.")
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--INPUT", str(item)]
    illumina_bead_pool_manifest_file_path = _core.existing_file(
        illumina_bead_pool_manifest_file, "illumina_bead_pool_manifest_file"
    )
    arguments += ["--ILLUMINA_BEAD_POOL_MANIFEST_FILE", str(illumina_bead_pool_manifest_file_path)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CompareGtcFiles"
    reserved = _core.names("ILLUMINA_BEAD_POOL_MANIFEST_FILE INPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_compare_gtc_files",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
        success_codes=(0, 1),
    )
    extras["has_findings"] = result["return_code"] != 0
    extras["files_match"] = result["return_code"] == 0
    return _core.finish(result, outputs=outs, **extras)


def gatk_create_extended_illumina_manifest(
    input_path: str,
    output_path: str,
    reference_sequence: str,
    report_file: str,
    *,
    bad_assays_file: str | None = None,
    cluster_file: str | None = None,
    dbsnp_file: str | None = None,
    target_build: str | None = None,
    supported_builds: Sequence[str] | None = None,
    supported_reference_files: Sequence[str] | None = None,
    supported_chain_files: Sequence[str] | None = None,
    flag_duplicates: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create the Extended Illumina Manifest (reference-specific, with mapped positions and
    alleles) that GtcToVcf needs.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Text version of the Illumina .bpm manifest.
    - ``output_path`` (``--OUTPUT``): New extended manifest.
    - ``reference_sequence`` (``--REFERENCE_SEQUENCE``): Reference FASTA of the target build
      (with .fai and .dict).
    - ``report_file`` (``--REPORT_FILE``): New report file.
    - ``bad_assays_file`` (``--BAD_ASSAYS_FILE``): New file with the unmappable assays.
    - ``cluster_file`` (``--CLUSTER_FILE``): Illumina cluster file (.egt).
    - ``dbsnp_file`` (``--DBSNP_FILE``): dbSNP VCF used to flag known variants.
    - ``target_build`` (``--TARGET_BUILD``): Genome build of the target reference, such as 37.
    - ``supported_builds`` (``--SUPPORTED_BUILD``): Other builds present in the manifest (same
      order as the reference and chain files).
    - ``supported_reference_files`` (``--SUPPORTED_REFERENCE_FILE``): Reference FASTA files of
      the supported builds.
    - ``supported_chain_files`` (``--SUPPORTED_CHAIN_FILE``): Chain files from the supported
      builds to the target build.
    - ``flag_duplicates`` (``--FLAG_DUPLICATES``): Flag assays that map to the same position.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CreateExtendedIlluminaManifest --help`` are accepted,
      a list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``row_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    reference_sequence_path = _core.existing_file(reference_sequence, "reference_sequence")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_sequence_path)]
    report_file_path = _core.new_file(report_file, "report_file")
    arguments += ["--REPORT_FILE", str(report_file_path)]
    expected.append(report_file_path)
    outs["report_file"] = report_file_path
    if bad_assays_file is not None:
        bad_assays_file_path = _core.new_file(bad_assays_file, "bad_assays_file")
        arguments += ["--BAD_ASSAYS_FILE", str(bad_assays_file_path)]
        expected.append(bad_assays_file_path)
        outs["bad_assays_file"] = bad_assays_file_path
    if cluster_file is not None:
        cluster_file_path = _core.existing_file(cluster_file, "cluster_file")
        arguments += ["--CLUSTER_FILE", str(cluster_file_path)]
    if dbsnp_file is not None:
        dbsnp_file_path = _core.existing_file(dbsnp_file, "dbsnp_file")
        arguments += ["--DBSNP_FILE", str(dbsnp_file_path)]
    if target_build is not None:
        target_build_text = _core.text(target_build, "target_build")
        arguments += ["--TARGET_BUILD", target_build_text]
    if supported_builds is not None:
        supported_builds_values = _core.texts(supported_builds, "supported_builds")
        for item in supported_builds_values:
            arguments += ["--SUPPORTED_BUILD", item]
    if supported_reference_files is not None:
        supported_reference_files_paths = _core.existing_files(
            supported_reference_files, "supported_reference_files"
        )
        for item in supported_reference_files_paths:
            arguments += ["--SUPPORTED_REFERENCE_FILE", str(item)]
    if supported_chain_files is not None:
        supported_chain_files_paths = _core.existing_files(
            supported_chain_files, "supported_chain_files"
        )
        for item in supported_chain_files_paths:
            arguments += ["--SUPPORTED_CHAIN_FILE", str(item)]
    if flag_duplicates is not None:
        arguments += [
            "--FLAG_DUPLICATES",
            "true" if _core.flag(flag_duplicates, "flag_duplicates") else "false",
        ]
    tool = "CreateExtendedIlluminaManifest"
    reserved = _core.names(
        "BAD_ASSAYS_FILE CLUSTER_FILE DBSNP_FILE FLAG_DUPLICATES INPUT OUTPUT "
        "REFERENCE_SEQUENCE REPORT_FILE SUPPORTED_BUILD SUPPORTED_CHAIN_FILE "
        "SUPPORTED_REFERENCE_FILE TARGET_BUILD"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_create_extended_illumina_manifest",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["row_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_gtc_to_vcf(
    input_path: str,
    output_path: str,
    extended_illumina_manifest: str,
    illumina_bead_pool_manifest_file: str,
    cluster_file: str,
    sample_alias: str,
    reference: str,
    *,
    expected_gender: str | None = None,
    gender_gtc: str | None = None,
    fingerprint_genotypes_vcf_file: str | None = None,
    analysis_version_number: int | None = None,
    pipeline_version: str | None = None,
    do_not_allow_calls_on_zeroed_out_assays: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert an Illumina GTC genotype file to a VCF with an extended manifest, bead pool manifest
    and cluster file.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Illumina GTC file.
    - ``output_path`` (``--OUTPUT``): New VCF.
    - ``extended_illumina_manifest`` (``--EXTENDED_ILLUMINA_MANIFEST``): Extended Illumina
      manifest (see gatk_create_extended_illumina_manifest).
    - ``illumina_bead_pool_manifest_file`` (``--ILLUMINA_BEAD_POOL_MANIFEST_FILE``): Illumina
      bead pool manifest (.bpm).
    - ``cluster_file`` (``--CLUSTER_FILE``): Illumina cluster file (.egt).
    - ``sample_alias`` (``--SAMPLE_ALIAS``): Sample name written to the VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``expected_gender`` (``--EXPECTED_GENDER``): Expected gender of the sample.
    - ``gender_gtc`` (``--GENDER_GTC``): GTC called with a gender-optimised cluster file.
    - ``fingerprint_genotypes_vcf_file`` (``--FINGERPRINT_GENOTYPES_VCF_FILE``): Fingerprint
      genotypes VCF used for the fingerprint check.
    - ``analysis_version_number`` (``--ANALYSIS_VERSION_NUMBER``): Analysis version number
      written to the VCF.
    - ``pipeline_version`` (``--PIPELINE_VERSION``): Pipeline version written to the VCF.
    - ``do_not_allow_calls_on_zeroed_out_assays``
      (``--DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS``): Do not call zeroed-out assays.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk GtcToVcf --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    extended_illumina_manifest_path = _core.existing_file(
        extended_illumina_manifest, "extended_illumina_manifest"
    )
    arguments += ["--EXTENDED_ILLUMINA_MANIFEST", str(extended_illumina_manifest_path)]
    illumina_bead_pool_manifest_file_path = _core.existing_file(
        illumina_bead_pool_manifest_file, "illumina_bead_pool_manifest_file"
    )
    arguments += ["--ILLUMINA_BEAD_POOL_MANIFEST_FILE", str(illumina_bead_pool_manifest_file_path)]
    cluster_file_path = _core.existing_file(cluster_file, "cluster_file")
    arguments += ["--CLUSTER_FILE", str(cluster_file_path)]
    sample_alias_text = _core.text(sample_alias, "sample_alias")
    arguments += ["--SAMPLE_ALIAS", sample_alias_text]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if expected_gender is not None:
        expected_gender_text = _core.text(expected_gender, "expected_gender")
        arguments += ["--EXPECTED_GENDER", expected_gender_text]
    if gender_gtc is not None:
        gender_gtc_path = _core.existing_file(gender_gtc, "gender_gtc")
        arguments += ["--GENDER_GTC", str(gender_gtc_path)]
    if fingerprint_genotypes_vcf_file is not None:
        fingerprint_genotypes_vcf_file_path = _core.existing_file(
            fingerprint_genotypes_vcf_file, "fingerprint_genotypes_vcf_file"
        )
        arguments += ["--FINGERPRINT_GENOTYPES_VCF_FILE", str(fingerprint_genotypes_vcf_file_path)]
    if analysis_version_number is not None:
        analysis_version_number_value = _core.integer(
            analysis_version_number, "analysis_version_number", minimum=0
        )
        arguments += ["--ANALYSIS_VERSION_NUMBER", str(analysis_version_number_value)]
    if pipeline_version is not None:
        pipeline_version_text = _core.text(pipeline_version, "pipeline_version")
        arguments += ["--PIPELINE_VERSION", pipeline_version_text]
    if do_not_allow_calls_on_zeroed_out_assays is not None:
        arguments += [
            "--DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS",
            "true"
            if _core.flag(
                do_not_allow_calls_on_zeroed_out_assays, "do_not_allow_calls_on_zeroed_out_assays"
            )
            else "false",
        ]
    tool = "GtcToVcf"
    reserved = _core.names(
        "ANALYSIS_VERSION_NUMBER CLUSTER_FILE "
        "DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS EXPECTED_GENDER "
        "EXTENDED_ILLUMINA_MANIFEST FINGERPRINT_GENOTYPES_VCF_FILE GENDER_GTC "
        "ILLUMINA_BEAD_POOL_MANIFEST_FILE INPUT OUTPUT PIPELINE_VERSION "
        "REFERENCE_SEQUENCE SAMPLE_ALIAS"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_gtc_to_vcf",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["record_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_merge_ped_into_vcf(
    original_vcf: str,
    ped_file: str,
    map_file: str,
    zcall_thresholds_file: str,
    zcall_version: str,
    output_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge a single-sample PED file from zCall into a single-sample VCF.

    Parameters (the GATK option is in parentheses):

    - ``original_vcf`` (``--ORIGINAL_VCF``): Single-sample VCF made by GtcToVcf.
    - ``ped_file`` (``--PED_FILE``): PED file from zCall.
    - ``map_file`` (``--MAP_FILE``): MAP file of the PED file.
    - ``zcall_thresholds_file`` (``--ZCALL_THRESHOLDS_FILE``): zCall thresholds file.
    - ``zcall_version`` (``--ZCALL_VERSION``): zCall version recorded in the VCF header.
    - ``output_path`` (``--OUTPUT``): New merged VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk MergePedIntoVcf --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    original_vcf_path = _core.existing_file(original_vcf, "original_vcf")
    arguments += ["--ORIGINAL_VCF", str(original_vcf_path)]
    ped_file_path = _core.existing_file(ped_file, "ped_file")
    arguments += ["--PED_FILE", str(ped_file_path)]
    map_file_path = _core.existing_file(map_file, "map_file")
    arguments += ["--MAP_FILE", str(map_file_path)]
    zcall_thresholds_file_path = _core.existing_file(zcall_thresholds_file, "zcall_thresholds_file")
    arguments += ["--ZCALL_THRESHOLDS_FILE", str(zcall_thresholds_file_path)]
    zcall_version_text = _core.text(zcall_version, "zcall_version")
    arguments += ["--ZCALL_VERSION", zcall_version_text]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "MergePedIntoVcf"
    reserved = _core.names(
        "MAP_FILE ORIGINAL_VCF OUTPUT PED_FILE REFERENCE_SEQUENCE "
        "ZCALL_THRESHOLDS_FILE ZCALL_VERSION"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_merge_ped_into_vcf",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["record_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_vcf_to_adpc(
    vcf_paths: Sequence[str],
    output_path: str,
    samples_file: str,
    num_markers_file: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert genotyping-array VCF files to an ADPC (adpc.bin) file for verifyIDintensity, with
    sample and marker-count files.

    Parameters (the GATK option is in parentheses):

    - ``vcf_paths`` (``--VCF``): Array VCF files with intensity data.
    - ``output_path`` (``--OUTPUT``): New adpc.bin file.
    - ``samples_file`` (``--SAMPLES_FILE``): New text file with the sample names (same order as
      adpc.bin).
    - ``num_markers_file`` (``--NUM_MARKERS_FILE``): New text file with the number of loci.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk VcfToAdpc --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``num_markers`` and ``samples``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    vcf_paths_paths = _core.existing_files(vcf_paths, "vcf_paths")
    for item in vcf_paths_paths:
        arguments += ["--VCF", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    samples_file_path = _core.new_file(samples_file, "samples_file")
    arguments += ["--SAMPLES_FILE", str(samples_file_path)]
    expected.append(samples_file_path)
    outs["samples_file"] = samples_file_path
    num_markers_file_path = _core.new_file(num_markers_file, "num_markers_file")
    arguments += ["--NUM_MARKERS_FILE", str(num_markers_file_path)]
    expected.append(num_markers_file_path)
    outs["num_markers_file"] = num_markers_file_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "VcfToAdpc"
    reserved = _core.names("NUM_MARKERS_FILE OUTPUT REFERENCE_SEQUENCE SAMPLES_FILE VCF")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_vcf_to_adpc",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["num_markers"] = _core.first_integer(outs["num_markers_file"])
    extras["samples"] = _core.read_text(outs["samples_file"]).split()
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_arrays_variant_calling_metrics(
    input_path: str,
    dbsnp: str,
    output_prefix: str,
    *,
    call_rate_pf_threshold: float | None = None,
    sequence_dictionary: str | None = None,
    num_processors: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect per-sample (detail) and aggregate (summary) metrics from a genotyping-array VCF:
    call rate, hom/het counts, dbSNP concordance.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Genotyping-array VCF.
    - ``dbsnp`` (``--DBSNP``): dbSNP VCF.
    - ``output_prefix`` (``--OUTPUT``): Prefix of the outputs
      (.arrays_variant_calling_detail_metrics, .arrays_variant_calling_summary_metrics ...).
    - ``call_rate_pf_threshold`` (``--CALL_RATE_PF_THRESHOLD``): Samples with a call rate below
      this fail the filter.
    - ``sequence_dictionary`` (``--SEQUENCE_DICTIONARY``): Sequence dictionary when the VCF has
      no contig lines.
    - ``num_processors`` (``--NUM_PROCESSORS``): Threads; 0 uses all cores.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectArraysVariantCallingMetrics --help`` are
      accepted, a list repeats the option, and options already set by a parameter above are
      refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics_files``, ``summary_metrics`` and ``detail_metrics``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    dbsnp_path = _core.existing_file(dbsnp, "dbsnp")
    arguments += ["--DBSNP", str(dbsnp_path)]
    output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
    arguments += ["--OUTPUT", str(output_prefix_path)]
    outs["output_prefix"] = output_prefix_path
    if call_rate_pf_threshold is not None:
        call_rate_pf_threshold_value = _core.number(
            call_rate_pf_threshold, "call_rate_pf_threshold", minimum=0.0, maximum=1.0
        )
        arguments += ["--CALL_RATE_PF_THRESHOLD", str(call_rate_pf_threshold_value)]
    if sequence_dictionary is not None:
        sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
        arguments += ["--SEQUENCE_DICTIONARY", str(sequence_dictionary_path)]
    if num_processors is not None:
        num_processors_value = _core.integer(num_processors, "num_processors")
        arguments += ["--NUM_PROCESSORS", str(num_processors_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CollectArraysVariantCallingMetrics"
    reserved = _core.names(
        "CALL_RATE_PF_THRESHOLD DBSNP INPUT NUM_PROCESSORS OUTPUT "
        "REFERENCE_SEQUENCE SEQUENCE_DICTIONARY"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_arrays_variant_calling_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    files = sorted(
        str(item)
        for item in output_prefix_path.parent.glob(
            output_prefix_path.name + ".arrays_variant_calling_*"
        )
    )
    extras["metrics_files"] = files
    for suffix, key in (
        (".arrays_variant_calling_summary_metrics", "summary_metrics"),
        (".arrays_variant_calling_detail_metrics", "detail_metrics"),
    ):
        chosen = [Path(item) for item in files if item.endswith(suffix)]
        if chosen:
            extras[key] = _core.picard_metrics(chosen[0])["metrics"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_add_flow_base_quality(
    input_paths: Sequence[str],
    output_path: str,
    *,
    maximal_quality_score: int | None = None,
    minimal_error_rate: float | None = None,
    replace_quality_mode: bool | None = None,
    keep_boundary_flows: bool | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a base-quality attribute to the reads of a flow-based SAM/BAM/CRAM file (Ultima-style
    flow data).

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_path`` (``--output``): New SAM/BAM/CRAM file.
    - ``maximal_quality_score`` (``--maximal-quality-score``): Largest quality written.
    - ``minimal_error_rate`` (``--minimal-error-rate``): Smallest error rate used.
    - ``replace_quality_mode`` (``--replace-quality-mode``): Replace the existing qualities
      instead of adding an attribute.
    - ``keep_boundary_flows`` (``--keep-boundary-flows``): Keep the first and last flows.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk AddFlowBaseQuality --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if maximal_quality_score is not None:
        maximal_quality_score_value = _core.integer(
            maximal_quality_score, "maximal_quality_score", minimum=0
        )
        arguments += ["--maximal-quality-score", str(maximal_quality_score_value)]
    if minimal_error_rate is not None:
        minimal_error_rate_value = _core.number(
            minimal_error_rate, "minimal_error_rate", minimum=0.0, maximum=1.0
        )
        arguments += ["--minimal-error-rate", str(minimal_error_rate_value)]
    if replace_quality_mode is not None:
        arguments += [
            "--replace-quality-mode",
            "true" if _core.flag(replace_quality_mode, "replace_quality_mode") else "false",
        ]
    if keep_boundary_flows is not None:
        arguments += [
            "--keep-boundary-flows",
            "true" if _core.flag(keep_boundary_flows, "keep_boundary_flows") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "AddFlowBaseQuality"
    reserved = _core.names(
        "input intervals keep-boundary-flows maximal-quality-score "
        "minimal-error-rate output reference replace-quality-mode"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_add_flow_base_quality",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_add_flow_snv_quality(
    input_paths: Sequence[str],
    output_path: str,
    *,
    snvq_mode: str | None = None,
    output_quality_attribute: str | None = None,
    max_phred_score: float | None = None,
    keep_supplementary_alignments: bool | None = None,
    include_qc_failed_read: bool | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add an SNV-quality attribute to the reads of a flow-based CRAM/BAM (Ultima-style flow data).

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_path`` (``--output``): New SAM/BAM/CRAM file.
    - ``snvq_mode`` (``--snvq-mode``): How the SNV quality is computed; GATK default Geometric.
    - ``output_quality_attribute`` (``--output-quality-attribute``): Name of the attribute that
      receives the SNV quality.
    - ``max_phred_score`` (``--max-phred-score``): Largest phred score written.
    - ``keep_supplementary_alignments`` (``--keep-supplementary-alignments``): Keep
      supplementary alignments.
    - ``include_qc_failed_read`` (``--include-qc-failed-read``): Keep reads that fail the vendor
      filter.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk AddFlowSNVQuality --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if snvq_mode is not None:
        snvq_mode_value = _core.choice(
            snvq_mode, ("Legacy", "Optimistic", "Pessimistic", "Geometric"), "snvq_mode"
        )
        arguments += ["--snvq-mode", snvq_mode_value]
    if output_quality_attribute is not None:
        output_quality_attribute_text = _core.text(
            output_quality_attribute, "output_quality_attribute"
        )
        arguments += ["--output-quality-attribute", output_quality_attribute_text]
    if max_phred_score is not None:
        max_phred_score_value = _core.number(max_phred_score, "max_phred_score")
        arguments += ["--max-phred-score", str(max_phred_score_value)]
    if keep_supplementary_alignments is not None:
        arguments += [
            "--keep-supplementary-alignments",
            "true"
            if _core.flag(keep_supplementary_alignments, "keep_supplementary_alignments")
            else "false",
        ]
    if include_qc_failed_read is not None:
        arguments += [
            "--include-qc-failed-read",
            "true" if _core.flag(include_qc_failed_read, "include_qc_failed_read") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "AddFlowSNVQuality"
    reserved = _core.names(
        "include-qc-failed-read input intervals keep-supplementary-alignments "
        "max-phred-score output output-quality-attribute reference snvq-mode"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_add_flow_snv_quality",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_calculate_average_combined_annotations(
    variant: str,
    summed_annotations_to_divide: Sequence[str],
    output_path: str,
    *,
    intervals: Sequence[str] | None = None,
    reference: str | None = None,
    variant_output_filtering: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Divide annotations that GenomicsDB summed over samples by the number of samples to give
    their average.

    Parameters (the GATK option is in parentheses):

    - ``variant`` (``--variant``): VCF with the summed annotations.
    - ``summed_annotations_to_divide`` (``--summed-annotation-to-divide``): Names of the
      annotations to divide by the sample count.
    - ``output_path`` (``--output``): New VCF.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``variant_output_filtering`` (``--variant-output-filtering``): Restrict the output
      variants to the intervals by this matching mode.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CalculateAverageCombinedAnnotations --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    variant_path = _core.existing_file(variant, "variant")
    arguments += ["--variant", str(variant_path)]
    summed_annotations_to_divide_values = _core.texts(
        summed_annotations_to_divide, "summed_annotations_to_divide"
    )
    for item in summed_annotations_to_divide_values:
        arguments += ["--summed-annotation-to-divide", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if variant_output_filtering is not None:
        variant_output_filtering_value = _core.choice(
            variant_output_filtering,
            ("STARTS_IN", "ENDS_IN", "OVERLAPS", "CONTAINED", "ANYWHERE"),
            "variant_output_filtering",
        )
        arguments += ["--variant-output-filtering", variant_output_filtering_value]
    tool = "CalculateAverageCombinedAnnotations"
    reserved = _core.names(
        "intervals output reference summed-annotation-to-divide variant variant-output-filtering"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_calculate_average_combined_annotations",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["record_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_flow_feature_mapper(
    input_paths: Sequence[str],
    output_path: str,
    *,
    mapping_feature: str | None = None,
    copy_attrs: Sequence[str] | None = None,
    copy_attr_prefix: str | None = None,
    min_score: float | None = None,
    max_score: float | None = None,
    limit_score: float | None = None,
    snv_identical_bases: int | None = None,
    snv_identical_bases_after: int | None = None,
    include_dup_reads: bool | None = None,
    include_qc_failed_reads: bool | None = None,
    keep_negatives: bool | None = None,
    exclude_nan_scores: bool | None = None,
    keep_supplementary_alignments: bool | None = None,
    emit_ref_confidence: str | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map features (initially SNVs) in a flow-based BAM and write them as a VCF.

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_path`` (``--output``): New VCF.
    - ``mapping_feature`` (``--mapping-feature``): Kind of feature to map.
    - ``copy_attrs`` (``--copy-attr``): Read attributes copied to the VCF records.
    - ``copy_attr_prefix`` (``--copy-attr-prefix``): Prefix of the copied attribute names.
    - ``min_score`` (``--min-score``): Minimum feature score kept.
    - ``max_score`` (``--max-score``): Maximum feature score kept.
    - ``limit_score`` (``--limit-score``): Scores are capped at this value.
    - ``snv_identical_bases`` (``--snv-identical-bases``): Identical bases required before the
      SNV.
    - ``snv_identical_bases_after`` (``--snv-identical-bases-after``): Identical bases required
      after the SNV.
    - ``include_dup_reads`` (``--include-dup-reads``): Use duplicate-marked reads.
    - ``include_qc_failed_reads`` (``--include-qc-failed-reads``): Use reads that fail the
      vendor filter.
    - ``keep_negatives`` (``--keep-negatives``): Keep features with negative scores.
    - ``exclude_nan_scores`` (``--exclude-nan-scores``): Drop features whose score is NaN.
    - ``keep_supplementary_alignments`` (``--keep-supplementary-alignments``): Keep
      supplementary alignments.
    - ``emit_ref_confidence`` (``--emit-ref-confidence``): Mode for emitting reference
      confidence.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk FlowFeatureMapper --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if mapping_feature is not None:
        mapping_feature_value = _core.choice(mapping_feature, ("SNV",), "mapping_feature")
        arguments += ["--mapping-feature", mapping_feature_value]
    if copy_attrs is not None:
        copy_attrs_values = _core.texts(copy_attrs, "copy_attrs")
        for item in copy_attrs_values:
            arguments += ["--copy-attr", item]
    if copy_attr_prefix is not None:
        copy_attr_prefix_text = _core.text(copy_attr_prefix, "copy_attr_prefix")
        arguments += ["--copy-attr-prefix", copy_attr_prefix_text]
    if min_score is not None:
        min_score_value = _core.number(min_score, "min_score")
        arguments += ["--min-score", str(min_score_value)]
    if max_score is not None:
        max_score_value = _core.number(max_score, "max_score")
        arguments += ["--max-score", str(max_score_value)]
    if limit_score is not None:
        limit_score_value = _core.number(limit_score, "limit_score")
        arguments += ["--limit-score", str(limit_score_value)]
    if snv_identical_bases is not None:
        snv_identical_bases_value = _core.integer(
            snv_identical_bases, "snv_identical_bases", minimum=0
        )
        arguments += ["--snv-identical-bases", str(snv_identical_bases_value)]
    if snv_identical_bases_after is not None:
        snv_identical_bases_after_value = _core.integer(
            snv_identical_bases_after, "snv_identical_bases_after", minimum=0
        )
        arguments += ["--snv-identical-bases-after", str(snv_identical_bases_after_value)]
    if include_dup_reads is not None:
        arguments += [
            "--include-dup-reads",
            "true" if _core.flag(include_dup_reads, "include_dup_reads") else "false",
        ]
    if include_qc_failed_reads is not None:
        arguments += [
            "--include-qc-failed-reads",
            "true" if _core.flag(include_qc_failed_reads, "include_qc_failed_reads") else "false",
        ]
    if keep_negatives is not None:
        arguments += [
            "--keep-negatives",
            "true" if _core.flag(keep_negatives, "keep_negatives") else "false",
        ]
    if exclude_nan_scores is not None:
        arguments += [
            "--exclude-nan-scores",
            "true" if _core.flag(exclude_nan_scores, "exclude_nan_scores") else "false",
        ]
    if keep_supplementary_alignments is not None:
        arguments += [
            "--keep-supplementary-alignments",
            "true"
            if _core.flag(keep_supplementary_alignments, "keep_supplementary_alignments")
            else "false",
        ]
    if emit_ref_confidence is not None:
        emit_ref_confidence_value = _core.choice(
            emit_ref_confidence, ("NONE", "BP_RESOLUTION", "GVCF"), "emit_ref_confidence"
        )
        arguments += ["--emit-ref-confidence", emit_ref_confidence_value]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "FlowFeatureMapper"
    reserved = _core.names(
        "copy-attr copy-attr-prefix emit-ref-confidence exclude-nan-scores "
        "include-dup-reads include-qc-failed-reads input intervals "
        "keep-negatives keep-supplementary-alignments limit-score "
        "mapping-feature max-score min-score output reference "
        "snv-identical-bases snv-identical-bases-after"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_flow_feature_mapper",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["record_count"] = _core.count_data_lines(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_flow_pair_hmm_align_reads_to_haplotypes(
    input_paths: Sequence[str],
    haplotypes: str,
    output_path: str,
    *,
    aligner: str | None = None,
    likelihood_calculation_engine: str | None = None,
    ref_haplotype: str | None = None,
    concise_output_format: bool | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Produce a read x haplotype matrix with the log-likelihood of every read against every
    haplotype using the flow-based PairHMM.

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``haplotypes`` (``--haplotypes``): Haplotypes as a FASTA/BAM file.
    - ``output_path`` (``--output``): File that receives the likelihood matrix (GATK takes it as
      a plain path).
    - ``aligner`` (``--aligner``): Likelihood implementation used for alignment.
    - ``likelihood_calculation_engine`` (``--likelihood-calculation-engine``): Likelihood
      calculation engine.
    - ``ref_haplotype`` (``--ref-haplotype``): Name of the reference haplotype.
    - ``concise_output_format`` (``--concise-output-format``): Write a concise matrix.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk FlowPairHMMAlignReadsToHaplotypes --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``matrix_written``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    haplotypes_path = _core.existing_file(haplotypes, "haplotypes")
    arguments += ["--haplotypes", str(haplotypes_path)]
    output_path_path = _core.new_prefix(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    outs["output_path"] = output_path_path
    if aligner is not None:
        aligner_value = _core.choice(aligner, ("PairHMM", "FlowBased", "FlowBasedHMM"), "aligner")
        arguments += ["--aligner", aligner_value]
    if likelihood_calculation_engine is not None:
        likelihood_calculation_engine_value = _core.choice(
            likelihood_calculation_engine,
            ("PairHMM", "FlowBased", "FlowBasedHMM"),
            "likelihood_calculation_engine",
        )
        arguments += ["--likelihood-calculation-engine", likelihood_calculation_engine_value]
    if ref_haplotype is not None:
        ref_haplotype_text = _core.text(ref_haplotype, "ref_haplotype")
        arguments += ["--ref-haplotype", ref_haplotype_text]
    if concise_output_format is not None:
        arguments += [
            "--concise-output-format",
            "true" if _core.flag(concise_output_format, "concise_output_format") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "FlowPairHMMAlignReadsToHaplotypes"
    reserved = _core.names(
        "aligner concise-output-format haplotypes input intervals "
        "likelihood-calculation-engine output ref-haplotype reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_flow_pair_hmm_align_reads_to_haplotypes",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["matrix_written"] = output_path_path.is_file()
    if not extras["matrix_written"]:
        raise BackendExecutionError(
            "gatk FlowPairHMMAlignReadsToHaplotypes did not write the matrix file.",
            code="GATK_OUTPUT_MISSING",
        )
    return _core.finish(result, outputs=outs, **extras)


def gatk_ground_truth_reads_builder(
    input_paths: Sequence[str],
    maternal_ref: str,
    paternal_ref: str,
    ancestral_translators_base_path: str,
    output_csv: str,
    *,
    subsampling_ratio: float | None = None,
    max_output_reads: int | None = None,
    min_haplotype_score: float | None = None,
    min_mq: float | None = None,
    include_supp_align: bool | None = None,
    fill_trimmed_reads: bool | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Produce a flexible and robust ground-truth set of reads for base-calling training from reads
    and the maternal and paternal haplotype references (CSV).

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``maternal_ref`` (``--maternal-ref``): Maternal haplotype reference FASTA.
    - ``paternal_ref`` (``--paternal-ref``): Paternal haplotype reference FASTA.
    - ``ancestral_translators_base_path`` (``--ancestral-translators-base-path``): Base path of
      the ancestral.<contig>.csv translation files.
    - ``output_csv`` (``--output-csv``): New main CSV (.csv or .csv.gz).
    - ``subsampling_ratio`` (``--subsampling-ratio``): Fraction of reads processed.
    - ``max_output_reads`` (``--max-output-reads``): Maximum number of reads written.
    - ``min_haplotype_score`` (``--min-haplotype-score``): Minimum haplotype score of a read.
    - ``min_mq`` (``--min-mq``): Minimum mapping quality.
    - ``include_supp_align`` (``--include-supp-align``): Include supplementary alignments.
    - ``fill_trimmed_reads`` (``--fill-trimmed-reads``): Fill trimmed ends of reads.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GroundTruthReadsBuilder --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    maternal_ref_path = _core.existing_file(maternal_ref, "maternal_ref")
    arguments += ["--maternal-ref", str(maternal_ref_path)]
    paternal_ref_path = _core.existing_file(paternal_ref, "paternal_ref")
    arguments += ["--paternal-ref", str(paternal_ref_path)]
    ancestral_translators_base_path_path = _core.new_prefix(
        ancestral_translators_base_path, "ancestral_translators_base_path"
    )
    arguments += ["--ancestral-translators-base-path", str(ancestral_translators_base_path_path)]
    outs["ancestral_translators_base_path"] = ancestral_translators_base_path_path
    output_csv_path = _core.new_file(output_csv, "output_csv")
    arguments += ["--output-csv", str(output_csv_path)]
    expected.append(output_csv_path)
    outs["output_csv"] = output_csv_path
    if subsampling_ratio is not None:
        subsampling_ratio_value = _core.number(
            subsampling_ratio, "subsampling_ratio", minimum=0.0, maximum=1.0
        )
        arguments += ["--subsampling-ratio", str(subsampling_ratio_value)]
    if max_output_reads is not None:
        max_output_reads_value = _core.integer(max_output_reads, "max_output_reads", minimum=1)
        arguments += ["--max-output-reads", str(max_output_reads_value)]
    if min_haplotype_score is not None:
        min_haplotype_score_value = _core.number(min_haplotype_score, "min_haplotype_score")
        arguments += ["--min-haplotype-score", str(min_haplotype_score_value)]
    if min_mq is not None:
        min_mq_value = _core.number(min_mq, "min_mq", minimum=0.0)
        arguments += ["--min-mq", str(min_mq_value)]
    if include_supp_align is not None:
        arguments += [
            "--include-supp-align",
            "true" if _core.flag(include_supp_align, "include_supp_align") else "false",
        ]
    if fill_trimmed_reads is not None:
        arguments += [
            "--fill-trimmed-reads",
            "true" if _core.flag(fill_trimmed_reads, "fill_trimmed_reads") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "GroundTruthReadsBuilder"
    reserved = _core.names(
        "ancestral-translators-base-path fill-trimmed-reads include-supp-align "
        "input intervals maternal-ref max-output-reads min-haplotype-score "
        "min-mq output-csv paternal-ref reference subsampling-ratio"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_ground_truth_reads_builder",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_ground_truth_scorer(
    input_paths: Sequence[str],
    output_csv: str,
    *,
    report_file: str | None = None,
    features_file: str | None = None,
    genome_prior: str | None = None,
    add_mean_call: bool | None = None,
    exclude_zero_flows: bool | None = None,
    omit_zeros_from_report: bool | None = None,
    quality_percentiles: str | None = None,
    use_softclipped_bases: bool | None = None,
    normalized_score_threshold: float | None = None,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score flow-based reads against a reference or ground truth and write the per-read CSV and an
    optional report.

    Flow-based tool for Ultima Genomics-style flow data (reads carry flow-order tags); GATK
    marks the flow tools experimental.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_csv`` (``--output-csv``): New main CSV (.csv or .csv.gz).
    - ``report_file`` (``--report-file``): New report file.
    - ``features_file`` (``--features-file``): VCF with features used to filter reads.
    - ``genome_prior`` (``--genome-prior``): CSV with genome priors (hmer frequencies).
    - ``add_mean_call`` (``--add-mean-call``): Add the mean call to the CSV.
    - ``exclude_zero_flows`` (``--exclude-zero-flows``): Exclude zero flows from the statistics.
    - ``omit_zeros_from_report`` (``--omit-zeros-from-report``): Omit zero flows from the
      report.
    - ``quality_percentiles`` (``--quality-percentiles``): Comma-separated quality percentiles
      for the report.
    - ``use_softclipped_bases`` (``--use-softclipped-bases``): Use soft-clipped bases.
    - ``normalized_score_threshold`` (``--normalized-score-threshold``): Threshold of the
      normalised score.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GroundTruthScorer --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    output_csv_path = _core.new_file(output_csv, "output_csv")
    arguments += ["--output-csv", str(output_csv_path)]
    expected.append(output_csv_path)
    outs["output_csv"] = output_csv_path
    if report_file is not None:
        report_file_path = _core.new_file(report_file, "report_file")
        arguments += ["--report-file", str(report_file_path)]
        expected.append(report_file_path)
        outs["report_file"] = report_file_path
    if features_file is not None:
        features_file_path = _core.existing_file(features_file, "features_file")
        arguments += ["--features-file", str(features_file_path)]
    if genome_prior is not None:
        genome_prior_path = _core.existing_file(genome_prior, "genome_prior")
        arguments += ["--genome-prior", str(genome_prior_path)]
    if add_mean_call is not None:
        arguments += [
            "--add-mean-call",
            "true" if _core.flag(add_mean_call, "add_mean_call") else "false",
        ]
    if exclude_zero_flows is not None:
        arguments += [
            "--exclude-zero-flows",
            "true" if _core.flag(exclude_zero_flows, "exclude_zero_flows") else "false",
        ]
    if omit_zeros_from_report is not None:
        arguments += [
            "--omit-zeros-from-report",
            "true" if _core.flag(omit_zeros_from_report, "omit_zeros_from_report") else "false",
        ]
    if quality_percentiles is not None:
        quality_percentiles_text = _core.text(quality_percentiles, "quality_percentiles")
        arguments += ["--quality-percentiles", quality_percentiles_text]
    if use_softclipped_bases is not None:
        arguments += [
            "--use-softclipped-bases",
            "true" if _core.flag(use_softclipped_bases, "use_softclipped_bases") else "false",
        ]
    if normalized_score_threshold is not None:
        normalized_score_threshold_value = _core.number(
            normalized_score_threshold, "normalized_score_threshold"
        )
        arguments += ["--normalized-score-threshold", str(normalized_score_threshold_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    tool = "GroundTruthScorer"
    reserved = _core.names(
        "add-mean-call exclude-zero-flows features-file genome-prior input "
        "intervals normalized-score-threshold omit-zeros-from-report output-csv "
        "quality-percentiles reference report-file use-softclipped-bases"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_ground_truth_scorer",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_pathseq_build_kmers(
    reference: str,
    output_path: str,
    *,
    kmer_size: int | None = None,
    kmer_mask: str | None = None,
    kmer_spacing: int | None = None,
    bloom_false_positive_probability: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the set of host reference k-mers (hash set, or Bloom filter when a false-positive
    probability is given) used by PathSeq to subtract host reads.

    Parameters (the GATK option is in parentheses):

    - ``reference`` (``--reference``): Host reference FASTA.
    - ``output_path`` (``--output``): Output path; GATK appends .hss (hash set) or .bfi (Bloom
      filter) when it is missing.
    - ``kmer_size`` (``--kmer-size``): k-mer size.
    - ``kmer_mask`` (``--kmer-mask``): Comma-separated 0-based base indices masked in each
      k-mer.
    - ``kmer_spacing`` (``--kmer-spacing``): Spacing between stored k-mers.
    - ``bloom_false_positive_probability`` (``--bloom-false-positive-probability``): Build a
      Bloom filter with this false-positive probability.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PathSeqBuildKmers --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``kmer_files`` (the .hss or .bfi file).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_prefix(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    outs["output_path"] = output_path_path
    if kmer_size is not None:
        kmer_size_value = _core.integer(kmer_size, "kmer_size", minimum=1)
        arguments += ["--kmer-size", str(kmer_size_value)]
    if kmer_mask is not None:
        kmer_mask_text = _core.text(kmer_mask, "kmer_mask")
        arguments += ["--kmer-mask", kmer_mask_text]
    if kmer_spacing is not None:
        kmer_spacing_value = _core.integer(kmer_spacing, "kmer_spacing", minimum=1)
        arguments += ["--kmer-spacing", str(kmer_spacing_value)]
    if bloom_false_positive_probability is not None:
        bloom_false_positive_probability_value = _core.number(
            bloom_false_positive_probability,
            "bloom_false_positive_probability",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--bloom-false-positive-probability",
            str(bloom_false_positive_probability_value),
        ]
    tool = "PathSeqBuildKmers"
    reserved = _core.names(
        "bloom-false-positive-probability kmer-mask kmer-size kmer-spacing output reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_pathseq_build_kmers",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    found = sorted(str(item) for item in output_path_path.parent.glob(output_path_path.name + "*"))
    if not found:
        raise BackendExecutionError(
            "gatk PathSeqBuildKmers wrote no k-mer file.", code="GATK_OUTPUT_MISSING"
        )
    extras["kmer_files"] = found
    return _core.finish(result, outputs=outs, **extras)


def gatk_pathseq_build_reference_taxonomy(
    reference: str,
    tax_dump: str,
    output_path: str,
    *,
    refseq_catalog: str | None = None,
    genbank_catalog: str | None = None,
    min_non_virus_contig_length: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the taxonomy data file (.db) of the microbe reference for PathSeq from the NCBI
    taxonomy dump and the reference FASTA.

    Parameters (the GATK option is in parentheses):

    - ``reference`` (``--reference``): Microbe reference FASTA.
    - ``tax_dump`` (``--tax-dump``): NCBI taxonomy dump tarball (taxdump.tar.gz).
    - ``output_path`` (``--output``): New taxonomy file, by convention with the extension .db.
    - ``refseq_catalog`` (``--refseq-catalog``): RefSeq release catalog, to map accessions to
      taxa.
    - ``genbank_catalog`` (``--genbank-catalog``): GenBank release catalog, to map accessions to
      taxa.
    - ``min_non_virus_contig_length`` (``--min-non-virus-contig-length``): Non-viral contigs
      shorter than this are dropped.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PathSeqBuildReferenceTaxonomy --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    tax_dump_path = _core.existing_file(tax_dump, "tax_dump")
    arguments += ["--tax-dump", str(tax_dump_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if refseq_catalog is not None:
        refseq_catalog_path = _core.existing_file(refseq_catalog, "refseq_catalog")
        arguments += ["--refseq-catalog", str(refseq_catalog_path)]
    if genbank_catalog is not None:
        genbank_catalog_path = _core.existing_file(genbank_catalog, "genbank_catalog")
        arguments += ["--genbank-catalog", str(genbank_catalog_path)]
    if min_non_virus_contig_length is not None:
        min_non_virus_contig_length_value = _core.integer(
            min_non_virus_contig_length, "min_non_virus_contig_length", minimum=0
        )
        arguments += ["--min-non-virus-contig-length", str(min_non_virus_contig_length_value)]
    tool = "PathSeqBuildReferenceTaxonomy"
    reserved = _core.names(
        "genbank-catalog min-non-virus-contig-length output reference refseq-catalog tax-dump"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_pathseq_build_reference_taxonomy",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)
