"""GATK command-line facades, part 2: copy number, read counts, coverage and pileups.

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

This part covers the somatic and germline copy-number pipeline (interval preparation, read and
allelic counts, denoising, segmentation, calling, plots, germline CNV), allele-specific counts,
callable loci, depth of coverage, expression counts, Mutect2 helper tables and base-quality
recalibration plots.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics import _gatk_p00_core as _core
from dnakit.exceptions import ConfigurationError

__all__ = [
    "gatk_analyze_covariates",
    "gatk_analyze_saturation_mutagenesis",
    "gatk_annotate_intervals",
    "gatk_ase_read_counter",
    "gatk_call_copy_ratio_segments",
    "gatk_callable_loci",
    "gatk_collect_allelic_counts",
    "gatk_collect_f1r2_counts",
    "gatk_collect_read_counts",
    "gatk_combine_segment_breakpoints",
    "gatk_count_reads_or_bases",
    "gatk_create_read_count_panel_of_normals",
    "gatk_denoise_read_counts",
    "gatk_depth_of_coverage",
    "gatk_determine_germline_contig_ploidy",
    "gatk_filter_intervals",
    "gatk_gene_expression_evaluation",
    "gatk_germline_cnv_caller",
    "gatk_get_normal_artifact_data",
    "gatk_get_pileup_summaries",
    "gatk_local_assembler",
    "gatk_merge_annotated_regions",
    "gatk_model_segments",
    "gatk_plot_denoised_copy_ratios",
    "gatk_plot_modeled_segments",
    "gatk_postprocess_germline_cnv_calls",
    "gatk_preprocess_intervals",
    "gatk_tag_germline_events",
]


def gatk_preprocess_intervals(
    reference: str,
    output_path: str,
    *,
    intervals: Sequence[str] | None = None,
    bin_length: int | None = None,
    padding: int | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Prepare equal-sized bins (or targets with padding) for coverage collection and write them as
    an interval list.

    Parameters (the GATK option is in parentheses):

    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New interval list.
    - ``intervals`` (``--intervals``): Loci or interval files to bin; None bins the whole
      reference.
    - ``bin_length`` (``--bin-length``): Bin length in bp; 0 keeps the targets unbinned (exome
      style).
    - ``padding`` (``--padding``): Bases added to both ends of each target or bin.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PreprocessIntervals --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if bin_length is not None:
        bin_length_value = _core.integer(bin_length, "bin_length", minimum=0)
        arguments += ["--bin-length", str(bin_length_value)]
    if padding is not None:
        padding_value = _core.integer(padding, "padding", minimum=0)
        arguments += ["--padding", str(padding_value)]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "PreprocessIntervals"
    reserved = _core.names("bin-length interval-merging-rule intervals output padding reference")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_preprocess_intervals",
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


def gatk_annotate_intervals(
    intervals: Sequence[str],
    reference: str,
    output_path: str,
    *,
    mappability_track: str | None = None,
    segmental_duplication_track: str | None = None,
    feature_query_lookahead: int | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate intervals with GC content and, when tracks are given, mappability and segmental-
    duplication content.

    Parameters (the GATK option is in parentheses):

    - ``intervals`` (``--intervals``): Interval lists or loci to annotate.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New annotated-intervals TSV.
    - ``mappability_track`` (``--mappability-track``): Mappability track (BED/bigWig readable by
      GATK).
    - ``segmental_duplication_track`` (``--segmental-duplication-track``): Segmental-duplication
      track (BED).
    - ``feature_query_lookahead`` (``--feature-query-lookahead``): Look-ahead window for feature
      queries.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk AnnotateIntervals --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (CONTIG, START, END, GC_CONTENT ...).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if mappability_track is not None:
        mappability_track_path = _core.existing_file(mappability_track, "mappability_track")
        arguments += ["--mappability-track", str(mappability_track_path)]
    if segmental_duplication_track is not None:
        segmental_duplication_track_path = _core.existing_file(
            segmental_duplication_track, "segmental_duplication_track"
        )
        arguments += ["--segmental-duplication-track", str(segmental_duplication_track_path)]
    if feature_query_lookahead is not None:
        feature_query_lookahead_value = _core.integer(
            feature_query_lookahead, "feature_query_lookahead", minimum=1
        )
        arguments += ["--feature-query-lookahead", str(feature_query_lookahead_value)]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "AnnotateIntervals"
    reserved = _core.names(
        "feature-query-lookahead interval-merging-rule intervals "
        "mappability-track output reference segmental-duplication-track"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_annotate_intervals",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        extras["rows"] = _core.tsv_rows(outs["output_path"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_read_counts(
    input_paths: Sequence[str],
    intervals: Sequence[str],
    output_path: str,
    *,
    output_format: str | None = None,
    reference: str | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect read counts at specified intervals (coverage-based copy-number input) into a TSV or
    HDF5 file.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``intervals`` (``--intervals``): Interval lists or loci (usually the output of
      PreprocessIntervals).
    - ``output_path`` (``--output``): New counts file.
    - ``output_format`` (``--format``): Output format; GATK default HDF5.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CollectReadCounts --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (first 200 rows, TSV output only).
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
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if output_format is not None:
        output_format_value = _core.choice(output_format, ("TSV", "HDF5"), "output_format")
        arguments += ["--format", output_format_value]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "CollectReadCounts"
    reserved = _core.names("format input interval-merging-rule intervals output reference")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_read_counts",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if str(outs["output_path"]).endswith(".tsv") or output_format == "TSV":
        extras["rows"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_allelic_counts(
    input_paths: Sequence[str],
    intervals: Sequence[str],
    reference: str,
    output_path: str,
    *,
    minimum_base_quality: int | None = None,
    max_depth_per_sample: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect reference and alternate allele counts at the sites of an interval list or VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``intervals`` (``--intervals``): Sites as an interval list, VCF or loci.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New allelic-counts TSV.
    - ``minimum_base_quality`` (``--minimum-base-quality``): Bases below this quality are not
      counted.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CollectAllelicCounts --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (CONTIG, POSITION, REF_COUNT, ALT_COUNT,
    REF_NUCLEOTIDE, ALT_NUCLEOTIDE; first 200).
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
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if minimum_base_quality is not None:
        minimum_base_quality_value = _core.integer(
            minimum_base_quality, "minimum_base_quality", minimum=0
        )
        arguments += ["--minimum-base-quality", str(minimum_base_quality_value)]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    tool = "CollectAllelicCounts"
    reserved = _core.names(
        "input intervals max-depth-per-sample minimum-base-quality output reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_allelic_counts",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["rows"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_filter_intervals(
    intervals: Sequence[str],
    output_path: str,
    *,
    annotated_intervals: str | None = None,
    input_paths: Sequence[str] | None = None,
    minimum_gc_content: float | None = None,
    maximum_gc_content: float | None = None,
    minimum_mappability: float | None = None,
    maximum_mappability: float | None = None,
    minimum_segmental_duplication_content: float | None = None,
    maximum_segmental_duplication_content: float | None = None,
    low_count_filter_count_threshold: int | None = None,
    low_count_filter_percentage_of_samples: float | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter intervals by annotations (GC content, mappability, segmental duplications) and by
    count statistics.

    Parameters (the GATK option is in parentheses):

    - ``intervals`` (``--intervals``): Interval lists to filter.
    - ``output_path`` (``--output``): New filtered interval list.
    - ``annotated_intervals`` (``--annotated-intervals``): Output of AnnotateIntervals, for
      annotation-based filters.
    - ``input_paths`` (``--input``): Read-count files (output of CollectReadCounts), for count-
      based filters.
    - ``minimum_gc_content`` (``--minimum-gc-content``): Keep intervals with at least this GC
      fraction.
    - ``maximum_gc_content`` (``--maximum-gc-content``): Keep intervals with at most this GC
      fraction.
    - ``minimum_mappability`` (``--minimum-mappability``): Keep intervals with at least this
      mappability.
    - ``maximum_mappability`` (``--maximum-mappability``): Keep intervals with at most this
      mappability.
    - ``minimum_segmental_duplication_content`` (``--minimum-segmental-duplication-content``):
      Minimum segmental-duplication fraction.
    - ``maximum_segmental_duplication_content`` (``--maximum-segmental-duplication-content``):
      Maximum segmental-duplication fraction.
    - ``low_count_filter_count_threshold`` (``--low-count-filter-count-threshold``): Intervals
      whose count is below this in the given percentage of samples are removed.
    - ``low_count_filter_percentage_of_samples`` (``--low-count-filter-percentage-of-samples``):
      Percentage of samples for the low-count filter.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk FilterIntervals --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``interval_count`` (intervals kept).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if annotated_intervals is not None:
        annotated_intervals_path = _core.existing_file(annotated_intervals, "annotated_intervals")
        arguments += ["--annotated-intervals", str(annotated_intervals_path)]
    if input_paths is not None:
        input_paths_paths = _core.existing_files(input_paths, "input_paths")
        for item in input_paths_paths:
            arguments += ["--input", str(item)]
    if minimum_gc_content is not None:
        minimum_gc_content_value = _core.number(
            minimum_gc_content, "minimum_gc_content", minimum=0.0, maximum=1.0
        )
        arguments += ["--minimum-gc-content", str(minimum_gc_content_value)]
    if maximum_gc_content is not None:
        maximum_gc_content_value = _core.number(
            maximum_gc_content, "maximum_gc_content", minimum=0.0, maximum=1.0
        )
        arguments += ["--maximum-gc-content", str(maximum_gc_content_value)]
    if minimum_mappability is not None:
        minimum_mappability_value = _core.number(
            minimum_mappability, "minimum_mappability", minimum=0.0, maximum=1.0
        )
        arguments += ["--minimum-mappability", str(minimum_mappability_value)]
    if maximum_mappability is not None:
        maximum_mappability_value = _core.number(
            maximum_mappability, "maximum_mappability", minimum=0.0, maximum=1.0
        )
        arguments += ["--maximum-mappability", str(maximum_mappability_value)]
    if minimum_segmental_duplication_content is not None:
        minimum_segmental_duplication_content_value = _core.number(
            minimum_segmental_duplication_content,
            "minimum_segmental_duplication_content",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--minimum-segmental-duplication-content",
            str(minimum_segmental_duplication_content_value),
        ]
    if maximum_segmental_duplication_content is not None:
        maximum_segmental_duplication_content_value = _core.number(
            maximum_segmental_duplication_content,
            "maximum_segmental_duplication_content",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--maximum-segmental-duplication-content",
            str(maximum_segmental_duplication_content_value),
        ]
    if low_count_filter_count_threshold is not None:
        low_count_filter_count_threshold_value = _core.integer(
            low_count_filter_count_threshold, "low_count_filter_count_threshold", minimum=0
        )
        arguments += [
            "--low-count-filter-count-threshold",
            str(low_count_filter_count_threshold_value),
        ]
    if low_count_filter_percentage_of_samples is not None:
        low_count_filter_percentage_of_samples_value = _core.number(
            low_count_filter_percentage_of_samples,
            "low_count_filter_percentage_of_samples",
            minimum=0.0,
            maximum=100.0,
        )
        arguments += [
            "--low-count-filter-percentage-of-samples",
            str(low_count_filter_percentage_of_samples_value),
        ]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "FilterIntervals"
    reserved = _core.names(
        "annotated-intervals input interval-merging-rule intervals "
        "low-count-filter-count-threshold "
        "low-count-filter-percentage-of-samples maximum-gc-content "
        "maximum-mappability maximum-segmental-duplication-content "
        "minimum-gc-content minimum-mappability "
        "minimum-segmental-duplication-content output"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_filter_intervals",
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


def gatk_create_read_count_panel_of_normals(
    input_paths: Sequence[str],
    output_path: str,
    *,
    annotated_intervals: str | None = None,
    number_of_eigensamples: int | None = None,
    minimum_interval_median_percentile: float | None = None,
    maximum_zeros_in_sample_percentage: float | None = None,
    maximum_zeros_in_interval_percentage: float | None = None,
    do_impute_zeros: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a panel of normals (PoN) for read-count denoising from the read counts of normal
    samples (HDF5).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Read-count files of the normal samples (output of
      CollectReadCounts).
    - ``output_path`` (``--output``): New PoN HDF5 file.
    - ``annotated_intervals`` (``--annotated-intervals``): Output of AnnotateIntervals, to
      correct for GC content.
    - ``number_of_eigensamples`` (``--number-of-eigensamples``): Number of principal components
      kept.
    - ``minimum_interval_median_percentile`` (``--minimum-interval-median-percentile``):
      Intervals below this median percentile are filtered.
    - ``maximum_zeros_in_sample_percentage`` (``--maximum-zeros-in-sample-percentage``): Samples
      with more zero-count intervals are filtered.
    - ``maximum_zeros_in_interval_percentage`` (``--maximum-zeros-in-interval-percentage``):
      Intervals with more zero-count samples are filtered.
    - ``do_impute_zeros`` (``--do-impute-zeros``): Impute zero counts.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CreateReadCountPanelOfNormals --help`` are accepted, a list
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
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if annotated_intervals is not None:
        annotated_intervals_path = _core.existing_file(annotated_intervals, "annotated_intervals")
        arguments += ["--annotated-intervals", str(annotated_intervals_path)]
    if number_of_eigensamples is not None:
        number_of_eigensamples_value = _core.integer(
            number_of_eigensamples, "number_of_eigensamples", minimum=0
        )
        arguments += ["--number-of-eigensamples", str(number_of_eigensamples_value)]
    if minimum_interval_median_percentile is not None:
        minimum_interval_median_percentile_value = _core.number(
            minimum_interval_median_percentile,
            "minimum_interval_median_percentile",
            minimum=0.0,
            maximum=100.0,
        )
        arguments += [
            "--minimum-interval-median-percentile",
            str(minimum_interval_median_percentile_value),
        ]
    if maximum_zeros_in_sample_percentage is not None:
        maximum_zeros_in_sample_percentage_value = _core.number(
            maximum_zeros_in_sample_percentage,
            "maximum_zeros_in_sample_percentage",
            minimum=0.0,
            maximum=100.0,
        )
        arguments += [
            "--maximum-zeros-in-sample-percentage",
            str(maximum_zeros_in_sample_percentage_value),
        ]
    if maximum_zeros_in_interval_percentage is not None:
        maximum_zeros_in_interval_percentage_value = _core.number(
            maximum_zeros_in_interval_percentage,
            "maximum_zeros_in_interval_percentage",
            minimum=0.0,
            maximum=100.0,
        )
        arguments += [
            "--maximum-zeros-in-interval-percentage",
            str(maximum_zeros_in_interval_percentage_value),
        ]
    if do_impute_zeros is not None:
        arguments += [
            "--do-impute-zeros",
            "true" if _core.flag(do_impute_zeros, "do_impute_zeros") else "false",
        ]
    tool = "CreateReadCountPanelOfNormals"
    reserved = _core.names(
        "annotated-intervals do-impute-zeros input "
        "maximum-zeros-in-interval-percentage "
        "maximum-zeros-in-sample-percentage minimum-interval-median-percentile "
        "number-of-eigensamples output"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_create_read_count_panel_of_normals",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_denoise_read_counts(
    input_path: str,
    standardized_copy_ratios: str,
    denoised_copy_ratios: str,
    *,
    count_panel_of_normals: str | None = None,
    annotated_intervals: str | None = None,
    number_of_eigensamples: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Denoise read counts with a panel of normals (or annotated intervals) to produce standardized
    and denoised copy ratios.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--input``): Read-count file of the case sample (TSV or HDF5).
    - ``standardized_copy_ratios`` (``--standardized-copy-ratios``): New standardized copy-ratio
      file.
    - ``denoised_copy_ratios`` (``--denoised-copy-ratios``): New denoised copy-ratio file.
    - ``count_panel_of_normals`` (``--count-panel-of-normals``): Panel of normals from
      CreateReadCountPanelOfNormals.
    - ``annotated_intervals`` (``--annotated-intervals``): Annotated intervals, used when no
      panel of normals is given.
    - ``number_of_eigensamples`` (``--number-of-eigensamples``): Principal components used for
      denoising.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk DenoiseReadCounts --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (denoised copy ratios, first 200).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--input", str(input_path_path)]
    standardized_copy_ratios_path = _core.new_file(
        standardized_copy_ratios, "standardized_copy_ratios"
    )
    arguments += ["--standardized-copy-ratios", str(standardized_copy_ratios_path)]
    expected.append(standardized_copy_ratios_path)
    outs["standardized_copy_ratios"] = standardized_copy_ratios_path
    denoised_copy_ratios_path = _core.new_file(denoised_copy_ratios, "denoised_copy_ratios")
    arguments += ["--denoised-copy-ratios", str(denoised_copy_ratios_path)]
    expected.append(denoised_copy_ratios_path)
    outs["denoised_copy_ratios"] = denoised_copy_ratios_path
    if count_panel_of_normals is not None:
        count_panel_of_normals_path = _core.existing_file(
            count_panel_of_normals, "count_panel_of_normals"
        )
        arguments += ["--count-panel-of-normals", str(count_panel_of_normals_path)]
    if annotated_intervals is not None:
        annotated_intervals_path = _core.existing_file(annotated_intervals, "annotated_intervals")
        arguments += ["--annotated-intervals", str(annotated_intervals_path)]
    if number_of_eigensamples is not None:
        number_of_eigensamples_value = _core.integer(
            number_of_eigensamples, "number_of_eigensamples", minimum=0
        )
        arguments += ["--number-of-eigensamples", str(number_of_eigensamples_value)]
    tool = "DenoiseReadCounts"
    reserved = _core.names(
        "annotated-intervals count-panel-of-normals denoised-copy-ratios input "
        "number-of-eigensamples standardized-copy-ratios"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_denoise_read_counts",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["rows"] = _core.tsv_rows(outs["denoised_copy_ratios"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_model_segments(
    output_dir: str,
    output_prefix: str,
    *,
    denoised_copy_ratios: str | None = None,
    allelic_counts: str | None = None,
    normal_allelic_counts: str | None = None,
    minimum_total_allele_count_case: int | None = None,
    minimum_total_allele_count_normal: int | None = None,
    genotyping_homozygous_log_ratio_threshold: float | None = None,
    genotyping_base_error_rate: float | None = None,
    maximum_number_of_segments_per_chromosome: int | None = None,
    window_sizes: Sequence[int] | None = None,
    number_of_changepoints_penalty_factor: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Model segmented copy ratios from denoised copy ratios and segmented minor-allele fractions
    from allelic counts.

    Parameters (the GATK option is in parentheses):

    - ``output_dir`` (``--output``): New folder for the segment files.
    - ``output_prefix`` (``--output-prefix``): Prefix of the output file names (usually the
      sample name).
    - ``denoised_copy_ratios`` (``--denoised-copy-ratios``): Denoised copy ratios (output of
      DenoiseReadCounts).
    - ``allelic_counts`` (``--allelic-counts``): Allelic counts of the case sample (output of
      CollectAllelicCounts).
    - ``normal_allelic_counts`` (``--normal-allelic-counts``): Allelic counts of the matched
      normal.
    - ``minimum_total_allele_count_case`` (``--minimum-total-allele-count-case``): Minimum total
      count for a het site in the case sample.
    - ``minimum_total_allele_count_normal`` (``--minimum-total-allele-count-normal``): Minimum
      total count for a het site in the normal.
    - ``genotyping_homozygous_log_ratio_threshold`` (``--genotyping-homozygous-log-ratio-
      threshold``): Log-ratio threshold below which a site is called homozygous.
    - ``genotyping_base_error_rate`` (``--genotyping-base-error-rate``): Base error rate used
      for genotyping.
    - ``maximum_number_of_segments_per_chromosome`` (``--maximum-number-of-segments-per-
      chromosome``): Upper limit of segments per chromosome.
    - ``window_sizes`` (``--window-size``): Window sizes for the multiscale changepoint search.
    - ``number_of_changepoints_penalty_factor`` (``--number-of-changepoints-penalty-factor``):
      Penalty on the number of changepoints.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk ModelSegments --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``segment_files`` and ``segments`` (rows of <prefix>.cr.seg, or
    <prefix>.modelFinal.seg when only allelic counts were given).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if denoised_copy_ratios is None and allelic_counts is None:
        raise ConfigurationError("Give denoised_copy_ratios, allelic_counts or both.")
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--output-prefix", output_prefix_text]
    if denoised_copy_ratios is not None:
        denoised_copy_ratios_path = _core.existing_file(
            denoised_copy_ratios, "denoised_copy_ratios"
        )
        arguments += ["--denoised-copy-ratios", str(denoised_copy_ratios_path)]
    if allelic_counts is not None:
        allelic_counts_path = _core.existing_file(allelic_counts, "allelic_counts")
        arguments += ["--allelic-counts", str(allelic_counts_path)]
    if normal_allelic_counts is not None:
        normal_allelic_counts_path = _core.existing_file(
            normal_allelic_counts, "normal_allelic_counts"
        )
        arguments += ["--normal-allelic-counts", str(normal_allelic_counts_path)]
    if minimum_total_allele_count_case is not None:
        minimum_total_allele_count_case_value = _core.integer(
            minimum_total_allele_count_case, "minimum_total_allele_count_case", minimum=0
        )
        arguments += [
            "--minimum-total-allele-count-case",
            str(minimum_total_allele_count_case_value),
        ]
    if minimum_total_allele_count_normal is not None:
        minimum_total_allele_count_normal_value = _core.integer(
            minimum_total_allele_count_normal, "minimum_total_allele_count_normal", minimum=0
        )
        arguments += [
            "--minimum-total-allele-count-normal",
            str(minimum_total_allele_count_normal_value),
        ]
    if genotyping_homozygous_log_ratio_threshold is not None:
        genotyping_homozygous_log_ratio_threshold_value = _core.number(
            genotyping_homozygous_log_ratio_threshold, "genotyping_homozygous_log_ratio_threshold"
        )
        arguments += [
            "--genotyping-homozygous-log-ratio-threshold",
            str(genotyping_homozygous_log_ratio_threshold_value),
        ]
    if genotyping_base_error_rate is not None:
        genotyping_base_error_rate_value = _core.number(
            genotyping_base_error_rate, "genotyping_base_error_rate", minimum=0.0, maximum=1.0
        )
        arguments += ["--genotyping-base-error-rate", str(genotyping_base_error_rate_value)]
    if maximum_number_of_segments_per_chromosome is not None:
        maximum_number_of_segments_per_chromosome_value = _core.integer(
            maximum_number_of_segments_per_chromosome,
            "maximum_number_of_segments_per_chromosome",
            minimum=1,
        )
        arguments += [
            "--maximum-number-of-segments-per-chromosome",
            str(maximum_number_of_segments_per_chromosome_value),
        ]
    if window_sizes is not None:
        window_sizes_values = _core.integers(window_sizes, "window_sizes", minimum=1)
        for item in window_sizes_values:
            arguments += ["--window-size", str(item)]
    if number_of_changepoints_penalty_factor is not None:
        number_of_changepoints_penalty_factor_value = _core.number(
            number_of_changepoints_penalty_factor,
            "number_of_changepoints_penalty_factor",
            minimum=0.0,
        )
        arguments += [
            "--number-of-changepoints-penalty-factor",
            str(number_of_changepoints_penalty_factor_value),
        ]
    tool = "ModelSegments"
    reserved = _core.names(
        "allelic-counts denoised-copy-ratios genotyping-base-error-rate "
        "genotyping-homozygous-log-ratio-threshold "
        "maximum-number-of-segments-per-chromosome "
        "minimum-total-allele-count-case minimum-total-allele-count-normal "
        "normal-allelic-counts number-of-changepoints-penalty-factor output "
        "output-prefix window-size"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_model_segments",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    files = sorted(str(item) for item in output_dir_path.glob(output_prefix_text + ".*"))
    extras["segment_files"] = files
    final = output_dir_path / (output_prefix_text + ".modelFinal.seg")
    copy_ratio = output_dir_path / (output_prefix_text + ".cr.seg")
    chosen = copy_ratio if copy_ratio.is_file() else final
    extras["segments"] = _core.tsv_rows(chosen, limit=200) if chosen.is_file() else []
    return _core.finish(result, outputs=outs, **extras)


def gatk_call_copy_ratio_segments(
    input_path: str,
    output_path: str,
    *,
    neutral_segment_copy_ratio_lower_bound: float | None = None,
    neutral_segment_copy_ratio_upper_bound: float | None = None,
    outlier_neutral_segment_copy_ratio_z_score_threshold: float | None = None,
    calling_copy_ratio_z_score_threshold: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call copy-ratio segments as amplified (+), deleted (-) or copy-number neutral (0).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--input``): Copy-ratio segment file (<prefix>.cr.seg from
      ModelSegments).
    - ``output_path`` (``--output``): New called-segment file.
    - ``neutral_segment_copy_ratio_lower_bound`` (``--neutral-segment-copy-ratio-lower-bound``):
      Lower bound of neutral copy ratios.
    - ``neutral_segment_copy_ratio_upper_bound`` (``--neutral-segment-copy-ratio-upper-bound``):
      Upper bound of neutral copy ratios.
    - ``outlier_neutral_segment_copy_ratio_z_score_threshold`` (``--outlier-neutral-segment-
      copy-ratio-z-score-threshold``): Z-score threshold that excludes outlier segments from the
      neutral set.
    - ``calling_copy_ratio_z_score_threshold`` (``--calling-copy-ratio-z-score-threshold``):
      Z-score threshold for calling amplifications and deletions.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CallCopyRatioSegments --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``segments`` (rows with a CALL column of +, - or 0; first 200).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--input", str(input_path_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if neutral_segment_copy_ratio_lower_bound is not None:
        neutral_segment_copy_ratio_lower_bound_value = _core.number(
            neutral_segment_copy_ratio_lower_bound,
            "neutral_segment_copy_ratio_lower_bound",
            minimum=0.0,
        )
        arguments += [
            "--neutral-segment-copy-ratio-lower-bound",
            str(neutral_segment_copy_ratio_lower_bound_value),
        ]
    if neutral_segment_copy_ratio_upper_bound is not None:
        neutral_segment_copy_ratio_upper_bound_value = _core.number(
            neutral_segment_copy_ratio_upper_bound,
            "neutral_segment_copy_ratio_upper_bound",
            minimum=0.0,
        )
        arguments += [
            "--neutral-segment-copy-ratio-upper-bound",
            str(neutral_segment_copy_ratio_upper_bound_value),
        ]
    if outlier_neutral_segment_copy_ratio_z_score_threshold is not None:
        outlier_neutral_segment_copy_ratio_z_score_threshold_value = _core.number(
            outlier_neutral_segment_copy_ratio_z_score_threshold,
            "outlier_neutral_segment_copy_ratio_z_score_threshold",
            minimum=0.0,
        )
        arguments += [
            "--outlier-neutral-segment-copy-ratio-z-score-threshold",
            str(outlier_neutral_segment_copy_ratio_z_score_threshold_value),
        ]
    if calling_copy_ratio_z_score_threshold is not None:
        calling_copy_ratio_z_score_threshold_value = _core.number(
            calling_copy_ratio_z_score_threshold,
            "calling_copy_ratio_z_score_threshold",
            minimum=0.0,
        )
        arguments += [
            "--calling-copy-ratio-z-score-threshold",
            str(calling_copy_ratio_z_score_threshold_value),
        ]
    tool = "CallCopyRatioSegments"
    reserved = _core.names(
        "calling-copy-ratio-z-score-threshold input "
        "neutral-segment-copy-ratio-lower-bound "
        "neutral-segment-copy-ratio-upper-bound "
        "outlier-neutral-segment-copy-ratio-z-score-threshold output"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_call_copy_ratio_segments",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["segments"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_combine_segment_breakpoints(
    segments: Sequence[str],
    columns_of_interest: Sequence[str],
    output_path: str,
    *,
    labels: Sequence[str] | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine the breakpoints of two segment files and annotate the resulting intervals with
    chosen columns from each file.

    Parameters (the GATK option is in parentheses):

    - ``segments`` (``--segments``): Exactly two segment files.
    - ``columns_of_interest`` (``--columns-of-interest``): Column names to carry into the
      output; a column found in both files gets the file label as suffix.
    - ``output_path`` (``--output``): New combined-segment file.
    - ``labels`` (``--labels``): Labels (suffixes) of the two files; GATK default 1 and 2.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CombineSegmentBreakpoints --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``segments`` (first 200 rows).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if isinstance(segments, (str, bytes)) or len(segments) != 2:
        raise ConfigurationError("segments must hold exactly two files.")
    segments_paths = _core.existing_files(segments, "segments")
    for item in segments_paths:
        arguments += ["--segments", str(item)]
    columns_of_interest_values = _core.texts(columns_of_interest, "columns_of_interest")
    for item in columns_of_interest_values:
        arguments += ["--columns-of-interest", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if labels is not None:
        labels_values = _core.texts(labels, "labels")
        for item in labels_values:
            arguments += ["--labels", item]
    tool = "CombineSegmentBreakpoints"
    reserved = _core.names("columns-of-interest labels output segments")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_combine_segment_breakpoints",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["segments"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_merge_annotated_regions(
    segments: str,
    output_path: str,
    *,
    reference: str | None = None,
    annotations_to_match: Sequence[str] | None = None,
    max_merge_distance: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge annotated genomic regions. Without annotations_to_match, touching or overlapping
    regions are merged (MergeAnnotatedRegions); with annotations_to_match, regions within
    max_merge_distance whose listed annotation values are identical are merged
    (MergeAnnotatedRegionsByAnnotation).

    Parameters (the GATK option is in parentheses):

    - ``segments`` (``--segments``): Segment file with annotation columns.
    - ``output_path`` (``--output``): New merged-region file.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files (required
      unless annotations_to_match is given).
    - ``annotations_to_match`` (``--annotations-to-match``): Annotation column names whose
      values must be identical for regions to merge.
    - ``max_merge_distance`` (``--max-merge-distance``): Maximum gap between regions that may be
      merged (with annotations_to_match only).
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk <tool> --help`` are accepted, a list repeats the option, and
      options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``regions`` (first 200 rows of the merged file).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if annotations_to_match is None and reference is None:
        raise ConfigurationError("reference is required when annotations_to_match is not given.")
    if annotations_to_match is None and max_merge_distance is not None:
        raise ConfigurationError(
            "max_merge_distance only applies together with annotations_to_match."
        )
    segments_path = _core.existing_file(segments, "segments")
    arguments += ["--segments", str(segments_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if annotations_to_match is not None:
        annotations_to_match_values = _core.texts(annotations_to_match, "annotations_to_match")
        for item in annotations_to_match_values:
            arguments += ["--annotations-to-match", item]
    if max_merge_distance is not None:
        max_merge_distance_value = _core.integer(
            max_merge_distance, "max_merge_distance", minimum=0
        )
        arguments += ["--max-merge-distance", str(max_merge_distance_value)]
    tool = (
        "MergeAnnotatedRegions"
        if annotations_to_match is None
        else "MergeAnnotatedRegionsByAnnotation"
    )
    reserved = _core.names("annotations-to-match max-merge-distance output reference segments")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_merge_annotated_regions",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["regions"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_plot_denoised_copy_ratios(
    standardized_copy_ratios: str,
    denoised_copy_ratios: str,
    sequence_dictionary: str,
    output_dir: str,
    output_prefix: str,
    *,
    minimum_contig_length: int | None = None,
    maximum_copy_ratio: float | None = None,
    point_size_copy_ratio: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create plots of standardized and denoised copy ratios and their median absolute deviations.

    The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK
    stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``standardized_copy_ratios`` (``--standardized-copy-ratios``): Standardized copy ratios
      (output of DenoiseReadCounts).
    - ``denoised_copy_ratios`` (``--denoised-copy-ratios``): Denoised copy ratios (output of
      DenoiseReadCounts).
    - ``sequence_dictionary`` (``--sequence-dictionary``): Sequence dictionary (.dict) of the
      reference.
    - ``output_dir`` (``--output``): New folder for the plots.
    - ``output_prefix`` (``--output-prefix``): Prefix of the output file names.
    - ``minimum_contig_length`` (``--minimum-contig-length``): Contigs shorter than this are not
      plotted.
    - ``maximum_copy_ratio`` (``--maximum-copy-ratio``): Upper limit of the copy-ratio axis.
    - ``point_size_copy_ratio`` (``--point-size-copy-ratio``): Point size of copy-ratio points.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PlotDenoisedCopyRatios --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``plot_files`` (PNG plots and the MAD text files).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    standardized_copy_ratios_path = _core.existing_file(
        standardized_copy_ratios, "standardized_copy_ratios"
    )
    arguments += ["--standardized-copy-ratios", str(standardized_copy_ratios_path)]
    denoised_copy_ratios_path = _core.existing_file(denoised_copy_ratios, "denoised_copy_ratios")
    arguments += ["--denoised-copy-ratios", str(denoised_copy_ratios_path)]
    sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
    arguments += ["--sequence-dictionary", str(sequence_dictionary_path)]
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--output-prefix", output_prefix_text]
    if minimum_contig_length is not None:
        minimum_contig_length_value = _core.integer(
            minimum_contig_length, "minimum_contig_length", minimum=0
        )
        arguments += ["--minimum-contig-length", str(minimum_contig_length_value)]
    if maximum_copy_ratio is not None:
        maximum_copy_ratio_value = _core.number(
            maximum_copy_ratio, "maximum_copy_ratio", minimum=0.0
        )
        arguments += ["--maximum-copy-ratio", str(maximum_copy_ratio_value)]
    if point_size_copy_ratio is not None:
        point_size_copy_ratio_value = _core.number(
            point_size_copy_ratio, "point_size_copy_ratio", minimum=0.0
        )
        arguments += ["--point-size-copy-ratio", str(point_size_copy_ratio_value)]
    tool = "PlotDenoisedCopyRatios"
    reserved = _core.names(
        "denoised-copy-ratios maximum-copy-ratio minimum-contig-length output "
        "output-prefix point-size-copy-ratio sequence-dictionary "
        "standardized-copy-ratios"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_plot_denoised_copy_ratios",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["plot_files"] = sorted(
        str(item) for item in output_dir_path.glob(output_prefix_text + ".*")
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_plot_modeled_segments(
    segments: str,
    sequence_dictionary: str,
    output_dir: str,
    output_prefix: str,
    *,
    denoised_copy_ratios: str | None = None,
    allelic_counts: str | None = None,
    minimum_contig_length: int | None = None,
    maximum_copy_ratio: float | None = None,
    point_size_copy_ratio: float | None = None,
    point_size_allele_fraction: float | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create plots of denoised and segmented copy-ratio and minor-allele-fraction estimates.

    The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK
    stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``segments`` (``--segments``): Modeled segments (<prefix>.modelFinal.seg from
      ModelSegments).
    - ``sequence_dictionary`` (``--sequence-dictionary``): Sequence dictionary (.dict) of the
      reference.
    - ``output_dir`` (``--output``): New folder for the plots.
    - ``output_prefix`` (``--output-prefix``): Prefix of the output file names.
    - ``denoised_copy_ratios`` (``--denoised-copy-ratios``): Denoised copy ratios.
    - ``allelic_counts`` (``--allelic-counts``): Allelic counts.
    - ``minimum_contig_length`` (``--minimum-contig-length``): Contigs shorter than this are not
      plotted.
    - ``maximum_copy_ratio`` (``--maximum-copy-ratio``): Upper limit of the copy-ratio axis.
    - ``point_size_copy_ratio`` (``--point-size-copy-ratio``): Point size of copy-ratio points.
    - ``point_size_allele_fraction`` (``--point-size-allele-fraction``): Point size of allele-
      fraction points.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PlotModeledSegments --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``plot_files``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    segments_path = _core.existing_file(segments, "segments")
    arguments += ["--segments", str(segments_path)]
    sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
    arguments += ["--sequence-dictionary", str(sequence_dictionary_path)]
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--output-prefix", output_prefix_text]
    if denoised_copy_ratios is not None:
        denoised_copy_ratios_path = _core.existing_file(
            denoised_copy_ratios, "denoised_copy_ratios"
        )
        arguments += ["--denoised-copy-ratios", str(denoised_copy_ratios_path)]
    if allelic_counts is not None:
        allelic_counts_path = _core.existing_file(allelic_counts, "allelic_counts")
        arguments += ["--allelic-counts", str(allelic_counts_path)]
    if minimum_contig_length is not None:
        minimum_contig_length_value = _core.integer(
            minimum_contig_length, "minimum_contig_length", minimum=0
        )
        arguments += ["--minimum-contig-length", str(minimum_contig_length_value)]
    if maximum_copy_ratio is not None:
        maximum_copy_ratio_value = _core.number(
            maximum_copy_ratio, "maximum_copy_ratio", minimum=0.0
        )
        arguments += ["--maximum-copy-ratio", str(maximum_copy_ratio_value)]
    if point_size_copy_ratio is not None:
        point_size_copy_ratio_value = _core.number(
            point_size_copy_ratio, "point_size_copy_ratio", minimum=0.0
        )
        arguments += ["--point-size-copy-ratio", str(point_size_copy_ratio_value)]
    if point_size_allele_fraction is not None:
        point_size_allele_fraction_value = _core.number(
            point_size_allele_fraction, "point_size_allele_fraction", minimum=0.0
        )
        arguments += ["--point-size-allele-fraction", str(point_size_allele_fraction_value)]
    tool = "PlotModeledSegments"
    reserved = _core.names(
        "allelic-counts denoised-copy-ratios maximum-copy-ratio "
        "minimum-contig-length output output-prefix point-size-allele-fraction "
        "point-size-copy-ratio segments sequence-dictionary"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_plot_modeled_segments",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["plot_files"] = sorted(
        str(item) for item in output_dir_path.glob(output_prefix_text + ".*")
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_tag_germline_events(
    segments: str,
    called_matched_normal_seg_file: str,
    reference: str,
    output_path: str,
    *,
    endpoint_padding: int | None = None,
    reciprocal_threshold: float | None = None,
    input_call_header: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tag segments of a called tumor segment file that overlap germline events found in the called
    matched normal.

    Parameters (the GATK option is in parentheses):

    - ``segments`` (``--segments``): Called tumor segment file.
    - ``called_matched_normal_seg_file`` (``--called-matched-normal-seg-file``): Called segment
      file of the matched normal.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New tagged segment file.
    - ``endpoint_padding`` (``--endpoint-padding``): Padding in bp when comparing segment
      endpoints.
    - ``reciprocal_threshold`` (``--reciprocal-threshold``): Minimum reciprocal overlap to call
      an event germline.
    - ``input_call_header`` (``--input-call-header``): Name of the column that holds the call.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk TagGermlineEvents --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``segments`` (first 200 rows).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    segments_path = _core.existing_file(segments, "segments")
    arguments += ["--segments", str(segments_path)]
    called_matched_normal_seg_file_path = _core.existing_file(
        called_matched_normal_seg_file, "called_matched_normal_seg_file"
    )
    arguments += ["--called-matched-normal-seg-file", str(called_matched_normal_seg_file_path)]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if endpoint_padding is not None:
        endpoint_padding_value = _core.integer(endpoint_padding, "endpoint_padding", minimum=0)
        arguments += ["--endpoint-padding", str(endpoint_padding_value)]
    if reciprocal_threshold is not None:
        reciprocal_threshold_value = _core.number(
            reciprocal_threshold, "reciprocal_threshold", minimum=0.0, maximum=1.0
        )
        arguments += ["--reciprocal-threshold", str(reciprocal_threshold_value)]
    if input_call_header is not None:
        input_call_header_text = _core.text(input_call_header, "input_call_header")
        arguments += ["--input-call-header", input_call_header_text]
    tool = "TagGermlineEvents"
    reserved = _core.names(
        "called-matched-normal-seg-file endpoint-padding input-call-header "
        "output reciprocal-threshold reference segments"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_tag_germline_events",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["segments"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_determine_germline_contig_ploidy(
    input_paths: Sequence[str],
    output_dir: str,
    output_prefix: str,
    *,
    contig_ploidy_priors: str | None = None,
    model: str | None = None,
    intervals: Sequence[str] | None = None,
    mapping_error_rate: float | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Determine the baseline contig ploidy of germline samples from read counts (cohort mode
    builds a ploidy model, case mode applies one).

    This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the
    GATK conda environment (gatkcondaenv.yml of the release zip).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Read-count files (output of CollectReadCounts) of all
      samples.
    - ``output_dir`` (``--output``): New folder for the ploidy model and calls.
    - ``output_prefix`` (``--output-prefix``): Prefix of the output file names.
    - ``contig_ploidy_priors`` (``--contig-ploidy-priors``): Contig-ploidy priors table
      (required for several samples, not allowed for one).
    - ``model`` (``--model``): Existing ploidy-model folder (required for a single sample).
    - ``intervals`` (``--intervals``): Intervals to analyse.
    - ``mapping_error_rate`` (``--mapping-error-rate``): Mapping error rate.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk DetermineGermlineContigPloidy --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``calls_directories`` (<prefix>-calls) and ``model_directories``
    (<prefix>-model).
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
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--output-prefix", output_prefix_text]
    if contig_ploidy_priors is not None:
        contig_ploidy_priors_path = _core.existing_file(
            contig_ploidy_priors, "contig_ploidy_priors"
        )
        arguments += ["--contig-ploidy-priors", str(contig_ploidy_priors_path)]
    if model is not None:
        model_path = _core.existing_directory(model, "model")
        arguments += ["--model", str(model_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if mapping_error_rate is not None:
        mapping_error_rate_value = _core.number(
            mapping_error_rate, "mapping_error_rate", minimum=0.0, maximum=1.0
        )
        arguments += ["--mapping-error-rate", str(mapping_error_rate_value)]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "DetermineGermlineContigPloidy"
    reserved = _core.names(
        "contig-ploidy-priors input interval-merging-rule intervals "
        "mapping-error-rate model output output-prefix"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_determine_germline_contig_ploidy",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["calls_directories"] = sorted(str(item) for item in output_dir_path.glob("*-calls"))
    extras["model_directories"] = sorted(str(item) for item in output_dir_path.glob("*-model"))
    return _core.finish(result, outputs=outs, **extras)


def gatk_germline_cnv_caller(
    run_mode: str,
    input_paths: Sequence[str],
    contig_ploidy_calls: str,
    output_dir: str,
    output_prefix: str,
    *,
    model: str | None = None,
    annotated_intervals: str | None = None,
    intervals: Sequence[str] | None = None,
    max_copy_number: int | None = None,
    p_alt: float | None = None,
    p_active: float | None = None,
    interval_merging_rule: str | None = "OVERLAPPING_ONLY",
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call copy-number variants in germline samples from read counts and the contig ploidy calls
    (COHORT mode learns a denoising model, CASE mode applies one).

    This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the
    GATK conda environment (gatkcondaenv.yml of the release zip).

    Parameters (the GATK option is in parentheses):

    - ``run_mode`` (``--run-mode``): COHORT builds a model from many samples; CASE calls samples
      with an existing model.
    - ``input_paths`` (``--input``): Read-count files (output of CollectReadCounts).
    - ``contig_ploidy_calls`` (``--contig-ploidy-calls``): Calls folder written by
      DetermineGermlineContigPloidy.
    - ``output_dir`` (``--output``): New folder for the calls and the model.
    - ``output_prefix`` (``--output-prefix``): Prefix of the output file names.
    - ``model`` (``--model``): Existing denoising-model folder (required in CASE mode).
    - ``annotated_intervals`` (``--annotated-intervals``): Output of AnnotateIntervals, to model
      GC bias.
    - ``intervals`` (``--intervals``): Intervals to analyse.
    - ``max_copy_number`` (``--max-copy-number``): Largest copy number called.
    - ``p_alt`` (``--p-alt``): Prior probability of an alternate copy-number state.
    - ``p_active`` (``--p-active``): Prior probability of an active copy-number class.
    - ``interval_merging_rule`` (``--interval-merging-rule``): How overlapping or adjacent
      intervals are merged. GATK's copy-number tools accept only OVERLAPPING_ONLY, which is
      therefore the default here (GATK's own default ALL stops these tools with an error); None
      omits the option.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GermlineCNVCaller --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``calls_directories`` and ``model_directories``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    run_mode_value = _core.choice(run_mode, ("COHORT", "CASE"), "run_mode")
    arguments += ["--run-mode", run_mode_value]
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    contig_ploidy_calls_path = _core.existing_directory(contig_ploidy_calls, "contig_ploidy_calls")
    arguments += ["--contig-ploidy-calls", str(contig_ploidy_calls_path)]
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    arguments += ["--output", str(output_dir_path)]
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--output-prefix", output_prefix_text]
    if model is not None:
        model_path = _core.existing_directory(model, "model")
        arguments += ["--model", str(model_path)]
    if annotated_intervals is not None:
        annotated_intervals_path = _core.existing_file(annotated_intervals, "annotated_intervals")
        arguments += ["--annotated-intervals", str(annotated_intervals_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if max_copy_number is not None:
        max_copy_number_value = _core.integer(max_copy_number, "max_copy_number", minimum=1)
        arguments += ["--max-copy-number", str(max_copy_number_value)]
    if p_alt is not None:
        p_alt_value = _core.number(p_alt, "p_alt", minimum=0.0, maximum=1.0)
        arguments += ["--p-alt", str(p_alt_value)]
    if p_active is not None:
        p_active_value = _core.number(p_active, "p_active", minimum=0.0, maximum=1.0)
        arguments += ["--p-active", str(p_active_value)]
    if interval_merging_rule is not None:
        interval_merging_rule_value = _core.choice(
            interval_merging_rule, ("ALL", "OVERLAPPING_ONLY"), "interval_merging_rule"
        )
        arguments += ["--interval-merging-rule", interval_merging_rule_value]
    tool = "GermlineCNVCaller"
    reserved = _core.names(
        "annotated-intervals contig-ploidy-calls input interval-merging-rule "
        "intervals max-copy-number model output output-prefix p-active p-alt "
        "run-mode"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_germline_cnv_caller",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["calls_directories"] = sorted(str(item) for item in output_dir_path.glob("*-calls"))
    extras["model_directories"] = sorted(str(item) for item in output_dir_path.glob("*-model"))
    return _core.finish(result, outputs=outs, **extras)


def gatk_postprocess_germline_cnv_calls(
    calls_shard_paths: Sequence[str],
    model_shard_paths: Sequence[str],
    contig_ploidy_calls: str,
    output_genotyped_intervals: str,
    output_genotyped_segments: str,
    output_denoised_copy_ratios: str,
    *,
    sample_index: int | None = None,
    allosomal_contigs: Sequence[str] | None = None,
    autosomal_ref_copy_number: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Postprocess GermlineCNVCaller output of one sample into an intervals VCF, a segments VCF and
    denoised copy ratios.

    This tool runs GATK's Python module gcnvkernel, so the gatk launcher must run inside the
    GATK conda environment (gatkcondaenv.yml of the release zip).

    Parameters (the GATK option is in parentheses):

    - ``calls_shard_paths`` (``--calls-shard-path``): Calls folders written by GermlineCNVCaller
      (one per shard).
    - ``model_shard_paths`` (``--model-shard-path``): Model folders written by GermlineCNVCaller
      (one per shard).
    - ``contig_ploidy_calls`` (``--contig-ploidy-calls``): Calls folder written by
      DetermineGermlineContigPloidy.
    - ``output_genotyped_intervals`` (``--output-genotyped-intervals``): New intervals VCF.
    - ``output_genotyped_segments`` (``--output-genotyped-segments``): New segments VCF.
    - ``output_denoised_copy_ratios`` (``--output-denoised-copy-ratios``): New denoised copy-
      ratio file.
    - ``sample_index`` (``--sample-index``): Index of the sample in the call set.
    - ``allosomal_contigs`` (``--allosomal-contig``): Contigs treated as allosomal (sex
      chromosomes).
    - ``autosomal_ref_copy_number`` (``--autosomal-ref-copy-number``): Reference copy number on
      autosomes.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk PostprocessGermlineCNVCalls --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``segment_record_count``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    calls_shard_paths_paths = tuple(
        _core.existing_directory(item, "calls_shard_paths") for item in calls_shard_paths
    )
    for item in calls_shard_paths_paths:
        arguments += ["--calls-shard-path", str(item)]
    model_shard_paths_paths = tuple(
        _core.existing_directory(item, "model_shard_paths") for item in model_shard_paths
    )
    for item in model_shard_paths_paths:
        arguments += ["--model-shard-path", str(item)]
    contig_ploidy_calls_path = _core.existing_directory(contig_ploidy_calls, "contig_ploidy_calls")
    arguments += ["--contig-ploidy-calls", str(contig_ploidy_calls_path)]
    output_genotyped_intervals_path = _core.new_file(
        output_genotyped_intervals, "output_genotyped_intervals"
    )
    arguments += ["--output-genotyped-intervals", str(output_genotyped_intervals_path)]
    expected.append(output_genotyped_intervals_path)
    outs["output_genotyped_intervals"] = output_genotyped_intervals_path
    output_genotyped_segments_path = _core.new_file(
        output_genotyped_segments, "output_genotyped_segments"
    )
    arguments += ["--output-genotyped-segments", str(output_genotyped_segments_path)]
    expected.append(output_genotyped_segments_path)
    outs["output_genotyped_segments"] = output_genotyped_segments_path
    output_denoised_copy_ratios_path = _core.new_file(
        output_denoised_copy_ratios, "output_denoised_copy_ratios"
    )
    arguments += ["--output-denoised-copy-ratios", str(output_denoised_copy_ratios_path)]
    expected.append(output_denoised_copy_ratios_path)
    outs["output_denoised_copy_ratios"] = output_denoised_copy_ratios_path
    if sample_index is not None:
        sample_index_value = _core.integer(sample_index, "sample_index", minimum=0)
        arguments += ["--sample-index", str(sample_index_value)]
    if allosomal_contigs is not None:
        allosomal_contigs_values = _core.texts(allosomal_contigs, "allosomal_contigs")
        for item in allosomal_contigs_values:
            arguments += ["--allosomal-contig", item]
    if autosomal_ref_copy_number is not None:
        autosomal_ref_copy_number_value = _core.integer(
            autosomal_ref_copy_number, "autosomal_ref_copy_number", minimum=0
        )
        arguments += ["--autosomal-ref-copy-number", str(autosomal_ref_copy_number_value)]
    tool = "PostprocessGermlineCNVCalls"
    reserved = _core.names(
        "allosomal-contig autosomal-ref-copy-number calls-shard-path "
        "contig-ploidy-calls model-shard-path output-denoised-copy-ratios "
        "output-genotyped-intervals output-genotyped-segments sample-index"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_postprocess_germline_cnv_calls",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["segment_record_count"] = _core.count_data_lines(outs["output_genotyped_segments"])
    return _core.finish(result, outputs=outs, **extras)


def gatk_ase_read_counter(
    input_paths: Sequence[str],
    variant_paths: Sequence[str],
    output_path: str,
    *,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    output_format: str | None = None,
    min_base_quality: int | None = None,
    min_mapping_quality: int | None = None,
    min_depth_of_non_filtered_base: int | None = None,
    max_depth_per_sample: int | None = None,
    count_overlap_reads_handling: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a table of filtered base counts at heterozygous sites for allele-specific
    expression.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``variant_paths`` (``--variant``): VCF files with the heterozygous sites.
    - ``output_path`` (``--output``): New count table.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``output_format`` (``--output-format``): Table format; GATK default RTABLE.
    - ``min_base_quality`` (``--min-base-quality``): Minimum base quality.
    - ``min_mapping_quality`` (``--min-mapping-quality``): Minimum mapping quality.
    - ``min_depth_of_non_filtered_base`` (``--min-depth-of-non-filtered-base``): Minimum depth
      of non-filtered bases at a site.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``count_overlap_reads_handling`` (``--count-overlap-reads-handling``): How overlapping
      read pairs are counted.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk ASEReadCounter --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (contig, position, refCount, altCount ...; first 200
    rows, tab-separated formats).
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
    variant_paths_paths = _core.existing_files(variant_paths, "variant_paths")
    for item in variant_paths_paths:
        arguments += ["--variant", str(item)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if output_format is not None:
        output_format_value = _core.choice(
            output_format, ("TABLE", "RTABLE", "CSV"), "output_format"
        )
        arguments += ["--output-format", output_format_value]
    if min_base_quality is not None:
        min_base_quality_value = _core.integer(min_base_quality, "min_base_quality", minimum=0)
        arguments += ["--min-base-quality", str(min_base_quality_value)]
    if min_mapping_quality is not None:
        min_mapping_quality_value = _core.integer(
            min_mapping_quality, "min_mapping_quality", minimum=0
        )
        arguments += ["--min-mapping-quality", str(min_mapping_quality_value)]
    if min_depth_of_non_filtered_base is not None:
        min_depth_of_non_filtered_base_value = _core.integer(
            min_depth_of_non_filtered_base, "min_depth_of_non_filtered_base"
        )
        arguments += ["--min-depth-of-non-filtered-base", str(min_depth_of_non_filtered_base_value)]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    if count_overlap_reads_handling is not None:
        count_overlap_reads_handling_value = _core.choice(
            count_overlap_reads_handling,
            ("COUNT_READS", "COUNT_FRAGMENTS", "COUNT_FRAGMENTS_REQUIRE_SAME_BASE"),
            "count_overlap_reads_handling",
        )
        arguments += ["--count-overlap-reads-handling", count_overlap_reads_handling_value]
    tool = "ASEReadCounter"
    reserved = _core.names(
        "count-overlap-reads-handling input intervals max-depth-per-sample "
        "min-base-quality min-depth-of-non-filtered-base min-mapping-quality "
        "output output-format reference variant"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_ase_read_counter",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["rows"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_analyze_saturation_mutagenesis(
    input_paths: Sequence[str],
    reference: str,
    orf: str,
    output_file_prefix: str,
    *,
    min_q: int | None = None,
    min_mapq: int | None = None,
    min_length: int | None = None,
    min_flanking_length: int | None = None,
    min_variant_obs: int | None = None,
    paired_mode: bool | None = None,
    codon_translation: str | None = None,
    write_rejected_reads: bool | None = None,
    find_large_deletions: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Process reads from a MITESeq or other saturation mutagenesis experiment into per-variant
    counts of codon and amino-acid changes.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Aligned BAM/SAM/CRAM file(s) of the experiment.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``orf`` (``--orf``): Open reading frame as 'start-end' (1-based) or several ranges
      separated by commas.
    - ``output_file_prefix`` (``--output-file-prefix``): Prefix of the output files.
    - ``min_q`` (``--min-q``): Minimum base quality of a variant base.
    - ``min_mapq`` (``--min-mapq``): Minimum mapping quality.
    - ``min_length`` (``--min-length``): Minimum read length.
    - ``min_flanking_length`` (``--min-flanking-length``): Minimum flanking match length around
      a variant.
    - ``min_variant_obs`` (``--min-variant-obs``): Minimum observations of a variant to report
      it.
    - ``paired_mode`` (``--paired-mode``): Treat reads as pairs.
    - ``codon_translation`` (``--codon-translation``): Codon translation table name or string.
    - ``write_rejected_reads`` (``--write-rejected-reads``): Write the rejected reads to a BAM.
    - ``find_large_deletions`` (``--find-large-deletions``): Look for large deletions.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk AnalyzeSaturationMutagenesis --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``output_files`` (every file that starts with the prefix, e.g.
    variantCounts, codonCounts, refCoverage, readCounts).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    orf_text = _core.text(orf, "orf")
    arguments += ["--orf", orf_text]
    output_file_prefix_path = _core.new_prefix(output_file_prefix, "output_file_prefix")
    arguments += ["--output-file-prefix", str(output_file_prefix_path)]
    outs["output_file_prefix"] = output_file_prefix_path
    if min_q is not None:
        min_q_value = _core.integer(min_q, "min_q", minimum=0)
        arguments += ["--min-q", str(min_q_value)]
    if min_mapq is not None:
        min_mapq_value = _core.integer(min_mapq, "min_mapq", minimum=0)
        arguments += ["--min-mapq", str(min_mapq_value)]
    if min_length is not None:
        min_length_value = _core.integer(min_length, "min_length", minimum=0)
        arguments += ["--min-length", str(min_length_value)]
    if min_flanking_length is not None:
        min_flanking_length_value = _core.integer(
            min_flanking_length, "min_flanking_length", minimum=0
        )
        arguments += ["--min-flanking-length", str(min_flanking_length_value)]
    if min_variant_obs is not None:
        min_variant_obs_value = _core.integer(min_variant_obs, "min_variant_obs", minimum=0)
        arguments += ["--min-variant-obs", str(min_variant_obs_value)]
    if paired_mode is not None:
        arguments += [
            "--paired-mode",
            "true" if _core.flag(paired_mode, "paired_mode") else "false",
        ]
    if codon_translation is not None:
        codon_translation_text = _core.text(codon_translation, "codon_translation")
        arguments += ["--codon-translation", codon_translation_text]
    if write_rejected_reads is not None:
        arguments += [
            "--write-rejected-reads",
            "true" if _core.flag(write_rejected_reads, "write_rejected_reads") else "false",
        ]
    if find_large_deletions is not None:
        arguments += [
            "--find-large-deletions",
            "true" if _core.flag(find_large_deletions, "find_large_deletions") else "false",
        ]
    tool = "AnalyzeSaturationMutagenesis"
    reserved = _core.names(
        "codon-translation find-large-deletions input min-flanking-length "
        "min-length min-mapq min-q min-variant-obs orf output-file-prefix "
        "paired-mode reference write-rejected-reads"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_analyze_saturation_mutagenesis",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["output_files"] = sorted(
        str(item)
        for item in output_file_prefix_path.parent.glob(output_file_prefix_path.name + "*")
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_callable_loci(
    input_paths: Sequence[str],
    reference: str,
    output_path: str,
    summary_path: str,
    *,
    intervals: Sequence[str] | None = None,
    min_depth: int | None = None,
    max_depth: int | None = None,
    min_mapping_quality: int | None = None,
    min_base_quality: int | None = None,
    min_depth_for_low_mapq: int | None = None,
    max_low_mapq: int | None = None,
    max_fraction_of_reads_with_low_mapq: float | None = None,
    output_format: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Determine the callable status of every locus (CALLABLE, LOW_COVERAGE, NO_COVERAGE,
    POOR_MAPPING_QUALITY, EXCESSIVE_COVERAGE, REF_N).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New BED (or per-base state) file.
    - ``summary_path`` (``--summary``): New summary table (bases per state).
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``min_depth`` (``--min-depth``): Minimum depth of non-filtered reads for CALLABLE.
    - ``max_depth`` (``--max-depth``): Depth above which a locus is EXCESSIVE_COVERAGE.
    - ``min_mapping_quality`` (``--min-mapping-quality``): Minimum mapping quality of a counted
      read.
    - ``min_base_quality`` (``--min-base-quality``): Minimum base quality of a counted base.
    - ``min_depth_for_low_mapq`` (``--min-depth-for-low-mapq``): Minimum depth for the low-
      mapping-quality state.
    - ``max_low_mapq`` (``--max-low-mapq``): Mapping quality at or below this counts as low.
    - ``max_fraction_of_reads_with_low_mapq`` (``--max-fraction-of-reads-with-low-mapq``):
      Maximum fraction of low-mapping-quality reads at a CALLABLE locus.
    - ``output_format`` (``--format``): BED intervals or one state per base.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CallableLoci --help`` are accepted, a list repeats the option,
      and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``summary`` (bases per state).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    summary_path_path = _core.new_file(summary_path, "summary_path")
    arguments += ["--summary", str(summary_path_path)]
    expected.append(summary_path_path)
    outs["summary_path"] = summary_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if min_depth is not None:
        min_depth_value = _core.integer(min_depth, "min_depth", minimum=0)
        arguments += ["--min-depth", str(min_depth_value)]
    if max_depth is not None:
        max_depth_value = _core.integer(max_depth, "max_depth", minimum=0)
        arguments += ["--max-depth", str(max_depth_value)]
    if min_mapping_quality is not None:
        min_mapping_quality_value = _core.integer(
            min_mapping_quality, "min_mapping_quality", minimum=0
        )
        arguments += ["--min-mapping-quality", str(min_mapping_quality_value)]
    if min_base_quality is not None:
        min_base_quality_value = _core.integer(min_base_quality, "min_base_quality", minimum=0)
        arguments += ["--min-base-quality", str(min_base_quality_value)]
    if min_depth_for_low_mapq is not None:
        min_depth_for_low_mapq_value = _core.integer(
            min_depth_for_low_mapq, "min_depth_for_low_mapq", minimum=0
        )
        arguments += ["--min-depth-for-low-mapq", str(min_depth_for_low_mapq_value)]
    if max_low_mapq is not None:
        max_low_mapq_value = _core.integer(max_low_mapq, "max_low_mapq", minimum=0)
        arguments += ["--max-low-mapq", str(max_low_mapq_value)]
    if max_fraction_of_reads_with_low_mapq is not None:
        max_fraction_of_reads_with_low_mapq_value = _core.number(
            max_fraction_of_reads_with_low_mapq,
            "max_fraction_of_reads_with_low_mapq",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--max-fraction-of-reads-with-low-mapq",
            str(max_fraction_of_reads_with_low_mapq_value),
        ]
    if output_format is not None:
        output_format_value = _core.choice(
            output_format, ("BED", "STATE_PER_BASE"), "output_format"
        )
        arguments += ["--format", output_format_value]
    tool = "CallableLoci"
    reserved = _core.names(
        "format input intervals max-depth max-fraction-of-reads-with-low-mapq "
        "max-low-mapq min-base-quality min-depth min-depth-for-low-mapq "
        "min-mapping-quality output reference summary"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_callable_loci",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    summary_rows = {}
    for line in _core.read_text(outs["summary_path"]).splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[1].isdigit():
            summary_rows[parts[0]] = int(parts[1])
    extras["summary"] = summary_rows
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_f1r2_counts(
    input_paths: Sequence[str],
    reference: str,
    output_path: str,
    *,
    intervals: Sequence[str] | None = None,
    f1r2_min_bq: int | None = None,
    f1r2_median_mq: int | None = None,
    f1r2_max_depth: int | None = None,
    max_depth_per_sample: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect F1R2 read counts per sequence context for the Mutect2 orientation-bias mixture model
    filter.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New tar.gz with the F1R2 counts.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``f1r2_min_bq`` (``--f1r2-min-bq``): Minimum base quality.
    - ``f1r2_median_mq`` (``--f1r2-median-mq``): Minimum median mapping quality.
    - ``f1r2_max_depth`` (``--f1r2-max-depth``): Largest depth recorded in the histograms.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CollectF1R2Counts --help`` are accepted, a list repeats the
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if f1r2_min_bq is not None:
        f1r2_min_bq_value = _core.integer(f1r2_min_bq, "f1r2_min_bq", minimum=0)
        arguments += ["--f1r2-min-bq", str(f1r2_min_bq_value)]
    if f1r2_median_mq is not None:
        f1r2_median_mq_value = _core.integer(f1r2_median_mq, "f1r2_median_mq", minimum=0)
        arguments += ["--f1r2-median-mq", str(f1r2_median_mq_value)]
    if f1r2_max_depth is not None:
        f1r2_max_depth_value = _core.integer(f1r2_max_depth, "f1r2_max_depth", minimum=1)
        arguments += ["--f1r2-max-depth", str(f1r2_max_depth_value)]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    tool = "CollectF1R2Counts"
    reserved = _core.names(
        "f1r2-max-depth f1r2-median-mq f1r2-min-bq input intervals "
        "max-depth-per-sample output reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_f1r2_counts",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_count_reads_or_bases(
    input_paths: Sequence[str],
    output_path: str,
    *,
    unit: str | None = None,
    intervals: Sequence[str] | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the reads (unit='reads', CountReads) or the bases (unit='bases', CountBases) of a
    SAM/BAM/CRAM file after GATK's default read filters.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_path`` (``--output``): New file that receives the count.
    - ``unit`` (``--None``): What to count; GATK CountReads (default) or CountBases.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk <tool> --help`` are accepted, a list repeats the option, and
      options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``count`` and ``unit``.
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
    if unit is not None:
        _core.choice(unit, ("reads", "bases"), "unit")
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    tool = "CountBases" if unit == "bases" else "CountReads"
    reserved = _core.names("input intervals output reference")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_count_reads_or_bases",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["count"] = _core.first_integer(outs["output_path"])
    extras["unit"] = unit or "reads"
    return _core.finish(result, outputs=outs, **extras)


def gatk_depth_of_coverage(
    input_paths: Sequence[str],
    reference: str,
    intervals: Sequence[str],
    output_prefix: str,
    *,
    partition_types: Sequence[str] | None = None,
    summary_coverage_thresholds: Sequence[int] | None = None,
    min_base_quality: int | None = None,
    count_type: str | None = None,
    output_format: str | None = None,
    calculate_coverage_over_genes: str | None = None,
    omit_depth_output_at_each_base: bool | None = None,
    omit_interval_statistics: bool | None = None,
    omit_locus_table: bool | None = None,
    omit_per_sample_statistics: bool | None = None,
    print_base_counts: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate coverage summary information for reads data: depth per locus, per interval, per
    gene and per sample, with coverage histograms.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``intervals`` (``--intervals``): Loci or interval files to cover.
    - ``output_prefix`` (``--output``): Prefix of the output files (.sample_summary,
      .sample_interval_summary, ...).
    - ``partition_types`` (``--partition-type``): How reads are grouped; GATK default sample.
    - ``summary_coverage_thresholds`` (``--summary-coverage-threshold``): Depth thresholds for
      the percentage-above-threshold columns.
    - ``min_base_quality`` (``--min-base-quality``): Minimum base quality.
    - ``count_type`` (``--count-type``): How overlapping pairs are counted.
    - ``output_format`` (``--output-format``): Table format; GATK default CSV.
    - ``calculate_coverage_over_genes`` (``--calculate-coverage-over-genes``): RefSeq gene-list
      file for per-gene coverage.
    - ``omit_depth_output_at_each_base`` (``--omit-depth-output-at-each-base``): Do not write
      the per-locus depth file.
    - ``omit_interval_statistics`` (``--omit-interval-statistics``): Do not write the interval
      files.
    - ``omit_locus_table`` (``--omit-locus-table``): Do not write the locus histogram table.
    - ``omit_per_sample_statistics`` (``--omit-per-sample-statistics``): Do not write the per-
      sample files.
    - ``print_base_counts`` (``--print-base-counts``): Write the A/C/G/T/N counts at each locus.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk DepthOfCoverage --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``output_files`` and ``sample_summary`` (rows of
    <prefix>.sample_summary: total, mean, granular quartiles and percent bases above the
    thresholds).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
    arguments += ["--output", str(output_prefix_path)]
    outs["output_prefix"] = output_prefix_path
    if partition_types is not None:
        partition_types_values = _core.choices(
            partition_types,
            (
                "readgroup",
                "sample",
                "library",
                "platform",
                "center",
                "sample_by_platform",
                "sample_by_center",
                "sample_by_platform_by_center",
            ),
            "partition_types",
        )
        for item in partition_types_values:
            arguments += ["--partition-type", item]
    if summary_coverage_thresholds is not None:
        summary_coverage_thresholds_values = _core.integers(
            summary_coverage_thresholds, "summary_coverage_thresholds", minimum=0
        )
        for item in summary_coverage_thresholds_values:
            arguments += ["--summary-coverage-threshold", str(item)]
    if min_base_quality is not None:
        min_base_quality_value = _core.integer(min_base_quality, "min_base_quality", minimum=0)
        arguments += ["--min-base-quality", str(min_base_quality_value)]
    if count_type is not None:
        count_type_value = _core.choice(
            count_type,
            ("COUNT_READS", "COUNT_FRAGMENTS", "COUNT_FRAGMENTS_REQUIRE_SAME_BASE"),
            "count_type",
        )
        arguments += ["--count-type", count_type_value]
    if output_format is not None:
        output_format_value = _core.choice(output_format, ("TABLE", "CSV"), "output_format")
        arguments += ["--output-format", output_format_value]
    if calculate_coverage_over_genes is not None:
        calculate_coverage_over_genes_path = _core.existing_file(
            calculate_coverage_over_genes, "calculate_coverage_over_genes"
        )
        arguments += ["--calculate-coverage-over-genes", str(calculate_coverage_over_genes_path)]
    if omit_depth_output_at_each_base is not None:
        arguments += [
            "--omit-depth-output-at-each-base",
            "true"
            if _core.flag(omit_depth_output_at_each_base, "omit_depth_output_at_each_base")
            else "false",
        ]
    if omit_interval_statistics is not None:
        arguments += [
            "--omit-interval-statistics",
            "true" if _core.flag(omit_interval_statistics, "omit_interval_statistics") else "false",
        ]
    if omit_locus_table is not None:
        arguments += [
            "--omit-locus-table",
            "true" if _core.flag(omit_locus_table, "omit_locus_table") else "false",
        ]
    if omit_per_sample_statistics is not None:
        arguments += [
            "--omit-per-sample-statistics",
            "true"
            if _core.flag(omit_per_sample_statistics, "omit_per_sample_statistics")
            else "false",
        ]
    if print_base_counts is not None:
        arguments += [
            "--print-base-counts",
            "true" if _core.flag(print_base_counts, "print_base_counts") else "false",
        ]
    tool = "DepthOfCoverage"
    reserved = _core.names(
        "calculate-coverage-over-genes count-type input intervals "
        "min-base-quality omit-depth-output-at-each-base "
        "omit-interval-statistics omit-locus-table omit-per-sample-statistics "
        "output output-format partition-type print-base-counts reference "
        "summary-coverage-threshold"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_depth_of_coverage",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    files = sorted(
        str(item) for item in output_prefix_path.parent.glob(output_prefix_path.name + "*")
    )
    extras["output_files"] = files
    summary_file = Path(str(output_prefix_path) + ".sample_summary")
    extras["sample_summary"] = _core.delimited_rows(summary_file) if summary_file.is_file() else []
    return _core.finish(result, outputs=outs, **extras)


def gatk_gene_expression_evaluation(
    input_paths: Sequence[str],
    gff_file: str,
    output_path: str,
    *,
    reference: str | None = None,
    intervals: Sequence[str] | None = None,
    grouping_types: Sequence[str] | None = None,
    overlap_types: Sequence[str] | None = None,
    feature_label_key: str | None = None,
    multi_map_method: str | None = None,
    multi_overlap_method: str | None = None,
    read_strands: str | None = None,
    unspliced: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate gene expression from RNA-seq reads aligned to the genome: read counts per gene,
    exon or other feature of a GFF.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``gff_file`` (``--gff-file``): GFF/GTF annotation file.
    - ``output_path`` (``--output``): New expression table.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``grouping_types`` (``--grouping-type``): Feature types that define the expression units;
      GATK default gene.
    - ``overlap_types`` (``--overlap-type``): Feature types that reads must overlap; GATK
      default exon.
    - ``feature_label_key`` (``--feature-label-key``): Attribute used as the feature label
      (NAME, ID ...).
    - ``multi_map_method`` (``--multi-map-method``): How multi-mapped reads are handled.
    - ``multi_overlap_method`` (``--multi-overlap-method``): How reads overlapping several
      features are shared.
    - ``read_strands`` (``--read-strands``): Strand relationship of read and feature.
    - ``unspliced`` (``--unspliced``): Count the read span instead of the spliced blocks.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GeneExpressionEvaluation --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (first 200 rows of the expression table).
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
    gff_file_path = _core.existing_file(gff_file, "gff_file")
    arguments += ["--gff-file", str(gff_file_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if grouping_types is not None:
        grouping_types_values = _core.texts(grouping_types, "grouping_types")
        for item in grouping_types_values:
            arguments += ["--grouping-type", item]
    if overlap_types is not None:
        overlap_types_values = _core.texts(overlap_types, "overlap_types")
        for item in overlap_types_values:
            arguments += ["--overlap-type", item]
    if feature_label_key is not None:
        feature_label_key_text = _core.text(feature_label_key, "feature_label_key")
        arguments += ["--feature-label-key", feature_label_key_text]
    if multi_map_method is not None:
        multi_map_method_value = _core.choice(
            multi_map_method, ("IGNORE", "EQUAL"), "multi_map_method"
        )
        arguments += ["--multi-map-method", multi_map_method_value]
    if multi_overlap_method is not None:
        multi_overlap_method_value = _core.choice(
            multi_overlap_method, ("EQUAL", "PROPORTIONAL"), "multi_overlap_method"
        )
        arguments += ["--multi-overlap-method", multi_overlap_method_value]
    if read_strands is not None:
        read_strands_value = _core.choice(
            read_strands,
            ("FORWARD_FORWARD", "FORWARD_REVERSE", "REVERSE_FORWARD", "REVERSE_REVERSE"),
            "read_strands",
        )
        arguments += ["--read-strands", read_strands_value]
    if unspliced is not None:
        arguments += ["--unspliced", "true" if _core.flag(unspliced, "unspliced") else "false"]
    tool = "GeneExpressionEvaluation"
    reserved = _core.names(
        "feature-label-key gff-file grouping-type input intervals "
        "multi-map-method multi-overlap-method output overlap-type read-strands "
        "reference unspliced"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_gene_expression_evaluation",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["rows"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_get_normal_artifact_data(
    input_paths: Sequence[str],
    normal_samples: Sequence[str],
    reference: str,
    output_path: str,
    *,
    intervals: Sequence[str] | None = None,
    error_prob: float | None = None,
    max_depth_per_sample: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect data for training the normal-artifact filter of Mutect2 (allele counts of normal
    samples at candidate sites).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``normal_samples`` (``--normal-sample``): Sample names of the normal samples.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New training-data table.
    - ``intervals`` (``--intervals``): Loci (such as chr1:1-1000) or interval files to restrict
      the analysis.
    - ``error_prob`` (``--error-prob``): Per-base error probability.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GetNormalArtifactData --help`` are accepted, a list repeats
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
    normal_samples_values = _core.texts(normal_samples, "normal_samples")
    for item in normal_samples_values:
        arguments += ["--normal-sample", item]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if error_prob is not None:
        error_prob_value = _core.number(error_prob, "error_prob", minimum=0.0, maximum=1.0)
        arguments += ["--error-prob", str(error_prob_value)]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    tool = "GetNormalArtifactData"
    reserved = _core.names(
        "error-prob input intervals max-depth-per-sample normal-sample output reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_get_normal_artifact_data",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_get_pileup_summaries(
    input_paths: Sequence[str],
    variant: str,
    intervals: Sequence[str],
    output_path: str,
    *,
    minimum_population_allele_frequency: float | None = None,
    maximum_population_allele_frequency: float | None = None,
    max_depth_per_sample: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tabulate pileup metrics (ref, alt and other-alt counts at common biallelic sites) for
    inferring cross-sample contamination.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``variant`` (``--variant``): VCF of common biallelic variants with population allele
      frequencies (AF).
    - ``intervals`` (``--intervals``): Loci or interval files, usually the same variant VCF.
    - ``output_path`` (``--output``): New pileup-summaries table.
    - ``minimum_population_allele_frequency`` (``--minimum-population-allele-frequency``):
      Minimum allele frequency of a site.
    - ``maximum_population_allele_frequency`` (``--maximum-population-allele-frequency``):
      Maximum allele frequency of a site.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GetPileupSummaries --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``rows`` (contig, position, ref_count, alt_count,
    other_alt_count, allele_frequency; first 200).
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
    variant_path = _core.existing_file(variant, "variant")
    arguments += ["--variant", str(variant_path)]
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if minimum_population_allele_frequency is not None:
        minimum_population_allele_frequency_value = _core.number(
            minimum_population_allele_frequency,
            "minimum_population_allele_frequency",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--minimum-population-allele-frequency",
            str(minimum_population_allele_frequency_value),
        ]
    if maximum_population_allele_frequency is not None:
        maximum_population_allele_frequency_value = _core.number(
            maximum_population_allele_frequency,
            "maximum_population_allele_frequency",
            minimum=0.0,
            maximum=1.0,
        )
        arguments += [
            "--maximum-population-allele-frequency",
            str(maximum_population_allele_frequency_value),
        ]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    tool = "GetPileupSummaries"
    reserved = _core.names(
        "input intervals max-depth-per-sample "
        "maximum-population-allele-frequency "
        "minimum-population-allele-frequency output reference variant"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_get_pileup_summaries",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["rows"] = _core.tsv_rows(outs["output_path"], limit=200)
    return _core.finish(result, outputs=outs, **extras)


def gatk_local_assembler(
    input_paths: Sequence[str],
    intervals: Sequence[str],
    assembly_name: str,
    *,
    fasta_file: str | None = None,
    gfa_file: str | None = None,
    reference: str | None = None,
    min_sv_size: int | None = None,
    no_scaffolding: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assemble small regions locally to discover structural variants (GATK BETA tool) and write
    the scaffolds or the assembly graph.

    GATK marks this tool BETA (work in progress).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``intervals`` (``--intervals``): Regions to assemble.
    - ``assembly_name`` (``--assembly-name``): Name used as a prefix of the traversal names.
    - ``fasta_file`` (``--fasta-file``): New FASTA with the scaffolds.
    - ``gfa_file`` (``--gfa-file``): New GFA with the assembly graph.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files; needed for
      CRAM input.
    - ``min_sv_size`` (``--min-sv-size``): Smallest structural variant of interest.
    - ``no_scaffolding`` (``--no-scaffolding``): Skip scaffolding.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk LocalAssembler --help`` are accepted, a list repeats the
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
    if fasta_file is None and gfa_file is None:
        raise ConfigurationError(
            "Give fasta_file, gfa_file or both; otherwise the assembly would not be saved."
        )
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    for item in input_paths_paths:
        arguments += ["--input", str(item)]
    intervals_values = _core.intervals(intervals, "intervals")
    for item in intervals_values:
        arguments += ["--intervals", item]
    assembly_name_text = _core.text(assembly_name, "assembly_name")
    arguments += ["--assembly-name", assembly_name_text]
    if fasta_file is not None:
        fasta_file_path = _core.new_file(fasta_file, "fasta_file")
        arguments += ["--fasta-file", str(fasta_file_path)]
        expected.append(fasta_file_path)
        outs["fasta_file"] = fasta_file_path
    if gfa_file is not None:
        gfa_file_path = _core.new_file(gfa_file, "gfa_file")
        arguments += ["--gfa-file", str(gfa_file_path)]
        expected.append(gfa_file_path)
        outs["gfa_file"] = gfa_file_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    if min_sv_size is not None:
        min_sv_size_value = _core.integer(min_sv_size, "min_sv_size", minimum=1)
        arguments += ["--min-sv-size", str(min_sv_size_value)]
    if no_scaffolding is not None:
        arguments += [
            "--no-scaffolding",
            "true" if _core.flag(no_scaffolding, "no_scaffolding") else "false",
        ]
    tool = "LocalAssembler"
    reserved = _core.names(
        "assembly-name fasta-file gfa-file input intervals min-sv-size no-scaffolding reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_local_assembler",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_analyze_covariates(
    *,
    before_report_file: str | None = None,
    after_report_file: str | None = None,
    bqsr_recal_file: str | None = None,
    plots_report_file: str | None = None,
    intermediate_csv_file: str | None = None,
    ignore_last_modification_times: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate and compare base quality score recalibration (BQSR) tables and plot them.

    The plots need R (Rscript) with the packages of the GATK conda environment; without it GATK
    stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``before_report_file`` (``--before-report-file``): BQSR first-pass report (before
      recalibration).
    - ``after_report_file`` (``--after-report-file``): BQSR second-pass report (after
      recalibration).
    - ``bqsr_recal_file`` (``--bqsr-recal-file``): Recalibration table used for on-the-fly
      recalibration.
    - ``plots_report_file`` (``--plots-report-file``): New PDF with the plots.
    - ``intermediate_csv_file`` (``--intermediate-csv-file``): New intermediate CSV.
    - ``ignore_last_modification_times`` (``--ignore-last-modification-times``): Do not warn
      about suspicious file modification order.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk AnalyzeCovariates --help`` are accepted, a list repeats the
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
    if before_report_file is None and after_report_file is None and bqsr_recal_file is None:
        raise ConfigurationError(
            "Give at least one of before_report_file, after_report_file or bqsr_recal_file."
        )
    if before_report_file is not None:
        before_report_file_path = _core.existing_file(before_report_file, "before_report_file")
        arguments += ["--before-report-file", str(before_report_file_path)]
    if after_report_file is not None:
        after_report_file_path = _core.existing_file(after_report_file, "after_report_file")
        arguments += ["--after-report-file", str(after_report_file_path)]
    if bqsr_recal_file is not None:
        bqsr_recal_file_path = _core.existing_file(bqsr_recal_file, "bqsr_recal_file")
        arguments += ["--bqsr-recal-file", str(bqsr_recal_file_path)]
    if plots_report_file is not None:
        plots_report_file_path = _core.new_file(plots_report_file, "plots_report_file")
        arguments += ["--plots-report-file", str(plots_report_file_path)]
        expected.append(plots_report_file_path)
        outs["plots_report_file"] = plots_report_file_path
    if intermediate_csv_file is not None:
        intermediate_csv_file_path = _core.new_file(intermediate_csv_file, "intermediate_csv_file")
        arguments += ["--intermediate-csv-file", str(intermediate_csv_file_path)]
        expected.append(intermediate_csv_file_path)
        outs["intermediate_csv_file"] = intermediate_csv_file_path
    if ignore_last_modification_times is not None:
        arguments += [
            "--ignore-last-modification-times",
            "true"
            if _core.flag(ignore_last_modification_times, "ignore_last_modification_times")
            else "false",
        ]
    tool = "AnalyzeCovariates"
    reserved = _core.names(
        "after-report-file before-report-file bqsr-recal-file "
        "ignore-last-modification-times intermediate-csv-file plots-report-file"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_analyze_covariates",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)
