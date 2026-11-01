"""GATK command-line facades, part 1: Illumina runs, alignment and read metrics, fingerprints.

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

This part covers the Illumina basecall tools, the Picard collectors that GATK bundles (alignment,
hybrid-selection, WGS, RNA-seq, RRBS, artifact, UMI and variant-calling metrics), fingerprint
tools, SAM/BAM validation, comparison and conversion tools.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics import _gatk_p00_core as _core
from dnakit.exceptions import BackendExecutionError, ConfigurationError

__all__ = [
    "gatk_calculate_fingerprint_metrics",
    "gatk_calculate_read_group_checksum",
    "gatk_check_duplicate_marking",
    "gatk_check_fingerprint",
    "gatk_check_pileup",
    "gatk_check_terminator_block",
    "gatk_collect_alignment_summary_metrics",
    "gatk_collect_base_distribution_by_cycle",
    "gatk_collect_hiseqx_pf_fail_metrics",
    "gatk_collect_hs_metrics",
    "gatk_collect_illumina_basecalling_metrics",
    "gatk_collect_illumina_lane_metrics",
    "gatk_collect_independent_replicate_metrics",
    "gatk_collect_jumping_library_metrics",
    "gatk_collect_multiple_metrics",
    "gatk_collect_oxog_metrics",
    "gatk_collect_quality_yield_metrics",
    "gatk_collect_rnaseq_metrics",
    "gatk_collect_rrbs_metrics",
    "gatk_collect_sequencing_artifact_metrics",
    "gatk_collect_targeted_pcr_metrics",
    "gatk_collect_umi_prevalence_metrics",
    "gatk_collect_variant_calling_metrics",
    "gatk_collect_wgs_metrics",
    "gatk_compare_base_qualities",
    "gatk_compare_duplicates_spark",
    "gatk_compare_sams",
    "gatk_convert_haplotype_database_to_vcf",
    "gatk_crosscheck_fingerprints",
    "gatk_dump_tabix_index",
    "gatk_extract_fingerprint",
    "gatk_extract_illumina_barcodes",
    "gatk_get_sample_name",
    "gatk_htsget_download",
    "gatk_identify_contaminant",
    "gatk_illumina_basecalls_to_reads",
    "gatk_lift_over_haplotype_map",
    "gatk_mark_illumina_adapters",
    "gatk_mean_quality_by_cycle",
    "gatk_quality_score_distribution",
    "gatk_split_cram",
    "gatk_validate_sam_file",
    "gatk_view_sam",
]


def gatk_collect_illumina_basecalling_metrics(
    basecalls_dir: str,
    lane: int,
    read_structure: str,
    *,
    output_path: str | None = None,
    barcodes_dir: str | None = None,
    barcode_file: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect Illumina basecalling metrics (clusters, PF clusters, per-barcode counts) for one
    lane of a sequencing run.

    Needs the BaseCalls directory of an Illumina run (BCL/CBCL, filter, locs or clocs files).

    Parameters (the GATK option is in parentheses):

    - ``basecalls_dir`` (``--BASECALLS_DIR``): Illumina basecalls directory
      (Data/Intensities/BaseCalls) that holds the run's cluster files.
    - ``lane`` (``--LANE``): Lane number to read.
    - ``read_structure`` (``--READ_STRUCTURE``): Read structure such as '76T8B76T' (cycle count
      plus T template, B barcode, S skip, M molecular index).
    - ``output_path`` (``--OUTPUT``): New metrics file; None writes the metrics next to the
      barcode files as GATK does by default.
    - ``barcodes_dir`` (``--BARCODES_DIR``): Folder with the *_barcode.txt files made by
      ExtractIlluminaBarcodes; None uses basecalls_dir.
    - ``barcode_file`` (``--INPUT``): Barcode data file listing the barcodes to expect from the
      run.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectIlluminaBasecallingMetrics --help`` are
      accepted, a list repeats the option, and options already set by a parameter above are
      refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (rows of the metrics table when output_path is given)
    and ``histogram``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    basecalls_dir_path = _core.existing_directory(basecalls_dir, "basecalls_dir")
    arguments += ["--BASECALLS_DIR", str(basecalls_dir_path)]
    lane_value = _core.integer(lane, "lane", minimum=1)
    arguments += ["--LANE", str(lane_value)]
    read_structure_text = _core.text(read_structure, "read_structure")
    arguments += ["--READ_STRUCTURE", read_structure_text]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if barcodes_dir is not None:
        barcodes_dir_path = _core.existing_directory(barcodes_dir, "barcodes_dir")
        arguments += ["--BARCODES_DIR", str(barcodes_dir_path)]
    if barcode_file is not None:
        barcode_file_path = _core.existing_file(barcode_file, "barcode_file")
        arguments += ["--INPUT", str(barcode_file_path)]
    tool = "CollectIlluminaBasecallingMetrics"
    reserved = _core.names("BARCODES_DIR BASECALLS_DIR INPUT LANE OUTPUT READ_STRUCTURE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_illumina_basecalling_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_illumina_lane_metrics(
    run_directory: str,
    output_directory: str,
    output_prefix: str,
    *,
    read_structure: str | None = None,
    file_extension: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect Illumina lane and phasing metrics (cluster density, PF clusters, phasing and
    prephasing) from a run directory.

    Needs the InterOp files of an Illumina run.

    Parameters (the GATK option is in parentheses):

    - ``run_directory`` (``--RUN_DIRECTORY``): Illumina run directory (the folder that holds
      Data/ and InterOp/).
    - ``output_directory`` (``--OUTPUT_DIRECTORY``): Existing folder that receives the metric
      files.
    - ``output_prefix`` (``--OUTPUT_PREFIX``): Prefix of the output file names; GATK adds
      .illumina_lane_metrics and .illumina_phasing_metrics.
    - ``read_structure`` (``--READ_STRUCTURE``): Read structure such as '76T8B76T'; None reads
      it from the run.
    - ``file_extension`` (``--FILE_EXTENSION``): Extension appended to every metrics file name.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectIlluminaLaneMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``metrics_files`` (the .illumina_lane_metrics and .illumina_phasing_metrics
    files written) and ``outputs``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    run_directory_path = _core.existing_directory(run_directory, "run_directory")
    arguments += ["--RUN_DIRECTORY", str(run_directory_path)]
    output_directory_path = _core.existing_directory(output_directory, "output_directory")
    arguments += ["--OUTPUT_DIRECTORY", str(output_directory_path)]
    output_prefix_text = _core.text(output_prefix, "output_prefix")
    arguments += ["--OUTPUT_PREFIX", output_prefix_text]
    if read_structure is not None:
        read_structure_text = _core.text(read_structure, "read_structure")
        arguments += ["--READ_STRUCTURE", read_structure_text]
    if file_extension is not None:
        file_extension_text = _core.text(file_extension, "file_extension")
        arguments += ["--FILE_EXTENSION", file_extension_text]
    tool = "CollectIlluminaLaneMetrics"
    reserved = _core.names(
        "FILE_EXTENSION OUTPUT_DIRECTORY OUTPUT_PREFIX READ_STRUCTURE RUN_DIRECTORY"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_illumina_lane_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    found = sorted(
        str(item) for item in output_directory_path.glob(output_prefix_text + ".illumina_*")
    )
    extras["metrics_files"] = found
    if not found:
        raise BackendExecutionError(
            "gatk CollectIlluminaLaneMetrics wrote no metrics file.", code="GATK_OUTPUT_MISSING"
        )
    outs["lane_metrics"] = Path(found[0])
    return _core.finish(result, outputs=outs, **extras)


def gatk_extract_illumina_barcodes(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    *,
    barcodes: Sequence[str] | None = None,
    barcode_file: str | None = None,
    metrics_file: str | None = None,
    output_dir: str | None = None,
    max_mismatches: int | None = None,
    min_mismatch_delta: int | None = None,
    max_no_calls: int | None = None,
    minimum_base_quality: int | None = None,
    minimum_quality: int | None = None,
    distance_mode: str | None = None,
    compress_outputs: bool | None = None,
    num_processors: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Determine the barcode of every read of an Illumina lane and write the *_barcode.txt files
    and per-barcode metrics.

    Needs the BaseCalls directory of an Illumina run.

    Parameters (the GATK option is in parentheses):

    - ``basecalls_dir`` (``--BASECALLS_DIR``): Illumina basecalls directory
      (Data/Intensities/BaseCalls).
    - ``lanes`` (``--LANE``): Lane numbers; reads with the same index in several lanes go to the
      same file.
    - ``read_structure`` (``--READ_STRUCTURE``): Read structure such as '76T8B76T'.
    - ``barcodes`` (``--BARCODE``): Barcode sequences (all the same length, unique); give either
      barcodes or barcode_file.
    - ``barcode_file`` (``--BARCODE_FILE``): Tab-delimited barcode file (barcode_sequence,
      barcode_name, library_name).
    - ``metrics_file`` (``--METRICS_FILE``): New per-barcode and per-lane metrics file.
    - ``output_dir`` (``--OUTPUT_DIR``): New folder for the *_barcode.txt files; None writes
      them into basecalls_dir.
    - ``max_mismatches`` (``--MAX_MISMATCHES``): Maximum mismatches for a barcode to be
      assigned.
    - ``min_mismatch_delta`` (``--MIN_MISMATCH_DELTA``): Minimum difference between the best and
      second-best barcode.
    - ``max_no_calls`` (``--MAX_NO_CALLS``): Maximum no-calls (N) allowed in a barcode read.
    - ``minimum_base_quality`` (``--MINIMUM_BASE_QUALITY``): Bases below this quality count as
      mismatches.
    - ``minimum_quality`` (``--MINIMUM_QUALITY``): Quality assigned to bases that are no-calls.
    - ``distance_mode`` (``--DISTANCE_MODE``): Distance metric used to compare barcode reads
      with the barcodes.
    - ``compress_outputs`` (``--COMPRESS_OUTPUTS``): Gzip the barcode files.
    - ``num_processors`` (``--NUM_PROCESSORS``): Threads; 0 or negative counts from the number
      of cores.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk ExtractIlluminaBarcodes --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (one row per barcode when metrics_file is given) and
    ``histogram``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if (barcodes is None) == (barcode_file is None):
        raise ConfigurationError("Give exactly one of barcodes or barcode_file.")
    basecalls_dir_path = _core.existing_directory(basecalls_dir, "basecalls_dir")
    arguments += ["--BASECALLS_DIR", str(basecalls_dir_path)]
    lanes_values = _core.integers(lanes, "lanes", minimum=1)
    for item in lanes_values:
        arguments += ["--LANE", str(item)]
    read_structure_text = _core.text(read_structure, "read_structure")
    arguments += ["--READ_STRUCTURE", read_structure_text]
    if barcodes is not None:
        barcodes_values = _core.texts(barcodes, "barcodes")
        for item in barcodes_values:
            arguments += ["--BARCODE", item]
    if barcode_file is not None:
        barcode_file_path = _core.existing_file(barcode_file, "barcode_file")
        arguments += ["--BARCODE_FILE", str(barcode_file_path)]
    if metrics_file is not None:
        metrics_file_path = _core.new_file(metrics_file, "metrics_file")
        arguments += ["--METRICS_FILE", str(metrics_file_path)]
        expected.append(metrics_file_path)
        outs["metrics_file"] = metrics_file_path
    if output_dir is not None:
        output_dir_path = _core.new_directory(output_dir, "output_dir")
        arguments += ["--OUTPUT_DIR", str(output_dir_path)]
        made_directories.append(output_dir_path)
        outs["output_dir"] = output_dir_path
    if max_mismatches is not None:
        max_mismatches_value = _core.integer(max_mismatches, "max_mismatches", minimum=0)
        arguments += ["--MAX_MISMATCHES", str(max_mismatches_value)]
    if min_mismatch_delta is not None:
        min_mismatch_delta_value = _core.integer(
            min_mismatch_delta, "min_mismatch_delta", minimum=0
        )
        arguments += ["--MIN_MISMATCH_DELTA", str(min_mismatch_delta_value)]
    if max_no_calls is not None:
        max_no_calls_value = _core.integer(max_no_calls, "max_no_calls", minimum=0)
        arguments += ["--MAX_NO_CALLS", str(max_no_calls_value)]
    if minimum_base_quality is not None:
        minimum_base_quality_value = _core.integer(
            minimum_base_quality, "minimum_base_quality", minimum=0
        )
        arguments += ["--MINIMUM_BASE_QUALITY", str(minimum_base_quality_value)]
    if minimum_quality is not None:
        minimum_quality_value = _core.integer(minimum_quality, "minimum_quality", minimum=0)
        arguments += ["--MINIMUM_QUALITY", str(minimum_quality_value)]
    if distance_mode is not None:
        distance_mode_value = _core.choice(
            distance_mode, ("HAMMING", "LENIENT_HAMMING", "FREE"), "distance_mode"
        )
        arguments += ["--DISTANCE_MODE", distance_mode_value]
    if compress_outputs is not None:
        arguments += [
            "--COMPRESS_OUTPUTS",
            "true" if _core.flag(compress_outputs, "compress_outputs") else "false",
        ]
    if num_processors is not None:
        num_processors_value = _core.integer(num_processors, "num_processors")
        arguments += ["--NUM_PROCESSORS", str(num_processors_value)]
    tool = "ExtractIlluminaBarcodes"
    reserved = _core.names(
        "BARCODE BARCODE_FILE BASECALLS_DIR COMPRESS_OUTPUTS DISTANCE_MODE LANE "
        "MAX_MISMATCHES MAX_NO_CALLS METRICS_FILE MINIMUM_BASE_QUALITY "
        "MINIMUM_QUALITY MIN_MISMATCH_DELTA NUM_PROCESSORS OUTPUT_DIR "
        "READ_STRUCTURE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_extract_illumina_barcodes",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "metrics_file" in outs and outs["metrics_file"].is_file():
        parsed = _core.picard_metrics(outs["metrics_file"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_illumina_basecalls_to_reads(
    basecalls_dir: str,
    lanes: Sequence[int],
    read_structure: str,
    run_barcode: str,
    parameters_file: str,
    output_format: str,
    *,
    output_prefix: str | None = None,
    sequencing_center: str | None = None,
    metrics_file: str | None = None,
    barcodes_dir: str | None = None,
    include_non_pf_reads: bool | None = None,
    max_mismatches: int | None = None,
    min_mismatch_delta: int | None = None,
    max_no_calls: int | None = None,
    minimum_quality: int | None = None,
    first_tile: int | None = None,
    tile_limit: int | None = None,
    num_processors: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert Illumina basecalls of one run to FASTQ files (output_format='fastq') or to unmapped
    SAM/BAM/CRAM files (output_format='sam').

    Needs the BaseCalls directory of an Illumina run.

    Parameters (the GATK option is in parentheses):

    - ``basecalls_dir`` (``--BASECALLS_DIR``): Illumina basecalls directory
      (Data/Intensities/BaseCalls).
    - ``lanes`` (``--LANE``): Lane numbers to convert.
    - ``read_structure`` (``--READ_STRUCTURE``): Read structure such as '76T8B76T'.
    - ``run_barcode`` (``--RUN_BARCODE``): Run barcode, used in read names and read groups.
    - ``parameters_file`` (``--LIBRARY_PARAMS/MULTIPLEX_PARAMS``): Tab-delimited file that names
      the output per barcode: MULTIPLEX_PARAMS for 'fastq', LIBRARY_PARAMS for 'sam' (columns
      OUTPUT, SAMPLE_ALIAS, LIBRARY_NAME, BARCODE_1...).
    - ``output_format`` (``--None``): 'fastq' runs IlluminaBasecallsToFastq, 'sam' runs
      IlluminaBasecallsToSam.
    - ``output_prefix`` (``--OUTPUT_PREFIX``): File-name prefix of the FASTQ files ('fastq'
      only).
    - ``sequencing_center`` (``--SEQUENCING_CENTER``): Sequencing centre written to the read
      groups ('sam' only).
    - ``metrics_file`` (``--METRICS_FILE``): New barcode metrics file.
    - ``barcodes_dir`` (``--BARCODES_DIR``): Folder with the *_barcode.txt files; None uses
      basecalls_dir.
    - ``include_non_pf_reads`` (``--INCLUDE_NON_PF_READS``): Keep reads that fail the vendor
      filter.
    - ``max_mismatches`` (``--MAX_MISMATCHES``): Maximum barcode mismatches.
    - ``min_mismatch_delta`` (``--MIN_MISMATCH_DELTA``): Minimum best/second-best barcode
      difference.
    - ``max_no_calls`` (``--MAX_NO_CALLS``): Maximum no-calls in a barcode read.
    - ``minimum_quality`` (``--MINIMUM_QUALITY``): Quality assigned to no-call bases.
    - ``first_tile`` (``--FIRST_TILE``): First tile to convert.
    - ``tile_limit`` (``--TILE_LIMIT``): Number of tiles to convert.
    - ``num_processors`` (``--NUM_PROCESSORS``): Threads; 0 uses all cores.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk <tool> --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``output_format``, ``metrics`` (barcode metrics when metrics_file
    is given). The FASTQ or BAM files are named by parameters_file.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    _core.choice(output_format, ("fastq", "sam"), "output_format")
    basecalls_dir_path = _core.existing_directory(basecalls_dir, "basecalls_dir")
    arguments += ["--BASECALLS_DIR", str(basecalls_dir_path)]
    lanes_values = _core.integers(lanes, "lanes", minimum=1)
    for item in lanes_values:
        arguments += ["--LANE", str(item)]
    read_structure_text = _core.text(read_structure, "read_structure")
    arguments += ["--READ_STRUCTURE", read_structure_text]
    run_barcode_text = _core.text(run_barcode, "run_barcode")
    arguments += ["--RUN_BARCODE", run_barcode_text]
    parameters_file_path = _core.existing_file(parameters_file, "parameters_file")
    arguments += [
        "--" + {"fastq": "MULTIPLEX_PARAMS", "sam": "LIBRARY_PARAMS"}[output_format],
        str(parameters_file_path),
    ]
    if output_prefix is not None and output_format not in ("fastq",):
        raise ConfigurationError("output_prefix only applies to output_format in ['fastq'].")
    if output_prefix is None and output_format in ("fastq",):
        raise ConfigurationError("output_prefix is required when output_format is in ['fastq'].")
    if output_prefix is not None:
        output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
        arguments += ["--OUTPUT_PREFIX", str(output_prefix_path)]
        outs["output_prefix"] = output_prefix_path
    if sequencing_center is not None and output_format not in ("sam",):
        raise ConfigurationError("sequencing_center only applies to output_format in ['sam'].")
    if sequencing_center is None and output_format in ("sam",):
        raise ConfigurationError("sequencing_center is required when output_format is in ['sam'].")
    if sequencing_center is not None:
        sequencing_center_text = _core.text(sequencing_center, "sequencing_center")
        arguments += ["--SEQUENCING_CENTER", sequencing_center_text]
    if metrics_file is not None:
        metrics_file_path = _core.new_file(metrics_file, "metrics_file")
        arguments += ["--METRICS_FILE", str(metrics_file_path)]
        expected.append(metrics_file_path)
        outs["metrics_file"] = metrics_file_path
    if barcodes_dir is not None:
        barcodes_dir_path = _core.existing_directory(barcodes_dir, "barcodes_dir")
        arguments += ["--BARCODES_DIR", str(barcodes_dir_path)]
    if include_non_pf_reads is not None:
        arguments += [
            "--INCLUDE_NON_PF_READS",
            "true" if _core.flag(include_non_pf_reads, "include_non_pf_reads") else "false",
        ]
    if max_mismatches is not None:
        max_mismatches_value = _core.integer(max_mismatches, "max_mismatches", minimum=0)
        arguments += ["--MAX_MISMATCHES", str(max_mismatches_value)]
    if min_mismatch_delta is not None:
        min_mismatch_delta_value = _core.integer(
            min_mismatch_delta, "min_mismatch_delta", minimum=0
        )
        arguments += ["--MIN_MISMATCH_DELTA", str(min_mismatch_delta_value)]
    if max_no_calls is not None:
        max_no_calls_value = _core.integer(max_no_calls, "max_no_calls", minimum=0)
        arguments += ["--MAX_NO_CALLS", str(max_no_calls_value)]
    if minimum_quality is not None:
        minimum_quality_value = _core.integer(minimum_quality, "minimum_quality", minimum=0)
        arguments += ["--MINIMUM_QUALITY", str(minimum_quality_value)]
    if first_tile is not None:
        first_tile_value = _core.integer(first_tile, "first_tile", minimum=0)
        arguments += ["--FIRST_TILE", str(first_tile_value)]
    if tile_limit is not None:
        tile_limit_value = _core.integer(tile_limit, "tile_limit", minimum=0)
        arguments += ["--TILE_LIMIT", str(tile_limit_value)]
    if num_processors is not None:
        num_processors_value = _core.integer(num_processors, "num_processors")
        arguments += ["--NUM_PROCESSORS", str(num_processors_value)]
    tool = {"fastq": "IlluminaBasecallsToFastq", "sam": "IlluminaBasecallsToSam"}[output_format]
    reserved = _core.names(
        "BARCODES_DIR BASECALLS_DIR FIRST_TILE INCLUDE_NON_PF_READS LANE "
        "LIBRARY_PARAMS MAX_MISMATCHES MAX_NO_CALLS METRICS_FILE "
        "MINIMUM_QUALITY MIN_MISMATCH_DELTA MULTIPLEX_PARAMS NUM_PROCESSORS "
        "OUTPUT_PREFIX READ_STRUCTURE RUN_BARCODE SEQUENCING_CENTER TILE_LIMIT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_illumina_basecalls_to_reads",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "metrics_file" in outs and outs["metrics_file"].is_file():
        parsed = _core.picard_metrics(outs["metrics_file"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_mark_illumina_adapters(
    input_path: str,
    metrics_file: str,
    *,
    output_path: str | None = None,
    adapters: Sequence[str] | None = None,
    five_prime_adapter: str | None = None,
    three_prime_adapter: str | None = None,
    max_error_rate_se: float | None = None,
    max_error_rate_pe: float | None = None,
    min_match_bases_se: int | None = None,
    min_match_bases_pe: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find Illumina adapter sequences in reads and rewrite the file with XT adapter-trimming tags,
    plus an adapter-clipping histogram.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Unaligned or aligned SAM/BAM/CRAM file.
    - ``metrics_file`` (``--METRICS``): New file for the histogram of clipped bases per read.
    - ``output_path`` (``--OUTPUT``): New BAM with the XT tags; None writes only the metrics.
    - ``adapters`` (``--ADAPTERS``): Adapter sets to search for.
    - ``five_prime_adapter`` (``--FIVE_PRIME_ADAPTER``): Custom 5' adapter sequence.
    - ``three_prime_adapter`` (``--THREE_PRIME_ADAPTER``): Custom 3' adapter sequence.
    - ``max_error_rate_se`` (``--MAX_ERROR_RATE_SE``): Maximum error rate for single-end
      matches.
    - ``max_error_rate_pe`` (``--MAX_ERROR_RATE_PE``): Maximum error rate for paired-end
      matches.
    - ``min_match_bases_se`` (``--MIN_MATCH_BASES_SE``): Minimum matching bases for single-end
      reads.
    - ``min_match_bases_pe`` (``--MIN_MATCH_BASES_PE``): Minimum matching bases for paired-end
      reads.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk MarkIlluminaAdapters --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` and ``histogram`` (bases clipped versus read count).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    metrics_file_path = _core.new_file(metrics_file, "metrics_file")
    arguments += ["--METRICS", str(metrics_file_path)]
    expected.append(metrics_file_path)
    outs["metrics_file"] = metrics_file_path
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if adapters is not None:
        adapters_values = _core.choices(
            adapters,
            (
                "PAIRED_END",
                "INDEXED",
                "SINGLE_END",
                "NEXTERA_V1",
                "NEXTERA_V2",
                "DUAL_INDEXED",
                "FLUIDIGM",
                "TRUSEQ_SMALLRNA",
                "ALTERNATIVE_SINGLE_END",
            ),
            "adapters",
        )
        for item in adapters_values:
            arguments += ["--ADAPTERS", item]
    if five_prime_adapter is not None:
        five_prime_adapter_text = _core.text(five_prime_adapter, "five_prime_adapter")
        arguments += ["--FIVE_PRIME_ADAPTER", five_prime_adapter_text]
    if three_prime_adapter is not None:
        three_prime_adapter_text = _core.text(three_prime_adapter, "three_prime_adapter")
        arguments += ["--THREE_PRIME_ADAPTER", three_prime_adapter_text]
    if max_error_rate_se is not None:
        max_error_rate_se_value = _core.number(
            max_error_rate_se, "max_error_rate_se", minimum=0.0, maximum=1.0
        )
        arguments += ["--MAX_ERROR_RATE_SE", str(max_error_rate_se_value)]
    if max_error_rate_pe is not None:
        max_error_rate_pe_value = _core.number(
            max_error_rate_pe, "max_error_rate_pe", minimum=0.0, maximum=1.0
        )
        arguments += ["--MAX_ERROR_RATE_PE", str(max_error_rate_pe_value)]
    if min_match_bases_se is not None:
        min_match_bases_se_value = _core.integer(
            min_match_bases_se, "min_match_bases_se", minimum=1
        )
        arguments += ["--MIN_MATCH_BASES_SE", str(min_match_bases_se_value)]
    if min_match_bases_pe is not None:
        min_match_bases_pe_value = _core.integer(
            min_match_bases_pe, "min_match_bases_pe", minimum=1
        )
        arguments += ["--MIN_MATCH_BASES_PE", str(min_match_bases_pe_value)]
    tool = "MarkIlluminaAdapters"
    reserved = _core.names(
        "ADAPTERS FIVE_PRIME_ADAPTER INPUT MAX_ERROR_RATE_PE MAX_ERROR_RATE_SE "
        "METRICS MIN_MATCH_BASES_PE MIN_MATCH_BASES_SE OUTPUT "
        "THREE_PRIME_ADAPTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_mark_illumina_adapters",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "metrics_file" in outs and outs["metrics_file"].is_file():
        parsed = _core.picard_metrics(outs["metrics_file"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_alignment_summary_metrics(
    input_path: str,
    output_path: str,
    *,
    reference: str | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    assume_sorted: bool | None = None,
    is_bisulfite_sequenced: bool | None = None,
    max_insert_size: int | None = None,
    expected_pair_orientations: Sequence[str] | None = None,
    adapter_sequences: Sequence[str] | None = None,
    histogram_file: str | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Produce a summary of alignment metrics (reads, PF reads, aligned bases, mismatch and indel
    rates, chimeras, adapters) from a SAM/BAM/CRAM file.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``is_bisulfite_sequenced`` (``--IS_BISULFITE_SEQUENCED``): Do not count C->T as a mismatch
      (bisulfite reads).
    - ``max_insert_size`` (``--MAX_INSERT_SIZE``): Larger insert sizes are treated as chimeras.
    - ``expected_pair_orientations`` (``--EXPECTED_PAIR_ORIENTATIONS``): Pair orientations that
      are not chimeric.
    - ``adapter_sequences`` (``--ADAPTER_SEQUENCE``): Adapter sequences to look for; None keeps
      the Illumina set of GATK.
    - ``histogram_file`` (``--HISTOGRAM_FILE``): New PDF read-length histogram (needs R).
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectAlignmentSummaryMetrics --help`` are accepted,
      a list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (one row per category and accumulation level, e.g.
    FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR) and ``histogram``.
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
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if is_bisulfite_sequenced is not None:
        arguments += [
            "--IS_BISULFITE_SEQUENCED",
            "true" if _core.flag(is_bisulfite_sequenced, "is_bisulfite_sequenced") else "false",
        ]
    if max_insert_size is not None:
        max_insert_size_value = _core.integer(max_insert_size, "max_insert_size", minimum=1)
        arguments += ["--MAX_INSERT_SIZE", str(max_insert_size_value)]
    if expected_pair_orientations is not None:
        expected_pair_orientations_values = _core.choices(
            expected_pair_orientations, ("FR", "RF", "TANDEM"), "expected_pair_orientations"
        )
        for item in expected_pair_orientations_values:
            arguments += ["--EXPECTED_PAIR_ORIENTATIONS", item]
    if adapter_sequences is not None:
        adapter_sequences_values = _core.texts(adapter_sequences, "adapter_sequences")
        for item in adapter_sequences_values:
            arguments += ["--ADAPTER_SEQUENCE", item]
    if histogram_file is not None:
        histogram_file_path = _core.new_file(histogram_file, "histogram_file")
        arguments += ["--HISTOGRAM_FILE", str(histogram_file_path)]
        expected.append(histogram_file_path)
        outs["histogram_file"] = histogram_file_path
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectAlignmentSummaryMetrics"
    reserved = _core.names(
        "ADAPTER_SEQUENCE ASSUME_SORTED EXPECTED_PAIR_ORIENTATIONS "
        "HISTOGRAM_FILE INPUT IS_BISULFITE_SEQUENCED MAX_INSERT_SIZE "
        "METRIC_ACCUMULATION_LEVEL OUTPUT REFERENCE_SEQUENCE STOP_AFTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_alignment_summary_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_base_distribution_by_cycle(
    input_path: str,
    output_path: str,
    chart_output: str,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    reference: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Chart the nucleotide distribution per sequencing cycle of a SAM/BAM/CRAM file and write the
    per-cycle percentages.

    The chart output needs R (Rscript) with the packages of the GATK conda environment; without
    it GATK stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New table of per-cycle base percentages.
    - ``chart_output`` (``--CHART_OUTPUT``): New PDF chart (needs R).
    - ``aligned_reads_only`` (``--ALIGNED_READS_ONLY``): Use only aligned reads.
    - ``pf_reads_only`` (``--PF_READS_ONLY``): Use only reads that pass the vendor filter.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectBaseDistributionByCycle --help`` are accepted,
      a list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (one row per read end and cycle with PCT_A, PCT_C,
    PCT_G, PCT_T, PCT_N).
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
    chart_output_path = _core.new_file(chart_output, "chart_output")
    arguments += ["--CHART_OUTPUT", str(chart_output_path)]
    expected.append(chart_output_path)
    outs["chart_output"] = chart_output_path
    if aligned_reads_only is not None:
        arguments += [
            "--ALIGNED_READS_ONLY",
            "true" if _core.flag(aligned_reads_only, "aligned_reads_only") else "false",
        ]
    if pf_reads_only is not None:
        arguments += [
            "--PF_READS_ONLY",
            "true" if _core.flag(pf_reads_only, "pf_reads_only") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectBaseDistributionByCycle"
    reserved = _core.names(
        "ALIGNED_READS_ONLY ASSUME_SORTED CHART_OUTPUT INPUT OUTPUT "
        "PF_READS_ONLY REFERENCE_SEQUENCE STOP_AFTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_base_distribution_by_cycle",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_mean_quality_by_cycle(
    input_path: str,
    output_path: str,
    chart_output: str,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    reference: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect the mean base quality for every sequencing cycle and chart it.

    The chart output needs R (Rscript) with the packages of the GATK conda environment; without
    it GATK stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New table of mean quality per cycle.
    - ``chart_output`` (``--CHART_OUTPUT``): New PDF chart (needs R).
    - ``aligned_reads_only`` (``--ALIGNED_READS_ONLY``): Use only aligned reads.
    - ``pf_reads_only`` (``--PF_READS_ONLY``): Use only reads that pass the vendor filter.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk MeanQualityByCycle --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` and ``histogram`` (cycle versus mean quality of all
    reads and of aligned reads).
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
    chart_output_path = _core.new_file(chart_output, "chart_output")
    arguments += ["--CHART_OUTPUT", str(chart_output_path)]
    expected.append(chart_output_path)
    outs["chart_output"] = chart_output_path
    if aligned_reads_only is not None:
        arguments += [
            "--ALIGNED_READS_ONLY",
            "true" if _core.flag(aligned_reads_only, "aligned_reads_only") else "false",
        ]
    if pf_reads_only is not None:
        arguments += [
            "--PF_READS_ONLY",
            "true" if _core.flag(pf_reads_only, "pf_reads_only") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "MeanQualityByCycle"
    reserved = _core.names(
        "ALIGNED_READS_ONLY ASSUME_SORTED CHART_OUTPUT INPUT OUTPUT "
        "PF_READS_ONLY REFERENCE_SEQUENCE STOP_AFTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_mean_quality_by_cycle",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_quality_score_distribution(
    input_path: str,
    output_path: str,
    chart_output: str,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    include_no_calls: bool | None = None,
    reference: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Chart the distribution of base quality scores (number of bases per quality value) of a
    SAM/BAM/CRAM file.

    The chart output needs R (Rscript) with the packages of the GATK conda environment; without
    it GATK stops with an error.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New table of quality score counts.
    - ``chart_output`` (``--CHART_OUTPUT``): New PDF chart (needs R).
    - ``aligned_reads_only`` (``--ALIGNED_READS_ONLY``): Use only aligned reads.
    - ``pf_reads_only`` (``--PF_READS_ONLY``): Use only reads that pass the vendor filter.
    - ``include_no_calls`` (``--INCLUDE_NO_CALLS``): Include bases that are N.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk QualityScoreDistribution --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` and ``histogram`` (quality score versus count of
    observations).
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
    chart_output_path = _core.new_file(chart_output, "chart_output")
    arguments += ["--CHART_OUTPUT", str(chart_output_path)]
    expected.append(chart_output_path)
    outs["chart_output"] = chart_output_path
    if aligned_reads_only is not None:
        arguments += [
            "--ALIGNED_READS_ONLY",
            "true" if _core.flag(aligned_reads_only, "aligned_reads_only") else "false",
        ]
    if pf_reads_only is not None:
        arguments += [
            "--PF_READS_ONLY",
            "true" if _core.flag(pf_reads_only, "pf_reads_only") else "false",
        ]
    if include_no_calls is not None:
        arguments += [
            "--INCLUDE_NO_CALLS",
            "true" if _core.flag(include_no_calls, "include_no_calls") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "QualityScoreDistribution"
    reserved = _core.names(
        "ALIGNED_READS_ONLY ASSUME_SORTED CHART_OUTPUT INCLUDE_NO_CALLS INPUT "
        "OUTPUT PF_READS_ONLY REFERENCE_SEQUENCE STOP_AFTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_quality_score_distribution",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_hiseqx_pf_fail_metrics(
    basecalls_dir: str,
    lane: int,
    output_path: str,
    *,
    n_cycles: int | None = None,
    prob_explicit_reads: float | None = None,
    num_processors: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify PF-failing reads of a HiSeq X BaseCalls directory into misaligned, empty,
    polyclonal and unexplained categories.

    Needs the BaseCalls directory of a HiSeq X run.

    Parameters (the GATK option is in parentheses):

    - ``basecalls_dir`` (``--BASECALLS_DIR``): Illumina HiSeq X BaseCalls directory.
    - ``lane`` (``--LANE``): Lane number.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``n_cycles`` (``--N_CYCLES``): Number of initial cycles examined.
    - ``prob_explicit_reads`` (``--PROB_EXPLICIT_READS``): Fraction of reads for which the
      classification is written explicitly.
    - ``num_processors`` (``--NUM_PROCESSORS``): Threads.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectHiSeqXPfFailMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics`` (counts per PF-fail category).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    basecalls_dir_path = _core.existing_directory(basecalls_dir, "basecalls_dir")
    arguments += ["--BASECALLS_DIR", str(basecalls_dir_path)]
    lane_value = _core.integer(lane, "lane", minimum=1)
    arguments += ["--LANE", str(lane_value)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if n_cycles is not None:
        n_cycles_value = _core.integer(n_cycles, "n_cycles", minimum=1)
        arguments += ["--N_CYCLES", str(n_cycles_value)]
    if prob_explicit_reads is not None:
        prob_explicit_reads_value = _core.number(
            prob_explicit_reads, "prob_explicit_reads", minimum=0.0, maximum=1.0
        )
        arguments += ["--PROB_EXPLICIT_READS", str(prob_explicit_reads_value)]
    if num_processors is not None:
        num_processors_value = _core.integer(num_processors, "num_processors")
        arguments += ["--NUM_PROCESSORS", str(num_processors_value)]
    tool = "CollectHiSeqXPfFailMetrics"
    reserved = _core.names("BASECALLS_DIR LANE NUM_PROCESSORS N_CYCLES OUTPUT PROB_EXPLICIT_READS")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_hiseqx_pf_fail_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_hs_metrics(
    input_path: str,
    output_path: str,
    bait_intervals: Sequence[str],
    target_intervals: Sequence[str],
    *,
    per_target_coverage: str | None = None,
    per_base_coverage: str | None = None,
    theoretical_sensitivity_output: str | None = None,
    reference: str | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    bait_set_name: str | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    clip_overlapping_reads: bool | None = None,
    coverage_cap: int | None = None,
    near_distance: int | None = None,
    include_indels: bool | None = None,
    sample_size: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect hybrid-selection (HS) metrics: on-bait and on-target bases, fold enrichment, target
    coverage and GC dropout.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``bait_intervals`` (``--BAIT_INTERVALS``): Interval-list files with the bait locations.
    - ``target_intervals`` (``--TARGET_INTERVALS``): Interval-list files with the target
      locations.
    - ``per_target_coverage`` (``--PER_TARGET_COVERAGE``): New per-target coverage table.
    - ``per_base_coverage`` (``--PER_BASE_COVERAGE``): New per-base coverage table.
    - ``theoretical_sensitivity_output`` (``--THEORETICAL_SENSITIVITY_OUTPUT``): New theoretical
      het-SNP sensitivity table.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``bait_set_name`` (``--BAIT_SET_NAME``): Name of the bait set; None derives it from the
      bait file.
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Reads below this mapping
      quality are ignored.
    - ``minimum_base_quality`` (``--MINIMUM_BASE_QUALITY``): Bases below this quality are
      ignored.
    - ``clip_overlapping_reads`` (``--CLIP_OVERLAPPING_READS``): Clip the overlap of read pairs.
    - ``coverage_cap`` (``--COVERAGE_CAP``): Coverage values above this are capped.
    - ``near_distance`` (``--NEAR_DISTANCE``): Distance around baits that counts as near-bait.
    - ``include_indels`` (``--INCLUDE_INDELS``): Count deletions as covered bases.
    - ``sample_size`` (``--SAMPLE_SIZE``): Sample size for the theoretical sensitivity.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectHsMetrics --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (one row per accumulation level with the HS metrics)
    and ``histogram``.
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
    bait_intervals_paths = _core.existing_files(bait_intervals, "bait_intervals")
    for item in bait_intervals_paths:
        arguments += ["--BAIT_INTERVALS", str(item)]
    target_intervals_paths = _core.existing_files(target_intervals, "target_intervals")
    for item in target_intervals_paths:
        arguments += ["--TARGET_INTERVALS", str(item)]
    if per_target_coverage is not None:
        per_target_coverage_path = _core.new_file(per_target_coverage, "per_target_coverage")
        arguments += ["--PER_TARGET_COVERAGE", str(per_target_coverage_path)]
        expected.append(per_target_coverage_path)
        outs["per_target_coverage"] = per_target_coverage_path
    if per_base_coverage is not None:
        per_base_coverage_path = _core.new_file(per_base_coverage, "per_base_coverage")
        arguments += ["--PER_BASE_COVERAGE", str(per_base_coverage_path)]
        expected.append(per_base_coverage_path)
        outs["per_base_coverage"] = per_base_coverage_path
    if theoretical_sensitivity_output is not None:
        theoretical_sensitivity_output_path = _core.new_file(
            theoretical_sensitivity_output, "theoretical_sensitivity_output"
        )
        arguments += ["--THEORETICAL_SENSITIVITY_OUTPUT", str(theoretical_sensitivity_output_path)]
        expected.append(theoretical_sensitivity_output_path)
        outs["theoretical_sensitivity_output"] = theoretical_sensitivity_output_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if bait_set_name is not None:
        bait_set_name_text = _core.text(bait_set_name, "bait_set_name")
        arguments += ["--BAIT_SET_NAME", bait_set_name_text]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if minimum_base_quality is not None:
        minimum_base_quality_value = _core.integer(
            minimum_base_quality, "minimum_base_quality", minimum=0
        )
        arguments += ["--MINIMUM_BASE_QUALITY", str(minimum_base_quality_value)]
    if clip_overlapping_reads is not None:
        arguments += [
            "--CLIP_OVERLAPPING_READS",
            "true" if _core.flag(clip_overlapping_reads, "clip_overlapping_reads") else "false",
        ]
    if coverage_cap is not None:
        coverage_cap_value = _core.integer(coverage_cap, "coverage_cap", minimum=1)
        arguments += ["--COVERAGE_CAP", str(coverage_cap_value)]
    if near_distance is not None:
        near_distance_value = _core.integer(near_distance, "near_distance", minimum=0)
        arguments += ["--NEAR_DISTANCE", str(near_distance_value)]
    if include_indels is not None:
        arguments += [
            "--INCLUDE_INDELS",
            "true" if _core.flag(include_indels, "include_indels") else "false",
        ]
    if sample_size is not None:
        sample_size_value = _core.integer(sample_size, "sample_size", minimum=1)
        arguments += ["--SAMPLE_SIZE", str(sample_size_value)]
    tool = "CollectHsMetrics"
    reserved = _core.names(
        "BAIT_INTERVALS BAIT_SET_NAME CLIP_OVERLAPPING_READS COVERAGE_CAP "
        "INCLUDE_INDELS INPUT METRIC_ACCUMULATION_LEVEL MINIMUM_BASE_QUALITY "
        "MINIMUM_MAPPING_QUALITY NEAR_DISTANCE OUTPUT PER_BASE_COVERAGE "
        "PER_TARGET_COVERAGE REFERENCE_SEQUENCE SAMPLE_SIZE TARGET_INTERVALS "
        "THEORETICAL_SENSITIVITY_OUTPUT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_hs_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_independent_replicate_metrics(
    input_path: str,
    vcf: str,
    output_path: str,
    *,
    matrix_output: str | None = None,
    sample: str | None = None,
    barcode_tag: str | None = None,
    barcode_bq: str | None = None,
    minimum_gq: int | None = None,
    minimum_bq: int | None = None,
    minimum_mq: int | None = None,
    minimum_barcode_bq: int | None = None,
    filter_unpaired_reads: bool | None = None,
    reference: str | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the fraction of duplicate-marked reads that are independent observations, using
    heterozygous sites of a VCF (experimental GATK tool).

    GATK marks this tool EXPERIMENTAL (use at your own risk).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Duplicate-marked BAM/CRAM with UMI tags.
    - ``vcf`` (``--VCF``): VCF of known heterozygous sites of the sample.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``matrix_output`` (``--MATRIX_OUTPUT``): New file for the duplicate-set classification
      matrix.
    - ``sample`` (``--SAMPLE``): Sample name in the VCF; None uses the only sample.
    - ``barcode_tag`` (``--BARCODE_TAG``): Read tag with the UMI.
    - ``barcode_bq`` (``--BARCODE_BQ``): Read tag with the UMI base qualities.
    - ``minimum_gq`` (``--MINIMUM_GQ``): Minimum genotype quality of a site.
    - ``minimum_bq`` (``--MINIMUM_BQ``): Minimum base quality.
    - ``minimum_mq`` (``--MINIMUM_MQ``): Minimum mapping quality.
    - ``minimum_barcode_bq`` (``--MINIMUM_BARCODE_BQ``): Minimum UMI base quality.
    - ``filter_unpaired_reads`` (``--FILTER_UNPAIRED_READS``): Ignore unpaired reads.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectIndependentReplicateMetrics --help`` are
      accepted, a list repeats the option, and options already set by a parameter above are
      refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    vcf_path = _core.existing_file(vcf, "vcf")
    arguments += ["--VCF", str(vcf_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if matrix_output is not None:
        matrix_output_path = _core.new_file(matrix_output, "matrix_output")
        arguments += ["--MATRIX_OUTPUT", str(matrix_output_path)]
        expected.append(matrix_output_path)
        outs["matrix_output"] = matrix_output_path
    if sample is not None:
        sample_text = _core.text(sample, "sample")
        arguments += ["--SAMPLE", sample_text]
    if barcode_tag is not None:
        barcode_tag_text = _core.text(barcode_tag, "barcode_tag")
        arguments += ["--BARCODE_TAG", barcode_tag_text]
    if barcode_bq is not None:
        barcode_bq_text = _core.text(barcode_bq, "barcode_bq")
        arguments += ["--BARCODE_BQ", barcode_bq_text]
    if minimum_gq is not None:
        minimum_gq_value = _core.integer(minimum_gq, "minimum_gq", minimum=0)
        arguments += ["--MINIMUM_GQ", str(minimum_gq_value)]
    if minimum_bq is not None:
        minimum_bq_value = _core.integer(minimum_bq, "minimum_bq", minimum=0)
        arguments += ["--MINIMUM_BQ", str(minimum_bq_value)]
    if minimum_mq is not None:
        minimum_mq_value = _core.integer(minimum_mq, "minimum_mq", minimum=0)
        arguments += ["--MINIMUM_MQ", str(minimum_mq_value)]
    if minimum_barcode_bq is not None:
        minimum_barcode_bq_value = _core.integer(
            minimum_barcode_bq, "minimum_barcode_bq", minimum=0
        )
        arguments += ["--MINIMUM_BARCODE_BQ", str(minimum_barcode_bq_value)]
    if filter_unpaired_reads is not None:
        arguments += [
            "--FILTER_UNPAIRED_READS",
            "true" if _core.flag(filter_unpaired_reads, "filter_unpaired_reads") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectIndependentReplicateMetrics"
    reserved = _core.names(
        "BARCODE_BQ BARCODE_TAG FILTER_UNPAIRED_READS INPUT MATRIX_OUTPUT "
        "MINIMUM_BARCODE_BQ MINIMUM_BQ MINIMUM_GQ MINIMUM_MQ OUTPUT "
        "REFERENCE_SEQUENCE SAMPLE STOP_AFTER VCF"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_independent_replicate_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_jumping_library_metrics(
    input_paths: Sequence[str],
    output_path: str,
    *,
    chimera_kb_min: int | None = None,
    minimum_mapping_quality: int | None = None,
    tail_limit: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect metrics of jumping (mate-pair) libraries: fragment, jump and non-jump pair counts
    and their orientations.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--INPUT``): Input SAM/BAM/CRAM files.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``chimera_kb_min`` (``--CHIMERA_KB_MIN``): Minimum insert size in kb to call a pair
      chimeric.
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Reads below this mapping
      quality are ignored.
    - ``tail_limit`` (``--TAIL_LIMIT``): Upper insert size kept in the histogram.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectJumpingLibraryMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` and ``histogram``.
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
    if chimera_kb_min is not None:
        chimera_kb_min_value = _core.integer(chimera_kb_min, "chimera_kb_min", minimum=1)
        arguments += ["--CHIMERA_KB_MIN", str(chimera_kb_min_value)]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if tail_limit is not None:
        tail_limit_value = _core.integer(tail_limit, "tail_limit", minimum=1)
        arguments += ["--TAIL_LIMIT", str(tail_limit_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CollectJumpingLibraryMetrics"
    reserved = _core.names(
        "CHIMERA_KB_MIN INPUT MINIMUM_MAPPING_QUALITY OUTPUT REFERENCE_SEQUENCE TAIL_LIMIT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_jumping_library_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_multiple_metrics(
    input_path: str,
    output_prefix: str,
    *,
    reference: str | None = None,
    programs: Sequence[str] | None = None,
    extra_arguments: Sequence[str] | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    db_snp: str | None = None,
    intervals_file: str | None = None,
    ref_flat: str | None = None,
    include_unpaired: bool | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run several Picard metric collectors (alignment summary, insert size, quality by cycle,
    quality distribution, base distribution, GC bias, sequencing artifacts, quality yield) in
    one pass over a file.

    All collectors except CollectQualityYieldMetrics and CollectSequencingArtifactMetrics write
    PDF charts and need R; pass programs=['null', 'CollectQualityYieldMetrics'] to run without
    R.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_prefix`` (``--OUTPUT``): File-name prefix of all outputs; each collector appends
      its own extension.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``programs`` (``--PROGRAM``): Collectors to run; None keeps the GATK default set
      (alignment summary, base distribution, insert size, mean quality, quality distribution).
      'null' first clears the default set.
    - ``extra_arguments`` (``--EXTRA_ARGUMENT``): Per-collector arguments written as
      'Collector::--OPTION value'.
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``db_snp`` (``--DB_SNP``): dbSNP VCF used by some collectors to exclude polymorphic sites.
    - ``intervals_file`` (``--INTERVALS``): Interval list that restricts the collectors that
      accept one.
    - ``ref_flat`` (``--REF_FLAT``): refFlat gene annotation (collectors that need it).
    - ``include_unpaired`` (``--INCLUDE_UNPAIRED``): Include unpaired reads in the insert size
      collector.
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectMultipleMetrics --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics_files`` (every file that starts with the prefix) and
    ``collected`` (the extensions written).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
    arguments += ["--OUTPUT", str(output_prefix_path)]
    outs["output_prefix"] = output_prefix_path
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if programs is not None:
        programs_values = _core.choices(
            programs,
            (
                "CollectAlignmentSummaryMetrics",
                "CollectInsertSizeMetrics",
                "QualityScoreDistribution",
                "MeanQualityByCycle",
                "CollectBaseDistributionByCycle",
                "CollectGcBiasMetrics",
                "CollectSequencingArtifactMetrics",
                "CollectQualityYieldMetrics",
                "null",
            ),
            "programs",
        )
        for item in programs_values:
            arguments += ["--PROGRAM", item]
    if extra_arguments is not None:
        extra_arguments_values = _core.texts(extra_arguments, "extra_arguments")
        for item in extra_arguments_values:
            arguments += ["--EXTRA_ARGUMENT", item]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if db_snp is not None:
        db_snp_path = _core.existing_file(db_snp, "db_snp")
        arguments += ["--DB_SNP", str(db_snp_path)]
    if intervals_file is not None:
        intervals_file_path = _core.existing_file(intervals_file, "intervals_file")
        arguments += ["--INTERVALS", str(intervals_file_path)]
    if ref_flat is not None:
        ref_flat_path = _core.existing_file(ref_flat, "ref_flat")
        arguments += ["--REF_FLAT", str(ref_flat_path)]
    if include_unpaired is not None:
        arguments += [
            "--INCLUDE_UNPAIRED",
            "true" if _core.flag(include_unpaired, "include_unpaired") else "false",
        ]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectMultipleMetrics"
    reserved = _core.names(
        "ASSUME_SORTED DB_SNP EXTRA_ARGUMENT INCLUDE_UNPAIRED INPUT INTERVALS "
        "METRIC_ACCUMULATION_LEVEL OUTPUT PROGRAM REFERENCE_SEQUENCE REF_FLAT "
        "STOP_AFTER"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_multiple_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    files = sorted(
        str(item) for item in output_prefix_path.parent.glob(output_prefix_path.name + ".*")
    )
    extras["metrics_files"] = files
    extras["collected"] = sorted(
        {Path(item).name[len(output_prefix_path.name) + 1 :].split(".")[0] for item in files}
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_oxog_metrics(
    input_path: str,
    output_path: str,
    reference: str,
    *,
    db_snp: str | None = None,
    intervals_file: str | None = None,
    contexts: Sequence[str] | None = None,
    context_size: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_quality_score: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_non_pf_reads: bool | None = None,
    use_oq: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect metrics that quantify oxidative (8-oxoguanine) artifacts: G>T error rates per
    sequence context.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``db_snp`` (``--DB_SNP``): dbSNP VCF whose sites are excluded.
    - ``intervals_file`` (``--INTERVALS``): Interval list to restrict the analysis.
    - ``contexts`` (``--CONTEXTS``): Sequence contexts to report, with N as wildcard.
    - ``context_size`` (``--CONTEXT_SIZE``): Number of bases on each side that define the
      context.
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Minimum mapping quality.
    - ``minimum_quality_score`` (``--MINIMUM_QUALITY_SCORE``): Minimum base quality.
    - ``minimum_insert_size`` (``--MINIMUM_INSERT_SIZE``): Minimum insert size.
    - ``maximum_insert_size`` (``--MAXIMUM_INSERT_SIZE``): Maximum insert size.
    - ``include_non_pf_reads`` (``--INCLUDE_NON_PF_READS``): Include reads that fail the vendor
      filter.
    - ``use_oq`` (``--USE_OQ``): Use the original base qualities (OQ tag) when present.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectOxoGMetrics --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics`` (one row per context, with the OxoG error rate).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if db_snp is not None:
        db_snp_path = _core.existing_file(db_snp, "db_snp")
        arguments += ["--DB_SNP", str(db_snp_path)]
    if intervals_file is not None:
        intervals_file_path = _core.existing_file(intervals_file, "intervals_file")
        arguments += ["--INTERVALS", str(intervals_file_path)]
    if contexts is not None:
        contexts_values = _core.texts(contexts, "contexts")
        for item in contexts_values:
            arguments += ["--CONTEXTS", item]
    if context_size is not None:
        context_size_value = _core.integer(context_size, "context_size", minimum=0)
        arguments += ["--CONTEXT_SIZE", str(context_size_value)]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if minimum_quality_score is not None:
        minimum_quality_score_value = _core.integer(
            minimum_quality_score, "minimum_quality_score", minimum=0
        )
        arguments += ["--MINIMUM_QUALITY_SCORE", str(minimum_quality_score_value)]
    if minimum_insert_size is not None:
        minimum_insert_size_value = _core.integer(
            minimum_insert_size, "minimum_insert_size", minimum=0
        )
        arguments += ["--MINIMUM_INSERT_SIZE", str(minimum_insert_size_value)]
    if maximum_insert_size is not None:
        maximum_insert_size_value = _core.integer(
            maximum_insert_size, "maximum_insert_size", minimum=0
        )
        arguments += ["--MAXIMUM_INSERT_SIZE", str(maximum_insert_size_value)]
    if include_non_pf_reads is not None:
        arguments += [
            "--INCLUDE_NON_PF_READS",
            "true" if _core.flag(include_non_pf_reads, "include_non_pf_reads") else "false",
        ]
    if use_oq is not None:
        arguments += ["--USE_OQ", "true" if _core.flag(use_oq, "use_oq") else "false"]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectOxoGMetrics"
    reserved = _core.names(
        "CONTEXTS CONTEXT_SIZE DB_SNP INCLUDE_NON_PF_READS INPUT INTERVALS "
        "MAXIMUM_INSERT_SIZE MINIMUM_INSERT_SIZE MINIMUM_MAPPING_QUALITY "
        "MINIMUM_QUALITY_SCORE OUTPUT REFERENCE_SEQUENCE STOP_AFTER USE_OQ"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_oxog_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_quality_yield_metrics(
    input_path: str,
    output_path: str,
    variant: str,
    *,
    use_original_qualities: bool | None = None,
    include_bq_histogram: bool | None = None,
    alternate_quality_attribute: str | None = None,
    include_secondary_alignments: bool | None = None,
    include_supplemental_alignments: bool | None = None,
    reference: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect yield metrics of reads that pass quality thresholds and vendor filters:
    variant='standard' (CollectQualityYieldMetrics), 'flow' (flow-based read files,
    CollectQualityYieldMetricsFlow) or 'snvq' (SNV quality, CollectQualityYieldMetricsSNVQ).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``variant`` (``--None``): Which collector to run.
    - ``use_original_qualities`` (``--USE_ORIGINAL_QUALITIES``): Use original qualities (OQ tag)
      when present ('standard' only).
    - ``include_bq_histogram`` (``--INCLUDE_BQ_HISTOGRAM``): Also report the base-quality
      histogram ('flow' and 'snvq').
    - ``alternate_quality_attribute`` (``--ALTERNATE_QUALITY_ATTRIBUTE``): Read attribute that
      holds the SNV quality ('snvq' only).
    - ``include_secondary_alignments`` (``--INCLUDE_SECONDARY_ALIGNMENTS``): Count secondary
      alignments.
    - ``include_supplemental_alignments`` (``--INCLUDE_SUPPLEMENTAL_ALIGNMENTS``): Count
      supplementary alignments.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk <tool> --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``variant``, ``metrics`` (total reads, PF reads, bases, Q20 and
    Q30 bases) and ``histogram``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    _core.choice(variant, ("standard", "flow", "snvq"), "variant")
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if use_original_qualities is not None and variant not in ("standard",):
        raise ConfigurationError("use_original_qualities only applies to variant in ['standard'].")
    if use_original_qualities is not None:
        arguments += [
            "--USE_ORIGINAL_QUALITIES",
            "true" if _core.flag(use_original_qualities, "use_original_qualities") else "false",
        ]
    if include_bq_histogram is not None and variant not in ("flow", "snvq"):
        raise ConfigurationError(
            "include_bq_histogram only applies to variant in ['flow', 'snvq']."
        )
    if include_bq_histogram is not None:
        arguments += [
            "--INCLUDE_BQ_HISTOGRAM",
            "true" if _core.flag(include_bq_histogram, "include_bq_histogram") else "false",
        ]
    if alternate_quality_attribute is not None and variant not in ("snvq",):
        raise ConfigurationError("alternate_quality_attribute only applies to variant in ['snvq'].")
    if alternate_quality_attribute is not None:
        alternate_quality_attribute_text = _core.text(
            alternate_quality_attribute, "alternate_quality_attribute"
        )
        arguments += ["--ALTERNATE_QUALITY_ATTRIBUTE", alternate_quality_attribute_text]
    if include_secondary_alignments is not None:
        arguments += [
            "--INCLUDE_SECONDARY_ALIGNMENTS",
            "true"
            if _core.flag(include_secondary_alignments, "include_secondary_alignments")
            else "false",
        ]
    if include_supplemental_alignments is not None:
        arguments += [
            "--INCLUDE_SUPPLEMENTAL_ALIGNMENTS",
            "true"
            if _core.flag(include_supplemental_alignments, "include_supplemental_alignments")
            else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = {
        "standard": "CollectQualityYieldMetrics",
        "flow": "CollectQualityYieldMetricsFlow",
        "snvq": "CollectQualityYieldMetricsSNVQ",
    }[variant]
    reserved = _core.names(
        "ALTERNATE_QUALITY_ATTRIBUTE ASSUME_SORTED INCLUDE_BQ_HISTOGRAM "
        "INCLUDE_SECONDARY_ALIGNMENTS INCLUDE_SUPPLEMENTAL_ALIGNMENTS INPUT "
        "OUTPUT REFERENCE_SEQUENCE STOP_AFTER USE_ORIGINAL_QUALITIES"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_quality_yield_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_rnaseq_metrics(
    input_path: str,
    output_path: str,
    ref_flat: str,
    strand_specificity: str,
    *,
    ribosomal_intervals: str | None = None,
    chart_output: str | None = None,
    end_bias_bases: int | None = None,
    minimum_length: int | None = None,
    rrna_fragment_percentage: float | None = None,
    ignore_sequences: Sequence[str] | None = None,
    reference: str | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Produce RNA-seq alignment metrics: bases in coding, UTR, intronic, intergenic and ribosomal
    regions, strand specificity and 5'/3' bias.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``ref_flat`` (``--REF_FLAT``): Gene annotations in refFlat format.
    - ``strand_specificity`` (``--STRAND_SPECIFICITY``): Library strandedness.
    - ``ribosomal_intervals`` (``--RIBOSOMAL_INTERVALS``): Interval list of ribosomal sequence.
    - ``chart_output`` (``--CHART_OUTPUT``): New PDF chart of coverage versus transcript
      position (needs R).
    - ``end_bias_bases`` (``--END_BIAS_BASES``): Bases from each transcript end used for the
      5'/3' bias.
    - ``minimum_length`` (``--MINIMUM_LENGTH``): Minimum transcript length for the coverage
      profile.
    - ``rrna_fragment_percentage`` (``--RRNA_FRAGMENT_PERCENTAGE``): Fragments overlapping this
      fraction of ribosomal sequence count as ribosomal.
    - ``ignore_sequences`` (``--IGNORE_SEQUENCE``): Sequences (for example chrM) to ignore.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectRnaSeqMetrics --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` and ``histogram`` (normalised coverage by transcript
    position).
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
    ref_flat_path = _core.existing_file(ref_flat, "ref_flat")
    arguments += ["--REF_FLAT", str(ref_flat_path)]
    strand_specificity_value = _core.choice(
        strand_specificity,
        ("NONE", "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"),
        "strand_specificity",
    )
    arguments += ["--STRAND_SPECIFICITY", strand_specificity_value]
    if ribosomal_intervals is not None:
        ribosomal_intervals_path = _core.existing_file(ribosomal_intervals, "ribosomal_intervals")
        arguments += ["--RIBOSOMAL_INTERVALS", str(ribosomal_intervals_path)]
    if chart_output is not None:
        chart_output_path = _core.new_file(chart_output, "chart_output")
        arguments += ["--CHART_OUTPUT", str(chart_output_path)]
        expected.append(chart_output_path)
        outs["chart_output"] = chart_output_path
    if end_bias_bases is not None:
        end_bias_bases_value = _core.integer(end_bias_bases, "end_bias_bases", minimum=1)
        arguments += ["--END_BIAS_BASES", str(end_bias_bases_value)]
    if minimum_length is not None:
        minimum_length_value = _core.integer(minimum_length, "minimum_length", minimum=1)
        arguments += ["--MINIMUM_LENGTH", str(minimum_length_value)]
    if rrna_fragment_percentage is not None:
        rrna_fragment_percentage_value = _core.number(
            rrna_fragment_percentage, "rrna_fragment_percentage", minimum=0.0, maximum=1.0
        )
        arguments += ["--RRNA_FRAGMENT_PERCENTAGE", str(rrna_fragment_percentage_value)]
    if ignore_sequences is not None:
        ignore_sequences_values = _core.texts(ignore_sequences, "ignore_sequences")
        for item in ignore_sequences_values:
            arguments += ["--IGNORE_SEQUENCE", item]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectRnaSeqMetrics"
    reserved = _core.names(
        "ASSUME_SORTED CHART_OUTPUT END_BIAS_BASES IGNORE_SEQUENCE INPUT "
        "METRIC_ACCUMULATION_LEVEL MINIMUM_LENGTH OUTPUT REFERENCE_SEQUENCE "
        "REF_FLAT RIBOSOMAL_INTERVALS RRNA_FRAGMENT_PERCENTAGE STOP_AFTER "
        "STRAND_SPECIFICITY"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_rnaseq_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_rrbs_metrics(
    input_path: str,
    metrics_file_prefix: str,
    reference: str,
    *,
    c_quality_threshold: int | None = None,
    next_base_quality_threshold: int | None = None,
    minimum_read_length: int | None = None,
    max_mismatch_rate: float | None = None,
    sequence_names: Sequence[str] | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    assume_sorted: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect metrics from reduced-representation bisulfite sequencing (RRBS): CpG conversion,
    coverage and read-start counts.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input BAM file of bisulfite reads.
    - ``metrics_file_prefix`` (``--METRICS_FILE_PREFIX``): Prefix of the output files
      (.rrbs_summary_metrics, .rrbs_detail_metrics and .rrbs_qc.pdf).
    - ``reference`` (``--REFERENCE``): Reference FASTA with its .fai and .dict files.
    - ``c_quality_threshold`` (``--C_QUALITY_THRESHOLD``): Minimum quality of a C base to count
      it.
    - ``next_base_quality_threshold`` (``--NEXT_BASE_QUALITY_THRESHOLD``): Minimum quality of
      the base after the C.
    - ``minimum_read_length`` (``--MINIMUM_READ_LENGTH``): Reads shorter than this are ignored.
    - ``max_mismatch_rate`` (``--MAX_MISMATCH_RATE``): Reads with a higher mismatch rate are
      ignored.
    - ``sequence_names`` (``--SEQUENCE_NAMES``): Sequences to analyse; None uses all.
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectRrbsMetrics --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics_files`` and ``metrics`` (rows of the summary metrics).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    metrics_file_prefix_path = _core.new_prefix(metrics_file_prefix, "metrics_file_prefix")
    arguments += ["--METRICS_FILE_PREFIX", str(metrics_file_prefix_path)]
    outs["metrics_file_prefix"] = metrics_file_prefix_path
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE", str(reference_path)]
    if c_quality_threshold is not None:
        c_quality_threshold_value = _core.integer(
            c_quality_threshold, "c_quality_threshold", minimum=0
        )
        arguments += ["--C_QUALITY_THRESHOLD", str(c_quality_threshold_value)]
    if next_base_quality_threshold is not None:
        next_base_quality_threshold_value = _core.integer(
            next_base_quality_threshold, "next_base_quality_threshold", minimum=0
        )
        arguments += ["--NEXT_BASE_QUALITY_THRESHOLD", str(next_base_quality_threshold_value)]
    if minimum_read_length is not None:
        minimum_read_length_value = _core.integer(
            minimum_read_length, "minimum_read_length", minimum=1
        )
        arguments += ["--MINIMUM_READ_LENGTH", str(minimum_read_length_value)]
    if max_mismatch_rate is not None:
        max_mismatch_rate_value = _core.number(
            max_mismatch_rate, "max_mismatch_rate", minimum=0.0, maximum=1.0
        )
        arguments += ["--MAX_MISMATCH_RATE", str(max_mismatch_rate_value)]
    if sequence_names is not None:
        sequence_names_values = _core.texts(sequence_names, "sequence_names")
        for item in sequence_names_values:
            arguments += ["--SEQUENCE_NAMES", item]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    tool = "CollectRrbsMetrics"
    reserved = _core.names(
        "ASSUME_SORTED C_QUALITY_THRESHOLD INPUT MAX_MISMATCH_RATE "
        "METRICS_FILE_PREFIX METRIC_ACCUMULATION_LEVEL MINIMUM_READ_LENGTH "
        "NEXT_BASE_QUALITY_THRESHOLD REFERENCE SEQUENCE_NAMES"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_rrbs_metrics",
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
        for item in metrics_file_prefix_path.parent.glob(metrics_file_prefix_path.name + ".rrbs_*")
    )
    extras["metrics_files"] = files
    summary = [Path(item) for item in files if item.endswith(".rrbs_summary_metrics")]
    if summary:
        parsed = _core.picard_metrics(summary[0])
        extras["metrics"] = parsed["metrics"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_sequencing_artifact_metrics(
    input_path: str,
    output_prefix: str,
    reference: str,
    *,
    db_snp: str | None = None,
    intervals_file: str | None = None,
    contexts_to_print: Sequence[str] | None = None,
    context_size: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_quality_score: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_duplicates: bool | None = None,
    include_non_pf_reads: bool | None = None,
    include_unpaired: bool | None = None,
    tandem_reads: bool | None = None,
    use_oq: bool | None = None,
    file_extension: str | None = None,
    assume_sorted: bool | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Quantify single-base sequencing artifacts (pre-adapter, bait-bias, error summaries by
    substitution and context).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_prefix`` (``--OUTPUT``): Prefix of the output files (.pre_adapter_detail_metrics,
      .bait_bias_detail_metrics, .error_summary_metrics ...).
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``db_snp`` (``--DB_SNP``): dbSNP VCF whose sites are excluded.
    - ``intervals_file`` (``--INTERVALS``): Interval list to restrict the analysis.
    - ``contexts_to_print`` (``--CONTEXTS_TO_PRINT``): Contexts written to the detail files.
    - ``context_size`` (``--CONTEXT_SIZE``): Number of bases on each side that define the
      context.
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Minimum mapping quality.
    - ``minimum_quality_score`` (``--MINIMUM_QUALITY_SCORE``): Minimum base quality.
    - ``minimum_insert_size`` (``--MINIMUM_INSERT_SIZE``): Minimum insert size.
    - ``maximum_insert_size`` (``--MAXIMUM_INSERT_SIZE``): Maximum insert size.
    - ``include_duplicates`` (``--INCLUDE_DUPLICATES``): Include duplicate-marked reads.
    - ``include_non_pf_reads`` (``--INCLUDE_NON_PF_READS``): Include reads that fail the vendor
      filter.
    - ``include_unpaired`` (``--INCLUDE_UNPAIRED``): Include unpaired reads.
    - ``tandem_reads`` (``--TANDEM_READS``): Pairs are in tandem orientation (mate-pair style).
    - ``use_oq`` (``--USE_OQ``): Use the original base qualities (OQ tag) when present.
    - ``file_extension`` (``--FILE_EXTENSION``): Extension appended to every output file name.
    - ``assume_sorted`` (``--ASSUME_SORTED``): Ignore the sort order stated in the header.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectSequencingArtifactMetrics --help`` are
      accepted, a list repeats the option, and options already set by a parameter above are
      refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics_files`` and ``metrics`` (rows of the
    .error_summary_metrics file).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
    arguments += ["--OUTPUT", str(output_prefix_path)]
    outs["output_prefix"] = output_prefix_path
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if db_snp is not None:
        db_snp_path = _core.existing_file(db_snp, "db_snp")
        arguments += ["--DB_SNP", str(db_snp_path)]
    if intervals_file is not None:
        intervals_file_path = _core.existing_file(intervals_file, "intervals_file")
        arguments += ["--INTERVALS", str(intervals_file_path)]
    if contexts_to_print is not None:
        contexts_to_print_values = _core.texts(contexts_to_print, "contexts_to_print")
        for item in contexts_to_print_values:
            arguments += ["--CONTEXTS_TO_PRINT", item]
    if context_size is not None:
        context_size_value = _core.integer(context_size, "context_size", minimum=0)
        arguments += ["--CONTEXT_SIZE", str(context_size_value)]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if minimum_quality_score is not None:
        minimum_quality_score_value = _core.integer(
            minimum_quality_score, "minimum_quality_score", minimum=0
        )
        arguments += ["--MINIMUM_QUALITY_SCORE", str(minimum_quality_score_value)]
    if minimum_insert_size is not None:
        minimum_insert_size_value = _core.integer(
            minimum_insert_size, "minimum_insert_size", minimum=0
        )
        arguments += ["--MINIMUM_INSERT_SIZE", str(minimum_insert_size_value)]
    if maximum_insert_size is not None:
        maximum_insert_size_value = _core.integer(
            maximum_insert_size, "maximum_insert_size", minimum=0
        )
        arguments += ["--MAXIMUM_INSERT_SIZE", str(maximum_insert_size_value)]
    if include_duplicates is not None:
        arguments += [
            "--INCLUDE_DUPLICATES",
            "true" if _core.flag(include_duplicates, "include_duplicates") else "false",
        ]
    if include_non_pf_reads is not None:
        arguments += [
            "--INCLUDE_NON_PF_READS",
            "true" if _core.flag(include_non_pf_reads, "include_non_pf_reads") else "false",
        ]
    if include_unpaired is not None:
        arguments += [
            "--INCLUDE_UNPAIRED",
            "true" if _core.flag(include_unpaired, "include_unpaired") else "false",
        ]
    if tandem_reads is not None:
        arguments += [
            "--TANDEM_READS",
            "true" if _core.flag(tandem_reads, "tandem_reads") else "false",
        ]
    if use_oq is not None:
        arguments += ["--USE_OQ", "true" if _core.flag(use_oq, "use_oq") else "false"]
    if file_extension is not None:
        file_extension_text = _core.text(file_extension, "file_extension")
        arguments += ["--FILE_EXTENSION", file_extension_text]
    if assume_sorted is not None:
        arguments += [
            "--ASSUME_SORTED",
            "true" if _core.flag(assume_sorted, "assume_sorted") else "false",
        ]
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectSequencingArtifactMetrics"
    reserved = _core.names(
        "ASSUME_SORTED CONTEXTS_TO_PRINT CONTEXT_SIZE DB_SNP FILE_EXTENSION "
        "INCLUDE_DUPLICATES INCLUDE_NON_PF_READS INCLUDE_UNPAIRED INPUT "
        "INTERVALS MAXIMUM_INSERT_SIZE MINIMUM_INSERT_SIZE "
        "MINIMUM_MAPPING_QUALITY MINIMUM_QUALITY_SCORE OUTPUT "
        "REFERENCE_SEQUENCE STOP_AFTER TANDEM_READS USE_OQ"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_sequencing_artifact_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    files = sorted(
        str(item) for item in output_prefix_path.parent.glob(output_prefix_path.name + ".*")
    )
    extras["metrics_files"] = files
    summary = [Path(item) for item in files if item.endswith(".error_summary_metrics")]
    if summary:
        extras["metrics"] = _core.picard_metrics(summary[0])["metrics"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_targeted_pcr_metrics(
    input_path: str,
    output_path: str,
    amplicon_intervals: str,
    target_intervals: Sequence[str],
    *,
    per_target_coverage: str | None = None,
    per_base_coverage: str | None = None,
    theoretical_sensitivity_output: str | None = None,
    custom_amplicon_set_name: str | None = None,
    reference: str | None = None,
    metric_accumulation_level: Sequence[str] | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    clip_overlapping_reads: bool | None = None,
    coverage_cap: int | None = None,
    near_distance: int | None = None,
    include_indels: bool | None = None,
    sample_size: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate PCR-related metrics from targeted amplicon sequencing: on-amplicon and on-target
    bases, coverage and sensitivity.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``amplicon_intervals`` (``--AMPLICON_INTERVALS``): Interval list of the amplicons.
    - ``target_intervals`` (``--TARGET_INTERVALS``): Interval lists of the targets.
    - ``per_target_coverage`` (``--PER_TARGET_COVERAGE``): New per-target coverage table.
    - ``per_base_coverage`` (``--PER_BASE_COVERAGE``): New per-base coverage table.
    - ``theoretical_sensitivity_output`` (``--THEORETICAL_SENSITIVITY_OUTPUT``): New theoretical
      het-SNP sensitivity table.
    - ``custom_amplicon_set_name`` (``--CUSTOM_AMPLICON_SET_NAME``): Name of the amplicon set.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``metric_accumulation_level`` (``--METRIC_ACCUMULATION_LEVEL``): Levels at which metrics
      are accumulated: ALL_READS, SAMPLE, LIBRARY, READ_GROUP (GATK default ALL_READS).
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Minimum mapping quality.
    - ``minimum_base_quality`` (``--MINIMUM_BASE_QUALITY``): Minimum base quality.
    - ``clip_overlapping_reads`` (``--CLIP_OVERLAPPING_READS``): Clip the overlap of read pairs.
    - ``coverage_cap`` (``--COVERAGE_CAP``): Coverage values above this are capped.
    - ``near_distance`` (``--NEAR_DISTANCE``): Distance around amplicons that counts as near-
      amplicon.
    - ``include_indels`` (``--INCLUDE_INDELS``): Count deletions as covered bases.
    - ``sample_size`` (``--SAMPLE_SIZE``): Sample size for the theoretical sensitivity.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectTargetedPcrMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics``.
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
    amplicon_intervals_path = _core.existing_file(amplicon_intervals, "amplicon_intervals")
    arguments += ["--AMPLICON_INTERVALS", str(amplicon_intervals_path)]
    target_intervals_paths = _core.existing_files(target_intervals, "target_intervals")
    for item in target_intervals_paths:
        arguments += ["--TARGET_INTERVALS", str(item)]
    if per_target_coverage is not None:
        per_target_coverage_path = _core.new_file(per_target_coverage, "per_target_coverage")
        arguments += ["--PER_TARGET_COVERAGE", str(per_target_coverage_path)]
        expected.append(per_target_coverage_path)
        outs["per_target_coverage"] = per_target_coverage_path
    if per_base_coverage is not None:
        per_base_coverage_path = _core.new_file(per_base_coverage, "per_base_coverage")
        arguments += ["--PER_BASE_COVERAGE", str(per_base_coverage_path)]
        expected.append(per_base_coverage_path)
        outs["per_base_coverage"] = per_base_coverage_path
    if theoretical_sensitivity_output is not None:
        theoretical_sensitivity_output_path = _core.new_file(
            theoretical_sensitivity_output, "theoretical_sensitivity_output"
        )
        arguments += ["--THEORETICAL_SENSITIVITY_OUTPUT", str(theoretical_sensitivity_output_path)]
        expected.append(theoretical_sensitivity_output_path)
        outs["theoretical_sensitivity_output"] = theoretical_sensitivity_output_path
    if custom_amplicon_set_name is not None:
        custom_amplicon_set_name_text = _core.text(
            custom_amplicon_set_name, "custom_amplicon_set_name"
        )
        arguments += ["--CUSTOM_AMPLICON_SET_NAME", custom_amplicon_set_name_text]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if metric_accumulation_level is not None:
        metric_accumulation_level_values = _core.choices(
            metric_accumulation_level,
            ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
            "metric_accumulation_level",
        )
        for item in metric_accumulation_level_values:
            arguments += ["--METRIC_ACCUMULATION_LEVEL", item]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if minimum_base_quality is not None:
        minimum_base_quality_value = _core.integer(
            minimum_base_quality, "minimum_base_quality", minimum=0
        )
        arguments += ["--MINIMUM_BASE_QUALITY", str(minimum_base_quality_value)]
    if clip_overlapping_reads is not None:
        arguments += [
            "--CLIP_OVERLAPPING_READS",
            "true" if _core.flag(clip_overlapping_reads, "clip_overlapping_reads") else "false",
        ]
    if coverage_cap is not None:
        coverage_cap_value = _core.integer(coverage_cap, "coverage_cap", minimum=1)
        arguments += ["--COVERAGE_CAP", str(coverage_cap_value)]
    if near_distance is not None:
        near_distance_value = _core.integer(near_distance, "near_distance", minimum=0)
        arguments += ["--NEAR_DISTANCE", str(near_distance_value)]
    if include_indels is not None:
        arguments += [
            "--INCLUDE_INDELS",
            "true" if _core.flag(include_indels, "include_indels") else "false",
        ]
    if sample_size is not None:
        sample_size_value = _core.integer(sample_size, "sample_size", minimum=1)
        arguments += ["--SAMPLE_SIZE", str(sample_size_value)]
    tool = "CollectTargetedPcrMetrics"
    reserved = _core.names(
        "AMPLICON_INTERVALS CLIP_OVERLAPPING_READS COVERAGE_CAP "
        "CUSTOM_AMPLICON_SET_NAME INCLUDE_INDELS INPUT "
        "METRIC_ACCUMULATION_LEVEL MINIMUM_BASE_QUALITY MINIMUM_MAPPING_QUALITY "
        "NEAR_DISTANCE OUTPUT PER_BASE_COVERAGE PER_TARGET_COVERAGE "
        "REFERENCE_SEQUENCE SAMPLE_SIZE TARGET_INTERVALS "
        "THEORETICAL_SENSITIVITY_OUTPUT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_targeted_pcr_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_umi_prevalence_metrics(
    input_path: str,
    output_path: str,
    *,
    barcode_tag: str | None = None,
    barcode_bq: str | None = None,
    minimum_barcode_bq: int | None = None,
    minimum_mq: int | None = None,
    filter_unpaired_reads: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tally the counts of UMIs in duplicate sets of a duplicate-marked BAM (how many duplicate
    sets have 1, 2, 3 ... distinct UMIs).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Duplicate-marked BAM/CRAM with UMI tags.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``barcode_tag`` (``--BARCODE_TAG``): Read tag with the UMI.
    - ``barcode_bq`` (``--BARCODE_BQ``): Read tag with the UMI base qualities.
    - ``minimum_barcode_bq`` (``--MINIMUM_BARCODE_BQ``): Minimum UMI base quality.
    - ``minimum_mq`` (``--MINIMUM_MQ``): Minimum mapping quality.
    - ``filter_unpaired_reads`` (``--FILTER_UNPAIRED_READS``): Ignore unpaired reads.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectUmiPrevalenceMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics``.
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
    if barcode_tag is not None:
        barcode_tag_text = _core.text(barcode_tag, "barcode_tag")
        arguments += ["--BARCODE_TAG", barcode_tag_text]
    if barcode_bq is not None:
        barcode_bq_text = _core.text(barcode_bq, "barcode_bq")
        arguments += ["--BARCODE_BQ", barcode_bq_text]
    if minimum_barcode_bq is not None:
        minimum_barcode_bq_value = _core.integer(
            minimum_barcode_bq, "minimum_barcode_bq", minimum=0
        )
        arguments += ["--MINIMUM_BARCODE_BQ", str(minimum_barcode_bq_value)]
    if minimum_mq is not None:
        minimum_mq_value = _core.integer(minimum_mq, "minimum_mq", minimum=0)
        arguments += ["--MINIMUM_MQ", str(minimum_mq_value)]
    if filter_unpaired_reads is not None:
        arguments += [
            "--FILTER_UNPAIRED_READS",
            "true" if _core.flag(filter_unpaired_reads, "filter_unpaired_reads") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CollectUmiPrevalenceMetrics"
    reserved = _core.names(
        "BARCODE_BQ BARCODE_TAG FILTER_UNPAIRED_READS INPUT MINIMUM_BARCODE_BQ "
        "MINIMUM_MQ OUTPUT REFERENCE_SEQUENCE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_umi_prevalence_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_variant_calling_metrics(
    input_path: str,
    dbsnp: str,
    output_prefix: str,
    *,
    target_intervals: str | None = None,
    sequence_dictionary: str | None = None,
    gvcf_input: bool | None = None,
    thread_count: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect per-sample and aggregate variant-calling metrics (SNP/indel counts, Ti/Tv, het/hom
    ratio, dbSNP concordance) from a VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input VCF (or GVCF with gvcf_input=True).
    - ``dbsnp`` (``--DBSNP``): dbSNP VCF used for the known/novel split.
    - ``output_prefix`` (``--OUTPUT``): Prefix of the outputs (.variant_calling_detail_metrics
      and .variant_calling_summary_metrics).
    - ``target_intervals`` (``--TARGET_INTERVALS``): Interval list to restrict the analysis.
    - ``sequence_dictionary`` (``--SEQUENCE_DICTIONARY``): Sequence dictionary when the VCF has
      no contig lines.
    - ``gvcf_input`` (``--GVCF_INPUT``): The input is a GVCF.
    - ``thread_count`` (``--THREAD_COUNT``): Threads.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectVariantCallingMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics_files``, ``summary_metrics`` and ``detail_metrics``
    (rows of the two tables).
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
    if target_intervals is not None:
        target_intervals_path = _core.existing_file(target_intervals, "target_intervals")
        arguments += ["--TARGET_INTERVALS", str(target_intervals_path)]
    if sequence_dictionary is not None:
        sequence_dictionary_path = _core.existing_file(sequence_dictionary, "sequence_dictionary")
        arguments += ["--SEQUENCE_DICTIONARY", str(sequence_dictionary_path)]
    if gvcf_input is not None:
        arguments += ["--GVCF_INPUT", "true" if _core.flag(gvcf_input, "gvcf_input") else "false"]
    if thread_count is not None:
        thread_count_value = _core.integer(thread_count, "thread_count", minimum=1)
        arguments += ["--THREAD_COUNT", str(thread_count_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CollectVariantCallingMetrics"
    reserved = _core.names(
        "DBSNP GVCF_INPUT INPUT OUTPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY "
        "TARGET_INTERVALS THREAD_COUNT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_variant_calling_metrics",
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
        for item in output_prefix_path.parent.glob(output_prefix_path.name + ".variant_calling_*")
    )
    extras["metrics_files"] = files
    for suffix, key in (
        (".variant_calling_summary_metrics", "summary_metrics"),
        (".variant_calling_detail_metrics", "detail_metrics"),
    ):
        chosen = [Path(item) for item in files if item.endswith(suffix)]
        if chosen:
            extras[key] = _core.picard_metrics(chosen[0])["metrics"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_collect_wgs_metrics(
    input_path: str,
    output_path: str,
    reference: str,
    *,
    intervals_file: str | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    coverage_cap: int | None = None,
    read_length: int | None = None,
    count_unpaired: bool | None = None,
    include_bq_histogram: bool | None = None,
    use_fast_algorithm: bool | None = None,
    sample_size: int | None = None,
    locus_accumulation_cap: int | None = None,
    theoretical_sensitivity_output: str | None = None,
    stop_after: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Collect coverage and performance metrics of a whole-genome sequencing experiment: mean
    coverage, fraction of bases at depth thresholds, excluded-base fractions.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``intervals_file`` (``--INTERVALS``): Interval list to restrict the analysis.
    - ``minimum_mapping_quality`` (``--MINIMUM_MAPPING_QUALITY``): Reads below this mapping
      quality are excluded.
    - ``minimum_base_quality`` (``--MINIMUM_BASE_QUALITY``): Bases below this quality are
      excluded.
    - ``coverage_cap`` (``--COVERAGE_CAP``): Coverage values above this are capped.
    - ``read_length`` (``--READ_LENGTH``): Read length for the theoretical sensitivity.
    - ``count_unpaired`` (``--COUNT_UNPAIRED``): Count unpaired reads and reads whose mate is
      unmapped.
    - ``include_bq_histogram`` (``--INCLUDE_BQ_HISTOGRAM``): Also report the base-quality
      histogram.
    - ``use_fast_algorithm`` (``--USE_FAST_ALGORITHM``): Use the faster WGS algorithm.
    - ``sample_size`` (``--SAMPLE_SIZE``): Sample size for the theoretical sensitivity.
    - ``locus_accumulation_cap`` (``--LOCUS_ACCUMULATION_CAP``): Maximum reads accumulated per
      locus.
    - ``theoretical_sensitivity_output`` (``--THEORETICAL_SENSITIVITY_OUTPUT``): New theoretical
      het-SNP sensitivity table.
    - ``stop_after`` (``--STOP_AFTER``): Stop after this many records; 0 reads everything (GATK
      default).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CollectWgsMetrics --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (GENOME_TERRITORY, MEAN_COVERAGE, PCT_* ...) and
    ``histogram`` (depth versus bases).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if intervals_file is not None:
        intervals_file_path = _core.existing_file(intervals_file, "intervals_file")
        arguments += ["--INTERVALS", str(intervals_file_path)]
    if minimum_mapping_quality is not None:
        minimum_mapping_quality_value = _core.integer(
            minimum_mapping_quality, "minimum_mapping_quality", minimum=0
        )
        arguments += ["--MINIMUM_MAPPING_QUALITY", str(minimum_mapping_quality_value)]
    if minimum_base_quality is not None:
        minimum_base_quality_value = _core.integer(
            minimum_base_quality, "minimum_base_quality", minimum=0
        )
        arguments += ["--MINIMUM_BASE_QUALITY", str(minimum_base_quality_value)]
    if coverage_cap is not None:
        coverage_cap_value = _core.integer(coverage_cap, "coverage_cap", minimum=1)
        arguments += ["--COVERAGE_CAP", str(coverage_cap_value)]
    if read_length is not None:
        read_length_value = _core.integer(read_length, "read_length", minimum=1)
        arguments += ["--READ_LENGTH", str(read_length_value)]
    if count_unpaired is not None:
        arguments += [
            "--COUNT_UNPAIRED",
            "true" if _core.flag(count_unpaired, "count_unpaired") else "false",
        ]
    if include_bq_histogram is not None:
        arguments += [
            "--INCLUDE_BQ_HISTOGRAM",
            "true" if _core.flag(include_bq_histogram, "include_bq_histogram") else "false",
        ]
    if use_fast_algorithm is not None:
        arguments += [
            "--USE_FAST_ALGORITHM",
            "true" if _core.flag(use_fast_algorithm, "use_fast_algorithm") else "false",
        ]
    if sample_size is not None:
        sample_size_value = _core.integer(sample_size, "sample_size", minimum=1)
        arguments += ["--SAMPLE_SIZE", str(sample_size_value)]
    if locus_accumulation_cap is not None:
        locus_accumulation_cap_value = _core.integer(
            locus_accumulation_cap, "locus_accumulation_cap", minimum=1
        )
        arguments += ["--LOCUS_ACCUMULATION_CAP", str(locus_accumulation_cap_value)]
    if theoretical_sensitivity_output is not None:
        theoretical_sensitivity_output_path = _core.new_file(
            theoretical_sensitivity_output, "theoretical_sensitivity_output"
        )
        arguments += ["--THEORETICAL_SENSITIVITY_OUTPUT", str(theoretical_sensitivity_output_path)]
        expected.append(theoretical_sensitivity_output_path)
        outs["theoretical_sensitivity_output"] = theoretical_sensitivity_output_path
    if stop_after is not None:
        stop_after_value = _core.integer(stop_after, "stop_after", minimum=0)
        arguments += ["--STOP_AFTER", str(stop_after_value)]
    tool = "CollectWgsMetrics"
    reserved = _core.names(
        "COUNT_UNPAIRED COVERAGE_CAP INCLUDE_BQ_HISTOGRAM INPUT INTERVALS "
        "LOCUS_ACCUMULATION_CAP MINIMUM_BASE_QUALITY MINIMUM_MAPPING_QUALITY "
        "OUTPUT READ_LENGTH REFERENCE_SEQUENCE SAMPLE_SIZE STOP_AFTER "
        "THEORETICAL_SENSITIVITY_OUTPUT USE_FAST_ALGORITHM"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_collect_wgs_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_calculate_fingerprint_metrics(
    haplotype_map: str,
    input_paths: Sequence[str],
    output_path: str,
    *,
    calculate_by: str | None = None,
    genotype_lod_threshold: float | None = None,
    number_of_sampling: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Calculate statistics on fingerprints (SAM/BAM/CRAM or VCF) against a haplotype map, checking
    their viability.

    Parameters (the GATK option is in parentheses):

    - ``haplotype_map`` (``--HAPLOTYPE_MAP``): Haplotype database (haplotype map) file.
    - ``input_paths`` (``--INPUT``): SAM/BAM/CRAM or VCF files (given as paths).
    - ``output_path`` (``--OUTPUT``): New metrics file.
    - ``calculate_by`` (``--CALCULATE_BY``): Aggregate the fingerprints by this level.
    - ``genotype_lod_threshold`` (``--GENOTYPE_LOD_THRESHOLD``): LOD threshold for a confident
      genotype.
    - ``number_of_sampling`` (``--NUMBER_OF_SAMPLING``): Number of samplings.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CalculateFingerprintMetrics --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``metrics``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths = [str(_core.existing_file(item, "input_paths")) for item in input_paths]
    haplotype_map_path = _core.existing_file(haplotype_map, "haplotype_map")
    arguments += ["--HAPLOTYPE_MAP", str(haplotype_map_path)]
    input_paths_values = _core.texts(input_paths, "input_paths")
    for item in input_paths_values:
        arguments += ["--INPUT", item]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if calculate_by is not None:
        calculate_by_value = _core.choice(
            calculate_by, ("READGROUP", "SAMPLE", "LIBRARY", "FILE"), "calculate_by"
        )
        arguments += ["--CALCULATE_BY", calculate_by_value]
    if genotype_lod_threshold is not None:
        genotype_lod_threshold_value = _core.number(
            genotype_lod_threshold, "genotype_lod_threshold"
        )
        arguments += ["--GENOTYPE_LOD_THRESHOLD", str(genotype_lod_threshold_value)]
    if number_of_sampling is not None:
        number_of_sampling_value = _core.integer(
            number_of_sampling, "number_of_sampling", minimum=1
        )
        arguments += ["--NUMBER_OF_SAMPLING", str(number_of_sampling_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CalculateFingerprintMetrics"
    reserved = _core.names(
        "CALCULATE_BY GENOTYPE_LOD_THRESHOLD HAPLOTYPE_MAP INPUT "
        "NUMBER_OF_SAMPLING OUTPUT REFERENCE_SEQUENCE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_calculate_fingerprint_metrics",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_calculate_read_group_checksum(
    input_path: str,
    output_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a hash code over the read groups (RG) of a SAM/BAM/CRAM file; files whose read-group
    content is the same give the same hash.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``output_path`` (``--OUTPUT``): New file for the hash code.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CalculateReadGroupChecksum --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``text`` and ``checksum`` (the hash code).
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
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CalculateReadGroupChecksum"
    reserved = _core.names("INPUT OUTPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_calculate_read_group_checksum",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        extras["text"] = _core.read_text(outs["output_path"]).strip()[:4000]
    if "text" in extras:
        extras["checksum"] = extras["text"].split()[0] if extras["text"] else None
    return _core.finish(result, outputs=outs, **extras)


def gatk_check_duplicate_marking(
    input_path: str,
    *,
    output_path: str | None = None,
    mode: str | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that reads of the same name (and mates) carry consistent duplicate markings.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM file.
    - ``output_path`` (``--OUTPUT``): New file listing the query names with inconsistent
      markings.
    - ``mode`` (``--MODE``): Which reads of the same name are compared.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CheckDuplicateMarking --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``has_findings`` (True when GATK reported inconsistent
    duplicate marks).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if mode is not None:
        mode_value = _core.choice(
            mode, ("ALL", "PRIMARY_ONLY", "PRIMARY_MAPPED_ONLY", "PRIMARY_PROPER_PAIR_ONLY"), "mode"
        )
        arguments += ["--MODE", mode_value]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CheckDuplicateMarking"
    reserved = _core.names("INPUT MODE OUTPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_check_duplicate_marking",
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
    return _core.finish(result, outputs=outs, **extras)


def gatk_check_fingerprint(
    input_path: str,
    genotypes: str,
    haplotype_map: str,
    *,
    output_prefix: str | None = None,
    summary_output: str | None = None,
    detail_output: str | None = None,
    expected_sample_alias: str | None = None,
    observed_sample_alias: str | None = None,
    genotype_lod_threshold: float | None = None,
    ignore_read_groups: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a fingerprint from a SAM/BAM/CRAM or VCF file and compare it with the provided
    genotypes (LOD of the same-sample hypothesis).

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM or VCF file.
    - ``genotypes`` (``--GENOTYPES``): VCF with the expected genotypes.
    - ``haplotype_map`` (``--HAPLOTYPE_MAP``): Haplotype database (haplotype map) file.
    - ``output_prefix`` (``--OUTPUT``): Prefix of the output files
      (.fingerprinting_summary_metrics and .fingerprinting_detail_metrics); alternative to
      summary_output and detail_output.
    - ``summary_output`` (``--SUMMARY_OUTPUT``): New summary metrics file (use with
      detail_output).
    - ``detail_output`` (``--DETAIL_OUTPUT``): New detail metrics file (use with
      summary_output).
    - ``expected_sample_alias`` (``--EXPECTED_SAMPLE_ALIAS``): Sample name in the genotypes file
      to compare with.
    - ``observed_sample_alias`` (``--OBSERVED_SAMPLE_ALIAS``): Sample name in the input VCF to
      use.
    - ``genotype_lod_threshold`` (``--GENOTYPE_LOD_THRESHOLD``): LOD threshold for a confident
      genotype.
    - ``ignore_read_groups`` (``--IGNORE_READ_GROUPS``): Pool all read groups.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CheckFingerprint --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``summary_metrics`` (LOD_EXPECTED_SAMPLE,
    LOD_EXPECTED_SAMPLE_TUMOR_NORMAL ... when the files exist).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if output_prefix is None and (summary_output is None or detail_output is None):
        raise ConfigurationError("Give output_prefix, or both summary_output and detail_output.")
    if output_prefix is not None and (summary_output is not None or detail_output is not None):
        raise ConfigurationError(
            "output_prefix cannot be combined with summary_output or detail_output."
        )
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    genotypes_path = _core.existing_file(genotypes, "genotypes")
    arguments += ["--GENOTYPES", str(genotypes_path)]
    haplotype_map_path = _core.existing_file(haplotype_map, "haplotype_map")
    arguments += ["--HAPLOTYPE_MAP", str(haplotype_map_path)]
    if output_prefix is not None:
        output_prefix_path = _core.new_prefix(output_prefix, "output_prefix")
        arguments += ["--OUTPUT", str(output_prefix_path)]
        outs["output_prefix"] = output_prefix_path
    if summary_output is not None:
        summary_output_path = _core.new_file(summary_output, "summary_output")
        arguments += ["--SUMMARY_OUTPUT", str(summary_output_path)]
        expected.append(summary_output_path)
        outs["summary_output"] = summary_output_path
    if detail_output is not None:
        detail_output_path = _core.new_file(detail_output, "detail_output")
        arguments += ["--DETAIL_OUTPUT", str(detail_output_path)]
        expected.append(detail_output_path)
        outs["detail_output"] = detail_output_path
    if expected_sample_alias is not None:
        expected_sample_alias_text = _core.text(expected_sample_alias, "expected_sample_alias")
        arguments += ["--EXPECTED_SAMPLE_ALIAS", expected_sample_alias_text]
    if observed_sample_alias is not None:
        observed_sample_alias_text = _core.text(observed_sample_alias, "observed_sample_alias")
        arguments += ["--OBSERVED_SAMPLE_ALIAS", observed_sample_alias_text]
    if genotype_lod_threshold is not None:
        genotype_lod_threshold_value = _core.number(
            genotype_lod_threshold, "genotype_lod_threshold"
        )
        arguments += ["--GENOTYPE_LOD_THRESHOLD", str(genotype_lod_threshold_value)]
    if ignore_read_groups is not None:
        arguments += [
            "--IGNORE_READ_GROUPS",
            "true" if _core.flag(ignore_read_groups, "ignore_read_groups") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CheckFingerprint"
    reserved = _core.names(
        "DETAIL_OUTPUT EXPECTED_SAMPLE_ALIAS GENOTYPES GENOTYPE_LOD_THRESHOLD "
        "HAPLOTYPE_MAP IGNORE_READ_GROUPS INPUT OBSERVED_SAMPLE_ALIAS OUTPUT "
        "REFERENCE_SEQUENCE SUMMARY_OUTPUT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_check_fingerprint",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    summary_file = None
    if output_prefix is not None:
        candidate = Path(str(output_prefix_path) + ".fingerprinting_summary_metrics")
        summary_file = candidate if candidate.is_file() else None
    elif summary_output is not None:
        summary_file = summary_output_path
    if summary_file is not None and summary_file.is_file():
        extras["summary_metrics"] = _core.picard_metrics(summary_file)["metrics"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_check_pileup(
    input_paths: Sequence[str],
    pileup: str,
    reference: str,
    *,
    output_path: str | None = None,
    intervals: Sequence[str] | None = None,
    ignore_overlaps: bool | None = None,
    continue_after_error: bool | None = None,
    max_depth_per_sample: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare GATK's internal pileup with a pileup made by samtools mpileup and report
    differences.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM files.
    - ``pileup`` (``--pileup``): Pileup file made by samtools mpileup.
    - ``reference`` (``--reference``): Reference FASTA with its .fai and .dict files.
    - ``output_path`` (``--output``): New report file; None prints to the log.
    - ``intervals`` (``--intervals``): Loci or interval files to restrict the comparison.
    - ``ignore_overlaps`` (``--ignore-overlaps``): Ignore overlapping bases of read pairs.
    - ``continue_after_error`` (``--continue-after-error``): Do not stop at the first
      difference.
    - ``max_depth_per_sample`` (``--max-depth-per-sample``): Downsample to this depth; 0 means
      no limit.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CheckPileup --help`` are accepted, a list repeats the option,
      and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``text`` (the report) when output_path is given.
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
    pileup_path = _core.existing_file(pileup, "pileup")
    arguments += ["--pileup", str(pileup_path)]
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--reference", str(reference_path)]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--output", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if ignore_overlaps is not None:
        arguments += [
            "--ignore-overlaps",
            "true" if _core.flag(ignore_overlaps, "ignore_overlaps") else "false",
        ]
    if continue_after_error is not None:
        arguments += [
            "--continue-after-error",
            "true" if _core.flag(continue_after_error, "continue_after_error") else "false",
        ]
    if max_depth_per_sample is not None:
        max_depth_per_sample_value = _core.integer(
            max_depth_per_sample, "max_depth_per_sample", minimum=0
        )
        arguments += ["--max-depth-per-sample", str(max_depth_per_sample_value)]
    tool = "CheckPileup"
    reserved = _core.names(
        "continue-after-error ignore-overlaps input intervals "
        "max-depth-per-sample output pileup reference"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_check_pileup",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        extras["text"] = _core.read_text(outs["output_path"]).strip()[:4000]
    return _core.finish(result, outputs=outs, **extras)


def gatk_check_terminator_block(
    input_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that the last block of a BGZF-compressed file (for example a BAM) is the well-formed
    terminator block.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Block-compressed file to check.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CheckTerminatorBlock --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``terminator_ok`` (False when GATK exits with 100, i.e. the last block is not a
    terminator) and ``has_findings``.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CheckTerminatorBlock"
    reserved = _core.names("INPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_check_terminator_block",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
        success_codes=(0, 100),
    )
    extras["has_findings"] = result["return_code"] != 0
    extras["terminator_ok"] = result["return_code"] == 0
    return _core.finish(result, outputs=outs, **extras)


def gatk_compare_base_qualities(
    input_paths: Sequence[str],
    *,
    output_path: str | None = None,
    throw_on_diff: bool | None = None,
    round_down_quantized: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the base qualities of two SAM/BAM/CRAM files that hold the same reads in the same
    order.

    Parameters (the GATK option is in parentheses):

    - ``input_paths``: Exactly two SAM/BAM/CRAM files with identical read names in identical
      order.
    - ``output_path`` (``--output``): New summary file.
    - ``throw_on_diff`` (``--throw-on-diff``): Fail when a difference is found.
    - ``round_down_quantized`` (``--round-down-quantized``): Round qualities down to the nearest
      quantized value.
    - ``reference`` (``--reference``): Reference FASTA (needed for CRAM).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CompareBaseQualities --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``text`` (the summary) when output_path is given.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if isinstance(input_paths, (str, bytes)) or len(input_paths) != 2:
        raise ConfigurationError("input_paths must hold exactly two files.")
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    positional += [str(item) for item in input_paths_paths]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--output", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if throw_on_diff is not None:
        arguments += [
            "--throw-on-diff",
            "true" if _core.flag(throw_on_diff, "throw_on_diff") else "false",
        ]
    if round_down_quantized is not None:
        arguments += [
            "--round-down-quantized",
            "true" if _core.flag(round_down_quantized, "round_down_quantized") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    tool = "CompareBaseQualities"
    reserved = _core.names("output reference round-down-quantized throw-on-diff")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_compare_base_qualities",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        extras["text"] = _core.read_text(outs["output_path"]).strip()[:4000]
    return _core.finish(result, outputs=outs, **extras)


def gatk_compare_duplicates_spark(
    input_paths: Sequence[str],
    input2: str,
    *,
    output_path: str | None = None,
    output2_path: str | None = None,
    throw_on_diff: bool | None = None,
    print_summary: bool | None = None,
    intervals: Sequence[str] | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Determine whether two potentially identical BAMs have the same duplicate reads (Spark tool,
    run locally).

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): First BAM/SAM/CRAM file(s).
    - ``input2`` (``--input2``): Second BAM/SAM/CRAM file.
    - ``output_path`` (``--output``): New BAM with the mismatching duplicate groups of the first
      input.
    - ``output2_path`` (``--output2``): New BAM with the mismatching duplicate groups of the
      second input.
    - ``throw_on_diff`` (``--throw-on-diff``): Fail when a difference is found.
    - ``print_summary`` (``--print-summary``): Print the summary to the log.
    - ``intervals`` (``--intervals``): Loci or interval files to restrict the comparison.
    - ``reference`` (``--reference``): Reference FASTA (needed for CRAM).
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk CompareDuplicatesSpark --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``summary_lines`` (the log lines about duplicate groups).
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
    input2_path = _core.existing_file(input2, "input2")
    arguments += ["--input2", str(input2_path)]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--output", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if output2_path is not None:
        output2_path_path = _core.new_file(output2_path, "output2_path")
        arguments += ["--output2", str(output2_path_path)]
        expected.append(output2_path_path)
        outs["output2_path"] = output2_path_path
    if throw_on_diff is not None:
        arguments += [
            "--throw-on-diff",
            "true" if _core.flag(throw_on_diff, "throw_on_diff") else "false",
        ]
    if print_summary is not None:
        arguments += [
            "--print-summary",
            "true" if _core.flag(print_summary, "print_summary") else "false",
        ]
    if intervals is not None:
        intervals_values = _core.intervals(intervals, "intervals")
        for item in intervals_values:
            arguments += ["--intervals", item]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    tool = "CompareDuplicatesSpark"
    reserved = _core.names(
        "input input2 intervals output output2 print-summary reference throw-on-diff"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_compare_duplicates_spark",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["summary_lines"] = [
        line
        for line in result["_log"].splitlines()
        if "duplicate" in line.lower() or "mismatch" in line.lower()
    ][:40]
    return _core.finish(result, outputs=outs, **extras)


def gatk_compare_sams(
    input_paths: Sequence[str],
    *,
    output_path: str | None = None,
    lenient_dup: bool | None = None,
    lenient_header: bool | None = None,
    lenient_low_mq_alignment: bool | None = None,
    lenient_unknown_mq_alignment: bool | None = None,
    low_mq_threshold: int | None = None,
    compare_mq: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two SAM/BAM/CRAM files (headers, alignments and duplicate marks) strictly or
    leniently and write the comparison.

    Parameters (the GATK option is in parentheses):

    - ``input_paths``: Exactly two SAM/BAM/CRAM files to compare.
    - ``output_path`` (``--OUTPUT``): New comparison table.
    - ``lenient_dup`` (``--LENIENT_DUP``): Allow differing duplicate marks.
    - ``lenient_header`` (``--LENIENT_HEADER``): Allow some header differences.
    - ``lenient_low_mq_alignment`` (``--LENIENT_LOW_MQ_ALIGNMENT``): Count reads below
      low_mq_threshold as equal even when aligned differently.
    - ``lenient_unknown_mq_alignment`` (``--LENIENT_UNKNOWN_MQ_ALIGNMENT``): Count reads with
      unknown mapping quality (255) as equal even when aligned differently.
    - ``low_mq_threshold`` (``--LOW_MQ_THRESHOLD``): Mapping quality below which lenient
      alignment mode applies.
    - ``compare_mq`` (``--COMPARE_MQ``): Also write a histogram of mapping-quality concordance.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CompareSAMs --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``files_match`` (GATK exit code 0) and ``has_findings`` (True when
    the files differ); ``text`` holds the comparison when output_path is given.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if isinstance(input_paths, (str, bytes)) or len(input_paths) != 2:
        raise ConfigurationError("input_paths must hold exactly two files.")
    input_paths_paths = _core.existing_files(input_paths, "input_paths")
    positional += [str(item) for item in input_paths_paths]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if lenient_dup is not None:
        arguments += [
            "--LENIENT_DUP",
            "true" if _core.flag(lenient_dup, "lenient_dup") else "false",
        ]
    if lenient_header is not None:
        arguments += [
            "--LENIENT_HEADER",
            "true" if _core.flag(lenient_header, "lenient_header") else "false",
        ]
    if lenient_low_mq_alignment is not None:
        arguments += [
            "--LENIENT_LOW_MQ_ALIGNMENT",
            "true" if _core.flag(lenient_low_mq_alignment, "lenient_low_mq_alignment") else "false",
        ]
    if lenient_unknown_mq_alignment is not None:
        arguments += [
            "--LENIENT_UNKNOWN_MQ_ALIGNMENT",
            "true"
            if _core.flag(lenient_unknown_mq_alignment, "lenient_unknown_mq_alignment")
            else "false",
        ]
    if low_mq_threshold is not None:
        low_mq_threshold_value = _core.integer(low_mq_threshold, "low_mq_threshold", minimum=0)
        arguments += ["--LOW_MQ_THRESHOLD", str(low_mq_threshold_value)]
    if compare_mq is not None:
        arguments += ["--COMPARE_MQ", "true" if _core.flag(compare_mq, "compare_mq") else "false"]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CompareSAMs"
    reserved = _core.names(
        "COMPARE_MQ LENIENT_DUP LENIENT_HEADER LENIENT_LOW_MQ_ALIGNMENT "
        "LENIENT_UNKNOWN_MQ_ALIGNMENT LOW_MQ_THRESHOLD OUTPUT "
        "REFERENCE_SEQUENCE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_compare_sams",
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
    if "output_path" in outs and outs["output_path"].is_file():
        extras["text"] = _core.read_text(outs["output_path"]).strip()[:4000]
    return _core.finish(result, outputs=outs, **extras)


def gatk_convert_haplotype_database_to_vcf(
    input_path: str,
    output_path: str,
    reference: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a haplotype database (haplotype map) file to a VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Haplotype database file.
    - ``output_path`` (``--OUTPUT``): New VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk ConvertHaplotypeDatabaseToVcf --help`` are accepted, a
      list repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``record_count`` (VCF records when the output is not
    compressed).
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
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "ConvertHaplotypeDatabaseToVcf"
    reserved = _core.names("INPUT OUTPUT REFERENCE_SEQUENCE")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_convert_haplotype_database_to_vcf",
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


def gatk_crosscheck_fingerprints(
    haplotype_map: str,
    input_paths: Sequence[str],
    *,
    output_path: str | None = None,
    second_input_paths: Sequence[str] | None = None,
    crosscheck_by: str | None = None,
    crosscheck_mode: str | None = None,
    lod_threshold: float | None = None,
    expect_all_groups_to_match: bool | None = None,
    output_errors_only: bool | None = None,
    matrix_output: str | None = None,
    num_threads: int | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that all data in the input files appear to have come from the same individual (or that
    groups that should match do).

    Parameters (the GATK option is in parentheses):

    - ``haplotype_map`` (``--HAPLOTYPE_MAP``): Haplotype database (haplotype map) file.
    - ``input_paths`` (``--INPUT``): SAM/BAM/CRAM or VCF files, or lists of files (given as
      paths).
    - ``output_path`` (``--OUTPUT``): New crosscheck metrics file.
    - ``second_input_paths`` (``--SECOND_INPUT``): Second set of files to crosscheck against the
      first.
    - ``crosscheck_by`` (``--CROSSCHECK_BY``): Level at which fingerprints are crosschecked.
    - ``crosscheck_mode`` (``--CROSSCHECK_MODE``): How INPUT and SECOND_INPUT are crosschecked.
    - ``lod_threshold`` (``--LOD_THRESHOLD``): LOD threshold between expected match and
      mismatch.
    - ``expect_all_groups_to_match`` (``--EXPECT_ALL_GROUPS_TO_MATCH``): Expect all groups to
      match each other.
    - ``output_errors_only`` (``--OUTPUT_ERRORS_ONLY``): Write only the failing comparisons.
    - ``matrix_output`` (``--MATRIX_OUTPUT``): New matrix of LOD scores.
    - ``num_threads`` (``--NUM_THREADS``): Threads.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk CrosscheckFingerprints --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``metrics`` (one row per compared pair with RESULT and LOD scores)
    and ``has_findings`` (True when GATK exits with 1 because of a mismatch).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_paths = [str(_core.existing_file(item, "input_paths")) for item in input_paths]
    if second_input_paths is not None:
        second_input_paths = [
            str(_core.existing_file(item, "second_input_paths")) for item in second_input_paths
        ]
    haplotype_map_path = _core.existing_file(haplotype_map, "haplotype_map")
    arguments += ["--HAPLOTYPE_MAP", str(haplotype_map_path)]
    input_paths_values = _core.texts(input_paths, "input_paths")
    for item in input_paths_values:
        arguments += ["--INPUT", item]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if second_input_paths is not None:
        second_input_paths_values = _core.texts(second_input_paths, "second_input_paths")
        for item in second_input_paths_values:
            arguments += ["--SECOND_INPUT", item]
    if crosscheck_by is not None:
        crosscheck_by_value = _core.choice(
            crosscheck_by, ("READGROUP", "SAMPLE", "LIBRARY", "FILE"), "crosscheck_by"
        )
        arguments += ["--CROSSCHECK_BY", crosscheck_by_value]
    if crosscheck_mode is not None:
        crosscheck_mode_value = _core.choice(
            crosscheck_mode, ("CHECK_SAME_SAMPLE", "CHECK_ALL_OTHERS"), "crosscheck_mode"
        )
        arguments += ["--CROSSCHECK_MODE", crosscheck_mode_value]
    if lod_threshold is not None:
        lod_threshold_value = _core.number(lod_threshold, "lod_threshold")
        arguments += ["--LOD_THRESHOLD", str(lod_threshold_value)]
    if expect_all_groups_to_match is not None:
        arguments += [
            "--EXPECT_ALL_GROUPS_TO_MATCH",
            "true"
            if _core.flag(expect_all_groups_to_match, "expect_all_groups_to_match")
            else "false",
        ]
    if output_errors_only is not None:
        arguments += [
            "--OUTPUT_ERRORS_ONLY",
            "true" if _core.flag(output_errors_only, "output_errors_only") else "false",
        ]
    if matrix_output is not None:
        matrix_output_path = _core.new_file(matrix_output, "matrix_output")
        arguments += ["--MATRIX_OUTPUT", str(matrix_output_path)]
        expected.append(matrix_output_path)
        outs["matrix_output"] = matrix_output_path
    if num_threads is not None:
        num_threads_value = _core.integer(num_threads, "num_threads", minimum=1)
        arguments += ["--NUM_THREADS", str(num_threads_value)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "CrosscheckFingerprints"
    reserved = _core.names(
        "CROSSCHECK_BY CROSSCHECK_MODE EXPECT_ALL_GROUPS_TO_MATCH HAPLOTYPE_MAP "
        "INPUT LOD_THRESHOLD MATRIX_OUTPUT NUM_THREADS OUTPUT "
        "OUTPUT_ERRORS_ONLY REFERENCE_SEQUENCE SECOND_INPUT"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_crosscheck_fingerprints",
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
    if "output_path" in outs and outs["output_path"].is_file():
        parsed = _core.picard_metrics(outs["output_path"])
        extras["metrics"] = parsed["metrics"]
        extras["histogram"] = parsed["histogram"]
    return _core.finish(result, outputs=outs, **extras)


def gatk_dump_tabix_index(
    tabix_index: str,
    output_path: str,
    *,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dump a tabix index file (.tbi) as text: sequence names, bins, chunks and linear index.

    Parameters (the GATK option is in parentheses):

    - ``tabix_index`` (``--tabix-index``): Tabix index file.
    - ``output_path`` (``--output``): New text file.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk DumpTabixIndex --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``text`` (first 4,000 characters of the dump).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    tabix_index_path = _core.existing_file(tabix_index, "tabix_index")
    arguments += ["--tabix-index", str(tabix_index_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--output", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    tool = "DumpTabixIndex"
    reserved = _core.names("output tabix-index")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_dump_tabix_index",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    if "output_path" in outs and outs["output_path"].is_file():
        extras["text"] = _core.read_text(outs["output_path"]).strip()[:4000]
    return _core.finish(result, outputs=outs, **extras)


def gatk_extract_fingerprint(
    input_path: str,
    haplotype_map: str,
    output_path: str,
    reference: str,
    *,
    contamination: float | None = None,
    sample_alias: str | None = None,
    locus_max_reads: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a fingerprint from a SAM/BAM/CRAM file against a haplotype map and write it as a
    VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``haplotype_map`` (``--HAPLOTYPE_MAP``): Haplotype database (haplotype map) file.
    - ``output_path`` (``--OUTPUT``): New fingerprint VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``contamination`` (``--CONTAMINATION``): Contamination fraction assumed when computing the
      fingerprint.
    - ``sample_alias`` (``--SAMPLE_ALIAS``): Sample name written to the VCF.
    - ``locus_max_reads`` (``--LOCUS_MAX_READS``): Maximum reads used per locus.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk ExtractFingerprint --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
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
    haplotype_map_path = _core.existing_file(haplotype_map, "haplotype_map")
    arguments += ["--HAPLOTYPE_MAP", str(haplotype_map_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if contamination is not None:
        contamination_value = _core.number(contamination, "contamination", minimum=0.0, maximum=1.0)
        arguments += ["--CONTAMINATION", str(contamination_value)]
    if sample_alias is not None:
        sample_alias_text = _core.text(sample_alias, "sample_alias")
        arguments += ["--SAMPLE_ALIAS", sample_alias_text]
    if locus_max_reads is not None:
        locus_max_reads_value = _core.integer(locus_max_reads, "locus_max_reads", minimum=1)
        arguments += ["--LOCUS_MAX_READS", str(locus_max_reads_value)]
    tool = "ExtractFingerprint"
    reserved = _core.names(
        "CONTAMINATION HAPLOTYPE_MAP INPUT LOCUS_MAX_READS OUTPUT REFERENCE_SEQUENCE SAMPLE_ALIAS"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_extract_fingerprint",
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


def gatk_get_sample_name(
    input_paths: Sequence[str],
    output_path: str,
    *,
    use_url_encoding: bool | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Emit the single sample name of a SAM/BAM/CRAM file to a file.

    Parameters (the GATK option is in parentheses):

    - ``input_paths`` (``--input``): Input SAM/BAM/CRAM file(s).
    - ``output_path`` (``--output``): New file that holds only the sample name.
    - ``use_url_encoding`` (``--use-url-encoding``): URL-encode spaces and special characters.
    - ``reference`` (``--reference``): Reference FASTA (needed for CRAM).
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk GetSampleName --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``sample_name``.
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
    if use_url_encoding is not None:
        arguments += [
            "--use-url-encoding",
            "true" if _core.flag(use_url_encoding, "use_url_encoding") else "false",
        ]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--reference", str(reference_path)]
    tool = "GetSampleName"
    reserved = _core.names("input output reference use-url-encoding")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_get_sample_name",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["sample_name"] = _core.read_text(outs["output_path"]).strip()
    return _core.finish(result, outputs=outs, **extras)


def gatk_identify_contaminant(
    input_path: str,
    haplotype_map: str,
    output_path: str,
    reference: str,
    *,
    contamination: float | None = None,
    extract_contaminated: bool | None = None,
    sample_alias: str | None = None,
    locus_max_reads: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute a fingerprint of the contaminating sample from a SAM/BAM/CRAM file, given a
    contamination estimate, and write it as a VCF.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``haplotype_map`` (``--HAPLOTYPE_MAP``): Haplotype database (haplotype map) file.
    - ``output_path`` (``--OUTPUT``): New fingerprint VCF.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files.
    - ``contamination`` (``--CONTAMINATION``): Contamination fraction of the sample.
    - ``extract_contaminated`` (``--EXTRACT_CONTAMINATED``): Extract the contaminated (majority)
      fingerprint instead of the contaminant.
    - ``sample_alias`` (``--SAMPLE_ALIAS``): Sample name written to the VCF.
    - ``locus_max_reads`` (``--LOCUS_MAX_READS``): Maximum reads used per locus.
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk IdentifyContaminant --help`` are accepted, a list
      repeats the option, and options already set by a parameter above are refused.
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
    haplotype_map_path = _core.existing_file(haplotype_map, "haplotype_map")
    arguments += ["--HAPLOTYPE_MAP", str(haplotype_map_path)]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--OUTPUT", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    reference_path = _core.existing_file(reference, "reference")
    arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    if contamination is not None:
        contamination_value = _core.number(contamination, "contamination", minimum=0.0, maximum=1.0)
        arguments += ["--CONTAMINATION", str(contamination_value)]
    if extract_contaminated is not None:
        arguments += [
            "--EXTRACT_CONTAMINATED",
            "true" if _core.flag(extract_contaminated, "extract_contaminated") else "false",
        ]
    if sample_alias is not None:
        sample_alias_text = _core.text(sample_alias, "sample_alias")
        arguments += ["--SAMPLE_ALIAS", sample_alias_text]
    if locus_max_reads is not None:
        locus_max_reads_value = _core.integer(locus_max_reads, "locus_max_reads", minimum=1)
        arguments += ["--LOCUS_MAX_READS", str(locus_max_reads_value)]
    tool = "IdentifyContaminant"
    reserved = _core.names(
        "CONTAMINATION EXTRACT_CONTAMINATED HAPLOTYPE_MAP INPUT LOCUS_MAX_READS "
        "OUTPUT REFERENCE_SEQUENCE SAMPLE_ALIAS"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_identify_contaminant",
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


def gatk_lift_over_haplotype_map(
    input_path: str,
    chain: str,
    sequence_dictionary: str,
    output_path: str,
    *,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lift over a haplotype database (haplotype map) from one reference build to another with a
    chain file.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Haplotype database on the old reference.
    - ``chain`` (``--CHAIN``): Chain file from the old to the new reference.
    - ``sequence_dictionary`` (``--SEQUENCE_DICTIONARY``): Sequence dictionary of the new
      reference.
    - ``output_path`` (``--OUTPUT``): New haplotype database.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk LiftOverHaplotypeMap --help`` are accepted, a list
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
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "LiftOverHaplotypeMap"
    reserved = _core.names("CHAIN INPUT OUTPUT REFERENCE_SEQUENCE SEQUENCE_DICTIONARY")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_lift_over_haplotype_map",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_validate_sam_file(
    input_path: str,
    *,
    mode: str | None = None,
    output_path: str | None = None,
    ignore: Sequence[str] | None = None,
    ignore_warnings: bool | None = None,
    max_output: int | None = None,
    is_bisulfite_sequenced: bool | None = None,
    skip_mate_validation: bool | None = None,
    index_validation_stringency: str | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate a SAM/BAM/CRAM file (headers, flags, mates, CIGAR, indices) and report errors and
    warnings.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``mode`` (``--MODE``): VERBOSE lists every problem, SUMMARY counts them.
    - ``output_path`` (``--OUTPUT``): New report file; None reports in the log.
    - ``ignore`` (``--IGNORE``): Validation error types to ignore, for example
      MISSING_READ_GROUP.
    - ``ignore_warnings`` (``--IGNORE_WARNINGS``): Report errors only.
    - ``max_output`` (``--MAX_OUTPUT``): Maximum lines reported in VERBOSE mode.
    - ``is_bisulfite_sequenced`` (``--IS_BISULFITE_SEQUENCED``): Do not count C->T as a
      mismatch.
    - ``skip_mate_validation`` (``--SKIP_MATE_VALIDATION``): Do not validate mate information.
    - ``index_validation_stringency`` (``--INDEX_VALIDATION_STRINGENCY``): How thoroughly an
      existing BAM index is validated.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk ValidateSamFile --help`` are accepted, a list repeats
      the option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``, ``valid`` (no error found), ``has_findings`` and ``problems`` (the
    ERROR and WARNING lines).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    if mode is not None:
        mode_value = _core.choice(mode, ("VERBOSE", "SUMMARY"), "mode")
        arguments += ["--MODE", mode_value]
    if output_path is not None:
        output_path_path = _core.new_file(output_path, "output_path")
        arguments += ["--OUTPUT", str(output_path_path)]
        expected.append(output_path_path)
        outs["output_path"] = output_path_path
    if ignore is not None:
        ignore_values = _core.texts(ignore, "ignore")
        for item in ignore_values:
            arguments += ["--IGNORE", item]
    if ignore_warnings is not None:
        arguments += [
            "--IGNORE_WARNINGS",
            "true" if _core.flag(ignore_warnings, "ignore_warnings") else "false",
        ]
    if max_output is not None:
        max_output_value = _core.integer(max_output, "max_output", minimum=1)
        arguments += ["--MAX_OUTPUT", str(max_output_value)]
    if is_bisulfite_sequenced is not None:
        arguments += [
            "--IS_BISULFITE_SEQUENCED",
            "true" if _core.flag(is_bisulfite_sequenced, "is_bisulfite_sequenced") else "false",
        ]
    if skip_mate_validation is not None:
        arguments += [
            "--SKIP_MATE_VALIDATION",
            "true" if _core.flag(skip_mate_validation, "skip_mate_validation") else "false",
        ]
    if index_validation_stringency is not None:
        index_validation_stringency_value = _core.choice(
            index_validation_stringency,
            ("EXHAUSTIVE", "LESS_EXHAUSTIVE", "NONE"),
            "index_validation_stringency",
        )
        arguments += ["--INDEX_VALIDATION_STRINGENCY", index_validation_stringency_value]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "ValidateSamFile"
    reserved = _core.names(
        "IGNORE IGNORE_WARNINGS INDEX_VALIDATION_STRINGENCY INPUT "
        "IS_BISULFITE_SEQUENCED MAX_OUTPUT MODE OUTPUT REFERENCE_SEQUENCE "
        "SKIP_MATE_VALIDATION"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_validate_sam_file",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
        success_codes=(0, 1, 4),
    )
    extras["has_findings"] = result["return_code"] != 0
    report = (
        _core.read_text(outs["output_path"])
        if "output_path" in outs and outs["output_path"].is_file()
        else result["_log"]
    )
    extras["valid"] = result["return_code"] == 0
    extras["problems"] = [
        line.strip() for line in report.splitlines() if line.startswith(("ERROR", "WARNING"))
    ][:200]
    return _core.finish(result, outputs=outs, **extras)


def gatk_view_sam(
    input_path: str,
    alignment_status: str,
    pf_status: str,
    *,
    header_only: bool | None = None,
    records_only: bool | None = None,
    interval_list: str | None = None,
    reference: str | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Print the records (and header) of a SAM/BAM/CRAM file, optionally only aligned or unaligned
    and only PF or non-PF reads.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--INPUT``): Input SAM/BAM/CRAM file.
    - ``alignment_status`` (``--ALIGNMENT_STATUS``): Which reads to print.
    - ``pf_status`` (``--PF_STATUS``): Which reads to print by vendor filter.
    - ``header_only`` (``--HEADER_ONLY``): Print the header only.
    - ``records_only`` (``--RECORDS_ONLY``): Print the records only.
    - ``interval_list`` (``--INTERVAL_LIST``): Interval list that restricts the records.
    - ``reference`` (``--REFERENCE_SEQUENCE``): Reference FASTA with its .fai and .dict files;
      None uses the GATK default (none).
    - ``options``: more long options of this tool by name, for example ``{'MAX_RECORDS_IN_RAM':
      ...}``; only options listed by ``gatk ViewSam --help`` are accepted, a list repeats the
      option, and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``header_lines``, ``record_count`` and ``records`` (at most 200 records); the
    full text is in the log of the run.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--INPUT", str(input_path_path)]
    alignment_status_value = _core.choice(
        alignment_status, ("Aligned", "Unaligned", "All"), "alignment_status"
    )
    arguments += ["--ALIGNMENT_STATUS", alignment_status_value]
    pf_status_value = _core.choice(pf_status, ("PF", "NonPF", "All"), "pf_status")
    arguments += ["--PF_STATUS", pf_status_value]
    if header_only is not None:
        arguments += [
            "--HEADER_ONLY",
            "true" if _core.flag(header_only, "header_only") else "false",
        ]
    if records_only is not None:
        arguments += [
            "--RECORDS_ONLY",
            "true" if _core.flag(records_only, "records_only") else "false",
        ]
    if interval_list is not None:
        interval_list_path = _core.existing_file(interval_list, "interval_list")
        arguments += ["--INTERVAL_LIST", str(interval_list_path)]
    if reference is not None:
        reference_path = _core.existing_file(reference, "reference")
        arguments += ["--REFERENCE_SEQUENCE", str(reference_path)]
    tool = "ViewSam"
    reserved = _core.names(
        "ALIGNMENT_STATUS HEADER_ONLY INPUT INTERVAL_LIST PF_STATUS RECORDS_ONLY REFERENCE_SEQUENCE"
    )
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_view_sam",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    lines = [
        line
        for line in result["_log"].splitlines()
        if line.startswith("@") or line.count("\t") >= 10
    ]
    extras["header_lines"] = [line for line in lines if line.startswith("@")][:200]
    extras["record_count"] = sum(1 for line in lines if not line.startswith("@"))
    extras["records"] = [line for line in lines if not line.startswith("@")][:200]
    return _core.finish(result, outputs=outs, **extras)


def gatk_htsget_download(
    url: str,
    record_id: str,
    output_path: str,
    *,
    data_format: str | None = None,
    data_class: str | None = None,
    intervals: str | None = None,
    check_md5: bool | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a file from a GA4GH htsget endpoint (a record id on a server URL) to a local file;
    needs network access.

    Parameters (the GATK option is in parentheses):

    - ``url`` (``--url``): URL of the htsget endpoint.
    - ``record_id`` (``--id``): Identifier of the record to request.
    - ``output_path`` (``--O``): New local file.
    - ``data_format`` (``--format``): Requested data format.
    - ``data_class`` (``--class``): Request the whole file or only the header.
    - ``intervals`` (``--intervals``): Region to request, such as chr1:1000-2000.
    - ``check_md5`` (``--check-md5``): Verify the MD5 sum of the downloaded data.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk HtsgetReader --help`` are accepted, a list repeats the option,
      and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs``. This function uses the network.
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    if not (url.startswith("http://") or url.startswith("https://")):
        raise ConfigurationError("url must be an http(s) URL of an htsget endpoint.")
    url_text = _core.text(url, "url")
    arguments += ["--url", url_text]
    record_id_text = _core.text(record_id, "record_id")
    arguments += ["--id", record_id_text]
    output_path_path = _core.new_file(output_path, "output_path")
    arguments += ["--O", str(output_path_path)]
    expected.append(output_path_path)
    outs["output_path"] = output_path_path
    if data_format is not None:
        data_format_value = _core.choice(data_format, ("BAM", "CRAM", "VCF", "BCF"), "data_format")
        arguments += ["--format", data_format_value]
    if data_class is not None:
        data_class_value = _core.choice(data_class, ("body", "header"), "data_class")
        arguments += ["--class", data_class_value]
    if intervals is not None:
        intervals_text = _core.interval(intervals, "intervals")
        arguments += ["--intervals", intervals_text]
    if check_md5 is not None:
        arguments += ["--check-md5", "true" if _core.flag(check_md5, "check_md5") else "false"]
    tool = "HtsgetReader"
    reserved = _core.names("O check-md5 class format id intervals url")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_htsget_download",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    return _core.finish(result, outputs=outs, **extras)


def gatk_split_cram(
    input_path: str,
    output_dir: str,
    *,
    shard_records: int | None = None,
    shard_max_output_count: int | None = None,
    options: Mapping[str, Any] | None = None,
    gatk_path: str | None = None,
    java_max_memory: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a CRAM file into smaller CRAM files of a given number of records, without
    recompressing the containers.

    Parameters (the GATK option is in parentheses):

    - ``input_path`` (``--input``): Input CRAM file.
    - ``output_dir`` (``--None``): New folder that receives the shards shard_0000.cram,
      shard_0001.cram ...
    - ``shard_records`` (``--shard-records``): Records per shard.
    - ``shard_max_output_count`` (``--shard-max-output-count``): Stop after this many shards; 0
      means no limit.
    - ``options``: more long options of this tool by name, for example ``{'verbosity': ...}``;
      only options listed by ``gatk SplitCRAM --help`` are accepted, a list repeats the option,
      and options already set by a parameter above are refused.
    - ``gatk_path``: path of the ``gatk`` launcher (default: ``DNAKIT_GATK``, then PATH).
    - ``java_max_memory``: Java heap such as ``'4g'``, passed as ``--java-options -Xmx4g``.

    Result keys: ``outputs`` and ``shards`` (the CRAM files written).
    """
    arguments: list[str] = []
    positional: list[str] = []
    expected: list[Path] = []
    made_directories: list[Path] = []
    outs: dict[str, Path] = {}
    extras: dict[str, Any] = {}
    input_path_path = _core.existing_file(input_path, "input_path")
    arguments += ["--input", str(input_path_path)]
    output_dir_path = _core.new_directory(output_dir, "output_dir")
    made_directories.append(output_dir_path)
    outs["output_dir"] = output_dir_path
    if shard_records is not None:
        shard_records_value = _core.integer(shard_records, "shard_records", minimum=1)
        arguments += ["--shard-records", str(shard_records_value)]
    if shard_max_output_count is not None:
        shard_max_output_count_value = _core.integer(
            shard_max_output_count, "shard_max_output_count", minimum=0
        )
        arguments += ["--shard-max-output-count", str(shard_max_output_count_value)]
    tool = "SplitCRAM"
    arguments += ["--output", str(output_dir_path / "shard_%04d.cram")]
    reserved = _core.names("input shard-max-output-count shard-records")
    arguments = [*positional, *arguments, *_core.option_arguments(tool, options, reserved=reserved)]
    cwd = _core.working_directory(expected, outs)
    for directory in made_directories:
        directory.mkdir()
    result = _core.run_tool(
        "gatk_split_cram",
        tool,
        arguments,
        cwd=cwd,
        expected=expected,
        gatk_path=gatk_path,
        java_max_memory=java_max_memory,
        verbose=verbose,
    )
    extras["shards"] = sorted(str(item) for item in output_dir_path.glob("shard_*.cram"))
    return _core.finish(result, outputs=outs, **extras)
