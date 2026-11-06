"""Shasta command-line facades: long-read de novo assembly with every pipeline stage selectable, and readers of its run folder.

The user installs Shasta (a release binary from https://github.com/paoloshasta/shasta/releases, or a source build; 0.14.0
is the tested version) and puts it on PATH as ``shasta``, sets ``DNAKIT_SHASTA`` or passes ``shasta_path``. The Shasta
*Python API* (``import shasta``, class ``Assembler``) is not wrapped: it is a developer interface of the source build whose
methods are single stages that work on the in-memory state of a previous stage, the release shared library needs Python 3.10
and Ubuntu 22.04 system libraries, and the command line does not expose the stages one by one. Each stage is instead reached
through the options of ``shasta --command assemble`` (``shasta_assemble``), and the files Shasta writes for the stages
(summary, histograms, read graph components, assembled FASTA and GFA, reads by segment) are read back by
``shasta_read_run_table`` and ``shasta_read_assembly_output``.

No shell is used, the executable is resolved from an explicit path, ``DNAKIT_SHASTA`` or PATH, option names are checked
against the table of Shasta 0.14.0 (``_shasta_p29_options.py``), and an option is passed only when it is given, so the
chosen ``config`` and Shasta's defaults apply otherwise. Results are JSON-serialisable dictionaries.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.assembly._shasta_p29_options import OPTIONS as _OPTIONS
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError

_PROGRAM = "shasta"
_ENV = "DNAKIT_SHASTA"
_HINT = (
    "Install Shasta 0.14 or later: download the release binary from https://github.com/paoloshasta/shasta/releases "
    "(for example shasta-Linux-0.14.0), make it executable, then put it on PATH as `shasta`, set DNAKIT_SHASTA or "
    "pass shasta_path."
)
_CONFIG_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.\-]{0,199}")
_FASTA_LIKE = re.compile(r"[A-Za-z0-9_.\-]+")
_MAX_ROWS = 1_000_000

#: run-folder tables readable by ``shasta_read_run_table``: name -> (file name, what it holds)
TABLES: dict[str, tuple[str, str]] = {
    "read_summary": (
        "ReadSummary.csv",
        "per read: lengths, marker count, Palindromic, Chimeric, alignment candidates, read graph neighbours",
    ),
    "read_low_hash_statistics": (
        "ReadLowHashStatistics.csv",
        "per read: MinHash/LowHash feature statistics",
    ),
    "read_length_histogram": (
        "ReadLengthHistogram.csv",
        "read length histogram with cumulative bases",
    ),
    "binned_read_length_histogram": (
        "Binned-ReadLengthHistogram.csv",
        "binned read length histogram",
    ),
    "duplicate_reads": (
        "DuplicateReads.csv",
        "reads with duplicate names and whether they were discarded",
    ),
    "low_hash_bucket_histogram": (
        "LowHashBucketHistogram.csv",
        "MinHash/LowHash bucket sizes per iteration",
    ),
    "suppressed_alignment_candidates": (
        "SuppressedAlignmentCandidates.csv",
        "alignment candidates suppressed during candidate search",
    ),
    "aligned_fraction_histogram": (
        "AlignedFractionHistogram.csv",
        "histogram of the aligned fraction of marker alignments",
    ),
    "alignment_marker_count_histogram": (
        "AlignmentMarkerCountHistogram.csv",
        "histogram of aligned marker counts",
    ),
    "alignment_skip_histogram": (
        "AlignmentSkipHistogram.csv",
        "histogram of the largest marker skip of alignments",
    ),
    "alignment_drift_histogram": (
        "AlignmentDriftHistogram.csv",
        "histogram of the largest marker drift of alignments",
    ),
    "alignment_trim_histogram": (
        "AlignmentTrimHistogram.csv",
        "histogram of the unaligned marker count at the alignment ends",
    ),
    "read_graph_components": ("ReadGraphComponents.csv", "connected components of the read graph"),
    "disjoint_sets_histogram": (
        "DisjointSetsHistogram.csv",
        "coverage histogram of the marker disjoint sets before marker graph vertex creation",
    ),
    "marker_graph_vertex_coverage_histogram": (
        "MarkerGraphVertexCoverageHistogram.csv",
        "marker graph vertex coverage histogram",
    ),
    "marker_graph_edge_coverage_histogram": (
        "MarkerGraphEdgeCoverageHistogram.csv",
        "marker graph edge coverage histogram",
    ),
    "anchor_coverage_histogram": (
        "AnchorCoverageHistogram.csv",
        "Mode 3 anchor coverage histogram",
    ),
    "reads_by_segment": (
        "ReadsBySegment.csv",
        "oriented reads contributing to each assembled segment (--write_reads_by_segment)",
    ),
    "assembly_chains": (
        "Assembly.csv",
        "Mode 3 assembly chains: bubble chain, sequence length, coverage",
    ),
    "components": ("Components.csv", "Mode 3 components: reads, segments, sequence, N50"),
    "oriented_reads_by_component": (
        "OrientedReadsByComponent.csv",
        "Mode 3 oriented reads of each component",
    ),
    "journeys": ("Journeys.csv", "Mode 3 journeys of the oriented reads on the anchor graph"),
}
_ASSEMBLY_OUTPUTS: dict[str, str] = {
    "fasta": "Assembly.fasta",
    "gfa": "Assembly.gfa",
    "gfa_both_strands": "Assembly-BothStrands.gfa",
    "gfa_no_sequence": "Assembly-BothStrands-NoSequence.gfa",
}


def _resolve(shasta_path: str | None) -> str | None:
    if shasta_path is not None:
        return shasta_path
    return os.environ.get(_ENV) or None


def _opt_int(
    name: str, value: int | None, flag: str, *, minimum: int | None = 0, maximum: int | None = None
) -> list[str]:
    if value is None:
        return []
    return [flag, str(_cli.integer(value, name, minimum=minimum, maximum=maximum))]


def _opt_float(
    name: str, value: float | None, flag: str, *, minimum: float = 0.0, maximum: float | None = None
) -> list[str]:
    if value is None:
        return []
    return [flag, f"{_cli.number(value, name, minimum=minimum, maximum=maximum):g}"]


def _opt_bool(name: str, value: bool | None, flag: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True, False or None.")
    return [flag] if value else []


def _opt_choice(name: str, value: int | None, flag: str, allowed: Sequence[int]) -> list[str]:
    if value is None:
        return []
    if isinstance(value, bool) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(str(a) for a in allowed)}.")
    return [flag, str(value)]


def _checked_options(options: Mapping[str, Any] | None, explicit_flags: set[str]) -> list[str]:
    """Translate ``{"Section.name": value}`` into flags, refusing unknown names and values of the wrong kind."""
    arguments: list[str] = []
    for key, value in (options or {}).items():
        if not isinstance(key, str) or key not in _OPTIONS:
            raise ConfigurationError(
                f"Unknown Shasta option {key!r}; the accepted names are Section.name as in `shasta --help` (0.14.0)."
            )
        flag = f"--{key}"
        if flag in explicit_flags:
            raise ConfigurationError(f"{key} is already set by a named parameter of this function.")
        kind = _OPTIONS[key]
        if kind == "flag":
            if not isinstance(value, bool):
                raise ConfigurationError(f"{key} is a switch; pass True or False.")
            if value:
                arguments.append(flag)
            continue
        if isinstance(value, bool) or value is None:
            raise ConfigurationError(f"{key} needs a value.")
        if kind == "int":
            if isinstance(value, float) and value.is_integer():
                value = int(value)
            arguments.extend([flag, str(_cli.integer(value, key))])
        elif kind == "float":
            arguments.extend([flag, f"{_cli.number(value, key):g}"])
        else:
            arguments.extend([flag, _cli.text(str(value), key)])
    return arguments


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _summary_fields(summary: dict[str, Any] | None) -> dict[str, Any]:
    if summary is None:
        return {}
    reads = summary.get("Reads used in this assembly", {})
    segments = summary.get("Assembled segments", {})
    return {
        "reads_used": reads.get("Number of reads"),
        "read_bases": reads.get("Number of read bases"),
        "read_n50": reads.get("Read N50"),
        "reads_flagged_palindromic": reads.get(
            "Number of reads flagged as palindromic by self alignment"
        ),
        "reads_flagged_chimeric": reads.get("Number of reads flagged as chimeric"),
        "assembled_segments": segments.get("Number of segments assembled"),
        "assembled_bases": segments.get("Total assembled segment length"),
        "segment_n50": segments.get("Assembled segments N50"),
        "summary_sections": sorted(summary),
    }


def shasta_assemble(
    reads: Sequence[str],
    output_directory: str,
    *,
    config: str = "Nanopore-May2022",
    assembly_mode: int | None = None,
    anchors_path: str | None = None,
    min_read_length: int | None = None,
    desired_coverage: str | int | None = None,
    palindromic_flagging: bool | None = None,
    kmer_generation_method: int | None = None,
    kmer_length: int | None = None,
    kmer_probability: float | None = None,
    minhash_m: int | None = None,
    minhash_hash_fraction: float | None = None,
    minhash_iteration_count: int | None = None,
    minhash_all_pairs: bool | None = None,
    align_method: int | None = None,
    read_graph_creation_method: int | None = None,
    strand_separation_method: int | None = None,
    max_chimeric_read_distance: int | None = None,
    flag_inconsistent_alignments: bool | None = None,
    marker_graph_min_coverage: int | None = None,
    marker_graph_max_coverage: int | None = None,
    marker_graph_min_edge_coverage: int | None = None,
    cleanup_duplicate_markers: bool | None = None,
    simplify_max_length: str | Sequence[int] | None = None,
    detangle_method: int | None = None,
    consensus_caller: str | None = None,
    write_reads_by_segment: bool | None = None,
    store_coverage_data: bool | None = None,
    memory_mode: str | None = None,
    memory_backing: str | None = None,
    options: Mapping[str, Any] | None = None,
    threads: int | None = None,
    shasta_path: str | None = None,
    timeout_seconds: float = 86_400.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """De novo assembly of long reads with Shasta, with every assembly stage selectable (``shasta --command assemble``).

    ``reads`` are FASTA/FASTQ files (``--input``) and ``output_directory`` is the new run folder (``--assemblyDirectory``;
    Shasta refuses an existing folder). ``config`` is a built-in configuration name (for example ``Nanopore-May2022``,
    ``Nanopore-Phased-R10-Slow-Nov2022``, ``Nanopore-ncm23-May2024`` or ``HiFi-Oct2021``) or the path of a configuration
    file (``--config``). Every parameter below is passed only when it is not ``None``, so the configuration decides
    otherwise. One call covers these Shasta stages, each selected by a parameter:

    - ``assembly_mode``: 0 (default marker-graph assembly), 2 (phased, needs ``strand_separation_method=2``) or 3 (Mode 3
      phased assembly from anchors, ``--Assembly.mode 3``; needs a raw-sequence configuration such as
      ``Nanopore-ncm23-May2024``); ``anchors_path`` is ``--anchors`` for Mode 3 assembly from a given anchor file.
    - Reads: ``min_read_length``, ``desired_coverage`` (bases, for example ``"120Gb"``), ``palindromic_flagging`` (False
      skips the palindromic read flagging, ``--Reads.palindromicReads.skipFlagging``).
    - Markers: ``kmer_generation_method`` (0 random, 1 excluding globally overenriched k-mers, 2 excluding k-mers
      overenriched in one read, 3 from file ``Kmers.file`` in ``options``, 4 excluding close repeats), ``kmer_length``,
      ``kmer_probability``.
    - Alignment candidates: ``minhash_m``, ``minhash_hash_fraction``, ``minhash_iteration_count`` (MinHash/LowHash search),
      ``minhash_all_pairs`` (skip MinHash, every pair is a candidate).
    - Alignment: ``align_method`` (0 old Shasta, 1 SeqAn, 3 banded SeqAn, 4, 5, 6 experimental; their tuning options
      ``Align.*`` go in ``options``).
    - Read graph: ``read_graph_creation_method`` (0 default criteria, 2 percentile criteria, 3, 4 Bayesian ranking of
      alignments; 4 needs ``strand_separation_method=0``), ``strand_separation_method`` (0 none, 1 approximate
      cross-strand edge flagging, 2 strict), ``max_chimeric_read_distance`` (chimeric read flagging, 0 turns it off),
      ``flag_inconsistent_alignments`` (triangle / least-squares test).
    - Marker graph: ``marker_graph_min_coverage`` / ``marker_graph_max_coverage`` (vertex creation),
      ``marker_graph_min_edge_coverage`` (edge creation), ``cleanup_duplicate_markers``, ``simplify_max_length`` (maximum
      lengths of the simplification iterations, text ``"10,100,1000"`` or a list of integers); approximate transitive
      reduction, strong-subgraph pruning and secondary edges are tuned with ``MarkerGraph.maxDistance``,
      ``MarkerGraph.pruneIterationCount`` and ``MarkerGraph.secondaryEdges.*`` in ``options``.
    - Assembly graph: ``detangle_method`` (0 none, 1 strict, 2 less strict), ``consensus_caller`` (repeat-count consensus:
      ``Modal`` or ``Bayesian:<name>``), ``write_reads_by_segment`` (``ReadsBySegment.csv``), ``store_coverage_data``
      (needs ``memory_mode="filesystem"``, which Shasta runs through sudo).

    ``options`` is a mapping ``{"Section.name": value}`` of any other Shasta 0.14.0 option (``Align.maxSkip``,
    ``Assembly.mode3.minAnchorCoverage``, ...); unknown names are refused. ``memory_mode`` (``anonymous``/``filesystem``) and
    ``memory_backing`` (``4K``/``2M``/``disk``) are Shasta's memory options. ``threads`` is ``--threads`` (Shasta's default
    uses all virtual processors). ``shasta_path`` is the executable, otherwise ``DNAKIT_SHASTA``, otherwise ``shasta`` on PATH.

    Returns ``status`` (``completed`` when the exit code is 0 and an assembled FASTA exists), ``assembly_directory``,
    ``assembly_fasta``, ``assembly_gfa``, ``summary`` (reads used, palindromic and chimeric reads, assembled segments, total
    bases, N50 from ``AssemblySummary.json``), ``outputs`` (every file of the run folder), ``command``, ``version`` and
    ``log_tail``.
    """
    inputs = [reads] if isinstance(reads, str) else list(reads)
    if not inputs:
        raise ConfigurationError("reads must hold at least one file.")
    read_files = [_cli.existing_file(item, "reads") for item in inputs]
    run_dir = Path(output_directory).expanduser().resolve()
    if run_dir.exists():
        raise ConfigurationError(
            f"output_directory must not exist yet (Shasta creates it): {run_dir}"
        )
    if not run_dir.parent.is_dir():
        raise ConfigurationError(
            f"The folder for output_directory does not exist: {run_dir.parent}"
        )
    config_file = Path(config).expanduser()
    if config_file.is_file():
        config_value = str(config_file.resolve())
    elif isinstance(config, str) and _CONFIG_NAME.fullmatch(config):
        config_value = config
    else:
        raise ConfigurationError(
            "config must be a built-in configuration name or an existing configuration file."
        )
    if desired_coverage is not None and (
        isinstance(desired_coverage, bool)
        or not re.fullmatch(r"[0-9]+(\.[0-9]+)?[A-Za-z]{0,2}", str(desired_coverage))
    ):
        raise ConfigurationError(
            "desired_coverage must be a number of bases such as 120000000 or '120Mb'."
        )
    if simplify_max_length is not None and not isinstance(simplify_max_length, str):
        simplify_max_length = ",".join(
            str(_cli.integer(item, "simplify_max_length", minimum=1))
            for item in simplify_max_length
        )
    if simplify_max_length is not None and not re.fullmatch(
        r"[0-9]+(,[0-9]+)*", simplify_max_length
    ):
        raise ConfigurationError("simplify_max_length must look like '10,100,1000'.")
    if consensus_caller is not None:
        _cli.text(consensus_caller, "consensus_caller")

    arguments = [
        "--input",
        *[str(path) for path in read_files],
        "--config",
        config_value,
        "--assemblyDirectory",
        str(run_dir),
    ]
    arguments += ["--command", "assemble"]
    if anchors_path is not None:
        arguments += ["--anchors", str(_cli.existing_file(anchors_path, "anchors_path"))]
    if memory_mode is not None:
        arguments += [
            "--memoryMode",
            _cli.choice(memory_mode, ("anonymous", "filesystem"), "memory_mode"),
        ]
    if memory_backing is not None:
        arguments += [
            "--memoryBacking",
            _cli.choice(memory_backing, ("4K", "2M", "disk"), "memory_backing"),
        ]
    arguments += _opt_int("threads", threads, "--threads", minimum=0, maximum=4096)
    named: list[str] = []
    named += _opt_choice("assembly_mode", assembly_mode, "--Assembly.mode", (0, 2, 3))
    named += _opt_int("min_read_length", min_read_length, "--Reads.minReadLength")
    named += [] if desired_coverage is None else ["--Reads.desiredCoverage", str(desired_coverage)]
    if palindromic_flagging is not None:
        if not isinstance(palindromic_flagging, bool):
            raise ConfigurationError("palindromic_flagging must be True, False or None.")
        named += [] if palindromic_flagging else ["--Reads.palindromicReads.skipFlagging"]
    named += _opt_choice(
        "kmer_generation_method",
        kmer_generation_method,
        "--Kmers.generationMethod",
        (0, 1, 2, 3, 4),
    )
    named += _opt_int("kmer_length", kmer_length, "--Kmers.k", minimum=1, maximum=32)
    named += _opt_float("kmer_probability", kmer_probability, "--Kmers.probability", maximum=1.0)
    named += _opt_int("minhash_m", minhash_m, "--MinHash.m", minimum=1)
    named += _opt_float(
        "minhash_hash_fraction", minhash_hash_fraction, "--MinHash.hashFraction", maximum=1.0
    )
    named += _opt_int(
        "minhash_iteration_count", minhash_iteration_count, "--MinHash.minHashIterationCount"
    )
    named += _opt_bool("minhash_all_pairs", minhash_all_pairs, "--MinHash.allPairs")
    named += _opt_choice("align_method", align_method, "--Align.alignMethod", (0, 1, 3, 4, 5, 6))
    named += _opt_choice(
        "read_graph_creation_method",
        read_graph_creation_method,
        "--ReadGraph.creationMethod",
        (0, 2, 3, 4),
    )
    named += _opt_choice(
        "strand_separation_method",
        strand_separation_method,
        "--ReadGraph.strandSeparationMethod",
        (0, 1, 2),
    )
    named += _opt_int(
        "max_chimeric_read_distance",
        max_chimeric_read_distance,
        "--ReadGraph.maxChimericReadDistance",
    )
    named += _opt_bool(
        "flag_inconsistent_alignments",
        flag_inconsistent_alignments,
        "--ReadGraph.flagInconsistentAlignments",
    )
    named += _opt_int(
        "marker_graph_min_coverage", marker_graph_min_coverage, "--MarkerGraph.minCoverage"
    )
    named += _opt_int(
        "marker_graph_max_coverage", marker_graph_max_coverage, "--MarkerGraph.maxCoverage"
    )
    named += _opt_int(
        "marker_graph_min_edge_coverage",
        marker_graph_min_edge_coverage,
        "--MarkerGraph.minEdgeCoverage",
    )
    named += _opt_bool(
        "cleanup_duplicate_markers",
        cleanup_duplicate_markers,
        "--MarkerGraph.cleanupDuplicateMarkers",
    )
    named += (
        []
        if simplify_max_length is None
        else ["--MarkerGraph.simplifyMaxLength", simplify_max_length]
    )
    named += _opt_choice("detangle_method", detangle_method, "--Assembly.detangleMethod", (0, 1, 2))
    named += [] if consensus_caller is None else ["--Assembly.consensusCaller", consensus_caller]
    named += _opt_bool(
        "write_reads_by_segment", write_reads_by_segment, "--Assembly.writeReadsByAssembledSegment"
    )
    named += _opt_bool("store_coverage_data", store_coverage_data, "--Assembly.storeCoverageData")
    explicit_flags = {item for item in named if item.startswith("--")}
    extra = _checked_options(options, explicit_flags)
    for key in ("Kmers.file", "Kmers.globalFrequencyOverrideDirectory"):
        if options and key in options and not Path(str(options[key])).expanduser().exists():
            raise ConfigurationError(f"{key} must be an existing path.")
    arguments += named + extra

    result = _cli.run(
        _PROGRAM,
        arguments,
        cwd=run_dir.parent,
        expected=[],
        explicit=_resolve(shasta_path),
        hint=_HINT,
        timeout_seconds=timeout_seconds,
    )
    produced = (
        sorted(path for path in run_dir.rglob("*") if path.is_file()) if run_dir.is_dir() else []
    )
    fasta = next(
        (
            run_dir / name
            for name in ("Assembly.fasta", "Assembly-Haploid.fasta")
            if (run_dir / name).is_file()
        ),
        None,
    )
    gfa = run_dir / "Assembly.gfa"
    result["operation"] = "shasta_assemble"
    result["assembly_directory"] = str(run_dir)
    result["assembly_fasta"] = None if fasta is None else str(fasta)
    result["assembly_gfa"] = str(gfa) if gfa.is_file() else None
    result["outputs"] = [str(path) for path in produced]
    result["summary"] = _summary_fields(_read_json(run_dir / "AssemblySummary.json"))
    result["status"] = "completed" if result["return_code"] == 0 and fasta is not None else "failed"
    return result


def _run_directory(value: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(
            f"run_directory must be an existing Shasta assembly folder: {path}"
        )
    return path


def _convert(value: str) -> Any:
    if value == "":
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        number = float(value)
    except ValueError:
        return value
    return number if math.isfinite(number) else value


def shasta_read_run_table(
    run_directory: str,
    table: str,
    *,
    max_rows: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one table Shasta wrote into its run folder, parsed into rows (no Shasta program needed).

    ``table`` selects the file and so the stage whose result is read: ``read_summary`` (palindromic and chimeric flag of
    every read, marker counts, alignment candidates, read graph neighbours), ``read_low_hash_statistics``,
    ``read_length_histogram``, ``low_hash_bucket_histogram`` and ``suppressed_alignment_candidates`` (MinHash candidate
    search), ``aligned_fraction_histogram``, ``alignment_marker_count_histogram``, ``alignment_skip_histogram``,
    ``alignment_drift_histogram``, ``alignment_trim_histogram`` (marker alignment), ``read_graph_components`` (connected
    components of the read graph), ``disjoint_sets_histogram``, ``marker_graph_vertex_coverage_histogram``,
    ``marker_graph_edge_coverage_histogram``, ``anchor_coverage_histogram`` (Mode 3), ``reads_by_segment`` (oriented reads
    that contribute to each assembled segment; the run needs ``write_reads_by_segment=True``) and the Mode 3 tables
    ``assembly_chains``, ``components``, ``oriented_reads_by_component`` and ``journeys``. A table that the run did not write
    raises ``ConfigurationError`` naming the file. At most ``max_rows`` rows are returned; ``row_count`` is the total.

    Returns ``table``, ``file``, ``columns``, ``row_count``, ``rows`` (a list of dictionaries with numbers converted) and
    ``truncated``.
    """
    if table not in TABLES:
        raise ConfigurationError(f"table must be one of {', '.join(sorted(TABLES))}.")
    max_rows = _cli.integer(max_rows, "max_rows", minimum=1, maximum=_MAX_ROWS)
    directory = _run_directory(run_directory)
    name, meaning = TABLES[table]
    path = directory / name
    if not path.is_file():
        raise ConfigurationError(
            f"{name} does not exist in {directory}; this run did not write the {table} table."
        )
    rows: list[dict[str, Any]] = []
    total = 0
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        columns = [column for column in header if column != ""]
        for record in reader:
            if not record:
                continue
            total += 1
            if len(rows) < max_rows:
                rows.append(
                    {
                        column: _convert(record[index]) if index < len(record) else None
                        for index, column in enumerate(columns)
                    }
                )
    return {
        "status": "completed",
        "operation": "shasta_read_run_table",
        "table": table,
        "meaning": meaning,
        "file": str(path),
        "columns": columns,
        "row_count": total,
        "rows": rows,
        "truncated": total > len(rows),
    }


def _n50(lengths: Sequence[int]) -> int:
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    running = 0
    for value in ordered:
        running += value
        if running >= half:
            return value
    return 0


def shasta_read_assembly_output(
    run_directory: str,
    *,
    output: str = "fasta",
    file_name: str | None = None,
    include_sequences: bool = False,
    max_records: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the assembled sequences or the assembly graph that Shasta exported (no Shasta program needed).

    ``output`` chooses the export: ``fasta`` (``Assembly.fasta``, the assembled segments), ``gfa`` (``Assembly.gfa``, the
    assembly graph in GFA 1.0, one strand), ``gfa_both_strands`` (``Assembly-BothStrands.gfa``, both strands) or
    ``gfa_no_sequence`` (``Assembly-BothStrands-NoSequence.gfa``, the graph with lengths only). ``file_name`` selects another
    file of the run folder (for example the Mode 2 files ``Assembly-Haploid.fasta`` or ``Assembly-Phased.gfa``); it must be
    a plain file name. Returns the records (name, length, GC fraction, and the sequence when ``include_sequences``), the
    GFA links, and totals: ``record_count``, ``total_length``, ``longest``, ``n50`` and, for GFA, ``link_count``. At most
    ``max_records`` records are listed.
    """
    kind = _cli.choice(output, tuple(_ASSEMBLY_OUTPUTS), "output")
    max_records = _cli.integer(max_records, "max_records", minimum=1, maximum=_MAX_ROWS)
    directory = _run_directory(run_directory)
    name = _ASSEMBLY_OUTPUTS[kind] if file_name is None else file_name
    if not _FASTA_LIKE.fullmatch(name):
        raise ConfigurationError("file_name must be a plain file name inside the run folder.")
    path = directory / name
    if not path.is_file():
        raise ConfigurationError(
            f"{name} does not exist in {directory}; this run did not write it."
        )
    is_fasta = path.name.endswith((".fasta", ".fa"))
    records: list[dict[str, Any]] = []
    lengths: list[int] = []
    links: list[dict[str, Any]] = []
    link_count = 0

    def add(record_name: str, length: int, gc: int, sequence: str | None) -> None:
        lengths.append(length)
        if len(records) < max_records:
            entry: dict[str, Any] = {
                "name": record_name,
                "length": length,
                "gc_fraction": None if length == 0 or sequence is None else round(gc / length, 6),
            }
            if include_sequences and sequence is not None:
                entry["sequence"] = sequence
            records.append(entry)

    with path.open(encoding="utf-8", errors="replace") as handle:
        if is_fasta:
            current: str | None = None
            parts: list[str] = []

            def flush() -> None:
                if current is not None:
                    sequence = "".join(parts)
                    add(
                        current.split()[0] if current.split() else current,
                        len(sequence),
                        sum(c in "GCgc" for c in sequence),
                        sequence,
                    )

            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    current, parts = line[1:], []
                elif current is not None:
                    parts.append(line.strip())
            flush()
        else:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S" and len(fields) >= 3:
                    sequence = fields[2]
                    tag = next((f[5:] for f in fields[3:] if f.startswith("LN:i:")), None)
                    if sequence == "*":
                        add(fields[1], int(tag) if tag is not None else 0, 0, None)
                    else:
                        add(fields[1], len(sequence), sum(c in "GCgc" for c in sequence), sequence)
                elif fields[0] == "L" and len(fields) >= 6:
                    link_count += 1
                    if len(links) < max_records:
                        links.append(
                            {
                                "from": fields[1],
                                "from_orient": fields[2],
                                "to": fields[3],
                                "to_orient": fields[4],
                                "overlap": fields[5],
                            }
                        )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "shasta_read_assembly_output",
        "output": kind,
        "file": str(path),
        "record_count": len(lengths),
        "total_length": sum(lengths),
        "longest": max(lengths, default=0),
        "n50": _n50(lengths),
        "records": records,
        "truncated": len(lengths) > len(records),
    }
    if not is_fasta:
        result["link_count"] = link_count
        result["links"] = links
    return result


__all__ = ["shasta_assemble", "shasta_read_assembly_output", "shasta_read_run_table"]

_register_public_verbose(globals(), tuple(__all__))
