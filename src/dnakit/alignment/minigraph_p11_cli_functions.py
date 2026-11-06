"""minigraph and its companion scripts as command-line facades (pack P11): fixed argument lists, no shell, bounded runs.

``map_sequences_to_graph``, ``graph_segment_coverage``, ``build_sequence_graph`` and ``call_graph_bubbles`` run the
``minigraph`` program (0.21 is the tested version). Every other function runs one sub-command of ``mgutils.js``,
``mgutils-es6.js`` or ``gafcall.js`` (the scripts in the ``misc/`` folder of the minigraph source tree, executed by the
``k8`` JavaScript shell). The user installs the tools (``conda install -c bioconda minigraph k8``, or builds them from
https://github.com/lh3/minigraph) and either puts them on PATH or passes ``minigraph_path``, ``k8_path`` and the script
path (``mgutils_path``, ``mgutils_es6_path``, ``gafcall_path``). Results go to a fresh ``output_directory`` (never
overwritten); programs that print to standard output are streamed into a file there and standard error goes to
``tool.log``. Results are JSON-serialisable dictionaries with the status, version, command, output files, small parsed
results and a bounded log tail. Optional parameters default to ``None`` and are then not passed, so the tool's own
defaults apply.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.bioinformatics import _p11_cli as _p
from dnakit.exceptions import ConfigurationError

_TESTED_VERSION = "0.21"
_CITATION = "https://doi.org/10.1186/s13059-020-02168-z"
_MINIGRAPH_HINT = (
    "Install minigraph 0.21 or later (for example `conda install -c bioconda minigraph`, or build it from "
    "https://github.com/lh3/minigraph), then put it on PATH or pass minigraph_path."
)
_SCRIPT_HINT = (
    "mgutils.js, mgutils-es6.js and gafcall.js are in the misc/ folder of the minigraph source tree "
    "(https://github.com/lh3/minigraph/tree/master/misc; `conda install -c bioconda minigraph k8` installs them). "
    "Put the script and the k8 JavaScript shell on PATH or pass the script path and k8_path."
)
_PRESETS = ("lr", "asm", "sr", "ggs")


# ---------------------------------------------------------------------------
# Running minigraph and the k8 scripts
# ---------------------------------------------------------------------------


def _pair(
    arguments: list[str],
    option: str,
    value: int | Sequence[int] | None,
    name: str,
    *,
    minimum: int = 0,
) -> None:
    """An option taking ``INT`` or ``INT,INT``."""
    if value is None:
        return
    values = [value] if isinstance(value, int) and not isinstance(value, bool) else list(value)
    if not 1 <= len(values) <= 2:
        raise ConfigurationError(f"{name} must be an integer or a pair of integers.")
    arguments.extend(
        (option, ",".join(str(_cli.integer(item, name, minimum=minimum)) for item in values))
    )


def _mapping_arguments(
    *,
    preset: str | None,
    base_alignment: bool,
    kmer_size: int | None,
    window_size: int | None,
    threads: int,
    min_mapq: int | None = None,
    min_aln_length: int | None = None,
    min_depth_aln_length: int | None = None,
    min_variant_length: int | None = None,
) -> list[str]:
    arguments: list[str] = []
    if preset is not None:
        arguments.extend(("-x", _cli.choice(preset, _PRESETS, "preset")))
    _p.switch(arguments, "-c", base_alignment, "base_alignment")
    _p.opt_int(arguments, "-k", kmer_size, "kmer_size", minimum=1, maximum=28)
    _p.opt_int(arguments, "-w", window_size, "window_size", minimum=1)
    _p.opt_int(arguments, "-q", min_mapq, "min_mapq", minimum=0)
    _p.opt_int(arguments, "-l", min_aln_length, "min_aln_length", minimum=1)
    _p.opt_int(arguments, "-d", min_depth_aln_length, "min_depth_aln_length", minimum=1)
    _p.opt_int(arguments, "-L", min_variant_length, "min_variant_length", minimum=1)
    arguments.extend(("-t", str(_cli.integer(threads, "threads", minimum=1, maximum=256))))
    return arguments


def _run_minigraph(
    operation: str,
    arguments: list[str],
    graph_file: str,
    query_files: Sequence[str],
    output_directory: str | Path,
    output_name: str,
    minigraph_path: str | None,
    timeout_seconds: float,
) -> dict[str, Any]:
    """Run ``minigraph <options> -o <output> <graph> <queries>`` (``-o`` lets minigraph write the file itself)."""
    executable = _cli.program("minigraph", explicit=minigraph_path, hint=_MINIGRAPH_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    output = directory / output_name
    command = [*arguments, "-o", str(output), graph_file, *query_files]
    _p.run_checked(operation, executable, command, directory, timeout_seconds=timeout_seconds)
    _p.require_outputs(operation, directory, {"output": output})
    return _p.result(
        operation,
        "minigraph",
        _p.version_of(executable),
        ["minigraph", *command],
        directory,
        {"output": output, "log": directory / _p.LOG_NAME},
        tool="minigraph",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
    )


def _script(
    operation: str,
    script_name: str,
    script_path: str | None,
    subcommand: str,
    build_arguments: Callable[[Path], list[str]],
    output_directory: str | Path,
    output_name: str,
    *,
    k8_path: str | None,
    timeout_seconds: float,
    with_version: bool = False,
) -> dict[str, Any]:
    """Run ``k8 <script> <subcommand> <arguments>`` with standard output streamed into ``output_name``.

    ``build_arguments`` receives the fresh output directory (for temporary input files that it must delete itself).
    """
    k8 = _cli.program("k8", explicit=k8_path, hint=_SCRIPT_HINT)
    script = _cli.program(script_name, explicit=script_path, hint=_SCRIPT_HINT)
    directory = _cli.new_directory(output_directory, "output_directory")
    output = directory / output_name
    arguments = build_arguments(directory)
    command = [script, subcommand, *arguments]
    _p.run_checked(
        operation, k8, command, directory, stdout_name=output_name, timeout_seconds=timeout_seconds
    )
    version = None
    if with_version:
        probe = directory / ".version"
        if (
            _p.stream(
                k8, [script, "version"], directory, stdout_name=".version", log_name=".version.log"
            )
            == 0
        ):
            version = probe.read_text(encoding="utf-8").strip() or None
        probe.unlink(missing_ok=True)
        (directory / ".version.log").unlink(missing_ok=True)
    return _p.result(
        operation,
        script_name,
        version,
        ["k8", script_name, subcommand, *arguments],
        directory,
        {"output": output, "log": directory / _p.LOG_NAME},
        tool="minigraph",
        tested_version=_TESTED_VERSION,
        citation=_CITATION,
        k8_version=_k8_version(k8, directory),
    )


def _k8_version(k8: str, directory: Path) -> str | None:
    """The k8 shell version (``k8 -v`` prints ``k8: 1.2-r137`` after the V8 line)."""
    probe = directory / ".k8_version"
    try:
        if (
            _p.stream(
                k8,
                ["-v"],
                directory,
                stdout_name=".k8_version",
                log_name=".k8_version.log",
                timeout_seconds=60.0,
            )
            == 0
        ):
            for line in probe.read_text(encoding="utf-8").splitlines():
                if line.startswith("k8:"):
                    return line.split(":", 1)[1].strip()
    finally:
        probe.unlink(missing_ok=True)
        (directory / ".k8_version.log").unlink(missing_ok=True)
    return None


def _output_name(value: str) -> str:
    return _p.file_stem(value, "output_name")


def _data_rows(result: dict[str, Any]) -> tuple[list[list[str]], int]:
    return _p.table(Path(result["outputs"]["output"]))


# ---------------------------------------------------------------------------
# minigraph itself
# ---------------------------------------------------------------------------


def map_sequences_to_graph(
    graph: str | Path,
    queries: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "mapping.gaf",
    preset: str | None = None,
    base_alignment: bool = False,
    kmer_size: int | None = None,
    window_size: int | None = None,
    repetitive_fraction: float | None = None,
    occurrence_range: int | Sequence[int] | None = None,
    divergence: float | None = None,
    chain_gap: int | None = None,
    bandwidth: int | Sequence[int] | None = None,
    min_minimizers: int | Sequence[int] | None = None,
    min_chain_score: int | Sequence[int] | None = None,
    min_secondary_ratio: float | None = None,
    max_secondary: int | None = None,
    skip_self_diagonal: bool = False,
    linear_chains: bool = False,
    vertex_coordinates: bool = False,
    threads: int = 4,
    minigraph_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map DNA sequences (reads or assemblies) to a sequence graph, or to a linear reference (minigraph).

    ``graph`` is a GFA file (or a FASTA file for a linear target); ``queries`` are FASTA/FASTQ files, optionally gzipped.
    ``preset`` is ``lr`` (noisy long reads, the default), ``asm`` (assembly-to-graph), ``sr`` (short reads) or ``ggs``.
    ``base_alignment`` is minigraph's ``-c`` (base-level alignment; recommended by the authors). The output is GAF
    (PAF when ``graph`` is linear) written to ``<output_directory>/<output_name>``; with ``vertex_coordinates`` the path
    coordinates are vertex coordinates. Result keys: ``records`` (number of alignment lines), ``mapped_queries`` (number
    of distinct query names), ``outputs``.
    """
    graph_file = _p.file_path(graph, "graph")
    query_files = _p.file_paths(queries, "queries")
    name = _output_name(output_name)
    arguments = _mapping_arguments(
        preset=preset,
        base_alignment=base_alignment,
        kmer_size=kmer_size,
        window_size=window_size,
        threads=threads,
    )
    _p.opt_num(
        arguments, "-f", repetitive_fraction, "repetitive_fraction", minimum=0.0, maximum=1.0
    )
    _pair(arguments, "-U", occurrence_range, "occurrence_range", minimum=1)
    _p.opt_num(arguments, "-j", divergence, "divergence", minimum=0.0, maximum=1.0)
    _p.opt_int(arguments, "-g", chain_gap, "chain_gap", minimum=1)
    _pair(arguments, "-r", bandwidth, "bandwidth", minimum=1)
    _pair(arguments, "-n", min_minimizers, "min_minimizers", minimum=1)
    _pair(arguments, "-m", min_chain_score, "min_chain_score", minimum=1)
    _p.opt_num(
        arguments, "-p", min_secondary_ratio, "min_secondary_ratio", minimum=0.0, maximum=1.0
    )
    _p.opt_int(arguments, "-N", max_secondary, "max_secondary", minimum=0)
    _p.switch(arguments, "-D", skip_self_diagonal, "skip_self_diagonal")
    _p.switch(arguments, "-S", linear_chains, "linear_chains")
    _p.switch(arguments, "--vc", vertex_coordinates, "vertex_coordinates")
    result = _run_minigraph(
        "map_sequences_to_graph",
        arguments,
        graph_file,
        query_files,
        output_directory,
        name,
        minigraph_path,
        timeout_seconds,
    )
    rows, total = _p.table(Path(result["outputs"]["output"]), limit=10_000_000)
    result["records"] = total
    result["mapped_queries"] = len({row[0] for row in rows})
    return result


def graph_segment_coverage(
    graph: str | Path,
    reads: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "coverage.gfa",
    preset: str | None = None,
    base_alignment: bool = True,
    min_aligned_length: int | None = None,
    min_mapq: int | None = None,
    kmer_size: int | None = None,
    window_size: int | None = None,
    threads: int = 4,
    minigraph_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remap reads (or assemblies) to a graph and report the coverage of every graph segment (minigraph ``--cov``).

    The output is the input graph in GFA format with a ``dc:f:`` (mean depth) tag on every covered segment.
    ``min_aligned_length`` and ``min_mapq`` are ``--min-cov-blen`` and ``--min-cov-mapq`` (alignments below them do not
    count). Result keys: ``segments``, ``covered_segments`` (segments carrying the depth tag), ``mean_depths`` (segment
    name to depth, at most 1000 segments).
    """
    graph_file = _p.file_path(graph, "graph")
    read_files = _p.file_paths(reads, "reads")
    name = _output_name(output_name)
    arguments = _mapping_arguments(
        preset=preset,
        base_alignment=base_alignment,
        kmer_size=kmer_size,
        window_size=window_size,
        threads=threads,
    )
    arguments.append("--cov")
    _p.opt_int(arguments, "--min-cov-blen", min_aligned_length, "min_aligned_length", minimum=0)
    _p.opt_int(arguments, "--min-cov-mapq", min_mapq, "min_mapq", minimum=0)
    result = _run_minigraph(
        "graph_segment_coverage",
        arguments,
        graph_file,
        read_files,
        output_directory,
        name,
        minigraph_path,
        timeout_seconds,
    )
    depths: dict[str, float] = {}
    segments = 0
    with Path(result["outputs"]["output"]).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("S\t"):
                continue
            segments += 1
            fields = line.rstrip("\n").split("\t")
            for tag in fields[3:]:
                if tag.startswith("dc:f:") and len(depths) < _p.ROW_LIMIT:
                    depths[fields[1]] = float(tag[5:])
    result["segments"] = segments
    result["covered_segments"] = len(depths) if len(depths) < _p.ROW_LIMIT else None
    result["mean_depths"] = depths
    return result


def build_sequence_graph(
    reference: str | Path,
    assemblies: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "graph.gfa",
    preset: str | None = "ggs",
    base_alignment: bool = True,
    min_mapq: int | None = None,
    min_aln_length: int | None = None,
    min_depth_aln_length: int | None = None,
    min_variant_length: int | None = None,
    kmer_size: int | None = None,
    window_size: int | None = None,
    threads: int = 4,
    minigraph_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a sequence graph by incrementally adding assemblies to a reference (minigraph ``-cxggs``).

    ``reference`` is a FASTA file or an existing GFA graph; the ``assemblies`` are added in the order given. ``preset``
    ``ggs`` (the default) is ``-x ggs``; ``None`` runs the bare ``--ggen`` mode without the preset parameters.
    ``min_mapq`` (``-q``), ``min_aln_length`` (``-l``), ``min_depth_aln_length`` (``-d``) and ``min_variant_length``
    (``-L``) keep minigraph's defaults (5, 100000, 20000 and 50) unless set. The output is GFA. Result keys: ``segments``,
    ``links``, ``paths``.
    """
    reference_file = _p.file_path(reference, "reference")
    assembly_files = _p.file_paths(assemblies, "assemblies")
    name = _output_name(output_name)
    arguments = _mapping_arguments(
        preset=preset,
        base_alignment=base_alignment,
        kmer_size=kmer_size,
        window_size=window_size,
        threads=threads,
        min_mapq=min_mapq,
        min_aln_length=min_aln_length,
        min_depth_aln_length=min_depth_aln_length,
        min_variant_length=min_variant_length,
    )
    if preset is None:
        arguments.append("--ggen")
    result = _run_minigraph(
        "build_sequence_graph",
        arguments,
        reference_file,
        assembly_files,
        output_directory,
        name,
        minigraph_path,
        timeout_seconds,
    )
    output = Path(result["outputs"]["output"])
    result["segments"] = _p.count_lines(output, "S\t")
    result["links"] = _p.count_lines(output, "L\t")
    result["paths"] = _p.count_lines(output, "W\t") + _p.count_lines(output, "P\t")
    return result


def call_graph_bubbles(
    graph: str | Path,
    assembly: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "calls.bed",
    preset: str | None = "asm",
    base_alignment: bool = True,
    min_mapq: int | None = None,
    min_aln_length: int | None = None,
    min_depth_aln_length: int | None = None,
    min_variant_length: int | None = None,
    kmer_size: int | None = None,
    window_size: int | None = None,
    threads: int = 4,
    minigraph_path: str | None = None,
    timeout_seconds: float = _cli.DEFAULT_TIMEOUT_SECONDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Call the path an assembly takes through every bubble of a graph, one BED row per bubble (minigraph ``--call``).

    This is the per-sample step of the minigraph pangenome pipeline (``minigraph -cxasm --call graph.gfa sample.fa``).
    The tool's own defaults for the alignment length thresholds (``-l`` 100000, ``-d`` 20000) are kept unless set; the
    minigraph documentation uses ``-l10000 -d10000`` for assemblies with short contigs. Result keys: ``bubbles`` (rows),
    ``rows`` (first 1000 rows split into columns).
    """
    graph_file = _p.file_path(graph, "graph")
    assembly_file = _p.file_path(assembly, "assembly")
    name = _output_name(output_name)
    arguments = _mapping_arguments(
        preset=preset,
        base_alignment=base_alignment,
        kmer_size=kmer_size,
        window_size=window_size,
        threads=threads,
        min_mapq=min_mapq,
        min_aln_length=min_aln_length,
        min_depth_aln_length=min_depth_aln_length,
        min_variant_length=min_variant_length,
    )
    arguments.append("--call")
    result = _run_minigraph(
        "call_graph_bubbles",
        arguments,
        graph_file,
        [assembly_file],
        output_directory,
        name,
        minigraph_path,
        timeout_seconds,
    )
    rows, total = _data_rows(result)
    result["bubbles"] = total
    result["rows"] = rows
    return result


# ---------------------------------------------------------------------------
# mgutils.js
# ---------------------------------------------------------------------------


def rename_fasta_sequences(
    fasta: str | Path,
    prefix: str,
    output_directory: str | Path,
    *,
    output_name: str = "renamed.fa",
    delimiter: str = "#",
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rename FASTA sequences to ``<prefix><delimiter><name>`` (the PanSN convention; ``mgutils.js renamefa``).

    An existing ``<something>#`` prefix of the sequence name is removed first. ``delimiter`` is mgutils' ``-d``
    (default ``#``). Result key: ``sequences`` (number of FASTA records written).
    """
    source = _p.file_path(fasta, "fasta")
    name = _output_name(output_name)
    prefix_text = _cli.text(prefix, "prefix")
    separator = _cli.text(delimiter, "delimiter")
    result = _script(
        "rename_fasta_sequences",
        "mgutils.js",
        mgutils_path,
        "renamefa",
        lambda directory: ["-d", separator, prefix_text, source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["sequences"] = _p.count_lines(Path(result["outputs"]["output"]), ">")
    return result


def blacklist_insert_alignments(
    paf: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "blacklist.bed",
    min_divergence: float = 0.01,
    identical_only: bool = False,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Blacklist reference regions that inserted sequences align to (``mgutils.js paf2bl``).

    ``paf`` holds alignments of the inserted sequences to the reference (with the ``tp`` and ``de`` tags written by
    minimap2). Alignments with mapping quality below 5 or length below 500 are ignored. By default regions with a
    divergence between ``min_divergence`` (``-d``, default 0.01) and 0.1 are written; ``identical_only`` (``-s``)
    writes regions with divergence up to 0.002 instead. The output is BED (reference, start, end). Result key:
    ``regions``.
    """
    source = _p.file_path(paf, "paf")
    name = _output_name(output_name)
    divergence = format(
        _cli.number(min_divergence, "min_divergence", minimum=0.0, maximum=1.0), ".12g"
    )
    only_identical = _p.flag(identical_only, "identical_only")

    def arguments(directory: Path) -> list[str]:
        values = ["-d", divergence]
        if only_identical:
            values.append("-s")
        return [*values, source]

    result = _script(
        "blacklist_insert_alignments",
        "mgutils.js",
        mgutils_path,
        "paf2bl",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["regions"] = _data_rows(result)[1]
    return result


def annotate_inserted_sequences(
    bed: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "anno.tsv",
    repeatmasker_out: str | Path | None = None,
    gaps: str | Path | None = None,
    low_complexity: str | Path | None = None,
    tandem_repeats: str | Path | None = None,
    paf: str | Path | None = None,
    centromere: str | Path | None = None,
    bubbles: str | Path | None = None,
    segdup: str | Path | None = None,
    min_feature_length: int | None = None,
    min_segdup_fraction: float | None = None,
    annotate_mobile: bool = False,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotate short inserted sequences (repeats, low-complexity regions, gaps, segmental duplications; ``mgutils.js anno``).

    ``bed`` lists the inserted sequences as ``contig start end length``. Evidence files (all optional): the RepeatMasker
    ``.out`` file (``-r``), the stretches of N reported by ``seqtk cutN``/``seqtk gap`` (``-g``), minimap2/sdust output
    for low-complexity regions (``-d``), eTRF output (``-e``), a PAF alignment against the reference (``-p``), dna-brnn
    centromere calls (``-c``), the ``gfatools bubble`` file (``-b``) and the segmental-duplication table (``-s``).
    ``min_feature_length`` is ``-l`` (default 30), ``min_segdup_fraction`` is ``-S`` (default 0.2) and
    ``annotate_mobile`` (``-m``) reports AluY and L1HS separately. Result key: ``annotated`` (rows written).
    """
    source = _p.file_path(bed, "bed")
    name = _output_name(output_name)
    arguments: list[str] = []
    _p.opt_file(arguments, "-r", repeatmasker_out, "repeatmasker_out")
    _p.opt_file(arguments, "-g", gaps, "gaps")
    _p.opt_file(arguments, "-d", low_complexity, "low_complexity")
    _p.opt_file(arguments, "-e", tandem_repeats, "tandem_repeats")
    _p.opt_file(arguments, "-p", paf, "paf")
    _p.opt_file(arguments, "-c", centromere, "centromere")
    _p.opt_file(arguments, "-b", bubbles, "bubbles")
    _p.opt_file(arguments, "-s", segdup, "segdup")
    _p.opt_int(arguments, "-l", min_feature_length, "min_feature_length", minimum=0)
    _p.opt_num(
        arguments, "-S", min_segdup_fraction, "min_segdup_fraction", minimum=0.0, maximum=1.0
    )
    _p.switch(arguments, "-m", annotate_mobile, "annotate_mobile")
    result = _script(
        "annotate_inserted_sequences",
        "mgutils.js",
        mgutils_path,
        "anno",
        lambda directory: [*arguments, source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["annotated"] = _data_rows(result)[1]
    return result


def summarize_insert_annotation(
    annotation: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "anno_table.tsv",
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the output of ``annotate_inserted_sequences`` as a table of repeat classes (``mgutils.js anno2tbl``).

    Gap rows and unplaced contigs (``chrUn``, ``_random``) are ignored. Each output row is the class key, its label and
    nine counts (inserts with 0, 1, 2 and 3-or-more annotated features, the summed insert length for each, and the summed
    annotated length). Result key: ``classes`` (rows, first 1000).
    """
    source = _p.file_path(annotation, "annotation")
    name = _output_name(output_name)
    result = _script(
        "summarize_insert_annotation",
        "mgutils.js",
        mgutils_path,
        "anno2tbl",
        lambda directory: [source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    rows, total = _data_rows(result)
    result["classes"] = [
        {"key": row[0], "label": row[1], "counts": [_p.number(item) for item in row[2:]]}
        for row in rows
    ]
    result["class_count"] = total
    return result


def convert_to_stable_gaf(
    graph: str | Path,
    alignments: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "stable.gaf",
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert graph alignments (GAF) to stable coordinates on the linear sequences of the graph (``mgutils.js stableGaf``).

    ``graph`` is the minigraph GFA with the ``SN``, ``SO`` and ``SR`` segment tags; ``alignments`` is the GAF written by
    minigraph. Alignments inside one reference-rank-0 segment path are rewritten in the coordinates of that sequence.
    Result key: ``records``.
    """
    graph_file = _p.file_path(graph, "graph")
    gaf = _p.file_path(alignments, "alignments")
    name = _output_name(output_name)
    result = _script(
        "convert_to_stable_gaf",
        "mgutils.js",
        mgutils_path,
        "stableGaf",
        lambda directory: [graph_file, gaf],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["records"] = _data_rows(result)[1]
    return result


def extract_gaf_segment_path(
    first_segment: str,
    second_segment: str,
    alignments: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "segment_path.bed",
    min_alignment_length: int = 100000,
    from_end: bool = False,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract the contig interval between two graph segments from GAF alignments (``mgutils.js extractseg``).

    ``first_segment`` and ``second_segment`` are segment names (for example ``s12``). For every alignment file the
    interval of each contig that lies between the two segments is written as BED (contig, start, end, ``*``, 0, strand).
    Alignments shorter than ``min_alignment_length`` (``-l``, default 100000) are ignored; ``from_end`` (``-e``) takes the
    interval from the end of the first segment to the start of the second. Result key: ``intervals``.
    """
    first = _cli.text(first_segment, "first_segment")
    second = _cli.text(second_segment, "second_segment")
    files = _p.file_paths(alignments, "alignments")
    name = _output_name(output_name)
    minimum = str(_cli.integer(min_alignment_length, "min_alignment_length", minimum=0))
    anchored = _p.flag(from_end, "from_end")

    def arguments(directory: Path) -> list[str]:
        values = ["-l", minimum]
        if anchored:
            values.append("-e")
        return [*values, first, second, *files]

    result = _script(
        "extract_gaf_segment_path",
        "mgutils.js",
        mgutils_path,
        "extractseg",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["intervals"] = _data_rows(result)[1]
    return result


def merge_bubble_calls(
    call_files: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "merged.txt",
    samples: Sequence[str] | None = None,
    annotation: str | Path | None = None,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge the per-sample ``--call`` BED files into one multi-sample table (``paste *.bed | mgutils.js merge``).

    ``call_files`` are the BED files written by ``call_graph_bubbles`` for the same graph (same bubbles in the same
    order), one per sample. ``samples`` are the sample names for the header (``-s``); ``annotation`` is a bubble
    annotation file whose twelfth column is added as ``ANNO`` (``-a``). The output has a VCF-like header, one row per
    bubble and one genotype column per sample. Result keys: ``bubbles`` (rows), ``sample_count``.
    """
    files = _p.file_paths(call_files, "call_files")
    name = _output_name(output_name)
    names = None if samples is None else [_cli.text(item, "samples") for item in samples]
    if names is not None and len(names) != len(files):
        raise ConfigurationError("samples must have one name per call file.")
    annotation_file = (
        None if annotation is None else str(_cli.existing_file(annotation, "annotation"))
    )
    columns = [Path(item).read_text(encoding="utf-8").splitlines() for item in files]
    if len({len(column) for column in columns}) != 1:
        raise ConfigurationError(
            "All call files must have the same number of rows (same graph and bubbles)."
        )

    def arguments(directory: Path) -> list[str]:
        pasted = directory / "pasted.bed"
        pasted.write_text(
            "".join("\t".join(parts) + "\n" for parts in zip(*columns, strict=True)),
            encoding="utf-8",
        )
        values: list[str] = []
        if annotation_file is not None:
            values.extend(("-a", annotation_file))
        if names is not None:
            listing = directory / "samples.txt"
            listing.write_text("".join(item + "\n" for item in names), encoding="utf-8")
            values.extend(("-s", str(listing)))
        return [*values, str(pasted)]

    result = _script(
        "merge_bubble_calls",
        "mgutils.js",
        mgutils_path,
        "merge",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    directory = Path(result["output_directory"])
    (directory / "pasted.bed").unlink(missing_ok=True)
    (directory / "samples.txt").unlink(missing_ok=True)
    result["bubbles"] = _data_rows(result)[1]
    result["sample_count"] = len(files)
    return result


def segment_allele_frequency(
    segment_bed: str | Path,
    merged_calls: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "segfreq.tsv",
    bubble_bed: str | Path | None = None,
    min_frequency: float | None = None,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the allele frequency of every graph segment from merged bubble calls (``mgutils.js segfreq``).

    ``segment_bed`` is the ``gfatools gfa2bed`` table of the graph, ``merged_calls`` the output of
    ``merge_bubble_calls`` and ``bubble_bed`` an optional ``gfatools bubble`` file. ``min_frequency`` is ``-f``
    (default 0.05). Result key: ``segments`` (rows).
    """
    segments = _p.file_path(segment_bed, "segment_bed")
    merged = _p.file_path(merged_calls, "merged_calls")
    bubbles = None if bubble_bed is None else _p.file_path(bubble_bed, "bubble_bed")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_num(options, "-f", min_frequency, "min_frequency", minimum=0.0, maximum=1.0)
    result = _script(
        "segment_allele_frequency",
        "mgutils.js",
        mgutils_path,
        "segfreq",
        lambda directory: [*options, segments, merged, *([] if bubbles is None else [bubbles])],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    result["segments"] = _data_rows(result)[1]
    return result


def count_gene_copies(
    alignments: str | Path,
    gene_bed: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "gene_copy.tsv",
    min_coverage: float | None = None,
    min_relative_coverage: float | None = None,
    mgutils_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the copies of each gene in an assembly from graph alignments of the gene sequences (``mgutils.js genecopy``).

    ``alignments`` is the GAF (with ``cg`` tags) of gene sequences mapped to the assembly graph, ``gene_bed`` the BED
    (name in column 4, strand in column 6) of the genes on the source reference. ``min_coverage`` is ``-c`` (default 0.8)
    and ``min_relative_coverage`` is ``-r`` (default 0.85). Result key: ``rows`` (output lines, first 1000).
    """
    gaf = _p.file_path(alignments, "alignments")
    genes = _p.file_path(gene_bed, "gene_bed")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_num(options, "-c", min_coverage, "min_coverage", minimum=0.0, maximum=1.0)
    _p.opt_num(
        options, "-r", min_relative_coverage, "min_relative_coverage", minimum=0.0, maximum=1.0
    )
    result = _script(
        "count_gene_copies",
        "mgutils.js",
        mgutils_path,
        "genecopy",
        lambda directory: [*options, gaf, genes],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    rows, total = _data_rows(result)
    result["rows"] = rows
    result["row_count"] = total
    return result


def convert_calls_to_vcf(
    merged_calls: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "calls.vcf",
    reference_index: int | None = None,
    max_alleles: int | None = None,
    samples: Sequence[str] | None = None,
    mgutils_es6_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert the merged bubble calls to VCF (``mgutils-es6.js merge2vcf``).

    ``reference_index`` is ``-r`` (which sample column is the reference, default 0), ``max_alleles`` is ``-a`` (default
    15) and ``samples`` the sample names (``-s``; one per call column). Result key: ``variants`` (VCF records).
    """
    merged = _p.file_path(merged_calls, "merged_calls")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_int(options, "-r", reference_index, "reference_index", minimum=0)
    _p.opt_int(options, "-a", max_alleles, "max_alleles", minimum=2)
    names = None if samples is None else [_cli.text(item, "samples") for item in samples]

    def arguments(directory: Path) -> list[str]:
        values = list(options)
        if names is not None:
            listing = directory / "samples.txt"
            listing.write_text("".join(item + "\n" for item in names), encoding="utf-8")
            values.extend(("-s", str(listing)))
        return [*values, merged]

    result = _script(
        "convert_calls_to_vcf",
        "mgutils-es6.js",
        mgutils_es6_path,
        "merge2vcf",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    (Path(result["output_directory"]) / "samples.txt").unlink(missing_ok=True)
    result["variants"] = _data_rows(result)[1]
    return result


def add_sample_names_to_calls(
    merged_calls: str | Path,
    samples: Sequence[str],
    output_directory: str | Path,
    *,
    output_name: str = "merged_named.txt",
    mgutils_es6_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Put sample names into the header of a merged bubble-call table (``mgutils-es6.js addsample``).

    Replaces the ``#CHROM`` header line of ``merged_calls`` with one that lists ``samples``. Result key:
    ``sample_count``.
    """
    merged = _p.file_path(merged_calls, "merged_calls")
    names = [_cli.text(item, "samples") for item in samples] if not isinstance(samples, str) else []
    if not names:
        raise ConfigurationError("samples must be a non-empty list of names.")
    name = _output_name(output_name)

    def arguments(directory: Path) -> list[str]:
        listing = directory / "samples.txt"
        listing.write_text("".join(item + "\n" for item in names), encoding="utf-8")
        return [merged, str(listing)]

    result = _script(
        "add_sample_names_to_calls",
        "mgutils-es6.js",
        mgutils_es6_path,
        "addsample",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
    )
    (Path(result["output_directory"]) / "samples.txt").unlink(missing_ok=True)
    result["sample_count"] = len(names)
    return result


# ---------------------------------------------------------------------------
# gafcall.js
# ---------------------------------------------------------------------------


def extract_long_indels(
    stable_gaf: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "extract.tsv",
    sample_name: str | None = None,
    min_indel_length: int | None = None,
    max_indels_per_10kb: int | None = None,
    min_mapq: int | None = None,
    min_mapq_ends: int | None = None,
    min_alignment_length_ends: int | None = None,
    min_alignment_length_middle: int | None = None,
    polya_penalty: int | None = None,
    centromere_bed: str | Path | None = None,
    gafcall_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract long INDELs and breakpoints from stable-coordinate GAF alignments of reads (``gafcall.js extract``).

    ``stable_gaf`` comes from ``convert_to_stable_gaf``. Optional parameters (gafcall defaults in brackets):
    ``sample_name`` (``-n``, foo), ``min_indel_length`` (``-l``, 100), ``max_indels_per_10kb`` (``-c``, 3),
    ``min_mapq`` (``-q``, 5), ``min_mapq_ends`` (``-Q``, 30), ``min_alignment_length_ends`` (``-e``, 2000),
    ``min_alignment_length_middle`` (``-m``, 50), ``polya_penalty`` (``-a``, 5) and a centromere BED (``-b``). Result key:
    ``records`` (evidence lines).
    """
    source = _p.file_path(stable_gaf, "stable_gaf")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_text(options, "-n", sample_name, "sample_name")
    _p.opt_int(options, "-l", min_indel_length, "min_indel_length", minimum=1)
    _p.opt_int(options, "-c", max_indels_per_10kb, "max_indels_per_10kb", minimum=0)
    _p.opt_int(options, "-q", min_mapq, "min_mapq", minimum=0)
    _p.opt_int(options, "-Q", min_mapq_ends, "min_mapq_ends", minimum=0)
    _p.opt_int(options, "-e", min_alignment_length_ends, "min_alignment_length_ends", minimum=0)
    _p.opt_int(options, "-m", min_alignment_length_middle, "min_alignment_length_middle", minimum=0)
    _p.opt_int(options, "-a", polya_penalty, "polya_penalty", minimum=0)
    _p.opt_file(options, "-b", centromere_bed, "centromere_bed")
    result = _script(
        "extract_long_indels",
        "gafcall.js",
        gafcall_path,
        "extract",
        lambda directory: [*options, source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        with_version=True,
    )
    result["records"] = _data_rows(result)[1]
    return result


def merge_sv_evidence(
    extract_files: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "calls.tsv",
    min_read_count: int | None = None,
    min_read_count_per_strand: int | None = None,
    window_size: int | None = None,
    max_length_difference: float | None = None,
    min_centromere_distance: int | None = None,
    min_retrotransposon_length: int | None = None,
    min_retrotransposon_reads: int | None = None,
    max_alleles: int | None = None,
    max_reads_compared: int | None = None,
    gafcall_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge read-level INDEL and breakpoint evidence into structural-variant calls (``gafcall.js merge``).

    ``extract_files`` are outputs of ``extract_long_indels`` (several files, for example of different read groups of the
    same sample, are concatenated and sorted by contig and position first, as ``sort -k1,1 -k2,2n`` does). Optional
    parameters (gafcall defaults in brackets): ``min_read_count`` (``-c``, 4), ``min_read_count_per_strand`` (``-s``, 2),
    ``window_size`` (``-w``, 100), ``max_length_difference`` (``-d``, 0.05), ``min_centromere_distance`` (``-e``, 500000),
    ``min_retrotransposon_length`` (``-r``, 0 disables), ``min_retrotransposon_reads`` (``-R``, 1), ``max_alleles``
    (``-A``, 100) and ``max_reads_compared`` (``-C``, 500). Result key: ``calls`` (rows).
    """
    files = _p.file_paths(extract_files, "extract_files")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_int(options, "-c", min_read_count, "min_read_count", minimum=1)
    _p.opt_int(options, "-s", min_read_count_per_strand, "min_read_count_per_strand", minimum=0)
    _p.opt_int(options, "-w", window_size, "window_size", minimum=1)
    _p.opt_num(
        options, "-d", max_length_difference, "max_length_difference", minimum=0.0, maximum=1.0
    )
    _p.opt_int(options, "-e", min_centromere_distance, "min_centromere_distance", minimum=0)
    _p.opt_int(options, "-r", min_retrotransposon_length, "min_retrotransposon_length", minimum=0)
    _p.opt_int(options, "-R", min_retrotransposon_reads, "min_retrotransposon_reads", minimum=0)
    _p.opt_int(options, "-A", max_alleles, "max_alleles", minimum=1)
    _p.opt_int(options, "-C", max_reads_compared, "max_reads_compared", minimum=1)
    lines = [
        line
        for item in files
        for line in Path(item).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    def sort_key(line: str) -> tuple[str, int]:
        fields = line.split("\t")
        try:
            return fields[0], int(fields[1])
        except (IndexError, ValueError):
            return fields[0], 0

    lines.sort(key=sort_key)

    def arguments(directory: Path) -> list[str]:
        merged_input = directory / "sorted_extract.tsv"
        merged_input.write_text("".join(line + "\n" for line in lines), encoding="utf-8")
        return [*options, str(merged_input)]

    result = _script(
        "merge_sv_evidence",
        "gafcall.js",
        gafcall_path,
        "merge",
        arguments,
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        with_version=True,
    )
    (Path(result["output_directory"]) / "sorted_extract.tsv").unlink(missing_ok=True)
    result["calls"] = _data_rows(result)[1]
    return result


def evaluate_sv_calls(
    vcf_files: Sequence[str | Path],
    output_directory: str | Path,
    *,
    output_name: str = "eval.txt",
    confident_regions: str | Path | None = None,
    min_sv_length: int | None = None,
    window_size: int | None = None,
    read_length_ratio: float | None = None,
    min_length_ratio: float | None = None,
    min_allele_frequency: float | None = None,
    ignore_filter: bool = False,
    check_genotype: bool = False,
    print_errors: bool = False,
    gafcall_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate structural-variant calls against a truth set, or compare several call sets (``gafcall.js eval``).

    With two VCF files, the first is the truth ("base") set and the second the call ("test") set; the result has the
    recall (``recall``: truth SVs missed and found) and precision-side counts (``precision``). With more than two files
    the result is the pairwise concordance matrix (``concordance``). Optional parameters (gafcall defaults in
    brackets): ``confident_regions`` BED (``-b``), ``min_sv_length`` (``-l``, 100), ``window_size`` (``-w``, 500),
    ``read_length_ratio`` (``-r``, 0.8), ``min_length_ratio`` (``-m``, 0.6), ``min_allele_frequency`` (``-v``, 0),
    ``ignore_filter`` (``-F``), ``check_genotype`` (``-G``) and ``print_errors`` (``-e``). The raw report is in
    ``outputs['output']``.
    """
    files = _p.file_paths(vcf_files, "vcf_files")
    if len(files) < 2:
        raise ConfigurationError("vcf_files needs at least two files: truth set then call set.")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_file(options, "-b", confident_regions, "confident_regions")
    _p.opt_int(options, "-l", min_sv_length, "min_sv_length", minimum=1)
    _p.opt_int(options, "-w", window_size, "window_size", minimum=0)
    _p.opt_num(options, "-r", read_length_ratio, "read_length_ratio", minimum=0.0, maximum=1.0)
    _p.opt_num(options, "-m", min_length_ratio, "min_length_ratio", minimum=0.0, maximum=1.0)
    _p.opt_num(
        options, "-v", min_allele_frequency, "min_allele_frequency", minimum=0.0, maximum=1.0
    )
    _p.switch(options, "-F", ignore_filter, "ignore_filter")
    _p.switch(options, "-G", check_genotype, "check_genotype")
    _p.switch(options, "-e", print_errors, "print_errors")
    result = _script(
        "evaluate_sv_calls",
        "gafcall.js",
        gafcall_path,
        "eval",
        lambda directory: [*options, *files],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        with_version=True,
    )
    lines = Path(result["outputs"]["output"]).read_text(encoding="utf-8").splitlines()
    if len(files) == 2:
        summary: dict[str, Any] = {}
        for line in lines:
            fields = line.split("\t")
            if fields and fields[0] in ("RN", "RP") and len(fields) >= 4:
                summary["recall" if fields[0] == "RN" else "precision"] = {
                    "total": _p.number(fields[1]),
                    "errors": _p.number(fields[2]),
                    "error_rate": _p.number(fields[3]),
                }
        result.update(summary)
    else:
        result["concordance"] = [line.split("\t") for line in lines if line.startswith("SN")]
    return result


def view_sv_calls(
    vcf: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "sv.tsv",
    min_length: int | None = None,
    regions: str | Path | None = None,
    ignore_filter: bool = False,
    check_genotype: bool = False,
    count_long: bool = False,
    gafcall_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Re-format the structural variants of a VCF as a table with their long-INDEL counts (``gafcall.js view``).

    One row per SV: contig, position, orientation, second contig, second position, SV type and SV length. With
    ``count_long`` (``-C``) the output is instead the number of SVs of at least 20 kb, 100 kb and 1 Mb and the number of
    translocations. ``min_length`` is ``-l`` (default 100), ``regions`` a BED of regions to keep (``-b``), and
    ``ignore_filter`` / ``check_genotype`` are ``-F`` / ``-G``. Result key: ``records`` (output lines).
    """
    source = _p.file_path(vcf, "vcf")
    name = _output_name(output_name)
    options: list[str] = []
    _p.opt_int(options, "-l", min_length, "min_length", minimum=1)
    _p.opt_file(options, "-b", regions, "regions")
    _p.switch(options, "-F", ignore_filter, "ignore_filter")
    _p.switch(options, "-G", check_genotype, "check_genotype")
    _p.switch(options, "-C", count_long, "count_long")
    result = _script(
        "view_sv_calls",
        "gafcall.js",
        gafcall_path,
        "view",
        lambda directory: [*options, source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        with_version=True,
    )
    rows, total = _data_rows(result)
    result["records"] = total
    result["rows"] = rows
    return result


def filter_sv_extract_by_type(
    filter_extract: str | Path,
    extract: str | Path,
    output_directory: str | Path,
    *,
    output_name: str = "joined.tsv",
    gafcall_path: str | None = None,
    k8_path: str | None = None,
    timeout_seconds: float = 3600.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the records of one SV-extract file whose read also carries the same SV type in another file (``gafcall.js join``).

    A record of ``extract`` is kept when its read name occurs in ``filter_extract`` with a matching type (insertion or
    duplication, deletion, ...) or when the type is undefined. Result key: ``records`` (lines kept).
    """
    filter_file = _p.file_path(filter_extract, "filter_extract")
    source = _p.file_path(extract, "extract")
    name = _output_name(output_name)
    result = _script(
        "filter_sv_extract_by_type",
        "gafcall.js",
        gafcall_path,
        "join",
        lambda directory: [filter_file, source],
        output_directory,
        name,
        k8_path=k8_path,
        timeout_seconds=timeout_seconds,
        with_version=True,
    )
    result["records"] = _data_rows(result)[1]
    return result


__all__ = [
    "add_sample_names_to_calls",
    "annotate_inserted_sequences",
    "blacklist_insert_alignments",
    "build_sequence_graph",
    "call_graph_bubbles",
    "convert_calls_to_vcf",
    "convert_to_stable_gaf",
    "count_gene_copies",
    "evaluate_sv_calls",
    "extract_gaf_segment_path",
    "extract_long_indels",
    "filter_sv_extract_by_type",
    "graph_segment_coverage",
    "map_sequences_to_graph",
    "merge_bubble_calls",
    "merge_sv_evidence",
    "rename_fasta_sequences",
    "segment_allele_frequency",
    "summarize_insert_annotation",
    "view_sv_calls",
]
_register_public_verbose(globals())
