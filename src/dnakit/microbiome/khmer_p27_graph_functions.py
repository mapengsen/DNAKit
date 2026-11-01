"""khmer k-mer table queries: hashes, abundances of a sequence, de Bruijn graph neighbourhoods and paths (pack P27).

Every function builds the khmer table it needs from one input: a saved table (``graph_path``), local reads
(``reads_path``, FASTA or FASTQ, or a list of files) or DNA strings (``sequences``). Counting questions use a
countgraph, topology questions a nodegraph; ``graph_type`` selects it. Defaults follow khmer's scripts:
k = 32, 4 tables, 1e6 buckets per table. A k-mer and its reverse complement are one node. khmer is imported
only when a function runs (``pip install "dnakit[bio-khmer]"``). The hash values returned are khmer's 2-bit
``forward_hash`` values of the canonical k-mer, not MurmurHash3 (use ``hash_kmer_murmur3`` for that).
"""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

from . import _khmer_p27_common as c

_TABLE_KINDS = ("countgraph", "nodegraph", "counttable", "nodetable")


def _graph(
    kind: str,
    graph_path: Any,
    reads_path: Any,
    sequences: Any,
    ksize: int,
    table_size: int,
    n_tables: int,
    max_reads: int,
    max_bases: int,
) -> tuple[Any, dict[str, Any]]:
    return c.make_graph(
        kind=kind,
        graph_path=graph_path,
        reads_path=reads_path,
        sequences=sequences,
        ksize=ksize,
        table_size=table_size,
        n_tables=n_tables,
        max_reads=max_reads,
        max_bases=max_bases,
    )


def khmer_graph_summary(
    graph_path: str | Path,
    *,
    graph_type: str = "countgraph",
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a saved khmer table and estimate how many distinct k-mers it holds (``n_unique_kmers``).

    ``graph_type`` is the class that wrote the file: ``countgraph`` (``khmer.load_countgraph``),
    ``nodegraph`` (``load_nodegraph``), ``counttable`` (``load_counttable``) or ``nodetable``
    (``load_nodetable``). khmer's own ``n_unique_kmers()`` counter is not saved with the table and reads 0
    after loading, so ``estimated_unique_kmers`` is computed from the occupancy of the first hash table by
    linear counting (``n_occupied`` is the number of occupied buckets there).
    """

    kind = c.choice(graph_type, _TABLE_KINDS, "graph_type")
    path = input_path(graph_path)
    khmer = require("khmer", c.KHMER_EXTRA)
    try:
        if kind in ("countgraph", "counttable"):
            info = khmer.extract_countgraph_info(str(path))
            flag = int(info.ht_type)
        else:
            flag = int(khmer.extract_nodegraph_info(str(path))[4])
        if flag != (1 if kind in ("countgraph", "counttable") else 2):
            raise ValueError("the file header is of another table family")
        table = getattr(khmer, f"load_{kind}")(str(path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer {kind}: {exc}") from exc
    sizes = [int(value) for value in table.hashsizes()]
    occupied = int(table.n_occupied())
    estimate = None
    if 0 < occupied < sizes[0]:
        estimate = round(-sizes[0] * math.log(1.0 - occupied / sizes[0]))
    elif occupied == 0:
        estimate = 0
    return {
        "status": "completed",
        "operation": "khmer_graph_summary",
        "source": str(path),
        "graph_type": kind,
        "ksize": int(table.ksize()),
        "hashsizes": sizes,
        "n_occupied": occupied,
        "estimated_unique_kmers": estimate,
        "provenance": c.provenance(
            loader=f"khmer.load_{kind}",
            convention="khmer's n_unique_kmers() counter is not stored in the file and reads 0 after "
            "loading, so the estimate is linear counting on the first table: "
            "-size * ln(1 - n_occupied / size); None when that table is full",
        ),
    }


def khmer_kmer_hashes(
    sequence: str,
    *,
    ksize: int = c.DEFAULT_KSIZE,
    as_set: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hash values of every k-mer of a DNA sequence, in order (khmer ``get_kmer_hashes``).

    ``as_set=False`` returns the ordered list of hash values (``Countgraph.get_kmer_hashes``);
    ``as_set=True`` returns the sorted distinct values (``get_kmer_hashes_as_hashset``). The values are
    khmer's 2-bit ``forward_hash`` of the canonical k-mer; ``kmers`` gives the k-mers they stand for.
    """

    k = c.ksize_of(ksize)
    text = c.dna(sequence, "sequence", k)
    use_set = c.flag(as_set, "as_set")
    khmer = require("khmer", c.KHMER_EXTRA)
    graph = khmer.Countgraph(k, 1000, 1)
    if use_set:
        hashes = sorted(int(v) for v in graph.get_kmer_hashes_as_hashset(text))
    else:
        hashes = [int(v) for v in graph.get_kmer_hashes(text)]
    return {
        "status": "completed",
        "operation": "khmer_kmer_hashes",
        "sequence_length": len(text),
        "ksize": k,
        "as_set": use_set,
        "n_hashes": len(hashes),
        "hashes": hashes,
        "kmers": [graph.reverse_hash(v) for v in hashes],
        "provenance": c.provenance(
            source_function="Countgraph.get_kmer_hashes_as_hashset"
            if use_set
            else "Countgraph.get_kmer_hashes",
            convention="input is uppercased; hashes are the canonical 2-bit forward hash",
        ),
    }


def khmer_sequence_abundance(
    *,
    sequences_to_query: Sequence[str] | None = None,
    query_path: str | Path | None = None,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    median_threshold: int | None = None,
    include_counts: bool = False,
    output_csv: str | Path | None = None,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """K-mer abundance statistics of query sequences against a countgraph (khmer ``count-median.py``).

    The countgraph comes from ``graph_path`` (saved countgraph), ``reads_path`` or ``sequences``. Query
    sequences are given as ``sequences_to_query`` (list of strings) and/or ``query_path`` (FASTA/FASTQ,
    N read as A, sequences shorter than k skipped, as ``count-median.py`` does). For every query the result
    holds ``min``/``max`` (``get_min_count``/``get_max_count``), ``median``, ``mean`` and ``sd``
    (``get_median_count``), optionally the per-k-mer ``counts`` (``get_kmer_counts``) and, when
    ``median_threshold`` is given, ``median_at_least`` (``median_at_least``). ``output_csv`` writes the
    ``count-median.py`` table (name, median, average, stddev, seqlen).
    """

    k_args = c.table_args(ksize, table_size, n_tables)
    threshold = (
        None
        if median_threshold is None
        else c.integer(median_threshold, "median_threshold", 0, 65535)
    )
    counts_on = c.flag(include_counts, "include_counts")
    if sequences_to_query is None and query_path is None:
        raise ConfigurationError("give sequences_to_query or query_path.")
    khmer = require("khmer", c.KHMER_EXTRA)
    graph, info = _graph(
        "countgraph",
        graph_path,
        reads_path,
        sequences,
        k_args[0],
        k_args[1],
        k_args[2],
        max_reads,
        max_bases,
    )
    k = int(graph.ksize())
    queries: list[tuple[str, str]] = []
    if sequences_to_query is not None:
        texts = c.sequence_list(sequences_to_query, "sequences_to_query")
        queries.extend((f"query{index}", text) for index, text in enumerate(texts))
    if query_path is not None:
        path = input_path(query_path)
        c.scan_reads(
            khmer,
            path,
            c.integer(max_reads, "max_reads", 1, 10**12),
            c.integer(max_bases, "max_bases", 1, 10**12),
        )
        for read in khmer.ReadParser(str(path)):
            queries.append((read.name, read.sequence.upper().replace("N", "A")))
    rows: list[dict[str, Any]] = []
    skipped_short = 0
    for name, text in queries:
        if len(text) < k:
            skipped_short += 1
            continue
        try:
            median, mean, sd = graph.get_median_count(text)
            row: dict[str, Any] = {
                "name": name,
                "length": len(text),
                "min": int(graph.get_min_count(text)),
                "max": int(graph.get_max_count(text)),
                "median": int(median),
                "mean": round(float(mean), 9),
                "sd": round(float(sd), 9),
            }
            if threshold is not None:
                row["median_at_least"] = bool(graph.median_at_least(text, threshold))
            if counts_on:
                row["counts"] = [int(v) for v in graph.get_kmer_counts(text)]
        except (ValueError, RuntimeError) as exc:
            raise ConfigurationError(f"khmer could not count {name}: {exc}") from exc
        rows.append(row)
    written = None
    if output_csv is not None:
        target = c.new_output(output_csv, "output_csv")
        with open(target, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["name", "median", "average", "stddev", "seqlen"])
            for row in rows:
                writer.writerow([row["name"], row["median"], row["mean"], row["sd"], row["length"]])
        written = str(target)
    return {
        "status": "completed",
        "operation": "khmer_sequence_abundance",
        "graph": info,
        "n_queries": len(rows),
        "skipped_shorter_than_k": skipped_short,
        "median_threshold": threshold,
        "results": rows,
        "output_csv": written,
        "provenance": c.provenance(
            source_function="Countgraph.get_median_count / get_min_count / get_max_count / get_kmer_counts / median_at_least",
            convention="counts are of canonical k-mers; a hash collision can only raise a count",
        ),
    }


def khmer_trim_by_abundance(
    sequence: str,
    cutoff: int,
    *,
    mode: str = "trim_below_abundance",
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim one sequence at the first k-mer whose countgraph abundance crosses ``cutoff``.

    ``mode="trim_below_abundance"`` calls ``Countgraph.trim_below_abundance`` and cuts at the first k-mer
    whose count is above the cutoff (high-abundance trimming); ``mode="trim_on_abundance"`` calls
    ``trim_on_abundance`` and cuts at the first k-mer below the cutoff (low-abundance trimming). The
    countgraph comes from ``graph_path``, ``reads_path`` or ``sequences``. Returns the trimmed sequence and
    the trim position (khmer returns the whole sequence length when nothing is trimmed).
    """

    method = c.choice(mode, ("trim_below_abundance", "trim_on_abundance"), "mode")
    limit = c.integer(cutoff, "cutoff", 0, 65535)
    k_args = c.table_args(ksize, table_size, n_tables)
    graph, info = _graph(
        "countgraph",
        graph_path,
        reads_path,
        sequences,
        k_args[0],
        k_args[1],
        k_args[2],
        max_reads,
        max_bases,
    )
    text = c.dna(sequence, "sequence", int(graph.ksize()))
    trimmed, position = getattr(graph, method)(text, limit)
    return {
        "status": "completed",
        "operation": "khmer_trim_by_abundance",
        "graph": info,
        "mode": method,
        "cutoff": limit,
        "input_length": len(text),
        "trimmed_sequence": trimmed,
        "trim_position": int(position),
        "provenance": c.provenance(source_function=f"Countgraph.{method}"),
    }


def khmer_spectral_error_positions(
    sequence: str,
    *,
    max_abundance: int = 1,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Positions of a sequence that lie in low-abundance k-mers, using ``find_spectral_error_positions``.

    A position is reported where khmer's k-mer spectrum shows an abundance of at most ``max_abundance``
    next to higher-abundance k-mers, which marks a probable sequencing error. The countgraph comes from
    ``graph_path``, ``reads_path`` or ``sequences``; positions are 0-based.
    """

    limit = c.integer(max_abundance, "max_abundance", 0, 65535)
    k_args = c.table_args(ksize, table_size, n_tables)
    graph, info = _graph(
        "countgraph",
        graph_path,
        reads_path,
        sequences,
        k_args[0],
        k_args[1],
        k_args[2],
        max_reads,
        max_bases,
    )
    text = c.dna(sequence, "sequence", int(graph.ksize()))
    positions = [int(p) for p in graph.find_spectral_error_positions(text, limit)]
    return {
        "status": "completed",
        "operation": "khmer_spectral_error_positions",
        "graph": info,
        "max_abundance": limit,
        "positions": positions,
        "provenance": c.provenance(
            source_function="Countgraph.find_spectral_error_positions",
            coordinate_system="0-based position in the query sequence",
        ),
    }


def khmer_merge_nodegraphs(
    first_graph_path: str | Path,
    second_graph_path: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Union of two saved nodegraphs (``Nodegraph.update``) written to a new file.

    Both graphs must have the same k-mer size, number of tables and table sizes. The output holds every
    k-mer present in either input.
    """

    first = input_path(first_graph_path)
    second = input_path(second_graph_path)
    target = c.new_output(output_path)
    khmer = require("khmer", c.KHMER_EXTRA)
    left = c.load_saved(khmer, first, "nodegraph")
    right = c.load_saved(khmer, second, "nodegraph")
    before = int(left.n_unique_kmers())
    try:
        left.update(right)
    except (ValueError, RuntimeError, OSError) as exc:
        raise ConfigurationError(f"the nodegraphs cannot be merged: {exc}") from exc
    left.save(str(target))
    return {
        "status": "completed",
        "operation": "khmer_merge_nodegraphs",
        "first": str(first),
        "second": str(second),
        "output_path": str(target),
        "ksize": int(left.ksize()),
        "unique_kmers_first": before,
        "unique_kmers_second": int(right.n_unique_kmers()),
        "unique_kmers_merged": int(left.n_unique_kmers()),
        "provenance": c.provenance(source_function="Nodegraph.update"),
    }


def _topology(
    kmer: str,
    graph_path: Any,
    reads_path: Any,
    sequences: Any,
    ksize: int,
    table_size: int,
    n_tables: int,
    graph_type: str,
    max_reads: int,
    max_bases: int,
) -> tuple[Any, dict[str, Any], str]:
    k_args = c.table_args(ksize, table_size, n_tables)
    graph, info = _graph(
        graph_type,
        graph_path,
        reads_path,
        sequences,
        k_args[0],
        k_args[1],
        k_args[2],
        max_reads,
        max_bases,
    )
    return graph, info, c.kmer_of(kmer, int(graph.ksize()))


def khmer_kmer_neighbors(
    kmer: str,
    *,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    graph_type: str = "nodegraph",
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """K-mers adjacent to ``kmer`` in the de Bruijn graph of the input (``neighbors``), with their degree.

    The input graph comes from ``graph_path`` (saved graph of ``graph_type``), ``reads_path`` or
    ``sequences``. ``degree`` is ``kmer_degree`` (number of immediate neighbours), so this one function
    answers both the neighbour list and the degree of a k-mer.
    """

    graph, info, text = _topology(
        kmer,
        graph_path,
        reads_path,
        sequences,
        ksize,
        table_size,
        n_tables,
        graph_type,
        max_reads,
        max_bases,
    )
    neighbours = [int(h) for h in graph.neighbors(text)]
    return {
        "status": "completed",
        "operation": "khmer_kmer_neighbors",
        "graph": info,
        "kmer": text,
        "degree": int(graph.kmer_degree(text)),
        "neighbors": c.kmers_of_hashes(graph, neighbours),
        "provenance": c.provenance(source_function="neighbors; kmer_degree"),
    }


def khmer_connected_component_size(
    kmer: str,
    *,
    max_size: int = 0,
    count_within_radius: int | None = None,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    graph_type: str = "nodegraph",
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of k-mers connected to ``kmer`` in the de Bruijn graph (``calc_connected_graph_size``).

    ``max_size`` stops the traversal once that many k-mers are counted (0 means no limit, as in khmer).
    With ``count_within_radius`` set, the number of k-mers within that graph radius is also returned
    (``count_kmers_within_radius``); the radius count is bounded by ``max_size`` when it is positive.
    """

    limit = c.integer(max_size, "max_size", 0, 10**12)
    radius = (
        None
        if count_within_radius is None
        else c.integer(count_within_radius, "count_within_radius", 0, 10**9)
    )
    graph, info, text = _topology(
        kmer,
        graph_path,
        reads_path,
        sequences,
        ksize,
        table_size,
        n_tables,
        graph_type,
        max_reads,
        max_bases,
    )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "khmer_connected_component_size",
        "graph": info,
        "kmer": text,
        "max_size": limit,
        "component_size": int(graph.calc_connected_graph_size(text, limit)),
        "count_within_radius": radius,
        "kmers_within_radius": None,
        "provenance": c.provenance(
            source_function="calc_connected_graph_size; count_kmers_within_radius"
        ),
    }
    if radius is not None:
        result["kmers_within_radius"] = int(graph.count_kmers_within_radius(text, radius, limit))
    return result


def khmer_high_degree_kmers(
    sequence: str,
    *,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    graph_type: str = "nodegraph",
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Branch points: k-mers of ``sequence`` with more than two neighbours in the graph (``find_high_degree_nodes``)."""

    k_args = c.table_args(ksize, table_size, n_tables)
    graph, info = _graph(
        graph_type,
        graph_path,
        reads_path,
        sequences,
        k_args[0],
        k_args[1],
        k_args[2],
        max_reads,
        max_bases,
    )
    text = c.dna(sequence, "sequence", int(graph.ksize()))
    hashes = sorted(int(h) for h in graph.find_high_degree_nodes(text))
    return {
        "status": "completed",
        "operation": "khmer_high_degree_kmers",
        "graph": info,
        "n_high_degree": len(hashes),
        "high_degree_kmers": c.kmers_of_hashes(graph, hashes),
        "provenance": c.provenance(source_function="find_high_degree_nodes"),
    }


def khmer_linear_path(
    kmer: str,
    *,
    output: str = "sequence",
    stop_sequence: str | None = None,
    graph_path: str | Path | None = None,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    graph_type: str = "nodegraph",
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Follow the unbranched path that starts at ``kmer`` in the de Bruijn graph.

    ``output="sequence"`` assembles the path into one DNA string (``assemble_linear_path``);
    ``output="kmers"`` returns the traversed k-mers and the branch k-mers met at its ends
    (``traverse_linear_path``). With ``output="kmers"``, ``stop_sequence`` supplies the branch points that
    end the traversal: the high-degree k-mers of that sequence (``find_high_degree_nodes``); without it
    the traversal follows the whole connected component. Assembly ends at a branch point, a dead end
    or after a loop.
    """

    mode = c.choice(output, ("sequence", "kmers"), "output")
    graph, info, text = _topology(
        kmer,
        graph_path,
        reads_path,
        sequences,
        ksize,
        table_size,
        n_tables,
        graph_type,
        max_reads,
        max_bases,
    )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "khmer_linear_path",
        "graph": info,
        "start_kmer": text,
        "output": mode,
        "provenance": c.provenance(
            source_function="assemble_linear_path" if mode == "sequence" else "traverse_linear_path"
        ),
    }
    if mode == "sequence":
        path = graph.assemble_linear_path(text)
        result["path_sequence"] = path
        result["path_length"] = len(path)
        return result
    khmer = require("khmer", c.KHMER_EXTRA)
    if stop_sequence is None:
        degree_nodes = khmer.HashSet(int(graph.ksize()))
    else:
        degree_nodes = graph.find_high_degree_nodes(
            c.dna(stop_sequence, "stop_sequence", int(graph.ksize()))
        )
    stop = khmer.Nodegraph(int(graph.ksize()), 1000, 1)
    size, connections, visited = graph.traverse_linear_path(text, degree_nodes, stop)
    result["path_size"] = int(size)
    result["branch_kmers"] = c.kmers_of_hashes(graph, sorted(int(h) for h in connections))
    result["path_kmers"] = c.kmers_of_hashes(graph, sorted(int(h) for h in visited))
    return result


__all__ = [
    "khmer_connected_component_size",
    "khmer_graph_summary",
    "khmer_high_degree_kmers",
    "khmer_kmer_hashes",
    "khmer_kmer_neighbors",
    "khmer_linear_path",
    "khmer_merge_nodegraphs",
    "khmer_sequence_abundance",
    "khmer_spectral_error_positions",
    "khmer_trim_by_abundance",
]

_register_public_verbose(globals(), tuple(__all__))
