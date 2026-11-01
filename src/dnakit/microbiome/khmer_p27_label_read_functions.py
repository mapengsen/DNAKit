"""khmer labelled graphs, HyperLogLog merging, read-to-graph alignment and read-file utilities (pack P27).

Labelled functions build a khmer ``GraphLabels`` object: a nodegraph whose tag k-mers carry integer labels.
Reads are given as local FASTA/FASTQ files (every read gets its record index in the file as label, khmer's
``consume_seqfile_and_tag_with_labels``) or as DNA strings with explicit labels. The read-file functions
reproduce khmer's ``extract-paired-reads.py``, ``split-paired-reads.py`` and ``readstats.py``. Defaults
follow khmer: k = 32, 4 tables, 1e6 buckets per table. khmer is imported only when a function runs
(``pip install "dnakit[bio-khmer]"``).
"""

from __future__ import annotations

import contextlib
import csv
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

from . import _khmer_p27_common as c


def _labels_for(sequences: Sequence[str], labels: Sequence[int] | None) -> list[int]:
    if labels is None:
        return list(range(len(sequences)))
    if len(labels) != len(sequences):
        raise ConfigurationError("labels must have one integer per sequence.")
    return [c.integer(value, "labels item", 0, 2**32 - 2) for value in labels]


def _labeled_graph(
    khmer: Any,
    *,
    reads_path: Any,
    sequences: Sequence[str] | None,
    labels: Sequence[int] | None,
    ksize: int,
    table_size: int,
    n_tables: int,
    max_reads: int,
    max_bases: int,
) -> tuple[Any, dict[str, Any]]:
    k_args = c.table_args(ksize, table_size, n_tables)
    if reads_path is None and sequences is None:
        raise ConfigurationError("give reads_path or sequences.")
    if labels is not None and sequences is None:
        raise ConfigurationError(
            "labels apply to sequences; reads_path files label by record index."
        )
    graph_labels = khmer.GraphLabels(*k_args)
    info: dict[str, Any] = {
        "ksize": k_args[0],
        "table_size": k_args[1],
        "n_tables": k_args[2],
        "reads": 0,
        "bases": 0,
    }
    limits = (
        c.integer(max_reads, "max_reads", 1, 10**12),
        c.integer(max_bases, "max_bases", 1, 10**12),
    )
    if reads_path is not None:
        items = [reads_path] if isinstance(reads_path, (str, Path)) else list(reads_path)
        for item in items:
            path = input_path(item)
            scanned = c.scan_reads(khmer, path, *limits)
            info["reads"] += scanned["reads"]
            info["bases"] += scanned["bases"]
            try:
                graph_labels.consume_seqfile_and_tag_with_labels(str(path))
            except (OSError, ValueError, RuntimeError) as exc:
                raise ConfigurationError(f"khmer could not label {path.name}: {exc}") from exc
    if sequences is not None:
        texts = c.sequence_list(sequences, "sequences", k_args[0])
        for text, label in zip(texts, _labels_for(texts, labels), strict=True):
            graph_labels.consume_sequence_and_tag_with_labels(text, label)
            info["reads"] += 1
            info["bases"] += len(text)
        if info["reads"] > limits[0] or info["bases"] > limits[1]:
            raise ConfigurationError(
                "sequence input exceeds max_reads or max_bases; raise the limit."
            )
    return graph_labels, info


def khmer_label_reads(
    *,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    labels: Sequence[int] | None = None,
    load_labels_path: str | Path | None = None,
    save_labels_path: str | Path | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert reads into a labelled graph and tag their k-mers with labels.

    Files in ``reads_path`` are processed with ``consume_seqfile_and_tag_with_labels`` (label = index of
    the record in its file); ``sequences`` with ``consume_sequence_and_tag_with_labels`` (label from
    ``labels``, default the position in the list). ``load_labels_path`` first loads a saved label and tag
    set (``load_labels_and_tags``); ``save_labels_path`` writes the result (``save_labels_and_tags``).
    Returns the number of labels and the labels carried by each tag k-mer.
    """

    khmer = require("khmer", c.KHMER_EXTRA)
    saved = None if save_labels_path is None else c.new_output(save_labels_path, "save_labels_path")
    if reads_path is None and sequences is None:
        if load_labels_path is None:
            raise ConfigurationError("give reads_path, sequences or load_labels_path.")
        k_args = c.table_args(ksize, table_size, n_tables)
        graph = khmer.GraphLabels(*k_args)
        info: dict[str, Any] = {"ksize": k_args[0], "reads": 0, "bases": 0}
    else:
        graph, info = _labeled_graph(
            khmer,
            reads_path=reads_path,
            sequences=sequences,
            labels=labels,
            ksize=ksize,
            table_size=table_size,
            n_tables=n_tables,
            max_reads=max_reads,
            max_bases=max_bases,
        )
    loaded = None
    if load_labels_path is not None:
        path = input_path(load_labels_path)
        try:
            graph.load_labels_and_tags(str(path))
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(
                f"{path.name} is not a khmer labels-and-tags file: {exc}"
            ) from exc
        loaded = str(path)
    if saved is not None:
        graph.save_labels_and_tags(str(saved))
    tag_labels = {
        str(tag): [int(label) for label in graph.get_tag_labels(str(tag))]
        for tag in sorted(graph.graph.get_tagset())[: c.MAX_LISTED]
    }
    return {
        "status": "completed",
        "operation": "khmer_label_reads",
        "graph": info,
        "loaded_labels": loaded,
        "n_labels": int(graph.n_labels()),
        "n_tags": int(graph.graph.n_tags()),
        "tag_labels": tag_labels,
        "saved_labels": None if saved is None else str(saved),
        "provenance": c.provenance(
            source_function="GraphLabels.consume_seqfile_and_tag_with_labels / consume_sequence_and_tag_with_labels; "
            "save_labels_and_tags; load_labels_and_tags"
        ),
    }


def khmer_label_partitioned_reads(
    partitioned_path: str | Path,
    *,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert partition-annotated reads into a labelled graph, using the partition ID as label.

    Calls ``consume_partitioned_fasta_and_tag_with_labels`` on a ``.part`` file written by
    ``khmer_partition_reads``. Returns the number of labels, tags and the labels of every tag.
    """

    path = input_path(partitioned_path)
    k_args = c.table_args(ksize, table_size, n_tables)
    khmer = require("khmer", c.KHMER_EXTRA)
    graph = khmer.GraphLabels(*k_args)
    try:
        graph.consume_partitioned_fasta_and_tag_with_labels(str(path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(
            f"{path.name} is not a partition-annotated read file: {exc}"
        ) from exc
    tag_labels = {
        str(tag): [int(label) for label in graph.get_tag_labels(str(tag))]
        for tag in sorted(graph.graph.get_tagset())[: c.MAX_LISTED]
    }
    return {
        "status": "completed",
        "operation": "khmer_label_partitioned_reads",
        "source": str(path),
        "ksize": k_args[0],
        "n_labels": int(graph.n_labels()),
        "n_tags": int(graph.graph.n_tags()),
        "tag_labels": tag_labels,
        "provenance": c.provenance(
            source_function="GraphLabels.consume_partitioned_fasta_and_tag_with_labels"
        ),
    }


def khmer_sweep_neighborhood(
    sequence: str,
    *,
    output: str = "labels",
    sweep_range: int = 0,
    break_on_stop_tags: bool = False,
    stop_big_traversals: bool = False,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    labels: Sequence[int] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Labels or tags found around a query sequence in a labelled graph.

    ``output="labels"`` calls ``sweep_label_neighborhood`` and returns the labels of tags near the
    query; ``output="tags"`` calls ``sweep_tag_neighborhood`` and returns the tag k-mers. ``sweep_range``
    is how far (in k-mers) the traversal goes from the query (0 keeps khmer's default behaviour of
    the nearest tags only). The graph is built as in ``khmer_label_reads``.
    """

    mode = c.choice(output, ("labels", "tags"), "output")
    distance = c.integer(sweep_range, "sweep_range", 0, 10**6)
    stops = c.flag(break_on_stop_tags, "break_on_stop_tags")
    big = c.flag(stop_big_traversals, "stop_big_traversals")
    khmer = require("khmer", c.KHMER_EXTRA)
    graph, info = _labeled_graph(
        khmer,
        reads_path=reads_path,
        sequences=sequences,
        labels=labels,
        ksize=ksize,
        table_size=table_size,
        n_tables=n_tables,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    text = c.dna(sequence, "sequence", int(graph.graph.ksize()))
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "khmer_sweep_neighborhood",
        "graph": info,
        "output": mode,
        "sweep_range": distance,
        "provenance": c.provenance(
            source_function="GraphLabels.sweep_label_neighborhood"
            if mode == "labels"
            else "GraphLabels.sweep_tag_neighborhood"
        ),
    }
    if mode == "labels":
        result["labels"] = sorted(
            int(v) for v in graph.sweep_label_neighborhood(text, distance, stops, big)
        )
    else:
        hashes = sorted(int(v) for v in graph.sweep_tag_neighborhood(text, distance, stops, big))
        result["tags"] = c.kmers_of_hashes(graph.graph, hashes)
    return result


def khmer_tag_labels(
    tag_kmer: str,
    *,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    labels: Sequence[int] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Labels carried by one tag k-mer of a labelled graph (``get_tag_labels``); empty if it is not a tag."""

    khmer = require("khmer", c.KHMER_EXTRA)
    graph, info = _labeled_graph(
        khmer,
        reads_path=reads_path,
        sequences=sequences,
        labels=labels,
        ksize=ksize,
        table_size=table_size,
        n_tables=n_tables,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    text = c.kmer_of(tag_kmer, int(graph.graph.ksize()), "tag_kmer")
    found = sorted(int(v) for v in graph.get_tag_labels(text))
    return {
        "status": "completed",
        "operation": "khmer_tag_labels",
        "graph": info,
        "tag_kmer": text,
        "labels": found,
        "provenance": c.provenance(source_function="GraphLabels.get_tag_labels"),
    }


def _contig_graph(
    khmer: Any,
    contigs: Sequence[str],
    labels: Sequence[int] | None,
    ksize: int,
    table_size: int,
    n_tables: int,
) -> tuple[Any, list[str], list[int], Any]:
    k_args = c.table_args(ksize, table_size, n_tables)
    texts = c.sequence_list(contigs, "contigs", k_args[0] + 1)
    ids = _labels_for(texts, None if labels is None else labels)
    if labels is None:
        ids = [index + 1 for index in range(len(texts))]
    graph = khmer.GraphLabels(*k_args)
    for text in texts:
        graph.graph.consume(text)
    hdn = graph.graph.find_high_degree_nodes(texts[0])
    for text in texts[1:]:
        hdn += graph.graph.find_high_degree_nodes(text)
    for text, label in zip(texts, ids, strict=True):
        graph.label_across_high_degree_nodes(text, hdn, label)
    return graph, texts, ids, hdn


def khmer_label_across_high_degree_nodes(
    contigs: Sequence[str],
    *,
    labels: Sequence[int] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Label the k-mers of contigs across the branch points they share (``label_across_high_degree_nodes``).

    All contigs are inserted into one graph; its high-degree k-mers (branch points, ``find_high_degree_nodes``)
    are found, and each contig labels itself and the tags around every branch point it crosses. Labels
    default to 1, 2, ... in list order. Returns the branch k-mers and the labels each one received.
    """

    khmer = require("khmer", c.KHMER_EXTRA)
    graph, texts, ids, hdn = _contig_graph(khmer, contigs, labels, ksize, table_size, n_tables)
    branch_labels = {
        str(graph.graph.reverse_hash(int(h))): sorted(int(v) for v in graph.get_tag_labels(int(h)))
        for h in sorted(int(x) for x in hdn)
    }
    return {
        "status": "completed",
        "operation": "khmer_label_across_high_degree_nodes",
        "n_contigs": len(texts),
        "contig_labels": ids,
        "n_high_degree": len(branch_labels),
        "high_degree_labels": branch_labels,
        "n_labels": int(graph.n_labels()),
        "provenance": c.provenance(source_function="GraphLabels.label_across_high_degree_nodes"),
    }


def khmer_assemble_labeled_path(
    start_kmer: str,
    contigs: Sequence[str],
    *,
    labels: Sequence[int] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assemble the paths through a branching graph that start at ``start_kmer``, resolved by labels.

    The contigs are inserted and labelled as in ``khmer_label_across_high_degree_nodes``; then
    ``assemble_labeled_path`` follows the graph from ``start_kmer`` and uses the labels to choose between
    branches, returning one sequence per consistent path.
    """

    khmer = require("khmer", c.KHMER_EXTRA)
    graph, _texts, ids, _hdn = _contig_graph(khmer, contigs, labels, ksize, table_size, n_tables)
    start = c.kmer_of(start_kmer, c.ksize_of(ksize), "start_kmer")
    paths = sorted(str(p) for p in graph.assemble_labeled_path(start))
    return {
        "status": "completed",
        "operation": "khmer_assemble_labeled_path",
        "start_kmer": start,
        "contig_labels": ids,
        "n_paths": len(paths),
        "paths": paths,
        "provenance": c.provenance(source_function="GraphLabels.assemble_labeled_path"),
    }


def _hll(
    khmer: Any, error_rate: float, ksize: int, paths: Sequence[Path], texts: Sequence[str]
) -> Any:
    counter = khmer.HLLCounter(error_rate, ksize)
    for path in paths:
        try:
            counter.consume_seqfile(str(path))
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(f"khmer could not read {path.name}: {exc}") from exc
    for text in texts:
        counter.consume_string(text)
    return counter


def khmer_merge_hll_counters(
    *,
    first_reads_path: str | Path | None = None,
    first_sequences: Sequence[str] | None = None,
    second_reads_path: str | Path | None = None,
    second_sequences: Sequence[str] | None = None,
    error_rate: float = 0.01,
    ksize: int = 20,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate distinct k-mers of two read sets and of their union by merging HyperLogLog counters (``HLLCounter.merge``).

    Two counters are filled, one per input (file and/or strings), then the second is merged into the
    first. Defaults follow khmer's ``unique-kmers.py``: error rate 0.01, k = 20. Both counters share the
    same error rate and k, as ``merge`` requires.
    """

    rate = c.number(error_rate, "error_rate", 0.0001, 0.367)
    k = c.integer(ksize, "ksize", 1, 32)
    sides: list[tuple[list[Path], list[str]]] = []
    for name, reads, texts in (
        ("first", first_reads_path, first_sequences),
        ("second", second_reads_path, second_sequences),
    ):
        if reads is None and texts is None:
            raise ConfigurationError(f"give {name}_reads_path or {name}_sequences.")
        paths = [] if reads is None else [input_path(reads)]
        strings = [] if texts is None else c.sequence_list(texts, f"{name}_sequences", k)
        sides.append((paths, strings))
    khmer = require("khmer", c.KHMER_EXTRA)
    first = _hll(khmer, rate, k, *sides[0])
    second = _hll(khmer, rate, k, *sides[1])
    first_estimate = int(first.estimate_cardinality())
    second_estimate = int(second.estimate_cardinality())
    first.merge(second)
    return {
        "status": "completed",
        "operation": "khmer_merge_hll_counters",
        "error_rate": rate,
        "ksize": k,
        "first_estimate": first_estimate,
        "second_estimate": second_estimate,
        "merged_estimate": int(first.estimate_cardinality()),
        "provenance": c.provenance(
            source_function="HLLCounter.consume_seqfile / consume_string / merge / estimate_cardinality",
            convention="estimates are HyperLogLog approximations with the stated relative error",
        ),
    }


def khmer_align_read(
    read: str,
    *,
    forward_only: bool = False,
    trusted_cov_cutoff: int = 2,
    bits_theta: float = 1.0,
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
    """Align a read to the k-mer graph of a countgraph with khmer's hidden Markov model (``ReadAligner``).

    ``forward_only=False`` calls ``ReadAligner.align``, which seeds on the read and extends in both
    directions; ``forward_only=True`` calls ``align_forward``, which extends forward only. The countgraph
    comes from ``graph_path``, ``reads_path`` or ``sequences``. ``trusted_cov_cutoff`` is the abundance at
    which a k-mer counts as trusted (khmer's sandbox default 2) and ``bits_theta`` the score threshold
    (default 1). Returns the alignment score, the aligned graph and read strings (gaps as ``-``) and
    whether the alignment was truncated.
    """

    forward = c.flag(forward_only, "forward_only")
    cutoff = c.integer(trusted_cov_cutoff, "trusted_cov_cutoff", 0, 65535)
    theta = c.number(bits_theta, "bits_theta", 0.0, 1e6)
    k_args = c.table_args(ksize, table_size, n_tables)
    graph, info = c.make_graph(
        kind="countgraph",
        graph_path=graph_path,
        reads_path=reads_path,
        sequences=sequences,
        ksize=k_args[0],
        table_size=k_args[1],
        n_tables=k_args[2],
        max_reads=max_reads,
        max_bases=max_bases,
    )
    text = c.dna(read, "read", int(graph.ksize()))
    khmer = require("khmer", c.KHMER_EXTRA)
    aligner = khmer.ReadAligner(graph, cutoff, theta)
    outcome = aligner.align_forward(text) if forward else aligner.align(text)
    score, graph_alignment, read_alignment, truncated = outcome[:4]
    return {
        "status": "completed",
        "operation": "khmer_align_read",
        "graph": info,
        "forward_only": forward,
        "trusted_cov_cutoff": cutoff,
        "bits_theta": theta,
        "score": float(score),
        "graph_alignment": graph_alignment,
        "read_alignment": read_alignment,
        "truncated": bool(truncated),
        "provenance": c.provenance(
            source_function="ReadAligner.align_forward" if forward else "ReadAligner.align",
            convention="khmer's default scoring matrix and transition probabilities",
        ),
    }


def khmer_extract_paired_reads(
    reads_path: str | Path,
    paired_output_path: str | Path,
    single_output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Separate a mixed read file into interleaved pairs and orphan reads (``extract-paired-reads.py``).

    Reads are paired by name (``/1`` and ``/2`` suffixes, or Casava ``1:``/``2:`` fields) with
    ``khmer.utils.broken_paired_reader``. Pairs are written interleaved to ``paired_output_path`` and
    unpaired reads to ``single_output_path``, in FASTA or FASTQ as the input is. Output files must not exist.
    khmer's script fails when no pair is found; so does this function.
    """

    source = input_path(reads_path)
    paired = c.new_output(paired_output_path, "paired_output_path")
    single = c.new_output(single_output_path, "single_output_path")
    if paired == single:
        raise ConfigurationError("paired_output_path and single_output_path must differ.")
    khmer = require("khmer", c.KHMER_EXTRA)
    utils = require("khmer.utils", c.KHMER_EXTRA)
    n_pairs = n_single = 0
    try:
        with open(paired, "wb") as pe, open(single, "wb") as se:
            for _index, is_pair, first, second in utils.broken_paired_reader(
                khmer.ReadParser(str(source))
            ):
                if is_pair:
                    utils.write_record_pair(first, second, pe)
                    n_pairs += 1
                else:
                    utils.write_record(first, se)
                    n_single += 1
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not read {source.name}: {exc}") from exc
    if n_pairs == 0:
        paired.unlink(missing_ok=True)
        single.unlink(missing_ok=True)
        raise ConfigurationError("no paired reads found; check the read names and the file format.")
    return {
        "status": "completed",
        "operation": "khmer_extract_paired_reads",
        "source": str(source),
        "paired_output_path": str(paired),
        "single_output_path": str(single),
        "n_reads": n_pairs * 2 + n_single,
        "n_pairs": n_pairs,
        "n_singletons": n_single,
        "provenance": c.provenance(source_function="khmer.utils.broken_paired_reader"),
    }


def khmer_split_paired_reads(
    reads_path: str | Path,
    first_output_path: str | Path,
    second_output_path: str | Path,
    *,
    orphan_output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split an interleaved paired-read file into a left-read and a right-read file (``split-paired-reads.py``).

    Without ``orphan_output_path`` an unpaired read is an error (khmer's ``require_paired``); with it,
    orphans are written there. Output files must not exist.
    """

    source = input_path(reads_path)
    left = c.new_output(first_output_path, "first_output_path")
    right = c.new_output(second_output_path, "second_output_path")
    orphan = (
        None
        if orphan_output_path is None
        else c.new_output(orphan_output_path, "orphan_output_path")
    )
    given = [left, right] + ([] if orphan is None else [orphan])
    if len(set(given)) != len(given):
        raise ConfigurationError("the output paths must all differ.")
    khmer = require("khmer", c.KHMER_EXTRA)
    utils = require("khmer.utils", c.KHMER_EXTRA)
    n_left = n_right = n_orphan = 0
    try:
        with contextlib.ExitStack() as stack:
            fp_left = stack.enter_context(open(left, "wb"))
            fp_right = stack.enter_context(open(right, "wb"))
            fp_orphan = None if orphan is None else stack.enter_context(open(orphan, "wb"))
            for _index, is_pair, first, second in utils.broken_paired_reader(
                khmer.ReadParser(str(source)), require_paired=orphan is None
            ):
                if is_pair:
                    utils.write_record(first, fp_left)
                    utils.write_record(second, fp_right)
                    n_left += 1
                    n_right += 1
                elif fp_orphan is not None:
                    utils.write_record(first, fp_orphan)
                    n_orphan += 1
    except (utils.UnpairedReadsError, OSError, ValueError) as exc:
        for path in given:
            path.unlink(missing_ok=True)
        if isinstance(exc, utils.UnpairedReadsError):
            raise ConfigurationError(
                f"the file holds unpaired reads and no orphan output was given: {exc}"
            ) from exc
        raise ConfigurationError(f"khmer could not read {source.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "khmer_split_paired_reads",
        "source": str(source),
        "first_output_path": str(left),
        "second_output_path": str(right),
        "orphan_output_path": None if orphan is None else str(orphan),
        "n_reads": n_left + n_right + n_orphan,
        "n_pairs": n_left,
        "n_orphans": n_orphan,
        "provenance": c.provenance(source_function="khmer.utils.broken_paired_reader"),
    }


def khmer_read_statistics(
    reads_paths: str | Path | Sequence[str | Path],
    *,
    output_csv: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of bases, sequences and mean sequence length of FASTA/FASTQ files (khmer ``readstats.py``).

    Returns one row per file plus the total, with the average length rounded to one decimal as the
    script prints it. A file that khmer cannot open or that holds no sequence is reported with zero sequences and the
    parser message in ``error``, and is left out of the total, as in the script. ``output_csv`` writes the script's ``--csv`` table (bp, seqs, avg_len, filename).
    """

    items = [reads_paths] if isinstance(reads_paths, (str, Path)) else list(reads_paths)
    if not items:
        raise ConfigurationError("reads_paths must not be empty.")
    files = [input_path(item) for item in items]
    target = None if output_csv is None else c.new_output(output_csv, "output_csv")
    khmer = require("khmer", c.KHMER_EXTRA)
    rows: list[dict[str, Any]] = []
    for path in files:
        bases = seqs = 0
        error = None
        try:
            for read in khmer.ReadParser(str(path)):
                bases += len(read.sequence)
                seqs += 1
        except (OSError, ValueError) as exc:
            error = str(exc)
        rows.append(
            {
                "filename": str(path),
                "bp": bases,
                "seqs": seqs,
                "avg_len": round(bases / seqs, 1) if seqs else None,
                "error": error,
            }
        )
    counted = [row for row in rows if row["seqs"]]
    total_bp = sum(row["bp"] for row in counted)
    total_seqs = sum(row["seqs"] for row in counted)
    if target is not None:
        with open(target, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["bp", "seqs", "avg_len", "filename"])
            for row in counted:
                writer.writerow([row["bp"], row["seqs"], f"{row['avg_len']:.1f}", row["filename"]])
    return {
        "status": "completed",
        "operation": "khmer_read_statistics",
        "files": rows,
        "total": {
            "bp": total_bp,
            "seqs": total_seqs,
            "avg_len": round(total_bp / total_seqs, 1) if total_seqs else None,
        },
        "output_csv": None if target is None else str(target),
        "provenance": c.provenance(source_function="khmer.ReadParser (readstats.py logic)"),
    }


__all__ = [
    "khmer_align_read",
    "khmer_assemble_labeled_path",
    "khmer_extract_paired_reads",
    "khmer_label_across_high_degree_nodes",
    "khmer_label_partitioned_reads",
    "khmer_label_reads",
    "khmer_merge_hll_counters",
    "khmer_read_statistics",
    "khmer_split_paired_reads",
    "khmer_sweep_neighborhood",
    "khmer_tag_labels",
]

_register_public_verbose(globals(), tuple(__all__))
