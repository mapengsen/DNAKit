"""khmer tags, partitions and stop tags: graph partitioning of reads into connected components (pack P27).

khmer marks a sparse set of "tag" k-mers in the de Bruijn graph of the reads, groups tags that are connected
into partitions, and annotates every read with the ID of its partition (khmer's ``do-partition.py``,
``extract-partitions.py``, ``find-knots.py`` and ``filter-stoptags.py`` workflows). "Stop tags" are
highly connected k-mers (knots) at which partitions are cut. Reads are given as local FASTA or FASTQ files
(consumed by khmer itself, so its read cleaning applies) or, for the small queries, as DNA strings. Defaults
follow khmer's scripts: k = 32, 4 tables, 1e6 buckets per table, subset size 1e5. khmer is imported only
when a function runs (``pip install "dnakit[bio-khmer]"``).
"""

from __future__ import annotations

import contextlib
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

from . import _khmer_p27_common as c

_DEFAULT_SUBSET_SIZE = 100_000
_EXCURSION_DISTANCE = 40
_EXCURSION_KMER_THRESHOLD = 200


def _paths(value: str | Path | Sequence[str | Path], name: str) -> list[Path]:
    items = [value] if isinstance(value, (str, Path)) else list(value)
    if not items:
        raise ConfigurationError(f"{name} must not be empty.")
    return [input_path(item) for item in items]


def _tagged(
    *,
    kind: str,
    reads_path: Any,
    sequences: Any,
    ksize: int,
    table_size: int,
    n_tables: int,
    max_reads: int,
    max_bases: int,
) -> tuple[Any, dict[str, Any]]:
    k_args = c.table_args(ksize, table_size, n_tables)
    return c.make_graph(
        kind=kind,
        graph_path=None,
        reads_path=reads_path,
        sequences=sequences,
        ksize=k_args[0],
        table_size=k_args[1],
        n_tables=k_args[2],
        tag=True,
        max_reads=max_reads,
        max_bases=max_bases,
    )


def _load_stop_tags(graph: Any, path: Path) -> None:
    try:
        graph.load_stop_tags(str(path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer stop-tag file: {exc}") from exc


def khmer_tag_reads(
    *,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    load_tagset_path: str | Path | None = None,
    extra_tag_kmers: Sequence[str] | None = None,
    save_tagset_path: str | Path | None = None,
    tagset_text_path: str | Path | None = None,
    graph_type: str = "nodegraph",
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert reads into a khmer graph and mark tag k-mers (``consume_seqfile_and_tag`` / ``consume_and_tag``).

    Tags are sparse k-mers placed along every read so that connected reads can be found later. The graph
    is built from ``reads_path`` (file or list of files) and/or ``sequences`` (strings). Optional steps,
    each a khmer call: ``load_tagset_path`` first loads a saved tag set (``load_tagset``);
    ``extra_tag_kmers`` marks more k-mers as tags (``add_tag``); ``save_tagset_path`` writes the tag set
    in khmer's binary format (``save_tagset``); ``tagset_text_path`` writes the tags as text, one k-mer per
    line (``print_tagset``). ``graph_type`` is ``nodegraph`` or ``countgraph``. Returns the tags and the
    number of reads and k-mers consumed.
    """

    kind = c.choice(graph_type, c.GRAPH_KINDS, "graph_type")
    k_args = c.table_args(ksize, table_size, n_tables)
    if (
        reads_path is None
        and sequences is None
        and load_tagset_path is None
        and not extra_tag_kmers
    ):
        raise ConfigurationError("give reads_path, sequences, load_tagset_path or extra_tag_kmers.")
    saved = None if save_tagset_path is None else c.new_output(save_tagset_path, "save_tagset_path")
    text = None if tagset_text_path is None else c.new_output(tagset_text_path, "tagset_text_path")
    khmer = require("khmer", c.KHMER_EXTRA)
    if reads_path is None and sequences is None:
        graph = c.new_graph(khmer, kind, *k_args)
        info: dict[str, Any] = {"source": "none", "ksize": k_args[0]}
    else:
        graph, info = _tagged(
            kind=kind,
            reads_path=reads_path,
            sequences=sequences,
            ksize=k_args[0],
            table_size=k_args[1],
            n_tables=k_args[2],
            max_reads=max_reads,
            max_bases=max_bases,
        )
    loaded = None
    if load_tagset_path is not None:
        tag_file = input_path(load_tagset_path)
        try:
            graph.load_tagset(str(tag_file), False)
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(
                f"{tag_file.name} is not a khmer tag set for this graph: {exc}"
            ) from exc
        loaded = str(tag_file)
    added = 0
    for kmer in [] if extra_tag_kmers is None else list(extra_tag_kmers):
        graph.add_tag(c.kmer_of(kmer, int(graph.ksize()), "extra_tag_kmers item"))
        added += 1
    tags = sorted(str(t) for t in graph.get_tagset())
    if saved is not None:
        graph.save_tagset(str(saved))
    if text is not None:
        graph.print_tagset(str(text))
    return {
        "status": "completed",
        "operation": "khmer_tag_reads",
        "graph": info,
        "loaded_tagset": loaded,
        "extra_tags_added": added,
        "n_tags": int(graph.n_tags()),
        "tags": tags[: c.MAX_LISTED],
        "n_unique_kmers": int(graph.n_unique_kmers()),
        "saved_tagset": None if saved is None else str(saved),
        "tagset_text": None if text is None else str(text),
        "provenance": c.provenance(
            source_function="consume_seqfile_and_tag / consume_and_tag; add_tag; load_tagset; save_tagset; print_tagset",
            convention="a tag is a canonical k-mer; the tag density is khmer's default (one tag about every 40 bases)",
        ),
    }


def khmer_tag_positions(
    sequence: str,
    *,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    graph_type: str = "nodegraph",
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Positions of the tag k-mers along ``sequence`` in a tagged graph (``get_tags_and_positions``).

    The graph is built and tagged from ``reads_path`` and/or ``sequences``; the query sequence need not be
    one of them. Positions are 0-based start offsets of the tag k-mer in ``sequence``.
    """

    graph, info = _tagged(
        kind=c.choice(graph_type, c.GRAPH_KINDS, "graph_type"),
        reads_path=reads_path,
        sequences=sequences,
        ksize=ksize,
        table_size=table_size,
        n_tables=n_tables,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    text = c.dna(sequence, "sequence", int(graph.ksize()))
    pairs = graph.get_tags_and_positions(text)
    return {
        "status": "completed",
        "operation": "khmer_tag_positions",
        "graph": info,
        "tags": [{"position": int(pos), "kmer": graph.reverse_hash(int(h))} for pos, h in pairs],
        "provenance": c.provenance(
            source_function="get_tags_and_positions",
            coordinate_system="0-based offset in the query sequence",
        ),
    }


def khmer_reachable_tags(
    kmer: str,
    *,
    assign_partition: bool = False,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    graph_type: str = "nodegraph",
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tags connected to ``kmer`` within khmer's tag range (``find_all_tags_list`` / ``find_all_tags``).

    ``find_all_tags_list`` gives the tag list. With ``assign_partition=True`` the tags found are also
    assigned one partition ID (``find_all_tags`` then ``assign_partition_id``), which is the ID of an
    existing partition when one of the tags already has one and a new ID otherwise. The graph is built and
    tagged from ``reads_path`` and/or ``sequences``.
    """

    assign = c.flag(assign_partition, "assign_partition")
    graph, info = _tagged(
        kind=c.choice(graph_type, c.GRAPH_KINDS, "graph_type"),
        reads_path=reads_path,
        sequences=sequences,
        ksize=ksize,
        table_size=table_size,
        n_tables=n_tables,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    text = c.kmer_of(kmer, int(graph.ksize()))
    tags = c.kmers_of_hashes(graph, sorted(int(h) for h in graph.find_all_tags_list(text)))
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "khmer_reachable_tags",
        "graph": info,
        "kmer": text,
        "n_tags": len(tags),
        "tags": tags,
        "partition_id": None,
        "provenance": c.provenance(
            source_function="find_all_tags_list; find_all_tags; assign_partition_id"
        ),
    }
    if assign:
        result["partition_id"] = int(graph.assign_partition_id(graph.find_all_tags(text)))
    return result


def khmer_unique_paths(
    sequence: str,
    *,
    min_length: int | None = None,
    min_unique_fraction: float = 1.0,
    reads_path: str | Path | Sequence[str | Path] | None = None,
    sequences: Sequence[str] | None = None,
    graph_type: str = "nodegraph",
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pieces of ``sequence`` whose k-mers are mostly absent from the graph (``extract_unique_paths``).

    A path is kept when it is at least ``min_length`` long (default: k) and at least
    ``min_unique_fraction`` of its k-mers are new to the graph built from ``reads_path`` and/or
    ``sequences``. khmer's own documentation marks this method for removal; it is kept because it is
    the khmer operation named "extract unique-coverage paths".
    """

    fraction = c.number(min_unique_fraction, "min_unique_fraction", 0.0, 1.0)
    kind = c.choice(graph_type, c.GRAPH_KINDS, "graph_type")
    k_args = c.table_args(ksize, table_size, n_tables)
    if reads_path is None and sequences is None:
        graph = c.new_graph(require("khmer", c.KHMER_EXTRA), kind, *k_args)
        info: dict[str, Any] = {"source": "empty graph", "ksize": k_args[0]}
    else:
        graph, info = c.make_graph(
            kind=kind,
            graph_path=None,
            reads_path=reads_path,
            sequences=sequences,
            ksize=k_args[0],
            table_size=k_args[1],
            n_tables=k_args[2],
            max_reads=max_reads,
            max_bases=max_bases,
        )
    k = int(graph.ksize())
    length = k if min_length is None else c.integer(min_length, "min_length", k, 10_000_000)
    text = c.dna(sequence, "sequence", k)
    paths = list(graph.extract_unique_paths(text, length, fraction))
    return {
        "status": "completed",
        "operation": "khmer_unique_paths",
        "graph": info,
        "min_length": length,
        "min_unique_fraction": fraction,
        "paths": paths,
        "provenance": c.provenance(source_function="extract_unique_paths"),
    }


def _subset_stats(graph: Any, subset: Any, index: int) -> dict[str, Any]:
    """Statistics of one subset. khmer's ``subset_count_partitions`` damages the subset map, so this runs
    only after the subset was saved and merged."""
    sizes, unassigned = graph.subset_partition_size_distribution(subset)
    return {
        "subset": index,
        "n_partitions": int(graph.subset_count_partitions(subset)[0]),
        "size_distribution": {int(s): int(n) for s, n in sizes},
        "n_unassigned_tags": int(unassigned),
    }


def _partition_graph(
    khmer: Any,
    graph: Any,
    *,
    subset_size: int,
    stop_big_traversals: bool,
    use_stop_tags: bool,
    abundance_range: tuple[int, int] | None,
    save_dir: Path | None,
    prefix: str,
) -> dict[str, Any]:
    """Partition a tagged graph subset by subset and merge the subsets; return per-subset statistics."""
    stats: list[dict[str, Any]] = []
    saved: list[str] = []
    if abundance_range is not None:
        subset = graph.do_subset_partition_with_abundance(abundance_range[0], abundance_range[1])
        if save_dir is not None:
            path = save_dir / f"{prefix}.subset.0.pmap"
            graph.save_subset_partitionmap(subset, str(path))
            saved.append(str(path))
        graph.merge_subset(subset)
        stats.append(_subset_stats(graph, subset, 0))
        return {"n_subsets": 1, "subsets": stats, "saved_partitionmaps": saved}
    bounds = [int(v) for v in graph.divide_tags_into_subsets(subset_size)]
    starts = bounds or [0]
    ends = [*bounds[1:], 0] if bounds else [0]
    if not bounds or len(bounds) == 1:
        starts, ends = [0], [0]
    for index, (start, end) in enumerate(zip(starts, ends, strict=True)):
        subset = graph.do_subset_partition(start, end, use_stop_tags, stop_big_traversals)
        if save_dir is not None:
            path = save_dir / f"{prefix}.subset.{index}.pmap"
            graph.save_subset_partitionmap(subset, str(path))
            saved.append(str(path))
        graph.merge_subset(subset)
        stats.append(_subset_stats(graph, subset, index))
    return {"n_subsets": len(starts), "subsets": stats, "saved_partitionmaps": saved}


def khmer_partition_reads(
    reads_path: str | Path | Sequence[str | Path],
    *,
    output_dir: str | Path | None = None,
    subset_size: int = _DEFAULT_SUBSET_SIZE,
    stop_big_traversals: bool = False,
    stoptags_path: str | Path | None = None,
    abundance_range: Sequence[int] | None = None,
    partitionmap_paths: Sequence[str | Path] | None = None,
    save_partitionmap_dir: str | Path | None = None,
    query_kmers: Sequence[str] | None = None,
    annotate: bool = True,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Partition reads into connected components of their de Bruijn graph and annotate them (``do-partition.py``).

    Steps, as in khmer: load the reads into a nodegraph and tag them (``consume_seqfile_and_tag``), divide
    the tags into subsets of ``subset_size`` (``divide_tags_into_subsets``), partition every subset
    (``do_subset_partition``, the ``partition-graph.py`` step) and merge the subsets (``merge_subset``, the
    ``merge-partitions.py`` step). With ``annotate=True`` every input file is rewritten into
    ``output_dir`` as ``<name>.part`` with the partition ID appended to each read name after a tab
    (``output_partitions``, the ``annotate-partitions.py`` step).

    Options, each one khmer call: ``stoptags_path`` cuts partitions at saved stop tags
    (``load_stop_tags``); ``stop_big_traversals`` truncates graph joins at big traversals;
    ``abundance_range=[min, max]`` partitions only k-mers whose count lies in the range, using a
    countgraph (``do_subset_partition_with_abundance``); ``save_partitionmap_dir`` saves each subset map
    (``save_subset_partitionmap``, files ``<name>.subset.N.pmap``); ``partitionmap_paths`` skips the
    partitioning and merges saved subset maps instead (``merge_subset_from_disk``; the reads are then only
    annotated); ``query_kmers`` returns the partition ID of tag k-mers (``get_partition_id``).
    ``partition_count`` is ``count_partitions()`` of the merged map (partitions, tags without partition; None when
    saved maps were merged, because that graph holds no tags, and the annotated files then give the count);
    ``per_subset`` holds the size distribution of every subset (``subset_partition_size_distribution``).
    """

    files = _paths(reads_path, "reads_path")
    size = c.integer(subset_size, "subset_size", 1, 10**12)
    big = c.flag(stop_big_traversals, "stop_big_traversals")
    annotate_on = c.flag(annotate, "annotate")
    window: tuple[int, int] | None = None
    if abundance_range is not None:
        if len(abundance_range) != 2:
            raise ConfigurationError("abundance_range must be [min_abundance, max_abundance].")
        window = (
            c.integer(abundance_range[0], "abundance_range min", 0, 65535),
            c.integer(abundance_range[1], "abundance_range max", 0, 65535),
        )
        if window[0] > window[1]:
            raise ConfigurationError("abundance_range min must not exceed max.")
    if partitionmap_paths is not None and (window is not None or stoptags_path is not None):
        raise ConfigurationError(
            "partitionmap_paths cannot be combined with abundance_range or stoptags_path."
        )
    out_dir = None
    if annotate_on:
        if output_dir is None:
            raise ConfigurationError("output_dir is required when annotate is True.")
        out_dir = Path(output_dir).expanduser().resolve()
        if not out_dir.is_dir():
            raise ConfigurationError(f"output_dir does not exist: {out_dir}")
        for source in files:
            if (out_dir / f"{source.name}.part").exists():
                raise ConfigurationError(
                    f"output already exists: {out_dir / (source.name + '.part')}"
                )
    save_dir = None
    if save_partitionmap_dir is not None:
        save_dir = Path(save_partitionmap_dir).expanduser().resolve()
        if not save_dir.is_dir():
            raise ConfigurationError(f"save_partitionmap_dir does not exist: {save_dir}")
    k_args = c.table_args(ksize, table_size, n_tables)
    khmer = require("khmer", c.KHMER_EXTRA)
    kind = "countgraph" if window is not None else "nodegraph"
    prefix = files[0].name
    run: dict[str, Any]
    if partitionmap_paths is not None:
        maps = _paths(list(partitionmap_paths), "partitionmap_paths")
        graph = khmer.Nodegraph(k_args[0], 1, 1)
        for pmap in maps:
            try:
                graph.merge_subset_from_disk(str(pmap))
            except (OSError, ValueError, RuntimeError) as exc:
                raise ConfigurationError(
                    f"{pmap.name} is not a khmer partition map: {exc}"
                ) from exc
        info: dict[str, Any] = {"source": "saved partition maps", "ksize": k_args[0]}
        run = {"n_subsets": len(maps), "subsets": [], "saved_partitionmaps": []}
    else:
        graph, info = c.make_graph(
            kind=kind,
            graph_path=None,
            reads_path=[str(p) for p in files],
            sequences=None,
            ksize=k_args[0],
            table_size=k_args[1],
            n_tables=k_args[2],
            tag=True,
            max_reads=max_reads,
            max_bases=max_bases,
        )
        use_stop = False
        if stoptags_path is not None:
            _load_stop_tags(graph, input_path(stoptags_path))
            use_stop = True
        run = _partition_graph(
            khmer,
            graph,
            subset_size=size,
            stop_big_traversals=big,
            use_stop_tags=use_stop,
            abundance_range=window,
            save_dir=save_dir,
            prefix=prefix,
        )
    counted: dict[str, int] | None = None
    if partitionmap_paths is None:
        n_partitions, n_unassigned = (int(v) for v in graph.count_partitions())
        counted = {"n_partitions": n_partitions, "n_unassigned_tags": n_unassigned}
    annotated: list[dict[str, Any]] = []
    if out_dir is not None:
        for source in files:
            target = out_dir / f"{source.name}.part"
            count = int(graph.output_partitions(str(source), str(target)))
            annotated.append({"input": str(source), "output": str(target), "n_partitions": count})
    queries: dict[str, int] = {}
    if query_kmers is not None:
        for kmer in query_kmers:
            text = c.kmer_of(kmer, k_args[0], "query_kmers item")
            queries[text] = int(graph.get_partition_id(text))
    return {
        "status": "completed",
        "operation": "khmer_partition_reads",
        "graph": info,
        "n_subsets": run["n_subsets"],
        "partition_count": counted,
        "per_subset": run["subsets"],
        "saved_partitionmaps": run["saved_partitionmaps"],
        "annotated": annotated,
        "partition_ids": queries,
        "provenance": c.provenance(
            source_function="consume_seqfile_and_tag; divide_tags_into_subsets; do_subset_partition; "
            "merge_subset; output_partitions",
            convention="partition ID 0 means unassigned; IDs are khmer's and are not stable between runs "
            "with different table sizes",
        ),
    }


def khmer_join_partitions(
    reads_path: str | Path | Sequence[str | Path],
    first_tag_kmer: str,
    second_tag_kmer: str,
    *,
    subset_size: int = _DEFAULT_SUBSET_SIZE,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Partition reads, then join the partitions of two tag k-mers (``join_partitions``).

    Both k-mers must be tags (see ``khmer_tag_reads``) of the graph built from ``reads_path``. The result
    gives the partition IDs before and after the join and the partition counts before and after.
    """

    files = _paths(reads_path, "reads_path")
    size = c.integer(subset_size, "subset_size", 1, 10**12)
    k_args = c.table_args(ksize, table_size, n_tables)
    khmer = require("khmer", c.KHMER_EXTRA)
    graph, info = c.make_graph(
        kind="nodegraph",
        graph_path=None,
        reads_path=[str(p) for p in files],
        sequences=None,
        ksize=k_args[0],
        table_size=k_args[1],
        n_tables=k_args[2],
        tag=True,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    first = c.kmer_of(first_tag_kmer, k_args[0], "first_tag_kmer")
    second = c.kmer_of(second_tag_kmer, k_args[0], "second_tag_kmer")
    tag_hashes = {int(graph.hash(str(t))) for t in graph.get_tagset()}
    for name, text in (("first_tag_kmer", first), ("second_tag_kmer", second)):
        if int(graph.hash(text)) not in tag_hashes:
            raise ConfigurationError(f"{name} is not a tag k-mer of this graph.")
    _partition_graph(
        khmer,
        graph,
        subset_size=size,
        stop_big_traversals=False,
        use_stop_tags=False,
        abundance_range=None,
        save_dir=None,
        prefix="",
    )
    before_ids = (int(graph.get_partition_id(first)), int(graph.get_partition_id(second)))
    before_count = [int(v) for v in graph.count_partitions()]
    graph.join_partitions(before_ids[0], before_ids[1])
    after_ids = (int(graph.get_partition_id(first)), int(graph.get_partition_id(second)))
    after_count = [int(v) for v in graph.count_partitions()]
    return {
        "status": "completed",
        "operation": "khmer_join_partitions",
        "graph": info,
        "partition_ids_before": list(before_ids),
        "partition_ids_after": list(after_ids),
        "partition_count_before": before_count,
        "partition_count_after": after_count,
        "provenance": c.provenance(
            source_function="join_partitions; get_partition_id; count_partitions"
        ),
    }


def khmer_load_partitioned_reads(
    partitioned_path: str | Path,
    *,
    query_kmers: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Insert partition-annotated reads (a ``.part`` file) into a new graph and rebuild its partitions.

    Calls ``consume_partitioned_fasta``: every read is counted into the graph and tagged, and its tag
    receives the partition ID written in the read name. ``tag_partition_ids`` maps each tag k-mer of the rebuilt graph to its partition ID; ``query_kmers`` returns
    the partition ID of given k-mers (``get_partition_id``; 0 unless the k-mer is one of those tags). Returns the read and k-mer counts and ``count_partitions()``.
    """

    path = input_path(partitioned_path)
    k_args = c.table_args(ksize, table_size, n_tables)
    khmer = require("khmer", c.KHMER_EXTRA)
    graph = khmer.Nodegraph(*k_args)
    try:
        reads, kmers = graph.consume_partitioned_fasta(str(path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(
            f"{path.name} is not a partition-annotated read file: {exc}"
        ) from exc
    n_partitions, n_unassigned = (int(v) for v in graph.count_partitions())
    ids: dict[str, int] = {}
    for kmer in [] if query_kmers is None else list(query_kmers):
        text = c.kmer_of(kmer, k_args[0], "query_kmers item")
        ids[text] = int(graph.get_partition_id(text))
    return {
        "status": "completed",
        "operation": "khmer_load_partitioned_reads",
        "source": str(path),
        "ksize": k_args[0],
        "n_reads": int(reads),
        "n_kmers_consumed": int(kmers),
        "n_partitions": n_partitions,
        "n_unassigned_tags": n_unassigned,
        "n_tags": int(graph.n_tags()),
        "tag_partition_ids": {
            str(tag): int(graph.get_partition_id(str(tag)))
            for tag in sorted(graph.get_tagset())[: c.MAX_LISTED]
        },
        "partition_ids": ids,
        "provenance": c.provenance(source_function="consume_partitioned_fasta; count_partitions"),
    }


def khmer_find_knots(
    reads_path: str | Path | Sequence[str | Path],
    *,
    stage: str = "knots",
    subset_size: int | None = None,
    excursion_distance: int = _EXCURSION_DISTANCE,
    excursion_kmer_threshold: int = _EXCURSION_KMER_THRESHOLD,
    excursion_kmer_count_threshold: int | None = None,
    load_stoptags_path: str | Path | None = None,
    extra_stop_kmers: Sequence[str] | None = None,
    save_stoptags_path: str | Path | None = None,
    stoptags_text_path: str | Path | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    table_size: int = c.DEFAULT_TABLE_SIZE,
    n_tables: int = c.DEFAULT_N_TABLES,
    counting_table_size: int = c.DEFAULT_TABLE_SIZE,
    max_reads: int = c.DEFAULT_MAX_READS,
    max_bases: int = c.DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find highly connected k-mers (knots) by repartitioning the largest partition and return them as stop tags.

    ``stage="initial"`` follows khmer's ``make-initial-stoptags.py``: partition only the first subset of
    ``subset_size`` tags (default 10000) and repartition its largest partition with an abundance threshold
    of 5. ``stage="knots"`` follows ``find-knots.py``: partition all subsets (default subset size 100000),
    repartition the largest partition of each (threshold 2) and finally of the merged map, which is
    ``repartition_largest_partition``. ``excursion_distance`` (40) and ``excursion_kmer_threshold`` (200) are
    khmer's defaults. Stop-tag handling, each a khmer call: ``load_stoptags_path`` starts from saved stop tags
    (``load_stop_tags``); ``extra_stop_kmers`` adds k-mers (``add_stop_tag``); ``save_stoptags_path`` writes
    the binary file (``save_stop_tags``); ``stoptags_text_path`` writes one k-mer per line
    (``print_stop_tags``).
    """

    mode = c.choice(stage, ("initial", "knots"), "stage")
    files = _paths(reads_path, "reads_path")
    size = c.integer(
        (10_000 if mode == "initial" else _DEFAULT_SUBSET_SIZE)
        if subset_size is None
        else subset_size,
        "subset_size",
        1,
        10**12,
    )
    distance = c.integer(excursion_distance, "excursion_distance", 1, 10**6)
    threshold = c.integer(excursion_kmer_threshold, "excursion_kmer_threshold", 1, 10**9)
    count_threshold = c.integer(
        (5 if mode == "initial" else 2)
        if excursion_kmer_count_threshold is None
        else excursion_kmer_count_threshold,
        "excursion_kmer_count_threshold",
        1,
        65535,
    )
    k_args = c.table_args(ksize, table_size, n_tables)
    counting_size = c.integer(counting_table_size, "counting_table_size", 1, 10**12)
    saved = (
        None
        if save_stoptags_path is None
        else c.new_output(save_stoptags_path, "save_stoptags_path")
    )
    text = (
        None
        if stoptags_text_path is None
        else c.new_output(stoptags_text_path, "stoptags_text_path")
    )
    khmer = require("khmer", c.KHMER_EXTRA)
    graph, info = c.make_graph(
        kind="nodegraph",
        graph_path=None,
        reads_path=[str(p) for p in files],
        sequences=None,
        ksize=k_args[0],
        table_size=k_args[1],
        n_tables=k_args[2],
        tag=True,
        max_reads=max_reads,
        max_bases=max_bases,
    )
    if load_stoptags_path is not None:
        _load_stop_tags(graph, input_path(load_stoptags_path))
    for kmer in [] if extra_stop_kmers is None else list(extra_stop_kmers):
        graph.add_stop_tag(c.kmer_of(kmer, k_args[0], "extra_stop_kmers item"))
    counting = khmer.Countgraph(k_args[0], counting_size, k_args[2])
    bounds = [int(v) for v in graph.divide_tags_into_subsets(size)]
    if mode == "initial":
        start, end = (0, 0) if len(bounds) <= 1 else (bounds[0], bounds[1])
        subset = graph.do_subset_partition(start, end)
        graph.repartition_largest_partition(subset, counting, distance, threshold, count_threshold)
    else:
        starts = [0] if len(bounds) <= 1 else bounds
        ends = [0] if len(bounds) <= 1 else [*bounds[1:], 0]
        for start, end in zip(starts, ends, strict=True):
            subset = graph.do_subset_partition(start, end)
            graph.repartition_largest_partition(
                subset, counting, distance, threshold, count_threshold
            )
            graph.merge_subset(subset)
        graph.repartition_largest_partition(None, counting, distance, threshold, count_threshold)
    stop_tags = sorted(str(t) for t in graph.get_stop_tags())
    if saved is not None:
        graph.save_stop_tags(str(saved))
    if text is not None:
        graph.print_stop_tags(str(text))
    return {
        "status": "completed",
        "operation": "khmer_find_knots",
        "graph": info,
        "stage": mode,
        "subset_size": size,
        "excursion_distance": distance,
        "excursion_kmer_threshold": threshold,
        "excursion_kmer_count_threshold": count_threshold,
        "n_stop_tags": len(stop_tags),
        "stop_tags": stop_tags[: c.MAX_LISTED],
        "saved_stoptags": None if saved is None else str(saved),
        "stoptags_text": None if text is None else str(text),
        "provenance": c.provenance(
            source_function="repartition_largest_partition; save_stop_tags; print_stop_tags; add_stop_tag; load_stop_tags",
            convention="the countgraph used for the excursions has the table size counting_table_size",
        ),
    }


def khmer_filter_stoptags(
    reads_path: str | Path,
    output_path: str | Path,
    *,
    stoptags_path: str | Path | None = None,
    stop_kmers: Sequence[str] | None = None,
    ksize: int = c.DEFAULT_KSIZE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim reads at stop tags and write the kept reads as FASTA (``filter-stoptags.py``, ``trim_on_stoptags``).

    Stop tags come from ``stoptags_path`` (``load_stop_tags``) and/or ``stop_kmers`` (``add_stop_tag``;
    k-mers must have length ``ksize``). As in khmer's script, reads containing N are dropped, a read is cut
    at the first stop tag and kept only if at least ``ksize`` bases remain. ``output_path`` must not exist.
    """

    source = input_path(reads_path)
    target = c.new_output(output_path)
    k = c.ksize_of(ksize)
    if stoptags_path is None and not stop_kmers:
        raise ConfigurationError("give stoptags_path or stop_kmers.")
    khmer = require("khmer", c.KHMER_EXTRA)
    graph = khmer.Nodegraph(k, 1, 1)
    if stoptags_path is not None:
        _load_stop_tags(graph, input_path(stoptags_path))
        if int(graph.ksize()) != k:
            raise ConfigurationError(f"the stop-tag file k-mer size does not match ksize={k}.")
    for kmer in [] if stop_kmers is None else list(stop_kmers):
        graph.add_stop_tag(c.kmer_of(kmer, k, "stop_kmers item"))
    total = kept = trimmed = dropped_n = 0
    try:
        with open(target, "w", encoding="utf-8") as handle:
            for read in khmer.ReadParser(str(source)):
                total += 1
                sequence = read.sequence
                if "N" in sequence:
                    dropped_n += 1
                    continue
                if not set(sequence.upper()) <= set("ACGT"):
                    continue
                cut, position = graph.trim_on_stoptags(sequence)
                if position >= k:
                    handle.write(f">{read.name}\n{cut}\n")
                    kept += 1
                    if position < len(sequence):
                        trimmed += 1
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not parse {source.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "khmer_filter_stoptags",
        "source": str(source),
        "output_path": str(target),
        "ksize": k,
        "n_reads": total,
        "n_kept": kept,
        "n_trimmed": trimmed,
        "n_dropped_with_n": dropped_n,
        "provenance": c.provenance(source_function="trim_on_stoptags"),
    }


def khmer_extract_partitions(
    partitioned_paths: str | Path | Sequence[str | Path],
    output_prefix: str,
    *,
    write_groups: bool = True,
    min_partition_size: int = 5,
    max_group_size: int = 1_000_000,
    write_unassigned: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group partition-annotated reads by partition size and write the size distribution (``extract-partitions.py``).

    Reads come from ``.part`` files written by ``khmer_partition_reads``. The distribution file
    ``<output_prefix>.dist`` has the columns: number of reads, count of partitions with that many reads,
    cumulative partitions, cumulative reads. With ``write_groups=True`` (khmer's default; ``False`` is the
    script's ``-n``) partitions with more than ``min_partition_size`` reads are packed, smallest first,
    into files ``<output_prefix>.groupNNNN.fa`` (or ``.fq``) of at most about ``max_group_size`` reads;
    ``write_unassigned`` also writes reads of partition 0 to ``<output_prefix>.unassigned``. Output files
    must not exist.
    """

    files = _paths(partitioned_paths, "partitioned_paths")
    minimum = c.integer(min_partition_size, "min_partition_size", 0, 10**12)
    limit = c.integer(max_group_size, "max_group_size", 1, 10**12)
    groups_on = c.flag(write_groups, "write_groups")
    unassigned_on = c.flag(write_unassigned, "write_unassigned")
    prefix = Path(output_prefix).expanduser()
    if not prefix.name or not prefix.resolve().parent.is_dir():
        raise ConfigurationError("output_prefix must name a file prefix inside an existing folder.")
    prefix = prefix.resolve()
    khmer = require("khmer", c.KHMER_EXTRA)
    records: list[tuple[Any, int]] = []
    counts: dict[int, int] = {}
    is_fastq: bool | None = None
    for path in files:
        try:
            for read in khmer.ReadParser(str(path)):
                try:
                    pid = int(read.name.rsplit("\t", 1)[1])
                except (IndexError, ValueError) as exc:
                    raise ConfigurationError(
                        f"{path.name} has a read without a tab-separated partition ID: {read.name!r}"
                    ) from exc
                fastq = hasattr(read, "quality") and read.quality is not None
                if is_fastq is None:
                    is_fastq = fastq
                elif is_fastq != fastq:
                    raise ConfigurationError("input files must all be FASTA or all be FASTQ.")
                counts[pid] = counts.get(pid, 0) + 1
                records.append((read, pid))
        except (OSError, ValueError) as exc:
            raise ConfigurationError(f"khmer could not parse {path.name}: {exc}") from exc
    suffix = "fq" if is_fastq else "fa"
    dist_path = Path(f"{prefix}.dist")
    unassigned_path = Path(f"{prefix}.unassigned.{suffix}")
    for target in (dist_path, unassigned_path if unassigned_on else None):
        if target is not None and target.exists():
            raise ConfigurationError(f"output already exists: {target}")
    histogram: dict[int, int] = {}
    for size in counts.values():
        histogram[size] = histogram.get(size, 0) + 1
    rows: list[list[int]] = []
    total = wtotal = 0
    for size, n_parts in sorted(histogram.items()):
        total += n_parts
        wtotal += size * n_parts
        rows.append([size, n_parts, total, wtotal])
    with open(dist_path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(" ".join(str(v) for v in row) + "\n")

    def write(read: Any, handle: Any) -> None:
        if hasattr(read, "quality") and read.quality is not None:
            handle.write(f"@{read.name}\n{read.sequence}\n+\n{read.quality}\n")
        else:
            handle.write(f">{read.name}\n{read.sequence}\n")

    n_unassigned = counts.get(0, 0)
    if unassigned_on:
        with open(unassigned_path, "w", encoding="utf-8") as handle:
            for read, pid in records:
                if pid == 0:
                    write(read, handle)
    group_files: list[str] = []
    group_of: dict[int, int] = {}
    if groups_on:
        ranked = sorted(
            ((pid, n) for pid, n in counts.items() if pid != 0), key=lambda item: item[1]
        )
        ranked = [item for item in ranked if item[1] > minimum]
        running = 0
        pending: set[int] = set()
        group_n = 0
        for pid, n_reads in ranked:
            pending.add(pid)
            running += n_reads
            if running > limit:
                for member in pending:
                    group_of[member] = group_n
                group_n += 1
                pending = set()
                running = 0
        if pending:
            for member in pending:
                group_of[member] = group_n
            group_n += 1
        names = [Path(f"{prefix}.group{index:04d}.{suffix}") for index in range(group_n)]
        for name in names:
            if name.exists():
                raise ConfigurationError(f"output already exists: {name}")
        with contextlib.ExitStack() as stack:
            handles = [stack.enter_context(open(name, "w", encoding="utf-8")) for name in names]
            for read, pid in records:
                if pid != 0 and pid in group_of:
                    write(read, handles[group_of[pid]])
        group_files = [str(name) for name in names]
    return {
        "status": "completed",
        "operation": "khmer_extract_partitions",
        "inputs": [str(p) for p in files],
        "n_reads": len(records),
        "n_partitions": len([pid for pid in counts if pid != 0]),
        "n_unassigned_reads": n_unassigned,
        "size_distribution": [
            {
                "reads": r[0],
                "n_partitions": r[1],
                "cumulative_partitions": r[2],
                "cumulative_reads": r[3],
            }
            for r in rows
        ],
        "distribution_path": str(dist_path),
        "group_files": group_files,
        "unassigned_path": str(unassigned_path) if unassigned_on else None,
        "provenance": c.provenance(
            source_function="extract-partitions.py logic (PartitionExtractor)",
            convention="partition 0 is unassigned and never grouped; a partition is kept when it has more than min_partition_size reads",
        ),
    }


__all__ = [
    "khmer_extract_partitions",
    "khmer_filter_stoptags",
    "khmer_find_knots",
    "khmer_join_partitions",
    "khmer_load_partitioned_reads",
    "khmer_partition_reads",
    "khmer_reachable_tags",
    "khmer_tag_positions",
    "khmer_tag_reads",
    "khmer_unique_paths",
]

_register_public_verbose(globals(), tuple(__all__))
