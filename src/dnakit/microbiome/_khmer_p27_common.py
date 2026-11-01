"""Shared helpers of the khmer P27 functions: argument checks, graph construction and result provenance.

The graph helpers build a khmer table from the same three kinds of input in every function: a saved
table file (``graph_path``), a local FASTA or FASTQ file (``reads_path``) or a list of DNA strings
(``sequences``). Files are consumed with khmer's own ``consume_seqfile*`` methods, as khmer's scripts do,
so khmer's read cleaning (upper-casing, N read as A, reads shorter than k skipped) applies. Strings are
upper-cased and must hold only A, C, G and T.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

KHMER_EXTRA = "bio-khmer"
DEFAULT_KSIZE = 32
DEFAULT_TABLE_SIZE = 1_000_000
DEFAULT_N_TABLES = 4
DEFAULT_MAX_READS = 1_000_000
DEFAULT_MAX_BASES = 100_000_000
MAX_LISTED = 100_000
GRAPH_KINDS = ("nodegraph", "countgraph")
_ACGT = re.compile(r"^[ACGT]+$")
_COUNT_LIMIT = 10**12


def version() -> str:
    try:
        return importlib.metadata.version("khmer")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def provenance(**notes: str) -> dict[str, Any]:
    return {"package": f"khmer {version()}", **notes}


def integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def number(value: object, name: str, minimum: float, maximum: float) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not minimum <= value <= maximum
    ):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    return float(value)


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def dna(value: object, name: str, minimum: int = 1) -> str:
    """Upper-case a DNA string; only A, C, G and T are accepted."""
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip().upper()
    if len(text) < minimum or len(text) > 10_000_000 or not _ACGT.match(text):
        raise ConfigurationError(
            f"{name} must be a DNA string of A, C, G and T with at least {minimum} bases."
        )
    return text


def kmer_of(value: object, ksize: int, name: str = "kmer") -> str:
    text = dna(value, name, ksize)
    if len(text) != ksize:
        raise ConfigurationError(f"{name} has length {len(text)}; the table k-mer size is {ksize}.")
    return text


def sequence_list(values: object, name: str, minimum: int = 1) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of DNA strings.")
    if not 1 <= len(values) <= MAX_LISTED:
        raise ConfigurationError(f"{name} must hold 1 to {MAX_LISTED} sequences.")
    return [dna(item, f"{name} item", minimum) for item in values]


def ksize_of(value: object) -> int:
    return integer(value, "ksize", 1, 32)


def table_args(ksize: object, table_size: object, n_tables: object) -> tuple[int, int, int]:
    return (
        ksize_of(ksize),
        integer(table_size, "table_size", 1, _COUNT_LIMIT),
        integer(n_tables, "n_tables", 1, 20),
    )


def new_output(value: str | Path, name: str = "output_path") -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def scan_reads(khmer: Any, path: Path, max_reads: int, max_bases: int) -> dict[str, int]:
    """Count reads and bases of a FASTA/FASTQ file through khmer's parser, enforcing the limits."""
    reads = 0
    bases = 0
    try:
        for read in khmer.ReadParser(str(path)):
            reads += 1
            bases += len(read.sequence)
            if reads > max_reads:
                raise ConfigurationError(
                    f"{path.name} exceeds max_reads={max_reads}; raise the limit."
                )
            if bases > max_bases:
                raise ConfigurationError(
                    f"{path.name} exceeds max_bases={max_bases}; raise the limit."
                )
    except ConfigurationError:
        raise
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not parse {path.name}: {exc}") from exc
    return {"reads": reads, "bases": bases}


def load_saved(khmer: Any, path: Path, kind: str) -> Any:
    """Load a saved nodegraph or countgraph after checking its header kind."""
    try:
        if kind == "countgraph":
            info = khmer.extract_countgraph_info(str(path))
            if int(info.ht_type) != 1:
                raise ValueError("the file header is not a countgraph")
            return khmer.load_countgraph(str(path))
        info = khmer.extract_nodegraph_info(str(path))
        if int(info[4]) != 2:
            raise ValueError("the file header is not a nodegraph")
        return khmer.load_nodegraph(str(path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer {kind}: {exc}") from exc


def new_graph(khmer: Any, kind: str, ksize: int, table_size: int, n_tables: int) -> Any:
    cls = khmer.Countgraph if kind == "countgraph" else khmer.Nodegraph
    return cls(ksize, table_size, n_tables)


def consume_into(
    khmer: Any,
    graph: Any,
    *,
    reads_paths: Sequence[Path],
    sequences: Sequence[str] | None,
    tag: bool,
    max_reads: int,
    max_bases: int,
) -> dict[str, int]:
    """Consume files (khmer's own methods) and strings into ``graph``; optionally tag the k-mers."""
    totals = {"reads": 0, "bases": 0, "kmers_consumed": 0}
    for path in reads_paths:
        scanned = scan_reads(khmer, path, max_reads, max_bases)
        totals["reads"] += scanned["reads"]
        totals["bases"] += scanned["bases"]
        try:
            result = (
                graph.consume_seqfile_and_tag(str(path))
                if tag
                else graph.consume_seqfile(str(path))
            )
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigurationError(f"khmer could not consume {path.name}: {exc}") from exc
        totals["kmers_consumed"] += int(result[1])
    if sequences is not None:
        ksize = int(graph.ksize())
        for text in sequences:
            if len(text) < ksize:
                continue
            totals["reads"] += 1
            totals["bases"] += len(text)
            if totals["reads"] > max_reads or totals["bases"] > max_bases:
                raise ConfigurationError(
                    "sequence input exceeds max_reads or max_bases; raise the limit."
                )
            if tag:
                result = graph.consume_and_tag(text)
                totals["kmers_consumed"] += (
                    int(result[0]) if isinstance(result, tuple) else int(result)
                )
            else:
                totals["kmers_consumed"] += int(graph.consume(text))
    return totals


def make_graph(
    *,
    kind: str,
    graph_path: str | Path | None,
    reads_path: str | Path | Sequence[str | Path] | None,
    sequences: Sequence[str] | None,
    ksize: int,
    table_size: int,
    n_tables: int,
    tag: bool = False,
    max_reads: int = DEFAULT_MAX_READS,
    max_bases: int = DEFAULT_MAX_BASES,
) -> tuple[Any, dict[str, Any]]:
    """Return ``(graph, info)`` from a saved table or from reads files and strings (or both)."""
    choice(kind, GRAPH_KINDS, "graph_type")
    flag(tag, "tag")
    limits = (
        integer(max_reads, "max_reads", 1, _COUNT_LIMIT),
        integer(max_bases, "max_bases", 1, _COUNT_LIMIT),
    )
    khmer = require("khmer", KHMER_EXTRA)
    paths: list[Path] = []
    if reads_path is not None:
        items = [reads_path] if isinstance(reads_path, (str, Path)) else list(reads_path)
        if not items:
            raise ConfigurationError("reads_path must not be an empty list.")
        paths = [input_path(item) for item in items]
    texts = None if sequences is None else sequence_list(sequences, "sequences")
    if graph_path is None and not paths and texts is None:
        raise ConfigurationError("give one of graph_path, reads_path or sequences.")
    if graph_path is not None and (paths or texts is not None):
        raise ConfigurationError("graph_path cannot be combined with reads_path or sequences.")
    if graph_path is not None:
        if tag:
            raise ConfigurationError(
                "tagging needs reads_path or sequences; a saved table holds no tags."
            )
        path = input_path(graph_path)
        graph = load_saved(khmer, path, kind)
        info: dict[str, Any] = {"source": str(path), "ksize": int(graph.ksize())}
        return graph, info
    graph = new_graph(khmer, kind, *table_args(ksize, table_size, n_tables))
    totals = consume_into(
        khmer,
        graph,
        reads_paths=paths,
        sequences=texts,
        tag=tag,
        max_reads=limits[0],
        max_bases=limits[1],
    )
    info = {
        "source": [str(path) for path in paths] if paths else "sequences",
        "ksize": int(graph.ksize()),
        "table_size": table_size,
        "n_tables": n_tables,
        **totals,
    }
    return graph, info


def kmers_of_hashes(graph: Any, hashes: Sequence[int]) -> list[str]:
    return [graph.reverse_hash(int(value)) for value in hashes]
