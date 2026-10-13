"""khmer k-mer table query for the table types that ``query_saved_countgraph`` cannot read.

``query_saved_countgraph`` loads a saved countgraph. A saved counttable, nodetable or nodegraph
has a different hash layout, so loading it as a countgraph gives wrong counts without an error.
:func:`query_saved_khmer_table` loads each type with its own khmer loader. khmer is imported only
when the function runs.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_KHMER_EXTRA = "bio-khmer"
_KSIZE_LIMIT = 32
_MAX_QUERY_KMERS = 100_000
_ACGT = frozenset("ACGT")
_TABLE_TYPES = ("counttable", "nodetable", "nodegraph")
_LOADERS = {
    "counttable": "load_counttable",
    "nodetable": "load_nodetable",
    "nodegraph": "load_nodegraph",
}
_HEADER_FLAG = {"counttable": 1, "nodetable": 2, "nodegraph": 2}


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _kmer(value: object) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError("each k-mer must be a non-empty DNA string.")
    text = value.upper()
    if set(text) - _ACGT:
        raise ConfigurationError("each k-mer must contain only A, C, G and T.")
    if len(text) > _KSIZE_LIMIT:
        raise ConfigurationError(f"each k-mer length must be at most {_KSIZE_LIMIT}.")
    return text


def _kmer_queries(value: object) -> list[str]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("kmers must be a non-empty list of DNA strings.")
    if len(value) > _MAX_QUERY_KMERS:
        raise ConfigurationError(f"kmers can contain at most {_MAX_QUERY_KMERS} items.")
    queries = [_kmer(item) for item in value]
    if len({len(item) for item in queries}) != 1:
        raise ConfigurationError("all k-mers must have the same length.")
    return queries


def query_saved_khmer_table(
    source: str | Path,
    kmers: Sequence[str],
    *,
    table_type: str = "counttable",
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up k-mers in a saved khmer counttable, nodetable or nodegraph.

    ``table_type`` must match the class that wrote the file, because khmer stores no marker for it:
    ``"counttable"`` calls ``khmer.load_counttable`` and returns counts; ``"nodetable"`` and
    ``"nodegraph"`` call ``khmer.load_nodetable`` and ``khmer.load_nodegraph`` and return 1 when
    the k-mer was added and 0 otherwise (a hash collision can only turn a 0 into a 1). Use
    ``query_saved_countgraph`` for a countgraph. All query k-mers must have the table's ``ksize``;
    a k-mer and its reverse complement share one entry. The result is sorted by k-mer.
    """

    if not isinstance(table_type, str) or table_type not in _TABLE_TYPES:
        raise ConfigurationError(f"table_type must be one of: {', '.join(_TABLE_TYPES)}.")
    path = input_path(source)
    queries = _kmer_queries(kmers)
    khmer = require("khmer", _KHMER_EXTRA)
    try:
        info = (
            khmer.extract_countgraph_info(str(path))
            if table_type == "counttable"
            else khmer.extract_nodegraph_info(str(path))
        )
        flag = int(info.ht_type) if table_type == "counttable" else int(info[4])
        if flag != _HEADER_FLAG[table_type]:
            raise ValueError("the file header is of another table family")
        table = getattr(khmer, _LOADERS[table_type])(str(path))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer {table_type}: {exc}") from exc
    k = int(table.ksize())
    if len(queries[0]) != k:
        raise ConfigurationError(f"query k-mers have length {len(queries[0])}; table ksize is {k}.")
    values = {kmer: int(table.get(kmer)) for kmer in sorted(set(queries))}
    return {
        "status": "completed",
        "operation": "query_saved_khmer_table",
        "source": str(path),
        "table_type": table_type,
        "ksize": k,
        "hashsizes": [int(value) for value in table.hashsizes()],
        "value_kind": "count" if table_type == "counttable" else "presence",
        "values": values,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "loader": f"khmer.{_LOADERS[table_type]}",
            "counting": "canonical k-mers; a k-mer and its reverse complement share one entry",
            "coordinate_system": "k-mer table; no genomic coordinates",
            "convention": "queries are uppercased; the table type must match the writing class",
        },
    }


__all__ = ["query_saved_khmer_table"]

_register_public_verbose(globals(), tuple(__all__))
