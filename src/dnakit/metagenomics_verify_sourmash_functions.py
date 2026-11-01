"""sourmash LCA taxonomy operations that write files or compare taxonomy tables.

These functions cover the sourmash command-line operations that the existing sourmash facades do not
expose: classifying and summarising several signature files with CSV output (``sourmash lca classify``
and ``sourmash lca summarize``), per-rank LCA counts (``sourmash lca rankinfo``), comparing two lineage
tables (``sourmash lca compare-csv``), and writing the anomaly report of an LCA index build (the
``--report`` file of ``sourmash lca index``). Each function calls the same sourmash functions as the
command it replaces, with explicit arguments instead of an argparse namespace. Nothing is downloaded.

Install the optional dependency with ``python -m pip install "dnakit[bio-sourmash]"``.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-sourmash"

__all__ = [
    "classify_lca_signatures",
    "compare_lineage_assignment_csvs",
    "lca_rank_info",
    "summarize_lca_signatures",
    "write_lca_index_report",
]


def _sourmash(module: str = "sourmash") -> Any:
    return require(module, _EXTRA)


def _existing(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {resolved}")
    return resolved


def _new(path: str | Path, name: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if resolved.exists():
        raise ConfigurationError(f"{name} must not exist yet: {resolved}")
    if not resolved.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {resolved.parent}")
    return resolved


def _paths(values: Sequence[str], name: str) -> list[Path]:
    if isinstance(values, str) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of file paths.")
    return [_existing(value, name) for value in values]


def _scaled_or_none(scaled: int | None) -> int | None:
    if scaled is not None and scaled < 1:
        raise ConfigurationError("scaled must be a positive integer.")
    return scaled


def _lineage(lineage: Sequence[Any]) -> list[dict[str, str]]:
    """Serialise a sourmash lineage; ranks without a name are placeholders and are dropped."""
    return [{"rank": pair.rank, "name": pair.name} for pair in lineage if pair.name]


def _lineage_pairs(rows: Sequence[Mapping[str, str]]) -> tuple[Any, ...]:
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    return tuple(lca_utils.LineagePair(str(row["rank"]), str(row["name"])) for row in rows)


def _at_scale(query: Any, scaled: int) -> Any:
    """Downsample a query signature to the database scaled value, as sourmash lca classify does."""
    if query.minhash.scaled == scaled:
        return query
    try:
        with query.update() as query:
            query.minhash = query.minhash.downsample(scaled=scaled)
    except ValueError as exc:
        raise ConfigurationError(
            f"Signature {query.name!r} has scaled={query.minhash.scaled}, coarser than the database scaled={scaled}."
        ) from exc
    return query


def classify_lca_signatures(
    query_files: Sequence[str],
    lca_databases: Sequence[str],
    output_csv: str,
    *,
    threshold: int = 5,
    majority: bool = False,
    scaled: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify every signature in several files against LCA databases and write the result as CSV.

    Same computation as ``sourmash lca classify`` (sourmash.lca.command_classify.classify_signature): each
    signature is downsampled to the databases' scaled value when needed, then assigned to the lowest common
    ancestor lineage of its k-mers (``majority=True`` uses the most frequent lineage). ``threshold`` is the
    minimum number of k-mers a lineage needs (sourmash default 5). The CSV has the columns ID, status and one
    column per taxonomic rank; status is ``found``, ``disagree`` or ``nomatch``. The returned dictionary holds
    the same rows. ``output_csv`` must not exist yet.
    """
    if threshold < 0:
        raise ConfigurationError("threshold must not be negative.")
    scaled = _scaled_or_none(scaled)
    queries = _paths(query_files, "query_files")
    databases = _paths(lca_databases, "lca_databases")
    out = _new(output_csv, "output_csv")
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    classify = _sourmash("sourmash.lca.command_classify")
    loader = _sourmash("sourmash")
    dblist, ksize, db_scaled = lca_utils.load_databases(
        [str(p) for p in databases], scaled, verbose=False
    )
    rows: list[dict[str, Any]] = []
    with open(out, "w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["ID", "status"] + list(lca_utils.taxlist()))
        for query_path in queries:
            for query in loader.load_file_as_signatures(str(query_path), ksize=ksize):
                query = _at_scale(query, db_scaled)
                lineage, status = classify.classify_signature(query, dblist, threshold, majority)
                writer.writerow([str(query), status] + lca_utils.zip_lineage(lineage))
                rows.append({"id": str(query), "status": status, "lineage": _lineage(lineage)})
    return {
        "status": "completed",
        "backend": "sourmash",
        "output_csv": str(out),
        "ksize": int(ksize),
        "scaled": int(db_scaled),
        "threshold": threshold,
        "majority": majority,
        "classification_count": len(rows),
        "classifications": rows,
    }


def summarize_lca_signatures(
    query_files: Sequence[str],
    lca_databases: Sequence[str],
    *,
    output_csv: str | None = None,
    threshold: int = 5,
    ignore_abundance: bool = False,
    scaled: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise the k-mers of every signature in several files over the LCA taxonomy.

    Same computation as ``sourmash lca summarize`` (sourmash.lca.command_summarize.summarize_main): each
    signature in the query files is counted on its own (load_singletons_and_count), its hashes are assigned
    to lineages with ``summarize`` and each lineage and its ancestors receive the number of hashes.
    Lineages with fewer than ``threshold`` hashes are not counted (sourmash default 5). ``ignore_abundance``
    counts each hash once. When ``output_csv`` is given, the rows are written in the sourmash summarize CSV
    format (columns count, one per rank, filename, sig_name, sig_md5, total_counts); the file must not exist.
    """
    if threshold < 0:
        raise ConfigurationError("threshold must not be negative.")
    scaled = _scaled_or_none(scaled)
    queries = _paths(query_files, "query_files")
    databases = _paths(lca_databases, "lca_databases")
    out = _new(output_csv, "output_csv") if output_csv is not None else None
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    summarizer = _sourmash("sourmash.lca.command_summarize")
    dblist, ksize, db_scaled = lca_utils.load_databases(
        [str(p) for p in databases], scaled, verbose=False
    )
    results: list[dict[str, Any]] = []
    csv_fp = open(out, "w", newline="") if out is not None else None
    write_header = True
    try:
        for filename, sig, hashvals in summarizer.load_singletons_and_count(
            [str(p) for p in queries], ksize, db_scaled, ignore_abundance
        ):
            lineage_counts = summarizer.summarize(hashvals, dblist, threshold, ignore_abundance)
            if not ignore_abundance:
                total = float(sum(hashvals.values()))
            else:
                total = float(len(hashvals))
            if csv_fp is not None:
                summarizer.output_csv(
                    lineage_counts, total, csv_fp, filename, sig, write_header=write_header
                )
                write_header = False
            results.append(
                {
                    "filename": str(filename),
                    "name": str(sig.name),
                    "md5": str(sig.md5sum()),
                    "total": total,
                    "counts": [
                        {"lineage": _lineage(lineage), "count": int(count)}
                        for lineage, count in lineage_counts.items()
                    ],
                }
            )
    finally:
        if csv_fp is not None:
            csv_fp.close()
    return {
        "status": "completed",
        "backend": "sourmash",
        "output_csv": str(out) if out is not None else None,
        "ksize": int(ksize),
        "scaled": int(db_scaled),
        "threshold": threshold,
        "ignore_abundance": ignore_abundance,
        "signature_count": len(results),
        "signatures": results,
    }


def lca_rank_info(
    lca_databases: Sequence[str],
    *,
    min_num: int = 0,
    scaled: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the k-mer hashes whose lowest common ancestor falls on each taxonomic rank.

    Same computation as ``sourmash lca rankinfo`` (sourmash.lca.command_rankinfo.rankinfo_main): the LCA
    counts of every lineage are collected with ``make_lca_counts``, summed by the rank of the lineage's last
    entry, and reported as a count and a percentage of all counted hashes for each rank in taxonomic order.
    ``min_num`` is the minimum number of hashes a lineage needs (sourmash default 0). Percentages are not
    rounded; the sourmash text output rounds them to one decimal.
    """
    if min_num < 0:
        raise ConfigurationError("min_num must not be negative.")
    scaled = _scaled_or_none(scaled)
    databases = _paths(lca_databases, "lca_databases")
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    rankinfo = _sourmash("sourmash.lca.command_rankinfo")
    dblist, ksize, db_scaled = lca_utils.load_databases(
        [str(p) for p in databases], scaled, verbose=False
    )
    counts = rankinfo.make_lca_counts(dblist, min_num)
    counts_by_rank: dict[str, int] = defaultdict(int)
    for lineage, count in counts.items():
        if lineage:
            counts_by_rank[lineage[-1].rank] += count
    total = float(sum(counts_by_rank.values()))
    ranks: list[dict[str, Any]] = []
    if total:
        for rank in lca_utils.taxlist():
            count = counts_by_rank.get(rank, 0)
            ranks.append({"rank": rank, "count": int(count), "percent": count / total * 100.0})
    return {
        "status": "completed",
        "backend": "sourmash",
        "ksize": int(ksize),
        "scaled": int(db_scaled),
        "min_num": min_num,
        "total": int(total),
        "note": None if total else "no hashvals with lineages found",
        "ranks": ranks,
    }


def compare_lineage_assignment_csvs(
    classify_csv: str,
    custom_csv: str,
    *,
    start_column: int = 2,
    tabs: bool = False,
    no_headers: bool = False,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare the lineage assigned to each identifier by two lineage tables.

    Same computation as ``sourmash lca compare-csv`` (sourmash.lca.command_compare_csv.compare_csv). The first
    table is a classify-style CSV (ID, status, then ranks; sourmash's load_taxonomy_assignments with
    start_column 3). The second is a custom table whose lineage starts at 1-based column ``start_column`` (at
    least 2) and is comma-separated, or tab-separated when ``tabs`` is true. For identifiers present in both
    tables with different lineages, the two lineages are compared: ``compatible`` when one lineage is an
    ancestor of the other, ``incompatible`` otherwise, with the lowest common ancestor reported. The result
    lists the identifiers missing from either table, the counts per verdict and the count of incompatible
    lineages for each rank.
    """
    if start_column < 2:
        raise ConfigurationError("start_column must be at least 2.")
    first = _existing(classify_csv, "classify_csv")
    second = _existing(custom_csv, "custom_csv")
    lca_index = _sourmash("sourmash.lca.command_index")
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    assignments0, num_rows0 = lca_index.load_taxonomy_assignments(
        str(first), start_column=3, force=force
    )
    delimiter = "\t" if tabs else ","
    assignments, num_rows = lca_index.load_taxonomy_assignments(
        str(second),
        delimiter=delimiter,
        start_column=start_column,
        use_headers=not no_headers,
        force=force,
    )
    missing_in_custom = sorted(set(assignments0) - set(assignments))
    missing_in_classify = sorted(set(assignments) - set(assignments0))
    common = sorted(set(assignments0) & set(assignments))
    differences: list[dict[str, Any]] = []
    incompatible_by_rank: dict[str, int] = defaultdict(int)
    for identifier in common:
        v0 = assignments0[identifier]
        v1 = assignments[identifier]
        if v0 == v1:
            continue
        tree = lca_utils.build_tree([v0])
        lca_utils.build_tree([v1], tree)
        lca, reason = lca_utils.find_lca(tree)
        verdict = "compatible" if reason == 0 else "incompatible"
        if verdict == "incompatible":
            rank = next(iter(lca_utils.taxlist()))
            if lca:
                rank = lca[-1].rank
            incompatible_by_rank[rank] += 1
        differences.append({"identifier": str(identifier), "verdict": verdict, "lca": _lineage(lca)})
    n_incompatible = sum(1 for row in differences if row["verdict"] == "incompatible")
    return {
        "status": "completed",
        "backend": "sourmash",
        "classify_csv": str(first),
        "custom_csv": str(second),
        "classify_rows": int(num_rows0),
        "custom_rows": int(num_rows),
        "classify_lineage_count": len(set(assignments0.values())),
        "custom_lineage_count": len(set(assignments.values())),
        "missing_in_custom": [str(x) for x in missing_in_custom],
        "missing_in_classify": [str(x) for x in missing_in_classify],
        "common_count": len(common),
        "different_count": len(differences),
        "compatible_count": len(differences) - n_incompatible,
        "incompatible_count": n_incompatible,
        "incompatible_by_rank": [
            {"rank": rank, "count": int(incompatible_by_rank.get(rank, 0))}
            for rank in lca_utils.taxlist()
        ],
        "differences": differences,
    }


def write_lca_index_report(
    duplicate_signatures: Sequence[str],
    identifiers_without_lineage: Sequence[str],
    identifiers_without_signature: Sequence[str],
    unused_identifiers: Sequence[str],
    unused_lineages: Sequence[Sequence[Mapping[str, str]]],
    output_report: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the anomaly report of an LCA index build as a text file.

    Same output as the ``--report`` file of ``sourmash lca index`` (sourmash.lca.command_index.generate_report):
    the number and the entries of duplicate signatures, identifiers with no lineage, identifiers with no
    signature, unused identifiers and unused lineages. Each lineage is a list of ``{"rank", "name"}`` mappings
    from the most general rank to the most specific. ``output_report`` must not exist yet.
    """
    out = _new(output_report, "output_report")
    lca_index = _sourmash("sourmash.lca.command_index")
    unused_pairs = [_lineage_pairs(lineage) for lineage in unused_lineages]
    lca_index.generate_report(
        [str(x) for x in duplicate_signatures],
        [str(x) for x in identifiers_without_lineage],
        [str(x) for x in identifiers_without_signature],
        unused_pairs,
        [str(x) for x in unused_identifiers],
        str(out),
    )
    return {
        "status": "completed",
        "backend": "sourmash",
        "report": str(out),
        "duplicate_count": len(duplicate_signatures),
        "no_lineage_count": len(identifiers_without_lineage),
        "no_signature_count": len(identifiers_without_signature),
        "unused_identifier_count": len(unused_identifiers),
        "unused_lineage_count": len(unused_lineages),
    }


_register_public_verbose(globals())
