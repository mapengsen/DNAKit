"""sourmash plot-data and needs-data operations as typed DNAKit functions.

These functions use sourmash's Python API on local files that the user supplies: a comparison matrix
written by ``sourmash compare``, LCA taxonomy databases, SBT or SQLite indexes, signature files, and
gather/taxonomy CSV files. Nothing is downloaded. Results are JSON-serialisable dictionaries; plots
return the image path together with the plotted data.

Install the optional dependency with ``python -m pip install "dnakit[bio-sourmash]"`` (plots also need
``dnakit[bio-plots]`` and ``dnakit[bio-analysis]``).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-sourmash"
_IMAGE_SUFFIXES = {".png", ".pdf", ".svg"}


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


def _lineage(lineage: Sequence[Any]) -> list[dict[str, str]]:
    """Serialise a sourmash lineage; ranks without a name are placeholders and are dropped."""
    return [{"rank": pair.rank, "name": pair.name} for pair in lineage if pair.name]


def _rank_key(lineage: list[dict[str, str]]) -> list[tuple[str, str]]:
    return [(item["rank"], item["name"]) for item in lineage]


def _lca_summary(path: Path, db: Any, ksize: int, scaled: int) -> dict[str, Any]:
    return {
        "path": str(path),
        "ksize": int(ksize),
        "scaled": int(scaled),
        "signature_count": len(db),
        "hash_count": len(db.hashvals),
    }


def _load_lca(
    lca_databases: Sequence[str], scaled: int | None
) -> tuple[list[Any], int, int, list[Path]]:
    paths = _paths(lca_databases, "lca_databases")
    lca_utils = _sourmash("sourmash.lca.lca_utils")
    dblist, ksize, db_scaled = lca_utils.load_databases(
        [str(p) for p in paths], scaled, verbose=False
    )
    return dblist, ksize, db_scaled, paths


def _signatures(path: Path, ksize: int | None = None) -> list[Any]:
    loader = _sourmash().load_file_as_signatures
    return list(loader(str(path), ksize=ksize))


def _single_signature(path: Path, ksize: int) -> Any:
    signatures = _signatures(path, ksize)
    if len(signatures) != 1:
        raise ConfigurationError(
            f"Expected exactly one signature with k={ksize} in {path}; found {len(signatures)}."
        )
    return signatures[0]


def _at_database_scale(query: Any, scaled: int) -> Any:
    """Downsample a finer query signature to the database scaled value, as sourmash's classify command does."""
    if query.minhash.scaled == scaled:
        return query
    try:
        query.minhash = query.minhash.downsample(scaled=scaled)
    except ValueError as exc:
        raise ConfigurationError(
            f"Signature {query.name!r} has scaled={query.minhash.scaled}, coarser than the database scaled={scaled}."
        ) from exc
    return query


def _index_summary(index: Any, path: Path, kind: str) -> dict[str, Any]:
    names: list[str] = []
    ksizes: set[int] = set()
    scaled: set[int] = set()
    for signature in index.signatures():
        names.append(str(signature.name))
        ksizes.add(int(signature.minhash.ksize))
        scaled.add(int(signature.minhash.scaled))
    return {
        "status": "completed",
        "backend": "sourmash",
        "path": str(path),
        "index_type": kind,
        "signature_count": len(names),
        "signature_names": names,
        "ksizes": sorted(ksizes),
        "scaled_values": sorted(scaled),
    }


def load_comparison_matrix(basefile: str, *, verbose: bool = True) -> dict[str, Any]:
    """Load a comparison matrix and its labels written by ``sourmash compare -o BASE``.

    Reads the ``.npy`` matrix at ``basefile`` and the ``<basefile>.labels.txt`` file next to it
    (sourmash.fig.load_matrix_and_labels). Returns the matrix as nested lists with the labels.
    """
    base = _existing(basefile, "basefile")
    _existing(str(base) + ".labels.txt", "labels file")
    numpy = require("numpy", "bio-analysis")
    matrix, labels = _sourmash("sourmash.fig").load_matrix_and_labels(str(base))
    matrix = numpy.asarray(matrix, dtype=float)
    return {
        "status": "completed",
        "backend": "sourmash",
        "basefile": str(base),
        "size": int(matrix.shape[0]),
        "labels": list(labels),
        "matrix": matrix.tolist(),
    }


def plot_comparison_matrix(
    basefile: str,
    output_image: str,
    *,
    show_labels: bool = True,
    vmin: float = 0.0,
    vmax: float = 1.0,
    force: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a dendrogram with the reordered distance heatmap for a ``sourmash compare`` matrix.

    Uses sourmash.fig.plot_composite_matrix (single-linkage clustering, YlGnBu colours). The image format
    follows the file suffix (.png, .pdf or .svg). Returns the image path, the leaf order and the reordered
    matrix. ``force`` rescales a matrix that is not within [0, 1] instead of raising an error.
    """
    base = _existing(basefile, "basefile")
    _existing(str(base) + ".labels.txt", "labels file")
    out = _new(output_image, "output_image")
    if out.suffix.lower() not in _IMAGE_SUFFIXES:
        raise ConfigurationError("output_image must end with .png, .pdf or .svg.")
    if vmin >= vmax:
        raise ConfigurationError("vmin must be smaller than vmax.")
    require("matplotlib", "bio-plots")
    require("scipy", "bio-analysis")
    pyplot = require("matplotlib.pyplot", "bio-plots")
    fig_module = _sourmash("sourmash.fig")
    matrix, labels = fig_module.load_matrix_and_labels(str(base))
    figure, ordered_labels, reordered = fig_module.plot_composite_matrix(
        matrix.copy(), labels, show_labels=show_labels, vmin=vmin, vmax=vmax, force=force
    )
    try:
        figure.savefig(out, dpi=150, bbox_inches="tight")
    finally:
        pyplot.close(figure)
    return {
        "status": "completed",
        "backend": "sourmash",
        "image": str(out),
        "labels": list(ordered_labels),
        "matrix": reordered.tolist(),
        "vmin": vmin,
        "vmax": vmax,
    }


def load_sqlite_index_summary(
    path: str, *, request_manifest: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Open a sourmash SQLite signature collection and summarise its signatures (load_sqlite_index)."""
    source = _existing(path, "path")
    index = _sourmash("sourmash.index.sqlite_index").load_sqlite_index(
        str(source), request_manifest=request_manifest
    )
    return _index_summary(index, source, type(index).__name__)


def load_sbt_index_summary(
    path: str, *, cache_size: int | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Open a sourmash SBT index file (``.sbt.zip``) and summarise its signatures (sbtmh.load_sbt_index)."""
    source = _existing(path, "path")
    index = _sourmash("sourmash.sbtmh").load_sbt_index(
        str(source), print_version_warning=False, cache_size=cache_size
    )
    return _index_summary(index, source, "SBT")


def load_lca_databases(
    lca_databases: Sequence[str], *, scaled: int | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Load LCA taxonomy databases written by ``sourmash lca index`` and describe each one.

    All databases must share one k-mer size. ``scaled`` (optional) downsamples the databases to a
    coarser resolution while loading (sourmash.lca.lca_utils.load_databases).
    """
    if scaled is not None and scaled < 1:
        raise ConfigurationError("scaled must be a positive integer.")
    dblist, ksize, db_scaled, paths = _load_lca(lca_databases, scaled)
    return {
        "status": "completed",
        "backend": "sourmash",
        "ksize": int(ksize),
        "scaled": int(db_scaled),
        "databases": [
            _lca_summary(p, db, ksize, db_scaled) for p, db in zip(paths, dblist, strict=False)
        ],
    }


def load_lca_database(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Load one LCA taxonomy database and report its k-mer size, scaled value and size (load_single_database)."""
    source = _existing(path, "path")
    db, ksize, scaled = _sourmash("sourmash.lca.lca_db").load_single_database(
        str(source), verbose=False
    )
    return {"status": "completed", "backend": "sourmash", **_lca_summary(source, db, ksize, scaled)}


def gather_lineage_assignments(
    signature: str, lca_databases: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """List the lineages each k-mer hash of a signature is assigned to across LCA databases.

    The signature file must hold exactly one signature with the databases' k-mer size. Uses
    sourmash.lca.lca_utils.gather_assignments. Hash values are returned as strings (JSON keys cannot be
    integers). Lineages are listed in sorted order per hash.
    """
    sig_path = _existing(signature, "signature")
    dblist, ksize, scaled, _ = _load_lca(lca_databases, None)
    sig = _at_database_scale(_single_signature(sig_path, ksize), scaled)
    assignments = _sourmash("sourmash.lca.lca_utils").gather_assignments(sig.minhash.hashes, dblist)
    rows = []
    for hashval in sorted(assignments):
        lineages = sorted(
            (_lineage(lineage) for lineage in assignments[hashval]),
            key=_rank_key,
        )
        rows.append({"hashval": str(hashval), "lineages": lineages})
    return {
        "status": "completed",
        "backend": "sourmash",
        "signature": str(sig_path),
        "ksize": int(ksize),
        "scaled": int(scaled),
        "assigned_hash_count": len(rows),
        "assignments": rows,
    }


def lca_rank_counts(
    lca_databases: Sequence[str], *, min_num: int = 0, verbose: bool = True
) -> dict[str, Any]:
    """Count, for every lineage, the k-mer hashes whose lowest common ancestor it is (sourmash.lca.command_rankinfo)."""
    if min_num < 0:
        raise ConfigurationError("min_num must not be negative.")
    dblist, ksize, scaled, _ = _load_lca(lca_databases, None)
    counts = _sourmash("sourmash.lca.command_rankinfo").make_lca_counts(dblist, min_num=min_num)
    rows = [
        {"lineage": _lineage(lineage), "count": int(count)} for lineage, count in counts.items()
    ]
    rows.sort(key=lambda row: (-row["count"], _rank_key(row["lineage"])))
    return {
        "status": "completed",
        "backend": "sourmash",
        "ksize": int(ksize),
        "scaled": int(scaled),
        "min_num": min_num,
        "counts": rows,
    }


def classify_signature_lineage(
    signature: str,
    lca_databases: Sequence[str],
    *,
    threshold: int = 5,
    majority: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Assign every signature in a file to the lowest common ancestor lineage of its k-mers (sourmash lca classify).

    ``threshold`` is the minimum number of k-mers a lineage needs before it is used (sourmash default 5).
    ``majority`` uses the most frequent lineage instead of the lowest common ancestor. Each result has
    status ``found``, ``disagree`` or ``nomatch``. A finer query is downsampled to the database scaled value.
    """
    if threshold < 0:
        raise ConfigurationError("threshold must not be negative.")
    sig_path = _existing(signature, "signature")
    dblist, ksize, scaled, _ = _load_lca(lca_databases, None)
    classify = _sourmash("sourmash.lca.command_classify")
    results = []
    for query in _signatures(sig_path, ksize):
        query = _at_database_scale(query, scaled)
        lineage, status = classify.classify_signature(query, dblist, threshold, majority)
        results.append({"name": str(query.name), "status": status, "lineage": _lineage(lineage)})
    return {
        "status": "completed",
        "backend": "sourmash",
        "ksize": int(ksize),
        "scaled": int(scaled),
        "threshold": threshold,
        "majority": majority,
        "classifications": results,
    }


def summarize_lineage_counts(
    signature: str,
    lca_databases: Sequence[str],
    *,
    threshold: int = 5,
    ignore_abundance: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a signature's k-mers over the lineage tree, counting each lineage and all its ancestors.

    The signature file must hold exactly one signature with the databases' k-mer size. Uses
    sourmash.lca.command_summarize.summarize. Lineages with fewer than ``threshold`` k-mers are not
    counted (sourmash default 5). ``ignore_abundance`` counts each hash once.
    """
    if threshold < 0:
        raise ConfigurationError("threshold must not be negative.")
    sig_path = _existing(signature, "signature")
    dblist, ksize, scaled, _ = _load_lca(lca_databases, None)
    query = _at_database_scale(_single_signature(sig_path, ksize), scaled)
    counts = _sourmash("sourmash.lca.command_summarize").summarize(
        query.minhash.hashes, dblist, threshold, ignore_abundance
    )
    rows = [
        {"lineage": _lineage(lineage), "count": int(count)} for lineage, count in counts.items()
    ]
    rows.sort(key=lambda row: (-row["count"], _rank_key(row["lineage"])))
    return {
        "status": "completed",
        "backend": "sourmash",
        "signature": str(sig_path),
        "ksize": int(ksize),
        "scaled": int(scaled),
        "threshold": threshold,
        "ignore_abundance": ignore_abundance,
        "counts": rows,
    }


def _taxonomy(taxonomy_files: Sequence[str], *, force: bool = False) -> Any:
    paths = _paths(taxonomy_files, "taxonomy_files")
    multi = _sourmash("sourmash.tax.tax_utils").MultiLineageDB
    return multi.load([str(p) for p in paths], force=force)


def _gather_summary(queries: Sequence[Any]) -> list[dict[str, Any]]:
    rows = []
    for query in queries:
        info = query.query_info
        rows.append(
            {
                "query_name": str(query.query_name),
                "query_md5": str(info.query_md5),
                "query_filename": str(info.query_filename),
                "ksize": int(info.ksize),
                "scaled": int(info.scaled),
                "match_count": len(query.raw_taxresults),
                "missed_identifiers": int(query.n_missed),
                "skipped_identifiers": int(query.n_skipped),
            }
        )
    return rows


def load_gather_results_summary(
    gather_csv: str,
    taxonomy_files: Sequence[str],
    *,
    fail_on_missing_taxonomy: bool = False,
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one ``sourmash gather`` CSV, attach lineages from taxonomy files and summarise each query.

    Uses sourmash.tax.tax_utils.load_gather_results. Taxonomy files are CSV or sqlite lineage files
    (sourmash lineage format). The result lists each query with its number of matches.
    """
    gather = _existing(gather_csv, "gather_csv")
    tax_assign = _taxonomy(taxonomy_files)
    tax_utils = _sourmash("sourmash.tax.tax_utils")
    results, header = tax_utils.load_gather_results(
        str(gather),
        tax_assign,
        fail_on_missing_taxonomy=fail_on_missing_taxonomy,
        keep_full_identifiers=keep_full_identifiers,
        keep_identifier_versions=keep_identifier_versions,
    )
    return {
        "status": "completed",
        "backend": "sourmash",
        "gather_csv": str(gather),
        "columns": list(header),
        "queries": _gather_summary(results.values()),
    }


def check_and_load_gather_csvs_summary(
    gather_csvs: Sequence[str],
    taxonomy_files: Sequence[str],
    *,
    fail_on_missing_taxonomy: bool = False,
    force: bool = False,
    keep_full_identifiers: bool = False,
    keep_identifier_versions: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load several gather CSVs against taxonomy files, checking for empty files and missing identifiers.

    Uses sourmash.tax.tax_utils.check_and_load_gather_csvs. ``force`` skips gather files that cannot be
    read instead of failing. Returns one summary row per query.
    """
    gathers = [str(_existing(path, "gather_csv")) for path in _paths(gather_csvs, "gather_csvs")]
    tax_assign = _taxonomy(taxonomy_files, force=force)
    queries = _sourmash("sourmash.tax.tax_utils").check_and_load_gather_csvs(
        gathers,
        tax_assign,
        fail_on_missing_taxonomy=fail_on_missing_taxonomy,
        force=force,
        keep_full_identifiers=keep_full_identifiers,
        keep_identifier_versions=keep_identifier_versions,
    )
    return {
        "status": "completed",
        "backend": "sourmash",
        "gather_csvs": gathers,
        "queries": _gather_summary(queries),
    }


__all__ = [
    "check_and_load_gather_csvs_summary",
    "classify_signature_lineage",
    "gather_lineage_assignments",
    "lca_rank_counts",
    "load_comparison_matrix",
    "load_gather_results_summary",
    "load_lca_database",
    "load_lca_databases",
    "load_sbt_index_summary",
    "load_sqlite_index_summary",
    "plot_comparison_matrix",
    "summarize_lineage_counts",
]
_register_public_verbose(globals())
