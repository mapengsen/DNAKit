"""Biopython motif and gene-data file functions of integration pack P18 (regulation area).

Motif functions read the output of motif finders and motif tables that ``Bio.motifs`` parses
(AlignACE, MEME XML, Cluster Buster, PFM tables, XMS) and compute the GC content of motifs.
Data-file functions read GEO SOFT files (``Bio.Geo``), Affymetrix CEL files
(``Bio.Affy.CelFile``), UniGene records (``Bio.UniGene``) and Cluster/TreeView expression
tables (``Bio.Cluster``), and save hierarchical or k-means clustering results as CDT, GTR, ATR
and KGG files.

Biopython is imported only when a function runs. Readers read local files and never touch the
network; the clustering function writes new files and never overwrites.
"""

from __future__ import annotations

import glob
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _biopython_p18_common as _c
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = _c.EXTRA
_MOTIF_FORMATS = (
    "alignace",
    "meme",
    "clusterbuster",
    "pfm-four-columns",
    "pfm-four-rows",
    "xms",
)
_LINKAGE = {"single": "s", "complete": "m", "centroid": "c", "average": "a"}
_DISTANCES = {
    "euclidean": "e",
    "city_block": "b",
    "pearson": "c",
    "absolute_pearson": "a",
    "uncentered": "u",
    "absolute_uncentered": "x",
    "spearman": "s",
    "kendall": "k",
}
_CENTERS = {"mean": "a", "median": "m"}
_DNA_SITE = frozenset("ACGTU")


def _optional_int(value: object) -> int | None:
    return None if value is None else int(value)  # type: ignore[call-overload]


def _object(value: object) -> Any:
    """Turn a Biopython record object into a plain dictionary of its public attributes."""

    if hasattr(value, "__dict__") and not isinstance(value, (str, bytes)):
        return {key: _object(item) for key, item in vars(value).items() if not key.startswith("_")}
    if isinstance(value, (list, tuple)):
        return [_object(item) for item in value[:100_000]]
    return _c.plain(value)


# ---------------------------------------------------------------- motifs


def _motif_summary(motif: Any, with_instances: bool, instance_limit: int) -> dict[str, Any]:
    counts = {str(letter): [float(v) for v in motif.counts[letter]] for letter in motif.alphabet}
    data: dict[str, Any] = {
        "name": str(motif.name),
        "length": int(motif.length),
        "alphabet": str(motif.alphabet),
        "consensus": str(motif.consensus),
        "counts": counts,
        "gc_content": float(motif.counts.gc_content),
    }
    for attribute in ("id", "alt_id", "evalue", "num_occurrences", "score"):
        if hasattr(motif, attribute):
            data[attribute] = _c.plain(getattr(motif, attribute))
    mask = getattr(motif, "mask", None)
    if mask is not None:
        data["mask"] = [int(v) for v in mask]
    sites = getattr(motif, "alignment", None)
    if sites is not None:
        data["instance_count"] = len(sites.sequences)
        if with_instances:
            data["instances"] = [
                {
                    "sequence": str(instance),
                    **{
                        key: _c.plain(value)
                        for key, value in vars(instance).items()
                        if not key.startswith("_")
                    },
                }
                for instance in sites.sequences[:instance_limit]
            ]
    return data


def read_motif_collection(
    path: str | Path,
    *,
    file_format: str = "meme",
    include_instances: bool = False,
    max_motifs: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every motif of an AlignACE, MEME, Cluster Buster, PFM or XMS file.

    ``file_format`` is ``alignace`` (AlignACE output), ``meme`` (MEME XML output),
    ``clusterbuster``, ``pfm-four-columns`` or ``pfm-four-rows`` (frequency tables with one
    row or one column per base) or ``xms``. Each motif has its name, length, consensus, base
    counts per position, GC content of the counts and, where the format has them, the
    identifier, E-value, number of occurrences, AlignACE score and mask. With
    ``include_instances`` the binding-site instances (sequence, and for MEME the sequence
    name, start, strand and p-value) are listed, at most 1000 per motif.
    """

    source = _c.source_file(path)
    fmt = _c.choice(file_format, _MOTIF_FORMATS, "file_format")
    with_instances = _c.flag(include_instances, "include_instances")
    limit = _c.integer(max_motifs, "max_motifs", minimum=1, maximum=100_000)
    motifs = require("Bio.motifs", _EXTRA)
    found: list[dict[str, Any]] = []
    try:
        with source.open(encoding="utf-8") as handle:
            for motif in motifs.parse(handle, fmt):
                if len(found) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_motifs={limit}.")
                found.append(_motif_summary(motif, with_instances, 1000))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid {fmt} motif file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_motif_collection",
        "source": str(source),
        "file_format": fmt,
        "motif_count": len(found),
        "motifs": found,
        "provenance": _c.provenance(
            parser="Bio.motifs.parse",
            coordinates="instance starts are 1-based as written by the motif program",
        ),
    }


def motif_gc_content(
    motifs: Mapping[str, Sequence[str] | Mapping[str, Sequence[float]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the GC fraction of motifs from their sites or their base-count matrices.

    ``motifs`` maps a motif name either to a list of equal-length DNA binding sites or to a
    count matrix ``{"A": [...], "C": [...], "G": [...], "T": [...]}``. The fraction is the
    total count of C and G over the total count of all bases across all positions
    (``GenericPositionMatrix.gc_content``). Position-specific scoring matrices have no GC
    content in Biopython, so log-odds matrices are not accepted.
    """

    if not isinstance(motifs, Mapping) or not motifs:
        raise ConfigurationError("motifs must be a non-empty mapping.")
    if len(motifs) > 10_000:
        raise ConfigurationError("motifs may hold at most 10000 entries.")
    bio_motifs = require("Bio.motifs", _EXTRA)
    seq_module = require("Bio.Seq", _EXTRA)
    results = []
    for name, value in motifs.items():
        _c.identifier(name, "motif name")
        if isinstance(value, Mapping):
            letters = sorted(str(key) for key in value)
            if letters not in (["A", "C", "G", "T"], ["A", "C", "G", "U"]):
                raise ConfigurationError(f"motif {name!r} counts need the four bases A, C, G, T.")
            lengths = {len(value[key]) for key in value}
            if len(lengths) != 1 or 0 in lengths:
                raise ConfigurationError(f"motif {name!r} counts need rows of equal length.")
            matrix = bio_motifs.matrix.FrequencyPositionMatrix(
                "".join(letters), {key: [float(v) for v in value[key]] for key in value}
            )
            length = lengths.pop()
            source_kind = "counts"
        else:
            if isinstance(value, str) or not isinstance(value, Sequence) or not value:
                raise ConfigurationError(f"motif {name!r} needs a non-empty list of sites.")
            sites = []
            for site in value:
                if not isinstance(site, str) or not site or set(site.upper()) - _DNA_SITE:
                    raise ConfigurationError(f"motif {name!r} has a site that is not A, C, G, T.")
                sites.append(seq_module.Seq(site.upper()))
            if len({len(site) for site in sites}) != 1:
                raise ConfigurationError(f"motif {name!r} sites must share one length.")
            motif = bio_motifs.create(sites)
            matrix = motif.counts
            length = len(sites[0])
            source_kind = "sites"
        results.append(
            {
                "name": name,
                "input": source_kind,
                "length": length,
                "gc_content": float(matrix.gc_content),
            }
        )
    return {
        "status": "completed",
        "operation": "motif_gc_content",
        "motif_count": len(results),
        "motifs": results,
        "provenance": _c.provenance(method="Bio.motifs.matrix.GenericPositionMatrix.gc_content"),
    }


# ---------------------------------------------------------------- gene-data files


def read_geo_soft(
    path: str | Path,
    *,
    max_records: int = 1_000,
    max_table_rows: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the entities (platform, sample, series, dataset) of a GEO SOFT file.

    Each record has its entity type (``SAMPLE``, ``SERIES``, ``PLATFORM`` ...), identifier,
    entity attributes (``!Sample_title`` and so on, without the leading marker), the column
    definitions and the table: the header row and up to ``max_table_rows`` data rows, with the
    total row count. ``Bio.Geo.parse`` is used; text that is not UTF-8 is read with
    replacement characters.
    """

    source = _c.source_file(path)
    record_limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    row_limit = _c.integer(max_table_rows, "max_table_rows", minimum=0, maximum=1_000_000)
    geo = require("Bio.Geo", _EXTRA)
    found: list[dict[str, Any]] = []
    try:
        with source.open(encoding="utf-8", errors="replace") as handle:
            for record in geo.parse(handle):
                if len(found) >= record_limit:
                    raise ConfigurationError(f"{source.name} exceeds max_records={record_limit}.")
                rows = record.table_rows
                found.append(
                    {
                        "entity_type": str(record.entity_type),
                        "entity_id": str(record.entity_id),
                        "entity_attributes": _c.plain(record.entity_attributes),
                        "col_defs": _c.plain(record.col_defs),
                        "table_header": _c.plain(rows[0]) if rows else [],
                        "table_row_count": max(len(rows) - 1, 0),
                        "table_rows": _c.plain([list(row) for row in rows[1 : 1 + row_limit]]),
                    }
                )
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid GEO SOFT file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_geo_soft",
        "source": str(source),
        "record_count": len(found),
        "records": found,
        "provenance": _c.provenance(parser="Bio.Geo.parse"),
    }


def read_affymetrix_cel(
    path: str | Path,
    *,
    include_values: bool = False,
    max_cells: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read an Affymetrix CEL file (version 3 text or version 4 binary) of probe intensities.

    The format version is detected from the file. Returned are the array size, the algorithm
    and its parameters, the grid corners, the DAT header, the numbers of masked and outlier
    cells and summary statistics (minimum, maximum, mean) of the intensities. With
    ``include_values`` the intensity, standard-deviation and pixel-count matrices are also
    returned as lists of rows, at most ``max_cells`` cells each.
    """

    source = _c.source_file(path)
    with_values = _c.flag(include_values, "include_values")
    cell_limit = _c.integer(max_cells, "max_cells", minimum=1, maximum=100_000_000)
    cel = require("Bio.Affy.CelFile", _EXTRA)
    with source.open("rb") as handle:
        binary = handle.read(1) == b"\x40"
    try:
        with source.open("rb") if binary else source.open("r", encoding="utf-8") as handle:
            record = cel.read(handle)
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid CEL file: {exc}") from exc
    values = record.intensities
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "read_affymetrix_cel",
        "source": str(source),
        "version": int(record.version),
        "rows": int(record.nrows),
        "columns": int(record.ncols),
        "cell_count": int(values.size),
        "algorithm": _c.plain(record.Algorithm),
        "algorithm_parameters": _c.plain(record.AlgorithmParameters),
        "dat_header": _c.plain(record.DatHeader),
        "grid_corners": {
            name: _c.plain(getattr(record, name))
            for name in ("GridCornerUL", "GridCornerUR", "GridCornerLR", "GridCornerLL")
        },
        "masked_cell_count": _optional_int(record.nmask),
        "outlier_count": _optional_int(record.noutliers),
        "intensity_minimum": float(values.min()),
        "intensity_maximum": float(values.max()),
        "intensity_mean": float(values.mean()),
        "provenance": _c.provenance(parser="Bio.Affy.CelFile.read"),
    }
    if with_values:
        if values.size > cell_limit:
            raise ConfigurationError(
                f"{source.name} has {values.size} cells; max_cells={cell_limit}."
            )
        result["intensities"] = values.tolist()
        result["standard_deviations"] = record.stdevs.tolist()
        result["pixel_counts"] = record.npix.tolist()
    return result


def read_unigene(
    path: str | Path, *, max_records: int = 1_000, verbose: bool = True
) -> dict[str, Any]:
    """Read the clusters of a UniGene ``.data`` flat file.

    Each record has the cluster identifier (for example ``Hs.2``), species, title, gene
    symbol, cytoband, chromosome, expression sources, gene and LocusLink identifiers, and the
    parsed protein-similarity (``protsim``), sequence, STS and transcript-map lines.
    """

    source = _c.source_file(path)
    limit = _c.integer(max_records, "max_records", minimum=1, maximum=_c.MAX_RECORDS)
    unigene = require("Bio.UniGene", _EXTRA)
    found: list[dict[str, Any]] = []
    try:
        with source.open(encoding="utf-8", errors="replace") as handle:
            for record in unigene.parse(handle):
                if len(found) >= limit:
                    raise ConfigurationError(f"{source.name} exceeds max_records={limit}.")
                found.append(_object(record))
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"{source.name} is not a valid UniGene file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_unigene",
        "source": str(source),
        "record_count": len(found),
        "records": found,
        "provenance": _c.provenance(parser="Bio.UniGene.parse"),
    }


def _expression_record(path: Path) -> Any:
    cluster = require("Bio.Cluster", _EXTRA)
    try:
        with path.open(encoding="utf-8") as handle:
            return cluster.read(handle)
    except Exception as exc:
        raise ConfigurationError(f"{path.name} is not a Cluster/TreeView table: {exc}") from exc


def _record_summary(record: Any, include_data: bool) -> dict[str, Any]:
    rows, columns = record.data.shape
    summary: dict[str, Any] = {
        "gene_count": int(rows),
        "experiment_count": int(columns),
        "gene_ids": _c.plain(record.geneid),
        "gene_names": _c.plain(record.genename),
        "experiment_ids": _c.plain(record.expid),
        "gene_weights": _c.plain(record.gweight),
        "experiment_weights": _c.plain(record.eweight),
        "gene_order": _c.plain(record.gorder),
        "experiment_order": _c.plain(record.eorder),
        "id_column_label": _c.plain(record.uniqid),
        "has_missing_values": record.mask is not None and bool((record.mask == 0).any()),
    }
    if include_data:
        summary["data"] = record.data.tolist()
        summary["mask"] = None if record.mask is None else record.mask.tolist()
    return summary


def read_cluster_expression_data(
    path: str | Path, *, include_data: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Read a Cluster/TreeView tab-delimited expression table (genes in rows, samples in columns).

    The first column holds the gene identifiers, optionally followed by a name column, a
    ``GWEIGHT`` column, the sample columns, and optional ``EWEIGHT``, ``GORDER`` and ``EORDER``
    rows (``Bio.Cluster.read``). Returned are the table size, identifiers, names, weights,
    orders and whether values are missing; with ``include_data`` the numeric matrix and the
    mask (1 = value present) are listed too.
    """

    source = _c.source_file(path)
    with_data = _c.flag(include_data, "include_data")
    record = _expression_record(source)
    return {
        "status": "completed",
        "operation": "read_cluster_expression_data",
        "source": str(source),
        **_record_summary(record, with_data),
        "provenance": _c.provenance(parser="Bio.Cluster.read"),
    }


def cluster_expression_data(
    path: str | Path,
    output_prefix: str | Path,
    *,
    algorithm: str = "hierarchical",
    linkage: str = "complete",
    distance: str = "euclidean",
    cluster_samples: bool = False,
    cluster_count: int = 2,
    passes: int = 1,
    center: str = "mean",
    initial_clusters: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cluster an expression table and save the result for Java TreeView (CDT, GTR, ATR, KGG).

    ``algorithm="hierarchical"`` clusters the genes with ``linkage`` ``single``, ``complete``
    (the Biopython default), ``centroid`` or ``average``; with ``cluster_samples`` the samples
    are clustered too. It writes ``<prefix>.cdt``, ``<prefix>.gtr`` and, for samples,
    ``<prefix>.atr``. ``algorithm="kmeans"`` runs k-means (``center="mean"``) or k-medians
    (``center="median"``) with ``cluster_count`` clusters and ``passes`` random restarts and
    writes ``<prefix>_K_G<k>.cdt`` and ``.kgg``; give ``initial_clusters`` (one cluster
    number per gene) for a deterministic run. ``distance`` is euclidean, city_block, pearson,
    absolute_pearson, uncentered, absolute_uncentered, spearman or kendall. The function stops
    if any file starting with the output prefix already exists.
    """

    source = _c.source_file(path)
    kind = _c.choice(algorithm, ("hierarchical", "kmeans"), "algorithm")
    method = _LINKAGE[_c.choice(linkage, tuple(_LINKAGE), "linkage")]
    metric = _DISTANCES[_c.choice(distance, tuple(_DISTANCES), "distance")]
    both = _c.flag(cluster_samples, "cluster_samples")
    groups = _c.integer(cluster_count, "cluster_count", minimum=2, maximum=10_000)
    restarts = _c.integer(passes, "passes", minimum=1, maximum=10_000)
    middle = _CENTERS[_c.choice(center, tuple(_CENTERS), "center")]
    if not isinstance(output_prefix, (str, Path)) or not str(output_prefix).strip():
        raise ConfigurationError("output_prefix must be a file path prefix.")
    prefix = Path(output_prefix).expanduser().resolve()
    if not prefix.parent.is_dir():
        raise ConfigurationError(f"output_prefix directory does not exist: {prefix.parent}")
    if glob.glob(glob.escape(str(prefix)) + "*"):
        raise ConfigurationError(f"Files starting with {prefix} already exist.")
    record = _expression_record(source)
    rows = record.data.shape[0]
    initial: Any = None
    if initial_clusters is not None:
        supplied: object = initial_clusters
        if isinstance(supplied, str) or not isinstance(supplied, Sequence):
            raise ConfigurationError("initial_clusters must be a list of integers.")
        if len(initial_clusters) != rows:
            raise ConfigurationError(f"initial_clusters needs {rows} values, one per gene.")
        initial = [
            _c.integer(v, "initial_clusters", minimum=0, maximum=groups - 1)
            for v in initial_clusters
        ]
        if len(set(initial)) != groups:
            raise ConfigurationError(
                "initial_clusters must use every cluster number 0..cluster_count-1."
            )
    extra: dict[str, Any] = {}
    try:
        if kind == "hierarchical":
            gene_tree = record.treecluster(method=method, dist=metric)
            sample_tree = (
                record.treecluster(transpose=True, method=method, dist=metric) if both else None
            )
            record.save(str(prefix), geneclusters=gene_tree, expclusters=sample_tree)
            extra = {
                "gene_tree_nodes": len(gene_tree),
                "sample_tree_nodes": None if sample_tree is None else len(sample_tree),
            }
        else:
            if rows < groups:
                raise ConfigurationError("cluster_count exceeds the number of genes.")
            labels, error, found = record.kcluster(
                nclusters=groups, npass=restarts, method=middle, dist=metric, initialid=initial
            )
            record.save(str(prefix), geneclusters=labels)
            extra = {
                "cluster_labels": [int(v) for v in labels],
                "within_cluster_error": float(error),
                "solutions_found": int(found),
            }
    except ConfigurationError:
        for leftover in glob.glob(glob.escape(str(prefix)) + "*"):
            Path(leftover).unlink(missing_ok=True)
        raise
    except Exception as exc:
        for leftover in glob.glob(glob.escape(str(prefix)) + "*"):
            Path(leftover).unlink(missing_ok=True)
        raise ConfigurationError(f"Clustering failed: {exc}") from exc
    files = sorted(glob.glob(glob.escape(str(prefix)) + "*"))
    return {
        "status": "completed",
        "operation": "cluster_expression_data",
        "source": str(source),
        "algorithm": kind,
        "gene_count": int(rows),
        "output_files": files,
        **extra,
        "provenance": _c.provenance(
            method="Bio.Cluster.Record.treecluster"
            if kind == "hierarchical"
            else "Bio.Cluster.Record.kcluster",
            writer="Bio.Cluster.Record.save",
        ),
    }


__all__ = [
    "cluster_expression_data",
    "motif_gc_content",
    "read_affymetrix_cel",
    "read_cluster_expression_data",
    "read_geo_soft",
    "read_motif_collection",
    "read_unigene",
]

_register_public_verbose(globals(), tuple(__all__))
