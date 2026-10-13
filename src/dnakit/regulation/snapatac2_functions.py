"""Fragment, peak, gene and motif functions over snapATAC2 for five-column fragment files.

Fragment inputs use the same format as ``read_fragment_qc()``: plain text with five tab-separated
columns (chromosome, start, end, barcode, read count), 0-based half-open coordinates. Each function
imports the fragments into a temporary snapATAC2 file that is removed on return. Gene annotations
are local GTF files, genomes are local FASTA files, and no data is downloaded.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _DEFAULT_MAX_CELLS,
    _DEFAULT_MAX_FRAGMENTS,
    _DEFAULT_MAX_PEAKS,
    _SNAPATAC2_EXTRA,
    _as_list,
    _chrom_sizes,
    _count_entries,
    _flag,
    _imported_fragments,
    _integer,
    _label,
    _peak_file,
    _snap_provenance,
    _upstream,
)

_MAX_REGIONS = 100_000
_MAX_GROUPS = 1_000
_MAX_MOTIFS = 1_000
_MAX_ENTRIES = 5_000_000
_MAX_SIZE = 10_000
_REGION = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")
_PEAK_COLUMNS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signal_value",
    "p_value",
    "q_value",
    "summit",
)


def _fragments(
    fragment_file: object,
    chrom_sizes: object,
    min_num_fragments: int,
    max_fragments: int,
    max_cells: int,
) -> tuple[Any, Any]:
    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    return snap, _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=_integer(min_num_fragments, "min_num_fragments", 1, 1_000_000),
        max_fragments=_integer(max_fragments, "max_fragments", 1, 1_000_000_000),
        max_cells=_integer(max_cells, "max_cells", 1, 10_000_000),
    )


def _regions(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of chr:start-end strings.")
    if len(values) > _MAX_REGIONS:
        raise ConfigurationError(f"{name} has more than {_MAX_REGIONS} regions.")
    checked = []
    for value in values:
        if not isinstance(value, str):
            raise ConfigurationError(f"{name} entries must be chr:start-end strings.")
        match = _REGION.match(value.strip())
        if not match or int(match.group(2)) >= int(match.group(3)):
            raise ConfigurationError(
                f"region {value!r} must look like chr1:100-200 with start < end."
            )
        checked.append(value.strip())
    return checked


def _gtf(path: object) -> Path:
    gtf = input_path(path)  # type: ignore[arg-type]
    if gtf.suffix.lower() not in {".gtf", ".gff", ".gff3", ".gz"}:
        raise ConfigurationError("gene annotation must be a GTF or GFF file.")
    return gtf


def _scalars(value: object) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and value != value:
        return None
    return value


def tss_enrichment(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    gene_gtf: str | Path,
    *,
    exclude_chroms: Sequence[str] = ("chrM", "M"),
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Transcription start site enrichment per cell, with snapATAC2 ``metrics.tsse``.

    TSS positions come from the transcripts in ``gene_gtf``. ``tsse`` is the cell's fragment
    density near TSSs relative to its flanks; ``library_tsse`` and ``frac_overlap_TSS`` summarise
    the whole library. Chromosomes in ``exclude_chroms`` are ignored.
    """

    chroms = [_label(name, "exclude_chroms entry") for name in exclude_chroms]
    gtf = _gtf(gene_gtf)
    snap, context = _fragments(
        fragment_file, chrom_sizes, min_num_fragments, max_fragments, max_cells
    )
    with context as (adata, _work):
        values = _upstream(
            lambda: snap.metrics.tsse(
                adata,
                gtf,
                exclude_chroms=chroms,
                inplace=False,
                n_jobs=1,
            )
        )
        barcodes = [str(item) for item in adata.obs_names]
    tsse = [_scalars(item) for item in _as_list(values["tsse"])]
    cells = [{"barcode": barcode, "tsse": tsse[index]} for index, barcode in enumerate(barcodes)]
    return {
        "status": "completed",
        "operation": "tss_enrichment",
        "cell_count": len(cells),
        "library_tsse": _scalars(values["library_tsse"]),
        "frac_overlap_TSS": _scalars(values["frac_overlap_TSS"]),
        "cells": cells,
        "provenance": _snap_provenance(
            tss="transcript starts from the GTF file, as snapATAC2 reads them",
            excluded_chromosomes=chroms,
        ),
    }


def fragment_in_peaks_fraction(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    peak_file: str | Path,
    *,
    normalized: bool = True,
    count_as_insertion: bool = False,
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_peaks: int = _DEFAULT_MAX_PEAKS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fraction of fragments per cell in a set of peaks (FRiP), with snapATAC2 ``metrics.frip``.

    ``peak_file`` is a BED file whose first three columns give the peaks. With ``normalized=True``
    the count is divided by each cell's total fragments. ``count_as_insertion=True`` counts Tn5
    insertion sites instead of fragments.
    """

    snap, context = _fragments(
        fragment_file, chrom_sizes, min_num_fragments, max_fragments, max_cells
    )
    sizes = _chrom_sizes(chrom_sizes)
    peaks = _peak_file(peak_file, sizes, _integer(max_peaks, "max_peaks", 1, 10_000_000))
    flag = _flag(normalized, "normalized")
    insertion = _flag(count_as_insertion, "count_as_insertion")
    with context as (adata, _work):
        values = _upstream(
            lambda: snap.metrics.frip(
                adata,
                {"peaks": peaks},
                normalized=flag,
                count_as_insertion=insertion,
                inplace=False,
                n_jobs=1,
            )
        )
        barcodes = [str(item) for item in adata.obs_names]
    fractions = [_scalars(item) for item in _as_list(values["peaks"])]
    cells = [
        {"barcode": barcode, "frip": fractions[index]} for index, barcode in enumerate(barcodes)
    ]
    return {
        "status": "completed",
        "operation": "fragment_in_peaks_fraction",
        "cell_count": len(cells),
        "normalized": flag,
        "count_as_insertion": insertion,
        "cells": cells,
        "provenance": _snap_provenance(
            frip="fragments or insertions in a peak, divided by the cell total when normalized",
        ),
    }


def fragment_size_distribution(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    max_recorded_size: int = 1000,
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dataset-wide fragment length histogram, with snapATAC2 ``metrics.frag_size_distr``.

    ``counts[i]`` is the number of fragments of length ``i`` base pairs, for lengths up to
    ``max_recorded_size``; longer fragments are not recorded.
    """

    limit = _integer(max_recorded_size, "max_recorded_size", 1, _MAX_SIZE)
    snap, context = _fragments(
        fragment_file, chrom_sizes, min_num_fragments, max_fragments, max_cells
    )
    with context as (adata, _work):
        values = _upstream(
            lambda: snap.metrics.frag_size_distr(
                adata,
                max_recorded_size=limit,
                inplace=False,
                n_jobs=1,
            )
        )
    counts = [int(item) for item in _as_list(values)]
    return {
        "status": "completed",
        "operation": "fragment_size_distribution",
        "max_recorded_size": limit,
        "fragment_count": sum(counts),
        "counts": counts,
        "provenance": _snap_provenance(
            size="fragment length in base pairs; counts[i] is the number of fragments of length i",
        ),
    }


def chromosome_fragment_counts(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    mode: str = "count",
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-cell fragment summary for each chromosome, with snapATAC2 ``metrics.summary_by_chrom``.

    ``mode`` is ``"count"`` (fragments), ``"sum"`` (sum of read counts) or ``"mean"``. Chromosomes
    are sorted by name, and each cell lists one value per chromosome in that order.
    """

    if mode not in {"count", "sum", "mean"}:
        raise ConfigurationError("mode must be one of count, sum, mean.")
    snap, context = _fragments(
        fragment_file, chrom_sizes, min_num_fragments, max_fragments, max_cells
    )
    with context as (adata, _work):
        values = _upstream(lambda: snap.metrics.summary_by_chrom(adata, mode=mode, n_jobs=1))
        barcodes = [str(item) for item in adata.obs_names]
    chromosomes = sorted(str(name) for name in values)
    columns = {name: [_scalars(item) for item in _as_list(values[name])] for name in chromosomes}
    cells = [
        {
            "barcode": barcode,
            "counts": [columns[name][index] for name in chromosomes],
        }
        for index, barcode in enumerate(barcodes)
    ]
    return {
        "status": "completed",
        "operation": "chromosome_fragment_counts",
        "mode": mode,
        "cell_count": len(cells),
        "chromosomes": chromosomes,
        "cells": cells,
        "provenance": _snap_provenance(mode=mode),
    }


def gene_activity_counts(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    gene_gtf: str | Path,
    *,
    upstream: int = 2000,
    downstream: int = 0,
    include_gene_body: bool = True,
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_entries: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cell-by-gene fragment counts from gene windows, with snapATAC2 ``pp.make_gene_matrix``.

    Each gene's window is ``upstream`` bases before its start, the gene body (when
    ``include_gene_body``) and ``downstream`` bases after its end, taken from ``gene_gtf``. Entries
    are ``[cell index, gene index, count]`` triplets sorted by row and column.
    """

    up = _integer(upstream, "upstream", 0, 10_000_000)
    down = _integer(downstream, "downstream", 0, 10_000_000)
    body = _flag(include_gene_body, "include_gene_body")
    gtf = _gtf(gene_gtf)
    entry_limit = _integer(max_entries, "max_entries", 1, _MAX_ENTRIES)
    snap, context = _fragments(
        fragment_file, chrom_sizes, min_num_fragments, max_fragments, max_cells
    )
    with context as (adata, work):
        gene_matrix = _upstream(
            lambda: snap.pp.make_gene_matrix(
                adata,
                gtf,
                inplace=False,
                file=work / "genes.h5ad",
                upstream=up,
                downstream=down,
                include_gene_body=body,
            )
        )
        try:
            barcodes = [str(item) for item in gene_matrix.obs_names]
            genes = [str(item) for item in gene_matrix.var_names]
            entries, total = _count_entries(gene_matrix.X[:], entry_limit)
        finally:
            gene_matrix.close()
    return {
        "status": "completed",
        "operation": "gene_activity_counts",
        "shape": [len(barcodes), len(genes)],
        "genes": genes,
        "nonzero_count": len(entries),
        "total_count": total,
        "entries": entries,
        "provenance": _snap_provenance(
            window=f"{up} bases upstream, {'gene body, ' if body else ''}{down} bases downstream",
        ),
    }


def regulatory_regions_to_genes(
    regions: Sequence[str],
    gene_gtf: str | Path,
    *,
    upstream: int = 250_000,
    downstream: int = 250_000,
    id_type: str = "gene_name",
    coding_gene_only: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Link regions to genes within a distance window (``tl.init_network_from_annotation``).

    A region is linked to a gene when the two are within ``upstream`` bases before the gene's start
    or ``downstream`` bases after its end. ``id_type`` names the gene identifier (``gene_name``,
    ``gene_id`` or ``transcript_id``). Edges are sorted by region, then gene.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    region_list = _regions(list(regions), "regions")
    gtf = _gtf(gene_gtf)
    up = _integer(upstream, "upstream", 0, 10_000_000)
    down = _integer(downstream, "downstream", 0, 10_000_000)
    kind = id_type
    if kind not in {"gene_name", "gene_id", "transcript_id"}:
        raise ConfigurationError("id_type must be gene_name, gene_id or transcript_id.")
    coding = _flag(coding_gene_only, "coding_gene_only")
    graph = _upstream(
        lambda: snap.tl.init_network_from_annotation(
            region_list,
            gtf,
            upstream=up,
            downstream=down,
            id_type=kind,
            coding_gene_only=coding,
        )
    )
    edges = []
    for source, target, data in graph.edge_index_map().values():
        left = graph[source]
        right = graph[target]
        if left.type == "region" and right.type == "gene":
            region_node, gene_node = left, right
        elif left.type == "gene" and right.type == "region":
            region_node, gene_node = right, left
        else:
            raise BackendExecutionError(
                "snapATAC2 returned an edge that does not join a region and a gene."
            )
        edges.append(
            {
                "region": str(region_node.id),
                "gene": str(gene_node.id),
                "distance": int(data.distance),
            }
        )
    edges.sort(key=lambda item: (item["region"], item["gene"]))
    return {
        "status": "completed",
        "operation": "regulatory_regions_to_genes",
        "region_count": len(region_list),
        "edge_count": len(edges),
        "edges": edges,
        "provenance": _snap_provenance(
            window=f"{up} bases upstream and {down} bases downstream of each gene",
            id_type=kind,
        ),
    }


def region_motif_enrichment(
    regions: Mapping[str, Sequence[str]],
    motifs: str | Path,
    genome: str | Path,
    *,
    background: Sequence[str] | None = None,
    method: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Motif enrichment in region groups against background (snapATAC2 ``tl.motif_enrichment``).

    ``regions`` maps group names to ``chr:start-end`` lists; ``motifs`` is a MEME file; ``genome``
    is a FASTA file supplying the region sequences. ``background`` (regions) defaults to the
    GC-matched set chosen by snapATAC2. ``method`` is ``"binomial"``, ``"hypergeometric"`` or None.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    if not isinstance(regions, Mapping) or not regions or len(regions) > _MAX_GROUPS:
        raise ConfigurationError(
            f"regions must map 1 to {_MAX_GROUPS} group names to region lists."
        )
    groups = {
        _label(name, "group name"): _regions(list(values), f"regions for {name!r}")
        for name, values in regions.items()
    }
    motif_path = input_path(motifs)
    fasta = input_path(genome)
    background_list = None if background is None else _regions(list(background), "background")
    if method not in {None, "binomial", "hypergeometric"}:
        raise ConfigurationError("method must be None, binomial or hypergeometric.")
    parsed = _upstream(lambda: snap.read_motifs(str(motif_path)))
    if not parsed:
        raise ConfigurationError(f"{motif_path.name} contains no motifs.")
    if len(parsed) > _MAX_MOTIFS:
        raise ConfigurationError(f"{motif_path.name} has more than {_MAX_MOTIFS} motifs.")
    tables = _upstream(
        lambda: snap.tl.motif_enrichment(
            parsed,
            groups,
            genome_fasta=fasta,
            background=background_list,
            method=method,
        )
    )
    results = []
    for group in sorted(tables):
        rows = []
        for row in tables[group].to_dicts():
            rows.append({str(key): _scalars(value) for key, value in row.items()})
        rows.sort(key=lambda item: str(item.get("id")))
        results.append({"group": str(group), "rows": rows})
    return {
        "status": "completed",
        "operation": "region_motif_enrichment",
        "group_count": len(results),
        "motif_count": len(parsed),
        "groups": results,
        "provenance": _snap_provenance(
            test=method or "snapATAC2 default",
            genome=fasta.name,
            motifs=motif_path.name,
        ),
    }


def merge_peak_sets(
    peaks_by_group: Mapping[str, Sequence[Mapping[str, object]]],
    chrom_sizes: Mapping[str, int],
    *,
    half_width: int = 250,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge group-specific peaks into one non-overlapping set, with snapATAC2 ``tl.merge_peaks``.

    Each group gives MACS2-style peak rows with ``chrom``, ``start``, ``end``, ``name``,
    ``score``, ``strand``, ``signal_value``, ``p_value``, ``q_value`` and ``summit``. ``summit`` is
    the offset of the summit from ``start``, as in narrowPeak column 10. Each merged peak spans
    ``half_width`` bases either side of the summit; ``groups`` lists the inputs that contain it.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    polars = require("polars", _SNAPATAC2_EXTRA)
    if not isinstance(peaks_by_group, Mapping) or not peaks_by_group:
        raise ConfigurationError("peaks_by_group must map group names to peak rows.")
    sizes = _chrom_sizes(chrom_sizes)
    width = _integer(half_width, "half_width", 0, 1_000_000)
    frames = {}
    for group, rows in peaks_by_group.items():
        name = _label(group, "group name")
        if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence) or not rows:
            raise ConfigurationError(f"peaks for group {name!r} must be a non-empty list of rows.")
        columns: dict[str, list[object]] = {column: [] for column in _PEAK_COLUMNS}
        for number, row in enumerate(rows, start=1):
            if not isinstance(row, Mapping):
                raise ConfigurationError(f"peak {number} of group {name!r} must be a mapping.")
            missing = [column for column in _PEAK_COLUMNS if column not in row]
            if missing:
                raise ConfigurationError(f"peak {number} of group {name!r} lacks {missing}.")
            chrom = _label(row["chrom"], "chrom")
            if chrom not in sizes:
                raise ConfigurationError(
                    f"peak {number} uses chromosome {chrom!r}, not in chrom_sizes."
                )
            start = _integer(row["start"], "start", 0, sizes[chrom])
            end = _integer(row["end"], "end", start + 1, sizes[chrom])
            summit = _integer(row["summit"], "summit", 0, end - start)
            for column, value in (
                ("chrom", chrom),
                ("start", start),
                ("end", end),
                ("name", str(row["name"])),
                ("score", _integer(row["score"], "score", 0, 65_535)),
                ("strand", str(row["strand"])),
                ("signal_value", float(row["signal_value"])),
                ("p_value", float(row["p_value"])),
                ("q_value", float(row["q_value"])),
                ("summit", summit),
            ):
                columns[column].append(value)
        frames[name] = polars.DataFrame(
            {
                "chrom": columns["chrom"],
                "start": polars.Series(columns["start"], dtype=polars.UInt64),
                "end": polars.Series(columns["end"], dtype=polars.UInt64),
                "name": columns["name"],
                "score": polars.Series(columns["score"], dtype=polars.UInt16),
                "strand": columns["strand"],
                "signal_value": polars.Series(columns["signal_value"], dtype=polars.Float64),
                "p_value": polars.Series(columns["p_value"], dtype=polars.Float64),
                "q_value": polars.Series(columns["q_value"], dtype=polars.Float64),
                "peak": polars.Series(columns["summit"], dtype=polars.UInt64),
            }
        )
    merged = _upstream(lambda: snap.tl.merge_peaks(frames, sizes, half_width=width))
    peaks = []
    for row in merged.to_dicts():
        text = str(row["Peaks"])
        chrom, span = text.rsplit(":", 1)
        start, end = (int(part) for part in span.split("-"))
        peaks.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "groups": sorted(name for name in frames if bool(row[name])),
            }
        )
    peaks.sort(key=lambda item: (item["chrom"], item["start"], item["end"]))
    return {
        "status": "completed",
        "operation": "merge_peak_sets",
        "group_count": len(frames),
        "peak_count": len(peaks),
        "half_width": width,
        "peaks": peaks,
        "provenance": _snap_provenance(
            merge="peaks widened by half_width around summits, with overlaps resolved by snapATAC2",
        ),
    }


__all__ = [
    "chromosome_fragment_counts",
    "fragment_in_peaks_fraction",
    "fragment_size_distribution",
    "gene_activity_counts",
    "merge_peak_sets",
    "region_motif_enrichment",
    "regulatory_regions_to_genes",
    "tss_enrichment",
]

_register_public_verbose(globals(), tuple(__all__))
