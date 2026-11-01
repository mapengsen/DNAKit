"""gtars region-set, gene-annotation, signal, tokenizer and region-database operations (integration
pack P24, part 2).

All region coordinates are 0-based half-open (BED). Region sets are read from BED files; results
that are region sets are returned inline (up to ``max_rows``) and, when ``output`` is given, written
to a new BED file. gtars is imported only when a function runs; a missing package raises
``BackendUnavailableError``."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.gtars_functions import (
    _BED_CONVENTION,
    _DEFAULT_ROWS,
    _bed_intervals,
    _choice,
    _flag,
    _integer,
    _lengths,
    _new_output,
    _plain,
    _provenance,
    _region_set,
    _rows_limit,
    _text,
    _upstream,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-gtars"
_PAIR_OPERATIONS = ("pintersect", "concat", "intersect_all")
_TOKENIZER_SOURCES = ("bed", "config", "pretrained")


def _models() -> Any:
    return require("gtars.models", _EXTRA)


def _load(path: str | Path, limit: int | None = None) -> tuple[Any, Path]:
    source = input_path(path)
    rows = _bed_intervals(source)
    if limit is not None and len(rows) > limit:
        raise ConfigurationError(f"Region set has {len(rows)} rows; max_rows={limit} is the limit.")
    return _region_set(_models(), source), source


def _region_row(region: Any) -> dict[str, Any]:
    row = {"chrom": str(region.chr), "start": int(region.start), "end": int(region.end)}
    rest = getattr(region, "rest", "")
    if rest:
        row["rest"] = str(rest)
    return row


def _emit(region_set: Any, output: str | Path | None, limit: int) -> dict[str, Any]:
    """Return the regions of a result set (capped by ``limit``) and write them when ``output`` is
    given."""

    written = None
    if output is not None:
        target = _new_output(output)
        _upstream(lambda: region_set.to_bed(str(target)), "write the region set")
        written = str(target)
    total = len(region_set)
    return {
        "region_count": total,
        "regions": [_region_row(region_set[i]) for i in range(min(total, limit))],
        "truncated": total > limit,
        "output": written,
    }


def _prov(note: str) -> dict[str, str]:
    return _provenance(_BED_CONVENTION + "; " + note)


def region_set_summary(
    regions: str | Path,
    *,
    include_widths: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a BED region set: counts, widths, total nucleotides, maximum end and per-chromosome
    statistics.

    Covers ``RegionSet.widths`` / ``region_widths`` (with ``include_widths=True``),
    ``mean_region_width``, ``get_nucleotide_length``, ``get_max_end_per_chr`` and
    ``chromosome_statistics``."""

    limit = _rows_limit(max_rows)
    widths_wanted = _flag(include_widths, "include_widths")
    region_set, _ = _load(regions, limit if widths_wanted else None)
    stats = _upstream(region_set.chromosome_statistics, "compute chromosome statistics")
    per_chrom = []
    for chrom in sorted(stats):
        item = stats[chrom]
        per_chrom.append(
            {
                "chromosome": str(chrom),
                "number_of_regions": int(item.number_of_regions),
                "minimum_region_length": int(item.minimum_region_length),
                "maximum_region_length": int(item.maximum_region_length),
                "mean_region_length": float(item.mean_region_length),
                "median_region_length": float(item.median_region_length),
                "start_nucleotide_position": int(item.start_nucleotide_position),
                "end_nucleotide_position": int(item.end_nucleotide_position),
            }
        )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "region_set_summary",
        "region_count": len(region_set),
        "mean_region_width": float(region_set.mean_region_width()),
        "total_nucleotides": int(region_set.get_nucleotide_length()),
        "max_end_per_chromosome": {
            k: int(v) for k, v in sorted(region_set.get_max_end_per_chr().items())
        },
        "chromosomes": per_chrom,
        "identifier": str(region_set.identifier),
    }
    if widths_wanted:
        result["widths"] = [int(item) for item in region_set.widths()]
    result["provenance"] = _prov("widths are end - start in file order")
    return result


def region_neighbor_distances(
    regions: str | Path,
    *,
    kind: str = "signed_gap",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distances between neighbouring regions on each chromosome.

    ``kind="signed_gap"`` returns the signed gaps between consecutive regions
    (``RegionSet.neighbor_distances``); ``kind="nearest"`` returns each region's distance to its
    nearest neighbour (``RegionSet.nearest_neighbors``). Chromosomes with too few regions are
    skipped, so the list is not aligned one to one with the input."""

    mode = _choice(kind, "kind", ("signed_gap", "nearest"))
    limit = _rows_limit(max_rows)
    region_set, _ = _load(regions, limit)
    call = region_set.neighbor_distances if mode == "signed_gap" else region_set.nearest_neighbors
    values = [int(item) for item in _upstream(call, "compute neighbour distances")]
    return {
        "status": "completed",
        "operation": "region_neighbor_distances",
        "kind": mode,
        "distance_count": len(values),
        "distances": values,
        "provenance": _prov("output is not aligned one to one with the input regions"),
    }


def region_distribution_bins(
    regions: str | Path,
    *,
    n_bins: int = 250,
    chrom_sizes: Mapping[str, int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count region midpoints in equal-width genomic bins (``RegionSet.distribution``).

    Without ``chrom_sizes`` the bin width comes from the largest end coordinate in the file;
    with it, each chromosome is cut into ``n_bins`` bins and regions beyond the stated size are
    skipped."""

    bins = _integer(n_bins, "n_bins", 1, 1_000_000)
    sizes = None if chrom_sizes is None else _lengths(chrom_sizes, "chrom_sizes")
    region_set, _ = _load(regions)
    rows = _upstream(
        lambda: region_set.distribution(n_bins=bins, chrom_sizes=sizes), "bin the regions"
    )
    table = [
        {
            "chrom": str(row["chr"]),
            "start": int(row["start"]),
            "end": int(row["end"]),
            "count": int(row["n"]),
            "bin": int(row["rid"]),
        }
        for row in rows
    ]
    return {
        "status": "completed",
        "operation": "region_distribution_bins",
        "n_bins": bins,
        "bin_count": len(table),
        "bins": table,
        "provenance": _prov("regions are assigned to bins by midpoint"),
    }


def region_promoters(
    regions: str | Path,
    upstream: int,
    downstream: int,
    *,
    output: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build promoter windows ``[start - upstream, start + downstream)`` around every region start.

    gtars measures from the region start and ignores strand; the windows are clamped at 0.
    """

    up = _integer(upstream, "upstream", 0, 10**9)
    down = _integer(downstream, "downstream", 0, 10**9)
    limit = _rows_limit(max_rows)
    region_set, _ = _load(regions, limit)
    promoters = _upstream(lambda: region_set.promoters(up, down), "build promoters")
    return {
        "status": "completed",
        "operation": "region_promoters",
        "upstream": up,
        "downstream": down,
        **_emit(promoters, output, limit),
        "provenance": _prov("window is relative to the region start, strand is not used"),
    }


def region_pairwise_operation(
    regions_a: str | Path,
    regions_b: str | Path,
    *,
    operation: str = "pintersect",
    output: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Combine two BED region sets.

    ``operation`` is ``pintersect`` (intersection of the i-th region of each set, by index),
    ``concat`` (all regions of both sets, nothing merged) or ``intersect_all`` (every overlapping
    pair, all-vs-all)."""

    mode = _choice(operation, "operation", _PAIR_OPERATIONS)
    limit = _rows_limit(max_rows)
    first, _ = _load(regions_a, limit)
    second, _ = _load(regions_b, limit)
    result = _upstream(lambda: getattr(first, mode)(second), f"run {mode}")
    return {
        "status": "completed",
        "operation": "region_pairwise_operation",
        "method": mode,
        **_emit(result, output, limit),
        "provenance": _prov("pintersect pairs regions by index position"),
    }


def region_overlap_coefficient(
    regions_a: str | Path,
    regions_b: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide overlap coefficient of two region sets: shared bases divided by the bases of the
    smaller set."""

    first, _ = _load(regions_a)
    second, _ = _load(regions_b)
    value = _upstream(lambda: first.overlap_coefficient(second), "compute the overlap coefficient")
    return {
        "status": "completed",
        "operation": "region_overlap_coefficient",
        "overlap_coefficient": float(value),
        "provenance": _prov("value in [0, 1]"),
    }


def region_subset_by_overlaps(
    regions_a: str | Path,
    regions_b: str | Path,
    *,
    output: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep the regions of the first set that overlap at least one region of the second set."""

    limit = _rows_limit(max_rows)
    first, _ = _load(regions_a, limit)
    second, _ = _load(regions_b, limit)
    kept = _upstream(lambda: first.subset_by_overlaps(second), "subset by overlaps")
    return {
        "status": "completed",
        "operation": "region_subset_by_overlaps",
        "input_region_count": len(first),
        **_emit(kept, output, limit),
        "provenance": _prov("regions of the first set are kept unchanged"),
    }


def region_overlap_indices(
    regions_a: str | Path,
    regions_b: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """For each region of the first set, the 0-based indices of the overlapping regions of the
    second set."""

    limit = _rows_limit(max_rows)
    first, _ = _load(regions_a, limit)
    second, _ = _load(regions_b, limit)
    found = _upstream(lambda: first.find_overlaps(second), "find overlaps")
    rows = [[int(i) for i in item] for item in found]
    return {
        "status": "completed",
        "operation": "region_overlap_indices",
        "region_count": len(rows),
        "overlapping_region_count": sum(1 for item in rows if item),
        "overlaps": rows,
        "provenance": _prov("indices refer to the second set in file order"),
    }


def region_gaps(
    regions: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    output: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the gaps between regions on each chromosome, from position 0 up to the chromosome
    size.

    Chromosomes of ``chrom_sizes`` without regions give one full-length gap; regions on other
    chromosomes are skipped."""

    sizes = _lengths(chrom_sizes, "chrom_sizes")
    limit = _rows_limit(max_rows)
    region_set, _ = _load(regions, limit)
    gaps = _upstream(lambda: region_set.gaps(sizes), "compute gaps")
    return {
        "status": "completed",
        "operation": "region_gaps",
        **_emit(gaps, output, limit),
        "provenance": _prov("leading and trailing gaps are included"),
    }


def region_disjoin(
    regions: str | Path,
    *,
    output: str | Path | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Break all regions into non-overlapping pieces: every region boundary becomes a piece
    boundary."""

    limit = _rows_limit(max_rows)
    region_set, _ = _load(regions, limit)
    pieces = _upstream(region_set.disjoin, "disjoin the regions")
    return {
        "status": "completed",
        "operation": "region_disjoin",
        **_emit(pieces, output, limit),
        "provenance": _prov("pieces tile the covered positions exactly"),
    }


# Gene annotation, TSS distances and signal matrices.


def _gene_sources(gtf: str | Path | None, annotation_binary: str | Path | None) -> tuple[str, Path]:
    if (gtf is None) == (annotation_binary is None):
        raise ConfigurationError("give exactly one of gtf and annotation_binary.")
    if gtf is not None:
        return "gtf", input_path(gtf)
    return "binary", input_path(annotation_binary)  # type: ignore[arg-type]


def _annotation(models: Any, kind: str, path: Path, coding: bool, convert: bool) -> Any:
    if kind == "gtf":
        return _upstream(
            lambda: models.GenomicDistAnnotation.from_gtf(
                str(path), filter_protein_coding=coding, convert_ensembl_ucsc=convert
            ),
            "read the GTF file",
        )
    return _upstream(
        lambda: models.GenomicDistAnnotation.load_bin(str(path)), "load the annotation binary"
    )


def gene_annotation_summary(
    *,
    gtf: str | Path | None = None,
    annotation_binary: str | Path | None = None,
    core_promoter_bp: int = 100,
    proximal_promoter_bp: int = 2000,
    protein_coding_only: bool = True,
    convert_ensembl_ucsc: bool = True,
    chrom_sizes: Mapping[str, int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the gtars gene model, partition list and TSS index from a GTF file (or a saved
    annotation binary).

    Returns the gene and exon counts, the partition names (core promoter, proximal promoter,
    exon, intron) and the TSS count. Covers ``GeneModel.from_gtf``,
    ``GenomicDistAnnotation.from_gtf`` / ``load_bin`` / ``gene_model`` / ``partition_list`` /
    ``tss_index`` and ``PartitionList.from_gene_model`` / ``from_gtf``."""

    kind, path = _gene_sources(gtf, annotation_binary)
    core = _integer(core_promoter_bp, "core_promoter_bp", 0, 100_000)
    proximal = _integer(proximal_promoter_bp, "proximal_promoter_bp", core, 1_000_000)
    coding = _flag(protein_coding_only, "protein_coding_only")
    convert = _flag(convert_ensembl_ucsc, "convert_ensembl_ucsc")
    sizes = None if chrom_sizes is None else _lengths(chrom_sizes, "chrom_sizes")
    models = _models()
    annotation = _annotation(models, kind, path, coding, convert)
    model = annotation.gene_model()
    partitions = _upstream(
        lambda: annotation.partition_list(core, proximal, sizes), "build the partitions"
    )
    index = _upstream(annotation.tss_index, "build the TSS index")
    return {
        "status": "completed",
        "operation": "gene_annotation_summary",
        "source": kind,
        "n_genes": int(model.n_genes),
        "n_exons": int(model.n_exons),
        "partition_names": [str(name) for name in partitions.partition_names()],
        "tss_count": len(index),
        "provenance": _provenance(
            f"core promoter {core} bp, proximal promoter {proximal} bp; "
            f"protein_coding_only={coding}"
        ),
    }


def tss_distances(
    regions: str | Path,
    *,
    tss_regions: str | Path | None = None,
    gtf: str | Path | None = None,
    annotation_binary: str | Path | None = None,
    signed: bool = False,
    protein_coding_only: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distance from each query region to the nearest TSS (or feature) with a gtars ``TssIndex``.

    The index comes from a BED file of TSS positions (``tss_regions``), a GTF file (strand-aware
    TSS) or a saved annotation binary. ``signed=False`` gives unsigned distances
    (``calc_tss_distances``); ``signed=True`` gives signed distances, null for chromosomes without
    features (``feature_distances``)."""

    flag = _flag(signed, "signed")
    coding = _flag(protein_coding_only, "protein_coding_only")
    limit = _rows_limit(max_rows)
    given = [item is not None for item in (tss_regions, gtf, annotation_binary)]
    if sum(given) != 1:
        raise ConfigurationError("give exactly one of tss_regions, gtf and annotation_binary.")
    models = _models()
    query, _ = _load(regions, limit)
    if tss_regions is not None:
        features, _ = _load(tss_regions)
        index = _upstream(lambda: models.TssIndex.from_regionset(features), "build the TSS index")
        origin = "regions"
    else:
        kind, path = _gene_sources(gtf, annotation_binary)
        index = _upstream(
            _annotation(models, kind, path, coding, True).tss_index, "build the TSS index"
        )
        origin = kind
    call = index.feature_distances if flag else index.calc_tss_distances
    values = _upstream(lambda: call(query), "compute TSS distances")
    rows = [{**_region_row(query[i]), "distance": _plain(values[i])} for i in range(len(values))]
    return {
        "status": "completed",
        "operation": "tss_distances",
        "index_source": origin,
        "signed": flag,
        "tss_count": len(index),
        "region_count": len(rows),
        "distances": rows,
        "provenance": _prov("distance is measured from the region midpoint to the TSS"),
    }


def signal_matrix_info(
    path: str | Path,
    *,
    binary: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a gtars signal matrix (TSV, or the packed binary with ``binary=True``) and report its
    shape and conditions."""

    source = input_path(path)
    packed = _flag(binary, "binary")
    models = _models()
    loader = models.SignalMatrix.load_bin if packed else models.SignalMatrix.from_tsv
    matrix = _upstream(lambda: loader(str(source)), "load the signal matrix")
    return {
        "status": "completed",
        "operation": "signal_matrix_info",
        "format": "binary" if packed else "tsv",
        "n_regions": int(matrix.n_regions),
        "n_conditions": int(matrix.n_conditions),
        "condition_names": [str(name) for name in matrix.condition_names],
        "provenance": _provenance("TSV first column holds labels chrom_start_end"),
    }


# Region tokenizers.


def _tokenizer(source: str, source_type: str) -> tuple[Any, str]:
    kind = _choice(source_type, "source_type", _TOKENIZER_SOURCES)
    tokenizers = require("gtars.tokenizers", _EXTRA)
    if kind == "pretrained":
        name = _text(source, "source", maximum=512)
        return _upstream(
            lambda: tokenizers.Tokenizer.from_pretrained(name), "load the pretrained tokenizer"
        ), name
    path = str(input_path(source))
    loader = tokenizers.Tokenizer.from_bed if kind == "bed" else tokenizers.Tokenizer.from_config
    return _upstream(lambda: loader(path), f"build the tokenizer from the {kind} file"), path


def region_tokenizer_vocabulary(
    source: str | Path,
    *,
    source_type: str = "bed",
    include_vocabulary: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a gtars region tokenizer and return its vocabulary and special tokens.

    ``source_type`` is ``bed`` (a BED universe, ``Tokenizer.from_bed``), ``config`` (a TOML
    tokenizer configuration, ``from_config``) or ``pretrained`` (a Hugging Face model name,
    ``from_pretrained``; needs network access). The vocabulary maps region tokens such as
    ``chr1:0-1000`` to integer ids."""

    keep = _flag(include_vocabulary, "include_vocabulary")
    tokenizer, label = _tokenizer(str(source), source_type)
    vocab = {str(k): int(v) for k, v in tokenizer.get_vocab().items()}
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "region_tokenizer_vocabulary",
        "source": label,
        "source_type": source_type,
        "vocab_size": int(tokenizer.vocab_size),
        "special_tokens": {k: v for k, v in sorted(tokenizer.special_tokens_map.items())},
    }
    if keep:
        result["vocabulary"] = dict(sorted(vocab.items(), key=lambda item: item[1]))
    result["provenance"] = _provenance(
        "token ids follow the universe order; special tokens come last"
    )
    return result


def region_token_ids(
    source: str | Path,
    values: Sequence[str | int],
    *,
    source_type: str = "bed",
    direction: str = "to_ids",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert between region tokens and integer ids with a gtars tokenizer.

    ``direction="to_ids"`` takes tokens such as ``chr1:0-1000`` and returns ids (unknown tokens map
    to the unknown id); ``direction="to_tokens"`` takes ids and returns tokens. Covers
    ``Tokenizer.encode`` / ``convert_tokens_to_ids`` and ``decode`` / ``convert_ids_to_tokens``."""

    way = _choice(direction, "direction", ("to_ids", "to_tokens"))
    if isinstance(values, (str, bytes)) or not values or len(values) > 1_000_000:
        raise ConfigurationError("values must be a non-empty list of at most 1000000 items.")
    tokenizer, label = _tokenizer(str(source), source_type)
    if way == "to_ids":
        tokens = [_text(item, "token", maximum=512) for item in values]
        converted = [
            int(item)
            for item in _upstream(lambda: tokenizer.convert_tokens_to_ids(tokens), "encode")
        ]
    else:
        ids = [_integer(item, "id", 0, 10**9) for item in values]
        converted = [
            str(item) for item in _upstream(lambda: tokenizer.convert_ids_to_tokens(ids), "decode")
        ]
    return {
        "status": "completed",
        "operation": "region_token_ids",
        "source": label,
        "direction": way,
        "inputs": list(values),
        "outputs": converted,
        "provenance": _provenance("one output per input, same order"),
    }


def region_special_tokens_mask(
    source: str | Path,
    tokens: Sequence[str],
    *,
    source_type: str = "bed",
    verbose: bool = True,
) -> dict[str, Any]:
    """Flag which of the given tokens are special tokens (``<cls>``, ``<pad>`` and so on) for a
    gtars tokenizer."""

    if isinstance(tokens, (str, bytes)) or not tokens or len(tokens) > 1_000_000:
        raise ConfigurationError("tokens must be a non-empty list of at most 1000000 strings.")
    items = [_text(item, "token", maximum=512) for item in tokens]
    tokenizer, label = _tokenizer(str(source), source_type)
    mask = _upstream(
        lambda: tokenizer.get_special_tokens_mask(items), "compute the special token mask"
    )
    return {
        "status": "completed",
        "operation": "region_special_tokens_mask",
        "source": label,
        "tokens": items,
        "mask": [bool(item) for item in mask],
        "provenance": _provenance("true marks a special token"),
    }


# LOLA region databases.


def _database(
    files: Sequence[str | Path] | None, folder: str | Path | None, names: Sequence[str] | None
) -> tuple[Any, str]:
    lola = require("gtars.lola", _EXTRA)
    if (files is None) == (folder is None):
        raise ConfigurationError("give exactly one of bed_files and folder.")
    if files is not None:
        if isinstance(files, (str, Path)) or not files:
            raise ConfigurationError("bed_files must be a non-empty list of BED paths.")
        paths = [str(input_path(item)) for item in files]
        labels = (
            [Path(item).name for item in paths]
            if names is None
            else [_text(n, "name", maximum=256) for n in names]
        )
        if len(labels) != len(paths) or len(set(labels)) != len(labels):
            raise ConfigurationError("names must be unique and match bed_files one to one.")
        return _upstream(
            lambda: lola.RegionDB.from_bed_files(paths, labels), "build the region database"
        ), "bed_files"
    directory = Path(folder).expanduser().resolve()  # type: ignore[arg-type]
    if not directory.is_dir():
        raise ConfigurationError(f"folder must be an existing LOLA database directory: {directory}")
    return _upstream(
        lambda: lola.RegionDB.from_folder(str(directory)), "read the LOLA folder"
    ), "folder"


def region_database_list_sets(
    *,
    bed_files: Sequence[str | Path] | None = None,
    folder: str | Path | None = None,
    names: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a LOLA region database from BED files (or a LOLA folder) and list its region sets.

    Covers ``RegionDB.from_bed_files``, ``RegionDB.from_folder`` and ``RegionDB.list_region_sets``.
    Set names default to the file names."""

    database, origin = _database(bed_files, folder, names)
    listing = [str(item) for item in database.list_region_sets()]
    return {
        "status": "completed",
        "operation": "region_database_list_sets",
        "source": origin,
        "set_count": int(database.num_region_sets),
        "sets": [{"index": i, "name": name} for i, name in enumerate(listing)],
        "provenance": _provenance("set index is the position in the database"),
    }


def region_database_get_sets(
    *,
    bed_files: Sequence[str | Path] | None = None,
    folder: str | Path | None = None,
    names: Sequence[str] | None = None,
    indices: Sequence[int] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the regions of selected region sets (by 0-based index, default all) of a LOLA region
    database."""

    limit = _rows_limit(max_rows)
    picked = None
    if indices is not None:
        if isinstance(indices, (str, bytes)) or not indices:
            raise ConfigurationError("indices must be a non-empty list of integers.")
        picked = [_integer(item, "index", 0, 10**9) for item in indices]
    database, origin = _database(bed_files, folder, names)
    total = int(database.num_region_sets)
    if picked is not None and any(item >= total for item in picked):
        raise ConfigurationError(f"indices must be below the number of region sets ({total}).")
    chosen = picked if picked is not None else list(range(total))
    labels = [str(item) for item in database.list_region_sets()]
    sets = []
    for index in chosen:
        group = _upstream(lambda i=index: database.get_region_sets([i]), "read the region set")
        region_set = group[0] if hasattr(group, "__getitem__") else next(iter(group))
        count = len(region_set)
        sets.append(
            {
                "index": index,
                "name": labels[index],
                "region_count": count,
                "regions": [_region_row(region_set[i]) for i in range(min(count, limit))],
                "truncated": count > limit,
            }
        )
    return {
        "status": "completed",
        "operation": "region_database_get_sets",
        "source": origin,
        "sets": sets,
        "provenance": _prov("regions of each set in file order"),
    }


__all__ = [
    "gene_annotation_summary",
    "region_database_get_sets",
    "region_database_list_sets",
    "region_disjoin",
    "region_distribution_bins",
    "region_gaps",
    "region_neighbor_distances",
    "region_overlap_coefficient",
    "region_overlap_indices",
    "region_pairwise_operation",
    "region_promoters",
    "region_set_summary",
    "region_special_tokens_mask",
    "region_subset_by_overlaps",
    "region_token_ids",
    "region_tokenizer_vocabulary",
    "signal_matrix_info",
    "tss_distances",
]

_register_public_verbose(globals(), tuple(__all__))
