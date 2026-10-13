"""Gtars wrappers for genome sequences, region sets, gene partitions, LOLA, tokens, and VRS alleles.

The upstream package exposes its functions from submodules (``gtars.genomic_distributions``,
``gtars.lola``, ``gtars.refget``, ``gtars.tokenizers``, ``gtars.vrs``, and ``gtars.vrs.hgvs``)
rather than from the top-level namespace, which is why the package inventory was empty.

Each public function reads local files, runs one upstream operation, and returns a
JSON-serializable dictionary with ``status``, ``operation``, result fields, and ``provenance``.
Region coordinates are 0-based half-open (BED) unless a docstring says otherwise.

Functions are importable as ``dnakit.annotation.gtars_functions.<name>``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_ALGORITHMS = ("sha512t24u", "md5")
_COMMENT_PREFIXES = ("#", "track", "browser")


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int = 512) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _choice(value: object, name: str, allowed: Iterable[str]) -> str:
    options = tuple(allowed)
    if not isinstance(value, str) or value not in options:
        raise ConfigurationError(f"{name} must be one of {', '.join(options)}.")
    return value


def _lengths(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value or len(value) > 100_000:
        raise ConfigurationError(
            f"{name} must be a non-empty mapping from chromosome name to length."
        )
    return {
        _text(chrom, f"{name} chromosome"): _integer(length, f"{name} length", 1, 10**12)
        for chrom, length in value.items()
    }


def _plain(value: Any) -> Any:
    """Convert numpy, gtars, and missing values to JSON-compatible Python objects."""

    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        return _plain(value.tolist())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _record_order(record: Mapping[str, Any]) -> tuple[tuple[int, Any], ...]:
    """Type-safe sort key so rows come back in one order regardless of hash seed."""

    keys: list[tuple[int, Any]] = []
    for value in record.values():
        if value is None:
            keys.append((3, ""))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            keys.append((0, float(value)))
        else:
            keys.append((2, str(value)))
    return tuple(keys)


def _ordered(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(records, key=_record_order)


def _checked(rows: list[Any], limit: int, what: str) -> list[Any]:
    if len(rows) > limit:
        raise ConfigurationError(f"{what} has {len(rows)} rows; max_rows={limit} is the limit.")
    return rows


def _new_output(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output file already exists and will not be overwritten: {path}")
    return path


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("gtars"), "convention": convention}


_BED_CONVENTION = "0-based half-open coordinates (BED)"


def _bed_intervals(path: Path, *, columns: int | None = None) -> list[tuple[str, int, int]]:
    """Parse BED lines for validation; upstream readers do the computation."""

    intervals: list[tuple[str, int, int]] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            if not text or text.startswith(_COMMENT_PREFIXES):
                continue
            fields = text.split("\t")
            if len(fields) < 3 or (columns is not None and len(fields) != columns):
                expected = f"exactly {columns}" if columns is not None else "at least three"
                raise ConfigurationError(f"{path}:{number} needs {expected} BED columns.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(f"{path}:{number} has non-integer coordinates.") from exc
            if start < 0 or end <= start:
                raise ConfigurationError(f"{path}:{number} must satisfy 0 <= start < end.")
            intervals.append((fields[0], start, end))
    if not intervals:
        raise ConfigurationError(f"{path} contains no intervals.")
    return intervals


def _region_set(gtars_models: Any, path: Path) -> Any:
    try:
        return gtars_models.RegionSet(str(path))
    except Exception as exc:  # gtars reports malformed files with runtime errors
        raise ConfigurationError(f"gtars could not read the region set {path}: {exc}") from exc


def _upstream(call: Any, what: str) -> Any:
    """Run one gtars call; its runtime errors name bad inputs, so report them as configuration."""

    try:
        return call()
    except Exception as exc:  # gtars raises RuntimeError for missing sequence or bad coordinates
        raise ConfigurationError(f"gtars could not {what}: {exc}") from exc


# Genome sequence and region-level composition.


def region_gc_content(
    regions: str | Path,
    genome_fasta: str | Path,
    *,
    ignore_unknown_chromosomes: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the GC fraction of each BED region, using the sequence in a FASTA genome."""

    region_path = input_path(regions)
    fasta_path = input_path(genome_fasta)
    limit = _rows_limit(max_rows)
    ignore = _flag(ignore_unknown_chromosomes, "ignore_unknown_chromosomes")
    _checked(_bed_intervals(region_path), limit, "Region set")
    models = require("gtars.models", "bio-gtars")
    genomic = require("gtars.genomic_distributions", "bio-gtars")
    region_set = _region_set(models, region_path)
    genome = models.GenomeAssembly(str(fasta_path))
    values = _upstream(
        lambda: genomic.calc_gc_content(region_set, genome, ignore_unk_chroms=ignore),
        "compute GC content for the regions",
    )
    labels = [(str(region.chr), int(region.start), int(region.end)) for region in region_set]
    if len(values) != len(labels):
        raise ConfigurationError("gtars returned a GC value for a different number of regions.")
    rows = _ordered(
        [
            {"chrom": chrom, "start": start, "end": end, "gc_fraction": _plain(value)}
            for (chrom, start, end), value in zip(labels, values, strict=True)
        ]
    )
    finite = [row["gc_fraction"] for row in rows if row["gc_fraction"] is not None]
    return {
        "status": "completed",
        "operation": "region_gc_content",
        "region_count": len(rows),
        "mean_gc_fraction": round(sum(finite) / len(finite), 6) if finite else None,
        "regions": rows,
        "provenance": _provenance(_BED_CONVENTION + "; GC fraction in [0, 1]; N bases are not GC"),
    }


def region_dinucleotide_frequencies(
    regions: str | Path,
    genome_fasta: str | Path,
    *,
    raw_counts: bool = False,
    ignore_unknown_chromosomes: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count or percentage-normalize the 16 dinucleotides inside each BED region."""

    region_path = input_path(regions)
    fasta_path = input_path(genome_fasta)
    limit = _rows_limit(max_rows)
    raw = _flag(raw_counts, "raw_counts")
    ignore = _flag(ignore_unknown_chromosomes, "ignore_unknown_chromosomes")
    _checked(_bed_intervals(region_path), limit, "Region set")
    models = require("gtars.models", "bio-gtars")
    genomic = require("gtars.genomic_distributions", "bio-gtars")
    region_set = _region_set(models, region_path)
    genome = models.GenomeAssembly(str(fasta_path))
    result = _upstream(
        lambda: genomic.calc_dinucl_freq(
            region_set, genome, raw_counts=raw, ignore_unk_chroms=ignore
        ),
        "count dinucleotides in the regions",
    )
    names = [str(name).upper() for name in result["dinucleotides"]]
    rows = [
        {
            "label": str(label),
            "frequencies": {name: _plain(value) for name, value in zip(names, values, strict=True)},
        }
        for label, values in zip(result["region_labels"], result["frequencies"], strict=True)
    ]
    pooled = {name: _plain(sum(float(row["frequencies"][name]) for row in rows)) for name in names}
    return {
        "status": "completed",
        "operation": "region_dinucleotide_frequencies",
        "dinucleotides": names,
        "region_count": len(rows),
        "value_unit": "count" if raw else "percent",
        "regions": _ordered(rows),
        "pooled": pooled,
        "provenance": _provenance(_BED_CONVENTION + "; labels are chrom_start_end"),
    }


# Gene partitions and signal summaries.


def partition_region_counts(
    regions: str | Path,
    gene_annotation: str | Path,
    *,
    core_promoter_bp: int = 100,
    proximal_promoter_bp: int = 2000,
    protein_coding_only: bool = True,
    chrom_sizes: Mapping[str, int] | None = None,
    bp_proportion: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Classify BED regions into core promoter, proximal promoter, exon, intron, and intergenic.

    Without ``chrom_sizes`` the result holds observed counts. With ``chrom_sizes`` it also holds
    expected counts under a uniform-genome model and chi-square p-values.
    """

    region_path = input_path(regions)
    gtf_path = input_path(gene_annotation)
    core = _integer(core_promoter_bp, "core_promoter_bp", 0, 100_000)
    proximal = _integer(proximal_promoter_bp, "proximal_promoter_bp", core, 1_000_000)
    coding = _flag(protein_coding_only, "protein_coding_only")
    proportion = _flag(bp_proportion, "bp_proportion")
    sizes = _lengths(chrom_sizes, "chrom_sizes") if chrom_sizes is not None else None
    _bed_intervals(region_path)
    models = require("gtars.models", "bio-gtars")
    genomic = require("gtars.genomic_distributions", "bio-gtars")
    region_set = _region_set(models, region_path)
    partitions = models.PartitionList.from_gtf(
        str(gtf_path),
        core_prom=core,
        prox_prom=proximal,
        filter_protein_coding=coding,
        convert_ensembl_ucsc=True,
        chrom_sizes=sizes,
    )
    if sizes is None:
        observed = genomic.calc_partitions(region_set, partitions, bp_proportion=proportion)
        rows = [
            {"partition": str(name), "observed": _plain(count)}
            for name, count in zip(observed["partition"], observed["count"], strict=True)
        ]
        total = int(observed["total"])
        mode = "observed"
    else:
        expected = genomic.calc_expected_partitions(
            region_set, partitions, sizes, bp_proportion=proportion
        )
        rows = [
            {
                "partition": str(name),
                "observed": _plain(obs),
                "expected": _plain(exp),
                "log10_odds": _plain(lor),
                "p_value": _plain(pval),
            }
            for name, obs, exp, lor, pval in zip(
                expected["partition"],
                expected["observed"],
                expected["expected"],
                expected["log10OE"],
                expected["pvalue"],
                strict=True,
            )
        ]
        total = sum(int(row["observed"] or 0) for row in rows)
        mode = "expected"
    return {
        "status": "completed",
        "operation": "partition_region_counts",
        "mode": mode,
        "unit": "base pairs" if proportion else "regions",
        "total": total,
        "partitions": rows,
        "provenance": _provenance(
            f"core promoter {core} bp, proximal promoter {proximal} bp; "
            f"protein_coding_only={coding}"
        ),
    }


def signal_matrix_summary(
    regions: str | Path,
    signal_table: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average per-condition signal over the query regions and summarize each condition.

    ``signal_table`` is a tab-separated file. Its first column holds labels ``chrom_start_end``
    and each further column holds one condition's value for that labeled region.
    """

    region_path = input_path(regions)
    table_path = input_path(signal_table)
    _bed_intervals(region_path)
    models = require("gtars.models", "bio-gtars")
    genomic = require("gtars.genomic_distributions", "bio-gtars")
    region_set = _region_set(models, region_path)
    try:
        matrix = models.SignalMatrix.from_tsv(str(table_path))
        summary = genomic.calc_summary_signal(region_set, matrix)
    except Exception as exc:  # gtars reports malformed matrices with runtime errors
        raise ConfigurationError(f"gtars could not summarize the signal table: {exc}") from exc
    conditions = [str(name) for name in summary["condition_names"]]
    rows = [
        {
            "label": str(label),
            "values": {name: _plain(value) for name, value in zip(conditions, values, strict=True)},
        }
        for label, values in zip(summary["region_labels"], summary["signal_matrix"], strict=True)
    ]
    stats = [
        {key: _plain(value) for key, value in item.items()} for item in summary["matrix_stats"]
    ]
    return {
        "status": "completed",
        "operation": "signal_matrix_summary",
        "condition_names": conditions,
        "region_count": len(rows),
        "regions": _ordered(rows),
        "matrix_stats": stats,
        "provenance": _provenance(_BED_CONVENTION + "; signal from overlapping labeled regions"),
    }


def consensus_regions(
    region_files: Sequence[str | Path],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count, for each distinct region, how many of the given BED files contain it."""

    paths = [input_path(item) for item in region_files]
    if len(paths) < 2:
        raise ConfigurationError("At least two region files are required for a consensus.")
    limit = _rows_limit(max_rows)
    models = require("gtars.models", "bio-gtars")
    genomic = require("gtars.genomic_distributions", "bio-gtars")
    sets = [_region_set(models, path) for path in paths]
    rows = [
        {
            "chrom": str(item["chr"]),
            "start": int(item["start"]),
            "end": int(item["end"]),
            "count": int(item["count"]),
        }
        for item in genomic.consensus(sets)
    ]
    _checked(rows, limit, "Consensus")
    return {
        "status": "completed",
        "operation": "consensus_regions",
        "file_count": len(paths),
        "region_count": len(rows),
        "regions": _ordered(rows),
        "provenance": _provenance(_BED_CONVENTION + "; count is the number of input files"),
    }


# LOLA enrichment and fragment tokens.


def lola_enrichment(
    user_regions: str | Path,
    universe: str | Path,
    database_files: Sequence[str | Path],
    *,
    min_overlap: int = 1,
    direction: str = "enrichment",
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether the user regions overlap each database region set more than expected.

    The user set and every database set must lie inside the universe. Database sets are named
    by file name, so file names in ``database_files`` must be unique.
    """

    user_path = input_path(user_regions)
    universe_path = input_path(universe)
    files = [input_path(item) for item in database_files]
    if not files:
        raise ConfigurationError("At least one database region file is required.")
    names = [path.name for path in files]
    if len(set(names)) != len(names):
        raise ConfigurationError("Database region files must have unique file names.")
    overlap = _integer(min_overlap, "min_overlap", 1, 10**12)
    mode = _choice(direction, "direction", ("enrichment", "depletion"))
    limit = _rows_limit(max_rows)
    user = [(chrom, start, end) for chrom, start, end in _bed_intervals(user_path)]
    background = [(chrom, start, end) for chrom, start, end in _bed_intervals(universe_path)]
    lola = require("gtars.lola", "bio-gtars")
    database = lola.RegionDB.from_bed_files([str(path) for path in files], names)
    result = lola.run_lola([user], background, database, min_overlap=overlap, direction=mode)
    columns = {key: list(values) for key, values in result.items()}
    total = len(columns.get("dbSet", []))
    if total > limit:
        raise ConfigurationError(f"LOLA result has {total} rows; max_rows={limit} is the limit.")
    rows = []
    for index in range(total):
        db_index = int(columns["dbSet"][index])
        rows.append(
            {
                "db_set": db_index,
                "db_name": str(columns["filename"][index] or names[db_index]),
                "p_value_log10": _plain(columns["pValueLog"][index]),
                "odds_ratio": _plain(columns["oddsRatio"][index]),
                "support": _plain(columns["support"][index]),
                "rank_p_value": _plain(columns["rnkPV"][index]),
                "b": _plain(columns["b"][index]),
                "c": _plain(columns["c"][index]),
                "d": _plain(columns["d"][index]),
                "db_size": _plain(columns["size"][index]),
                "q_value": _plain(columns["qValue"][index]),
            }
        )
    return {
        "status": "completed",
        "operation": "lola_enrichment",
        "direction": mode,
        "database_count": len(files),
        "results": _ordered(rows),
        "provenance": _provenance(
            _BED_CONVENTION + f"; LOLA min_overlap={overlap}; universe-restricted contingency"
        ),
    }


def tokenize_fragments(
    fragments: str | Path,
    universe: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map fragment file records onto BED3 universe regions, one token list per barcode."""

    fragment_path = input_path(fragments)
    universe_path = input_path(universe)
    limit = _rows_limit(max_rows)
    universe_regions = _bed_intervals(universe_path, columns=3)
    _checked(universe_regions, limit, "Universe")
    tokenizers = require("gtars.tokenizers", "bio-gtars")
    try:
        tokenizer = tokenizers.Tokenizer.from_bed(str(universe_path))
        mapping = tokenizers.tokenize_fragment_file(str(fragment_path), tokenizer)
    except Exception as exc:  # gtars reports malformed fragment files with runtime errors
        raise ConfigurationError(f"gtars could not tokenize the fragments: {exc}") from exc
    barcodes = [
        {
            "barcode": str(barcode),
            "token_ids": sorted(int(token) for token in tokens),
            "token_count": len(tokens),
        }
        for barcode, tokens in mapping.items()
    ]
    _checked(barcodes, limit, "Barcodes")
    return {
        "status": "completed",
        "operation": "tokenize_fragments",
        "universe_size": len(universe_regions),
        "barcode_count": len(barcodes),
        "barcodes": _ordered(barcodes),
        "provenance": _provenance(
            "fragment records map to universe regions by overlap; token ids are universe indexes"
        ),
    }


# FASTA sequences, digests, and identifiers.


def digest_fasta_sequences(
    fasta: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute GA4GH sha512t24u and MD5 digests for every sequence in a FASTA file."""

    source = input_path(fasta)
    limit = _rows_limit(max_rows)
    refget = require("gtars.refget", "bio-gtars")
    try:
        collection = refget.digest_fasta(str(source))
    except Exception as exc:  # gtars reports unreadable FASTA with runtime errors
        raise ConfigurationError(f"gtars could not digest {source.name}: {exc}") from exc
    rows = []
    for record in collection.sequences:
        metadata = record.metadata
        rows.append(
            {
                "name": str(metadata.name),
                "length": int(metadata.length),
                "sha512t24u": str(metadata.sha512t24u),
                "md5": str(metadata.md5),
            }
        )
    _checked(rows, limit, "FASTA file")
    return {
        "status": "completed",
        "operation": "digest_fasta_sequences",
        "sequence_count": len(rows),
        "collection_digest": str(collection.digest),
        "sequences": _ordered(rows),
        "provenance": _provenance("digests follow the GA4GH refget sequence specification"),
    }


def fasta_index_records(
    fasta: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the samtools-style FASTA index (name and length) for each sequence in a FASTA file."""

    source = input_path(fasta)
    limit = _rows_limit(max_rows)
    refget = require("gtars.refget", "bio-gtars")
    try:
        records = refget.compute_fai(str(source))
    except Exception as exc:  # uncompressed FASTA only; gtars reports other failures
        raise ConfigurationError(f"gtars could not index {source.name}: {exc}") from exc
    rows = [{"name": str(item.name), "length": int(item.length)} for item in records]
    _checked(rows, limit, "FASTA index")
    return {
        "status": "completed",
        "operation": "fasta_index_records",
        "sequence_count": len(rows),
        "records": _ordered(rows),
        "provenance": _provenance("uncompressed FASTA; lengths are residue counts"),
    }


def sequence_digest(
    sequence: str,
    *,
    algorithm: str = "sha512t24u",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the sha512t24u (GA4GH refget) or MD5 digest of one DNA sequence."""

    text = sequence_text(sequence)
    name = _choice(algorithm, "algorithm", _ALGORITHMS)
    refget = require("gtars.refget", "bio-gtars")
    digest = refget.sha512t24u_digest(text) if name == "sha512t24u" else refget.md5_digest(text)
    return {
        "status": "completed",
        "operation": "sequence_digest",
        "algorithm": name,
        "digest": str(digest),
        "length": len(text),
        "provenance": _provenance("sequence is uppercased before digesting, as in digest_fasta"),
    }


# Variant normalization and identifiers.


def normalize_variant_allele(
    sequence: str,
    start: int,
    ref_allele: str,
    alt_allele: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim shared bases from an allele pair and expand it over repeats (VRS normalization)."""

    text = sequence_text(sequence)
    position = _integer(start, "start", 0, len(text))
    ref = sequence_text(ref_allele)
    alt = sequence_text(alt_allele)
    if text[position : position + len(ref)] != ref:
        raise ConfigurationError("ref_allele does not match the reference at start.")
    vrs = require("gtars.vrs", "bio-gtars")
    try:
        normalized = vrs.normalize_allele(text, position, ref, alt)
    except Exception as exc:  # gtars rejects invalid allele pairs with runtime errors
        raise ConfigurationError(f"gtars could not normalize the allele: {exc}") from exc
    return {
        "status": "completed",
        "operation": "normalize_variant_allele",
        "input_start": position,
        "start": int(normalized["start"]),
        "end": int(normalized["end"]),
        "allele": str(normalized["allele"]),
        "provenance": _provenance(_BED_CONVENTION + "; start is 0-based on the reference"),
    }


def vrs_allele_id(
    sequence_digest_value: str,
    start: int,
    end: int,
    alt: str,
    *,
    prefixed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the GA4GH VRS identifier (``ga4gh:VA.<digest>``) of one allele on a sequence.

    With ``prefixed=False`` the bare VRS digest without the ``ga4gh:VA.`` prefix is returned.
    """

    digest = _text(sequence_digest_value, "sequence_digest", maximum=256)
    first = _integer(start, "start", 0, 10**12)
    last = _integer(end, "end", first, 10**12)
    allele = sequence_text(alt)
    flag = _flag(prefixed, "prefixed")
    vrs = require("gtars.vrs", "bio-gtars")
    identifier = (
        vrs.vrs_id(digest, first, last, allele)
        if flag
        else vrs.vrs_digest(digest, first, last, allele)
    )
    return {
        "status": "completed",
        "operation": "vrs_allele_id",
        "identifier": str(identifier),
        "prefixed": flag,
        "sequence_digest": digest,
        "start": first,
        "end": last,
        "alt": allele,
        "provenance": _provenance("VRS SequenceLocation with 0-based start and exclusive end"),
    }


def parse_hgvs_variant(
    text: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse one HGVS expression such as ``NC_000001.11:g.12345A>G`` into its components."""

    expression = _text(text, "text", maximum=2048)
    hgvs = require("gtars.vrs.hgvs", "bio-gtars")
    try:
        parsed = hgvs.parse_hgvs(expression)
    except Exception as exc:  # HgvsError for unsupported or malformed expressions
        raise ConfigurationError(f"gtars could not parse the HGVS expression: {exc}") from exc
    variant = _plain(parsed.to_dict())
    return {
        "status": "completed",
        "operation": "parse_hgvs_variant",
        "accession": variant.get("accession"),
        "reference_type": variant.get("reference_type"),
        "variant": variant,
        "provenance": _provenance("HGVS parsed to its accession, reference type, and edits"),
    }


__all__ = [
    "consensus_regions",
    "digest_fasta_sequences",
    "fasta_index_records",
    "lola_enrichment",
    "normalize_variant_allele",
    "parse_hgvs_variant",
    "partition_region_counts",
    "region_dinucleotide_frequencies",
    "region_gc_content",
    "sequence_digest",
    "signal_matrix_summary",
    "tokenize_fragments",
    "vrs_allele_id",
]

_register_public_verbose(globals(), tuple(__all__))
