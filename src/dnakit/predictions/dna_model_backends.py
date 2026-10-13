"""Optional adapters for DNA-level data utilities of DNA AI model packages.

Each public function calls one upstream package through ``require()``:

- ``alphagenome``: interval and variant parsing, interval merging, and GENCODE
  gene and transcript extraction. These run offline.
- ``grelu``: reverse complement, dinucleotide shuffling, resizing and single-base
  mutation of DNA strings.
- ``selene-sdk``: FASTA-backed sequence retrieval, contig lengths and one-hot encoding.
- ``enformer-pytorch``: one-hot encoding and FASTA windows built from BED intervals.
- ``genomic-benchmarks``: local dataset summaries, dataset names and cache checks.
- ``evo2`` and ``kipoi``: local Evo 2 generation and kipoi model descriptions.

Model inference and hosted prediction APIs are already in :mod:`dnakit.predictions`
and :mod:`dnakit.databases`, so they are not repeated here. Packages are imported
only when a function runs. A missing package raises ``BackendUnavailableError``
with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.metadata as metadata
import io
import math
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_ALPHAGENOME = "alphagenome"
_GRELU = "grelu"
_SELENE = "selene-sdk"
_ENFORMER = "enformer-pytorch"
_GENOMIC_BENCHMARKS = "genomic-benchmarks"
_EVO2 = "evo2"
_KIPOI = "kipoi"

_EXTRAS = {
    _ALPHAGENOME: "bio-alphagenome",
    _GRELU: "bio-grelu",
    _SELENE: "bio-selene",
    _ENFORMER: "neural-enformer",
    _GENOMIC_BENCHMARKS: "bio-genomic-benchmarks",
    _EVO2: "neural-evo2",
    _KIPOI: "bio-kipoi",
}
_DNA = frozenset("ACGTN")
_STRANDS = ("+", "-", ".")
_CHANNELS = ["A", "C", "G", "T"]
_VARIANT_FORMATS = ("default", "gtex", "open_targets", "open_targets_bigquery", "gnomad")
_GENCODE_FEATURES = frozenset({"gene", "transcript", "exon", "CDS", "start_codon", "stop_codon"})
_TRANSCRIPT_FEATURES = frozenset({"transcript", "exon", "CDS", "start_codon", "stop_codon"})
_GTF_COLUMNS = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Feature",
    "Frame",
    "gene_id",
    "gene_name",
    "gene_type",
    "transcript_id",
    "transcript_type",
    "protein_id",
]
_GTF_ATTRIBUTE = re.compile(r'([A-Za-z_][\w.-]*) "([^"]*)"')
_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,127}$")
# Copied from evo2/utils.py MODEL_NAMES (evo2 0.6.0); evo2 cannot be imported here.
_EVO2_MODELS = (
    "evo2_40b",
    "evo2_7b",
    "evo2_20b",
    "evo2_40b_base",
    "evo2_7b_base",
    "evo2_1b_base",
    "evo2_7b_262k",
    "evo2_7b_microviridae",
)
_MAX_TEXT = 256
_MAX_SEQUENCE = 1_000_000
_MAX_SEQUENCES = 10_000
_MAX_ONE_HOT = 100_000
_MAX_COORDINATE = 10**12
_MAX_RECORDS = 5_000_000
_DEFAULT_RECORDS = 200_000
_ENFORMER_CONTEXT = 196_608
_MAX_CONTEXT = 1_000_000
_MAX_WINDOWS = 10_000
_MAX_REGIONS = 100_000


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _real(value: object, name: str, *, low: float, high: float, low_open: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    below = value <= low if low_open else value < low
    if below or value > high:
        raise ConfigurationError(f"{name} must be in {'(' if low_open else '['}{low}, {high}].")
    return float(value)


def _text(value: object, name: str, *, maximum: int = _MAX_TEXT) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if "\x00" in value:
        raise ConfigurationError(f"{name} must not contain NUL characters.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum} characters.")
    return value.strip()


def _dna(value: object, name: str, *, maximum: int = _MAX_SEQUENCE) -> str:
    text = _text(value, name, maximum=maximum).upper()
    if set(text) - _DNA:
        raise ConfigurationError(f"{name} may contain only A, C, G, T and N.")
    return text


def _dna_batch(value: object, name: str) -> list[str]:
    if isinstance(value, str):
        return [_dna(value, name)]
    if isinstance(value, (bytes, Mapping)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a DNA string or a non-empty list of strings.")
    if len(value) > _MAX_SEQUENCES:
        raise ConfigurationError(f"{name} may contain at most {_MAX_SEQUENCES} sequences.")
    return [_dna(item, f"{name}[{index}]") for index, item in enumerate(value)]


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _clean(value: object) -> Any:
    """Convert numpy, torch and pandas values to JSON-safe Python values."""

    if hasattr(value, "tolist") and not isinstance(value, (str, bytes)):
        value = value.tolist()
    if isinstance(value, Mapping):
        return {str(key): _clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _provenance(package: str, **notes: str) -> dict[str, Any]:
    try:
        version = metadata.version(package)
    except metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": package, "version": version, **notes}


def _rejected(package: str, exc: Exception) -> ConfigurationError:
    return ConfigurationError(f"{package} rejected the input: {exc}")


def _names(value: object) -> list[str]:
    """Input or target names of a kipoi schema: a dict, a list, or one array object."""

    if value is None:
        return []
    if isinstance(value, Mapping):
        return [str(key) for key in value]
    if isinstance(value, (list, tuple)):
        return [str(getattr(item, "name", item)) for item in value]
    name = getattr(value, "name", None)
    return [] if name is None else [str(name)]


def _interval_record(interval: Any) -> dict[str, Any]:
    return {
        "chromosome": str(interval.chromosome),
        "start": int(interval.start),
        "end": int(interval.end),
        "strand": str(interval.strand),
        "width": int(interval.width),
    }


def _gtf_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in _GTF_ATTRIBUTE.findall(text)}


def _gtf_frame(path: Path, features: frozenset[str], max_records: int) -> Any:
    """Read the selected GENCODE rows into a DataFrame with 0-based half-open coordinates."""

    pandas = require("pandas", _EXTRAS[_ALPHAGENOME])
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[dict[str, Any]] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ConfigurationError(f"GTF line {line_number} has fewer than 9 columns.")
            feature = fields[2]
            if feature not in features:
                continue
            if len(rows) >= max_records:
                raise ConfigurationError(
                    f"GTF has more than max_records={max_records} selected rows; "
                    "increase the explicit limit."
                )
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError as exc:
                raise ConfigurationError(
                    f"GTF line {line_number} has invalid coordinates."
                ) from exc
            if fields[6] not in _STRANDS or start < 0 or end <= start:
                raise ConfigurationError(f"GTF line {line_number} has an invalid strand or span.")
            frame: int | None = None
            if feature in ("CDS", "stop_codon"):
                try:
                    frame = int(fields[7])
                except ValueError as exc:
                    raise ConfigurationError(
                        f"GTF line {line_number}: {feature} rows need a numeric frame."
                    ) from exc
            attributes = _gtf_attributes(fields[8])
            rows.append(
                {
                    "Chromosome": fields[0],
                    "Start": start,
                    "End": end,
                    "Strand": fields[6],
                    "Feature": feature,
                    "Frame": frame,
                    "gene_id": attributes.get("gene_id"),
                    "gene_name": attributes.get("gene_name"),
                    "gene_type": attributes.get("gene_type"),
                    "transcript_id": attributes.get("transcript_id"),
                    "transcript_type": attributes.get("transcript_type"),
                    "protein_id": attributes.get("protein_id"),
                }
            )
    return pandas.DataFrame(rows, columns=_GTF_COLUMNS)


def _region(item: object, index: int) -> tuple[str, int, int, str]:
    if not isinstance(item, Mapping):
        raise ConfigurationError(f"regions[{index}] must be a mapping with chrom, start and end.")
    name = f"regions[{index}]"
    chrom = _text(item.get("chrom"), f"{name}.chrom", maximum=128)
    start = _integer(item.get("start"), f"{name}.start", minimum=0, maximum=_MAX_COORDINATE)
    end = _integer(item.get("end"), f"{name}.end", minimum=1, maximum=_MAX_COORDINATE)
    if end <= start:
        raise ConfigurationError(f"{name}.end must be greater than start.")
    strand = _choice(item.get("strand", "+"), f"{name}.strand", _STRANDS)
    return chrom, start, end, strand


# alphagenome: offline interval, variant and GENCODE helpers


def parse_alphagenome_interval(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse ``chr:start-end[:strand]`` with alphagenome's Interval parser."""

    value = _text(text, "interval")
    genome = require("alphagenome.data.genome", _EXTRAS[_ALPHAGENOME])
    try:
        interval = genome.Interval.from_str(value)
    except (ValueError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "parse_alphagenome_interval",
        "interval": _interval_record(interval),
        "provenance": _provenance(
            _ALPHAGENOME,
            coordinate_convention="start is used as given (0-based); end is exclusive",
        ),
    }


def parse_alphagenome_variant(
    text: str, *, variant_format: str = "default", verbose: bool = True
) -> dict[str, Any]:
    """Parse a variant string such as ``chr1:10:A>G`` with alphagenome."""

    value = _text(text, "variant")
    fmt = _choice(variant_format, "variant_format", _VARIANT_FORMATS)
    genome = require("alphagenome.data.genome", _EXTRAS[_ALPHAGENOME])
    try:
        variant = genome.Variant.from_str(value, variant_format=genome.VariantFormat(fmt))
    except (ValueError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "parse_alphagenome_variant",
        "variant": {
            "chromosome": str(variant.chromosome),
            "position": int(variant.position),
            "reference": str(variant.reference_bases),
            "alternate": str(variant.alternate_bases),
            "start": int(variant.start),
            "end": int(variant.end),
            "is_snv": bool(variant.is_snv),
            "is_indel": bool(variant.is_indel),
            "is_insertion": bool(variant.is_insertion),
            "is_deletion": bool(variant.is_deletion),
            "is_structural": bool(variant.is_structural),
        },
        "provenance": _provenance(
            _ALPHAGENOME,
            coordinate_convention="position is 1-based as in VCF; start is position - 1",
            variant_format=fmt,
        ),
    }


def merge_alphagenome_intervals(
    intervals: Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Merge overlapping interval strings with alphagenome's merge helper."""

    if isinstance(intervals, (str, bytes)) or not isinstance(intervals, Sequence) or not intervals:
        raise ConfigurationError("intervals must be a non-empty list of interval strings.")
    if len(intervals) > _MAX_SEQUENCES:
        raise ConfigurationError(f"intervals may contain at most {_MAX_SEQUENCES} entries.")
    genome = require("alphagenome.data.genome", _EXTRAS[_ALPHAGENOME])
    parsed = []
    for index, text in enumerate(intervals):
        value = _text(text, f"intervals[{index}]")
        try:
            parsed.append(genome.Interval.from_str(value))
        except (ValueError, TypeError, IndexError) as exc:
            raise _rejected(_ALPHAGENOME, exc) from exc
    merged = genome.merge_overlapping_intervals(parsed)
    return {
        "status": "completed",
        "operation": "merge_alphagenome_intervals",
        "input_count": len(parsed),
        "merged_count": len(merged),
        "intervals": [_interval_record(item) for item in merged],
        "provenance": _provenance(_ALPHAGENOME, coordinate_convention="0-based half-open"),
    }


def read_gencode_gene_intervals(
    gtf_path: str | Path,
    *,
    gene_symbols: Sequence[str] | None = None,
    gene_ids: Sequence[str] | None = None,
    max_records: int = _DEFAULT_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return gene intervals from a GENCODE GTF by symbol or Ensembl gene ID."""

    path = input_path(gtf_path)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    if (gene_symbols is None) == (gene_ids is None):
        raise ConfigurationError("Pass exactly one of gene_symbols or gene_ids.")
    key = "gene_symbols" if gene_symbols is not None else "gene_ids"
    values = gene_symbols if gene_symbols is not None else gene_ids
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{key} must be a non-empty list of strings.")
    names = [_text(item, f"{key}[{index}]") for index, item in enumerate(values)]
    gene_annotation = require("alphagenome.data.gene_annotation", _EXTRAS[_ALPHAGENOME])
    frame = _gtf_frame(path, frozenset({"gene"}), limit)
    try:
        intervals = gene_annotation.get_gene_intervals(
            frame,
            gene_symbols=names if key == "gene_symbols" else None,
            gene_ids=names if key == "gene_ids" else None,
        )
    except (ValueError, KeyError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "read_gencode_gene_intervals",
        "source": str(path),
        "requested": names,
        "genes": [
            {**_interval_record(interval), "name": str(interval.name)} for interval in intervals
        ],
        "provenance": _provenance(
            _ALPHAGENOME,
            coordinate_convention="GTF 1-based start converted to 0-based; end is exclusive",
        ),
    }


def extract_gencode_transcripts(
    gtf_path: str | Path,
    region: str,
    *,
    max_records: int = _DEFAULT_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return GENCODE transcripts overlapping ``region`` with alphagenome's TranscriptExtractor."""

    path = input_path(gtf_path)
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    query_text = _text(region, "region")
    genome = require("alphagenome.data.genome", _EXTRAS[_ALPHAGENOME])
    transcript_module = require("alphagenome.data.transcript", _EXTRAS[_ALPHAGENOME])
    try:
        query = genome.Interval.from_str(query_text)
    except (ValueError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    frame = _gtf_frame(path, _TRANSCRIPT_FEATURES, limit)
    try:
        extractor = transcript_module.TranscriptExtractor(frame)
        transcripts = extractor.extract(query)
    except (ValueError, KeyError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    records = []
    for item in transcripts:
        records.append(
            {
                "transcript_id": item.transcript_id,
                "gene_id": item.gene_id,
                "protein_id": item.protein_id,
                "chromosome": str(item.chromosome),
                "strand": "-" if item.is_negative_strand else "+",
                "transcript_type": item.info.get("transcript_type"),
                "is_coding": bool(item.is_coding),
                "exons": [[int(exon.start), int(exon.end)] for exon in item.exons],
                "cds": [[int(part.start), int(part.end)] for part in (item.cds or [])],
            }
        )
    return {
        "status": "completed",
        "operation": "extract_gencode_transcripts",
        "source": str(path),
        "region": {
            "chromosome": str(query.chromosome),
            "start": int(query.start),
            "end": int(query.end),
        },
        "transcript_count": len(records),
        "transcripts": records,
        "provenance": _provenance(
            _ALPHAGENOME,
            coordinate_convention="0-based half-open; region uses the same convention",
        ),
    }


# grelu: DNA string utilities (torch-based upstream functions)


def reverse_complement_with_grelu(
    sequences: str | Sequence[str], *, verbose: bool = True
) -> dict[str, Any]:
    """Reverse-complement one DNA string or a list of DNA strings with grelu."""

    values = _dna_batch(sequences, "sequences")
    utils = require("grelu.sequence.utils", _EXTRAS[_GRELU])
    try:
        result = utils.reverse_complement(values)
    except (ValueError, TypeError, KeyError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "reverse_complement_with_grelu",
        "sequences": [str(item) for item in result],
        "provenance": _provenance(_GRELU, alphabet="ACGTN"),
    }


def shuffle_dinucleotides(
    sequence: str,
    *,
    n_shuffles: int = 1,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return dinucleotide-preserving shuffles of one DNA string, using grelu."""

    value = _dna(sequence, "sequence")
    count = _integer(n_shuffles, "n_shuffles", minimum=1, maximum=1_000)
    if seed is not None:
        seed = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    utils = require("grelu.sequence.utils", _EXTRAS[_GRELU])
    try:
        shuffles = utils.dinuc_shuffle(value, n_shuffles=count, seed=seed)
    except ValueError as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "shuffle_dinucleotides",
        "sequence": value,
        "n_shuffles": count,
        "shuffles": [str(item) for item in shuffles],
        "provenance": _provenance(_GRELU, method="grelu dinuc_shuffle (tangermeme backend)"),
    }


def resize_dna_sequences(
    sequences: str | Sequence[str],
    seq_len: int,
    *,
    end: str = "both",
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim or N-pad DNA strings to ``seq_len`` with grelu."""

    values = _dna_batch(sequences, "sequences")
    length = _integer(seq_len, "seq_len", minimum=1, maximum=_MAX_SEQUENCE)
    side = _choice(end, "end", ("left", "right", "both"))
    utils = require("grelu.sequence.utils", _EXTRAS[_GRELU])
    try:
        result = utils.resize(values, length, end=side)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "resize_dna_sequences",
        "seq_len": length,
        "end": side,
        "sequences": [str(item) for item in result],
        "provenance": _provenance(_GRELU, padding_base="N", trimming="symmetric when end is both"),
    }


def mutate_dna_sequence(
    sequence: str, allele: str, *, position: int, verbose: bool = True
) -> dict[str, Any]:
    """Replace one base at a 0-based ``position`` with ``allele``, using grelu."""

    value = _dna(sequence, "sequence")
    base = _dna(allele, "allele", maximum=1)
    index = _integer(position, "position", minimum=0, maximum=len(value) - 1)
    mutate_module = require("grelu.sequence.mutate", _EXTRAS[_GRELU])
    try:
        mutated = mutate_module.mutate(value, base, pos=index)
    except (ValueError, TypeError) as exc:
        raise _rejected(_GRELU, exc) from exc
    return {
        "status": "completed",
        "operation": "mutate_dna_sequence",
        "reference": value,
        "sequence": str(mutated),
        "position": index,
        "reference_base": value[index],
        "allele": base,
        "provenance": _provenance(_GRELU, position_convention="0-based"),
    }


# selene-sdk: FASTA-backed sequences and encodings


def fetch_genome_sequences(
    fasta_path: str | Path,
    regions: Sequence[Mapping[str, Any]],
    *,
    pad: bool = False,
    max_regions: int = 1_000,
    max_region_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch DNA for 0-based regions from an indexed FASTA, using selene-sdk.

    Each region is a mapping with ``chrom``, ``start``, ``end`` and an optional
    ``strand``. pyfaidx writes a ``.fai`` index next to the FASTA file.
    """

    path = input_path(fasta_path)
    if isinstance(regions, (str, bytes, Mapping)) or not isinstance(regions, Sequence):
        raise ConfigurationError("regions must be a list of mappings.")
    limit = _integer(max_regions, "max_regions", minimum=1, maximum=_MAX_REGIONS)
    span = _integer(max_region_length, "max_region_length", minimum=1, maximum=_MAX_SEQUENCE)
    if not regions or len(regions) > limit:
        raise ConfigurationError(f"regions must contain between 1 and {limit} entries.")
    parsed = [_region(item, index) for index, item in enumerate(regions)]
    for index, (_chrom, start, end, _strand) in enumerate(parsed):
        if end - start > span:
            raise ConfigurationError(f"regions[{index}] is longer than max_region_length={span}.")
    selene = require("selene_sdk.sequences", _EXTRAS[_SELENE])
    genome = selene.Genome(str(path))
    known = set(genome.get_chrs())
    results = []
    for chrom, start, end, strand in parsed:
        if chrom not in known:
            raise ConfigurationError(f"unknown chromosome {chrom!r} in {path.name}.")
        if not pad and not genome.coords_in_bounds(chrom, start, end):
            raise ConfigurationError(
                f"{chrom}:{start}-{end} is outside the FASTA; pass pad=True to pad with N."
            )
        sequence = genome.get_sequence_from_coords(chrom, start, end, strand=strand, pad=pad)
        results.append(
            {"chrom": chrom, "start": start, "end": end, "strand": strand, "sequence": sequence}
        )
    return {
        "status": "completed",
        "operation": "fetch_genome_sequences",
        "source": str(path),
        "regions": results,
        "provenance": _provenance(
            _SELENE,
            coordinate_convention="0-based half-open",
            padding="N" if pad else "none",
        ),
    }


def read_contig_lengths(
    fasta_path: str | Path, *, max_contigs: int = 100_000, verbose: bool = True
) -> dict[str, Any]:
    """List contig names and lengths of a FASTA using selene-sdk's Genome index."""

    path = input_path(fasta_path)
    limit = _integer(max_contigs, "max_contigs", minimum=1, maximum=_MAX_REGIONS)
    selene = require("selene_sdk.sequences", _EXTRAS[_SELENE])
    lengths = selene.Genome(str(path)).get_chr_lens()
    if len(lengths) > limit:
        raise ConfigurationError(f"FASTA has more than max_contigs={limit} contigs.")
    return {
        "status": "completed",
        "operation": "read_contig_lengths",
        "source": str(path),
        "contig_count": len(lengths),
        "contigs": [{"name": str(name), "length": int(length)} for name, length in lengths],
        "provenance": _provenance(_SELENE),
    }


def encode_one_hot_selene(
    sequence: str, *, max_length: int = 10_000, verbose: bool = True
) -> dict[str, Any]:
    """One-hot encode DNA with selene-sdk; N is uniform 0.25 in every channel."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_ONE_HOT)
    value = _dna(sequence, "sequence", maximum=limit)
    selene = require("selene_sdk.sequences", _EXTRAS[_SELENE])
    encoding = selene.Genome.sequence_to_encoding(value)
    return {
        "status": "completed",
        "operation": "encode_one_hot_selene",
        "length": len(value),
        "channels": list(_CHANNELS),
        "one_hot": _clean(encoding),
        "provenance": _provenance(_SELENE, unknown_base_encoding="0.25 in each channel"),
    }


# enformer-pytorch: one-hot encoding and BED-driven FASTA windows


def encode_one_hot_enformer(
    sequence: str, *, max_length: int = 10_000, verbose: bool = True
) -> dict[str, Any]:
    """One-hot encode DNA with enformer-pytorch; N is all zeros."""

    limit = _integer(max_length, "max_length", minimum=1, maximum=_MAX_ONE_HOT)
    value = _dna(sequence, "sequence", maximum=limit)
    data = require("enformer_pytorch.data", _EXTRAS[_ENFORMER])
    encoding = data.str_to_one_hot(value)
    return {
        "status": "completed",
        "operation": "encode_one_hot_enformer",
        "length": len(value),
        "channels": list(_CHANNELS),
        "one_hot": _clean(encoding),
        "provenance": _provenance(_ENFORMER, unknown_base_encoding="all zeros"),
    }


def read_enformer_bed_windows(
    bed_path: str | Path,
    fasta_path: str | Path,
    *,
    context_length: int = _ENFORMER_CONTEXT,
    max_windows: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize enformer-pytorch's one-hot windows for BED intervals in a FASTA.

    Windows shorter than ``context_length`` are extended symmetrically with padding,
    as in ``FastaInterval``. Longer windows keep their BED length. Only counts are
    returned, so large one-hot matrices are not copied into the result.
    """

    bed = input_path(bed_path)
    fasta = input_path(fasta_path)
    context = _integer(context_length, "context_length", minimum=1, maximum=_MAX_CONTEXT)
    limit = _integer(max_windows, "max_windows", minimum=1, maximum=_MAX_WINDOWS)
    data = require("enformer_pytorch.data", _EXTRAS[_ENFORMER])
    try:
        dataset = data.GenomeIntervalDataset(
            bed_file=str(bed),
            fasta_file=str(fasta),
            context_length=context,
        )
    except (ValueError, KeyError, OSError) as exc:
        raise _rejected(_ENFORMER, exc) from exc
    count = len(dataset)
    if count == 0 or count > limit:
        raise ConfigurationError(f"BED must contain between 1 and {limit} intervals; got {count}.")
    windows = []
    for index in range(count):
        chrom, start, end = dataset.df.row(index)[:3]
        try:
            one_hot = dataset[index]
        except (ValueError, KeyError) as exc:
            raise _rejected(_ENFORMER, exc) from exc
        padding = int((one_hot == 0.25).all(dim=-1).sum().item())
        unknown = int((one_hot == 0.0).all(dim=-1).sum().item())
        length = int(one_hot.shape[0])
        windows.append(
            {
                "index": index,
                "chromosome": str(chrom),
                "start": int(start),
                "end": int(end),
                "length": length,
                "padding_positions": padding,
                "unknown_positions": unknown,
                "acgt_positions": length - padding - unknown,
            }
        )
    return {
        "status": "completed",
        "operation": "read_enformer_bed_windows",
        "source": str(bed),
        "fasta": str(fasta),
        "context_length": context,
        "window_count": count,
        "windows": windows,
        "provenance": _provenance(
            _ENFORMER,
            padding="'.' rows are 0.25 in every channel",
            unknown_base="N rows are all zero",
        ),
    }


# genomic-benchmarks: local dataset metadata and cache state


def list_genomic_benchmark_datasets(
    *, dataset_path: str | Path | None = None, verbose: bool = True
) -> dict[str, Any]:
    """List dataset names from the package cache index, or dataset folders under a directory."""

    gb = require("genomic_benchmarks.data_check", _EXTRAS[_GENOMIC_BENCHMARKS])
    if dataset_path is None:
        names = sorted(str(name) for name in gb.list_datasets())
        source = "package_index"
        location = None
    else:
        folder = Path(_text(str(dataset_path), "dataset_path", maximum=4096)).expanduser()
        folder = folder.resolve()
        if not folder.is_dir():
            raise ConfigurationError(f"dataset_path is not a directory: {folder}")
        names = sorted(str(name) for name in gb.list_datasets(folder, local_repo=True))
        source = "local_directory"
        location = str(folder)
    return {
        "status": "completed",
        "operation": "list_genomic_benchmark_datasets",
        "source": source,
        "location": location,
        "count": len(names),
        "datasets": names,
        "provenance": _provenance(_GENOMIC_BENCHMARKS, network="not used for either source"),
    }


def describe_genomic_benchmark(
    dataset_path: str | Path,
    *,
    version: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize class and split counts of a local genomic-benchmarks dataset folder."""

    folder = Path(_text(str(dataset_path), "dataset_path", maximum=4096)).expanduser().resolve()
    if not folder.is_dir():
        raise ConfigurationError(f"dataset_path is not a directory: {folder}")
    expected = None if version is None else _integer(version, "version", minimum=0, maximum=10**6)
    gb = require("genomic_benchmarks.data_check", _EXTRAS[_GENOMIC_BENCHMARKS])
    printed = io.StringIO()
    try:
        with contextlib.redirect_stdout(printed), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = gb.info(str(folder), version=expected, local_repo=True)
    except (FileNotFoundError, ValueError, KeyError, OSError) as exc:
        raise ConfigurationError(f"could not describe {folder.name}: {exc}") from exc
    counts = {
        str(label): {str(split): int(value) for split, value in row.items()}
        for label, row in table.to_dict(orient="index").items()
    }
    return {
        "status": "completed",
        "operation": "describe_genomic_benchmark",
        "dataset": folder.name,
        "classes": counts,
        "total_sequences": int(table.to_numpy().sum()),
        "messages": [line for line in printed.getvalue().splitlines() if line.strip()],
        "provenance": _provenance(_GENOMIC_BENCHMARKS, version_checked=str(expected is not None)),
    }


def check_genomic_benchmark_cache(
    dataset: str,
    *,
    cache_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report whether a named genomic-benchmarks dataset exists in the local cache folder."""

    name = _text(dataset, "dataset", maximum=128)
    if not _NAME.match(name):
        raise ConfigurationError(
            "dataset must be a plain benchmark name such as 'human_enhancers_cohn'."
        )
    gb = require("genomic_benchmarks.data_check", _EXTRAS[_GENOMIC_BENCHMARKS])
    paths = require("genomic_benchmarks.utils.paths", _EXTRAS[_GENOMIC_BENCHMARKS])
    if cache_path is None:
        location = Path(paths.CACHE_PATH)
        downloaded = bool(gb.is_downloaded(name))
    else:
        location = Path(_text(str(cache_path), "cache_path", maximum=4096)).expanduser().resolve()
        downloaded = bool(gb.is_downloaded(name, cache_path=location))
    return {
        "status": "completed",
        "operation": "check_genomic_benchmark_cache",
        "dataset": name,
        "cache_path": str(location),
        "downloaded": downloaded,
        "provenance": _provenance(_GENOMIC_BENCHMARKS, network="not used"),
    }


# evo2: local generation (fake-tested only in this environment)


def generate_evo2_sequences(
    prompts: Sequence[str],
    *,
    weights_path: str | Path,
    model_name: str = "evo2_7b",
    n_tokens: int = 100,
    temperature: float = 1.0,
    top_k: int = 4,
    top_p: float = 1.0,
    max_prompts: int = 8,
    max_prompt_length: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate DNA continuations of prompts with a local Evo 2 checkpoint.

    ``weights_path`` must be an existing checkpoint file. The function never
    downloads weights, so the Hugging Face fallback of the upstream constructor
    is not used.
    """

    if isinstance(prompts, (str, bytes)) or not isinstance(prompts, Sequence) or not prompts:
        raise ConfigurationError("prompts must be a non-empty list of DNA strings.")
    limit = _integer(max_prompts, "max_prompts", minimum=1, maximum=64)
    if len(prompts) > limit:
        raise ConfigurationError(f"prompts may contain at most max_prompts={limit} entries.")
    values = [
        _dna(item, f"prompts[{index}]", maximum=max_prompt_length)
        for index, item in enumerate(prompts)
    ]
    weights = input_path(weights_path)
    model = _choice(model_name, "model_name", _EVO2_MODELS)
    tokens = _integer(n_tokens, "n_tokens", minimum=1, maximum=10_000)
    heat = _real(temperature, "temperature", low=0.0, high=10.0, low_open=True)
    top = _integer(top_k, "top_k", minimum=1, maximum=1_000)
    nucleus = _real(top_p, "top_p", low=0.0, high=1.0, low_open=True)
    evo2 = require("evo2", _EXTRAS[_EVO2])
    try:
        wrapper = evo2.Evo2(model_name=model, local_path=str(weights))
        sequences, scores = wrapper.generate(
            prompt_seqs=values,
            n_tokens=tokens,
            temperature=heat,
            top_k=top,
            top_p=nucleus,
            batched=True,
            verbose=0,
        )
    except Exception as exc:
        raise BackendExecutionError(
            "Evo 2 generation failed.",
            code="EVO2_GENERATION_FAILED",
            context={"error_type": type(exc).__name__},
        ) from exc
    return {
        "status": "completed",
        "operation": "generate_evo2_sequences",
        "model_name": model,
        "prompt_count": len(values),
        "n_tokens": tokens,
        "sequences": [str(item) for item in sequences],
        "scores": [None if score is None else _clean(float(score)) for score in scores],
        "provenance": _provenance(_EVO2, weights="explicit local file only"),
    }


# kipoi: model zoo metadata (network and git source; fake-tested here)


def list_kipoi_models(
    *,
    sources: Sequence[str] | None = None,
    limit: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List model descriptions from kipoi sources.

    Calling kipoi reads its git-backed model sources. The first call can clone a
    repository, and the default source list includes every configured source.
    """

    rows_limit = _integer(limit, "limit", minimum=1, maximum=10_000)
    names = None
    if sources is not None:
        if isinstance(sources, (str, bytes)) or not isinstance(sources, Sequence) or not sources:
            raise ConfigurationError("sources must be a non-empty list of source names.")
        names = [_text(item, f"sources[{index}]") for index, item in enumerate(sources)]
    kipoi = require("kipoi", _EXTRAS[_KIPOI])
    try:
        table = kipoi.list_models(sources=names)
    except (ValueError, KeyError, OSError) as exc:
        raise _rejected(_KIPOI, exc) from exc
    columns = [
        column
        for column in ("model", "version", "type", "license", "inputs", "targets")
        if column in table.columns
    ]
    models = _clean(table.head(rows_limit)[columns].to_dict(orient="records"))
    return {
        "status": "completed",
        "operation": "list_kipoi_models",
        "sources": names if names is not None else "all configured sources",
        "total_models": len(table),
        "returned": len(models),
        "models": models,
        "provenance": _provenance(_KIPOI, source_access="git-backed kipoi sources"),
    }


def describe_kipoi_model(
    model: str, *, source: str = "kipoi", verbose: bool = True
) -> dict[str, Any]:
    """Return the documented metadata and input and target names of one kipoi model."""

    name = _text(model, "model", maximum=512)
    source_name = _text(source, "source", maximum=128)
    kipoi = require("kipoi", _EXTRAS[_KIPOI])
    try:
        descr = kipoi.get_model_descr(name, source=source_name)
    except (ValueError, KeyError, OSError) as exc:
        raise _rejected(_KIPOI, exc) from exc
    info = descr.info
    return {
        "status": "completed",
        "operation": "describe_kipoi_model",
        "model": name,
        "source": source_name,
        "description": {
            "type": _clean(getattr(descr, "type", None)),
            "version": _clean(info.version),
            "doc": _clean(info.doc),
            "license": _clean(info.license),
            "cite_as": _clean(getattr(info, "cite_as", None)),
            "inputs": _names(descr.schema.inputs),
            "targets": _names(descr.schema.targets),
        },
        "provenance": _provenance(_KIPOI, source_access="git-backed kipoi sources"),
    }


__all__ = [
    "check_genomic_benchmark_cache",
    "describe_genomic_benchmark",
    "describe_kipoi_model",
    "encode_one_hot_enformer",
    "encode_one_hot_selene",
    "extract_gencode_transcripts",
    "fetch_genome_sequences",
    "generate_evo2_sequences",
    "list_genomic_benchmark_datasets",
    "list_kipoi_models",
    "merge_alphagenome_intervals",
    "mutate_dna_sequence",
    "parse_alphagenome_interval",
    "parse_alphagenome_variant",
    "read_contig_lengths",
    "read_enformer_bed_windows",
    "read_gencode_gene_intervals",
    "resize_dna_sequences",
    "reverse_complement_with_grelu",
    "shuffle_dinucleotides",
]

_register_public_verbose(globals(), tuple(__all__))
