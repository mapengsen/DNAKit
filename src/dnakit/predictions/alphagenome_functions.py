"""Offline AlphaGenome utilities for DNA sequences, intervals, genes and ISM.

Every public function imports ``alphagenome`` only when it runs and uses local inputs:
sequences, interval strings, GENCODE GTF files, Ensembl ID lists, fold BED files and
ISM scores. Model inference, API clients, gRPC services and plotting are not wrapped.
"""

from __future__ import annotations

import gzip
import itertools
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _dna,
    _integer,
    _interval_record,
    _provenance,
    _rejected,
    _text,
)

_ALPHAGENOME = "alphagenome"
_EXTRAS = "bio-alphagenome"
_GTF_FEATURES = frozenset({"gene", "transcript"})
_GTF_COLUMNS = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Feature",
    "gene_id",
    "gene_name",
    "gene_type",
    "transcript_id",
    "transcript_type",
    "transcript_support_level",
    "tag",
]
_GTF_ATTRIBUTE = re.compile(r'([A-Za-z_][\w.-]*) "([^"]*)"')
_STRANDS = ("+", "-", ".")
_BASES = frozenset("ACGTN")
_OPERATIONS = ("intersect", "union")
_TRANSCRIPT_METHODS = (
    "protein_coding",
    "transcript_type",
    "longest",
    "mane_select",
    "support_level",
)
_SUPPORT_LEVELS = ("1", "2", "3", "4", "5")
_FOLD_VERSIONS = ("all_folds", "fold_0", "fold_1", "fold_2", "fold_3")
_SUBSETS = ("train", "valid", "test")
_ORGANISMS = ("homo_sapiens", "mus_musculus")
_MAX_TEXT = 256
_MAX_ISM_WIDTH = 2_000
_MAX_INTERVALS = 10_000
_MAX_IDS = 100_000
_MAX_RECORDS = 5_000_000
_DEFAULT_RECORDS = 1_000_000
_MAX_COORDINATE = 10**12


class _WindowExtractor:
    """Serve reference bases from one local window to alphagenome's left-aligner."""

    def __init__(self, chromosome: str, start: int, sequence: str) -> None:
        self._chromosome = chromosome
        self._start = start
        self._sequence = sequence

    def extract(self, interval: Any) -> str:
        end = self._start + len(self._sequence)
        if str(interval.chromosome) != self._chromosome or not (
            self._start <= interval.start <= interval.end <= end
        ):
            raise ValueError("reference window does not cover the requested interval")
        return self._sequence[interval.start - self._start : interval.end - self._start]


def _variant_record(variant: Any) -> dict[str, Any]:
    return {
        "chromosome": str(variant.chromosome),
        "position": int(variant.position),
        "reference": str(variant.reference_bases),
        "alternate": str(variant.alternate_bases),
    }


def _parse_interval(genome: Any, text: str, name: str) -> Any:
    value = _text(text, name)
    try:
        return genome.Interval.from_str(value)
    except (ValueError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc


def _interval_list(genome: Any, values: object, name: str) -> list[Any]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of interval strings.")
    if len(values) > _MAX_INTERVALS:
        raise ConfigurationError(f"{name} may contain at most {_MAX_INTERVALS} intervals.")
    return [_parse_interval(genome, item, f"{name}[{index}]") for index, item in enumerate(values)]


def _require_disjoint(intervals: list[Any], name: str) -> None:
    """Alphagenome's set operations assume non-overlapping intervals within each side."""

    by_chromosome: dict[str, list[Any]] = {}
    for interval in intervals:
        by_chromosome.setdefault(str(interval.chromosome), []).append(interval)
    for chromosome, items in by_chromosome.items():
        ordered = sorted(items, key=lambda item: (item.start, item.end))
        for previous, current in itertools.pairwise(ordered):
            if current.start < previous.end:
                raise ConfigurationError(
                    f"{name} intervals must not overlap on {chromosome}: "
                    f"{previous.start}-{previous.end} and {current.start}-{current.end}."
                )


def _nullable(value: Any) -> str | None:
    if value is None or value != value:
        return None
    return str(value)


def _annotation_frame(path: Path, max_records: int) -> Any:
    """Read selected GENCODE rows with 0-based half-open coordinates and the fields
    used by alphagenome's transcript filters (tag and transcript_support_level)."""

    pandas = require("pandas", _EXTRAS)
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[dict[str, Any]] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ConfigurationError(f"GTF line {line_number} has fewer than 9 columns.")
            if fields[2] not in _GTF_FEATURES:
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
            values: dict[str, list[str]] = {}
            for key, value in _GTF_ATTRIBUTE.findall(fields[8]):
                values.setdefault(key, []).append(value)
            attributes = {key: " ".join(items) for key, items in values.items()}
            rows.append(
                {
                    "Chromosome": fields[0],
                    "Start": start,
                    "End": end,
                    "Strand": fields[6],
                    "Feature": fields[2],
                    "gene_id": attributes.get("gene_id"),
                    "gene_name": attributes.get("gene_name"),
                    "gene_type": attributes.get("gene_type"),
                    "transcript_id": attributes.get("transcript_id"),
                    "transcript_type": attributes.get("transcript_type"),
                    "transcript_support_level": attributes.get("transcript_support_level"),
                    "tag": attributes.get("tag"),
                }
            )
    return pandas.DataFrame(rows, columns=_GTF_COLUMNS)


def _gtf_path_and_limit(gtf_path: str | Path, max_records: int) -> tuple[Path, int]:
    path = input_path(gtf_path)
    return path, _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)


def normalize_alphagenome_variant(
    variant: str,
    *,
    reference_sequence: str,
    reference_start: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Left-align a variant with alphagenome's normalize_variant on a local reference window.

    ``reference_sequence`` holds the reference bases from 0-based ``reference_start``
    on the variant's chromosome. The window must include one base before the variant
    and the reference span of the variant.
    """

    value = _text(variant, "variant")
    window = _dna(reference_sequence, "reference_sequence")
    offset = _integer(reference_start, "reference_start", minimum=0, maximum=_MAX_COORDINATE)
    genome = require("alphagenome.data.genome", _EXTRAS)
    try:
        parsed = genome.Variant.from_str(value)
    except (ValueError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    extractor = _WindowExtractor(str(parsed.chromosome), offset, window)
    try:
        normalized = genome.normalize_variant(parsed, extractor)
    except (ValueError, IndexError, KeyError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    input_record = _variant_record(parsed)
    output_record = _variant_record(normalized)
    return {
        "status": "completed",
        "operation": "normalize_alphagenome_variant",
        "input": input_record,
        "normalized": output_record,
        "changed": input_record != output_record,
        "provenance": _provenance(
            _ALPHAGENOME,
            method="left-alignment (doi:10.1093/bioinformatics/btv112)",
            position_convention="1-based, as in VCF",
        ),
    }


def combine_alphagenome_intervals(
    left: Sequence[str],
    right: Sequence[str],
    *,
    operation: str = "intersect",
    verbose: bool = True,
) -> dict[str, Any]:
    """Intersect or union two interval sets with alphagenome's interval algebra.

    Each side must be non-overlapping. Output intervals are sorted by chromosome and
    coordinates; adjacent intervals are kept separate.
    """

    mode = _choice(operation, "operation", _OPERATIONS)
    genome = require("alphagenome.data.genome", _EXTRAS)
    lhs = _interval_list(genome, left, "left")
    rhs = _interval_list(genome, right, "right")
    _require_disjoint(lhs, "left")
    _require_disjoint(rhs, "right")
    function = genome.intersect_intervals if mode == "intersect" else genome.union_intervals
    try:
        combined = list(function(lhs, rhs))
    except (ValueError, TypeError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    records = sorted(
        (_interval_record(item) for item in combined),
        key=lambda item: (item["chromosome"], item["start"], item["end"]),
    )
    return {
        "status": "completed",
        "operation": "combine_alphagenome_intervals",
        "mode": mode,
        "left_count": len(lhs),
        "right_count": len(rhs),
        "count": len(records),
        "intervals": records,
        "provenance": _provenance(_ALPHAGENOME, coordinate_convention="0-based half-open"),
    }


def extract_alphagenome_tss(
    gtf_path: str | Path,
    *,
    feature: str = "transcript",
    max_records: int = _DEFAULT_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return transcription start sites from a GENCODE GTF as zero-width points.

    ``feature`` selects ``transcript`` or ``gene`` rows. TSSs use the 5' end of
    the feature: the start on ``+`` and the exclusive end on ``-``.
    """

    path, limit = _gtf_path_and_limit(gtf_path, max_records)
    choice = _choice(feature, "feature", ("transcript", "gene"))
    frame = _annotation_frame(path, limit)
    annotation = require("alphagenome.data.gene_annotation", _EXTRAS)
    try:
        tss = annotation.extract_tss(frame, feature=choice)
    except (ValueError, KeyError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    records = [
        {
            "chromosome": str(row.Chromosome),
            "position": int(row.Start),
            "strand": str(row.Strand),
            "gene_id": _nullable(row.gene_id),
            "gene_name": _nullable(row.gene_name),
            "transcript_id": _nullable(row.transcript_id),
        }
        for row in tss.itertuples(index=False)
    ]
    records.sort(
        key=lambda item: (
            item["chromosome"],
            item["position"],
            item["strand"],
            item["gene_id"] or "",
            item["transcript_id"] or "",
        )
    )
    return {
        "status": "completed",
        "operation": "extract_alphagenome_tss",
        "source": str(path),
        "feature": choice,
        "count": len(records),
        "tss": records,
        "provenance": _provenance(
            _ALPHAGENOME,
            coordinate_convention="TSS is a 0-based point; position is the point start",
        ),
    }


def filter_alphagenome_transcripts(
    gtf_path: str | Path,
    *,
    method: str = "protein_coding",
    transcript_types: Sequence[str] | None = None,
    support_levels: Sequence[str | int] | None = None,
    include_gene_entries: bool = False,
    max_records: int = _DEFAULT_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Filter a GENCODE GTF to protein-coding, typed, longest, MANE or TSL-selected rows.

    ``method`` selects one of alphagenome's filters: ``protein_coding``,
    ``transcript_type`` (needs ``transcript_types``), ``longest`` (longest transcript
    per gene), ``mane_select`` (needs a MANE tag, human only) or ``support_level``
    (needs ``support_levels`` from 1 to 5).
    """

    path, limit = _gtf_path_and_limit(gtf_path, max_records)
    choice = _choice(method, "method", _TRANSCRIPT_METHODS)
    if not isinstance(include_gene_entries, bool):
        raise ConfigurationError("include_gene_entries must be boolean.")
    frame = _annotation_frame(path, limit)
    annotation = require("alphagenome.data.gene_annotation", _EXTRAS)
    try:
        if choice == "protein_coding":
            result = annotation.filter_protein_coding(
                frame, include_gene_entries=include_gene_entries
            )
        elif choice == "transcript_type":
            result = annotation.filter_transcript_type(
                frame, _transcript_types(annotation, transcript_types)
            )
        elif choice == "longest":
            result = annotation.filter_to_longest_transcript(frame)
        elif choice == "mane_select":
            result = annotation.filter_to_mane_select_transcript(frame)
        else:
            result = annotation.filter_transcript_support_level(
                frame, _support_levels(support_levels)
            )
    except (ValueError, KeyError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    feature_counts = result["Feature"].value_counts()
    return {
        "status": "completed",
        "operation": "filter_alphagenome_transcripts",
        "source": str(path),
        "method": choice,
        "row_count": len(result),
        "feature_counts": {
            str(key): int(feature_counts[key]) for key in sorted(feature_counts.index)
        },
        "gene_ids": sorted({str(item) for item in result["gene_id"].dropna()}),
        "transcript_ids": sorted({str(item) for item in result["transcript_id"].dropna()}),
        "provenance": _provenance(_ALPHAGENOME, annotation="GENCODE GTF"),
    }


def _transcript_types(annotation: Any, values: object) -> tuple[Any, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("transcript_types must be a non-empty list of type names.")
    known = {item.value: item for item in annotation.TranscriptType}
    selected = []
    for item in values:
        if item not in known:
            raise ConfigurationError(f"Unknown transcript type: {item!r}.")
        selected.append(known[item])
    return tuple(selected)


def _support_levels(values: object) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("support_levels must be a non-empty list of levels 1 to 5.")
    levels = []
    for item in values:
        text = str(item) if isinstance(item, (int, str)) and not isinstance(item, bool) else ""
        if text not in _SUPPORT_LEVELS:
            raise ConfigurationError(f"Unknown transcript support level: {item!r}.")
        levels.append(text)
    return levels


def upgrade_alphagenome_annotation_ids(
    old_ids: Sequence[str],
    new_ids: Sequence[str],
    *,
    patchless: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map old Ensembl IDs to new IDs with alphagenome's version-aware matcher.

    Versions are dropped for matching, so ``ENST00010.1`` maps to ``ENST00010.3``.
    Unmatched old IDs map to ``None``. Set ``patchless`` when old IDs lack a version.
    """

    old = _id_list(old_ids, "old_ids")
    new = _id_list(new_ids, "new_ids")
    if not isinstance(patchless, bool):
        raise ConfigurationError("patchless must be boolean.")
    pandas = require("pandas", _EXTRAS)
    annotation = require("alphagenome.data.gene_annotation", _EXTRAS)
    try:
        mapped = annotation.upgrade_annotation_ids(
            pandas.Series(old, dtype=object),
            pandas.Series(new, dtype=object),
            patchless=patchless,
        )
    except (ValueError, AssertionError, KeyError, AttributeError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    mappings = [
        {"old_id": old_id, "new_id": _nullable(new_id)}
        for old_id, new_id in zip(old, mapped.tolist(), strict=True)
    ]
    matched = sum(1 for item in mappings if item["new_id"] is not None)
    return {
        "status": "completed",
        "operation": "upgrade_alphagenome_annotation_ids",
        "count": len(mappings),
        "matched_count": matched,
        "mappings": mappings,
        "provenance": _provenance(_ALPHAGENOME, patchless=str(patchless).lower()),
    }


def _id_list(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of identifiers.")
    if len(values) > _MAX_IDS:
        raise ConfigurationError(f"{name} may contain at most {_MAX_IDS} identifiers.")
    return [_text(item, f"{name}[{index}]") for index, item in enumerate(values)]


def enumerate_alphagenome_ism_variants(
    interval: str,
    sequence: str,
    *,
    vocabulary: str = "ACGT",
    skip_n: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List every single-base alternative for an interval, as alphagenome's ISM does.

    Variants are ordered by position, then by the ``vocabulary`` order; reference bases
    are skipped. ``position`` is 1-based, as in VCF.
    """

    genome = require("alphagenome.data.genome", _EXTRAS)
    ism = require("alphagenome.interpretation.ism", _EXTRAS)
    region = _parse_interval(genome, interval, "interval")
    bases = _dna(sequence, "sequence", maximum=_MAX_ISM_WIDTH)
    alphabet = _vocabulary(vocabulary)
    if len(bases) != region.width:
        raise ConfigurationError("sequence length must equal the interval width.")
    if not isinstance(skip_n, bool):
        raise ConfigurationError("skip_n must be boolean.")
    try:
        variants = ism.ism_variants(region, bases, vocabulary=alphabet, skip_n=skip_n)
    except ValueError as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    records = [_variant_record(item) for item in variants]
    return {
        "status": "completed",
        "operation": "enumerate_alphagenome_ism_variants",
        "interval": _interval_record(region),
        "vocabulary": alphabet,
        "count": len(records),
        "variants": records,
        "provenance": _provenance(
            _ALPHAGENOME,
            position_convention="1-based, as in VCF; interval start is 0-based",
        ),
    }


def _vocabulary(value: object) -> str:
    text = _text(value, "vocabulary", maximum=_MAX_TEXT).upper()
    if len(set(text)) != len(text) or set(text) - _BASES:
        raise ConfigurationError("vocabulary must contain unique letters from A, C, G, T and N.")
    return text


def build_alphagenome_ism_matrix(
    interval: str,
    sequence: str,
    scores: Sequence[float],
    *,
    vocabulary: str = "ACGT",
    multiply_by_sequence: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build alphagenome's ISM (position, base) matrix from one score per ISM variant.

    ``scores`` follows the order of ``enumerate_alphagenome_ism_variants`` for the same
    interval and sequence, and the sequence must contain only A, C, G and T. The matrix
    is centred per position; with ``multiply_by_sequence`` only reference-base entries
    stay non-zero.
    """

    genome = require("alphagenome.data.genome", _EXTRAS)
    ism = require("alphagenome.interpretation.ism", _EXTRAS)
    region = _parse_interval(genome, interval, "interval")
    bases = _dna(sequence, "sequence", maximum=_MAX_ISM_WIDTH)
    alphabet = _vocabulary(vocabulary)
    if len(alphabet) != 4 or set(bases) - set("ACGT"):
        raise ConfigurationError("ISM matrices need a 4-letter vocabulary and an A/C/G/T sequence.")
    if len(bases) != region.width:
        raise ConfigurationError("sequence length must equal the interval width.")
    if not isinstance(multiply_by_sequence, bool):
        raise ConfigurationError("multiply_by_sequence must be boolean.")
    values = _scores(scores)
    variants = ism.ism_variants(region, bases, vocabulary=alphabet, skip_n=False)
    if len(values) != len(variants):
        raise ConfigurationError(
            f"scores must have one value per ISM variant: expected {len(variants)}, "
            f"got {len(values)}."
        )
    try:
        matrix = ism.ism_matrix(
            values,
            variants,
            interval=region,
            multiply_by_sequence=multiply_by_sequence,
            vocabulary=alphabet,
        )
    except ValueError as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    # Adding 0.0 turns negative zeros into 0.0 so the JSON output is stable.
    rows = [[float(item) + 0.0 for item in row] for row in matrix.tolist()]
    return {
        "status": "completed",
        "operation": "build_alphagenome_ism_matrix",
        "interval": _interval_record(region),
        "columns": list(alphabet),
        "matrix": rows,
        "provenance": _provenance(
            _ALPHAGENOME,
            centring="score minus the per-position mean over the vocabulary",
        ),
    }


def _scores(values: object) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("scores must be a non-empty list of numbers.")
    if len(values) > 3 * _MAX_ISM_WIDTH:
        raise ConfigurationError("scores may contain at most 3 * 2000 values.")
    result = []
    for index, item in enumerate(values):
        if isinstance(item, bool) or not isinstance(item, (int, float)):
            raise ConfigurationError(f"scores[{index}] must be a finite number.")
        number = float(item)
        if number != number or number in (float("inf"), float("-inf")):
            raise ConfigurationError(f"scores[{index}] must be a finite number.")
        result.append(number)
    return result


def validate_alphagenome_sequence_length(
    sequence_length: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check that a sequence length is one alphagenome's model accepts.

    Raises ``ConfigurationError`` for lengths outside the supported set (16 KB,
    100 KB, 500 KB and 1 MB). No model is loaded and no request is sent.
    """

    length = _integer(sequence_length, "sequence_length", minimum=1, maximum=_MAX_COORDINATE)
    validator = require("alphagenome.models.dna_client", _EXTRAS)
    model_client = require("alphagenome.models.v1.dna_client", _EXTRAS)
    try:
        validator.validate_sequence_length(length)
    except ValueError as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    supported = sorted(int(item) for item in model_client.SUPPORTED_SEQUENCE_LENGTHS.values())
    return {
        "status": "completed",
        "operation": "validate_alphagenome_sequence_length",
        "sequence_length": length,
        "supported_lengths": supported,
        "provenance": _provenance(_ALPHAGENOME, scope="length check only; no model call"),
    }


def get_alphagenome_fold_intervals(
    regions_path: str | Path,
    *,
    model_version: str = "fold_0",
    subset: str = "test",
    organism: str = "homo_sapiens",
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the training, validation or test intervals of an alphagenome data fold.

    ``regions_path`` is a local tab-separated file with chromosome, start, end and fold
    columns and no header. The default regions download from GitHub is not used here.
    """

    path = input_path(regions_path)
    version_name = _choice(model_version, "model_version", _FOLD_VERSIONS)
    subset_name = _choice(subset, "subset", _SUBSETS)
    organism_name = _choice(organism, "organism", _ORGANISMS)
    folds = require("alphagenome.data.fold_intervals", _EXTRAS)
    model = require("alphagenome.models.v1.dna_model", _EXTRAS)
    try:
        version = model.ModelVersion[version_name.upper()]
        species = model.Organism[organism_name.upper()]
        subset_value = folds.Subset[subset_name.upper()]
        frame = folds.get_fold_intervals(
            version, species, subset_value, example_regions_path=str(path)
        )
        fold_names = folds.get_fold_names(version, subset_value)
    except (ValueError, KeyError, OSError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    records = sorted(
        (
            {
                "chromosome": str(row.chromosome),
                "start": int(row.start),
                "end": int(row.end),
                "fold": str(row.fold),
            }
            for row in frame.itertuples(index=False)
        ),
        key=lambda item: (item["chromosome"], item["start"], item["end"], item["fold"]),
    )
    return {
        "status": "completed",
        "operation": "get_alphagenome_fold_intervals",
        "source": str(path),
        "model_version": version_name,
        "subset": subset_name,
        "organism": organism_name,
        "fold_names": sorted(str(name) for name in fold_names),
        "count": len(records),
        "intervals": records,
        "provenance": _provenance(_ALPHAGENOME, coordinate_convention="0-based half-open"),
    }


__all__ = [
    "build_alphagenome_ism_matrix",
    "combine_alphagenome_intervals",
    "enumerate_alphagenome_ism_variants",
    "extract_alphagenome_tss",
    "filter_alphagenome_transcripts",
    "get_alphagenome_fold_intervals",
    "normalize_alphagenome_variant",
    "upgrade_alphagenome_annotation_ids",
    "validate_alphagenome_sequence_length",
]

_register_public_verbose(globals(), tuple(__all__))
