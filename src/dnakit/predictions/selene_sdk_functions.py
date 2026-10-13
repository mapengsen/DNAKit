"""Offline selene-sdk helpers for one-hot DNA, mutants, VCF variants and attribution scores.

Each public function imports ``selene_sdk`` only when it runs. Inputs are local: one-hot
rows as nested lists in A, C, G, T order, DNA strings, a local FASTA with its VCF, and
score or variant tables on disk. Model prediction, training, configuration files, logging
and plotting are not wrapped here.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _clean,
    _dna,
    _integer,
    _provenance,
    _rejected,
    _text,
)

_SELENE = "selene-sdk"
_EXTRAS = "bio-selene"
_BASES = ("A", "C", "G", "T")
_UNKNOWN_ROW = 0.25
_SCALINGS = ("identity", "probability", "max_effect")
_MAX_ROWS = 100_000
_MAX_MUTATIONS = 1_000
_MAX_MUTANT_SEQUENCE = 10_000
_MAX_MUTANTS = 1_000_000
_MAX_VCF_BYTES = 50_000_000
_MAX_VARIANTS = 1_000_000
_MAX_INTERVALS = 1_000_000
_MAX_INTERVAL_LENGTH = 10**12
_MAX_SCORE_FILE_BYTES = 100_000_000
_MAX_LABELS = 200_000
_MAX_CHROMOSOMES = 10_000
_MAX_ALLELE_TEXT = 256


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _is_sequence(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, Mapping))


def _rows(values: object, name: str, *, limit: int = _MAX_ROWS) -> list[list[float]]:
    """Validate a non-empty list of 4-value rows and return them as floats."""

    if not _is_sequence(values) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of rows.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} may contain at most {limit} rows.")
    rows: list[list[float]] = []
    for index, row in enumerate(values):
        label = f"{name}[{index}]"
        if not _is_sequence(row) or len(row) != 4:
            raise ConfigurationError(f"{label} must hold four numbers for A, C, G and T.")
        rows.append([_number(item, label) for item in row])
    return rows


def _encoding(values: object, name: str) -> list[list[float]]:
    """Validate one-hot rows: a single 1 and three 0s, or a uniform 0.25 (unknown) row."""

    rows = _rows(values, name)
    for index, row in enumerate(rows):
        _row_base(row, f"{name}[{index}]")
    return rows


def _row_base(row: list[float], label: str) -> int:
    """Index of the base in a one-hot row, or -1 for a uniform 0.25 row."""

    if all(math.isclose(value, _UNKNOWN_ROW, abs_tol=1e-9) for value in row):
        return -1
    ones = [index for index, value in enumerate(row) if value == 1.0]
    zeros = [value for index, value in enumerate(row) if index not in ones]
    if len(ones) == 1 and all(value == 0.0 for value in zeros):
        return ones[0]
    raise ConfigurationError(f"{label} must be one-hot (one 1 and three 0s) or uniform 0.25.")


def _array(rows: list[list[float]]) -> Any:
    numpy = require("numpy", _EXTRAS)
    return numpy.asarray(rows, dtype=float)


def _genome_class() -> Any:
    return require("selene_sdk.sequences.genome", _EXTRAS).Genome


def _complement(genome: Any) -> dict[str, str]:
    return dict(genome.COMPLEMENTARY_BASE_DICT)


def _pairs(values: object, name: str, *, length: int) -> list[tuple[int, str]]:
    """Validate [position, base] pairs against a sequence of ``length`` positions."""

    if not _is_sequence(values):
        raise ConfigurationError(f"{name} must be a list of [position, base] pairs.")
    if len(values) > _MAX_MUTATIONS:
        raise ConfigurationError(f"{name} may contain at most {_MAX_MUTATIONS} pairs.")
    pairs: list[tuple[int, str]] = []
    seen: set[int] = set()
    for index, item in enumerate(values):
        label = f"{name}[{index}]"
        if not _is_sequence(item) or len(item) != 2:
            raise ConfigurationError(f"{label} must be a [position, base] pair.")
        position = _integer(item[0], f"{label}[0]", minimum=0, maximum=length - 1)
        if position in seen:
            raise ConfigurationError(f"{label} repeats position {position}.")
        seen.add(position)
        base = item[1].upper() if isinstance(item[1], str) else item[1]
        pairs.append((position, _choice(base, f"{label}[1]", _BASES)))
    return pairs


def _combination_count(counts: Sequence[int], choose: int) -> int:
    """Sum over all choose-subsets of the product of counts (elementary symmetric sum)."""

    table = [1] + [0] * choose
    for count in counts:
        for k in range(choose, 0, -1):
            table[k] += table[k - 1] * count
    return table[choose]


def decode_one_hot_sequence(
    encoding: Sequence[Sequence[float]],
    *,
    unk_base: str = "N",
    verbose: bool = True,
) -> dict[str, Any]:
    """Decode (L, 4) one-hot rows in A, C, G, T order to DNA; uniform 0.25 rows give unk_base."""

    rows = _encoding(encoding, "encoding")
    unknown = _choice(unk_base, "unk_base", (*_BASES, "N"))
    sequence_module = require("selene_sdk.sequences.sequence", _EXTRAS)
    genome = _genome_class()
    sequence = sequence_module.encoding_to_sequence(_array(rows), list(genome.BASES_ARR), unknown)
    return {
        "status": "completed",
        "operation": "decode_one_hot_sequence",
        "length": len(sequence),
        "sequence": sequence,
        "unk_base": unknown,
        "provenance": _provenance(_SELENE, channel_order="A,C,G,T", unknown_row="uniform 0.25"),
    }


def reverse_complement_one_hot(
    encoding: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement (L, 4) one-hot rows: permute channels A<->T and C<->G, then flip rows."""

    rows = _encoding(encoding, "encoding")
    model_predict = require("selene_sdk.predict.model_predict", _EXTRAS)
    genome = _genome_class()
    encoded = model_predict.get_reverse_complement_encoding(
        _array(rows), list(genome.BASES_ARR), _complement(genome)
    )
    return {
        "status": "completed",
        "operation": "reverse_complement_one_hot",
        "length": len(rows),
        "encoding": _clean(encoded),
        "provenance": _provenance(_SELENE, complement="A-T, C-G, N-N"),
    }


def reverse_complement_allele(allele: str, *, verbose: bool = True) -> dict[str, Any]:
    """Reverse-complement a variant allele; '-' and '*' (deletion or missing) stay '*'."""

    text = _text(allele, "allele", maximum=_MAX_ALLELE_TEXT).upper()
    if text not in ("-", "*") and (set(text) - set("ACGTN") or "-" in text or "*" in text):
        raise ConfigurationError(
            "allele must use A, C, G, T or N only, or be a single '-' or '*' allele."
        )
    model_predict = require("selene_sdk.predict.model_predict", _EXTRAS)
    result = model_predict.get_reverse_complement(text, _complement(_genome_class()))
    return {
        "status": "completed",
        "operation": "reverse_complement_allele",
        "allele": text,
        "reverse_complement": str(result),
        "provenance": _provenance(_SELENE, convention="'-' and '*' map to '*'"),
    }


def apply_one_hot_mutations(
    encoding: Sequence[Sequence[float]],
    mutations: Sequence[Sequence[Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Set given 0-based positions of (L, 4) one-hot rows to new bases (A, C, G or T)."""

    rows = _encoding(encoding, "encoding")
    pairs = _pairs(mutations, "mutations", length=len(rows))
    model_predict = require("selene_sdk.predict.model_predict", _EXTRAS)
    mutated = model_predict.mutate_sequence(_array(rows), pairs, _genome_class())
    return {
        "status": "completed",
        "operation": "apply_one_hot_mutations",
        "length": len(rows),
        "mutation_count": len(pairs),
        "mutations": [{"position": position, "base": base} for position, base in sorted(pairs)],
        "encoding": _clean(mutated),
        "provenance": _provenance(_SELENE, position_convention="0-based"),
    }


def list_in_silico_mutants(
    sequence: str,
    *,
    mutate_n_bases: int = 1,
    start_position: int = 0,
    end_position: int | None = None,
    max_mutants: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """All mutants changing exactly ``mutate_n_bases`` bases in a 0-based window of a DNA string."""

    value = _dna(sequence, "sequence", maximum=_MAX_MUTANT_SEQUENCE)
    width = _integer(mutate_n_bases, "mutate_n_bases", minimum=1, maximum=4)
    start = _integer(start_position, "start_position", minimum=0, maximum=len(value) - 1)
    if end_position is None:
        end = len(value)
    else:
        end = _integer(end_position, "end_position", minimum=1, maximum=len(value))
    if end <= start:
        raise ConfigurationError("end_position must be greater than start_position.")
    if end - start < width:
        raise ConfigurationError("the window is shorter than mutate_n_bases.")
    limit = _integer(max_mutants, "max_mutants", minimum=1, maximum=_MAX_MUTANTS)
    # selene-sdk offers the three other bases at each ACGT position and all four at N.
    options = [4 if base == "N" else 3 for base in value[start:end]]
    expected = _combination_count(options, width)
    if expected > limit:
        raise ConfigurationError(
            f"the window would give {expected} mutants; raise max_mutants or narrow the window."
        )
    model_predict = require("selene_sdk.predict.model_predict", _EXTRAS)
    mutants = model_predict.in_silico_mutagenesis_sequences(
        value,
        mutate_n_bases=width,
        reference_sequence=_genome_class(),
        start_position=start,
        end_position=end,
    )
    return {
        "status": "completed",
        "operation": "list_in_silico_mutants",
        "length": len(value),
        "start_position": start,
        "end_position": end,
        "mutate_n_bases": width,
        "count": len(mutants),
        "mutants": [
            [[int(position), str(base)] for position, base in mutant] for mutant in mutants
        ],
        "provenance": _provenance(_SELENE, order="combinations of positions, then base order"),
    }


def read_selene_vcf(
    vcf_path: str,
    fasta_path: str,
    *,
    max_variants: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read VCF records with selene-sdk and check chromosome, bounds and REF against the FASTA."""

    vcf = input_path(vcf_path)
    fasta = input_path(fasta_path)
    if vcf.stat().st_size > _MAX_VCF_BYTES:
        raise ConfigurationError(f"VCF file is larger than {_MAX_VCF_BYTES} bytes.")
    limit = _integer(max_variants, "max_variants", minimum=1, maximum=_MAX_VARIANTS)
    model_predict = require("selene_sdk.predict.model_predict", _EXTRAS)
    genome = _genome_class()(str(fasta))
    chromosomes = set(genome.get_chrs())
    try:
        records = model_predict.read_vcf_file(str(vcf), reference_sequence=genome)
    except (ValueError, KeyError, IndexError, OSError) as exc:
        raise _rejected(_SELENE, exc) from exc
    if len(records) > limit:
        raise ConfigurationError(f"VCF has more than max_variants={limit} records.")
    variants = []
    for chrom, position, name, ref, alt, strand in records:
        if chrom not in chromosomes:
            raise ConfigurationError(f"VCF chromosome {chrom!r} is not in {fasta.name}.")
        start = int(position) - 1
        # An insertion (empty REF) is checked at its anchor base, since empty intervals are invalid.
        in_bounds = bool(genome.coords_in_bounds(chrom, start, start + max(len(ref), 1)))
        matches = None
        if ref and in_bounds:
            found = genome.get_sequence_from_coords(chrom, start, start + len(ref))
            matches = str(found).upper() == ref.upper()
        variants.append(
            {
                "chromosome": str(chrom),
                "position": int(position),
                "id": str(name),
                "ref": str(ref),
                "alt": str(alt),
                "strand": str(strand),
                "in_bounds": in_bounds,
                "ref_matches_fasta": matches,
            }
        )
    return {
        "status": "completed",
        "operation": "read_selene_vcf",
        "fasta": fasta.name,
        "variant_count": len(variants),
        "variants": variants,
        "provenance": _provenance(
            _SELENE, coordinates="1-based VCF POS", ref_check="reported per variant"
        ),
    }


def unpack_bit_packed_sequence(
    packed_rows: Sequence[Sequence[int]],
    sequence_length: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unpack selene's bit-packed one-hot rows (uint8, 4 channels) to the first L positions."""

    rows = _rows(packed_rows, "packed_rows")
    for index, row in enumerate(rows):
        for value in row:
            if not value.is_integer() or not 0 <= value <= 255:
                raise ConfigurationError(f"packed_rows[{index}] must hold integers from 0 to 255.")
    length = _integer(sequence_length, "sequence_length", minimum=1, maximum=_MAX_ROWS * 8)
    if len(rows) < (length + 7) // 8:
        raise ConfigurationError("packed_rows holds fewer than ceil(sequence_length / 8) rows.")
    numpy = require("numpy", _EXTRAS)
    sequence_module = require("selene_sdk.samplers.dataloader", _EXTRAS)
    packed = numpy.asarray(rows, dtype=numpy.uint8)
    unpacked = sequence_module.unpackbits_sequence(packed, length)
    return {
        "status": "completed",
        "operation": "unpack_bit_packed_sequence",
        "length": length,
        "encoding": _clean(unpacked),
        "provenance": _provenance(
            _SELENE, unknown_row="rows with all four bits set become uniform 0.25"
        ),
    }


def interval_length_probabilities(
    interval_lengths: Sequence[int],
    *,
    indices: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sampling probabilities of intervals proportional to length; zero-length ones drop out."""

    if not _is_sequence(interval_lengths) or not interval_lengths:
        raise ConfigurationError("interval_lengths must be a non-empty list of integers.")
    if len(interval_lengths) > _MAX_INTERVALS:
        raise ConfigurationError(f"interval_lengths may contain at most {_MAX_INTERVALS} items.")
    lengths = [
        _integer(item, f"interval_lengths[{index}]", minimum=0, maximum=_MAX_INTERVAL_LENGTH)
        for index, item in enumerate(interval_lengths)
    ]
    if indices is None:
        selected = list(range(len(lengths)))
    else:
        if not _is_sequence(indices) or not indices:
            raise ConfigurationError("indices must be a non-empty list of interval indices.")
        selected = [
            _integer(item, f"indices[{index}]", minimum=0, maximum=len(lengths) - 1)
            for index, item in enumerate(indices)
        ]
        if len(set(selected)) != len(selected):
            raise ConfigurationError("indices must not repeat.")
    if sum(lengths[index] for index in selected) <= 0:
        raise ConfigurationError("at least one selected interval must have a positive length.")
    utils = require("selene_sdk.utils", _EXTRAS)
    kept, weights = utils.get_indices_and_probabilities(lengths, selected)
    kept_list = [int(index) for index in kept]
    kept_set = set(kept_list)
    return {
        "status": "completed",
        "operation": "interval_length_probabilities",
        "interval_count": len(lengths),
        "indices": kept_list,
        "probabilities": [float(weight) for weight in weights],
        "dropped_indices": [index for index in selected if index not in kept_set],
        "provenance": _provenance(_SELENE, weighting="proportional to interval length"),
    }


def rescale_attribution_scores(
    score_matrix: Sequence[Sequence[float]],
    *,
    base_scaling: str = "identity",
    position_scaling: str = "identity",
    verbose: bool = True,
) -> dict[str, Any]:
    """Rescale an (L, 4) attribution matrix by base and by position, as selene-sdk does."""

    rows = _rows(score_matrix, "score_matrix")
    base_mode = _choice(base_scaling, "base_scaling", _SCALINGS)
    position_mode = _choice(position_scaling, "position_scaling", _SCALINGS)
    interpret = require("selene_sdk.interpret", _EXTRAS)
    numpy = require("numpy", _EXTRAS)
    scores = _array(rows)
    with warnings.catch_warnings():
        # Zero sums are reported below as a ConfigurationError, not as numpy warnings.
        warnings.simplefilter("ignore", RuntimeWarning)
        rescaled = interpret.rescale_score_matrix(
            scores, base_scaling=base_mode, position_scaling=position_mode
        )
    bad = [index for index, bad_row in enumerate(~numpy.isfinite(rescaled).all(axis=1)) if bad_row]
    if bad:
        raise ConfigurationError(
            f"rescaling divides by zero at position(s) {bad[:5]}; check the score matrix."
        )
    return {
        "status": "completed",
        "operation": "rescale_attribution_scores",
        "length": len(rows),
        "base_scaling": base_mode,
        "position_scaling": position_mode,
        "scores": _clean(rescaled),
        "provenance": _provenance(_SELENE, note="base 'probability' is the identity in selene-sdk"),
    }


def load_variant_diff_scores(
    scores_path: str,
    *,
    max_rows: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a variant-by-feature score table: five variant columns, then one column per feature."""

    path = input_path(scores_path)
    if path.stat().st_size > _MAX_SCORE_FILE_BYTES:
        raise ConfigurationError(f"score file is larger than {_MAX_SCORE_FILE_BYTES} bytes.")
    limit = _integer(max_rows, "max_rows", minimum=1, maximum=1_000_000)
    interpret = require("selene_sdk.interpret", _EXTRAS)
    try:
        scores, labels, features = interpret.load_variant_abs_diff_scores(str(path))
    except (ValueError, IndexError) as exc:
        raise _rejected(_SELENE, exc) from exc
    if not features:
        raise ConfigurationError(
            "the score file has no feature columns after the five variant columns."
        )
    if len(labels) > limit:
        raise ConfigurationError(f"the score file has more than max_rows={limit} variants.")
    variants = []
    for label in labels:
        chrom, position, name, ref, alt = label
        if not str(position).isdigit():
            raise ConfigurationError(f"variant position {position!r} is not an integer.")
        variants.append(
            {
                "chromosome": str(chrom),
                "position": int(position),
                "name": str(name),
                "ref": str(ref),
                "alt": str(alt),
            }
        )
    return {
        "status": "completed",
        "operation": "load_variant_diff_scores",
        "variant_count": len(variants),
        "feature_count": len(features),
        "features": [str(feature) for feature in features],
        "variants": variants,
        "scores": _clean(scores),
        "provenance": _provenance(_SELENE, columns="5 variant columns then features"),
    }


def order_variant_labels(
    labels: Sequence[Sequence[Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Order variants by chromosome (1-22, X, Y, M, other), then position, ref and alt."""

    if not _is_sequence(labels) or not labels:
        raise ConfigurationError("labels must be a non-empty list of [chrom, pos, name, ref, alt].")
    if len(labels) > _MAX_LABELS:
        raise ConfigurationError(f"labels may contain at most {_MAX_LABELS} items.")
    records = []
    for index, item in enumerate(labels):
        label = f"labels[{index}]"
        if not _is_sequence(item) or len(item) != 5:
            raise ConfigurationError(f"{label} must be [chrom, pos, name, ref, alt].")
        chrom = _text(item[0], f"{label}.chrom", maximum=128)
        if not chrom.startswith("chr") or len(chrom) < 4:
            raise ConfigurationError(f"{label}.chrom must start with 'chr'.")
        position = _integer(item[1], f"{label}.pos", minimum=0, maximum=10**12)
        name, ref, alt = (str(part) for part in item[2:])
        if "\x00" in name + ref + alt or len(ref + alt) > 2_000:
            raise ConfigurationError(f"{label} has an invalid name, ref or alt.")
        records.append((index, chrom, position, name, ref, alt))
    interpret = require("selene_sdk.interpret", _EXTRAS)
    # Pre-sort so that ties resolve the same way for any input order.
    presorted = sorted(
        records,
        key=lambda rec: (
            interpret.sort_standard_chrs(rec[1]),
            rec[1],
            rec[2],
            rec[4],
            rec[5],
            rec[3],
            rec[0],
        ),
    )
    tuples = [(rec[1], rec[2], rec[3], rec[4], rec[5]) for rec in presorted]
    ordered, positions = interpret.ordered_variants_and_indices(tuples)
    original = [presorted[pos][0] for pos in positions]
    variants = [
        {"index": original[k], "chromosome": chrom, "position": int(pos), "ref": ref, "alt": alt}
        for k, (chrom, pos, ref, alt) in enumerate(ordered)
    ]
    return {
        "status": "completed",
        "operation": "order_variant_labels",
        "count": len(variants),
        "order": original,
        "variants": variants,
        "provenance": _provenance(
            _SELENE, sort="selene-sdk chromosome key, then position, ref, alt"
        ),
    }


def sort_standard_chromosomes(
    chromosomes: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sort unique 'chr' names with selene-sdk's key: chr1-22, chrX, chrY, chrM, then others."""

    if not _is_sequence(chromosomes) or not chromosomes:
        raise ConfigurationError("chromosomes must be a non-empty list of 'chr' names.")
    if len(chromosomes) > _MAX_CHROMOSOMES:
        raise ConfigurationError(f"chromosomes may contain at most {_MAX_CHROMOSOMES} names.")
    names = [
        _text(item, f"chromosomes[{index}]", maximum=128) for index, item in enumerate(chromosomes)
    ]
    for name in names:
        if not name.startswith("chr") or len(name) < 4:
            raise ConfigurationError(f"chromosome {name!r} must start with 'chr'.")
    if len(set(names)) != len(names):
        raise ConfigurationError("chromosomes must not repeat.")
    interpret = require("selene_sdk.interpret", _EXTRAS)
    ranks = {name: int(interpret.sort_standard_chrs(name)) for name in names}
    ordered = sorted(names, key=lambda name: (ranks[name], name))
    return {
        "status": "completed",
        "operation": "sort_standard_chromosomes",
        "count": len(ordered),
        "chromosomes": [{"chromosome": name, "rank": ranks[name]} for name in ordered],
        "provenance": _provenance(_SELENE, key="selene-sdk sort_standard_chrs, ties by name"),
    }


__all__ = [
    "apply_one_hot_mutations",
    "decode_one_hot_sequence",
    "interval_length_probabilities",
    "list_in_silico_mutants",
    "load_variant_diff_scores",
    "order_variant_labels",
    "read_selene_vcf",
    "rescale_attribution_scores",
    "reverse_complement_allele",
    "reverse_complement_one_hot",
    "sort_standard_chromosomes",
    "unpack_bit_packed_sequence",
]

_register_public_verbose(globals(), tuple(__all__))
