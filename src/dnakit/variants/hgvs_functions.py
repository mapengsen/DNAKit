"""Offline hgvs wrappers: indel normalization, CIGAR building, repeat blocks and accession checks.

These functions take the reference sequence or alignment blocks as arguments, so they run
without seqrepo, UTA or PostgreSQL. Projection and sequence-fetching functions of hgvs
(``full_house``, ``variant_context_w_alignment``) need those services and are not wrapped.

Coordinates are 0-based and half-open for sequence intervals, matching hgvs's internal
normalization routines. Packages are imported only when a function runs; a missing hgvs
raises ``BackendUnavailableError`` with the ``pip install "dnakit[bio-hgvs]"`` hint.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from importlib import metadata
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_HGVS_EXTRA = "bio-hgvs"
_MAX_SEQUENCE = 100_000
_MAX_BLOCKS = 10_000
_BASES = re.compile(r"[ACGTN]*")
_CIGAR = re.compile(r"(\d+[=DIMNX])+")


def _provenance() -> dict[str, str]:
    return {
        "package": f"hgvs {metadata.version('hgvs')}",
        "scope": "sequence and alignment arithmetic only; seqrepo and UTA are not used",
    }


def _module(name: str) -> Any:
    """Import hgvs.parser first: importing hgvs submodules directly hits a circular import."""

    require("hgvs.parser", _HGVS_EXTRA)
    return require(name, _HGVS_EXTRA)


def _bases(value: object, name: str, *, allow_empty: bool = True) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string of A, C, G, T or N.")
    text = value.upper()
    if (not text and not allow_empty) or len(text) > _MAX_SEQUENCE or not _BASES.fullmatch(text):
        raise ConfigurationError(f"{name} must contain only A, C, G, T or N.")
    return text


def _int(value: object, name: str, *, minimum: int = 0, maximum: int = _MAX_SEQUENCE) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _mapping(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping.")
    return dict(value)


def normalize_indel_alleles(
    reference: str,
    start: int,
    stop: int,
    alleles: Sequence[str],
    *,
    left: bool = True,
    bound: int | None = None,
    ref_step: int = 1,
    shuffle: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shift an indel to its leftmost (or rightmost) equivalent position in a reference string."""

    sequence = _bases(reference, "reference", allow_empty=False)
    left_edge = _int(start, "start", maximum=len(sequence))
    right_edge = _int(stop, "stop", maximum=len(sequence))
    if right_edge < left_edge:
        raise ConfigurationError("stop must be at least start.")
    if not isinstance(alleles, Sequence) or isinstance(alleles, str) or len(alleles) < 2:
        raise ConfigurationError("alleles must be a sequence of at least two allele strings.")
    cleaned = tuple(_bases(allele, "alleles item") for allele in alleles)
    if not isinstance(left, bool) or not isinstance(shuffle, bool):
        raise ConfigurationError("left and shuffle must be True or False.")
    step = _int(ref_step, "ref_step", minimum=1, maximum=len(sequence) or 1)
    if bound is None:
        limit = 0 if left else len(sequence)
    else:
        limit = _int(bound, "bound", maximum=len(sequence))
    module = _module("hgvs.utils.norm")
    normalized = module.normalize_alleles(
        sequence, left_edge, right_edge, cleaned, limit, step, left, shuffle=shuffle
    )
    return {
        "status": "completed",
        "operation": "normalize_indel_alleles",
        "direction": "left" if left else "right",
        "input": {"start": left_edge, "stop": right_edge, "alleles": list(cleaned)},
        "start": int(normalized.start),
        "stop": int(normalized.stop),
        "alleles": list(normalized.alleles),
        "changed": (int(normalized.start), int(normalized.stop), tuple(normalized.alleles))
        != (left_edge, right_edge, cleaned),
        "provenance": _provenance(),
    }


def trim_common_prefix_strings(
    strings: Sequence[str], *, min_length: int = 0, verbose: bool = True
) -> dict[str, Any]:
    """Remove the prefix shared by a list of allele strings and report the bases removed."""

    if not isinstance(strings, Sequence) or isinstance(strings, str) or len(strings) < 2:
        raise ConfigurationError("strings must be a sequence of at least two allele strings.")
    cleaned = [_bases(item, "strings item") for item in strings]
    floor = _int(min_length, "min_length", maximum=_MAX_SEQUENCE)
    module = _module("hgvs.utils.norm")
    trimmed, remaining = module.trim_common_prefixes(cleaned, min_len=floor)
    return {
        "status": "completed",
        "operation": "trim_common_prefix_strings",
        "input": cleaned,
        "trimmed_length": int(trimmed),
        "strings": list(remaining),
        "provenance": _provenance(),
    }


def build_transcript_cigar(
    exons: Sequence[Mapping[str, Any]],
    strand: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join per-exon CIGAR strings into one transcript CIGAR, with N for each intron."""

    if strand not in (1, -1) or isinstance(strand, bool):
        raise ConfigurationError("strand must be 1 or -1.")
    if not isinstance(exons, Sequence) or isinstance(exons, str):
        raise ConfigurationError("exons must be a list of exon mappings in transcript order.")
    records = []
    for index, exon in enumerate(exons):
        item = _mapping(exon, f"exons[{index}]")
        missing = sorted({"cigar", "alt_start_i", "alt_end_i"} - item.keys())
        if missing:
            raise ConfigurationError(f"exons[{index}] is missing fields: {', '.join(missing)}.")
        cigar = item["cigar"]
        if not isinstance(cigar, str) or not _CIGAR.fullmatch(cigar):
            raise ConfigurationError(f"exons[{index}].cigar must look like '10=' or '3M2D5='.")
        start = _int(item["alt_start_i"], f"exons[{index}].alt_start_i", maximum=2_000_000_000)
        end = _int(item["alt_end_i"], f"exons[{index}].alt_end_i", maximum=2_000_000_000)
        records.append({"cigar": cigar, "alt_start_i": start, "alt_end_i": end})
    module = _module("hgvs.utils")
    cigar = module.build_tx_cigar(records, strand) if records else None
    return {
        "status": "completed",
        "operation": "build_transcript_cigar",
        "strand": strand,
        "exon_count": len(records),
        "cigar": cigar,
        "provenance": _provenance(),
    }


_LEVELS = {1: "VALID", 2: "WARNING", 3: "ERROR"}


def validate_hgvs_type_accession(
    variant_type: str, accession: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Check that an accession such as NM_000088.3 fits an HGVS type (c, g, m, n, p or r)."""

    module = _module("hgvs.utils.validation")
    if not isinstance(variant_type, str) or variant_type not in module.valid_pairs:
        raise ConfigurationError("variant_type must be one of: c, g, m, n, p, r.")
    if not isinstance(accession, str) or not accession or any(c.isspace() for c in accession):
        raise ConfigurationError("accession must be a non-empty string without whitespace.")
    level, message = module.validate_type_ac_pair(variant_type, accession)
    return {
        "status": "completed",
        "operation": "validate_hgvs_type_accession",
        "variant_type": variant_type,
        "accession": accession,
        "level": _LEVELS[int(level.value)],
        "message": message,
        "provenance": _provenance(),
    }


def _block_row(unit: Any) -> dict[str, Any]:
    return {
        "repeat_unit": str(unit.repeat_unit),
        "repeat_count": int(unit.repeat_count),
        "block_size": int(unit.block_size),
        "block": str(unit.block),
    }


def detect_repeat_blocks(
    sequence: str, *, reverse: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Find tandem-repeat blocks in a sequence, scanning left to right or right to left."""

    text = _bases(sequence, "sequence", allow_empty=False)
    if not isinstance(reverse, bool):
        raise ConfigurationError("reverse must be True or False.")
    module = _module("hgvs.repeats")
    units = module.detect_repetitive_block_lengths(text, reverse=reverse)
    return {
        "status": "completed",
        "operation": "detect_repeat_blocks",
        "sequence": text,
        "reverse": reverse,
        "count": len(units),
        "blocks": [_block_row(unit) for unit in units],
        "provenance": _provenance(),
    }


def assemble_repeat_blocks(
    sequence: str,
    blocks: Sequence[Mapping[str, Any]],
    *,
    reverse: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format detected repeat blocks as HGVS-style repeat notation, such as AC[3]G[1]T[4]."""

    text = _bases(sequence, "sequence", allow_empty=False)
    if not isinstance(blocks, Sequence) or isinstance(blocks, str) or not blocks:
        raise ConfigurationError("blocks must be a non-empty list of repeat-block mappings.")
    if len(blocks) > _MAX_BLOCKS:
        raise ConfigurationError(f"blocks may hold at most {_MAX_BLOCKS} entries.")
    module = _module("hgvs.repeats")
    units = []
    for index, block in enumerate(blocks):
        item = _mapping(block, f"blocks[{index}]")
        missing = sorted({"repeat_unit", "repeat_count", "block_size", "block"} - item.keys())
        if missing:
            raise ConfigurationError(f"blocks[{index}] is missing fields: {', '.join(missing)}.")
        units.append(
            module.RepeatUnit(
                repeat_count=_int(item["repeat_count"], f"blocks[{index}].repeat_count", minimum=1),
                repeat_unit=_bases(
                    item["repeat_unit"], f"blocks[{index}].repeat_unit", allow_empty=False
                ),
                block_size=_int(item["block_size"], f"blocks[{index}].block_size", minimum=1),
                block=_bases(item["block"], f"blocks[{index}].block", allow_empty=False),
            )
        )
    if not isinstance(reverse, bool):
        raise ConfigurationError("reverse must be True or False.")
    return {
        "status": "completed",
        "operation": "assemble_repeat_blocks",
        "sequence": text,
        "reverse": reverse,
        "formatted": module.assemble_repeat_string(text, units, reverse=reverse),
        "provenance": _provenance(),
    }


__all__ = [
    "assemble_repeat_blocks",
    "build_transcript_cigar",
    "detect_repeat_blocks",
    "normalize_indel_alleles",
    "trim_common_prefix_strings",
    "validate_hgvs_type_accession",
]

_register_public_verbose(globals(), tuple(__all__))
