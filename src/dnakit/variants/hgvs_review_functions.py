"""Offline hgvs helpers for allele suffix trimming, interbase spans and repeat notation.

- ``trim_common_suffix_strings()`` removes the suffix shared by allele strings.
- ``hgvs_interbase_span()`` parses one HGVS string without reference data and returns its
  0-based, half-open interval.
- ``repeat_notation_for_alleles()`` formats repeat blocks from two sequences as HGVS-style
  repeat notation, using hgvs' own repeat assembly.

Nothing here needs seqrepo, UTA or PostgreSQL. ``hgvs`` is imported only when a function runs;
a missing package raises ``BackendUnavailableError`` with the ``bio-hgvs`` install hint.
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
_MAX_ALLELES = 1_000
_MAX_BLOCKS = 10_000
_BASES = re.compile(r"[ACGTN]*")
_VARIANT = re.compile(r"[A-Za-z0-9_.()-]+:[gcnmrp]\..+")


def _provenance() -> dict[str, str]:
    return {
        "package": f"hgvs {metadata.version('hgvs')}",
        "scope": "string and parser arithmetic only; seqrepo, UTA and PostgreSQL are not used",
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


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _alleles(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a sequence of allele strings.")
    if len(values) < 2 or len(values) > _MAX_ALLELES:
        raise ConfigurationError(f"{name} must hold between 2 and {_MAX_ALLELES} alleles.")
    return [_bases(item, f"{name} item", allow_empty=False) for item in values]


def _blocks(values: object, name: str) -> list[dict[str, Any]]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a sequence of repeat-block mappings.")
    if len(values) > _MAX_BLOCKS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_BLOCKS} entries.")
    blocks: list[dict[str, Any]] = []
    for index, block in enumerate(values):
        if not isinstance(block, Mapping):
            raise ConfigurationError(f"{name}[{index}] must be a mapping.")
        missing = sorted({"repeat_unit", "repeat_count", "block_size", "block"} - block.keys())
        if missing:
            raise ConfigurationError(f"{name}[{index}] is missing fields: {', '.join(missing)}.")
        blocks.append(dict(block))
    return blocks


def _repeat_units(module: Any, blocks: list[dict[str, Any]], name: str) -> list[Any]:
    units = []
    for index, item in enumerate(blocks):
        units.append(
            module.RepeatUnit(
                repeat_count=_int(item["repeat_count"], f"{name}[{index}].repeat_count", minimum=1),
                repeat_unit=_bases(
                    item["repeat_unit"], f"{name}[{index}].repeat_unit", allow_empty=False
                ),
                block_size=_int(item["block_size"], f"{name}[{index}].block_size", minimum=1),
                block=_bases(item["block"], f"{name}[{index}].block", allow_empty=False),
            )
        )
    return units


def trim_common_suffix_strings(
    strings: Sequence[str], *, min_length: int = 0, verbose: bool = True
) -> dict[str, Any]:
    """Remove the suffix shared by a list of allele strings and report the bases removed."""

    cleaned = _alleles(strings, "strings")
    floor = _int(min_length, "min_length", maximum=_MAX_SEQUENCE)
    module = _module("hgvs.utils.norm")
    trimmed, remaining = module.trim_common_suffixes(cleaned, min_len=floor)
    return {
        "status": "completed",
        "operation": "trim_common_suffix_strings",
        "input": cleaned,
        "trimmed_length": int(trimmed),
        "strings": list(remaining),
        "provenance": _provenance(),
    }


def hgvs_interbase_span(
    variant: str, *, outer_confidence: bool = True, verbose: bool = True
) -> dict[str, Any]:
    """Parse one HGVS variant and return its 0-based, half-open (interbase) start and end."""

    if not isinstance(variant, str) or not _VARIANT.fullmatch(variant.strip()):
        raise ConfigurationError(
            "variant must look like NC_000001.11:g.100_102del (accession:type.location)."
        )
    _flag(outer_confidence, "outer_confidence")
    text = variant.strip()
    parser_module = _module("hgvs.parser")
    utility = _module("hgvs.utils.position")
    try:
        parsed = parser_module.Parser().parse_hgvs_variant(text)
    except Exception as exc:  # hgvs raises its own parse-error hierarchy
        raise ConfigurationError(f"variant is not valid HGVS: {exc}") from exc
    start, end = utility.get_start_end_interbase(
        parsed.posedit.pos, outer_confidence=outer_confidence
    )
    return {
        "status": "completed",
        "operation": "hgvs_interbase_span",
        "input": text,
        "start": int(start),
        "end": int(end),
        "coordinates": "0-based, half-open (interbase)",
        "provenance": _provenance(),
    }


def repeat_notation_for_alleles(
    sequence: str,
    other_sequence: str,
    blocks: Sequence[Mapping[str, Any]],
    other_blocks: Sequence[Mapping[str, Any]] = (),
    *,
    reverse: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format repeat notation for one allele, such as AC[3], from its repeat blocks.

    ``blocks`` and ``other_blocks`` are the repeat-block rows returned by
    ``detect_repeat_blocks()``. ``formatted`` is ``None`` when hgvs finds no notation.
    """

    primary_text = _bases(sequence, "sequence")
    other_text = _bases(other_sequence, "other_sequence")
    primary = _blocks(blocks, "blocks")
    other = _blocks(other_blocks, "other_blocks")
    _flag(reverse, "reverse")
    module = _module("hgvs.repeats")
    formatted = module.get_repeat_str(
        primary_text,
        other_text,
        _repeat_units(module, primary, "blocks"),
        _repeat_units(module, other, "other_blocks"),
        reverse=reverse,
    )
    return {
        "status": "completed",
        "operation": "repeat_notation_for_alleles",
        "sequence": primary_text,
        "other_sequence": other_text,
        "reverse": reverse,
        "formatted": formatted,
        "provenance": _provenance(),
    }


__all__ = [
    "hgvs_interbase_span",
    "repeat_notation_for_alleles",
    "trim_common_suffix_strings",
]

_register_public_verbose(globals(), tuple(__all__))
