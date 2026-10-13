"""Offline pyfaidx helpers for BED lines, UCSC regions and chromosome-name maps.

The functions below call the upstream pyfaidx helpers on explicit strings or
name lists. They never open FASTA indexes, touch the network, or read files.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.sequence.backends import _provenance, _text

_PYFAIDX_EXTRA = "bio-pyfaidx"
_MAX_NAMES = 100_000
_MAX_NAME = 512
_MAX_LINE = 65_536


def _name_list(value: object, name: str, *, maximum: int) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of chromosome names.")
    if len(value) > maximum:
        raise ConfigurationError(f"{name} has more than {maximum} entries.")
    return [_text(item, name, maximum=_MAX_NAME) for item in value]


def translate_chromosome_names(
    names: Sequence[str],
    *,
    from_names: Sequence[str],
    to_names: Sequence[str],
    verbose: bool = True,
) -> dict[str, Any]:
    """Map chromosome names through explicit old-to-new pairs using pyfaidx.translate_chr_name."""

    source = _name_list(names, "names", maximum=_MAX_NAMES)
    old = _name_list(from_names, "from_names", maximum=_MAX_NAMES)
    new = _name_list(to_names, "to_names", maximum=_MAX_NAMES)
    if len(old) != len(new):
        raise ConfigurationError("from_names and to_names must have the same length.")
    if len(set(old)) != len(old):
        raise ConfigurationError("from_names must not repeat a name.")
    pyfaidx = require("pyfaidx", _PYFAIDX_EXTRA)
    mapper = pyfaidx.translate_chr_name(old, new)
    translations: list[dict[str, str]] = []
    unmapped: set[str] = set()
    for name in source:
        try:
            translated = mapper(name)
        except KeyError:
            unmapped.add(name)
            continue
        translations.append({"name": name, "translated": str(translated)})
    return {
        "status": "completed",
        "operation": "translate_chromosome_names",
        "requested_count": len(source),
        "translated_count": len(translations),
        "translations": translations,
        "unmapped": sorted(unmapped),
        "provenance": _provenance(
            "pyfaidx",
            method="pyfaidx.translate_chr_name lookup",
            unmapped="names without an explicit pair are reported, never guessed",
        ),
    }


def split_bed_entry(line: str, *, verbose: bool = True) -> dict[str, Any]:
    """Split one BED line into chromosome, 0-based start and end using pyfaidx.bed_split."""

    text = _text(line, "line", maximum=_MAX_LINE)
    pyfaidx = require("pyfaidx", _PYFAIDX_EXTRA)
    try:
        chromosome, start, end = pyfaidx.bed_split(text)
    except (pyfaidx.BedError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"pyfaidx could not parse the BED line: {exc}") from exc
    if chromosome is None:
        return {
            "status": "completed",
            "operation": "split_bed_entry",
            "is_comment": True,
            "chromosome": None,
            "start": None,
            "end": None,
            "length": None,
            "provenance": _provenance("pyfaidx", comment_rule="lines starting with '#'"),
        }
    if start < 0 or end < start:
        raise ConfigurationError(
            "BED start must be non-negative and end must not precede start; "
            f"got start={start}, end={end}."
        )
    return {
        "status": "completed",
        "operation": "split_bed_entry",
        "is_comment": False,
        "chromosome": str(chromosome),
        "start": start,
        "end": end,
        "length": end - start,
        "provenance": _provenance(
            "pyfaidx",
            coordinate_system="0-based start, 0-based exclusive end, as stored in BED",
        ),
    }


def split_ucsc_region(region: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a UCSC region such as chr1:101-200 into a 0-based half-open interval."""

    text = _text(region, "region", maximum=_MAX_NAME + 64)
    pyfaidx = require("pyfaidx", _PYFAIDX_EXTRA)
    chromosome, start, end = pyfaidx.ucsc_split(text)
    if start is None:
        # pyfaidx returns (name, None, None) for a bare name and for malformed intervals.
        if ":" in text:
            raise ConfigurationError(
                f"region {text!r} must be chromosome:start-end with integer coordinates; "
                "names that contain a colon need an explicit interval."
            )
        return {
            "status": "completed",
            "operation": "split_ucsc_region",
            "chromosome": str(chromosome),
            "start": None,
            "end": None,
            "length": None,
            "is_whole_chromosome": True,
            "provenance": _provenance("pyfaidx", coordinate_system="whole sequence"),
        }
    if start < 0:
        raise ConfigurationError("UCSC start coordinates are 1-based and must be at least 1.")
    if end <= start:
        raise ConfigurationError(f"UCSC region {text!r} ends before it starts.")
    return {
        "status": "completed",
        "operation": "split_ucsc_region",
        "chromosome": str(chromosome),
        "start": start,
        "end": end,
        "length": end - start,
        "is_whole_chromosome": False,
        "provenance": _provenance(
            "pyfaidx",
            coordinate_system="1-based inclusive input converted to 0-based half-open",
            ambiguity="names containing a colon are rejected unless they carry an interval",
        ),
    }


__all__ = [
    "split_bed_entry",
    "split_ucsc_region",
    "translate_chromosome_names",
]

_register_public_verbose(globals(), tuple(__all__))
