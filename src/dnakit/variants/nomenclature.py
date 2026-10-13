"""Optional variant-nomenclature adapters for varcode and hgvs.

The functions in this module are explicit adapters over optional packages:

- ``varcode`` builds variant objects. ``describe_variant()`` reports the minimal (trimmed)
  representation and the kind of a variant from its alleles alone. ``load_vcf_variants()``
  loads VCF records through ``varcode.load_vcf``. Effect prediction needs Ensembl data that
  pyensembl downloads, so it is not wrapped here.
- ``hgvs`` parses, formats and syntax-checks HGVS strings (c., g., n., m., r., p.).
  Only the parser and the formatter are used. The reference-based parts of hgvs (sequence
  fetching through seqrepo, projection through UTA, validation and normalization) are not
  wrapped, because they need reference data that is not available offline.

Coordinates are 1-based. None of these functions reads a reference sequence.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import re
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_VARCODE_EXTRA = "bio-varcode"
_HGVS_EXTRA = "bio-hgvs"
_MAX_HGVS_LENGTH = 1_000
_MAX_HGVS_BATCH = 10_000
_HGVS_PARSER: Any = None
_ALLELE_PATTERN = re.compile(r"[ACGTN]*")
_MAX_ALLELE_LENGTH = 10_000
_MAX_CONTIG_LENGTH = 64
_MAX_POSITION = 2_000_000_000
_DEFAULT_MAX_VARIANTS = 100_000
_MAX_VARIANTS_LIMIT = 1_000_000
_DEFAULT_GENOME = "GRCh38"


def _contig(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > _MAX_CONTIG_LENGTH
        or any(character.isspace() or character == "\x00" for character in value)
    ):
        raise ConfigurationError("contig must be a non-empty name without whitespace.")
    return value


def _start(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_POSITION:
        raise ConfigurationError(f"start must be a 1-based integer in [1, {_MAX_POSITION}].")
    return value


def _allele(value: object, name: str) -> str:
    if not isinstance(value, str) or len(value) > _MAX_ALLELE_LENGTH:
        raise ConfigurationError(f"{name} must be a string of at most {_MAX_ALLELE_LENGTH} bases.")
    text = value.upper()
    if not _ALLELE_PATTERN.fullmatch(text):
        raise ConfigurationError(f"{name} may contain only A, C, G, T and N.")
    return text


def _kind(variant: Any) -> str:
    if variant.is_snv:
        return "snv"
    if variant.is_insertion:
        return "insertion"
    if variant.is_deletion:
        return "deletion"
    if len(variant.trimmed_ref) == len(variant.trimmed_alt):
        return "substitution"
    return "complex"


def _describe(variant: Any) -> dict[str, Any]:
    # varcode trims shared prefixes in its constructor; the original_* fields keep the input form.
    return {
        "contig": str(variant.contig),
        "start": int(variant.original_start),
        "ref": str(variant.original_ref),
        "alt": str(variant.original_alt),
        "trimmed_start": int(variant.trimmed_base1_start),
        "trimmed_end": int(variant.trimmed_base1_end),
        "trimmed_ref": str(variant.trimmed_ref),
        "trimmed_alt": str(variant.trimmed_alt),
        "kind": _kind(variant),
        "is_transition": bool(variant.is_transition) if variant.is_snv else None,
    }


def _provenance(**notes: str) -> dict[str, str]:
    return {"package": f"varcode {metadata.version('varcode')}", **notes}


def describe_variant(
    contig: str, start: int, ref: str, alt: str, *, verbose: bool = True
) -> dict[str, Any]:
    """Describe one variant by its trimmed alleles and kind, without a reference sequence."""

    name = _contig(contig)
    position = _start(start)
    first = _allele(ref, "ref")
    second = _allele(alt, "alt")
    if not first and not second:
        raise ConfigurationError("ref and alt cannot both be empty.")
    if first == second:
        raise ConfigurationError("ref and alt must differ.")
    module = require("varcode", _VARCODE_EXTRA)
    variant = module.Variant(name, position, first, second)
    return {
        "status": "completed",
        "operation": "describe_variant",
        **_describe(variant),
        "provenance": _provenance(
            reference="none: alleles are trimmed without reading a reference sequence"
        ),
    }


def load_vcf_variants(
    path: object,
    *,
    genome: str = _DEFAULT_GENOME,
    only_passing: bool = True,
    max_variants: int = _DEFAULT_MAX_VARIANTS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load VCF records as varcode variants, with trimmed alleles and kinds for each record."""

    if isinstance(max_variants, bool) or not isinstance(max_variants, int):
        raise ConfigurationError(f"max_variants must be an integer in [1, {_MAX_VARIANTS_LIMIT}].")
    if not 1 <= max_variants <= _MAX_VARIANTS_LIMIT:
        raise ConfigurationError(f"max_variants must be an integer in [1, {_MAX_VARIANTS_LIMIT}].")
    if not isinstance(only_passing, bool):
        raise ConfigurationError("only_passing must be True or False.")
    if not isinstance(genome, str) or not genome or any(c.isspace() for c in genome):
        raise ConfigurationError("genome must be a non-empty reference name such as 'GRCh38'.")
    source: Path = input_path(path)
    module = require("varcode", _VARCODE_EXTRA)
    # varcode truncates silently at max_variants, so request one extra record to detect overflow.
    collection = module.load_vcf(
        str(source),
        genome=genome,
        only_passing=only_passing,
        max_variants=max_variants + 1,
        distinct=False,
        include_info=False,
    )
    if len(collection) > max_variants:
        raise ConfigurationError(f"VCF has more than max_variants={max_variants} records.")
    return {
        "status": "completed",
        "operation": "load_vcf_variants",
        "count": len(collection),
        "variants": [_describe(variant) for variant in collection],
        "provenance": _provenance(
            reference=f"genome name {genome} is used only to infer contig names",
            filter="only_passing keeps PASS and missing FILTER records" if only_passing else "none",
        ),
    }


def _hgvs_text(value: object, name: str = "text") -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > _MAX_HGVS_LENGTH
        or any(character.isspace() or character == "\x00" for character in value)
    ):
        raise ConfigurationError(
            f"{name} must be a non-empty HGVS string of at most {_MAX_HGVS_LENGTH} characters "
            "without whitespace."
        )
    return value


def _hgvs_parser() -> Any:
    """Return the hgvs parser module and one shared Parser instance."""

    global _HGVS_PARSER
    module = require("hgvs.parser", _HGVS_EXTRA)
    if _HGVS_PARSER is None:
        _HGVS_PARSER = module.Parser()
    return module


def _position(item: Any) -> dict[str, int | None] | None:
    if item is None:
        return None
    base = getattr(item, "base", None)
    offset = getattr(item, "offset", None)
    return {
        "base": int(base) if base is not None else None,
        "offset": int(offset) if offset is not None else None,
    }


def _parse_checked(text: str) -> Any:
    """Parse an HGVS string offline and reject reversed intervals."""

    module = _hgvs_parser()
    try:
        variant = _HGVS_PARSER.parse_hgvs_variant(text)
    except module.HGVSParseError as exc:
        raise ConfigurationError(f"Invalid HGVS string {text!r}: {exc}") from exc
    position = getattr(variant.posedit, "pos", None)
    start, end = getattr(position, "start", None), getattr(position, "end", None)
    if start is not None and end is not None:
        first = _position(start) or {}
        last = _position(end) or {}
        if (
            first.get("base") is not None
            and last.get("base") is not None
            and (first["base"], first["offset"] or 0) > (last["base"], last["offset"] or 0)
        ):
            raise ConfigurationError(
                f"Invalid HGVS string {text!r}: interval start is after its end."
            )
    return variant


def _hgvs_record(variant: Any) -> dict[str, Any]:
    posedit = variant.posedit
    edit = getattr(posedit, "edit", None)
    position = getattr(posedit, "pos", None)
    start = getattr(position, "start", None)
    end = getattr(position, "end", None)
    return {
        "reference": variant.ac,
        "gene": variant.gene,
        "coordinate_type": variant.type,
        "start": _position(start),
        "end": _position(end),
        "edit_type": getattr(edit, "type", None),
        "reference_bases": getattr(edit, "ref", None) or None,
        "alternate_bases": getattr(edit, "alt", None) or None,
        "formatted": variant.format(),
    }


def parse_hgvs_variant(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse one HGVS string (c., g., n., m., r. or p.) without reference data."""

    value = _hgvs_text(text)
    variant = _parse_checked(value)
    return {
        "status": "completed",
        "operation": "parse_hgvs_variant",
        "input": value,
        **_hgvs_record(variant),
        "provenance": _hgvs_provenance(),
    }


def format_hgvs_variant(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Format one HGVS string into hgvs's canonical text without reference data."""

    value = _hgvs_text(text)
    formatted = _parse_checked(value).format()
    return {
        "status": "completed",
        "operation": "format_hgvs_variant",
        "input": value,
        "formatted": formatted,
        "changed": formatted != value,
        "provenance": _hgvs_provenance(),
    }


def check_hgvs_syntax(
    texts: list[str], *, max_variants: int = _MAX_HGVS_BATCH, verbose: bool = True
) -> dict[str, Any]:
    """Check a batch of HGVS strings and report valid and invalid entries one by one."""

    if isinstance(max_variants, bool) or not isinstance(max_variants, int):
        raise ConfigurationError(f"max_variants must be an integer in [1, {_MAX_HGVS_BATCH}].")
    if not 1 <= max_variants <= _MAX_HGVS_BATCH:
        raise ConfigurationError(f"max_variants must be an integer in [1, {_MAX_HGVS_BATCH}].")
    if not isinstance(texts, list) or not texts:
        raise ConfigurationError("texts must be a non-empty list of HGVS strings.")
    if len(texts) > max_variants:
        raise ConfigurationError(
            f"texts has {len(texts)} entries, above max_variants={max_variants}."
        )
    _hgvs_parser()
    items: list[dict[str, Any]] = []
    for text in texts:
        if not isinstance(text, str):
            raise ConfigurationError("every entry in texts must be a string.")
        try:
            value = _hgvs_text(text)
            record = _hgvs_record(_parse_checked(value))
            items.append(
                {"input": value, "valid": True, "formatted": record["formatted"], "error": None}
            )
        except ConfigurationError as exc:
            items.append({"input": text, "valid": False, "formatted": None, "error": str(exc)})
    valid = sum(1 for item in items if item["valid"])
    return {
        "status": "completed",
        "operation": "check_hgvs_syntax",
        "count": len(items),
        "valid": valid,
        "invalid": len(items) - valid,
        "items": items,
        "provenance": _hgvs_provenance(),
    }


def _hgvs_provenance() -> dict[str, str]:
    return {
        "package": f"hgvs {metadata.version('hgvs')}",
        "scope": "parsing and formatting only; reference data (seqrepo, UTA) is not used",
    }


__all__ = [
    "check_hgvs_syntax",
    "describe_variant",
    "format_hgvs_variant",
    "load_vcf_variants",
    "parse_hgvs_variant",
]

_register_public_verbose(globals(), tuple(__all__))
