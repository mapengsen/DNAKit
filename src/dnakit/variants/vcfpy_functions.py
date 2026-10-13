"""Optional vcfpy adapters for VCF ALT alleles, field values and header lines.

- ``classify_vcf_alt()`` wraps ``vcfpy.parser.process_alt``: sequence, symbolic, breakend and
  single-breakend ALT alleles.
- ``classify_vcf_substitution()`` wraps ``vcfpy.parser.process_sub``: SNV, MNV, INS, DEL and INDEL
  from REF and ALT lengths and bases.
- ``parse_vcf_breakend()`` wraps ``vcfpy.parser.parse_breakend`` for bracket breakend notation.
- ``parse_vcf_field_value()`` and ``format_vcf_field_value()`` wrap vcfpy's INFO/FORMAT value
  parser and writer for a declared Number and Type.
- ``remove_vcf_header_lines()`` wraps ``vcfpy.header_without_lines`` on the header of a VCF file.

Records are not read. Each function handles one allele, one field value or one header.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``bio-vcfpy`` install hint.
"""

from __future__ import annotations

import importlib
import importlib.metadata
import math
import re
import warnings
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-vcfpy"
_VALUE_TYPES = ("Integer", "Float", "Flag", "Character", "String")
_SPECIAL_NUMBERS = ("A", "R", "G", ".")
_ALT_KINDS = {
    "Substitution": "substitution",
    "SymbolicAllele": "symbolic",
    "BreakEnd": "breakend",
    "SingleBreakEnd": "single_breakend",
}
_BASES = re.compile(r"^[ACGTNacgtn]+$")
_ALT_CHARS = re.compile(r"^[A-Za-z0-9<>\[\]:.*_-]+$")
_MAX_ALLELE = 100_000
_MAX_HEADER_REMOVALS = 10_000


def _vcfpy() -> Any:
    return require("vcfpy", _EXTRA)


def _submodule(name: str) -> Any:
    require("vcfpy", _EXTRA)
    return importlib.import_module(name)


def _provenance(convention: str) -> dict[str, str]:
    return {
        "package": "vcfpy",
        "version": importlib.metadata.version("vcfpy"),
        "convention": convention,
    }


def _call(label: str, function: Callable[..., Any], *args: Any) -> Any:
    vcfpy = _vcfpy()
    try:
        return function(*args)
    except (vcfpy.exceptions.VCFPyException, ValueError, TypeError, KeyError, IndexError) as exc:
        raise ConfigurationError(f"{label} failed: {exc}") from exc


def _bases(value: object, name: str) -> str:
    if not isinstance(value, str) or not _BASES.match(value):
        raise ConfigurationError(f"{name} must be a non-empty string of A, C, G, T or N.")
    if len(value) > _MAX_ALLELE:
        raise ConfigurationError(f"{name} exceeds {_MAX_ALLELE} bases.")
    return value.upper()


def _alt(value: object) -> str:
    if (
        not isinstance(value, str)
        or value == "."
        or not _ALT_CHARS.match(value)
        or len(value) > _MAX_ALLELE
    ):
        raise ConfigurationError(f"alt must be a VCF ALT allele; got {value!r}.")
    return value


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _number(value: object) -> int | str:
    if isinstance(value, bool):
        raise ConfigurationError("number must be a non-negative integer or A, R, G or '.'.")
    if isinstance(value, int):
        if value < 0:
            raise ConfigurationError("number must be a non-negative integer.")
        return value
    if value in _SPECIAL_NUMBERS:
        return str(value)
    raise ConfigurationError("number must be a non-negative integer or A, R, G or '.'.")


def _json(value: object) -> Any:
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, (list, tuple)):
        return [_json(item) for item in value]
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        return float(value) if math.isfinite(value) else None
    return str(value)


def _field_info(value_type: str, number: int | str, field_id: str) -> Any:
    header = _submodule("vcfpy.header")
    return header.FieldInfo(value_type, number, None, field_id or None)


def _atoms_match(value_type: str, items: Sequence[object]) -> None:
    for item in items:
        if item is None:
            continue
        if value_type in ("String", "Character"):
            ok = isinstance(item, str)
        elif value_type == "Integer":
            ok = isinstance(item, int) and not isinstance(item, bool)
        else:
            ok = isinstance(item, (int, float)) and not isinstance(item, bool)
        if not ok:
            raise ConfigurationError(f"value {item!r} does not match the {value_type} type.")


def classify_vcf_alt(ref: str, alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Classify one ALT allele against its REF bases with vcfpy.

    ``kind`` is ``substitution``, ``symbolic``, ``breakend`` or ``single_breakend``. Substitutions
    also report ``subtype`` (SNV, MNV, INS, DEL or INDEL) and the inserted or retained ``value``.
    """

    ref_text = _bases(ref, "ref")
    alt_text = _alt(alt)
    parser = _submodule("vcfpy.parser")
    header = _submodule("vcfpy.header")
    empty = header.Header(lines=[], samples=header.SamplesInfos([]))
    obj = _call("process_alt", parser.process_alt, empty, ref_text, alt_text)
    kind = _ALT_KINDS.get(type(obj).__name__)
    if kind is None:
        raise ConfigurationError(f"unrecognized ALT allele class {type(obj).__name__}.")
    breakend = None
    if kind in ("breakend", "single_breakend"):
        breakend = {key: _json(value) for key, value in sorted(vars(obj).items()) if key != "type"}
    return {
        "status": "completed",
        "operation": "classify_vcf_alt",
        "ref": ref_text,
        "alt": alt_text,
        "kind": kind,
        "subtype": _json(obj.type) if kind == "substitution" else None,
        "value": _json(getattr(obj, "value", None))
        if kind in ("substitution", "symbolic")
        else None,
        "breakend": breakend,
        "provenance": _provenance("VCF 4.x ALT classes; the upstream header argument is unused"),
    }


def parse_vcf_breakend(alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a bracket breakend ALT such as ``G]17:198982]`` into its mate and orientation.

    ``orientation`` and ``mate_orientation`` are ``+`` or ``-``. Angle-bracketed contigs such as
    ``<ctg1>:5`` set ``within_main_assembly`` to false.
    """

    text = _alt(alt)
    if "[" not in text and "]" not in text:
        raise ConfigurationError("a breakend ALT uses bracket notation such as 'G]17:198982]'.")
    parser = _submodule("vcfpy.parser")
    mate_chrom, mate_pos, orientation, mate_orientation, sequence, within = _call(
        "parse_breakend", parser.parse_breakend, text
    )
    return {
        "status": "completed",
        "operation": "parse_vcf_breakend",
        "alt": text,
        "mate_chrom": mate_chrom,
        "mate_pos": int(mate_pos),
        "orientation": _json(orientation),
        "mate_orientation": _json(mate_orientation),
        "sequence": sequence,
        "within_main_assembly": bool(within),
        "provenance": _provenance("VCF 4.x breakend notation; POS is 1-based"),
    }


def classify_vcf_substitution(ref: str, alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Classify a REF/ALT pair by allele lengths as SNV, MNV, INS, DEL or INDEL.

    The class depends only on lengths and on whether the first bases match, so identical REF and
    ALT of length one are reported as SNV.
    """

    ref_text = _bases(ref, "ref")
    alt_text = _bases(alt, "alt")
    parser = _submodule("vcfpy.parser")
    obj = _call("process_sub", parser.process_sub, ref_text, alt_text)
    return {
        "status": "completed",
        "operation": "classify_vcf_substitution",
        "ref": ref_text,
        "alt": alt_text,
        "subtype": _json(obj.type),
        "value": _json(obj.value),
        "provenance": _provenance("length-based REF/ALT class; bases are upper-cased"),
    }


def parse_vcf_field_value(
    text: str,
    *,
    value_type: str,
    number: int | str,
    field_id: str = "",
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse one INFO or FORMAT value according to its declared Number and Type.

    ``number`` is an integer or one of ``A``, ``R``, ``G`` or ``.``. A Number of 1 returns a
    scalar; other Numbers return a list. ``.`` is missing. Values that fail to convert raise an
    error instead of staying as strings.
    """

    selected = _choice(value_type, "value_type", _VALUE_TYPES)
    declared = _number(number)
    if not isinstance(text, str) or any(c in text for c in "\t\r\n\x00"):
        raise ConfigurationError("text must be a single VCF field value.")
    parser = _submodule("vcfpy.parser")
    exceptions = _submodule("vcfpy.exceptions")
    info = _field_info(selected, declared, field_id)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = _call("parse_field_value", parser.parse_field_value, info, text)
    if any(issubclass(item.category, exceptions.CannotConvertValue) for item in caught):
        raise ConfigurationError(f"{text!r} cannot be converted to {selected}.")
    return {
        "status": "completed",
        "operation": "parse_vcf_field_value",
        "value_type": selected,
        "number": declared,
        "value": _json(value),
        "provenance": _provenance("VCF field value conversion by declared Number and Type"),
    }


def format_vcf_field_value(
    value: object,
    *,
    value_type: str,
    number: int | str,
    section: str = "INFO",
    field_id: str = "",
    verbose: bool = True,
) -> dict[str, Any]:
    """Format one INFO or FORMAT value as VCF text, the inverse of ``parse_vcf_field_value``.

    A Number of 1 takes one value (or None for missing). Other Numbers take a list. Flag values
    are not formatted here; a Flag INFO key is written as a bare key.
    """

    selected = _choice(value_type, "value_type", _VALUE_TYPES)
    if selected == "Flag":
        raise ConfigurationError("Flag values are written as bare keys and are not formatted here.")
    declared = _number(number)
    layout = _choice(section, "section", ("INFO", "FORMAT"))
    writer = _submodule("vcfpy.writer")
    if declared == 1:
        if isinstance(value, (list, tuple)):
            raise ConfigurationError("a Number of 1 takes one value, not a list.")
        items: list[object] = [] if value is None else [value]
    else:
        if value is None:
            items = []
        elif isinstance(value, (list, tuple)):
            items = list(value)
        else:
            raise ConfigurationError("this Number takes a list of values.")
    _atoms_match(selected, items)
    info = _field_info(selected, declared, field_id)
    payload: Any = value if declared == 1 else (items or None)
    text = _call("format_value", writer.format_value, info, payload, layout)
    return {
        "status": "completed",
        "operation": "format_vcf_field_value",
        "value_type": selected,
        "number": declared,
        "section": layout,
        "text": text,
        "provenance": _provenance("VCF field text; floats keep the Python repr"),
    }


def remove_vcf_header_lines(
    source: str | Path,
    remove: Sequence[Sequence[str | None]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the header of a VCF file without the lines named by ``(key, ID)`` pairs.

    For lines that have no ID mapping, give the whole value in place of the ID.
    """

    path = input_path(str(source))
    if isinstance(remove, (str, bytes)) or not isinstance(remove, Sequence):
        raise ConfigurationError("remove must be a list of (key, ID) pairs.")
    if len(remove) > _MAX_HEADER_REMOVALS:
        raise ConfigurationError(f"remove exceeds the limit of {_MAX_HEADER_REMOVALS} pairs.")
    pairs: list[tuple[str, str | None]] = []
    for item in remove:
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)) or len(item) != 2:
            raise ConfigurationError("each removal must be a (key, ID) pair.")
        key, identifier = item
        if not isinstance(key, str) or not key:
            raise ConfigurationError("each removal key must be a non-empty string.")
        if identifier is not None and not isinstance(identifier, str):
            raise ConfigurationError("each removal ID must be a string or None.")
        pairs.append((key, identifier))
    vcfpy = _vcfpy()
    reader = _call("Reader.from_path", vcfpy.Reader.from_path, str(path))
    try:
        header = reader.header
    finally:
        reader.close()
    trimmed = _call("header_without_lines", vcfpy.header_without_lines, header, pairs)
    return {
        "status": "completed",
        "operation": "remove_vcf_header_lines",
        "source": str(path),
        "header_lines": [line.serialize() for line in trimmed.lines],
        "samples": list(trimmed.samples.names),
        "removed_lines": len(header.lines) - len(trimmed.lines),
        "provenance": _provenance(
            "header lines matched by (key, ID); the #CHROM line is not included"
        ),
    }


__all__ = [
    "classify_vcf_alt",
    "classify_vcf_substitution",
    "format_vcf_field_value",
    "parse_vcf_breakend",
    "parse_vcf_field_value",
    "remove_vcf_header_lines",
]

_register_public_verbose(globals(), tuple(__all__))
