"""SBOL3 design-identity and copy operations backed by pySBOL3 (``sbol3``).

pySBOL3 is imported only inside each function. Every function returns a
JSON-safe dictionary with ``status``, ``operation``, its result fields and
``provenance``.
"""

from __future__ import annotations

import importlib.metadata
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MAX_FILE_BYTES = 50_000_000
_MAX_TOP_LEVELS = 10_000
_MAX_TEXT = 256
_MAX_IDENTITY = 2048
_DISPLAY_ID_RULE = "alphanumeric or underscore characters only, and not beginning with a digit"
_OUTPUT_FORMATS = {
    ".ttl": "TURTLE",
    ".nt": "NTRIPLES",
    ".rdf": "RDF_XML",
    ".xml": "RDF_XML",
    ".jsonld": "JSONLD",
}


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version("sbol3")
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": "sbol3", "version": version, **notes}


def _text(value: object, name: str, *, maximum: int = _MAX_TEXT) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum or "\x00" in value:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    return value


def _source(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("source must be a local SBOL file path.")
    path = input_path(value)
    if path.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"Input file exceeds {_MAX_FILE_BYTES} bytes: {path.name}")
    return path


def _namespace(value: object) -> str:
    text = _text(value, "into_namespace").strip().rstrip("/")
    parsed = urlparse(text)
    if not parsed.scheme or not parsed.netloc or any(ch.isspace() for ch in text):
        raise ConfigurationError(
            "into_namespace must be an absolute URI such as https://example.org/lab."
        )
    return text


def _output(value: object) -> tuple[Path, str]:
    text = _text(str(value) if isinstance(value, Path) else value, "output", maximum=4096)
    path = Path(text).expanduser().resolve()
    format_name = _OUTPUT_FORMATS.get(path.suffix.lower())
    if format_name is None:
        raise ConfigurationError("output must end with .ttl, .nt, .rdf, .xml or .jsonld.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"output already exists; choose a new file: {path}")
    return path, format_name


def _read_document(sbol3: Any, path: Path) -> Any:
    document = sbol3.Document()
    try:
        document.read(str(path))
    except (OSError, SyntaxError, TypeError, ValueError, RuntimeError, sbol3.SBOLError) as exc:
        raise ConfigurationError(f"pySBOL3 could not read {path.name}: {exc}") from exc
    return document


def copy_sbol3_design(
    source: str | Path,
    *,
    into_namespace: str,
    output: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Copy every top-level SBOL3 object of a local file into a new namespace.

    The clones keep their display IDs and paths. References between the copied
    objects are remapped to the new namespace, regardless of the order in which
    the objects were read. ``output`` optionally writes the copies to a new file;
    its suffix selects the RDF serialization.
    """

    sbol3 = require("sbol3", "bio-pysbol3")
    path = _source(source)
    namespace = _namespace(into_namespace)
    target = _output(output) if output is not None else None
    document = _read_document(sbol3, path)
    originals = sorted(
        (item for item in document.objects if isinstance(item, sbol3.TopLevel)),
        key=lambda item: str(item.identity),
    )
    if not originals:
        raise ConfigurationError(f"{path.name} contains no top-level SBOL objects.")
    if len(originals) > _MAX_TOP_LEVELS:
        raise ConfigurationError(f"SBOL file has more than {_MAX_TOP_LEVELS} top-level objects.")
    source_identities = [str(item.identity) for item in originals]
    clones = sbol3.copy(originals)
    # sbol3.copy(into_namespace=...) can leave a reference to a clone that was
    # renamed later in the same pass. Repeating the reference update over the
    # whole group fixes this without depending on the input order.
    sbol3.Document.change_object_namespace(clones, namespace, update_references=clones)
    copied = sbol3.Document()
    copied.add(clones)
    output_format = None
    if target is not None:
        target_path, output_format = target
        copied.write(str(target_path), getattr(sbol3, output_format))
    rows = [
        {
            "type": type(clone).__name__,
            "source_identity": source_identity,
            "identity": str(clone.identity),
            "display_id": clone.display_id,
        }
        for source_identity, clone in zip(source_identities, clones, strict=True)
    ]
    rows.sort(key=lambda row: row["identity"])
    return {
        "status": "completed",
        "operation": "copy_sbol3_design",
        "source": str(path),
        "into_namespace": namespace,
        "copied_count": len(rows),
        "copied": rows,
        "output": str(target[0]) if target is not None else None,
        "output_format": output_format,
        "provenance": _provenance(
            standard="SBOL3", reference_repair="second update_references pass after sbol3.copy"
        ),
    }


def to_sbol3_display_id(name: str, *, verbose: bool = True) -> dict[str, Any]:
    """Convert a free-text name into a valid SBOL3 display ID with pySBOL3."""

    sbol3 = require("sbol3", "bio-pysbol3")
    text = _text(name, "name")
    display_id = str(sbol3.string_to_display_id(text))
    return {
        "status": "completed",
        "operation": "to_sbol3_display_id",
        "name": text,
        "display_id": display_id,
        "valid": bool(sbol3.is_valid_display_id(display_id)),
        "rule": _DISPLAY_ID_RULE,
        "provenance": _provenance(),
    }


def check_sbol3_display_id(display_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Report whether a string is a valid SBOL3 display ID."""

    sbol3 = require("sbol3", "bio-pysbol3")
    if not isinstance(display_id, str) or len(display_id) > _MAX_TEXT or "\x00" in display_id:
        raise ConfigurationError(f"display_id must be a string of at most {_MAX_TEXT} characters.")
    return {
        "status": "completed",
        "operation": "check_sbol3_display_id",
        "display_id": display_id,
        "valid": bool(sbol3.is_valid_display_id(display_id)),
        "rule": _DISPLAY_ID_RULE,
        "provenance": _provenance(),
    }


def extract_sbol3_display_id(identity: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the display ID at the end of an SBOL3 URL identity.

    Identities that are not URLs have no display ID and return ``None``. A URL
    whose last path segment is not a valid display ID raises ``ConfigurationError``.
    """

    sbol3 = require("sbol3", "bio-pysbol3")
    if not isinstance(identity, str) or len(identity) > _MAX_IDENTITY:
        raise ConfigurationError(
            f"identity must be a string of at most {_MAX_IDENTITY} characters."
        )
    try:
        display_id = sbol3.extract_display_id(identity)
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    return {
        "status": "completed",
        "operation": "extract_sbol3_display_id",
        "identity": identity,
        "display_id": display_id,
        "found": display_id is not None,
        "provenance": _provenance(),
    }


__all__ = [
    "check_sbol3_display_id",
    "copy_sbol3_design",
    "extract_sbol3_display_id",
    "to_sbol3_display_id",
]

_register_public_verbose(globals(), tuple(__all__))
