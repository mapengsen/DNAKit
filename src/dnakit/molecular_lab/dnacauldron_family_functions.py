"""DnaCauldron sequence readers that the family check found to have no existing wrapper (H2).

Each function imports DnaCauldron (and Biopython) only when it runs. A missing package
raises ``BackendUnavailableError`` with the ``dnakit[bio-dnacauldron]`` install hint.
Every function returns a JSON-safe dictionary with ``status``, ``operation``, the result
fields and ``provenance``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnacauldron"
_MAX_FILES = 200
_MAX_FILE_BYTES = 20_000_000
_MAX_RECORDS = 5000
_MAX_BASES = 2_000_000
_MAX_TEXT_BYTES = 1_000_000


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnacauldron")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnacauldron", "version": installed, **notes}


def _record_row(record: Any) -> dict[str, Any]:
    sequence = str(record.seq).upper()
    source = getattr(record, "source_file", None)
    return {
        "id": str(record.id),
        "file": Path(str(source)).name if source else None,
        "length": len(sequence),
        "topology": str(record.annotations.get("topology", "linear")),
        "sequence": sequence,
    }


def _source_files(sources: object, folder: object) -> list[Path]:
    paths: list[Path] = []
    if sources:
        if isinstance(sources, (str, bytes, Path)) or not isinstance(sources, Sequence):
            raise ConfigurationError("sources must be a list of local file paths.")
        for item in sources:
            if not isinstance(item, (str, Path)):
                raise ConfigurationError("sources must be a list of local file paths.")
            paths.append(input_path(item))
    if folder is not None:
        if not isinstance(folder, (str, Path)):
            raise ConfigurationError("folder must be a local directory path.")
        root = Path(folder).expanduser().resolve()
        if not root.is_dir():
            raise ConfigurationError(f"folder does not exist: {root}")
        paths.extend(sorted(path for path in root.rglob("*") if path.is_file()))
    if not paths:
        raise ConfigurationError("Provide at least one file in sources or a folder.")
    if len(paths) > _MAX_FILES:
        raise ConfigurationError(f"At most {_MAX_FILES} files can be read at once.")
    for path in paths:
        if path.stat().st_size > _MAX_FILE_BYTES:
            raise ConfigurationError(f"Input file exceeds {_MAX_FILE_BYTES} bytes: {path.name}")
    return paths


def read_dnacauldron_records(
    sources: Sequence[str | Path] = (),
    *,
    folder: str | Path | None = None,
    use_file_names_as_ids: bool = False,
    max_records: int = 500,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read every record of local FASTA, GenBank or SnapGene files with DnaCauldron.

    Unlike ``read_dnacauldron_parts``, a file may hold several records, the format is
    detected from the file content, ``.zip`` archives of sequence files are opened, and a
    whole ``folder`` can be read (its files, sorted by path). ``sources`` lists files.
    A record without an ID takes its file name, with a number when a file holds several
    records. ``use_file_names_as_ids`` gives a record the name of its file when the file
    holds one record. ``files`` reports the detected format of each file (``zip`` for an
    archive); ``records`` keeps the order of the files and of the records within them.
    """

    module = require("dnacauldron.biotools.sequence_io", _EXTRA)
    names_as_ids = _flag(use_file_names_as_ids, "use_file_names_as_ids")
    limit = _integer(max_records, "max_records", minimum=1, maximum=_MAX_RECORDS)
    paths = _source_files(sources, folder)
    records: list[Any] = []
    files: list[dict[str, Any]] = []
    for path in paths:
        is_zip = path.name.lower().endswith("zip")
        try:
            loaded = module.load_records_from_files(
                files=[str(path)], use_file_names_as_ids=names_as_ids
            )
            detected = "zip" if is_zip else str(module.load_records_from_file(str(path))[1])
        except (OSError, TypeError, ValueError, RuntimeError, KeyError) as exc:
            raise ConfigurationError(f"DnaCauldron could not read {path.name}: {exc}") from exc
        files.append({"file": path.name, "format": detected, "record_count": len(loaded)})
        records.extend(loaded)
        if len(records) > limit:
            raise ConfigurationError(f"The files hold more than {limit} records.")
    rows = [_record_row(record) for record in records]
    if sum(row["length"] for row in rows) > _MAX_BASES:
        raise ConfigurationError(f"Records exceed {_MAX_BASES} bases in total.")
    return {
        "status": "completed",
        "operation": "read_dnacauldron_records",
        "file_count": len(files),
        "record_count": len(rows),
        "files": files,
        "records": rows,
        "provenance": _provenance(
            method="dnacauldron.biotools.sequence_io.load_records_from_files",
            names_as_ids=names_as_ids,
        ),
    }


def parse_dnacauldron_sequence_text(text: str, *, verbose: bool = True) -> dict[str, Any]:
    """Detect the format of sequence text and read its records with DnaCauldron.

    The text is a bare ``ACGT`` string, FASTA or GenBank. ``format`` is ``"ATGC"`` for a
    bare string (one record without a name), ``"fasta"`` or ``"genbank"``. Unlike
    ``parse_sequence_text``, a bare sequence is accepted and the detected format is
    reported. Text in another format raises ``ConfigurationError``.
    """

    module = require("dnacauldron.biotools.sequence_io", _EXTRA)
    if not isinstance(text, str) or not text.strip():
        raise ConfigurationError("text must be a non-empty bare, FASTA or GenBank sequence.")
    if len(text.encode("utf-8")) > _MAX_TEXT_BYTES:
        raise ConfigurationError(f"text exceeds {_MAX_TEXT_BYTES} bytes.")
    try:
        parsed = module.string_to_records(text)
    except (OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaCauldron could not read the text: {exc}") from exc
    if not (isinstance(parsed, tuple) and len(parsed) == 2):
        raise ConfigurationError("DnaCauldron could not detect the format of the text.")
    records, detected = parsed
    rows = []
    for record in records:
        sequence = str(record.seq).upper()
        rows.append(
            {
                "id": str(record.id),
                "length": len(sequence),
                "topology": str(record.annotations.get("topology", "linear")),
                "sequence": sequence,
            }
        )
    return {
        "status": "completed",
        "operation": "parse_dnacauldron_sequence_text",
        "format": str(detected),
        "record_count": len(rows),
        "records": rows,
        "provenance": _provenance(method="dnacauldron.biotools.sequence_io.string_to_records"),
    }


__all__ = [
    "parse_dnacauldron_sequence_text",
    "read_dnacauldron_records",
]

_register_public_verbose(globals(), tuple(__all__))
