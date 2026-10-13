"""Offline genomepy helper for chromosome-name mapping between installed genome folders.

The mapping reads the README.txt and assembly_report.txt that genomepy saved when it installed a
genome. The wrapper checks that both files exist before it calls genomepy, so it never downloads
the NCBI assembly report. Readme bookkeeping, cache inventory and configuration helpers are not
wrapped here; the ledger records them as internal.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-genomepy]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "bio-genomepy"
_DISTRIBUTION = "genomepy"
_TARGETS = ("ucsc", "ncbi", "ensembl")
_GENOME_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}")
_NO_NAME = "na"


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version(_DISTRIBUTION)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"{_DISTRIBUTION} {version}", **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": provenance}


def _upstream(operation: str, call: Callable[[], Any]) -> Any:
    """Run one upstream call and report unexpected failures as backend errors."""
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"genomepy {operation} failed: {exc}") from exc


def _genome_name(value: object) -> str:
    text = _text(value, "genome_name", limit=128)
    if _GENOME_NAME.fullmatch(text) is None:
        raise ConfigurationError(
            "genome_name must be a folder name of letters, digits, '.', '_' or '-'."
        )
    return text


def _genomes_directory(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("genomes_dir must be a path to an existing directory.")
    path = Path(value).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"genomes_dir is not an existing directory: {path}")
    return path


def map_genome_contig_names(
    genome_name: str,
    target_provider: str,
    *,
    genomes_dir: str | Path,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map contig names of an installed genome to the naming scheme of another provider.

    Reads the installed genome's README.txt and assembly_report.txt only. A source and target
    provider that match, or a target whose names the report does not carry, give no mapping
    (``available`` is false). An ``na`` name in the report is returned as null.
    """
    name = _genome_name(genome_name)
    target = _text(target_provider, "target_provider", limit=16).lower()
    if target not in _TARGETS:
        raise ConfigurationError("target_provider must be 'ucsc', 'ncbi' or 'ensembl'.")
    root = _genomes_directory(genomes_dir)
    genome_dir = root / name
    readme = genome_dir / "README.txt"
    report = genome_dir / "assembly_report.txt"
    if not readme.is_file():
        raise BackendUnavailableError(
            f"No README.txt for genome {name!r} in {root}.",
            hint="Install the genome with genomepy, or pass the genomes directory that holds it.",
        )
    if not report.is_file():
        raise BackendUnavailableError(
            f"assembly_report.txt for genome {name!r} is not saved in {genome_dir}.",
            hint="genomepy saves this report when it installs a genome; this wrapper never "
            "downloads it.",
        )
    files = require("genomepy.files", _EXTRA)
    metadata, _ = _upstream("read_readme", lambda: files.read_readme(str(readme)))
    source = str(metadata.get("provider") or "").lower()
    if not source:
        raise BackendExecutionError(f"README.txt of genome {name!r} has no provider field.")
    provenance = _provenance(
        input="README.txt and assembly_report.txt saved in the genome folder",
        network="none; cached files only",
    )
    if source == target:
        return _document(
            "map_genome_contig_names",
            provenance,
            genome=name,
            source_provider=source,
            target_provider=target,
            available=False,
            reason="source and target provider are the same",
            record_count=0,
            records=[],
        )
    providers = require("genomepy.providers", _EXTRA)
    frame = _upstream(
        "map_locations",
        lambda: providers.map_locations(name, target, genomes_dir=str(root)),
    )
    if frame is None:
        return _document(
            "map_genome_contig_names",
            provenance,
            genome=name,
            source_provider=source,
            target_provider=target,
            available=False,
            reason="the assembly report has no names for the target provider",
            record_count=0,
            records=[],
        )
    column = str(frame.columns[0])
    values = frame[column].tolist()
    pairs = sorted({(str(key), str(value)) for key, value in zip(frame.index, values, strict=True)})
    records = [
        {"source_name": key, "target_name": None if value == _NO_NAME else value}
        for key, value in pairs
    ]
    return _document(
        "map_genome_contig_names",
        provenance,
        genome=name,
        source_provider=source,
        target_provider=target,
        available=True,
        reason=None,
        record_count=len(records),
        records=records,
    )


__all__ = ["map_genome_contig_names"]

_register_public_verbose(globals(), tuple(__all__))
