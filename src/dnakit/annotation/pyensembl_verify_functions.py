"""Install a gene annotation from local files into a pyensembl cache (the ``pyensembl install`` action).

``pyensembl.shell.run`` with the ``install`` action builds a ``Genome`` from explicit GTF and FASTA paths
(``--gtf``, ``--transcript-fasta``, ``--protein-fasta``, ``--genome-fasta-path``) and then calls
``genome.download`` and ``genome.index``. ``pyensembl_install_local_annotation`` does the same with the
pyensembl API and never reaches the network: every source must be a local file. The existing DNAKit
pyensembl wrappers only read installed caches.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_PYENSEMBL_EXTRA = "bio-pyensembl"
_MAX_FASTAS = 10


def _local_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing local file: {path}")
    return path


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 128:
        raise ConfigurationError(f"{name} must be a non-empty name of at most 128 characters.")
    return value.strip()


def _files(value: object, name: str) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of local file paths.")
    if len(value) > _MAX_FASTAS:
        raise ConfigurationError(f"{name} lists more than {_MAX_FASTAS} files.")
    return [str(_local_file(item, name)) for item in value]


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"


def pyensembl_install_local_annotation(
    reference_name: str,
    annotation_name: str,
    gtf_path: str | Path,
    *,
    annotation_version: int | None = None,
    transcript_fasta_paths: Sequence[str | Path] = (),
    protein_fasta_paths: Sequence[str | Path] = (),
    genome_fasta_path: str | Path | None = None,
    cache_directory_path: str | Path | None = None,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the pyensembl index for a local GTF annotation (and optional local FASTA sources).

    ``reference_name`` and ``annotation_name`` name the genome in pyensembl's cache. ``gtf_path`` is the
    annotation; ``transcript_fasta_paths``, ``protein_fasta_paths`` and ``genome_fasta_path`` are optional
    local sequence files. ``cache_directory_path`` (optional) replaces pyensembl's default cache folder,
    which is otherwise ``$PYENSEMBL_CACHE_DIR`` or ``~/.cache/pyensembl``. ``overwrite=True`` rebuilds an
    existing index. The result counts the indexed genes and transcripts and gives the install location.
    """
    pyensembl = require("pyensembl", _PYENSEMBL_EXTRA)
    reference = _text(reference_name, "reference_name")
    annotation = _text(annotation_name, "annotation_name")
    gtf = _local_file(gtf_path, "gtf_path")
    if annotation_version is not None:
        if isinstance(annotation_version, bool) or not isinstance(annotation_version, int) or annotation_version < 1:
            raise ConfigurationError("annotation_version must be a positive integer.")
    transcripts = _files(list(transcript_fasta_paths), "transcript_fasta_paths")
    proteins = _files(list(protein_fasta_paths), "protein_fasta_paths")
    genome_fasta = None if genome_fasta_path is None else str(_local_file(genome_fasta_path, "genome_fasta_path"))
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be true or false.")
    cache_dir = None
    if cache_directory_path is not None:
        cache_dir = Path(cache_directory_path).expanduser().resolve()
        if not cache_dir.is_dir():
            raise ConfigurationError(f"cache_directory_path must be an existing folder: {cache_dir}")
    try:
        genome = pyensembl.Genome(
            reference_name=reference,
            annotation_name=annotation,
            annotation_version=annotation_version,
            gtf_path_or_url=str(gtf),
            transcript_fasta_paths_or_urls=transcripts or None,
            protein_fasta_paths_or_urls=proteins or None,
            genome_fasta_path_or_url=genome_fasta,
            cache_directory_path=str(cache_dir) if cache_dir is not None else None,
        )
        genome.download(overwrite=overwrite, show_progress=False)
        genome.index(overwrite=overwrite, show_progress=False)
        installed = bool(genome.installed())
        gene_ids = list(genome.gene_ids())
        transcript_ids = list(genome.transcript_ids())
        cache_location = str(genome.download_cache.cache_directory_path)
    except (ValueError, TypeError, OSError) as exc:
        raise ConfigurationError(f"pyensembl could not install the annotation: {exc}") from exc
    return {
        "status": "installed" if installed else "not_installed",
        "operation": "pyensembl_install_local_annotation",
        "reference_name": reference,
        "annotation_name": annotation,
        "annotation_version": annotation_version,
        "gtf": gtf.name,
        "cache_directory": cache_location,
        "genes": len(gene_ids),
        "transcripts": len(transcript_ids),
        "provenance": {
            "package": f"pyensembl {_version('pyensembl')}",
            "convention": "Genome(local paths) then download and index, as pyensembl install does; no network",
        },
    }



def pyensembl_install_ensembl_release(
    release: int | str,
    species: str = "human",
    *,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download and index one Ensembl release for a species (the ``pyensembl install --release`` action).

    This needs network access to the Ensembl FTP server: pyensembl downloads the GTF annotation and the
    transcript and protein FASTA files, then builds the gene index in its cache. ``release`` is a numbered
    release (for example 110) or an annotation date in YYYY_MM form. ``species`` is pyensembl's species
    name, for example ``caenorhabditis_elegans`` or ``human``. The result records the release, the species
    and the Ensembl FTP server, and counts the indexed genes and transcripts.
    """
    pyensembl = require("pyensembl", _PYENSEMBL_EXTRA)
    if isinstance(release, bool) or not isinstance(release, (int, str)):
        raise ConfigurationError("release must be an Ensembl release number or a YYYY_MM date.")
    if isinstance(release, int) and not 1 <= release <= 1000:
        raise ConfigurationError("release must be between 1 and 1000.")
    species_name = _text(species, "species")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be true or false.")
    try:
        ensembl = pyensembl.EnsemblRelease(release=release, species=species_name)
        ensembl.download(overwrite=overwrite, show_progress=False)
        ensembl.index(overwrite=overwrite, show_progress=False)
        installed = bool(ensembl.installed())
        gene_count = len(list(ensembl.gene_ids()))
        transcript_count = len(list(ensembl.transcript_ids()))
        cache_location = str(ensembl.download_cache.cache_directory_path)
        gtf_url = str(ensembl.gtf_url)
    except (ValueError, TypeError, OSError) as exc:
        raise ConfigurationError(f"pyensembl could not install Ensembl release {release}: {exc}") from exc
    return {
        "status": "installed" if installed else "not_installed",
        "operation": "pyensembl_install_ensembl_release",
        "ensembl_release": release,
        "species": species_name,
        "gtf_url": gtf_url,
        "cache_directory": cache_location,
        "genes": gene_count,
        "transcripts": transcript_count,
        "provenance": {
            "package": f"pyensembl {_version('pyensembl')}",
            "convention": "EnsemblRelease download and index from the Ensembl FTP server (network)",
        },
    }

__all__ = ["pyensembl_install_local_annotation", "pyensembl_install_ensembl_release"]

_register_public_verbose(globals())
