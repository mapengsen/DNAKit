"""pyensembl operations not covered by the existing DNAKit pyensembl helpers (Fig. 2 coverage gap).

Listing functions call pyensembl's dated-release listings, which read Ensembl's new-platform directory listings
over the network; pyensembl reuses cached listings unless a refresh is requested. The cache functions open an
annotation that is already installed in the local pyensembl cache and never download anything: when the annotation
is missing they raise BackendUnavailableError whose hint names the pyensembl install command. Downloads, index
builds and the pyensembl command-line entry points are not wrapped (see the disposition file).
"""

from __future__ import annotations

import datetime as _dt
import math
import re
from collections.abc import Callable, Mapping
from datetime import timedelta
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "bio-pyensembl"
_ACCESSION = re.compile(r"GC[AF]_[0-9]{9}\.[0-9]+")
_PROVIDER = re.compile(r"[a-z][a-z_]{0,31}")
_SERVICE_LISTING = "Ensembl new-platform directory listings (pyensembl dated releases)"
_SERVICE_CACHE = "local pyensembl annotation cache (read only)"


def _module(name: str) -> Any:
    return require(name, _EXTRA)


def _retrieved() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _text(value: object, name: str, *, limit: int = 128) -> str:
    if not isinstance(value, str) or not value.strip() or len(value.strip()) > limit:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return value.strip()


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _accession(value: object) -> str:
    text = _text(value, "assembly_accession", limit=32)
    if _ACCESSION.fullmatch(text) is None:
        raise ConfigurationError(
            "assembly_accession must be a versioned GCA or GCF accession, for example GCA_000001405.29."
        )
    return text


def _provider(value: object) -> str:
    text = _text(value, "provider", limit=32)
    if _PROVIDER.fullmatch(text) is None:
        raise ConfigurationError(
            "provider must be a lower-case provider directory name, for example ensembl."
        )
    return text


def _release(value: object) -> int | str:
    """A numbered Ensembl release (positive integer) or an annotation date such as 2026_04."""
    if isinstance(value, bool):
        raise ConfigurationError(
            "release must be a positive integer or an annotation date such as 2026_04."
        )
    if isinstance(value, int) and value >= 1:
        return value
    if isinstance(value, str) and value.strip():
        return value.strip()
    raise ConfigurationError(
        "release must be a positive integer or an annotation date such as 2026_04."
    )


def _nonnegative_number(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0
    ):
        raise ConfigurationError(f"{name} must be a finite number of zero or more.")
    return float(value)


def _positive_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _upstream(name: str, call: Callable[[], Any]) -> Any:
    """Run one pyensembl call; bad names and values become ConfigurationError, failed reads BackendExecutionError."""
    try:
        return call()
    except (ConfigurationError, BackendUnavailableError):
        raise
    except (
        ValueError
    ) as exc:  # pyensembl reports unknown species, references and releases as ValueError
        raise ConfigurationError(f"pyensembl {name} rejected the input: {exc}") from exc
    except OSError as exc:  # an unreachable listing, or a listing that names no dates
        raise BackendExecutionError(f"pyensembl {name} could not read the listing: {exc}") from exc


def _listing_result(
    operation: str, request: Mapping[str, Any], data: Mapping[str, Any]
) -> dict[str, Any]:
    pyensembl = _module("pyensembl")
    return {
        "status": "complete",
        "operation": operation,
        "service": _SERVICE_LISTING,
        "backend": "pyensembl",
        "version": getattr(pyensembl, "__version__", None),
        "request": dict(request),
        "retrieved_utc": _retrieved(),
        "data": dict(data),
    }


def list_ensembl_annotation_dates(
    assembly_accession: str,
    *,
    provider: str = "ensembl",
    refresh: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the annotation dates Ensembl publishes for one assembly accession (pyensembl available_annotation_dates).

    Needs network access the first time a listing is read; later calls reuse the cached listing unless
    refresh is true. Dates are returned oldest first in YYYY_MM form.
    """
    accession = _accession(assembly_accession)
    provider_name = _provider(provider)
    refresh_flag = _flag(refresh, "refresh")
    dated = _module("pyensembl.dated_releases")
    dates = _upstream(
        "available_annotation_dates",
        lambda: dated.available_annotation_dates(accession, provider_name, refresh=refresh_flag),
    )
    listed = [str(date) for date in dates]
    return _listing_result(
        "list_ensembl_annotation_dates",
        {"assembly_accession": accession, "provider": provider_name, "refresh": refresh_flag},
        {"dates": listed, "count": len(listed)},
    )


def list_ensembl_dated_releases(
    species: str = "human",
    *,
    refresh: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the annotation dates published for a species' current assembly (pyensembl available_dated_releases).

    The species name is resolved in pyensembl's registry first; a species without dated releases is refused
    before any network access. Needs network access unless a cached listing exists and refresh is false.
    """
    name = _text(species, "species")
    refresh_flag = _flag(refresh, "refresh")
    species_module = _module("pyensembl.species")
    resolved = _upstream("find_species_by_name", lambda: species_module.check_species_object(name))
    if resolved.dated_releases is None:
        raise ConfigurationError(f"{resolved.latin_name} has no dated Ensembl releases.")
    accession, provider = resolved.dated_releases
    dated = _module("pyensembl.dated_releases")
    dates = _upstream(
        "available_dated_releases",
        lambda: dated.available_dated_releases(name, refresh=refresh_flag),
    )
    listed = [str(date) for date in dates]
    return _listing_result(
        "list_ensembl_dated_releases",
        {"species": name, "refresh": refresh_flag},
        {
            "latin_name": resolved.latin_name,
            "assembly_accession": accession,
            "provider": provider,
            "dates": listed,
            "count": len(listed),
        },
    )


def fetch_all_ensembl_dated_releases(
    *,
    expire_after_days: float = 1.0,
    timeout_seconds: int = 10,
    verbose: bool = True,
) -> dict[str, Any]:
    """Annotation dates of every species with dated releases, keyed by latin name (pyensembl fetch_all_dated_releases).

    Cached listings younger than expire_after_days are used as they are. Human is fetched first; if its listing
    cannot be read, every species falls back to its cached dates. A species whose dates are unknown maps to None.
    Needs network access for listings that are missing or expired.
    """
    expire_days = _nonnegative_number(expire_after_days, "expire_after_days")
    timeout = _positive_int(timeout_seconds, "timeout_seconds")
    dated = _module("pyensembl.dated_releases")
    found = _upstream(
        "fetch_all_dated_releases",
        lambda: dated.fetch_all_dated_releases(
            expire_after=timedelta(days=expire_days), timeout=timeout
        ),
    )
    if found is None:
        raise BackendExecutionError(
            "pyensembl could not fetch any dated-release listing and none are cached."
        )
    by_species = {
        str(latin): (None if dates is None else [str(date) for date in dates])
        for latin, dates in sorted(found.items())
    }
    unknown = sorted(latin for latin, dates in by_species.items() if dates is None)
    return _listing_result(
        "fetch_all_ensembl_dated_releases",
        {"expire_after_days": expire_days, "timeout_seconds": timeout},
        {
            "dates_by_species": by_species,
            "species_count": len(by_species),
            "species_with_dates": len(by_species) - len(unknown),
            "species_without_dates": unknown,
        },
    )


def _install_hint(genome: Any) -> str:
    return (
        "Install the annotation with the pyensembl command line first, for example: "
        f"pyensembl install --release {genome.release} --species {genome.species.latin_name}; "
        "then call this function again."
    )


def _installed_release(operation: str, genome: Any, request: Mapping[str, Any]) -> dict[str, Any]:
    """Describe an EnsemblRelease that is fully installed locally; never downloads anything."""
    report = genome.inspect_data()  # reads the cache only: no download, no index build
    if not report["installed"]:
        raise BackendUnavailableError(
            f"The {genome.species.latin_name} {genome.reference_name} release {genome.release} is not "
            f"installed in the local pyensembl cache (annotation: {report['annotation'] or 'not configured'}, "
            f"reference DNA: {report['reference_dna'] or 'not configured'}).",
            hint=_install_hint(genome),
        )
    files = {
        role: {"status": inspection.status, "path": str(inspection.path)}
        for role, inspection in report["files"].items()
    }
    return {
        "status": "complete",
        "operation": operation,
        "service": _SERVICE_CACHE,
        "backend": "pyensembl",
        "version": getattr(_module("pyensembl"), "__version__", None),
        "request": dict(request),
        "retrieved_utc": _retrieved(),
        "data": {
            "species": genome.species.latin_name,
            "reference_name": genome.reference_name,
            "release": str(genome.release),
            "annotation_name": genome.annotation_name,
            "annotation_version": str(genome.annotation_version),
            "cache_directory": str(genome.download_cache.cache_directory_path),
            "installed": True,
            "annotation_status": report["annotation"],
            "reference_dna_status": report["reference_dna"],
            "gene_count": len(genome.gene_ids()),
            "contig_count": len(genome.contigs()),
            "files": files,
        },
    }


def cached_ensembl_release(
    release: int | str,
    *,
    species: str = "human",
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe an installed Ensembl release from the local pyensembl cache (pyensembl cached_release).

    Reads only the local cache: gene and contig counts come from the installed annotation index. If the
    release is not installed, BackendUnavailableError is raised with the pyensembl install command as the hint.
    """
    release_value = _release(release)
    species_name = _text(species, "species")
    ensembl_release = _module("pyensembl.ensembl_release")
    genome = _upstream(
        "cached_release",
        lambda: ensembl_release.cached_release(release_value, species=species_name),
    )
    return _installed_release(
        "cached_ensembl_release",
        genome,
        {"release": release_value, "species": species_name},
    )


def ensembl_release_for_reference_name(
    reference_name: str,
    *,
    allow_older_downloaded_release: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Describe the installed Ensembl release for a reference name such as GRCh38 (pyensembl genome_for_reference_name).

    With allow_older_downloaded_release true, the newest installed release is chosen; otherwise the newest
    release with downloaded files. Nothing is downloaded: if the chosen release is not installed,
    BackendUnavailableError is raised with the pyensembl install command as the hint.
    """
    reference = _text(reference_name, "reference_name")
    allow_older = _flag(allow_older_downloaded_release, "allow_older_downloaded_release")
    reference_module = _module("pyensembl.reference_name")
    genome = _upstream(
        "genome_for_reference_name",
        lambda: reference_module.genome_for_reference_name(
            reference, allow_older_downloaded_release=allow_older
        ),
    )
    return _installed_release(
        "ensembl_release_for_reference_name",
        genome,
        {"reference_name": reference, "allow_older_downloaded_release": allow_older},
    )


__all__ = [
    "cached_ensembl_release",
    "ensembl_release_for_reference_name",
    "fetch_all_ensembl_dated_releases",
    "list_ensembl_annotation_dates",
    "list_ensembl_dated_releases",
]

_register_public_verbose(globals())
