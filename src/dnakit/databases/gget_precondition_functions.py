"""gget operations that need an optional package or an account: CELLxGENE expression queries and COSMIC references.

``gget_cellxgene_query`` runs gget's CELLxGENE Census query for explicit genes and metadata filters and writes the
result to an explicit ``.h5ad`` file (expression) or ``.csv`` file (metadata only). The optional ``cellxgene-census``
package is checked first, and a missing package raises ``BackendUnavailableError`` with the install hint.

``gget_cosmic_download_reference`` and ``gget_cosmic_select_reference`` download COSMIC reference archives with the
user's own COSMIC account. Requests go over HTTPS to ``cancer.sanger.ac.uk`` only, the credentials are sent only to
that host and are never written to a result, the archive is extracted with member checks, and no prompt is shown
(gget asks for input when the credentials are missing; these functions require them as arguments).

Results are JSON-serialisable dictionaries that record the service, the access date and the gget version.
"""

from __future__ import annotations

import base64
import gzip
import json
import re
import shutil
import tarfile
import urllib.request
from collections.abc import Sequence
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.databases.gget_gap_functions import (
    _document,
    _flag,
    _int,
    _new_folder,
    _new_path,
    _text,
    _upstream,
    _upstream_module,
)
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_CELLXGENE_EXTRA = "bio-cellxgene"
_CELLXGENE_SERVICE = "CELLxGENE Discover Census"
_CELLXGENE_SPECIES = (
    "homo_sapiens",
    "mus_musculus",
    "macaca_mulatta",
    "callithrix_jacchus",
    "pan_troglodytes",
)
_GENE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")
_LABEL = re.compile(r"[A-Za-z0-9][A-Za-z0-9 ._()/,+-]{0,127}")
_CENSUS_VERSION = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,31}")
_MAX_IDS = 500
_MAX_LABELS = 50
_MAX_LISTED = 200

_COSMIC_SERVICE = "COSMIC (Sanger Institute)"
_COSMIC_HOST = "cancer.sanger.ac.uk"
_COSMIC_LINK = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path={path}&bucket=downloads"
_EMAIL = re.compile(r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+")
# project -> (folder under the release folder, product name, contained file stem)
_COSMIC_PROJECTS = {
    "cancer": ("cmc", "CancerMutationCensus_AllData_Tsv", "CancerMutationCensus_AllData"),
    "cell_line": (
        "cell_lines",
        "CellLinesProject_GenomeScreensMutant_Tsv",
        "CellLinesProject_GenomeScreensMutant",
    ),
    "census": ("cosmic", "Cosmic_MutantCensus_Tsv", "Cosmic_MutantCensus"),
    "resistance": ("cosmic", "Cosmic_ResistanceMutations_Tsv", "Cosmic_ResistanceMutations"),
    "genome_screen": ("cosmic", "Cosmic_GenomeScreensMutant_Tsv", "Cosmic_GenomeScreensMutant"),
    "targeted_screen": (
        "cosmic",
        "Cosmic_CompleteTargetedScreensMutant_Tsv",
        "Cosmic_CompleteTargetedScreensMutant",
    ),
}
_MAX_LISTED_FILES = 200
_TIMEOUT_SECONDS = 600


# --- validation ------------------------------------------------------------------------------


def _names(values: object, name: str, pattern: re.Pattern[str], limit: int) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of strings.")
    if not values or len(values) > limit:
        raise ConfigurationError(f"{name} must hold between 1 and {limit} entries.")
    result: list[str] = []
    for item in values:
        if not isinstance(item, str) or pattern.fullmatch(item.strip()) is None:
            raise ConfigurationError(f"{name} has an unsupported entry: {item!r}.")
        result.append(item.strip())
    return result


def _optional_labels(values: object, name: str) -> list[str] | None:
    if values is None:
        return None
    return _names(values, name, _LABEL, _MAX_LABELS)


def _email(value: object) -> str:
    text = _text(value, "email", limit=254)
    if _EMAIL.fullmatch(text) is None:
        raise ConfigurationError("email is not a valid address.")
    return text


def _password(value: object) -> str:
    if not isinstance(value, str) or not value or len(value) > 512:
        raise ConfigurationError("password must be a non-empty string of at most 512 characters.")
    return value


def _https_cosmic(link: object) -> str:
    text = _text(link, "download_link", limit=2048)
    parts = urlsplit(text)
    if parts.scheme != "https" or parts.hostname != _COSMIC_HOST:
        raise ConfigurationError(f"download_link must be an https address on {_COSMIC_HOST}.")
    return text


def _safe_member(name: str) -> bool:
    path = Path(name)
    return not path.is_absolute() and ".." not in path.parts and name not in ("", ".")


# --- CELLxGENE -------------------------------------------------------------------------------


def gget_cellxgene_query(
    output_path: str,
    *,
    species: str = "homo_sapiens",
    genes: list[str] | None = None,
    ensembl: bool = False,
    tissue: list[str] | None = None,
    cell_type: list[str] | None = None,
    disease: list[str] | None = None,
    sex: list[str] | None = None,
    assay: list[str] | None = None,
    column_names: list[str] | None = None,
    is_primary_data: bool = True,
    meta_only: bool = False,
    census_version: str = "stable",
    verbose: bool = True,
) -> dict[str, Any]:
    """Query CELLxGENE Census expression for genes and cell metadata filters (gget.cellxgene).

    ``genes`` are gene symbols, or Ensembl IDs when ``ensembl`` is True. The metadata filters (``tissue``,
    ``cell_type``, ``disease``, ``sex``, ``assay``) are lists of labels as used in the Census schema. With
    ``meta_only`` True only the cell metadata table is written, as CSV; otherwise the expression AnnData is
    written as an ``.h5ad`` file. ``output_path`` must not exist yet. The Census is queried over the network,
    and the service version is recorded in the result.
    """
    target = _new_path(output_path, "output_path")
    flag_meta = _flag(meta_only, "meta_only")
    expected = ".csv" if flag_meta else ".h5ad"
    if target.suffix.lower() != expected:
        raise ConfigurationError(f"output_path must end in {expected} for this query.")
    if species not in _CELLXGENE_SPECIES:
        raise ConfigurationError(f"species must be one of: {', '.join(_CELLXGENE_SPECIES)}.")
    if not isinstance(census_version, str) or _CENSUS_VERSION.fullmatch(census_version) is None:
        raise ConfigurationError("census_version has an unsupported format.")
    gene_list = None if genes is None else _names(genes, "genes", _GENE, _MAX_IDS)
    columns = (
        None if column_names is None else _names(column_names, "column_names", _LABEL, _MAX_LABELS)
    )
    filters = {
        "tissue": _optional_labels(tissue, "tissue"),
        "cell_type": _optional_labels(cell_type, "cell_type"),
        "disease": _optional_labels(disease, "disease"),
        "sex": _optional_labels(sex, "sex"),
        "assay": _optional_labels(assay, "assay"),
    }
    use_ensembl = _flag(ensembl, "ensembl")
    primary = _flag(is_primary_data, "is_primary_data")
    flag_verbose = _flag(verbose, "verbose")

    require("cellxgene_census", _CELLXGENE_EXTRA)
    cellxgene = _upstream_module("gget_cellxgene").cellxgene
    result = _upstream(
        "cellxgene",
        lambda: cellxgene(
            species=species,
            gene=gene_list,
            ensembl=use_ensembl,
            column_names=columns,
            meta_only=flag_meta,
            tissue=filters["tissue"],
            cell_type=filters["cell_type"],
            disease=filters["disease"],
            sex=filters["sex"],
            assay=filters["assay"],
            is_primary_data=primary,
            census_version=census_version,
            verbose=flag_verbose,
            out=str(target),
        ),
    )
    if result is None:
        raise BackendExecutionError(
            "gget cellxgene returned no data; check the filters and the Census version."
        )
    if flag_meta:
        counts = {"rows": len(result)}
        obs_columns = [str(column) for column in result.columns]
        var_names: list[str] = []
    else:
        counts = {"cells": int(result.n_obs), "genes": int(result.n_vars)}
        obs_columns = [str(column) for column in result.obs.columns]
        var_names = [str(name) for name in list(result.var_names)[:_MAX_LISTED]]
    return _document(
        "cellxgene",
        _CELLXGENE_SERVICE,
        species=species,
        census_version=census_version,
        genes=gene_list,
        ensembl=use_ensembl,
        filters={key: value for key, value in filters.items() if value is not None},
        is_primary_data=primary,
        meta_only=flag_meta,
        output_path=str(target),
        counts=counts,
        obs_columns=obs_columns[:_MAX_LISTED],
        var_names=var_names,
    )


# --- COSMIC ----------------------------------------------------------------------------------


def _open_url(request: urllib.request.Request) -> Any:
    return urllib.request.urlopen(request, timeout=_TIMEOUT_SECONDS)


def _extract_archive(tar_path: Path, folder: Path) -> list[Path]:
    with tarfile.open(tar_path, "r") as archive:
        members = archive.getmembers()
        for member in members:
            if not _safe_member(member.name) or member.issym() or member.islnk():
                raise BackendExecutionError(
                    f"Unsafe member in the COSMIC archive: {member.name!r}."
                )
        if hasattr(tarfile, "data_filter"):
            archive.extractall(folder, members=members, filter="data")
        else:
            archive.extractall(folder, members=members)
    return sorted(path for path in folder.rglob("*") if path.is_file() and path != tar_path)


def gget_cosmic_download_reference(
    download_link: str,
    output_folder: str,
    *,
    email: str,
    password: str,
    unzip: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download one COSMIC reference archive with the user's account and extract it (gget.cosmic.download_reference).

    The link must be an https address on cancer.sanger.ac.uk. The first request sends the account as HTTP basic
    authentication and receives a signed download address; the archive is then fetched from that address without
    the credentials. The archive is saved as ``output_folder/<name>.tar`` and extracted into
    ``output_folder/<name>/``. With ``unzip`` True each extracted ``.gz`` file is decompressed next to it. The
    password is not stored or returned.
    """
    link = _https_cosmic(download_link)
    address = _email(email)
    secret = _password(password)
    flag_unzip = _flag(unzip, "unzip")
    _flag(verbose, "verbose")
    folder = _new_folder(output_folder, "output_folder")
    archive_name = Path(parse_qs(urlsplit(link).query).get("path", [""])[0]).name
    if not archive_name.endswith(".tar") or not _safe_member(archive_name):
        raise ConfigurationError("download_link must name a .tar archive in its path parameter.")
    token = base64.b64encode(f"{address}:{secret}\n".encode()).decode("ascii")
    request = urllib.request.Request(link, headers={"Authorization": f"Basic {token}"})
    try:
        with _open_url(request) as response:
            payload = json.loads(response.read(65536).decode("utf-8"))
    except (OSError, ValueError) as exc:
        raise BackendExecutionError(
            "COSMIC did not return a download address; check the account and the release."
        ) from exc
    signed = payload.get("url") if isinstance(payload, dict) else None
    if not isinstance(signed, str) or urlsplit(signed).scheme != "https":
        raise BackendExecutionError("COSMIC did not return a valid https download address.")
    tar_path = folder / archive_name
    try:
        with _open_url(urllib.request.Request(signed)) as response, tar_path.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except OSError as exc:
        raise BackendExecutionError(f"COSMIC archive download failed: {exc}") from exc
    extracted_folder = folder / archive_name[: -len(".tar")]
    extracted_folder.mkdir(exist_ok=True)
    files = _extract_archive(tar_path, extracted_folder)
    unzipped: list[str] = []
    if flag_unzip:
        for item in files:
            if item.suffix == ".gz":
                target = item.with_suffix("")
                with gzip.open(item, "rb") as source, target.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
                unzipped.append(str(target))
    return _document(
        "download_reference",
        _COSMIC_SERVICE,
        download_link=link,
        tar_path=str(tar_path),
        extracted_folder=str(extracted_folder),
        files=[str(item) for item in files[:_MAX_LISTED_FILES]],
        unzipped=unzipped,
    )


def gget_cosmic_select_reference(
    output_folder: str,
    cosmic_project: str,
    grch_version: int,
    cosmic_version: int,
    *,
    email: str,
    password: str,
    unzip: bool = True,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download one COSMIC project release for a genome build, or reuse it (gget.cosmic.select_reference).

    ``cosmic_project`` is one of ``cancer`` (CancerMutationCensus, GRCh37 only), ``cell_line``, ``census``,
    ``resistance``, ``genome_screen`` or ``targeted_screen``. ``grch_version`` is 37 or 38 and ``cosmic_version``
    is the release number. When the expected file already exists in ``output_folder`` it is reused unless
    ``overwrite`` is True. Otherwise the archive is downloaded with the account given by ``email`` and
    ``password`` (no prompt is shown), and the file path is returned.
    """
    project = _text(cosmic_project, "cosmic_project", limit=32)
    if project not in _COSMIC_PROJECTS:
        raise ConfigurationError(f"cosmic_project must be one of: {', '.join(_COSMIC_PROJECTS)}.")
    grch = _int(grch_version, "grch_version", 37, 38)
    release = _int(cosmic_version, "cosmic_version", 1, 200)
    if project == "cancer" and grch != 37:
        raise ConfigurationError(
            "CancerMutationCensus data is only available for GRCh37; use grch_version=37."
        )
    flag_unzip = _flag(unzip, "unzip")
    flag_overwrite = _flag(overwrite, "overwrite")
    _flag(verbose, "verbose")
    folder = _new_folder(output_folder, "output_folder")
    subdir, product, stem = _COSMIC_PROJECTS[project]
    prefix = "GRCh" if project == "cancer" and release in (99, 100) else "grch"
    suffix = "_v2" if project == "cancer" and release == 100 else ""
    name = f"{product}_v{release}_GRCh{grch}{suffix}"
    link = _COSMIC_LINK.format(path=f"{prefix}{grch}/{subdir}/v{release}/{name}.tar")
    contained = f"{stem}_v{release}_GRCh{grch}.tsv"
    expected = folder / name / (contained if flag_unzip else f"{contained}.gz")
    if expected.exists() and not flag_overwrite:
        return _document(
            "select_reference",
            _COSMIC_SERVICE,
            cosmic_project=project,
            grch_version=grch,
            cosmic_version=release,
            file_path=str(expected),
            downloaded=False,
            overwritten=False,
        )
    fetched = gget_cosmic_download_reference(
        link, str(folder), email=email, password=password, unzip=flag_unzip, verbose=verbose
    )
    if not expected.is_file():
        raise BackendExecutionError(f"The COSMIC archive did not contain {expected.name}.")
    return _document(
        "select_reference",
        _COSMIC_SERVICE,
        cosmic_project=project,
        grch_version=grch,
        cosmic_version=release,
        download_link=link,
        file_path=str(expected),
        downloaded=True,
        overwritten=flag_overwrite,
        tar_path=fetched["tar_path"],
        file_count=len(fetched["files"]),
    )


_register_public_verbose(globals())
