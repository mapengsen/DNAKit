"""NCBI Datasets genome-reference download adapter."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
import zipfile
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, TypeAlias, cast
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from dnakit._version import __version__
from dnakit.core import (
    Citation,
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    OriginClass,
    Provenance,
    ReferenceInfo,
)
from dnakit.exceptions import ConfigurationError, DownloadError

from .models import DownloadProgress, GenomeAssembly, GenomeDownloadResult

DEFAULT_NCBI_API_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DEFAULT_TIMEOUT = 60.0
DEFAULT_CHUNK_SIZE = 1024 * 1024
DEFAULT_MAX_DOWNLOAD_BYTES = 20_000_000_000
_ACCESSION_PATTERN = re.compile(r"^GC[AF]_\d+\.\d+$", flags=re.IGNORECASE)

# These aliases are intentionally versioned.  Use an explicit accession when exact
# reproducibility matters; taxon names are resolved against NCBI's current report.
_FIXED_ASSEMBLY_ALIASES: dict[str, str] = {
    "hg38": "GCF_000001405.40",
    "grch38": "GCF_000001405.40",
    "grch38.p14": "GCF_000001405.40",
    "mm39": "GCF_000001635.27",
    "grcm39": "GCF_000001635.27",
    "tair10": "GCF_000001735.4",
}
_TAXON_ALIASES: dict[str, str] = {
    "arabidopsis": "Arabidopsis thaliana",
    "arabidopsis thaliana": "Arabidopsis thaliana",
    "arabidopsis_thaliana": "Arabidopsis thaliana",
    "human": "human",
    "homo sapiens": "Homo sapiens",
    "mouse": "mouse",
    "mice": "mouse",
    "mus musculus": "Mus musculus",
}
ProgressCallback: TypeAlias = Callable[[DownloadProgress], None]


def supported_genome_aliases() -> dict[str, str]:
    """Return the fixed aliases supported without a metadata lookup."""

    return dict(_FIXED_ASSEMBLY_ALIASES)


def _validate_positive(value: int | float, name: str) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ConfigurationError(f"{name} must be positive.", code="INVALID_DOWNLOAD_CONFIG")


def _new_digest(name: str) -> Any:
    if name == "md5":
        return hashlib.md5(usedforsecurity=False)
    return hashlib.new(name)


def _validate_common(
    *,
    timeout: float,
    chunk_size: int,
    max_download_bytes: int,
    api_base_url: str,
    api_key: str | None,
    progress: ProgressCallback | None,
) -> None:
    _validate_positive(timeout, "timeout")
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size < 1:
        raise ConfigurationError("chunk_size must be positive.", code="INVALID_DOWNLOAD_CONFIG")
    if (
        isinstance(max_download_bytes, bool)
        or not isinstance(max_download_bytes, int)
        or max_download_bytes < 1
    ):
        raise ConfigurationError(
            "max_download_bytes must be positive.", code="INVALID_DOWNLOAD_CONFIG"
        )
    if not isinstance(api_base_url, str) or not api_base_url.startswith("https://"):
        raise ConfigurationError(
            "api_base_url must be an HTTPS URL.", code="INVALID_DOWNLOAD_CONFIG"
        )
    if api_key is not None and (not isinstance(api_key, str) or not api_key.strip()):
        raise ConfigurationError(
            "api_key must be None or non-empty text.", code="INVALID_DOWNLOAD_CONFIG"
        )
    if progress is not None and not callable(progress):
        raise ConfigurationError(
            "progress must be callable or None.", code="INVALID_DOWNLOAD_CONFIG"
        )


def _request_headers(api_key: str | None) -> dict[str, str]:
    headers = {
        "Accept": "application/json",
        "User-Agent": f"DNAKit/{__version__} NCBI-reference-downloader",
    }
    if api_key is not None:
        headers["api-key"] = api_key
    return headers


def _request_json(
    url: str,
    *,
    timeout: float,
    api_key: str | None,
    max_bytes: int = 20_000_000,
) -> dict[str, Any]:
    request = Request(url, headers=_request_headers(api_key), method="GET")
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = _read_limited(response, max_bytes=max_bytes)
    except HTTPError as exc:
        raise DownloadError(
            "NCBI Datasets returned an HTTP error.",
            code="NCBI_HTTP_ERROR",
            context={"url": url, "status": exc.code},
        ) from exc
    except (OSError, TimeoutError, URLError) as exc:
        raise DownloadError(
            "Could not reach NCBI Datasets.",
            code="NCBI_NETWORK_ERROR",
            context={"url": url},
        ) from exc
    try:
        decoded = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise DownloadError(
            "NCBI Datasets returned invalid JSON.",
            code="NCBI_RESPONSE_ERROR",
            context={"url": url},
        ) from exc
    if not isinstance(decoded, dict):
        raise DownloadError(
            "NCBI Datasets returned an unexpected JSON shape.",
            code="NCBI_RESPONSE_ERROR",
            context={"url": url},
        )
    return cast(dict[str, Any], decoded)


def _read_limited(stream: BinaryIO, *, max_bytes: int) -> bytes:
    chunks: list[bytes] = []
    count = 0
    while True:
        chunk = stream.read(min(DEFAULT_CHUNK_SIZE, max_bytes - count + 1))
        if not chunk:
            break
        count += len(chunk)
        if count > max_bytes:
            raise DownloadError(
                "NCBI metadata response exceeds the configured limit.",
                code="NCBI_RESPONSE_SIZE_LIMIT",
                context={"max_bytes": max_bytes},
            )
        chunks.append(chunk)
    return b"".join(chunks)


def _join_url(api_base_url: str, path: str, params: Iterable[tuple[str, str]]) -> str:
    query = urlencode(tuple(params), doseq=True)
    return f"{api_base_url.rstrip('/')}/{path.lstrip('/')}" + (f"?{query}" if query else "")


def _resolve_taxon_alias(query: str) -> str:
    return _TAXON_ALIASES.get(query.casefold(), query)


def resolve_genome_assembly(
    query: str,
    *,
    timeout: float = DEFAULT_TIMEOUT,
    api_key: str | None = None,
    api_base_url: str = DEFAULT_NCBI_API_BASE_URL,
) -> GenomeAssembly:
    """Resolve an accession, fixed alias, or taxon name to one NCBI assembly.

    Taxon resolution selects one current RefSeq reference assembly.  If NCBI
    returns multiple candidates, callers must use an explicit accession.
    """

    if not isinstance(query, str) or not query.strip():
        raise ConfigurationError("query must be non-empty text.", code="INVALID_GENOME_QUERY")
    _validate_common(
        timeout=timeout,
        chunk_size=DEFAULT_CHUNK_SIZE,
        max_download_bytes=DEFAULT_MAX_DOWNLOAD_BYTES,
        api_base_url=api_base_url,
        api_key=api_key,
        progress=None,
    )
    normalized = query.strip()
    if _ACCESSION_PATTERN.fullmatch(normalized):
        return GenomeAssembly(query=normalized, accession=normalized.upper())
    fixed_accession = _FIXED_ASSEMBLY_ALIASES.get(normalized.casefold())
    if fixed_accession is not None:
        return GenomeAssembly(query=normalized, accession=fixed_accession)

    taxon = _resolve_taxon_alias(normalized)
    url = _join_url(
        api_base_url,
        f"/genome/taxon/{quote(taxon, safe='')}/dataset_report",
        (
            ("filters.reference_only", "true"),
            ("filters.assembly_source", "refseq"),
            ("filters.assembly_version", "current"),
            ("filters.exclude_paired_reports", "true"),
            ("filters.exclude_atypical", "true"),
            ("tax_exact_match", "true"),
            ("returned_content", "ASSM_ACC"),
            ("page_size", "1000"),
        ),
    )
    payload = _request_json(url, timeout=timeout, api_key=api_key)
    reports = payload.get("reports")
    if not isinstance(reports, list):
        raise DownloadError(
            "NCBI Datasets taxon response did not contain reports.",
            code="NCBI_RESPONSE_ERROR",
            context={"query": query},
        )
    accessions: list[str] = []
    for item in reports:
        if not isinstance(item, Mapping):
            continue
        candidate = item.get("accession")
        if isinstance(candidate, str):
            accessions.append(candidate)
    if not accessions:
        raise ConfigurationError(
            f"No current RefSeq reference assembly was found for {query!r}.",
            code="GENOME_ASSEMBLY_NOT_FOUND",
            context={"query": query},
            hint="Use an explicit GCA_/GCF_ assembly accession.",
        )
    if len(accessions) != 1:
        raise ConfigurationError(
            f"Taxon {query!r} resolved to multiple reference assemblies.",
            code="AMBIGUOUS_GENOME_ASSEMBLY",
            context={"query": query, "accessions": accessions},
            hint="Use an explicit GCA_/GCF_ assembly accession.",
        )
    accession = accessions[0]
    if _ACCESSION_PATTERN.fullmatch(accession) is None:
        raise DownloadError(
            "NCBI returned an invalid assembly accession.",
            code="NCBI_RESPONSE_ERROR",
            context={"query": query, "accession": accession},
        )
    source_database = None
    selected_report = next(
        (
            item
            for item in reports
            if isinstance(item, Mapping) and item.get("accession") == accession
        ),
        None,
    )
    if isinstance(selected_report, Mapping) and isinstance(
        selected_report.get("source_database"), str
    ):
        source_database = selected_report["source_database"]
    return GenomeAssembly(
        query=normalized,
        accession=accession,
        source_database=source_database,
    )


def _content_length(response: Any) -> int | None:
    headers = getattr(response, "headers", None)
    if headers is None:
        return None
    value = headers.get("Content-Length")
    if value is None:
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 0 else None


def _download_package(
    url: str,
    target: Path,
    *,
    timeout: float,
    api_key: str | None,
    chunk_size: int,
    max_download_bytes: int,
    assembly: GenomeAssembly,
    progress: ProgressCallback | None,
) -> tuple[int, str]:
    request = Request(url, headers={**_request_headers(api_key), "Accept": "application/zip"})
    try:
        with urlopen(request, timeout=timeout) as response:
            total_bytes = _content_length(response)
            if total_bytes is not None and total_bytes > max_download_bytes:
                raise DownloadError(
                    "The NCBI genome package exceeds max_download_bytes.",
                    code="DOWNLOAD_SIZE_LIMIT",
                    context={
                        "accession": assembly.accession,
                        "content_length": total_bytes,
                        "max_download_bytes": max_download_bytes,
                    },
                )
            completed = 0
            digest = hashlib.sha256()
            with target.open("wb") as output:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    completed += len(chunk)
                    if completed > max_download_bytes:
                        raise DownloadError(
                            "The NCBI genome package exceeds max_download_bytes.",
                            code="DOWNLOAD_SIZE_LIMIT",
                            context={
                                "accession": assembly.accession,
                                "downloaded_bytes": completed,
                                "max_download_bytes": max_download_bytes,
                            },
                        )
                    output.write(chunk)
                    digest.update(chunk)
                    if progress is not None:
                        progress(
                            DownloadProgress(
                                phase="download",
                                query=assembly.query,
                                accession=assembly.accession,
                                bytes_completed=completed,
                                total_bytes=total_bytes,
                            )
                        )
            if completed < 1:
                raise DownloadError(
                    "NCBI returned an empty genome package.",
                    code="EMPTY_DOWNLOAD",
                    context={"accession": assembly.accession},
                )
            return completed, digest.hexdigest()
    except DownloadError:
        raise
    except HTTPError as exc:
        raise DownloadError(
            "NCBI Datasets returned an HTTP error while downloading the genome.",
            code="NCBI_HTTP_ERROR",
            context={"url": url, "status": exc.code},
        ) from exc
    except (OSError, TimeoutError, URLError) as exc:
        raise DownloadError(
            "The NCBI genome package could not be downloaded.",
            code="NCBI_NETWORK_ERROR",
            context={"url": url},
        ) from exc


def _safe_member_name(name: str) -> str:
    if not isinstance(name, str) or not name or "\x00" in name or "\\" in name:
        raise DownloadError("NCBI package contains an unsafe member name.", code="UNSAFE_ARCHIVE")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise DownloadError("NCBI package contains an unsafe member path.", code="UNSAFE_ARCHIVE")
    return "/".join(path.parts)


def _validate_member(info: zipfile.ZipInfo) -> str:
    name = _safe_member_name(info.filename)
    mode = (info.external_attr >> 16) & 0o170000
    if mode == 0o120000:
        raise DownloadError("NCBI package contains a symbolic link.", code="UNSAFE_ARCHIVE")
    if info.is_dir():
        raise DownloadError(
            "Expected NCBI package file but found a directory.", code="INVALID_ARCHIVE"
        )
    return name


def _find_member(infos: list[zipfile.ZipInfo], filename: str) -> zipfile.ZipInfo:
    matches = [info for info in infos if PurePosixPath(info.filename).name == filename]
    if len(matches) != 1:
        raise DownloadError(
            f"NCBI package must contain exactly one {filename} file.",
            code="INVALID_ARCHIVE",
            context={"matches": len(matches)},
        )
    _validate_member(matches[0])
    return matches[0]


def _copy_member(
    archive: zipfile.ZipFile,
    info: zipfile.ZipInfo,
    target: Path,
    *,
    chunk_size: int,
    assembly: GenomeAssembly,
    progress: ProgressCallback | None,
    hash_algorithms: tuple[str, ...] = (),
) -> dict[str, str]:
    digests = {name: _new_digest(name) for name in hash_algorithms}
    completed = 0
    try:
        with archive.open(info, "r") as source, target.open("wb") as output:
            while True:
                chunk = source.read(chunk_size)
                if not chunk:
                    break
                output.write(chunk)
                completed += len(chunk)
                for digest in digests.values():
                    digest.update(chunk)
                if progress is not None:
                    progress(
                        DownloadProgress(
                            phase="extract",
                            query=assembly.query,
                            accession=assembly.accession,
                            bytes_completed=completed,
                            total_bytes=info.file_size,
                            current_file=PurePosixPath(info.filename).name,
                        )
                    )
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise DownloadError(
            "Could not extract the NCBI genome package.",
            code="ARCHIVE_EXTRACTION_ERROR",
            context={"member": info.filename},
        ) from exc
    return {name: digest.hexdigest() for name, digest in digests.items()}


def _read_member_text(archive: zipfile.ZipFile, info: zipfile.ZipInfo) -> str:
    try:
        with archive.open(info, "r") as source:
            payload = source.read(20_000_001)
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise DownloadError(
            "Could not read an NCBI package metadata file.",
            code="ARCHIVE_EXTRACTION_ERROR",
            context={"member": info.filename},
        ) from exc
    if len(payload) > 20_000_000:
        raise DownloadError(
            "NCBI package metadata exceeds the configured limit.",
            code="ARCHIVE_METADATA_SIZE_LIMIT",
            context={"member": info.filename},
        )
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise DownloadError(
            "NCBI package metadata is not UTF-8.",
            code="ARCHIVE_METADATA_ERROR",
            context={"member": info.filename},
        ) from exc


def _extract_report_details(text: str, accession: str) -> tuple[str | None, str | None]:
    for line in text.splitlines():
        if not line.strip():
            continue
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            raise DownloadError(
                "NCBI assembly_data_report.jsonl is invalid JSON.",
                code="ARCHIVE_METADATA_ERROR",
            ) from exc
        if not isinstance(payload, Mapping):
            continue
        reported_accession = payload.get("accession")
        assembly = payload.get("assembly")
        organism = payload.get("organism")
        if isinstance(assembly, Mapping):
            reported_accession = assembly.get("accession", reported_accession)
        if reported_accession is not None and reported_accession != accession:
            continue
        assembly_name = None
        if isinstance(assembly, Mapping):
            for key in ("name", "assembly_name"):
                if isinstance(assembly.get(key), str) and assembly[key].strip():
                    assembly_name = assembly[key]
                    break
        for key in ("assembly_name", "assminfo_name"):
            if assembly_name is None and isinstance(payload.get(key), str):
                assembly_name = payload[key]
        organism_name = None
        if isinstance(organism, Mapping):
            for key in ("organism_name", "name"):
                if isinstance(organism.get(key), str) and organism[key].strip():
                    organism_name = organism[key]
                    break
        for key in ("organism_name", "organism"):
            if organism_name is None and isinstance(payload.get(key), str):
                organism_name = payload[key]
        return organism_name, assembly_name
    return None, None


def _expected_md5(text: str, member_name: str) -> str:
    normalized_member = member_name.lstrip("./")
    basename = PurePosixPath(member_name).name
    for line in text.splitlines():
        fields = line.strip().split(maxsplit=1)
        if len(fields) != 2:
            continue
        digest, path = fields
        normalized_path = path.lstrip("*").lstrip("./")
        if (
            normalized_path == normalized_member or PurePosixPath(normalized_path).name == basename
        ) and re.fullmatch(r"[0-9a-fA-F]{32}", digest):
            return digest.lower()
    raise DownloadError(
        "NCBI md5sum.txt does not contain the genomic FASTA entry.",
        code="CHECKSUM_MISSING",
        context={"member": member_name},
    )


def _install_staged(staged_to_target: Iterable[tuple[Path, Path]], *, overwrite: bool) -> None:
    pairs = tuple(staged_to_target)
    targets = tuple(target for _, target in pairs)
    if len(targets) != len(set(targets)):
        raise ConfigurationError("Download output paths collide.", code="OUTPUT_COLLISION")
    existing = [target for target in targets if target.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing download outputs: {existing}")
    if any(target.exists() and not target.is_file() for target in targets):
        raise IsADirectoryError("Download output targets must be files.")
    installed: list[Path] = []
    backups: list[tuple[Path, Path]] = []
    try:
        if overwrite:
            backup_paths = [
                target.with_name(f".{target.name}.dnakit-backup-{index}")
                for index, target in enumerate(targets)
                if target.exists()
            ]
            if any(path.exists() for path in backup_paths):
                raise FileExistsError("A stale DNAKit download backup already exists.")
        for index, (staged, target) in enumerate(pairs):
            if overwrite:
                if target.exists():
                    backup = target.with_name(f".{target.name}.dnakit-backup-{index}")
                    os.replace(target, backup)
                    backups.append((target, backup))
                os.replace(staged, target)
            else:
                os.link(staged, target)
                staged.unlink()
            installed.append(target)
    except BaseException:
        for target in installed:
            target.unlink(missing_ok=True)
        for target, backup in reversed(backups):
            if backup.exists():
                os.replace(backup, target)
        raise
    for _, backup in backups:
        backup.unlink(missing_ok=True)


def download_genome(
    query: str,
    output_dir: str | os.PathLike[str],
    *,
    overwrite: bool = False,
    keep_package: bool = False,
    timeout: float = DEFAULT_TIMEOUT,
    api_key: str | None = None,
    api_base_url: str = DEFAULT_NCBI_API_BASE_URL,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    max_download_bytes: int = DEFAULT_MAX_DOWNLOAD_BYTES,
    progress: ProgressCallback | None = None,
) -> GenomeDownloadResult:
    """Download one NCBI genomic FASTA and its provenance files.

    ``query`` may be an explicit GCA_/GCF_ assembly accession, a fixed alias
    such as ``hg38`` or ``mm39``, or a taxon name.  Taxon names are resolved
    to one current RefSeq reference assembly and fail when ambiguous.
    """

    _validate_common(
        timeout=timeout,
        chunk_size=chunk_size,
        max_download_bytes=max_download_bytes,
        api_base_url=api_base_url,
        api_key=api_key,
        progress=progress,
    )
    if not isinstance(overwrite, bool) or not isinstance(keep_package, bool):
        raise ConfigurationError(
            "overwrite and keep_package must be booleans.", code="INVALID_DOWNLOAD_CONFIG"
        )
    assembly = resolve_genome_assembly(
        query,
        timeout=timeout,
        api_key=api_key,
        api_base_url=api_base_url,
    )
    if progress is not None:
        progress(
            DownloadProgress(
                phase="resolve",
                query=assembly.query,
                accession=assembly.accession,
                bytes_completed=0,
            )
        )
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(f"output_dir must be a directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    url = _join_url(
        api_base_url,
        f"/genome/accession/{quote(assembly.accession, safe='')}/download",
        (("include_annotation_type", "GENOME_FASTA"), ("hydrated", "FULLY_HYDRATED")),
    )

    with tempfile.TemporaryDirectory(prefix=".dnakit-ncbi-", dir=str(output)) as raw_work:
        work = Path(raw_work)
        package = work / "ncbi_dataset.zip"
        downloaded_bytes, package_sha256 = _download_package(
            url,
            package,
            timeout=timeout,
            api_key=api_key,
            chunk_size=chunk_size,
            max_download_bytes=max_download_bytes,
            assembly=assembly,
            progress=progress,
        )
        try:
            with zipfile.ZipFile(package, "r") as archive:
                infos = archive.infolist()
                for info in infos:
                    _safe_member_name(info.filename)
                fasta_info = next(
                    (
                        info
                        for info in infos
                        if PurePosixPath(info.filename).name.endswith("_genomic.fna")
                    ),
                    None,
                )
                if fasta_info is None:
                    raise DownloadError(
                        "NCBI package does not contain a genomic FASTA file.",
                        code="INVALID_ARCHIVE",
                    )
                fasta_member_name = _validate_member(fasta_info)
                report_info = _find_member(infos, "assembly_data_report.jsonl")
                checksum_info = _find_member(infos, "md5sum.txt")
                report_text = _read_member_text(archive, report_info)
                checksum_text = _read_member_text(archive, checksum_info)
                expected_md5 = _expected_md5(checksum_text, fasta_member_name)
                staged_fasta = work / PurePosixPath(fasta_member_name).name
                digests = _copy_member(
                    archive,
                    fasta_info,
                    staged_fasta,
                    chunk_size=chunk_size,
                    assembly=assembly,
                    progress=progress,
                    hash_algorithms=("md5", "sha256"),
                )
                if digests["md5"] != expected_md5:
                    raise DownloadError(
                        "The downloaded genomic FASTA failed the NCBI MD5 check.",
                        code="CHECKSUM_MISMATCH",
                        context={"expected_md5": expected_md5, "actual_md5": digests["md5"]},
                    )
                organism, assembly_name = _extract_report_details(report_text, assembly.accession)
                staged_report = work / f"{assembly.accession}_assembly_data_report.jsonl"
                staged_report.write_text(report_text, encoding="utf-8")
                staged_checksum = work / f"{assembly.accession}_md5sum.txt"
                staged_checksum.write_text(checksum_text, encoding="utf-8")
                staged_package = package
                fasta_target = output / staged_fasta.name
                report_target = output / staged_report.name
                checksum_target = output / staged_checksum.name
                package_target = output / f"{assembly.accession}_ncbi_dataset.zip"
                pairs: list[tuple[Path, Path]] = [
                    (staged_fasta, fasta_target),
                    (staged_report, report_target),
                    (staged_checksum, checksum_target),
                ]
                if keep_package:
                    pairs.append((staged_package, package_target))
                _install_staged(pairs, overwrite=overwrite)
        except DownloadError:
            raise
        except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
            raise DownloadError(
                "The NCBI genome package is invalid or could not be extracted.",
                code="INVALID_ARCHIVE",
                context={"accession": assembly.accession},
            ) from exc

    return GenomeDownloadResult(
        query=assembly.query,
        accession=assembly.accession,
        organism=organism,
        assembly_name=assembly_name,
        output_directory=str(output),
        fasta_path=str(fasta_target),
        metadata_path=str(report_target),
        checksum_path=str(checksum_target),
        package_path=str(package_target) if keep_package else None,
        download_url=url,
        downloaded_bytes=downloaded_bytes,
        package_sha256=package_sha256,
        fasta_sha256=digests["sha256"],
        fasta_md5=digests["md5"],
        provenance=Provenance(
            implementation=ImplementationInfo(
                label=ImplementationLabel.ADAPTER,
                execution_mode=ExecutionMode.EXTERNAL,
                origin_class=OriginClass.INTEGRATION,
                citations=(
                    Citation(
                        "ncbi-datasets",
                        title="NCBI Datasets Genome Package",
                        url=(
                            "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/"
                            "reference-docs/data-packages/genome/"
                        ),
                    ),
                ),
            ),
            reference=ReferenceInfo(
                name="NCBI genome assembly",
                version=assembly.accession,
                checksum=digests["sha256"],
                filters={
                    "query": assembly.query,
                    "include_annotation_type": "GENOME_FASTA",
                    "hydrated": "FULLY_HYDRATED",
                },
            ),
        ),
    )


__all__ = [
    "DEFAULT_CHUNK_SIZE",
    "DEFAULT_MAX_DOWNLOAD_BYTES",
    "DEFAULT_NCBI_API_BASE_URL",
    "DEFAULT_TIMEOUT",
    "ProgressCallback",
    "download_genome",
    "resolve_genome_assembly",
    "supported_genome_aliases",
]
