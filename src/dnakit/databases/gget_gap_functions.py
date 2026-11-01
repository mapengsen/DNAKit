"""gget operations not covered by gget_functions or gene_backends: Ensembl reference downloads, NCBI Virus and
GenBank data, cBioPortal study data and heatmap, and two alignment programs (MUSCLE 5 and DIAMOND).

Network functions call one gget upstream function with the arguments named here and return a JSON-serialisable
dictionary that records the service, the access date and the gget version. The remote service is contacted at
call time; nothing is cached by DNAKit. Output locations are passed explicitly, and an output file that already
exists is refused. Temporary files that gget itself creates are gget's, not DNAKit's.

External programs: ``gget_muscle_align`` runs MUSCLE 5 (``muscle -align`` or ``-super5``) and
``gget_diamond_align`` runs DIAMOND (``makedb`` then ``blastp``/``blastx`` with a fixed ``outfmt 6`` column list).
Each executable is taken from ``muscle_path``/``diamond_path``, else from PATH. No shell is used. A missing
program raises ``BackendUnavailableError`` with an install hint. gget ships Linux builds of both programs under
``gget/bins``; they can be passed as the explicit path.

Secrets: API keys are passed to the upstream call only and are never returned or logged.
"""

from __future__ import annotations

import datetime as _dt
import importlib
import importlib.metadata
import math
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import require
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    BackendVersionError,
    ConfigurationError,
)

_GGET = "bio-gget"
_MUSCLE_HINT = (
    "Install MUSCLE 5 (conda install -c bioconda muscle=5, or a release binary from "
    "https://github.com/rcedgar/muscle/releases); put it on PATH or pass muscle_path. "
    "gget ships a Linux build under gget/bins/Linux/muscle."
)
_DIAMOND_HINT = (
    "Install DIAMOND 2.1 or later (conda install -c bioconda diamond, or a release binary from "
    "https://github.com/bbuchfink/diamond/releases); put it on PATH or pass diamond_path. "
    "gget ships a Linux build under gget/bins/Linux/diamond."
)
_SENSITIVITY = (
    "fast",
    "mid-sensitive",
    "sensitive",
    "more-sensitive",
    "very-sensitive",
    "ultra-sensitive",
)
_DIAMOND_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "qlen",
    "slen",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
_INTEGER_COLUMNS = frozenset(
    {"qlen", "slen", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"}
)
_REFERENCE_KINDS = ("gtf", "cdna", "dna", "cds", "cdrna", "pep")
_IDENTIFIER = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/+-]{0,63}")
_SPECIES = re.compile(r"[a-z][a-z0-9_]{1,80}")
_GENE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:-]{0,63}")
_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
_VIRUS = re.compile(r"[A-Za-z0-9][A-Za-z0-9 ._()/+-]{0,127}")
_MAX_TEXT = 10_000
_MAX_IDS = 10_000
_MAX_GENES = 200
_MAX_ROWS = 100_000
_MAX_TAIL = 4_000
_MAX_LENGTH = 1_000_000_000
_MAX_THREADS = 256
_DIAMOND_SERVICE = "DIAMOND local binary"
_MUSCLE_SERVICE = "MUSCLE 5 local binary"


# --- validation, conversion and provenance -----------------------------------------------------


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return value.strip()


def _pattern(value: object, name: str, pattern: re.Pattern[str]) -> str:
    text = _text(value, name)
    if pattern.fullmatch(text) is None:
        raise ConfigurationError(f"{name} has an unsupported format: {text!r}.")
    return text


def _identifiers(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list of identifiers.")
    if not values or len(values) > _MAX_IDS:
        raise ConfigurationError(f"{name} must hold between 1 and {_MAX_IDS} identifiers.")
    return [_pattern(item, f"{name}[]", _IDENTIFIER) for item in values]


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _optional_int(value: object, name: str, minimum: int, maximum: int) -> int | None:
    return None if value is None else _int(value, name, minimum, maximum)


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _optional_flag(value: object, name: str) -> bool | None:
    return None if value is None else _flag(value, name)


def _optional_date(value: object, name: str) -> str | None:
    return None if value is None else _pattern(value, name, _DATE)


def _plain(value: Any) -> Any:
    """Make upstream results JSON-serialisable: no NumPy scalars, NaN or frames."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_plain(item) for item in value]
    if hasattr(value, "to_dict"):
        return _plain(value.to_dict(orient="records"))
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    return str(value)


def _gget_version() -> str:
    try:
        return importlib.metadata.version("gget")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _document(operation: str, service: str, **fields: Any) -> dict[str, Any]:
    return {
        "operation": operation,
        "status": "completed",
        "backend": "gget",
        "package": f"gget {_gget_version()}",
        "service": service,
        "accessed": _dt.date.today().isoformat(),
        **{key: _plain(value) for key, value in fields.items()},
    }


def _upstream(name: str, call: Any) -> Any:
    """Run one gget call; upstream failures become BackendExecutionError, configuration errors pass through."""
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:  # gget raises a mix of requests, ValueError and SystemExit-free errors
        raise BackendExecutionError(f"gget {name} failed: {exc}") from exc


def _upstream_module(name: str) -> Any:
    """Import one gget submodule (for example ``gget_virus``) after checking that gget is installed."""
    require("gget", _GGET)
    return importlib.import_module(f"gget.{name}")


def _new_path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"{name} already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {target.parent}")
    return target


def _new_folder(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a folder path.")
    target = Path(value).expanduser().resolve()
    if target.exists() and not target.is_dir():
        raise ConfigurationError(f"{name} exists and is not a folder: {target}")
    target.mkdir(parents=True, exist_ok=True)
    return target


def _choice(value: object) -> str:
    if value not in _SENSITIVITY:
        raise ConfigurationError(f"sensitivity must be one of: {', '.join(_SENSITIVITY)}.")
    return str(value)


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    target = Path(value).expanduser().resolve()
    if not target.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {target}")
    return target


# --- network: cBioPortal ---------------------------------------------------------------------


def gget_cbio_plot(
    study_ids: Sequence[str],
    genes: Sequence[str],
    *,
    figure_dir: str,
    figure_filename: str,
    data_dir: str,
    stratification: str = "tissue",
    variation_type: str = "mutation_occurrences",
    confirm_download: bool = False,
    dpi: int = 100,
    figure_title: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download cBioPortal study data and draw a heatmap of gene alterations per stratum (gget.cbio_plot).

    ``data_dir`` receives the downloaded study files and ``figure_dir/figure_filename`` the PNG. Both are explicit,
    so nothing is written to the working directory. ``stratification`` and ``variation_type`` are passed to gget.
    """
    studies = _identifiers(study_ids, "study_ids")
    if (
        isinstance(genes, (str, bytes))
        or not isinstance(genes, Sequence)
        or not genes
        or len(genes) > _MAX_GENES
    ):
        raise ConfigurationError(f"genes must hold between 1 and {_MAX_GENES} gene symbols.")
    symbols = [_pattern(gene, "genes[]", _GENE) for gene in genes]
    strat = _text(stratification, "stratification", limit=64)
    variation = _text(variation_type, "variation_type", limit=64)
    name = _text(figure_filename, "figure_filename", limit=128)
    if "/" in name or "\\" in name or name in (".", ".."):
        raise ConfigurationError("figure_filename must be a plain file name.")
    figures = _new_folder(figure_dir, "figure_dir")
    data = _new_folder(data_dir, "data_dir")
    dots = _int(dpi, "dpi", 30, 1_200)
    title = None if figure_title is None else _text(figure_title, "figure_title", limit=256)
    gget = require("gget", _GGET)
    ok = _upstream(
        "cbio_plot",
        lambda: gget.cbio_plot(
            studies,
            symbols,
            stratification=strat,
            variation_type=variation,
            data_dir=str(data),
            figure_dir=str(figures),
            figure_filename=name,
            verbose=_flag(verbose, "verbose"),
            confirm_download=_flag(confirm_download, "confirm_download"),
            dpi=dots,
            show=False,
            figure_title=title,
        ),
    )
    if not ok:
        raise BackendExecutionError(
            "gget cbio_plot reported failure; see the gget log for the cause."
        )
    figure = figures / name
    return _document(
        "cbio_plot",
        "cBioPortal (www.cbioportal.org) via gget",
        study_ids=studies,
        genes=symbols,
        figure=str(figure) if figure.exists() else None,
        data_dir=str(data),
    )


def gget_cbio_download(
    study_ids: Sequence[str],
    *,
    out_dir: str,
    confirm_download: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download cBioPortal study archives (mutations, copy-number, structural variants, clinical) into ``out_dir``."""
    studies = _identifiers(study_ids, "study_ids")
    target = _new_folder(out_dir, "out_dir")
    cbio = _upstream_module("gget_cbio")
    ok = _upstream(
        "download_cbioportal_data",
        lambda: cbio.download_cbioportal_data(
            studies,
            verbose=False,  # DNAKit reports progress itself; keep gget quiet
            out_dir=str(target),
            confirm_download=_flag(confirm_download, "confirm_download"),
        ),
    )
    if not ok:
        raise BackendExecutionError(
            "gget reported that not all cBioPortal study files were downloaded."
        )
    return _document(
        "download_cbioportal_data",
        "cBioPortal study downloads via gget",
        study_ids=studies,
        out_dir=str(target),
    )


# --- network: Ensembl reference files ----------------------------------------------------------


def gget_ref(
    species: str,
    *,
    which: Sequence[str] | str = "all",
    release: int | None = None,
    ftp: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the Ensembl FTP links for a species' reference files (gget.ref), without downloading them.

    ``species`` is ``genus_species`` (for example ``homo_sapiens``) or a shortcut such as ``human``.
    ``which`` is ``all`` or a subset of gtf, cdna, dna, cds, cdrna and pep. ``ftp=True`` returns only the links.
    """
    name = _pattern(species, "species", _SPECIES)
    if isinstance(which, str):
        kinds: list[str] = [_text(which, "which", limit=16)]
    else:
        kinds = [_text(item, "which[]", limit=16) for item in which]
    if kinds != ["all"] and any(kind not in _REFERENCE_KINDS for kind in kinds):
        raise ConfigurationError(
            f"which must be 'all' or a subset of: {', '.join(_REFERENCE_KINDS)}."
        )
    gget = require("gget", _GGET)
    rel = _optional_int(release, "release", 1, 1_000)
    links = _upstream(
        "ref",
        lambda: gget.ref(
            name,
            which=kinds[0] if len(kinds) == 1 else kinds,
            release=rel,
            ftp=_flag(ftp, "ftp"),
            save=False,
            verbose=_flag(verbose, "verbose"),
        ),
    )
    return _document(
        "ref",
        "Ensembl FTP (ftp.ensembl.org) via gget",
        species=name,
        which=kinds,
        release=rel,
        ftp_only=_flag(ftp, "ftp"),
        links=links,
    )


# --- network: NCBI Virus and GenBank -----------------------------------------------------------


def gget_virus(
    virus: str,
    *,
    outfolder: str,
    is_accession: bool = False,
    host: str | None = None,
    min_seq_length: int | None = None,
    max_seq_length: int | None = None,
    nuc_completeness: str | None = None,
    annotated: bool | None = None,
    geographic_location: str | None = None,
    min_release_date: str | None = None,
    max_release_date: str | None = None,
    lineage: str | None = None,
    genbank_metadata: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a virus genome dataset (sequences and metadata) from NCBI Virus (gget.virus) into ``outfolder``.

    ``virus`` is a taxon name such as ``"Zika virus"`` or, with ``is_accession=True``, an accession list. Filters
    are passed to gget unchanged; unset filters stay unset.
    """
    taxon = _text(virus, "virus", limit=256)
    if not is_accession and _VIRUS.fullmatch(taxon) is None:
        raise ConfigurationError(
            "virus must be a taxon name (letters, digits and common punctuation)."
        )
    folder = _new_folder(outfolder, "outfolder")
    lengths = (
        _optional_int(min_seq_length, "min_seq_length", 0, _MAX_LENGTH),
        _optional_int(max_seq_length, "max_seq_length", 0, _MAX_LENGTH),
    )
    if lengths[0] is not None and lengths[1] is not None and lengths[0] > lengths[1]:
        raise ConfigurationError("min_seq_length must not exceed max_seq_length.")
    gget = require("gget", _GGET)
    result = _upstream(
        "virus",
        lambda: gget.virus(
            taxon,
            is_accession=_flag(is_accession, "is_accession"),
            outfolder=str(folder),
            host=host,
            min_seq_length=lengths[0],
            max_seq_length=lengths[1],
            nuc_completeness=nuc_completeness,
            annotated=_optional_flag(annotated, "annotated"),
            geographic_location=geographic_location,
            min_release_date=_optional_date(min_release_date, "min_release_date"),
            max_release_date=_optional_date(max_release_date, "max_release_date"),
            lineage=lineage,
            genbank_metadata=_flag(genbank_metadata, "genbank_metadata"),
            verbose=_flag(verbose, "verbose"),
        ),
    )
    return _document(
        "virus",
        "NCBI Virus (www.ncbi.nlm.nih.gov/labs/virus) via gget",
        virus=taxon,
        outfolder=str(folder),
        result=result,
    )


def gget_virus_metadata(
    virus: str,
    *,
    accession: bool = False,
    host: str | None = None,
    geographic_location: str | None = None,
    annotated: bool | None = None,
    complete_only: bool | None = None,
    min_release_date: str | None = None,
    refseq_only: bool | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch virus genome metadata from the NCBI Datasets API without downloading sequences (gget.gget_virus)."""
    taxon = _text(virus, "virus", limit=256)
    upstream = _upstream_module("gget_virus")
    rows = _upstream(
        "fetch_virus_metadata",
        lambda: upstream.fetch_virus_metadata(
            taxon,
            accession=_flag(accession, "accession"),
            host=host,
            geographic_location=geographic_location,
            annotated=_optional_flag(annotated, "annotated"),
            complete_only=_optional_flag(complete_only, "complete_only"),
            min_release_date=_optional_date(min_release_date, "min_release_date"),
            refseq_only=_optional_flag(refseq_only, "refseq_only"),
        ),
    )
    return _document(
        "fetch_virus_metadata",
        "NCBI Datasets API (api.ncbi.nlm.nih.gov/datasets) via gget",
        virus=taxon,
        records=rows,
        count=len(rows) if isinstance(rows, (list, tuple)) else None,
    )


def gget_virus_metadata_chunked(
    virus: str,
    *,
    accession: bool = False,
    host: str | None = None,
    geographic_location: str | None = None,
    annotated: bool | None = None,
    complete_only: bool = False,
    min_release_date: str | None = None,
    max_release_date: str | None = None,
    refseq_only: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch metadata for a very large virus dataset in release-date chunks (gget.gget_virus, chunked strategy)."""
    taxon = _text(virus, "virus", limit=256)
    upstream = _upstream_module("gget_virus")
    rows = _upstream(
        "fetch_virus_metadata_chunked",
        lambda: upstream.fetch_virus_metadata_chunked(
            taxon,
            accession=_flag(accession, "accession"),
            host=host,
            geographic_location=geographic_location,
            annotated=_optional_flag(annotated, "annotated"),
            complete_only=_flag(complete_only, "complete_only"),
            min_release_date=_optional_date(min_release_date, "min_release_date"),
            max_release_date=_optional_date(max_release_date, "max_release_date"),
            refseq_only=_flag(refseq_only, "refseq_only"),
        ),
    )
    return _document(
        "fetch_virus_metadata_chunked",
        "NCBI Datasets API (api.ncbi.nlm.nih.gov/datasets) via gget",
        virus=taxon,
        records=rows,
        count=len(rows) if isinstance(rows, (list, tuple)) else None,
    )


def gget_sars_cov2_download(
    *,
    outdir: str,
    host: str | None = None,
    complete_only: bool | None = None,
    annotated: bool | None = None,
    lineage: str | None = None,
    accession: str | None = None,
    use_accession: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download SARS-CoV-2 sequences from NCBI's cached data packages (gget.gget_virus, optimized route).

    Returns the ZIP path with the filters NCBI applied and the filters it could not apply.
    """
    folder = _new_folder(outdir, "outdir")
    upstream = _upstream_module("gget_virus")
    zip_path, applied, missing = _upstream(
        "download_sars_cov2_optimized",
        lambda: upstream.download_sars_cov2_optimized(
            host=host,
            complete_only=_optional_flag(complete_only, "complete_only"),
            annotated=_optional_flag(annotated, "annotated"),
            outdir=str(folder),
            lineage=lineage,
            accession=accession,
            use_accession=_flag(use_accession, "use_accession"),
        ),
    )
    return _document(
        "download_sars_cov2_optimized",
        "NCBI Virus cached data packages via gget",
        outdir=str(folder),
        zip_path=zip_path,
        applied_filters=applied,
        missing_filters=missing,
    )


def gget_alphainfluenza_download(
    *,
    outdir: str,
    host: str | None = None,
    complete_only: bool | None = None,
    annotated: bool | None = None,
    accession: str | None = None,
    use_accession: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download influenza A sequences from NCBI's cached data packages (gget.gget_virus, optimized route)."""
    folder = _new_folder(outdir, "outdir")
    upstream = _upstream_module("gget_virus")
    zip_path, applied, missing = _upstream(
        "download_alphainfluenza_optimized",
        lambda: upstream.download_alphainfluenza_optimized(
            host=host,
            complete_only=_optional_flag(complete_only, "complete_only"),
            annotated=_optional_flag(annotated, "annotated"),
            outdir=str(folder),
            accession=accession,
            use_accession=_flag(use_accession, "use_accession"),
        ),
    )
    return _document(
        "download_alphainfluenza_optimized",
        "NCBI Virus cached data packages via gget",
        outdir=str(folder),
        zip_path=zip_path,
        applied_filters=applied,
        missing_filters=missing,
    )


def gget_virus_sequences_by_accession(
    accessions: Sequence[str],
    *,
    outdir: str,
    batch_size: int = 200,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download virus genome sequences for a list of accession numbers from NCBI (gget.gget_virus, batched)."""
    ids = _identifiers(accessions, "accessions")
    folder = _new_folder(outdir, "outdir")
    size = _int(batch_size, "batch_size", 1, 1_000)
    key = None if api_key is None else _text(api_key, "api_key", limit=256)
    upstream = _upstream_module("gget_virus")
    fasta_path = _upstream(
        "download_sequences_by_accessions",
        lambda: upstream.download_sequences_by_accessions(
            ids, outdir=str(folder), batch_size=size, api_key=key
        ),
    )
    return _document(
        "download_sequences_by_accessions",
        "NCBI Virus sequence download (E-utilities) via gget",
        accessions=len(ids),
        outdir=str(folder),
        fasta_path=fasta_path,
    )


def gget_genbank_metadata(
    accessions: Sequence[str],
    *,
    xml_path: str,
    csv_path: str,
    batch_size: int = 200,
    delay: float = 0.5,
    failed_log_path: str | None = None,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fetch GenBank records for accessions through NCBI E-utilities and write the XML and a metadata CSV.

    ``xml_path`` and ``csv_path`` must not exist yet. ``delay`` seconds separate E-utilities requests.
    """
    ids = _identifiers(accessions, "accessions")
    xml_target = _new_path(xml_path, "xml_path")
    csv_target = _new_path(csv_path, "csv_path")
    failed = None if failed_log_path is None else _new_path(failed_log_path, "failed_log_path")
    size = _int(batch_size, "batch_size", 1, 1_000)
    if (
        not isinstance(delay, (int, float))
        or isinstance(delay, bool)
        or not 0 <= float(delay) <= 60
    ):
        raise ConfigurationError("delay must be a number of seconds from 0 to 60.")
    key = None if api_key is None else _text(api_key, "api_key", limit=256)
    upstream = _upstream_module("gget_virus")
    summary = _upstream(
        "fetch_genbank_metadata",
        lambda: upstream.fetch_genbank_metadata(
            ids,
            str(xml_target),
            str(csv_target),
            batch_size=size,
            delay=float(delay),
            failed_log_path=None if failed is None else str(failed),
            api_key=key,
        ),
    )
    return _document(
        "fetch_genbank_metadata",
        "NCBI E-utilities (efetch/esummary) via gget",
        accessions=len(ids),
        xml_path=str(xml_target),
        csv_path=str(csv_target),
        summary=summary,
    )


# --- external programs ---------------------------------------------------------------------------


def _program(explicit: str | None, name: str, hint: str) -> str:
    if explicit is not None:
        candidate = Path(_text(explicit, f"{name}_path", limit=4096)).expanduser()
        if not candidate.is_file():
            raise BackendUnavailableError(f"{name}_path does not exist: {candidate}", hint=hint)
        return str(candidate)
    found = shutil.which(name)
    if found is None:
        raise BackendUnavailableError(f"{name} was not found on PATH.", hint=hint)
    return found


def _run(program: str, args: list[str], *, backend: str, cwd: Path) -> tuple[int, str]:
    result = execute_bounded_command(program, args, backend_id=backend, cwd=cwd, verbose=False)
    return result.return_code, (result.output or "")[-_MAX_TAIL:]


def _version_text(program: str, args: list[str], *, backend: str, cwd: Path) -> str:
    code, output = _run(program, args, backend=backend, cwd=cwd)
    if code != 0:
        raise BackendExecutionError(f"{backend} version check failed: {output.strip()[:300]}")
    return output.strip().splitlines()[0] if output.strip() else "unknown"


def _fasta_summary(path: Path) -> dict[str, Any]:
    """Record count, aligned length (all characters, gaps included) and identifiers of an aligned FASTA."""
    records: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if text.startswith(">"):
                header = text[1:].split()
                records.append({"id": header[0] if header else "", "length": 0})
            elif records:
                records[-1]["length"] += len(text)
    lengths = {record["length"] for record in records}
    return {
        "records": len(records),
        "aligned_length": lengths.pop() if len(lengths) == 1 else None,
        "ids": [record["id"] for record in records[:_MAX_ROWS]],
    }


def gget_muscle_align(
    fasta: str,
    output_path: str,
    *,
    super5: bool = False,
    muscle_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align nucleotide or amino-acid sequences with MUSCLE 5 and write an aligned FASTA (gget.muscle).

    ``super5=True`` uses the Super5 algorithm for large inputs (a few hundred sequences). The command is
    ``muscle -align <fasta> -output <output_path>`` (or ``-super5``), the same as gget's.
    """
    source = _existing_file(fasta, "fasta")
    target = _new_path(output_path, "output_path")
    program = _program(muscle_path, "muscle", _MUSCLE_HINT)
    version = _version_text(program, ["-version"], backend="muscle", cwd=target.parent)
    if re.search(r"\b5\.\d", version) is None:
        raise BackendVersionError(f"MUSCLE 5 is required (found: {version}).")
    option = "-super5" if _flag(super5, "super5") else "-align"
    code, tail = _run(
        program, [option, str(source), "-output", str(target)], backend="muscle", cwd=target.parent
    )
    if code != 0 or not target.is_file():
        raise BackendExecutionError(f"MUSCLE exited with code {code}: {tail.strip()[-300:]}")
    summary = _fasta_summary(target)
    return _document(
        "muscle_align",
        _MUSCLE_SERVICE,
        program={"name": "muscle", "version": version, "path": program},
        command=["muscle", option, str(source), "-output", str(target)],
        output=str(target),
        records=summary["records"],
        aligned_length=summary["aligned_length"],
        ids=summary["ids"],
        log_tail=tail,
    )


def _diamond_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(_DIAMOND_COLUMNS):
                raise BackendExecutionError(
                    "DIAMOND wrote a row with an unexpected number of columns."
                )
            row: dict[str, Any] = {}
            for name, value in zip(_DIAMOND_COLUMNS, fields, strict=True):
                if name in _INTEGER_COLUMNS:
                    row[name] = int(value)
                elif name == "qseqid" or name == "sseqid":
                    row[name] = value
                else:
                    row[name] = float(value)
            rows.append(row)
            if len(rows) >= _MAX_ROWS:
                break
    return rows


def gget_diamond_align(
    query: str,
    reference: str,
    output_path: str,
    *,
    translated: bool = False,
    sensitivity: str = "very-sensitive",
    threads: int = 1,
    diamond_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Align query sequences to reference sequences with DIAMOND and write the tabular (outfmt 6) hits (gget.diamond).

    ``translated=False`` runs ``blastp`` (protein query against protein reference); ``translated=True`` runs
    ``blastx`` (nucleotide query against protein reference). The reference database is built into a temporary
    folder next to ``output_path`` and removed afterwards. The result lists every hit with the fixed columns.
    """
    query_path = _existing_file(query, "query")
    reference_path = _existing_file(reference, "reference")
    target = _new_path(output_path, "output_path")
    mode = _choice(sensitivity)
    workers = _int(threads, "threads", 1, _MAX_THREADS)
    program = _program(diamond_path, "diamond", _DIAMOND_HINT)
    version = _version_text(program, ["version"], backend="diamond", cwd=target.parent)
    program_name = "blastx" if _flag(translated, "translated") else "blastp"
    with tempfile.TemporaryDirectory(prefix="dnakit-diamond-", dir=target.parent) as scratch:
        database = Path(scratch) / "reference"
        code, tail = _run(
            program,
            [
                "makedb",
                "--quiet",
                "--in",
                str(reference_path),
                "--db",
                str(database),
                "--threads",
                str(workers),
            ],
            backend="diamond",
            cwd=Path(scratch),
        )
        if code != 0:
            raise BackendExecutionError(
                f"DIAMOND makedb exited with code {code}: {tail.strip()[-300:]}"
            )
        args = [
            program_name,
            "--outfmt",
            "6",
            *_DIAMOND_COLUMNS,
            "--quiet",
            "--query",
            str(query_path),
            "--db",
            f"{database}.dmnd",
            "--out",
            str(target),
            f"--{mode}",
            "--threads",
            str(workers),
        ]
        code, tail = _run(program, args, backend="diamond", cwd=Path(scratch))
    if code != 0:
        raise BackendExecutionError(
            f"DIAMOND {program_name} exited with code {code}: {tail.strip()[-300:]}"
        )
    rows = _diamond_rows(target) if target.is_file() else []
    return _document(
        "diamond_align",
        _DIAMOND_SERVICE,
        program={"name": "diamond", "version": version, "path": program},
        mode=program_name,
        sensitivity=mode,
        output=str(target),
        columns=list(_DIAMOND_COLUMNS),
        hits=len(rows),
        rows=rows,
        log_tail=tail,
    )


__all__ = [
    "gget_alphainfluenza_download",
    "gget_cbio_download",
    "gget_cbio_plot",
    "gget_diamond_align",
    "gget_genbank_metadata",
    "gget_muscle_align",
    "gget_ref",
    "gget_sars_cov2_download",
    "gget_virus",
    "gget_virus_metadata",
    "gget_virus_metadata_chunked",
    "gget_virus_sequences_by_accession",
]
_register_public_verbose(globals())
