"""genomepy operations not covered by the existing DNAKit genomepy helpers (Fig. 2 coverage gap).

Catalogue, assembly-report and annotation wrappers call genomepy provider functions at run time, so
they need network access; tests replace those functions. URL reachability probes and generic
transport helpers are not wrapped (see the disposition file).

Aligner index builds run as explicit argument vectors without a shell. genomepy's own
``run_index_cmd`` passes a string to a shell; the wrappers here accept only the aligner's index
program as the first argument.
"""

from __future__ import annotations

import contextlib
import importlib.metadata
import io
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, require, run
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_EXTRA = "genomepy"
_INDEX_PROGRAMS = {
    "bowtie2": "bowtie2-build",
    "bwa": "bwa",
    "gmap": "gmap_build",
    "hisat2": "hisat2-build",
    "minimap2": "minimap2",
    "star": "STAR",
}
_UCSC_TOOLS_HINT = (
    "Install the UCSC conversion tools (for example conda install -c bioconda ucsc-gtftogenepred "
    "ucsc-genepredtobed ucsc-bedtobigbed) and put them on PATH."
)
_LOG_TAIL = 2000
_SPECIES = ("human", "mouse")


def _module(name: str) -> Any:
    return require(name, _EXTRA)


def _version() -> str | None:
    try:
        return importlib.metadata.version("genomepy")
    except importlib.metadata.PackageNotFoundError:
        return None


def _text(value: object, name: str, *, limit: int = 512) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return value.strip()


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a file path.")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} does not exist: {path}")
    return path


def _fresh_path(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new path: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _positive_or_none(value: object, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(f"{name} must be a positive integer or None.")
    return value


def _plain(value: Any) -> Any:
    """Convert pandas, numpy and nested containers to JSON-compatible Python values."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    item = getattr(value, "item", None)
    if callable(item) and not isinstance(value, (str, bytes)):
        try:
            return item()
        except (TypeError, ValueError):
            return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _frame_rows(frame: Any) -> list[dict[str, Any]]:
    return [_plain(record) for record in frame.to_dict(orient="records")]


def _result(operation: str, outputs: Mapping[str, object] | None = None, **extra: Any) -> dict:
    return {
        "status": "complete",
        "operation": operation,
        "backend": "genomepy",
        "version": _version(),
        "outputs": {key: str(value) for key, value in (outputs or {}).items()},
        **extra,
    }


def _require_ucsc_tools(module: Any) -> None:
    """genomepy's check_ucsc_tools raises OSError; report it as a missing backend instead."""
    try:
        module.check_ucsc_tools()
    except OSError as exc:
        raise BackendUnavailableError(
            f"UCSC conversion tools are not installed: {exc}", hint=_UCSC_TOOLS_HINT
        ) from exc


def genomepy_list_available_genomes(
    *,
    provider: str | None = None,
    include_size: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List genomes offered by the online providers (``genomepy.list_available_genomes``).

    Each row is ``[name, provider, description...]``. Needs network access.
    """
    functions = _module("genomepy.functions")
    rows = [_plain(row) for row in functions.list_available_genomes(provider, bool(include_size))]
    return _result("list_available_genomes", rows=rows, row_count=len(rows))


def genomepy_search_genomes(
    term: str,
    *,
    provider: str | None = None,
    exact: bool = False,
    include_size: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search genome names and descriptions across providers (``genomepy.search``).

    Needs network access.
    """
    query = _text(term, "term") if isinstance(term, str) else str(term)
    providers = _module("genomepy.providers")
    rows = [
        _plain(row) for row in providers.search(query, provider, bool(exact), bool(include_size))
    ]
    return _result("search_genomes", rows=rows, row_count=len(rows))


def genomepy_head_annotations(
    name: str,
    *,
    provider: str | None = None,
    lines: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Show the first annotation lines of a genome (``genomepy.head_annotations``).

    genomepy downloads and converts the annotation with the UCSC tools, so those are required.
    The lines genomepy prints are captured and returned. Needs network access.
    """
    genome = _text(name, "name")
    count = _positive_or_none(lines, "lines")
    if count is None:
        raise ConfigurationError("lines must be a positive integer.")
    functions = _module("genomepy.functions")
    _require_ucsc_tools(_module("genomepy.annotation.utils"))
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        functions.head_annotations(genome, provider, n=count)
    printed = [line for line in buffer.getvalue().splitlines() if line.strip()]
    return _result("head_annotations", genome=genome, lines=printed)


def genomepy_generate_annotation(
    template: str,
    target: str,
    *,
    overwrite: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert GTF or BED annotation into the other format with UCSC tools (``generate_annot``)."""
    source = _existing_file(template, "template")
    output = Path(target).expanduser().resolve()
    if output.exists() and not overwrite:
        raise ConfigurationError(f"target already exists; pass overwrite=True: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    annotation = _module("genomepy.annotation.utils")
    _require_ucsc_tools(annotation)
    try:
        annotation.generate_annot(str(source), str(output), overwrite=bool(overwrite))
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    if not output.is_file():
        raise BackendExecutionError(f"generate_annot did not write {output}.")
    return _result("generate_annotation", {"annotation": output})


def genomepy_download_assembly_report(
    accession: str,
    *,
    output_path: str | None = None,
    quiet: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Retrieve an NCBI assembly report (chromosome names, roles, lengths) for an accession.

    With ``output_path`` the report is written as a new TSV and the path is returned; otherwise the
    rows are returned. ``accession`` must start with GCA or GCF. Needs network access.
    """
    acc = _text(accession, "accession")
    ncbi = _module("genomepy.providers.ncbi")
    target = None if output_path is None else _fresh_path(output_path, "output_path")
    if target is None:
        frame = ncbi.download_assembly_report(acc, quiet=bool(quiet))
        if frame is None:
            raise BackendExecutionError(f"No NCBI assembly report for {acc}.")
        rows = _frame_rows(frame)
        return _result("assembly_report", accession=acc, rows=rows, row_count=len(rows))
    ncbi.download_assembly_report(acc, fname=str(target), quiet=bool(quiet))
    if not target.is_file():
        raise BackendExecutionError(f"No NCBI assembly report for {acc}.")
    return _result("assembly_report", {"report": target}, accession=acc)


def genomepy_genome_size(accession: str, *, verbose: bool = True) -> dict[str, Any]:
    """Total sequence length of an NCBI assembly, summed from its assembly report.

    Used for the ncbi, gencode and ucsc providers (they share one implementation). The size is
    ``None`` when NCBI has no report for the accession. Needs network access.
    """
    acc = _text(accession, "accession")
    size = _module("genomepy.providers.ncbi").get_genome_size(acc)
    found = isinstance(size, int) and size >= 0
    return _result("genome_size", accession=acc, genome_size=int(size) if found else None)


def genomepy_ncbi_genomes(assembly_url: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse the NCBI assembly summary files into a genome catalogue (``ncbi.get_genomes``)."""
    url = _text(assembly_url, "assembly_url")
    genomes = _module("genomepy.providers.ncbi").get_genomes(url)
    return _result("ncbi_genomes", genomes=_plain(genomes), genome_count=len(genomes))


def genomepy_ensembl_genomes(rest_url: str, *, verbose: bool = True) -> dict[str, Any]:
    """Genome catalogue from the Ensembl REST API (``ensembl.get_genomes``). Needs network."""
    url = _text(rest_url, "rest_url")
    genomes = _module("genomepy.providers.ensembl").get_genomes(url)
    return _result("ensembl_genomes", genomes=_plain(genomes), genome_count=len(genomes))


def genomepy_ucsc_genomes(rest_url: str, *, verbose: bool = True) -> dict[str, Any]:
    """Genome catalogue from the UCSC REST API (``ucsc.get_genomes``). Needs network access."""
    url = _text(rest_url, "rest_url")
    genomes = _module("genomepy.providers.ucsc").get_genomes(url)
    return _result("ucsc_genomes", genomes=_plain(genomes), genome_count=len(genomes))


def genomepy_gencode_genomes(ftp_link: str, *, verbose: bool = True) -> dict[str, Any]:
    """Latest GENCODE release of each major assembly from the GENCODE FTP site. Needs network."""
    link = _text(ftp_link, "ftp_link")
    genomes = _module("genomepy.providers.gencode").get_genomes(link)
    return _result("gencode_genomes", genomes=_plain(genomes), genome_count=len(genomes))


def genomepy_gencode_releases(
    listing: Sequence[str],
    species: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """GENCODE release numbers from a listing of release directories, highest first.

    Pure function over the listing (no network). ``species`` is ``human`` or ``mouse``; mouse
    release numbers get the ``M`` prefix that GENCODE uses.
    """
    if isinstance(listing, (str, bytes)) or not all(isinstance(item, str) for item in listing):
        raise ConfigurationError("listing must be a list of directory names.")
    if species not in _SPECIES:
        raise ConfigurationError(f"species must be one of {list(_SPECIES)}.")
    releases = _module("genomepy.providers.gencode").get_releases(list(listing), species)
    return _result("gencode_releases", releases=list(releases), release_count=len(releases))


def genomepy_gencode_to_ucsc_names(
    genomes: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map GENCODE chromosome names to UCSC names for the given catalogue (``get_gencode2ucsc``)."""
    if not isinstance(genomes, Mapping) or not genomes:
        raise ConfigurationError("genomes must be a non-empty catalogue mapping.")
    mapping = _module("genomepy.providers.gencode").get_gencode2ucsc(dict(genomes))
    return _result("gencode_to_ucsc_names", mapping=_plain(mapping), entry_count=len(mapping))


def _annotation_result(operation: str, directory: Path) -> dict[str, Any]:
    files = sorted(artifacts(directory).values())
    return _result(operation, {"directory": directory, **_numbered(files, "file")})


def _numbered(paths: Sequence[Path], prefix: str) -> dict[str, Path]:
    return {f"{prefix}_{index}": path for index, path in enumerate(paths)}


def genomepy_download_annotation(
    genomes_dir: str,
    annotation_url: str,
    localname: str,
    *,
    n: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download an annotation file and convert it to the genePred and BED outputs of genomepy.

    Shared by the base and gencode providers. Requires the UCSC tools; needs network access.
    """
    base_dir = Path(genomes_dir).expanduser().resolve()
    name = _text(localname, "localname")
    url = _text(annotation_url, "annotation_url")
    annotation = _module("genomepy.annotation.utils")
    _require_ucsc_tools(annotation)
    providers = _module("genomepy.providers.base")
    providers.download_annotation(str(base_dir), url, name, n=_positive_or_none(n, "n"))
    directory = base_dir / name
    if not directory.is_dir():
        raise BackendExecutionError(f"genomepy did not create {directory}.")
    return _annotation_result("download_annotation", directory)


def genomepy_ucsc_download_annotation(
    name: str,
    annotation: str,
    genomes_dir: str,
    localname: str,
    *,
    n: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download a UCSC annotation (extended genePred from the UCSC MySQL server). Needs network."""
    genome = _text(name, "name")
    kind = _text(annotation, "annotation")
    base_dir = Path(genomes_dir).expanduser().resolve()
    local = _text(localname, "localname")
    _require_ucsc_tools(_module("genomepy.annotation.utils"))
    ucsc = _module("genomepy.providers.ucsc")
    ucsc.download_annotation(genome, kind, str(base_dir), local, n=_positive_or_none(n, "n"))
    directory = base_dir / local
    if not directory.is_dir():
        raise BackendExecutionError(f"genomepy did not create {directory}.")
    return _annotation_result("ucsc_download_annotation", directory)


def genomepy_download_annotation_head(
    annotation_url: str,
    annotation_file: str,
    *,
    n: int = 5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Download only the first ``n`` lines of a (gzipped) annotation file. Needs network access."""
    url = _text(annotation_url, "annotation_url")
    target = _fresh_path(annotation_file, "annotation_file")
    count = _positive_or_none(n, "n")
    if count is None:
        raise ConfigurationError("n must be a positive integer.")
    _module("genomepy.providers.base").download_head(url, str(target), n=count)
    if not target.is_file():
        raise BackendExecutionError(f"No annotation lines were written to {target}.")
    return _result("annotation_head", {"head": target}, lines=count)


def genomepy_find_annotation_urls(url: str, *, verbose: bool = True) -> dict[str, Any]:
    """Look for GTF or GFF3 annotation files next to a genome URL. Needs network access."""
    location = _text(url, "url")
    links = _module("genomepy.providers.url").search_url_for_annotations(location)
    return _result("find_annotation_urls", urls=[str(link) for link in links], url_count=len(links))


def _run_index(
    aligner: str,
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None,
) -> dict[str, Any]:
    if (
        isinstance(command, (str, bytes))
        or not command
        or not all(isinstance(item, str) for item in command)
    ):
        raise ConfigurationError("command must be a non-empty list of strings (no shell string).")
    expected = _INDEX_PROGRAMS[aligner]
    argv = list(command)
    if Path(argv[0]).name != expected:
        raise ConfigurationError(f"The {aligner} index command must start with {expected!r}.")
    directory = Path(working_directory).expanduser().resolve()
    if not directory.is_dir():
        raise ConfigurationError(f"working_directory does not exist: {directory}")
    program = executable_path or shutil.which(expected)
    if program is None:
        raise BackendUnavailableError(
            f"The {aligner} index program {expected!r} was not found on PATH.",
            hint=f"Install {aligner} (for example conda install -c bioconda {aligner}) and put it "
            "on PATH, or pass executable_path.",
        )
    before = set(artifacts(directory))
    log = run(program, argv[1:], directory, ExecutionConfig())
    created = sorted(path for key, path in artifacts(directory).items() if key not in before)
    if not created:
        raise BackendExecutionError(f"{expected} wrote no index files.", context={"log": log})
    return _result(
        f"{aligner}_index_command",
        {f"index_{i}": path for i, path in enumerate(created)},
        backend=aligner,
        executable=str(Path(program).resolve()),
        log_tail=log[-_LOG_TAIL:],
    )


def genomepy_bowtie2_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a bowtie2-build index command given as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "bowtie2", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_bwa_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a bwa index command given as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "bwa", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_gmap_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a gmap_build index command given as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "gmap", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_hisat2_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a hisat2-build index command given as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "hisat2", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_minimap2_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a minimap2 index command given as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "minimap2", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_star_index_command(
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a STAR genome-generation command as an argument list (no shell; genomepy plugin)."""
    return _run_index(
        "star", command, working_directory=working_directory, executable_path=executable_path
    )


def genomepy_run_index_cmd(
    name: str,
    command: Sequence[str],
    *,
    working_directory: str,
    executable_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run one named aligner's index command (``genomepy.utils.run_index_cmd``, no shell).

    ``name`` is one of bowtie2, bwa, gmap, hisat2, minimap2 or star; the command must start with
    that aligner's index program.
    """
    if name not in _INDEX_PROGRAMS:
        raise ConfigurationError(f"name must be one of {sorted(_INDEX_PROGRAMS)}.")
    return _run_index(
        name, command, working_directory=working_directory, executable_path=executable_path
    )


__all__ = [
    "genomepy_bowtie2_index_command",
    "genomepy_bwa_index_command",
    "genomepy_download_annotation",
    "genomepy_download_annotation_head",
    "genomepy_download_assembly_report",
    "genomepy_ensembl_genomes",
    "genomepy_find_annotation_urls",
    "genomepy_gencode_genomes",
    "genomepy_gencode_releases",
    "genomepy_gencode_to_ucsc_names",
    "genomepy_generate_annotation",
    "genomepy_genome_size",
    "genomepy_gmap_index_command",
    "genomepy_head_annotations",
    "genomepy_hisat2_index_command",
    "genomepy_list_available_genomes",
    "genomepy_minimap2_index_command",
    "genomepy_ncbi_genomes",
    "genomepy_run_index_cmd",
    "genomepy_search_genomes",
    "genomepy_star_index_command",
    "genomepy_ucsc_download_annotation",
    "genomepy_ucsc_genomes",
]

_register_public_verbose(globals())
