"""Optional service clients and a local Toil job graph for gene, sequence, and API lookups.

- ``compute_toil_sequence_stats()`` runs two chained Toil jobs over a batch of DNA
  sequences (reverse complement, then GC fraction) on Toil's single-machine batch system.
  The sequence operations are Biopython calls made inside the jobs.
- ``search_smartapi_registry()`` and ``describe_smartapi_api()`` use the
  ``smartapi-mcp`` registry query, ranking, and OpenAPI spec loader.
- ``get_biomcp_gene()``, ``search_biomcp_genes()``, and ``get_biomcp_variant()``
  call the ``biomcp`` command from ``biomcp-cli`` with a fixed argument list.
  biomcp has no Python API, so these functions use its JSON CLI output.
- ``list_galaxy_tools()`` and ``show_galaxy_tool()`` read a Galaxy tool catalogue
  through ``bioblend``, using the redirect-safe client from ``galaxy-mcp``.
- ``analyze_fasta_with_clawbio()`` runs the packaged ``analyze-fasta`` script
  from ``clawbio`` on one bounded nucleotide FASTA file.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import contextlib
import importlib.metadata as metadata
import importlib.util
import json
import math
import os
import re
import shutil
import sys
import tempfile
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import input_path, named_sequence_items, require, sequence_text
from dnakit.core import DNA
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    InputFormatError,
    QueryError,
)
from dnakit.tools.serialization import serialize_for_agent

_TOIL_EXTRA = "bio-toil"
_SMARTAPI_EXTRA = "bio-smartapi"
_BIOMCP_EXTRA = "bio-biomcp"
_GALAXY_EXTRA = "bio-galaxy"
_CLAWBIO_EXTRA = "bio-clawbio"
_BIOMCP_TIMEOUT_SECONDS = 60.0
_BIOMCP_MAX_OUTPUT_BYTES = 1_000_000
_HTTP_METHODS = ("get", "put", "post", "delete", "options", "head", "patch", "trace")
_GENE_SYMBOL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,31}")
_VARIANT_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.:>()+*=\[\],-]{0,127}")
_SMARTAPI_ID = re.compile(r"[0-9a-f]{32}")
_GALAXY_TOOL_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_./:-]{0,255}")
_API_KEY = re.compile(r"[!-~]{1,256}")


def _bounded_int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _matching(value: object, name: str, pattern: re.Pattern[str]) -> str:
    text = _text(value, name, maximum=256)
    if not pattern.fullmatch(text):
        raise ConfigurationError(f"{name} has an unsupported format: {text!r}.")
    return text


def _install_hint(extra: str) -> str:
    return f'Install the optional dependency with: python -m pip install "dnakit[{extra}]"'


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {metadata.version(distribution)}"
    except metadata.PackageNotFoundError:
        return f"{distribution} unknown"


def _number_or_none(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


@contextlib.contextmanager
def _interpreter_bin_on_path() -> Iterator[None]:
    """Expose the running interpreter's console scripts to Toil's worker command."""

    previous = os.environ.get("PATH")
    bin_dir = str(Path(sys.executable).absolute().parent)
    os.environ["PATH"] = bin_dir if not previous else os.pathsep.join([bin_dir, previous])
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop("PATH", None)
        else:
            os.environ["PATH"] = previous


def _toil_reverse_complement_all(job: Any, items: list[tuple[str, str]]) -> list[tuple[str, str]]:
    from Bio.Seq import reverse_complement

    return [(name, str(reverse_complement(sequence))) for name, sequence in items]


def _toil_gc_fraction_all(job: Any, reversed_items: list[tuple[str, str]]) -> list[dict[str, Any]]:
    from Bio.SeqUtils import gc_fraction

    return [
        {
            "id": name,
            "length": len(sequence),
            "reverse_complement": sequence,
            "gc_fraction": float(gc_fraction(sequence, ambiguous="ignore")),
        }
        for name, sequence in reversed_items
    ]


def _toil_sequence_root(job: Any, items: list[tuple[str, str]]) -> Any:
    """Chain reverse-complement and GC-fraction jobs over the whole batch."""

    reverse = job.addChildJobFn(_toil_reverse_complement_all, items)
    gc = reverse.addFollowOnJobFn(_toil_gc_fraction_all, reverse.rv())
    return gc.rv()


def _toil_items(
    sequences: object, max_sequences: int, max_total_bases: int
) -> list[tuple[str, str]]:
    count_limit = _bounded_int(max_sequences, "max_sequences", 1, 1000)
    base_limit = _bounded_int(max_total_bases, "max_total_bases", 1, 100_000_000)
    if not isinstance(sequences, (DNA, Mapping)):
        raise ConfigurationError("sequences must be DNA([...]) or an ID-to-sequence mapping.")
    pairs = named_sequence_items(sequences)
    if not pairs:
        raise ConfigurationError("At least one sequence is required.")
    if len(pairs) > count_limit:
        raise ConfigurationError(f"Too many sequences: {len(pairs)} > max_sequences={count_limit}.")
    items: list[tuple[str, str]] = []
    total = 0
    for name, value in pairs:
        identifier = _text(name, "sequence name", maximum=200)
        sequence = sequence_text(value)
        total += len(sequence)
        if total > base_limit:
            raise ConfigurationError(f"Total bases exceed max_total_bases={base_limit}.")
        items.append((identifier, sequence))
    return items


def compute_toil_sequence_stats(
    sequences: Mapping[str, str] | DNA,
    *,
    max_sequences: int = 100,
    max_total_bases: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run reverse complement and GC fraction as two chained Toil jobs over a batch."""

    items = _toil_items(sequences, max_sequences, max_total_bases)
    toil_job = require("toil.job", _TOIL_EXTRA)
    toil_common = require("toil.common", _TOIL_EXTRA)
    require("Bio.Seq", _TOIL_EXTRA)
    with _interpreter_bin_on_path(), tempfile.TemporaryDirectory(prefix="dnakit-toil-") as workdir:
        options = toil_job.Job.Runner.getDefaultOptions(str(Path(workdir) / "jobstore"))
        options.workDir = workdir
        options.logLevel = "ERROR"
        options.clean = "always"
        options.batchSystem = "single_machine"
        root = toil_job.Job.wrapJobFn(_toil_sequence_root, items)
        with toil_common.Toil(options) as toil:
            rows = toil.start(root)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "compute_toil_sequence_stats",
            "sequence_count": len(rows),
            "sequences": rows,
            "provenance": {
                "package": _package("toil"),
                "dependency": _package("biopython"),
                "batch_system": "single_machine",
                "job_chain": "root -> reverse_complement (child) -> gc_fraction (follow-on)",
                "gc_convention": "G+C divided by all characters, including ambiguous codes",
                "input_convention": "sequences are uppercased before processing",
            },
        },
        verbose=verbose,
    )


def _run_coroutine(coroutine: Any) -> Any:
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coroutine)
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        return pool.submit(asyncio.run, coroutine).result()


def search_smartapi_registry(
    query: str, *, limit: int = 10, verbose: bool = True
) -> dict[str, Any]:
    """Search the SmartAPI registry by keyword and rank the matching APIs."""

    text = _text(query, "query", maximum=200)
    count = _bounded_int(limit, "limit", 1, 100)
    smartapi = require("smartapi_mcp.smartapi", _SMARTAPI_EXTRA)
    biothings = require("smartapi_mcp.biothings", _SMARTAPI_EXTRA)
    # The registry's relevance depends on case ("MyGene" and "mygene" return different hits);
    # rank_apis is case-insensitive, so the lowercase form is sent.
    registry_query = text.lower()
    try:
        entries = _run_coroutine(smartapi.get_smartapi_registry(q=registry_query))
    except Exception as exc:
        raise QueryError(
            "SmartAPI registry query failed.",
            context={"error_type": type(exc).__name__},
        ) from exc
    registry = biothings.build_registry_from_entries(entries)
    ranked = biothings.rank_apis(registry, keyword=text, limit=count)
    matches = [{"smartapi_id": registry[item["name"]].smartapi_id, **item} for item in ranked]
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "search_smartapi_registry",
            "query": text,
            "registry_query": registry_query,
            "registry_hits": len(entries),
            "matches": matches,
            "provenance": {
                "package": _package("smartapi-mcp"),
                "source": "https://smart-api.info/api/query",
                "convention": "ranking by smartapi_mcp rank_apis; zero-score APIs are omitted",
            },
        },
        verbose=verbose,
    )


def describe_smartapi_api(
    smartapi_id: str, *, max_operations: int = 200, verbose: bool = True
) -> dict[str, Any]:
    """Load one SmartAPI OpenAPI spec and list its HTTP operations."""

    identifier = _matching(smartapi_id, "smartapi_id", _SMARTAPI_ID)
    limit = _bounded_int(max_operations, "max_operations", 1, 1000)
    smartapi = require("smartapi_mcp.smartapi", _SMARTAPI_EXTRA)
    try:
        spec = smartapi.load_api_spec(identifier)
        server_url = smartapi.get_base_server_url(spec)
    except Exception as exc:
        raise QueryError(
            "SmartAPI OpenAPI spec could not be loaded.",
            context={"smartapi_id": identifier, "error_type": type(exc).__name__},
        ) from exc
    operations: list[dict[str, Any]] = []
    paths = spec.get("paths", {})
    for path in sorted(paths):
        item = paths[path]
        if not isinstance(item, Mapping):
            continue
        for method in _HTTP_METHODS:
            operation = item.get(method)
            if isinstance(operation, Mapping):
                operations.append(
                    {
                        "method": method.upper(),
                        "path": path,
                        "operation_id": operation.get("operationId"),
                        "summary": operation.get("summary"),
                    }
                )
    info = spec.get("info", {})
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "describe_smartapi_api",
            "smartapi_id": identifier,
            "title": info.get("title"),
            "version": info.get("version"),
            "openapi": spec.get("openapi"),
            "server_url": server_url,
            "operation_count": len(operations),
            "truncated": len(operations) > limit,
            "operations": operations[:limit],
            "provenance": {
                "package": _package("smartapi-mcp"),
                "source": "https://smart-api.info/api/metadata/{smartapi_id}",
                "convention": "sorted by path, then by method in the order GET, PUT, POST, DELETE",
            },
        },
        verbose=verbose,
    )


def _biomcp_executable() -> str | None:
    found = shutil.which("biomcp")
    if found:
        return str(Path(found).absolute())
    candidate = Path(sys.executable).absolute().parent / "biomcp"
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return None


def _biomcp_json(arguments: list[str], *, verbose: bool) -> dict[str, Any]:
    executable = _biomcp_executable()
    if executable is None:
        raise BackendUnavailableError(
            "The biomcp executable was not found.",
            code="BACKEND_UNAVAILABLE",
            hint=_install_hint(_BIOMCP_EXTRA),
        )
    with tempfile.TemporaryDirectory(prefix="dnakit-biomcp-") as directory:
        result = execute_bounded_command(
            executable,
            ["-j", *arguments],
            backend_id="biomcp",
            cwd=directory,
            timeout_seconds=_BIOMCP_TIMEOUT_SECONDS,
            max_output_bytes=_BIOMCP_MAX_OUTPUT_BYTES,
            verbose=verbose,
        )
    try:
        payload = json.loads(result.output)
    except json.JSONDecodeError as exc:
        raise BackendExecutionError(
            "biomcp did not return JSON.",
            context={"return_code": result.return_code},
        ) from exc
    if not isinstance(payload, dict):
        raise BackendExecutionError("biomcp returned an unexpected JSON value.")
    error = payload.get("error")
    if isinstance(error, Mapping):
        lines = [line.strip() for line in str(error.get("message", "")).splitlines()]
        message = next((line for line in lines if line), "biomcp query failed.")
        code = error.get("code")
        query_code = f"BIOMCP_{str(code).upper()}" if code else None
        raise QueryError(message, code=query_code, context={"biomcp_code": code})
    if result.return_code != 0:
        raise BackendExecutionError(
            "biomcp failed without a JSON error.",
            context={"return_code": result.return_code},
        )
    return payload


def get_biomcp_gene(symbol: str, *, verbose: bool = True) -> dict[str, Any]:
    """Look up one human gene symbol through the biomcp CLI."""

    name = _matching(symbol, "symbol", _GENE_SYMBOL)
    payload = _biomcp_json(["get", "gene", name], verbose=verbose)
    coordinates = payload.get("genomic_coordinates") or {}
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "get_biomcp_gene",
            "symbol": payload.get("symbol"),
            "name": payload.get("name"),
            "entrez_id": payload.get("entrez_id"),
            "ensembl_id": payload.get("ensembl_id"),
            "uniprot_id": payload.get("uniprot_id"),
            "omim_id": payload.get("omim_id"),
            "gene_type": payload.get("gene_type"),
            "location": payload.get("location"),
            "coordinate": coordinates.get("coordinate"),
            "genome_build": coordinates.get("genome_build"),
            "aliases": list(payload.get("aliases") or []),
            "clinical_diseases": list(payload.get("clinical_diseases") or []),
            "summary": payload.get("summary"),
            "provenance": {
                "package": _package("biomcp-cli"),
                "source": "biomcp get gene (MyGene.info and NCBI Gene)",
                "convention": "coordinates and summary are passed through, not recomputed",
            },
        },
        verbose=verbose,
    )


def search_biomcp_genes(query: str, *, limit: int = 10, verbose: bool = True) -> dict[str, Any]:
    """Search gene records by free text through the biomcp CLI."""

    text = _text(query, "query", maximum=200)
    if text.startswith("-"):
        raise ConfigurationError("query cannot start with '-'.")
    count = _bounded_int(limit, "limit", 1, 50)
    payload = _biomcp_json(["search", "gene", "-q", text, "--limit", str(count)], verbose=verbose)
    results = payload.get("results") or []
    rows = [
        {
            "symbol": item.get("symbol"),
            "name": item.get("name"),
            "entrez_id": item.get("entrez_id"),
            "uniprot_id": item.get("uniprot_id"),
            "omim_id": item.get("omim_id"),
            "coordinate": (item.get("genomic_coordinates") or {}).get("coordinate"),
        }
        for item in results
        if isinstance(item, Mapping)
    ]
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "search_biomcp_genes",
            "query": text,
            "total": (payload.get("pagination") or {}).get("total", len(rows)),
            "genes": rows,
            "provenance": {
                "package": _package("biomcp-cli"),
                "source": "biomcp search gene (MyGene.info and NCBI Gene)",
                "convention": "results are returned in the order given by biomcp",
            },
        },
        verbose=verbose,
    )


def get_biomcp_variant(variant_id: str, *, verbose: bool = True) -> dict[str, Any]:
    """Look up one variant by rsID or genomic HGVS through the biomcp CLI."""

    identifier = _matching(variant_id, "variant_id", _VARIANT_ID)
    payload = _biomcp_json(["get", "variant", identifier], verbose=verbose)
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "get_biomcp_variant",
            "id": payload.get("id"),
            "rsid": payload.get("rsid"),
            "gene": payload.get("gene"),
            "hgvs_c": payload.get("hgvs_c"),
            "hgvs_p": payload.get("hgvs_p"),
            "consequence": payload.get("consequence"),
            "significance": payload.get("significance"),
            "transcript": payload.get("transcript"),
            "genome_build": payload.get("genome_build"),
            "cadd_score": _number_or_none(payload.get("cadd_score")),
            "sift_pred": payload.get("sift_pred"),
            "polyphen_pred": payload.get("polyphen_pred"),
            "provenance": {
                "package": _package("biomcp-cli"),
                "source": "biomcp get variant (MyVariant.info and related annotation sources)",
                "convention": "annotations are passed through from biomcp, not recomputed",
            },
        },
        verbose=verbose,
    )


def _galaxy_client(url: str) -> tuple[Any, str]:
    client = require("galaxy_mcp.client", _GALAXY_EXTRA)
    text = _text(url, "url", maximum=2048)
    try:
        normalized = client.normalize_galaxy_url(text)
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc
    return client, normalized


def _galaxy_key(api_key: object) -> str | None:
    if api_key is None:
        return None
    return _matching(api_key, "api_key", _API_KEY)


def _galaxy_tool_row(item: Mapping[str, Any]) -> dict[str, str]:
    return {
        "id": str(item.get("id", "")),
        "name": str(item.get("name", "")),
        "version": str(item.get("version", "")),
        "description": str(item.get("description", "")),
    }


def list_galaxy_tools(
    url: str,
    *,
    api_key: str | None = None,
    limit: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """List tools from a Galaxy server; the key is optional for public servers."""

    count = _bounded_int(limit, "limit", 1, 500)
    client, server = _galaxy_client(url)
    key = _galaxy_key(api_key)
    gi = client.GalaxyInstance(url=server, key=key)
    try:
        tools = gi.tools.get_tools()
    except Exception as exc:
        raise QueryError(
            "Galaxy did not return a tool list.",
            context={"url": server, "error_type": type(exc).__name__},
        ) from exc
    rows = sorted(
        (_galaxy_tool_row(item) for item in tools if isinstance(item, Mapping)),
        key=lambda row: row["id"],
    )
    selected = rows[:count]
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "list_galaxy_tools",
            "url": server,
            "total_tools": len(rows),
            "returned": len(selected),
            "tools": selected,
            "provenance": {
                "package": _package("galaxy-mcp"),
                "dependency": _package("bioblend"),
                "convention": "sorted by id; the full list is downloaded before limit applies",
                "credentials": "api_key is sent in the x-api-key header and is never returned",
            },
        },
        verbose=verbose,
    )


def show_galaxy_tool(
    url: str,
    *,
    tool_id: str,
    api_key: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Show one Galaxy tool's metadata; the key is optional for public servers."""

    identifier = _matching(tool_id, "tool_id", _GALAXY_TOOL_ID)
    if ".." in identifier:
        raise ConfigurationError("tool_id cannot contain '..'.")
    client, server = _galaxy_client(url)
    key = _galaxy_key(api_key)
    gi = client.GalaxyInstance(url=server, key=key)
    try:
        item = gi.tools.show_tool(identifier)
    except Exception as exc:
        raise QueryError(
            "Galaxy did not return the requested tool.",
            context={"url": server, "tool_id": identifier, "error_type": type(exc).__name__},
        ) from exc
    if not isinstance(item, Mapping):
        raise BackendExecutionError("Galaxy returned an unexpected tool record.")
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "show_galaxy_tool",
            "url": server,
            "tool": {**_galaxy_tool_row(item), "model_class": str(item.get("model_class", ""))},
            "provenance": {
                "package": _package("galaxy-mcp"),
                "dependency": _package("bioblend"),
                "credentials": "api_key is sent in the x-api-key header and is never returned",
            },
        },
        verbose=verbose,
    )


def _load_clawbio_analyze_fasta() -> Any:
    clawbio = require("clawbio", _CLAWBIO_EXTRA)
    script = (
        Path(clawbio.__file__).resolve().parent / "skills" / "analyze-fasta" / "analyze_fasta.py"
    )
    if not script.is_file():
        raise BackendUnavailableError(
            "The installed clawbio package does not provide the analyze-fasta skill.",
            hint=_install_hint(_CLAWBIO_EXTRA),
        )
    spec = importlib.util.spec_from_file_location("_dnakit_clawbio_analyze_fasta", script)
    if spec is None or spec.loader is None:
        raise BackendUnavailableError("clawbio analyze-fasta script could not be loaded.")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        raise BackendUnavailableError(
            f"clawbio analyze-fasta needs an optional dependency: {exc}",
            hint=_install_hint(_CLAWBIO_EXTRA),
        ) from exc
    return module


def analyze_fasta_with_clawbio(
    path: str | Path,
    *,
    max_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize one nucleotide FASTA file with the clawbio analyze-fasta skill."""

    limit = _bounded_int(max_bytes, "max_bytes", 1, 100_000_000)
    source = input_path(path)
    size = source.stat().st_size
    if size > limit:
        raise ConfigurationError(f"FASTA file has {size} bytes; max_bytes is {limit}.")
    module = _load_clawbio_analyze_fasta()
    try:
        report = module.analyze_fasta(str(source))
    except (OSError, ValueError) as exc:
        raise InputFormatError(f"clawbio could not parse the FASTA file: {exc}") from exc
    if not isinstance(report, Mapping):
        raise BackendExecutionError("clawbio returned an unexpected report.")
    if "error" in report:
        raise InputFormatError(f"clawbio rejected the FASTA file: {report['error']}")
    if report.get("sequence_type") != "nucleotide":
        raise InputFormatError(
            "clawbio classified the first record as protein; nucleotide FASTA is required."
        )
    sequences = [
        {
            "index": index,
            "id": str(item.get("id", "")),
            "length_bp": item.get("length_bp"),
            "gc_percent": _number_or_none(item.get("gc_content")),
            "composition": {
                str(key): int(value) for key, value in (item.get("composition") or {}).items()
            },
        }
        for index, item in enumerate(report.get("sequences", []))
        if isinstance(item, Mapping)
    ]
    summary = report.get("summary") or {}
    return serialize_for_agent(
        {
            "status": "completed",
            "operation": "analyze_fasta_with_clawbio",
            "file": source.name,
            "sequence_type": "nucleotide",
            "sequences": sequences,
            "summary": {
                "total_sequences": summary.get("total_sequences"),
                "total_residues": summary.get("total_residues"),
                "min_length": summary.get("min_length"),
                "max_length": summary.get("max_length"),
                "n50": summary.get("n50"),
                "avg_gc_percent": _number_or_none(summary.get("avg_gc_content")),
            },
            "provenance": {
                "package": _package("clawbio"),
                "dependency": _package("biopython"),
                "source": "clawbio/skills/analyze-fasta/analyze_fasta.py",
                "convention": "GC is a percentage; sequence type is decided from the first record",
            },
        },
        verbose=verbose,
    )


__all__ = [
    "analyze_fasta_with_clawbio",
    "compute_toil_sequence_stats",
    "describe_smartapi_api",
    "get_biomcp_gene",
    "get_biomcp_variant",
    "list_galaxy_tools",
    "search_biomcp_genes",
    "search_smartapi_registry",
    "show_galaxy_tool",
]

_register_public_verbose(globals(), tuple(__all__))
