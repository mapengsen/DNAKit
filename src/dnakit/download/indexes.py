"""Explicit, bounded adapters for common nucleotide index builders."""

from __future__ import annotations

import hashlib
import os
import tempfile
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.search._shared import adapter_provenance

from .files import _install_staged_files, _write_manifest, resolved_config
from .models import DownloadConfig, IndexArtifact, IndexBuildResult

_TOOLS = frozenset(
    {"makeblastdb", "bwa", "bowtie2-build", "minimap2", "samtools-fai", "samtools-dict"}
)


@dataclass(frozen=True, slots=True)
class IndexProgress:
    """Start/finish progress for one external index build."""

    tool: str
    stage: str


def _known_outputs(tool: str, prefix: Path) -> tuple[Path, ...]:
    if tool == "makeblastdb":
        return tuple(
            Path(f"{prefix}{suffix}")
            for suffix in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto")
        )
    if tool == "bwa":
        return tuple(
            Path(f"{prefix}{suffix}") for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")
        )
    if tool == "bowtie2-build":
        return tuple(
            Path(f"{prefix}{suffix}")
            for suffix in (
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
                ".1.bt2l",
                ".2.bt2l",
                ".3.bt2l",
                ".4.bt2l",
                ".rev.1.bt2l",
                ".rev.2.bt2l",
            )
        )
    if tool == "minimap2":
        return (Path(f"{prefix}.mmi"),)
    if tool == "samtools-fai":
        return (Path(f"{prefix}.fai"),)
    return (Path(f"{prefix}.dict"),)


def _arguments(tool: str, fasta: Path, prefix: Path, threads: int) -> tuple[str, ...]:
    if tool == "makeblastdb":
        return ("-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix), "-parse_seqids")
    if tool == "bwa":
        return ("index", "-p", str(prefix), str(fasta))
    if tool == "bowtie2-build":
        return ("--threads", str(threads), str(fasta), str(prefix))
    if tool == "minimap2":
        return ("-d", f"{prefix}.mmi", "-t", str(threads), str(fasta))
    if tool == "samtools-fai":
        return ("faidx", "-o", f"{prefix}.fai", str(fasta))
    return ("dict", "-o", f"{prefix}.dict", str(fasta))


def _actual_outputs(tool: str, known: Iterable[Path]) -> tuple[Path, ...]:
    existing = tuple(path for path in known if path.is_file())
    if tool == "bowtie2-build":
        standard = tuple(path for path in existing if path.suffix == ".bt2")
        large = tuple(path for path in existing if path.suffix == ".bt2l")
        return standard or large
    return existing


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1_048_576):
            digest.update(chunk)
    return digest.hexdigest()


def build_index(
    fasta_path: str | os.PathLike[str],
    *,
    tool: str,
    executable_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    prefix: str | None = None,
    threads: int = 1,
    timeout_seconds: float = 3_600.0,
    config: DownloadConfig | None = None,
    progress: Callable[[IndexProgress], None] | None = None,
) -> IndexBuildResult:
    """Build BLAST/BWA/Bowtie2/minimap2/samtools indexes using an explicit executable."""

    resolved = resolved_config(config)
    tool_value = tool.strip() if isinstance(tool, str) else ""
    if tool_value not in _TOOLS:
        raise ConfigurationError(
            "Unsupported index tool.",
            code="INVALID_INDEX_TOOL",
            context={"allowed": sorted(_TOOLS)},
        )
    try:
        fasta = Path(fasta_path).expanduser().resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError("fasta_path must be an existing file.") from exc
    if not fasta.is_file():
        raise ConfigurationError("fasta_path must be an existing file.")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 256:
        raise ConfigurationError("threads must be in [1, 256].")
    if progress is not None and not callable(progress):
        raise ConfigurationError("progress must be callable or None.")
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and not output.is_dir():
        raise NotADirectoryError(f"output_dir must be a directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    prefix_value = fasta.stem if prefix is None else prefix
    if (
        not isinstance(prefix_value, str)
        or not prefix_value
        or prefix_value in {".", ".."}
        or "/" in prefix_value
        or "\\" in prefix_value
        or "\x00" in prefix_value
    ):
        raise ConfigurationError("prefix must be one safe filename component.")
    final_prefix = output / prefix_value
    manifest = output / f"{prefix_value}.{tool_value}.manifest.json"
    possible_targets = (*_known_outputs(tool_value, final_prefix), manifest)
    if any(path.exists() and not path.is_file() for path in possible_targets):
        raise IsADirectoryError("Every existing index output must be a file.")
    if not resolved.overwrite and any(path.exists() for path in possible_targets):
        raise FileExistsError("Refusing to overwrite existing index outputs.")
    if progress is not None:
        progress(IndexProgress(tool_value, "start"))
    with tempfile.TemporaryDirectory(prefix=".dnakit-index-", dir=str(output)) as work_raw:
        work = Path(work_raw)
        staged_prefix = work / prefix_value
        known = _known_outputs(tool_value, staged_prefix)
        command = execute_bounded_command(
            Path(executable_path),
            _arguments(tool_value, fasta, staged_prefix, threads),
            backend_id=tool_value,
            cwd=work,
            timeout_seconds=timeout_seconds,
            max_output_bytes=10_000_000,
            monitored_output_paths=known,
            max_monitored_output_bytes=min(resolved.max_total_bytes, 10_000_000_000),
        )
        if command.return_code != 0:
            raise BackendExecutionError(
                "Index builder returned a non-zero status.",
                code="INDEX_BUILD_FAILED",
                context={"tool": tool_value, "return_code": command.return_code},
            )
        staged = _actual_outputs(tool_value, known)
        minimum = 6 if tool_value == "bowtie2-build" else 1
        if len(staged) < minimum:
            raise BackendExecutionError(
                "Index builder did not create its expected artifacts.",
                code="INDEX_OUTPUT_MISSING",
                context={"tool": tool_value},
            )
        targets = tuple(
            output / path.name.replace(prefix_value, final_prefix.name, 1) for path in staged
        )
        artifacts = tuple(
            IndexArtifact(
                str(target),
                source.stat().st_size,
                _sha256_path(source),
            )
            for source, target in zip(staged, targets, strict=True)
        )
        provenance = adapter_provenance(
            tool_value,
            citation_url="https://www.ncbi.nlm.nih.gov/books/NBK569856/",
            filters={
                "tool": tool_value,
                "fasta_sha256": _sha256_path(fasta),
                "threads": threads,
            },
        )
        payload = {
            "tool": tool_value,
            "fasta_path": str(fasta),
            "output_prefix": str(final_prefix),
            "artifacts": [
                {"path": item.path, "byte_size": item.byte_size, "sha256": item.sha256}
                for item in artifacts
            ],
            "elapsed_seconds": command.elapsed_seconds,
            "provenance": provenance.to_dict(),
        }
        staged_manifest = work / manifest.name
        _write_manifest(payload, staged_manifest, overwrite=False)
        _install_staged_files(
            (*staged, staged_manifest),
            (*targets, manifest),
            overwrite=resolved.overwrite,
        )
    if progress is not None:
        progress(IndexProgress(tool_value, "complete"))
    return IndexBuildResult(
        tool_value,
        str(fasta),
        str(final_prefix),
        artifacts,
        command.output,
        command.elapsed_seconds,
        str(manifest),
        provenance,
    )


__all__ = ["IndexProgress", "build_index"]
