"""Minimal DNAKit command-line application."""

from __future__ import annotations

import json
import os
import platform
import re
import sys
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Annotated, NoReturn, Protocol

import typer
from rich.console import Console
from rich.progress import BarColumn, DownloadColumn, SpinnerColumn, TaskProgressColumn, TextColumn
from rich.table import Table

from dnakit._progress import Progress
from dnakit._version import __version__
from dnakit.cli.database import app as database_app
from dnakit.cli.gene_annotation import annotate_genes_command
from dnakit.cli.skills import app as skills_app
from dnakit.cli.visualization import plot_distribution_command
from dnakit.core import DNAAlphabet, DNARecord, DNASequence
from dnakit.datasets import DeduplicationConfig, SplitConfig, deduplicate, split
from dnakit.descriptors import (
    all_descriptors,
    base_composition,
    exact_repeat_fraction,
    gc_at_content,
    linguistic_complexity,
)
from dnakit.exceptions import DNAKitError
from dnakit.fingerprints import fracminhash, kmer_fingerprint, minhash
from dnakit.io import ReadConfig, WriteConfig, read, read_set, write
from dnakit.patterns import scan_motif, scan_orfs
from dnakit.references import DownloadProgress, GenomeDownloadResult, download_genome
from dnakit.similarity import compare
from dnakit.standardize import (
    AmbiguityPolicy,
    NormalizationConfig,
    UPolicy,
    ValidationConfig,
)
from dnakit.standardize import normalize as normalize_sequence
from dnakit.standardize import validate as validate_sequence
from dnakit.visualization import SaveConfig, build_html_report, save_html_report

if TYPE_CHECKING:
    from dnakit.predictions import DataMapping, LLMSchemaConfig

app = typer.Typer(
    help="Deterministic tools for DNA sequence analysis.",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
)
console = Console()
error_console = Console(stderr=True)

app.command("annotate-genes")(annotate_genes_command)
app.command("plot-distribution")(plot_distribution_command)
app.add_typer(skills_app, name="skills")
app.add_typer(database_app, name="database")


class _JSONResult(Protocol):
    def to_dict(self) -> dict[str, object]: ...


def _version_callback(value: bool) -> None:
    """Print the package version before command dispatch."""
    if value:
        typer.echo(__version__)
        raise typer.Exit


@app.callback()
def main(
    ctx: typer.Context,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            callback=_version_callback,
            help="Show the DNAKit version and exit.",
            is_eager=True,
        ),
    ] = False,
    workers: Annotated[
        int,
        typer.Option(
            "--workers", help="CPU workers: -1 for all available CPUs; 1 for serial processing."
        ),
    ] = -1,
) -> None:
    """Run DNAKit commands."""
    from dnakit._parallel import worker_scope

    try:
        scope = worker_scope(workers)
        scope.__enter__()
        ctx.call_on_close(lambda: scope.__exit__(None, None, None))
    except (DNAKitError, ValueError, TypeError) as exc:
        _abort(exc)


@app.command("info")
def info() -> None:
    """Show the local DNAKit runtime information."""
    table = Table(title="DNAKit runtime", show_header=False)
    table.add_column("Field", style="bold")
    table.add_column("Value")
    table.add_row("DNAKit", __version__)
    table.add_row("Python", platform.python_version())
    table.add_row("Executable", sys.executable)
    table.add_row("Platform", platform.platform())
    console.print(table)


@app.command("backends")
def backends() -> None:
    """List registered optional backends without probing external programs."""
    from dnakit.backends import backend_registry

    table = Table(title="DNAKit backends")
    table.add_column("Backend")
    table.add_column("Status")
    registrations = backend_registry.registrations()
    if not registrations:
        table.add_row("-", "No optional backends are registered in this build.")
    for registration in registrations:
        table.add_row(registration.backend_id, f"registered ({registration.source}); not probed")
    console.print(table)


def _download_genome_with_progress(
    query: str,
    output_dir: Path,
    *,
    overwrite: bool,
    keep_package: bool,
    api_key: str | None,
    enabled: bool,
) -> GenomeDownloadResult:
    if not enabled:
        return download_genome(
            query,
            output_dir,
            overwrite=overwrite,
            keep_package=keep_package,
            api_key=api_key,
            show_progress=False,
        )
    with Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        DownloadColumn(),
        console=error_console,
        transient=True,
    ) as display:
        task = display.add_task("Resolving assembly", total=None)
        previous_phase: tuple[str, str | None] | None = None

        def on_progress(event: DownloadProgress) -> None:
            nonlocal previous_phase
            if event.phase == "resolve":
                description = f"Resolving {event.query}"
            elif event.phase == "download":
                description = f"Downloading {event.accession}"
            else:
                description = f"Extracting {event.current_file or 'FASTA'}"
            phase = (event.phase, event.current_file)
            if phase != previous_phase:
                display.reset(task, description=description, total=event.total_bytes)
                previous_phase = phase
            display.update(
                task,
                description=description,
                total=event.total_bytes,
                completed=event.bytes_completed,
            )

        return download_genome(
            query,
            output_dir,
            overwrite=overwrite,
            keep_package=keep_package,
            api_key=api_key,
            progress=on_progress,
        )


@app.command("download-genome")
def download_genome_command(
    query: Annotated[
        str,
        typer.Argument(help="NCBI assembly accession, fixed alias, or species name."),
    ],
    output_dir: Annotated[
        Path,
        typer.Argument(help="Directory for the FASTA and NCBI provenance files."),
    ],
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
    keep_package: Annotated[
        bool,
        typer.Option("--keep-package", help="Keep the downloaded NCBI ZIP package."),
    ] = False,
    api_key: Annotated[
        str | None,
        typer.Option("--api-key", help="Optional NCBI Datasets API key."),
    ] = None,
    progress: Annotated[
        bool,
        typer.Option("--progress/--no-progress", help="Show download progress on stderr."),
    ] = True,
) -> None:
    """Download one NCBI genomic FASTA and verify its checksum."""

    try:
        result = _download_genome_with_progress(
            query,
            output_dir,
            overwrite=overwrite,
            keep_package=keep_package,
            api_key=api_key,
            enabled=progress,
        )
        _json_output(result.to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


def _json_output(value: object) -> None:
    typer.echo(json.dumps(value, ensure_ascii=False, sort_keys=True))


def _abort(error: BaseException) -> NoReturn:
    error_console.print(f"[red]{error}[/red]")
    raise typer.Exit(code=2)


def _tracked(records: Iterable[DNARecord], *, enabled: bool) -> Iterator[DNARecord]:
    if not enabled:
        yield from records
        return
    with Progress(
        SpinnerColumn(),
        TextColumn("Processed {task.completed} records"),
        console=error_console,
        transient=True,
    ) as progress:
        task = progress.add_task("records", total=None)
        for record in records:
            yield record
            progress.advance(task)


def _normalization_config(
    *,
    alphabet: str,
    keep_ambiguous: bool,
    keep_u: bool,
    keep_other: bool,
    u_policy: str | None,
    ambiguity_policy: str | None,
) -> NormalizationConfig:
    try:
        resolved_alphabet = None if alphabet == "auto" else DNAAlphabet(alphabet)
        return NormalizationConfig(
            alphabet=resolved_alphabet,
            keep_ambiguous=keep_ambiguous,
            keep_u=keep_u,
            keep_other=keep_other,
            u_policy=None if u_policy is None else UPolicy(u_policy),
            ambiguity_policy=(
                None if ambiguity_policy is None else AmbiguityPolicy(ambiguity_policy)
            ),
        )
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("normalize")
def normalize_command(
    sequence: Annotated[str, typer.Argument(help="Raw DNA sequence text.")],
    alphabet: Annotated[
        str, typer.Option(help="Output alphabet: auto, strict, or iupac.")
    ] = "auto",
    keep_ambiguous: Annotated[
        bool,
        typer.Option(
            "--keep-ambiguous/--delete-ambiguous",
            help="Keep IUPAC ambiguity symbols such as N, R, and Y.",
        ),
    ] = True,
    keep_u: Annotated[
        bool,
        typer.Option("--keep-u/--delete-u", help="Keep uracil U instead of deleting it."),
    ] = False,
    keep_other: Annotated[
        bool,
        typer.Option(
            "--keep-other/--delete-other",
            help="Keep characters outside the DNA/IUPAC alphabet.",
        ),
    ] = False,
    u_policy: Annotated[
        str | None,
        typer.Option(
            "--u-policy",
            help="Advanced U override: delete, replace, error, warn, or keep.",
        ),
    ] = None,
    ambiguity_policy: Annotated[
        str | None,
        typer.Option(
            "--ambiguity-policy",
            help="Advanced IUPAC override: delete, ignore, error, mask, or probability.",
        ),
    ] = None,
    include_raw: Annotated[
        bool,
        typer.Option(
            "--include-raw/--omit-raw",
            help="Include sensitive raw sequence text in JSON output.",
        ),
    ] = False,
) -> None:
    """Normalize one raw sequence and print an auditable JSON result."""

    try:
        result = normalize_sequence(
            sequence,
            config=_normalization_config(
                alphabet=alphabet,
                keep_ambiguous=keep_ambiguous,
                keep_u=keep_u,
                keep_other=keep_other,
                u_policy=u_policy,
                ambiguity_policy=ambiguity_policy,
            ),
        )
        _json_output(result.to_dict(include_raw_content=include_raw))
        if not result.is_valid:
            raise typer.Exit(code=2)
    except typer.Exit:
        raise
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("validate")
def validate_command(
    sequence: Annotated[str, typer.Argument(help="Raw DNA sequence text.")],
    allow_empty: Annotated[
        bool, typer.Option("--allow-empty/--reject-empty", help="Allow an empty sequence.")
    ] = False,
    sequence_length: Annotated[
        int | None,
        typer.Option(
            "--sequence-length",
            help="Require this exact normalized nucleotide symbol length.",
        ),
    ] = None,
) -> None:
    """Normalize and validate one raw sequence, returning JSON."""

    try:
        normalized = normalize_sequence(sequence)
        if normalized.sequence is None:
            _json_output(normalized.to_dict())
            raise typer.Exit(code=2)
        report = validate_sequence(
            normalized.sequence,
            config=ValidationConfig(
                allow_empty=allow_empty,
                min_length=0 if allow_empty else 1,
                sequence_length=sequence_length,
            ),
        )
        _json_output(report.to_dict())
        if not report.is_valid:
            raise typer.Exit(code=2)
    except typer.Exit:
        raise
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


def _normalized_sequence(raw: str) -> DNASequence:
    result = normalize_sequence(raw)
    if result.sequence is None:
        raise ValueError("Sequence normalization did not produce a valid sequence.")
    return result.sequence


@app.command("describe")
def describe_command(
    sequence: Annotated[str, typer.Argument(help="Raw DNA sequence text.")],
    all_fields: Annotated[
        bool,
        typer.Option(
            "--all/--compact",
            help="Return the fixed 240-field schema or the legacy compact report.",
        ),
    ] = True,
) -> None:
    """Calculate the complete or compact descriptor report for one sequence."""

    try:
        value = _normalized_sequence(sequence)
        if all_fields:
            _json_output(all_descriptors(value).to_dict())
            return
        _json_output(
            {
                "base_composition": base_composition(value, ambiguity_policy="ignore").to_dict(),
                "gc_at": gc_at_content(value, ambiguity_policy="ignore").to_dict(),
                "complexity": linguistic_complexity(value, ambiguity_policy="ignore").to_dict(),
                "repeat": exact_repeat_fraction(value, ambiguity_policy="ignore").to_dict(),
            }
        )
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("fingerprint")
def fingerprint_command(
    sequence: Annotated[str, typer.Argument(help="Raw DNA sequence text.")],
    kind: Annotated[
        str, typer.Option(help="Fingerprint type: kmer, minhash, or fracminhash.")
    ] = "kmer",
    k: Annotated[int, typer.Option(help="k-mer word length.")] = 3,
    size: Annotated[int, typer.Option(help="MinHash count or FracMinHash scaled value.")] = 1_000,
) -> None:
    """Build one native exact or sketch fingerprint."""

    try:
        value = _normalized_sequence(sequence)
        if kind == "kmer":
            result: _JSONResult = kmer_fingerprint(value, k=k)
        elif kind == "minhash":
            result = minhash(value, k=k, num_hashes=size)
        elif kind == "fracminhash":
            result = fracminhash(value, k=k, scaled=size)
        else:
            raise ValueError("kind must be kmer, minhash, or fracminhash")
        _json_output(result.to_dict())
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("search")
def search_command(
    target: Annotated[str, typer.Argument(help="Raw target DNA sequence.")],
    pattern: Annotated[str, typer.Argument(help="Exact/IUPAC/regex motif.")],
    mode: Annotated[str, typer.Option(help="Motif mode: exact, iupac, or regex.")] = "iupac",
) -> None:
    """Search one bounded motif on both strands."""

    try:
        result = scan_motif(
            _normalized_sequence(target),
            pattern,
            mode=mode,  # type: ignore[arg-type]
        )
        _json_output(result.to_dict())
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("orfs")
def orfs_command(
    sequence: Annotated[str, typer.Argument(help="Raw DNA sequence text.")],
    min_length: Annotated[int, typer.Option("--min-length", help="Minimum nucleotide length.")] = 0,
    complete_only: Annotated[
        bool, typer.Option("--complete-only/--allow-incomplete", help="Require a stop codon.")
    ] = True,
) -> None:
    """Scan six reading frames for start-anchored ORFs."""

    try:
        result = scan_orfs(
            _normalized_sequence(sequence),
            min_length=min_length,
            require_complete=complete_only,
        )
        _json_output(result.to_dict())
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("compare")
def compare_command(
    left: Annotated[str, typer.Argument(help="Left raw DNA sequence.")],
    right: Annotated[str, typer.Argument(help="Right raw DNA sequence.")],
    method: Annotated[
        str, typer.Option(help="exact, hamming, levenshtein, or kmer_jaccard.")
    ] = "exact",
    k: Annotated[int, typer.Option(help="k value for k-mer methods.")] = 3,
) -> None:
    """Compare two normalized sequences with one deterministic method."""

    try:
        left_sequence = _normalized_sequence(left)
        right_sequence = _normalized_sequence(right)
        result = (
            compare(
                left_sequence,
                right_sequence,
                method=method,
                k=k,
            )
            if method.startswith("kmer_")
            else compare(
                left_sequence,
                right_sequence,
                method=method,
            )
        )
        _json_output(result.to_dict())
    except (DNAKitError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("report")
def report_command(
    input_path: Annotated[Path, typer.Argument(help="Input sequence file.")],
    output_path: Annotated[Path, typer.Argument(help="Self-contained HTML output file.")],
    input_format: Annotated[
        str | None, typer.Option("--input-format", help="Override the inferred input format.")
    ] = None,
    title: Annotated[str, typer.Option(help="Report title.")] = "DNAKit report",
    max_records: Annotated[
        int, typer.Option("--max-records", help="Strict input and report record limit.")
    ] = 1_000,
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Build a bounded, read-only HTML report without external resources."""

    try:
        records = read_set(
            input_path,
            format=input_format,
            config=ReadConfig(max_records=max_records),
        )
        artifact = build_html_report(
            records,
            title=title,
            max_records=max_records,
            results={
                "dataset_summary": {
                    "record_count": len(records),
                    "total_symbol_length": sum(record.sequence.symbol_length for record in records),
                }
            },
        )
        result = save_html_report(
            artifact,
            output_path,
            config=SaveConfig(overwrite=overwrite),
        )
        _json_output(result.to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("workflow")
def workflow_command(
    config_path: Annotated[Path, typer.Argument(help="Workflow JSON/YAML configuration.")],
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Validate the workflow without creating outputs."),
    ] = False,
    resume: Annotated[
        bool,
        typer.Option("--resume", help="Resume only checksum-verified output steps."),
    ] = False,
    progress: Annotated[
        bool,
        typer.Option("--progress/--no-progress", help="Show workflow progress on stderr."),
    ] = True,
) -> None:
    """Run a strict, bounded local workflow without commands or network access."""

    from dnakit.cli.workflow import main as workflow_main

    arguments = ["run", str(config_path)]
    if dry_run:
        arguments.append("--dry-run")
    if resume:
        arguments.append("--resume")
    if not progress:
        arguments.append("--no-progress")
    code = workflow_main(arguments)
    if code:
        raise typer.Exit(code=code)


@app.command("convert")
def convert_command(
    input_path: Annotated[Path, typer.Argument(help="Input sequence file.")],
    output_path: Annotated[Path, typer.Argument(help="Output sequence file.")],
    input_format: Annotated[
        str | None, typer.Option("--input-format", help="Override the inferred input format.")
    ] = None,
    output_format: Annotated[
        str | None, typer.Option("--output-format", help="Override the inferred output format.")
    ] = None,
    overwrite: Annotated[
        bool, typer.Option("--overwrite", help="Replace an existing output atomically.")
    ] = False,
    progress: Annotated[
        bool, typer.Option("--progress/--no-progress", help="Show record progress on stderr.")
    ] = True,
) -> None:
    """Stream records between supported FASTA/FASTQ/table formats."""

    try:
        with read(input_path, format=input_format) as records:
            result = write(
                _tracked(records, enabled=progress),
                output_path,
                format=output_format,
                config=WriteConfig(overwrite=overwrite),
            )
        _json_output(result.to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("deduplicate")
def deduplicate_command(
    input_path: Annotated[Path, typer.Argument(help="Input sequence file.")],
    output_path: Annotated[Path, typer.Argument(help="Non-redundant output file.")],
    equivalence: Annotated[
        str,
        typer.Option(help="Duplicate equivalence: exact or reverse_complement."),
    ] = "exact",
    input_format: Annotated[str | None, typer.Option("--input-format")] = None,
    output_format: Annotated[str | None, typer.Option("--output-format")] = None,
    representative: Annotated[
        str, typer.Option(help="Representative policy: first, last, or best_quality.")
    ] = "first",
    conflict_field: Annotated[str | None, typer.Option(help="Metadata field to audit.")] = None,
    conflict_policy: Annotated[
        str,
        typer.Option(help="Conflict policy: error, drop_group, keep_representative, or keep_all."),
    ] = "error",
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Materialize and deduplicate records with an auditable mapping."""

    try:
        records = read_set(input_path, format=input_format)
        result = deduplicate(
            records,
            equivalence=equivalence,  # type: ignore[arg-type]
            config=DeduplicationConfig(
                representative_policy=representative,  # type: ignore[arg-type]
                conflict_field=conflict_field,
                conflict_policy=conflict_policy,  # type: ignore[arg-type]
            ),
        )
        write(
            result.records,
            output_path,
            format=output_format,
            config=WriteConfig(overwrite=overwrite),
        )
        _json_output(result.to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


_SAFE_SPLIT_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*", flags=re.ASCII)


def _parse_ratios(value: str) -> dict[str, float]:
    if not isinstance(value, str):
        raise ValueError("ratios must look like train=0.8,valid=0.1,test=0.1")
    ratios: dict[str, float] = {}
    for item in value.split(","):
        pair = item.split("=", maxsplit=1)
        if len(pair) != 2:
            raise ValueError("ratios must look like train=0.8,valid=0.1,test=0.1")
        name = pair[0].strip()
        if _SAFE_SPLIT_NAME.fullmatch(name) is None:
            raise ValueError(
                "split names must contain only ASCII letters, digits, '_' or '-', "
                "and must start with a letter or digit"
            )
        if name in ratios:
            raise ValueError(f"duplicate split name: {name!r}")
        try:
            ratios[name] = float(pair[1])
        except ValueError as exc:
            raise ValueError("ratios must look like train=0.8,valid=0.1,test=0.1") from exc
    return ratios


def _safe_output_suffix(value: str) -> str:
    normalized = value.strip().lower().lstrip(".")
    allowed = {
        "fa",
        "fasta",
        "fna",
        "ffn",
        "fq",
        "fastq",
        "csv",
        "tsv",
        "json",
        "jsonl",
    }
    if normalized not in allowed:
        raise ValueError("output-format must be FASTA, FASTQ, CSV, TSV, JSON, or JSONL")
    return normalized


def _write_manifest(path: Path, payload: object, *, overwrite: bool) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary_name)
    os.close(descriptor)
    try:
        temporary_path.write_text(
            json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if overwrite:
            os.replace(temporary_path, path)
        else:
            os.link(temporary_path, path)
            temporary_path.unlink()
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def _commit_staged_files(files: Sequence[tuple[Path, Path]], *, overwrite: bool) -> None:
    """Install a staged split as one rollback-capable local transaction."""

    targets = tuple(target for _, target in files)
    if len(targets) != len(set(targets)):
        raise ValueError("split output paths collide with each other or the manifest")
    existing = [str(path) for path in targets if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing split outputs: {existing}")
    invalid_targets = [str(path) for path in targets if path.exists() and not path.is_file()]
    if invalid_targets:
        raise IsADirectoryError(f"Split output targets must be files: {invalid_targets}")

    installed: list[Path] = []
    backups: list[tuple[Path, Path]] = []
    try:
        if overwrite:
            for index, (_, target) in enumerate(files):
                if target.exists():
                    backup = files[0][0].parent / f".backup-{index}"
                    os.replace(target, backup)
                    backups.append((target, backup))
            for staged, target in files:
                os.replace(staged, target)
                installed.append(target)
        else:
            for staged, target in files:
                os.link(staged, target)
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


@app.command("split")
def split_command(
    input_path: Annotated[Path, typer.Argument(help="Input sequence file.")],
    output_dir: Annotated[Path, typer.Argument(help="Directory for named subset files.")],
    ratios: Annotated[
        str,
        typer.Option(
            help="Named fractions for ratio methods, e.g. train=0.8,test=0.2; omit with kfold."
        ),
    ] = "train=0.8,valid=0.1,test=0.1",
    method: Annotated[
        str,
        typer.Option(help="Split method: random, hash, kfold, stratified, group, or similarity."),
    ] = "random",
    metadata_key: Annotated[
        str | None, typer.Option(help="Metadata key used by stratified or group splitting.")
    ] = None,
    seed: Annotated[
        int, typer.Option(help="Deterministic random seed; keep fixed across folds.")
    ] = 42,
    n_splits: Annotated[
        int | None, typer.Option(help="Number of folds for method=kfold (default 10).")
    ] = None,
    fold_id: Annotated[
        int | None,
        typer.Option(help="Validation fold for method=kfold: 0 to n-splits-1 (default 0)."),
    ] = None,
    shuffle: Annotated[
        bool, typer.Option("--shuffle/--no-shuffle", help="Shuffle before partitioning.")
    ] = True,
    input_format: Annotated[str | None, typer.Option("--input-format")] = None,
    output_format: Annotated[str, typer.Option("--output-format")] = "fasta",
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Create named dataset subsets and an assignments JSON manifest."""

    try:
        records = read_set(input_path, format=input_format)
        result = split(
            records,
            config=SplitConfig(
                method=method,  # type: ignore[arg-type]
                ratios=_parse_ratios(ratios),
                seed=seed,
                metadata_key=metadata_key,
                n_splits=n_splits,
                fold_id=fold_id,
                shuffle=shuffle,
            ),
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        suffix = _safe_output_suffix(output_format)
        targets = {subset.name: output_dir / f"{subset.name}.{suffix}" for subset in result.subsets}
        manifest_path = output_dir / "assignments.json"
        final_paths = (*targets.values(), manifest_path)
        if len(final_paths) != len(set(final_paths)):
            raise ValueError("a split output path collides with assignments.json")
        existing = [str(path) for path in final_paths if path.exists()]
        if existing and not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing split outputs: {existing}")
        with tempfile.TemporaryDirectory(prefix=".dnakit-split-", dir=output_dir) as staging:
            staging_dir = Path(staging)
            staged_files: list[tuple[Path, Path]] = []
            for subset in result.subsets:
                staged_path = staging_dir / f"{subset.name}.{suffix}"
                write(
                    subset.records,
                    staged_path,
                    format=suffix,
                )
                staged_files.append((staged_path, targets[subset.name]))
            staged_manifest = staging_dir / "assignments.json"
            _write_manifest(staged_manifest, result.to_dict(), overwrite=False)
            staged_files.append((staged_manifest, manifest_path))
            _commit_staged_files(staged_files, overwrite=overwrite)
        _json_output(result.to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


def _comma_fields(value: str | None) -> tuple[str, ...]:
    if value is None:
        return ()
    fields = tuple(item.strip() for item in value.split(",") if item.strip())
    if not fields:
        raise ValueError("A comma-separated field list must not be empty.")
    return fields


def _custom_data_mapping(
    *,
    input_format: str,
    sequence_field: str | None,
    label_fields: str | None,
    id_field: str | None,
    split_field: str | None,
    record_path: str | None,
    delimiter: str | None,
    sheet_name: str | None,
    table_name: str | None,
) -> DataMapping | None:
    from dnakit.predictions import DataMapping

    explicit = (
        any(
            value is not None
            for value in (
                sequence_field,
                label_fields,
                id_field,
                split_field,
                record_path,
                delimiter,
                sheet_name,
                table_name,
            )
        )
        or input_format != "auto"
    )
    if not explicit:
        return None
    return DataMapping(
        format=input_format,  # type: ignore[arg-type]
        sequence_field=sequence_field or "sequence",
        label_fields=_comma_fields(label_fields),
        id_field=id_field,
        split_field=split_field,
        record_path=record_path,
        delimiter=delimiter,
        sheet_name=sheet_name,
        table_name=table_name,
    )


def _llm_schema_config(
    *,
    use_llm: bool,
    model: str,
    base_url: str,
    api_key_env: str,
    share_values: bool,
) -> LLMSchemaConfig | None:
    from dnakit.predictions import LLMSchemaConfig

    if not use_llm:
        return None
    return LLMSchemaConfig(
        model=model,
        base_url=base_url,
        api_key_env=api_key_env,
        include_sample_values=share_values,
    )


@app.command("property-models")
def property_models_command() -> None:
    """List DNA models available for full fine-tuning and linear probing."""

    from dnakit.predictions import available_custom_training_models

    _json_output([model.to_dict() for model in available_custom_training_models()])


@app.command("inspect-property-data")
def inspect_property_data_command(
    input_path: Annotated[Path, typer.Argument(help="Local user dataset to inspect.")],
    input_format: Annotated[
        str,
        typer.Option(
            "--input-format",
            help="auto, csv, tsv, json, jsonl, yaml, xml, excel, parquet, sqlite, or FASTA-like.",
        ),
    ] = "auto",
    sequence_field: Annotated[str | None, typer.Option(help="Sequence field or dot path.")] = None,
    label_fields: Annotated[
        str | None, typer.Option(help="Comma-separated label fields or dot paths.")
    ] = None,
    id_field: Annotated[str | None, typer.Option(help="Optional record-ID field.")] = None,
    split_field: Annotated[
        str | None, typer.Option(help="Optional train/validation/test field.")
    ] = None,
    record_path: Annotated[
        str | None, typer.Option(help="Optional structured-record container path.")
    ] = None,
    delimiter: Annotated[str | None, typer.Option(help="Optional one-character delimiter.")] = None,
    sheet_name: Annotated[str | None, typer.Option(help="Optional Excel sheet name.")] = None,
    table_name: Annotated[str | None, typer.Option(help="Optional SQLite table name.")] = None,
    use_llm: Annotated[
        bool, typer.Option("--use-llm/--no-llm", help="Use bounded LLM schema analysis.")
    ] = False,
    llm_model: Annotated[str, typer.Option(help="OpenAI-compatible model name.")] = (
        "deepseek-v4-flash"
    ),
    llm_base_url: Annotated[str, typer.Option(help="OpenAI-compatible API base URL.")] = (
        "https://api.deepseek.com"
    ),
    llm_api_key_env: Annotated[
        str, typer.Option(help="Environment variable containing the LLM API key.")
    ] = "DEEPSEEK_API_KEY",
    llm_share_values: Annotated[
        bool,
        typer.Option(
            "--llm-share-values/--llm-structure-only",
            help="Share bounded sample values or only field/type structure.",
        ),
    ] = True,
) -> None:
    """Inspect arbitrary local data and print the selected declarative mapping."""

    from dnakit.predictions import inspect_property_data

    try:
        mapping = _custom_data_mapping(
            input_format=input_format,
            sequence_field=sequence_field,
            label_fields=label_fields,
            id_field=id_field,
            split_field=split_field,
            record_path=record_path,
            delimiter=delimiter,
            sheet_name=sheet_name,
            table_name=table_name,
        )
        llm = _llm_schema_config(
            use_llm=use_llm,
            model=llm_model,
            base_url=llm_base_url,
            api_key_env=llm_api_key_env,
            share_values=llm_share_values,
        )
        _json_output(inspect_property_data(input_path, mapping=mapping, llm=llm).to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("train-property")
def train_property_command(
    input_path: Annotated[Path, typer.Argument(help="Local supervised dataset.")],
    output_dir: Annotated[Path, typer.Argument(help="New model artifact directory.")],
    task_type: Annotated[
        str, typer.Option(help="Required task type: classification or regression.")
    ],
    model: Annotated[str, typer.Option(help="Registered DNA foundation model.")] = "dnabert2",
    finetune_method: Annotated[str, typer.Option(help="full (default) or linear_probe.")] = "full",
    classification_mode: Annotated[
        str,
        typer.Option(help="auto, binary, multiclass, multilabel, or multioutput."),
    ] = "auto",
    input_format: Annotated[str, typer.Option("--input-format")] = "auto",
    sequence_field: Annotated[str | None, typer.Option(help="Sequence field or dot path.")] = None,
    label_fields: Annotated[
        str | None, typer.Option(help="Comma-separated target fields or dot paths.")
    ] = None,
    id_field: Annotated[str | None, typer.Option(help="Optional record-ID field.")] = None,
    split_field: Annotated[str | None, typer.Option(help="Optional split field.")] = None,
    record_path: Annotated[
        str | None, typer.Option(help="Structured-record container path.")
    ] = None,
    delimiter: Annotated[str | None, typer.Option(help="Optional one-character delimiter.")] = None,
    sheet_name: Annotated[str | None, typer.Option(help="Optional Excel sheet name.")] = None,
    table_name: Annotated[str | None, typer.Option(help="Optional SQLite table name.")] = None,
    use_llm: Annotated[bool, typer.Option("--use-llm/--no-llm")] = False,
    llm_model: Annotated[str, typer.Option(help="OpenAI-compatible schema model.")] = (
        "deepseek-v4-flash"
    ),
    llm_base_url: Annotated[str, typer.Option(help="OpenAI-compatible API base URL.")] = (
        "https://api.deepseek.com"
    ),
    llm_api_key_env: Annotated[
        str, typer.Option(help="Environment variable containing the LLM API key.")
    ] = "DEEPSEEK_API_KEY",
    llm_share_values: Annotated[
        bool, typer.Option("--llm-share-values/--llm-structure-only")
    ] = True,
    train_ratio: Annotated[float, typer.Option(help="Training fraction.")] = 0.8,
    validation_ratio: Annotated[float, typer.Option(help="Validation fraction.")] = 0.1,
    test_ratio: Annotated[float, typer.Option(help="Test fraction.")] = 0.1,
    seed: Annotated[int, typer.Option(help="Random split/training seed.")] = 42,
    epochs: Annotated[int, typer.Option(help="Maximum training epochs.")] = 10,
    batch_size: Annotated[int, typer.Option(help="Training batch size.")] = 8,
    learning_rate: Annotated[
        float | None, typer.Option(help="Optional learning-rate override.")
    ] = None,
    checkpoint_dir: Annotated[Path | None, typer.Option(help="Checkpoint cache root.")] = None,
    checkpoint_path: Annotated[Path | None, typer.Option(help="Existing local checkpoint.")] = None,
    max_length: Annotated[int | None, typer.Option(help="Tokenizer/model input limit.")] = None,
    device: Annotated[str, typer.Option(help="auto, cpu, cuda, or cuda:N.")] = "auto",
    allow_remote_code: Annotated[
        bool,
        typer.Option(
            "--allow-remote-code/--deny-remote-code",
            help="Opt in to reviewed checkpoint-provided Python code.",
        ),
    ] = False,
    overwrite: Annotated[bool, typer.Option("--overwrite")] = False,
) -> None:
    """Train a custom DNA property classifier or regressor and emit full metrics."""

    from dnakit.predictions import CustomPropertyTrainingConfig, train_custom_property_model

    try:
        mapping = _custom_data_mapping(
            input_format=input_format,
            sequence_field=sequence_field,
            label_fields=label_fields,
            id_field=id_field,
            split_field=split_field,
            record_path=record_path,
            delimiter=delimiter,
            sheet_name=sheet_name,
            table_name=table_name,
        )
        llm = _llm_schema_config(
            use_llm=use_llm,
            model=llm_model,
            base_url=llm_base_url,
            api_key_env=llm_api_key_env,
            share_values=llm_share_values,
        )
        config = CustomPropertyTrainingConfig(
            task_type=task_type,  # type: ignore[arg-type]
            model=model,
            finetune_method=finetune_method,  # type: ignore[arg-type]
            classification_mode=classification_mode,  # type: ignore[arg-type]
            mapping=mapping,
            llm=llm,
            output_dir=output_dir,
            train_ratio=train_ratio,
            validation_ratio=validation_ratio,
            test_ratio=test_ratio,
            seed=seed,
            epochs=epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            checkpoint_dir=checkpoint_dir,
            checkpoint_path=checkpoint_path,
            max_length=max_length,
            device=device,
            allow_remote_code=allow_remote_code,
            overwrite=overwrite,
        )
        _json_output(train_custom_property_model(input_path, config=config).to_dict())
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


@app.command("predict-property")
def predict_property_command(
    input_path: Annotated[Path, typer.Argument(help="Local sequence dataset to predict.")],
    model_dir: Annotated[Path, typer.Argument(help="Saved custom-property model directory.")],
    input_format: Annotated[str, typer.Option("--input-format")] = "auto",
    sequence_field: Annotated[str | None, typer.Option(help="Sequence field or dot path.")] = None,
    id_field: Annotated[str | None, typer.Option(help="Optional record-ID field.")] = None,
    record_path: Annotated[
        str | None, typer.Option(help="Structured-record container path.")
    ] = None,
    delimiter: Annotated[str | None, typer.Option(help="Optional one-character delimiter.")] = None,
    sheet_name: Annotated[str | None, typer.Option(help="Optional Excel sheet name.")] = None,
    table_name: Annotated[str | None, typer.Option(help="Optional SQLite table name.")] = None,
    batch_size: Annotated[int, typer.Option(help="Inference batch size.")] = 8,
    device: Annotated[str, typer.Option(help="auto, cpu, cuda, or cuda:N.")] = "auto",
    threshold: Annotated[float | None, typer.Option(help="Multilabel threshold override.")] = None,
    allow_remote_code: Annotated[
        bool, typer.Option("--allow-remote-code/--deny-remote-code")
    ] = False,
) -> None:
    """Predict properties with a saved user-defined model."""

    from dnakit.predictions import predict_custom_properties

    try:
        mapping = _custom_data_mapping(
            input_format=input_format,
            sequence_field=sequence_field,
            label_fields=None,
            id_field=id_field,
            split_field=None,
            record_path=record_path,
            delimiter=delimiter,
            sheet_name=sheet_name,
            table_name=table_name,
        )
        _json_output(
            predict_custom_properties(
                input_path,
                model_dir=model_dir,
                mapping=mapping,
                batch_size=batch_size,
                device=device,
                threshold=threshold,
                allow_remote_code=allow_remote_code,
            ).to_dict()
        )
    except (DNAKitError, OSError, TypeError, ValueError) as exc:
        _abort(exc)


def run() -> None:
    """Run the command-line application."""
    app()
