"""Run reproducible, bounded microbenchmarks for native DNAKit operations."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.metadata
import inspect
import json
import os
import platform
import random
import statistics
import sys
import tempfile
import time
import tracemalloc
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Final

from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

import dnakit
from dnakit import DNASequence
from dnakit.descriptors import gc_at_content
from dnakit.fingerprints import kmer_fingerprint, minhash
from dnakit.ops import reverse_complement
from dnakit.similarity import subsequence_search
from dnakit.standardize import normalize

REPORT_SCHEMA_VERSION: Final = "dnakit.benchmark.v1"
MAX_SIZE: Final = 200_000
MAX_REPEATS: Final = 50
MAX_WARMUPS: Final = 20
MAX_CASES: Final = 500
MAX_TOTAL_NUCLEOTIDE_RUNS: Final = 20_000_000
DEFAULT_TASKS: Final = (
    "construct",
    "normalize",
    "gc_content",
    "kmer_fingerprint",
    "minhash",
    "subsequence_search",
    "reverse_complement",
)
DEFAULT_IMPLEMENTATIONS: Final = ("dnakit",)
SUPPORTED_IMPLEMENTATIONS: Final = frozenset({"dnakit", "biopython"})
_BIOPYTHON_TASKS: Final = frozenset({"construct", "gc_content", "reverse_complement"})

BenchmarkCallable = Callable[[], object]


def _installed_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _environment() -> dict[str, Any]:
    return {
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_compiler": platform.python_compiler(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or None,
        "cpu_count": os.cpu_count(),
        "dnakit": dnakit.__version__,
        "rich": _installed_version("rich"),
        "biopython": _installed_version("biopython"),
    }


def _random_dna(length: int, seed: int) -> str:
    generator = random.Random(seed)
    return "".join(generator.choices("ACGT", k=length))


def _build_dnakit_task(name: str, sequence_text: str, seed: int) -> BenchmarkCallable:
    sequence = DNASequence(sequence_text)
    if name == "construct":
        return lambda: DNASequence(sequence_text)
    if name == "normalize":
        raw = sequence_text.lower()
        return lambda: normalize(raw)
    if name == "gc_content":
        return lambda: gc_at_content(sequence)
    if name == "kmer_fingerprint":
        return lambda: kmer_fingerprint(sequence, k=4, canonical=True)
    if name == "minhash":
        return lambda: minhash(sequence, k=15, num_hashes=256, seed=seed)
    if name == "subsequence_search":
        query_length = min(21, len(sequence_text))
        query = DNASequence(sequence_text[:query_length])
        return lambda: subsequence_search(query, sequence, max_matches=len(sequence_text) + 1)
    if name == "reverse_complement":
        return lambda: reverse_complement(sequence)
    raise ValueError(f"Unknown benchmark task: {name}")


def _build_biopython_task(name: str, sequence_text: str) -> BenchmarkCallable:
    try:
        from Bio.Seq import Seq
        from Bio.SeqUtils import gc_fraction
    except ImportError as exc:
        raise ValueError("Biopython is not installed for the requested comparison.") from exc
    sequence = Seq(sequence_text)
    if name == "construct":
        return lambda: Seq(sequence_text)
    if name == "gc_content":
        return lambda: gc_fraction(sequence)  # type: ignore[no-untyped-call]
    if name == "reverse_complement":
        return sequence.reverse_complement
    raise ValueError(f"Biopython has no configured equivalent for benchmark task: {name}")


def _build_task(implementation: str, name: str, sequence_text: str, seed: int) -> BenchmarkCallable:
    if implementation == "dnakit":
        return _build_dnakit_task(name, sequence_text, seed)
    if implementation == "biopython":
        return _build_biopython_task(name, sequence_text)
    raise ValueError(f"Unknown benchmark implementation: {implementation}")


def _source_object(implementation: str, task: str) -> Any | None:
    if implementation == "dnakit":
        return {
            "construct": DNASequence,
            "normalize": normalize,
            "gc_content": gc_at_content,
            "kmer_fingerprint": kmer_fingerprint,
            "minhash": minhash,
            "subsequence_search": subsequence_search,
            "reverse_complement": reverse_complement,
        }.get(task)
    if implementation == "biopython":
        try:
            from Bio.Seq import Seq
            from Bio.SeqUtils import gc_fraction
        except ImportError:
            return None
        return {
            "construct": Seq,
            "gc_content": gc_fraction,
            "reverse_complement": Seq.reverse_complement,
        }.get(task)
    return None


def _source_metrics(implementations: Sequence[str], tasks: Sequence[str]) -> dict[str, Any]:
    """Count non-empty, non-comment source lines for configured public callables."""

    metrics: dict[str, Any] = {}
    for implementation in implementations:
        entries: list[dict[str, object]] = []
        total_lines = 0
        for task in tasks:
            source_object = _source_object(implementation, task)
            if source_object is None:
                continue
            try:
                lines, start_line = inspect.getsourcelines(source_object)
                source_file = inspect.getsourcefile(source_object)
            except (OSError, TypeError):
                continue
            count = sum(bool(line.strip()) and not line.lstrip().startswith("#") for line in lines)
            total_lines += count
            entries.append(
                {
                    "task": task,
                    "qualified_name": (
                        f"{getattr(source_object, '__module__', '')}."
                        f"{getattr(source_object, '__qualname__', type(source_object).__name__)}"
                    ).strip("."),
                    "source_file": source_file,
                    "start_line": start_line,
                    "nonblank_noncomment_lines": count,
                }
            )
        metrics[implementation] = {
            "definition": (
                "inspect.getsourcelines public callable; counts non-empty lines whose "
                "first non-whitespace character is not #; not a quality metric"
            ),
            "entries": entries,
            "total_nonblank_noncomment_lines": total_lines,
        }
    return metrics


def _summary(values: Sequence[int]) -> dict[str, float | int]:
    if not values:
        raise ValueError("At least one measurement is required.")
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "mean": statistics.fmean(values),
        "maximum": max(values),
        "standard_deviation": statistics.pstdev(values),
    }


def _measure(operation: BenchmarkCallable) -> tuple[int, int]:
    tracemalloc.start()
    try:
        started = time.perf_counter_ns()
        result = operation()
        elapsed = time.perf_counter_ns() - started
        _current, peak = tracemalloc.get_traced_memory()
        # Keep the return value alive until after peak-memory sampling.
        del result
    finally:
        tracemalloc.stop()
    return elapsed, peak


def _validate_request(
    sizes: Sequence[int],
    repeats: int,
    warmups: int,
    seed: int,
    tasks: Sequence[str],
    implementations: Sequence[str],
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[str, ...]]:
    if not sizes:
        raise ValueError("sizes must contain at least one positive integer.")
    resolved_sizes = tuple(sizes)
    if any(isinstance(size, bool) or not isinstance(size, int) or size <= 0 for size in sizes):
        raise ValueError("Every size must be a positive integer.")
    if max(resolved_sizes) > MAX_SIZE:
        raise ValueError(f"A size exceeds the hard safety ceiling of {MAX_SIZE} nt.")
    if isinstance(repeats, bool) or not isinstance(repeats, int) or not 1 <= repeats <= MAX_REPEATS:
        raise ValueError(f"repeats must be an integer in [1, {MAX_REPEATS}].")
    if isinstance(warmups, bool) or not isinstance(warmups, int) or not 0 <= warmups <= MAX_WARMUPS:
        raise ValueError(f"warmups must be an integer in [0, {MAX_WARMUPS}].")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise ValueError("seed must be an integer in [0, 2**64).")
    if not tasks:
        raise ValueError("tasks must contain at least one benchmark name.")
    resolved_tasks = tuple(dict.fromkeys(tasks))
    unknown = sorted(set(resolved_tasks) - set(DEFAULT_TASKS))
    if unknown:
        raise ValueError(f"Unknown benchmark tasks: {', '.join(unknown)}")
    if not implementations:
        raise ValueError("implementations must contain at least one implementation name.")
    resolved_implementations = tuple(dict.fromkeys(implementations))
    unknown_implementations = sorted(set(resolved_implementations) - SUPPORTED_IMPLEMENTATIONS)
    if unknown_implementations:
        raise ValueError(f"Unknown benchmark implementations: {', '.join(unknown_implementations)}")
    if "biopython" in resolved_implementations and _installed_version("biopython") is None:
        raise ValueError("Biopython is not installed for the requested comparison.")
    executable_pairs = sum(
        1
        for implementation in resolved_implementations
        for task in resolved_tasks
        if implementation == "dnakit" or task in _BIOPYTHON_TASKS
    )
    case_count = len(resolved_sizes) * executable_pairs
    if case_count > MAX_CASES:
        raise ValueError(f"The request exceeds the hard ceiling of {MAX_CASES} cases.")
    nucleotide_runs = sum(resolved_sizes) * executable_pairs * (repeats + warmups)
    if nucleotide_runs > MAX_TOTAL_NUCLEOTIDE_RUNS:
        raise ValueError(
            "The request exceeds the aggregate safety ceiling of "
            f"{MAX_TOTAL_NUCLEOTIDE_RUNS} nucleotide-runs."
        )
    return resolved_sizes, resolved_tasks, resolved_implementations


def run_benchmarks(
    *,
    sizes: Sequence[int],
    repeats: int,
    warmups: int,
    seed: int,
    tasks: Sequence[str] = DEFAULT_TASKS,
    implementations: Sequence[str] = DEFAULT_IMPLEMENTATIONS,
    show_progress: bool = True,
) -> dict[str, Any]:
    """Return an in-memory benchmark report without writing files."""

    resolved_sizes, resolved_tasks, resolved_implementations = _validate_request(
        sizes, repeats, warmups, seed, tasks, implementations
    )
    executable_pairs = tuple(
        (implementation, task)
        for implementation in resolved_implementations
        for task in resolved_tasks
        if implementation == "dnakit" or task in _BIOPYTHON_TASKS
    )
    total_steps = len(resolved_sizes) * len(executable_pairs) * (repeats + warmups)
    progress = Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=Console(stderr=True),
        disable=not show_progress,
    )
    cases: list[dict[str, Any]] = []
    with progress:
        progress_task = progress.add_task("DNAKit benchmark", total=total_steps)
        for size_index, size in enumerate(resolved_sizes):
            sequence_text = _random_dna(size, seed + size_index)
            input_digest = hashlib.sha256(sequence_text.encode("ascii")).hexdigest()
            for implementation, task_name in executable_pairs:
                operation = _build_task(implementation, task_name, sequence_text, seed)
                for _ in range(warmups):
                    operation()
                    progress.advance(progress_task)
                durations: list[int] = []
                peaks: list[int] = []
                for _ in range(repeats):
                    duration, peak = _measure(operation)
                    durations.append(duration)
                    peaks.append(peak)
                    progress.advance(progress_task)
                cases.append(
                    {
                        "implementation": implementation,
                        "task": task_name,
                        "size_nt": size,
                        "input_sha256": input_digest,
                        "duration_ns": _summary(durations),
                        "peak_tracemalloc_bytes": _summary(peaks),
                        "samples": [
                            {"duration_ns": duration, "peak_tracemalloc_bytes": peak}
                            for duration, peak in zip(durations, peaks, strict=True)
                        ],
                    }
                )
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark_kind": "local_microbenchmark_not_a_cross_machine_claim",
        "timer": "time.perf_counter_ns",
        "memory_meter": "tracemalloc_peak_python_allocations",
        "environment": _environment(),
        "source_code_metrics": _source_metrics(resolved_implementations, resolved_tasks),
        "parameters": {
            "sizes": list(resolved_sizes),
            "repeats": repeats,
            "warmups": warmups,
            "seed": seed,
            "tasks": list(resolved_tasks),
            "implementations": list(resolved_implementations),
            "comparison_task_policy": {
                "dnakit": list(resolved_tasks),
                "biopython": [task for task in resolved_tasks if task in _BIOPYTHON_TASKS],
            },
            "safety_limits": {
                "max_size": MAX_SIZE,
                "max_repeats": MAX_REPEATS,
                "max_warmups": MAX_WARMUPS,
                "max_cases": MAX_CASES,
                "max_total_nucleotide_runs": MAX_TOTAL_NUCLEOTIDE_RUNS,
            },
        },
        "cases": cases,
        "interpretation": (
            "Compare repeated runs only under controlled hardware and software conditions. "
            "tracemalloc excludes native allocations outside Python's allocator. Source-line "
            "counts describe selected public callable bodies and are not a quality metric."
        ),
    }


def write_report(report: dict[str, Any], output: Path, *, force: bool = False) -> None:
    """Atomically write a report, refusing accidental replacement by default."""

    resolved = output.expanduser().resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{resolved.name}.", suffix=".tmp", dir=resolved.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(payload)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if force:
            os.replace(temporary_name, resolved)
        else:
            try:
                os.link(temporary_name, resolved)
            except FileExistsError as exc:
                raise FileExistsError(
                    f"Output already exists: {resolved}. Pass --force to replace it."
                ) from exc
            os.unlink(temporary_name)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(temporary_name)
        raise


def _csv_ints(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("sizes must be comma-separated integers") from exc
    if not values:
        raise argparse.ArgumentTypeError("sizes cannot be empty")
    return values


def _csv_strings(value: str) -> tuple[str, ...]:
    values = tuple(item.strip() for item in value.split(",") if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("tasks cannot be empty")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=_csv_ints, default=(100, 1_000, 10_000))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--tasks", type=_csv_strings, default=DEFAULT_TASKS)
    parser.add_argument(
        "--implementations",
        type=_csv_strings,
        default=DEFAULT_IMPLEMENTATIONS,
        help="Comma-separated implementations: dnakit,biopython",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/results/local_benchmark_report.json"),
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        report = run_benchmarks(
            sizes=arguments.sizes,
            repeats=arguments.repeats,
            warmups=arguments.warmups,
            seed=arguments.seed,
            tasks=arguments.tasks,
            implementations=arguments.implementations,
        )
        write_report(report, arguments.output, force=arguments.force)
    except (FileExistsError, ValueError) as exc:
        parser.error(str(exc))
    output = arguments.output.expanduser().resolve()
    Console(stderr=True).print(f"[green]Benchmark report written:[/green] {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
