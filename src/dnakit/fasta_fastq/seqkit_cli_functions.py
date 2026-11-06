"""SeqKit command-line facades: fixed argument whitelists, no shell, bounded execution.

Each public function runs one SeqKit subcommand (2.14.0 is the tested version) through
``dnakit.backends.external.execute_bounded_command``. The user installs SeqKit (conda, apt or a
release binary) and either puts it on PATH or passes ``seqkit_path``. Outputs are written to
explicit paths; standard output is never the only result. SeqKit's own ``version`` subcommand is
not used because it checks for updates over the network; the version is read from ``--help``.
Results are JSON-serialisable dictionaries with status, backend, version, outputs, parsed rows and
a bounded log tail.
"""

from __future__ import annotations

import math
import os
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_PROGRAM = "seqkit"
_TESTED_VERSION = "2.14.0"
_CITATION = "https://doi.org/10.1002/imt2.191"
_INSTALL_HINT = (
    "Install SeqKit 2.14.0 or later (for example `conda install -c bioconda seqkit`, or a release "
    "binary from https://github.com/shenwei356/seqkit/releases), then put it on PATH or pass "
    "seqkit_path."
)
_TIMEOUT_SECONDS = 3600.0
_VERSION_TIMEOUT_SECONDS = 60.0
_MAX_OUTPUT_BYTES = 5_000_000
_MAX_ARTIFACT_BYTES = 2_000_000_000
_LOG_TAIL_CHARS = 4000
_VERSION_LINE = re.compile(r"^Version:\s*(\S+)", re.MULTILINE)
_BUNDLE_SUMMARY = re.compile(r"Written (\d+) BAM records to (\d+) loci and (\d+) bundles")
_REGION = re.compile(r"^-?\d+:-?\d+$")
_ALIGN_PARAMS = re.compile(r"^-?\d+,-?\d+,-?\d+,-?\d+$")
_BAM_FIELD = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")
_FILL = re.compile(r"^[A-Za-z-]+$")
_READ_FIELDS = ("ReadLen", "MeanQual", "GC", "GCSkew")
_FQ_ENCODINGS = ("sanger", "solexa", "illumina-1.3+", "illumina-1.5+", "illumina-1.8+")


@dataclass(frozen=True)
class _Operation:
    """A SeqKit subcommand and the only long flags this facade may pass to it."""

    subcommand: str
    flags: frozenset[str]


# Every operation name is the coverage-gap operation slug from matches/seqkit.csv.
_OPERATIONS: Mapping[str, _Operation] = {
    "flanking_sequences_of_features": _Operation(
        "subseq",
        frozenset({"--bed", "--gtf", "--feature", "--up-stream", "--down-stream", "--only-flank"}),
    ),
    "fastx_extended_statistics_n50_quality": _Operation(
        "stats", frozenset({"--all", "--tabular", "--fq-encoding", "--N"})
    ),
    "sequence_set_digest": _Operation(
        "sum", frozenset({"--include-id", "--remove-gaps", "--rna2dna"})
    ),
    "read_feature_histogram_monitor": _Operation(
        "watch",
        frozenset(
            {
                "--fields",
                "--bins",
                "--log",
                "--validate-seq",
                "--img",
                "--qual-ascii-base",
                "--dump",
                "--quiet-mode",
            }
        ),
    ),
    "retrieve_fastq_records_by_fasta": _Operation(
        "fa2fq", frozenset({"--fasta-file", "--only-positive-strand"})
    ),
    "mean_read_quality": _Operation(
        "fx2tab", frozenset({"--name", "--avg-qual", "--qual-ascii-base"})
    ),
    "local_alignment_search_short_queries": _Operation(
        "fish", frozenset({"--query-fastx", "--min-qual", "--aln-params", "--all", "--stranded"})
    ),
    "filter_records_by_sequence_motif": _Operation(
        "grep",
        frozenset(
            {
                "--by-seq",
                "--pattern",
                "--degenerate",
                "--max-mismatch",
                "--only-positive-strand",
                "--ignore-case",
                "--invert-match",
            }
        ),
    ),
    "count_matching_records": _Operation(
        "grep",
        frozenset(
            {
                "--count",
                "--by-name",
                "--by-seq",
                "--use-regexp",
                "--pattern",
                "--invert-match",
                "--ignore-case",
            }
        ),
    ),
    "duplicate_records": _Operation("duplicate", frozenset({"--times"})),
    "select_records_by_position": _Operation("range", frozenset({"--range"})),
    "first_genome_by_name_prefix": _Operation("head-genome", frozenset({"--mini-common-words"})),
    "split_by_region_subsequence": _Operation("split", frozenset({"--by-region", "--out-dir"})),
    "concatenate_same_id_sequences": _Operation(
        "concat", frozenset({"--full", "--fill", "--separator"})
    ),
    "make_duplicate_ids_unique": _Operation(
        "rename", frozenset({"--by-name", "--start-num", "--separator", "--rename-1st-rec"})
    ),
    "replace_sequence_regex": _Operation(
        "replace", frozenset({"--by-seq", "--pattern", "--replacement", "--ignore-case"})
    ),
    "bam_alignment_statistics": _Operation(
        "bam", frozenset({"--stat", "--idx-stat", "--quiet-mode"})
    ),
    "bam_feature_histogram": _Operation(
        "bam", frozenset({"--field", "--bins", "--log", "--dump", "--img", "--quiet-mode"})
    ),
    "bam_top_records_by_field": _Operation(
        "bam", frozenset({"--field", "--top-bam", "--top-size", "--quiet-mode"})
    ),
    "bam_bundle_partition": _Operation("bam", frozenset({"--bundle", "--quiet-mode"})),
    "merge_sliding_windows_to_bed": _Operation(
        "merge-slides", frozenset({"--min-overlap", "--max-gap"})
    ),
}


@dataclass(frozen=True)
class _Completed:
    """Captured output of one bounded SeqKit run."""

    text: str
    version: str | None
    executable: str
    command: tuple[str, ...]


# ---------------------------------------------------------------------------
# Validation and path handling
# ---------------------------------------------------------------------------


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty string.")
    return value


def _int(
    value: object, name: str, *, minimum: int | None = None, maximum: int | None = None
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _input(value: object, name: str = "input_path") -> Path:
    path = Path(_text(value, name)).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}")
    return path


def _inputs(values: object, name: str = "input_paths") -> tuple[Path, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError(f"{name} must be a non-empty list of file paths.")
    return tuple(_input(value, name) for value in values)


def _output(value: object, name: str = "output_path") -> Path:
    """A new file: its parent must exist and the file itself must not, so nothing is overwritten."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _new_directory(value: object, name: str = "output_dir") -> Path:
    """A new directory that SeqKit or this facade creates; it must not exist yet."""
    path = Path(_text(value, name)).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output directory already exists: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Parent of the output directory does not exist: {path.parent}")
    return path


def _patterns(values: object) -> tuple[str, ...]:
    if isinstance(values, str) or not isinstance(values, Sequence) or not values:
        raise ConfigurationError("patterns must be a non-empty list of strings.")
    for value in values:
        if not isinstance(value, str) or not value or "," in value:
            raise ConfigurationError(
                "each pattern must be a non-empty string without commas (SeqKit splits "
                "pattern lists on commas)."
            )
    return tuple(values)


def _check_options(operation: str, options: Sequence[str]) -> None:
    spec = _OPERATIONS[operation]
    for token in options:
        name = token.split("=", 1)[0]
        if not name.startswith("--") or name not in spec.flags:
            raise ConfigurationError(
                f"seqkit {spec.subcommand} option {name!r} is not allowed by this facade."
            )


def _resolve_executable(seqkit_path: str | None) -> str:
    """Find SeqKit on PATH, or use an explicit path. Never a hard-coded location."""
    if seqkit_path is None or seqkit_path == "":
        found = shutil.which(_PROGRAM)
        if found is None:
            raise BackendUnavailableError(
                "External executable 'seqkit' was not found on PATH.",
                hint=_INSTALL_HINT,
                code="SEQKIT_UNAVAILABLE",
            )
        return str(Path(found).absolute())
    if os.sep not in seqkit_path:
        found = shutil.which(seqkit_path)
        if found is None:
            raise BackendUnavailableError(
                f"External executable {seqkit_path!r} was not found on PATH.",
                hint=_INSTALL_HINT,
                code="SEQKIT_UNAVAILABLE",
            )
        return str(Path(found).absolute())
    candidate = Path(seqkit_path).expanduser()
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise BackendUnavailableError(
            f"SeqKit executable is missing or not executable: {seqkit_path}",
            hint=_INSTALL_HINT,
            code="SEQKIT_UNAVAILABLE",
        )
    return str(candidate.absolute())


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


def _version(executable: str, cwd: Path, verbose: bool) -> str | None:
    result = execute_bounded_command(
        executable,
        ("--help",),
        backend_id=_PROGRAM,
        cwd=cwd,
        timeout_seconds=_VERSION_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_OUTPUT_BYTES,
        verbose=verbose,
    )
    match = _VERSION_LINE.search(result.output)
    return match.group(1) if match else None


def _run(
    operation: str,
    options: Sequence[str],
    inputs: Sequence[Path],
    *,
    cwd: Path,
    output: Path | None = None,
    artifacts: Sequence[Path] = (),
    seqkit_path: str | None,
    verbose: bool,
) -> _Completed:
    """Run one whitelisted SeqKit subcommand. Positional inputs are absolute paths; no shell."""
    _check_options(operation, options)
    subcommand = _OPERATIONS[operation].subcommand
    executable = _resolve_executable(seqkit_path)
    version = _version(executable, cwd, verbose)
    monitored: tuple[Path, ...] = (output,) if output is not None else tuple(artifacts)
    arguments: tuple[str, ...] = (subcommand, "--quiet", *options, *(str(p) for p in inputs))
    if output is not None:
        arguments = (*arguments, "-o", str(output))
    result = execute_bounded_command(
        executable,
        arguments,
        backend_id=_PROGRAM,
        cwd=cwd,
        timeout_seconds=_TIMEOUT_SECONDS,
        max_output_bytes=_MAX_OUTPUT_BYTES,
        monitored_output_paths=monitored,
        max_monitored_output_bytes=_MAX_ARTIFACT_BYTES,
        verbose=verbose,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"SeqKit {subcommand} exited with code {result.return_code}.",
            code="SEQKIT_EXECUTION_FAILED",
            context={"operation": operation, "log_tail": _tail(result.output)},
        )
    missing = [str(path) for path in monitored if not path.is_file()]
    if missing:
        raise BackendExecutionError(
            f"SeqKit {subcommand} did not write the expected files: {missing}",
            code="SEQKIT_OUTPUT_MISSING",
            context={"operation": operation},
        )
    return _Completed(result.output, version, executable, (_PROGRAM, *arguments))


def _tail(text: str) -> str:
    return text[-_LOG_TAIL_CHARS:]


def _result(
    operation: str,
    run: _Completed,
    *,
    outputs: Mapping[str, Path],
    records: list[dict[str, object]] | None = None,
    **extra: object,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "backend": _PROGRAM,
        "operation": operation,
        "version": run.version,
        "command": list(run.command),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "log_tail": _tail(run.text),
        "provenance": {
            "tool": "SeqKit",
            "version": run.version,
            "tested_version": _TESTED_VERSION,
            "executable": run.executable,
            "citation": _CITATION,
        },
    }
    if records is not None:
        result["records"] = records
    result.update(extra)
    return result


# ---------------------------------------------------------------------------
# Parsing helpers (SeqKit writes tables to -o files or, for some modes, to stderr)
# ---------------------------------------------------------------------------


def _number(value: str) -> object:
    try:
        return int(value)
    except ValueError:
        pass
    try:
        number = float(value)
    except ValueError:
        return value
    return number if math.isfinite(number) else value


def _row(
    names: Sequence[str], fields: Sequence[str], text_keys: frozenset[str]
) -> dict[str, object]:
    return {
        name: value if name in text_keys else _number(value)
        for name, value in zip(names, fields, strict=False)
    }


def _file_rows(
    path: Path, names: Sequence[str], *, text_keys: frozenset[str] = frozenset(), header: bool
) -> list[dict[str, object]]:
    lines = [
        line for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line
    ]
    if header:
        if not lines:
            return []
        names = lines[0].rstrip("\t").split("\t")
        lines = lines[1:]
    return [_row(names, line.rstrip("\t").split("\t"), text_keys) for line in lines]


def _table_lines(text: str, first_column: str) -> list[str]:
    """The header line and the following rows of a tab-separated block SeqKit prints to stderr/stdout."""
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if "\t" in line and line.split("\t", 1)[0] == first_column:
            width = len(line.rstrip("\t").split("\t"))
            block = [line]
            for body in lines[index + 1 :]:
                if len(body.rstrip("\t").split("\t")) != width:
                    break
                block.append(body)
            return block
    return []


def _log_table(
    text: str, first_column: str, *, text_keys: frozenset[str] = frozenset()
) -> list[dict[str, object]]:
    """Parse a tab-separated block that SeqKit prints to stderr (or stdout) after a header line."""
    block = _table_lines(text, first_column)
    if not block:
        return []
    names = block[0].rstrip("\t").split("\t")
    return [_row(names, line.rstrip("\t").split("\t"), text_keys) for line in block[1:]]


def _record_count(path: Path) -> int:
    """Number of FASTA or (4-line) FASTQ records in a SeqKit output file."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline()
        if not first:
            return 0
        if first.startswith(">"):
            return 1 + sum(1 for line in handle if line.startswith(">"))
        if first.startswith("@"):
            return (1 + sum(1 for _ in handle) + 3) // 4
    return 0


def _directory_files(directory: Path) -> list[dict[str, object]]:
    return [
        {"file": path.name, "record_count": _record_count(path)}
        for path in sorted(directory.iterdir())
        if path.is_file()
    ]


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def seqkit_flanking_sequences(
    sequence_path: str,
    feature_path: str,
    output_path: str,
    *,
    feature_format: str = "bed",
    up_stream: int = 0,
    down_stream: int = 0,
    only_flank: bool = False,
    feature_types: Sequence[str] | None = None,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract BED or GTF features with up-/down-stream flanks (``seqkit subseq --bed/--gtf``).

    ``only_flank`` returns only the flank (SeqKit allows one of ``up_stream``/``down_stream`` in
    that mode). ``feature_types`` filters GTF features (SeqKit ``--feature``). SeqKit writes a
    ``<sequence>.seqkit.fai`` index next to the FASTA; it is reported in ``index_file``.
    """
    source = _input(sequence_path, "sequence_path")
    features = _input(feature_path, "feature_path")
    target = _output(output_path)
    fmt = _choice(feature_format, ("bed", "gtf"), "feature_format")
    up = _int(up_stream, "up_stream", minimum=0)
    down = _int(down_stream, "down_stream", minimum=0)
    flank_only = _flag(only_flank, "only_flank")
    if flank_only and (up > 0) == (down > 0):
        raise ConfigurationError(
            "only_flank needs exactly one of up_stream or down_stream above 0."
        )
    if feature_types is not None and fmt != "gtf":
        raise ConfigurationError("feature_types only applies to feature_format='gtf'.")
    options = [f"--{fmt}={features}"]
    if feature_types:
        options += [f"--feature={_text(name, 'feature_types')}" for name in feature_types]
    if up:
        options.append(f"--up-stream={up}")
    if down:
        options.append(f"--down-stream={down}")
    if flank_only:
        options.append("--only-flank")
    run = _run(
        "flanking_sequences_of_features",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    index = source.with_name(source.name + ".seqkit.fai")
    return _result(
        "flanking_sequences_of_features",
        run,
        outputs={"output": target},
        index_file=str(index) if index.exists() else None,
        record_count=_record_count(target),
    )


def seqkit_extended_stats(
    input_paths: Sequence[str],
    output_path: str,
    *,
    fastq_encoding: str = "sanger",
    n_values: Sequence[int] = (),
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extended statistics per file: length quartiles, N50 and Nx, Q20/Q30, AvgQual, GC (``seqkit stats -a -T``).

    ``n_values`` adds further Nx columns (0-100); N50/L50 are already part of ``-a``.
    """
    sources = _inputs(input_paths)
    target = _output(output_path)
    encoding = _choice(fastq_encoding, _FQ_ENCODINGS, "fastq_encoding")
    options = ["--all", "--tabular", f"--fq-encoding={encoding}"]
    for value in n_values:
        options.append(f"--N={_int(value, 'n_values', minimum=0, maximum=100)}")
    run = _run(
        "fastx_extended_statistics_n50_quality",
        options,
        sources,
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _file_rows(target, (), text_keys=frozenset({"file", "format", "type"}), header=True)
    return _result(
        "fastx_extended_statistics_n50_quality", run, outputs={"output": target}, records=rows
    )


def seqkit_sequence_digest(
    input_paths: Sequence[str],
    output_path: str,
    *,
    include_ids: bool = False,
    remove_gaps: bool = False,
    rna_to_dna: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Order-aware message digest of the sequences of each file (``seqkit sum``).

    The digest is ``seqkit.v0.1_...``; it is not a SHA-256 of the file. ``include_ids`` adds the IDs,
    ``remove_gaps`` strips gap characters, ``rna_to_dna`` converts U to T before hashing.
    """
    sources = _inputs(input_paths)
    target = _output(output_path)
    options: list[str] = []
    if _flag(include_ids, "include_ids"):
        options.append("--include-id")
    if _flag(remove_gaps, "remove_gaps"):
        options.append("--remove-gaps")
    if _flag(rna_to_dna, "rna_to_dna"):
        options.append("--rna2dna")
    run = _run(
        "sequence_set_digest",
        options,
        sources,
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _file_rows(
        target, ("digest", "file"), text_keys=frozenset({"digest", "file"}), header=False
    )
    return _result("sequence_set_digest", run, outputs={"output": target}, records=rows)


def seqkit_read_feature_histogram(
    input_path: str,
    field: str,
    *,
    bins: int = -1,
    log_transform: bool = False,
    quality_ascii_base: int = 33,
    validate_seq: bool = False,
    image_path: str | None = None,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of one read feature of a FASTA/FASTQ file, reported once at end of file (``seqkit watch -y``).

    ``field`` is one of ReadLen, MeanQual, GC, GCSkew. The bin table is parsed from SeqKit's
    stderr. ``image_path`` (a new file, with a PDF/PNG/... extension) also saves the plot.
    Online monitoring (``--print-freq``, ``--delay``, pass-through) is not offered: the report is
    one-shot, so the facade cannot follow a growing file.
    """
    source = _input(input_path)
    chosen = _choice(field, _READ_FIELDS, "field")
    image = _output(image_path, "image_path") if image_path is not None else None
    options = [
        f"--fields={chosen}",
        f"--bins={_int(bins, 'bins', minimum=-1)}",
        f"--qual-ascii-base={_int(quality_ascii_base, 'quality_ascii_base', minimum=0, maximum=127)}",
        "--dump",
        "--quiet-mode",
    ]
    if _flag(log_transform, "log_transform"):
        options.append("--log")
    if _flag(validate_seq, "validate_seq"):
        options.append("--validate-seq")
    if image is not None:
        options.append(f"--img={image}")
    run = _run(
        "read_feature_histogram_monitor",
        options,
        (source,),
        cwd=(image or source).parent,
        artifacts=(image,) if image is not None else (),
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _log_table(run.text, "Bin")
    outputs = {"image": image} if image is not None else {}
    return _result(
        "read_feature_histogram_monitor", run, outputs=outputs, records=rows, field=chosen
    )


def seqkit_fastq_records_for_fasta(
    fastq_path: str,
    fasta_path: str,
    output_path: str,
    *,
    only_positive_strand: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """FASTQ records that have a FASTA record with the same ID and the same sequence (``seqkit fa2fq``).

    SeqKit 2.14.0 matches the full sequence on both strands (reverse complement included) unless
    ``only_positive_strand`` is True. A FASTQ record whose ID is absent from the FASTA, or whose
    sequence differs, is not returned.
    """
    source = _input(fastq_path, "fastq_path")
    fasta = _input(fasta_path, "fasta_path")
    target = _output(output_path)
    options = [f"--fasta-file={fasta}"]
    if _flag(only_positive_strand, "only_positive_strand"):
        options.append("--only-positive-strand")
    run = _run(
        "retrieve_fastq_records_by_fasta",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "retrieve_fastq_records_by_fasta",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_mean_read_quality(
    fastq_path: str,
    output_path: str,
    *,
    quality_ascii_base: int = 33,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average quality of each FASTQ read, one row per read (``seqkit fx2tab --name --avg-qual``).

    ``quality_ascii_base`` 33 is Phred+33 (SeqKit default). The value is SeqKit's own average.
    """
    source = _input(fastq_path, "fastq_path")
    target = _output(output_path)
    options = [
        "--name",
        "--avg-qual",
        f"--qual-ascii-base={_int(quality_ascii_base, 'quality_ascii_base', minimum=0, maximum=127)}",
    ]
    run = _run(
        "mean_read_quality",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _file_rows(target, ("id", "avg_qual"), text_keys=frozenset({"id"}), header=False)
    return _result("mean_read_quality", run, outputs={"output": target}, records=rows)


def seqkit_fish_search(
    target_path: str,
    query_path: str,
    output_path: str,
    *,
    min_mapping_quality: float = 5.0,
    search_all: bool = False,
    stranded: bool = False,
    alignment_params: str | None = None,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate short query sequences in longer targets by local alignment (``seqkit fish``).

    ``alignment_params`` is "match,mismatch,gap_open,gap_extend" (SeqKit default 4,-4,-2,-1).
    ``search_all`` reports every hit, not only the best one per query; ``stranded`` searches the
    plus strand only. Hits below ``min_mapping_quality`` are dropped by SeqKit. SeqKit 2.14.0 ignores
    ``-o`` for fish, so this facade writes the table to ``output_path`` itself. SeqKit 2.14.0 panics
    (exit 2, reported as an error) when a target contains N, so use N-free targets.
    """
    target = _input(target_path, "target_path")
    query = _input(query_path, "query_path")
    out = _output(output_path)
    if not isinstance(min_mapping_quality, (int, float)) or isinstance(min_mapping_quality, bool):
        raise ConfigurationError("min_mapping_quality must be a number.")
    options = [f"--query-fastx={query}", f"--min-qual={float(min_mapping_quality)!r}"]
    if alignment_params is not None:
        if not _ALIGN_PARAMS.match(alignment_params):
            raise ConfigurationError("alignment_params must look like '4,-4,-2,-1'.")
        options.append(f"--aln-params={alignment_params}")
    if _flag(search_all, "search_all"):
        options.append("--all")
    if _flag(stranded, "stranded"):
        options.append("--stranded")
    # SeqKit 2.14.0 fish ignores -o and prints its table to stdout, so the table is written here.
    run = _run(
        "local_alignment_search_short_queries",
        options,
        (target,),
        cwd=out.parent,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    block = _table_lines(run.text, "Ref")
    out.write_text("".join(f"{line}\n" for line in block), encoding="utf-8")
    rows = _log_table(run.text, "Ref", text_keys=frozenset({"Ref", "Query", "Strand"}))
    return _result(
        "local_alignment_search_short_queries", run, outputs={"output": out}, records=rows
    )


def seqkit_grep_motif(
    input_path: str,
    output_path: str,
    *,
    motifs: Sequence[str],
    degenerate: bool = False,
    max_mismatch: int = 0,
    only_positive_strand: bool = False,
    ignore_case: bool = False,
    invert_match: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep records whose sequence contains any of the motifs (``seqkit grep --by-seq``).

    Literal matching by default; ``degenerate`` reads IUPAC codes. ``max_mismatch`` allows that many
    mismatches. Matches are on both strands unless ``only_positive_strand``. ``invert_match`` keeps the
    records without a match.
    """
    source = _input(input_path)
    target = _output(output_path)
    options = ["--by-seq"]
    options += [f"--pattern={pattern}" for pattern in _patterns(motifs)]
    if _flag(degenerate, "degenerate"):
        options.append("--degenerate")
    mismatch = _int(max_mismatch, "max_mismatch", minimum=0)
    if mismatch:
        options.append(f"--max-mismatch={mismatch}")
    if _flag(only_positive_strand, "only_positive_strand"):
        options.append("--only-positive-strand")
    if _flag(ignore_case, "ignore_case"):
        options.append("--ignore-case")
    if _flag(invert_match, "invert_match"):
        options.append("--invert-match")
    run = _run(
        "filter_records_by_sequence_motif",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "filter_records_by_sequence_motif",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_count_matching_records(
    input_path: str,
    *,
    patterns: Sequence[str],
    by_name: bool = False,
    by_sequence: bool = False,
    use_regexp: bool = False,
    invert_match: bool = False,
    ignore_case: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count records whose ID (default), full name, or sequence matches any pattern (``seqkit grep --count``).

    Patterns match IDs exactly unless ``use_regexp``. ``by_sequence`` searches sequences as motifs.
    With ``invert_match`` the count is of non-matching records. The count is read from stdout.
    """
    source = _input(input_path)
    options = ["--count"]
    options += [f"--pattern={pattern}" for pattern in _patterns(patterns)]
    if _flag(by_name, "by_name"):
        options.append("--by-name")
    if _flag(by_sequence, "by_sequence"):
        options.append("--by-seq")
    if _flag(use_regexp, "use_regexp"):
        options.append("--use-regexp")
    if _flag(invert_match, "invert_match"):
        options.append("--invert-match")
    if _flag(ignore_case, "ignore_case"):
        options.append("--ignore-case")
    run = _run(
        "count_matching_records",
        options,
        (source,),
        cwd=source.parent,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    numbers = [line.strip() for line in run.text.splitlines() if line.strip().isdigit()]
    if not numbers:
        raise BackendExecutionError(
            "SeqKit grep --count printed no count.",
            code="SEQKIT_OUTPUT_UNREADABLE",
            context={"log_tail": _tail(run.text)},
        )
    return _result("count_matching_records", run, outputs={}, count=int(numbers[-1]))


def seqkit_duplicate_records(
    input_path: str,
    output_path: str,
    *,
    times: int,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write every record ``times`` times in place (``seqkit duplicate -n``); copies keep their IDs."""
    source = _input(input_path)
    target = _output(output_path)
    count = _int(times, "times", minimum=1)
    run = _run(
        "duplicate_records",
        [f"--times={count}"],
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "duplicate_records", run, outputs={"output": target}, record_count=_record_count(target)
    )


def seqkit_range_records(
    input_path: str,
    output_path: str,
    *,
    start: int,
    end: int,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Records from position ``start`` to ``end`` (1-based, inclusive; negative counts from the end).

    ``seqkit range -r start:end``, e.g. (1, 12) is the first 12 records and (-3, -1) the last three.
    """
    source = _input(input_path)
    target = _output(output_path)
    first = _int(start, "start")
    last = _int(end, "end")
    run = _run(
        "select_records_by_position",
        [f"--range={first}:{last}"],
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "select_records_by_position",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_head_genome(
    input_path: str,
    output_path: str,
    *,
    minimal_common_words: int = 1,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Records of the first genome, where a genome is a run of records sharing name prefix words (``seqkit head-genome``).

    ``minimal_common_words`` is the number of shared leading words that define one genome.
    """
    source = _input(input_path)
    target = _output(output_path)
    words = _int(minimal_common_words, "minimal_common_words", minimum=1)
    run = _run(
        "first_genome_by_name_prefix",
        [f"--mini-common-words={words}"],
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "first_genome_by_name_prefix",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_split_by_region(
    input_path: str,
    output_dir: str,
    *,
    region: str,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split records into files by the subsequence of a region, e.g. a barcode (``seqkit split -r``).

    ``region`` uses SeqKit syntax: "1:5" for the first five bases, "-12:-1" for the last twelve.
    ``output_dir`` must not exist; SeqKit creates it and one file per distinct subsequence.
    """
    source = _input(input_path)
    directory = _new_directory(output_dir)
    if not isinstance(region, str) or not _REGION.match(region):
        raise ConfigurationError("region must look like '1:5' or '-12:-1'.")
    run = _run(
        "split_by_region_subsequence",
        [f"--by-region={region}", f"--out-dir={directory}"],
        (source,),
        cwd=directory.parent,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    if not directory.is_dir():
        raise BackendExecutionError(
            "SeqKit split did not create the output directory.",
            code="SEQKIT_OUTPUT_MISSING",
            context={"operation": "split_by_region_subsequence"},
        )
    return _result(
        "split_by_region_subsequence",
        run,
        outputs={"directory": directory},
        records=_directory_files(directory),
    )


def seqkit_concat_by_id(
    input_paths: Sequence[str],
    output_path: str,
    *,
    full_join: bool = False,
    fill: str | None = None,
    separator: str = "|",
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Concatenate sequences that share an ID across files, in file order (``seqkit concat``).

    ``full_join`` keeps IDs missing from some files and fills the gaps with ``fill`` (letters only).
    ``separator`` joins the descriptions of the merged records.
    """
    sources = _inputs(input_paths)
    target = _output(output_path)
    options: list[str] = []
    if _flag(full_join, "full_join"):
        options.append("--full")
    if fill is not None:
        if not full_join:
            raise ConfigurationError("fill only applies with full_join=True.")
        if not isinstance(fill, str) or not _FILL.match(fill):
            raise ConfigurationError("fill must be letters or '-'.")
        options.append(f"--fill={fill}")
    if not isinstance(separator, str):
        raise ConfigurationError("separator must be a string.")
    options.append(f"--separator={separator}")
    run = _run(
        "concatenate_same_id_sequences",
        options,
        sources,
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "concatenate_same_id_sequences",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_rename_duplicate_ids(
    input_path: str,
    output_path: str,
    *,
    by_name: bool = False,
    start_number: int = 2,
    separator: str = "_",
    rename_first: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Suffix repeated IDs with a counter, e.g. id, id_2, id_3 (``seqkit rename``).

    The first occurrence keeps its ID unless ``rename_first``. ``by_name`` checks full headers
    instead of IDs; ``start_number`` is the first counter value (SeqKit default 2).
    """
    source = _input(input_path)
    target = _output(output_path)
    options = [
        f"--start-num={_int(start_number, 'start_number', minimum=1)}",
        f"--separator={_text(separator, 'separator')}",
    ]
    if _flag(by_name, "by_name"):
        options.append("--by-name")
    if _flag(rename_first, "rename_first"):
        options.append("--rename-1st-rec")
    run = _run(
        "make_duplicate_ids_unique",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "make_duplicate_ids_unique",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_replace_sequence_regex(
    input_path: str,
    output_path: str,
    *,
    pattern: str,
    replacement: str,
    ignore_case: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Regex substitution applied to the sequences (``seqkit replace --by-seq``).

    ``replacement`` may use capture groups ($1). FASTQ quality strings must keep their length, so
    the replacement must not lengthen a FASTQ record.
    """
    source = _input(input_path)
    target = _output(output_path)
    if not isinstance(replacement, str):
        raise ConfigurationError("replacement must be a string.")
    options = [
        "--by-seq",
        f"--pattern={_text(pattern, 'pattern')}",
        f"--replacement={replacement}",
    ]
    if _flag(ignore_case, "ignore_case"):
        options.append("--ignore-case")
    run = _run(
        "replace_sequence_regex",
        options,
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "replace_sequence_regex",
        run,
        outputs={"output": target},
        record_count=_record_count(target),
    )


def seqkit_bam_alignment_stats(
    bam_path: str,
    *,
    index_based: bool = False,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mapped, unmapped, primary, secondary and supplementary counts of a BAM (``seqkit bam -s``).

    ``index_based`` uses the BAM index (``-i``, needs a ``.bai``) instead of scanning all records.
    The table is parsed from SeqKit's stderr and includes the multimapping percentage in full-scan mode.
    """
    source = _input(bam_path, "bam_path")
    mode = "--idx-stat" if _flag(index_based, "index_based") else "--stat"
    run = _run(
        "bam_alignment_statistics",
        [mode, "--quiet-mode"],
        (source,),
        cwd=source.parent,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    first = "AlnPerc" if mode == "--idx-stat" else "PrimAlnPerc"
    rows = _log_table(run.text, first, text_keys=frozenset({"File"}))
    return _result(
        "bam_alignment_statistics", run, outputs={}, records=rows, index_based=mode == "--idx-stat"
    )


def seqkit_bam_feature_histogram(
    bam_path: str,
    field: str,
    *,
    bins: int = -1,
    log_transform: bool = False,
    image_path: str | None = None,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of one per-record BAM feature, e.g. ReadLen, MapQual or MeanQual (``seqkit bam -f -y``).

    Known features are listed by ``seqkit bam -H``. Some features need tags the BAM may not carry;
    for example Acc panics in SeqKit 2.14.0 when the BAM has no NM tag, so the facade does not pass it.
    ``image_path`` (a new file) also saves the plot.
    """
    source = _input(bam_path, "bam_path")
    if not isinstance(field, str) or not _BAM_FIELD.match(field):
        raise ConfigurationError("field must be a SeqKit BAM feature name such as ReadLen.")
    image = _output(image_path, "image_path") if image_path is not None else None
    options = [
        f"--field={field}",
        f"--bins={_int(bins, 'bins', minimum=-1)}",
        "--dump",
        "--quiet-mode",
    ]
    if _flag(log_transform, "log_transform"):
        options.append("--log")
    if image is not None:
        options.append(f"--img={image}")
    run = _run(
        "bam_feature_histogram",
        options,
        (source,),
        cwd=(image or source).parent,
        artifacts=(image,) if image is not None else (),
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _log_table(run.text, "Bin")
    outputs = {"image": image} if image is not None else {}
    return _result("bam_feature_histogram", run, outputs=outputs, records=rows, field=field)


def seqkit_bam_top_records(
    bam_path: str,
    output_bam: str,
    *,
    field: str,
    top_size: int = 100,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Save the ``top_size`` BAM records with the largest value of a feature to a new BAM (``seqkit bam -@ -?``)."""
    source = _input(bam_path, "bam_path")
    target = _output(output_bam, "output_bam")
    if not isinstance(field, str) or not _BAM_FIELD.match(field):
        raise ConfigurationError("field must be a SeqKit BAM feature name such as ReadLen.")
    size = _int(top_size, "top_size", minimum=1)
    run = _run(
        "bam_top_records_by_field",
        [f"--field={field}", f"--top-bam={target}", f"--top-size={size}", "--quiet-mode"],
        (source,),
        cwd=target.parent,
        artifacts=(target,),
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    return _result(
        "bam_top_records_by_field", run, outputs={"bam": target}, field=field, top_size=size
    )


def seqkit_bam_bundles(
    bam_path: str,
    output_dir: str,
    *,
    min_reads: int = -1,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Partition a BAM into loci (``min_reads=-1``) or bundles of at least ``min_reads`` reads (``seqkit bam -N``).

    ``output_dir`` must not exist; SeqKit writes ``<bam name>_bundles/`` into it, and the bundle table is parsed.
    """
    source = _input(bam_path, "bam_path")
    directory = _new_directory(output_dir)
    bundle = _int(min_reads, "min_reads", minimum=-1)
    directory.mkdir()
    run = _run(
        "bam_bundle_partition",
        [f"--bundle={bundle}", "--quiet-mode"],
        (source,),
        cwd=directory,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _log_table(run.text, "Bundle", text_keys=frozenset({"Chrom"}))
    summary = _BUNDLE_SUMMARY.search(run.text)
    counts: dict[str, object] = {}
    if summary is not None:
        counts = {
            "record_count": int(summary.group(1)),
            "locus_count": int(summary.group(2)),
            "bundle_count": int(summary.group(3)),
        }
    bundle_dir = directory / f"{source.name}_bundles"
    return _result(
        "bam_bundle_partition",
        run,
        outputs={"directory": directory, "bundles": bundle_dir},
        records=rows,
        min_reads=bundle,
        **counts,
    )


def seqkit_merge_slides_to_bed(
    window_names_path: str,
    output_path: str,
    *,
    min_overlap: int = 1,
    max_gap: int = 0,
    seqkit_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge overlapping sliding-window IDs into BED intervals (``seqkit merge-slides``).

    The input is a text file with one window ID per line, such as ``seqkit seq -n`` prints for the
    output of ``seqkit sliding`` (IDs like ``chr1_sliding:1-1000``). Output is BED: chromosome, 0-based
    start, end. ``min_overlap`` is the minimum overlap of adjacent windows (SeqKit default 1);
    ``max_gap`` limits the distance between window starts (0 means no limit).
    """
    source = _input(window_names_path, "window_names_path")
    target = _output(output_path)
    run = _run(
        "merge_sliding_windows_to_bed",
        [
            f"--min-overlap={_int(min_overlap, 'min_overlap', minimum=1)}",
            f"--max-gap={_int(max_gap, 'max_gap', minimum=0)}",
        ],
        (source,),
        cwd=target.parent,
        output=target,
        seqkit_path=seqkit_path,
        verbose=verbose,
    )
    rows = _file_rows(
        target, ("chrom", "start", "end"), text_keys=frozenset({"chrom"}), header=False
    )
    return _result("merge_sliding_windows_to_bed", run, outputs={"output": target}, records=rows)


__all__ = [
    "seqkit_bam_alignment_stats",
    "seqkit_bam_bundles",
    "seqkit_bam_feature_histogram",
    "seqkit_bam_top_records",
    "seqkit_concat_by_id",
    "seqkit_count_matching_records",
    "seqkit_duplicate_records",
    "seqkit_extended_stats",
    "seqkit_fastq_records_for_fasta",
    "seqkit_fish_search",
    "seqkit_flanking_sequences",
    "seqkit_grep_motif",
    "seqkit_head_genome",
    "seqkit_mean_read_quality",
    "seqkit_merge_slides_to_bed",
    "seqkit_range_records",
    "seqkit_read_feature_histogram",
    "seqkit_rename_duplicate_ids",
    "seqkit_replace_sequence_regex",
    "seqkit_sequence_digest",
    "seqkit_split_by_region",
]

_register_public_verbose(globals(), tuple(__all__))
