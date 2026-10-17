"""Generic adapters around external command-line programs.

Programs are always started from an argument *list* (never through a shell), with a timeout, after being resolved
with ``shutil.which``.  When the requested program is not installed, ``run`` returns ``status: "not_implemented"``
naming the missing program.  The conversion and parsing helpers are pure text transformations.

A *command* is given as a program name or argv list (``["samtools", "view", "-h", "x.bam"]``), a command-line string
(split with ``shlex``, no shell involved) or a mapping ``{"program": str, "args": [...], "flags": {"--threads": 4,
"-v": True, "--exclude": ["a", "b"]}}``.  ``options["subcommand"]`` (a string or a list) is inserted right after
the program, before any flag, for programs with sub-commands (``samtools view``).
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import BackendTimeoutError, ConfigurationError, InputFormatError
from dnakit.sequence import normalize_dna_sequence

_MODULE = "external_adapters"
_MAX_TEXT = 50_000_000
_MAX_LOG_LINES = 200_000


def _opts(options: Mapping[str, object] | None) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping.", code="INVALID_OPTION")
    return dict(options or {})


def _token(value: object, what: str) -> str:
    if isinstance(value, bool) or not isinstance(value, (str, int, float, os.PathLike)):
        raise InputFormatError(f"{what} must be a string, number or path, got {type(value).__name__}.", code="INVALID_COMMAND")
    text = os.fspath(value) if isinstance(value, os.PathLike) else str(value)
    if "\x00" in text:
        raise InputFormatError(f"{what} contains a NUL byte.", code="INVALID_COMMAND")
    return text


def _flag_tokens(flags: object) -> list[str]:
    if not isinstance(flags, Mapping):
        raise InputFormatError("flags must be a mapping of flag -> value.", code="INVALID_COMMAND")
    tokens: list[str] = []
    for flag, value in flags.items():
        name = _token(flag, "flag")
        if not name.startswith("-") or name in {"-", "--"}:
            raise InputFormatError(f"Flag {name!r} must start with '-'.", code="INVALID_COMMAND")
        if value is None or value is False:
            continue
        if value is True:
            tokens.append(name)
        elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            for item in value:
                tokens.extend((name, _token(item, f"value of {name}")))
        else:
            tokens.extend((name, _token(value, f"value of {name}")))
    return tokens


def _subcommand_tokens(settings: Mapping[str, object]) -> list[str]:
    """Tokens of ``options['subcommand']`` (one string or a list of strings); empty when absent."""

    sub = settings.get("subcommand")
    if sub is None:
        return []
    items = [sub] if isinstance(sub, (str, os.PathLike)) else list(sub) if isinstance(sub, Sequence) and not isinstance(sub, bytes) else None  # type: ignore[call-overload]
    if items is None or not items:
        raise ConfigurationError("options.subcommand must be a non-empty string or list of strings.", code="INVALID_OPTION", hint="For example 'view' or ['index', 'build'].")
    tokens = [_token(item, "subcommand") for item in items]
    if any(not token.strip() for token in tokens):
        raise ConfigurationError("options.subcommand must not contain empty strings.", code="INVALID_OPTION", hint="For example 'view' or ['index', 'build'].")
    return tokens


def _argv(command: object, settings: Mapping[str, object]) -> list[str]:
    """Build the argv list from a program name, argv list, command-line string or mapping (+ options flags/args).

    Order: program, ``options['subcommand']``, the mapping's flags, the mapping's args, ``options['flags']``,
    ``options['args']``.  Without a subcommand this is the previous order.
    """

    extra_flags, extra_args = settings.get("flags"), settings.get("args")
    subcommand = _subcommand_tokens(settings)
    if isinstance(command, Mapping):
        if "program" not in command:
            raise InputFormatError("A command mapping needs a 'program' key.", code="INVALID_COMMAND")
        argv = [_token(command["program"], "program"), *subcommand]
        argv += _flag_tokens(command.get("flags", {}))
        argv += [_token(item, "argument") for item in command.get("args", [])]
    elif isinstance(command, str):
        try:
            argv = shlex.split(command)
        except ValueError as exc:
            raise InputFormatError(f"Cannot split the command line: {exc}", code="INVALID_COMMAND") from exc
    elif isinstance(command, Sequence) and not isinstance(command, bytes):
        argv = [_token(item, "argument") for item in command]
    elif isinstance(command, os.PathLike):
        argv = [_token(command, "program")]
    else:
        raise InputFormatError("input must be a program name, argv list, command line or mapping.", code="INVALID_COMMAND")
    if subcommand and not isinstance(command, Mapping) and argv:
        argv[1:1] = subcommand
    if extra_flags:
        argv += _flag_tokens(extra_flags)
    if extra_args:
        argv += [_token(item, "argument") for item in extra_args]  # type: ignore[union-attr]
    if not argv or not argv[0].strip():
        raise InputFormatError("The command needs a program name.", code="INVALID_COMMAND")
    return argv


def _text_of(value: object) -> str:
    """Text from a string, a file path, or a mapping holding ``stdout``/``text``/``path``."""

    if isinstance(value, Mapping):
        for key in ("stdout", "text", "output"):
            if isinstance(value.get(key), str):
                return str(value[key])
        if "path" in value:
            return _text_of(Path(str(value["path"])))
        if isinstance(value.get("result"), Mapping):
            return _text_of(value["result"])
        raise InputFormatError("The mapping needs a 'stdout', 'text' or 'path' key.", code="INVALID_INPUT")
    if isinstance(value, os.PathLike):
        path = Path(value)
        if not path.is_file():
            raise InputFormatError(f"File not found: {path}", code="FILE_NOT_FOUND")
        if path.stat().st_size > _MAX_TEXT:
            raise InputFormatError(f"File exceeds {_MAX_TEXT} bytes.", code="FILE_TOO_LARGE")
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        if len(value) > _MAX_TEXT:
            raise InputFormatError("Text is too large.", code="FILE_TOO_LARGE")
        return value
    raise InputFormatError(f"Unsupported input type {type(value).__name__}.", code="INVALID_INPUT")


def _convert(cell: str) -> object:
    for kind in (int, float):
        try:
            return kind(cell)
        except ValueError:
            continue
    return cell


def _parse_records(text: str, settings: Mapping[str, object]) -> tuple[str, list[object]]:
    """Parse tool output into ``(format, records)``; formats: tsv, csv, json, jsonl, fasta, lines."""

    fmt = str(settings.get("format", "auto"))
    stripped = text.lstrip()
    if fmt == "auto":
        if stripped.startswith(">"):
            fmt = "fasta"
        elif stripped[:1] in {"{", "["}:
            fmt = "json"
        elif "\t" in text:
            fmt = "tsv"
        else:
            fmt = "lines"
    if fmt not in {"tsv", "csv", "json", "jsonl", "fasta", "lines"}:
        raise ConfigurationError("options.format must be auto, tsv, csv, json, jsonl, fasta or lines.", code="INVALID_OPTION")
    if fmt == "lines":
        return fmt, [line for line in text.splitlines() if line.strip()]
    if fmt == "json":
        try:
            value = json.loads(text)
        except ValueError as exc:
            raise InputFormatError(f"Invalid JSON output: {exc}", code="INVALID_OUTPUT") from exc
        return fmt, value if isinstance(value, list) else [value]
    if fmt == "jsonl":
        try:
            return fmt, [json.loads(line) for line in text.splitlines() if line.strip()]
        except ValueError as exc:
            raise InputFormatError(f"Invalid JSON Lines output: {exc}", code="INVALID_OUTPUT") from exc
    if fmt == "fasta":
        records: list[object] = []
        current: dict[str, object] | None = None
        chunks: list[str] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    current["sequence"] = "".join(chunks)
                    records.append(current)
                header = line[1:].strip()
                identifier, _, description = header.partition(" ")
                current, chunks = {"id": identifier, "description": description.strip()}, []
            elif current is None:
                raise InputFormatError("FASTA output has sequence before the first header.", code="INVALID_OUTPUT")
            else:
                chunks.append(line)
        if current is not None:
            current["sequence"] = "".join(chunks)
            records.append(current)
        return fmt, records
    delimiter = str(settings.get("delimiter", "\t" if fmt == "tsv" else ","))
    prefix = settings.get("comment_prefix", "#")
    lines = [line for line in text.splitlines() if line.strip() and not (prefix and line.startswith(str(prefix)))]
    rows = list(csv.reader(lines, delimiter=delimiter))
    convert = settings.get("convert_types", True) is not False
    columns = settings.get("columns")
    if columns is not None:
        names = [str(name) for name in columns]  # type: ignore[union-attr]
        body = rows
    elif settings.get("header", True) is not False:
        names, body = (rows[0], rows[1:]) if rows else ([], [])
    else:
        names, body = [], rows
    records = []
    for index, row in enumerate(body, start=1):
        keys = names or [f"column_{i}" for i in range(1, len(row) + 1)]
        if names and len(row) != len(names):
            raise InputFormatError(f"Row {index} has {len(row)} fields, expected {len(names)}.", code="INVALID_OUTPUT")
        records.append({key: (_convert(cell) if convert else cell) for key, cell in zip(keys, row)})
    return fmt, records


def build_command(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build an argv list (never a shell string) for an external program from a program, flags and arguments.

    ``input``: a program name, an argv list, a command-line string (split with ``shlex``) or a mapping
    ``{"program", "args": [...], "flags": {...}}``.  Flags are emitted in mapping order: ``true`` gives the bare flag,
    ``false``/``null`` omits it, a list repeats the flag per item, anything else gives ``flag value``; flag names must
    start with ``-``.  Options ``flags`` and ``args`` are appended after the input's own.  The mapping form emits
    ``program, flags, args``, which is wrong for programs with sub-commands (``samtools -b view x.bam``): give the
    sub-command in ``options['subcommand']`` (a string such as ``"view"`` or a list such as ``["index", "build"]``); it is
    emitted right after the program and before the flags, so ``{"program": "samtools", "args": ["x.bam"], "flags":
    {"-b": True}}`` with ``subcommand="view"`` gives ``samtools view -b x.bam``.  Without it the order is unchanged.  The
    subcommand is inserted after the first token for every input form (also an argv list or command-line string).
    Result: ``argv``, ``program`` and ``command_line`` (``shlex.join``, for display only).  Nothing is executed.
    """

    settings = _opts(options)
    argv = _argv(input, settings)
    return candidate_ok(_MODULE, "build_command", input, {"argv": argv, "program": argv[0], "command_line": shlex.join(argv)}, options=options)


def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Check whether external programs are installed (``shutil.which``), optionally reading their version.

    ``input``: a program name/path, or a list of them.  Options: ``version_args`` (list, e.g. ``["--version"]``;
    when given, ``[program, *version_args]`` is run with ``timeout`` seconds, default 10, and the first line of its
    combined output is returned as ``version``; by default nothing is executed).  A missing program is reported with
    ``available: false`` (the check itself succeeds).  Result: ``all_available`` and one row per program.
    """

    settings = _opts(options)
    programs = [input] if isinstance(input, (str, os.PathLike)) else list(input) if isinstance(input, Sequence) else None  # type: ignore[call-overload]
    if not programs:
        raise InputFormatError("input must be a program name or a list of program names.", code="INVALID_INPUT")
    version_args = settings.get("version_args")
    timeout = float(settings.get("timeout", 10))
    rows = []
    for program in programs:
        name = _token(program, "program")
        path = shutil.which(name)
        row: dict[str, object] = {"program": name, "available": path is not None, "path": path, "version": None}
        if path is not None and version_args is not None:
            try:
                completed = subprocess.run(
                    [path, *[_token(item, "version argument") for item in version_args]],  # type: ignore[union-attr]
                    capture_output=True, text=True, timeout=timeout, check=False, stdin=subprocess.DEVNULL,
                )
                lines = (completed.stdout + completed.stderr).strip().splitlines()
                row["version"] = lines[0] if lines else ""
            except subprocess.TimeoutExpired:
                row["version"] = None
                row["version_error"] = f"timed out after {timeout} s"
        rows.append(row)
    return candidate_ok(_MODULE, "check", input, {"all_available": all(r["available"] for r in rows), "programs": rows}, options=options)


def convert_input(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert DNA sequence records into the text format an external tool expects (FASTA, FASTQ or TSV).

    ``input``: FASTA text, a list of sequences, ``(id, sequence)`` pairs or mappings with ``id``/``sequence``.
    Sequences are normalized (uppercase, ``U`` to ``T``, IUPAC letters only).  Options: ``target_format``
    ("fasta" default, "fastq", "tsv"), ``line_width`` (FASTA wrapping, default 80, 0 = no wrapping), ``quality``
    (single character used for every base in FASTQ, default "I").  Records without an id are named ``seq1``, ``seq2``...
    Result: ``format``, ``text`` and ``record_count``.
    """

    settings = _opts(options)
    target = str(settings.get("target_format", "fasta"))
    if target not in {"fasta", "fastq", "tsv"}:
        raise ConfigurationError("options.target_format must be fasta, fastq or tsv.", code="INVALID_OPTION")
    width = int(settings.get("line_width", 80))
    quality = str(settings.get("quality", "I"))
    if width < 0 or len(quality) != 1:
        raise ConfigurationError("line_width must be >= 0 and quality a single character.", code="INVALID_OPTION")
    records: list[tuple[str, str]] = []

    def add(identifier: object, sequence: object) -> None:
        try:
            normalized = normalize_dna_sequence(str(sequence))
        except ValueError as exc:
            raise InputFormatError(str(exc), code="INVALID_SEQUENCE") from exc
        records.append((str(identifier) if identifier not in (None, "") else f"seq{len(records) + 1}", normalized))

    if isinstance(input, (str, os.PathLike)) and (isinstance(input, os.PathLike) or str(input).lstrip().startswith(">")):
        _, parsed = _parse_records(_text_of(input), {"format": "fasta"})
        for record in parsed:
            add(record["id"], record["sequence"])  # type: ignore[index]
    elif isinstance(input, str):
        add(None, input)
    elif isinstance(input, Sequence):
        for item in input:
            if isinstance(item, Mapping):
                add(item.get("id"), item.get("sequence"))
            elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                add(item[0], item[1])
            else:
                add(None, item)
    else:
        raise InputFormatError("input must be FASTA text, a sequence or a list of records.", code="INVALID_INPUT")
    if not records:
        raise InputFormatError("No sequence records were found.", code="INVALID_INPUT")
    lines: list[str] = []
    for identifier, sequence in records:
        if target == "fasta":
            lines.append(f">{identifier}")
            lines.extend([sequence[i:i + width] for i in range(0, len(sequence), width)] if width else [sequence])
        elif target == "fastq":
            lines.extend((f"@{identifier}", sequence, "+", quality * len(sequence)))
        else:
            lines.append(f"{identifier}\t{sequence}")
    return candidate_ok(_MODULE, "convert_input", input, {"format": target, "text": "\n".join(lines) + "\n", "record_count": len(records)}, options=options)


def convert_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert external-tool output (TSV, CSV, JSON, JSON Lines, FASTA, lines) into JSON, JSON Lines, TSV or CSV text.

    ``input``: output text, a file path (``pathlib.Path``) or a run result mapping with ``stdout``.  The source is
    parsed as in ``parse_output`` (same options: ``format``, ``delimiter``, ``header``, ``columns``,
    ``comment_prefix``, ``convert_types``) and written in ``target_format`` ("json" default, "jsonl", "tsv", "csv").
    Tabular targets need dictionary records; their columns are the union of keys in first-seen order.
    Result: ``format``, ``text``, ``record_count``.
    """

    settings = _opts(options)
    target = str(settings.get("target_format", "json"))
    if target not in {"json", "jsonl", "tsv", "csv"}:
        raise ConfigurationError("options.target_format must be json, jsonl, tsv or csv.", code="INVALID_OPTION")
    _, records = _parse_records(_text_of(input), settings)
    if target == "json":
        text = json.dumps(records, indent=2, ensure_ascii=False) + "\n"
    elif target == "jsonl":
        text = "".join(json.dumps(record, ensure_ascii=False) + "\n" for record in records)
    else:
        if any(not isinstance(record, Mapping) for record in records):
            raise InputFormatError("Tabular targets need dictionary records.", code="INVALID_OUTPUT")
        columns = list(dict.fromkeys(key for record in records for key in record))  # type: ignore[union-attr]
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t" if target == "tsv" else ",", lineterminator="\n")
        writer.writerow(columns)
        for record in records:
            writer.writerow(["" if record.get(key) is None else record[key] for key in columns])  # type: ignore[union-attr]
        text = buffer.getvalue()
    return candidate_ok(_MODULE, "convert_output", input, {"format": target, "text": text, "record_count": len(records)}, options=options)


_LEVELS = (
    ("error", re.compile(r"\b(?:fatal|critical|error|exception|traceback|segmentation fault)\b", re.IGNORECASE)),
    ("warning", re.compile(r"\b(?:warn|warning)\b", re.IGNORECASE)),
    ("info", re.compile(r"\binfo\b", re.IGNORECASE)),
    ("debug", re.compile(r"\bdebug\b", re.IGNORECASE)),
)


def make_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarize one external-program run (a result of ``run``) as a Markdown or JSON report.

    ``input``: the ``run`` result mapping (``argv``, ``returncode``, ``stdout``, ``stderr``, ``duration_s``), or the
    full ``run`` envelope.  The report lists the command, exit code, duration, ``parse_log`` counts for stderr and
    the first ``max_lines`` (default 20) lines of stdout and stderr.  Options: ``format`` ("markdown" default or
    "json"), ``title`` (default "External program run").
    """

    settings = _opts(options)
    run = input["result"] if isinstance(input, Mapping) and isinstance(input.get("result"), Mapping) else input
    if not isinstance(run, Mapping) or "argv" not in run or "returncode" not in run:
        raise InputFormatError("input must be a run result with 'argv' and 'returncode'.", code="INVALID_INPUT")
    output_format = str(settings.get("format", "markdown"))
    if output_format not in {"markdown", "json"}:
        raise ConfigurationError("options.format must be markdown or json.", code="INVALID_OPTION")
    max_lines = int(settings.get("max_lines", 20))
    stdout, stderr = str(run.get("stdout", "")), str(run.get("stderr", ""))
    log = _log_summary(stderr + "\n" + stdout if stderr else stdout)
    summary = {
        "command": shlex.join([str(item) for item in run["argv"]]), "returncode": run["returncode"],
        "success": run["returncode"] == 0, "duration_s": run.get("duration_s"), "log": log,
        "stdout_head": stdout.splitlines()[:max_lines], "stderr_head": stderr.splitlines()[:max_lines],
    }
    if output_format == "json":
        return candidate_ok(_MODULE, "make_report", input, {"format": "json", "report": summary}, options=options)
    lines = [f"# {settings.get('title', 'External program run')}", "", f"- command: `{summary['command']}`", f"- return code: {summary['returncode']}",
             f"- success: {'yes' if summary['success'] else 'no'}", f"- duration: {summary['duration_s']} s",
             f"- log: {log['levels']['error']} error, {log['levels']['warning']} warning line(s)", "", "## stdout (head)", ""]
    lines += ["    " + line for line in summary["stdout_head"]] or ["(empty)"]
    lines += ["", "## stderr (head)", ""] + (["    " + line for line in summary["stderr_head"]] or ["(empty)"])
    return candidate_ok(_MODULE, "make_report", input, {"format": "markdown", "report": "\n".join(lines) + "\n"}, options=options)


def _log_summary(text: str, *, sample: int = 100) -> dict[str, object]:
    levels = {name: 0 for name, _ in _LEVELS}
    samples: dict[str, list[dict[str, object]]] = {"error": [], "warning": []}
    count = 0
    for number, line in enumerate(text.splitlines(), start=1):
        if number > _MAX_LOG_LINES:
            raise InputFormatError(f"Log exceeds {_MAX_LOG_LINES} lines.", code="FILE_TOO_LARGE")
        count = number
        for name, pattern in _LEVELS:  # the first (most severe) matching level classifies the line
            if pattern.search(line):
                levels[name] += 1
                if name in samples and len(samples[name]) < sample:
                    samples[name].append({"line": number, "text": line.strip()[:500]})
                break
    return {"line_count": count, "levels": levels, "errors": samples["error"], "warnings": samples["warning"], "has_errors": levels["error"] > 0}


def parse_log(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Classify the lines of an external-program log by severity (heuristic keyword rule).

    ``input``: log text, a file path (``pathlib.Path``) or a run result mapping (stderr then stdout are scanned
    if ``stderr`` is present, else ``stdout``).  A line is classified by the first whole word that matches, in order:
    error (fatal, critical, error, exception, traceback, segmentation fault), warning (warn, warning), info, debug,
    case-insensitively; lines without keyword are unclassified.  Result: ``line_count``, ``levels`` (counts),
    ``errors``/``warnings`` (up to 100 sample ``{line, text}``) and ``has_errors``.
    """

    _opts(options)
    if isinstance(input, Mapping) and "stderr" in input:
        text = str(input.get("stderr", "")) + "\n" + str(input.get("stdout", ""))
    else:
        text = _text_of(input)
    return candidate_ok(_MODULE, "parse_log", input, _log_summary(text), options=options)


def parse_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Parse external-tool output into structured records.

    ``input``: output text, a file path (``pathlib.Path``) or a run result mapping with ``stdout``.  Options:
    ``format`` ("auto" default: FASTA if it starts with ``>``, JSON if it starts with ``{``/``[``, TSV if it contains
    tabs, else lines; or tsv, csv, json, jsonl, fasta, lines), ``delimiter``, ``header`` (default true: first row names
    the columns), ``columns`` (explicit column names, used instead of a header), ``comment_prefix`` (default "#";
    such lines are skipped in TSV/CSV), ``convert_types`` (default true: integer/float cells become numbers).
    FASTA records are ``{id, description, sequence}``.  Result: ``format``, ``records``, ``record_count``.
    """

    settings = _opts(options)
    fmt, records = _parse_records(_text_of(input), settings)
    return candidate_ok(_MODULE, "parse_output", input, {"format": fmt, "records": records, "record_count": len(records)}, options=options)


def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run an external program (argv list, no shell, with a timeout) and capture its output.

    ``input``: a command as described in the module docstring.  The program is resolved with ``shutil.which``; if it
    is not installed the result is ``status: "not_implemented"`` naming it.  Options: ``timeout`` seconds (default 60,
    maximum 3600; on expiry the process is killed and ``BackendTimeoutError`` is raised), ``cwd`` (existing directory),
    ``env`` (mapping of extra environment variables), ``stdin`` (text sent to the program), ``max_output_bytes``
    (default 1,000,000; longer stdout/stderr are truncated and flagged), ``flags``/``args`` (appended) and
    ``subcommand`` (inserted right after the program; the argv is built exactly as in ``build_command``).  A non-zero exit code is reported in the result (``returncode``, ``success``), not raised.
    Result: ``argv`` (with the resolved program path), ``returncode``, ``success``, ``stdout``, ``stderr``,
    ``duration_s`` and ``truncated``.
    """

    settings = _opts(options)
    argv = _argv(input, settings)
    timeout = float(settings.get("timeout", 60))
    if not 0 < timeout <= 3600:
        raise ConfigurationError("options.timeout must be in (0, 3600] seconds.", code="INVALID_OPTION")
    cap = int(settings.get("max_output_bytes", 1_000_000))
    if cap < 1:
        raise ConfigurationError("options.max_output_bytes must be positive.", code="INVALID_OPTION")
    program = shutil.which(argv[0])
    if program is None:
        return candidate_not_implemented(
            _MODULE, "run", input, options=options, requires=f"the external program {argv[0]!r} on PATH",
            message=f"The external program {argv[0]!r} was not found on PATH, so it was not run.",
        )
    cwd = settings.get("cwd")
    if cwd is not None and not Path(str(cwd)).is_dir():
        raise ConfigurationError(f"options.cwd is not a directory: {cwd}", code="INVALID_OPTION")
    environment = None
    if settings.get("env") is not None:
        extra = settings["env"]
        if not isinstance(extra, Mapping):
            raise ConfigurationError("options.env must be a mapping.", code="INVALID_OPTION")
        environment = {**os.environ, **{_token(k, "env name"): _token(v, "env value") for k, v in extra.items()}}
    stdin_text = settings.get("stdin")
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [program, *argv[1:]], capture_output=True, timeout=timeout, check=False, cwd=cwd, env=environment,
            input=None if stdin_text is None else str(stdin_text).encode("utf-8"),
            stdin=subprocess.DEVNULL if stdin_text is None else None,
        )
    except subprocess.TimeoutExpired as exc:
        raise BackendTimeoutError(f"{argv[0]} did not finish within {timeout} s and was killed.", code="EXTERNAL_PROGRAM_TIMEOUT") from exc
    except OSError as exc:
        raise ConfigurationError(f"{argv[0]} could not be started: {exc}", code="EXTERNAL_PROGRAM_FAILED_TO_START") from exc
    duration = round(time.monotonic() - started, 6)
    truncated = len(completed.stdout) > cap or len(completed.stderr) > cap
    result = {
        "argv": [program, *argv[1:]], "returncode": completed.returncode, "success": completed.returncode == 0,
        "stdout": completed.stdout[:cap].decode("utf-8", errors="replace"),
        "stderr": completed.stderr[:cap].decode("utf-8", errors="replace"),
        "duration_s": duration, "truncated": truncated,
    }
    warnings = [f"output truncated to {cap} bytes per stream"] if truncated else []
    return candidate_ok(_MODULE, "run", input, result, options=options, warnings=warnings, implementation="subprocess")


def validate_output(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Check that an external program produced usable output.

    ``input``: output text, a file path (``pathlib.Path``) or a run result mapping.  Checks (each reported in
    ``checks`` as ``{name, passed, detail}``): ``not_empty`` (always), ``returncode`` equals ``expected_returncode``
    (default 0, only for run results), ``parses`` (the text parses with ``options.format``, when given), ``min_records``,
    ``required_columns`` (for tabular/JSON dictionary records) and ``must_contain`` (list of substrings).  Options:
    ``format``, ``min_records``, ``required_columns``, ``must_contain``, ``expected_returncode`` (and the parse options of
    ``parse_output``).  Result: ``valid`` (all checks passed) and ``checks``.
    """

    settings = _opts(options)
    text = _text_of(input)
    checks: list[dict[str, object]] = [{"name": "not_empty", "passed": bool(text.strip()), "detail": f"{len(text)} characters"}]
    if isinstance(input, Mapping) and "returncode" in input:
        expected = int(settings.get("expected_returncode", 0))
        checks.append({"name": "returncode", "passed": input["returncode"] == expected, "detail": f"got {input['returncode']}, expected {expected}"})
    records: list[object] | None = None
    if "format" in settings:
        try:
            _, records = _parse_records(text, settings)
            checks.append({"name": "parses", "passed": True, "detail": f"{len(records)} records"})
        except InputFormatError as exc:
            checks.append({"name": "parses", "passed": False, "detail": exc.message})
    if "min_records" in settings:
        count = None if records is None else len(records)
        checks.append({"name": "min_records", "passed": count is not None and count >= int(settings["min_records"]), "detail": f"{count} records, need {settings['min_records']}"})  # type: ignore[call-overload]
    if "required_columns" in settings:
        keys = set().union(*(r.keys() for r in records if isinstance(r, Mapping))) if records else set()
        missing = [str(c) for c in settings["required_columns"] if c not in keys]  # type: ignore[union-attr]
        checks.append({"name": "required_columns", "passed": records is not None and not missing, "detail": f"missing: {missing}"})
    for needle in settings.get("must_contain", []):  # type: ignore[union-attr]
        checks.append({"name": "must_contain", "passed": str(needle) in text, "detail": f"{needle!r}"})
    return candidate_ok(_MODULE, "validate_output", input, {"valid": all(c["passed"] for c in checks), "checks": checks}, options=options)


__all__ = ['build_command', 'check', 'convert_input', 'convert_output', 'make_report', 'parse_log', 'parse_output', 'run', 'validate_output']

__candidate_adapter_operations__ = set(['build_command', 'check', 'convert_input', 'convert_output', 'make_report', 'parse_log', 'parse_output', 'run', 'validate_output'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
