# External adapters

This category documents **9** functions of the `external_adapters` module, and all 9 compute a real result. `run` starts an installed program from an argument list (never through a shell, with a timeout); it returns `status: "not_implemented"` only when the requested program is not installed. `check` looks programs up on `PATH`; the other seven (`build_command`, `convert_input`, `convert_output`, `make_report`, `parse_log`, `parse_output`, `validate_output`) are pure text transformations. None of them wraps a scientific tool by itself. Example outputs are real; `input_summary` is abbreviated to `{...}`.

## 1) build_command

- **Purpose:** Build an argv list (never a shell string) for an external program from a program, flags and arguments. Nothing is executed.

- **Method:** The program comes first, then `options.subcommand` (when given), then the flags, then the arguments. Flags are emitted in mapping order: `true` gives the bare flag, `false`/`null` omits it, a list repeats the flag once per item, anything else gives `flag value`; flag names must start with `-` and no token may contain a NUL byte or be a boolean. `options.flags` and `options.args` are appended after those of the input. Because the flags precede the arguments, `{"program": "samtools", "args": ["view", "x.bam"], "flags": {"-b": true}}` would give `samtools -b view x.bam`, which is wrong for a program that expects a sub-command first. Give the sub-command as `options.subcommand` instead: it is emitted right after the program and before the flags (`samtools view -b x.bam`). For the argv-list and command-line-string forms the subcommand is inserted after the first token. Without `subcommand` the order is unchanged.

- **API:** `dnakit.external_adapters.build_command(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a program name, an argv list, a command-line string (split with `shlex`, no shell involved), a path-like, or a mapping `{"program": str, "args": [...], "flags": {...}}`. `options`: `subcommand` (a string such as `"view"`, or a list such as `["index", "build"]`; each item a non-empty string, otherwise `ConfigurationError` `INVALID_OPTION`), `flags` and `args` (same forms as above, appended). Invalid commands raise `InputFormatError` (`INVALID_COMMAND`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `argv`, `program` and `command_line` (`shlex.join`, for display only).

- **Example code:**

```python
import json
from dnakit.external_adapters import build_command

result = build_command(
    {"program": "fastqc", "args": ["reads.fq"], "flags": {"--quiet": True, "--threads": 4, "--nogroup": None, "--adapters": ["a.txt", "b.txt"]}}
)
print(json.dumps(result, indent=2))

# a program with a sub-command: options.subcommand goes right after the program, before the flags
view = build_command(
    {"program": "samtools", "args": ["x.bam"], "flags": {"-b": True}}, options={"subcommand": "view"}
)
print(json.dumps(view["result"], indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "build_command",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "argv": [
      "fastqc",
      "--quiet",
      "--threads",
      "4",
      "--adapters",
      "a.txt",
      "--adapters",
      "b.txt",
      "reads.fq"
    ],
    "program": "fastqc",
    "command_line": "fastqc --quiet --threads 4 --adapters a.txt --adapters b.txt reads.fq"
  }
}
{
  "argv": [
    "samtools",
    "view",
    "-b",
    "x.bam"
  ],
  "program": "samtools",
  "command_line": "samtools view -b x.bam"
}
```

## 2) check

- **Purpose:** Check whether external programs are installed, optionally reading their version.

- **Method:** Each program is resolved with `shutil.which`. Only when `version_args` is given is `[program, *version_args]` run (no shell, no stdin, with a timeout) and the first line of its combined stdout and stderr returned as `version`. A missing program is reported with `available: false`; the call itself succeeds.

- **API:** `dnakit.external_adapters.check(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a program name/path, or a list of them. `options`: `version_args` (list such as `["--version"]`; default none, so nothing is executed), `timeout` (seconds, default 10; on expiry `version` is `null` and `version_error` is set).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `all_available` and `programs`, one row per program with `program`, `available`, `path` and `version`. The `path` and `version` shown depend on the machine.

- **Example code:**

```python
import json
from dnakit.external_adapters import check

result = check(["python3", "no_such_tool_xyz"], options={"version_args": ["--version"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "check",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "all_available": false,
    "programs": [
      {
        "program": "python3",
        "available": true,
        "path": "/opt/conda/bin/python3",
        "version": "Python 3.8.12"
      },
      {
        "program": "no_such_tool_xyz",
        "available": false,
        "path": null,
        "version": null
      }
    ]
  }
}
```

## 3) convert_input

- **Purpose:** Convert DNA sequence records into the text format an external tool expects (FASTA, FASTQ or TSV).

- **Method:** Sequences are normalized with `dnakit.sequence.normalize_dna_sequence` (uppercase, `U` to `T`, IUPAC letters only; anything else raises `InputFormatError` `INVALID_SEQUENCE`). Records without an id are named `seq1`, `seq2`, ... FASTQ uses one quality character for every base (no real qualities exist), so it is a format conversion, not quality data. FASTA headers keep only the id (descriptions are dropped).

- **API:** `dnakit.external_adapters.convert_input(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text (starting with `>`), a path-like to a FASTA file, a single sequence string, or a list of sequences, `(id, sequence)` pairs or mappings with `id`/`sequence`. `options`: `target_format` (`"fasta"` default, `"fastq"`, `"tsv"`), `line_width` (FASTA wrapping, default 80, 0 = no wrapping), `quality` (one character used for every base in FASTQ, default `"I"`).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format`, `text` and `record_count`.

- **Example code:**

```python
import json
from dnakit.external_adapters import convert_input

result = convert_input(">a desc\nacgu\n>b\nGGAA\n", options={"target_format": "fastq"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "convert_input",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fastq",
    "text": "@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n",
    "record_count": 2
  }
}
```

## 4) convert_output

- **Purpose:** Convert external-tool output (TSV, CSV, JSON, JSON Lines, FASTA, plain lines) into JSON, JSON Lines, TSV or CSV text.

- **Method:** The source is parsed exactly as in `parse_output` and then written in `target_format`. Tabular targets need dictionary records; their columns are the union of the keys in first-seen order and `null` becomes an empty cell. JSON is written with 2-space indent.

- **API:** `dnakit.external_adapters.convert_output(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a text string, a `pathlib.Path` of an existing file, or a run-result mapping holding `stdout` (or `text`/`output`/`path`). `options`: `target_format` (`"json"` default, `"jsonl"`, `"tsv"`, `"csv"`) and the parse options of `parse_output`: `format` (`"auto"` default: FASTA if the text starts with `>`, JSON if it starts with `{` or `[`, TSV if it contains tabs, else lines; or `tsv`, `csv`, `json`, `jsonl`, `fasta`, `lines`), `delimiter`, `header` (default true: the first row names the columns), `columns` (explicit column names, used instead of a header), `comment_prefix` (default `"#"`; such lines are skipped in TSV/CSV), `convert_types` (default true: integer and float cells become numbers).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format` (the target), `text` and `record_count`.

- **Example code:**

```python
import json
from dnakit.external_adapters import convert_output

result = convert_output("#comment\nid\tcount\tname\nr1\t3\talpha\nr2\t4.5\tbeta\n", options={"target_format": "csv"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "convert_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "csv",
    "text": "id,count,name\nr1,3,alpha\nr2,4.5,beta\n",
    "record_count": 2
  }
}
```

## 5) make_report

- **Purpose:** Summarize one external-program run (a result of `run`) as a Markdown or JSON report.

- **Method:** The report lists the command (`shlex.join` of `argv`), return code, success (return code 0), duration, the severity counts from the `parse_log` keyword rule applied to stderr followed by stdout, and the first `max_lines` lines of stdout and stderr. Nothing is executed.

- **API:** `dnakit.external_adapters.make_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the `run` result mapping (needs `argv` and `returncode`; optional `stdout`, `stderr`, `duration_s`) or the whole `run` envelope. `options`: `format` (`"markdown"` default or `"json"`), `title` (default `"External program run"`), `max_lines` (default 20).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format` and `report` (the Markdown text, or with `format: "json"` a mapping with `command`, `returncode`, `success`, `duration_s`, `log`, `stdout_head`, `stderr_head`).

- **Example code:**

```python
import json
from dnakit.external_adapters import make_report

run_result = {
    "argv": ["/usr/bin/samtools", "flagstat", "x.bam"],
    "returncode": 1,
    "stdout": "10 + 0 in total\n",
    "stderr": "[W::hts_open] warning: low quality\n[E::main] error: truncated file\n",
    "duration_s": 0.5,
}
result = make_report(run_result, options={"title": "samtools flagstat"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "make_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "markdown",
    "report": "# samtools flagstat\n\n- command: `/usr/bin/samtools flagstat x.bam`\n- return code: 1\n- success: no\n- duration: 0.5 s\n- log: 1 error, 1 warning line(s)\n\n## stdout (head)\n\n    10 + 0 in total\n\n## stderr (head)\n\n    [W::hts_open] warning: low quality\n    [E::main] error: truncated file\n"
  }
}
```

## 6) parse_log

- **Purpose:** Classify the lines of an external-program log by severity with a heuristic keyword rule.

- **Method:** A line is classified by the first matching whole word, tested in this order and case-insensitively: error (`fatal`, `critical`, `error`, `exception`, `traceback`, `segmentation fault`), warning (`warn`, `warning`), `info`, `debug`; lines without keyword are unclassified (counted only in `line_count`). The rule is a heuristic, not a parser of any specific tool's log format. Logs above 200,000 lines raise `InputFormatError`.

- **API:** `dnakit.external_adapters.parse_log(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: log text, a `pathlib.Path` of a file, or a run-result mapping (stderr then stdout are scanned if `stderr` is present, else `stdout`). `options` is unused.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `line_count`, `levels` (counts for error, warning, info, debug), `errors` and `warnings` (up to 100 samples `{line, text}` each, text cut at 500 characters) and `has_errors`.

- **Example code:**

```python
import json
from dnakit.external_adapters import parse_log

result = parse_log("INFO start\nWARNING low coverage\nERROR failed to open file\nplain line\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "parse_log",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "line_count": 4,
    "levels": {
      "error": 1,
      "warning": 1,
      "info": 1,
      "debug": 0
    },
    "errors": [
      {
        "line": 3,
        "text": "ERROR failed to open file"
      }
    ],
    "warnings": [
      {
        "line": 2,
        "text": "WARNING low coverage"
      }
    ],
    "has_errors": true
  }
}
```

## 7) parse_output

- **Purpose:** Parse external-tool output into structured records.

- **Method:** TSV/CSV rows become dictionaries keyed by the header (or `columns`, or `column_1`, `column_2`, ... without a header); a row with a different number of fields raises `InputFormatError` (`INVALID_OUTPUT`). With `convert_types` a cell is converted to an integer, then to a float, otherwise kept as text. FASTA records are `{id, description, sequence}` (sequence lines joined). JSON gives the list (a single value is wrapped in a list), JSON Lines one record per line, `lines` the non-empty lines. Input larger than 50,000,000 bytes is refused.

- **API:** `dnakit.external_adapters.parse_output(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a text string, a `pathlib.Path` of an existing file, or a run-result mapping holding `stdout` (or `text`/`output`/`path`). `options`: `format` (`"auto"` default: FASTA if the text starts with `>`, JSON if it starts with `{` or `[`, TSV if it contains tabs, else lines; or `tsv`, `csv`, `json`, `jsonl`, `fasta`, `lines`), `delimiter`, `header` (default true: the first row names the columns), `columns` (explicit column names, used instead of a header), `comment_prefix` (default `"#"`; such lines are skipped in TSV/CSV), `convert_types` (default true: integer and float cells become numbers).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `format` (the format actually used), `records` and `record_count`.

- **Example code:**

```python
import json
from dnakit.external_adapters import parse_output

result = parse_output("id\tcount\tname\nr1\t3\talpha\nr2\t4.5\tbeta\n")
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "parse_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "records": [
      {
        "id": "r1",
        "count": 3,
        "name": "alpha"
      },
      {
        "id": "r2",
        "count": 4.5,
        "name": "beta"
      }
    ],
    "record_count": 2
  }
}
```

## 8) run

- **Purpose:** Run an external program (argv list, no shell, with a timeout) and capture its output.

- **Method:** The command is turned into an argv list as in `build_command` and the program is resolved with `shutil.which`; it is started with `subprocess.run` without a shell. If the program is not installed the call returns `status: "not_implemented"` with `required_dependency: "the external program '<name>' on PATH"` and a `message` saying it was not found on PATH. A non-zero exit code is reported in the result, not raised. On timeout the process is killed and `BackendTimeoutError` (`EXTERNAL_PROGRAM_TIMEOUT`) is raised; a program that cannot be started raises `ConfigurationError`. stdout and stderr are each cut at `max_output_bytes` and flagged with `truncated` plus a warning.

- **API:** `dnakit.external_adapters.run(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a program name, argv list, command-line string, path-like or mapping (see `build_command`). `options`: `timeout` (seconds, default 60, must be in (0, 3600]), `cwd` (existing directory), `env` (mapping of extra environment variables), `stdin` (text sent to the program), `max_output_bytes` (default 1,000,000), `flags` and `args` (appended) and `subcommand` (inserted right after the program, exactly as in `build_command`, so `run` executes the argv that `build_command` shows).

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result` (`provenance.implementation` is `subprocess`); `result` has `argv` (with the resolved program path), `returncode`, `success`, `stdout`, `stderr`, `duration_s` and `truncated`. If the program is missing, `status` is `not_implemented` instead and `result` is absent. In the example the resolved path of `echo` depends on the machine and the measured `duration_s` is replaced by the placeholder `"<seconds>"`.

- **Example code:**

```python
import json
from dnakit.external_adapters import run

result = run(["echo", "hello"], options={"timeout": 10})
result["result"]["duration_s"] = "<seconds>"
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "run",
  "provenance": {
    "implementation": "subprocess",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "argv": [
      "/usr/bin/echo",
      "hello"
    ],
    "returncode": 0,
    "success": true,
    "stdout": "hello\n",
    "stderr": "",
    "duration_s": "<seconds>",
    "truncated": false
  }
}
```

## 9) validate_output

- **Purpose:** Check that an external program produced usable output.

- **Method:** Each check is reported in `checks` as `{name, passed, detail}`: `not_empty` (always), `returncode` equals `expected_returncode` (only for run results), `parses` (the text parses with the given `format`), `min_records`, `required_columns` (for tabular or JSON dictionary records) and one `must_contain` check per substring. `min_records` and `required_columns` fail when the output was not parsed (no `format`). Nothing is executed.

- **API:** `dnakit.external_adapters.validate_output(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a text string, a `pathlib.Path` of an existing file, or a run-result mapping holding `stdout` (or `text`/`output`/`path`). `options`: `format` (enables the `parses` check), `min_records`, `required_columns`, `must_contain` (list of substrings), `expected_returncode` (default 0) and the parse options of `parse_output`.

- **Output:** Envelope with `status: "ok"`, `module`, `operation`, `provenance`, `input_summary`, `warnings` and `result`; `result` has `valid` (all checks passed) and `checks`.

- **Example code:**

```python
import json
from dnakit.external_adapters import validate_output

result = validate_output(
    "id\tcount\nr1\t3\nr2\t4\n",
    options={"format": "tsv", "min_records": 3, "required_columns": ["id", "count", "name"], "must_contain": ["r1"]},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "external_adapters",
  "operation": "validate_output",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "checks": [
      {
        "name": "not_empty",
        "passed": true,
        "detail": "19 characters"
      },
      {
        "name": "parses",
        "passed": true,
        "detail": "2 records"
      },
      {
        "name": "min_records",
        "passed": false,
        "detail": "2 records, need 3"
      },
      {
        "name": "required_columns",
        "passed": false,
        "detail": "missing: ['name']"
      },
      {
        "name": "must_contain",
        "passed": true,
        "detail": "'r1'"
      }
    ]
  }
}
```
