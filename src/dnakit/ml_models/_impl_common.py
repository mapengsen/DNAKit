"""Shared input parsing and output writing for the ``ml_models`` candidate operations."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_ROWS = 1_000_000
MAX_FILE_BYTES = 100_000_000
MAX_INLINE_CHARS = 2_000_000
_TABLE_SUFFIXES = {".csv": ",", ".tsv": "\t"}


class Unavailable(Exception):
    """Raised by an operation that cannot run here; the public wrapper reports ``not_implemented``."""

    def __init__(self, requires: str, message: str) -> None:
        super().__init__(message)
        self.requires = requires
        self.message = message


def bad_input(
    message: str, *, code: str = "INVALID_ML_INPUT", hint: str | None = None
) -> InputFormatError:
    """Build an ``InputFormatError`` with a stable code."""

    return InputFormatError(message, code=code, hint=hint)


def bad_option(
    message: str, *, code: str = "INVALID_ML_OPTION", hint: str | None = None
) -> ConfigurationError:
    """Build a ``ConfigurationError`` with a stable code."""

    return ConfigurationError(message, code=code, hint=hint)


def get_float(
    options: Mapping[str, object],
    key: str,
    default: float,
    *,
    low: float | None = None,
    high: float | None = None,
    inclusive_low: bool = True,
    inclusive_high: bool = True,
) -> float:
    """Read a finite float option within optional bounds."""

    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise bad_option(f"option '{key}' must be a finite number.")
    value = float(value)
    if low is not None and (value < low or (value == low and not inclusive_low)):
        raise bad_option(f"option '{key}' must be {'>=' if inclusive_low else '>'} {low}.")
    if high is not None and (value > high or (value == high and not inclusive_high)):
        raise bad_option(f"option '{key}' must be {'<=' if inclusive_high else '<'} {high}.")
    return value


def get_int(
    options: Mapping[str, object],
    key: str,
    default: int,
    *,
    low: int | None = None,
    high: int | None = None,
) -> int:
    """Read an integer option within optional bounds."""

    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int):
        raise bad_option(f"option '{key}' must be an integer.")
    if low is not None and value < low:
        raise bad_option(f"option '{key}' must be >= {low}.")
    if high is not None and value > high:
        raise bad_option(f"option '{key}' must be <= {high}.")
    return value


def get_str(options: Mapping[str, object], key: str, default: str | None) -> str | None:
    """Read a non-empty string option."""

    value = options.get(key, default)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise bad_option(f"option '{key}' must be a non-empty string.")
    return value


def get_choice(
    options: Mapping[str, object], key: str, default: str, choices: Sequence[str]
) -> str:
    """Read a string option restricted to ``choices`` (case-insensitive)."""

    value = options.get(key, default)
    if not isinstance(value, str) or value.strip().lower() not in choices:
        raise bad_option(f"option '{key}' must be one of {list(choices)}.")
    return value.strip().lower()


def to_float(value: object, *, what: str) -> float:
    """Convert a number or numeric string to a finite float."""

    if isinstance(value, bool):
        return float(value)
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise bad_input(f"{what} contains a non-numeric value: {value!r}.") from None
    if not math.isfinite(number):
        raise bad_input(f"{what} contains a non-finite value: {value!r}.")
    return number


def float_list(values: object, *, what: str) -> list[float]:
    """Convert a flat sequence into a list of finite floats."""

    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise bad_input(f"{what} must be a list of numbers.")
    if len(values) > MAX_ROWS:
        raise bad_input(f"{what} exceeds the limit of {MAX_ROWS} values.")
    return [to_float(item, what=what) for item in values]


def _read_text(path: Path) -> str:
    try:
        if path.stat().st_size > MAX_FILE_BYTES:
            raise bad_input(
                f"file is larger than {MAX_FILE_BYTES} bytes: {path}", code="ML_INPUT_TOO_LARGE"
            )
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise bad_input(
            f"could not read UTF-8 text file {path}: {exc}", code="ML_INPUT_UNREADABLE"
        ) from exc


def _is_path_like(value: object) -> bool:
    if isinstance(value, Path):
        return True
    return (
        isinstance(value, str)
        and "\n" not in value
        and 0 < len(value) < 4096
        and os.path.isfile(value)
    )


def read_fasta(text: str) -> list[dict[str, object]]:
    """Parse FASTA text into ``{"id", "sequence"}`` rows."""

    rows: list[dict[str, object]] = []
    name: str | None = None
    parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                rows.append({"id": name, "sequence": "".join(parts)})
            name, parts = line[1:].split()[0] if line[1:].split() else "", []
        elif name is None:
            raise bad_input("FASTA data must start with a '>' header line.")
        else:
            parts.append(line)
    if name is not None:
        rows.append({"id": name, "sequence": "".join(parts)})
    return rows


def load_rows(value: object, *, what: str = "input") -> list[dict[str, object]]:
    """Return table rows from records, a column mapping, or a CSV/TSV/JSON/JSONL/FASTA file path.

    Accepted forms: a list of mappings; a mapping with a ``records``/``rows``/``data`` list;
    a mapping of equally long column lists; or the path of a ``.csv``, ``.tsv``, ``.json``,
    ``.jsonl``, ``.fa``/``.fasta`` file (UTF-8, at most 100 MB, at most 1,000,000 rows).
    """

    if _is_path_like(value):
        path = Path(str(value)).expanduser()
        suffix = path.suffix.lower()
        text = _read_text(path)
        if suffix in _TABLE_SUFFIXES:
            reader = csv.DictReader(io.StringIO(text), delimiter=_TABLE_SUFFIXES[suffix])
            rows = []
            for row in reader:
                rows.append({str(k): v for k, v in row.items() if k is not None})
                if len(rows) > MAX_ROWS:
                    raise bad_input(f"{what} exceeds {MAX_ROWS} rows.", code="ML_INPUT_TOO_LARGE")
            return rows
        if suffix == ".jsonl":
            rows = []
            for number, line in enumerate(text.splitlines(), 1):
                if line.strip():
                    try:
                        item = json.loads(line)
                    except json.JSONDecodeError as exc:
                        raise bad_input(f"invalid JSON on line {number} of {path.name}.") from exc
                    if not isinstance(item, Mapping):
                        raise bad_input(f"line {number} of {path.name} is not a JSON object.")
                    rows.append(dict(item))
            return rows
        if suffix == ".json":
            try:
                return load_rows(json.loads(text), what=what)
            except json.JSONDecodeError as exc:
                raise bad_input(f"invalid JSON file {path.name}: {exc}") from exc
        if suffix in {".fa", ".fasta", ".fna"}:
            return read_fasta(text)
        raise bad_input(
            f"unsupported file type '{suffix}' for {what}.",
            code="UNSUPPORTED_ML_INPUT_FORMAT",
            hint="Use .csv, .tsv, .json, .jsonl or .fasta, or pass records directly.",
        )
    if isinstance(value, Mapping):
        for key in ("records", "rows", "data"):
            inner = value.get(key)
            if isinstance(inner, Sequence) and not isinstance(inner, (str, bytes)):
                return load_rows(inner, what=what)
        columns = [
            v for v in value.values() if isinstance(v, Sequence) and not isinstance(v, (str, bytes))
        ]
        if value and len(columns) == len(value) and len({len(c) for c in columns}) == 1:
            size = len(columns[0])
            if size > MAX_ROWS:
                raise bad_input(f"{what} exceeds {MAX_ROWS} rows.", code="ML_INPUT_TOO_LARGE")
            return [{str(k): value[k][i] for k in value} for i in range(size)]
        raise bad_input(f"{what} mapping must contain a records list or equally long column lists.")
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) > MAX_ROWS:
            raise bad_input(f"{what} exceeds {MAX_ROWS} rows.", code="ML_INPUT_TOO_LARGE")
        rows = []
        for item in value:
            if not isinstance(item, Mapping):
                raise bad_input(f"{what} must be a list of mappings (records).")
            rows.append({str(k): v for k, v in item.items()})
        return rows
    raise bad_input(
        f"{what} must be records, a column mapping, or a file path.",
        hint="Pass a list of dicts, a dict of equally long lists, or a CSV/TSV/JSON/JSONL/FASTA path.",
    )


def render_rows(
    rows: Sequence[Mapping[str, object]], fmt: str, columns: Sequence[str] | None = None
) -> str:
    """Serialise rows as csv, tsv, json, jsonl or fasta text."""

    names = list(columns) if columns else list(dict.fromkeys(k for row in rows for k in row))
    if fmt in {"csv", "tsv"}:
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="," if fmt == "csv" else "\t", lineterminator="\n")
        writer.writerow(names)
        for row in rows:
            writer.writerow(["" if row.get(name) is None else row.get(name) for name in names])
        return buffer.getvalue()
    selected = [{name: row.get(name) for name in names} for row in rows]
    if fmt == "json":
        return json.dumps(selected, indent=2, ensure_ascii=False, default=str) + "\n"
    if fmt == "jsonl":
        return "".join(json.dumps(r, ensure_ascii=False, default=str) + "\n" for r in selected)
    if fmt == "fasta":
        lines = []
        for index, row in enumerate(rows):
            if "sequence" not in row:
                raise bad_input("fasta output needs a 'sequence' column in every record.")
            lines.append(f">{row.get('id', f'seq{index + 1}')}\n{row['sequence']}\n")
        return "".join(lines)
    raise bad_option(
        f"unsupported output format '{fmt}'.", hint="Use csv, tsv, json, jsonl or fasta."
    )


def output_format(options: Mapping[str, object], default: str, allowed: Sequence[str]) -> str:
    """Resolve the output format from ``options['format']`` or the ``output`` suffix."""

    fmt = options.get("format")
    if fmt is None:
        output = options.get("output")
        suffix = Path(str(output)).suffix.lower().lstrip(".") if output else ""
        fmt = suffix if suffix in allowed else default
    if not isinstance(fmt, str) or fmt.lower() not in allowed:
        raise bad_option(f"option 'format' must be one of {list(allowed)}.")
    return fmt.lower()


def write_or_return(text: str, options: Mapping[str, object]) -> dict[str, object]:
    """Write ``text`` to ``options['output']`` (atomic, no silent overwrite) or return it inline.

    Option ``overwrite`` (default False) allows replacing an existing file. Without an
    ``output`` path nothing is written and the text is returned under ``content``.
    """

    output = options.get("output")
    data = text.encode("utf-8")
    info: dict[str, object] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    if output is None:
        if len(text) > MAX_INLINE_CHARS:
            raise bad_option(
                "content is too large to return inline; set options['output'].",
                code="ML_OUTPUT_TOO_LARGE",
            )
        info["content"] = text
        info["output_path"] = None
        return info
    if not isinstance(output, (str, os.PathLike)) or not str(output).strip() or "\0" in str(output):
        raise bad_option("option 'output' must be a file path.")
    path = Path(output).expanduser()
    if path.exists() and not options.get("overwrite", False):
        raise bad_option(
            f"output file already exists: {path}",
            code="ML_OUTPUT_EXISTS",
            hint="Set options['overwrite']=True or choose another path.",
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(dir=str(path.parent), prefix=".dnakit-", suffix=".tmp")
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
        os.replace(temp_name, path)
    except BaseException:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
        raise
    info["output_path"] = str(path)
    return info


def sequences_from(
    value: object, options: Mapping[str, object], *, limit: int = 10_000
) -> list[tuple[str, str]]:
    """Return ``(id, sequence)`` pairs from a sequence string, list of strings, records or FASTA/table path.

    Option ``sequence_column`` (default ``sequence``) and ``id_column`` (default ``id``) pick
    the columns of record input. A string that is not a file path is one sequence (or FASTA
    text when it starts with ``>``).
    """

    seq_col = str(options.get("sequence_column", "sequence"))
    id_col = str(options.get("id_column", "id"))
    pairs: list[tuple[str, str]] = []
    if isinstance(value, str) and not _is_path_like(value):
        text = value.strip()
        if text.startswith(">"):
            pairs = [(str(r["id"]), str(r["sequence"])) for r in read_fasta(text)]
        elif text:
            pairs = [("seq1", "".join(text.split()))]
    elif (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and value
        and all(isinstance(item, str) for item in value)
    ):
        pairs = [(f"seq{i + 1}", "".join(str(item).split())) for i, item in enumerate(value)]
    else:
        for index, row in enumerate(load_rows(value, what="sequence input")):
            if seq_col not in row or not isinstance(row[seq_col], str):
                raise bad_input(f"record {index} has no string column '{seq_col}'.")
            pairs.append(
                (str(row.get(id_col, f"seq{index + 1}")), "".join(str(row[seq_col]).split()))
            )
    if not pairs:
        raise bad_input("no sequences were provided.")
    if len(pairs) > limit:
        raise bad_input(
            f"at most {limit} sequences are accepted per call.", code="ML_INPUT_TOO_LARGE"
        )
    for name, sequence in pairs:
        if not sequence:
            raise bad_input(f"sequence '{name}' is empty.")
    return pairs
