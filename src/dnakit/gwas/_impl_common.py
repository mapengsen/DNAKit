"""Shared input parsing, option validation and output helpers for the statistical-genetics
candidates.

Used by ``dnakit.gwas``, ``dnakit.population``, ``dnakit.breeding``,
``dnakit.forensic_reproductive`` and ``dnakit.experimental_design``.  Pure standard library.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
from collections.abc import Mapping, Sequence
from numbers import Real
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, DNAKitError, InputFormatError

MAX_TABLE_ROWS = 500_000
MAX_MATRIX_CELLS = 20_000_000
_MISSING_TOKENS = {"", "na", "nan", "n/a", ".", "null", "none"}


# --------------------------------------------------------------------------- options
def get_options(options: Mapping[str, object] | None) -> dict[str, Any]:
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ConfigurationError(
            "options must be a mapping of option names to values.",
            code="INVALID_OPTIONS",
            hint="Pass options={'key': value}.",
        )
    return dict(options)


def opt_number(
    opts: Mapping[str, Any],
    key: str,
    default: float | None,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    integer: bool = False,
    exclusive: bool = False,
) -> Any:
    value = opts.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ConfigurationError(
            f"option '{key}' must be a finite number, got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options['{key}'] to a number.",
        )
    number = float(value)
    if integer:
        if number != int(number):
            raise ConfigurationError(
                f"option '{key}' must be an integer, got {value!r}.",
                code="INVALID_OPTION",
                hint=f"Set options['{key}'] to a whole number.",
            )
        number = int(number)
    low_bad = minimum is not None and (number <= minimum if exclusive else number < minimum)
    high_bad = maximum is not None and (number >= maximum if exclusive else number > maximum)
    if low_bad or high_bad:
        bounds = (
            f"{'>' if exclusive else '>='} {minimum}"
            if low_bad
            else f"{'<' if exclusive else '<='} {maximum}"
        )
        raise ConfigurationError(
            f"option '{key}' must be {bounds}, got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Choose a valid value for options['{key}'].",
        )
    return number


def opt_choice(opts: Mapping[str, Any], key: str, default: str, choices: Sequence[str]) -> str:
    value = opts.get(key, default)
    if not isinstance(value, str) or value.strip().lower() not in choices:
        raise ConfigurationError(
            f"option '{key}' must be one of {list(choices)}, got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Use one of: {', '.join(choices)}.",
        )
    return value.strip().lower()


def opt_bool(opts: Mapping[str, Any], key: str, default: bool) -> bool:
    value = opts.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"option '{key}' must be true or false, got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options['{key}'] to True or False.",
        )
    return value


def require_option(opts: Mapping[str, Any], key: str, why: str) -> Any:
    if opts.get(key) is None:
        raise ConfigurationError(
            f"option '{key}' is required: {why}",
            code="MISSING_OPTION",
            hint=f"Pass options={{'{key}': ...}}.",
        )
    return opts[key]


# --------------------------------------------------------------------------- scalar parsing
def to_number(value: object) -> float | None:
    """Numeric value of a cell, or None when it is missing; raises ValueError for non-numeric
    text.
    """

    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError("boolean is not a number")
    if isinstance(value, Real):
        number = float(value)
        return None if math.isnan(number) else number
    text = str(value).strip()
    if text.lower() in _MISSING_TOKENS:
        return None
    return float(text)


def _convert_cell(text: str) -> object:
    stripped = text.strip()
    if stripped.lower() in _MISSING_TOKENS:
        return None
    try:
        number = float(stripped)
    except ValueError:
        return stripped
    if re.fullmatch(r"[+-]?\d+", stripped):
        return int(stripped)
    return number


# --------------------------------------------------------------------------- tables
def read_text_source(source: object, *, what: str = "input") -> str:
    """Text of a path / path string (existing file) / literal text; supports ``.gz`` files."""

    if isinstance(source, bytes):
        return source.decode("utf-8", errors="replace")
    path: Path | None = None
    if isinstance(source, Path):
        path = source
    elif (
        isinstance(source, str) and "\n" not in source and "\t" not in source and len(source) < 4096
    ):
        candidate = Path(source)
        try:
            if candidate.is_file():
                path = candidate
        except OSError:
            path = None
    if path is not None:
        if not path.is_file():
            raise InputFormatError(
                f"{what} file not found: {path}",
                code="FILE_NOT_FOUND",
                hint="Check the path or pass the table text/records directly.",
            )
        data = path.read_bytes()
        if path.suffix == ".gz" or data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        return data.decode("utf-8", errors="replace")
    if isinstance(source, str):
        return source
    raise InputFormatError(
        f"{what} must be text, a file path or a table of records, got {type(source).__name__}.",
        code="INVALID_INPUT_TYPE",
        hint="Pass records, a column dict, a path or delimited text.",
    )


def parse_delimited(
    text: str, *, delimiter: str | None = None, max_rows: int = MAX_TABLE_ROWS
) -> list[dict[str, object]]:
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("##")]
    if len(lines) < 1:
        raise InputFormatError(
            "the table text is empty.", code="EMPTY_TABLE", hint="Provide a header line and rows."
        )
    header_line = lines[0].lstrip("#")
    if delimiter is None:
        if "\t" in header_line:
            delimiter = "\t"
        elif "," in header_line:
            delimiter = ","
        elif ";" in header_line:
            delimiter = ";"
        else:
            delimiter = " "
    if delimiter == " ":
        rows = [ln.lstrip("#").split() if i == 0 else ln.split() for i, ln in enumerate(lines)]
    else:
        rows = list(
            csv.reader(io.StringIO("\n".join([header_line, *lines[1:]])), delimiter=delimiter)
        )
    header = [h.strip() for h in rows[0]]
    if len(set(header)) != len(header) or not all(header):
        raise InputFormatError(
            "the table header has empty or duplicated column names.",
            code="INVALID_HEADER",
            hint="Give every column a unique name.",
        )
    body = rows[1:]
    if len(body) > max_rows:
        raise InputFormatError(
            f"the table has {len(body)} rows, more than the limit of {max_rows}.",
            code="TABLE_TOO_LARGE",
            hint="Raise options['max_rows'] or pre-filter the table.",
        )
    records = []
    for number, row in enumerate(body, start=2):
        if len(row) != len(header):
            raise InputFormatError(
                f"line {number} has {len(row)} fields but the header has {len(header)}.",
                code="RAGGED_TABLE",
                hint="Every row needs the same number of fields as the header.",
            )
        records.append({name: _convert_cell(cell) for name, cell in zip(header, row, strict=False)})
    return records


def load_table(
    source: object,
    *,
    what: str = "table",
    max_rows: int = MAX_TABLE_ROWS,
    delimiter: str | None = None,
) -> list[dict[str, object]]:
    """Records from a list of mappings, a column mapping, a ``{'records': [...]}`` wrapper, a
    path or text.
    """

    if isinstance(source, Mapping):
        for key in ("records", "rows", "table", "data"):
            if (
                key in source
                and isinstance(source[key], (Sequence, Mapping))
                and not isinstance(source[key], str)
            ):
                return load_table(source[key], what=what, max_rows=max_rows, delimiter=delimiter)
        columns = {str(k): v for k, v in source.items()}
        if columns and all(
            isinstance(v, Sequence) and not isinstance(v, (str, bytes)) for v in columns.values()
        ):
            lengths = {len(v) for v in columns.values()}
            if len(lengths) != 1:
                raise InputFormatError(
                    f"{what}: columns have different lengths {sorted(lengths)}.",
                    code="RAGGED_TABLE",
                    hint="Give every column the same number of values.",
                )
            length = lengths.pop()
            if length > max_rows:
                raise InputFormatError(
                    f"{what} has {length} rows, more than the limit of {max_rows}.",
                    code="TABLE_TOO_LARGE",
                    hint="Raise options['max_rows'] or pre-filter the table.",
                )
            return [{k: columns[k][i] for k in columns} for i in range(length)]
        raise InputFormatError(
            f"{what}: a mapping must hold equal-length columns or a 'records' list.",
            code="INVALID_TABLE",
            hint="Use a list of dict records or a dict of equal-length lists.",
        )
    if isinstance(source, Sequence) and not isinstance(source, (str, bytes)):
        if len(source) > max_rows:
            raise InputFormatError(
                f"{what} has {len(source)} rows, more than the limit of {max_rows}.",
                code="TABLE_TOO_LARGE",
                hint="Raise options['max_rows'] or pre-filter the table.",
            )
        if not all(isinstance(r, Mapping) for r in source):
            raise InputFormatError(
                f"{what}: every record must be a mapping of column names to values.",
                code="INVALID_TABLE",
                hint="Use a list of dict records.",
            )
        return [{str(k): v for k, v in r.items()} for r in source]
    return parse_delimited(
        read_text_source(source, what=what), delimiter=delimiter, max_rows=max_rows
    )


def serialize_records(
    records: Sequence[Mapping[str, object]],
    columns: Sequence[str],
    delimiter: str,
    missing: str = "NA",
) -> str:
    out = io.StringIO()
    writer = csv.writer(out, delimiter=delimiter, lineterminator="\n")
    writer.writerow(columns)
    for rec in records:
        row = []
        for col in columns:
            value = rec.get(col)
            if value is None or (isinstance(value, float) and math.isnan(value)):
                row.append(missing)
            elif isinstance(value, float):
                row.append(repr(value))
            else:
                row.append(value)
        writer.writerow(row)
    return out.getvalue()


def write_text_file(path: object, text: str, *, overwrite: bool = False) -> str:
    """Write ``text`` to an explicitly requested path (gzip when it ends in .gz); refuses to
    overwrite by default.
    """

    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise ConfigurationError(
            "the output path must be a non-empty string or Path.",
            code="INVALID_OUTPUT_PATH",
            hint="Pass options={'path': 'out.tsv'}.",
        )
    target = Path(path)
    if target.exists() and not overwrite:
        raise ConfigurationError(
            f"{target} already exists.",
            code="OUTPUT_EXISTS",
            hint="Pass options['overwrite']=True or choose another path.",
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    data = text.encode("utf-8")
    target.write_bytes(gzip.compress(data) if target.suffix == ".gz" else data)
    return str(target)


# --------------------------------------------------------------------------- genotype matrices
def parse_genotype_matrix(
    source: object, *, what: str = "genotypes", min_samples: int = 1, min_variants: int = 1
) -> dict[str, Any]:
    """Validate a samples x variants dosage matrix.

    Accepts a list of rows, or a mapping with the matrix under
    ``genotypes``/``matrix``/``dosages`` and optional ``sample_ids``, ``variant_ids``,
    ``chrom``, ``pos``, ``populations`` lists.  Entries are dosages of the alternate/effect
    allele in [0, 2] (fractional imputed dosages allowed) or ``None`` for missing.
    """

    extra: dict[str, Any] = {}
    matrix: object = source
    if isinstance(source, Mapping):
        for key in ("genotypes", "matrix", "dosages"):
            if key in source:
                matrix = source[key]
                break
        else:
            raise InputFormatError(
                f"{what}: the mapping needs a 'genotypes' (or 'matrix'/'dosages') key.",
                code="MISSING_GENOTYPES",
                hint="Use {'genotypes': [[0,1,2],...], 'sample_ids': [...]}.",
            )
        extra = {k: v for k, v in source.items() if k not in ("genotypes", "matrix", "dosages")}
    if not isinstance(matrix, Sequence) or isinstance(matrix, (str, bytes)) or not matrix:
        raise InputFormatError(
            f"{what} must be a non-empty list of rows (samples x variants).",
            code="INVALID_GENOTYPES",
            hint="Use a list of equal-length lists of dosages 0/1/2.",
        )
    rows: list[list[float | None]] = []
    width = None
    for r, row in enumerate(matrix):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            raise InputFormatError(
                f"{what}: row {r} is not a list of dosages.",
                code="INVALID_GENOTYPES",
                hint="Use a list of equal-length lists of dosages 0/1/2.",
            )
        if width is None:
            width = len(row)
        if len(row) != width:
            raise InputFormatError(
                f"{what}: row {r} has {len(row)} values, expected {width}.",
                code="RAGGED_GENOTYPES",
                hint="Every sample needs one dosage per variant.",
            )
        parsed: list[float | None] = []
        for c, cell in enumerate(row):
            try:
                value = to_number(cell)
            except ValueError as exc:
                raise InputFormatError(
                    f"{what}: entry ({r}, {c}) = {cell!r} is not numeric.",
                    code="INVALID_GENOTYPES",
                    hint="Dosages must be numbers in [0, 2] or null.",
                ) from exc
            if value is not None and not 0.0 <= value <= 2.0:
                raise InputFormatError(
                    f"{what}: dosage {value} at ({r}, {c}) is outside [0, 2].",
                    code="INVALID_GENOTYPES",
                    hint="Use alternate-allele dosages 0, 1, 2 (or imputed values in [0, 2]).",
                )
            parsed.append(value)
        rows.append(parsed)
    n, m = len(rows), int(width or 0)
    if n < min_samples or m < min_variants:
        raise InputFormatError(
            f"{what} needs at least {min_samples} samples and {min_variants} variants, got {n} x "
            f"{m}.",
            code="GENOTYPES_TOO_SMALL",
            hint="Provide more samples or variants.",
        )
    if n * m > MAX_MATRIX_CELLS:
        raise InputFormatError(
            f"{what} has {n * m} cells, more than the limit of {MAX_MATRIX_CELLS}.",
            code="GENOTYPES_TOO_LARGE",
            hint="Subset samples or variants first.",
        )

    def labels(key: str, count: int, prefix: str) -> list[str]:
        value = extra.get(key)
        if value is None:
            return [f"{prefix}{i + 1}" for i in range(count)]
        if (
            not isinstance(value, Sequence)
            or isinstance(value, (str, bytes))
            or len(value) != count
        ):
            raise InputFormatError(
                f"{what}: '{key}' must be a list with {count} entries.",
                code="INVALID_LABELS",
                hint=f"Give one entry per {'sample' if key.startswith('sample') else 'variant'}.",
            )
        text = [str(v) for v in value]
        if key in ("sample_ids", "variant_ids") and len(set(text)) != count:
            raise InputFormatError(
                f"{what}: '{key}' contains duplicates.",
                code="DUPLICATE_IDS",
                hint="Identifiers must be unique.",
            )
        return text

    result: dict[str, Any] = {
        "matrix": rows,
        "samples": labels("sample_ids", n, "sample"),
        "variants": labels("variant_ids", m, "variant"),
        "n_samples": n,
        "n_variants": m,
    }
    for key in ("chrom", "pos", "populations", "ref", "alt"):
        if key in extra:
            value = extra[key]
            expected = n if key == "populations" else m
            if (
                not isinstance(value, Sequence)
                or isinstance(value, (str, bytes))
                or len(value) != expected
            ):
                raise InputFormatError(
                    f"{what}: '{key}' must be a list with {expected} entries.",
                    code="INVALID_LABELS",
                    hint="Give one entry per "
                    + ("sample" if key == "populations" else "variant")
                    + ".",
                )
            result[key] = list(value)
    return result


def column_values(matrix: Sequence[Sequence[float | None]], j: int) -> list[float | None]:
    return [row[j] for row in matrix]


def dumps_json(value: object) -> str:
    return json.dumps(value, sort_keys=False, ensure_ascii=False, allow_nan=True)


# --------------------------------------------------------------------------- SVG
def svg_escape(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def svg_wrap(width: int, height: int, body: Sequence[str], title: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 '
        f'{width} {height}">'
        f'<title>{svg_escape(title)}</title><rect width="{width}" height="{height}" fill="white"/>'
        + "".join(body)
        + "</svg>"
    )


def nice_ticks(lo: float, hi: float, count: int = 5) -> list[float]:
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / count
    magnitude = 10 ** math.floor(math.log10(raw))
    step = min((s for s in (1, 2, 2.5, 5, 10) if s * magnitude >= raw), default=10) * magnitude
    first = math.ceil(lo / step) * step
    ticks = []
    value = first
    while value <= hi + 1e-9 * step:
        ticks.append(round(value, 12))
        value += step
    return ticks


class CandidateUnavailable(Exception):
    """Raised when the operation cannot run here (missing program, library or data)."""

    def __init__(self, requires: str, message: str) -> None:
        super().__init__(message)
        self.requires = requires
        self.message = message


def finish(
    module: str, operation: str, input: object, options: Mapping[str, object] | None, impl: Any
) -> dict[str, object]:
    """Run ``impl(input, opts) -> (result, warnings)`` and wrap it in the candidate envelope."""

    from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok

    opts = get_options(options)
    try:
        outcome = impl(input, opts)
        result, warnings = outcome if isinstance(outcome, tuple) else (outcome, [])
    except DNAKitError:
        raise
    except CandidateUnavailable as exc:
        return candidate_not_implemented(
            module, operation, input, options=opts, requires=exc.requires, message=exc.message
        )
    except (ArithmeticError, IndexError, KeyError, TypeError, ValueError) as exc:
        raise InputFormatError(
            f"{module}.{operation}: the input could not be processed "
            f"({type(exc).__name__}: {exc}).",
            code="INVALID_INPUT",
            hint="Check the input structure against the function documentation.",
        ) from exc
    return candidate_ok(module, operation, input, result, options=opts, warnings=warnings)
