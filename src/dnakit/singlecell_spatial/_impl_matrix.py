"""Shared count-matrix input parsing and option helpers for single-cell, spatial and bulk RNA functions.

Accepted forms of a count matrix (``parse_matrix``):

* a mapping ``{gene: [count, count, ...]}`` or ``{gene: {cell: count}}`` (rows are genes);
* a mapping ``{"genes": [...], "cells": [...], "counts": [[...], ...]}`` (``counts`` may also be
  called ``values`` or ``matrix``; rows are genes);
* a list of rows (list of numeric lists, or list of ``{column: value}`` mappings);
* CSV or TSV text (or a path to a ``.csv``/``.tsv``/``.txt`` file, optionally ``.gz``) whose first
  line is a header (first field = row-label header, other fields = column names) and whose first
  column holds the row names;
* a 2-D numpy array or a pandas DataFrame (anything exposing ``tolist()`` / ``columns`` + ``index``);
* a 10x Genomics directory holding ``matrix.mtx[.gz]``, ``features.tsv[.gz]`` (or ``genes.tsv``)
  and ``barcodes.tsv[.gz]``.

Which axis holds the genes is decided by ``options["orientation"]``: ``"genes_by_cells"`` (default;
rows are genes, columns are cells/samples) or ``"cells_by_genes"`` (rows are cells/samples). The
rule applies to every form except the 10x directory (always genes by cells). Names that are not
supplied are generated (``gene_1``..., ``cell_1``...) or taken from ``options["genes"]`` and
``options["cells"]``.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError

DEFAULT_MAX_ELEMENTS = 5_000_000
MAX_FILE_BYTES = 256 * 1024 * 1024
ORIENTATIONS = ("genes_by_cells", "cells_by_genes")


class MissingDependency(Exception):
    """An optional third-party package needed by the operation is not installed."""

    def __init__(self, package: str) -> None:
        super().__init__(package)
        self.package = package


def numpy_module():
    """Import numpy lazily; report the operation as not implemented when it is missing."""

    try:
        import numpy
    except ImportError as exc:  # pragma: no cover - numpy is present in the test environment
        raise MissingDependency("numpy") from exc
    return numpy


def scipy_optimize_module():
    """Import scipy.optimize lazily (used by NNLS deconvolution)."""

    try:
        from scipy import optimize
    except ImportError as exc:  # pragma: no cover - scipy is present in the test environment
        raise MissingDependency("scipy") from exc
    return optimize


# ---------------------------------------------------------------------------------------------
# option helpers
# ---------------------------------------------------------------------------------------------


def as_options(options: object) -> dict[str, object]:
    if options is None:
        return {}
    if not isinstance(options, Mapping):
        raise ConfigurationError(
            "options must be a mapping of option names to values.",
            code="INVALID_OPTIONS",
            hint="Pass a dict.",
        )
    return {str(key): value for key, value in options.items()}


def opt_int(
    options: Mapping[str, object],
    key: str,
    default: int | None,
    minimum: int = 0,
    maximum: int = 10**9,
) -> int | None:
    value = options.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"options[{key!r}] must be an integer in [{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options[{key!r}] to an integer between {minimum} and {maximum}.",
        )
    return value


def opt_float(
    options: Mapping[str, object],
    key: str,
    default: float | None,
    minimum: float = -math.inf,
    maximum: float = math.inf,
    *,
    exclusive_min: bool = False,
) -> float | None:
    value = options.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"options[{key!r}] must be a finite number; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options[{key!r}] to a number.",
        )
    number = float(value)
    if number < minimum or number > maximum or (exclusive_min and number == minimum):
        raise ConfigurationError(
            f"options[{key!r}] must be within {'(' if exclusive_min else '['}{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Choose a value inside the documented range for {key!r}.",
        )
    return number


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"options[{key!r}] must be true or false; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options[{key!r}] to a boolean.",
        )
    return value


def opt_choice(
    options: Mapping[str, object], key: str, default: str, allowed: Sequence[str]
) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(
            f"options[{key!r}] must be one of {', '.join(allowed)}; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Use one of: {', '.join(allowed)}.",
        )
    return value


def opt_str(options: Mapping[str, object], key: str, default: str | None) -> str | None:
    value = options.get(key, default)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(
            f"options[{key!r}] must be a non-empty string; got {value!r}.",
            code="INVALID_OPTION",
            hint=f"Set options[{key!r}] to text.",
        )
    return value


def require_option(options: Mapping[str, object], key: str, what: str) -> object:
    if options.get(key) is None:
        raise ConfigurationError(
            f"options[{key!r}] is required: {what}.",
            code="MISSING_OPTION",
            hint=f"Pass options={{{key!r}: ...}}.",
        )
    return options[key]


# ---------------------------------------------------------------------------------------------
# matrix container and parsing
# ---------------------------------------------------------------------------------------------


@dataclass
class CountMatrix:
    """A dense matrix whose rows are genes and whose columns are cells (or samples/spots)."""

    genes: list[str]
    cells: list[str]
    values: list[list[float]]  # len(genes) rows, each len(cells)

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_cells(self) -> int:
        return len(self.cells)

    def cell_vectors(self) -> list[list[float]]:
        """Return one expression vector (over genes) per cell."""

        if not self.values:
            return [[] for _ in self.cells]
        return [list(column) for column in zip(*self.values, strict=False)]

    def column_totals(self) -> list[float]:
        return (
            [math.fsum(column) for column in zip(*self.values, strict=False)]
            if self.values
            else [0.0] * self.n_cells
        )

    def gene_index(self) -> dict[str, int]:
        return {gene: index for index, gene in enumerate(self.genes)}

    def cell_index(self) -> dict[str, int]:
        return {cell: index for index, cell in enumerate(self.cells)}

    def subset_genes(self, indices: Sequence[int]) -> CountMatrix:
        return CountMatrix(
            [self.genes[i] for i in indices], list(self.cells), [self.values[i] for i in indices]
        )

    def subset_cells(self, indices: Sequence[int]) -> CountMatrix:
        return CountMatrix(
            list(self.genes),
            [self.cells[j] for j in indices],
            [[row[j] for j in indices] for row in self.values],
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "orientation": "genes_by_cells",
            "genes": list(self.genes),
            "cells": list(self.cells),
            "values": self.values,
        }


def _number(token: object, where: str) -> float:
    if isinstance(token, bool):
        raise InputFormatError(
            f"{where}: boolean values are not counts.",
            code="INVALID_MATRIX_VALUE",
            hint="Use numbers.",
        )
    if isinstance(token, (int, float)):
        number = float(token)
    else:
        text = str(token).strip()
        try:
            number = float(text)
        except ValueError:
            raise InputFormatError(
                f"{where}: {text!r} is not a number.",
                code="INVALID_MATRIX_VALUE",
                hint="The matrix must contain only numeric values (header row and first column hold the names).",
            ) from None
    if not math.isfinite(number):
        raise InputFormatError(
            f"{where}: non-finite value {token!r}.",
            code="INVALID_MATRIX_VALUE",
            hint="Remove NaN/inf entries.",
        )
    return number


def _read_text(path: Path) -> str:
    if path.stat().st_size > MAX_FILE_BYTES:
        raise InputFormatError(
            f"{path.name} is larger than {MAX_FILE_BYTES} bytes.",
            code="INPUT_TOO_LARGE",
            hint="Subset the matrix first.",
        )
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read(MAX_FILE_BYTES + 1)
    return path.read_text(encoding="utf-8")


def _as_existing_path(value: object) -> Path | None:
    if isinstance(value, Path):
        return value
    if (
        isinstance(value, str)
        and "\n" not in value
        and 0 < len(value) < 2048
        and "\x00" not in value
    ):
        try:
            if os.path.exists(value):
                return Path(value)
        except OSError:
            return None
    return None


def text_of(value: object, what: str) -> str:
    """Return the text of a string input or of the file a path points to."""

    path = _as_existing_path(value)
    if path is not None:
        if not path.is_file():
            raise InputFormatError(
                f"{what}: {path} is not a file.",
                code="INPUT_NOT_FILE",
                hint="Pass a file path or the text itself.",
            )
        return _read_text(path)
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, str):
        return value
    raise InputFormatError(
        f"{what}: expected text or a file path.",
        code="INPUT_FORMAT",
        hint="Pass CSV/TSV text or a path.",
    )


def _delimited_rows(text: str, delimiter: str | None) -> list[list[str]]:
    lines = [line for line in text.splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        raise InputFormatError(
            "The matrix text is empty.",
            code="EMPTY_INPUT",
            hint="Provide a header line and at least one row.",
        )
    if delimiter is None:
        delimiter = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    if delimiter is None:
        return [line.split() for line in lines]
    return [row for row in csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter)]


def _from_delimited(
    text: str, delimiter: str | None
) -> tuple[list[str], list[str], list[list[float]]]:
    rows = _delimited_rows(text, delimiter)
    header = [item.strip() for item in rows[0]]
    columns = header[1:]
    if not columns:
        raise InputFormatError(
            "The header line needs a row-label column followed by column names.",
            code="MATRIX_HEADER",
            hint="Use 'name<TAB>col1<TAB>col2'.",
        )
    names: list[str] = []
    data: list[list[float]] = []
    for number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise InputFormatError(
                f"Line {number} has {len(row)} fields; the header has {len(header)}.",
                code="MATRIX_NOT_RECTANGULAR",
                hint="Every row needs a name and one value per column.",
            )
        names.append(row[0].strip())
        data.append([_number(item, f"line {number}") for item in row[1:]])
    if not data:
        raise InputFormatError(
            "The matrix has a header but no rows.",
            code="EMPTY_INPUT",
            hint="Add at least one data row.",
        )
    return names, columns, data


def _read_mtx(directory: Path, feature_column: str) -> CountMatrix:
    def find(*names: str) -> Path:
        for name in names:
            for candidate in (directory / name, directory / (name + ".gz")):
                if candidate.is_file():
                    return candidate
        raise InputFormatError(
            f"{directory} lacks {' / '.join(names)}.",
            code="MTX_FILE_MISSING",
            hint="A 10x directory needs matrix.mtx, features.tsv (or genes.tsv) and barcodes.tsv (optionally .gz).",
        )

    matrix_path = find("matrix.mtx")
    feature_path = find("features.tsv", "genes.tsv")
    barcode_path = find("barcodes.tsv")
    features = [line.split("\t") for line in _read_text(feature_path).splitlines() if line.strip()]
    barcodes = [
        line.split("\t")[0].strip()
        for line in _read_text(barcode_path).splitlines()
        if line.strip()
    ]
    use_name = feature_column == "name" and all(len(row) >= 2 for row in features)
    genes = [row[1].strip() if use_name else row[0].strip() for row in features]
    lines = iter(_read_text(matrix_path).splitlines())
    header = next(lines, "")
    if not header.startswith("%%MatrixMarket") or "coordinate" not in header:
        raise InputFormatError(
            "matrix.mtx must be a MatrixMarket coordinate file.",
            code="MTX_FORMAT",
            hint="Check the %%MatrixMarket header.",
        )
    dims = None
    for line in lines:
        if line.startswith("%") or not line.strip():
            continue
        dims = line.split()
        break
    if dims is None or len(dims) != 3:
        raise InputFormatError(
            "matrix.mtx lacks the 'rows cols entries' line.",
            code="MTX_FORMAT",
            hint="Re-export the matrix.",
        )
    n_rows, n_cols, n_entries = (int(item) for item in dims)
    if n_rows != len(genes) or n_cols != len(barcodes):
        raise InputFormatError(
            f"matrix.mtx is {n_rows}x{n_cols} but there are {len(genes)} features and {len(barcodes)} barcodes.",
            code="MTX_DIMENSIONS",
            hint="Use the features/barcodes files that belong to this matrix.",
        )
    if n_rows * n_cols > DEFAULT_MAX_ELEMENTS:
        raise InputFormatError(
            f"The matrix has {n_rows * n_cols} elements, above the {DEFAULT_MAX_ELEMENTS} limit.",
            code="INPUT_TOO_LARGE",
            hint="Subset the genes/cells before loading.",
        )
    values = [[0.0] * n_cols for _ in range(n_rows)]
    seen = 0
    for line in lines:
        if not line.strip() or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) < 3:
            raise InputFormatError(
                f"Bad MatrixMarket entry: {line!r}",
                code="MTX_FORMAT",
                hint="Entries are 'row col value'.",
            )
        row, col = int(parts[0]) - 1, int(parts[1]) - 1
        if not (0 <= row < n_rows and 0 <= col < n_cols):
            raise InputFormatError(
                f"MatrixMarket index out of range: {line!r}",
                code="MTX_FORMAT",
                hint="Indices are 1-based.",
            )
        values[row][col] += _number(parts[2], "matrix.mtx")
        seen += 1
    if seen != n_entries:
        raise InputFormatError(
            f"matrix.mtx declares {n_entries} entries but has {seen}.",
            code="MTX_FORMAT",
            hint="The file is truncated.",
        )
    return CountMatrix(genes, barcodes, values)


def _names(
    options: Mapping[str, object], key: str, expected: int, prefix: str, given: Sequence[str] | None
) -> list[str]:
    explicit = options.get(key)
    if explicit is not None:
        if isinstance(explicit, (str, bytes)) or not isinstance(explicit, Sequence):
            raise ConfigurationError(
                f"options[{key!r}] must be a list of names.",
                code="INVALID_OPTION",
                hint="Pass a list of strings.",
            )
        names = [str(item) for item in explicit]
        if len(names) != expected:
            raise InputFormatError(
                f"options[{key!r}] has {len(names)} names but the matrix has {expected}.",
                code="NAME_COUNT_MISMATCH",
                hint="Provide exactly one name per row/column.",
            )
        return names
    if given is not None:
        return list(given)
    return [f"{prefix}_{index}" for index in range(1, expected + 1)]


def _column_names(
    options: Mapping[str, object], mappings: Sequence[Mapping[object, object]]
) -> list[str]:
    """Columns of a sparse ``{row: {column: value}}`` input: the names given in options, else first-appearance order."""

    key = "cells" if options.get("orientation", "genes_by_cells") == "genes_by_cells" else "genes"
    explicit = options.get(key)
    if (
        explicit is not None
        and not isinstance(explicit, (str, bytes))
        and isinstance(explicit, Sequence)
    ):
        return [str(name) for name in explicit]
    columns: list[str] = []
    seen: set[str] = set()
    for item in mappings:
        for column in item:
            if str(column) not in seen:
                seen.add(str(column))
                columns.append(str(column))
    return columns


def _table_from_object(
    value: object, options: Mapping[str, object]
) -> tuple[list[str] | None, list[str] | None, list[list[float]]]:
    """Return (row names, column names, rows) for any supported non-text container."""

    if (
        hasattr(value, "columns") and hasattr(value, "index") and hasattr(value, "values")
    ):  # pandas DataFrame
        rows = [[_number(item, "DataFrame") for item in row] for row in value.values.tolist()]
        return [str(i) for i in value.index], [str(c) for c in value.columns], rows
    if (
        hasattr(value, "shape") and hasattr(value, "tolist") and not isinstance(value, (str, bytes))
    ):  # numpy array
        listed = value.tolist()
        if not isinstance(listed, list) or (listed and not isinstance(listed[0], list)):
            raise InputFormatError(
                "The array must be two-dimensional.",
                code="MATRIX_DIMENSIONS",
                hint="Pass a 2-D array.",
            )
        return None, None, [[_number(item, "array") for item in row] for row in listed]
    if isinstance(value, Mapping):
        keys = {str(key) for key in value}
        for counts_key in ("counts", "values", "matrix"):
            if "genes" in keys and counts_key in keys and isinstance(value[counts_key], Sequence):
                body = value[counts_key]
                rows = [
                    [_number(item, f"{counts_key}[{i}]") for item in row]
                    for i, row in enumerate(body)
                ]
                cells = value.get("cells")
                return (
                    [str(g) for g in value["genes"]],
                    ([str(c) for c in cells] if cells is not None else None),
                    rows,
                )
        items = list(value.items())
        if not items:
            raise InputFormatError(
                "The matrix mapping is empty.",
                code="EMPTY_INPUT",
                hint="Provide gene -> counts entries.",
            )
        if all(isinstance(item, Mapping) for _, item in items):
            columns = _column_names(options, [item for _, item in items])
            keyed = [{str(k): v for k, v in item.items()} for _, item in items]
            rows = [
                [_number(item.get(column, 0), str(name)) for column in columns]
                for (name, _), item in zip(items, keyed, strict=False)
            ]
            return [str(name) for name, _ in items], columns, rows
        rows = []
        for name, item in items:
            if isinstance(item, (str, bytes)) or not isinstance(item, Sequence):
                raise InputFormatError(
                    f"Row {name!r} must be a list of numbers or a {{column: number}} mapping.",
                    code="MATRIX_ROW",
                    hint="Use {gene: [counts...]} or {gene: {cell: count}}.",
                )
            rows.append([_number(x, str(name)) for x in item])
        return [str(name) for name, _ in items], None, rows
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if not value:
            raise InputFormatError(
                "The matrix has no rows.", code="EMPTY_INPUT", hint="Provide at least one row."
            )
        if all(isinstance(row, Mapping) for row in value):
            columns = _column_names(options, list(value))
            keyed_rows = [{str(k): v for k, v in row.items()} for row in value]
            rows = [
                [_number(row.get(column, 0), "row") for column in columns] for row in keyed_rows
            ]
            return None, columns, rows
        rows = []
        for index, row in enumerate(value):
            if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
                raise InputFormatError(
                    f"Row {index} is not a list of numbers.",
                    code="MATRIX_ROW",
                    hint="Pass a list of numeric lists.",
                )
            rows.append([_number(x, f"row {index}") for x in row])
        return None, None, rows
    raise InputFormatError(
        f"Unsupported matrix input of type {type(value).__name__}.",
        code="UNSUPPORTED_MATRIX_INPUT",
        hint="Use a mapping, list of rows, CSV/TSV text or path, numpy array, DataFrame or 10x directory.",
    )


def parse_matrix(
    value: object,
    options: Mapping[str, object] | None = None,
    *,
    allow_negative: bool = False,
    check_names: bool = True,
    default_orientation: str = "genes_by_cells",
    max_elements: int | None = None,
) -> CountMatrix:
    """Parse any supported matrix form (see the module docstring) into a genes-by-cells matrix."""

    options = options or {}
    if value is None:
        raise InputFormatError(
            "input is required.", code="EMPTY_INPUT", hint="Provide a count matrix."
        )
    orientation = opt_choice(options, "orientation", default_orientation, ORIENTATIONS)
    limit = (
        max_elements
        if max_elements is not None
        else opt_int(options, "max_matrix_elements", DEFAULT_MAX_ELEMENTS, 1, 10**9)
    )

    path = _as_existing_path(value)
    row_names: list[str] | None
    col_names: list[str] | None
    if path is not None and path.is_dir():
        matrix = _read_mtx(path, opt_choice(options, "feature_column", "name", ("name", "id")))
        row_names = col_names = None
        rows = matrix.values
        orientation = "genes_by_cells"
        row_names, col_names = matrix.genes, matrix.cells
    elif path is not None or isinstance(value, (str, bytes)):
        delimiter = options.get("delimiter")
        if delimiter is not None and (not isinstance(delimiter, str) or len(delimiter) != 1):
            raise ConfigurationError(
                "options['delimiter'] must be one character.",
                code="INVALID_OPTION",
                hint="Use ',' or a tab.",
            )
        row_names, col_names, rows = _from_delimited(text_of(value, "matrix"), delimiter)
    else:
        row_names, col_names, rows = _table_from_object(value, options)

    if not rows or not rows[0]:
        raise InputFormatError(
            "The matrix has no values.",
            code="EMPTY_INPUT",
            hint="Provide at least one row and column.",
        )
    width = len(rows[0])
    for index, row in enumerate(rows):
        if len(row) != width:
            raise InputFormatError(
                f"Row {index} has {len(row)} values; expected {width}.",
                code="MATRIX_NOT_RECTANGULAR",
                hint="Every row must have the same number of values.",
            )
    if len(rows) * width > limit:
        raise InputFormatError(
            f"The matrix has {len(rows) * width} elements, above max_matrix_elements={limit}.",
            code="INPUT_TOO_LARGE",
            hint="Subset genes/cells or raise options['max_matrix_elements'].",
        )
    if not allow_negative and any(x < 0 for row in rows for x in row):
        raise InputFormatError(
            "The matrix contains negative values.",
            code="NEGATIVE_VALUES",
            hint="Counts must be non-negative.",
        )

    if orientation == "genes_by_cells":
        genes = _names(options, "genes", len(rows), "gene", row_names)
        cells = _names(options, "cells", width, "cell", col_names)
        values = [list(row) for row in rows]
    else:
        cells = _names(options, "cells", len(rows), "cell", row_names)
        genes = _names(options, "genes", width, "gene", col_names)
        values = [list(column) for column in zip(*rows, strict=False)]
    if check_names:
        for what, names in (("gene", genes), ("cell", cells)):
            if len(set(names)) != len(names):
                duplicate = next(name for name in names if names.count(name) > 1)
                raise InputFormatError(
                    f"Duplicate {what} name {duplicate!r}.",
                    code="DUPLICATE_NAMES",
                    hint=f"Make {what} names unique (or aggregate duplicates before calling).",
                )
    return CountMatrix(genes, cells, values)


# ---------------------------------------------------------------------------------------------
# labels, coordinates, shared small helpers
# ---------------------------------------------------------------------------------------------


def labels_for(cells: Sequence[str], raw: object, what: str) -> list[str]:
    """Resolve per-cell labels from a cell->label mapping or a list aligned with the cells."""

    if isinstance(raw, Mapping):
        mapping = {str(key): str(item) for key, item in raw.items()}
        missing = [cell for cell in cells if cell not in mapping]
        if missing:
            raise InputFormatError(
                f"{what} has no label for {len(missing)} cells (first: {missing[0]!r}).",
                code="MISSING_LABELS",
                hint="Provide a label for every cell.",
            )
        return [mapping[cell] for cell in cells]
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
        if len(raw) != len(cells):
            raise InputFormatError(
                f"{what} has {len(raw)} labels for {len(cells)} cells.",
                code="LABEL_COUNT_MISMATCH",
                hint="Give exactly one label per cell, in matrix column order.",
            )
        return [str(item) for item in raw]
    raise ConfigurationError(
        f"{what} must be a mapping cell -> label or a list of labels.",
        code="INVALID_OPTION",
        hint="Pass a dict or list.",
    )


def option_labels(cm: CountMatrix, options: Mapping[str, object], key: str, what: str) -> list[str]:
    return labels_for(cm.cells, require_option(options, key, what), f"options[{key!r}]")


def parse_coordinates(
    value: object, cells: Sequence[str], options: Mapping[str, object]
) -> list[tuple[float, float]]:
    """Resolve 2-D spatial coordinates for ``cells``.

    Accepted: a mapping ``cell -> [x, y]``; a list/array of ``[x, y]`` aligned with the cells; CSV/TSV
    text or a file with a header whose first column is the barcode and that has ``x``/``y`` columns
    (or ``pxl_col_in_fullres``/``pxl_row_in_fullres``, ``array_col``/``array_row``); or a headerless
    10x Visium ``tissue_positions_list.csv`` (barcode, in_tissue, array_row, array_col, pxl_row, pxl_col).
    ``options['x_column']`` / ``options['y_column']`` pick other column names.
    """

    if value is None:
        raise ConfigurationError(
            "options['coordinates'] is required for spatial operations.",
            code="MISSING_OPTION",
            hint="Pass a mapping cell -> [x, y], a list of [x, y], or a tissue-positions CSV.",
        )
    mapping: dict[str, tuple[float, float]] = {}
    path = _as_existing_path(value)
    if path is not None or isinstance(value, (str, bytes)):
        rows = _delimited_rows(text_of(value, "coordinates"), None)
        first = [item.strip() for item in rows[0]]
        headerless = len(first) >= 6 and first[1] in {"0", "1"}
        if headerless:
            body = rows
            x_i, y_i = 5, 4  # pxl_col (x), pxl_row (y)
        else:
            body = rows[1:]
            lowered = [item.casefold() for item in first]
            x_name = str(options.get("x_column", "")).casefold()
            y_name = str(options.get("y_column", "")).casefold()

            def pick(preferred: str, candidates: Sequence[str]) -> int:
                for name in ([preferred] if preferred else []) + list(candidates):
                    if name and name in lowered:
                        return lowered.index(name)
                raise InputFormatError(
                    f"No coordinate column found in header {first}.",
                    code="COORDINATE_COLUMNS",
                    hint="Name the columns x/y (or pxl_col_in_fullres/pxl_row_in_fullres) or set options['x_column'/'y_column'].",
                )

            x_i = pick(x_name, ("x", "pxl_col_in_fullres", "imagecol", "array_col"))
            y_i = pick(y_name, ("y", "pxl_row_in_fullres", "imagerow", "array_row"))
        for row in body:
            if len(row) <= max(x_i, y_i):
                raise InputFormatError(
                    "A coordinate row is too short.",
                    code="COORDINATE_ROW",
                    hint="Check the delimiter.",
                )
            mapping[row[0].strip()] = (
                _number(row[x_i], "coordinates"),
                _number(row[y_i], "coordinates"),
            )
    elif isinstance(value, Mapping):
        for key, item in value.items():
            if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) < 2:
                raise InputFormatError(
                    f"Coordinates of {key!r} must be [x, y].",
                    code="COORDINATE_ROW",
                    hint="Use two numbers per cell.",
                )
            mapping[str(key)] = (_number(item[0], str(key)), _number(item[1], str(key)))
    else:
        _, _, rows_list = _table_from_object(value, options)
        if len(rows_list) != len(cells) or any(len(row) < 2 for row in rows_list):
            raise InputFormatError(
                f"coordinates has {len(rows_list)} rows for {len(cells)} cells.",
                code="COORDINATE_COUNT",
                hint="Provide one [x, y] per cell in matrix column order.",
            )
        return [(row[0], row[1]) for row in rows_list]
    missing = [cell for cell in cells if cell not in mapping]
    if missing:
        raise InputFormatError(
            f"{len(missing)} cells have no coordinates (first: {missing[0]!r}).",
            code="MISSING_COORDINATES",
            hint="Provide coordinates for every cell/spot in the matrix.",
        )
    return [mapping[cell] for cell in cells]


def write_matrix_tsv(cm: CountMatrix, options: Mapping[str, object]) -> str | None:
    """Write ``cm`` as TSV only when ``options['output_path']`` is given (never overwrites)."""

    target = opt_str(options, "output_path", None)
    if target is None:
        return None
    path = Path(target).expanduser()
    if path.exists():
        raise ConfigurationError(
            f"Output file already exists: {path}",
            code="OUTPUT_EXISTS",
            hint="Choose a new output_path.",
        )
    if not path.parent.is_dir():
        raise ConfigurationError(
            f"Output directory does not exist: {path.parent}",
            code="OUTPUT_DIR",
            hint="Create the directory first.",
        )
    lines = ["gene\t" + "\t".join(cm.cells)]
    for gene, row in zip(cm.genes, cm.values, strict=False):
        lines.append(gene + "\t" + "\t".join(repr(x) for x in row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def matrix_result(
    cm: CountMatrix, options: Mapping[str, object], extra: Mapping[str, object] | None = None
) -> dict[str, object]:
    """Standard result block for functions that return a transformed matrix."""

    result: dict[str, object] = dict(extra or {})
    result["n_genes"] = cm.n_genes
    result["n_cells"] = cm.n_cells
    if opt_bool(options, "return_matrix", True):
        result["matrix"] = cm.as_dict()
    result["output_path"] = write_matrix_tsv(cm, options)
    return result


# ---------------------------------------------------------------------------------------------
# envelope runner
# ---------------------------------------------------------------------------------------------

Implementation = Callable[
    [object, dict[str, object]], "tuple[object, list[str]] | dict[str, object]"
]


def run(
    module: str, operation: str, input: object, options: object, implementation: Implementation
) -> dict[str, object]:
    """Run ``implementation(input, options)`` and wrap it in the DNAKit candidate envelope."""

    opts = as_options(options)
    try:
        outcome = implementation(input, opts)
    except MissingDependency as exc:
        return candidate_not_implemented(
            module,
            operation,
            input,
            options=opts,
            requires=f"pip install {exc.package}",
            message=f"{module}.{operation} needs the optional package {exc.package}, which is not installed.",
        )
    if isinstance(
        outcome, dict
    ):  # the implementation already built an envelope (e.g. external program missing)
        return outcome
    result, warnings = outcome
    return candidate_ok(module, operation, input, result, options=opts, warnings=warnings)
