"""Shared input checks and result conversion for the scikit-bio P20 functions (private)."""

from __future__ import annotations

import importlib.metadata
import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

EXTRA = "bio-skibio"
PACKAGE = "scikit-bio"
MAX_ROWS = 5000
MAX_COLUMNS = 20000
MAX_FILE_BYTES = 256 * 1024 * 1024


def provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": PACKAGE, "package_version": installed, **notes}


def skbio() -> Any:
    return require("skbio", EXTRA)


def flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def number(
    value: object, name: str, *, minimum: float | None = None, maximum: float | None = None
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    if (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
        raise ConfigurationError(f"{name} is outside the allowed range [{minimum}, {maximum}].")
    return float(value)


def text(value: object, name: str, *, maximum: int = 10_000_000) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    return value


def input_path(value: object, name: str = "path") -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if not path.is_file():
        raise ConfigurationError(f"{name} does not exist: {path}")
    if path.stat().st_size > MAX_FILE_BYTES:
        raise ConfigurationError(f"{name} is larger than {MAX_FILE_BYTES} bytes.")
    return path


def output_path(value: object, name: str = "output_path") -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if path.is_dir():
        raise ConfigurationError(f"{name} is a directory: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def matrix(value: object, name: str, *, square: bool = False, nonnegative: bool = False) -> Any:
    """Validate a rectangular numeric matrix and return a float numpy array."""

    import numpy as np

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of numeric rows.")
    if len(value) > MAX_ROWS:
        raise ConfigurationError(f"{name} may have at most {MAX_ROWS} rows.")
    try:
        array = np.array(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular numeric matrix.") from exc
    if array.ndim != 2 or array.shape[1] == 0 or array.shape[1] > MAX_COLUMNS:
        raise ConfigurationError(
            f"{name} must be a rectangular matrix with 1 to {MAX_COLUMNS} columns."
        )
    if not np.isfinite(array).all():
        raise ConfigurationError(f"{name} must contain only finite numbers.")
    if square and array.shape[0] != array.shape[1]:
        raise ConfigurationError(f"{name} must be square.")
    if nonnegative and (array < 0).any():
        raise ConfigurationError(f"{name} cannot contain negative values.")
    return array


def labels(value: object | None, expected: int, name: str, prefix: str) -> list[str]:
    if value is None:
        return [f"{prefix}{i + 1}" for i in range(expected)]
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of strings.")
    result = [str(item) for item in value]
    if len(result) != expected or len(set(result)) != expected or not all(result):
        raise ConfigurationError(f"{name} must hold exactly {expected} unique non-empty strings.")
    return result


def dataframe(columns: object, sample_ids: Sequence[str], name: str) -> Any:
    """Build a pandas DataFrame from ``{column: [value per sample]}`` indexed by sample id."""

    import pandas as pd

    if not isinstance(columns, Mapping) or not columns:
        raise ConfigurationError(
            f"{name} must be a mapping of column name to one value per sample."
        )
    data: dict[str, list[Any]] = {}
    for key, values in columns.items():
        if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
            raise ConfigurationError(f"{name}[{key!r}] must be a list with one value per sample.")
        if len(values) != len(sample_ids):
            raise ConfigurationError(f"{name}[{key!r}] must have {len(sample_ids)} values.")
        data[str(key)] = list(values)
    return pd.DataFrame(data, index=list(sample_ids))


def distance_matrix(values: object, ids: object | None, name: str = "distance_matrix") -> Any:
    """Build a scikit-bio DistanceMatrix from a square symmetric list matrix."""

    package = skbio()
    array = matrix(values, name, square=True, nonnegative=True)
    names = labels(ids, array.shape[0], "ids", "S")
    try:
        return package.DistanceMatrix(array, ids=names)
    except Exception as exc:  # scikit-bio raises several DistanceMatrixError variants
        raise ConfigurationError(f"{name} is not a valid distance matrix: {exc}") from exc


def jsonable(value: Any) -> Any:
    """Convert numpy, pandas and scikit-bio values to plain JSON-serialisable structures."""

    import numpy as np
    import pandas as pd

    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number_ = float(value)
        return number_ if math.isfinite(number_) else None
    if isinstance(value, pd.DataFrame):
        return {
            "index": [jsonable(i) for i in value.index],
            "columns": [str(c) for c in value.columns],
            "data": [[jsonable(v) for v in row] for row in value.to_numpy().tolist()],
        }
    if isinstance(value, pd.Series):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return [jsonable(v) for v in value.tolist()]
    if isinstance(value, Mapping):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return sorted((jsonable(v) for v in value), key=str)
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return str(value)


def ordination_to_dict(result: Any) -> dict[str, Any]:
    """Summarise an OrdinationResults as plain lists."""

    return {
        "method": getattr(result, "short_method_name", None),
        "long_method_name": getattr(result, "long_method_name", None),
        "eigvals": jsonable(result.eigvals),
        "proportion_explained": jsonable(result.proportion_explained),
        "samples": jsonable(result.samples),
        "features": jsonable(result.features),
        "biplot_scores": jsonable(result.biplot_scores),
        "sample_constraints": jsonable(result.sample_constraints),
    }


def table_frame(
    counts: object,
    sample_ids: object | None,
    feature_ids: object | None,
    *,
    nonnegative: bool = True,
) -> Any:
    """Samples-by-features pandas DataFrame from a list matrix."""

    import pandas as pd

    array = matrix(counts, "counts", nonnegative=nonnegative)
    samples = labels(sample_ids, array.shape[0], "sample_ids", "S")
    features = labels(feature_ids, array.shape[1], "feature_ids", "F")
    return pd.DataFrame(array, index=samples, columns=features)


def grouping_list(value: object, expected: int, name: str = "grouping") -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != expected:
        raise ConfigurationError(
            f"{name} must be a list with one group label per sample ({expected})."
        )
    return [str(item) for item in value]


def frame_result(frame: Any) -> dict[str, Any]:
    """Row-wise records of a pandas DataFrame (index levels become columns)."""

    flat = frame.reset_index()
    return {
        "columns": [str(c) for c in flat.columns],
        "rows": [[jsonable(v) for v in row] for row in flat.to_numpy().tolist()],
    }
