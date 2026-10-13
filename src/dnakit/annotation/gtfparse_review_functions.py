"""GTF attribute-string expansion and missing transcript or gene rows through gtfparse (review
group G18g).

Functions read local GTF files or attribute strings only and write nothing. Coordinates are 1-based
inclusive as stored in GTF. gtfparse is imported only when a function runs; a missing package raises
BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _checked,
    _gtf_limit,
    _gtf_rows,
    _package,
    _plain,
    _rows_limit,
    _text,
    _text_list,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def expand_gtfparse_attribute_strings(
    attribute_strings: Sequence[str | None],
    *,
    quote_char: str = "'",
    usecols: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expand GTF attribute strings into one value list per key with gtfparse."""

    limit = _rows_limit(max_rows)
    if isinstance(attribute_strings, (str, bytes)) or not isinstance(attribute_strings, Iterable):
        raise ConfigurationError("attribute_strings must be a list of strings or None.")
    rows = list(attribute_strings)
    _checked(rows, limit, "attribute strings")
    for item in rows:
        if item is not None and not isinstance(item, str):
            raise ConfigurationError("each attribute_strings item must be a string or None.")
    quote = _text(quote_char, "quote_char", maximum=1)
    keys = None if usecols is None else list(_text_list(usecols, "usecols", maximum=10_000))
    attribute_parsing = require("gtfparse.attribute_parsing", "bio-gtfparse")
    try:
        parsed = attribute_parsing.expand_attribute_strings(rows, quote_char=quote, usecols=keys)
    except Exception as exc:  # gtfparse raises for values it cannot split
        raise ConfigurationError(f"gtfparse could not expand the attributes: {exc}") from exc
    columns = {
        str(name): [_plain(value) for value in values] for name, values in sorted(parsed.items())
    }
    return {
        "status": "completed",
        "operation": "expand_gtfparse_attribute_strings",
        "row_count": len(rows),
        "column_count": len(columns),
        "columns": columns,
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "one entry per input row in each key list; a repeated key in one row is "
            "joined with commas; a key absent from a row holds an empty string; None rows have "
            "no attributes",
        },
    }


def create_missing_gtfparse_features(
    source: str | Path,
    *,
    feature: str = "transcript",
    key: str = "transcript_id",
    extra_columns: Sequence[str] = (),
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create rows for a missing feature (transcript or gene) from exon rows with gtfparse."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    feature_name = _text(feature, "feature", maximum=64)
    key_name = _text(key, "key", maximum=64)
    extras = list(_text_list(extra_columns, "extra_columns", maximum=64))
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    creator = require("gtfparse.create_missing_features", "bio-gtfparse")
    try:
        frame = gtfparse.parse_gtf_and_expand_attributes(str(path))
    except Exception as exc:  # gtfparse raises ParsingError and pandas errors
        raise ConfigurationError(f"gtfparse could not read {path.name}: {exc}") from exc
    original = len(frame)
    try:
        complete = creator.create_missing_features(
            frame,
            unique_keys={feature_name: key_name},
            extra_columns={feature_name: extras},
        )
    except (KeyError, ValueError) as exc:  # gtfparse names the missing column in the message
        raise ConfigurationError(f"gtfparse could not create {feature_name} rows: {exc}") from exc
    records = _gtf_rows(complete.iloc[original:], limit, "created rows")
    return {
        "status": "completed",
        "operation": "create_missing_gtfparse_features",
        "feature": feature_name,
        "key": key_name,
        "created_count": len(records),
        "records": records,
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "one row per distinct key from the exon rows; start is the minimum and "
            "end the maximum of member rows; no rows are created when the feature already exists; "
            "missing values are None",
        },
    }


__all__ = [
    "create_missing_gtfparse_features",
    "expand_gtfparse_attribute_strings",
]

_register_public_verbose(globals(), tuple(__all__))
