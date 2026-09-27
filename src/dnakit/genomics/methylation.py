"""Strict modkit bedMethyl input and descriptive, coverage-weighted summaries.

Format reference: https://nanoporetech.github.io/modkit/intro_pileup.html
This independent parser does not redistribute modkit. Column 11 is a percentage,
despite the abbreviated README calling it a fraction. No statistical testing or
disease prediction is performed here.
"""

from __future__ import annotations

import gzip
import math
import re
import zlib
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._progress import track
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import ConfigurationError, InputFormatError

_COUNTS = (
    "valid_coverage",
    "modified_count",
    "canonical_count",
    "other_modified_count",
    "deletion_count",
    "failed_count",
    "different_base_count",
    "no_call_count",
)
_INTEGERS = ("start", "end", "score", "thick_start", "thick_end", *_COUNTS)
_CODE = re.compile(r"(?:[A-Za-z]|[0-9]+)\Z")
_MOTIF = re.compile(r"[ACGTRYSWKMBDHVN]+\Z")


def _positive(value: int, name: str) -> None:
    if type(value) is not int or value <= 0:
        raise ConfigurationError(f"{name} must be a positive integer.")


def _invalid(message: str, line: int) -> InputFormatError:
    return InputFormatError(
        f"bedMethyl record at line/index {line}: {message}",
        code="INVALID_BEDMETHYL",
        context={"line": line},
    )


def _validate(record: Mapping[str, Any], line: int) -> None:
    for key in _INTEGERS:
        if type(record.get(key)) is not int or record[key] < 0:
            raise _invalid(f"{key} must be a non-negative integer.", line)
    chrom = record.get("chromosome")
    if not isinstance(chrom, str) or not chrom or any(c.isspace() for c in chrom):
        raise _invalid("chromosome must be non-empty text without whitespace.", line)
    code = record.get("modification_code")
    if not isinstance(code, str) or _CODE.fullmatch(code) is None:
        raise _invalid("modification code must be a letter or numeric ChEBI code.", line)
    if record.get("strand") not in ("+", "-", "."):
        raise _invalid("strand must be +, -, or . (combined strands).", line)
    if not (record["start"] <= record["thick_start"] <= record["thick_end"] <= record["end"]):
        raise _invalid("compatibility coordinates must lie within the BED interval.", line)
    if record["start"] >= record["end"]:
        raise _invalid("coordinates must define a non-empty 0-based half-open interval.", line)
    coverage = record["valid_coverage"]
    if coverage != sum(record[key] for key in _COUNTS[1:4]):
        raise _invalid("valid coverage must equal modified + canonical + other modified.", line)
    if record["score"] != coverage:
        raise _invalid("score must equal valid coverage for modkit bedMethyl.", line)
    percent = record.get("percent_modified")
    if (
        isinstance(percent, bool)
        or not isinstance(percent, (int, float))
        or not 0 <= percent <= 100
        or not math.isfinite(percent)
    ):
        raise _invalid("percent_modified must be finite and between 0 and 100.", line)
    expected = 100 * record["modified_count"] / coverage if coverage else 0.0
    # Modkit writes percentages rounded to two decimal places. Counts, rather
    # than that rounded field, remain authoritative for every summary.
    if not math.isclose(percent, expected, rel_tol=0, abs_tol=0.005001):
        raise _invalid("percentage disagrees with counts (expected 0..100 units).", line)
    motif, offset = record.get("motif"), record.get("motif_offset")
    if motif is None:
        if offset is not None:
            raise _invalid("motif_offset requires a motif.", line)
    elif (
        not isinstance(motif, str)
        or _MOTIF.fullmatch(motif) is None
        or type(offset) is not int
        or not 0 <= offset < len(motif)
    ):
        raise _invalid("motif must use IUPAC DNA and contain the 0-based motif offset.", line)


def _parse(text: str, line: int) -> dict[str, Any]:
    fields = text.split()  # Older modkit mixes tabs (BED9) and spaces (counts).
    if len(fields) != 18:
        raise _invalid(f"expected 18 columns, found {len(fields)}.", line)
    integers = {
        1: "start",
        2: "end",
        4: "score",
        6: "thick_start",
        7: "thick_end",
        9: "valid_coverage",
        **{i: key for i, key in enumerate(_COUNTS[1:], 11)},
    }
    record: dict[str, Any] = {
        "chromosome": fields[0],
        "name": fields[3],
        "strand": fields[5],
        "color": fields[8],
    }
    for index, key in integers.items():
        if re.fullmatch(r"[0-9]+", fields[index]) is None:
            raise _invalid(f"{key} must be a non-negative integer.", line)
        try:
            record[key] = int(fields[index])
        except ValueError as exc:
            raise _invalid(f"{key} integer is too large.", line) from exc
    try:
        record["percent_modified"] = float(fields[10])
    except ValueError as exc:
        raise _invalid("invalid percentage.", line) from exc
    colors = fields[8].split(",")
    if (
        len(colors) != 3
        or any(re.fullmatch(r"[0-9]{1,3}", c) is None for c in colors)
        or any(int(c) > 255 for c in colors)
    ):
        raise _invalid("color must contain three RGB integers between 0 and 255.", line)
    parts = fields[3].split(",")
    if len(parts) not in (1, 3):
        raise _invalid("name must be a modification code or code,motif,offset.", line)
    record["modification_code"] = parts[0]
    record["motif"] = parts[1] if len(parts) == 3 else None
    record["motif_offset"] = None
    if len(parts) == 3:
        if re.fullmatch(r"[0-9]{1,6}", parts[2]) is None:
            raise _invalid("invalid motif offset.", line)
        record["motif_offset"] = int(parts[2])
    _validate(record, line)
    record["fraction_modified"] = (
        record["modified_count"] / record["valid_coverage"] if record["valid_coverage"] else None
    )
    return record


def read_bedmethyl(
    path: str | Path,
    *,
    max_records: int = 1_000_000,
    max_bytes: int = 268_435_456,
    max_line_bytes: int = 65_536,
    verbose: bool = True,
) -> list[dict[str, Any]]:
    """Read modkit's 18-column BED, detecting gzip/BGZF by its magic bytes.

    All limits apply to decompressed input. Exceeding a limit raises, never
    silently truncates. Empty lines, # comments, track and browser headers are
    ignored; malformed records are rejected with a physical line number.
    Coordinates are 0-based, end-exclusive. Column 11 is percent (0..100);
    fraction_modified is calculated from counts and is None at zero coverage.
    bedRMod has a different denominator and is deliberately unsupported.
    """
    for name, value in (
        ("max_records", max_records),
        ("max_bytes", max_bytes),
        ("max_line_bytes", max_line_bytes),
    ):
        _positive(value, name)
    if not isinstance(path, (str, Path)):
        raise ConfigurationError("path must be a filesystem path string or Path.")
    source = Path(path)
    records: list[dict[str, Any]] = []

    def lines() -> Iterable[tuple[int, bytes]]:
        with source.open("rb") as raw:
            compressed = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            stream = gzip.GzipFile(fileobj=raw) if compressed else raw
            try:
                consumed = 0
                line_number = 0
                while True:
                    data = stream.readline(min(max_line_bytes, max_bytes - consumed) + 1)
                    if not data:
                        break
                    line_number += 1
                    consumed += len(data)
                    if len(data) > max_line_bytes or consumed > max_bytes:
                        raise _invalid(
                            "decompressed byte or line-size limit exceeded.", line_number
                        )
                    yield line_number, data
            finally:
                if compressed:
                    stream.close()

    try:
        for line_number, data in track(lines(), description="Read bedMethyl", disable=not verbose):
            try:
                text = data.decode("utf-8").strip()
            except UnicodeDecodeError as exc:
                raise _invalid("input must be UTF-8 text.", line_number) from exc
            if text.startswith("#") and "bedrmod" in text.lower():
                raise _invalid("bedRMod is not the supported bedMethyl format.", line_number)
            if not text or text.startswith(("#", "track ", "browser ")):
                continue
            if len(records) >= max_records:
                raise _invalid("max_records exceeded; no partial result returned.", line_number)
            records.append(_parse(text, line_number))
    except (gzip.BadGzipFile, EOFError, zlib.error) as exc:
        raise InputFormatError("Invalid or truncated gzip bedMethyl input.") from exc
    return records


def _filter(values: Sequence[str] | None, name: str) -> frozenset[str] | None:
    if values is None:
        return None
    if (
        not isinstance(values, Sequence)
        or isinstance(values, (str, bytes))
        or any(not isinstance(v, str) or not v for v in values)
    ):
        raise ConfigurationError(f"{name} must be a sequence of non-empty strings.")
    return frozenset(values)


def summarize_methylation(
    records: Iterable[Mapping[str, Any]],
    *,
    modification_codes: Sequence[str] | None = None,
    chromosomes: Sequence[str] | None = None,
    strands: Sequence[str] | None = None,
    min_coverage: int = 0,
    max_records: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize validated records from read_bedmethyl using sum(mod)/sum(valid).

    Groups always preserve chromosome, modification code, strand and motif.
    In particular, 5mC (m), 5hmC (h), and combined strands (.) remain separate.
    Duplicated/overlapping records are counted as supplied; record_count is not
    a count of unique CpGs. No cross-sample tests or disease inference is made.
    Empty filter sequences explicitly select nothing. Every input is validated
    before filtering, and undefined zero-coverage fractions remain None.
    """
    _positive(max_records, "max_records")
    if isinstance(records, (str, bytes, Mapping)):
        raise ConfigurationError("records must be an iterable of bedMethyl record mappings.")
    try:
        iterator = iter(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be an iterable of bedMethyl record mappings."
        ) from exc
    if type(min_coverage) is not int or min_coverage < 0:
        raise ConfigurationError("min_coverage must be a non-negative integer.")
    codes = _filter(modification_codes, "modification_codes")
    chroms = _filter(chromosomes, "chromosomes")
    selected_strands = _filter(strands, "strands")
    if selected_strands is not None and not selected_strands <= {"+", "-", "."}:
        raise ConfigurationError("strands may contain only +, -, and .")
    grouped: dict[tuple[Any, ...], dict[str, Any]] = {}
    total = selected = 0
    for total, record in enumerate(
        track(iterator, description="Summarize methylation", disable=not verbose), 1
    ):
        if total > max_records:
            raise _invalid("max_records exceeded; no partial summary returned.", total)
        if not isinstance(record, Mapping):
            raise _invalid("expected a record mapping from read_bedmethyl.", total)
        _validate(record, total)
        if (
            (codes is not None and record["modification_code"] not in codes)
            or (chroms is not None and record["chromosome"] not in chroms)
            or (selected_strands is not None and record["strand"] not in selected_strands)
            or record["valid_coverage"] < min_coverage
        ):
            continue
        selected += 1
        fields = ("chromosome", "modification_code", "strand", "motif", "motif_offset")
        key = tuple(record.get(field) for field in fields)
        if key not in grouped:
            grouped[key] = {
                **dict(zip(fields, key, strict=True)),
                "record_count": 0,
                "zero_coverage_records": 0,
                **dict.fromkeys(_COUNTS, 0),
            }
        group = grouped[key]
        group["record_count"] += 1
        group["zero_coverage_records"] += record["valid_coverage"] == 0
        for field in _COUNTS:
            group[field] += record[field]
    groups = list(grouped.values())
    for group in groups:
        fraction = (
            group["modified_count"] / group["valid_coverage"] if group["valid_coverage"] else None
        )
        group["fraction_modified"] = fraction
        group["percent_modified"] = 100 * fraction if fraction is not None else None
    return {
        "coordinate_system": "0-based-half-open",
        "records_total": total,
        "records_selected": selected,
        "records_filtered": total - selected,
        "aggregation": "sum(modified_count)/sum(valid_coverage) within each group",
        "groups": groups,
    }


__all__ = ["read_bedmethyl", "summarize_methylation"]
_register_public_verbose(globals(), tuple(__all__))
