"""Review wrappers for gtars LOLA universe helpers, region-token files, and VRS locations.

Each public function reads local BED or token files, runs one upstream gtars operation, and
returns a JSON-serializable dictionary with ``status``, ``operation``, result fields, and
``provenance``. Region coordinates are 0-based half-open (BED). Results are sorted so repeated
runs return identical lists.

Functions are importable as ``dnakit.annotation.gtars_review_functions.<name>``.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_COMMENT_PREFIXES = ("#", "track", "browser")
_BED_CONVENTION = "0-based half-open coordinates (BED)"
_MAX_TOKEN = 2**32 - 1


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int = 512) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _new_output(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"Output file already exists and will not be overwritten: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {path.parent}")
    return path


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _provenance(convention: str) -> dict[str, str]:
    return {"package": _package("gtars"), "convention": convention}


def _regions(path: Path) -> list[tuple[str, int, int]]:
    """Parse BED3 intervals for validation; upstream LOLA helpers do the computation."""

    regions: list[tuple[str, int, int]] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            if not text or text.startswith(_COMMENT_PREFIXES):
                continue
            fields = text.split("\t")
            if len(fields) < 3:
                raise ConfigurationError(f"{path}:{number} needs at least three BED columns.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ConfigurationError(f"{path}:{number} has non-integer coordinates.") from exc
            if start < 0 or end <= start:
                raise ConfigurationError(f"{path}:{number} must satisfy 0 <= start < end.")
            regions.append((fields[0], start, end))
    if not regions:
        raise ConfigurationError(f"{path} contains no intervals.")
    return regions


def _user_sets(user_regions: Sequence[str | Path], limit: int) -> list[list[tuple[str, int, int]]]:
    if isinstance(user_regions, (str, bytes, Path)) or not isinstance(user_regions, Iterable):
        raise ConfigurationError("user_regions must be a list of BED files, one per user set.")
    paths = [input_path(item) for item in user_regions]
    if not paths:
        raise ConfigurationError("At least one user set BED file is required.")
    sets = [_regions(path) for path in paths]
    total = sum(len(item) for item in sets)
    if total > limit:
        raise ConfigurationError(f"User sets hold {total} regions; max_rows={limit} is the limit.")
    return sets


def _region_record(region: Sequence[Any]) -> dict[str, Any]:
    chrom, start, end = region
    return {"chrom": str(chrom), "start": int(start), "end": int(end)}


def _sorted_regions(regions: Iterable[Sequence[Any]]) -> list[dict[str, Any]]:
    records = [_region_record(region) for region in regions]
    return sorted(records, key=lambda item: (item["chrom"], item["start"], item["end"]))


# LOLA universe helpers.


def lola_restricted_universe(
    user_regions: Sequence[str | Path],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the restricted LOLA universe from user region BED files, splitting overlaps.

    Each BED file is one user set. The result covers only the bases the user sets touch, cut
    into non-overlapping pieces, so it can serve as the universe for ``lola_enrichment``.
    """

    limit = _rows_limit(max_rows)
    sets = _user_sets(user_regions, limit)
    lola = require("gtars.lola", "bio-gtars")
    universe = lola.build_restricted_universe(sets)
    regions = _sorted_regions(universe)
    _check_count(regions, limit, "Restricted universe")
    return {
        "status": "completed",
        "operation": "lola_restricted_universe",
        "user_set_count": len(sets),
        "region_count": len(regions),
        "regions": regions,
        "provenance": _provenance(_BED_CONVENTION + "; overlapping user regions are split"),
    }


def lola_check_universe(
    user_regions: Sequence[str | Path],
    universe: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check how well a universe BED file covers each user set, with LOLA's warnings."""

    limit = _rows_limit(max_rows)
    sets = _user_sets(user_regions, limit)
    universe_path = input_path(universe)
    background = _regions(universe_path)
    _check_count(background, limit, "Universe")
    lola = require("gtars.lola", "bio-gtars")
    try:
        report = lola.check_universe(sets, background)
    except Exception as exc:  # gtars rejects malformed region tuples with runtime errors
        raise ConfigurationError(f"gtars could not check the universe: {exc}") from exc
    columns = {key: list(values) for key, values in report.items() if key != "warnings"}
    per_set = []
    for index in range(len(sets)):
        per_set.append(
            {
                "user_set": index,
                "total_regions": int(columns["totalRegions"][index]),
                "regions_in_universe": int(columns["regionsInUniverse"][index]),
                "coverage": float(columns["coverage"][index]),
                "many_to_many": int(columns["manyToMany"][index]),
            }
        )
    return {
        "status": "completed",
        "operation": "lola_check_universe",
        "user_set_count": len(sets),
        "user_sets": per_set,
        "warnings": sorted(str(item) for item in report.get("warnings", [])),
        "provenance": _provenance(_BED_CONVENTION + "; coverage is the fraction of regions hit"),
    }


def lola_redefine_user_sets(
    user_regions: Sequence[str | Path],
    universe: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Redefine each user set as the universe regions it overlaps, removing many-to-many noise."""

    limit = _rows_limit(max_rows)
    sets = _user_sets(user_regions, limit)
    universe_path = input_path(universe)
    background = _regions(universe_path)
    _check_count(background, limit, "Universe")
    lola = require("gtars.lola", "bio-gtars")
    try:
        redefined = lola.redefine_user_sets(sets, background)
    except Exception as exc:  # gtars rejects malformed region tuples with runtime errors
        raise ConfigurationError(f"gtars could not redefine the user sets: {exc}") from exc
    user_sets = [
        {"user_set": index, "regions": _sorted_regions(regions)}
        for index, regions in enumerate(redefined)
    ]
    return {
        "status": "completed",
        "operation": "lola_redefine_user_sets",
        "user_set_count": len(user_sets),
        "user_sets": user_sets,
        "provenance": _provenance(
            _BED_CONVENTION + "; regions are the overlapping universe regions"
        ),
    }


def _check_count(rows: Sequence[Any], limit: int, what: str) -> None:
    if len(rows) > limit:
        raise ConfigurationError(f"{what} has {len(rows)} rows; max_rows={limit} is the limit.")


# Region-token files.


def read_region_tokens(
    path: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GTOK binary file of region token ids, in the order they were stored."""

    source = input_path(path)
    limit = _rows_limit(max_rows)
    utils = require("gtars.utils", "bio-gtars")
    try:
        tokens = [int(token) for token in utils.read_tokens_from_gtok(str(source))]
    except Exception as exc:  # gtars reports malformed GTOK files with runtime errors
        raise ConfigurationError(f"gtars could not read the GTOK file: {exc}") from exc
    _check_count(tokens, limit, "Token file")
    return {
        "status": "completed",
        "operation": "read_region_tokens",
        "token_count": len(tokens),
        "tokens": tokens,
        "provenance": _provenance("token ids are universe indexes, in stored order"),
    }


def write_region_tokens(
    tokens: Sequence[int],
    destination: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write region token ids to a new GTOK binary file for later reading."""

    if isinstance(tokens, (str, bytes, Mapping)) or not isinstance(tokens, Iterable):
        raise ConfigurationError("tokens must be a list of non-negative integers.")
    values = list(tokens)
    if not values or len(values) > _MAX_ROWS:
        raise ConfigurationError(f"tokens must hold between 1 and {_MAX_ROWS} ids.")
    checked = [_integer(item, "token", 0, _MAX_TOKEN) for item in values]
    target = _new_output(destination)
    utils = require("gtars.utils", "bio-gtars")
    try:
        utils.write_tokens_to_gtok(str(target), checked)
    except Exception as exc:  # gtars reports unwritable paths with runtime errors
        raise ConfigurationError(f"gtars could not write the GTOK file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_region_tokens",
        "token_count": len(checked),
        "destination": str(target),
        "provenance": _provenance("token ids are stored as given, in order, in GTOK format"),
    }


# VRS sequence locations.


def vrs_location_digest(
    sequence_digest_value: str,
    start: int,
    end: int,
    *,
    prefixed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the GA4GH VRS SequenceLocation digest (``ga4gh:SL.<digest>``) of one interval.

    With ``prefixed=False`` the bare digest without the ``ga4gh:SL.`` prefix is returned.
    """

    digest = _text(sequence_digest_value, "sequence_digest", maximum=256)
    first = _integer(start, "start", 0, 10**12)
    last = _integer(end, "end", first + 1, 10**12)
    flag = _flag(prefixed, "prefixed")
    vrs = require("gtars.vrs", "bio-gtars")
    try:
        bare = str(vrs.location_digest(digest, first, last))
    except Exception as exc:  # gtars rejects unknown digests with runtime errors
        raise ConfigurationError(f"gtars could not digest the location: {exc}") from exc
    identifier = f"ga4gh:SL.{bare}" if flag else bare
    return {
        "status": "completed",
        "operation": "vrs_location_digest",
        "identifier": identifier,
        "digest": bare,
        "prefixed": flag,
        "sequence_digest": digest,
        "start": first,
        "end": last,
        "provenance": _provenance("VRS SequenceLocation with 0-based start and exclusive end"),
    }


__all__ = [
    "lola_check_universe",
    "lola_redefine_user_sets",
    "lola_restricted_universe",
    "read_region_tokens",
    "vrs_location_digest",
    "write_region_tokens",
]

_register_public_verbose(globals(), tuple(__all__))
