"""Sourmash wrappers for variants that the existing sourmash wrappers do not reach.

These functions cover signature archives (zip and SQLite), signature JSON text, filtered signature
loading, abundance-weighted comparison and search, sequence bloom tree search and num-based
sketches. Sketches are plain mappings with ``ksize``, ``scaled`` (or ``num``) and ``hashes``;
abundance-weighted functions also read an optional ``abundances`` list that is parallel to
``hashes``.

sourmash is imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-sourmash]"`` hint.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.microbiome.sourmash_functions import (
    _HASH_LIMIT,
    _KSIZE_LIMIT,
    _MATRIX_LIMIT,
    _SOURMASH_EXTRA,
    _integer,
    _real,
    _sketch_fields,
    _version,
)

_ARCHIVE_SUFFIXES = (".sig.zip", ".zip", ".sqldb")
_MAX_SIGNATURES = 100_000
_MAX_JSON_CHARS = 200_000_000
_SEARCH_MODES = ("jaccard", "containment", "max_containment", "abundance")
_LOAD_EXPECTATIONS = ("many", "one", "query")


def _package() -> str:
    return f"sourmash {_version('sourmash')}"


def _sketch_list(value: object, name: str, minimum: int, maximum: int) -> Sequence[Any]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of sketch mappings.")
    if not minimum <= len(value) <= maximum:
        raise ConfigurationError(f"{name} must contain between {minimum} and {maximum} sketches.")
    return value


def _abundance_sketch(
    sourmash: Any, sketch: object, name: str, *, default_name: str
) -> tuple[Any, bool]:
    """Build a signature; abundances are kept when the sketch gives them."""

    identifier, ksize, scaled, _hashes = _sketch_fields(sketch, name)
    assert isinstance(sketch, Mapping)
    abundances = sketch.get("abundances")
    track = abundances is not None
    try:
        minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=track)
        if track:
            raw = sketch["hashes"]
            if (
                isinstance(abundances, (str, bytes))
                or not isinstance(abundances, Sequence)
                or len(abundances) != len(raw)
            ):
                raise ConfigurationError(f"{name} abundances must parallel its hashes.")
            table: dict[int, int] = {}
            for hash_value, count in zip(raw, abundances, strict=True):
                if isinstance(count, bool) or not isinstance(count, int) or count < 1:
                    raise ConfigurationError(f"{name} abundances must be positive integers.")
                if table.setdefault(int(hash_value), count) != count:
                    raise ConfigurationError(f"{name} repeats a hash with different abundances.")
            minhash.set_abundances(table)
            kept = len(minhash)
            expected = len(table)
        else:
            minhash.add_many(_hashes)
            kept = len(minhash)
            expected = len(set(_hashes))
    except (AttributeError, OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not build a MinHash: {exc}") from exc
    if kept != expected:
        raise ConfigurationError(f"{name} has hashes outside the scaled range of its sketch.")
    signature = sourmash.signature.SourmashSignature(minhash, name=identifier or default_name)
    return signature, track


def _signature_record(item: Any) -> dict[str, Any]:
    minhash = item.minhash
    if minhash.moltype != "DNA":
        raise ConfigurationError(
            f"signature {item.name!r} uses {minhash.moltype}; only DNA signatures are supported."
        )
    if int(minhash.scaled) == 0:
        raise ConfigurationError(
            f"signature {item.name!r} is num-based; only scaled signatures are supported."
        )
    pairs = sorted((int(key), int(count)) for key, count in minhash.hashes.items())
    record: dict[str, Any] = {
        "id": str(item.name),
        "ksize": int(minhash.ksize),
        "scaled": int(minhash.scaled),
        "hash_count": len(pairs),
        "hashes": [key for key, _count in pairs],
        "md5sum": str(item.md5sum()),
    }
    if minhash.track_abundance:
        record["abundances"] = [count for _key, count in pairs]
    return record


def write_sourmash_signature_archive(
    sketches: Sequence[Mapping[str, Any]],
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write sketches to a sourmash zip archive or SQLite database.

    The format follows the output name: ``.zip`` and ``.sig.zip`` give a zip collection and
    ``.sqldb`` gives a SQLite signature index. Both are written with
    ``sourmash.save_load.SaveSignaturesToLocation`` and can be read back with
    ``read_sourmash_signatures``. Abundances are kept when a sketch has an ``abundances`` list.
    The output must not exist yet.
    """

    items = _sketch_list(sketches, "sketches", 1, _MAX_SIGNATURES)
    path = Path(output).expanduser().resolve()
    suffix = next((s for s in _ARCHIVE_SUFFIXES if path.name.endswith(s)), None)
    if suffix is None:
        raise ConfigurationError(f"output must end with one of: {', '.join(_ARCHIVE_SUFFIXES)}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"output already exists and will not be overwritten: {path}")
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    save_load = require("sourmash.save_load", _SOURMASH_EXTRA)
    signatures = [
        _abundance_sketch(sourmash, item, f"sketches[{index}]", default_name=f"sketch_{index}")[0]
        for index, item in enumerate(items)
    ]
    try:
        with save_load.SaveSignaturesToLocation(str(path)) as saver:
            for signature in signatures:
                saver.add(signature)
    except (OSError, ValueError, RuntimeError) as exc:
        path.unlink(missing_ok=True)
        raise ConfigurationError(f"sourmash could not write {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "write_sourmash_signature_archive",
        "output": str(path),
        "format": "sqlite" if suffix == ".sqldb" else "zip",
        "signatures": len(signatures),
        "provenance": {
            "package": _package(),
            "format": "sourmash zip collection (.zip, .sig.zip) or SQLite index (.sqldb)",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
        },
    }


def parse_sourmash_signature_json(
    json_text: str,
    *,
    ksize: int | None = None,
    require_single: bool = False,
    max_signatures: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse sourmash signature JSON text into sketches, without reading a file.

    ``require_single=False`` calls ``sourmash.signature.load_signatures_from_json`` and returns
    every signature, optionally only those with the given ``ksize``. ``require_single=True`` calls
    ``load_one_signature_from_json``, which fails unless exactly one signature matches.
    """

    if not isinstance(json_text, str) or not json_text.strip():
        raise ConfigurationError("json_text must be non-empty sourmash signature JSON text.")
    if len(json_text) > _MAX_JSON_CHARS:
        raise ConfigurationError(f"json_text is longer than {_MAX_JSON_CHARS} characters.")
    if ksize is not None:
        ksize = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    limit = _integer(max_signatures, "max_signatures", 1, _MAX_SIGNATURES)
    if not isinstance(require_single, bool):
        raise ConfigurationError("require_single must be a bool.")
    signature_module = require("sourmash.signature", _SOURMASH_EXTRA)
    try:
        if require_single:
            loaded = [signature_module.load_one_signature_from_json(json_text, ksize=ksize)]
        else:
            loaded = []
            for item in signature_module.load_signatures_from_json(json_text, ksize=ksize):
                if len(loaded) >= limit:
                    raise ConfigurationError(
                        f"JSON exceeds max_signatures={limit}; increase the explicit limit."
                    )
                loaded.append(item)
    except ConfigurationError:
        raise
    except Exception as exc:  # sourmash raises generic exceptions for malformed JSON
        raise ConfigurationError(f"sourmash could not parse the signature JSON: {exc}") from exc
    if not loaded:
        raise ConfigurationError("the JSON text holds no signature that matches the selection.")
    sketches = sorted(
        (_signature_record(item) for item in loaded),
        key=lambda record: (record["id"], record["ksize"], record["md5sum"]),
    )
    return {
        "status": "completed",
        "operation": "parse_sourmash_signature_json",
        "signatures": len(sketches),
        "require_single": require_single,
        "ksize": ksize,
        "sketches": sketches,
        "provenance": {
            "package": _package(),
            "source_function": (
                "sourmash.signature.load_one_signature_from_json"
                if require_single
                else "sourmash.signature.load_signatures_from_json"
            ),
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "output is sorted by id, ksize and md5sum; abundances are kept when "
            "the signature tracks them",
        },
    }


def select_sourmash_signatures(
    sources: Sequence[str | Path],
    *,
    ksize: int | None = None,
    name_pattern: str | None = None,
    md5_prefix: str | None = None,
    expect: str = "many",
    max_signatures: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load DNA signatures from sourmash files with k-mer size, name and md5 filters.

    ``expect="many"`` calls ``sourmash.sourmash_args.load_many_signatures`` on every source.
    ``expect="one"`` calls ``load_one_signature`` and needs exactly one matching signature in each
    source. ``expect="query"`` calls ``load_query_signature`` on a single source and needs exactly
    one signature after the ``ksize`` and ``md5_prefix`` selection. ``name_pattern`` is a regular
    expression (case-insensitive) searched in the name, file name and md5 of each
    signature; ``many`` and ``one`` use it.
    """

    if expect not in _LOAD_EXPECTATIONS:
        raise ConfigurationError(f"expect must be one of: {', '.join(_LOAD_EXPECTATIONS)}.")
    if isinstance(sources, (str, bytes, Path)) or not isinstance(sources, Sequence) or not sources:
        raise ConfigurationError("sources must be a non-empty list of signature file paths.")
    paths = [str(input_path(item)) for item in sources]
    if expect == "query" and len(paths) != 1:
        raise ConfigurationError("expect='query' takes exactly one source.")
    if ksize is not None:
        ksize = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    limit = _integer(max_signatures, "max_signatures", 1, _MAX_SIGNATURES)
    pattern: Any = None
    if name_pattern is not None:
        if not isinstance(name_pattern, str) or not name_pattern:
            raise ConfigurationError("name_pattern must be a non-empty regular expression.")
        if expect == "query":
            raise ConfigurationError("name_pattern is not used with expect='query'.")
        try:
            compiled = re.compile(name_pattern, re.IGNORECASE)
        except re.error as exc:
            raise ConfigurationError(
                f"name_pattern is not a valid regular expression: {exc}"
            ) from exc

        def pattern(values: Sequence[str]) -> bool:
            return any(compiled.search(value) for value in values)

    if md5_prefix is not None:
        if expect != "query":
            raise ConfigurationError("md5_prefix is only used with expect='query'.")
        if not isinstance(md5_prefix, str) or not re.fullmatch(r"[0-9a-fA-F]{1,32}", md5_prefix):
            raise ConfigurationError("md5_prefix must be 1 to 32 hexadecimal characters.")
    args = require("sourmash.sourmash_args", _SOURMASH_EXTRA)
    loaded: list[tuple[Any, str]] = []
    try:
        if expect == "many":
            progress = args.SignatureLoadingProgress()
            for item, location in args.load_many_signatures(
                paths, progress, ksize=ksize, moltype="DNA", pattern=pattern
            ):
                if len(loaded) >= limit:
                    raise ConfigurationError(
                        f"sources exceed max_signatures={limit}; increase the explicit limit."
                    )
                loaded.append((item, str(location)))
        elif expect == "one":
            for path in paths:
                loaded.append(
                    (
                        args.load_one_signature(
                            path, select_moltype="DNA", ksize=ksize, pattern=pattern
                        ),
                        path,
                    )
                )
        else:
            loaded.append((args.load_query_signature(paths[0], ksize, "DNA", md5_prefix), paths[0]))
    except ConfigurationError:
        raise
    except SystemExit as exc:
        raise ConfigurationError(
            "sourmash could not select signatures: no file, no match, or more than one match "
            f"where one was needed (exit code {exc.code})."
        ) from exc
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not select signatures: {exc}") from exc
    sketches = []
    for item, location in loaded:
        record = _signature_record(item)
        record["source"] = location
        sketches.append(record)
    sketches.sort(key=lambda record: (record["source"], record["id"], record["md5sum"]))
    return {
        "status": "completed",
        "operation": "select_sourmash_signatures",
        "expect": expect,
        "ksize": ksize,
        "name_pattern": name_pattern,
        "md5_prefix": md5_prefix,
        "signatures": len(sketches),
        "sketches": sketches,
        "provenance": {
            "package": _package(),
            "source_function": {
                "many": "sourmash.sourmash_args.load_many_signatures",
                "one": "sourmash.sourmash_args.load_one_signature",
                "query": "sourmash.sourmash_args.load_query_signature",
            }[expect],
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "only DNA scaled signatures; output is sorted by source, id and md5sum",
        },
    }


def compare_minhash_abundance_matrix(
    sketches: Sequence[Mapping[str, Any]],
    *,
    ignore_abundance: bool = False,
    return_ani: bool = False,
    n_jobs: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise matrix of scaled sketches that can use abundance and run in several processes.

    Calls ``sourmash.compare.compare_all_pairs``, which uses ``compare_serial`` for ``n_jobs=1``
    and ``compare_parallel`` otherwise. Sketches with an ``abundances`` list give the angular
    (cosine) similarity unless ``ignore_abundance`` is true; sketches without abundances give the
    Jaccard index. ``return_ani=True`` returns the estimated average nucleotide identity instead
    of the similarity.
    """

    items = _sketch_list(sketches, "sketches", 2, _MATRIX_LIMIT)
    for name, flag in (("ignore_abundance", ignore_abundance), ("return_ani", return_ani)):
        if not isinstance(flag, bool):
            raise ConfigurationError(f"{name} must be a bool.")
    jobs = _integer(n_jobs, "n_jobs", 1, 16)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    compare_module = require("sourmash.compare", _SOURMASH_EXTRA)
    built = [
        _abundance_sketch(sourmash, item, f"sketches[{index}]", default_name=f"sketch_{index}")
        for index, item in enumerate(items)
    ]
    signatures = [signature for signature, _track in built]
    tracked = {track for _signature, track in built}
    if len(tracked) != 1:
        raise ConfigurationError("give abundances for all sketches or for none.")
    first = signatures[0].minhash
    if any(
        (item.minhash.ksize, item.minhash.scaled) != (first.ksize, first.scaled)
        for item in signatures
    ):
        raise ConfigurationError("sketches must share ksize and scaled to be compared.")
    try:
        matrix = compare_module.compare_all_pairs(
            signatures,
            ignore_abundance,
            downsample=False,
            n_jobs=None if jobs == 1 else jobs,
            return_ani=return_ani,
        )
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash comparison failed: {exc}") from exc
    weighted = tracked == {True} and not ignore_abundance
    return {
        "status": "completed",
        "operation": "compare_minhash_abundance_matrix",
        "measure": "ani" if return_ani else ("angular_similarity" if weighted else "jaccard"),
        "abundance_weighted": weighted,
        "n_jobs": jobs,
        "ids": [signature.name for signature in signatures],
        "ksize": int(first.ksize),
        "scaled": int(first.scaled),
        "matrix": [[float(value) for value in row] for row in matrix.tolist()],
        "provenance": {
            "package": _package(),
            "source_function": "sourmash.compare.compare_all_pairs"
            + (" (compare_parallel)" if jobs > 1 else " (compare_serial)"),
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "matrix[i][j] compares sketch i with sketch j; the matrix is symmetric",
        },
    }


def search_minhash_database_ranked(
    query: Mapping[str, Any],
    database: Sequence[Mapping[str, Any]],
    *,
    mode: str = "jaccard",
    threshold: float = 0.08,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rank database sketches by similarity to a query sketch, as the sourmash search command does.

    ``mode="jaccard"``, ``"containment"`` and ``"max_containment"`` call
    ``sourmash.search.search_databases_with_flat_query`` with the matching score.
    ``mode="abundance"`` calls ``search_databases_with_abund_query`` and needs ``abundances`` on the
    query and all database sketches; its score is the angular similarity. Matches below
    ``threshold`` are dropped and the rest are sorted by score, then id.
    """

    if mode not in _SEARCH_MODES:
        raise ConfigurationError(f"mode must be one of: {', '.join(_SEARCH_MODES)}.")
    cutoff = _real(threshold, "threshold")
    if not 0.0 <= cutoff <= 1.0:
        raise ConfigurationError("threshold must be between 0 and 1.")
    items = _sketch_list(database, "database", 1, _MAX_SIGNATURES)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    search = require("sourmash.search", _SOURMASH_EXTRA)
    index_module = require("sourmash.index", _SOURMASH_EXTRA)
    query_signature, query_tracked = _abundance_sketch(
        sourmash, query, "query", default_name="query"
    )
    built = [
        _abundance_sketch(sourmash, item, f"database[{index}]", default_name=f"database_{index}")
        for index, item in enumerate(items)
    ]
    if mode == "abundance" and not (query_tracked and all(track for _s, track in built)):
        raise ConfigurationError("mode='abundance' needs abundances on the query and every sketch.")
    if mode != "abundance" and (query_tracked or any(track for _s, track in built)):
        raise ConfigurationError(
            "abundances are only used with mode='abundance'; remove them for the other modes."
        )
    signatures = [signature for signature, _track in built]
    query_minhash = query_signature.minhash
    if any(item.minhash.ksize != query_minhash.ksize for item in signatures):
        raise ConfigurationError("database sketches must share the query ksize.")
    try:
        index = index_module.LinearIndex(signatures)
        if mode == "abundance":
            results = search.search_databases_with_abund_query(
                query_signature, [index], threshold=cutoff
            )
        else:
            results = search.search_databases_with_flat_query(
                query_signature,
                [index],
                threshold=cutoff,
                do_containment=mode == "containment",
                do_max_containment=mode == "max_containment",
            )
    except (AssertionError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash search failed: {exc}") from exc
    matches = [
        {
            "id": str(result.match_name),
            "md5": str(result.match_md5),
            "score": float(result.similarity),
        }
        for result in results
    ]
    matches.sort(key=lambda item: (-item["score"], item["id"], item["md5"]))
    return {
        "status": "completed",
        "operation": "search_minhash_database_ranked",
        "query_id": str(query_signature.name),
        "ksize": int(query_minhash.ksize),
        "mode": mode,
        "threshold": cutoff,
        "database_size": len(signatures),
        "match_count": len(matches),
        "matches": matches,
        "provenance": {
            "package": _package(),
            "source_function": "sourmash.search.search_databases_with_abund_query"
            if mode == "abundance"
            else "sourmash.search.search_databases_with_flat_query",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": {
                "jaccard": "score is the Jaccard index",
                "containment": "score is the fraction of the query found in the match",
                "max_containment": "score is the larger of the two containments",
                "abundance": "score is the angular similarity of abundance-weighted sketches",
            }[mode],
        },
    }


def search_sourmash_sbt_index(
    index: str | Path,
    query: Mapping[str, Any],
    *,
    threshold: float = 0.08,
    verbose: bool = True,
) -> dict[str, Any]:
    """Search a saved sourmash sequence bloom tree (SBT) index with a query sketch.

    Loads the index with ``sourmash.sbtmh.load_sbt_index`` and calls
    ``sourmash.sbtmh.search_sbt_index``. The score is the Jaccard index, so the matches equal
    those of ``search_minhash_database_ranked(mode="jaccard")`` over the same sketches. The index
    is a ``.sbt.zip`` file or ``.sbt.json`` file written by ``sourmash index``; the query must
    share the ``ksize`` and ``scaled`` of the indexed sketches.
    """

    path = input_path(index)
    cutoff = _real(threshold, "threshold")
    if not 0.0 <= cutoff <= 1.0:
        raise ConfigurationError("threshold must be between 0 and 1.")
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    sbtmh = require("sourmash.sbtmh", _SOURMASH_EXTRA)
    query_signature, query_tracked = _abundance_sketch(
        sourmash, query, "query", default_name="query"
    )
    if query_tracked:
        raise ConfigurationError("the SBT search ignores abundance; remove abundances.")
    try:
        tree = sbtmh.load_sbt_index(str(path))
        results = list(sbtmh.search_sbt_index(tree, query_signature, cutoff))
    except ConfigurationError:
        raise
    except Exception as exc:  # sourmash raises generic exceptions for bad index files
        raise ConfigurationError(f"sourmash could not search {path.name}: {exc}") from exc
    matches = [
        {"id": str(match.name), "md5": str(match.md5sum()), "jaccard": float(score)}
        for match, score in results
    ]
    matches.sort(key=lambda item: (-item["jaccard"], item["id"], item["md5"]))
    return {
        "status": "completed",
        "operation": "search_sourmash_sbt_index",
        "index": str(path),
        "query_id": str(query_signature.name),
        "ksize": int(query_signature.minhash.ksize),
        "threshold": cutoff,
        "match_count": len(matches),
        "matches": matches,
        "provenance": {
            "package": _package(),
            "source_function": "sourmash.sbtmh.search_sbt_index",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "matches are sorted by Jaccard index, then id and md5",
        },
    }


def downsample_minhash_num_sketch(
    sketch: Mapping[str, Any],
    *,
    num: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Downsample a num-based MinHash sketch to fewer hashes with sourmash.

    The sketch has ``ksize``, ``num`` and ``hashes``: the ``num`` smallest hash values of the
    sequence. Calls ``sourmash.minhash.flatten_and_downsample_num``, which keeps the smallest
    ``num`` hashes. The new ``num`` cannot exceed the sketch's own ``num``.
    """

    if not isinstance(sketch, Mapping):
        raise ConfigurationError("sketch must be a mapping with ksize, num and hashes.")
    ksize = _integer(sketch.get("ksize"), "sketch ksize", 1, _KSIZE_LIMIT)
    source_num = _integer(sketch.get("num"), "sketch num", 1, 10_000_000)
    hashes = sketch.get("hashes")
    if isinstance(hashes, (str, bytes)) or not isinstance(hashes, Sequence) or not hashes:
        raise ConfigurationError("sketch hashes must be a non-empty list of integers.")
    values: list[int] = []
    for value in hashes:
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _HASH_LIMIT:
            raise ConfigurationError("sketch hashes must be unsigned 64-bit integers.")
        values.append(value)
    if len(set(values)) > source_num:
        raise ConfigurationError("sketch has more distinct hashes than its num.")
    target = _integer(num, "num", 1, source_num)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    minhash_module = require("sourmash.minhash", _SOURMASH_EXTRA)
    try:
        original = sourmash.MinHash(n=source_num, ksize=ksize)
        original.add_many(sorted(set(values)))
        result = minhash_module.flatten_and_downsample_num(original, target)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not downsample the sketch: {exc}") from exc
    kept = sorted(int(value) for value in result.hashes)
    identifier = sketch.get("id")
    return {
        "status": "completed",
        "operation": "downsample_minhash_num_sketch",
        "id": None if identifier is None else str(identifier),
        "ksize": ksize,
        "input_num": source_num,
        "num": int(result.num),
        "input_hash_count": len(set(values)),
        "hash_count": len(kept),
        "hashes": kept,
        "provenance": {
            "package": _package(),
            "algorithm": "flatten_and_downsample_num: keep the smallest num hashes",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "abundance is not kept",
        },
    }


__all__ = [
    "compare_minhash_abundance_matrix",
    "downsample_minhash_num_sketch",
    "parse_sourmash_signature_json",
    "search_minhash_database_ranked",
    "search_sourmash_sbt_index",
    "select_sourmash_signatures",
    "write_sourmash_signature_archive",
]

_register_public_verbose(globals(), tuple(__all__))
