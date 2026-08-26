"""Deterministic in-memory sketch index with integrity-checked JSON persistence."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

from dnakit.core import DNA, DNARecord, DNASet
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.fingerprints import SketchResult, minhash

from .results import NearestNeighborHit, NearestNeighborResult
from .sketch import sketch_similarity

_INDEX_SCHEMA = "dnakit.sketch-index.v1"
DEFAULT_MAX_INDEX_HASHES = 10_000_000
DEFAULT_MAX_INDEX_FILE_BYTES = 512 * 1024 * 1024


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise InputFormatError(
                "Sketch index contains a duplicate JSON key.",
                code="INVALID_SKETCH_INDEX",
                context={"key": key},
            )
        result[key] = value
    return result


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


@dataclass(frozen=True)
class SketchIndex:
    """Versioned exact scan index over compatible bottom-k sketches."""

    ids: tuple[str, ...]
    sketches: tuple[SketchResult, ...]
    k: int
    num_hashes: int
    canonical: bool
    seed: int
    schema_version: str = _INDEX_SCHEMA

    def __post_init__(self) -> None:
        if self.schema_version != _INDEX_SCHEMA:
            raise ConfigurationError("Unsupported sketch index schema.")
        if not isinstance(self.ids, tuple) or not isinstance(self.sketches, tuple):
            raise ConfigurationError("Sketch index IDs and sketches must be immutable tuples.")
        if len(self.ids) != len(self.sketches) or len(set(self.ids)) != len(self.ids):
            raise ConfigurationError("Sketch index IDs must be unique and align with sketches.")
        if any(not isinstance(item, str) or not item.strip() for item in self.ids):
            raise ConfigurationError("Sketch index IDs must be non-empty strings.")
        if any(not isinstance(item, SketchResult) for item in self.sketches):
            raise ConfigurationError("Sketch index entries must be SketchResult objects.")
        for name in ("k", "num_hashes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ConfigurationError(f"Sketch index {name} must be a positive integer.")
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("Sketch index canonical must be boolean.")
        if (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, int)
            or not 0 <= self.seed < 2**64
        ):
            raise ConfigurationError("Sketch index seed must be an unsigned 64-bit integer.")
        total_hashes = sum(len(item.hashes) for item in self.sketches)
        if total_hashes > DEFAULT_MAX_INDEX_HASHES:
            raise ConfigurationError(
                "Sketch index exceeds the hard total-hash limit.",
                code="INDEX_HASH_LIMIT_EXCEEDED",
                context={
                    "total_hashes": total_hashes,
                    "hard_limit": DEFAULT_MAX_INDEX_HASHES,
                },
            )
        for item in self.sketches:
            if (
                item.selection != "bottom_k"
                or item.k != self.k
                or item.num_hashes != self.num_hashes
                or item.canonical != self.canonical
                or item.seed != self.seed
            ):
                raise ConfigurationError("Sketch index contains an incompatible entry.")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def build_sketch_index(
    records: DNA | DNASet | Iterable[DNARecord],
    *,
    k: int = 15,
    num_hashes: int = 1_000,
    canonical: bool = True,
    seed: int = 0,
    max_records: int = 100_000,
) -> SketchIndex:
    """Build a bounded, deterministic bottom-k index in input order."""

    if isinstance(max_records, bool) or not isinstance(max_records, int) or max_records <= 0:
        raise ConfigurationError("max_records must be a positive integer.")
    if isinstance(k, bool) or not isinstance(k, int) or k <= 0:
        raise ConfigurationError("k must be a positive integer.")
    if isinstance(num_hashes, bool) or not isinstance(num_hashes, int) or num_hashes <= 0:
        raise ConfigurationError("num_hashes must be a positive integer.")
    if not isinstance(canonical, bool):
        raise ConfigurationError("canonical must be boolean.")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**64:
        raise ConfigurationError("seed must be an unsigned 64-bit integer.")
    try:
        source = records.records if isinstance(records, (DNA, DNASet)) else iter(records)
    except TypeError as exc:
        raise ConfigurationError("records must be DNASet or an iterable of DNARecord.") from exc
    materialized: list[DNARecord] = []
    for index, record in enumerate(source):
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "Index inputs must be DNARecord objects.",
                code="INVALID_INDEX_RECORD",
                context={"record_index": index},
            )
        if index >= max_records:
            raise ConfigurationError(
                "Index input exceeds max_records.",
                code="INDEX_RECORD_LIMIT_EXCEEDED",
                context={"record_count_lower_bound": max_records + 1},
            )
        materialized.append(record)
    ids = tuple(record.id for record in materialized)
    if len(set(ids)) != len(ids):
        raise ConfigurationError("Index record IDs must be unique.", code="DUPLICATE_INDEX_ID")
    sketches_list: list[SketchResult] = []
    total_hashes = 0
    for record in materialized:
        sketch = minhash(record, k=k, num_hashes=num_hashes, canonical=canonical, seed=seed)
        total_hashes += len(sketch.hashes)
        if total_hashes > DEFAULT_MAX_INDEX_HASHES:
            raise ConfigurationError(
                "Sketch index exceeds the hard total-hash limit.",
                code="INDEX_HASH_LIMIT_EXCEEDED",
                context={"total_hashes_lower_bound": total_hashes},
            )
        sketches_list.append(sketch)
    sketches = tuple(sketches_list)
    return SketchIndex(ids, sketches, k, num_hashes, canonical, seed)


def nearest_neighbors(
    query: DNA | DNARecord | SketchResult,
    index: SketchIndex,
    *,
    top_k: int = 10,
    min_similarity: float = 0.0,
) -> NearestNeighborResult:
    """Return stable top-k Jaccard hits from an in-memory sketch index."""

    if not isinstance(index, SketchIndex):
        raise ConfigurationError("index must be SketchIndex.", code="INVALID_SKETCH_INDEX")
    if isinstance(top_k, bool) or not isinstance(top_k, int) or top_k <= 0:
        raise ConfigurationError("top_k must be a positive integer.")
    if (
        isinstance(min_similarity, bool)
        or not isinstance(min_similarity, (int, float))
        or not 0 <= min_similarity <= 1
    ):
        raise ConfigurationError("min_similarity must be in [0, 1].")
    sketch = (
        query
        if isinstance(query, SketchResult)
        else minhash(
            query,
            k=index.k,
            num_hashes=index.num_hashes,
            canonical=index.canonical,
            seed=index.seed,
        )
    )
    scored: list[NearestNeighborHit] = []
    for position, (record_id, target) in enumerate(zip(index.ids, index.sketches, strict=True)):
        comparison = sketch_similarity(sketch, target)
        if comparison.value >= min_similarity:
            scored.append(
                NearestNeighborHit(
                    record_id=record_id,
                    index=position,
                    similarity=comparison.value,
                    shared_hash_count=comparison.shared_hash_count,
                )
            )
    scored.sort(key=lambda item: (-item.similarity, -item.shared_hash_count, item.index))
    return NearestNeighborResult(
        name="nearest_neighbors",
        method="exact-scan-minhash-jaccard",
        algorithm_version="dnakit-sketch-nearest-neighbor-v1",
        query_id=sketch.sequence_id,
        index_size=len(index.ids),
        hits=tuple(scored[:top_k]),
        top_k=top_k,
        min_similarity=float(min_similarity),
        parameters=FrozenDict(
            {
                "index_schema": index.schema_version,
                "k": index.k,
                "num_hashes": index.num_hashes,
                "canonical": index.canonical,
                "seed": index.seed,
                "tie_break": "similarity-desc,shared-desc,input-index-asc",
            }
        ),
    )


def save_sketch_index(
    index: SketchIndex,
    path: str | os.PathLike[str],
    *,
    overwrite: bool = False,
) -> str:
    """Atomically persist a self-checksummed index and return its SHA-256."""

    if not isinstance(index, SketchIndex):
        raise ConfigurationError("index must be SketchIndex.", code="INVALID_SKETCH_INDEX")
    if not isinstance(overwrite, bool):
        raise ConfigurationError("overwrite must be boolean.")
    target = Path(path).expanduser().absolute()
    if target.exists() and not overwrite:
        raise ConfigurationError("Index target exists.", code="INDEX_EXISTS")
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = index.to_dict()
    digest = hashlib.sha256(_canonical_json(payload)).hexdigest()
    envelope = {"schema_version": _INDEX_SCHEMA, "sha256": digest, "payload": payload}
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(_canonical_json(envelope) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary_name, target)
        else:
            try:
                os.link(temporary_name, target)
            except FileExistsError as exc:
                raise ConfigurationError("Index target exists.", code="INDEX_EXISTS") from exc
            Path(temporary_name).unlink()
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return digest


def _sketch_from_payload(payload: object) -> SketchResult:
    if not isinstance(payload, dict):
        raise InputFormatError("Invalid sketch index entry.", code="INVALID_SKETCH_INDEX")
    try:
        return SketchResult(
            name=payload["name"],
            method=payload["method"],
            schema_version=payload["schema_version"],
            sequence_id=payload.get("sequence_id"),
            symbol_length=payload["symbol_length"],
            k=payload["k"],
            canonical=payload["canonical"],
            seed=payload["seed"],
            selection=payload["selection"],
            hashes=tuple(payload["hashes"]),
            num_hashes=payload.get("num_hashes"),
            scaled=payload.get("scaled"),
            threshold=payload.get("threshold"),
            observation_count=payload["observation_count"],
            unique_hash_count=payload["unique_hash_count"],
            max_hashes=payload["max_hashes"],
            parameters=FrozenDict(payload["parameters"]),
        )
    except (KeyError, TypeError, ValueError, ConfigurationError) as exc:
        raise InputFormatError("Invalid sketch index entry.", code="INVALID_SKETCH_INDEX") from exc


def load_sketch_index(path: str | os.PathLike[str]) -> SketchIndex:
    """Load and verify a DNAKit sketch index."""

    target = Path(path).expanduser().absolute()
    if not target.is_file() or target.is_symlink():
        raise InputFormatError("Index path must be a regular file.", code="INVALID_SKETCH_INDEX")
    if target.stat().st_size > DEFAULT_MAX_INDEX_FILE_BYTES:
        raise InputFormatError(
            "Sketch index file exceeds the hard byte limit.",
            code="SKETCH_INDEX_FILE_LIMIT",
            context={
                "byte_size": target.stat().st_size,
                "hard_limit": DEFAULT_MAX_INDEX_FILE_BYTES,
            },
        )
    try:
        with target.open("rb") as handle:
            raw = handle.read(DEFAULT_MAX_INDEX_FILE_BYTES + 1)
        if len(raw) > DEFAULT_MAX_INDEX_FILE_BYTES:
            raise InputFormatError(
                "Sketch index grew beyond the hard byte limit while being read.",
                code="SKETCH_INDEX_FILE_LIMIT",
                context={"hard_limit": DEFAULT_MAX_INDEX_FILE_BYTES},
            )
        envelope = json.loads(raw.decode("utf-8"), object_pairs_hook=_unique_json_object)
        if not isinstance(envelope, dict) or envelope.get("schema_version") != _INDEX_SCHEMA:
            raise ValueError("unsupported envelope schema")
        payload = envelope["payload"]
        expected = envelope["sha256"]
    except InputFormatError:
        raise
    except (
        OSError,
        UnicodeError,
        json.JSONDecodeError,
        KeyError,
        TypeError,
        ValueError,
        RecursionError,
    ) as exc:
        raise InputFormatError("Cannot decode sketch index.", code="INVALID_SKETCH_INDEX") from exc
    actual = hashlib.sha256(_canonical_json(payload)).hexdigest()
    if not isinstance(expected, str) or actual != expected:
        raise InputFormatError("Sketch index checksum mismatch.", code="SKETCH_INDEX_CHECKSUM")
    if not isinstance(payload, dict) or payload.get("schema_version") != _INDEX_SCHEMA:
        raise InputFormatError("Unsupported sketch index schema.", code="INVALID_SKETCH_INDEX")
    try:
        return SketchIndex(
            ids=tuple(payload["ids"]),
            sketches=tuple(_sketch_from_payload(item) for item in payload["sketches"]),
            k=payload["k"],
            num_hashes=payload["num_hashes"],
            canonical=payload["canonical"],
            seed=payload["seed"],
            schema_version=payload["schema_version"],
        )
    except (KeyError, TypeError, ValueError, ConfigurationError) as exc:
        raise InputFormatError(
            "Invalid sketch index payload.", code="INVALID_SKETCH_INDEX"
        ) from exc


__all__ = [
    "DEFAULT_MAX_INDEX_FILE_BYTES",
    "DEFAULT_MAX_INDEX_HASHES",
    "SketchIndex",
    "build_sketch_index",
    "load_sketch_index",
    "nearest_neighbors",
    "save_sketch_index",
]
