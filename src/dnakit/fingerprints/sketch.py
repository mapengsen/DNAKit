"""Deterministic bottom-k MinHash and FracMinHash sketches."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable

from dnakit.core import Gap, Topology
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError

from ._shared import SequenceInput, sequence_and_id, validate_bool, validate_positive_int
from .results import SketchResult

DEFAULT_MAX_SKETCH_HASHES = 1_000_000
DEFAULT_MAX_UNIQUE_HASHES = 1_000_000
_UINT64_SPACE = 1 << 64
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _hash64(word: str, seed: int) -> int:
    payload = seed.to_bytes(8, "little", signed=False) + word.encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def _canonical(word: str) -> str:
    reverse = word.translate(_COMPLEMENT)[::-1]
    return min(word, reverse)


def _validate_common(
    *,
    k: int,
    canonical: bool,
    seed: int,
    max_hashes: int,
    max_unique_hashes: int,
) -> None:
    validate_positive_int(k, "k")
    validate_bool(canonical, "canonical")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < _UINT64_SPACE:
        raise ConfigurationError(
            "seed must be an integer in [0, 2**64).",
            code="INVALID_SKETCH_SEED",
        )
    validate_positive_int(max_hashes, "max_hashes")
    if max_hashes > DEFAULT_MAX_SKETCH_HASHES:
        raise ConfigurationError(
            "max_hashes exceeds the hard safety ceiling.",
            code="SKETCH_HASH_LIMIT_INVALID",
            context={"max_hashes": max_hashes, "hard_limit": DEFAULT_MAX_SKETCH_HASHES},
        )
    validate_positive_int(max_unique_hashes, "max_unique_hashes")
    if max_unique_hashes > DEFAULT_MAX_UNIQUE_HASHES:
        raise ConfigurationError(
            "max_unique_hashes exceeds the hard safety ceiling.",
            code="SKETCH_UNIQUE_HASH_LIMIT_INVALID",
            context={
                "max_unique_hashes": max_unique_hashes,
                "hard_limit": DEFAULT_MAX_UNIQUE_HASHES,
            },
        )


def _unique_hashes(
    value: SequenceInput,
    *,
    k: int,
    canonical: bool,
    seed: int,
    max_unique_hashes: int,
) -> tuple[str | None, int, int, tuple[int, ...]]:
    sequence, sequence_id = sequence_and_id(value)
    if sequence.topology is Topology.CIRCULAR:
        raise ConfigurationError(
            "Native sketches currently require a linear sequence.",
            code="SKETCH_CIRCULAR_UNSUPPORTED",
            hint="Rotate a circular sequence to an explicit origin before sketching.",
        )
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    if gaps:
        raise UnsupportedGapOperationError(
            "Sketches do not silently omit or cross explicit Gap objects.",
            code="SKETCH_GAP_NOT_ALLOWED",
            context={"gap_count": len(gaps)},
        )
    if any(symbol not in "ACGT" for symbol in sequence.symbols):
        raise ConfigurationError(
            "Native sketches require canonical A/C/G/T symbols.",
            code="SKETCH_AMBIGUITY_NOT_ALLOWED",
        )
    words: Iterable[str] = (
        sequence.symbols[index : index + k]
        for index in range(max(0, sequence.symbol_length - k + 1))
    )
    if canonical:
        words = (_canonical(word) for word in words)
    unique_hashes: set[int] = set()
    for word in words:
        unique_hashes.add(_hash64(word, seed))
        if len(unique_hashes) > max_unique_hashes:
            raise ConfigurationError(
                "Sketch calculation exceeds max_unique_hashes.",
                code="SKETCH_UNIQUE_HASH_LIMIT_EXCEEDED",
                context={"unique_hash_count_lower_bound": max_unique_hashes + 1},
            )
    hashes = tuple(sorted(unique_hashes))
    observation_count = max(0, sequence.symbol_length - k + 1)
    return sequence_id, sequence.symbol_length, observation_count, hashes


def minhash(
    value: SequenceInput,
    *,
    k: int,
    num_hashes: int = 1_000,
    canonical: bool = True,
    seed: int = 0,
    max_hashes: int = DEFAULT_MAX_SKETCH_HASHES,
    max_unique_hashes: int = DEFAULT_MAX_UNIQUE_HASHES,
) -> SketchResult:
    """Return the smallest distinct 64-bit SHA-256-derived k-mer hashes."""

    _validate_common(
        k=k,
        canonical=canonical,
        seed=seed,
        max_hashes=max_hashes,
        max_unique_hashes=max_unique_hashes,
    )
    validate_positive_int(num_hashes, "num_hashes")
    if num_hashes > max_hashes:
        raise ConfigurationError(
            "num_hashes cannot exceed max_hashes.",
            code="SKETCH_HASH_LIMIT_EXCEEDED",
        )
    sequence_id, symbol_length, observation_count, hashes = _unique_hashes(
        value,
        k=k,
        canonical=canonical,
        seed=seed,
        max_unique_hashes=max_unique_hashes,
    )
    return SketchResult(
        name="minhash",
        method="bottom-k-sha256-64",
        schema_version="dnakit.sketch.minhash.v1",
        sequence_id=sequence_id,
        symbol_length=symbol_length,
        k=k,
        canonical=canonical,
        seed=seed,
        selection="bottom_k",
        hashes=hashes[:num_hashes],
        num_hashes=num_hashes,
        scaled=None,
        threshold=None,
        observation_count=observation_count,
        unique_hash_count=len(hashes),
        max_hashes=max_hashes,
        parameters=FrozenDict(
            {
                "hash": "sha256-first-64-big-endian",
                "deduplicate_kmers": True,
                "max_unique_hashes": max_unique_hashes,
                "topology": "linear",
            }
        ),
    )


def fracminhash(
    value: SequenceInput,
    *,
    k: int,
    scaled: int = 1_000,
    canonical: bool = True,
    seed: int = 0,
    max_hashes: int = DEFAULT_MAX_SKETCH_HASHES,
    max_unique_hashes: int = DEFAULT_MAX_UNIQUE_HASHES,
) -> SketchResult:
    """Retain hashes below ``floor(2**64 / scaled)`` with a hard output limit."""

    _validate_common(
        k=k,
        canonical=canonical,
        seed=seed,
        max_hashes=max_hashes,
        max_unique_hashes=max_unique_hashes,
    )
    validate_positive_int(scaled, "scaled")
    if scaled > _UINT64_SPACE:
        raise ConfigurationError(
            "scaled must not exceed 2**64.",
            code="INVALID_FRACMINHASH_SCALED",
        )
    sequence_id, symbol_length, observation_count, hashes = _unique_hashes(
        value,
        k=k,
        canonical=canonical,
        seed=seed,
        max_unique_hashes=max_unique_hashes,
    )
    threshold = _UINT64_SPACE // scaled
    retained = tuple(item for item in hashes if item < threshold)
    if len(retained) > max_hashes:
        raise ConfigurationError(
            "FracMinHash output exceeds max_hashes.",
            code="SKETCH_HASH_LIMIT_EXCEEDED",
            context={"retained_hashes": len(retained), "max_hashes": max_hashes},
            hint="Increase scaled or explicitly raise max_hashes within the hard ceiling.",
        )
    return SketchResult(
        name="fracminhash",
        method="threshold-sha256-64",
        schema_version="dnakit.sketch.fracminhash.v1",
        sequence_id=sequence_id,
        symbol_length=symbol_length,
        k=k,
        canonical=canonical,
        seed=seed,
        selection="scaled",
        hashes=retained,
        num_hashes=None,
        scaled=scaled,
        threshold=threshold,
        observation_count=observation_count,
        unique_hash_count=len(hashes),
        max_hashes=max_hashes,
        parameters=FrozenDict(
            {
                "hash": "sha256-first-64-big-endian",
                "deduplicate_kmers": True,
                "max_unique_hashes": max_unique_hashes,
                "topology": "linear",
            }
        ),
    )


__all__ = [
    "DEFAULT_MAX_SKETCH_HASHES",
    "DEFAULT_MAX_UNIQUE_HASHES",
    "fracminhash",
    "minhash",
]
