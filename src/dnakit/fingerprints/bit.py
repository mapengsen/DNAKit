"""Fixed-length hashed k-mer and named-panel bit fingerprints."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from typing import Literal

from dnakit.core import Gap
from dnakit.core._json import FrozenDict
from dnakit.descriptors import kmer_statistics
from dnakit.exceptions import ConfigurationError
from dnakit.patterns import scan_motif

from ._shared import (
    FingerprintAmbiguityPolicy,
    FingerprintRepresentation,
    SequenceInput,
    coerce_enum,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from .results import BitFingerprintResult, BitFingerprintValues

PanelMode = Literal["exact", "iupac"]
DEFAULT_HASHED_KMER_BITS = 2_048
MAX_BIT_FINGERPRINT_DIMENSION = 1_000_000
DEFAULT_MAX_PANEL_SIZE = 10_000
_UINT64_SPACE = 1 << 64


def _hash_bit(word: str, *, seed: int, n_bits: int) -> int:
    payload = seed.to_bytes(8, "little", signed=False) + word.encode("ascii")
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return value % n_bits


def _materialize_bits(
    feature_names: tuple[str, ...],
    set_indices: set[int],
    representation: FingerprintRepresentation,
) -> BitFingerprintValues:
    if representation is FingerprintRepresentation.DENSE:
        return tuple(int(index in set_indices) for index in range(len(feature_names)))
    return FrozenDict({feature_names[index]: 1 for index in sorted(set_indices)})


def hashed_kmer_fingerprint(
    value: SequenceInput,
    *,
    k: int,
    n_bits: int = DEFAULT_HASHED_KMER_BITS,
    canonical: bool = True,
    seed: int = 0,
    representation: FingerprintRepresentation | str = FingerprintRepresentation.DENSE,
    ambiguity_policy: FingerprintAmbiguityPolicy | str = FingerprintAmbiguityPolicy.ERROR,
    overlapping: bool = True,
    cross_gaps: bool = False,
) -> BitFingerprintResult:
    """Hash observed k-mers into a deterministic fixed-length binary vector."""

    validate_positive_int(k, "k")
    validate_positive_int(n_bits, "n_bits")
    if n_bits > MAX_BIT_FINGERPRINT_DIMENSION:
        raise ConfigurationError(
            "n_bits exceeds the bit fingerprint dimension limit.",
            code="BIT_FINGERPRINT_DIMENSION_LIMIT",
            context={
                "n_bits": n_bits,
                "max_dimension": MAX_BIT_FINGERPRINT_DIMENSION,
            },
        )
    validate_bool(canonical, "canonical")
    validate_bool(overlapping, "overlapping")
    validate_bool(cross_gaps, "cross_gaps")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < _UINT64_SPACE:
        raise ConfigurationError(
            "seed must be an integer in [0, 2**64).",
            code="INVALID_BIT_FINGERPRINT_SEED",
        )
    resolved_representation = coerce_enum(
        representation,
        FingerprintRepresentation,
        "fingerprint representation",
    )
    resolved_ambiguity = coerce_enum(
        ambiguity_policy,
        FingerprintAmbiguityPolicy,
        "fingerprint ambiguity policy",
    )
    sequence, sequence_id = sequence_and_id(value)
    statistics = kmer_statistics(
        value,
        k,
        overlapping=overlapping,
        canonical=canonical,
        ambiguity_policy=resolved_ambiguity.value,
        cross_gaps=cross_gaps,
    )
    set_indices = {_hash_bit(word, seed=seed, n_bits=n_bits) for word in statistics.counts}
    feature_names = tuple(f"bit:{index}" for index in range(n_bits))
    return BitFingerprintResult(
        name="hashed_kmer_fingerprint",
        method="sha256-kmer-modulo-bitset",
        schema_version="dnakit.hashed-kmer-bit.v1",
        sequence_id=sequence_id,
        symbol_length=sequence.symbol_length,
        gap_count=statistics.gap_count,
        unknown_gap_count=statistics.unknown_gap_count,
        representation=resolved_representation,
        feature_names=feature_names,
        values=_materialize_bits(feature_names, set_indices, resolved_representation),
        parameters=FrozenDict(
            {
                "k": k,
                "n_bits": n_bits,
                "canonical": canonical,
                "seed": seed,
                "hash": "sha256-first-64-big-endian-modulo-n-bits",
                "ambiguity_policy": resolved_ambiguity.value,
                "overlapping": overlapping,
                "cross_gaps": cross_gaps,
            }
        ),
        observation_count=statistics.denominator,
        ignored_ambiguity_count=statistics.ignored_ambiguity_count,
        max_dimension=MAX_BIT_FINGERPRINT_DIMENSION,
    )


def panel_fingerprint(
    value: SequenceInput,
    panel: Mapping[str, str],
    *,
    mode: PanelMode = "iupac",
    overlapping: bool = True,
    representation: FingerprintRepresentation | str = FingerprintRepresentation.DENSE,
    max_panel_size: int = DEFAULT_MAX_PANEL_SIZE,
    max_matches_per_pattern: int = 100_000,
) -> BitFingerprintResult:
    """Encode the presence of each named exact/IUPAC pattern in a panel."""

    if not isinstance(panel, Mapping) or not panel:
        raise ConfigurationError("panel must be a non-empty name-to-pattern mapping.")
    validate_positive_int(max_panel_size, "max_panel_size")
    if max_panel_size > DEFAULT_MAX_PANEL_SIZE:
        raise ConfigurationError(
            "max_panel_size exceeds the hard safety ceiling.",
            code="PANEL_SIZE_LIMIT_INVALID",
            context={
                "max_panel_size": max_panel_size,
                "hard_limit": DEFAULT_MAX_PANEL_SIZE,
            },
        )
    if len(panel) > max_panel_size:
        raise ConfigurationError(
            "panel exceeds max_panel_size.",
            code="PANEL_SIZE_LIMIT_EXCEEDED",
            context={"panel_size": len(panel), "max_panel_size": max_panel_size},
        )
    if any(not isinstance(name, str) or not name.strip() for name in panel):
        raise ConfigurationError("panel names must be non-empty strings.")
    if any(not isinstance(pattern, str) or not pattern.strip() for pattern in panel.values()):
        raise ConfigurationError("panel patterns must be non-empty strings.")
    if mode not in ("exact", "iupac"):
        raise ConfigurationError(
            "panel mode must be 'exact' or 'iupac'.",
            context={"mode": mode},
        )
    validate_bool(overlapping, "overlapping")
    validate_positive_int(max_matches_per_pattern, "max_matches_per_pattern")
    resolved_representation = coerce_enum(
        representation,
        FingerprintRepresentation,
        "fingerprint representation",
    )
    sequence, sequence_id = sequence_and_id(value)
    names = tuple(sorted(panel))
    results = tuple(
        scan_motif(
            value,
            panel[name],
            mode=mode,
            name=name,
            strand="both",
            overlapping=overlapping,
            merge_strands=True,
            max_matches=max_matches_per_pattern,
        )
        for name in names
    )
    set_indices = {index for index, result in enumerate(results) if result.hits}
    feature_names = tuple(f"panel:{name}" for name in names)
    gaps = tuple(part for part in sequence.parts if isinstance(part, Gap))
    return BitFingerprintResult(
        name="panel_fingerprint",
        method=f"named-{mode}-pattern-presence",
        schema_version="dnakit.panel-presence-bit.v1",
        sequence_id=sequence_id,
        symbol_length=sequence.symbol_length,
        gap_count=len(gaps),
        unknown_gap_count=sum(gap.length is None for gap in gaps),
        representation=resolved_representation,
        feature_names=feature_names,
        values=_materialize_bits(feature_names, set_indices, resolved_representation),
        parameters=FrozenDict(
            {
                "panel": {name: panel[name] for name in names},
                "mode": mode,
                "strand": "both",
                "merge_strands": True,
                "overlapping": overlapping,
            }
        ),
        observation_count=sum(len(result.hits) for result in results),
        ignored_ambiguity_count=0,
        max_dimension=max_panel_size,
    )


__all__ = [
    "DEFAULT_HASHED_KMER_BITS",
    "DEFAULT_MAX_PANEL_SIZE",
    "MAX_BIT_FINGERPRINT_DIMENSION",
    "PanelMode",
    "hashed_kmer_fingerprint",
    "panel_fingerprint",
]
