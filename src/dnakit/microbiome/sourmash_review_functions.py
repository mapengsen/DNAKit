"""sourmash mutation-model, threshold and k-mer LCA helpers for DNA, one wrapper per upstream call.

These wrappers cover the DNA-level sourmash functions that the existing ``sourmash_functions``
module does not expose: the k-mer mutation model behind the mutation-rate estimate, the conversion
of a bp overlap threshold into hashes, and the lowest common ancestor of k-mer assignments. They
take plain numbers, DNA lengths or local lineage strings, and nothing is read from disk or from
the network.

sourmash is imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-sourmash]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_SOURMASH_EXTRA = "bio-sourmash"
_KSIZE_LIMIT = 100
_SCALED_LIMIT = 1_000_000_000
_COUNT_LIMIT = 10**12
_HASH_LIMIT = 2**64 - 1
_HASH_COUNT_LIMIT = 100_000
_LINEAGE_RANKS = 8


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _rate(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not 0.0 <= value <= 1.0
    ):
        raise ConfigurationError(f"{name} must be a number in [0, 1].")
    return float(value)


def _mutation_inputs(
    ksize: object,
    mutation_rate: object,
    sequence_len_bp: object,
    n_unique_kmers: object,
) -> tuple[int, float, int]:
    """Validate k, the per-base mutation rate and the number of distinct k-mers."""

    k = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    rate = _rate(mutation_rate, "mutation_rate")
    if (sequence_len_bp is None) == (n_unique_kmers is None):
        raise ConfigurationError("give exactly one of sequence_len_bp or n_unique_kmers.")
    if n_unique_kmers is not None:
        return k, rate, _integer(n_unique_kmers, "n_unique_kmers", 1, _COUNT_LIMIT)
    length = _integer(sequence_len_bp, "sequence_len_bp", k, _COUNT_LIMIT)
    return k, rate, length - (k - 1)


def _lineage_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty ';'-separated lineage string.")
    if "," in value:
        raise ConfigurationError(f"{name} must use ';' as the separator, not ','.")
    parts = [part.strip() for part in value.split(";")]
    if len(parts) > _LINEAGE_RANKS:
        raise ConfigurationError(f"{name} has more than {_LINEAGE_RANKS} ranks.")
    if any(not part for part in parts):
        raise ConfigurationError(f"{name} contains an empty taxon name.")
    return ";".join(parts)


def kmer_mutation_probability(
    mutation_rate: float,
    ksize: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability that a k-mer contains at least one mutated base, from a per-base rate.

    Calls ``sourmash.distance_utils.r1_to_q``, which computes ``1 - (1 - mutation_rate) ** ksize``.
    The bases are assumed to mutate independently. ``mutation_rate`` is a fraction in [0, 1].
    """

    k = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    rate = _rate(mutation_rate, "mutation_rate")
    distance_utils = require("sourmash.distance_utils", _SOURMASH_EXTRA)
    probability = float(distance_utils.r1_to_q(k, rate))
    return {
        "status": "completed",
        "operation": "kmer_mutation_probability",
        "mutation_rate": rate,
        "ksize": k,
        "kmer_mutation_probability": probability,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "sourmash.distance_utils.r1_to_q: 1 - (1 - mutation_rate) ** ksize",
            "coordinate_system": "per-base mutation rate; no genomic coordinates",
            "convention": "the probability is per k-mer and assumes independent base mutations",
        },
    }


def expected_mutated_kmers(
    ksize: int,
    mutation_rate: float,
    *,
    sequence_len_bp: int | None = None,
    n_unique_kmers: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected number of distinct k-mers of a sequence that a mutation rate changes.

    Calls ``sourmash.distance_utils.exp_n_mutated`` on ``n_unique_kmers`` (distinct k-mers, which
    is ``sequence_len_bp - ksize + 1`` for a sequence without repeats). Give exactly one of
    ``sequence_len_bp`` and ``n_unique_kmers``. The result is an expectation, not a count.
    """

    k, rate, n = _mutation_inputs(ksize, mutation_rate, sequence_len_bp, n_unique_kmers)
    distance_utils = require("sourmash.distance_utils", _SOURMASH_EXTRA)
    expected = float(distance_utils.exp_n_mutated(n, k, rate))
    return {
        "status": "completed",
        "operation": "expected_mutated_kmers",
        "ksize": k,
        "mutation_rate": rate,
        "n_unique_kmers": n,
        "kmer_mutation_probability": float(distance_utils.r1_to_q(k, rate)),
        "expected_mutated_kmers": expected,
        "expected_unmutated_kmers": float(n) - expected,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "sourmash.distance_utils.exp_n_mutated: n_unique_kmers times the "
            "per-k-mer mutation probability",
            "coordinate_system": "distinct k-mers; no genomic coordinates",
            "convention": "sequence_len_bp - ksize + 1 distinct k-mers; repeats are not collapsed",
        },
    }


def variance_mutated_kmers(
    ksize: int,
    mutation_rate: float,
    *,
    sequence_len_bp: int | None = None,
    n_unique_kmers: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Variance of the number of mutated distinct k-mers, from the same mutation model.

    Calls ``sourmash.distance_utils.var_n_mutated``. Give exactly one of ``sequence_len_bp`` and
    ``n_unique_kmers``. Very small rates can make the closed form negative in floating point;
    sourmash raises for that case and this wrapper reports it as a configuration error.
    """

    k, rate, n = _mutation_inputs(ksize, mutation_rate, sequence_len_bp, n_unique_kmers)
    distance_utils = require("sourmash.distance_utils", _SOURMASH_EXTRA)
    try:
        variance = float(distance_utils.var_n_mutated(n, k, rate))
    except ValueError as exc:
        raise ConfigurationError(f"sourmash could not compute the variance: {exc}") from exc
    return {
        "status": "completed",
        "operation": "variance_mutated_kmers",
        "ksize": k,
        "mutation_rate": rate,
        "n_unique_kmers": n,
        "variance_mutated_kmers": variance,
        "standard_deviation_mutated_kmers": math.sqrt(variance),
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "sourmash.distance_utils.var_n_mutated (closed-form variance of the "
            "mutated k-mer count)",
            "coordinate_system": "distinct k-mers; no genomic coordinates",
            "convention": "variance is for the count of mutated distinct k-mers, not the "
            "sketch estimate",
        },
    }


def probability_no_shared_kmers(
    mutation_rate: float,
    ksize: int,
    scaled: int,
    *,
    sequence_len_bp: int | None = None,
    n_unique_kmers: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability that a scaled sketch of a sequence and of its mutated copy share no hash.

    Calls ``sourmash.distance_utils.get_exp_probability_nothing_common``. The estimate is
    ``exp((n - E[mutated]) * ln(1 - 1/scaled))``. A rate of 0 gives 0 and a rate of 1 gives 1.
    ``scaled`` must be at least 2, because a scale of 1 keeps every k-mer and the logarithm is
    undefined. Give exactly one of ``sequence_len_bp`` and ``n_unique_kmers``.
    """

    k, rate, n = _mutation_inputs(ksize, mutation_rate, sequence_len_bp, n_unique_kmers)
    s = _integer(scaled, "scaled", 2, _SCALED_LIMIT)
    distance_utils = require("sourmash.distance_utils", _SOURMASH_EXTRA)
    probability = float(
        distance_utils.get_exp_probability_nothing_common(rate, k, s, n_unique_kmers=n)
    )
    if not math.isfinite(probability) or not 0.0 <= probability <= 1.0:
        raise ConfigurationError("sourmash returned a probability outside [0, 1].")
    return {
        "status": "completed",
        "operation": "probability_no_shared_kmers",
        "mutation_rate": rate,
        "ksize": k,
        "scaled": s,
        "n_unique_kmers": n,
        "probability_nothing_in_common": probability,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "sourmash.distance_utils.get_exp_probability_nothing_common",
            "coordinate_system": "scaled k-mer sketch; no genomic coordinates",
            "convention": "approximate expectation; a large value means the sketches of the "
            "sequence and its mutated copy may share nothing",
        },
    }


def bp_threshold_to_hashes(
    threshold_bp: int,
    scaled: int,
    query_size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a bp overlap threshold into a hash count and a containment fraction of the query.

    Calls ``sourmash.search.calc_threshold_from_bp``. ``threshold_hashes`` is
    ``threshold_bp / scaled`` and ``containment_threshold`` is that count divided by
    ``query_size``, the number of hashes in the query sketch. A threshold above the query's own bp
    cannot be reached, and sourmash raises for it; this wrapper reports that as a configuration
    error.
    """

    bp = _integer(threshold_bp, "threshold_bp", 0, _COUNT_LIMIT)
    s = _integer(scaled, "scaled", 1, _SCALED_LIMIT)
    size = _integer(query_size, "query_size", 1, _COUNT_LIMIT)
    search = require("sourmash.search", _SOURMASH_EXTRA)
    try:
        fraction, hashes = search.calc_threshold_from_bp(bp, s, size)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not convert the threshold: {exc}") from exc
    return {
        "status": "completed",
        "operation": "bp_threshold_to_hashes",
        "threshold_bp": bp,
        "scaled": s,
        "query_size": size,
        "threshold_hashes": float(hashes),
        "containment_threshold": float(fraction),
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "sourmash.search.calc_threshold_from_bp: threshold_bp / scaled hashes, "
            "divided by the query size",
            "coordinate_system": "k-mer hashes; overlap in bp, no genomic coordinates",
            "convention": "threshold_hashes may be fractional; containment_threshold is in [0, 1]",
        },
    }


def _hash_key(key: object) -> object:
    """Turn a digit string into an int so hash keys survive a JSON round trip."""

    return int(key) if isinstance(key, str) and key.isdigit() else key


def count_kmer_lca_assignments(
    assignments: Mapping[str, Sequence[str]],
    *,
    hashval_counts: Mapping[str, int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count, for each k-mer hash, the lowest common ancestor of its taxonomic assignments.

    ``assignments`` maps a 64-bit hash to one or more ``;``-separated lineage strings. Calls
    ``sourmash.lca.lca_utils.count_lca_for_assignments``. ``hashval_counts``, if given, must map
    exactly the same hashes to positive integer weights. An empty lineage in the output is the
    root, where the assignments disagree at the superkingdom rank. Counts are sorted by count
    (largest first), then by lineage.
    """

    if not isinstance(assignments, Mapping) or not assignments:
        raise ConfigurationError("assignments must be a non-empty mapping of hash to lineages.")
    if len(assignments) > _HASH_COUNT_LIMIT:
        raise ConfigurationError(f"at most {_HASH_COUNT_LIMIT} hashes can be counted at once.")
    # JSON callers send hashes as digit strings; integers are accepted as well.
    assignments = {_hash_key(key): value for key, value in assignments.items()}
    if isinstance(hashval_counts, Mapping):
        hashval_counts = {_hash_key(key): value for key, value in hashval_counts.items()}
    checked: dict[int, list[str]] = {}
    for key, lineages in assignments.items():
        identifier = _integer(key, "assignment hash", 0, _HASH_LIMIT)
        if isinstance(lineages, (str, bytes, Mapping)) or not isinstance(lineages, Sequence):
            raise ConfigurationError(f"lineages for hash {identifier} must be a list of strings.")
        if not lineages:
            raise ConfigurationError(f"hash {identifier} has no lineage assignments.")
        checked[identifier] = [_lineage_text(item, f"hash {identifier}") for item in lineages]
    weights: dict[int, int] | None = None
    if hashval_counts is not None:
        if not isinstance(hashval_counts, Mapping) or set(hashval_counts) != set(checked):
            raise ConfigurationError("hashval_counts must map exactly the hashes in assignments.")
        weights = {
            key: _integer(hashval_counts[key], f"count for hash {key}", 1, _COUNT_LIMIT)
            for key in sorted(checked)
        }
    lca_utils = require("sourmash.lca.lca_utils", _SOURMASH_EXTRA)
    try:
        parsed = {
            key: [lca_utils.make_lineage(text) for text in checked[key]] for key in sorted(checked)
        }
        counter: Counter[Any] = lca_utils.count_lca_for_assignments(parsed, weights)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not resolve the assignments: {exc}") from exc
    counts = []
    for lineage, count in counter.items():
        names = [str(item.name) for item in lineage]
        counts.append(
            {
                "lineage": ";".join(names),
                "depth": len(names),
                "lowest_rank": str(lineage[-1].rank) if lineage else None,
                "count": int(count),
            }
        )
    counts.sort(key=lambda item: (-item["count"], item["lineage"]))
    return {
        "status": "completed",
        "operation": "count_kmer_lca_assignments",
        "hash_count": len(checked),
        "weighted": weights is not None,
        "total_count": sum(item["count"] for item in counts),
        "counts": counts,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "lowest common ancestor per hash over its lineages, with "
            "sourmash.lca.lca_utils.count_lca_for_assignments",
            "coordinate_system": "taxonomic ranks; hashes are k-mer hashes, not coordinates",
            "convention": "ranks are superkingdom to strain; an empty lineage is the root",
        },
    }


__all__ = [
    "bp_threshold_to_hashes",
    "count_kmer_lca_assignments",
    "expected_mutated_kmers",
    "kmer_mutation_probability",
    "probability_no_shared_kmers",
    "variance_mutated_kmers",
]

_register_public_verbose(globals(), tuple(__all__))
