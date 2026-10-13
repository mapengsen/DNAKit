"""Sourmash sketch and k-mer helpers for microbial DNA, one wrapper per upstream operation.

The functions in this module call sourmash directly on sketches that are passed in as plain
mappings with ``ksize``, ``scaled`` and ``hashes`` (and optional ``id``), such as the items of
``dnakit.microbiome.sketch_fasta_minhash()["sketches"]``. Nothing is downloaded and no file is
read. Sketches are FracMinHash sketches over canonical DNA k-mers; coordinates are not available.

sourmash is imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-sourmash]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_SOURMASH_EXTRA = "bio-sourmash"
_HASH_LIMIT = 2**64 - 1
_KSIZE_LIMIT = 100
_SCALED_LIMIT = 1_000_000_000
_SET_SIZE_LIMIT = 10**12
_MATRIX_LIMIT = 500
_ACGT = frozenset("ACGT")
_COMPARE_METRICS = {
    "jaccard": "compare_serial",
    "containment": "compare_serial_containment",
    "max_containment": "compare_serial_max_containment",
    "avg_containment": "compare_serial_avg_containment",
}
_MEASURES = ("jaccard", "containment")


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _real(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"{name} must be a finite number.")
    return float(value)


def _dna(value: object, name: str, max_length: int) -> str:
    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = value.upper()
    if set(text) - _ACGT:
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    if len(text) > max_length:
        raise ConfigurationError(f"{name} length must be at most {max_length}.")
    return text


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _sketch_fields(sketch: object, name: str) -> tuple[str | None, int, int, list[int]]:
    if not isinstance(sketch, Mapping):
        raise ConfigurationError(
            f"{name} must be a sketch mapping, such as one item of "
            "sketch_fasta_minhash()['sketches']."
        )
    ksize = _integer(sketch.get("ksize"), f"{name} ksize", 1, _KSIZE_LIMIT)
    scaled = _integer(sketch.get("scaled"), f"{name} scaled", 1, _SCALED_LIMIT)
    hashes = sketch.get("hashes")
    if isinstance(hashes, (str, bytes)) or not isinstance(hashes, Sequence) or not hashes:
        raise ConfigurationError(f"{name} hashes must be a non-empty list of integers.")
    values: list[int] = []
    for value in hashes:
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _HASH_LIMIT:
            raise ConfigurationError(f"{name} hashes must be unsigned 64-bit integers.")
        values.append(value)
    identifier = sketch.get("id")
    return (
        None if identifier is None else str(identifier),
        ksize,
        scaled,
        sorted(set(values)),
    )


def _minhash(sourmash: Any, ksize: int, scaled: int, hashes: Sequence[int]) -> Any:
    """Rebuild a sourmash MinHash from explicit hashes; hashes above the scaled cut are rejected."""

    try:
        minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        minhash.add_many(hashes)
    except (AttributeError, OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not build a MinHash: {exc}") from exc
    if len(minhash) != len(set(hashes)):
        raise ConfigurationError("some hashes lie outside the scaled range of their sketch.")
    return minhash


def _signature(sourmash: Any, minhash: Any, name: str) -> Any:
    return sourmash.signature.SourmashSignature(minhash, name=name)


def compare_minhash_matrix(
    sketches: Sequence[Mapping[str, Any]],
    *,
    metric: str = "jaccard",
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise similarity matrix of scaled MinHash sketches, computed with sourmash.

    ``metric="jaccard"`` calls ``sourmash.compare.compare_serial``. The containment metrics call
    ``compare_serial_containment``, ``compare_serial_max_containment`` or
    ``compare_serial_avg_containment``. Entry ``matrix[i][j]`` compares sketch ``i`` with sketch
    ``j``. Containment is directional: it is the fraction of sketch ``i`` found in sketch ``j``.
    """

    if metric not in _COMPARE_METRICS:
        raise ConfigurationError(f"metric must be one of: {', '.join(_COMPARE_METRICS)}.")
    if (
        isinstance(sketches, (str, bytes, Mapping))
        or not isinstance(sketches, Sequence)
        or len(sketches) < 2
    ):
        raise ConfigurationError("sketches must be a list of at least two sketch mappings.")
    if len(sketches) > _MATRIX_LIMIT:
        raise ConfigurationError(f"at most {_MATRIX_LIMIT} sketches can be compared at once.")
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    compare_module = require("sourmash.compare", _SOURMASH_EXTRA)
    fields = [_sketch_fields(item, f"sketches[{index}]") for index, item in enumerate(sketches)]
    ksize, scaled = fields[0][1], fields[0][2]
    if any((item[1], item[2]) != (ksize, scaled) for item in fields):
        raise ConfigurationError("sketches must share ksize and scaled to be compared.")
    signatures = [
        _signature(
            sourmash,
            _minhash(sourmash, ksize, scaled, hashes),
            identifier or f"sketch_{index}",
        )
        for index, (identifier, _k, _s, hashes) in enumerate(fields)
    ]
    try:
        if metric == "jaccard":
            matrix = compare_module.compare_serial(signatures, True, downsample=False)
        else:
            compare = getattr(compare_module, _COMPARE_METRICS[metric])
            # sourmash stores containments[i][j] = sketch j contained in sketch i; transpose so
            # that the wrapper reports the fraction of sketch i found in sketch j.
            matrix = compare(signatures, downsample=False).T
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash comparison failed: {exc}") from exc
    rows = [[float(value) for value in row] for row in matrix.tolist()]
    return {
        "status": "completed",
        "operation": "compare_minhash_matrix",
        "metric": metric,
        "ids": [signature.name for signature in signatures],
        "ksize": ksize,
        "scaled": scaled,
        "matrix": rows,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": (
                "set Jaccard from compare_serial"
                if metric == "jaccard"
                else f"directional containment from {_COMPARE_METRICS[metric]}"
            ),
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": (
                "matrix[i][j] compares sketch i with sketch j; containment is the fraction of "
                "sketch i found in sketch j"
            ),
        },
    }


def estimate_mutation_rate(
    value: float,
    *,
    measure: str = "jaccard",
    ksize: int = 21,
    scaled: int = 1_000,
    n_unique_kmers: int | None = None,
    sequence_len_bp: int | None = None,
    confidence: float = 0.95,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the per-base mutation rate from a sketch Jaccard index or containment.

    ``measure="jaccard"`` calls ``sourmash.minhash.jaccard_to_distance``; ``measure="containment"``
    calls ``sourmash.minhash.containment_to_distance`` with a confidence interval. The estimate
    needs the sequence length in bp or the number of distinct k-mers. It is an approximation, and
    the result reports whether the sketch size is likely too small for it.
    """

    if measure not in _MEASURES:
        raise ConfigurationError(f"measure must be one of: {', '.join(_MEASURES)}.")
    point = _real(value, "value")
    if not 0.0 < point < 1.0:
        raise ConfigurationError("value must be strictly between 0 and 1.")
    k = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    s = _integer(scaled, "scaled", 1, _SCALED_LIMIT)
    if n_unique_kmers is None and sequence_len_bp is None:
        raise ConfigurationError("give sequence_len_bp or n_unique_kmers to estimate the rate.")
    if n_unique_kmers is not None:
        n_unique_kmers = _integer(n_unique_kmers, "n_unique_kmers", 1, _SET_SIZE_LIMIT)
    if sequence_len_bp is not None:
        sequence_len_bp = _integer(sequence_len_bp, "sequence_len_bp", 1, _SET_SIZE_LIMIT)
    level = _real(confidence, "confidence")
    if not 0.0 < level < 1.0:
        raise ConfigurationError("confidence must be strictly between 0 and 1.")
    minhash = require("sourmash.minhash", _SOURMASH_EXTRA)
    try:
        if measure == "jaccard":
            result = minhash.jaccard_to_distance(
                point,
                k,
                s,
                n_unique_kmers=n_unique_kmers,
                sequence_len_bp=sequence_len_bp,
            )
            low = high = None
            jaccard_error = float(result.jaccard_error)
        else:
            result = minhash.containment_to_distance(
                point,
                k,
                s,
                n_unique_kmers=n_unique_kmers,
                sequence_len_bp=sequence_len_bp,
                confidence=level,
                estimate_ci=True,
            )
            low = None if result.dist_low is None else float(result.dist_low)
            high = None if result.dist_high is None else float(result.dist_high)
            jaccard_error = None
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"sourmash could not estimate the rate: {exc}") from exc
    distance = float(result.dist)
    return {
        "status": "completed",
        "operation": "estimate_mutation_rate",
        "measure": measure,
        "input_value": point,
        "ksize": k,
        "scaled": s,
        "distance": distance,
        "ani": 1.0 - distance,
        "distance_low": low,
        "distance_high": high,
        "confidence": level if measure == "containment" else None,
        "jaccard_error": jaccard_error,
        "p_nothing_in_common": float(result.p_nothing_in_common),
        "size_is_inaccurate": bool(result.size_is_inaccurate),
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": (
                "point estimate from the Jaccard index (sourmash.minhash.jaccard_to_distance)"
                if measure == "jaccard"
                else "point estimate and interval from containment "
                "(sourmash.minhash.containment_to_distance)"
            ),
            "coordinate_system": "per-base mutation rate; no genomic coordinates",
            "convention": "ani is 1 - distance; sizes are approximate and need sequence_len_bp "
            "or n_unique_kmers",
        },
    }


def set_size_error_probability(
    set_size: int,
    scaled: int,
    *,
    relative_error: float = 0.05,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability that a scaled sketch estimates a k-mer set size within a relative tolerance.

    Calls ``sourmash.minhash.set_size_exact_prob``. The sketch size is binomial with parameters
    ``set_size`` and ``1 / scaled``. The returned value is P(|sketch_size * scaled - set_size| <=
    relative_error * set_size), the probability of an estimate inside the tolerance. The upstream
    docstring describes the complement (a deviation), but the computed value is this one.
    """

    n = _integer(set_size, "set_size", 1, _SET_SIZE_LIMIT)
    s = _integer(scaled, "scaled", 1, _SCALED_LIMIT)
    tolerance = _real(relative_error, "relative_error")
    if not 0.0 < tolerance < 1.0:
        raise ConfigurationError("relative_error must be strictly between 0 and 1.")
    minhash = require("sourmash.minhash", _SOURMASH_EXTRA)
    try:
        probability = float(minhash.set_size_exact_prob(n, s, relative_error=tolerance))
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"sourmash could not compute the probability: {exc}") from exc
    return {
        "status": "completed",
        "operation": "set_size_error_probability",
        "set_size": n,
        "scaled": s,
        "relative_error": tolerance,
        "expected_sketch_size": n / s,
        "probability_within_tolerance": probability,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "exact binomial probability (sourmash.minhash.set_size_exact_prob)",
            "coordinate_system": "k-mer set size; no genomic coordinates",
            "convention": "probability that |sketch_size * scaled - set_size| <= "
            "relative_error * set_size, boundary included",
        },
    }


def hash_canonical_kmer(
    kmer: str,
    *,
    seed: int = 42,
    verbose: bool = True,
) -> dict[str, Any]:
    """MurmurHash of the canonical form of a DNA k-mer, as sourmash hashes it in a sketch.

    The canonical form is the lexicographically smaller of the k-mer and its reverse complement.
    The hash is ``sourmash.minhash.hash_murmur(canonical, seed)``, which is the value a scaled=1
    sketch stores for the k-mer. Lowercase input is uppercased first.
    """

    text = _dna(kmer, "kmer", _KSIZE_LIMIT)
    value = _integer(seed, "seed", 0, _HASH_LIMIT)
    reverse = _reverse_complement(text)
    canonical = min(text, reverse)
    minhash = require("sourmash.minhash", _SOURMASH_EXTRA)
    hashed = int(minhash.hash_murmur(canonical, value))
    if reverse == text:
        strand = "self_reverse_complement"
    elif canonical == text:
        strand = "forward"
    else:
        strand = "reverse_complement"
    return {
        "status": "completed",
        "operation": "hash_canonical_kmer",
        "kmer": text,
        "canonical": canonical,
        "strand": strand,
        "ksize": len(text),
        "seed": value,
        "hash": hashed,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "hash": "MurmurHash64 of the canonical k-mer string (sourmash.minhash.hash_murmur)",
            "coordinate_system": "k-mer hash; no genomic coordinates",
            "convention": "the canonical form is the smaller of the k-mer and its reverse "
            "complement as strings",
        },
    }


def downsample_minhash_sketch(
    sketch: Mapping[str, Any],
    *,
    scaled: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Downsample a scaled sketch to a coarser scale with sourmash.

    Calls ``sourmash.minhash.flatten_and_downsample_scaled``. The target ``scaled`` must be at
    least the sketch's own scale; sketches cannot be upsampled. Hashes above sourmash's scaled cut
    are dropped. The cut is computed in floating point, so a hash within a few hundred units of
    ``2**64 / scaled`` can be decided either way; the output is sourmash's own decision.
    """

    identifier, ksize, source_scaled, hashes = _sketch_fields(sketch, "sketch")
    target = _integer(scaled, "scaled", 1, _SCALED_LIMIT)
    if target < source_scaled:
        raise ConfigurationError(
            f"scaled={target} is finer than the sketch scaled={source_scaled}; "
            "sketches cannot be upsampled."
        )
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    minhash = require("sourmash.minhash", _SOURMASH_EXTRA)
    original = _minhash(sourmash, ksize, source_scaled, hashes)
    try:
        result = minhash.flatten_and_downsample_scaled(original, target)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"sourmash could not downsample the sketch: {exc}") from exc
    kept = sorted(int(value) for value in result.hashes)
    return {
        "status": "completed",
        "operation": "downsample_minhash_sketch",
        "id": identifier,
        "ksize": ksize,
        "input_scaled": source_scaled,
        "scaled": int(result.scaled),
        "input_hash_count": len(hashes),
        "hash_count": len(kept),
        "hashes": kept,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "flatten_and_downsample_scaled: keep hashes at or below the new cut",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "the output scale is the requested scale; abundance is not kept",
        },
    }


def intersect_minhash_sketches(
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hashes shared by two scaled sketches, after both are downsampled to the coarser scale.

    Calls ``sourmash.minhash.flatten_and_intersect_scaled``. The two sketches must share ``ksize``;
    their scales may differ. The result is at the coarser of the two scales.
    """

    left_id, left_k, left_s, left_hashes = _sketch_fields(left, "left")
    right_id, right_k, right_s, right_hashes = _sketch_fields(right, "right")
    if left_k != right_k:
        raise ConfigurationError("sketches must share ksize to be intersected.")
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    minhash = require("sourmash.minhash", _SOURMASH_EXTRA)
    left_mh = _minhash(sourmash, left_k, left_s, left_hashes)
    right_mh = _minhash(sourmash, right_k, right_s, right_hashes)
    try:
        shared = minhash.flatten_and_intersect_scaled(left_mh, right_mh)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ConfigurationError(f"sourmash could not intersect the sketches: {exc}") from exc
    hashes = sorted(int(value) for value in shared.hashes)
    return {
        "status": "completed",
        "operation": "intersect_minhash_sketches",
        "left_id": left_id,
        "right_id": right_id,
        "ksize": left_k,
        "scaled": int(shared.scaled),
        "hash_count": len(hashes),
        "hashes": hashes,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "flatten_and_intersect_scaled: downsample both, then intersect hashes",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "the result is at the coarser of the two scales",
        },
    }


def search_minhash_database(
    query: Mapping[str, Any],
    database: Sequence[Mapping[str, Any]],
    *,
    threshold_bp: int = 5_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find database sketches that share at least ``threshold_bp`` of k-mers with a query sketch.

    Calls ``sourmash.search.prefetch_database`` on a ``sourmash.index.LinearIndex`` built from the
    database sketches. Overlap in bp is the number of shared hashes times the coarser scale. Matches
    are sorted by overlap in bp (largest first), then by id, so the output is deterministic.
    A threshold above the query's own bp is unattainable, and sourmash raises for it.
    """

    query_id, query_k, query_s, query_hashes = _sketch_fields(query, "query")
    threshold = _integer(threshold_bp, "threshold_bp", 1, 10**12)
    if isinstance(database, (str, bytes, Mapping)) or not isinstance(database, Sequence):
        raise ConfigurationError("database must be a list of sketch mappings.")
    if not database:
        raise ConfigurationError("database must contain at least one sketch.")
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    search = require("sourmash.search", _SOURMASH_EXTRA)
    index_module = require("sourmash.index", _SOURMASH_EXTRA)
    signatures = []
    for index, item in enumerate(database):
        identifier, item_k, item_s, item_hashes = _sketch_fields(item, f"database[{index}]")
        if item_k != query_k:
            raise ConfigurationError("database sketches must share the query ksize.")
        minhash = _minhash(sourmash, item_k, item_s, item_hashes)
        signatures.append(_signature(sourmash, minhash, identifier or f"database_{index}"))
    query_minhash = _minhash(sourmash, query_k, query_s, query_hashes)
    query_signature = _signature(sourmash, query_minhash, query_id or "query")
    try:
        results = list(
            search.prefetch_database(
                query_signature,
                index_module.LinearIndex(signatures),
                threshold,
            )
        )
    except (AssertionError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash search failed: {exc}") from exc
    matches = [
        {
            "id": str(result.name),
            "md5": str(result.md5),
            "intersect_bp": int(result.intersect_bp),
            "jaccard": float(result.jaccard),
            "fraction_of_query_in_match": float(result.f_match_query),
            "fraction_of_match_in_query": float(result.f_query_match),
            "max_containment": float(result.max_containment),
        }
        for result in results
    ]
    matches.sort(key=lambda item: (-item["intersect_bp"], item["id"], item["md5"]))
    return {
        "status": "completed",
        "operation": "search_minhash_database",
        "query_id": query_id,
        "ksize": query_k,
        "threshold_bp": threshold,
        "database_size": len(signatures),
        "match_count": len(matches),
        "matches": matches,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "prefetch: containment search over a LinearIndex (sourmash.search)",
            "coordinate_system": "k-mer hashes; overlap is reported in bp, not coordinates",
            "convention": "intersect_bp is shared hashes times the coarser scale; "
            "fraction_of_query_in_match is the query share found in the match, "
            "fraction_of_match_in_query is the reverse",
        },
    }


def _lineage_string(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty ';'-separated lineage string.")
    if "," in value:
        raise ConfigurationError(f"{name} must use ';' as the separator, not ','.")
    parts = [part.strip() for part in value.split(";")]
    if len(parts) > 8:
        raise ConfigurationError(f"{name} has more than 8 ranks (superkingdom to strain).")
    if any(not part for part in parts):
        raise ConfigurationError(f"{name} contains an empty taxon name.")
    return ";".join(parts)


def lowest_common_lineage(
    lineages: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Lowest common ancestor of taxonomic lineages, using sourmash's lineage tree helpers.

    Each lineage is a ``;``-separated string from superkingdom downwards, such as
    ``"d__Bacteria;p__Firmicutes;c__Bacilli"``. Calls ``sourmash.lca.lca_utils.make_lineage``,
    ``build_tree`` and ``find_lca``. ``distinct_children_below_lca`` is 0 when every lineage agrees
    to its full depth, and otherwise the number of different taxa at the next rank.
    """

    if isinstance(lineages, (str, bytes, Mapping)) or not isinstance(lineages, Sequence):
        raise ConfigurationError("lineages must be a list of lineage strings.")
    if not lineages:
        raise ConfigurationError("lineages must contain at least one lineage.")
    texts = [_lineage_string(item, f"lineages[{i}]") for i, item in enumerate(lineages)]
    lca_utils = require("sourmash.lca.lca_utils", _SOURMASH_EXTRA)
    try:
        assignments = [lca_utils.make_lineage(text) for text in texts]
        tree = lca_utils.build_tree(assignments)
        lca, children = lca_utils.find_lca(tree)
    except (AttributeError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not resolve the lineages: {exc}") from exc
    ranks = [{"rank": str(item.rank), "name": str(item.name)} for item in lca]
    return {
        "status": "completed",
        "operation": "lowest_common_lineage",
        "assignment_count": len(texts),
        "lineage": ranks,
        "lineage_string": ";".join(item["name"] for item in ranks),
        "lowest_rank": ranks[-1]["rank"] if ranks else None,
        "distinct_children_below_lca": int(children),
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "algorithm": "lineage tree built with sourmash.lca.lca_utils.build_tree and find_lca",
            "coordinate_system": "taxonomic ranks; no genomic coordinates",
            "convention": "ranks are superkingdom, phylum, class, order, family, genus, species, "
            "strain, in that order",
        },
    }


__all__ = [
    "compare_minhash_matrix",
    "downsample_minhash_sketch",
    "estimate_mutation_rate",
    "hash_canonical_kmer",
    "intersect_minhash_sketches",
    "lowest_common_lineage",
    "search_minhash_database",
    "set_size_error_probability",
]

_register_public_verbose(globals(), tuple(__all__))
