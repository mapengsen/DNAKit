"""Immutable structured results for native similarity operations."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Literal, TypeAlias, cast

from dnakit.core import Provenance, Strand
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError

DistanceKind: TypeAlias = Literal["hamming", "levenshtein"]
EditOperation: TypeAlias = Literal["match", "substitute", "insert", "delete"]
ValueKind: TypeAlias = Literal["similarity", "distance"]
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class _SerializableResult:
    def to_dict(self) -> dict[str, Any]:
        """Return a plain JSON-compatible representation."""

        return cast(dict[str, Any], to_json_compatible(self))


def _non_negative_int(value: object, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer.")


def _optional_identifier(value: str | None, name: str) -> None:
    if value is not None and (not isinstance(value, str) or not value.strip()):
        raise ConfigurationError(f"{name} must be a non-empty string or None.")


@dataclass(frozen=True)
class SequenceMatch(_SerializableResult):
    """One exact hit in target 0-based, half-open coordinates."""

    target_index: int
    target_id: str | None
    start: int
    end: int
    strand: Strand

    def __post_init__(self) -> None:
        _non_negative_int(self.target_index, "target_index")
        _optional_identifier(self.target_id, "target_id")
        _non_negative_int(self.start, "start")
        _non_negative_int(self.end, "end")
        if self.end < self.start:
            raise ConfigurationError("SequenceMatch end cannot be smaller than start.")
        if not isinstance(self.strand, Strand):
            raise ConfigurationError("SequenceMatch strand must be Strand.")


@dataclass(frozen=True)
class ApproximateMatch(SequenceMatch):
    """One approximate substring hit with its weighted edit distance."""

    distance: float

    def __post_init__(self) -> None:
        super().__post_init__()
        if (
            isinstance(self.distance, bool)
            or not isinstance(self.distance, (int, float))
            or not math.isfinite(self.distance)
            or self.distance < 0
        ):
            raise ConfigurationError("Approximate match distance must be non-negative and finite.")


@dataclass(frozen=True)
class SearchResult(_SerializableResult):
    """Exact full-sequence or exact-subsequence search result."""

    name: str
    method: str
    query_id: str | None
    query_length: int
    target_count: int
    matches: tuple[SequenceMatch, ...]
    overlapping: bool
    reverse_complement: bool
    merge_strands: bool
    full_length: bool
    max_targets: int
    max_matches: int
    iupac_matching: Literal["literal"] = "literal"
    coordinate_system: Literal["0-based-half-open"] = "0-based-half-open"
    circular_wrap: bool = False

    def __post_init__(self) -> None:
        for name in ("name", "method"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"SearchResult {name} must be non-empty.")
        _optional_identifier(self.query_id, "query_id")
        _non_negative_int(self.query_length, "query_length")
        _non_negative_int(self.target_count, "target_count")
        if any(not isinstance(match, SequenceMatch) for match in self.matches):
            raise ConfigurationError("SearchResult matches must all be SequenceMatch objects.")
        if any(match.target_index >= self.target_count for match in self.matches):
            raise ConfigurationError("A search match target_index is outside target_count.")
        for name in ("overlapping", "reverse_complement", "merge_strands", "full_length"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"SearchResult {name} must be boolean.")
        for name in ("max_targets", "max_matches"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ConfigurationError(f"SearchResult {name} must be a positive integer.")
        if self.target_count > self.max_targets:
            raise ConfigurationError("SearchResult target_count cannot exceed max_targets.")
        if len(self.matches) > self.max_matches:
            raise ConfigurationError("SearchResult matches cannot exceed max_matches.")
        if self.iupac_matching != "literal":
            raise ConfigurationError("Only literal IUPAC matching is implemented in the MVP.")
        if self.coordinate_system != "0-based-half-open":
            raise ConfigurationError("Internal search coordinates must be 0-based half-open.")
        if self.circular_wrap:
            raise ConfigurationError("Circular origin-wrapping search is not implemented in MVP.")

    @property
    def found(self) -> bool:
        return bool(self.matches)


@dataclass(frozen=True)
class ApproximateSearchResult(_SerializableResult):
    """Audited collection of bounded approximate substring matches."""

    name: str
    method: str
    algorithm_version: str
    query_id: str | None
    query_length: int
    target_count: int
    matches: tuple[ApproximateMatch, ...]
    max_distance: float
    reverse_complement: bool
    max_targets: int
    max_matches: int
    max_cells: int
    dp_cells: int
    parameters: FrozenDict

    def __post_init__(self) -> None:
        for name in ("name", "method", "algorithm_version"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ConfigurationError(f"ApproximateSearchResult {name} must be non-empty.")
        _optional_identifier(self.query_id, "query_id")
        for name in ("query_length", "target_count", "dp_cells"):
            _non_negative_int(getattr(self, name), name)
        for name in ("max_targets", "max_matches", "max_cells"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ConfigurationError(f"{name} must be a positive integer.")
        if any(not isinstance(item, ApproximateMatch) for item in self.matches):
            raise ConfigurationError("matches must contain ApproximateMatch objects.")
        if any(item.target_index >= self.target_count for item in self.matches):
            raise ConfigurationError("Approximate match target index is outside target_count.")
        if any(item.distance > self.max_distance for item in self.matches):
            raise ConfigurationError("Approximate match exceeds max_distance.")
        if self.dp_cells > self.max_cells or len(self.matches) > self.max_matches:
            raise ConfigurationError("Approximate result exceeds a declared resource limit.")
        if not isinstance(self.reverse_complement, bool):
            raise ConfigurationError("reverse_complement must be boolean.")

    @property
    def found(self) -> bool:
        return bool(self.matches)


@dataclass(frozen=True)
class Mismatch(_SerializableResult):
    """One literal symbol mismatch at a shared sequence position."""

    position: int
    left_symbol: str
    right_symbol: str

    def __post_init__(self) -> None:
        _non_negative_int(self.position, "position")
        if len(self.left_symbol) != 1 or len(self.right_symbol) != 1:
            raise ConfigurationError("Mismatch symbols must each contain one character.")


@dataclass(frozen=True)
class EditStep(_SerializableResult):
    """One deterministic Levenshtein traceback step."""

    operation: EditOperation
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    left_symbol: str | None
    right_symbol: str | None

    def __post_init__(self) -> None:
        if self.operation not in ("match", "substitute", "insert", "delete"):
            raise ConfigurationError("Unknown edit-path operation.")
        for name in ("left_start", "left_end", "right_start", "right_end"):
            _non_negative_int(getattr(self, name), name)
        if self.left_end < self.left_start or self.right_end < self.right_start:
            raise ConfigurationError("EditStep coordinates must be ordered.")
        for name in ("left_symbol", "right_symbol"):
            symbol = getattr(self, name)
            if symbol is not None and (not isinstance(symbol, str) or len(symbol) != 1):
                raise ConfigurationError(f"EditStep {name} must be one character or None.")


@dataclass(frozen=True)
class DistanceResult(_SerializableResult):
    """Hamming or weighted Levenshtein distance with optional details."""

    name: str
    method: DistanceKind
    left_id: str | None
    right_id: str | None
    left_length: int
    right_length: int
    distance: float
    mismatches: tuple[Mismatch, ...]
    edit_path: tuple[EditStep, ...] | None
    costs: FrozenDict
    max_distance: float | None
    exceeded_max_distance: bool
    iupac_matching: Literal["literal"] = "literal"
    max_cells: int | None = None
    dp_cells: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigurationError("DistanceResult name must be non-empty.")
        if self.method not in ("hamming", "levenshtein"):
            raise ConfigurationError("Unknown distance method.")
        _optional_identifier(self.left_id, "left_id")
        _optional_identifier(self.right_id, "right_id")
        _non_negative_int(self.left_length, "left_length")
        _non_negative_int(self.right_length, "right_length")
        if (
            isinstance(self.distance, bool)
            or not isinstance(self.distance, (int, float))
            or not math.isfinite(self.distance)
            or self.distance < 0
        ):
            raise ConfigurationError("distance must be finite and non-negative.")
        if any(not isinstance(item, Mismatch) for item in self.mismatches):
            raise ConfigurationError("mismatches must contain Mismatch objects.")
        if self.edit_path is not None and any(
            not isinstance(item, EditStep) for item in self.edit_path
        ):
            raise ConfigurationError("edit_path must contain EditStep objects or be None.")
        if self.max_distance is not None and (
            isinstance(self.max_distance, bool)
            or not isinstance(self.max_distance, (int, float))
            or not math.isfinite(self.max_distance)
            or self.max_distance < 0
        ):
            raise ConfigurationError("max_distance must be finite and non-negative or None.")
        if not isinstance(self.exceeded_max_distance, bool):
            raise ConfigurationError("exceeded_max_distance must be boolean.")
        if self.max_cells is not None and (
            isinstance(self.max_cells, bool)
            or not isinstance(self.max_cells, int)
            or self.max_cells <= 0
        ):
            raise ConfigurationError("max_cells must be a positive integer or None.")
        if self.dp_cells is not None and (
            isinstance(self.dp_cells, bool)
            or not isinstance(self.dp_cells, int)
            or self.dp_cells <= 0
        ):
            raise ConfigurationError("dp_cells must be a positive integer or None.")
        if (
            self.max_cells is not None
            and self.dp_cells is not None
            and self.dp_cells > self.max_cells
        ):
            raise ConfigurationError("dp_cells cannot exceed max_cells in a completed result.")


@dataclass(frozen=True)
class SimilarityResult(_SerializableResult):
    """Scalar sequence, k-mer, or fingerprint comparison."""

    name: str
    method: str
    value: float
    value_kind: ValueKind
    left_id: str | None
    right_id: str | None
    left_dimension: int
    right_dimension: int
    parameters: FrozenDict
    components: FrozenDict
    zero_vector_policy: str | None
    iupac_matching: Literal["literal"] | None

    def __post_init__(self) -> None:
        for name in ("name", "method"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"SimilarityResult {name} must be non-empty.")
        if self.value_kind not in ("similarity", "distance"):
            raise ConfigurationError("value_kind must be 'similarity' or 'distance'.")
        if (
            isinstance(self.value, bool)
            or not isinstance(self.value, (int, float))
            or not math.isfinite(self.value)
        ):
            raise ConfigurationError("SimilarityResult value must be finite.")
        if self.value_kind == "distance" and self.value < 0:
            raise ConfigurationError("A distance value cannot be negative.")
        if self.value_kind == "similarity" and not -1 <= self.value <= 1:
            raise ConfigurationError("A similarity value must be in [-1, 1].")
        _optional_identifier(self.left_id, "left_id")
        _optional_identifier(self.right_id, "right_id")
        _non_negative_int(self.left_dimension, "left_dimension")
        _non_negative_int(self.right_dimension, "right_dimension")
        if self.zero_vector_policy is not None and (
            not isinstance(self.zero_vector_policy, str) or not self.zero_vector_policy.strip()
        ):
            raise ConfigurationError("zero_vector_policy must be non-empty or None.")
        if self.iupac_matching not in (None, "literal"):
            raise ConfigurationError("Only literal IUPAC matching is implemented in the MVP.")


@dataclass(frozen=True)
class SimilarityMatrixResult(_SerializableResult):
    """Bounded in-memory dense pairwise matrix."""

    name: str
    method: str
    value_kind: ValueKind
    labels: tuple[str, ...]
    values: tuple[tuple[float, ...], ...]
    symmetric: bool
    max_items: int
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigurationError("SimilarityMatrixResult name must be non-empty.")
        if not isinstance(self.method, str) or not self.method.strip():
            raise ConfigurationError("SimilarityMatrixResult method must be non-empty.")
        if self.value_kind not in ("similarity", "distance"):
            raise ConfigurationError("Matrix value_kind is invalid.")
        if any(not isinstance(label, str) or not label for label in self.labels):
            raise ConfigurationError("Matrix labels must be non-empty strings.")
        _non_negative_int(self.max_items, "max_items")
        size = len(self.labels)
        if len(self.values) != size or any(len(row) != size for row in self.values):
            raise ConfigurationError("Matrix values must be square and match labels.")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or (self.value_kind == "distance" and value < 0)
            or (self.value_kind == "similarity" and not -1 <= value <= 1)
            for row in self.values
            for value in row
        ):
            raise ConfigurationError("Matrix contains an invalid value.")
        if self.symmetric and any(
            not math.isclose(self.values[i][j], self.values[j][i])
            for i in range(size)
            for j in range(i)
        ):
            raise ConfigurationError("Matrix declared symmetric but values are asymmetric.")

    @property
    def item_count(self) -> int:
        return len(self.labels)


def _validate_dashing_common(
    algorithm_version: str,
    input_sha256: tuple[str, ...],
    provenance: Provenance,
    raw_output_path: str | None,
    raw_output_sha256: str,
) -> None:
    if not isinstance(algorithm_version, str) or not algorithm_version.strip():
        raise ConfigurationError("Dashing algorithm_version must be non-empty.")
    if not input_sha256 or any(not _SHA256.fullmatch(value) for value in input_sha256):
        raise ConfigurationError("Dashing input_sha256 must contain lowercase SHA-256 digests.")
    if not isinstance(provenance, Provenance):
        raise ConfigurationError("Dashing provenance must be Provenance.")
    _optional_identifier(raw_output_path, "raw_output_path")
    if not _SHA256.fullmatch(raw_output_sha256):
        raise ConfigurationError("Dashing raw_output_sha256 must be a lowercase SHA-256 digest.")


@dataclass(frozen=True)
class DashingJaccardMatrixResult(_SerializableResult):
    """Validated symmetric Jaccard matrix emitted by a user-supplied Dashing CLI."""

    name: str
    method: Literal["exact-kmer-set", "hll-sketch"]
    algorithm_version: str
    labels: tuple[str, ...]
    values: tuple[tuple[float, ...], ...]
    input_sha256: tuple[str, ...]
    parameters: FrozenDict
    provenance: Provenance
    raw_output_path: str | None
    raw_output_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigurationError("Dashing matrix name must be non-empty.")
        if self.method not in {"exact-kmer-set", "hll-sketch"}:
            raise ConfigurationError("Dashing matrix method is invalid.")
        if len(self.labels) < 2 or any(
            not isinstance(label, str) or not label for label in self.labels
        ):
            raise ConfigurationError("Dashing matrix requires at least two non-empty labels.")
        size = len(self.labels)
        if len(self.values) != size or any(len(row) != size for row in self.values):
            raise ConfigurationError("Dashing matrix values must be square and match labels.")
        if len(self.input_sha256) != size:
            raise ConfigurationError("Dashing matrix checksums must match labels.")
        if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0 <= value <= 1
            for row in self.values
            for value in row
        ):
            raise ConfigurationError("Dashing matrix contains an invalid Jaccard value.")
        if any(not math.isclose(self.values[index][index], 1.0) for index in range(size)):
            raise ConfigurationError("Dashing Jaccard matrix diagonal must equal one.")
        if any(
            not math.isclose(self.values[left][right], self.values[right][left])
            for left in range(size)
            for right in range(left)
        ):
            raise ConfigurationError("Dashing Jaccard matrix must be symmetric.")
        _validate_dashing_common(
            self.algorithm_version,
            self.input_sha256,
            self.provenance,
            self.raw_output_path,
            self.raw_output_sha256,
        )

    @property
    def item_count(self) -> int:
        return len(self.labels)


@dataclass(frozen=True)
class DashingNeighborHit(_SerializableResult):
    """One deterministic neighbor selected from a validated Dashing matrix."""

    label: str
    index: int
    jaccard: float

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label:
            raise ConfigurationError("Dashing neighbor label must be non-empty.")
        _non_negative_int(self.index, "index")
        if (
            isinstance(self.jaccard, bool)
            or not isinstance(self.jaccard, (int, float))
            or not math.isfinite(self.jaccard)
            or not 0 <= self.jaccard <= 1
        ):
            raise ConfigurationError("Dashing neighbor Jaccard must be finite in [0, 1].")


@dataclass(frozen=True)
class DashingNeighborRow(_SerializableResult):
    """Top-k neighbors for one matrix row."""

    query_label: str
    query_index: int
    hits: tuple[DashingNeighborHit, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.query_label, str) or not self.query_label:
            raise ConfigurationError("Dashing query label must be non-empty.")
        _non_negative_int(self.query_index, "query_index")
        if any(not isinstance(hit, DashingNeighborHit) for hit in self.hits):
            raise ConfigurationError("Dashing hits must contain DashingNeighborHit objects.")
        if any(hit.index == self.query_index for hit in self.hits):
            raise ConfigurationError("Dashing Top-k rows cannot contain the query itself.")
        expected = tuple(sorted(self.hits, key=lambda hit: (-hit.jaccard, hit.index)))
        if self.hits != expected:
            raise ConfigurationError("Dashing Top-k hits must use deterministic score/index order.")


@dataclass(frozen=True)
class DashingTopKResult(_SerializableResult):
    """Per-item Top-k neighbors derived from a validated Dashing Jaccard matrix."""

    name: str
    method: Literal["exact-kmer-set", "hll-sketch"]
    algorithm_version: str
    labels: tuple[str, ...]
    rows: tuple[DashingNeighborRow, ...]
    top_k: int
    input_sha256: tuple[str, ...]
    parameters: FrozenDict
    provenance: Provenance
    raw_output_path: str | None
    raw_output_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ConfigurationError("Dashing Top-k name must be non-empty.")
        if self.method not in {"exact-kmer-set", "hll-sketch"}:
            raise ConfigurationError("Dashing Top-k method is invalid.")
        size = len(self.labels)
        if size < 2 or any(not isinstance(label, str) or not label for label in self.labels):
            raise ConfigurationError("Dashing Top-k requires at least two non-empty labels.")
        if (
            isinstance(self.top_k, bool)
            or not isinstance(self.top_k, int)
            or not 1 <= self.top_k < size
        ):
            raise ConfigurationError("Dashing top_k must be in [1, item_count - 1].")
        if len(self.rows) != size or any(
            not isinstance(row, DashingNeighborRow) for row in self.rows
        ):
            raise ConfigurationError("Dashing Top-k rows must match labels.")
        if tuple(row.query_index for row in self.rows) != tuple(range(size)):
            raise ConfigurationError("Dashing Top-k rows must cover query indexes in order.")
        if any(row.query_label != self.labels[row.query_index] for row in self.rows):
            raise ConfigurationError("Dashing Top-k query labels do not match labels.")
        if any(len(row.hits) != self.top_k for row in self.rows):
            raise ConfigurationError("Every Dashing Top-k row must contain top_k hits.")
        if any(
            hit.index >= size or hit.label != self.labels[hit.index]
            for row in self.rows
            for hit in row.hits
        ):
            raise ConfigurationError("A Dashing Top-k hit is outside the result labels.")
        if len(self.input_sha256) != size:
            raise ConfigurationError("Dashing Top-k checksums must match labels.")
        _validate_dashing_common(
            self.algorithm_version,
            self.input_sha256,
            self.provenance,
            self.raw_output_path,
            self.raw_output_sha256,
        )


@dataclass(frozen=True)
class SketchSimilarityResult(_SerializableResult):
    """Set-based comparison of two compatible sketches."""

    name: str
    method: str
    value: float
    metric: Literal["jaccard", "containment"]
    left_id: str | None
    right_id: str | None
    left_hash_count: int
    right_hash_count: int
    shared_hash_count: int
    union_hash_count: int
    denominator: int
    passed_min_shared_hashes: bool
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if not self.name or not self.method:
            raise ConfigurationError("Sketch similarity name and method must be non-empty.")
        if self.metric not in {"jaccard", "containment"}:
            raise ConfigurationError("Sketch similarity metric is invalid.")
        if (
            isinstance(self.value, bool)
            or not isinstance(self.value, (int, float))
            or not math.isfinite(self.value)
            or not 0 <= self.value <= 1
        ):
            raise ConfigurationError("Sketch similarity value must be finite in [0, 1].")
        _optional_identifier(self.left_id, "left_id")
        _optional_identifier(self.right_id, "right_id")
        for name in (
            "left_hash_count",
            "right_hash_count",
            "shared_hash_count",
            "union_hash_count",
            "denominator",
        ):
            _non_negative_int(getattr(self, name), name)
        if self.shared_hash_count > min(self.left_hash_count, self.right_hash_count):
            raise ConfigurationError("Shared sketch hashes exceed an input hash count.")
        if self.union_hash_count < max(self.left_hash_count, self.right_hash_count):
            raise ConfigurationError("Sketch union is smaller than an input hash count.")
        if not isinstance(self.passed_min_shared_hashes, bool):
            raise ConfigurationError("passed_min_shared_hashes must be boolean.")


@dataclass(frozen=True)
class NearestNeighborHit(_SerializableResult):
    record_id: str
    index: int
    similarity: float
    shared_hash_count: int

    def __post_init__(self) -> None:
        if not isinstance(self.record_id, str) or not self.record_id.strip():
            raise ConfigurationError("Nearest-neighbor ID must be non-empty.")
        _non_negative_int(self.index, "index")
        _non_negative_int(self.shared_hash_count, "shared_hash_count")
        if (
            isinstance(self.similarity, bool)
            or not isinstance(self.similarity, (int, float))
            or not math.isfinite(self.similarity)
            or not 0 <= self.similarity <= 1
        ):
            raise ConfigurationError("Nearest-neighbor similarity must be in [0, 1].")


@dataclass(frozen=True)
class NearestNeighborResult(_SerializableResult):
    name: str
    method: str
    algorithm_version: str
    query_id: str | None
    index_size: int
    hits: tuple[NearestNeighborHit, ...]
    top_k: int
    min_similarity: float
    parameters: FrozenDict

    def __post_init__(self) -> None:
        for name in ("name", "method", "algorithm_version"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise ConfigurationError(f"NearestNeighborResult {name} must be non-empty.")
        _optional_identifier(self.query_id, "query_id")
        _non_negative_int(self.index_size, "index_size")
        if isinstance(self.top_k, bool) or not isinstance(self.top_k, int) or self.top_k <= 0:
            raise ConfigurationError("top_k must be positive.")
        if len(self.hits) > self.top_k or any(
            not isinstance(item, NearestNeighborHit) for item in self.hits
        ):
            raise ConfigurationError("Nearest-neighbor hits violate the declared top_k.")
        if any(item.index >= self.index_size for item in self.hits):
            raise ConfigurationError("Nearest-neighbor hit index is outside the index.")
        if (
            isinstance(self.min_similarity, bool)
            or not isinstance(self.min_similarity, (int, float))
            or not 0 <= self.min_similarity <= 1
        ):
            raise ConfigurationError("min_similarity must be in [0, 1].")


__all__ = [
    "ApproximateMatch",
    "ApproximateSearchResult",
    "DashingJaccardMatrixResult",
    "DashingNeighborHit",
    "DashingNeighborRow",
    "DashingTopKResult",
    "DistanceResult",
    "EditStep",
    "Mismatch",
    "NearestNeighborHit",
    "NearestNeighborResult",
    "SearchResult",
    "SequenceMatch",
    "SimilarityMatrixResult",
    "SimilarityResult",
    "SketchSimilarityResult",
]
