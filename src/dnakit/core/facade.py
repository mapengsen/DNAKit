"""One user-facing entry point for one or many DNA records."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, replace
from itertools import islice
from typing import TYPE_CHECKING, TypeAlias, cast

from dnakit.core._json import FrozenDict
from dnakit.core.collection import DNASet
from dnakit.core.coordinates import CompoundLocation, Interval, Location, UnresolvedLocation
from dnakit.core.enums import DNAAlphabet, Strandedness, Topology
from dnakit.core.feature import DNAFeature
from dnakit.core.gap import Gap
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import (
    ConfigurationError,
    FeatureError,
    InvalidAlphabetError,
    SequenceError,
)

if TYPE_CHECKING:
    from dnakit.standardize.config import NormalizationConfig
    from dnakit.standardize.results import NormalizationResult


FeatureInput: TypeAlias = DNAFeature | Mapping[str, object]
_MAX_DEFAULT_RECORDS = 100_000
_RECORD_SPEC_KEYS = frozenset(
    {
        "sequence",
        "parts",
        "id",
        "description",
        "features",
        "metadata",
        "letter_annotations",
        "alphabet",
        "topology",
        "strandedness",
        "normalization_config",
    }
)
_FEATURE_SPEC_KEYS = frozenset(
    {
        "type",
        "location",
        "start",
        "end",
        "id",
        "strand",
        "label",
        "score",
        "phase",
        "qualifiers",
        "source",
    }
)


@dataclass(frozen=True, slots=True)
class _RecordOptions:
    id: str | None = None
    description: str | None = None
    features: tuple[DNAFeature, ...] | None = None
    metadata: Mapping[str, object] | None = None
    letter_annotations: Mapping[str, tuple[int | float, ...]] | None = None
    alphabet: DNAAlphabet | str | None = None
    topology: Topology | str | None = None
    strandedness: Strandedness | str | None = None
    normalization_config: NormalizationConfig | None = None


@dataclass(frozen=True, init=False)
class DNA:
    """A stable facade over one or many DNA records.

    Raw text is normalized with a retained audit. A mapping can describe one
    record, while an iterable can contain raw sequences, records, or record
    mappings. Existing ``DNASequence``, ``DNARecord`` and ``DNASet`` values are
    accepted so advanced and legacy code can enter the same public path.
    """

    _dataset: DNASet
    normalizations: tuple[NormalizationResult | None, ...]

    __hash__ = None  # type: ignore[assignment]

    def __init__(
        self,
        data: object,
        *,
        id: str | None = None,
        description: str | None = None,
        features: FeatureInput | Iterable[FeatureInput] | None = None,
        metadata: Mapping[str, object] | None = None,
        letter_annotations: Mapping[str, Iterable[int | float]] | None = None,
        alphabet: DNAAlphabet | str | None = None,
        topology: Topology | str | None = None,
        strandedness: Strandedness | str | None = None,
        normalization_config: NormalizationConfig | None = None,
        name: str | None = None,
        source: str | None = None,
        version: str | None = None,
        collection_metadata: Mapping[str, object] | None = None,
        max_records: int = _MAX_DEFAULT_RECORDS,
    ) -> None:
        """Build one consistent object from one sequence or a bounded collection."""

        _validate_max_records(max_records)
        defaults = _RecordOptions(
            id=id,
            description=description,
            features=None if features is None else _coerce_features(features),
            metadata=metadata,
            letter_annotations=_materialize_annotations(letter_annotations),
            alphabet=alphabet,
            topology=topology,
            strandedness=strandedness,
            normalization_config=normalization_config,
        )
        items, inherited_normalizations, inherited_dataset = _materialize_data(
            data, max_records=max_records
        )
        if len(items) > 1 and defaults.id is not None:
            raise ConfigurationError(
                "A shared id cannot identify multiple DNA records.",
                code="MULTIPLE_DNA_SHARED_ID",
                hint="Put a distinct 'id' in each record mapping or omit it for generated IDs.",
            )

        records: list[DNARecord] = []
        audits: list[NormalizationResult | None] = []
        for index, item in enumerate(items):
            inherited_audit = (
                inherited_normalizations[index] if index < len(inherited_normalizations) else None
            )
            record, audit = _build_record(
                item,
                index=index,
                defaults=defaults,
                inherited_audit=inherited_audit,
            )
            records.append(record)
            audits.append(audit)

        resolved_name = name if name is not None else getattr(inherited_dataset, "name", None)
        resolved_source = (
            source if source is not None else getattr(inherited_dataset, "source", None)
        )
        resolved_version = (
            version if version is not None else getattr(inherited_dataset, "version", None)
        )
        resolved_collection_metadata = (
            collection_metadata
            if collection_metadata is not None
            else getattr(inherited_dataset, "metadata", None)
        )
        dataset = DNASet(
            records,
            name=resolved_name,
            source=resolved_source,
            version=resolved_version,
            metadata=resolved_collection_metadata,
        )
        object.__setattr__(self, "_dataset", dataset)
        object.__setattr__(self, "normalizations", tuple(audits))

    def __len__(self) -> int:
        """Return the number of contained DNA records."""

        return len(self._dataset)

    def __iter__(self) -> Iterator[DNARecord]:
        """Iterate over compatible internal records for dataset algorithms."""

        return iter(self._dataset)

    def __getitem__(self, index: int | slice) -> DNA:
        """Select one or more records while staying in the public DNA type."""

        indices = range(len(self))[index]
        if isinstance(indices, int):
            return self._derive(((indices, self.records[indices]),))
        return self._derive((position, self.records[position]) for position in indices)

    @property
    def records(self) -> tuple[DNARecord, ...]:
        """The immutable compatible records used by advanced APIs."""

        return self._dataset.records

    @property
    def record_count(self) -> int:
        """Number of contained records."""

        return len(self)

    @property
    def is_single(self) -> bool:
        """Whether this object contains exactly one DNA record."""

        return len(self) == 1

    @property
    def record(self) -> DNARecord:
        """Return the only record, rejecting empty or multi-record data."""

        if not self.is_single:
            raise SequenceError(
                "This operation requires exactly one DNA record.",
                code="SINGLE_DNA_REQUIRED",
                context={"record_count": len(self)},
                hint="Select one record or use a dataset-level operation.",
            )
        return self.records[0]

    @property
    def sequence(self) -> DNASequence:
        """Return the sequence of a single-record object."""

        return self.record.sequence

    @property
    def id(self) -> str:
        """Return the identifier of a single-record object."""

        return self.record.id

    @property
    def ids(self) -> tuple[str, ...]:
        """Return every record identifier in stable input order."""

        return self._dataset.ids

    @property
    def description(self) -> str:
        """Return the description of a single-record object."""

        return self.record.description

    @property
    def features(self) -> tuple[DNAFeature, ...]:
        """Return the features of a single-record object."""

        return self.record.features

    @property
    def metadata(self) -> Mapping[str, object]:
        """Return record metadata for one record, otherwise collection metadata."""

        return self.record.metadata if self.is_single else self._dataset.metadata

    @property
    def collection_metadata(self) -> FrozenDict:
        """Return metadata attached to the collection itself."""

        return self._dataset.metadata

    @property
    def name(self) -> str | None:
        return self._dataset.name

    @property
    def source(self) -> str | None:
        return self._dataset.source

    @property
    def version(self) -> str | None:
        return self._dataset.version

    @property
    def normalization(self) -> NormalizationResult | None:
        """Return the retained normalization audit for one record."""

        _ = self.record
        return self.normalizations[0]

    @property
    def symbols(self) -> str:
        """Convenience access to normalized symbols for one record."""

        return self.sequence.symbols

    @property
    def symbol_length(self) -> int:
        return self.sequence.symbol_length

    @property
    def parts(self) -> tuple[str | Gap, ...]:
        return self.sequence.parts

    @property
    def is_gapped(self) -> bool:
        return self.sequence.is_gapped

    @property
    def has_unknown_length(self) -> bool:
        return self.sequence.has_unknown_length

    @property
    def canonical_base_count(self) -> int:
        return self.sequence.canonical_base_count

    @property
    def ambiguity_count(self) -> int:
        return self.sequence.ambiguity_count

    @property
    def coordinate_span(self) -> int | None:
        return self.sequence.coordinate_span

    @property
    def alphabet(self) -> DNAAlphabet:
        return self.sequence.alphabet

    @property
    def topology(self) -> Topology:
        return self.sequence.topology

    @property
    def strandedness(self) -> Strandedness:
        return self.sequence.strandedness

    @property
    def dataset(self) -> DNASet:
        """Return the compatible immutable dataset for advanced APIs."""

        return self._dataset

    def _derive(self, indexed_records: Iterable[tuple[int, DNARecord]]) -> DNA:
        """Build a facade result while retaining audits for unchanged source sequences."""

        selected = tuple(indexed_records)
        indices = tuple(index for index, _ in selected)
        if any(
            isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(self)
            for index in indices
        ):
            raise IndexError("A derived DNA record index is outside the source object.")
        if len(indices) != len(set(indices)):
            raise ValueError("A derived DNA object cannot select one source record twice.")
        dataset = DNASet(
            (record for _, record in selected),
            name=self.name,
            source=self.source,
            version=self.version,
            metadata=self.collection_metadata,
        )
        derived = object.__new__(DNA)
        object.__setattr__(derived, "_dataset", dataset)
        object.__setattr__(
            derived, "normalizations", tuple(self.normalizations[i] for i in indices)
        )
        return derived


SingleDNAInput: TypeAlias = DNA | DNASequence | DNARecord


def resolve_single_dna(value: object) -> tuple[DNASequence, DNARecord | None]:
    """Resolve one facade/record/sequence without exposing constructor choices."""

    if isinstance(value, DNA):
        record = value.record
        return record.sequence, record
    if isinstance(value, DNARecord):
        return value.sequence, value
    if isinstance(value, DNASequence):
        return value, None
    raise TypeError("value must be DNA, DNARecord, or DNASequence.")


def _validate_max_records(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            "max_records must be a positive integer.",
            code="INVALID_DNA_RECORD_LIMIT",
            context={"max_records": value},
        )


def _materialize_data(
    data: object,
    *,
    max_records: int,
) -> tuple[tuple[object, ...], tuple[NormalizationResult | None, ...], DNASet | None]:
    if isinstance(data, DNA):
        return data.records, data.normalizations, data.dataset
    if isinstance(data, DNASet):
        return data.records, (None,) * len(data), data
    if isinstance(data, (str, bytes, DNASequence, DNARecord, Mapping)):
        return (data,), (), None
    if not isinstance(data, Iterable):
        raise TypeError(
            "DNA data must be sequence text, a record mapping, a DNA object, or an iterable."
        )

    items = tuple(islice(data, max_records + 1))
    if len(items) > max_records:
        raise ConfigurationError(
            "DNA input exceeds max_records.",
            code="DNA_RECORD_LIMIT_EXCEEDED",
            context={"max_records": max_records},
        )
    if any(isinstance(item, Gap) for item in items):
        if any(not isinstance(item, (str, Gap)) for item in items):
            raise TypeError("A gapped sequence may contain only text fragments and Gap objects.")
        return (items,), (), None
    return cast(tuple[object, ...], items), (), None


def _build_record(
    item: object,
    *,
    index: int,
    defaults: _RecordOptions,
    inherited_audit: NormalizationResult | None,
) -> tuple[DNARecord, NormalizationResult | None]:
    raw, options = _record_spec(item, defaults=defaults)
    source_record = raw if isinstance(raw, DNARecord) else None
    sequence_input = source_record.sequence if source_record is not None else raw

    sequence_changed = any(
        value is not None
        for value in (
            options.alphabet,
            options.topology,
            options.strandedness,
            options.normalization_config,
        )
    )
    if source_record is not None and not sequence_changed:
        sequence = source_record.sequence
        audit = inherited_audit
    else:
        sequence, audit = _normalize_sequence(
            sequence_input,
            alphabet=options.alphabet,
            topology=options.topology,
            strandedness=options.strandedness,
            config=options.normalization_config,
        )

    record_id = options.id or (source_record.id if source_record is not None else None)
    if record_id is None:
        record_id = f"sequence_{index + 1}"
    description = (
        options.description
        if options.description is not None
        else source_record.description
        if source_record is not None
        else ""
    )
    features = (
        options.features
        if options.features is not None
        else source_record.features
        if source_record is not None
        else ()
    )
    metadata = (
        options.metadata
        if options.metadata is not None
        else source_record.metadata
        if source_record is not None
        else None
    )
    annotations = (
        options.letter_annotations
        if options.letter_annotations is not None
        else source_record.letter_annotations
        if source_record is not None
        else None
    )
    return (
        DNARecord(
            sequence,
            record_id,
            description=description,
            features=features,
            metadata=metadata,
            letter_annotations=annotations,
        ),
        audit,
    )


def _record_spec(item: object, *, defaults: _RecordOptions) -> tuple[object, _RecordOptions]:
    from dnakit.standardize.config import NormalizationConfig

    if not isinstance(item, Mapping):
        return item, defaults
    unknown = set(item) - _RECORD_SPEC_KEYS
    if unknown:
        raise ConfigurationError(
            "A DNA record mapping contains unknown fields.",
            code="UNKNOWN_DNA_RECORD_FIELD",
            context={"fields": sorted(str(key) for key in unknown)},
        )
    has_sequence = "sequence" in item
    has_parts = "parts" in item
    if has_sequence == has_parts:
        raise ConfigurationError(
            "A DNA record mapping must contain exactly one of 'sequence' or 'parts'.",
            code="DNA_RECORD_SEQUENCE_REQUIRED",
        )
    raw = item["sequence"] if has_sequence else item["parts"]

    local_features = (
        _coerce_features(cast(object, item["features"]))
        if "features" in item
        else defaults.features
    )
    local_metadata = _optional_mapping(item.get("metadata"), name="metadata")
    if "metadata" not in item:
        local_metadata = defaults.metadata
    elif defaults.metadata is not None and local_metadata is not None:
        local_metadata = {**defaults.metadata, **local_metadata}
    local_annotations = (
        _materialize_annotations(
            cast(Mapping[str, Iterable[int | float]] | None, item.get("letter_annotations"))
        )
        if "letter_annotations" in item
        else defaults.letter_annotations
    )
    return raw, _RecordOptions(
        id=cast(str | None, item.get("id", defaults.id)),
        description=cast(str | None, item.get("description", defaults.description)),
        features=local_features,
        metadata=local_metadata,
        letter_annotations=local_annotations,
        alphabet=cast(DNAAlphabet | str | None, item.get("alphabet", defaults.alphabet)),
        topology=cast(Topology | str | None, item.get("topology", defaults.topology)),
        strandedness=cast(
            Strandedness | str | None,
            item.get("strandedness", defaults.strandedness),
        ),
        normalization_config=cast(
            NormalizationConfig | None,
            item.get("normalization_config", defaults.normalization_config),
        ),
    )


def _normalize_sequence(
    raw: object,
    *,
    alphabet: DNAAlphabet | str | None,
    topology: Topology | str | None,
    strandedness: Strandedness | str | None,
    config: NormalizationConfig | None,
) -> tuple[DNASequence, NormalizationResult]:
    from dnakit.standardize.config import NormalizationConfig
    from dnakit.standardize.normalize import normalize

    if not isinstance(raw, (str, bytes, DNASequence, Iterable)) or isinstance(raw, Mapping):
        raise TypeError("Each DNA item must be text, DNASequence, DNARecord, or string/Gap parts.")
    if config is not None and not isinstance(config, NormalizationConfig):
        raise ConfigurationError(
            "normalization_config must be NormalizationConfig or None.",
            code="INVALID_DNA_NORMALIZATION_CONFIG",
        )
    resolved_alphabet = _coerce_alphabet(alphabet)
    if config is not None and resolved_alphabet is not None:
        if config.alphabet is not None and config.alphabet is not resolved_alphabet:
            raise ConfigurationError(
                "alphabet conflicts with normalization_config.alphabet.",
                code="DNA_ALPHABET_CONFIG_CONFLICT",
            )
        config = replace(config, alphabet=resolved_alphabet)
    elif config is None and resolved_alphabet is not None:
        config = NormalizationConfig(alphabet=resolved_alphabet)

    normalized = normalize(
        cast(str | bytes | Iterable[str | Gap] | DNASequence, raw), config=config
    )
    if normalized.sequence is None:
        raise InvalidAlphabetError(
            "DNA input could not be converted into a valid sequence.",
            code="DNA_NORMALIZATION_FAILED",
            context={"issue_codes": [issue.code for issue in normalized.issues]},
            hint="Inspect the retained normalization policy or correct the input symbols.",
        )
    source_sequence = raw if isinstance(raw, DNASequence) else None
    resolved_topology = (
        topology
        if topology is not None
        else source_sequence.topology
        if source_sequence is not None
        else Topology.LINEAR
    )
    resolved_strandedness = (
        strandedness
        if strandedness is not None
        else source_sequence.strandedness
        if source_sequence is not None
        else Strandedness.SINGLE
    )
    sequence = DNASequence(
        normalized.sequence.parts,
        alphabet=normalized.sequence.alphabet,
        topology=resolved_topology,
        strandedness=resolved_strandedness,
    )
    return sequence, replace(normalized, sequence=sequence)


def _coerce_alphabet(value: DNAAlphabet | str | None) -> DNAAlphabet | None:
    if value is None:
        return None
    try:
        return value if isinstance(value, DNAAlphabet) else DNAAlphabet(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Unknown DNA alphabet.",
            code="INVALID_DNA_ALPHABET",
            context={"alphabet": value},
            hint="Choose 'strict' or 'iupac'.",
        ) from exc


def _coerce_features(value: object) -> tuple[DNAFeature, ...]:
    if value is None:
        return ()
    if isinstance(value, (DNAFeature, Mapping)):
        items = (value,)
    elif isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        items = tuple(value)
    else:
        raise FeatureError(
            "features must be a feature, a feature mapping, or an iterable of them.",
            code="INVALID_DNA_FEATURES",
        )
    return tuple(_coerce_feature(item) for item in items)


def _coerce_feature(value: object) -> DNAFeature:
    if isinstance(value, DNAFeature):
        return value
    if not isinstance(value, Mapping):
        raise FeatureError(
            "Every feature must be DNAFeature or a mapping.",
            code="INVALID_DNA_FEATURE",
            context={"type": type(value).__name__},
        )
    unknown = set(value) - _FEATURE_SPEC_KEYS
    if unknown:
        raise FeatureError(
            "A feature mapping contains unknown fields.",
            code="UNKNOWN_DNA_FEATURE_FIELD",
            context={"fields": sorted(str(key) for key in unknown)},
        )
    if "type" not in value:
        raise FeatureError("A feature mapping requires 'type'.", code="DNA_FEATURE_TYPE_REQUIRED")
    has_location = "location" in value
    has_interval = "start" in value or "end" in value
    if has_location and has_interval:
        raise FeatureError(
            "Use either location or start/end for a feature, not both.",
            code="DNA_FEATURE_LOCATION_CONFLICT",
        )
    if has_location:
        location = value["location"]
        if not isinstance(location, (Interval, CompoundLocation, UnresolvedLocation)):
            raise FeatureError(
                "feature location must be Interval, CompoundLocation, or UnresolvedLocation.",
                code="INVALID_DNA_FEATURE_LOCATION",
            )
    else:
        if "start" not in value or "end" not in value:
            raise FeatureError(
                "A simple feature mapping requires both start and end.",
                code="DNA_FEATURE_INTERVAL_REQUIRED",
            )
        location = Interval(cast(int, value["start"]), cast(int, value["end"]))
    qualifiers = _optional_mapping(value.get("qualifiers"), name="qualifiers")
    return DNAFeature(
        cast(str, value["type"]),
        cast(Location, location),
        id=cast(str | None, value.get("id")),
        strand=cast(str, value.get("strand", "unknown")),
        label=cast(str | None, value.get("label")),
        score=cast(float | None, value.get("score")),
        phase=cast(int | None, value.get("phase")),
        qualifiers=qualifiers,
        source=cast(str | None, value.get("source")),
    )


def _optional_mapping(value: object, *, name: str) -> Mapping[str, object] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            f"{name} must be a mapping or None.",
            code="INVALID_DNA_MAPPING",
            context={"field": name, "type": type(value).__name__},
        )
    return cast(Mapping[str, object], value)


def _materialize_annotations(
    value: Mapping[str, Iterable[int | float]] | None,
) -> Mapping[str, tuple[int | float, ...]] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            "letter_annotations must be a mapping or None.",
            code="INVALID_DNA_LETTER_ANNOTATIONS",
        )
    return {name: tuple(values) for name, values in value.items()}


__all__ = ["DNA", "FeatureInput", "SingleDNAInput", "resolve_single_dna"]
