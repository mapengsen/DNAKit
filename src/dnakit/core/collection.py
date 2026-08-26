"""In-memory immutable DNA record collection."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import overload

from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import DuplicateIDError, SequenceError

IDFactory = Callable[[int, DNASequence], str]


def _default_id_factory(index: int, sequence: DNASequence) -> str:
    del sequence
    return f"sequence_{index + 1}"


@dataclass(frozen=True, init=False)
class DNASet:
    """A repeatable, materialized collection of :class:`DNARecord` objects."""

    records: tuple[DNARecord, ...]
    name: str | None
    source: str | None
    version: str | None
    metadata: FrozenDict

    __hash__ = None  # type: ignore[assignment]

    def __init__(
        self,
        records: Iterable[DNARecord],
        *,
        name: str | None = None,
        source: str | None = None,
        version: str | None = None,
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        record_tuple = tuple(records)
        if any(not isinstance(record, DNARecord) for record in record_tuple):
            raise SequenceError(
                "DNASet records must all be DNARecord objects.",
                code="INVALID_DNASET_RECORD",
            )
        for field_name, value in (("name", name), ("source", source), ("version", version)):
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise SequenceError(
                    f"DNASet {field_name} must be None or a non-empty string.",
                    code="INVALID_DNASET_METADATA",
                    context={field_name: value},
                )
        object.__setattr__(self, "records", record_tuple)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))

    @classmethod
    def from_records(
        cls,
        records: Iterable[DNARecord],
        *,
        name: str | None = None,
        source: str | None = None,
        version: str | None = None,
        metadata: Mapping[str, object] | None = None,
    ) -> DNASet:
        return cls(
            records,
            name=name,
            source=source,
            version=version,
            metadata=metadata,
        )

    @classmethod
    def from_sequences(
        cls,
        sequences: Iterable[DNASequence],
        *,
        id_factory: IDFactory | None = None,
        name: str | None = None,
        source: str | None = None,
        version: str | None = None,
        metadata: Mapping[str, object] | None = None,
    ) -> DNASet:
        factory = _default_id_factory if id_factory is None else id_factory
        records: list[DNARecord] = []
        for index, sequence in enumerate(sequences):
            if not isinstance(sequence, DNASequence):
                raise SequenceError(
                    "DNASet.from_sequences() accepts only DNASequence objects.",
                    code="INVALID_DNASET_SEQUENCE",
                    context={"index": index, "value_type": type(sequence).__name__},
                )
            record_id = factory(index, sequence)
            records.append(DNARecord(sequence, record_id))
        return cls(
            records,
            name=name,
            source=source,
            version=version,
            metadata=metadata,
        )

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Iterator[DNARecord]:
        return iter(self.records)

    @overload
    def __getitem__(self, index: int) -> DNARecord: ...

    @overload
    def __getitem__(self, index: slice) -> DNASet: ...

    def __getitem__(self, index: int | slice) -> DNARecord | DNASet:
        if isinstance(index, slice):
            return DNASet(
                self.records[index],
                name=self.name,
                source=self.source,
                version=self.version,
                metadata=self.metadata,
            )
        return self.records[index]

    @property
    def ids(self) -> tuple[str, ...]:
        return tuple(record.id for record in self.records)

    def get(self, record_id: str) -> DNARecord:
        """Return one record by ID, rejecting missing and duplicate matches."""

        matches = tuple(record for record in self.records if record.id == record_id)
        if not matches:
            raise KeyError(record_id)
        if len(matches) > 1:
            raise DuplicateIDError(
                "Record ID is not unique in this DNASet.",
                context={"record_id": record_id, "count": len(matches)},
                hint="Run validate() and resolve duplicate IDs before keyed lookup.",
            )
        return matches[0]

    def select(self, indices: Iterable[int]) -> DNASet:
        """Return records at the requested input-order indices."""

        return DNASet(
            (self.records[index] for index in indices),
            name=self.name,
            source=self.source,
            version=self.version,
            metadata=self.metadata,
        )


__all__ = ["DNASet", "IDFactory"]
