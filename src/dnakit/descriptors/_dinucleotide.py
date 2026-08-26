"""Validated user-supplied dinucleotide-property tables.

DNAKit deliberately ships property names and units only. Numerical
coefficients must be supplied by the caller in an explicitly loaded JSON file.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

DINUCLEOTIDES = (
    "AA",
    "AC",
    "AG",
    "AT",
    "CA",
    "CC",
    "CG",
    "CT",
    "GA",
    "GC",
    "GG",
    "GT",
    "TA",
    "TC",
    "TG",
    "TT",
)
DINUCLEOTIDE_TABLE_SCHEMA_VERSION = "dnakit.dinucleotide_property_table.v1"
DEFAULT_MAX_DINUCLEOTIDE_TABLE_BYTES = 1_000_000
_MAX_METADATA_LENGTH = 1_000
_MAX_ABSOLUTE_VALUE = 1.0e12


@dataclass(frozen=True, slots=True)
class DinucleotidePropertySpec:
    """Stable output identity for one optional dinucleotide-property scale."""

    key: str
    display_name: str
    unit: str


DINUCLEOTIDE_PROPERTY_SPECS = (
    DinucleotidePropertySpec("twist", "Twist", "degree"),
    DinucleotidePropertySpec("tilt", "Tilt", "degree"),
    DinucleotidePropertySpec("roll", "Roll", "degree"),
    DinucleotidePropertySpec("shift", "Shift", "angstrom"),
    DinucleotidePropertySpec("slide", "Slide", "angstrom"),
    DinucleotidePropertySpec("rise", "Rise", "angstrom"),
    DinucleotidePropertySpec("bend", "Bend", "degree"),
    DinucleotidePropertySpec("inclination", "Inclination", "degree"),
    DinucleotidePropertySpec("direction", "Direction", "degree"),
    DinucleotidePropertySpec("propeller_twist", "Propeller twist", "degree"),
    DinucleotidePropertySpec("major_groove_width", "Major groove width", "angstrom"),
    DinucleotidePropertySpec("minor_groove_width", "Minor groove width", "angstrom"),
    DinucleotidePropertySpec("persistence_length", "Persistence length", "nanometer"),
    DinucleotidePropertySpec("stacking_energy", "Stacking energy", "kcal/mol"),
    DinucleotidePropertySpec("free_energy", "Free energy", "kcal/mol"),
)


@dataclass(frozen=True, slots=True)
class DinucleotideProperty:
    """One fully validated 16-value property supplied by the caller."""

    key: str
    display_name: str
    unit: str
    values: FrozenDict

    def __post_init__(self) -> None:
        spec = next((item for item in DINUCLEOTIDE_PROPERTY_SPECS if item.key == self.key), None)
        if spec is None or self.display_name != spec.display_name or self.unit != spec.unit:
            raise ConfigurationError(
                "Dinucleotide property identity does not match the fixed schema.",
                code="INVALID_DINUCLEOTIDE_TABLE_PROPERTY",
                context={"property": self.key},
            )
        object.__setattr__(self, "values", _validated_coefficients(self.values, spec))


@dataclass(frozen=True, slots=True)
class DinucleotidePropertyTable:
    """An immutable, checksummed table loaded from a user-owned JSON file."""

    name: str
    version: str
    source: str
    sha256: str
    properties: tuple[DinucleotideProperty, ...]

    def __post_init__(self) -> None:
        for field_name in ("name", "version", "source"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip() or len(value) > _MAX_METADATA_LENGTH:
                raise ConfigurationError(
                    f"Dinucleotide table {field_name} must be non-empty bounded text.",
                    code="INVALID_DINUCLEOTIDE_TABLE_METADATA",
                )
        if (
            not isinstance(self.sha256, str)
            or len(self.sha256) != 64
            or any(character not in "0123456789abcdef" for character in self.sha256)
        ):
            raise ConfigurationError(
                "Dinucleotide table sha256 must contain 64 lowercase hexadecimal digits.",
                code="INVALID_DINUCLEOTIDE_TABLE_CHECKSUM",
            )
        if not isinstance(self.properties, tuple) or not all(
            isinstance(item, DinucleotideProperty) for item in self.properties
        ):
            raise ConfigurationError(
                "Dinucleotide table properties must be a tuple of validated properties.",
                code="INVALID_DINUCLEOTIDE_TABLE_PROPERTIES",
            )
        expected = tuple(spec.key for spec in DINUCLEOTIDE_PROPERTY_SPECS)
        observed = tuple(item.key for item in self.properties)
        if observed != expected:
            raise ConfigurationError(
                "Dinucleotide table properties must follow the required 15-property order.",
                code="INVALID_DINUCLEOTIDE_TABLE_PROPERTIES",
                context={"expected": expected, "observed": observed},
            )

    def property(self, key: str) -> DinucleotideProperty:
        """Return one required property by its stable key."""

        for item in self.properties:
            if item.key == key:
                return item
        raise KeyError(key)  # pragma: no cover - construction enforces all keys.


def _object_without_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ConfigurationError(
                "Dinucleotide table JSON contains a duplicate object key.",
                code="DUPLICATE_DINUCLEOTIDE_TABLE_KEY",
                context={"key": key},
            )
        result[key] = value
    return result


def _bounded_metadata(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > _MAX_METADATA_LENGTH:
        raise ConfigurationError(
            f"Dinucleotide table {name} must be non-empty bounded text.",
            code="INVALID_DINUCLEOTIDE_TABLE_METADATA",
            context={"field": name},
        )
    return value.strip()


def _validated_coefficients(
    raw_values: object,
    spec: DinucleotidePropertySpec,
) -> FrozenDict:
    if not isinstance(raw_values, Mapping) or set(raw_values) != set(DINUCLEOTIDES):
        raise ConfigurationError(
            "Each property must define exactly the 16 canonical dinucleotides.",
            code="INCOMPLETE_DINUCLEOTIDE_TABLE_PROPERTY",
            context={"property": spec.key, "required_dinucleotides": DINUCLEOTIDES},
        )
    resolved: dict[str, float] = {}
    for word in DINUCLEOTIDES:
        item = raw_values[word]
        if (
            isinstance(item, bool)
            or not isinstance(item, (int, float))
            or not math.isfinite(item)
            or abs(float(item)) > _MAX_ABSOLUTE_VALUE
        ):
            raise ConfigurationError(
                "Dinucleotide coefficients must be finite bounded numbers.",
                code="INVALID_DINUCLEOTIDE_TABLE_VALUE",
                context={"property": spec.key, "dinucleotide": word},
            )
        resolved[word] = float(item)
    return FrozenDict(resolved)


def _property_values(value: object, spec: DinucleotidePropertySpec) -> FrozenDict:
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            "Each dinucleotide property must be a JSON object.",
            code="INVALID_DINUCLEOTIDE_TABLE_PROPERTY",
            context={"property": spec.key},
        )
    if set(value) != {"unit", "values"}:
        raise ConfigurationError(
            "Each dinucleotide property must contain exactly 'unit' and 'values'.",
            code="INVALID_DINUCLEOTIDE_TABLE_PROPERTY",
            context={"property": spec.key},
        )
    if value["unit"] != spec.unit:
        raise ConfigurationError(
            "Dinucleotide property unit does not match the fixed descriptor schema.",
            code="DINUCLEOTIDE_TABLE_UNIT_MISMATCH",
            context={
                "property": spec.key,
                "expected_unit": spec.unit,
                "observed_unit": value["unit"],
            },
        )
    return _validated_coefficients(value["values"], spec)


def load_dinucleotide_property_table(
    path: str | Path,
    *,
    max_bytes: int = DEFAULT_MAX_DINUCLEOTIDE_TABLE_BYTES,
) -> DinucleotidePropertyTable:
    """Load a bounded JSON table supplied and licensed by the caller.

    Loading is explicit and local. DNAKit does not fetch, generate, cache, or
    redistribute the table and records a SHA-256 digest for every calculation.
    """

    if (
        isinstance(max_bytes, bool)
        or not isinstance(max_bytes, int)
        or not 1 <= max_bytes <= 10_000_000
    ):
        raise ConfigurationError("max_bytes must be an integer in [1, 10000000].")
    try:
        resolved_path = Path(path).expanduser().resolve(strict=True)
        size = resolved_path.stat().st_size
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Dinucleotide property table is unavailable.",
            code="DINUCLEOTIDE_TABLE_UNAVAILABLE",
            context={"path": str(path)},
        ) from exc
    if not resolved_path.is_file():
        raise ConfigurationError(
            "Dinucleotide property table must be a regular file.",
            code="DINUCLEOTIDE_TABLE_UNAVAILABLE",
            context={"path": str(resolved_path)},
        )
    if size > max_bytes:
        raise ConfigurationError(
            "Dinucleotide property table exceeds max_bytes.",
            code="DINUCLEOTIDE_TABLE_SIZE_LIMIT",
            context={"size": size, "max_bytes": max_bytes},
        )
    try:
        with resolved_path.open("rb") as handle:
            payload = handle.read(max_bytes + 1)
    except OSError as exc:
        raise ConfigurationError(
            "Dinucleotide property table could not be read.",
            code="DINUCLEOTIDE_TABLE_UNAVAILABLE",
            context={"path": str(resolved_path)},
        ) from exc
    if len(payload) > max_bytes:
        raise ConfigurationError(
            "Dinucleotide property table exceeds max_bytes.",
            code="DINUCLEOTIDE_TABLE_SIZE_LIMIT",
            context={"size_at_least": len(payload), "max_bytes": max_bytes},
        )
    try:
        decoded = payload.decode("utf-8")
        raw = json.loads(decoded, object_pairs_hook=_object_without_duplicates)
    except ConfigurationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as exc:
        raise ConfigurationError(
            "Dinucleotide property table must be valid UTF-8 JSON.",
            code="INVALID_DINUCLEOTIDE_TABLE_JSON",
            context={"path": str(resolved_path), "error_type": type(exc).__name__},
        ) from exc
    if not isinstance(raw, Mapping) or set(raw) != {
        "schema_version",
        "name",
        "version",
        "source",
        "properties",
    }:
        raise ConfigurationError(
            "Dinucleotide table root fields do not match the required schema.",
            code="INVALID_DINUCLEOTIDE_TABLE_SCHEMA",
        )
    if raw["schema_version"] != DINUCLEOTIDE_TABLE_SCHEMA_VERSION:
        raise ConfigurationError(
            "Unsupported dinucleotide table schema_version.",
            code="UNSUPPORTED_DINUCLEOTIDE_TABLE_SCHEMA",
            context={"schema_version": raw["schema_version"]},
        )
    raw_properties = raw["properties"]
    if not isinstance(raw_properties, Mapping):
        raise ConfigurationError(
            "Dinucleotide table properties must be a JSON object.",
            code="INVALID_DINUCLEOTIDE_TABLE_PROPERTIES",
        )
    expected_keys = {spec.key for spec in DINUCLEOTIDE_PROPERTY_SPECS}
    if set(raw_properties) != expected_keys:
        raise ConfigurationError(
            "Dinucleotide table must contain exactly the 15 required properties.",
            code="INVALID_DINUCLEOTIDE_TABLE_PROPERTIES",
            context={
                "missing": tuple(sorted(expected_keys - set(raw_properties))),
                "extra": tuple(sorted(set(raw_properties) - expected_keys)),
            },
        )
    properties = tuple(
        DinucleotideProperty(
            key=spec.key,
            display_name=spec.display_name,
            unit=spec.unit,
            values=_property_values(raw_properties[spec.key], spec),
        )
        for spec in DINUCLEOTIDE_PROPERTY_SPECS
    )
    return DinucleotidePropertyTable(
        name=_bounded_metadata(raw["name"], "name"),
        version=_bounded_metadata(raw["version"], "version"),
        source=_bounded_metadata(raw["source"], "source"),
        sha256=hashlib.sha256(payload).hexdigest(),
        properties=properties,
    )


__all__ = [
    "DEFAULT_MAX_DINUCLEOTIDE_TABLE_BYTES",
    "DINUCLEOTIDES",
    "DINUCLEOTIDE_PROPERTY_SPECS",
    "DINUCLEOTIDE_TABLE_SCHEMA_VERSION",
    "DinucleotideProperty",
    "DinucleotidePropertySpec",
    "DinucleotidePropertyTable",
    "load_dinucleotide_property_table",
]
