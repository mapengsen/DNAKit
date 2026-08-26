"""Load JSON/TOML configuration with explicit precedence and validation."""

from __future__ import annotations

import json
import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import TYPE_CHECKING, Literal, cast

import yaml  # type: ignore[import-untyped]

if TYPE_CHECKING:
    import tomli as toml_reader
else:
    try:
        import tomllib as toml_reader
    except ImportError:  # pragma: no cover - Python 3.10 compatibility branch.
        import tomli as toml_reader

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.exceptions import ConfigurationError, InputFormatError

ConfigFormat = Literal["json", "toml", "yaml"]
_MAX_CONFIG_BYTES = 1_000_000
_MAX_CONFIG_DEPTH = 32
_MAX_CONFIG_NODES = 100_000
_MAX_CONFIG_KEYS = 100_000


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    values: dict[str, object] = {}
    for key, value in pairs:
        if key in values:
            raise InputFormatError(
                "Configuration contains a duplicate key.",
                code="DUPLICATE_CONFIG_KEY",
                context={"key": key},
            )
        values[key] = value
    return values


class _UniqueKeySafeLoader(yaml.SafeLoader):  # type: ignore[misc]
    """Safe YAML loader that rejects duplicate keys."""


def _construct_unique_mapping(
    loader: yaml.SafeLoader, node: yaml.nodes.MappingNode, deep: bool = False
) -> dict[object, object]:
    loader.flatten_mapping(node)
    values: dict[object, object] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            duplicate = key in values
        except TypeError as exc:
            raise InputFormatError(
                "Configuration keys must be hashable strings.", code="INVALID_CONFIG_KEY"
            ) from exc
        if duplicate:
            raise InputFormatError(
                "Configuration contains a duplicate key.",
                code="DUPLICATE_CONFIG_KEY",
                context={"key": str(key)},
            )
        values[key] = loader.construct_object(value_node, deep=deep)
    return values


_UniqueKeySafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


def _validate_structure(value: object) -> None:
    nodes = [0]
    active: set[int] = set()

    def visit(item: object, depth: int) -> None:
        nodes[0] += 1
        if nodes[0] > _MAX_CONFIG_NODES or depth > _MAX_CONFIG_DEPTH:
            raise InputFormatError(
                "Configuration exceeds structural limits.",
                code="CONFIG_STRUCTURE_LIMIT",
                context={"max_nodes": _MAX_CONFIG_NODES, "max_depth": _MAX_CONFIG_DEPTH},
            )
        if isinstance(item, (Mapping, list, tuple)):
            identity = id(item)
            if identity in active:
                raise InputFormatError(
                    "Configuration contains a recursive reference.",
                    code="CONFIG_RECURSIVE_REFERENCE",
                )
            active.add(identity)
            try:
                if isinstance(item, Mapping):
                    for key in item:
                        visit(key, depth + 1)
                        visit(item[key], depth + 1)
                else:
                    for child in item:
                        visit(child, depth + 1)
            finally:
                active.remove(identity)

    visit(value, 0)


@dataclass(frozen=True, slots=True)
class ConfigLayer:
    """One applied configuration source in precedence order."""

    source: str
    values: FrozenDict

    def __init__(self, source: str, values: Mapping[str, object]) -> None:
        if not isinstance(source, str) or not source.strip():
            raise ConfigurationError("Config layer source must be non-empty.")
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "values", freeze_mapping(values))


@dataclass(frozen=True, slots=True)
class LoadedConfig:
    """Fully resolved configuration plus its ordered source audit."""

    values: FrozenDict
    layers: tuple[ConfigLayer, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.values, FrozenDict):
            raise ConfigurationError("LoadedConfig values must be FrozenDict.")
        if (
            not isinstance(self.layers, tuple)
            or not self.layers
            or any(not isinstance(layer, ConfigLayer) for layer in self.layers)
        ):
            raise ConfigurationError("LoadedConfig layers must contain ConfigLayer objects.")

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], to_json_compatible(self))


def _plain_mapping(value: object, *, source: str) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise InputFormatError(
            "Configuration root must be an object/table.",
            code="INVALID_CONFIG_ROOT",
            context={"source": source},
        )
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str) or not key:
            raise InputFormatError(
                "Configuration keys must be non-empty strings.",
                code="INVALID_CONFIG_KEY",
                context={"source": source},
            )
        result[key] = item
    freeze_mapping(result)
    return result


def _read_file(path: Path, format: ConfigFormat | None) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Configuration file does not exist: {path}")
    if path.is_symlink():
        raise InputFormatError(
            "Configuration must be a regular non-symlink file.",
            code="INVALID_CONFIG_PATH",
            context={"path": str(path)},
        )
    byte_size = path.stat().st_size
    if byte_size > _MAX_CONFIG_BYTES:
        raise InputFormatError(
            "Configuration exceeds the byte limit.",
            code="CONFIG_SIZE_LIMIT",
            context={"byte_size": byte_size, "max_bytes": _MAX_CONFIG_BYTES},
        )
    resolved_format = format
    if resolved_format is None:
        suffix = path.suffix.lower()
        resolved_format = (
            "toml"
            if suffix == ".toml"
            else "json"
            if suffix == ".json"
            else "yaml"
            if suffix in {".yaml", ".yml"}
            else None
        )
    if resolved_format not in ("json", "toml", "yaml"):
        raise ConfigurationError(
            "Configuration format must be json, toml, or yaml.",
            code="UNSUPPORTED_CONFIG_FORMAT",
            context={"path": str(path), "format": format},
        )
    try:
        with path.open("rb") as handle:
            raw = handle.read(_MAX_CONFIG_BYTES + 1)
            if len(raw) > _MAX_CONFIG_BYTES:
                raise InputFormatError(
                    "Configuration grew beyond the byte limit while being read.",
                    code="CONFIG_SIZE_LIMIT",
                    context={"max_bytes": _MAX_CONFIG_BYTES},
                )
            if resolved_format == "json":
                parsed: object = json.loads(
                    raw.decode("utf-8"), object_pairs_hook=_unique_json_object
                )
            elif resolved_format == "toml":
                parsed = toml_reader.loads(raw.decode("utf-8"))
            else:
                text = raw.decode("utf-8")
                if any(
                    isinstance(token, yaml.tokens.AliasToken)
                    for token in yaml.scan(text, Loader=_UniqueKeySafeLoader)
                ):
                    raise InputFormatError(
                        "Configuration YAML aliases are disabled.",
                        code="CONFIG_YAML_ALIAS_DISABLED",
                    )
                parsed = yaml.load(text, Loader=_UniqueKeySafeLoader)
    except (
        UnicodeError,
        json.JSONDecodeError,
        toml_reader.TOMLDecodeError,
        yaml.YAMLError,
        ValueError,
        RecursionError,
    ) as exc:
        raise InputFormatError(
            "Configuration file could not be parsed.",
            code="INVALID_CONFIG_FILE",
            context={"path": str(path), "format": resolved_format},
        ) from exc
    _validate_structure(parsed)
    return _plain_mapping(parsed, source=str(path))


def _environment_values(
    environment: Mapping[str, str],
    *,
    prefix: str,
    allowed_keys: frozenset[str],
) -> dict[str, object]:
    if not isinstance(environment, Mapping) or any(
        not isinstance(name, str) or not isinstance(value, str)
        for name, value in environment.items()
    ):
        raise ConfigurationError(
            "environment must map string names to string values.",
            code="INVALID_CONFIG_ENVIRONMENT",
        )
    result: dict[str, object] = {}
    for name in sorted(environment):
        if not name.startswith(prefix):
            continue
        key = name[len(prefix) :].lower()
        if not key or key not in allowed_keys:
            raise ConfigurationError(
                "Unknown DNAKit environment configuration key.",
                code="UNKNOWN_CONFIG_KEY",
                context={"environment_variable": name},
            )
        raw = environment[name]
        if len(raw.encode("utf-8")) > _MAX_CONFIG_BYTES:
            raise InputFormatError(
                "Environment configuration value exceeds the byte limit.",
                code="CONFIG_ENVIRONMENT_SIZE_LIMIT",
                context={"environment_variable": name, "max_bytes": _MAX_CONFIG_BYTES},
            )
        try:
            result[key] = json.loads(raw)
        except json.JSONDecodeError:
            result[key] = raw
        except RecursionError as exc:
            raise InputFormatError(
                "Environment configuration value exceeds structural limits.",
                code="CONFIG_STRUCTURE_LIMIT",
                context={"environment_variable": name, "max_depth": _MAX_CONFIG_DEPTH},
            ) from exc
        _validate_structure(result[key])
    return result


def _validated_layer(
    source: str,
    values: Mapping[str, object],
    *,
    allowed_keys: frozenset[str],
) -> ConfigLayer:
    _validate_structure(values)
    unknown = sorted(set(values) - allowed_keys)
    if unknown:
        raise ConfigurationError(
            "Configuration contains unknown keys.",
            code="UNKNOWN_CONFIG_KEY",
            context={"source": source, "keys": unknown},
        )
    return ConfigLayer(source, values)


def load_config(
    *,
    defaults: Mapping[str, object],
    allowed_keys: Iterable[str],
    path: str | os.PathLike[str] | None = None,
    format: ConfigFormat | None = None,
    environment: Mapping[str, str] | None = None,
    env_prefix: str = "DNAKIT_",
    cli: Mapping[str, object] | None = None,
) -> LoadedConfig:
    """Resolve defaults < file < environment < explicit CLI values."""

    if not isinstance(defaults, Mapping):
        raise ConfigurationError("defaults must be a mapping.", code="INVALID_CONFIG_DEFAULTS")
    if isinstance(allowed_keys, (str, bytes)):
        raise ConfigurationError("allowed_keys must be an iterable of configuration names.")
    try:
        allowed_items = tuple(islice(iter(allowed_keys), _MAX_CONFIG_KEYS + 1))
    except TypeError as exc:
        raise ConfigurationError(
            "allowed_keys must be an iterable of configuration names."
        ) from exc
    if not allowed_items or any(not isinstance(key, str) or not key for key in allowed_items):
        raise ConfigurationError("allowed_keys must contain non-empty strings.")
    if len(allowed_items) > _MAX_CONFIG_KEYS:
        raise ConfigurationError(
            "allowed_keys exceeds the configured limit.",
            code="CONFIG_KEY_LIMIT",
            context={"max_keys": _MAX_CONFIG_KEYS},
        )
    allowed = frozenset(allowed_items)
    if not isinstance(env_prefix, str) or not env_prefix:
        raise ConfigurationError("env_prefix must be non-empty.")
    if format is not None and format not in ("json", "toml", "yaml"):
        raise ConfigurationError(
            "Configuration format must be json, toml, or yaml.",
            code="UNSUPPORTED_CONFIG_FORMAT",
            context={"format": format},
        )
    if environment is not None and not isinstance(environment, Mapping):
        raise ConfigurationError(
            "environment must be a mapping or None.", code="INVALID_CONFIG_ENVIRONMENT"
        )
    if cli is not None and not isinstance(cli, Mapping):
        raise ConfigurationError("cli must be a mapping or None.", code="INVALID_CONFIG_CLI")
    layers: list[ConfigLayer] = [_validated_layer("defaults", defaults, allowed_keys=allowed)]
    if path is not None:
        try:
            config_path = Path(path)
        except TypeError as exc:
            raise ConfigurationError(
                "path must be a filesystem path or None.", code="INVALID_CONFIG_PATH"
            ) from exc
        layers.append(
            _validated_layer(
                str(config_path), _read_file(config_path, format), allowed_keys=allowed
            )
        )
    env_values = _environment_values(
        os.environ if environment is None else environment,
        prefix=env_prefix,
        allowed_keys=allowed,
    )
    if env_values:
        layers.append(_validated_layer("environment", env_values, allowed_keys=allowed))
    if cli:
        layers.append(_validated_layer("cli", cli, allowed_keys=allowed))
    resolved: dict[str, object] = {}
    for layer in layers:
        resolved.update(layer.values)
    return LoadedConfig(freeze_mapping(resolved), tuple(layers))


__all__ = ["ConfigLayer", "LoadedConfig", "load_config"]
