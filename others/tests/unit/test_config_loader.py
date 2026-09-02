"""Tests for strict layered workflow configuration."""

import json
import sys
from importlib import import_module
from pathlib import Path

import pytest

from dnakit.config import load_config
from dnakit.exceptions import ConfigurationError, InputFormatError


def test_toml_reader_matches_supported_python_runtime() -> None:
    reader = vars(import_module("dnakit.config.loader"))["toml_reader"]
    if sys.version_info >= (3, 11):
        assert reader.__name__ == "tomllib"
    else:
        assert reader.__name__ == "tomli"


def test_config_precedence_and_source_audit(tmp_path: Path) -> None:
    path = tmp_path / "run.toml"
    path.write_text('jobs = 2\nmode = "file"\n', encoding="utf-8")

    loaded = load_config(
        defaults={"jobs": 1, "mode": "default", "seed": None},
        allowed_keys={"jobs", "mode", "seed"},
        path=path,
        environment={"DNAKIT_JOBS": "3", "IGNORED": "x"},
        cli={"mode": "cli"},
    )

    assert loaded.values == {"jobs": 3, "mode": "cli", "seed": None}
    assert [layer.source for layer in loaded.layers] == [
        "defaults",
        str(path),
        "environment",
        "cli",
    ]
    json.dumps(loaded.to_dict(), sort_keys=True)


def test_json_config_and_environment_scalar_decoding(tmp_path: Path) -> None:
    path = tmp_path / "run.json"
    path.write_text('{"progress": false}', encoding="utf-8")

    loaded = load_config(
        defaults={"progress": True, "label": "default"},
        allowed_keys={"progress", "label"},
        path=path,
        environment={"DNAKIT_LABEL": "plain-text"},
    )

    assert loaded.values == {"progress": False, "label": "plain-text"}


def test_yaml_config_is_safe_and_supported(tmp_path: Path) -> None:
    path = tmp_path / "run.yaml"
    path.write_text("jobs: 4\nmode: yaml\n", encoding="utf-8")

    loaded = load_config(
        defaults={"jobs": 1},
        allowed_keys={"jobs", "mode"},
        path=path,
        environment={},
    )

    assert loaded.values == {"jobs": 4, "mode": "yaml"}


def test_unknown_and_invalid_configuration_are_rejected(tmp_path: Path) -> None:
    unknown = tmp_path / "unknown.json"
    unknown.write_text('{"typo": 1}', encoding="utf-8")
    invalid = tmp_path / "invalid.toml"
    invalid.write_text("bad = [", encoding="utf-8")

    with pytest.raises(ConfigurationError) as file_error:
        load_config(defaults={}, allowed_keys={"known"}, path=unknown)
    assert file_error.value.code == "UNKNOWN_CONFIG_KEY"
    with pytest.raises(ConfigurationError) as env_error:
        load_config(
            defaults={},
            allowed_keys={"known"},
            environment={"DNAKIT_TYPO": "1"},
        )
    assert env_error.value.code == "UNKNOWN_CONFIG_KEY"
    with pytest.raises(InputFormatError) as parse_error:
        load_config(defaults={}, allowed_keys={"known"}, path=invalid)
    assert parse_error.value.code == "INVALID_CONFIG_FILE"


def test_invalid_yaml_and_public_argument_types_are_normalized(tmp_path: Path) -> None:
    invalid_yaml = tmp_path / "invalid.yaml"
    invalid_yaml.write_text("key: [", encoding="utf-8")

    with pytest.raises(InputFormatError) as yaml_error:
        load_config(defaults={}, allowed_keys={"key"}, path=invalid_yaml)
    assert yaml_error.value.code == "INVALID_CONFIG_FILE"
    with pytest.raises(ConfigurationError) as format_error:
        load_config(defaults={}, allowed_keys={"key"}, format="ini")  # type: ignore[arg-type]
    assert format_error.value.code == "UNSUPPORTED_CONFIG_FORMAT"
    with pytest.raises(ConfigurationError) as environment_error:
        load_config(
            defaults={},
            allowed_keys={"key"},
            environment={"DNAKIT_KEY": 1},  # type: ignore[dict-item]
        )
    assert environment_error.value.code == "INVALID_CONFIG_ENVIRONMENT"


@pytest.mark.parametrize(
    ("name", "content", "code"),
    [
        ("alias.yaml", "a: &a [*a]\n", "CONFIG_YAML_ALIAS_DISABLED"),
        ("duplicate.yaml", "key: one\nkey: two\n", "DUPLICATE_CONFIG_KEY"),
        ("duplicate.json", '{"key": 1, "key": 2}', "DUPLICATE_CONFIG_KEY"),
    ],
)
def test_config_rejects_aliases_and_duplicate_keys(
    tmp_path: Path, name: str, content: str, code: str
) -> None:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")

    with pytest.raises(InputFormatError) as error:
        load_config(defaults={}, allowed_keys={"key", "a"}, path=path, environment={})
    assert error.value.code == code


def test_config_rejects_byte_and_depth_limits(tmp_path: Path) -> None:
    large = tmp_path / "large.json"
    large.write_text('{"key":"' + "x" * 1_000_001 + '"}', encoding="utf-8")
    deep = tmp_path / "deep.json"
    value: object = 1
    for _ in range(40):
        value = [value]
    deep.write_text(json.dumps({"key": value}), encoding="utf-8")

    with pytest.raises(InputFormatError) as size_error:
        load_config(defaults={}, allowed_keys={"key"}, path=large, environment={})
    assert size_error.value.code == "CONFIG_SIZE_LIMIT"
    with pytest.raises(InputFormatError) as depth_error:
        load_config(defaults={}, allowed_keys={"key"}, path=deep, environment={})
    assert depth_error.value.code == "CONFIG_STRUCTURE_LIMIT"


@pytest.mark.parametrize("source", ["defaults", "environment", "cli"])
def test_config_rejects_deep_non_file_layers(source: str) -> None:
    value: object = 1
    for _ in range(40):
        value = [value]
    kwargs: dict[str, object] = {"defaults": {}, "allowed_keys": {"key"}, "environment": {}}
    if source == "defaults":
        kwargs["defaults"] = {"key": value}
    elif source == "environment":
        kwargs["environment"] = {"DNAKIT_KEY": json.dumps(value)}
    else:
        kwargs["cli"] = {"key": value}

    with pytest.raises(InputFormatError) as error:
        load_config(**kwargs)  # type: ignore[arg-type]
    assert error.value.code == "CONFIG_STRUCTURE_LIMIT"


def test_config_allows_shared_but_non_recursive_values() -> None:
    shared = [1, 2]
    loaded = load_config(
        defaults={"key": [shared, shared]},
        allowed_keys={"key"},
        environment={},
    )
    assert loaded.values["key"] == ((1, 2), (1, 2))


def test_config_rejects_oversized_environment_value() -> None:
    with pytest.raises(InputFormatError) as error:
        load_config(
            defaults={},
            allowed_keys={"key"},
            environment={"DNAKIT_KEY": "x" * 1_000_001},
        )
    assert error.value.code == "CONFIG_ENVIRONMENT_SIZE_LIMIT"
