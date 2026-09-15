"""DNAKit-owned dispatch and validation for directly bundled database adapters."""

from __future__ import annotations

import copy
import importlib
import json
import os
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from typing import Any

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.search._http import redact_url


@lru_cache(maxsize=1)
def definitions() -> dict[str, Any]:
    entries: dict[str, Any] = json.loads(
        Path(__file__).with_name("providers.json").read_text(encoding="utf-8")
    )
    from ._routes import EXTRA_FUNCTIONS

    extra = {
        "pdb_get_structure": (
            {
                "pdb_id": {"type": "string"},
                "format": {"type": "string", "default": "cif"},
                "max_bytes": {"type": "integer", "default": 10_000_000},
            },
            ["pdb_id"],
            "PDB_get_structure",
        ),
        "orphanet_get_disease_genes": (
            {"orpha_code": {"type": "string"}},
            ["orpha_code"],
            "Orphanet_get_disease_genes",
        ),
        "cadd_get_scores": (
            {"variant": {"type": "string"}, "assembly": {"type": "string", "default": "hg38"}},
            ["variant"],
            "CADD_get_scores",
        ),
        "string_get_functional_enrichment": (
            {
                "identifiers": {"type": "string"},
                "species": {"type": "integer", "default": 9606},
                "category": {"type": "string"},
            },
            ["identifiers"],
            "STRING_get_functional_enrichment",
        ),
    }
    for name, (properties, required, alias) in extra.items():
        schema = {"type": "object", "properties": properties, "required": required}
        entries[name] = {
            "name": name,
            "upstream_name": alias,
            "native_function": True,
            "effect": "network",
            "config": {
                "name": name,
                "description": EXTRA_FUNCTIONS[name].__doc__,
                "parameter": schema,
            },
            "input_schema": schema,
            "implementation": "DNAKit workflow extension",
        }
    return entries


def describe(operation: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the function's actual upstream API schema and implementation provenance."""
    entries = definitions()
    if operation in entries:
        return copy.deepcopy(entries[operation])
    found = [entry for entry in entries.values() if entry["upstream_name"] == operation]
    if len(found) == 1:
        return copy.deepcopy(found[0])
    raise ConfigurationError(
        f"Unknown DNAKit database operation: {operation}", code="UNKNOWN_DATABASE_OPERATION"
    )


def search(query: str, *, limit: int = 10, verbose: bool = True) -> list[dict[str, str]]:
    """Find bundled scientific API operations locally by task, name or description."""
    import re

    if (
        not isinstance(query, str)
        or not query.strip()
        or not isinstance(limit, int)
        or isinstance(limit, bool)
        or not 1 <= limit <= 100
    ):
        raise ConfigurationError("Provide a nonempty query and limit between 1 and 100.")
    terms = set(re.findall(r"\w+", query.casefold()))
    matches = []
    for name, entry in definitions().items():
        description = entry["config"].get("description", "")
        corpus = (name + " " + entry["upstream_name"] + " " + description).casefold()
        score = sum(term in corpus for term in terms)
        if score:
            matches.append(
                (
                    score,
                    {
                        "name": name,
                        "description": description,
                        "upstream_name": entry["upstream_name"],
                    },
                )
            )
    return [item for _, item in sorted(matches, key=lambda row: (-row[0], row[1]["name"]))[:limit]]


def list_operations(*, verbose: bool = True) -> list[dict[str, str]]:
    """List native operation names and recorded upstream aliases without loading providers."""
    return [
        {"name": name, "upstream_name": entry["upstream_name"]}
        for name, entry in definitions().items()
    ]


def execute(
    request: str | Mapping[str, object],
    arguments: Mapping[str, object] | None = None,
    *,
    verbose: bool = True,
) -> object:
    """Execute an operation name plus arguments, or a JSON name/arguments request."""
    if isinstance(request, str):
        return call(request, arguments)
    if not isinstance(request, Mapping) or not isinstance(request.get("name"), str):
        raise ConfigurationError("request must contain an operation name and arguments.")
    if arguments is not None or set(request) - {"name", "arguments"}:
        raise ConfigurationError("Use either a request mapping or separate arguments.")
    operation, payload = request["name"], request.get("arguments")
    assert isinstance(operation, str)
    if payload is not None and not isinstance(payload, Mapping):
        raise ConfigurationError("arguments must be a mapping.")
    return call(operation, payload)


def _redact(value: object, secrets: tuple[str, ...]) -> object:
    if isinstance(value, Mapping):
        return {
            str(key): (
                "REDACTED"
                if str(key).lower()
                in {
                    "api_key",
                    "apikey",
                    "api-key",
                    "authorization",
                    "access_token",
                    "token",
                    "password",
                }
                else _redact(item, secrets)
            )
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_redact(item, secrets) for item in value]
    if isinstance(value, str):
        if value.startswith(("https://", "http://")) and "?" in value:
            value = redact_url(value)
        for secret in secrets:
            value = value.replace(secret, "REDACTED")
    return value


def call(
    operation: str,
    arguments: Mapping[str, object] | None = None,
    *,
    max_result_bytes: int = 20_000_000,
    verbose: bool = True,
) -> object:
    """Call a bundled DNAKit database adapter without installing ToolUniverse.

    An operation may use its native snake_case name or recorded upstream name.
    Input JSON is checked against the pinned service schema; returned provider
    envelopes preserve the service's fields, pagination and evidence details.
    Failure raises an exception, including provider error envelopes. Models,
    reference data, API credentials and specialist libraries remain external.
    """
    from jsonschema import Draft7Validator

    entry = describe(operation)
    if arguments is not None and not isinstance(arguments, Mapping):
        raise ConfigurationError("arguments must be a mapping.", code="INVALID_DATABASE_ARGUMENTS")
    if (
        isinstance(max_result_bytes, bool)
        or not isinstance(max_result_bytes, int)
        or not 1 <= max_result_bytes <= 100_000_000
    ):
        raise ConfigurationError("max_result_bytes must be between 1 and 100,000,000.")
    values = copy.deepcopy(dict(arguments or {}))
    from ._extensions import finish, prepare

    extension_schema = entry.get("input_schema")
    if extension_schema is not None:
        check = copy.deepcopy(extension_schema)
        for keyword in ("required", "oneOf", "anyOf", "allOf"):
            check.pop(keyword, None)
        check.setdefault("additionalProperties", False)
        error = next(Draft7Validator(check).iter_errors(values), None)
        if error:
            raise ConfigurationError(
                f"Invalid arguments for {entry['name']}; inspect its input_schema.",
                code="INVALID_DATABASE_ARGUMENTS",
            )
    from ._routes import route

    routed = route(entry["name"], values)
    config = entry["config"]
    extension_options: dict[str, Any] = {}
    if routed is None:
        values, extension_options = prepare(entry, values)
        schema = copy.deepcopy(config.get("parameter") or {"type": "object"})
        properties = schema.get("properties") or {}
        for name, item in properties.items():
            if name not in values and "default" in item:
                values[name] = copy.deepcopy(item["default"])
        fixed = (config.get("fields") or {}).get("operation")
        if fixed:
            values.setdefault("operation", fixed)
            if "operation" not in properties:
                properties["operation"] = {"type": "string", "const": fixed}
        schema["required"] = sorted(
            set(schema.get("required") or [])
            | {key for key, item in properties.items() if item.get("required") is True}
        )
        schema["properties"] = properties
        schema.setdefault("type", "object")
        schema.setdefault("additionalProperties", False)
        error = next(Draft7Validator(schema).iter_errors(values), None)
        if error:
            path = ".".join(map(str, error.path)) or "arguments"
            raise ConfigurationError(
                f"Invalid {path} for {entry['name']}; inspect its input schema.",
                code="INVALID_DATABASE_ARGUMENTS",
            )
    required_keys = config.get("required_api_keys") or []
    secrets = tuple(
        value
        for key, value in os.environ.items()
        if len(value) >= 6
        and (key in required_keys or key.endswith(("_API_KEY", "_API_TOKEN", "_JWT")))
    )
    secrets += tuple(
        str(value)
        for key, value in values.items()
        if key.lower() in {"api_key", "token", "password", "access_token"} and len(str(value)) >= 6
    )
    try:
        if routed is not None:
            result = routed()
        elif entry.get("native_function"):
            from ._routes import EXTRA_FUNCTIONS

            result = EXTRA_FUNCTIONS[entry["name"]](**values)
        else:
            module = importlib.import_module("dnakit.databases._providers." + entry["module"])
            implementation = getattr(module, entry["class"])(copy.deepcopy(config))
            result = implementation.run(values)
    except ImportError as exc:
        dependency = getattr(exc, "name", None) or "an optional provider dependency"
        raise BackendUnavailableError(
            f"{entry['name']} requires {dependency}; "
            "prepare the scientific dependencies listed by the skill.",
            code="DATABASE_DEPENDENCY_UNAVAILABLE",
        ) from exc
    except Exception as exc:
        message = _redact(str(exc), secrets)
        raise BackendExecutionError(
            f"{entry['name']}: {message}", code="DATABASE_PROVIDER_ERROR"
        ) from exc
    if result is None:
        raise BackendExecutionError(
            f"{entry['name']} returned no result.", code="EMPTY_DATABASE_RESULT"
        )
    if isinstance(result, Mapping) and (
        result.get("error")
        or str(result.get("status", "")).lower() in {"error", "failed", "failure"}
    ):
        raise BackendExecutionError(
            str(_redact(result.get("error") or "Provider reported failure", secrets)),
            code="DATABASE_PROVIDER_ERROR",
        )
    if entry["name"] == "ch_ipatlas_enrichment_analysis" and isinstance(result, dict):
        result["status"] = "requires_submission"
        result["analysis_completed"] = False
        result["scope"] = "Submission instructions only; no enrichment statistics computed."
    result = finish(result, extension_options)
    result = _redact(result, secrets)
    try:
        encoded = json.dumps(result, ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise BackendExecutionError(
            "Provider returned a non-JSON result.", code="INVALID_DATABASE_RESULT"
        ) from exc
    if len(encoded) > max_result_bytes:
        raise BackendExecutionError(
            "Provider result exceeds max_result_bytes.", code="DATABASE_RESULT_LIMIT"
        )
    return result
