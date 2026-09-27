"""Versioned, text-only local descriptions for native and upstream MCP tools."""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from dnakit.exceptions import ConfigurationError

if TYPE_CHECKING:
    from .catalog import AgentToolSpec

_DESCRIPTION_PATH = Path(__file__).with_name("_data") / "descriptions.json"
_TEXT_FIELDS = frozenset({"description", "full_description", "parameters", "examples"})
_EXAMPLES = "\n\nDNAKit usage examples:\n"


def description_path(root: Path) -> Path:
    return root / "src" / "dnakit" / "tools" / "_data" / "descriptions.json"


def _json(value: object) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False)


def read_descriptions(path: Path | None = None) -> dict[str, Any]:
    path = _DESCRIPTION_PATH if path is None else path
    if not path.exists():
        return {"schema_version": 1, "native": {}, "mcp": {}}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        if (
            not isinstance(value, dict)
            or value.get("schema_version") != 1
            or not isinstance(value.get("native"), dict)
            or not isinstance(value.get("mcp"), dict)
        ):
            raise ValueError("Invalid description metadata")
        return value
    except (OSError, ValueError) as exc:
        raise ConfigurationError(
            "Local tool description metadata is unreadable or invalid.",
            code="TOOL_DESCRIPTION_INVALID",
            context={"path": str(path)},
        ) from exc


def schema_contract(schema: object) -> object:
    """Remove schema annotations while preserving property names and defaults."""
    if not isinstance(schema, dict):
        return deepcopy(schema)
    result: dict[str, Any] = {}
    for key, value in schema.items():
        if key in {"description", "examples"}:
            continue
        if key in {
            "properties",
            "patternProperties",
            "$defs",
            "definitions",
            "dependentSchemas",
        } and isinstance(value, dict):
            result[key] = {name: schema_contract(child) for name, child in value.items()}
        elif key in {"items", "additionalProperties", "contains", "not", "if", "then", "else"}:
            result[key] = schema_contract(value)
        elif key in {"allOf", "anyOf", "oneOf", "prefixItems"} and isinstance(value, list):
            result[key] = [schema_contract(child) for child in value]
        else:
            result[key] = deepcopy(value)
    return result


def description_contract(tool: dict[str, Any]) -> str:
    fields: dict[str, Any] = {
        key: tool[key]
        for key in (
            "name",
            "function",
            "signature",
            "effect",
            "requires_confirmation",
            "open_world",
            "agent_compatible",
            "incompatibilities",
            "python_only_optional_parameters",
            "annotations",
        )
        if key in tool
    }
    for key in ("input_schema", "output_schema", "inputSchema", "outputSchema"):
        if key in tool:
            fields[key] = schema_contract(tool[key])
    return hashlib.sha256(_json(fields).encode("utf-8")).hexdigest()


def validate_changes(changes: object, schema: dict[str, Any], *, native: bool) -> dict[str, Any]:
    allowed = _TEXT_FIELDS if native else _TEXT_FIELDS - {"full_description"}
    if not isinstance(changes, dict) or not changes or set(changes) - allowed:
        raise ConfigurationError(
            "Only tool text, parameter descriptions and examples may change.",
            code="TOOL_DESCRIPTION_CONTRACT_CHANGE",
        )
    for field in ("description", "full_description"):
        if field in changes and (
            not isinstance(changes[field], str)
            or not changes[field].strip()
            or len(changes[field]) > 32_000
        ):
            raise ConfigurationError(
                "Description must be bounded text.", code="TOOL_DESCRIPTION_INVALID"
            )
    parameters = changes.get("parameters", {})
    if not isinstance(parameters, dict) or any(
        name not in schema.get("properties", {})
        or not isinstance(text, str)
        or not text.strip()
        or len(text) > 8000
        for name, text in parameters.items()
    ):
        raise ConfigurationError(
            "Parameter changes must name existing parameters and contain description text only.",
            code="TOOL_DESCRIPTION_CONTRACT_CHANGE",
        )
    if "examples" in changes and (
        not isinstance(changes["examples"], list) or len(changes["examples"]) > 20
    ):
        raise ConfigurationError(
            "examples must be a list of at most 20 values.", code="TOOL_DESCRIPTION_INVALID"
        )
    try:
        if len(_json(changes).encode("utf-8")) > 100_000:
            raise ValueError("Too large")
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(
            "Descriptions must be bounded finite JSON.", code="TOOL_DESCRIPTION_INVALID"
        ) from exc
    return deepcopy(changes)


def patch_description(
    tool: dict[str, Any], changes: dict[str, Any], *, native: bool
) -> dict[str, Any]:
    schema_key = "input_schema" if native else "inputSchema"
    result = deepcopy(tool)
    changes = validate_changes(changes, result.get(schema_key, {}), native=native)
    field = "full_description" if native else "description"
    old_text, separator, old_examples = result.get(field, "").partition(_EXAMPLES)
    result[field] = old_text
    for key in ("description", "full_description"):
        if key in changes:
            result[key] = changes[key]
    if native:
        if "full_description" in changes and "description" not in changes:
            result["description"] = changes["full_description"].split("\n\n", 1)[0][:500]
        elif "description" in changes and "full_description" not in changes:
            paragraphs = result["full_description"].split("\n\n", 1)
            result["full_description"] = changes["description"] + (
                "\n\n" + paragraphs[1] if len(paragraphs) == 2 else ""
            )
    for name, text in changes.get("parameters", {}).items():
        result[schema_key]["properties"][name]["description"] = text
    if "examples" in changes:
        old_examples = _json(changes["examples"])
        separator = _EXAMPLES
    if separator:
        result[field] = result.get(field, "").partition(_EXAMPLES)[0] + _EXAMPLES + old_examples
    if native or "description" in result[schema_key]:
        result[schema_key]["description"] = result[field]
    return result


def _active(entry: object, contract: str) -> dict[str, Any] | None:
    if not isinstance(entry, dict) or not isinstance(entry.get("versions"), list):
        return None
    for version in reversed(entry["versions"]):
        if isinstance(version, dict) and version.get("id") == entry.get("active"):
            if version.get("contract") == contract and isinstance(version.get("changes"), dict):
                return deepcopy(cast(dict[str, Any], version["changes"]))
            return None
    return None


def append_description_version(
    store: dict[str, Any],
    *,
    kind: str,
    name: str,
    server_id: str | None,
    contract: str,
    changes: dict[str, Any],
    change_id: str,
) -> dict[str, Any]:
    result = deepcopy(store)
    entries = result["native"] if kind == "native" else result["mcp"].setdefault(server_id, {})
    previous = entries.get(name, {})
    merged = _active(previous, contract) or {}
    merged_parameters = {**merged.get("parameters", {}), **changes.get("parameters", {})}
    merged.update(deepcopy(changes))
    if kind == "native":
        if (
            "description" in changes
            and "full_description" not in changes
            and "full_description" in merged
        ):
            paragraphs = merged["full_description"].split("\n\n", 1)
            merged["full_description"] = changes["description"] + (
                "\n\n" + paragraphs[1] if len(paragraphs) == 2 else ""
            )
        elif "full_description" in changes and "description" not in changes:
            merged["description"] = changes["full_description"].split("\n\n", 1)[0][:500]
    if merged_parameters:
        merged["parameters"] = merged_parameters
    entries[name] = {
        "active": change_id,
        "versions": [
            *previous.get("versions", []),
            {
                "id": change_id,
                "parent": previous.get("active"),
                "contract": contract,
                "changes": merged,
            },
        ],
    }
    return result


def apply_native_description(spec: AgentToolSpec) -> AgentToolSpec:
    """Return an independent spec; callable, argument types and defaults are untouched."""
    tool = spec.to_dict()
    changes = _active(read_descriptions()["native"].get(spec.name), description_contract(tool))
    if changes is None:
        return spec
    updated = patch_description(tool, changes, native=True)
    properties = updated["input_schema"].get("properties", {})
    plan = replace(
        spec.call_plan,
        input_schema=deepcopy(updated["input_schema"]),
        parameters=tuple(
            replace(parameter, schema=deepcopy(properties.get(parameter.name, parameter.schema)))
            for parameter in spec.call_plan.parameters
        ),
    )
    return replace(
        spec,
        description=updated["description"],
        full_description=updated["full_description"],
        call_plan=plan,
        search_text=_json(
            {
                "description": updated["description"],
                "full_description": updated["full_description"],
                "parameters": {
                    name: schema.get("description", "") for name, schema in properties.items()
                },
            }
        ),
    )


def overlay_mcp_tools(server_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Overlay local listing text only; live MCP invocation schemas stay upstream-owned."""
    entries = read_descriptions()["mcp"].get(server_id, {})
    result = deepcopy(payload)
    for index, tool in enumerate(result.get("tools", [])):
        changes = _active(entries.get(tool.get("name")), description_contract(tool))
        if changes is not None:
            result["tools"][index] = patch_description(tool, changes, native=False)
    return result
