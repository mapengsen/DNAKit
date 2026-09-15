"""Port the pinned, skill-referenced scientific adapters into DNAKit.

This imports implementation source, not an installed ToolUniverse runtime.
Preserve Apache notices and original file hashes. Public callables and resource
schemas are generated from the pinned definitions, with lazy provider imports.
"""

from __future__ import annotations

import ast
import hashlib
import json
import keyword
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "temp/all-skills-migration"
SOURCE = ROOT / "temp/ToolUniverse/src/tooluniverse"
DEST = ROOT / "src/dnakit/databases/_providers"
PREFIX = "dnakit.databases._providers"
REVISION = "7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7"


def native_name(name: str) -> str:
    value = re.sub(r"(?<=[a-z0-9])(?=[A-Z])", "_", name).lower()
    value = re.sub(r"[^a-z0-9_]", "_", value)
    return value + "_" if keyword.iskeyword(value) else value


def annotation(schema: dict) -> str:
    kind = schema.get("type")
    if isinstance(kind, list):
        return " | ".join(dict.fromkeys(annotation({"type": item}) for item in kind))
    return {
        "string": "str",
        "integer": "int",
        "number": "float",
        "boolean": "bool",
        "object": "dict[str, _Object]",
        "null": "None",
    }.get(
        kind,
        ("list[" + annotation(schema.get("items") or {}) + "]" if kind == "array" else "_Object"),
    )


def write_functions(entries: dict) -> None:
    lines = [
        '"""Generated DNAKit database APIs; implementations are loaded only when called."""',
        "from __future__ import annotations",
        "from builtins import object as _Object",
        "",
        "from ._runtime import call as _call",
        "",
    ]
    for name, entry in entries.items():
        schema = entry.get("input_schema") or entry["config"].get("parameter") or {}
        required = set(schema.get("required") or [])
        properties = schema.get("properties") or {}
        required.update(key for key, value in properties.items() if value.get("required") is True)
        required -= {key for key, value in properties.items() if "default" in value}
        if (entry["config"].get("fields") or {}).get("operation"):
            required.discard("operation")
        aliases = {
            key: (re.sub(r"\W", "_", key) + ("_" if keyword.iskeyword(key) else ""))
            for key in properties
        }
        if any(not alias.isidentifier() for alias in aliases.values()) or len(
            set(aliases.values())
        ) != len(aliases):
            raise ValueError(f"Invalid Python argument names: {name}")
        args = [
            f"    {aliases[key]}: {annotation(value)}"
            + ("" if key in required else " | None = None")
            + ","
            for key, value in properties.items()
        ]
        args.append("    verbose: bool = True,")
        lines.extend(["", f"def {name}(" + ("*," if args else ""), *args, ") -> _Object:"])
        description = str(entry["config"].get("description", name)).replace('"""', "'''")
        lines.append(
            "    "
            + repr(
                description
                + "\n\nBundled adapter provenance: mims-harvard/ToolUniverse, Apache-2.0."
            )
        )
        lines.append("    _dnakit_arguments: dict[str, _Object] = {}")
        for key in properties:
            if key not in required:
                lines.append(f"    if {aliases[key]} is not None:")
            indent = "        " if key not in required else "    "
            lines.append(f"{indent}_dnakit_arguments[{key!r}] = {aliases[key]}")
        lines.append(f"    return _call({name!r}, _dnakit_arguments, verbose=verbose)")
    lines.extend(["", "__all__ = " + repr(sorted(entries)), ""])
    (DEST.parent / "_generated.py").write_text("\n".join(lines), encoding="utf-8")
    (DEST.parent / "__init__.py").write_text(
        '"""DNAKit scientific database APIs with local, attributed provider implementations."""\n\nfrom ._generated import *  # noqa: F403\nfrom ._generated import __all__ as _functions\nfrom ._routes import (\n    cadd_get_scores,\n    orphanet_get_disease_genes,\n    pdb_get_structure,\n    string_get_functional_enrichment,\n)\nfrom ._runtime import call, describe, execute, list_operations, search\n\n__all__ = [\n    "call",\n    "describe",\n    "search",\n    "execute",\n    "list_operations",\n    *_functions,\n    "pdb_get_structure",\n    "orphanet_get_disease_genes",\n    "cadd_get_scores",\n    "string_get_functional_enrichment",\n]\n\nfrom dnakit._verbosity import register_public as _register_public_verbose\n\n_register_public_verbose(globals())\n'
    )


def main() -> None:
    import subprocess

    WORK.mkdir(parents=True, exist_ok=True)

    revision = subprocess.check_output(
        ["git", "-C", str(SOURCE.parents[1]), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != REVISION:
        raise ValueError("Provider source revision differs from reviewed registry")
    plan = json.loads((ROOT / "registry/native-providers.json").read_text())
    classes = plan["classes"]
    for node in ast.parse((SOURCE / "_lazy_registry_static.py").read_text()).body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            try:
                mapping = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                continue
            if isinstance(mapping, dict):
                for name, module in mapping.items():
                    if (
                        isinstance(name, str)
                        and isinstance(module, str)
                        and (SOURCE / (module + ".py")).exists()
                    ):
                        classes.setdefault(name, {"module": module, "class": name})
    configs = {
        name: config
        for name, config in plan["configs"].items()
        if config.get("type") not in {"secret", "api_key", "SpecialTool"}
    }
    missing = {
        name: cfg.get("type") for name, cfg in configs.items() if cfg.get("type") not in classes
    }
    if missing:
        raise ValueError(f"Unresolved referenced provider classes: {missing}")
    public = {name: native_name(name) for name in configs}
    if len(set(public.values())) != len(public):
        raise ValueError("Native function names collide")
    pending = {classes[cfg["type"]]["module"] for cfg in configs.values()}
    pending.update({"base_tool", "exceptions", "utils", "http_utils"})
    custom = {"tool_registry", "logging_config"}
    copied = set()
    provenance = []
    DEST.mkdir(parents=True, exist_ok=True)
    while pending:
        module = pending.pop()
        if module in copied or module in custom:
            continue
        path = SOURCE.joinpath(*module.split(".")).with_suffix(".py")
        if not path.exists():
            raise ValueError(f"Missing adapter dependency: {module}")
        original = path.read_bytes()
        text = original.decode("utf-8")
        tree = ast.parse(text)
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                dep = node.module
                if node.level:
                    package = module.split(".")[:-1]
                    prefix = package[: len(package) - node.level + 1]
                    dep = ".".join([*prefix, dep])
                elif dep.startswith("tooluniverse."):
                    dep = dep.removeprefix("tooluniverse.")
                else:
                    continue
                pending.add(dep)
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("tooluniverse."):
                        pending.add(alias.name.removeprefix("tooluniverse."))
        # Redirect both absolute imports and importlib.resources package names.
        text = re.sub(r"\btooluniverse\.", PREFIX + ".", text)
        text = text.replace('"tooluniverse"', '"' + PREFIX + '"')
        text = text.replace("'tooluniverse'", "'" + PREFIX + "'")
        text = text.replace("TOOLUNIVERSE_CACHE_DIR", "DNAKIT_CACHE_DIR")
        text = text.replace("ToolUniverse/", "DNAKit/")
        header = (
            f"# Adapted from mims-harvard/ToolUniverse, Apache-2.0, revision {REVISION}.\n"
            f"# Original: src/tooluniverse/{module.replace('.', '/')}.py\n"
            "# DNAKit modification: local provider namespace and DNAKit cache/user-agent names.\n"
            "# See LICENSE and UPSTREAM.json in this directory.\n"
        )
        target = DEST.joinpath(*module.split(".")).with_suffix(".py")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(header + text, encoding="utf-8")
        provenance.append(
            {
                "path": str(path.relative_to(SOURCE)),
                "sha256": hashlib.sha256(original).hexdigest(),
                "adapted_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
            }
        )
        copied.add(module)
    # Defaults and GraphQL schema resources are used by the actual adapters.
    data = DEST / "data"
    data.mkdir(exist_ok=True)
    for path in (SOURCE / "data").rglob("*"):
        if not path.is_file() or (
            "defaults" not in path.name and path.suffix not in {".graphql", ".gql"}
        ):
            continue
        target = data / path.relative_to(SOURCE / "data")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
        provenance.append(
            {
                "path": str(path.relative_to(SOURCE)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    for directory in [
        DEST,
        *[p for p in DEST.rglob("*") if p.is_dir() and "__pycache__" not in p.parts],
    ]:
        (directory / "__init__.py").write_text(
            '"""Pinned scientific adapter implementation; loaded on demand by DNAKit."""\n'
        )
    (DEST / "LICENSE").write_bytes((SOURCE.parents[1] / "LICENSE").read_bytes())
    (DEST / "UPSTREAM.json").write_text(
        json.dumps(
            {
                "source": "https://github.com/mims-harvard/ToolUniverse",
                "revision": REVISION,
                "license": "Apache-2.0",
                "attribution": "ToolUniverse contributors, mims-harvard",
                "scope": "Only adapters referenced by bundled skills and their implementation dependencies; no ToolUniverse orchestrator.",
                "files": sorted(provenance, key=lambda item: item["path"]),
            },
            indent=2,
        )
        + "\n"
    )
    entries = {}
    for original_name, config in sorted(configs.items()):
        name = public[original_name]
        entries[name] = {
            "name": name,
            "upstream_name": original_name,
            "config": config,
            **classes[config["type"]],
        }
    import copy

    extensions = json.loads((ROOT / "registry/native-extensions.json").read_text())
    for name, entry in entries.items():
        extension = extensions.get(name, {})
        entry["extension"] = extension
        schema = copy.deepcopy(entry["config"].get("parameter") or {"type": "object"})
        properties = schema.setdefault("properties", {})
        required = set(schema.get("required") or [])
        for alias, canonical in extension.get("aliases", {}).items():
            properties[alias] = copy.deepcopy(properties[canonical])
            properties[alias].pop("default", None)
            properties[alias].pop("required", None)
            properties[alias]["description"] = "DNAKit optional alias for " + canonical
            required.discard(canonical)
            properties[canonical].pop("required", None)
        properties.update(extension.get("parameters", {}))
        required.difference_update(extension.get("optional", []))
        schema["required"] = sorted(required)
        entry["input_schema"] = schema
    effects = json.loads((ROOT / "registry/provider-effects.json").read_text())
    for name, entry in entries.items():
        entry["effect"] = effects[name]
    (DEST.parent / "providers.json").write_text(
        json.dumps(entries, indent=2, ensure_ascii=False) + "\n"
    )
    write_functions(entries)
    (WORK / "provider-name-map.json").write_text(json.dumps(public, indent=2) + "\n")
    print(f"Ported {len(entries)} referenced API definitions and {len(copied)} dependency modules")


if __name__ == "__main__":
    main()
