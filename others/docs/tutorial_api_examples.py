"""Source-backed call examples for functional tutorials; never execute providers."""

from __future__ import annotations

import ast
import inspect
import json
import pprint
from pathlib import Path
from typing import Literal, get_args, get_origin, get_type_hints

from jsonschema import Draft7Validator
from jsonschema.exceptions import ValidationError

OVERRIDES = json.loads(Path(__file__).with_suffix(".json").read_text(encoding="utf-8"))


def validate_arguments(record: dict, arguments: dict) -> None:
    """Check examples against both the Python signature and registered schema."""
    inspect.signature(record["fn"]).bind(**arguments)
    entry = record["entry"]
    if entry:
        Draft7Validator(entry["input_schema"]).validate(arguments)
    else:
        hints = get_type_hints(inspect.unwrap(record["fn"]))
        for name, value in arguments.items():
            annotation = hints.get(name)
            if get_origin(annotation) is Literal and value not in get_args(annotation):
                raise ValueError(f"{record['name']}: invalid {name}={value!r}")


def arguments_for(record: dict) -> dict:
    name = record["name"]
    if name in OVERRIDES:
        arguments = OVERRIDES[name]
        validate_arguments(record, arguments)
        return arguments
    candidates = (record["entry"] or {}).get("config", {}).get("test_examples", [])
    errors = []
    for arguments in candidates:
        if not isinstance(arguments, dict):
            continue
        try:
            validate_arguments(record, arguments)
        except (TypeError, ValueError, ValidationError) as exc:
            # Invalid upstream examples must never appear as valid Python calls.
            errors.append(str(exc))
        else:
            return arguments
    raise ValueError(f"No valid documentation example for {name}: {errors}")


def render_example(record: dict) -> str:
    arguments = {**arguments_for(record), "verbose": False}
    lines = ["import dnakit", ""]
    call_args = []
    for name, value in arguments.items():
        literal = pprint.pformat(value, width=84, sort_dicts=False)
        if "\n" in literal or len(literal) > 78:
            variable = f"example_{name}"
            lines += [f"{variable} = {literal}", ""]
            literal = variable
        call_args.append(f"    {name}={literal},")
    lines += [f"result = {record['name']}(", *call_args, ")"]
    if record["name"] in {"dnakit.workflows.run", "dnakit.workflows.run_support"}:
        lines += [
            'if result["returncode"] != 0:',
            '    raise RuntimeError(result["output"])',
            'print(result["output"])',
        ]
    else:
        lines += ["print(result)"]
    code = "\n".join(lines)
    ast.parse(code)
    return code
