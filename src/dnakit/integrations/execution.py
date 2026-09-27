"""Execute reviewed argument templates without shells, downloads or implicit installs."""

from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

from dnakit._progress import Progress
from dnakit.bioinformatics._common import ExecutionConfig, artifacts, run, run_to_file, workspace
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

from .catalog import describe_analysis


def _parameter(name: str, value: Any, schema: dict[str, Any]) -> str:
    kind = schema["type"]
    valid = False
    if kind == "integer":
        valid = isinstance(value, int) and not isinstance(value, bool)
    elif kind == "number":
        valid = isinstance(value, (int, float)) and not isinstance(value, bool)
        try:
            valid = valid and math.isfinite(value)
        except OverflowError:
            valid = False
    elif kind == "string":
        valid = isinstance(value, str) and bool(value) and not value.startswith("-")
        valid = valid and not any(c in value for c in "\x00\r\n") and len(value) <= 4096
    if not valid:
        raise ConfigurationError(f"Invalid {kind} parameter: {name}.")
    if "enum" in schema and value not in schema["enum"]:
        raise ConfigurationError(f"{name} must be one of {schema['enum']}.")
    if "minimum" in schema and value < schema["minimum"]:
        raise ConfigurationError(f"{name} is below {schema['minimum']}.")
    if "maximum" in schema and value > schema["maximum"]:
        raise ConfigurationError(f"{name} exceeds {schema['maximum']}.")
    if "pattern" in schema and re.fullmatch(schema["pattern"], value) is None:
        raise ConfigurationError(f"{name} does not match the required format.")
    return str(value)


def _input(name: str, value: str, schema: dict[str, Any]) -> Path:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ConfigurationError(f"Input {name} must be a non-empty path string.")
    path = Path(value).expanduser().resolve()
    if schema["type"] == "prefix":
        missing = [
            str(path) + suffix
            for suffix in schema["suffixes"]
            if not Path(str(path) + suffix).is_file()
        ]
        if missing:
            raise ConfigurationError(f"Incomplete {name} file set: {missing}.")
    elif schema["type"] == "directory":
        if not path.is_dir():
            raise ConfigurationError(f"Input directory {name} does not exist: {path}.")
        missing = [child for child in schema.get("contains", []) if not (path / child).is_file()]
        if missing:
            raise ConfigurationError(f"Input directory {name} is missing {missing}.")
    elif not path.is_file():
        raise ConfigurationError(f"Input file {name} does not exist: {path}.")
    for suffix in schema.get("required_suffixes", []):
        if not Path(str(path) + suffix).is_file():
            raise ConfigurationError(f"Input {name} requires an existing {suffix} index: {path}.")
    return path


def plan_analysis(
    analysis: str,
    inputs: dict[str, str],
    output_directory: str,
    *,
    parameters: dict[str, Any] | None = None,
    threads: int = 1,
    executable_paths: dict[str, str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate an analysis and return exact argument vectors without writing files.

    Input prefixes refer to complete file sets, e.g. BED/BIM/FAM. Paths are
    resolved before changing the working directory. Missing executables remain
    explicit in the plan; running requires a separately installed backend.
    """
    definition = describe_analysis(analysis, verbose=False)
    if not isinstance(inputs, dict) or set(inputs) != set(definition["inputs"]):
        raise ConfigurationError(f"inputs must contain exactly {sorted(definition['inputs'])}.")
    if parameters is not None and not isinstance(parameters, dict):
        raise ConfigurationError("parameters must be a dictionary.")
    parameters = dict(parameters or {})
    if set(parameters) - set(definition["parameters"]):
        raise ConfigurationError(
            f"Unknown parameters: {sorted(set(parameters) - set(definition['parameters']))}."
        )
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 1024:
        raise ConfigurationError("threads must be an integer in [1, 1024].")
    if not isinstance(output_directory, str) or not output_directory or "\x00" in output_directory:
        raise ConfigurationError("output_directory must be a non-empty path string.")
    directory = Path(output_directory).expanduser().resolve()
    if directory.exists():
        raise ConfigurationError(f"Output directory already exists: {directory}.")
    resolved_inputs = {
        name: _input(name, inputs[name], schema) for name, schema in definition["inputs"].items()
    }
    values = {name: str(path) for name, path in resolved_inputs.items()}
    values.update(out=str(directory), threads=str(threads))
    selected = {}
    for name, schema in definition["parameters"].items():
        value = parameters.get(name, schema.get("default"))
        values[name] = _parameter(name, value, schema)
        selected[name] = value
    programs = {step["program"] for step in definition["steps"]}
    if executable_paths is not None and not isinstance(executable_paths, dict):
        raise ConfigurationError("executable_paths must be a dictionary.")
    overrides = executable_paths or {}
    if set(overrides) - programs:
        raise ConfigurationError("executable_paths contains an unrelated program.")
    steps = []
    for step in definition["steps"]:
        program = overrides.get(step["program"], step["program"])
        if not isinstance(program, str) or not program or "\x00" in program:
            raise ConfigurationError("Executable paths must be non-empty strings.")
        resolved = shutil.which(program)
        steps.append(
            {
                "program": step["program"],
                "executable": resolved,
                "arguments": [argument.format_map(values) for argument in step["arguments"]],
                "stdout": str(directory / step["stdout"]) if step.get("stdout") else None,
                "expected_outputs": list(step["expected_outputs"]),
                "allow_empty_outputs": list(step.get("allow_empty_outputs", [])),
            }
        )
    return {
        "analysis": analysis,
        "project_id": definition.get("project_id"),
        "status": "planned",
        "output_directory": str(directory),
        "inputs": {key: str(value) for key, value in resolved_inputs.items()},
        "parameters": selected,
        "threads": threads,
        "steps": steps,
        "sources": definition["sources"],
        "installation": definition["installation"],
        "validation": definition["validation"],
    }


def run_analysis(
    analysis: str,
    inputs: dict[str, str],
    output_directory: str,
    *,
    parameters: dict[str, Any] | None = None,
    threads: int = 1,
    executable_paths: dict[str, str] | None = None,
    timeout_seconds: float = 3600.0,
    max_artifact_bytes: int = 1_000_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run an explicitly selected scientific backend in a fresh output directory.

    No packages, databases, model weights or credentials are installed or fetched.
    Results include original sources, exact commands and persistent output paths.
    Process success verifies produced artifacts, not scientific validity.
    """
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or not 0 < timeout_seconds <= 86_400
    ):
        raise ConfigurationError("timeout_seconds must be a finite number in (0, 86400].")
    config = ExecutionConfig(
        threads=threads, timeout_seconds=timeout_seconds, max_artifact_bytes=max_artifact_bytes
    )
    plan = plan_analysis(
        analysis,
        inputs,
        output_directory,
        parameters=parameters,
        threads=threads,
        executable_paths=executable_paths,
        verbose=False,
    )
    missing = [step["program"] for step in plan["steps"] if not step["executable"]]
    if missing:
        raise BackendUnavailableError(
            f"Install external programs before running: {', '.join(missing)}.",
            hint=plan["installation"],
        )
    directory = workspace(plan["output_directory"])
    manifest = directory / "dnakit-run.json"
    plan["status"] = "running"
    manifest.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    logs = []
    try:
        with Progress(disable=not verbose) as progress:
            task = progress.add_task(analysis, total=len(plan["steps"]))
            for step in plan["steps"]:
                if step["stdout"]:
                    log = run_to_file(
                        step["executable"], step["arguments"], Path(step["stdout"]), config
                    )
                else:
                    monitored = [
                        directory / path
                        for path in step["expected_outputs"]
                        if not any(c in path for c in "*?[")
                    ]
                    log = run(
                        step["executable"], step["arguments"], directory, config, outputs=monitored
                    )
                logs.append(log)
                for pattern in step["expected_outputs"]:
                    matches = [path for path in directory.glob(pattern) if path.is_file()]
                    if not matches:
                        raise BackendExecutionError(
                            f"{step['program']} produced no expected {pattern} artifact."
                        )
                    for path in matches:
                        if (
                            not path.resolve().is_relative_to(directory)
                            or path.stat().st_size > max_artifact_bytes
                        ):
                            raise BackendExecutionError(
                                f"Invalid or oversized artifact: {path.name}."
                            )
                        if path.stat().st_size == 0 and pattern not in step["allow_empty_outputs"]:
                            raise BackendExecutionError(
                                f"Expected non-empty artifact: {path.name}."
                            )
                progress.advance(task)
        produced = artifacts(directory)
        for name, path in produced.items():
            if path == manifest:
                continue
            if (
                not path.resolve().is_relative_to(directory)
                or path.stat().st_size > max_artifact_bytes
            ):
                raise BackendExecutionError(f"Invalid or oversized artifact: {name}.")
        plan["status"] = "completed"
        plan["artifacts"] = {name: str(path) for name, path in produced.items()}
        plan["log"] = "\n".join(logs)
    except Exception as exc:
        plan.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        if isinstance(exc, BackendExecutionError):
            plan["backend_error"] = dict(exc.context)
        raise
    finally:
        manifest.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return plan
