"""Resolve verified resources in both source checkouts and installed wheels."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Any

from dnakit.backends.external import execute_bounded_command
from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.skills.catalog import _locations, _resource_bytes, get_skill


def scripts(skill_id: str, *, verbose: bool = True) -> list[str]:
    """List bundled Python, R and shell resources; consult their usage before execution."""
    entry = get_skill(skill_id)
    return [
        item["path"]
        for item in entry["files"]
        if Path(item["path"]).suffix.lower() in {".py", ".r", ".sh"}
        and not any(part.startswith("test") for part in Path(item["path"]).parts)
    ]


def resource_path(skill_id: str, resource: str, *, verbose: bool = True) -> Path:
    """Locate a checksum-verified bundled file, including resources of dependent skills."""
    entry = get_skill(skill_id)
    _resource_bytes(entry, resource)
    root, _ = _locations()
    return Path(root / entry["id"] / resource)


def _command(
    skill_id: str, script: str, arguments: list[str] | None
) -> tuple[dict[str, Any], list[str]]:
    entry = get_skill(skill_id)
    # Verify the whole bundle: scripts may import helpers or read local assets.
    for item in entry["files"]:
        _resource_bytes(entry, item["path"])
    path = resource_path(entry["id"], script)
    interpreter = {
        ".py": sys.executable,
        ".r": shutil.which("Rscript"),
        ".sh": shutil.which("bash"),
    }.get(path.suffix.lower())
    if path.suffix.lower() not in {".py", ".r", ".sh"}:
        raise ConfigurationError("Only bundled Python, R and shell scripts can be run.")
    if not interpreter:
        raise BackendUnavailableError(f"Interpreter for {path.suffix} is unavailable.")
    if arguments is not None and (
        not isinstance(arguments, list)
        or any(not isinstance(arg, str) or "\0" in arg for arg in arguments)
    ):
        raise ConfigurationError("arguments must be a list of strings without NUL bytes.")
    return entry, [interpreter, str(path), *(arguments or [])]


def run(
    skill_id: str,
    script: str,
    *,
    arguments: list[str] | None = None,
    workdir: str | Path | None = None,
    timeout_seconds: float = 300.0,
    max_output_bytes: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a bundled workflow with its original CLI options and true scientific backends.

    The caller chooses an explicit script and working directory. Python uses the
    current environment; R, shell tools, trained models and data retain their own
    requirements. This executes code and may write analysis outputs or access APIs.
    Nonzero exit status is returned, not interpreted as scientific success.
    """
    entry, command = _command(skill_id, script, arguments)
    result = execute_bounded_command(
        command[0],
        command[1:],
        backend_id="skill:" + entry["id"],
        cwd=workdir or Path.cwd(),
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        verbose=verbose,
    )
    return {
        "skill": entry["id"],
        "script": script,
        "returncode": result.return_code,
        "output": result.output,
        "upstream_url": entry["upstream_url"],
        "upstream_revision": entry["upstream_revision"],
        "license": entry["license"],
    }


def _support_command(name: str, arguments: list[str]) -> list[str]:
    import hashlib
    import json

    allowed = {
        "clinical-variant-reporter",
        "clinical-trial-finder",
        "omics-target-evidence-mapper",
        "gwas-lookup",
        "fine-mapping",
    }
    if name not in allowed:
        raise ConfigurationError("Unknown bundled support workflow.")
    root = Path(__file__).parent / "support"
    manifest = json.loads((root / "UPSTREAM.json").read_text())
    prefix = name.replace("-", "_")
    for entry in manifest["files"]:
        if not entry["path"].startswith(prefix + "/"):
            continue
        path = root / entry["path"]
        if (
            path.is_symlink()
            or not path.resolve().is_relative_to(root.resolve())
            or hashlib.sha256(path.read_bytes()).hexdigest() != entry["adapted_sha256"]
        ):
            raise ConfigurationError("Support workflow checksum mismatch.")
    return [sys.executable, str(root / prefix / (prefix + ".py")), *arguments]


def run_support(
    name: str,
    *,
    arguments: list[str] | None = None,
    workdir: str | Path | None = None,
    timeout_seconds: float = 300,
    max_output_bytes: int = 1_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a licensed support stage of a composite skill with its documented CLI options.

    Names: clinical-variant-reporter, clinical-trial-finder, omics-target-evidence-mapper,
    gwas-lookup, fine-mapping.
    Source implementations are bundled in DNAKit; scientific dependencies remain external.
    """
    command = _support_command(name, arguments or [])
    result = execute_bounded_command(
        command[0],
        command[1:],
        backend_id="workflow:" + name,
        cwd=workdir or Path.cwd(),
        timeout_seconds=timeout_seconds,
        max_output_bytes=max_output_bytes,
        verbose=verbose,
    )
    return {"name": name, "returncode": result.return_code, "output": result.output}
