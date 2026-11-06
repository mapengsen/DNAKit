"""Local, deterministic workflow-security functions.

Every function works on workflow descriptions, on files named by the caller or on the local Python
environment.  None of them executes a user command.  Operations that would need a job scheduler, a cloud
account or a running process (``benchmark``, ``cancel``, ``clear``, ``monitor``, ``retry``,
``submit_cloud_job``, ``submit_hpc_job``) return ``status: "not_implemented"``.

A *workflow specification* is a mapping ``{"name": ..., "steps": [{"id", "command", "depends_on", "inputs",
"outputs", "env", "cpus", "memory_mb", "runtime_s", "disk_mb"}, ...]}``; ``steps`` may also be a mapping
``id -> step``.  It can be given as a mapping, as JSON/YAML text or as the path of a JSON/YAML file.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import math
import os
import platform
import random
import re
import shutil
import sys
import tempfile
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit._metadata_pseudonymize import pseudonymize_envelope
from dnakit._version import __version__
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.workflow_security import _impl

_MODULE = "workflow_security"


def _opts(options: Mapping[str, object] | None) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError("options must be a mapping.", code="INVALID_OPTION")
    return dict(options or {})


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _timestamp(options: Mapping[str, object]) -> str:
    value = options.get("timestamp")
    if value is None:
        return _now()
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            "options.timestamp must be a non-empty string.", code="INVALID_OPTION"
        )
    return value


def _analysis(
    input: object, *, security: bool = True
) -> tuple[dict[str, object], dict[str, object], str]:
    document, source = _impl.load_document(input)
    return document, _impl.analyze_workflow(document, security=security), source


def _names(input: object, what: str) -> list[str]:
    """Names from a string (comma/whitespace separated) or a sequence of strings."""

    if isinstance(input, str):
        names = [item for item in re.split(r"[,\s]+", input) if item]
    elif isinstance(input, Sequence) and not isinstance(input, (bytes, bytearray)):
        if any(not isinstance(item, str) for item in input):
            raise InputFormatError(f"{what} must be strings.", code="INVALID_INPUT")
        names = [item.strip() for item in input if item.strip()]
    else:
        raise InputFormatError(
            f"{what} must be a string or a list of strings.", code="INVALID_INPUT"
        )
    if not names:
        raise InputFormatError(f"No {what} were given.", code="INVALID_INPUT")
    return names


def _version_tuple(version: str) -> tuple[int, ...]:
    parts: list[int] = []
    for component in version.split("."):
        match = re.match(r"\d+", component)
        if not match:
            break
        parts.append(int(match.group(0)))
    return tuple(parts)


def _compare(left: tuple[int, ...], operator: str, right: tuple[int, ...]) -> bool:
    width = max(len(left), len(right))
    a, b = left + (0,) * (width - len(left)), right + (0,) * (width - len(right))
    return {"==": a == b, "!=": a != b, ">=": a >= b, "<=": a <= b, ">": a > b, "<": a < b}[
        operator
    ]


_REQUIREMENT = re.compile(
    r"^\s*([A-Za-z0-9][A-Za-z0-9._-]*)\s*(?:(==|!=|>=|<=|>|<)\s*([0-9][A-Za-z0-9.]*))?\s*$"
)


def benchmark(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: benchmarking would have to execute the workflow or command being measured."""

    return candidate_not_implemented(
        _MODULE,
        "benchmark",
        input,
        options=options,
        requires="an execution backend allowed to run the workflow under test",
        message="Benchmarking needs to execute the workload; workflow_security never executes user commands.",
    )


def cache(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute the content-addressed cache key of a step description and report whether a cache entry exists.

    ``input`` is any JSON-serializable step description (mapping, list or string).  The key is the SHA-256 of the
    canonical JSON (sorted keys, compact separators) of ``{"step": input, "files": {path: sha256}, "salt": salt}``.
    Options: ``input_files`` (list of paths whose content digests are part of the key, default none),
    ``salt`` (string, default ""), ``cache_dir`` (directory; when given, ``hit`` is true if ``<cache_dir>/<key>``
    or ``<cache_dir>/<key>.json`` exists, otherwise ``hit`` is null).  Nothing is written.
    """

    settings = _opts(options)
    files: dict[str, str] = {}
    for name in settings.get("input_files", []):  # type: ignore[union-attr]
        path = Path(str(name))
        if not path.is_file():
            raise InputFormatError(f"Input file does not exist: {name}", code="FILE_NOT_FOUND")
        files[str(name)] = _impl.sha256_file(path)
    key = _impl.sha256_text(
        _impl.canonical_json({"step": input, "files": files, "salt": str(settings.get("salt", ""))})
    )
    hit = None
    cache_dir = _impl.option_path(settings, "cache_dir")
    if cache_dir is not None:
        hit = (cache_dir / key).exists() or (cache_dir / f"{key}.json").exists()
    return candidate_ok(
        _MODULE,
        "cache",
        input,
        {"cache_key": key, "hit": hit, "input_files": files},
        options=options,
    )


def cancel(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: cancelling needs a running job or scheduler handle, which this library does not keep."""

    return candidate_not_implemented(
        _MODULE,
        "cancel",
        input,
        options=options,
        requires="a job scheduler or process supervisor that owns the running job",
        message="There is no job registry in this process, so there is nothing that cancel could stop.",
    )


def check_dependencies(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Check that Python packages are installed (and satisfy an optional version constraint) without importing them.

    ``input``: a comma/newline separated string, a list of requirement strings such as ``"numpy>=1.24"``
    (operators ``== != >= <= > <``, one clause), or a mapping ``name -> constraint``.  Versions are compared as
    dotted integer components (suffixes such as ``rc1`` are ignored).  A missing package is reported, not raised.
    The result has ``all_satisfied`` and one row per requirement with ``installed``, ``version`` and ``satisfied``.
    """

    _opts(options)
    if isinstance(input, Mapping):
        requirements = [f"{name}{spec or ''}" for name, spec in input.items()]
    else:
        requirements = _names(input, "requirements")
    rows = []
    for requirement in requirements:
        match = _REQUIREMENT.match(requirement)
        if match is None:
            raise InputFormatError(
                f"Unsupported requirement {requirement!r}.",
                code="INVALID_REQUIREMENT",
                hint="Use 'name' or 'name<op>version' with <op> in == != >= <= > <.",
            )
        name, operator, wanted = match.groups()
        version: str | None
        try:
            version = importlib.metadata.version(name)
            installed = True
        except importlib.metadata.PackageNotFoundError:
            version = None
            installed = importlib.util.find_spec(name.replace("-", "_").split(".")[0]) is not None
        satisfied = installed
        if installed and operator:
            satisfied = version is not None and _compare(
                _version_tuple(version), operator, _version_tuple(wanted)
            )
        rows.append(
            {
                "requirement": requirement.strip(),
                "name": name,
                "installed": installed,
                "version": version,
                "satisfied": satisfied,
            }
        )
    result = {"all_satisfied": all(row["satisfied"] for row in rows), "requirements": rows}
    return candidate_ok(_MODULE, "check_dependencies", input, result, options=options)


def check_external_tools(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Look up external programs on PATH with ``shutil.which`` (the programs are not executed).

    ``input``: a program name, a comma-separated string or list of names/paths, or a workflow specification (the
    program of each step command, i.e. its first argv token, is checked).  The result has ``all_available`` and one
    row per distinct tool with ``available`` and the resolved ``path``.
    """

    _opts(options)
    if isinstance(input, Mapping) or (
        isinstance(input, str) and ("\n" in input or input.strip()[:1] in {"{", "["})
    ):
        document, _ = _impl.load_document(input)
        steps, _ = _impl.normalize_steps(document)
        tools = []
        for step in steps:
            tokens = _impl.command_tokens(step.get("command"))
            if tokens:
                tools.append(tokens[0])
        if not tools:
            raise InputFormatError(
                "The workflow has no step commands to check.", code="INVALID_INPUT"
            )
    else:
        tools = _names(input, "tool names")
    rows = []
    for tool in dict.fromkeys(tools):
        path = shutil.which(tool)
        rows.append({"tool": tool, "available": path is not None, "path": path})
    result = {"all_available": all(row["available"] for row in rows), "tools": rows}
    return candidate_ok(_MODULE, "check_external_tools", input, result, options=options)


def check_numerical_reproducibility(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Compare repeated numeric results elementwise and report whether they agree within a tolerance.

    ``input``: a list of at least two runs (each a number or a flat/nested list of numbers) or a mapping
    ``{"runs": [...]}``.  Every run is compared with the first.  Two values ``a`` and ``b`` agree when
    ``|a-b| <= max(rtol * max(|a|, |b|), atol)``.  Options: ``rtol`` (default 1e-9), ``atol`` (default 0.0),
    ``equal_nan`` (default false).  Result: ``reproducible``, ``bitwise_identical``, ``max_abs_diff``,
    ``max_rel_diff`` (relative to ``max(|a|, |b|)``) and one row per compared run.
    """

    settings = _opts(options)
    rtol, atol = float(settings.get("rtol", 1e-9)), float(settings.get("atol", 0.0))
    equal_nan = _impl.option_bool(settings, "equal_nan", False)
    if rtol < 0 or atol < 0:
        raise ConfigurationError("rtol and atol must be non-negative.", code="INVALID_OPTION")
    runs = input.get("runs") if isinstance(input, Mapping) else input
    if not isinstance(runs, Sequence) or isinstance(runs, (str, bytes)) or len(runs) < 2:
        raise InputFormatError("Provide at least two runs to compare.", code="INVALID_INPUT")

    def flatten(value: object) -> list[float]:
        if isinstance(value, bool):
            raise InputFormatError("Runs must contain numbers, not booleans.", code="INVALID_INPUT")
        if isinstance(value, (int, float)):
            return [float(value)]
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            return [number for item in value for number in flatten(item)]
        raise InputFormatError(
            f"Runs must contain only numbers, got {type(value).__name__}.", code="INVALID_INPUT"
        )

    flat = [flatten(run) for run in runs]
    rows = []
    for index, other in enumerate(flat[1:], start=1):
        base = flat[0]
        if len(other) != len(base):
            rows.append(
                {
                    "run": index,
                    "length_match": False,
                    "agrees": False,
                    "max_abs_diff": None,
                    "max_rel_diff": None,
                    "bitwise_identical": False,
                }
            )
            continue
        max_abs = max_rel = 0.0
        agrees = identical = True
        for a, b in zip(base, other, strict=False):
            if math.isnan(a) or math.isnan(b):
                same = equal_nan and math.isnan(a) and math.isnan(b)
                agrees &= same
                identical &= math.isnan(a) and math.isnan(b)
                continue
            difference = abs(a - b)
            scale = max(abs(a), abs(b))
            identical &= a == b
            agrees &= difference <= max(rtol * scale, atol)
            max_abs = max(max_abs, difference)
            max_rel = max(max_rel, difference / scale if scale else 0.0)
        rows.append(
            {
                "run": index,
                "length_match": True,
                "agrees": bool(agrees),
                "max_abs_diff": max_abs,
                "max_rel_diff": max_rel,
                "bitwise_identical": bool(identical),
            }
        )
    finite = [row["max_abs_diff"] for row in rows if row["max_abs_diff"] is not None]
    result = {
        "reproducible": all(row["agrees"] for row in rows),
        "bitwise_identical": all(row["bitwise_identical"] for row in rows),
        "max_abs_diff": max(finite) if finite else None,
        "max_rel_diff": max(
            (row["max_rel_diff"] for row in rows if row["max_rel_diff"] is not None), default=None
        ),
        "rtol": rtol,
        "atol": atol,
        "runs": rows,
    }
    return candidate_ok(_MODULE, "check_numerical_reproducibility", input, result, options=options)


def checkpoint(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a checksummed checkpoint document of a JSON-serializable workflow state, optionally writing it.

    The document is ``{"format": "dnakit-checkpoint-v1", "step", "created_at", "state", "state_sha256"}`` where
    ``state_sha256`` is the SHA-256 of the canonical JSON of ``state``.  Options: ``step`` (label, default
    "checkpoint"), ``timestamp`` (ISO string, default now), ``path`` (file to write atomically; nothing is written
    without it), ``overwrite`` (default false).  Use ``resume`` to verify and read a checkpoint back.
    """

    settings = _opts(options)
    state_json = _impl.canonical_json(input)
    document = {
        "format": _impl.CHECKPOINT_FORMAT,
        "step": str(settings.get("step", "checkpoint")),
        "created_at": _timestamp(settings),
        "state": input,
        "state_sha256": _impl.sha256_text(state_json),
    }
    target = _impl.option_path(settings, "path")
    written = None
    if target is not None:
        _impl.write_atomic(
            target,
            (_impl.canonical_json(document) + "\n").encode("utf-8"),
            overwrite=_impl.option_bool(settings, "overwrite", False),
        )
        written = str(target)
    return candidate_ok(
        _MODULE, "checkpoint", input, {"checkpoint": document, "path": written}, options=options
    )


def cleanup_temp_workspace(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Inventory (and, only when asked, delete) a temporary workspace directory.

    ``input`` is the directory path.  Safety rules: it must be a real directory (not a symlink) strictly inside
    ``allowed_root`` (default: the system temporary directory) and its name must start with ``prefix`` (default
    "dnakit-").  Options: ``dry_run`` (default true: only count files and bytes), ``allowed_root``, ``prefix``.
    Deletion is ``shutil.rmtree`` of that directory; symlinks inside are removed, never followed.
    """

    settings = _opts(options)
    if not isinstance(input, (str, os.PathLike)):
        raise InputFormatError("input must be the workspace directory path.", code="INVALID_INPUT")
    path = Path(input)
    if path.is_symlink():
        raise ConfigurationError("Refusing to clean a symbolic link.", code="UNSAFE_CLEANUP_TARGET")
    if not path.is_dir():
        raise InputFormatError(f"Not a directory: {path}", code="DIRECTORY_NOT_FOUND")
    root = Path(str(settings.get("allowed_root") or tempfile.gettempdir())).resolve()
    prefix = str(settings.get("prefix", "dnakit-"))
    resolved = path.resolve()
    if resolved == root or root not in resolved.parents:
        raise ConfigurationError(
            f"{resolved} is not strictly inside the allowed root {root}.",
            code="UNSAFE_CLEANUP_TARGET",
            hint="Set options.allowed_root to the directory that owns the workspace.",
        )
    if not resolved.name.startswith(prefix):
        raise ConfigurationError(
            f"Workspace name {resolved.name!r} does not start with {prefix!r}.",
            code="UNSAFE_CLEANUP_TARGET",
            hint="Set options.prefix to the prefix used when the workspace was created.",
        )
    files = size = 0
    for current, _, names in os.walk(resolved, followlinks=False):
        for name in names:
            files += 1
            candidate = Path(current, name)
            if not candidate.is_symlink():
                size += candidate.stat().st_size
    dry_run = _impl.option_bool(settings, "dry_run", True)
    if not dry_run:
        shutil.rmtree(resolved)
    result = {
        "path": str(resolved),
        "files": files,
        "bytes": size,
        "dry_run": dry_run,
        "deleted": not dry_run,
    }
    return candidate_ok(_MODULE, "cleanup_temp_workspace", input, result, options=options)


def clear(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: the name does not say what to clear and this library keeps no workflow state to reset."""

    return candidate_not_implemented(
        _MODULE,
        "clear",
        input,
        options=options,
        requires="a concrete state store (cache directory, job queue) to clear",
        message="clear has no defined target; use cleanup_temp_workspace for a temporary directory.",
    )


def create(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a canonical, validated workflow specification from a list of steps or a partial specification.

    ``input``: a mapping with ``steps`` (and optional ``name``), or a list whose items are step mappings or plain
    command strings.  Steps without an id are named ``step_1``, ``step_2``, ... by position.  The specification is
    validated like ``validate`` (structure, dependencies, cycles, secrets, dangerous commands, paths); any error
    raises ``InputFormatError``.  Result: ``format``, ``name``, ``steps`` (each with a ``depends_on`` list),
    ``order`` (topological order) and ``sha256`` of the canonical JSON.  Option: ``name`` (default "workflow").
    """

    settings = _opts(options)
    if isinstance(input, Mapping):
        raw_steps = input.get("steps")
        name = str(input.get("name", settings.get("name", "workflow")))
    else:
        raw_steps, name = input, str(settings.get("name", "workflow"))
    if isinstance(raw_steps, Mapping):
        document: dict[str, object] = {"name": name, "steps": raw_steps}
    elif isinstance(raw_steps, Sequence) and not isinstance(raw_steps, (str, bytes)):
        items = []
        for index, step in enumerate(raw_steps, start=1):
            if isinstance(step, str):
                step = {"command": step}
            if isinstance(step, Mapping):
                step = {"id": f"step_{index}", **step}
            items.append(step)
        document = {"name": name, "steps": items}
    else:
        raise InputFormatError(
            "input must define a list or mapping of steps.", code="INVALID_INPUT"
        )
    analysis = _impl.analyze_workflow(document)
    _impl.raise_if_invalid(analysis, what="creation")
    spec = {
        "format": _impl.SPEC_FORMAT,
        "name": name,
        "steps": analysis["steps"],
        "order": analysis["order"],
    }
    spec["sha256"] = _impl.sha256_text(_impl.canonical_json({k: v for k, v in spec.items()}))
    return candidate_ok(
        _MODULE,
        "create",
        input,
        spec,
        options=options,
        warnings=[str(item["message"]) for item in analysis["warnings"]],
    )  # type: ignore[union-attr]


def _fernet(options: Mapping[str, object], *, allow_generate: bool) -> tuple[object, str | None]:
    from cryptography.fernet import Fernet

    key = options.get("key")
    env_name = options.get("key_env")
    generated = None
    if key is None and env_name is not None:
        key = os.environ.get(str(env_name))
        if not key:
            raise ConfigurationError(
                f"Environment variable {env_name!r} is not set.", code="KEY_MISSING"
            )
    if key is None and allow_generate and _impl.option_bool(options, "generate_key", False):
        generated = Fernet.generate_key().decode("ascii")
        key = generated
    if key is None:
        raise ConfigurationError(
            "An encryption key is required.",
            code="KEY_MISSING",
            hint="Pass options.key (urlsafe base64, 32 bytes), options.key_env, or options.generate_key=true when encrypting.",
        )
    try:
        return Fernet(key if isinstance(key, bytes) else str(key).encode("ascii")), generated
    except (ValueError, TypeError, UnicodeEncodeError) as exc:
        raise ConfigurationError(
            "The key is not a valid Fernet key (32 url-safe base64 bytes).", code="INVALID_KEY"
        ) from exc


def _crypt_file(operation: str, input: object, options: Mapping[str, object] | None) -> object:
    settings = _opts(options)
    try:
        import cryptography.fernet  # noqa: F401
    except ImportError:
        return candidate_not_implemented(
            _MODULE,
            operation,
            input,
            options=options,
            requires="pip install cryptography",
            message="File encryption uses cryptography.fernet (AES-128-CBC + HMAC-SHA256), which is not installed.",
        )
    if not isinstance(input, (str, os.PathLike)):
        raise InputFormatError("input must be the path of the file.", code="INVALID_INPUT")
    source = Path(input)
    if not source.is_file():
        raise InputFormatError(f"File not found: {source}", code="FILE_NOT_FOUND")
    limit = int(settings.get("max_bytes", 256 * 1024 * 1024))
    if not 1 <= limit <= 2 * 1024**3 or source.stat().st_size > limit:
        raise ConfigurationError(
            f"File exceeds max_bytes ({limit}); Fernet processes the whole file in memory.",
            code="FILE_TOO_LARGE",
        )
    encrypting = operation == "encrypt_genomic_file"
    output = _impl.option_path(settings, "output")
    if output is None:
        if encrypting:
            output = source.with_name(source.name + ".enc")
        else:
            output = (
                source.with_suffix("")
                if source.suffix == ".enc"
                else source.with_name(source.name + ".dec")
            )
    if output.resolve() == source.resolve():
        raise ConfigurationError(
            "The output path must differ from the input path.", code="OUTPUT_EXISTS"
        )
    cipher, generated = _fernet(settings, allow_generate=encrypting)
    data = source.read_bytes()
    if encrypting:
        payload = cipher.encrypt(data)  # type: ignore[attr-defined]
    else:
        from cryptography.fernet import InvalidToken

        try:
            payload = cipher.decrypt(data)  # type: ignore[attr-defined]
        except InvalidToken as exc:
            raise InputFormatError(
                "Decryption failed: wrong key or the file is not a Fernet token.",
                code="DECRYPTION_FAILED",
            ) from exc
    _impl.write_atomic(output, payload, overwrite=_impl.option_bool(settings, "overwrite", False))
    result: dict[str, object] = {
        "output": str(output),
        "bytes_in": len(data),
        "bytes_out": len(payload),
        "plaintext_sha256": hashlib.sha256(data if encrypting else payload).hexdigest(),
        "algorithm": "fernet-aes128-cbc-hmac-sha256",
    }
    warnings = []
    if generated is not None:
        result["generated_key"] = generated
        warnings.append(
            "A key was generated; store it safely, the data cannot be decrypted without it."
        )
    return candidate_ok(_MODULE, operation, input, result, options=options, warnings=warnings)


def decrypt_genomic_file(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Decrypt a file produced by ``encrypt_genomic_file`` (Fernet token) with the given key.

    ``input`` is the path of the encrypted file.  Options: ``key`` (urlsafe-base64 Fernet key) or ``key_env``
    (name of an environment variable holding it), ``output`` (default: input without ``.enc``, else input +
    ".dec"), ``overwrite`` (default false), ``max_bytes`` (default 256 MiB).  A wrong key raises
    ``InputFormatError``.  Requires the optional package ``cryptography``.
    """

    return _crypt_file("decrypt_genomic_file", input, options)


def encrypt_genomic_file(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Encrypt a file with authenticated symmetric encryption (Fernet: AES-128-CBC + HMAC-SHA256).

    ``input`` is the path of the plaintext file.  Options: ``key`` (urlsafe-base64 Fernet key) or ``key_env``
    (environment variable name) or ``generate_key`` (true: create a new key and return it as ``generated_key``),
    ``output`` (default input + ".enc"), ``overwrite`` (default false), ``max_bytes`` (default 256 MiB, the whole
    file is processed in memory).  The result holds sizes and the SHA-256 of the plaintext.  Requires the optional
    package ``cryptography``.
    """

    return _crypt_file("encrypt_genomic_file", input, options)


def estimate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate wall time and peak resources of a workflow specification from per-step resource declarations.

    Steps declare ``cpus``, ``memory_mb``, ``runtime_s``, ``disk_mb`` and ``depends_on``.  The estimate is an
    as-soon-as-possible schedule with unlimited workers: ``wall_time_s`` is the critical-path length, peaks are taken
    over steps running at the same time, ``total_cpu_seconds`` = sum(cpus * runtime_s), ``total_disk_mb`` = sum of
    ``disk_mb``.  Options: ``default_cpus`` (1), ``default_memory_mb`` (0), ``default_runtime_s`` (0),
    ``default_disk_mb`` (0) for steps that omit a value, ``max_cpus`` (adds a warning when exceeded).
    An invalid dependency graph raises ``InputFormatError``.
    """

    settings = _opts(options)
    _, analysis, _ = _analysis(input, security=False)
    _impl.raise_if_invalid(analysis, what="estimation")
    result = _impl.estimate_resources(
        analysis["steps"],
        analysis["order"],  # type: ignore[arg-type]
        default_cpus=float(settings.get("default_cpus", 1)),
        default_memory_mb=float(settings.get("default_memory_mb", 0)),
        default_runtime_s=float(settings.get("default_runtime_s", 0)),
        default_disk_mb=float(settings.get("default_disk_mb", 0)),
    )
    warnings = []
    if "max_cpus" in settings and result["peak_cpus"] > float(settings["max_cpus"]):  # type: ignore[operator]
        warnings.append(f"peak cpus {result['peak_cpus']} exceed max_cpus {settings['max_cpus']}")
    return candidate_ok(_MODULE, "estimate", input, result, options=options, warnings=warnings)


def monitor(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: monitoring needs a running job or scheduler to observe."""

    return candidate_not_implemented(
        _MODULE,
        "monitor",
        input,
        options=options,
        requires="a job scheduler or process supervisor exposing job status",
        message="There is no running job handle to monitor in this library.",
    )


def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a workflow specification or configuration from a mapping, JSON/YAML text or a JSON/YAML file.

    JSON is tried first, then the YAML safe loader (no tags, no code execution); files larger than 1 MB are
    refused.  The result holds the parsed ``document``, its ``source`` ("inline" or the file path), the SHA-256 of
    its canonical JSON and the number of top-level keys.
    """

    _opts(options)
    document, source = _impl.load_document(input)
    result = {
        "document": document,
        "source": source,
        "sha256": _impl.sha256_text(_impl.canonical_json(document)),
        "top_level_keys": sorted(document),
    }
    return candidate_ok(_MODULE, "read", input, result, options=options)


def record_provenance(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Create a provenance record of a processing step with SHA-256 digests of its input and output files.

    ``input`` is a mapping such as ``{"tool", "version", "parameters", "command", "inputs": [paths],
    "outputs": [paths]}``; ``inputs``/``outputs`` are hashed (``exists``, ``size``, ``sha256``; missing files are
    recorded with ``exists: false`` and a warning).  The record also holds the Python, platform and DNAKit
    versions, a timestamp and ``record_sha256`` (SHA-256 of the canonical JSON of everything else).  Text
    containing credential-looking strings is refused (``ConfigurationError``).  Options: ``timestamp`` (ISO string,
    default now), ``output`` (JSON file to write atomically; nothing is written without it), ``overwrite`` (false).
    """

    settings = _opts(options)
    if not isinstance(input, Mapping):
        raise InputFormatError("input must be a mapping describing the step.", code="INVALID_INPUT")
    for location, text in _impl._walk_strings(dict(input)):
        if _impl.scan_secrets(text):
            raise ConfigurationError(
                f"Credential-looking text found at {location}; remove it before recording provenance.",
                code="SECRET_IN_PROVENANCE",
            )
    warnings: list[str] = []

    def describe(paths: object) -> list[dict[str, object]]:
        values = (
            [paths]
            if isinstance(paths, str)
            else (list(paths) if isinstance(paths, Sequence) else [])
        )  # type: ignore[call-overload]
        if len(values) > 1000:
            raise ConfigurationError("At most 1000 files can be hashed.", code="TOO_MANY_FILES")
        rows = []
        for value in values:
            path = Path(str(value))
            if path.is_file():
                rows.append(
                    {
                        "path": str(value),
                        "exists": True,
                        "size": path.stat().st_size,
                        "sha256": _impl.sha256_file(path),
                    }
                )
            else:
                warnings.append(f"file not found: {value}")
                rows.append({"path": str(value), "exists": False, "size": None, "sha256": None})
        return rows

    record: dict[str, object] = {
        "format": "dnakit-provenance-v1",
        "created_at": _timestamp(settings),
        "step": {key: value for key, value in input.items() if key not in {"inputs", "outputs"}},
        "inputs": describe(input.get("inputs", [])),
        "outputs": describe(input.get("outputs", [])),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "dnakit": __version__,
        },
    }
    record["record_sha256"] = _impl.sha256_text(_impl.canonical_json(record))
    target = _impl.option_path(settings, "output")
    if target is not None:
        _impl.write_atomic(
            target,
            (_impl.canonical_json(record) + "\n").encode("utf-8"),
            overwrite=_impl.option_bool(settings, "overwrite", False),
        )
    return candidate_ok(
        _MODULE, "record_provenance", input, record, options=options, warnings=warnings
    )


def pseudonymize_metadata_fields(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Pseudonymise sensitive metadata fields (redact, drop or keyed-hash); key-name and pattern based, NOT de-identification.

    Input: a mapping, a list of mappings, or JSON text of those; anything else raises InputFormatError. Nested mappings
    and lists are searched. A key is sensitive when its name, lower-cased with everything except letters and digits
    removed (``Patient_ID`` -> ``patientid``), is in a default list (name, first/last/full name, patient, patient id/name,
    mrn, medical record number, ssn, social security number, email, phone, telephone, mobile, address, street, city, zip,
    postal code, dob, date of birth, birth date, national id, passport, insurance, account number, ip address, device id,
    photo, biometric and a few variants).
    Options: ``mode`` ('redact' default: value -> "[REDACTED]"; 'drop': key removed; 'hash': HMAC-SHA256 token),
    ``salt`` (required for 'hash', string of >= 8 characters, else ConfigurationError), ``digest_chars`` (hash token
    length, 8-64, default 16), ``keys`` (extra sensitive key names), ``keep`` (key names never treated as sensitive),
    ``scan_text`` (True: e-mail, phone and SSN-like substrings of every string value become "[REDACTED:email]" etc.).
    The input is not modified. The result has ``metadata``, ``mode``, ``n_redacted``, ``redacted_paths`` (key paths only,
    never values) and ``warnings``; ``input_summary`` has only ``record_count`` and ``n_keys``. This is NOT a guarantee
    of de-identification: dates, ages over 89, small geographic units and names in free text are not handled.
    """
    return pseudonymize_envelope(
        "workflow_security", "pseudonymize_metadata_fields", input, options
    )


def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarize a workflow specification (validity, execution order, findings, resource estimate) as Markdown or JSON.

    Built from the same checks as ``validate`` and ``estimate`` (dependency analysis, secret/command/path scan,
    as-soon-as-possible resource schedule).  Options: ``format`` ("markdown" default or "json").  The result has
    ``format``, ``valid``, ``report`` (text or mapping).
    """

    settings = _opts(options)
    output_format = str(settings.get("format", "markdown"))
    if output_format not in {"markdown", "json"}:
        raise ConfigurationError(
            "options.format must be 'markdown' or 'json'.", code="INVALID_OPTION"
        )
    document, analysis, source = _analysis(input)
    errors, order = analysis["errors"], analysis["order"]
    estimate_result = None
    if order is not None and not any(
        item["code"] == "INVALID_RESOURCE" for item in analysis["findings"]
    ):  # type: ignore[union-attr]
        estimate_result = _impl.estimate_resources(
            analysis["steps"],
            order,
            default_cpus=1,
            default_memory_mb=0,
            default_runtime_s=0,
            default_disk_mb=0,
        )  # type: ignore[arg-type]
    name = str(document.get("name", "workflow"))
    summary = {
        "name": name,
        "source": source,
        "valid": not errors,
        "step_count": len(analysis["steps"]),
        "order": order,
        "findings": analysis["findings"],
        "estimate": estimate_result,
    }  # type: ignore[arg-type]
    if output_format == "json":
        return candidate_ok(
            _MODULE,
            "report",
            input,
            {"format": "json", "valid": not errors, "report": summary},
            options=options,
        )
    lines = [
        f"# Workflow report: {name}",
        "",
        f"- source: {source}",
        f"- steps: {summary['step_count']}",
        f"- valid: {'yes' if not errors else 'no'}",
        "",
        "## Execution order",
        "",
    ]
    lines += [f"{index}. {item}" for index, item in enumerate(order or [], start=1)] or [
        "(no valid order)"
    ]
    lines += ["", "## Findings", ""]
    lines += [
        f"- [{item['severity']}] {item['code']}"
        + (f" ({item['step']})" if item["step"] else "")
        + f": {item['message']}"
        for item in analysis["findings"]
    ] or ["(none)"]  # type: ignore[union-attr]
    if estimate_result is not None:
        lines += [
            "",
            "## Resource estimate",
            "",
            f"- wall time: {estimate_result['wall_time_s']} s",
            f"- total cpu time: {estimate_result['total_cpu_seconds']} cpu-s",
            f"- peak cpus: {estimate_result['peak_cpus']}",
            f"- peak memory: {estimate_result['peak_memory_mb']} MB",
            f"- critical path: {' -> '.join(estimate_result['critical_path'])}",
        ]  # type: ignore[arg-type]
    return candidate_ok(
        _MODULE,
        "report",
        input,
        {"format": "markdown", "valid": not errors, "report": "\n".join(lines) + "\n"},
        options=options,
    )


def resume(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Resume from a verified checkpoint, or resume a DNAKit workflow run from its manifest.

    If ``input`` (mapping, JSON text or file path) is a ``dnakit-checkpoint-v1`` document, its ``state_sha256`` is
    re-computed; the result holds ``valid``, ``step`` and ``state`` and a mismatch raises ``InputFormatError``.
    Otherwise ``input`` must be the path of a DNAKit workflow file and ``dnakit.workflows.run_workflow(path,
    resume=True)`` is called (only built-in deterministic operations run, never commands); the result is its
    ``to_dict()``.
    """

    _opts(options)
    document: dict[str, object] | None = None
    inline = isinstance(input, str) and ("\n" in input or input.strip()[:1] in {"{", "["})
    if isinstance(input, Mapping):
        document = dict(input)
    elif inline or (isinstance(input, (str, os.PathLike)) and Path(input).is_file()):
        try:
            document, _ = _impl.load_document(input)
        except InputFormatError:
            document = None  # not parseable here: let run_workflow report the problem
    if document is not None and document.get("format") == _impl.CHECKPOINT_FORMAT:
        state = document.get("state")
        if _impl.sha256_text(_impl.canonical_json(state)) != document.get("state_sha256"):
            raise InputFormatError(
                "Checkpoint state does not match its recorded SHA-256.", code="CHECKPOINT_CORRUPT"
            )
        result = {
            "kind": "checkpoint",
            "valid": True,
            "step": document.get("step"),
            "created_at": document.get("created_at"),
            "state": state,
        }
        return candidate_ok(_MODULE, "resume", input, result, options=options)
    if inline or not isinstance(input, (str, os.PathLike)):
        raise InputFormatError(
            "input must be a checkpoint document or the path of a workflow file.",
            code="INVALID_INPUT",
        )
    from dnakit.workflows import run_workflow

    run_result = run_workflow(input, resume=True, verbose=False)
    return candidate_ok(
        _MODULE,
        "resume",
        input,
        {"kind": "workflow", **run_result.to_dict()},
        options=options,
        implementation="dnakit.workflows.run_workflow",
    )


def retry(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: retrying means executing a failed job again, and no user command is ever executed here."""

    return candidate_not_implemented(
        _MODULE,
        "retry",
        input,
        options=options,
        requires="an execution backend that owns the failed job",
        message="Use resume to continue a DNAKit workflow or checkpoint; re-running arbitrary jobs is not supported.",
    )


def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run a DNAKit workflow file with the built-in deterministic operations (no commands, no network).

    ``input`` is the path of a JSON/YAML workflow in the ``dnakit-workflow-v1`` schema (see
    ``dnakit.workflows.load_workflow``); it is executed by ``dnakit.workflows.run_workflow`` and writes only into
    the workflow's own ``output_dir``.  Option: ``dry_run`` (default false) validates and plans without writing.
    The result is the run result as a dictionary.
    """

    settings = _opts(options)
    if not isinstance(input, (str, os.PathLike)):
        raise InputFormatError("input must be the path of a workflow file.", code="INVALID_INPUT")
    from dnakit.workflows import run_workflow

    run_result = run_workflow(
        input, dry_run=_impl.option_bool(settings, "dry_run", False), verbose=False
    )
    return candidate_ok(
        _MODULE,
        "run",
        input,
        run_result.to_dict(),
        options=options,
        implementation="dnakit.workflows.run_workflow",
    )


def save_config(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Save a configuration mapping as canonical pretty JSON, refusing configurations that contain secrets.

    ``input`` is a JSON-serializable mapping.  Options: ``output`` (required file path), ``overwrite`` (default
    false), ``allow_secrets`` (default false: credential-looking text, e.g. AWS keys, private-key blocks or
    ``password=...`` values, raises ``ConfigurationError``).  The file is written atomically (sorted keys, 2-space
    indent, trailing newline); the result holds ``path``, ``bytes`` and ``sha256``.
    """

    settings = _opts(options)
    target = _impl.option_path(settings, "output")
    if target is None:
        raise ConfigurationError(
            "options.output is required.",
            code="OUTPUT_REQUIRED",
            hint="Pass the file path to write.",
        )
    if not isinstance(input, Mapping):
        raise InputFormatError("input must be a mapping.", code="INVALID_INPUT")
    if not _impl.option_bool(settings, "allow_secrets", False):
        for location, text in _impl._walk_strings(dict(input)):
            hits = _impl.scan_secrets(text)
            if hits:
                raise ConfigurationError(
                    f"Credential-looking text ({hits[0]['rule']}, {hits[0]['masked']}) found at {location}.",
                    code="SECRET_IN_CONFIG",
                    hint="Reference secrets through environment variables, or pass allow_secrets=true.",
                )
    import json

    _impl.canonical_json(input)  # rejects NaN and non-JSON values with a DNAKit error
    payload = (json.dumps(input, sort_keys=True, indent=2, ensure_ascii=False) + "\n").encode(
        "utf-8"
    )
    _impl.write_atomic(target, payload, overwrite=_impl.option_bool(settings, "overwrite", False))
    return candidate_ok(
        _MODULE,
        "save_config",
        input,
        {"path": str(target), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()},
        options=options,
    )


def set_random_seed(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Seed Python's ``random`` (and NumPy / PyTorch when available) with one integer for reproducible runs.

    The seed is ``options["seed"]`` if given, else ``input`` (an integer in [0, 2**32 - 1]).  ``random.seed`` is always
    applied; ``numpy.random.seed`` when NumPy is installed; ``torch.manual_seed`` only if PyTorch is already
    imported.  This changes global interpreter state by design.  Result: ``seed`` and the list ``seeded``;
    ``PYTHONHASHSEED`` cannot be changed after start-up and is only reported.
    """

    settings = _opts(options)
    seed = settings.get("seed", input)
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**32:
        raise InputFormatError(
            "The seed must be an integer in [0, 2**32 - 1].", code="INVALID_SEED"
        )
    random.seed(seed)
    seeded = ["random"]
    if importlib.util.find_spec("numpy") is not None:
        import numpy

        numpy.random.seed(seed)
        seeded.append("numpy")
    if "torch" in sys.modules:
        sys.modules["torch"].manual_seed(seed)
        seeded.append("torch")
    result = {"seed": seed, "seeded": seeded, "python_hash_seed": os.environ.get("PYTHONHASHSEED")}
    return candidate_ok(_MODULE, "set_random_seed", input, result, options=options)


def submit_cloud_job(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: submitting a cloud job needs a cloud account, credentials and a provider SDK."""

    return candidate_not_implemented(
        _MODULE,
        "submit_cloud_job",
        input,
        options=options,
        requires="cloud provider credentials and SDK (AWS Batch, Google Batch, Azure Batch, ...)",
        message="No cloud backend is configured, and workflow_security does not execute user commands remotely either.",
    )


def submit_hpc_job(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: submitting an HPC job means running sbatch/qsub on a cluster with a user-provided script."""

    return candidate_not_implemented(
        _MODULE,
        "submit_hpc_job",
        input,
        options=options,
        requires="a cluster scheduler (sbatch, qsub, bsub) on PATH",
        message="Job submission executes a user script through a scheduler; workflow_security never executes user commands.",
    )


def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a workflow specification: structure, dependency graph (unknown ids, cycles) and a security scan.

    ``input``: mapping, JSON/YAML text or file path (see the module description).  Checks: non-empty steps with
    unique ids matching ``[A-Za-z0-9][A-Za-z0-9_.-]{0,63}``; ``command`` a non-empty string or list of strings;
    ``cpus`` > 0 and ``memory_mb``/``runtime_s``/``disk_mb`` >= 0; ``depends_on`` references existing steps, no
    self-dependency, no cycle (the cycle path is reported).  Security scan (disable with ``security: false``):
    dangerous commands (``rm -rf /``, ``curl|sh``, fork bomb, ``dd of=/dev/..``, ``mkfs``, ``sudo``, ``chmod 777``),
    shell metacharacters in string commands (warning), literal credentials (AWS key ids, private-key blocks, GitHub
    and Slack tokens, bearer tokens, ``password=...``-style assignments, secret-named ``env`` values; secrets are
    masked in messages), and paths with ``..``, sensitive locations such as ``/etc`` or ``~/.ssh`` (errors) or absolute
    paths (warning).  Result: ``valid`` (no error), ``errors``, ``warnings``, ``order`` (topological order or null).
    """

    settings = _opts(options)
    document, analysis, source = _analysis(
        input, security=_impl.option_bool(settings, "security", True)
    )
    result = {
        "valid": not analysis["errors"],
        "errors": analysis["errors"],
        "warnings": analysis["warnings"],
        "order": analysis["order"],
        "step_count": len(analysis["steps"]),  # type: ignore[arg-type]
        "sha256": _impl.sha256_text(_impl.canonical_json(document)),
        "source": source,
    }
    return candidate_ok(_MODULE, "validate", input, result, options=options)


def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialize a workflow specification to JSON or YAML text, optionally writing it to a file.

    ``input``: a mapping (or JSON/YAML text/file path that parses to one).  Options: ``format`` ("json" default,
    "yaml"), ``output`` (file to write atomically; nothing is written without it), ``overwrite`` (default false).
    JSON is canonical (sorted keys, 2-space indent); YAML uses the safe dumper.  Result: ``text``, ``format``,
    ``sha256`` of the text and ``path``.
    """

    settings = _opts(options)
    output_format = str(settings.get("format", "json"))
    if output_format not in {"json", "yaml"}:
        raise ConfigurationError("options.format must be 'json' or 'yaml'.", code="INVALID_OPTION")
    document, _ = _impl.load_document(input)
    _impl.canonical_json(document)
    if output_format == "json":
        import json

        text = json.dumps(document, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    else:
        import yaml

        text = yaml.safe_dump(document, sort_keys=True, allow_unicode=True)
    target = _impl.option_path(settings, "output")
    written = None
    if target is not None:
        _impl.write_atomic(
            target, text.encode("utf-8"), overwrite=_impl.option_bool(settings, "overwrite", False)
        )
        written = str(target)
    return candidate_ok(
        _MODULE,
        "write",
        input,
        {"format": output_format, "text": text, "sha256": _impl.sha256_text(text), "path": written},
        options=options,
    )


__all__ = [
    "benchmark",
    "cache",
    "cancel",
    "check_dependencies",
    "check_external_tools",
    "check_numerical_reproducibility",
    "checkpoint",
    "cleanup_temp_workspace",
    "clear",
    "create",
    "decrypt_genomic_file",
    "encrypt_genomic_file",
    "estimate",
    "monitor",
    "pseudonymize_metadata_fields",
    "read",
    "record_provenance",
    "report",
    "resume",
    "retry",
    "run",
    "save_config",
    "set_random_seed",
    "submit_cloud_job",
    "submit_hpc_job",
    "validate",
    "write",
]

__candidate_adapter_operations__ = set(
    [
        "benchmark",
        "cache",
        "cancel",
        "check_dependencies",
        "check_external_tools",
        "check_numerical_reproducibility",
        "checkpoint",
        "cleanup_temp_workspace",
        "clear",
        "create",
        "decrypt_genomic_file",
        "encrypt_genomic_file",
        "estimate",
        "monitor",
        "pseudonymize_metadata_fields",
        "read",
        "record_provenance",
        "report",
        "resume",
        "retry",
        "run",
        "save_config",
        "set_random_seed",
        "submit_cloud_job",
        "submit_hpc_job",
        "validate",
        "write",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
