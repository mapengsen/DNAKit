"""Local, deterministic helpers behind the ``dnakit.workflow_security`` functions.

Nothing here executes a command, opens a network connection or modifies a file
unless the caller passes an explicit output path.  A *workflow specification* is
a plain mapping::

    {
        "name": "wf",
        "steps": [
            {
                "id": "qc",
                "command": ["fastqc", "reads.fq"],
                "depends_on": [],
                "inputs": ["reads.fq"],
                "outputs": ["qc/"],
                "env": {"THREADS": "4"},
                "cpus": 4,
                "memory_mb": 2048,
                "runtime_s": 600,
                "disk_mb": 100,
            },
            ...,
        ],
    }

``steps`` may also be a mapping ``id -> step``.  ``needs`` is accepted as an alias of ``depends_on``.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shlex
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_DOCUMENT_BYTES = 1_000_000
MAX_STEPS = 10_000
STEP_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
CHECKPOINT_FORMAT = "dnakit-checkpoint-v1"
SPEC_FORMAT = "dnakit-workflow-security-spec-v1"

# --- documents --------------------------------------------------------------------------------


def canonical_json(value: object) -> str:
    """Compact, key-sorted JSON used for every digest in this module."""

    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        )
    except (TypeError, ValueError) as exc:
        raise InputFormatError(
            "The value must be JSON-serializable (no NaN/Infinity, sets or custom objects).",
            code="NOT_JSON_SERIALIZABLE",
        ) from exc


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_text(text: str, *, source: str) -> object:
    """Parse JSON, falling back to the YAML safe loader (no tags, no code)."""

    if len(text.encode("utf-8")) > MAX_DOCUMENT_BYTES:
        raise InputFormatError(
            f"Document {source} exceeds {MAX_DOCUMENT_BYTES} bytes.", code="DOCUMENT_TOO_LARGE"
        )
    try:
        return json.loads(text)
    except ValueError:
        pass
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - PyYAML is a DNAKit dependency
        raise InputFormatError(
            "Document is not JSON and PyYAML is unavailable.", code="YAML_UNAVAILABLE"
        ) from exc
    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise InputFormatError(
            f"Document {source} is neither valid JSON nor valid YAML.",
            code="INVALID_DOCUMENT",
            hint="Pass a mapping, JSON text, YAML text or a path to a JSON/YAML file.",
        ) from exc


def load_document(value: object) -> tuple[dict[str, object], str]:
    """Return ``(mapping, source)`` from a mapping, JSON/YAML text, or a file path."""

    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}, "inline"
    if isinstance(value, os.PathLike) or (
        isinstance(value, str) and "\n" not in value and value.strip()[:1] not in {"{", "["}
    ):
        path = Path(value)  # type: ignore[arg-type]
        if not path.is_file():
            raise InputFormatError(
                f"Not a mapping, JSON/YAML text or an existing file: {str(value)[:80]!r}.",
                code="DOCUMENT_NOT_FOUND",
            )
        if path.stat().st_size > MAX_DOCUMENT_BYTES:
            raise InputFormatError(
                f"File {path} exceeds {MAX_DOCUMENT_BYTES} bytes.", code="DOCUMENT_TOO_LARGE"
            )
        parsed = parse_text(path.read_text(encoding="utf-8"), source=str(path))
        source = str(path)
    elif isinstance(value, str):
        parsed = parse_text(value, source="inline text")
        source = "inline"
    else:
        raise InputFormatError(
            f"Unsupported document input type {type(value).__name__}.",
            code="INVALID_DOCUMENT",
            hint="Pass a mapping, JSON/YAML text or a file path.",
        )
    if not isinstance(parsed, dict):
        raise InputFormatError(
            "The document must contain a mapping at its top level.", code="INVALID_DOCUMENT"
        )
    return parsed, source


def write_atomic(path: Path, payload: bytes, *, overwrite: bool) -> None:
    """Write ``payload`` through a temporary file in the target directory."""

    if not path.parent.is_dir():
        raise ConfigurationError(
            f"Output directory does not exist: {path.parent}", code="OUTPUT_DIR_MISSING"
        )
    if path.exists() and not overwrite:
        raise ConfigurationError(
            f"Refusing to overwrite existing file: {path}",
            code="OUTPUT_EXISTS",
            hint="Pass options={'overwrite': True} to replace it.",
        )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def option_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(f"options.{key} must be a boolean.", code="INVALID_OPTION")
    return value


def option_path(options: Mapping[str, object], key: str) -> Path | None:
    value = options.get(key)
    if value is None:
        return None
    if not isinstance(value, (str, os.PathLike)) or not str(value).strip():
        raise ConfigurationError(f"options.{key} must be a non-empty path.", code="INVALID_OPTION")
    return Path(value)


# --- scanning rules ---------------------------------------------------------------------------

# (code, regex) pairs for literal credentials; matches are masked in every report.
SECRET_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("AWS_ACCESS_KEY_ID", re.compile(r"\b(?:AKIA|ASIA)[0-9A-Z]{16}\b")),
    ("PRIVATE_KEY_BLOCK", re.compile(r"-----BEGIN (?:[A-Z]+ )?PRIVATE KEY-----")),
    ("GITHUB_TOKEN", re.compile(r"\bgh[pousr]_[A-Za-z0-9]{36,}\b")),
    ("SLACK_TOKEN", re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{10,}\b")),
    ("BEARER_TOKEN", re.compile(r"(?i)\bBearer\s+[A-Za-z0-9._~+/=-]{20,}")),
    (
        "CREDENTIAL_ASSIGNMENT",
        re.compile(
            r"(?i)\b[A-Za-z0-9_.-]*(?:password|passwd|secret|token|api[_-]?key|access[_-]?key)\b"
            r"\s*[=:]\s*(?![$<{(])['\"]?([^\s'\"]{4,})"
        ),
    ),
)
_SECRET_NAME = re.compile(
    r"(?i)(password|passwd|secret|token|api[_-]?key|access[_-]?key|credential)"
)
_REFERENCE_VALUE = re.compile(r"^\s*(?:\$\{?[A-Za-z_]|<.*>|\{\{.*\}\}|\*+$|$)")

# Commands that are rejected outright (error) ...
DANGEROUS_COMMAND_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "RECURSIVE_DELETE_ROOT",
        re.compile(
            r"\brm\s+(?:-[A-Za-z]*[rR][A-Za-z]*\s+)+(?:--no-preserve-root\s+)?(?:/|~|\$HOME)(?:\s|$|\*)"
        ),
    ),
    ("PIPE_TO_SHELL", re.compile(r"\b(?:curl|wget)\b[^|;&]*\|\s*(?:sudo\s+)?(?:ba|z|da)?sh\b")),
    ("FORK_BOMB", re.compile(r":\s*\(\s*\)\s*\{")),
    ("DISK_OVERWRITE", re.compile(r"\b(?:dd\s+[^;|&]*\bof=/dev/|mkfs(?:\.\w+)?\s|>\s*/dev/sd)")),
    ("PRIVILEGE_ESCALATION", re.compile(r"(?:^|[\s;&|])(?:sudo|su)\s")),
    ("WORLD_WRITABLE", re.compile(r"\bchmod\s+(?:-R\s+)?(?:0?777|a\+rwx)\b")),
)
# ... and shell metacharacters that only warrant a warning (list-form commands never reach a shell).
_SHELL_META = re.compile(r"[;&|`<>\n]|\$\(")
SENSITIVE_PATH_PREFIXES = (
    "/etc/",
    "/root/",
    "/proc/",
    "/sys/",
    "/dev/",
    "/boot/",
    "/var/run/",
    "/var/lib/docker",
)
SENSITIVE_PATH_PARTS = (".ssh", ".aws", ".gnupg", ".kube", ".docker")


def _mask(text: str) -> str:
    return text[:2] + "****" if len(text) >= 12 else "****"


def scan_secrets(text: str) -> list[dict[str, str]]:
    """Return ``[{rule, masked}]`` for credential-looking substrings; the secret is never returned."""

    found: list[dict[str, str]] = []
    for rule, pattern in SECRET_PATTERNS:
        for match in pattern.finditer(text):
            secret = match.group(1) if pattern.groups else match.group(0)
            found.append({"rule": rule, "masked": _mask(secret)})
    return found


def command_tokens(command: object) -> list[str]:
    """Split a command given as a list or a string into argv tokens (no execution)."""

    if isinstance(command, str):
        try:
            return shlex.split(command)
        except ValueError:
            return command.split()
    if isinstance(command, Sequence) and not isinstance(command, (bytes, bytearray)):
        return [str(item) for item in command]
    return []


def scan_path(path: str) -> list[tuple[str, str, str]]:
    """Return ``(severity, code, message)`` findings for one path string."""

    findings: list[tuple[str, str, str]] = []
    if "\x00" in path:
        return [("error", "NUL_IN_PATH", "path contains a NUL byte")]
    parts = Path(path.replace("\\", "/")).parts
    if ".." in parts:
        findings.append(("error", "PATH_TRAVERSAL", f"path {path!r} contains '..'"))
    normalized = path.replace("\\", "/")
    if normalized.startswith("~") or normalized.startswith("/"):
        sensitive = normalized.startswith(SENSITIVE_PATH_PREFIXES) or normalized == "/"
        sensitive = sensitive or any(part in SENSITIVE_PATH_PARTS for part in parts)
        if sensitive:
            findings.append(
                ("error", "SENSITIVE_PATH", f"path {path!r} points to a sensitive location")
            )
        else:
            findings.append(
                (
                    "warning",
                    "ABSOLUTE_PATH",
                    f"path {path!r} is absolute; prefer paths relative to the workspace",
                )
            )
    elif any(part in SENSITIVE_PATH_PARTS for part in parts):
        findings.append(
            ("error", "SENSITIVE_PATH", f"path {path!r} points to a sensitive location")
        )
    return findings


def _walk_strings(value: object, prefix: str = "") -> Iterable[tuple[str, str]]:
    """Yield ``(location, text)`` for every string in a nested JSON-like value."""

    if isinstance(value, str):
        yield prefix, value
    elif isinstance(value, Mapping):
        for key, item in value.items():
            yield from _walk_strings(item, f"{prefix}.{key}" if prefix else str(key))
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        for index, item in enumerate(value):
            yield from _walk_strings(item, f"{prefix}[{index}]")


# --- workflow analysis ------------------------------------------------------------------------


def _finding(severity: str, code: str, message: str, step: str | None = None) -> dict[str, object]:
    return {"severity": severity, "code": code, "step": step, "message": message}


def normalize_steps(
    document: Mapping[str, object],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Return ``(steps, findings)``; each step keeps its fields and gets a list ``depends_on``."""

    findings: list[dict[str, object]] = []
    raw = document.get("steps")
    items: list[tuple[str | None, object]]
    if isinstance(raw, Mapping):
        items = [(str(key), value) for key, value in raw.items()]
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        items = [(None, value) for value in raw]
    else:
        return [], [
            _finding("error", "NO_STEPS", "the workflow needs a non-empty 'steps' list or mapping")
        ]
    if not items:
        return [], [
            _finding("error", "NO_STEPS", "the workflow needs a non-empty 'steps' list or mapping")
        ]
    if len(items) > MAX_STEPS:
        raise InputFormatError(f"At most {MAX_STEPS} steps are supported.", code="TOO_MANY_STEPS")
    steps: list[dict[str, object]] = []
    seen: set[str] = set()
    for index, (key, value) in enumerate(items):
        if not isinstance(value, Mapping):
            findings.append(
                _finding("error", "STEP_NOT_MAPPING", f"step {index} must be a mapping")
            )
            continue
        step = {str(k): v for k, v in value.items()}
        identifier = step.get("id", key)
        if not isinstance(identifier, str) or not STEP_ID.match(identifier):
            findings.append(
                _finding(
                    "error",
                    "INVALID_STEP_ID",
                    f"step {index} needs an id matching {STEP_ID.pattern}",
                )
            )
            continue
        if identifier in seen:
            findings.append(
                _finding(
                    "error",
                    "DUPLICATE_STEP_ID",
                    f"step id {identifier!r} is used more than once",
                    identifier,
                )
            )
            continue
        seen.add(identifier)
        step["id"] = identifier
        dependencies = step.get("depends_on", step.get("needs", []))
        if isinstance(dependencies, str):
            dependencies = [dependencies]
        if not isinstance(dependencies, Sequence) or any(
            not isinstance(item, str) for item in dependencies
        ):
            findings.append(
                _finding(
                    "error",
                    "INVALID_DEPENDENCIES",
                    "depends_on must be a list of step ids",
                    identifier,
                )
            )
            dependencies = []
        step["depends_on"] = list(dependencies)
        steps.append(step)
    return steps, findings


def topological_order(
    steps: Sequence[Mapping[str, object]],
) -> tuple[list[str] | None, list[dict[str, object]]]:
    """Kahn ordering (ties broken by declaration order) plus dependency/cycle findings.

    Returns ``(order, findings)``; ``order`` is ``None`` when a dependency is unknown or a cycle exists.
    """

    findings: list[dict[str, object]] = []
    ids = [str(step["id"]) for step in steps]
    known = set(ids)
    graph: dict[str, list[str]] = {identifier: [] for identifier in ids}
    for step in steps:
        identifier = str(step["id"])
        for dependency in step["depends_on"]:  # type: ignore[union-attr]
            if dependency == identifier:
                findings.append(
                    _finding(
                        "error", "SELF_DEPENDENCY", "a step cannot depend on itself", identifier
                    )
                )
            elif dependency not in known:
                findings.append(
                    _finding(
                        "error",
                        "UNKNOWN_DEPENDENCY",
                        f"depends on unknown step {dependency!r}",
                        identifier,
                    )
                )
            else:
                graph[identifier].append(dependency)
    if findings:
        return None, findings
    remaining = {identifier: len(set(deps)) for identifier, deps in graph.items()}
    children: dict[str, list[str]] = {identifier: [] for identifier in ids}
    for identifier, deps in graph.items():
        for dependency in set(deps):
            children[dependency].append(identifier)
    position = {identifier: index for index, identifier in enumerate(ids)}
    ready = [identifier for identifier in ids if remaining[identifier] == 0]
    order: list[str] = []
    while ready:
        ready.sort(key=position.__getitem__)
        current = ready.pop(0)
        order.append(current)
        for child in children[current]:
            remaining[child] -= 1
            if remaining[child] == 0:
                ready.append(child)
    if len(order) != len(ids):
        cycle = _find_cycle(graph)
        findings.append(_finding("error", "CYCLE", "dependency cycle: " + " -> ".join(cycle)))
        return None, findings
    return order, findings


def _find_cycle(graph: Mapping[str, Sequence[str]]) -> list[str]:
    """One concrete cycle ``[a, b, ..., a]`` (iterative DFS; ``graph[x]`` lists the dependencies of x)."""

    state: dict[str, int] = {}
    for root in graph:
        if state.get(root):
            continue
        stack = [(root, iter(graph[root]))]
        path = [root]
        state[root] = 1
        while stack:
            node, iterator = stack[-1]
            advanced = False
            for neighbour in iterator:
                if state.get(neighbour) == 1:
                    return path[path.index(neighbour) :] + [neighbour]
                if not state.get(neighbour):
                    state[neighbour] = 1
                    stack.append((neighbour, iter(graph[neighbour])))
                    path.append(neighbour)
                    advanced = True
                    break
            if not advanced:
                state[node] = 2
                stack.pop()
                path.pop()
    return []


def _is_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def analyze_workflow(document: Mapping[str, object], *, security: bool = True) -> dict[str, object]:
    """Validate structure, resources, dependency graph and (optionally) scan commands, secrets and paths."""

    steps, findings = normalize_steps(document)
    order, graph_findings = topological_order(steps) if steps else (None, [])
    findings.extend(graph_findings)
    for step in steps:
        identifier = str(step["id"])
        command = step.get("command")
        if command is not None:
            if isinstance(command, str):
                valid = bool(command.strip())
            else:
                valid = (
                    isinstance(command, Sequence)
                    and not isinstance(command, (bytes, bytearray))
                    and len(command) > 0
                    and all(isinstance(item, str) for item in command)
                    and bool(command[0].strip())
                )
            if not valid:
                findings.append(
                    _finding(
                        "error",
                        "INVALID_COMMAND",
                        "command must be a non-empty string or list of strings",
                        identifier,
                    )
                )
        for field in ("cpus", "memory_mb", "runtime_s", "disk_mb"):
            if field not in step:
                continue
            number = step[field]
            if not _is_number(number) or number < 0 or (field == "cpus" and number == 0):  # type: ignore[operator]
                requirement = "positive" if field == "cpus" else "non-negative"
                findings.append(
                    _finding(
                        "error",
                        "INVALID_RESOURCE",
                        f"{field} must be a {requirement} finite number",
                        identifier,
                    )
                )
        for field in ("inputs", "outputs"):
            value = step.get(field, [])
            if isinstance(value, str):
                value = [value]
            if not isinstance(value, Sequence) or any(not isinstance(item, str) for item in value):
                findings.append(
                    _finding(
                        "error",
                        "INVALID_PATH_LIST",
                        f"{field} must be a list of path strings",
                        identifier,
                    )
                )
        if "env" in step and not isinstance(step["env"], Mapping):
            findings.append(_finding("error", "INVALID_ENV", "env must be a mapping", identifier))
    if security:
        findings.extend(_security_findings(document, steps))
    return {
        "steps": steps,
        "order": order,
        "findings": findings,
        "errors": [item for item in findings if item["severity"] == "error"],
        "warnings": [item for item in findings if item["severity"] == "warning"],
    }


def _security_findings(
    document: Mapping[str, object], steps: Sequence[Mapping[str, object]]
) -> list[dict[str, object]]:
    findings: list[dict[str, object]] = []
    for step in steps:
        identifier = str(step["id"])
        command = step.get("command")
        if command is not None:
            joined = command if isinstance(command, str) else " ".join(command_tokens(command))
            for code, pattern in DANGEROUS_COMMAND_PATTERNS:
                if pattern.search(joined):
                    findings.append(
                        _finding(
                            "error", code, "command matches a dangerous-command rule", identifier
                        )
                    )
            if isinstance(command, str) and _SHELL_META.search(command):
                findings.append(
                    _finding(
                        "warning",
                        "SHELL_METACHARACTERS",
                        "string command contains shell metacharacters; use a list-form argv",
                        identifier,
                    )
                )
            elif not isinstance(command, str) and any(
                _SHELL_META.search(token) for token in command_tokens(command)[:1]
            ):
                findings.append(
                    _finding(
                        "warning",
                        "SHELL_METACHARACTERS",
                        "the program name contains shell metacharacters",
                        identifier,
                    )
                )
        for field in ("inputs", "outputs"):
            value = step.get(field, [])
            for path in (
                [value]
                if isinstance(value, str)
                else (value if isinstance(value, Sequence) else [])
            ):
                if isinstance(path, str):
                    for severity, code, message in scan_path(path):
                        findings.append(_finding(severity, code, message, identifier))
        env = step.get("env")
        if isinstance(env, Mapping):
            for name, value in env.items():
                if (
                    _SECRET_NAME.search(str(name))
                    and isinstance(value, str)
                    and not _REFERENCE_VALUE.match(value)
                ):
                    findings.append(
                        _finding(
                            "error",
                            "SECRET_IN_ENV",
                            f"env {name!r} holds a literal secret-like value ({_mask(value)})",
                            identifier,
                        )
                    )
    for location, text in _walk_strings(document):
        step_id = _step_of(location, steps)
        for hit in scan_secrets(text):
            findings.append(
                _finding(
                    "error",
                    "SECRET_DETECTED",
                    f"{hit['rule']} at {location}: {hit['masked']}",
                    step_id,
                )
            )
    # Deduplicate (a literal env secret is reported by both the env rule and the generic scan).
    unique: dict[tuple[object, ...], dict[str, object]] = {}
    for item in findings:
        unique.setdefault((item["code"], item["step"], item["message"]), item)
    return list(unique.values())


def _step_of(location: str, steps: Sequence[Mapping[str, object]]) -> str | None:
    match = re.match(r"steps\[(\d+)\]", location)
    if match and int(match.group(1)) < len(steps):
        return str(steps[int(match.group(1))]["id"])
    return None


def raise_if_invalid(analysis: Mapping[str, object], *, what: str) -> None:
    errors = analysis["errors"]
    if errors:
        raise InputFormatError(
            f"The workflow is not valid for {what}: {errors[0]['message']}",  # type: ignore[index]
            code="INVALID_WORKFLOW_SPEC",
            context={"errors": [dict(item) for item in errors]},  # type: ignore[union-attr]
            hint="Run workflow_security.validate for the complete list of findings.",
        )


# --- resource estimation ----------------------------------------------------------------------


def estimate_resources(
    steps: Sequence[Mapping[str, object]],
    order: Sequence[str],
    *,
    default_cpus: float,
    default_memory_mb: float,
    default_runtime_s: float,
    default_disk_mb: float,
) -> dict[str, object]:
    """As-soon-as-possible schedule with unlimited workers.

    A step starts when all its dependencies have ended.  ``wall_time_s`` is the critical-path length;
    peaks are taken over the steps that run simultaneously (half-open ``[start, end)`` intervals of
    positive duration).  ``total_cpu_seconds`` is the sum of ``cpus * runtime_s``.
    """

    by_id = {str(step["id"]): step for step in steps}
    start: dict[str, float] = {}
    end: dict[str, float] = {}
    via: dict[str, str | None] = {}
    for identifier in order:
        step = by_id[identifier]
        dependencies = list(step["depends_on"])  # type: ignore[call-overload]
        begin = max((end[item] for item in dependencies), default=0.0)
        via[identifier] = max(dependencies, key=lambda item: end[item]) if dependencies else None
        runtime = float(step.get("runtime_s", default_runtime_s))  # type: ignore[arg-type]
        start[identifier], end[identifier] = begin, begin + runtime
    wall = max(end.values(), default=0.0)
    tail = max(order, key=lambda item: end[item]) if order else None
    critical: list[str] = []
    while tail is not None:
        critical.append(tail)
        tail = via[tail]
    critical.reverse()

    def resource(step: Mapping[str, object], name: str, default: float) -> float:
        return float(step.get(name, default))  # type: ignore[arg-type]

    events: list[tuple[float, int, float, float]] = []
    for identifier in order:
        if end[identifier] > start[identifier]:
            step = by_id[identifier]
            cpus, memory = (
                resource(step, "cpus", default_cpus),
                resource(step, "memory_mb", default_memory_mb),
            )
            events.append((start[identifier], 1, cpus, memory))
            events.append((end[identifier], 0, -cpus, -memory))
    events.sort(
        key=lambda item: (item[0], item[1])
    )  # at equal times, ends (0) are processed before starts (1)
    cpus_now = memory_now = peak_cpus = peak_memory = 0.0
    concurrent = peak_concurrent = 0
    for _, kind, cpu_delta, memory_delta in events:
        cpus_now += cpu_delta
        memory_now += memory_delta
        concurrent += 1 if kind else -1
        peak_cpus, peak_memory = max(peak_cpus, cpus_now), max(peak_memory, memory_now)
        peak_concurrent = max(peak_concurrent, concurrent)
    return {
        "wall_time_s": wall,
        "critical_path": critical,
        "total_cpu_seconds": sum(
            resource(by_id[item], "cpus", default_cpus)
            * resource(by_id[item], "runtime_s", default_runtime_s)
            for item in order
        ),
        "peak_cpus": peak_cpus,
        "peak_memory_mb": peak_memory,
        "peak_concurrent_steps": peak_concurrent,
        "total_disk_mb": sum(resource(by_id[item], "disk_mb", default_disk_mb) for item in order),
        "schedule": [{"id": item, "start_s": start[item], "end_s": end[item]} for item in order],
        "steps_without_runtime": [item for item in order if "runtime_s" not in by_id[item]],
    }
