"""Pipeline-definition tools for WDL, CWL, Snakemake, nf-core and Galaxy.

These functions are not DNA algorithms. They parse and validate workflow
definitions, list their tasks, rules, processes or steps, plan a Snakemake
dry run, and read Galaxy server metadata. No sequence is read or computed,
and no pipeline step is executed.

- ``parse_wdl_document()`` and ``check_wdl_inputs()`` use miniwdl (``WDL``).
- ``read_cwl_workflow()`` loads a CWL workflow with cwltool and schema-salad.
- ``read_snakefile()`` and ``plan_snakemake_jobs()`` use the Snakemake API.
- ``list_nfcore_processes()`` reads an nf-core checkout offline with nf_core.
- ``list_galaxy_tools()``, ``list_galaxy_histories()`` and
  ``list_galaxy_workflows()`` call a Galaxy server through bioblend.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata as metadata
import math
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendError, ConfigurationError

_EXTRAS = {
    "WDL": "bio-wdl",
    "cwltool": "bio-cwl",
    "snakemake": "bio-snakemake",
    "nf_core": "bio-nf-core",
    "bioblend": "bio-bioblend",
}
_DEFAULT_MAX_BYTES = 1_000_000
_MAX_BYTES_LIMIT = 50_000_000
_DEFAULT_MAX_ITEMS = 1_000
_MAX_ITEMS_LIMIT = 100_000
_MAX_NAMES = 1_000
_MAX_COMPONENTS = 10_000
_DEFAULT_MAX_JOBS = 10_000
_MAX_JOBS_LIMIT = 1_000_000


def _limit(value: object, name: str, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [1, {maximum}].")
    return value


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def _names(values: object, name: str) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Iterable):
        raise ConfigurationError(f"{name} must be a list of strings.")
    items = sorted({_text(item, name) for item in values})
    if len(items) > _MAX_NAMES:
        raise ConfigurationError(f"{name} accepts at most {_MAX_NAMES} entries.")
    return items


def _document(source: str | Path, max_bytes: object) -> Path:
    path = input_path(source)
    limit = _limit(max_bytes, "max_bytes", _MAX_BYTES_LIMIT)
    size = path.stat().st_size
    if size > limit:
        raise ConfigurationError(f"{path.name} has {size} bytes; max_bytes is {limit}.")
    return path


def _directory(source: str | Path, name: str) -> Path:
    path = Path(source).expanduser().resolve()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing directory: {path}")
    return path


def _plain(value: object) -> Any:
    """Convert upstream containers and scalars to JSON-serializable Python values."""

    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if value is None or isinstance(value, (bool, int, str)):
        return value
    return str(value)


def _provenance(distribution: str, **notes: str) -> dict[str, str]:
    return {"package": f"{distribution} {metadata.version(distribution)}", **notes}


# miniwdl ---------------------------------------------------------------------


def _wdl_document(path: Path) -> tuple[Any, Any]:
    wdl = require("WDL", _EXTRAS["WDL"])
    try:
        document = wdl.load(str(path))
    except (
        wdl.Error.SyntaxError,
        wdl.Error.ValidationError,
        wdl.Error.MultipleValidationErrors,
        wdl.Error.ImportError,
    ) as exc:
        raise ConfigurationError(f"miniwdl rejected {path.name}: {exc}") from exc
    return wdl, document


def _wdl_bindings(bindings: Iterable[Any], required: set[str]) -> list[dict[str, Any]]:
    return [
        {"name": item.name, "type": str(item.value.type), "required": item.name in required}
        for item in sorted(bindings, key=lambda item: item.name)
    ]


def _wdl_task(task: Any) -> dict[str, Any]:
    required = {item.name for item in task.required_inputs}
    outputs = sorted(task.effective_outputs, key=lambda item: item.name)
    return {
        "name": task.name,
        "inputs": _wdl_bindings(task.available_inputs, required),
        "outputs": [{"name": item.name, "type": str(item.value)} for item in outputs],
    }


def _wdl_calls(wdl: Any, nodes: Iterable[Any]) -> list[Any]:
    calls: list[Any] = []
    for node in nodes:
        if isinstance(node, wdl.Tree.Call):
            calls.append(node)
        else:
            calls.extend(_wdl_calls(wdl, getattr(node, "body", ())))
    return calls


def _wdl_workflow(wdl: Any, workflow: Any) -> dict[str, Any]:
    required = {item.name for item in workflow.required_inputs}
    calls = sorted(_wdl_calls(wdl, workflow.body), key=lambda call: call.name)
    outputs = sorted(workflow.effective_outputs, key=lambda item: item.name)
    return {
        "name": workflow.name,
        "inputs": _wdl_bindings(workflow.available_inputs, required),
        "outputs": [{"name": item.name, "type": str(item.value)} for item in outputs],
        "calls": [{"name": call.name, "task": call.callee.name} for call in calls],
    }


def parse_wdl_document(
    source: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a WDL document and list its tasks, inputs, outputs and workflow calls."""

    path = _document(source, max_bytes)
    wdl, document = _wdl_document(path)
    workflow = None if document.workflow is None else _wdl_workflow(wdl, document.workflow)
    return {
        "status": "completed",
        "operation": "parse_wdl_document",
        "source": str(path),
        "wdl_version": _plain(document.effective_wdl_version),
        "imports": sorted(item.namespace for item in document.imports),
        "tasks": [_wdl_task(task) for task in sorted(document.tasks, key=lambda t: t.name)],
        "workflow": workflow,
        "provenance": _provenance(
            "miniwdl",
            checks="miniwdl parses and type-checks the document; no task is run",
            imports="imports are followed by miniwdl and must exist on disk",
        ),
    }


def _input_mapping(inputs: object) -> dict[str, object]:
    if not isinstance(inputs, Mapping) or len(inputs) > _MAX_NAMES:
        raise ConfigurationError(f"inputs must be a mapping with at most {_MAX_NAMES} entries.")
    return {_text(key, "input name"): value for key, value in inputs.items()}


def check_wdl_inputs(
    source: str | Path,
    inputs: Mapping[str, object],
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Check a JSON-style input mapping against the inputs of a WDL workflow."""

    path = _document(source, max_bytes)
    values = _input_mapping(inputs)
    wdl, document = _wdl_document(path)
    workflow = document.workflow
    if workflow is None:
        raise ConfigurationError("check_wdl_inputs needs a document that contains a workflow.")
    declared = {item.name: item.value for item in workflow.available_inputs}
    required = {item.name for item in workflow.required_inputs}
    provided = set(values)
    missing = sorted(required - provided)
    unknown = sorted(provided - set(declared))
    invalid: list[dict[str, str]] = []
    for name in sorted(provided & set(declared)):
        try:
            wdl.Value.from_json(declared[name].type, values[name])
        except (wdl.Error.InputError, TypeError) as exc:
            invalid.append({"name": name, "type": str(declared[name].type), "message": str(exc)})
    return {
        "status": "completed",
        "operation": "check_wdl_inputs",
        "source": str(path),
        "workflow": workflow.name,
        "valid": not (missing or unknown or invalid),
        "missing": missing,
        "unknown": unknown,
        "invalid": invalid,
        "provenance": _provenance(
            "miniwdl",
            checks="input names, required inputs and JSON-to-WDL value types",
            files="File values are not checked for existence",
        ),
    }


# cwltool ---------------------------------------------------------------------


def _cwl_name(identifier: object) -> str:
    return str(identifier).split("#", 1)[-1]


def _cwl_sources(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [_cwl_name(value)]
    if isinstance(value, list):
        return sorted(_cwl_name(item) for item in value)
    return []


def _cwl_location(uri: str, base: Path) -> str:
    parsed = urlsplit(uri)
    if parsed.scheme != "file":
        return uri
    file_path = Path(unquote(parsed.path))
    try:
        location = file_path.relative_to(base).as_posix()
    except ValueError:
        location = file_path.as_posix()
    return f"{location}#{parsed.fragment}" if parsed.fragment else location


def _cwl_step(step: Any, base: Path) -> dict[str, Any]:
    run = step.embedded_tool.tool
    inputs = sorted(
        (
            {
                "name": _cwl_name(item["id"]).rsplit("/", 1)[-1],
                "sources": _cwl_sources(item.get("source")),
            }
            for item in step.tool.get("inputs", [])
        ),
        key=lambda item: item["name"],
    )
    outputs = sorted(
        _cwl_name(item["id"]).rsplit("/", 1)[-1] for item in step.tool.get("outputs", [])
    )
    return {
        "name": _cwl_name(step.id),
        "run": _cwl_location(str(run.get("id", "")), base),
        "run_class": _plain(run.get("class")),
        "inputs": inputs,
        "outputs": outputs,
    }


def read_cwl_workflow(
    source: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a CWL Workflow with cwltool and list its inputs, steps and outputs."""

    path = _document(source, max_bytes)
    require("cwltool", _EXTRAS["cwltool"])
    from cwltool.context import LoadingContext
    from cwltool.load_tool import load_tool
    from cwltool.workflow import default_make_tool
    from schema_salad.exceptions import SchemaSaladException

    try:
        process = load_tool(str(path), LoadingContext({"construct_tool_object": default_make_tool}))
    except SchemaSaladException as exc:
        raise ConfigurationError(f"cwltool rejected {path.name}: {exc}") from exc
    document = process.tool
    if document.get("class") != "Workflow":
        raise ConfigurationError(
            f"read_cwl_workflow needs class: Workflow; {path.name} is {document.get('class')}."
        )
    base = path.parent
    steps = sorted(process.steps, key=lambda step: _cwl_name(step.id))
    inputs = sorted(
        (
            {"name": _cwl_name(item["id"]), "type": _plain(item.get("type"))}
            for item in document.get("inputs", [])
        ),
        key=lambda item: item["name"],
    )
    outputs = sorted(
        (
            {
                "name": _cwl_name(item["id"]),
                "type": _plain(item.get("type")),
                "sources": _cwl_sources(item.get("outputSource")),
            }
            for item in document.get("outputs", [])
        ),
        key=lambda item: item["name"],
    )
    return {
        "status": "completed",
        "operation": "read_cwl_workflow",
        "source": str(path),
        "class": "Workflow",
        "inputs": inputs,
        "steps": [_cwl_step(step, base) for step in steps],
        "outputs": outputs,
        "provenance": _provenance(
            "cwltool",
            checks="schema-salad validation while loading; no step is executed",
            order="steps are sorted by id because cwltool's order depends on PYTHONHASHSEED",
        ),
    }


# Snakemake -------------------------------------------------------------------


def _snakemake_api() -> Any:
    from snakemake.api import SnakemakeApi
    from snakemake.settings.enums import Quietness
    from snakemake.settings.types import OutputSettings

    return SnakemakeApi(OutputSettings(quiet={Quietness.ALL}))


def _snakemake_workflow(api: Any, path: Path, workdir: Path) -> Any:
    from snakemake.settings.types import (
        ConfigSettings,
        DeploymentSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    return api.workflow(
        resource_settings=ResourceSettings(),
        config_settings=ConfigSettings(),
        storage_settings=StorageSettings(),
        workflow_settings=WorkflowSettings(),
        deployment_settings=DeploymentSettings(),
        snakefile=path,
        workdir=workdir,
    )


def _snakemake_errors() -> tuple[type[BaseException], ...]:
    from snakemake.exceptions import RuleException, WorkflowError

    return (SyntaxError, WorkflowError, RuleException)


def read_snakefile(
    source: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a Snakefile and list its rules with their inputs and outputs."""

    path = _document(source, max_bytes)
    require("snakemake", _EXTRAS["snakemake"])
    errors = _snakemake_errors()
    # Snakemake writes .snakemake/ metadata into the workdir, so parsing uses a scratch dir.
    with tempfile.TemporaryDirectory(prefix="dnakit-snakemake-") as scratch:
        try:
            with _snakemake_api() as api:
                # Private upstream attribute; the bio-snakemake range pins it.
                workflow = _snakemake_workflow(api, path, Path(scratch))._workflow
                rules = sorted(workflow.rules, key=lambda rule: rule.name)
                rule_rows = [
                    {
                        "name": rule.name,
                        "inputs": [str(item) for item in rule.input],
                        "outputs": [str(item) for item in rule.output],
                    }
                    for rule in rules
                ]
                default_target = _plain(workflow.default_target)
        except errors as exc:
            raise ConfigurationError(f"Snakemake rejected {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_snakefile",
        "source": str(path),
        "default_target": default_target,
        "rules": rule_rows,
        "provenance": _provenance(
            "snakemake",
            checks="parsing and include resolution only; no rule is run",
            order="rules are sorted by name; includes resolve relative to the Snakefile",
        ),
    }


def plan_snakemake_jobs(
    source: str | Path,
    *,
    targets: Iterable[str] = (),
    workdir: str | Path | None = None,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    max_jobs: int = _DEFAULT_MAX_JOBS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Dry-run a Snakefile and list the jobs Snakemake would run for the targets."""

    path = _document(source, max_bytes)
    selected = _names(targets, "targets")
    limit = _limit(max_jobs, "max_jobs", _MAX_JOBS_LIMIT)
    base = path.parent if workdir is None else _directory(workdir, "workdir")
    require("snakemake", _EXTRAS["snakemake"])
    from snakemake.settings.types import DAGSettings

    errors = _snakemake_errors()
    try:
        with _snakemake_api() as api:
            workflow_api = _snakemake_workflow(api, path, base)
            # DAGApi stores its settings on the upstream workflow when it is created.
            workflow_api.dag(DAGSettings(targets=frozenset(selected)))
            workflow = workflow_api._workflow  # private upstream attribute, pinned by bio-snakemake
            workflow._prepare_dag(
                forceall=False, ignore_incomplete=True, lock_warn_only=True, nolock=True
            )
            workflow._build_dag()
            jobs = list(workflow.dag.needrun_jobs(exclude_finished=True))
            default_target = _plain(workflow.default_target)
            rows = [
                {
                    "rule": job.rule.name,
                    "inputs": [str(item) for item in job.input],
                    "outputs": [str(item) for item in job.output],
                    "wildcards": {
                        key: str(value) for key, value in sorted(job.wildcards_dict.items())
                    },
                }
                for job in jobs
            ]
    except errors as exc:
        raise ConfigurationError(f"Snakemake could not plan {path.name}: {exc}") from exc
    if len(rows) > limit:
        raise ConfigurationError(f"the plan has {len(rows)} jobs; max_jobs is {limit}.")
    rows.sort(key=lambda row: (row["rule"], row["outputs"], row["inputs"]))
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["rule"]] = counts.get(row["rule"], 0) + 1
    return {
        "status": "completed",
        "operation": "plan_snakemake_jobs",
        "source": str(path),
        "workdir": str(base),
        "targets": selected,
        "default_target": default_target,
        "job_count": len(rows),
        "jobs_by_rule": dict(sorted(counts.items())),
        "jobs": rows,
        "provenance": _provenance(
            "snakemake",
            mode="dry run: no rule is executed and no lock is taken",
            files="inputs are checked in the workdir, which gets a .snakemake/ folder with a log",
            order="jobs are sorted by rule, outputs and inputs",
        ),
    }


# nf-core ---------------------------------------------------------------------


def list_nfcore_processes(
    source: str | Path,
    *,
    max_components: int = 500,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the Nextflow processes in an nf-core pipeline checkout, offline."""

    base = _directory(source, "source")
    component_limit = _limit(max_components, "max_components", _MAX_COMPONENTS)
    size_limit = _limit(max_bytes, "max_bytes", _MAX_BYTES_LIMIT)
    require("nf_core", _EXTRAS["nf_core"])
    from nf_core.components.nfcore_component import NFCoreComponent
    from nf_core.modules.modules_utils import scan_modules_dir

    modules = base / "modules"
    if not modules.is_dir():
        raise ConfigurationError(
            f"{base} has no modules/ directory; not an nf-core pipeline checkout."
        )
    found: list[tuple[str, str, Path, bool]] = []
    for kind in ("nf-core", "local"):
        root = modules / kind
        for name in scan_modules_dir(root):
            found.append((kind, name, root / name, True))
        if root.is_dir():
            for file in root.glob("*.nf"):
                found.append((kind, file.stem, file, False))
    if len(found) > component_limit:
        raise ConfigurationError(
            f"found {len(found)} components; max_components is {component_limit}."
        )
    rows: list[dict[str, Any]] = []
    for kind, name, location, folder in found:
        main = location / "main.nf" if folder else location
        size = main.stat().st_size
        if size > size_limit:
            raise ConfigurationError(f"{main.name} has {size} bytes; max_bytes is {size_limit}.")
        component = NFCoreComponent(
            name, None, location, "pipeline", base, "modules", remote_component=folder
        )
        rows.append(
            {
                "name": name,
                "kind": kind,
                "process": component.process_name or None,
                "path": main.relative_to(base).as_posix(),
            }
        )
    rows.sort(key=lambda row: (row["kind"], row["name"], row["path"]))
    return {
        "status": "completed",
        "operation": "list_nfcore_processes",
        "source": str(base),
        "component_count": len(rows),
        "process_count": sum(1 for row in rows if row["process"] is not None),
        "components": rows,
        "provenance": _provenance(
            "nf-core",
            mode="offline: main.nf files are read from the checkout; nothing is cloned or fetched",
            scope="process names and paths only; containers and channel shapes are not read",
            order="components are sorted by kind, name and path",
        ),
    }


# Galaxy (bioblend) -----------------------------------------------------------


def _galaxy_url(value: object) -> str:
    text = _text(value, "url")
    parts = urlsplit(text)
    if parts.scheme not in ("http", "https") or not parts.hostname:
        raise ConfigurationError("url must be an http:// or https:// Galaxy server address.")
    if parts.username or parts.password or parts.query or parts.fragment:
        raise ConfigurationError(
            "url must not contain credentials, a query or a fragment; pass api_key."
        )
    return text.rstrip("/")


def _galaxy_client(url: object, api_key: object) -> tuple[Any, Any, str]:
    bioblend = require("bioblend", _EXTRAS["bioblend"])
    from bioblend.galaxy import GalaxyInstance

    address = _galaxy_url(url)
    key = None if api_key is None else _text(api_key, "api_key")
    return bioblend, GalaxyInstance(url=address, key=key), address


def _galaxy_rows(
    bioblend: Any, operation: str, request: Any, max_items: int
) -> list[dict[str, Any]]:
    requests = require("requests", _EXTRAS["bioblend"])
    try:
        rows = request()
    except (bioblend.ConnectionError, requests.RequestException) as exc:
        raise BackendError(
            f"Galaxy {operation} request failed: {exc}",
            hint="Check the server URL and api_key.",
        ) from exc
    if not isinstance(rows, list):
        raise BackendError(f"Galaxy {operation} response was not a list.")
    if len(rows) > max_items:
        raise ConfigurationError(
            f"Galaxy returned {len(rows)} {operation}; max_items is {max_items}."
        )
    return rows


def _galaxy_result(
    operation: str, address: str, method: str, items: list[dict[str, Any]]
) -> dict[str, Any]:
    items.sort(key=lambda item: str(item.get("id")))
    return {
        "status": "completed",
        "operation": operation,
        "server": address,
        "count": len(items),
        "items": items,
        "provenance": _provenance(
            "bioblend",
            request=f"GET through bioblend {method}",
            secrets="api_key is sent in the x-api-key header and is never included in the output",
            order="items are sorted by id",
        ),
    }


def list_galaxy_tools(
    url: str,
    *,
    api_key: str | None = None,
    name: str | None = None,
    max_items: int = _DEFAULT_MAX_ITEMS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the tools installed on a Galaxy server through bioblend."""

    bioblend, gi, address = _galaxy_client(url, api_key)
    selector = None if name is None else _text(name, "name")
    limit = _limit(max_items, "max_items", _MAX_ITEMS_LIMIT)
    rows = _galaxy_rows(bioblend, "tools", lambda: gi.tools.get_tools(name=selector), limit)
    items = [
        {
            "id": _plain(row.get("id")),
            "name": _plain(row.get("name")),
            "version": _plain(row.get("version")),
            "section": _plain(row.get("panel_section_name")),
        }
        for row in rows
    ]
    return _galaxy_result("list_galaxy_tools", address, "tools.get_tools", items)


def list_galaxy_histories(
    url: str,
    *,
    api_key: str | None = None,
    name: str | None = None,
    max_items: int = _DEFAULT_MAX_ITEMS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the active histories of the Galaxy user through bioblend."""

    bioblend, gi, address = _galaxy_client(url, api_key)
    selector = None if name is None else _text(name, "name")
    limit = _limit(max_items, "max_items", _MAX_ITEMS_LIMIT)
    rows = _galaxy_rows(
        bioblend, "histories", lambda: gi.histories.get_histories(name=selector), limit
    )
    items = [
        {
            "id": _plain(row.get("id")),
            "name": _plain(row.get("name")),
            "count": _plain(row.get("count")),
            "deleted": _plain(row.get("deleted")),
            "published": _plain(row.get("published")),
        }
        for row in rows
    ]
    return _galaxy_result("list_galaxy_histories", address, "histories.get_histories", items)


def list_galaxy_workflows(
    url: str,
    *,
    api_key: str | None = None,
    name: str | None = None,
    max_items: int = _DEFAULT_MAX_ITEMS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the stored Galaxy workflows, the pipeline definitions, through bioblend."""

    bioblend, gi, address = _galaxy_client(url, api_key)
    selector = None if name is None else _text(name, "name")
    limit = _limit(max_items, "max_items", _MAX_ITEMS_LIMIT)
    rows = _galaxy_rows(
        bioblend, "workflows", lambda: gi.workflows.get_workflows(name=selector), limit
    )
    items = [
        {
            "id": _plain(row.get("id")),
            "name": _plain(row.get("name")),
            "published": _plain(row.get("published")),
        }
        for row in rows
    ]
    return _galaxy_result("list_galaxy_workflows", address, "workflows.get_workflows", items)


__all__ = [
    "check_wdl_inputs",
    "list_galaxy_histories",
    "list_galaxy_tools",
    "list_galaxy_workflows",
    "list_nfcore_processes",
    "parse_wdl_document",
    "plan_snakemake_jobs",
    "read_cwl_workflow",
    "read_snakefile",
]

_register_public_verbose(globals(), tuple(__all__))
