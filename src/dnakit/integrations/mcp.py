"""On-demand connections to verified, separately installed genomics MCP servers."""

from __future__ import annotations

import asyncio
import importlib
import json
import logging
import math
import os
import re
import shutil
from collections.abc import Awaitable, Callable, Iterator
from contextlib import contextmanager, suppress
from contextvars import ContextVar
from pathlib import Path
from typing import Any, cast
from urllib.parse import urlsplit
from urllib.request import getproxies, proxy_bypass

from dnakit._progress import Progress
from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)

_CATALOG = Path(__file__).with_name("mcp_servers.json")
_ENV_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")
_RUNTIME_ENVIRONMENT = frozenset({"PYTHONPATH", "VIRTUAL_ENV", "CONDA_PREFIX", "LANG", "LC_ALL"})
_BRIDGE_LOG_CONTEXT: ContextVar[bool] = ContextVar("dnakit_mcp_bridge_logs", default=False)


class _BridgeLogFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        # SDK parsing exceptions can contain raw server output and credentials.
        # Context-local suppression preserves diagnostics of unrelated clients.
        return not _BRIDGE_LOG_CONTEXT.get()


_BRIDGE_LOG_FILTER = _BridgeLogFilter()


@contextmanager
def _quiet_sdk_logs() -> Iterator[None]:
    for name in ("mcp.client.stdio", "mcp.client.session", "fastmcp.client.client"):
        logging.getLogger(name).addFilter(_BRIDGE_LOG_FILTER)
    token = _BRIDGE_LOG_CONTEXT.set(True)
    try:
        yield
    finally:
        _BRIDGE_LOG_CONTEXT.reset(token)


def _definitions() -> list[dict[str, Any]]:
    return cast(list[dict[str, Any]], json.loads(_CATALOG.read_text(encoding="utf-8"))["servers"])


def list_mcp_servers(*, verbose: bool = True) -> list[dict[str, Any]]:
    """Return source-backed connection recipes without importing SDKs or starting services.

    Catalog inclusion confirms an upstream startup recipe, not local availability,
    server health or scientific validation. No credentials are read by this call.
    """
    return _definitions()


def _definition(server_id: str) -> dict[str, Any]:
    if not isinstance(server_id, str):
        raise ConfigurationError("server_id must be a catalog identifier.")
    for server in _definitions():
        if server["id"] == server_id:
            return server
    raise ConfigurationError("Unknown MCP server; use list_mcp_servers().")


def _sdk() -> tuple[Any, Any, Any]:
    try:
        client = importlib.import_module("fastmcp").Client
        transports = importlib.import_module("fastmcp.client.transports")
        return client, transports.StdioTransport, transports.StreamableHttpTransport
    except ImportError:
        raise BackendUnavailableError(
            "The optional FastMCP client SDK is not installed.",
            hint="Install the DNAKit agent extra in this environment: pip install 'dnakit[agent]'.",
        ) from None


def _http_proxy_support(url: str) -> None:
    parsed = urlsplit(url)
    if parsed.hostname and proxy_bypass(parsed.hostname):
        return
    proxies = getproxies()
    proxy = proxies.get(parsed.scheme) or proxies.get("all", "")
    if proxy.lower().startswith("socks"):
        try:
            importlib.import_module("socksio")
        except ImportError:
            raise BackendUnavailableError(
                "HTTP MCP requires SOCKS support for the configured environment proxy.",
                hint="Install or upgrade 'dnakit[agent]', or explicitly install 'httpx[socks]'.",
            ) from None


def _limits(timeout_seconds: float, max_response_bytes: int, max_tools: int) -> None:
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or not math.isfinite(timeout_seconds)
        or timeout_seconds <= 0
    ):
        raise ConfigurationError("timeout_seconds must be a finite positive number.")
    for name, value in (("max_response_bytes", max_response_bytes), ("max_tools", max_tools)):
        if type(value) is not int or value <= 0:
            raise ConfigurationError(f"{name} must be a positive integer.")


def _json_size(value: Any, limit: int, *, configuration: bool = False) -> None:
    error = ConfigurationError if configuration else BackendExecutionError
    try:
        size = 0
        for chunk in json.JSONEncoder(ensure_ascii=False, allow_nan=False).iterencode(value):
            size += len(chunk.encode("utf-8"))
            if size > limit:
                raise error("MCP payload exceeds the configured response/input byte limit.")
    except (TypeError, ValueError, RecursionError):
        raise error("MCP data must be finite JSON-compatible values.") from None
    if size > limit:
        raise error("MCP payload exceeds the configured response/input byte limit.")


def _environment(
    server: dict[str, Any], names: dict[str, str] | None
) -> tuple[dict[str, str], list[str]]:
    allowed = set(server["required_environment"]) | set(server["optional_environment"])
    allowed |= _RUNTIME_ENVIRONMENT
    if names is not None and not isinstance(names, dict):
        raise ConfigurationError(
            "environment_variables maps server names to source environment variable names."
        )
    mapping = {name: name for name in allowed if name not in _RUNTIME_ENVIRONMENT}
    for target, source in (names or {}).items():
        if (
            target not in allowed
            or not isinstance(source, str)
            or _ENV_NAME.fullmatch(source) is None
        ):
            raise ConfigurationError("Invalid or unsupported named environment variable mapping.")
        mapping[target] = source
    values: dict[str, str] = {}
    secrets: list[str] = []
    for target, source in mapping.items():
        value = os.environ.get(source)
        if value:
            if "\x00" in value:
                raise ConfigurationError("Environment values cannot contain NUL characters.")
            values[target] = value
            if target in server["secret_environment"]:
                secrets.append(value)
        elif target in server["required_environment"] or target in (names or {}):
            raise ConfigurationError(f"Required environment variable is not set: {source}.")
    for first, second in server.get("environment_pairs", []):
        if (first in values) != (second in values):
            raise ConfigurationError(f"Configure both {first} and {second}, or neither.")
    # Verified transport selectors are constants, never credentials supplied by callers.
    values.update(server.get("fixed_environment", {}))
    return values, sorted(set(secrets), key=len, reverse=True)


def _redact(value: Any, secrets: list[str]) -> Any:
    if isinstance(value, str):
        for secret in secrets:
            value = value.replace(secret, "[REDACTED]")
        return value
    if isinstance(value, list):
        return [_redact(item, secrets) for item in value]
    if isinstance(value, dict):
        return {_redact(key, secrets): _redact(item, secrets) for key, item in value.items()}
    return value


def _redact_result(value: dict[str, Any], secrets: list[str]) -> dict[str, Any]:
    return cast(dict[str, Any], _redact(value, secrets))


def _connection(
    server: dict[str, Any],
    transport: str,
    executable_path: str | None,
    workdir: str | None,
    environment_variables: dict[str, str] | None,
) -> tuple[dict[str, Any], list[str]]:
    if transport not in {"stdio", "http"}:
        raise ConfigurationError("transport must be stdio or http.")
    supported = server.get("transports", ["stdio", "http"])
    if transport not in supported:
        raise ConfigurationError(
            f"The selected MCP server does not support {transport.upper()} transport.",
            hint="Select one of its catalog transports explicitly: " + ", ".join(supported),
        )
    if transport == "http":
        if not server.get("http_url"):
            raise ConfigurationError(
                "No author-verified HTTP endpoint is recorded for this server."
            )
        if executable_path is not None or workdir is not None or environment_variables:
            raise ConfigurationError("Local process configuration is only supported with stdio.")
        return {"url": server["http_url"]}, []
    program = executable_path if executable_path is not None else server["command"]
    if not isinstance(program, str) or not program or any(c in program for c in "\x00\r\n"):
        raise ConfigurationError("executable_path must name an installed executable.")
    command = shutil.which(program)
    if command is None:
        raise BackendUnavailableError(
            "The selected MCP server executable is not installed.", hint=server["installation"]
        )
    if Path(command).name.lower() in {"uv", "uvx", "bunx", "npx", "npm", "pip", "pip3", "pipx"}:
        raise ConfigurationError("Package installers/runners cannot be used as MCP executables.")
    cwd = None
    if workdir is not None:
        if not isinstance(workdir, str) or not workdir or "\x00" in workdir:
            raise ConfigurationError("workdir must be an existing directory.")
        cwd = Path(workdir).expanduser().resolve()
        if not cwd.is_dir():
            raise ConfigurationError("workdir must be an existing directory.")
    if server["required_files"]:
        if cwd is None:
            raise ConfigurationError("This source-based server requires an explicit workdir.")
        if any(not (cwd / name).is_file() for name in server["required_files"]):
            raise ConfigurationError("workdir is missing required upstream server files.")
    env, secrets = _environment(server, environment_variables)
    return {
        "command": command,
        "args": list(server["arguments"]),
        "env": env,
        "cwd": str(cwd) if cwd is not None else None,
    }, secrets


async def _discard_log(*args: Any, **kwargs: Any) -> None:
    """Do not forward untrusted server log/progress messages or credentials."""


async def _tools(client: Any, max_tools: int, max_response_bytes: int) -> list[dict[str, Any]]:
    tools: list[dict[str, Any]] = []
    cursor = None
    seen: set[str] = set()
    for _ in range(100):
        page = await client.list_tools_mcp(cursor=cursor)
        tools.extend(
            tool.model_dump(mode="json", by_alias=True, exclude_none=True) for tool in page.tools
        )
        if len(tools) > max_tools:
            raise BackendExecutionError("MCP tool list exceeds max_tools.")
        _json_size(tools, max_response_bytes)
        cursor = page.nextCursor
        if not cursor:
            names = [tool["name"] for tool in tools]
            if len(names) != len(set(names)):
                raise BackendExecutionError("MCP server returned duplicate tool names.")
            return tools
        if cursor in seen:
            raise BackendExecutionError("MCP server returned a repeated pagination cursor.")
        seen.add(cursor)
    raise BackendExecutionError("MCP tool pagination exceeds 100 pages.")


def _validate_arguments(arguments: dict[str, Any], schema: dict[str, Any]) -> None:
    from jsonschema import Draft202012Validator
    from jsonschema.exceptions import SchemaError, ValidationError
    from jsonschema.validators import validator_for
    from referencing import Registry
    from referencing.exceptions import Unresolvable

    try:
        validator = validator_for(schema, default=Draft202012Validator)
        validator.check_schema(schema)
        # Explicit empty registry prevents schema references from fetching URLs.
        validator(schema, registry=Registry()).validate(arguments)
    except ValidationError as exc:
        raise ConfigurationError(
            "Arguments do not satisfy the MCP tool inputSchema.",
            code="INVALID_MCP_ARGUMENTS",
            context={"argument_path": list(exc.path)},
        ) from None
    except (SchemaError, Unresolvable):
        raise BackendExecutionError(
            "MCP tool has an invalid or externally referenced schema."
        ) from None


def _sync(
    operation: Callable[[], Awaitable[dict[str, Any]]], timeout_seconds: float
) -> dict[str, Any]:
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        pass
    else:
        raise ConfigurationError(
            "This synchronous MCP API cannot run inside an active event loop.",
            hint="From async code use await asyncio.to_thread(the_mcp_function, ...).",
        )

    async def bounded() -> dict[str, Any]:
        return await asyncio.wait_for(operation(), timeout_seconds)

    try:
        return asyncio.run(bounded())
    except (ConfigurationError, BackendExecutionError, BackendUnavailableError):
        raise
    except (TimeoutError, asyncio.TimeoutError):
        raise BackendTimeoutError("MCP connection or invocation timed out.") from None
    except Exception as exc:
        # SDK errors may embed process environment, tool input, or HTTP headers.
        raise BackendExecutionError(
            "MCP connection or invocation failed.", context={"error_type": type(exc).__name__}
        ) from None


def _request(
    server_id: str,
    *,
    tool_name: str | None,
    arguments: dict[str, Any] | None,
    transport: str,
    executable_path: str | None,
    workdir: str | None,
    environment_variables: dict[str, str] | None,
    timeout_seconds: float,
    max_response_bytes: int,
    max_tools: int,
    verbose: bool,
    workflow_calls: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    _limits(timeout_seconds, max_response_bytes, max_tools)
    server = _definition(server_id)
    params, secrets = _connection(
        server, transport, executable_path, workdir, environment_variables
    )
    client_type, stdio_type, http_type = _sdk()
    if transport == "http":
        _http_proxy_support(params["url"])

    async def operation() -> dict[str, Any]:
        with (
            open(os.devnull, "w", encoding="utf-8") as quiet,
            Progress(disable=not verbose) as progress,
            _quiet_sdk_logs(),
        ):
            total = 1 + 2 * len(workflow_calls) if workflow_calls else (3 if tool_name else 2)
            task = progress.add_task("External MCP", total=total)
            connection = (
                stdio_type(**params, keep_alive=False, log_file=quiet)
                if transport == "stdio"
                else http_type(**params)
            )
            client = client_type(
                connection,
                timeout=timeout_seconds,
                init_timeout=timeout_seconds,
                log_handler=_discard_log,
                progress_handler=_discard_log,
            )
            async with client:
                progress.advance(task)
                if workflow_calls is not None:
                    steps: list[dict[str, Any]] = []
                    for index, call in enumerate(workflow_calls):
                        # Stateful servers may change their advertised tools after loading.
                        available = await _tools(client, max_tools, max_response_bytes)
                        progress.advance(task)
                        step = await _invoke(
                            client,
                            available,
                            call["tool_name"],
                            call["arguments"],
                            timeout_seconds,
                            max_response_bytes,
                            application_error_prefixes=server.get("application_error_prefixes", []),
                            application_error_keys=server.get("application_error_keys", []),
                        )
                        steps.append({**step, "step_index": index})
                        _json_size(steps, max_response_bytes)
                        progress.advance(task)
                        if step["isError"]:
                            break
                    failed = bool(steps[-1]["isError"])
                    return _redact_result(
                        {
                            "server_id": server_id,
                            "transport": transport,
                            "status": "remote_error" if failed else "completed",
                            "isError": failed,
                            "requested_steps": len(workflow_calls),
                            "executed_steps": len(steps),
                            "steps": steps,
                        },
                        secrets,
                    )
                tools = await _tools(client, max_tools, max_response_bytes)
                progress.advance(task)
                if tool_name is None:
                    return _redact_result(
                        {
                            "server_id": server_id,
                            "transport": transport,
                            "tool_count": len(tools),
                            "tools": tools,
                        },
                        secrets,
                    )
                payload = await _invoke(
                    client,
                    tools,
                    tool_name,
                    arguments or {},
                    timeout_seconds,
                    max_response_bytes,
                    application_error_prefixes=server.get("application_error_prefixes", []),
                    application_error_keys=server.get("application_error_keys", []),
                )
                progress.advance(task)
                return _redact_result(
                    {
                        **payload,
                        "server_id": server_id,
                        "transport": transport,
                    },
                    secrets,
                )

    return _sync(operation, timeout_seconds)


async def _invoke(
    client: Any,
    tools: list[dict[str, Any]],
    tool_name: str,
    arguments: dict[str, Any],
    timeout_seconds: float,
    max_response_bytes: int,
    *,
    application_error_prefixes: list[str] | None = None,
    application_error_keys: list[str] | None = None,
) -> dict[str, Any]:
    tool = next((tool for tool in tools if tool["name"] == tool_name), None)
    if tool is None:
        raise ConfigurationError("Tool is not advertised by the selected MCP server.")
    _validate_arguments(arguments, tool["inputSchema"])
    result = await client.call_tool_mcp(tool_name, arguments, timeout=timeout_seconds)
    payload = result.model_dump(mode="json", by_alias=True, exclude_none=True)
    _json_size(payload, max_response_bytes)
    prefixes = [
        prefix.replace("{tool_name}", tool_name) for prefix in application_error_prefixes or []
    ]
    text_error = any(
        block.get("type") == "text"
        and any(block.get("text", "").startswith(prefix) for prefix in prefixes)
        for block in payload.get("content", [])
    )
    documents = [payload.get("structuredContent")]
    if application_error_keys:
        for block in payload.get("content", []):
            if block.get("type") == "text":
                with suppress(TypeError, ValueError):
                    documents.append(json.loads(block.get("text", "")))
    json_error = any(
        isinstance(document, dict)
        and any(document.get(key) for key in application_error_keys or [])
        for document in documents
    )
    failed = bool(result.isError) or text_error or json_error
    if (text_error or json_error) and not result.isError:
        # A few source-reviewed servers encode known failures as ordinary text.
        # Preserve the original protocol flag while reporting the identified failure.
        payload["upstream_isError"] = False
        payload["error_detection"] = (
            "verified_upstream_json_error" if json_error else "verified_upstream_text_prefix"
        )
    return {
        **payload,
        "tool_name": tool_name,
        "isError": failed,
        "status": "remote_error" if failed else "completed",
    }


def list_mcp_tools(
    server_id: str,
    *,
    transport: str = "stdio",
    executable_path: str | None = None,
    workdir: str | None = None,
    environment_variables: dict[str, str] | None = None,
    timeout_seconds: float = 60.0,
    max_response_bytes: int = 4_194_304,
    max_tools: int = 5000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Connect to one catalog server and retrieve its live tool names and schemas.

    Stdio launches the separately installed executable; source servers require
    workdir. No installer or shell is used. environment_variables maps supported
    server variable names to existing environment variable names, never values.
    Only explicitly requested runtime variables (e.g. PYTHONPATH) are forwarded.
    HTTP is supported only for author-documented catalog endpoints. Response
    byte limits are checked after SDK decoding, not at the transport wire layer.
    Validated local description revisions apply only to this listing; upstream
    invocation schemas remain authoritative.
    """
    from dnakit.tools.descriptions import overlay_mcp_tools

    payload = _request(
        server_id,
        tool_name=None,
        arguments=None,
        transport=transport,
        executable_path=executable_path,
        workdir=workdir,
        environment_variables=environment_variables,
        timeout_seconds=timeout_seconds,
        max_response_bytes=max_response_bytes,
        max_tools=max_tools,
        verbose=verbose,
    )
    payload = overlay_mcp_tools(server_id, payload)
    _json_size(payload, max_response_bytes)
    return payload


def call_mcp_tool(
    server_id: str,
    tool_name: str,
    arguments: dict[str, Any],
    *,
    transport: str = "stdio",
    executable_path: str | None = None,
    workdir: str | None = None,
    environment_variables: dict[str, str] | None = None,
    timeout_seconds: float = 60.0,
    max_response_bytes: int = 4_194_304,
    max_tools: int = 5000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Validate against a server's live inputSchema, then invoke exactly one tool.

    Tool effects belong to the caller-selected upstream tool: external services
    may read/write files, contact networks or load models. Returned content and
    structuredContent retain the MCP isError flag; remote errors are explicitly
    marked remote_error. Known forwarded credential values are redacted. Each
    call starts a fresh session, so upstream in-memory state is not retained.
    """
    if not isinstance(tool_name, str) or not tool_name or len(tool_name) > 512:
        raise ConfigurationError("tool_name must be a non-empty MCP tool name.")
    if not isinstance(arguments, dict) or any(not isinstance(key, str) for key in arguments):
        raise ConfigurationError("arguments must be a JSON object with string keys.")
    _limits(timeout_seconds, max_response_bytes, max_tools)
    _json_size(arguments, max_response_bytes, configuration=True)
    return _request(
        server_id,
        tool_name=tool_name,
        arguments=arguments,
        transport=transport,
        executable_path=executable_path,
        workdir=workdir,
        environment_variables=environment_variables,
        timeout_seconds=timeout_seconds,
        max_response_bytes=max_response_bytes,
        max_tools=max_tools,
        verbose=verbose,
    )


def call_mcp_workflow(
    server_id: str,
    calls: list[dict[str, Any]],
    *,
    transport: str = "stdio",
    executable_path: str | None = None,
    workdir: str | None = None,
    environment_variables: dict[str, str] | None = None,
    timeout_seconds: float = 60.0,
    max_response_bytes: int = 4_194_304,
    max_tools: int = 5000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run explicit tool calls sequentially within one retained MCP session.

    Each call is {"tool_name": str, "arguments": dict}; 1 to 100 calls are allowed.
    Schemas are retrieved and validated before every step. Server in-memory state
    persists between steps, supporting model loading followed by prediction. Calls
    use explicit arguments; output substitution and arbitrary scripts are not run.
    A remote isError stops the workflow and returns all executed steps, including
    the error. Configuration/connection errors raise; earlier effects cannot be
    rolled back. timeout_seconds bounds the complete session, and the response
    byte limit applies to accumulated decoded step results.
    """
    _limits(timeout_seconds, max_response_bytes, max_tools)
    if not isinstance(calls, list) or not 1 <= len(calls) <= 100:
        raise ConfigurationError("calls must contain between 1 and 100 tool calls.")
    for call in calls:
        if not isinstance(call, dict) or set(call) != {"tool_name", "arguments"}:
            raise ConfigurationError("Each workflow call requires tool_name and arguments only.")
        if (
            not isinstance(call["tool_name"], str)
            or not call["tool_name"]
            or len(call["tool_name"]) > 512
        ):
            raise ConfigurationError("tool_name must be a non-empty MCP tool name.")
        if not isinstance(call["arguments"], dict) or any(
            not isinstance(key, str) for key in call["arguments"]
        ):
            raise ConfigurationError("arguments must be a JSON object with string keys.")
    _json_size(calls, max_response_bytes, configuration=True)
    return _request(
        server_id,
        tool_name=None,
        arguments=None,
        workflow_calls=calls,
        transport=transport,
        executable_path=executable_path,
        workdir=workdir,
        environment_variables=environment_variables,
        timeout_seconds=timeout_seconds,
        max_response_bytes=max_response_bytes,
        max_tools=max_tools,
        verbose=verbose,
    )


__all__ = ["call_mcp_tool", "call_mcp_workflow", "list_mcp_servers", "list_mcp_tools"]
_register_public_verbose(globals(), tuple(__all__))
