"""The Chat tab: tool-calling loop against a local stand-in for an OpenAI-compatible API."""

from __future__ import annotations

import json
import threading
from collections.abc import Callable, Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pytest

from dnakit.gui import GuiSettings

pytest.importorskip("gradio")

from dnakit.gui import _chat
from dnakit.gui._backend import GuiBackend

KEY = "sk-test-0123456789"


class FakePlatform:
    """Serves ``/v1/models`` and scripted ``/v1/chat/completions`` replies on 127.0.0.1."""

    def __init__(self) -> None:
        self.replies: list[tuple[int, list[dict[str, Any]] | dict[str, Any]]] = []
        self.requests: list[dict[str, Any]] = []
        self.headers: list[dict[str, str]] = []
        platform = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args: Any) -> None:
                pass

            def _send(self, status: int, body: bytes, content_type: str) -> None:
                self.send_response(status)
                self.send_header("Content-Type", content_type)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self) -> None:
                platform.headers.append(dict(self.headers))
                models = [{"id": "fake-chat"}, {"id": "text-embedding-fake"}, {"id": "fake-pro"}]
                self._send(200, json.dumps({"data": models}).encode(), "application/json")

            def do_POST(self) -> None:
                length = int(self.headers.get("Content-Length", "0"))
                platform.requests.append(json.loads(self.rfile.read(length)))
                platform.headers.append(dict(self.headers))
                status, reply = platform.replies.pop(0)
                if isinstance(reply, dict):
                    self._send(status, json.dumps(reply).encode(), "application/json")
                    return
                lines = [f"data: {json.dumps(chunk)}\n\n" for chunk in reply] + ["data: [DONE]\n\n"]
                self._send(status, "".join(lines).encode(), "text/event-stream")

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.base_url = f"http://127.0.0.1:{self.server.server_address[1]}/v1"
        threading.Thread(target=self.server.serve_forever, daemon=True).start()

    def close(self) -> None:
        self.server.shutdown()
        self.server.server_close()


def delta(**fields: Any) -> dict[str, Any]:
    return {"choices": [{"index": 0, "delta": fields}]}


def tool_call(index: int, call_id: str, name: str, *fragments: str) -> list[dict[str, Any]]:
    """Stream one tool call the way platforms do: id and name first, arguments in pieces."""

    first = {"index": index, "id": call_id, "type": "function", "function": {"name": name}}
    chunks = [delta(tool_calls=[first])]
    chunks += [
        delta(tool_calls=[{"index": index, "function": {"arguments": piece}}])
        for piece in fragments
    ]
    return chunks


@pytest.fixture
def platform() -> Iterator[FakePlatform]:
    fake = FakePlatform()
    yield fake
    fake.close()


@pytest.fixture
def make_agent(tmp_path: Path) -> Iterator[Callable[[str], Any]]:
    backends: list[GuiBackend] = []

    def make(mode: str) -> Any:
        folder = tmp_path / mode
        settings = (
            GuiSettings.public_demo(folder) if mode == "public" else GuiSettings.local(folder)
        )
        backends.append(GuiBackend(settings))
        return _chat.ChatAgent(backends[-1])

    yield make
    for backend in backends:
        backend.worker.stop()


def ask(agent: Any, platform: FakePlatform, text: str, **options: Any) -> tuple[Any, ...]:
    """Send one message through a custom platform and return the final outputs of the tab."""

    outputs = list(
        agent.respond(
            text,
            options.get("shown", []),
            options.get("history", []),
            "custom",
            platform.base_url,
            options.get("key", KEY),
            options.get("model", "fake-chat"),
            options.get("allow", False),
            None,
        )
    )
    assert outputs[0][4] == "", "the first update empties the message box"
    return tuple(outputs[-1])


def test_chat_runs_the_tools_the_model_asks_for(
    make_agent: Callable[[str], Any], platform: FakePlatform
) -> None:
    agent = make_agent("local")
    call = '{"tool_name": "dnakit.descriptors.cpg_features", "arguments": {"value": "ACGCGT"}}'
    platform.replies = [
        (
            200,
            [
                delta(reasoning_content="I need the CpG count."),
                *tool_call(0, "call_a", "search_dnakit_tools", '{"query": ', '"cpg island"}'),
                *tool_call(1, "call_b", "call_dnakit_tool", call[:30], call[30:]),
            ],
        ),
        (200, [delta(content="The sequence has "), delta(content="2 CpG sites."), {"usage": {}}]),
    ]

    _, rendered, shown, history, _, _, _, status = ask(agent, platform, "Count CpG in ACGCGT")

    first, second = platform.requests
    assert first["model"] == "fake-chat" and first["stream"] is True
    assert [tool["function"]["name"] for tool in first["tools"]][:3] == [
        "search_dnakit_tools",
        "describe_dnakit_tool",
        "call_dnakit_tool",
    ]
    assert first["messages"][0]["role"] == "system"
    assert "example.fasta" in first["messages"][0]["content"]
    assert platform.headers[0]["Authorization"] == f"Bearer {KEY}"
    # The second request carries the assistant's tool calls (with its reasoning) and results.
    assistant, searched, called = second["messages"][2:5]
    assert assistant["reasoning_content"] == "I need the CpG count."
    assert [c["id"] for c in assistant["tool_calls"]] == ["call_a", "call_b"]
    assert searched["tool_call_id"] == "call_a" and called["tool_call_id"] == "call_b"
    found = json.loads(searched["content"])
    assert found["ok"] and all(tool["runs_here"] for tool in found["tools"][:1])
    assert json.loads(called["content"])["result"]["cpg_count"] == 2
    assert shown[-1] == {"role": "assistant", "content": "The sequence has 2 CpG sites."}
    assert [m["metadata"]["title"] for m in shown if "metadata" in m] == [
        "💭 Thinking",
        "🔎 Search functions: cpg island",
        "🔧 Run: dnakit.descriptors.cpg_features",
    ]
    assert rendered == shown
    assert [m["role"] for m in history] == ["user", "assistant", "tool", "tool", "assistant"]
    assert status.startswith("`fake-chat` · 2 tool calls")

    # The next message sends the conversation so far, without the earlier reasoning.
    platform.replies = [(200, [delta(content="Yes.")])]
    ask(agent, platform, "Sure?", shown=shown, history=history)
    resent = platform.requests[-1]["messages"]
    assert [m["role"] for m in resent] == [
        "system",
        "user",
        "assistant",
        "tool",
        "tool",
        "assistant",
        "user",
    ]
    assert all("reasoning_content" not in message for message in resent)


def test_chat_keeps_side_effects_behind_the_checkbox_and_out_of_the_web_app(
    make_agent: Callable[[str], Any], tmp_path: Path
) -> None:
    local, public = make_agent("local"), make_agent("public")
    name = "dnakit.replication_repair.plot"
    arguments = {"tool_name": name, "arguments": {"input": "ACGTTGCAACGTGGCC"}}
    folder = tmp_path / "local"

    refused, _ = local.run_tool("call_dnakit_tool", arguments, folder, False)
    assert not refused["ok"] and "tick 'Allow tools with side effects'" in refused["error"]
    drawn, svg = local.run_tool("call_dnakit_tool", arguments, folder, True)
    assert drawn["ok"] and svg is not None and svg.startswith("<svg")
    assert "<svg" not in json.dumps(drawn) and "[SVG figure" in json.dumps(drawn)

    blocked, _ = public.run_tool("call_dnakit_tool", arguments, tmp_path / "public", True)
    assert not blocked["ok"] and "Not executed in the web app" in blocked["error"]
    outside = {"tool_name": "dnakit.fasta_fastq.count_records", "arguments": {"path": "/etc/hosts"}}
    confined, _ = public.run_tool("call_dnakit_tool", outside, tmp_path / "public", False)
    assert not confined["ok"] and "session folder" in confined["error"]
    described, _ = public.run_tool("describe_dnakit_tool", {"tool_name": name}, folder, False)
    assert described["tool"]["runs_here"] is False and "output_schema" not in described["tool"]
    unknown, _ = public.run_tool("call_dnakit_tool", {"tool_name": "dnakit.nope"}, folder, False)
    assert not unknown["ok"] and "AGENT_TOOL_NOT_FOUND" in unknown["error"]


def test_chat_shows_figures_and_reads_skills(
    make_agent: Callable[[str], Any], tmp_path: Path
) -> None:
    import gradio as gr

    agent = make_agent("local")
    skills, _ = agent.run_tool(
        "search_dnakit_skills", {"query": "variant annotation"}, Path(), False
    )
    assert skills["ok"] and skills["skills"]
    page, _ = agent.run_tool(
        "read_dnakit_skill", {"skill_id": skills["skills"][0]["id"]}, Path(), False
    )
    assert page["ok"] and len(page["content"]) <= _chat.SKILL_PAGE_CHARS

    figure = tmp_path / "figure.svg"
    figure.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="4" height="4"/>')
    shown = [{"role": "user", "content": "plot"}, {"role": "assistant", "figure": str(figure)}]
    rendered = _chat.render(shown)
    assert rendered[0] == shown[0] and isinstance(rendered[1]["content"], gr.Image)


def test_chat_reports_platform_errors_without_the_key(
    make_agent: Callable[[str], Any], platform: FakePlatform
) -> None:
    agent = make_agent("local")
    platform.replies = [(401, {"error": {"message": f"Incorrect API key provided: {KEY}"}})]

    _, _, shown, history, _, _, _, status = ask(agent, platform, "hello")

    assert "HTTP 401" in shown[-1]["content"] and "Check the API key" in shown[-1]["content"]
    assert KEY not in shown[-1]["content"]
    assert [m["role"] for m in history] == ["user"]
    assert "0 tool calls" in status


def test_chat_ends_with_an_answer_when_the_tool_rounds_are_used_up(
    make_agent: Callable[[str], Any], platform: FakePlatform, monkeypatch: pytest.MonkeyPatch
) -> None:
    agent = make_agent("local")
    monkeypatch.setattr(_chat, "MAX_TOOL_ROUNDS", 2)
    listing = tool_call(0, "call_a", "list_session_files", "{}")
    platform.replies = [(200, listing), (200, [delta(content="Here is what I have.")])]

    _, _, shown, history, _, _, _, _ = ask(agent, platform, "keep going")

    last = platform.requests[-1]
    assert last["tool_choice"] == "none"
    assert last["messages"][-1] == {"role": "user", "content": _chat.LAST_ROUND_NOTE}
    assert shown[-1]["content"] == "Here is what I have."
    assert history[-1] == {"role": "assistant", "content": "Here is what I have."}


def test_model_list_and_platform_addresses(
    make_agent: Callable[[str], Any], platform: FakePlatform
) -> None:
    import gradio as gr

    agent = make_agent("local")
    update, note = agent.load_models("custom", platform.base_url, KEY, None)
    assert update["choices"] == ["fake-chat", "fake-pro"] and update["value"] == "fake-chat"
    assert "2 models loaded" in note
    assert agent.load_models("deepseek", "", "", None)[0] == gr.update()

    assert _chat.resolve_base_url("deepseek", "ignored", public=True) == "https://api.deepseek.com"
    local = _chat.resolve_base_url("custom", platform.base_url + "/chat/completions/", public=False)
    assert local == platform.base_url
    for address in (platform.base_url, "https://127.0.0.1/v1", "https://localhost/v1", "ftp://x"):
        with pytest.raises(_chat.ChatError):
            _chat.resolve_base_url("custom", address, public=True)
    with pytest.raises(_chat.ChatError):
        _chat.resolve_base_url("custom", "https://user:pw@example.com/v1", public=False)

    with pytest.raises(_chat.ChatError, match="Could not reach 127"):
        _chat.list_models("http://127.0.0.1:9/v1", KEY)  # nothing listens on the discard port

    public = make_agent("public")
    with pytest.raises(gr.Error, match="API key first"):
        next(public.respond("hi", [], [], "custom", platform.base_url, "", "m", False, None))


def test_long_conversations_drop_old_tool_results_first() -> None:
    messages: list[dict[str, Any]] = [
        {"role": "user", "content": "first"},
        {"role": "assistant", "content": "", "tool_calls": [{"id": "a"}]},
        {"role": "tool", "tool_call_id": "a", "content": "x" * 5000},
        {"role": "assistant", "content": "done"},
        {"role": "user", "content": "second"},
    ]
    _chat.fit_history(messages, limit=1000)
    assert messages[2]["content"] == _chat.OMITTED_RESULT and len(messages) == 5

    messages[3]["content"] = "y" * 5000
    _chat.fit_history(messages, limit=1000)
    assert messages == [{"role": "user", "content": "second"}]


def test_search_lists_built_in_functions_first_and_remembers_missing_packages(
    make_agent: Callable[[str], Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    agent = make_agent("public")
    name = "dnakit.sequence.pyfastx_count_sequences"
    query = {"query": "count sequences in a fasta file", "limit": 20}

    before, _ = agent.run_tool("search_dnakit_tools", query, tmp_path, False)
    kinds = [tool["implementation"] for tool in before["tools"] if tool["runs_here"]]
    assert kinds == sorted(kinds, key=lambda kind: kind != "built in")
    runnable = [tool["runs_here"] for tool in before["tools"]]
    assert runnable == sorted(runnable, reverse=True)
    wrapper = next(tool for tool in before["tools"] if tool["name"] == name)
    assert wrapper["runs_here"] and wrapper["implementation"] == "optional package"

    error = (
        "BackendUnavailableError: [BACKEND_UNAVAILABLE] pyfastx could not be imported: "
        "No module named 'pyfastx'"
    )
    monkeypatch.setattr(agent.backend.worker, "call", lambda *args, **kwargs: (False, error))
    call = {"tool_name": name, "arguments": {"path": "example.fasta"}}
    failed, _ = agent.run_tool("call_dnakit_tool", call, tmp_path, False)
    assert not failed["ok"]

    after, _ = agent.run_tool("search_dnakit_tools", query, tmp_path, False)
    wrapper = next(tool for tool in after["tools"] if tool["name"] == name)
    assert not wrapper["runs_here"] and "optional package pyfastx" in wrapper["note"]
    described, _ = agent.run_tool("describe_dnakit_tool", {"tool_name": name}, tmp_path, False)
    assert described["tool"]["runs_here"] is False and "pyfastx" in described["tool"]["note"]
