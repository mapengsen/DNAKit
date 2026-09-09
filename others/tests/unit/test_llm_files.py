"""Offline tests for explicit, bounded local-file LLM input."""

from __future__ import annotations

import asyncio
import io
import json
import os
import runpy
import threading
import traceback
import urllib.request
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.llm import LLMConfig, acall_llm, call_llm, read_llm_file
from dnakit.llm import files as file_module


def _mock_api(monkeypatch: pytest.MonkeyPatch) -> list[urllib.request.Request]:
    requests: list[urllib.request.Request] = []

    def capture(request: urllib.request.Request, *, timeout: float) -> io.BytesIO:
        requests.append(request)
        return io.BytesIO(
            json.dumps(
                {
                    "choices": [
                        {
                            "message": {"content": '{"summary": "Processed file"}'},
                            "finish_reason": "stop",
                        }
                    ],
                }
            ).encode()
        )

    monkeypatch.setattr(urllib.request, "build_opener", lambda *args: SimpleNamespace(open=capture))
    return requests


@pytest.mark.parametrize("provider", ["deepseek", "gemini", "qwen", "kimi", "openrouter", "openai"])
def test_file_content_not_path_is_sent(
    provider: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "LICENSE"
    path.write_text("A sample license.\nPreserve this line.\n", encoding="utf-8")
    (tmp_path / "neighbor-secret.txt").write_text("Never read adjacent files", encoding="utf-8")
    requests = _mock_api(monkeypatch)
    result = call_llm(
        "Summarize the file as JSON.",
        file_path=path,
        config=LLMConfig(provider=provider, model="selected-model", api_key="test-key"),
        output_format="json",
    )
    assert result.data == {"summary": "Processed file"}
    assert len(requests) == 1
    assert isinstance(requests[0].data, bytes)
    payload = json.loads(requests[0].data)
    serialized = payload["messages"][-1]["content"].split("\n", 1)[1]
    assert json.loads(serialized) == {
        "file_name": "LICENSE",
        "content": "A sample license.\nPreserve this line.\n",
    }
    assert str(tmp_path) not in json.dumps(payload)
    assert "Never read adjacent files" not in json.dumps(payload)
    assert "test-key" not in json.dumps(payload)


@pytest.mark.parametrize(
    ("filename", "encoding"),
    [
        ("LICENSE", "utf-8"),
        ("notes.md", "utf-8-sig"),
        ("data.csv", "gb18030"),
        ("records.json", "utf-16"),
        ("script.py", "utf-32"),
    ],
)
def test_local_preview_text_formats_and_encodings(
    filename: str, encoding: str, tmp_path: Path
) -> None:
    path = tmp_path / filename
    content = "中文数据\nsecond line\tvalue\r\n"
    path.write_bytes(content.encode(encoding))
    selected_encoding = "utf-8-sig" if encoding == "utf-8" else encoding
    assert read_llm_file(path, encoding=selected_encoding) == content


def test_relative_and_space_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    path = tmp_path / "file with spaces.txt"
    path.write_text("contents", encoding="utf-8")
    assert read_llm_file("file with spaces.txt", max_bytes=8) == "contents"
    with pytest.raises(InputFormatError, match="LLM_FILE_LIMIT"):
        read_llm_file("file with spaces.txt", max_bytes=7)


@pytest.mark.parametrize("data", [b"", b" \n\t", b"\xef\xbb\xbf"])
def test_empty_files_are_rejected(data: bytes, tmp_path: Path) -> None:
    path = tmp_path / "empty"
    path.write_bytes(data)
    with pytest.raises(InputFormatError, match="EMPTY_LLM_FILE"):
        read_llm_file(path)


@pytest.mark.parametrize(
    "data",
    [
        b"%PDF-1.7\nASCII-only fake PDF",
        b"PK\x03\x04fake-office-file",
        b"\x89PNG\r\n\x1a\n",
        b"GIF89a",
        b"\x1f\x8bcompressed",
        b"binary\x00data",
        b"log\x1b[31mred",
    ],
)
def test_binary_files_are_rejected_even_with_permissive_encoding(
    data: bytes, tmp_path: Path
) -> None:
    path = tmp_path / "binary"
    path.write_bytes(data)
    with pytest.raises(InputFormatError, match="UNSUPPORTED_LLM_FILE_FORMAT"):
        read_llm_file(path, encoding="latin-1")


@pytest.mark.parametrize(
    "file_path", ["", "   ", "https://example.com/file", "file:///tmp/file", "bad\x00path", 3]
)
def test_invalid_file_paths(file_path: Any) -> None:
    with pytest.raises(InputFormatError, match="INVALID_LLM_FILE_PATH"):
        read_llm_file(file_path)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_bytes": 0},
        {"max_bytes": True},
        {"max_bytes": 1.5},
        {"encoding": "no-such-encoding"},
        {"encoding": "base64_codec"},
        {"encoding": None},
    ],
)
def test_invalid_reader_options(kwargs: dict[str, Any], tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="INVALID_LLM_FILE_CONFIG"):
        read_llm_file(tmp_path / "not-opened", **kwargs)


def test_directory_missing_and_bad_encoding_fail_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests = _mock_api(monkeypatch)
    configured = LLMConfig(model="m", api_key="test-key")
    with pytest.raises(InputFormatError, match="INVALID_LLM_FILE_TYPE"):
        call_llm("Task", file_path=tmp_path, config=configured)
    with pytest.raises(InputFormatError, match="LLM_FILE_NOT_FOUND"):
        call_llm("Task", file_path=tmp_path / "missing", config=configured)
    path = tmp_path / "encoded.txt"
    path.write_bytes("中文".encode("gb18030"))
    with pytest.raises(InputFormatError, match="LLM_FILE_ENCODING_ERROR"):
        call_llm("Task", file_path=path, config=configured)
    assert not requests


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFOs are unavailable on this platform")
def test_fifo_is_rejected_without_waiting_for_writer(tmp_path: Path) -> None:
    fifo = tmp_path / "pipe"
    os.mkfifo(fifo)
    with pytest.raises(InputFormatError, match="INVALID_LLM_FILE_TYPE"):
        read_llm_file(fifo)


def test_local_read_errors_do_not_echo_file_paths(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    path = tmp_path / "private-file-name"

    def unreadable(*args: object, **kwargs: object) -> Any:
        raise PermissionError("private-file-name and other sensitive details")

    with monkeypatch.context() as patch:
        patch.setattr(Path, "stat", unreadable)
        with pytest.raises(InputFormatError, match="LLM_FILE_READ_ERROR") as caught:
            read_llm_file(path)
    assert "private-file-name" not in "".join(traceback.format_exception(caught.value))


def test_file_and_request_limits_are_independent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests = _mock_api(monkeypatch)
    path = tmp_path / "large.txt"
    path.write_text("a" * 1000, encoding="utf-8")
    with pytest.raises(InputFormatError, match="LLM_FILE_LIMIT"):
        call_llm(
            "Task", file_path=path, max_file_bytes=999, config=LLMConfig(model="m", api_key="k")
        )
    with pytest.raises(InputFormatError, match="LLM_REQUEST_LIMIT"):
        call_llm(
            "Task",
            file_path=path,
            max_file_bytes=1000,
            config=LLMConfig(model="m", api_key="k", max_request_bytes=1000),
        )
    assert not requests


@pytest.mark.parametrize("input_data", ["text", {}, [], False, 0])
def test_data_and_file_are_mutually_exclusive(input_data: object, tmp_path: Path) -> None:
    # Nonexistent source confirms the conflict is detected before opening a file.
    with pytest.raises(ConfigurationError, match="either input_data or file_path"):
        call_llm(
            "Task",
            input_data=input_data,
            file_path=tmp_path / "missing",
            config=LLMConfig(model="m"),
        )


def test_async_file_read_and_call_run_outside_event_loop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    requests = _mock_api(monkeypatch)
    path = tmp_path / "utf16.txt"
    path.write_bytes("中文".encode("utf-16"))
    main_thread = threading.get_ident()
    reading_threads: list[int] = []
    original_opener = file_module._nonblocking_open

    def capture_open(path: str, flags: int) -> int:
        reading_threads.append(threading.get_ident())
        return original_opener(path, flags)

    monkeypatch.setattr(file_module, "_nonblocking_open", capture_open)
    result = asyncio.run(
        acall_llm(
            "Summarize",
            file_path=path,
            file_encoding="utf-16",
            max_file_bytes=100,
            config=LLMConfig(provider="gemini", model="m", api_key="test-key"),
            output_format="json",
        )
    )
    assert result.data == {"summary": "Processed file"}
    assert len(requests) == 1
    assert reading_threads and all(item != main_thread for item in reading_threads)


def test_real_repository_license_is_read_locally() -> None:
    source = Path(__file__).resolve().parents[3] / "LICENSE"
    content = read_llm_file(source)
    assert content.startswith("MIT License")
    assert content == source.read_bytes().decode("utf-8-sig")


def test_file_example_submits_only_one_request(monkeypatch: pytest.MonkeyPatch) -> None:
    requests = _mock_api(monkeypatch)
    monkeypatch.setenv("DNAKIT_LLM_API_KEY", "example-test-key")
    monkeypatch.setenv("DNAKIT_LLM_MODEL", "selected-model")
    source = Path(__file__).resolve().parents[2] / "examples" / "custom_llm_file.py"
    runpy.run_path(str(source), run_name="__main__")
    assert len(requests) == 1
    assert isinstance(requests[0].data, bytes)
    assert "MIT License" in requests[0].data.decode()
