"""Offline checks for generic web access, real link extraction and search failures."""

from __future__ import annotations

import io
import json
from email.message import Message
from typing import Any
from unittest.mock import Mock
from urllib.request import Request

import pytest

from dnakit.exceptions import ConfigurationError, QueryError
from dnakit.search import SearchConfig, read_webpage, web, web_search
from dnakit.tools import default_tool_registry


class _Response(io.BytesIO):
    def __init__(self, body: bytes, url: str, content_type: str, status: int) -> None:
        super().__init__(body)
        self.url = url
        self.status = status
        self.headers = Message()
        self.headers["Content-Type"] = content_type

    def geturl(self) -> str:
        return self.url


def _serve(
    monkeypatch: pytest.MonkeyPatch,
    body: str,
    *,
    url: str = "https://files.example/repo/tree/main/",
    content_type: str = "text/html; charset=utf-8",
    status: int = 200,
) -> Mock:
    opener = Mock()
    opener.open.side_effect = lambda *args, **kwargs: _Response(
        body.encode(), url, content_type, status
    )
    monkeypatch.setattr(web, "build_opener", Mock(return_value=opener))
    return opener


def test_read_page_resolves_observed_links_and_filters_scripts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _serve(
        monkeypatch,
        """<title>Repository</title><script>invented.parquet</script>
      <style>.hidden{}</style><p>H2AFZ train and test</p>
      <a href="../H2AFZ/train.parquet?download=true&amp;x=1">train</a>
      <a href="https://files.example/test.parquet">test</a>
      <a href="javascript:alert(1)">ignore</a><a href="file:///tmp/x">ignore</a>""",
    )
    result = read_webpage("https://files.example/start", verbose=False)
    assert result.request_url == "https://files.example/repo/tree/main/"
    assert [row["url"] for row in result.records] == [
        "https://files.example/repo/tree/H2AFZ/train.parquet?download=true&x=1",
        "https://files.example/test.parquet",
    ]
    assert "invented" not in result.metadata["text"]
    assert result.metadata["title"] == "Repository"


def test_page_limits_have_continuations_and_match(monkeypatch: pytest.MonkeyPatch) -> None:
    _serve(
        monkeypatch,
        "<p>H2AFZ files</p><p>other</p>"
        + "".join(f'<p><a href="/H2AFZ/{i}.csv">H2AFZ {i}</a></p>' for i in range(4)),
    )
    first = read_webpage("https://files.example", match="H2AFZ", max_characters=5, max_links=2)
    assert first.metadata["text_truncated"]
    assert first.metadata["next_offset"] == 5
    assert first.metadata["next_link_offset"] == 2
    assert len(first.records) == 2
    second = read_webpage("https://files.example", match="H2AFZ", offset=5, link_offset=2)
    assert not second.metadata["links_truncated"]
    assert len(second.records) == 2
    assert "other" not in second.metadata["text"]


def test_json_api_is_read_without_inventing_download_paths(monkeypatch: pytest.MonkeyPatch) -> None:
    _serve(
        monkeypatch,
        '{"name":"train.parquet","url":"https:\\/\\/files.example/train"}',
        content_type="application/json",
    )
    result = read_webpage("https://files.example/api")
    assert result.records[0]["url"] == "https://files.example/train"
    assert '"name":"train.parquet"' in result.metadata["text"]


@pytest.mark.parametrize(
    "url",
    [
        "http://files.example",
        "file:///tmp/a",
        "https://user:pass@files.example",
        "https://files.example:bad",
        "https://files.example/a\nb",
    ],
)
def test_invalid_urls_fail_before_network(monkeypatch: pytest.MonkeyPatch, url: str) -> None:
    opener = _serve(monkeypatch, "unused")
    with pytest.raises(ConfigurationError, match="INVALID_WEB_URL"):
        read_webpage(url)
    opener.open.assert_not_called()


def test_redirects_are_validated_before_following() -> None:
    handler = web._WebRedirects()
    request = Request("https://files.example/start")
    with pytest.raises(ConfigurationError, match="INVALID_WEB_URL"):
        handler.redirect_request(request, None, 302, "Found", {}, "http://files.example/file")
    request.add_header("Authorization", "Bearer test-only")
    with pytest.raises(QueryError, match="WEB_REDIRECT"):
        handler.redirect_request(request, None, 302, "Found", {}, "https://other.example/")


def test_binary_and_oversized_responses_are_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    _serve(monkeypatch, "PAR1 binary", content_type="application/octet-stream")
    with pytest.raises(QueryError, match="WEB_BINARY_CONTENT"):
        read_webpage("https://files.example/train.parquet")
    _serve(monkeypatch, "0123456789")
    with pytest.raises(QueryError, match="QUERY_RESPONSE_SIZE_LIMIT"):
        read_webpage("https://files.example", config=SearchConfig(max_response_bytes=5))


def test_duckduckgo_unwraps_actual_search_links(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("TAVILY_API_KEY", raising=False)
    _serve(
        monkeypatch,
        '<a class="result__a" href="//duckduckgo.com/l/?uddg=https%3A%2F%2Ffiles.example%2Frepo">'
        "Dataset revised</a>",
        url="https://html.duckduckgo.com/html/",
    )
    result = web_search("dataset revised", verbose=False)
    assert result.provider == "duckduckgo"
    assert result.records[0]["url"] == "https://files.example/repo"
    assert result.metadata["verified_downloads"] is False


@pytest.mark.parametrize(
    "body,status,code",
    [
        ('<form id="challenge-form">', 202, "WEB_SEARCH_BLOCKED"),
        ("unexpected markup", 200, "WEB_SEARCH_RESPONSE"),
    ],
)
def test_search_failures_do_not_masquerade_as_empty_results(
    monkeypatch: pytest.MonkeyPatch, body: str, status: int, code: str
) -> None:
    _serve(monkeypatch, body, status=status)
    with pytest.raises(QueryError, match=code):
        web_search("dataset", provider="duckduckgo")


def test_tavily_uses_environment_key_without_disclosing_it(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("TAVILY_API_KEY", "test-only-search-key")
    opener = _serve(
        monkeypatch,
        json.dumps(
            {
                "results": [
                    {
                        "url": "https://files.example/repo",
                        "title": "Dataset",
                        "content": "Description",
                    }
                ]
            }
        ),
        content_type="application/json",
        url="https://api.tavily.com/search",
    )
    result = web_search("dataset", limit=3)
    request: Request = opener.open.call_args.args[0]
    assert request.get_header("Authorization") == "Bearer test-only-search-key"
    assert json.loads(request.data)["max_results"] == 3
    assert "test-only-search-key" not in str(result.to_dict())
    assert result.provider == "tavily"


def test_web_tools_are_callable_from_default_catalog(monkeypatch: pytest.MonkeyPatch) -> None:
    _serve(monkeypatch, '<a href="/file.csv">download</a>')
    registry = default_tool_registry(verbose=False)
    name = "dnakit.search.read_webpage"
    assert registry.describe(name)["effect"] == "network"
    result: Any = registry.execute(name, {"url": "https://files.example"})
    assert result["records"][0]["url"] == "https://files.example/file.csv"
    assert registry.describe("dnakit.search.web_search")["effect"] == "network"
