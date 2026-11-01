"""Offline fixed-route, JSON, size/error, pagination and actual ZIP-byte contracts."""

from __future__ import annotations

import inspect
import io
import json
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlsplit

import pytest

from dnakit.databases import _network_p44 as transport
from dnakit.databases import public_apis_p44_functions as api
from dnakit.download import files
from dnakit.exceptions import ConfigurationError, QueryError


class Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url
        self.headers = {"Date": "Tue, 20 Oct 2026 00:00:00 GMT", "Content-Length": str(len(data))}

    def geturl(self):
        return self.url


def args_for(function, tmp_path):
    values = {}
    for name, param in inspect.signature(function).parameters.items():
        if param.default is not inspect.Parameter.empty:
            continue
        if name == "query":
            values[name] = {}
        elif name == "output_path":
            values[name] = str(tmp_path / "package.zip")
        else:
            values[name] = "id with/slash?value"
    values["verbose"] = False
    return values


@pytest.mark.parametrize("name", api.__all__)
@pytest.mark.parametrize("payload", [{"records": [{"id": "x"}], "next_page_token": "NEXT"}, [], {}])
def test_every_explicit_operation_preserves_upstream_payload(monkeypatch, tmp_path, name, payload):
    function = getattr(api, name)
    kwargs = args_for(function, tmp_path)
    seen = []

    def opener(request, **kw):
        seen.append(request)
        return Response(json.dumps(payload).encode(), request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    monkeypatch.setattr(
        files, "urlopen", lambda request, **kw: Response(b"PK\x03\x04ZIP fixture", request.full_url)
    )
    result = function(**kwargs)
    if "output_path" in kwargs:
        assert (tmp_path / "package.zip").read_bytes() == b"PK\x03\x04ZIP fixture"
        assert result["file"]["byte_size"] == 15
    else:
        assert result["data"] == payload
        assert seen[0].full_url.startswith(
            ("https://rest.ensembl.org/", "https://api.ncbi.nlm.nih.gov/datasets/v2/")
        )
        if "id with/slash?value" in kwargs.values():
            assert "%2F" in seen[0].full_url
        if "query" in kwargs:
            assert seen[0].method == "POST"
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("name", api.__all__)
def test_unknown_arguments_rejected_before_request(monkeypatch, tmp_path, name):
    function = getattr(api, name)
    kwargs = args_for(function, tmp_path)
    kwargs["query" if "query" in kwargs else "parameters"] = {"unregistered_parameter": "x"}
    with pytest.raises(ConfigurationError):
        function(**kwargs)


@pytest.mark.parametrize("status_data", [b'{"error":"not found"}', b"not JSON"])
def test_bad_response_not_reported_completed(monkeypatch, status_data):
    monkeypatch.setattr(
        transport, "urlopen", lambda request, **kw: Response(status_data, request.full_url)
    )
    with pytest.raises(QueryError):
        api.ensembl_info_analysis("human", verbose=False)


def test_http_error_and_size_limit(monkeypatch):
    def fail(request, **kw):
        raise HTTPError(request.full_url, 503, "busy", {}, None)

    monkeypatch.setattr(transport, "urlopen", fail)
    with pytest.raises(QueryError):
        api.ensembl_info_analysis("human", verbose=False)
    monkeypatch.setattr(
        transport, "urlopen", lambda request, **kw: Response(b'[{"too":"large"}]', request.full_url)
    )
    with pytest.raises(QueryError):
        api.ensembl_info_analysis("human", max_bytes=3, verbose=False)


def test_empty_204_style_body_is_valid_no_hits(monkeypatch):
    monkeypatch.setattr(transport, "urlopen", lambda request, **kw: Response(b"", request.full_url))
    assert transport._request("https://search.rcsb.org/rcsbsearch/v2/query")["data"] is None


@pytest.mark.parametrize(
    "by,fragment",
    [("accession", "/accession/"), ("assembly_name", "/assembly_name/"), ("wgs", "/wgs/")],
)
def test_genome_identifier_selectors(monkeypatch, by, fragment):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b"{}", request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    api.datasets_genome_report("x", by=by, verbose=False)
    assert fragment in seen[0]


@pytest.mark.parametrize("by", ["id", "accession", "symbol", "locus_tag"])
def test_gene_product_identifier_selectors(monkeypatch, by):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b"{}", request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    api.datasets_gene_product_report("TP53", by=by, taxon="9606", verbose=False)
    assert "/gene/" + by + "/" in seen[0]
    if by == "symbol":
        assert "/taxon/9606/" in seen[0]


def test_page_parameters_use_official_query_names(monkeypatch):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b'{"next_page_token":"more"}', request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = api.datasets_genome_report(
        "GCF_000001405.40", parameters={"page_size": 1, "page_token": "t"}, verbose=False
    )
    assert parse_qs(urlsplit(seen[0]).query) == {"page_size": ["1"], "page_token": ["t"]}
    assert result["data"]["next_page_token"] == "more"
