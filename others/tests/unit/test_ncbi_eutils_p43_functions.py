"""Offline E-utilities contracts; no accession data are fetched in these tests."""

from __future__ import annotations

import io
import json
from urllib.parse import parse_qs, urlsplit

import pytest

from dnakit.databases import _network_p44 as transport
from dnakit.databases import ncbi_eutils_p43_functions as api
from dnakit.exceptions import ConfigurationError, QueryError


class Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url
        self.headers = {"Date": "Tue, 20 Oct 2026 00:00:00 GMT"}

    def geturl(self):
        return self.url


@pytest.fixture(autouse=True)
def no_wait(monkeypatch):
    monkeypatch.setattr(api.time, "sleep", lambda value: None)


@pytest.mark.parametrize(
    "database,rettype", [(db, rt) for db, values in api._FETCH.items() for rt in values]
)
def test_fetch_modes_preserve_response(monkeypatch, database, rettype):
    seen = []
    payload = b'<?xml version="1.0"?><records><id>1</id></records>'

    def opener(request, **kwargs):
        seen.append(request.full_url)
        return Response(payload, request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = api.ncbi_fetch_dna_records(
        database, ["1", "2"], rettype=rettype, retmode="xml", verbose=False
    )
    assert result["data"] == payload.decode()
    assert parse_qs(urlsplit(seen[0]).query)["id"] == ["1,2"]
    assert parse_qs(urlsplit(seen[0]).query)["rettype"] == [rettype]


@pytest.mark.parametrize("linkname,pair", list(api._LINKS.items()))
def test_link_routes_are_consistent(monkeypatch, linkname, pair):
    seen = []

    def opener(request, **kwargs):
        seen.append(request.full_url)
        return Response(b'{"linksets": []}', request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = api.ncbi_link_dna_records(linkname, ["123"], verbose=False)
    query = parse_qs(urlsplit(seen[0]).query)
    assert (query["dbfrom"][0], query["db"][0]) == pair
    assert query["linkname"] == [linkname]
    assert result["data"] == {"linksets": []}


@pytest.mark.parametrize("database", sorted(api._SEARCH))
def test_search_preserves_query_warnings_and_zero_hits(monkeypatch, database):
    data = {"esearchresult": {"count": "0", "idlist": [], "warninglist": {"phrasesignored": ["x"]}}}
    monkeypatch.setattr(
        transport,
        "urlopen",
        lambda request, **kwargs: Response(json.dumps(data).encode(), request.full_url),
    )
    result = api.ncbi_search_dna_database(database, "x", retmax=0, verbose=False)
    assert result["data"] == data


@pytest.mark.parametrize("database", sorted(api._SUMMARY))
def test_summary(monkeypatch, database):
    monkeypatch.setattr(
        transport,
        "urlopen",
        lambda request, **kwargs: Response(b'{"result":{"uids":[]}}', request.full_url),
    )
    assert (
        api.ncbi_summarize_dna_records(database, ["1"], verbose=False)["data"]["result"]["uids"]
        == []
    )


@pytest.mark.parametrize("ids", [[], "123", [""]])
def test_bad_ids_rejected(ids):
    with pytest.raises(ConfigurationError):
        api.ncbi_link_dna_records("snp_gene", ids, verbose=False)


def test_xml_service_error_is_not_success(monkeypatch):
    monkeypatch.setattr(
        transport,
        "urlopen",
        lambda request, **kwargs: Response(
            b"<eFetchResult><ERROR>bad ID</ERROR></eFetchResult>", request.full_url
        ),
    )
    with pytest.raises(QueryError):
        api.ncbi_fetch_dna_records("nuccore", ["bad"], verbose=False)


def test_credentials_are_redacted(monkeypatch):
    monkeypatch.setattr(
        transport,
        "urlopen",
        lambda request, **kwargs: Response(b'{"linksets":[]}', request.full_url),
    )
    result = api.ncbi_link_dna_records(
        "snp_gene", ["1"], api_key="secret", email="private@example.org", verbose=False
    )
    assert "secret" not in json.dumps(result) and "private" not in json.dumps(result)
