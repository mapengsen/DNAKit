"""Offline ENA metadata routes and streamed record bytes, including empty/error edges."""

from __future__ import annotations

import io
import json

import pytest

from dnakit.databases import _network_p44 as transport
from dnakit.download import ena_p44_functions as api
from dnakit.download import files
from dnakit.exceptions import ConfigurationError


class Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url
        self.headers = {"Content-Length": str(len(data)), "Date": "Tue, 20 Oct 2026 00:00:00 GMT"}

    def geturl(self):
        return self.url


@pytest.mark.parametrize("format", ["embl", "fasta", "xml"])
def test_record_download_really_writes_bytes(monkeypatch, tmp_path, format):
    seen = []
    body = b">ENA_record\nACGT\n"

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(body, request.full_url)

    monkeypatch.setattr(files, "urlopen", opener)
    result = api.ena_download_record(
        "ABCD01", str(tmp_path / "records"), format=format, verbose=False
    )
    assert (tmp_path / "records").read_bytes() == body
    assert f"/{format}/ABCD01" in seen[0]
    assert result["file"]["byte_size"] == len(body)


@pytest.mark.parametrize(
    "function,kwargs,fragment",
    [
        (api.ena_download_record, {"accession": "ABCD01", "full_contig_set": True}, "/embl/ABCD01"),
        (
            api.ena_download_search,
            {"query": "tax_tree(9606)", "result": "sequence", "format": "fasta"},
            "/fasta/search",
        ),
        (
            api.ena_download_search,
            {"query": "x", "result": "read_run", "format": "xml"},
            "/xml/search",
        ),
        (
            api.ena_download_search,
            {"query": "x", "result": "wgs_set", "format": "embl"},
            "/embl/search",
        ),
        (
            api.ena_download_linked_records,
            {"link_type": "study", "accession": "PRJEB1", "result": "sequence"},
            "/embl/links/study",
        ),
        (
            api.ena_download_keyword_results,
            {"query": "BRCA1", "domain": "embl", "result": "sequence"},
            "/embl/textsearch",
        ),
        (api.ena_download_cram_reference, {"checksum": "a" * 32}, "/cram/md5/"),
        (
            api.ena_download_cram_reference,
            {"checksum": "a" * 40, "algorithm": "sha1"},
            "/cram/sha1/",
        ),
    ],
)
def test_all_download_routes(monkeypatch, tmp_path, function, kwargs, fragment):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b"ACGT", request.full_url)

    monkeypatch.setattr(files, "urlopen", opener)
    result = function(output_path=str(tmp_path / "data"), verbose=False, **kwargs)
    assert fragment in seen[0]
    assert (tmp_path / "data").read_bytes() == b"ACGT"
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize(
    "function,kwargs,fragment",
    [
        (
            api.ena_count_keyword_results,
            {"query": "x", "domain": "embl", "result": "sequence"},
            "/embl/textsearch/count",
        ),
        (
            api.ena_linked_accessions,
            {"link_type": "sample", "accession": "SAMEA1", "result": "sequence"},
            "/links/sample",
        ),
        (api.ena_search_cross_references, {"source": "ENA", "target": "Ensembl"}, "/tsv/search"),
    ],
)
def test_text_metadata_preserves_service_output(monkeypatch, function, kwargs, fragment):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b"count\n0\n", request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = function(**kwargs, verbose=False)
    assert result["data"] == "count\n0\n" and fragment in seen[0]


@pytest.mark.parametrize(
    "kwargs",
    [{"checksum": "x" * 32}, {"checksum": "a" * 31}, {"checksum": "a" * 32, "algorithm": "sha256"}],
)
def test_bad_checksum_input_never_requests(tmp_path, kwargs):
    with pytest.raises(ConfigurationError):
        api.ena_download_cram_reference(output_path=str(tmp_path / "data"), verbose=False, **kwargs)


def test_existing_target_protected(monkeypatch, tmp_path):
    target = tmp_path / "data"
    target.write_text("original")
    with pytest.raises(FileExistsError):
        api.ena_download_record("ABCD01", str(target), verbose=False)
    assert target.read_text() == "original"
