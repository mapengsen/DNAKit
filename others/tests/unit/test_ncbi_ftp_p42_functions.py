"""Offline tests of actual streamed bytes, catalog routing and integrity failures."""

from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path

import pytest

from dnakit.download import files
from dnakit.download import ncbi_ftp_p42_functions as api
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import _http


class Response(io.BytesIO):
    def __init__(self, data: bytes, url: str) -> None:
        super().__init__(data)
        self.headers = {
            "Content-Length": str(len(data)),
            "Last-Modified": "Tue, 20 Oct 2026 00:00:00 GMT",
        }
        self.url = url

    def geturl(self) -> str:
        return self.url


@pytest.mark.parametrize("collection", list(api.NCBI_FTP_COLLECTIONS))
def test_each_curated_dataset_downloads_real_bytes(monkeypatch, tmp_path, collection):
    body = b">record\nACGT\n"
    urls = []

    def opener(request, **kwargs):
        urls.append(request.full_url)
        return Response(body, request.full_url)

    monkeypatch.setattr(files, "urlopen", opener)
    result = api.download_ncbi_ftp_collection(
        collection,
        ["release/input.fa"],
        str(tmp_path),
        expected_md5={"release/input.fa": hashlib.md5(body).hexdigest()},
        verbose=False,
    )
    assert (tmp_path / "input.fa").read_bytes() == body
    assert result["files"][0]["checksum_verified"] is True
    assert Path(result["manifest_path"]).is_file()
    assert urls == [
        "https://ftp.ncbi.nlm.nih.gov/"
        + api.NCBI_FTP_COLLECTIONS[collection]["root"]
        + "release/input.fa"
    ]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize(
    "path",
    [
        "../bad",
        "%2e%2e/bad",
        "/bad",
        "a/../../x",
        "a\\x",
        "a?x",
        "a#x",
        "https://example.org/x",
        "directory/",
    ],
)
def test_path_escape_rejected_without_network(path, tmp_path):
    with pytest.raises(ConfigurationError):
        api.download_ncbi_ftp_collection(
            next(iter(api.NCBI_FTP_COLLECTIONS)), [path], str(tmp_path), verbose=False
        )


def test_listing_keeps_only_collection_children(monkeypatch):
    body = b'<a href="../">parent</a><a href="release/">release</a><a href="x.fa.gz">file</a><a href="https://evil.invalid/x">offsite</a>'
    monkeypatch.setattr(
        _http, "urlopen", lambda request, **kwargs: Response(body, request.full_url)
    )
    result = api.list_ncbi_ftp_collection("cds_nucleotide_sequences_fasta", verbose=False)
    assert [f["path"] for f in result["files"]] == ["release/", "x.fa.gz"]
    assert result["files"][0]["is_directory"] is True


def test_failed_checksum_installs_no_file(monkeypatch, tmp_path):
    monkeypatch.setattr(
        files, "urlopen", lambda request, **kwargs: Response(b"ACGT", request.full_url)
    )
    with pytest.raises(DownloadError):
        api.download_ncbi_ftp_collection(
            "cds_nucleotide_sequences_fasta",
            ["x.fna"],
            str(tmp_path),
            expected_md5={"x.fna": "0" * 32},
            verbose=False,
        )
    assert not (tmp_path / "x.fna").exists()


def test_duplicate_basenames_and_unknown_catalog_rejected(tmp_path):
    with pytest.raises(ConfigurationError):
        api.download_ncbi_ftp_collection(
            "cds_nucleotide_sequences_fasta", ["a/x", "b/x"], str(tmp_path), verbose=False
        )
    with pytest.raises(ConfigurationError):
        api.list_ncbi_ftp_collection("unknown", verbose=False)
