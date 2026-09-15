"""F1: offline fault injection through public DNAKit download entry points."""

from __future__ import annotations

import hashlib
import io
import json
import time
import zipfile
from email.message import Message
from http.client import IncompleteRead
from pathlib import Path
from typing import Any
from urllib.error import HTTPError

import pytest

import dnakit.download.files as remote
import dnakit.references.ncbi as ncbi
from dnakit.download import DownloadConfig, RemoteFile, download_file
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.references import download_genome, resolve_genome_assembly

ACCESSION = "GCF_000001405.40"


class Response(io.BytesIO):
    def __init__(self, data: bytes, length: int | None = None) -> None:
        super().__init__(data)
        self.headers = {} if length is None else {"Content-Length": str(length)}


def package(
    *, accession: Any = ACCESSION, report_accession: Any = ACCESSION, checksum: Any = None
) -> Any:
    payload = io.BytesIO()
    sequence = b">chr1\nAAGC\n"
    member = f"ncbi_dataset/data/{accession}/{accession}_Test_genomic.fna"
    with zipfile.ZipFile(payload, "w") as archive:
        archive.writestr(member, sequence)
        archive.writestr(
            "md5sum.txt", f"{checksum or hashlib.md5(sequence).hexdigest()}  {member}\n"
        )
        archive.writestr(
            "ncbi_dataset/data/assembly_data_report.jsonl",
            json.dumps({"accession": report_accession}),
        )
    return payload.getvalue()


@pytest.fixture(autouse=True)
def no_retry_delay(monkeypatch: Any) -> Any:
    monkeypatch.setattr(time, "sleep", lambda seconds: None)


@pytest.mark.parametrize("adapter", ["genome", "file"])
@pytest.mark.parametrize(
    "fault,attempts", [(404, 1), (429, 3), (500, 3), ("timeout", 3), ("disconnect", 3)]
)
def test_network_failure_is_bounded_and_preserves_output(
    tmp_path: Any, monkeypatch: Any, adapter: Any, fault: Any, attempts: Any
) -> None:
    target = tmp_path / "existing.txt"
    target.write_text("previous complete result")
    calls: list[int] = []

    def fail(*args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        if isinstance(fault, int):
            raise HTTPError("https://example.test/data", fault, "test fault", Message(), None)
        if fault == "timeout":
            raise TimeoutError("test timeout")
        raise IncompleteRead(b"partial", 100)

    module = ncbi if adapter == "genome" else remote
    monkeypatch.setattr(module, "urlopen", fail)
    with pytest.raises(DownloadError) as error:
        if adapter == "genome":
            download_genome("hg38", tmp_path, overwrite=True, verbose=False)
        else:
            download_file(
                RemoteFile("https://example.test/data"),
                target,
                config=DownloadConfig(overwrite=True),
                verbose=False,
            )
    assert len(calls) == attempts
    if attempts > 1:
        assert error.value.context["attempts"] == attempts
    assert target.read_text() == "previous complete result"
    assert list(tmp_path.iterdir()) == [target]


@pytest.mark.parametrize("adapter", ["genome", "file"])
def test_retry_restarts_and_recovers(tmp_path: Any, monkeypatch: Any, adapter: Any) -> None:
    data = package() if adapter == "genome" else b"ACGT"
    calls: list[int] = []

    def recover(*args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        if len(calls) == 1:
            raise TimeoutError("temporary")
        return Response(data, len(data))

    monkeypatch.setattr(ncbi if adapter == "genome" else remote, "urlopen", recover)
    if adapter == "genome":
        result = download_genome("hg38", tmp_path, verbose=False)
        assert Path(result.fasta_path).read_bytes() == b">chr1\nAAGC\n"
    else:
        download_file(RemoteFile("https://example.test/data"), tmp_path / "dna.txt", verbose=False)
        assert (tmp_path / "dna.txt").read_bytes() == data
    assert len(calls) == 2


@pytest.mark.parametrize("adapter", ["genome", "file"])
def test_truncated_response_never_becomes_success(
    tmp_path: Any, monkeypatch: Any, adapter: Any
) -> None:
    data = package() if adapter == "genome" else b"ACGT"
    monkeypatch.setattr(
        ncbi if adapter == "genome" else remote,
        "urlopen",
        lambda *args, **kwargs: Response(data, len(data) + 1),
    )
    with pytest.raises(DownloadError, match="DOWNLOAD_LENGTH_MISMATCH"):
        if adapter == "genome":
            download_genome("hg38", tmp_path, verbose=False)
        else:
            download_file(
                RemoteFile("https://example.test/data"), tmp_path / "dna.txt", verbose=False
            )
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize(
    "data,code",
    [
        (b"truncated zip", "INVALID_ARCHIVE"),
        (package(checksum="0" * 32), "CHECKSUM_MISMATCH"),
        (package(report_accession="GCF_000001405.39"), "ASSEMBLY_MISMATCH"),
        (package(accession="GCF_000001405.39"), "INVALID_ARCHIVE"),
    ],
)
def test_invalid_genome_package_preserves_all_previous_files(
    tmp_path: Any, monkeypatch: Any, data: Any, code: Any
) -> None:
    old = tmp_path / f"{ACCESSION}_Test_genomic.fna"
    old.write_text("previous genome")
    monkeypatch.setattr(ncbi, "urlopen", lambda *args, **kwargs: Response(data))
    with pytest.raises(DownloadError, match=code):
        download_genome("hg38", tmp_path, overwrite=True, verbose=False)
    assert old.read_text() == "previous genome"
    assert list(tmp_path.iterdir()) == [old]


def test_invalid_genome_query_does_not_fall_back(monkeypatch: Any) -> None:
    monkeypatch.setattr(ncbi, "urlopen", lambda *args, **kwargs: Response(b'{"reports": []}'))
    with pytest.raises(ConfigurationError, match="GENOME_ASSEMBLY_NOT_FOUND"):
        resolve_genome_assembly("hg999999", verbose=False)


@pytest.mark.parametrize("adapter", ["genome", "file"])
def test_existing_output_is_not_overwritten(tmp_path: Any, monkeypatch: Any, adapter: Any) -> None:
    target = tmp_path / (f"{ACCESSION}_Test_genomic.fna" if adapter == "genome" else "dna.txt")
    target.write_text("old")
    monkeypatch.setattr(ncbi, "urlopen", lambda *args, **kwargs: Response(package()))
    with pytest.raises(FileExistsError):
        if adapter == "genome":
            download_genome("hg38", tmp_path, verbose=False)
        else:
            download_file(RemoteFile("https://example.test/data"), target, verbose=False)
    assert target.read_text() == "old"


@pytest.mark.parametrize("adapter", ["genome", "file"])
def test_unwritable_directory_has_clear_failure(
    tmp_path: Any, monkeypatch: Any, adapter: Any
) -> None:
    def denied(*args: Any, **kwargs: Any) -> Any:
        raise PermissionError("Read-only output directory")

    monkeypatch.setattr(Path, "mkdir", denied)
    with pytest.raises(PermissionError, match="Read-only"):
        if adapter == "genome":
            download_genome("hg38", tmp_path, verbose=False)
        else:
            download_file(
                RemoteFile("https://example.test/data"), tmp_path / "dna.txt", verbose=False
            )
    assert not list(tmp_path.iterdir())


def test_midstream_disconnect_discards_partial_bytes_before_retry(
    tmp_path: Any, monkeypatch: Any
) -> None:
    class Interrupted(Response):
        def read(self, size: Any = -1) -> Any:
            if self.tell():
                raise IncompleteRead(b"partial", 10)
            return super().read(2)

    responses = iter([Interrupted(b"ACGT"), Response(b"TGCA", 4)])
    monkeypatch.setattr(remote, "urlopen", lambda *args, **kwargs: next(responses))
    target = tmp_path / "dna.txt"
    download_file(RemoteFile("https://example.test/dna"), target, verbose=False)
    assert target.read_bytes() == b"TGCA"
    assert list(tmp_path.iterdir()) == [target]


def test_shared_metadata_query_retries_transient_get(monkeypatch: Any) -> None:
    import dnakit.search._http as http
    from dnakit.search import SearchConfig

    calls: list[int] = []

    def response(*args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        if len(calls) < 3:
            raise TimeoutError("transient")
        return Response(b'{"records": []}')

    monkeypatch.setattr(http, "urlopen", response)
    assert http.request_json("https://example.test/query", SearchConfig(), provider="test") == {
        "records": []
    }
    assert len(calls) == 3


def test_post_is_not_retried_automatically(monkeypatch: Any) -> None:
    import dnakit.search._http as http
    from dnakit.exceptions import QueryError
    from dnakit.search import SearchConfig

    calls: list[int] = []

    def failed(*args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        raise TimeoutError("ambiguous server state")

    monkeypatch.setattr(http, "urlopen", failed)
    with pytest.raises(QueryError):
        http.request_text(
            "https://example.test/jobs", SearchConfig(), provider="test", method="POST", data=b"job"
        )
    assert len(calls) == 1
