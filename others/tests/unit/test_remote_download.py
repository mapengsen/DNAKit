"""Tests for checksummed public-data downloads and local index adapters."""

from __future__ import annotations

import hashlib
import importlib
import io
import json
from pathlib import Path
from typing import Any
from urllib.error import URLError

import pytest

import dnakit.download.files as download_files
from dnakit.download import (
    DownloadConfig,
    DownloadProgress,
    IndexProgress,
    RemoteFile,
    annotation,
    build_index,
    dataset,
    download_file,
    encode_files,
    expression,
    reads,
    sequence,
    variants,
    virus_package,
)
from dnakit.download import (
    metadata as export_metadata,
)
from dnakit.download import (
    ucsc_files as download_ucsc_files,
)
from dnakit.exceptions import DownloadError
from dnakit.search import QueryResult
from dnakit.search._shared import adapter_provenance

sequence_module = importlib.import_module("dnakit.download.sequence")


class _FakeResponse:
    def __init__(self, payload: bytes, url: str) -> None:
        self._stream = io.BytesIO(payload)
        self._url = url
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        self._stream.close()

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def geturl(self) -> str:
        return self._url


def test_download_file_is_atomic_and_checks_md5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"ACGT\n"
    calls: list[Any] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        return _FakeResponse(payload, request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    events: list[DownloadProgress] = []
    target = tmp_path / "data.fa"
    result = download_file(
        RemoteFile(
            "https://example.test/data.fa",
            expected_md5=hashlib.md5(payload, usedforsecurity=False).hexdigest(),
        ),
        target,
        progress=events.append,
    )

    assert target.read_bytes() == payload
    assert result.checksum_verified is True
    assert result.sha256 == hashlib.sha256(payload).hexdigest()
    assert events and events[-1].bytes_completed == len(payload)
    assert calls[0].headers["User-agent"].startswith("DNAKit/")


def test_checksum_failure_leaves_no_target(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payload = b"wrong"
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda request, *, timeout: _FakeResponse(payload, request.full_url),
    )
    target = tmp_path / "data.bin"
    with pytest.raises(DownloadError) as error:
        download_file(
            RemoteFile("https://example.test/data.bin", expected_md5="0" * 32),
            target,
        )
    assert error.value.code == "CHECKSUM_MISMATCH"
    assert not target.exists()
    assert list(tmp_path.glob("*.part")) == []


def test_dataset_writes_integrity_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payloads = {"a.txt": b"a", "b.txt": b"bb"}

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        name = Path(request.full_url).name
        return _FakeResponse(payloads[name], request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    result = dataset(
        (
            RemoteFile("https://example.test/a.txt"),
            RemoteFile("https://example.test/b.txt"),
        ),
        tmp_path / "dataset",
        kind="test_data",
    )

    manifest = json.loads(Path(result.manifest_path).read_text(encoding="utf-8"))
    assert result.downloaded_bytes == 3
    assert len(manifest["files"]) == 2
    assert manifest["files"][1]["sha256"] == hashlib.sha256(b"bb").hexdigest()


def test_ncbi_annotation_package_uses_selected_formats(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[Any] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        return _FakeResponse(b"zip", request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    result = annotation(
        "GCF_000001405.40",
        tmp_path / "annotation",
        formats=("GENOME_GFF", "SEQUENCE_REPORT"),
        config=DownloadConfig(api_key="secret"),
        api_base_url="https://example.test/v2",
    )

    assert Path(result.files[0].path).name.endswith(".zip")
    assert "include_annotation_type=GENOME_GFF" in calls[0].full_url
    assert "include_annotation_type=SEQUENCE_REPORT" in calls[0].full_url
    assert calls[0].headers["Api-key"] == "secret"


def test_ena_reads_downloads_paired_files_with_advertised_checksums(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = b"read-one"
    second = b"read-two"
    md5s = ";".join(
        hashlib.md5(value, usedforsecurity=False).hexdigest() for value in (first, second)
    )
    query = QueryResult(
        "reads",
        "ENA",
        "https://example.test/search",
        (
            {
                "fastq_ftp": "ftp.sra.ebi.ac.uk/run_1.fastq.gz;ftp.sra.ebi.ac.uk/run_2.fastq.gz",
                "fastq_md5": md5s,
                "fastq_bytes": f"{len(first)};{len(second)}",
            },
        ),
        adapter_provenance(
            "European Nucleotide Archive",
            citation_url="https://www.ebi.ac.uk/ena/portal/api/",
        ),
        total_count=1,
    )

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        payload = first if request.full_url.endswith("run_1.fastq.gz") else second
        return _FakeResponse(payload, request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    result = reads(query, tmp_path / "reads")

    assert len(result.files) == 2
    assert all(item.checksum_verified for item in result.files)
    assert all(item.url.startswith("https://ftp.sra.ebi.ac.uk") for item in result.files)


def test_coordinate_sequence_download_writes_fasta_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = QueryResult(
        "sequence",
        "Ensembl",
        "https://example.test/sequence",
        ({"id": "region1", "seq": "ACGT", "requested_region": "1:0-4"},),
        adapter_provenance("Ensembl REST", citation_url="https://rest.ensembl.org/"),
        total_count=1,
    )
    monkeypatch.setattr(sequence_module, "search_sequence", lambda *_args, **_kwargs: query)
    target = tmp_path / "region.fa"
    result = sequence("human", "1:0-4", target)

    assert target.read_text(encoding="ascii") == ">region1 requested=1:0-4\nACGT\n"
    assert Path(result.manifest_path).is_file()
    assert result.files[0].sha256 == hashlib.sha256(target.read_bytes()).hexdigest()


def test_explicit_makeblastdb_index_adapter(tmp_path: Path) -> None:
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">ref\nACGT\n", encoding="ascii")
    executable = tmp_path / "fake-makeblastdb"
    executable.write_text(
        "#!/bin/sh\n"
        "out=''\n"
        'while [ "$#" -gt 0 ]; do\n'
        "  if [ \"$1\" = '-out' ]; then out=$2; shift 2; else shift; fi\n"
        "done\n"
        "printf 'index' > \"${out}.nhr\"\n",
        encoding="ascii",
    )
    executable.chmod(0o700)
    events: list[IndexProgress] = []
    result = build_index(
        fasta,
        tool="makeblastdb",
        executable_path=executable,
        output_dir=tmp_path / "index",
        progress=events.append,
    )

    assert len(result.artifacts) == 1
    assert Path(result.artifacts[0].path).read_text(encoding="ascii") == "index"
    assert [event.stage for event in events] == ["start", "complete"]


def test_overwrite_dataset_rolls_back_when_later_download_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "dataset"
    output.mkdir()
    first = output / "a.txt"
    second = output / "b.txt"
    manifest = output / "test_data_manifest.json"
    first.write_text("old-a", encoding="ascii")
    second.write_text("old-b", encoding="ascii")
    manifest.write_text("old-manifest", encoding="ascii")
    call_count = 0

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise URLError("network failed")
        return _FakeResponse(b"new-a", request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    with pytest.raises(DownloadError):
        dataset(
            (
                RemoteFile("https://example.test/a.txt"),
                RemoteFile("https://example.test/b.txt"),
            ),
            output,
            kind="test_data",
            config=DownloadConfig(overwrite=True),
        )

    assert first.read_text(encoding="ascii") == "old-a"
    assert second.read_text(encoding="ascii") == "old-b"
    assert manifest.read_text(encoding="ascii") == "old-manifest"
    assert not tuple(output.glob(".dnakit-*"))


def test_sequence_preflights_manifest_before_remote_query(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "region.fa"
    manifest = tmp_path / "region.fa.manifest.json"
    manifest.write_text("existing", encoding="ascii")
    called = False

    def unexpected_query(*_args: object, **_kwargs: object) -> QueryResult:
        nonlocal called
        called = True
        raise AssertionError("remote query must not run")

    monkeypatch.setattr(sequence_module, "search_sequence", unexpected_query)
    with pytest.raises(FileExistsError):
        sequence("human", "1:0-4", target)
    assert called is False
    assert not target.exists()


def test_query_metadata_exports_csv_and_xml_with_manifests(tmp_path: Path) -> None:
    query = QueryResult(
        "sample",
        "test-provider",
        "https://example.test/query",
        ({"id": "S1", "nested": {"tissue": "liver"}},),
        adapter_provenance("Test Provider", citation_url="https://example.test"),
        total_count=1,
    )

    csv_result = export_metadata(query, tmp_path / "samples.csv")
    xml_result = export_metadata(query, tmp_path / "samples.xml")

    assert (tmp_path / "samples.csv").read_text(encoding="utf-8").splitlines()[0] == "id,nested"
    assert '""tissue"": ""liver""' in (tmp_path / "samples.csv").read_text(encoding="utf-8")
    assert b'<field name="id">S1</field>' in (tmp_path / "samples.xml").read_bytes()
    assert Path(csv_result.manifest_path).is_file()
    assert Path(xml_result.manifest_path).is_file()


def test_ucsc_catalog_download_uses_bounded_https_resources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = QueryResult(
        "files",
        "UCSC",
        "https://api.genome.ucsc.edu/list/files?genome=hg38",
        ({"url": "goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", "sizeBytes": 5},),
        adapter_provenance(
            "UCSC Genome Browser REST API",
            citation_url="https://genome.ucsc.edu/goldenPath/help/api.html",
        ),
        total_count=1,
    )
    calls: list[Any] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        return _FakeResponse(b"chain", request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    result = download_ucsc_files(query, tmp_path / "ucsc")

    assert Path(result.files[0].path).read_bytes() == b"chain"
    assert Path(result.files[0].path).name.endswith("hg38ToHg19.over.chain.gz")
    assert calls[0].full_url.startswith("https://hgdownload.soe.ucsc.edu/goldenPath/")


def test_virus_dbsnp_and_geo_matrix_download_urls(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[Any] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        return _FakeResponse(b"data", request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    virus_result = virus_package(
        "NC_045512.2",
        tmp_path / "virus",
        api_base_url="https://example.test/v2",
    )
    dbsnp_result = variants(
        tmp_path / "dbsnp",
        source="dbsnp",
        include_index=False,
    )
    matrix_result = expression(
        "GSE100",
        tmp_path / "geo",
        format="matrix",
    )

    assert "/virus/accession/NC_045512.2/genome/download" in calls[0].full_url
    assert any("GCF_000001405.40.gz" in call.full_url for call in calls)
    assert calls[-1].full_url.endswith("/GSE100/matrix/GSE100_series_matrix.txt.gz")
    assert Path(virus_result.manifest_path).is_file()
    assert Path(dbsnp_result.manifest_path).is_file()
    assert Path(matrix_result.manifest_path).is_file()


def test_encode_file_download_verifies_provider_md5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"bigwig"
    query = QueryResult(
        "regulation",
        "ENCODE",
        "https://www.encodeproject.org/search/?type=File",
        (
            {
                "accession": "ENCFF000AAA",
                "href": "/files/ENCFF000AAA/@@download/ENCFF000AAA.bigWig",
                "file_size": len(payload),
                "md5sum": hashlib.md5(payload, usedforsecurity=False).hexdigest(),
                "file_format": "bigWig",
            },
        ),
        adapter_provenance(
            "ENCODE Portal", citation_url="https://www.encodeproject.org/help/rest-api/"
        ),
        total_count=1,
    )
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda request, *, timeout: _FakeResponse(payload, request.full_url),
    )

    result = encode_files(query, tmp_path / "encode")

    assert result.files[0].checksum_verified is True
    assert Path(result.files[0].path).name == "ENCFF000AAA.bigWig"
