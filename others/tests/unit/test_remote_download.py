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
    download_url,
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
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.search import QueryProgress, QueryResult
from dnakit.search._shared import adapter_provenance
from dnakit.tools import default_tool_registry

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


@pytest.mark.parametrize(
    "url,filename,expected_name",
    [
        ("https://example.test/data.fa.gz", None, "data.fa.gz"),
        ("https://example.test/my%20data.fa?version=2", None, "my data.fa"),
        ("https://example.test/%E5%BA%8F%E5%88%97.fa", None, "序列.fa"),
        ("https://example.test/?id=123", "chosen.fa", "chosen.fa"),
    ],
)
def test_download_url_saves_to_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    url: str,
    filename: str | None,
    expected_name: str,
) -> None:
    payload = b">seq\nACGT\n"
    calls: list[Any] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append((request, timeout))
        return _FakeResponse(payload, request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    output = tmp_path / "new" / "downloads"
    events: list[DownloadProgress] = []
    checksum = hashlib.sha256(payload).hexdigest()
    result = download_url(
        url,
        output,
        filename=filename,
        expected_sha256=checksum,
        config=DownloadConfig(chunk_size=2, timeout=12),
        headers={"Accept": "application/octet-stream"},
        progress=events.append,
    )

    target = output / expected_name
    assert target.read_bytes() == payload
    assert result.path == str(target.resolve())
    assert result.byte_size == len(payload)
    assert result.sha256 == checksum
    assert result.checksum_verified is True
    assert list(output.iterdir()) == [target]
    assert calls[0][0].full_url == url
    assert calls[0][0].headers["Accept"] == "application/octet-stream"
    assert calls[0][1] == 12
    assert events[-1].bytes_completed == len(payload)
    assert events[-1].target == str(target.resolve())


@pytest.mark.parametrize(
    "url,filename",
    [
        ("https://example.test/", None),
        ("https://example.test/directory/", None),
        ("https://example.test/directory/.", None),
        ("https://example.test/%2e%2e", None),
        ("https://example.test/..%2foutside.fa", None),
        ("https://example.test/%5coutside.fa", None),
        ("https://example.test/data.fa", "../outside.fa"),
        ("https://example.test/data.fa", "/outside.fa"),
        ("https://example.test/data.fa", "..\\outside.fa"),
        ("https://example.test/data.fa", "C:outside.fa"),
    ],
)
def test_download_url_rejects_missing_or_unsafe_names_before_downloading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    url: str,
    filename: str | None,
) -> None:
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda *_args, **_kwargs: pytest.fail("Invalid filename must not make a request"),
    )
    output = tmp_path / "downloads"
    with pytest.raises(ConfigurationError, match="INVALID_DOWNLOAD_FILENAME"):
        download_url(url, output, filename=filename, verbose=False)
    assert not output.exists()


def test_download_url_preserves_existing_file_until_verified_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "data.fa"
    target.write_bytes(b"previous complete file")
    calls: list[Any] = []
    payload = b"new file"

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        return _FakeResponse(payload, request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    url = "https://example.test/data.fa"
    with pytest.raises(FileExistsError):
        download_url(url, tmp_path, verbose=False)
    assert not calls

    config = DownloadConfig(overwrite=True, show_progress=False)
    with pytest.raises(DownloadError, match="CHECKSUM_MISMATCH"):
        download_url(url, tmp_path, config=config, expected_md5="0" * 32)
    assert target.read_bytes() == b"previous complete file"
    assert list(tmp_path.iterdir()) == [target]

    result = download_url(
        url,
        tmp_path,
        config=config,
        expected_md5=hashlib.md5(payload, usedforsecurity=False).hexdigest(),
    )
    assert target.read_bytes() == payload
    assert result.checksum_verified is True
    assert list(tmp_path.iterdir()) == [target]


def test_download_url_rejects_file_as_output_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "existing.txt"
    output.write_bytes(b"keep")
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda *_args, **_kwargs: pytest.fail("Invalid output directory must not make a request"),
    )
    with pytest.raises(NotADirectoryError):
        download_url("https://example.test/data.fa", output, verbose=False)
    assert output.read_bytes() == b"keep"


@pytest.mark.parametrize("show_progress,verbose", [(True, True), (False, True), (True, False)])
def test_download_url_progress_is_not_duplicated(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    show_progress: bool,
    verbose: bool,
) -> None:
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda request, *, timeout: _FakeResponse(b"ACGT", request.full_url),
    )
    download_url(
        "https://example.test/data.fa",
        tmp_path,
        config=DownloadConfig(show_progress=show_progress),
        verbose=verbose,
    )
    output = capsys.readouterr()
    assert output.out == ""
    if show_progress and verbose:
        assert output.err.count("Starting download:") == 1
        assert output.err.count("Download complete:") == 1
        assert "Downloading 1/1" in output.err
    else:
        assert output.err == ""


def test_download_url_is_discoverable_and_callable_by_agent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    registry = default_tool_registry()
    tool = registry.search("download_url")[0]
    assert tool.name == "dnakit.download.download_url"
    assert tool.agent_compatible is True
    assert tool.call_plan.input_schema["required"] == ["url", "output_dir"]
    monkeypatch.setattr(
        download_files,
        "urlopen",
        lambda request, *, timeout: _FakeResponse(b"ACGT", request.full_url),
    )
    result = registry.execute(
        tool.name,
        {
            "url": "https://example.test/data.fa",
            "output_dir": str(tmp_path),
            "config": {"show_progress": False},
        },
        allow_side_effects=True,
    )
    assert isinstance(result, dict)
    assert result["path"] == str((tmp_path / "data.fa").resolve())
    assert Path(result["path"]).read_bytes() == b"ACGT"
    assert result["byte_size"] == 4


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


def test_checksum_failure_leaves_no_target(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
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
    output = capsys.readouterr().err
    assert "Download failed or interrupted:" in output
    assert "Download complete:" not in output


@pytest.mark.parametrize("show_progress", [True, False])
@pytest.mark.parametrize("verbose", [True, False])
def test_dataset_writes_integrity_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    show_progress: bool,
    verbose: bool,
) -> None:
    payloads = {"a.txt": b"a", "b.txt": b"bb"}

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        name = Path(request.full_url).name
        return _FakeResponse(payloads[name], request.full_url)

    monkeypatch.setattr(download_files, "urlopen", fake_urlopen)
    config = DownloadConfig(show_progress=show_progress)
    result = dataset(
        (
            RemoteFile("https://example.test/a.txt"),
            RemoteFile("https://example.test/b.txt"),
        ),
        tmp_path / "dataset",
        kind="test_data",
        config=config,
        verbose=verbose,
    )

    assert config.show_progress is show_progress
    manifest = json.loads(Path(result.manifest_path).read_text(encoding="utf-8"))
    assert result.downloaded_bytes == 3
    assert len(manifest["files"]) == 2
    assert manifest["files"][1]["sha256"] == hashlib.sha256(b"bb").hexdigest()
    output = capsys.readouterr()
    assert output.out == ""
    if show_progress and verbose:
        assert output.err.count("Starting download:") == 1
        assert output.err.count("Download complete:") == 1
        assert "Downloading 2/2" in output.err
        assert ".dnakit-dataset-" not in output.err
    else:
        assert output.err == ""


@pytest.mark.parametrize("show_progress", [True, False])
def test_file_progress_without_content_length(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    show_progress: bool,
) -> None:
    response = _FakeResponse(b"ACGT", "https://example.test/seq.fa")
    response.headers.clear()
    monkeypatch.setattr(download_files, "urlopen", lambda *_args, **_kwargs: response)
    result = download_file(
        RemoteFile("https://example.test/seq.fa"),
        tmp_path / "seq.fa",
        config=DownloadConfig(show_progress=show_progress),
    )
    assert Path(result.path).read_bytes() == b"ACGT"
    output = capsys.readouterr()
    assert output.out == ""
    if show_progress:
        assert "Starting download:" in output.err
        assert "Download complete:" in output.err
        assert "100%" in output.err
    else:
        assert output.err == ""


def test_download_config_rejects_non_boolean_progress() -> None:
    with pytest.raises(ConfigurationError, match="show_progress"):
        DownloadConfig(show_progress="false")  # type: ignore[arg-type]


@pytest.mark.parametrize("show_progress", [True, False])
def test_sequence_query_progress(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    show_progress: bool,
) -> None:
    query = QueryResult(
        "sequence",
        "Ensembl",
        "https://example.test/sequence",
        ({"id": "region1", "seq": "ACGT"},),
        adapter_provenance("Ensembl REST", citation_url="https://rest.ensembl.org/"),
    )

    def search(*_args: Any, **kwargs: Any) -> QueryResult:
        if kwargs["progress"] is not None:
            kwargs["progress"](QueryProgress(1, 1, "1:0-4"))
        return query

    monkeypatch.setattr(sequence_module, "search_sequence", search)
    sequence(
        "human", "1:0-4", tmp_path / "seq.fa", config=DownloadConfig(show_progress=show_progress)
    )
    output = capsys.readouterr()
    assert output.out == ""
    if show_progress:
        assert "Fetching sequences" in output.err
        assert "Download complete:" in output.err
    else:
        assert output.err == ""


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
        if request.full_url.endswith("/b.txt"):
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
