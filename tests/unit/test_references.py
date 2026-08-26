"""Tests for the NCBI reference-genome downloader."""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path
from typing import Any

import pytest

import dnakit.references.ncbi as ncbi
from dnakit.exceptions import ConfigurationError, DownloadError
from dnakit.references import download_genome, resolve_genome_assembly


class _FakeResponse:
    def __init__(self, payload: bytes, *, content_length: int | None = None) -> None:
        self._stream = io.BytesIO(payload)
        self.headers = {}
        if content_length is not None:
            self.headers["Content-Length"] = str(content_length)

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        self._stream.close()

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)


def _package_payload(
    *, accession: str = "GCF_000001405.40", checksum_override: str | None = None
) -> bytes:
    fasta_member = f"ncbi_dataset/data/{accession}/{accession}_Test_genomic.fna"
    fasta = b">NC_000001.1 Test chromosome\nACGTACGT\n"
    md5 = hashlib.md5(fasta, usedforsecurity=False).hexdigest()
    checksum = checksum_override or md5
    report = {
        "accession": accession,
        "assembly": {"accession": accession, "name": "TestAssembly"},
        "organism": {"organism_name": "Test species"},
    }
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.md", "test package\n")
        archive.writestr("md5sum.txt", f"{checksum}  {fasta_member}\n")
        archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl", json.dumps(report) + "\n")
        archive.writestr(fasta_member, fasta)
    return stream.getvalue()


def test_fixed_alias_resolves_without_network() -> None:
    result = resolve_genome_assembly("hg38", api_base_url="https://example.test/api")

    assert result.accession == "GCF_000001405.40"
    assert result.query == "hg38"


def test_download_extracts_fasta_and_verifies_md5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = _package_payload()
    events: list[object] = []

    def fake_urlopen(_request: Any, *, timeout: float) -> _FakeResponse:
        assert timeout == 12.0
        return _FakeResponse(payload, content_length=len(payload))

    monkeypatch.setattr(ncbi, "urlopen", fake_urlopen)
    result = download_genome(
        "hg38",
        tmp_path / "reference",
        api_base_url="https://example.test/api",
        timeout=12.0,
        progress=events.append,
    )

    fasta_path = Path(result.fasta_path)
    assert fasta_path.name == "GCF_000001405.40_Test_genomic.fna"
    assert fasta_path.read_text(encoding="utf-8").startswith(">NC_000001.1")
    assert result.organism == "Test species"
    assert result.assembly_name == "TestAssembly"
    assert (
        result.fasta_md5 == hashlib.md5(fasta_path.read_bytes(), usedforsecurity=False).hexdigest()
    )
    assert Path(result.metadata_path).is_file()
    assert Path(result.checksum_path).is_file()
    assert result.package_path is None
    assert any(getattr(event, "phase", None) == "download" for event in events)
    assert any(getattr(event, "phase", None) == "extract" for event in events)


def test_taxon_resolution_requires_one_current_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = _package_payload()
    calls: list[str] = []

    def fake_urlopen(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request.full_url)
        if "dataset_report" in request.full_url:
            return _FakeResponse(
                json.dumps(
                    {
                        "reports": [
                            {
                                "accession": "GCF_000001405.40",
                                "source_database": "SOURCE_DATABASE_REFSEQ",
                            }
                        ]
                    }
                ).encode("utf-8")
            )
        return _FakeResponse(payload, content_length=len(payload))

    monkeypatch.setattr(ncbi, "urlopen", fake_urlopen)
    result = download_genome(
        "human",
        tmp_path / "reference",
        api_base_url="https://example.test/api",
        progress=None,
    )

    assert result.accession == "GCF_000001405.40"
    assert len(calls) == 2
    assert "filters.reference_only=true" in calls[0]
    assert "tax_exact_match=true" in calls[0]


def test_taxon_resolution_rejects_ambiguous_results(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_urlopen(_request: Any, *, timeout: float) -> _FakeResponse:
        return _FakeResponse(
            json.dumps(
                {
                    "reports": [
                        {"accession": "GCF_000000001.1"},
                        {"accession": "GCF_000000002.1"},
                    ]
                }
            ).encode("utf-8")
        )

    monkeypatch.setattr(ncbi, "urlopen", fake_urlopen)

    with pytest.raises(ConfigurationError) as error:
        resolve_genome_assembly("ambiguous species", api_base_url="https://example.test/api")
    assert error.value.code == "AMBIGUOUS_GENOME_ASSEMBLY"


def test_download_rejects_bad_checksum(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payload = _package_payload(checksum_override="0" * 32)
    monkeypatch.setattr(
        ncbi,
        "urlopen",
        lambda _request, *, timeout: _FakeResponse(payload, content_length=len(payload)),
    )

    with pytest.raises(DownloadError) as error:
        download_genome(
            "hg38",
            tmp_path / "reference",
            api_base_url="https://example.test/api",
        )
    assert error.value.code == "CHECKSUM_MISMATCH"
    reference_dir = tmp_path / "reference"
    if reference_dir.exists():
        assert list(reference_dir.glob("*")) == []


def test_download_enforces_package_size_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = _package_payload()
    monkeypatch.setattr(
        ncbi,
        "urlopen",
        lambda _request, *, timeout: _FakeResponse(payload, content_length=len(payload)),
    )

    with pytest.raises(DownloadError) as error:
        download_genome(
            "hg38",
            tmp_path / "reference",
            api_base_url="https://example.test/api",
            max_download_bytes=len(payload) - 1,
        )
    assert error.value.code == "DOWNLOAD_SIZE_LIMIT"
