"""Tests for bounded public biological-database query adapters."""

from __future__ import annotations

import io
import json
from collections.abc import Callable, Mapping
from typing import Any
from urllib.parse import parse_qs, urlsplit

import pytest

import dnakit
import dnakit.search._http as remote_http
import dnakit.search.blast as blast_module
from dnakit.exceptions import ConfigurationError, QueryError
from dnakit.search import (
    BlastJob,
    QueryProgress,
    SearchConfig,
    annotation,
    assembly,
    blast_results,
    encode_search,
    ensembl_lookup,
    entrez,
    gene,
    identify,
    reads,
    sequence,
    submit_blast,
    taxonomy,
    transcripts,
    ucsc_chromosomes,
    ucsc_files,
    ucsc_sequence,
    ucsc_track_data,
    ucsc_tracks,
    variant,
    virus,
)


class _FakeResponse:
    def __init__(self, payload: bytes, *, url: str = "https://example.test/result") -> None:
        self._stream = io.BytesIO(payload)
        self._url = url
        self.headers: dict[str, str] = {"Content-Length": str(len(payload))}

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, *_args: object) -> None:
        self._stream.close()

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def geturl(self) -> str:
        return self._url


def _json_opener(payloads: list[object], calls: list[Any]) -> Callable[..., _FakeResponse]:
    encoded = [json.dumps(payload).encode("utf-8") for payload in payloads]

    def open_fake(request: Any, *, timeout: float) -> _FakeResponse:
        calls.append(request)
        assert timeout > 0
        return _FakeResponse(encoded.pop(0), url=request.full_url)

    return open_fake


def test_public_namespaces_are_loaded() -> None:
    assert dnakit.search.taxonomy is taxonomy
    assert callable(dnakit.download.genome)


def test_search_config_rejects_unbounded_values() -> None:
    with pytest.raises(ConfigurationError):
        SearchConfig(max_records=10_001)
    with pytest.raises(ConfigurationError):
        SearchConfig(max_response_bytes=0)


def test_taxonomy_and_assembly_build_current_ncbi_dataset_queries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                {"reports": [{"taxonomy": {"tax_id": 9606}}], "total_count": 1},
                {"reports": [{"accession": "GCF_000001405.40"}], "total_count": 1},
            ],
            calls,
        ),
    )
    config = SearchConfig(api_key="secret")
    taxon = taxonomy(
        "Homo sapiens",
        config=config,
        api_base_url="https://example.test/v2",
    )
    genomes = assembly(
        "human",
        reference_only=True,
        config=config,
        api_base_url="https://example.test/v2",
    )

    assert taxon.records[0]["taxonomy"] == {"tax_id": 9606}
    assert taxon.total_count == 1
    assert calls[0].headers["Api-key"] == "secret"
    assert "Homo%20sapiens" in calls[0].full_url
    assembly_query = parse_qs(urlsplit(calls[1].full_url).query)
    assert assembly_query["filters.reference_only"] == ["true"]
    assert assembly_query["filters.assembly_version"] == ["current"]
    assert "secret" not in genomes.request_url


def test_gene_symbol_product_report_uses_taxon(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener([{"reports": [{"gene": {"gene_id": "675"}}]}], calls),
    )

    result = gene(
        "BRCA2",
        taxon="human",
        report="product",
        api_base_url="https://example.test/v2",
    )

    assert result.records[0]["gene"] == {"gene_id": "675"}
    assert "/gene/symbol/BRCA2/taxon/human/product_report" in calls[0].full_url


def test_entrez_runs_search_then_summary_and_redacts_credentials(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                {"esearchresult": {"count": "2", "idlist": ["11", "22"]}},
                {
                    "result": {
                        "uids": ["11", "22"],
                        "11": {"uid": "11", "title": "one"},
                        "22": {"uid": "22", "title": "two"},
                    }
                },
            ],
            calls,
        ),
    )
    result = entrez(
        "pubmed",
        "BRCA2",
        config=SearchConfig(api_key="key-value", email="person@example.org"),
        api_base_url="https://example.test/eutils",
    )

    assert [record["uid"] for record in result.records] == ["11", "22"]
    assert result.total_count == 2
    assert "key-value" not in result.request_url
    assert "person%40example.org" not in result.request_url
    assert "REDACTED" in result.request_url
    assert calls[0].full_url.endswith("api_key=key-value")


def test_coordinate_sequence_converts_zero_based_regions_and_reports_progress(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                {"id": "first", "seq": "ACGT"},
                {"id": "second", "seq": "TGCA"},
            ],
            calls,
        ),
    )
    progress: list[QueryProgress] = []

    result = sequence(
        "human",
        ("1:0-4", "X:10-14"),
        strand=-1,
        upstream=2,
        progress=progress.append,
        api_base_url="https://example.test",
    )

    assert len(result.records) == 2
    assert result.records[0]["requested_start_0based"] == 0
    assert "/1:1..4:-1?" in calls[0].full_url
    assert parse_qs(urlsplit(calls[0].full_url).query)["expand_5prime"] == ["2"]
    assert [item.completed for item in progress] == [1, 2]


def test_annotation_and_variant_add_zero_based_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                [{"id": "ENSG1", "start": 11, "end": 20}],
                {
                    "name": "rs1",
                    "mappings": [{"start": 101, "end": 101, "assembly_name": "GRCh38"}],
                },
            ],
            calls,
        ),
    )

    features = annotation(
        "human", "1:10-20", features=("gene",), api_base_url="https://example.test"
    )
    snp = variant("rs1", api_base_url="https://example.test")

    assert features.records[0]["start_0based"] == 10
    mappings = snp.records[0]["mappings"]
    assert isinstance(mappings, tuple)
    first_mapping = mappings[0]
    assert isinstance(first_mapping, Mapping)
    assert first_mapping["start_0based"] == 100
    assert parse_qs(urlsplit(calls[1].full_url).query)["pops"] == ["true"]


def test_ena_and_encode_queries_return_provider_records(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                [{"accession": "SRR1", "fastq_ftp": "ftp.sra.ebi.ac.uk/a.fastq.gz"}],
                {"total": 1, "@graph": [{"accession": "ENCSR1"}]},
            ],
            calls,
        ),
    )

    runs = reads(
        'run_accession="SRR1"',
        fields=("accession", "fastq_ftp"),
        api_base_url="https://example.test/ena",
    )
    encode = encode_search(
        search_term="CTCF",
        filters={"assay_title": "ChIP-seq"},
        api_base_url="https://example.test",
    )

    assert runs.records[0]["accession"] == "SRR1"
    assert encode.records[0]["accession"] == "ENCSR1"
    assert parse_qs(urlsplit(calls[1].full_url).query)["assay_title"] == ["ChIP-seq"]


def test_blast_submit_status_and_normalized_results(monkeypatch: pytest.MonkeyPatch) -> None:
    config = SearchConfig(email="person@example.org")
    responses = iter(
        [
            "RID = ABCDEFGH1234\nRTOE = 30\n",
            "Status=READY\nThereAreHits=yes\n",
            "RID = ABCDEFGH1234\nRTOE = 30\n",
        ]
    )
    monkeypatch.setattr(blast_module, "request_text", lambda *_args, **_kwargs: next(responses))
    job = submit_blast(
        "ACGTACGT",
        config=config,
        api_url="https://example.test/Blast.cgi",
    )
    status = blast_module.blast_status(
        job,
        config=config,
        api_url="https://example.test/Blast.cgi",
    )

    assert job.rid == "ABCDEFGH1234"
    assert status.status == "ready"
    assert status.has_hits is True
    identified = identify(
        "ACGT",
        config=config,
        api_url="https://example.test/Blast.cgi",
    )
    assert isinstance(identified, BlastJob)
    assert identified.rid == "ABCDEFGH1234"

    payload = {
        "BlastOutput2": [
            {
                "report": {
                    "results": {
                        "search": {
                            "query_len": 100,
                            "hits": [
                                {
                                    "description": [
                                        {
                                            "accession": "NC_1",
                                            "title": "test",
                                            "taxid": 9606,
                                            "sciname": "Homo sapiens",
                                        }
                                    ],
                                    "hsps": [
                                        {
                                            "identity": 72,
                                            "align_len": 80,
                                            "bit_score": 100.0,
                                            "evalue": 1e-20,
                                        }
                                    ],
                                }
                            ],
                        }
                    }
                }
            }
        ]
    }
    monkeypatch.setattr(blast_module, "request_json", lambda *_args, **_kwargs: payload)
    result = blast_results(job, config=config, api_url="https://example.test/Blast.cgi")

    assert result.records[0]["identity"] == pytest.approx(0.9)
    assert result.records[0]["query_coverage"] == pytest.approx(0.8)
    assert result.metadata["novelty_score"] == pytest.approx(0.1)


def test_response_limits_are_enforced(monkeypatch: pytest.MonkeyPatch) -> None:
    payload = b"{" + b" " * 100 + b"}"
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        lambda request, *, timeout: _FakeResponse(payload, url=request.full_url),
    )
    with pytest.raises(QueryError) as error:
        taxonomy(
            "human",
            config=SearchConfig(max_response_bytes=10),
            api_base_url="https://example.test",
        )
    assert error.value.code == "QUERY_RESPONSE_SIZE_LIMIT"


def test_ensembl_lookup_and_transcripts_expand_nested_coordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    payload = {
        "id": "ENSG1",
        "start": 11,
        "end": 30,
        "canonical_transcript": "ENST1.1",
        "Transcript": [
            {
                "id": "ENST1",
                "start": 12,
                "end": 29,
                "Exon": [{"id": "ENSE1", "start": 12, "end": 15}],
            }
        ],
    }
    monkeypatch.setattr(remote_http, "urlopen", _json_opener([payload, payload], calls))

    by_id = ensembl_lookup("ENSG1", api_base_url="https://example.test")
    by_symbol = transcripts("BRCA2", api_base_url="https://example.test")

    assert by_id.records[0]["start_0based"] == 10
    nested = by_symbol.records[0]["Transcript"]
    assert isinstance(nested, tuple)
    assert isinstance(nested[0], Mapping)
    assert nested[0]["start_0based"] == 11
    assert "/lookup/id/ENSG1" in calls[0].full_url
    assert "/lookup/symbol/human/BRCA2" in calls[1].full_url


def test_ucsc_sequence_tracks_chromosomes_and_file_catalog(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener(
            [
                {
                    "genome": "hg38",
                    "chromCount": 2,
                    "dataTime": "2026-01-01",
                    "chromosomes": {"chrM": 16569, "chr1": 248956422},
                },
                {"genome": "hg38", "chrom": "chrM", "start": 0, "end": 4, "dna": "GATC"},
                {
                    "dataTime": "2026-01-01",
                    "hg38": {
                        "cpgIslandExt": {
                            "shortLabel": "CpG Islands",
                            "type": "bed 4 +",
                        }
                    },
                },
                {
                    "track": "cpgIslandExt",
                    "trackType": "bed 4 +",
                    "itemsReturned": 1,
                    "cpgIslandExt": [{"chrom": "chr1", "chromStart": 100, "chromEnd": 200}],
                },
                {
                    "itemsReturned": 2,
                    "maxItemsLimit": False,
                    "urlList": [
                        {"url": "goldenPath/hg38/bigZips/hg38.fa.gz", "sizeBytes": 10},
                        {
                            "url": "goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",
                            "sizeBytes": 20,
                        },
                    ],
                },
            ],
            calls,
        ),
    )

    chroms = ucsc_chromosomes("hg38", api_base_url="https://example.test")
    bases = ucsc_sequence("hg38", "chrM:0-4", api_base_url="https://example.test")
    track_list = ucsc_tracks("hg38", contains="CpG", api_base_url="https://example.test")
    features = ucsc_track_data(
        "hg38",
        "cpgIslandExt",
        "chr1:0-1000",
        api_base_url="https://example.test",
    )
    files = ucsc_files(
        "hg38",
        pattern="*chain*",
        limit=2,
        api_base_url="https://example.test",
    )

    assert chroms.records[0]["name"] == "chr1"
    assert bases.records[0]["seq"] == "GATC"
    assert track_list.records[0]["track"] == "cpgIslandExt"
    assert features.records[0]["start_0based"] == 100
    file_url = files.records[0]["url"]
    assert isinstance(file_url, str)
    assert file_url.endswith("chain.gz")
    assert parse_qs(urlsplit(calls[1].full_url).query)["revComp"] == ["false"]


def test_ncbi_virus_query_uses_accession_endpoint(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[Any] = []
    monkeypatch.setattr(
        remote_http,
        "urlopen",
        _json_opener([{"reports": [{"accession": "NC_045512.2"}], "total_count": 1}], calls),
    )

    result = virus(
        "NC_045512.2",
        refseq_only=True,
        page_size=1,
        api_base_url="https://example.test/v2",
    )

    assert result.records[0]["accession"] == "NC_045512.2"
    assert "/virus/accession/NC_045512.2/dataset_report" in calls[0].full_url
    assert parse_qs(urlsplit(calls[0].full_url).query)["refseq_only"] == ["true"]
