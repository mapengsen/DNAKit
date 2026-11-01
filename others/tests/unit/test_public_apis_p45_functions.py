"""Offline UCSC/ENCODE/RCSB request contracts and limits; no live service calls."""

from __future__ import annotations

import io
import json
from urllib.parse import parse_qs, urlsplit

import pytest

from dnakit.databases import _network_p44 as transport
from dnakit.databases import public_apis_p45_functions as api
from dnakit.exceptions import ConfigurationError


class Response(io.BytesIO):
    def __init__(self, data, url):
        super().__init__(data)
        self.url = url
        self.headers = {"Date": "Tue, 20 Oct 2026 00:00:00 GMT"}

    def geturl(self):
        return self.url


@pytest.mark.parametrize(
    "function,kwargs,path",
    [
        (api.ucsc_track_schema, {"genome": "hg38", "track": "knownGene"}, "/list/schema"),
        (
            api.ucsc_search_assembly,
            {"genome": "hg38", "term": "BRCA1", "categories": "trackDb"},
            "/search",
        ),
        (api.encode_object_metadata, {"accession": "ENCFF000AAA"}, "/ENCFF000AAA/"),
        (api.encode_object_metadata, {"accession": "/files/ENCFF000AAA/"}, "/files/ENCFF000AAA/"),
        (api.encode_experiment_matrix, {"filters": {"assay_title": "ATAC-seq"}}, "/matrix/"),
        (api.encode_body_system_summary, {}, "/summary/"),
    ],
)
@pytest.mark.parametrize("payload", [{}, {"results": []}, {"results": [{"id": "x"}]}])
def test_fixed_json_routes_preserve_data(monkeypatch, function, kwargs, path, payload):
    seen = []

    def opener(request, **kw):
        seen.append(request)
        return Response(json.dumps(payload).encode(), request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = function(**kwargs, verbose=False)
    assert urlsplit(seen[0].full_url).path == path
    assert result["data"] == payload
    json.dumps(result, allow_nan=False)


def test_encode_report_field_list_is_repeated(monkeypatch):
    seen = []

    def opener(request, **kw):
        seen.append(request.full_url)
        return Response(b"accession\tassembly\nA\thg38\n", request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    result = api.encode_tabular_report(
        ["accession", "assembly"], filters={"type": "File"}, verbose=False
    )
    assert parse_qs(urlsplit(seen[0]).query)["field"] == ["accession", "assembly"]
    assert result["data"].startswith("accession\tassembly")


def test_encode_files_list_is_not_claimed_download(monkeypatch):
    monkeypatch.setattr(
        transport,
        "urlopen",
        lambda request, **kw: Response(
            b"https://example.org/metadata.tsv\nhttps://example.org/reads.bam\n", request.full_url
        ),
    )
    result = api.encode_download_urls(filters={"type": "File", "assembly": "hg38"}, verbose=False)
    assert len(result["urls"]) == 2
    assert "file" not in result and "outputs" not in result


@pytest.mark.parametrize("sequence_type", ["dna", "rna"])
def test_rcsb_sequence_query_has_original_defaults(monkeypatch, sequence_type):
    seen = []

    def opener(request, **kw):
        seen.append(request)
        return Response(b'{"total_count":0,"result_set":[]}', request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    sequence = ("ACGT" if sequence_type == "dna" else "ACGU") * 7
    result = api.rcsb_nucleotide_sequence_search(
        sequence, sequence_type=sequence_type, verbose=False
    )
    body = json.loads(seen[0].data)
    assert body["query"]["parameters"] == {
        "value": sequence,
        "sequence_type": sequence_type,
        "identity_cutoff": 0.0,
        "evalue_cutoff": 0.1,
    }
    assert body["return_type"] == "polymer_entity"
    assert result["data"]["total_count"] == 0


@pytest.mark.parametrize("pattern_type", ["simple", "prosite", "regex"])
@pytest.mark.parametrize("sequence_type", ["dna", "rna"])
def test_rcsb_motif_request(monkeypatch, pattern_type, sequence_type):
    seen = []

    def opener(request, **kw):
        seen.append(request)
        return Response(b"{}", request.full_url)

    monkeypatch.setattr(transport, "urlopen", opener)
    api.rcsb_nucleotide_motif_search(
        "ACG", pattern_type=pattern_type, sequence_type=sequence_type, verbose=False
    )
    body = json.loads(seen[0].data)
    assert body["query"]["service"] == "seqmotif"
    assert body["query"]["parameters"]["sequence_type"] == sequence_type
    assert body["query"]["parameters"]["pattern_type"] == pattern_type


@pytest.mark.parametrize(
    "kwargs",
    [
        {"sequence": "ACG"},
        {"sequence": "A" * 25, "identity_cutoff": float("nan")},
        {"sequence": "A" * 25, "evalue_cutoff": 0},
        {"sequence": "A" * 25, "sequence_type": "protein"},
        {"sequence": "A" * 25, "rows": 0},
    ],
)
def test_rcsb_invalid_inputs_do_not_request(kwargs):
    with pytest.raises(ConfigurationError):
        api.rcsb_nucleotide_sequence_search(**kwargs, verbose=False)


def test_unknown_encode_filters_rejected():
    with pytest.raises(ConfigurationError):
        api.encode_experiment_matrix(filters={"arbitrary_api_parameter": "x"}, verbose=False)
