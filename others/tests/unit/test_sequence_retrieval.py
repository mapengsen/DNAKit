"""Native sequence workflows: actual HTTP boundaries, files, and Agent calls."""

from __future__ import annotations

import builtins
import hashlib
import io
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlsplit

import pytest

import dnakit
import dnakit.search._http as remote_http
from dnakit.download import DownloadConfig
from dnakit.exceptions import ConfigurationError, DownloadError, QueryError
from dnakit.search import SearchConfig
from dnakit.tools import default_tool_registry

FASTA = ">NC_000913.3 fixture\nACGTACGTACGT\n"
GENBANK = """LOCUS       TEST                      12 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  Synthetic test fixture.
ACCESSION   NC_000913
VERSION     NC_000913.3
KEYWORDS    .
SOURCE      Escherichia coli
  ORGANISM  Escherichia coli
            Bacteria.
FEATURES             Location/Qualifiers
     source          1..12
                     /organism="Escherichia coli"
     gene            1..12
                     /gene="test"
ORIGIN
        1 acgtacgtacgt
//
"""


class Response:
    def __init__(self, payload: bytes, url: str):
        self.stream, self.url = io.BytesIO(payload), url

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        self.stream.close()

    def read(self, size=-1):
        return self.stream.read(size)

    def geturl(self):
        return self.url


def stub(monkeypatch, payloads):
    calls = []

    def open_fake(request, *, timeout):
        assert timeout > 0
        calls.append(request.full_url)
        payload = payloads.pop(0)
        data = payload.encode() if isinstance(payload, str) else json.dumps(payload).encode()
        return Response(data, request.full_url)

    monkeypatch.setattr(remote_http, "urlopen", open_fake)
    return calls


def test_structured_search_and_legacy_query_are_compatible(monkeypatch):
    calls = stub(
        monkeypatch,
        [
            {"esearchresult": {"idlist": [], "count": "0"}},
            {"esearchresult": {"idlist": [], "count": "0"}},
        ],
    )
    result = dnakit.search.ncbi.entrez(
        "nuccore",
        "BRCA1 OR BRCA2",
        organism="Homo sapiens",
        gene="BRCA1",
        strain="isolate A",
        seq_type="mrna",
        source="refseq",
        retmax=3,
        config=SearchConfig(api_key="private-key", email="private@example.test"),
    )
    query = parse_qs(urlsplit(calls[0]).query)
    assert query["term"] == [
        '(BRCA1 OR BRCA2) AND ("Homo sapiens"[Organism]) AND ("BRCA1"[Gene Name]) '
        'AND ("isolate A"[All Fields]) AND ("biomol mrna"[Properties]) AND (refseq[filter])'
    ]
    assert query["retmax"] == ["3"]
    assert result.metadata["uids"] == () and result.total_count == 0
    assert "private-key" not in json.dumps(result.to_dict())
    dnakit.search.entrez("pubmed", "DNA repair OR replication")
    assert parse_qs(urlsplit(calls[1]).query)["term"] == ["DNA repair OR replication"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"database": "pubmed", "organism": "human"},
        {"seq_type": "protein"},
        {"seq_type": "refseq", "source": "genbank"},
        {"organism": 'human" OR virus'},
        {"organism": "human", "seq_type": "typo"},
        {"organism": "human", "source": "typo"},
    ],
)
def test_invalid_search_fails_before_http(monkeypatch, kwargs):
    calls = stub(monkeypatch, [])
    with pytest.raises(ConfigurationError):
        dnakit.search.entrez(**{"database": "nuccore", **kwargs})
    assert calls == []


def test_accessions_preserve_uid_order_and_duplicates(monkeypatch):
    calls = stub(
        monkeypatch,
        [
            {
                "result": {
                    "uids": ["2", "1"],
                    "1": {"accessionversion": "NC_000913.3"},
                    "2": {"accessionversion": "NZ_CP000001.1"},
                }
            }
        ],
    )
    assert dnakit.search.ncbi.accessions([1, 2, "1"]) == (
        "NC_000913.3",
        "NZ_CP000001.1",
        "NC_000913.3",
    )
    assert parse_qs(urlsplit(calls[0]).query)["id"] == ["1,2"]
    assert dnakit.search.ncbi.accessions([]) == ()
    assert len(calls) == 1


def test_accessions_batch_at_200(monkeypatch):
    ids = [str(i) for i in range(1, 202)]

    def payload(batch):
        return {"result": {uid: {"accessionversion": f"AB{uid}.1"} for uid in batch}}

    calls = stub(monkeypatch, [payload(ids[:200]), payload(ids[200:])])
    assert len(dnakit.search.accessions(ids)) == 201
    assert len(parse_qs(urlsplit(calls[0]).query)["id"][0].split(",")) == 200
    assert parse_qs(urlsplit(calls[1]).query)["id"] == ["201"]


def test_unresolved_uid_is_an_error(monkeypatch):
    stub(monkeypatch, [{"result": {"1": {"error": "UID not found"}}}])
    with pytest.raises(QueryError):
        dnakit.search.accessions([1])


@pytest.mark.parametrize("uids", ["123", [True], ["1,2"], ["x"], [-1]])
def test_invalid_uids_do_not_reach_network(monkeypatch, uids):
    calls = stub(monkeypatch, [])
    with pytest.raises(ConfigurationError):
        dnakit.search.accessions(uids)
    assert not calls


def test_download_saves_sequence_and_redacted_manifest(monkeypatch, tmp_path):
    calls = stub(monkeypatch, [FASTA])
    target = tmp_path / "sequence.fasta"
    result = dnakit.download.ncbi.sequence(
        "NC_000913.3",
        target,
        seq_start=1,
        seq_stop=12,
        strand=2,
        search_config=SearchConfig(api_key="secret-key", email="secret@example.test"),
    )
    params = parse_qs(urlsplit(calls[0]).query)
    assert params["id"] == ["NC_000913.3"]
    assert params["seq_start"] == ["1"] and params["seq_stop"] == ["12"]
    assert params["strand"] == ["2"] and params["rettype"] == ["fasta"]
    assert target.read_text() == FASTA
    manifest = Path(result.manifest_path).read_text()
    assert "secret-key" not in manifest and "secret@example.test" not in manifest
    assert result.files[0].sha256 == hashlib.sha256(FASTA.encode()).hexdigest()
    assert json.loads(manifest)["metadata"]["coordinate_system"] == "1-based inclusive"
    with pytest.raises(FileExistsError):
        dnakit.download.ncbi.sequence("NC_000913.3", target)
    assert len(calls) == 1  # Protected before another request.


@pytest.mark.parametrize(
    "payload", ["", "Error: invalid ID", "<html>unavailable</html>", ">id\n", ">id\nACGT\n>empty\n"]
)
def test_provider_errors_never_install_a_sequence(monkeypatch, tmp_path, payload):
    stub(monkeypatch, [payload])
    with pytest.raises(DownloadError):
        dnakit.download.ncbi.sequence("NC_000913.3", tmp_path / "invalid.fasta")
    assert list(tmp_path.iterdir()) == []


def test_response_limit_preserves_existing_file(monkeypatch, tmp_path):
    stub(monkeypatch, [FASTA])
    target = tmp_path / "old.fasta"
    target.write_text("original")
    with pytest.raises(QueryError, match="max_response_bytes"):
        dnakit.download.ncbi.sequence(
            "NC_000913.3",
            target,
            config=DownloadConfig(overwrite=True),
            search_config=SearchConfig(max_response_bytes=5),
        )
    assert target.read_text() == "original"
    assert not target.with_name("old.fasta.manifest.json").exists()


@pytest.mark.parametrize(
    "text,format",
    [
        (FASTA.replace("NC_000913.3", "NC_000913.2"), "fasta"),
        (GENBANK.replace("12 bp", "13 bp"), "genbank"),
        (GENBANK.replace("NC_000913.3", "NC_045512.2"), "genbank"),
    ],
)
def test_identity_and_length_mismatch_preserve_existing_file(monkeypatch, tmp_path, text, format):
    stub(monkeypatch, [text])
    target = tmp_path / "existing"
    target.write_text("original")
    with pytest.raises(DownloadError):
        dnakit.download.ncbi.sequence(
            "NC_000913.3", target, format=format, config=DownloadConfig(overwrite=True)
        )
    assert target.read_text() == "original"
    assert not (tmp_path / "existing.manifest.json").exists()


def test_invalid_entrez_query_is_not_an_empty_success(monkeypatch):
    stub(monkeypatch, [{"esearchresult": {"idlist": [], "errorlist": {"FieldNotFound": ["typo"]}}}])
    with pytest.raises(QueryError):
        dnakit.search.entrez("nuccore", "x[typo]")


def test_ncbi_annotation_and_embl_conversion(monkeypatch, tmp_path):
    SeqIO = pytest.importorskip("Bio.SeqIO")
    calls = stub(monkeypatch, [GENBANK, GENBANK])
    gb = dnakit.download.ncbi.sequence("NC_000913.3", tmp_path / "seq.gb", format="genbank")
    embl = dnakit.download.ncbi.sequence("NC_000913.3", tmp_path / "seq.embl", format="embl")
    assert parse_qs(urlsplit(calls[0]).query)["rettype"] == ["gbwithparts"]
    original = SeqIO.read(gb.files[0].path, "genbank")
    converted = SeqIO.read(embl.files[0].path, "embl")
    assert str(original.seq) == str(converted.seq) == "ACGTACGTACGT"
    assert original.id == converted.id
    assert [f.type for f in converted.features] == ["source", "gene"]
    assert embl.metadata["converted"] is True


def test_missing_conversion_dependency_is_actionable(monkeypatch, tmp_path):
    calls = stub(monkeypatch, [])
    real_import = builtins.__import__

    def without_bio(name, *args, **kwargs):
        if name == "Bio":
            raise ImportError("fixture")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_bio)
    with pytest.raises(ConfigurationError, match="bio-formats"):
        dnakit.download.ncbi.sequence("NC_000913.3", tmp_path / "seq.embl", format="embl")
    assert not calls


def test_ena_native_and_converted_formats(monkeypatch, tmp_path):
    SeqIO = pytest.importorskip("Bio.SeqIO")
    output = io.StringIO()
    SeqIO.convert(io.StringIO(GENBANK), "genbank", output, "embl")
    ena_embl = output.getvalue().replace("NC_000913", "U00096")
    calls = stub(monkeypatch, [FASTA.replace("NC_000913", "U00096"), ena_embl, ena_embl])
    dnakit.download.ena.sequence("U00096.3", tmp_path / "ena.fasta")
    dnakit.download.ena.sequence("U00096.3", tmp_path / "ena.embl", format="embl")
    gb = dnakit.download.ena.sequence("U00096.3", tmp_path / "ena.gb", format="genbank")
    assert calls[0].endswith("/fasta/U00096.3")
    assert calls[1].endswith("/embl/U00096.3")
    assert str(SeqIO.read(gb.files[0].path, "genbank").seq) == "ACGTACGTACGT"
    assert gb.metadata["converted"] is True


@pytest.mark.parametrize(
    "accession", ["NC_000913.3", "NZ_CP000001.1", "GCF_000005845.2", "U00096.3,AB123.1"]
)
def test_ena_rejects_refseq_and_invalid_ids(monkeypatch, tmp_path, accession):
    calls = stub(monkeypatch, [])
    with pytest.raises(ConfigurationError):
        dnakit.download.ena.sequence(accession, tmp_path / "ena.fasta")
    assert not calls


def test_protein_uses_protein_database_and_genpept(monkeypatch, tmp_path):
    calls = stub(monkeypatch, [GENBANK.replace("NC_000913.3", "NP_009225.1")])
    dnakit.download.ncbi.sequence(
        "NP_009225.1", tmp_path / "seq.gp", database="protein", format="genbank"
    )
    params = parse_qs(urlsplit(calls[0]).query)
    assert params["db"] == ["protein"] and params["rettype"] == ["gp"]
    assert "strand" not in params
    with pytest.raises(ConfigurationError):
        dnakit.download.ncbi.sequence(
            "NP_009225.1", tmp_path / "seq.embl", database="protein", format="embl"
        )


def test_agent_can_describe_and_execute_native_download(monkeypatch, tmp_path):
    stub(monkeypatch, [FASTA])
    registry = default_tool_registry()
    search = registry.describe("dnakit.search.entrez")
    assert search["agent_compatible"]
    properties = search["input_schema"]["properties"]
    assert {"organism", "gene", "strain", "seq_type", "source"} <= properties.keys()
    assert "term" not in search["input_schema"].get("required", [])
    for name in (
        "dnakit.search.accessions",
        "dnakit.download.ncbi_sequence",
        "dnakit.download.ena_sequence",
    ):
        assert registry.describe(name)["agent_compatible"]
    # Exercise the same coercion/serialization path used by MCP execution.
    result = registry.execute(
        "dnakit.download.ncbi_sequence",
        {
            "accession": "NC_000913.3",
            "output_path": str(tmp_path / "agent.fasta"),
        },
        allow_side_effects=True,
    )
    assert result["kind"] == "accession_sequence"
    assert (tmp_path / "agent.fasta").read_text() == FASTA


def test_skill_examples_run_independently(monkeypatch, tmp_path):
    SeqIO = pytest.importorskip("Bio.SeqIO")
    output = io.StringIO()
    SeqIO.convert(io.StringIO(GENBANK), "genbank", output, "embl")

    def open_fake(request: Any, *, timeout: float):
        url = request.full_url
        parts, params = urlsplit(url), parse_qs(urlsplit(url).query)
        if parts.path.endswith("esearch.fcgi"):
            value = {"esearchresult": {"idlist": ["1"], "count": "1"}}
        elif parts.path.endswith("esummary.fcgi"):
            value = {"result": {"uids": ["1"], "1": {"accessionversion": "NC_000913.3"}}}
        elif parts.path.endswith("efetch.fcgi"):
            value = FASTA if params["rettype"] == ["fasta"] else GENBANK
            value = value.replace("NC_000913.3", params["id"][0])
        elif "/fasta/" in parts.path:
            value = FASTA.replace("NC_000913.3", parts.path.rsplit("/", 1)[1])
        elif "/embl/" in parts.path:
            value = output.getvalue().replace("NC_000913", "U00096")
        else:
            value = [{"accession": "U00096", "description": "fixture"}]
        payload = value.encode() if isinstance(value, str) else json.dumps(value).encode()
        return Response(payload, url)

    monkeypatch.setattr(remote_http, "urlopen", open_fake)
    root = Path(__file__).resolve().parents[3]
    examples = (root / "others/skills/adaptations/sequence-retrieval/examples.md").read_text()
    for index, code in enumerate(re.findall(r"```python\n(.*?)```", examples, re.S)):
        workspace = tmp_path / str(index)
        workspace.mkdir()
        monkeypatch.chdir(workspace)
        exec(compile(code, f"sequence-retrieval/examples.md:example-{index + 1}", "exec"), {})
        assert list(workspace.glob("*.manifest.json")) or index == 3
