"""Contract tests for optional service clients and the local Toil job graph."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import pytest

from dnakit.core import DNA
from dnakit.exceptions import (
    BackendExecutionError,
    BackendUnavailableError,
    ConfigurationError,
    InputFormatError,
    QueryError,
)
from dnakit.integrations import service_backends as backends
from dnakit.integrations.service_backends import (
    analyze_fasta_with_clawbio,
    compute_toil_sequence_stats,
    describe_smartapi_api,
    get_biomcp_gene,
    get_biomcp_variant,
    list_galaxy_tools,
    search_biomcp_genes,
    search_smartapi_registry,
    show_galaxy_tool,
)

LIVE = pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="set DNAKIT_LIVE_TESTS=1 to run live network tests",
)
POSIX_ONLY = pytest.mark.skipif(os.name != "posix", reason="fake CLI is a POSIX shell script")
MYGENE_ID = "59dce17363dce279d389100834e43648"
GO_ID = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
PROTEIN_ID = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"

BRCA1_GENE = {
    "symbol": "BRCA1",
    "name": "BRCA1 DNA repair associated",
    "entrez_id": "672",
    "ensembl_id": "ENSG00000012048",
    "uniprot_id": "P38398",
    "omim_id": "113705",
    "gene_type": "protein-coding",
    "location": "17",
    "genomic_coordinates": {
        "coordinate": "17:43044292-43170245 (strand: -1)",
        "genome_build": "GRCh38",
        "source": "MyGene.info",
    },
    "aliases": ["BRCAI", "BRCC1"],
    "clinical_diseases": ["breast cancer"],
    "summary": "This gene encodes a 190 kD nuclear phosphoprotein.",
}
BRCA1_SEARCH = {
    "pagination": {"offset": 0, "limit": 2, "returned": 1, "total": 1, "has_more": False},
    "count": 1,
    "results": [
        {
            "symbol": "BRCA1",
            "name": "BRCA1 DNA repair associated",
            "entrez_id": "672",
            "genomic_coordinates": {"coordinate": "17:43044292-43170245", "genome_build": "GRCh38"},
            "uniprot_id": "P38398",
            "omim_id": "113705",
        }
    ],
}
RS334_VARIANT = {
    "id": "chr11:g.5248232T>A",
    "rsid": "rs334",
    "gene": "HBB",
    "hgvs_c": "c.20A>T",
    "hgvs_p": "p.Glu7Val",
    "consequence": "missense_variant",
    "significance": "Pathogenic",
    "transcript": "NM_000518.4",
    "genome_build": "GRCh37",
    "cadd_score": 7.277,
    "sift_pred": "Deleterious",
    "polyphen_pred": "Benign",
}
CHR7_VARIANT = {
    "id": "chr7:g.140753336A>T",
    "rsid": "rs113488022",
    "gene": "BRAF",
    "hgvs_p": "p.Val600Glu",
    "consequence": None,
    "significance": "Pathogenic",
    "genome_build": "GRCh38",
}
BIOMCP_NOT_FOUND = {
    "error": {
        "code": "not_found",
        "message": "gene 'NOTAGENEXYZ' not found.\n\nTry searching: biomcp search gene -q X",
    },
    "_meta": {"not_found": True},
}


def _smartapi_hit(identifier: str, title: str, description: str, tags: list[str]) -> dict[str, Any]:
    return {
        "_id": identifier,
        "info": {"title": title, "description": description},
        "tags": [{"name": tag} for tag in tags],
    }


SMARTAPI_HITS = [
    _smartapi_hit(MYGENE_ID, "MyGene.info API", "gene annotation service", ["biothings"]),
    _smartapi_hit(GO_ID, "Gene Ontology API", "gene ontology terms", []),
    _smartapi_hit(PROTEIN_ID, "Protein Atlas API", "protein expression", []),
]

SMARTAPI_SPEC = {
    "openapi": "3.0.3",
    "info": {"title": "Fixture Gene API", "version": "1.0"},
    "servers": [{"url": "https://example.org/v1"}],
    "paths": {
        "/query": {
            "get": {"operationId": "queryGenes", "summary": "Query genes"},
            "post": {"operationId": "postQuery"},
        },
        "/gene/{id}": {"get": {"operationId": "getGene", "summary": "Get one gene"}},
        "/meta": {"parameters": [], "get": {"operationId": "getMeta"}},
    },
}

GALAXY_TOOLS = [
    {
        "id": "bwa_mem",
        "name": "Map with BWA-MEM",
        "version": "0.7.17",
        "description": "aligns reads",
        "model_class": "Tool",
    },
    {
        "id": "Cut1",
        "name": "Cut",
        "version": "1.0.2",
        "description": "columns from a table",
        "model_class": "Tool",
    },
    {"id": "Sort1", "name": "Sort", "version": "1.0.0", "description": "", "model_class": "Tool"},
]


def _mock_httpx(monkeypatch: pytest.MonkeyPatch, handler: Any) -> None:
    import httpx

    real_client = httpx.AsyncClient

    def factory(*args: Any, **kwargs: Any) -> Any:
        kwargs["transport"] = httpx.MockTransport(handler)
        return real_client(*args, **kwargs)

    monkeypatch.setattr(httpx, "AsyncClient", factory)


def _fake_biomcp(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    response: str,
    *,
    status: int = 0,
) -> Path:
    """Put a POSIX stub named biomcp first on PATH; return the file that records its arguments."""

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    response_file = tmp_path / "response.json"
    response_file.write_text(response, encoding="utf-8")
    args_file = tmp_path / "args.txt"
    script = bin_dir / "biomcp"
    script.write_text(
        f'#!/bin/sh\nprintf "%s\\n" "$@" > "{args_file}"\ncat "{response_file}"\nexit {status}\n',
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return args_file


@pytest.fixture
def galaxy_server(monkeypatch: pytest.MonkeyPatch):
    """A localhost stub of the Galaxy tool endpoints, recording requests."""

    seen: list[dict[str, Any]] = []
    state = {"tools_status": 200}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:
            return

        def do_GET(self) -> None:
            parsed = urlsplit(self.path)
            seen.append({"path": parsed.path, "key": self.headers.get("x-api-key")})
            if parsed.path == "/api/tools":
                status, body = state["tools_status"], GALAXY_TOOLS
            elif parsed.path == "/api/tools/Cut1":
                status, body = 200, GALAXY_TOOLS[1]
            else:
                status, body = 404, {"err_msg": "not found"}
            data = json.dumps(body).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost")
    monkeypatch.setenv("no_proxy", "127.0.0.1,localhost")
    try:
        yield {"url": f"http://127.0.0.1:{server.server_address[1]}/", "seen": seen, "state": state}
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


# Toil: two chained jobs over DNA strings, checked against Biopython and hand values.


def test_toil_chain_matches_biopython_and_hand_values() -> None:
    pytest.importorskip("toil.job")
    pytest.importorskip("toil.common")
    seq_module = pytest.importorskip("Bio.Seq")
    seq_utils = pytest.importorskip("Bio.SeqUtils")
    result = compute_toil_sequence_stats(
        {"first": "ATGCGCAATT", "second": "GGGCCCAAAT", "third": "GCNN"}
    )
    rows = result["sequences"]
    assert result["sequence_count"] == 3
    assert [row["id"] for row in rows] == ["first", "second", "third"]
    # Hand values: ATGCGCAATT reverses to TTAACGCGTA and complements to AATTGCGCAT;
    # GC is 4 of 10 (0.4) and 6 of 10 (0.6). GCNN has 2 of 4 GC characters (0.5).
    assert [row["reverse_complement"] for row in rows] == ["AATTGCGCAT", "ATTTGGGCCC", "NNGC"]
    assert [row["gc_fraction"] for row in rows] == pytest.approx([0.4, 0.6, 0.5])
    assert [row["length"] for row in rows] == [10, 10, 4]
    for row, source in zip(rows, ["ATGCGCAATT", "GGGCCCAAAT", "GCNN"], strict=True):
        assert row["reverse_complement"] == str(seq_module.reverse_complement(source))
        assert row["gc_fraction"] == pytest.approx(
            float(seq_utils.gc_fraction(source, ambiguous="ignore"))
        )
    assert result["provenance"]["package"].startswith("toil ")
    # ACGTTG reverses to GTTGCA and complements to CAACGT; GC is 3 of 6 (0.5).
    dna_result = compute_toil_sequence_stats(DNA([{"id": "dna_one", "sequence": "ACGTTG"}]))
    assert dna_result["sequences"] == [
        {"id": "dna_one", "length": 6, "reverse_complement": "CAACGT", "gc_fraction": 0.5}
    ]


def test_toil_rejects_invalid_sequence_input() -> None:
    for bad in [{}, {"a": ""}, {"a": "ACGTX"}, {"a": "AC GT"}, {"a\x00": "ACGT"}, 5, "ACGT"]:
        with pytest.raises(ConfigurationError):
            compute_toil_sequence_stats(bad)


def test_toil_enforces_sequence_and_base_limits() -> None:
    with pytest.raises(ConfigurationError, match="max_sequences"):
        compute_toil_sequence_stats({"a": "ACGT", "b": "ACGT"}, max_sequences=1)
    with pytest.raises(ConfigurationError, match="max_total_bases"):
        compute_toil_sequence_stats({"a": "ACGT"}, max_total_bases=3)
    with pytest.raises(ConfigurationError):
        compute_toil_sequence_stats({"a": "ACGT"}, max_sequences=True)


def test_toil_missing_package_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("toil", "toil.job", "toil.common"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(BackendUnavailableError, match="bio-toil"):
        compute_toil_sequence_stats({"a": "ACGT"})


# SmartAPI: real smartapi-mcp parsing and ranking, with the HTTP layer mocked.


def test_smartapi_search_ranks_mocked_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    import httpx

    pytest.importorskip("smartapi_mcp.smartapi")
    biothings = pytest.importorskip("smartapi_mcp.biothings")
    requested: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        requested.append(request.url.params["q"])
        return httpx.Response(200, json={"hits": SMARTAPI_HITS})

    _mock_httpx(monkeypatch, handler)
    result = search_smartapi_registry("Gene Annotation", limit=5)
    assert requested == ["gene annotation"]
    assert result["query"] == "Gene Annotation"
    assert result["registry_query"] == "gene annotation"
    assert result["registry_hits"] == 3
    # Protein Atlas shares no query term and is omitted. MyGene matches both terms and
    # is a core API, so it outranks Gene Ontology, which matches only "gene" in its name.
    assert [match["smartapi_id"] for match in result["matches"]] == [MYGENE_ID, GO_ID]
    entries = [
        {
            "_id": MYGENE_ID,
            "title": "MyGene.info API",
            "description": "gene annotation service",
            "tags": ["biothings"],
        },
        {
            "_id": GO_ID,
            "title": "Gene Ontology API",
            "description": "gene ontology terms",
            "tags": [],
        },
        {
            "_id": PROTEIN_ID,
            "title": "Protein Atlas API",
            "description": "protein expression",
            "tags": [],
        },
    ]
    direct = biothings.rank_apis(
        biothings.build_registry_from_entries(entries), keyword="gene annotation", limit=5
    )
    assert [match["name"] for match in result["matches"]] == [item["name"] for item in direct]
    assert [match["score"] for match in result["matches"]] == pytest.approx(
        [item["score"] for item in direct]
    )


def test_smartapi_search_reports_registry_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    import httpx

    pytest.importorskip("smartapi_mcp.smartapi")

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(503, json={"error": "unavailable"})

    _mock_httpx(monkeypatch, handler)
    with pytest.raises(QueryError) as excinfo:
        search_smartapi_registry("gene")
    assert excinfo.value.context["error_type"] == "HTTPStatusError"


def test_smartapi_search_validates_before_network() -> None:
    with pytest.raises(ConfigurationError):
        search_smartapi_registry("")
    with pytest.raises(ConfigurationError):
        search_smartapi_registry("gene", limit=0)


def test_smartapi_missing_package_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "smartapi_mcp.smartapi", None)
    monkeypatch.setitem(sys.modules, "smartapi_mcp.biothings", None)
    with pytest.raises(BackendUnavailableError, match="bio-smartapi"):
        search_smartapi_registry("gene")


def test_smartapi_describe_lists_operations_from_spec(monkeypatch: pytest.MonkeyPatch) -> None:
    smartapi_module = pytest.importorskip("smartapi_mcp.smartapi")
    urls: list[str] = []

    def fake_fetch(url: str, **_: Any) -> dict[str, Any]:
        urls.append(url)
        return SMARTAPI_SPEC

    monkeypatch.setattr(smartapi_module, "fetch_spec", fake_fetch)
    result = describe_smartapi_api(GO_ID)
    assert urls == [f"https://smart-api.info/api/metadata/{GO_ID}"]
    assert result["title"] == "Fixture Gene API"
    assert result["version"] == "1.0"
    assert result["openapi"] == "3.0.3"
    assert result["server_url"] == "https://example.org/v1"
    # Four operations in total: GET /gene/{id}, GET /meta, GET /query, POST /query.
    assert result["operation_count"] == 4
    assert result["truncated"] is False
    assert [(op["method"], op["path"]) for op in result["operations"]] == [
        ("GET", "/gene/{id}"),
        ("GET", "/meta"),
        ("GET", "/query"),
        ("POST", "/query"),
    ]
    assert result["operations"][2]["operation_id"] == "queryGenes"

    truncated = describe_smartapi_api(GO_ID, max_operations=2)
    assert truncated["truncated"] is True
    assert len(truncated["operations"]) == 2


def test_smartapi_describe_reports_spec_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    smartapi_module = pytest.importorskip("smartapi_mcp.smartapi")

    def failing_fetch(url: str, **_: Any) -> dict[str, Any]:
        raise ValueError("spec too large")

    monkeypatch.setattr(smartapi_module, "fetch_spec", failing_fetch)
    with pytest.raises(QueryError, match="could not be loaded"):
        describe_smartapi_api(GO_ID)


def test_smartapi_describe_validates_identifier() -> None:
    with pytest.raises(ConfigurationError):
        describe_smartapi_api("not-a-smartapi-id")
    with pytest.raises(ConfigurationError):
        describe_smartapi_api(GO_ID.upper())


# biomcp: a fixed-argument CLI call through the bounded executor, with a POSIX stub on PATH.


@POSIX_ONLY
def test_biomcp_gene_uses_fixed_arguments_and_parses_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args_file = _fake_biomcp(tmp_path, monkeypatch, json.dumps(BRCA1_GENE))
    result = get_biomcp_gene("BRCA1")
    assert args_file.read_text(encoding="utf-8").split() == ["-j", "get", "gene", "BRCA1"]
    assert result["symbol"] == "BRCA1"
    assert result["entrez_id"] == "672"
    assert result["location"] == "17"
    assert result["coordinate"] == "17:43044292-43170245 (strand: -1)"
    assert result["genome_build"] == "GRCh38"
    assert result["aliases"] == ["BRCAI", "BRCC1"]
    assert result["provenance"]["package"].startswith("biomcp-cli ")


@POSIX_ONLY
def test_biomcp_search_and_variant_parse_cli_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args_file = _fake_biomcp(tmp_path, monkeypatch, json.dumps(BRCA1_SEARCH))
    search = search_biomcp_genes("BRCA1", limit=2)
    assert args_file.read_text(encoding="utf-8").split() == [
        "-j",
        "search",
        "gene",
        "-q",
        "BRCA1",
        "--limit",
        "2",
    ]
    assert search["total"] == 1
    assert search["genes"] == [
        {
            "symbol": "BRCA1",
            "name": "BRCA1 DNA repair associated",
            "entrez_id": "672",
            "uniprot_id": "P38398",
            "omim_id": "113705",
            "coordinate": "17:43044292-43170245",
        }
    ]

    args_file = _fake_biomcp(tmp_path, monkeypatch, json.dumps(RS334_VARIANT))
    variant = get_biomcp_variant("rs334")
    assert args_file.read_text(encoding="utf-8").split() == ["-j", "get", "variant", "rs334"]
    assert variant["gene"] == "HBB"
    assert variant["hgvs_p"] == "p.Glu7Val"
    assert variant["cadd_score"] == pytest.approx(7.277)
    assert variant["significance"] == "Pathogenic"


@POSIX_ONLY
def test_biomcp_variant_keeps_hgvs_characters_as_one_argument(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args_file = _fake_biomcp(tmp_path, monkeypatch, json.dumps(CHR7_VARIANT))
    variant = get_biomcp_variant("chr7:g.140753336A>T")
    assert args_file.read_text(encoding="utf-8").split() == [
        "-j",
        "get",
        "variant",
        "chr7:g.140753336A>T",
    ]
    assert variant["rsid"] == "rs113488022"
    assert variant["consequence"] is None


@POSIX_ONLY
def test_biomcp_reports_cli_error_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _fake_biomcp(tmp_path, monkeypatch, json.dumps(BIOMCP_NOT_FOUND), status=1)
    with pytest.raises(QueryError, match="not found") as excinfo:
        get_biomcp_gene("NOTAGENEXYZ")
    assert excinfo.value.code == "BIOMCP_NOT_FOUND"
    assert str(excinfo.value).count("\n") == 0


@POSIX_ONLY
def test_biomcp_rejects_non_json_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _fake_biomcp(tmp_path, monkeypatch, "plain text, not JSON")
    with pytest.raises(BackendExecutionError, match="did not return JSON"):
        get_biomcp_gene("BRCA1")


@POSIX_ONLY
def test_biomcp_validates_input_before_running(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    args_file = _fake_biomcp(tmp_path, monkeypatch, json.dumps(BRCA1_GENE))
    for bad_symbol in ["BRCA1; rm -rf /", "-j", "", "BRCA1 extra"]:
        with pytest.raises(ConfigurationError):
            get_biomcp_gene(bad_symbol)
    with pytest.raises(ConfigurationError):
        search_biomcp_genes("-leading-dash")
    with pytest.raises(ConfigurationError):
        get_biomcp_variant("rs1 --help")
    assert not args_file.exists()


def test_biomcp_missing_executable_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(backends, "_biomcp_executable", lambda: None)
    with pytest.raises(BackendUnavailableError, match="bio-biomcp"):
        get_biomcp_gene("BRCA1")


@LIVE
def test_live_biomcp_matches_direct_cli() -> None:
    executable = backends._biomcp_executable()
    if executable is None:
        pytest.skip("biomcp executable is not installed")
    direct = json.loads(
        subprocess.run(
            [executable, "-j", "get", "gene", "BRCA1"],
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
        ).stdout
    )
    wrapped = get_biomcp_gene("BRCA1")
    assert wrapped["entrez_id"] == direct["entrez_id"] == "672"
    assert wrapped["coordinate"] == direct["genomic_coordinates"]["coordinate"]
    variant = get_biomcp_variant("rs334")
    assert variant["hgvs_p"] == "p.Glu7Val"
    assert variant["rsid"] == "rs334"


# Galaxy: bioblend and galaxy-mcp's redirect-safe client against a localhost stub server.


def test_galaxy_lists_sorted_tools_and_sends_key_in_header(galaxy_server: dict[str, Any]) -> None:
    pytest.importorskip("galaxy_mcp.client")
    result = list_galaxy_tools(galaxy_server["url"], api_key="TESTKEY-123", limit=2)
    assert result["url"] == galaxy_server["url"]
    assert result["total_tools"] == 3
    assert result["returned"] == 2
    # Sorted by id with uppercase first: Cut1, Sort1, then bwa_mem.
    assert [tool["id"] for tool in result["tools"]] == ["Cut1", "Sort1"]
    assert result["tools"][0]["name"] == "Cut"
    tools_requests = [item for item in galaxy_server["seen"] if item["path"] == "/api/tools"]
    assert tools_requests[-1]["key"] == "TESTKEY-123"
    assert "TESTKEY-123" not in json.dumps(result)


def test_galaxy_anonymous_request_sends_no_key(galaxy_server: dict[str, Any]) -> None:
    pytest.importorskip("galaxy_mcp.client")
    result = list_galaxy_tools(galaxy_server["url"], limit=1)
    assert result["returned"] == 1
    assert all(item["key"] is None for item in galaxy_server["seen"])


def test_galaxy_show_tool_reads_detail_endpoint(galaxy_server: dict[str, Any]) -> None:
    pytest.importorskip("galaxy_mcp.client")
    result = show_galaxy_tool(galaxy_server["url"], tool_id="Cut1")
    assert result["tool"] == {
        "id": "Cut1",
        "name": "Cut",
        "version": "1.0.2",
        "description": "columns from a table",
        "model_class": "Tool",
    }
    assert any(item["path"] == "/api/tools/Cut1" for item in galaxy_server["seen"])


def test_galaxy_reports_server_failure(galaxy_server: dict[str, Any]) -> None:
    pytest.importorskip("galaxy_mcp.client")
    galaxy_server["state"]["tools_status"] = 500
    with pytest.raises(QueryError, match="did not return a tool list"):
        list_galaxy_tools(galaxy_server["url"])


def test_galaxy_rejects_bad_urls_before_any_request(galaxy_server: dict[str, Any]) -> None:
    pytest.importorskip("galaxy_mcp.client")
    bad_urls = [
        "ftp://example.org/",
        "http://user:secret@example.org/",
        "http://example.org/?key=1",
    ]
    for url in bad_urls:
        with pytest.raises(ConfigurationError):
            list_galaxy_tools(url)
    with pytest.raises(ConfigurationError):
        list_galaxy_tools(galaxy_server["url"], api_key="bad key with spaces")
    with pytest.raises(ConfigurationError):
        show_galaxy_tool(galaxy_server["url"], tool_id="../etc")
    assert galaxy_server["seen"] == []


def test_galaxy_missing_package_names_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "galaxy_mcp", None)
    monkeypatch.setitem(sys.modules, "galaxy_mcp.client", None)
    with pytest.raises(BackendUnavailableError, match="bio-galaxy"):
        list_galaxy_tools("https://usegalaxy.org")


@LIVE
def test_live_galaxy_public_server_lists_tools() -> None:
    pytest.importorskip("galaxy_mcp.client")
    result = list_galaxy_tools("https://usegalaxy.org", limit=2)
    assert result["url"] == "https://usegalaxy.org/"
    assert result["total_tools"] > 1000
    assert result["returned"] == 2
    cut = show_galaxy_tool("https://usegalaxy.org", tool_id="Cut1")
    assert cut["tool"]["id"] == "Cut1"


# clawbio: the packaged analyze-fasta script on a small local FASTA, checked by hand.

FASTA_TEXT = ">seq1 test\nATGCGCAA\nTT\n>seq2\nGGGCCCAAAT\n"


def test_clawbio_summary_matches_hand_counts(tmp_path: Path) -> None:
    clawbio = pytest.importorskip("clawbio")
    pytest.importorskip("Bio.SeqIO")
    source = tmp_path / "sequences.fa"
    source.write_text(FASTA_TEXT, encoding="utf-8")
    result = analyze_fasta_with_clawbio(source)
    # seq1 = ATGCGCAATT: 10 bp, G+C = 4 (40.0 %), A3 T3 G2 C2 N0.
    # seq2 = GGGCCCAAAT: 10 bp, G+C = 6 (60.0 %), A3 T1 G3 C3 N0.
    assert result["sequences"] == [
        {
            "index": 0,
            "id": "seq1",
            "length_bp": 10,
            "gc_percent": 40.0,
            "composition": {"A": 3, "T": 3, "G": 2, "C": 2, "N": 0},
        },
        {
            "index": 1,
            "id": "seq2",
            "length_bp": 10,
            "gc_percent": 60.0,
            "composition": {"A": 3, "T": 1, "G": 3, "C": 3, "N": 0},
        },
    ]
    assert result["summary"] == {
        "total_sequences": 2,
        "total_residues": 20,
        "min_length": 10,
        "max_length": 10,
        "n50": 10,
        "avg_gc_percent": 50.0,
    }
    assert "analysis_date" not in result and "file_path" not in result

    script = (
        Path(clawbio.__file__).resolve().parent / "skills" / "analyze-fasta" / "analyze_fasta.py"
    )
    spec = importlib.util.spec_from_file_location("direct_analyze_fasta", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    direct = module.analyze_fasta(str(source))
    assert [row["gc_percent"] for row in result["sequences"]] == [
        item["gc_content"] for item in direct["sequences"]
    ]


def test_clawbio_rejects_protein_short_and_oversized_input(tmp_path: Path) -> None:
    pytest.importorskip("clawbio")
    pytest.importorskip("Bio.SeqIO")
    protein = tmp_path / "protein.faa"
    protein.write_text(">p1\nMKTAYIAKQRQISFVKSHFSRQ\n", encoding="utf-8")
    with pytest.raises(InputFormatError, match="protein"):
        analyze_fasta_with_clawbio(protein)
    short = tmp_path / "short.fa"
    short.write_text(">s\nACGT\n", encoding="utf-8")
    with pytest.raises(InputFormatError):
        analyze_fasta_with_clawbio(short)
    with pytest.raises(ConfigurationError, match="max_bytes"):
        analyze_fasta_with_clawbio(protein, max_bytes=5)
    with pytest.raises(ConfigurationError):
        analyze_fasta_with_clawbio(tmp_path / "missing.fa")


def test_clawbio_missing_package_names_extra(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    source = tmp_path / "sequences.fa"
    source.write_text(FASTA_TEXT, encoding="utf-8")
    monkeypatch.setitem(sys.modules, "clawbio", None)
    with pytest.raises(BackendUnavailableError, match="bio-clawbio"):
        analyze_fasta_with_clawbio(source)


@LIVE
def test_live_smartapi_search_and_describe() -> None:
    pytest.importorskip("smartapi_mcp.smartapi")
    result = search_smartapi_registry("mygene", limit=3)
    assert MYGENE_ID in [match["smartapi_id"] for match in result["matches"]]
    spec = describe_smartapi_api(MYGENE_ID)
    assert spec["title"] == "MyGene.info API"
    assert spec["server_url"] == "https://mygene.info/v3"
    assert spec["operation_count"] == 6
