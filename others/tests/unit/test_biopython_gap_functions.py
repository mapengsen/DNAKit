"""Tests for the Biopython network and plot wrappers in dnakit.sequence.biopython_gap_functions.

Remote services are never contacted: each upstream Biopython callable is replaced with a fake
that records its arguments and returns a canned response.
"""

from __future__ import annotations

import io
import socket
import sys
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("Bio")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequence import biopython_gap_functions as bg
from dnakit.sequence.biopython_gap_functions import (
    blast_qblast_search,
    gc_skew_window_plot,
    kegg_convert_identifiers,
    kegg_database_info,
    kegg_find_entries,
    kegg_get_entries,
    kegg_link_entries,
    kegg_list_entries,
    ncbi_efetch,
    ncbi_einfo,
    ncbi_elink,
    ncbi_epost,
    ncbi_esearch,
    ncbi_esummary,
    ncbiwww_qblast_search,
    phylo_draw_ascii_text,
    phylo_draw_newick_png,
    togows_convert_data,
    togows_count_matches,
    togows_fetch_entry,
    togows_iterate_ids,
    togows_search_ids,
)

ESEARCH_XML = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eSearchResult PUBLIC "-//NLM//DTD esearch 20060628//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/esearch.dtd">\n'
    b"<eSearchResult><Count>2</Count><RetMax>2</RetMax><RetStart>0</RetStart>"
    b"<IdList><Id>123</Id><Id>456</Id></IdList><TranslationSet/>"
    b"<QueryTranslation>BRCA1[All Fields]</QueryTranslation></eSearchResult>\n"
)
EINFO_XML = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD einfo 20190110//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20190110/einfo.dtd">\n'
    b"<eInfoResult><DbList><DbName>pubmed</DbName><DbName>nucleotide</DbName></DbList>"
    b"</eInfoResult>\n"
)
EPOST_XML = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE ePostResult PUBLIC "-//NLM//DTD epost 20060628//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20060628/epost.dtd">\n'
    b"<ePostResult><QueryKey>1</QueryKey><WebEnv>MCID_test_env</WebEnv></ePostResult>\n"
)
ESUMMARY_XML = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eSummaryResult PUBLIC "-//NLM//DTD eSummaryResult, 29 October 2004//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/eSummary_041029.dtd">\n'
    b"<eSummaryResult><DocSum><Id>123</Id>"
    b'<Item Name="Title" Type="String">Test record</Item></DocSum></eSummaryResult>\n'
)
ELINK_XML = (
    b'<?xml version="1.0" encoding="UTF-8" ?>\n'
    b'<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD elink 20101123//EN" '
    b'"https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">\n'
    b"<eLinkResult><LinkSet><DbFrom>nucleotide</DbFrom><IdList><Id>123</Id></IdList>"
    b"<LinkSetDb><DbTo>pubmed</DbTo><LinkName>nuccore_pubmed</LinkName>"
    b"<Link><Id>999</Id></Link></LinkSetDb></LinkSet></eLinkResult>\n"
)
FASTA_TEXT = ">seq1 test\nACGTACGT\n"

# A minimal NCBI BLAST XML document with one hit and one HSP for a 8 nt DNA query.
BLAST_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"
  "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_reference>Test reference.</BlastOutput_reference>
  <BlastOutput_version>BLASTN 2.14.0+</BlastOutput_version>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query_dna</BlastOutput_query-def>
  <BlastOutput_query-len>8</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_gap-open>5</Parameters_gap-open>
      <Parameters_gap-extend>2</Parameters_gap-extend>
      <Parameters_filter>L;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>query_dna</Iteration_query-def>
      <Iteration_query-len>8</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|1|ref|NM_000001.1|</Hit_id>
          <Hit_def>Homo sapiens test gene</Hit_def>
          <Hit_accession>NM_000001</Hit_accession>
          <Hit_len>20</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>14.6</Hsp_bit-score>
              <Hsp_score>8</Hsp_score>
              <Hsp_evalue>0.01</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>8</Hsp_query-to>
              <Hsp_hit-from>5</Hsp_hit-from>
              <Hsp_hit-to>12</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>8</Hsp_identity>
              <Hsp_positive>8</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>8</Hsp_align-len>
              <Hsp_qseq>ACGTACGT</Hsp_qseq>
              <Hsp_hseq>ACGTACGT</Hsp_hseq>
              <Hsp_midline>||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1</Statistics_db-num>
          <Statistics_db-len>20</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>1.37</Statistics_lambda>
          <Statistics_entropy>1.0</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


class _Recorder:
    """Fake upstream callable: records the keyword arguments and returns a canned handle."""

    def __init__(self, payload: bytes | str, *, seen_email: list[Any] | None = None, entrez=None):
        self.payload = payload
        self.calls: list[dict[str, Any]] = []
        self.seen_email = seen_email
        self.entrez = entrez

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append({"args": args, "kwargs": kwargs})
        if self.entrez is not None and self.seen_email is not None:
            self.seen_email.append(self.entrez.email)
        if isinstance(self.payload, str):
            return io.StringIO(self.payload)
        return io.BytesIO(self.payload)


def _entrez(monkeypatch, name: str, payload: bytes) -> _Recorder:
    from Bio import Entrez

    fake = _Recorder(payload)
    monkeypatch.setattr(Entrez, name, fake)
    return fake


# --- NCBI Entrez ---------------------------------------------------------------


def test_ncbi_esearch_passes_arguments_and_parses_ids(monkeypatch):
    fake = _entrez(monkeypatch, "esearch", ESEARCH_XML)

    result = ncbi_esearch("nucleotide", "BRCA1", retmax=5)

    assert fake.calls == [
        {"args": (), "kwargs": {"db": "nucleotide", "term": "BRCA1", "retmax": 5}}
    ]
    assert result["status"] == "completed"
    assert result["service"].startswith("NCBI Entrez")
    assert result["upstream"] == "Bio.Entrez.esearch"
    assert result["parameters"] == {"db": "nucleotide", "term": "BRCA1", "retmax": 5}
    assert result["output_path"] is None
    assert result["retrieved_utc"].endswith("Z")
    assert result["biopython_version"]
    assert result["data"]["count"] == 2
    assert result["data"]["ids"] == ["123", "456"]
    assert result["data"]["QueryTranslation"] == "BRCA1[All Fields]"


def test_ncbi_esearch_default_retmax_is_twenty(monkeypatch):
    fake = _entrez(monkeypatch, "esearch", ESEARCH_XML)

    ncbi_esearch("pubmed", "dna repair")

    assert fake.calls[0]["kwargs"]["retmax"] == 20


def test_ncbi_efetch_returns_text_and_writes_new_file(monkeypatch, tmp_path: Path):
    fake = _entrez(monkeypatch, "efetch", FASTA_TEXT.encode())
    target = tmp_path / "records.fasta"

    result = ncbi_efetch("nucleotide", ["123", "456"], output_path=str(target))

    assert fake.calls[0]["kwargs"] == {
        "db": "nucleotide",
        "id": "123,456",
        "rettype": "fasta",
        "retmode": "text",
    }
    assert target.read_text(encoding="utf-8") == FASTA_TEXT
    assert result["output_path"] == str(target.resolve())
    assert result["data"]["text"] == FASTA_TEXT
    assert result["data"]["characters"] == len(FASTA_TEXT)


def test_ncbi_efetch_without_output_path_returns_text_only(monkeypatch):
    _entrez(monkeypatch, "efetch", FASTA_TEXT.encode())

    result = ncbi_efetch("nucleotide", ["123"])

    assert result["output_path"] is None
    assert result["data"]["text"] == FASTA_TEXT


def test_ncbi_efetch_refuses_existing_output(monkeypatch, tmp_path: Path):
    _entrez(monkeypatch, "efetch", FASTA_TEXT.encode())
    target = tmp_path / "exists.fasta"
    target.write_text("keep me", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="must not exist yet"):
        ncbi_efetch("nucleotide", ["123"], output_path=str(target))
    assert target.read_text(encoding="utf-8") == "keep me"


def test_ncbi_einfo_lists_databases_and_describes_one(monkeypatch):
    fake = _entrez(monkeypatch, "einfo", EINFO_XML)

    listing = ncbi_einfo()
    described = ncbi_einfo(db="nucleotide")

    assert fake.calls[0]["kwargs"] == {}
    assert fake.calls[1]["kwargs"] == {"db": "nucleotide"}
    assert listing["data"]["DbList"] == ["pubmed", "nucleotide"]
    assert described["parameters"] == {"db": "nucleotide"}


def test_ncbi_esummary_elink_epost_pass_identifiers(monkeypatch):
    summary = _entrez(monkeypatch, "esummary", ESUMMARY_XML)
    link = _entrez(monkeypatch, "elink", ELINK_XML)
    post = _entrez(monkeypatch, "epost", EPOST_XML)

    summary_result = ncbi_esummary("nucleotide", ["123", "456"])
    link_result = ncbi_elink("nucleotide", "pubmed", ["123"], linkname="nuccore_pubmed")
    post_result = ncbi_epost("nucleotide", ["123", "456"])

    assert summary.calls[0]["kwargs"] == {"db": "nucleotide", "id": "123,456"}
    assert link.calls[0]["kwargs"] == {
        "dbfrom": "nucleotide",
        "db": "pubmed",
        "id": "123",
        "linkname": "nuccore_pubmed",
    }
    assert post.calls[0]["kwargs"] == {"db": "nucleotide", "id": "123,456"}
    assert summary_result["data"][0]["Id"] == "123"
    assert summary_result["data"][0]["Title"] == "Test record"
    linkset = link_result["data"][0]
    assert linkset["DbFrom"] == "nucleotide"
    assert linkset["LinkSetDb"][0]["DbTo"] == "pubmed"
    assert linkset["LinkSetDb"][0]["Link"][0]["Id"] == "999"
    assert post_result["data"] == {"QueryKey": "1", "WebEnv": "MCID_test_env"}
    assert post_result["upstream"] == "Bio.Entrez.epost"


def test_entrez_email_is_set_only_during_the_call_and_not_recorded(monkeypatch):
    from Bio import Entrez

    seen: list[Any] = []
    fake = _Recorder(ESEARCH_XML, seen_email=seen, entrez=Entrez)
    monkeypatch.setattr(Entrez, "esearch", fake)
    monkeypatch.setattr(Entrez, "email", None)

    result = ncbi_esearch("pubmed", "x", email="reviewer@example.org")

    assert seen == ["reviewer@example.org"]
    assert Entrez.email is None
    assert "reviewer@example.org" not in repr(result["parameters"])
    assert "email" not in result["parameters"]


def test_entrez_without_email_leaves_module_setting_alone(monkeypatch):
    from Bio import Entrez

    monkeypatch.setattr(Entrez, "email", "previous@example.org")
    seen: list[Any] = []
    monkeypatch.setattr(Entrez, "esearch", _Recorder(ESEARCH_XML, seen_email=seen, entrez=Entrez))

    ncbi_esearch("pubmed", "x")

    assert seen == ["previous@example.org"]
    assert Entrez.email == "previous@example.org"


def test_entrez_rejects_bad_input_before_any_call(monkeypatch):
    fake = _entrez(monkeypatch, "esearch", ESEARCH_XML)

    with pytest.raises(ConfigurationError):
        ncbi_esearch("pub med", "x")
    with pytest.raises(ConfigurationError):
        ncbi_esearch("pubmed", "")
    with pytest.raises(ConfigurationError):
        ncbi_esearch("pubmed", "x", retmax=True)
    with pytest.raises(ConfigurationError):
        ncbi_esearch("pubmed", "x", email="not-an-address")
    with pytest.raises(ConfigurationError):
        ncbi_esummary("nucleotide", [])
    with pytest.raises(ConfigurationError):
        ncbi_esummary("nucleotide", ["1,2"])
    with pytest.raises(ConfigurationError):
        ncbi_esummary("nucleotide", [str(i) for i in range(bg._MAX_IDS + 1)])
    assert fake.calls == []


def test_entrez_verbose_must_be_boolean(monkeypatch):
    _entrez(monkeypatch, "esearch", ESEARCH_XML)

    with pytest.raises(ConfigurationError, match="verbose must be boolean"):
        ncbi_esearch("pubmed", "x", verbose="yes")


def test_missing_biopython_reports_install_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "Bio.Entrez", None)

    with pytest.raises(BackendUnavailableError) as raised:
        ncbi_esearch("pubmed", "x")
    assert "bio-formats" in (raised.value.hint or "")


# --- NCBI BLAST QBLAST ----------------------------------------------------------------------------


def _blast_fake(monkeypatch, module_name: str) -> _Recorder:
    module = __import__(module_name, fromlist=["qblast"])
    fake = _Recorder(BLAST_XML)
    monkeypatch.setattr(module, "qblast", fake)
    return fake


@pytest.mark.parametrize(
    ("function", "module_name", "upstream"),
    [
        (ncbiwww_qblast_search, "Bio.Blast.NCBIWWW", "Bio.Blast.NCBIWWW.qblast"),
        (blast_qblast_search, "Bio.Blast", "Bio.Blast.qblast"),
    ],
)
def test_qblast_wrappers_pass_arguments_and_summarise_hits(
    monkeypatch, function, module_name, upstream
):
    fake = _blast_fake(monkeypatch, module_name)

    result = function(
        "blastn",
        "nt",
        "ACGTACGT",
        hitlist_size=10,
        expect=0.5,
        entrez_query="Homo sapiens[organism]",
    )

    assert len(fake.calls) == 1
    kwargs = fake.calls[0]["kwargs"]
    assert kwargs == {
        "program": "blastn",
        "database": "nt",
        "sequence": "ACGTACGT",
        "hitlist_size": 10,
        "expect": 0.5,
        "entrez_query": "Homo sapiens[organism]",
        "format_type": "XML",
    }
    assert result["upstream"] == upstream
    assert result["parameters"]["sequence_length"] == 8
    assert "sequence" not in result["parameters"]
    data = result["data"]
    assert data["query_id"] == "Query_1"
    assert data["query"] == "query_dna"
    assert data["database"] == "nt"
    assert data["alignment_count"] == 1
    alignment = data["alignments"][0]
    assert alignment["hit_def"] == "Homo sapiens test gene"
    assert alignment["length"] == 20
    hsp = alignment["hsps"][0]
    assert hsp["score"] == 8
    assert hsp["identities"] == 8
    assert hsp["align_length"] == 8
    assert (hsp["query_start"], hsp["query_end"]) == (1, 8)
    assert (hsp["subject_start"], hsp["subject_end"]) == (5, 12)
    assert hsp["expect"] == pytest.approx(0.01)


def test_qblast_writes_raw_xml_to_new_file(monkeypatch, tmp_path: Path):
    _blast_fake(monkeypatch, "Bio.Blast.NCBIWWW")
    target = tmp_path / "hits.xml"

    result = ncbiwww_qblast_search("blastn", "nt", "ACGTACGT", output_xml=str(target))

    assert result["output_path"] == str(target.resolve())
    assert target.read_text(encoding="utf-8").startswith("<?xml")


def test_qblast_email_is_set_on_the_upstream_module_only_while_running(monkeypatch):
    from Bio.Blast import NCBIWWW

    seen: list[Any] = []

    class _EmailRecorder(_Recorder):
        def __call__(self, *args: Any, **kwargs: Any) -> Any:
            seen.append(NCBIWWW.email)
            return super().__call__(*args, **kwargs)

    monkeypatch.setattr(NCBIWWW, "qblast", _EmailRecorder(BLAST_XML))
    monkeypatch.setattr(NCBIWWW, "email", None)

    result = ncbiwww_qblast_search("blastn", "nt", "ACGTACGT", email="reviewer@example.org")

    assert seen == ["reviewer@example.org"]
    assert NCBIWWW.email is None
    assert "email" not in result["parameters"]


@pytest.mark.parametrize("function", [ncbiwww_qblast_search, blast_qblast_search])
def test_qblast_rejects_unknown_program_and_bad_sequence(monkeypatch, function):
    fake = _blast_fake(
        monkeypatch, "Bio.Blast.NCBIWWW" if function is ncbiwww_qblast_search else "Bio.Blast"
    )

    with pytest.raises(ConfigurationError, match="program must be one of"):
        function("psiblast", "nt", "ACGT")
    with pytest.raises(ConfigurationError):
        function("blastn", "nt", "   ")
    with pytest.raises(ConfigurationError):
        function("blastn", "nt", "ACGT", hitlist_size=0)
    with pytest.raises(ConfigurationError):
        function("blastn", "nt", "ACGT", expect=0)
    assert fake.calls == []


def test_qblast_refuses_existing_xml_output(monkeypatch, tmp_path: Path):
    fake = _blast_fake(monkeypatch, "Bio.Blast.NCBIWWW")
    target = tmp_path / "hits.xml"
    target.write_text("old", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="must not exist yet"):
        ncbiwww_qblast_search("blastn", "nt", "ACGT", output_xml=str(target))
    assert fake.calls == []


# --- KEGG REST ------------------------------------------------------------------------------------


def _kegg(monkeypatch, name: str, payload: str) -> _Recorder:
    from Bio.KEGG import REST

    fake = _Recorder(payload)
    monkeypatch.setattr(REST, name, fake)
    return fake


def test_kegg_database_info_and_list(monkeypatch):
    info = _kegg(monkeypatch, "kegg_info", "pathway          KEGG Pathway Database\nRelease 1\n")
    listing = _kegg(
        monkeypatch,
        "kegg_list",
        "path:hsa00010\tGlycolysis / Gluconeogenesis\npath:hsa00020\tCitrate cycle (TCA cycle)\n",
    )

    described = kegg_database_info("pathway")
    organism = kegg_list_entries("pathway", organism="hsa")
    entries = kegg_list_entries(["hsa:10458", "ece:Z5100"])

    assert info.calls[0]["args"] == ("pathway",)
    assert described["upstream"] == "Bio.KEGG.REST.kegg_info"
    assert described["service"].startswith("KEGG REST")
    assert described["data"]["line_count"] == 2
    assert listing.calls[0]["args"] == ("pathway", "hsa")
    assert organism["data"]["rows"][0] == ["path:hsa00010", "Glycolysis / Gluconeogenesis"]
    assert organism["data"]["row_count"] == 2
    assert listing.calls[1]["args"] == (["hsa:10458", "ece:Z5100"], None)
    assert entries["parameters"]["database"] == ["hsa:10458", "ece:Z5100"]


def test_kegg_find_convert_and_link_pass_options(monkeypatch):
    find = _kegg(monkeypatch, "kegg_find", "hsa:10458\tDNA repair gene\n")
    conv = _kegg(monkeypatch, "kegg_conv", "hsa:10458\tncbi-geneid:10458\n")
    link = _kegg(monkeypatch, "kegg_link", "hsa:10458\tpath:hsa00010\n")

    found = kegg_find_entries("genes", "DNA repair", option=None)
    converted = kegg_convert_identifiers("ncbi-geneid", ["hsa:10458"])
    linked = kegg_link_entries("pathway", "hsa:10458", option=None)

    assert find.calls[0]["args"] == ("genes", "DNA repair", None)
    assert found["data"]["rows"] == [["hsa:10458", "DNA repair gene"]]
    assert conv.calls[0]["args"] == ("ncbi-geneid", ["hsa:10458"], None)
    assert converted["data"]["rows"] == [["hsa:10458", "ncbi-geneid:10458"]]
    assert converted["upstream"] == "Bio.KEGG.REST.kegg_conv"
    assert link.calls[0]["args"] == ("pathway", "hsa:10458", None)
    assert linked["data"]["row_count"] == 1
    assert linked["upstream"] == "Bio.KEGG.REST.kegg_link"


def test_kegg_get_passes_option_and_bounds_entry_count(monkeypatch):
    fake = _kegg(monkeypatch, "kegg_get", ">hsa:10458\nATGCATGC\n")

    result = kegg_get_entries(["hsa:10458", "ece:Z5100"], option="ntseq")
    single = kegg_get_entries("hsa:10458+ece:Z5100")

    assert fake.calls[0]["args"] == (["hsa:10458", "ece:Z5100"], "ntseq")
    assert result["data"]["text"] == ">hsa:10458\nATGCATGC\n"
    assert result["parameters"] == {"entries": ["hsa:10458", "ece:Z5100"], "option": "ntseq"}
    assert fake.calls[1]["args"] == ("hsa:10458+ece:Z5100", None)
    assert single["parameters"]["entries"] == "hsa:10458+ece:Z5100"


def test_kegg_rejects_bad_arguments_before_calling_upstream(monkeypatch):
    fake = _kegg(monkeypatch, "kegg_get", "")

    with pytest.raises(ConfigurationError, match="1 to 10"):
        kegg_get_entries([f"hsa:{i}" for i in range(11)])
    with pytest.raises(ConfigurationError, match="image is not supported"):
        kegg_get_entries("hsa:10458", option="image")
    with pytest.raises(ConfigurationError):
        kegg_database_info("path way")
    with pytest.raises(ConfigurationError):
        kegg_find_entries("genes", "")
    with pytest.raises(ConfigurationError):
        kegg_find_entries("genes", "x", option="image")
    with pytest.raises(ConfigurationError):
        kegg_convert_identifiers("ncbi-geneid", ["hsa:10458"], option="json")
    assert fake.calls == []


# --- TogoWS ---------------------------------------------------------------------------------------


def _togows(monkeypatch, name: str, payload: Any) -> _Recorder:
    from Bio import TogoWS

    fake = _Recorder(payload) if isinstance(payload, str) else payload
    monkeypatch.setattr(TogoWS, name, fake)
    return fake


def test_togows_fetch_entry_passes_format_and_field(monkeypatch):
    fake = _togows(monkeypatch, "entry", '{"id": 7157}')

    result = togows_fetch_entry("ncbi-geneid", "7157", output_format="json")
    plain = togows_fetch_entry("ncbi-geneid", "7157", field="symbol")

    assert fake.calls[0]["kwargs"] == {
        "db": "ncbi-geneid",
        "id": "7157",
        "format": "json",
        "field": None,
    }
    assert result["upstream"] == "Bio.TogoWS.entry"
    assert result["parameters"] == {
        "database": "ncbi-geneid",
        "identifier": "7157",
        "output_format": "json",
    }
    assert result["data"]["text"] == '{"id": 7157}'
    assert fake.calls[1]["kwargs"]["field"] == "symbol"
    assert plain["parameters"] == {
        "database": "ncbi-geneid",
        "identifier": "7157",
        "field": "symbol",
    }


def test_togows_search_ids_returns_identifiers_per_page(monkeypatch):
    fake = _togows(monkeypatch, "search", "101\n102\n103\n")

    result = togows_search_ids("pubmed", "diabetes+human", offset=1, limit=3)

    assert fake.calls[0]["kwargs"] == {
        "db": "pubmed",
        "query": "diabetes+human",
        "offset": 1,
        "limit": 3,
        "format": None,
    }
    assert result["upstream"] == "Bio.TogoWS.search"
    assert result["data"] == {"ids": ["101", "102", "103"], "count": 3}


def test_togows_search_ids_in_other_format_returns_text(monkeypatch):
    fake = _togows(monkeypatch, "search", '{"hits": []}')

    result = togows_search_ids("pubmed", "x", output_format="json")

    assert fake.calls[0]["kwargs"]["format"] == "json"
    assert result["data"] == {"text": '{"hits": []}'}


def test_togows_count_matches_returns_integer(monkeypatch):
    calls: list[dict[str, Any]] = []

    def _count(**kwargs):
        calls.append(kwargs)
        return 42

    _togows(monkeypatch, "search_count", _count)

    result = togows_count_matches("pubmed", "diabetes")

    assert calls == [{"db": "pubmed", "query": "diabetes"}]
    assert result["upstream"] == "Bio.TogoWS.search_count"
    assert result["data"] == {"count": 42}


def test_togows_iterate_ids_is_bounded(monkeypatch):
    calls: list[dict[str, Any]] = []

    def _iterate(**kwargs):
        calls.append(kwargs)
        yield from ["1", "2", "3", "4", "5"]

    _togows(monkeypatch, "search_iter", _iterate)

    result = togows_iterate_ids("pubmed", "diabetes", limit=3)

    assert calls == [{"db": "pubmed", "query": "diabetes", "limit": 3}]
    assert result["upstream"] == "Bio.TogoWS.search_iter"
    assert result["data"] == {"ids": ["1", "2", "3"], "returned_count": 3}
    with pytest.raises(ConfigurationError):
        togows_iterate_ids("pubmed", "diabetes", limit=0)
    with pytest.raises(ConfigurationError):
        togows_iterate_ids("pubmed", "diabetes", limit=10_000)


def test_togows_convert_data_passes_formats(monkeypatch):
    fake = _togows(monkeypatch, "convert", ">X\nACGT\n")

    result = togows_convert_data("ID   X; SV 1.\nSQ   ACGT\n//\n", "genbank", "fasta")

    assert fake.calls[0]["args"][1:] == ("genbank", "fasta")
    assert result["upstream"] == "Bio.TogoWS.convert"
    assert result["data"]["text"] == ">X\nACGT\n"
    assert result["parameters"]["in_format"] == "genbank"
    assert result["parameters"]["out_format"] == "fasta"


def test_togows_rejects_bad_identifiers(monkeypatch):
    fake = _togows(monkeypatch, "entry", "x")

    with pytest.raises(ConfigurationError):
        togows_fetch_entry("ncbi-geneid", "7157 9")
    with pytest.raises(ConfigurationError):
        togows_fetch_entry("ncbi geneid", "7157")
    with pytest.raises(ConfigurationError):
        togows_search_ids("pubmed", "x", offset=0)
    assert fake.calls == []


# --- Phylo trees and GC skew ----------------------------------------------------------------------

TREE = "((A:0.1,B:0.2):0.3,C:0.5);"


def test_phylo_draw_newick_png_writes_image_and_tree_data(tmp_path: Path):
    target = tmp_path / "tree.png"

    result = phylo_draw_newick_png(TREE, str(target))

    assert result["status"] == "completed"
    assert result["format"] == "png"
    assert result["image_path"] == str(target.resolve())
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    data = result["data"]
    assert data["tip_names"] == ["A", "B", "C"]
    assert data["terminal_count"] == 3
    assert data["clade_count"] == 5
    assert data["total_branch_length"] == pytest.approx(1.1)
    # Bio.Phylo.write adds a root branch length, so the Newick text is Biopython's normalised form.
    assert data["newick"] == "((A:0.1,B:0.2):0.3,C:0.5):0;"


def test_phylo_draw_newick_png_refuses_existing_or_non_png(tmp_path: Path):
    target = tmp_path / "tree.png"
    target.write_bytes(b"old")

    with pytest.raises(ConfigurationError, match="must not exist yet"):
        phylo_draw_newick_png(TREE, str(target))
    with pytest.raises(ConfigurationError, match=r"must end with \.png"):
        phylo_draw_newick_png(TREE, str(tmp_path / "tree.svg"))
    assert target.read_bytes() == b"old"


def test_phylo_draw_ascii_text_writes_text_in_tip_order(tmp_path: Path):
    target = tmp_path / "tree.txt"

    result = phylo_draw_ascii_text(TREE, str(target), column_width=60)

    assert result["status"] == "completed"
    assert result["output_path"] == str(target.resolve())
    assert "image_path" not in result
    text = target.read_text(encoding="utf-8")
    assert text == result["data"]["text"]
    assert result["data"]["tip_names"] == ["A", "B", "C"]
    assert result["data"]["column_width"] == 60
    lines = text.splitlines()
    rows = [next(i for i, line in enumerate(lines) if tip in line) for tip in ("A", "B", "C")]
    assert rows == sorted(rows)
    assert len(set(rows)) == 3


def test_phylo_draw_ascii_text_refuses_existing_output(tmp_path: Path):
    target = tmp_path / "tree.txt"
    target.write_text("old", encoding="utf-8")

    with pytest.raises(ConfigurationError, match="must not exist yet"):
        phylo_draw_ascii_text(TREE, str(target))
    assert target.read_text(encoding="utf-8") == "old"


def test_gc_skew_window_plot_values_on_known_input(tmp_path: Path):
    target = tmp_path / "skew.png"

    result = gc_skew_window_plot("GGGGCCCCAAAATTTT", str(target), window=4)

    assert result["status"] == "completed"
    assert result["image_path"] == str(target.resolve())
    assert target.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    data = result["data"]
    assert data["sequence_length"] == 16
    assert data["window"] == 4
    assert data["window_count"] == 4
    assert data["window_starts"] == [0, 4, 8, 12]
    assert data["gc_skew"] == pytest.approx([1.0, -1.0, 0.0, 0.0])
    assert data["accumulated_gc_skew"] == pytest.approx([1.0, 0.0, 0.0, 0.0])


def test_gc_skew_window_plot_lowercase_and_single_window(tmp_path: Path):
    lower = gc_skew_window_plot("ggggcccc", str(tmp_path / "lower.png"), window=4)
    single = gc_skew_window_plot("GGCC", str(tmp_path / "single.png"))

    assert lower["data"]["gc_skew"] == pytest.approx([1.0, -1.0])
    assert single["data"]["window_count"] == 1
    assert single["data"]["gc_skew"] == [0.0]


def test_gc_skew_window_plot_rejects_bad_input(tmp_path: Path):
    with pytest.raises(ConfigurationError, match="IUPAC DNA"):
        gc_skew_window_plot("ACGTX", str(tmp_path / "bad.png"))
    with pytest.raises(ConfigurationError):
        gc_skew_window_plot("ACGT", str(tmp_path / "bad.png"), window=0)
    with pytest.raises(ConfigurationError):
        gc_skew_window_plot("ACGT", str(tmp_path / "bad.png"), dpi=10)
    assert not (tmp_path / "bad.png").exists()


# --- public surface -------------------------------------------------------------------------------

EXPECTED_PUBLIC = {
    "ncbi_einfo",
    "ncbi_esearch",
    "ncbi_efetch",
    "ncbi_esummary",
    "ncbi_elink",
    "ncbi_epost",
    "ncbiwww_qblast_search",
    "blast_qblast_search",
    "kegg_database_info",
    "kegg_list_entries",
    "kegg_find_entries",
    "kegg_convert_identifiers",
    "kegg_link_entries",
    "kegg_get_entries",
    "togows_fetch_entry",
    "togows_search_ids",
    "togows_count_matches",
    "togows_iterate_ids",
    "togows_convert_data",
    "phylo_draw_newick_png",
    "phylo_draw_ascii_text",
    "gc_skew_window_plot",
}


def test_public_names_match_the_rows_and_each_takes_verbose_last():
    import inspect

    assert set(bg.__all__) == EXPECTED_PUBLIC
    for name in bg.__all__:
        parameters = list(inspect.signature(getattr(bg, name)).parameters.values())
        assert parameters[-1].name == "verbose", name
        assert parameters[-1].kind is inspect.Parameter.KEYWORD_ONLY, name


def test_upstream_qualnames_exist_in_installed_biopython():
    import importlib

    qualnames = [
        "Bio.Blast.NCBIWWW.qblast",
        "Bio.Blast.qblast",
        "Bio.Entrez.einfo",
        "Bio.Entrez.esearch",
        "Bio.Entrez.efetch",
        "Bio.Entrez.esummary",
        "Bio.Entrez.elink",
        "Bio.Entrez.epost",
        "Bio.KEGG.REST.kegg_info",
        "Bio.KEGG.REST.kegg_list",
        "Bio.KEGG.REST.kegg_find",
        "Bio.KEGG.REST.kegg_conv",
        "Bio.KEGG.REST.kegg_link",
        "Bio.KEGG.REST.kegg_get",
        "Bio.TogoWS.entry",
        "Bio.TogoWS.search",
        "Bio.TogoWS.search_count",
        "Bio.TogoWS.search_iter",
        "Bio.TogoWS.convert",
        "Bio.Phylo.draw",
        "Bio.Phylo.draw_ascii",
        "Bio.SeqUtils.xGC_skew",
        "Bio.SeqUtils.GC_skew",
    ]
    for qualname in qualnames:
        module_name, _, attribute = qualname.rpartition(".")
        module = importlib.import_module(module_name)
        assert callable(getattr(module, attribute)), qualname
