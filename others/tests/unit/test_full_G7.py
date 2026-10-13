"""Unit tests for the G7 database wrappers: gget, mygene, myvariant and biothings_client.

Default runs never touch the network. Network-only functions replace the upstream client
call or the HTTP method of the upstream client; the real upstream code then runs on the
canned payload. Set DNAKIT_LIVE_TESTS=1 to run one optional live check per package.
"""

from __future__ import annotations

import inspect
import os
from pathlib import Path

import pytest

from dnakit.databases import biothings_client_functions as biothings_fns
from dnakit.databases import gget_functions as gget_fns
from dnakit.databases import mygene_functions as mygene_fns
from dnakit.databases import myvariant_functions as myvariant_fns
from dnakit.exceptions import BackendExecutionError, ConfigurationError

LIVE = os.environ.get("DNAKIT_LIVE_TESTS") == "1"
live_only = pytest.mark.skipif(not LIVE, reason="set DNAKIT_LIVE_TESTS=1 for live network checks")

COSMIC_CANCER_HEADER = (
    "GENE_NAME\tACCESSION_NUMBER\tLEGACY_MUTATION_ID\tMutation CDS\tMutation AA\t"
    "GENOMIC_MUTATION_ID\tPrimary site\n"
)
COSMIC_CANCER_ROWS = [
    "EGFR\tNM_005228\tCOSM6223\tc.2573T>G\tp.L858R\tCOSV55365\tlung\n",
    "BRAF\tNM_004333\tCOSM476\tc.1799T>A\tp.V600E\tCOSV55410\tskin\n",
    "EGFR\tNM_005228\tCOSM6224\tc.2236_2250del\tp.E746_A750del\tCOSV55366\tlung\n",
]


@pytest.fixture
def gget_mod():
    return pytest.importorskip("gget")


@pytest.fixture
def mygene_mod():
    return pytest.importorskip("mygene")


@pytest.fixture
def myvariant_mod():
    return pytest.importorskip("myvariant")


@pytest.fixture
def biothings_mod():
    return pytest.importorskip("biothings_client")


def _recorder(monkeypatch, owner, name, payload):
    """Replace one upstream callable; return the list that records every call."""
    calls = []

    def fake(*args, **kwargs):
        calls.append((args, kwargs))
        return payload

    monkeypatch.setattr(owner, name, fake)
    return calls


def _cosmic_fixture(tmp_path: Path, name: str = "CancerMutationCensus_AllData_mini.tsv") -> Path:
    path = tmp_path / name
    path.write_text(COSMIC_CANCER_HEADER + "".join(COSMIC_CANCER_ROWS), encoding="utf-8")
    return path


# --- gget: offline local COSMIC -------------------------------------------------------


def test_gget_cosmic_local_matches_direct_upstream_call(gget_mod, tmp_path):
    path = _cosmic_fixture(tmp_path)
    result = gget_fns.gget_cosmic_local("EGFR", path, verbose=False)
    direct = gget_mod.cosmic(
        "EGFR", cosmic_tsv_path=str(path), cosmic_project="cancer", json=True, verbose=False
    )
    assert result["status"] == "complete"
    assert result["operation"] == "gget_cosmic_local"
    assert result["cosmic_project"] == "cancer"
    assert result["count"] == 2
    assert result["records"] == direct
    # Hand check: the two EGFR rows, in file order, with spaces in column names replaced.
    assert [row["Mutation_AA"] for row in result["records"]] == ["p.L858R", "p.E746_A750del"]
    assert result["provenance"]["network"] is False
    assert result["provenance"]["package"] == "gget"


def test_gget_cosmic_local_no_match_is_empty_not_error(gget_mod, tmp_path):
    path = _cosmic_fixture(tmp_path)
    result = gget_fns.gget_cosmic_local("KRAS", path, verbose=False)
    assert result["count"] == 0
    assert result["records"] == []


def test_gget_cosmic_local_census_columns_and_explicit_project(gget_mod, tmp_path):
    path = tmp_path / "census_subset.tsv"
    path.write_text(
        "GENE_SYMBOL\tTRANSCRIPT_ACCESSION\tCOSMIC_GENE_ID\tCOSMIC_SAMPLE_ID\t"
        "COSMIC_PHENOTYPE_ID\tGENOMIC_MUTATION_ID\tLEGACY_MUTATION_ID\tSAMPLE_NAME\t"
        "MUTATION_CDS\tMUTATION_AA\tMUTATION_ID\tCOSMIC_STUDY_ID\n"
        "TP53\tENST00000269305\t1\t10\t20\tCOSV1\tCOSM1\tS1\tc.817C>T\tp.R273C\tM1\t5\n",
        encoding="utf-8",
    )
    result = gget_fns.gget_cosmic_local("TP53", path, cosmic_project="census", verbose=False)
    assert result["count"] == 1
    assert result["records"][0]["MUTATION_AA"] == "p.R273C"
    assert result["cosmic_project"] == "census"


def test_gget_cosmic_local_validates_inputs(gget_mod, tmp_path):
    path = _cosmic_fixture(tmp_path)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_cosmic_local("   ", path, verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_cosmic_local("EGFR", tmp_path / "missing.tsv", verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_cosmic_local("EGFR", path, cosmic_project="unknown", verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_cosmic_local("EGFR", path, limit=0, verbose=False)


# --- gget: network-only, upstream replaced ----------------------------------------------


def test_gget_blast_passes_validated_arguments(gget_mod, monkeypatch):
    calls = _recorder(
        monkeypatch,
        gget_mod,
        "blast",
        [{"Accession": "NM_000546", "Percent Identity": 99.5, "E-Value": 0.0}],
    )
    result = gget_fns.gget_blast("acgtacgt", program="blastn", limit=5, expect=0.01, verbose=False)
    args, kwargs = calls[0]
    assert args == ("ACGTACGT",)
    assert kwargs["program"] == "blastn" and kwargs["limit"] == 5 and kwargs["expect"] == 0.01
    assert kwargs["json"] is True
    assert result["query_length"] == 8
    assert result["count"] == 1
    assert result["hits"][0]["Accession"] == "NM_000546"
    with pytest.raises(ConfigurationError):
        gget_fns.gget_blast("ACGT!", verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_blast("ACGT", program="blastz", verbose=False)


def test_gget_blat_returns_hits_and_checks_seqtype(gget_mod, monkeypatch):
    calls = _recorder(monkeypatch, gget_mod, "blat", [{"chrom": "chr7", "start": 55019016}])
    result = gget_fns.gget_blat("ACGTTGCA", seqtype="DNA", assembly="human", verbose=False)
    assert calls[0][0] == ("ACGTTGCA",)
    assert calls[0][1]["seqtype"] == "DNA"
    assert result["count"] == 1 and result["hits"][0]["chrom"] == "chr7"
    with pytest.raises(ConfigurationError):
        gget_fns.gget_blat("ACGT", seqtype="RNA", verbose=False)


def test_gget_opentargets_checks_ensembl_gene_and_resource(gget_mod, monkeypatch):
    calls = _recorder(monkeypatch, gget_mod, "opentargets", [{"drug.name": "BELVARAFENIB"}])
    result = gget_fns.gget_opentargets("ENSG00000157764", resource="drugs", limit=2, verbose=False)
    args, kwargs = calls[0]
    assert args == ("ENSG00000157764",)
    assert kwargs["resource"] == "drugs" and kwargs["limit"] == 2
    assert result["records"] == [{"drug.name": "BELVARAFENIB"}]
    with pytest.raises(ConfigurationError):
        gget_fns.gget_opentargets("BRAF", verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_opentargets("ENSG00000157764", resource="proteins", verbose=False)


def test_gget_archs4_checks_which_species_and_ensembl_flag(gget_mod, monkeypatch):
    calls = _recorder(monkeypatch, gget_mod, "archs4", [{"gene": "SOX2", "correlation": 0.9}])
    result = gget_fns.gget_archs4("SOX2", which="correlation", gene_count=10, verbose=False)
    assert calls[0][0] == ("SOX2",)
    assert calls[0][1]["gene_count"] == 10 and calls[0][1]["ensembl"] is False
    assert result["count"] == 1 and result["species"] == "human"
    with pytest.raises(ConfigurationError):
        gget_fns.gget_archs4("SOX2", species="zebrafish", verbose=False)
    with pytest.raises(ConfigurationError):
        gget_fns.gget_archs4("SOX2", ensembl=True, verbose=False)


def test_gget_bgee_orthologs_need_one_gene(gget_mod, monkeypatch):
    calls = _recorder(monkeypatch, gget_mod, "bgee", [{"geneId": "ENSMUSG00000020122"}])
    result = gget_fns.gget_bgee("ENSG00000139618", record_type="orthologs", verbose=False)
    assert calls[0][0] == ("ENSG00000139618",)
    assert calls[0][1]["type"] == "orthologs"
    assert result["gene_ids"] == ["ENSG00000139618"]
    with pytest.raises(ConfigurationError):
        gget_fns.gget_bgee(["ENSG00000139618", "ENSG00000157764"], verbose=False)


def test_gget_enrichr_sends_sorted_gene_list_and_library(gget_mod, monkeypatch):
    calls = _recorder(monkeypatch, gget_mod, "enrichr", [{"Term": "p53 pathway", "P-value": 0.01}])
    genes = ["TP53", "MDM2", "CDKN1A"]
    result = gget_fns.gget_enrichr(genes, database="KEGG_2021_Human", verbose=False)
    args, kwargs = calls[0]
    assert args == (genes,)
    assert kwargs["database"] == "KEGG_2021_Human"
    assert result["gene_count"] == 3 and result["count"] == 1
    assert result["records"][0]["Term"] == "p53 pathway"
    with pytest.raises(ConfigurationError):
        gget_fns.gget_enrichr([], database="KEGG_2021_Human", verbose=False)


@pytest.mark.parametrize(
    ("name", "kwargs"),
    [
        ("specificity", {}),
        ("gene_expression", {"analysis_level": "Across_tissues", "analysis_type": "Sex"}),
        ("psi_block", {"analysis_level": "Across_tissues", "analysis_type": "Sex"}),
    ],
)
def test_gget_8cube_functions_forward_genes(gget_mod, monkeypatch, name, kwargs):
    calls = _recorder(monkeypatch, gget_mod, name, [{"gene": "SOX2", "value": 1.5}])
    wrapper = getattr(gget_fns, f"gget_{name}")
    result = wrapper(["SOX2", "PAX6"], verbose=False, **kwargs)
    args, call_kwargs = calls[0]
    assert args == (["SOX2", "PAX6"],)
    assert call_kwargs["json"] is True
    assert result["gene_count"] == 2 and result["count"] == 1
    with pytest.raises(ConfigurationError):
        wrapper([], verbose=False, **kwargs)


# --- mygene: MyGene.info batch gene lookup (client HTTP replaced) -----------------------


def _fake_post(monkeypatch, mygene_mod, rows):
    """Replace the HTTP POST of the biothings client; the real query code still runs."""
    from biothings_client.client.base import BiothingClient

    seen = []

    def fake(self, url, params=None, verbose=True):
        seen.append((url, dict(params or {})))
        return (False, rows)

    monkeypatch.setattr(BiothingClient, "_post", fake)
    return seen


def test_mygene_getgenes_keeps_input_order_and_reports_missing(mygene_mod, monkeypatch):
    rows = [
        {"query": "7157", "_id": "7157", "symbol": "TP53", "taxid": 9606},
        {"query": "999999999", "notfound": True},
        {"query": "1017", "_id": "1017", "symbol": "CDK2", "taxid": 9606},
    ]
    seen = _fake_post(monkeypatch, mygene_mod, rows)
    result = mygene_fns.mygene_getgenes(["1017", "999999999", "7157"], verbose=False)
    assert result["requested"] == 3 and result["found"] == 2
    assert [record["symbol"] for record in result["records"]] == ["CDK2", "TP53"]
    assert result["not_found"] == ["999999999"]
    assert seen[0][0].endswith("/gene/")
    assert seen[0][1]["species"] == "human"


def test_mygene_getgenes_validates_ids(mygene_mod):
    with pytest.raises(ConfigurationError):
        mygene_fns.mygene_getgenes([], verbose=False)
    with pytest.raises(ConfigurationError):
        mygene_fns.mygene_getgenes(["bad id with space"], verbose=False)


# --- myvariant: MyVariant.info batch variant lookup (client HTTP replaced) ---------------


def test_myvariant_querymany_orders_hits_and_lists_unmatched(myvariant_mod, monkeypatch):
    rows = [
        {"query": "rs4244285", "_id": "chr10:g.94781859G>C", "_score": 17.0},
        {"query": "rs58991260", "_id": "chr1:g.218458480G>A", "_score": 16.0},
        {"query": "rs4244285", "_id": "chr10:g.94781859G>A", "_score": 15.0},
        {"query": "rs0000000", "notfound": True},
    ]
    seen = _fake_post(monkeypatch, myvariant_mod, rows)
    terms = ["rs58991260", "rs4244285", "rs0000000"]
    result = myvariant_fns.myvariant_querymany(terms, assembly="hg19", verbose=False)
    assert result["requested"] == 3
    assert result["matched_terms"] == 2
    assert result["hit_count"] == 3
    assert [row["query"] for row in result["records"]] == [
        "rs58991260",
        "rs4244285",
        "rs4244285",
    ]
    assert result["not_found"] == ["rs0000000"]
    assert result["provenance"]["assembly"] == "hg19"
    assert seen[0][1]["assembly"] == "hg19"
    assert seen[0][1]["scopes"] == "dbsnp.rsid"


def test_myvariant_querymany_rejects_bad_terms(myvariant_mod):
    with pytest.raises(ConfigurationError):
        myvariant_fns.myvariant_querymany([], verbose=False)
    with pytest.raises(ConfigurationError):
        myvariant_fns.myvariant_querymany(["rs1 2"], verbose=False)
    with pytest.raises(ConfigurationError):
        myvariant_fns.myvariant_querymany(["rs1"], assembly="hg36", verbose=False)


# --- biothings_client: homolog lookup through the gene client -------------------------


class _FakeGeneClient:
    """Two-call fake of MyGene's getgenes: homologene annotations, then the homolog records."""

    def __init__(self, homologene, records):
        self.homologene = homologene
        self.records = records
        self.calls = []

    def getgenes(self, gene_list, **kwargs):
        ids = [str(item) for item in gene_list]
        self.calls.append((ids, dict(kwargs)))
        if kwargs.get("fields") == "homologene":
            return [{"_id": gid, "homologene": self.homologene[gid]} for gid in ids]
        return [self.records[gid] for gid in ids if gid in self.records]


def test_biothings_get_homologs_uses_real_join_logic(biothings_mod, monkeypatch):
    from biothings_client.utils import gene as gene_utils

    fake = _FakeGeneClient(
        homologene={"7157": {"genes": [[10090, 22059]]}},
        records={"22059": {"_id": "22059", "symbol": "Trp53"}},
    )
    monkeypatch.setattr(gene_utils, "get_client", lambda *args, **kwargs: fake)
    result = biothings_fns.biothings_client_get_homologs(["7157"], fields="symbol", verbose=False)
    # Hand calculation: TP53 (7157) has one homolog, mouse Trp53 (22059, taxon 10090).
    assert result["results"] == [
        {
            "gene": "7157",
            "homologs": [
                {
                    "taxid": 10090,
                    "genes": [{"_id": "22059", "symbol": "Trp53"}],
                }
            ],
        }
    ]
    assert result["requested"] == 1
    assert fake.calls[0][1]["fields"] == "homologene"


def test_biothings_get_homologs_species_filter_and_validation(biothings_mod, monkeypatch):
    from biothings_client.utils import gene as gene_utils

    fake = _FakeGeneClient(
        homologene={"7157": {"genes": [[10090, 22059]]}},
        records={"22059": {"_id": "22059", "symbol": "Trp53"}},
    )
    monkeypatch.setattr(gene_utils, "get_client", lambda *args, **kwargs: fake)
    result = biothings_fns.biothings_client_get_homologs(["7157"], species="9606", verbose=False)
    assert result["results"] == [{"gene": "7157", "homologs": []}]
    with pytest.raises(ConfigurationError):
        biothings_fns.biothings_client_get_homologs([], verbose=False)
    with pytest.raises(ConfigurationError):
        biothings_fns.biothings_client_get_homologs(["7157"], species="human", verbose=False)


# --- public surface -------------------------------------------------------------------


ALL_WRAPPERS = [
    gget_fns.gget_archs4,
    gget_fns.gget_bgee,
    gget_fns.gget_blast,
    gget_fns.gget_blat,
    gget_fns.gget_cosmic_local,
    gget_fns.gget_enrichr,
    gget_fns.gget_gene_expression,
    gget_fns.gget_opentargets,
    gget_fns.gget_psi_block,
    gget_fns.gget_specificity,
    mygene_fns.mygene_getgenes,
    myvariant_fns.myvariant_querymany,
    biothings_fns.biothings_client_get_homologs,
]


@pytest.mark.parametrize("function", ALL_WRAPPERS, ids=lambda fn: fn.__name__)
def test_each_wrapper_has_keyword_only_verbose(function):
    parameter = inspect.signature(function).parameters["verbose"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is True


def test_module_exports_match_wrappers():
    for module, names in (
        (gget_fns, gget_fns.__all__),
        (mygene_fns, mygene_fns.__all__),
        (myvariant_fns, myvariant_fns.__all__),
        (biothings_fns, biothings_fns.__all__),
    ):
        assert names, module.__name__
        for name in names:
            assert callable(getattr(module, name))


# --- optional live checks (one per package) ------------------------------------------


@live_only
def test_live_gget_opentargets_drugs_for_braf(gget_mod):
    result = gget_fns.gget_opentargets("ENSG00000157764", resource="drugs", limit=2)
    assert result["count"] == 2
    assert result["records"]


@live_only
def test_live_mygene_getgenes_cdk2(mygene_mod):
    result = mygene_fns.mygene_getgenes(["1017"], fields="symbol", verbose=False)
    assert result["records"][0]["symbol"] == "CDK2"


@live_only
def test_live_myvariant_querymany_rsid(myvariant_mod):
    result = myvariant_fns.myvariant_querymany(["rs58991260"], verbose=False)
    assert result["matched_terms"] == 1


@live_only
def test_live_biothings_homologs_tp53_mouse(biothings_mod):
    result = biothings_fns.biothings_client_get_homologs(["7157"], species="10090", verbose=False)
    symbols = [g["symbol"] for h in result["results"][0]["homologs"] for g in h["genes"]]
    assert symbols == ["Trp53"]


def test_upstream_failure_is_reported_as_backend_error(gget_mod, monkeypatch):
    def unavailable(*args, **kwargs):
        raise RuntimeError("service unavailable")

    monkeypatch.setattr(gget_mod, "opentargets", unavailable)
    with pytest.raises(BackendExecutionError, match="service unavailable"):
        gget_fns.gget_opentargets("ENSG00000157764", verbose=False)
