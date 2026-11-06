"""Tests for the genomepy catalogue transforms, compared with the upstream genomepy steps.

The upstream UCSC steps read the UCSC MySQL server through ``query_ucsc``. These tests replace that
function with a fake that returns fixed rows, so no network or database connection is opened.
"""

from __future__ import annotations

import copy
import inspect
import math

import pytest

pytest.importorskip("genomepy")
ucsc = pytest.importorskip("genomepy.providers.ucsc")
ensembl = pytest.importorskip("genomepy.providers.ensembl")
gencode = pytest.importorskip("genomepy.providers.gencode")

from dnakit.databases import genomepy_verify_functions as gv
from dnakit.exceptions import ConfigurationError

PUBLIC = [
    gv.genomepy_ucsc_assembly_accessions,
    gv.genomepy_ucsc_sourcename_accessions,
    gv.genomepy_ucsc_annotation_links,
    gv.genomepy_ensembl_add_grch37,
    gv.genomepy_gencode_add_grch37,
]

EQUIVALENCES = [
    ("hg38", "GCA_000001405.28_GRCh38.p13", 3),
    ("hg38", "GCF_000001405.39_GRCh38.p13", 2),  # lower match, dropped
    ("mm10", "GCF_000001635.27_GRCm39", 5),
    ("ce11", "GCF_000002985.6_WBcel235", 4),  # tie with GCA: only GCA kept
    ("ce11", "GCA_000002985.3_WBcel235", 4),
    ("dm6", "GCF_000001215.4_Release_6_plus_ISO1_MT", 7),
    ("danRer11", "GCF_000002035.6_GRCz11", 3),
    ("rn7", "GCF_015227675.2_mRatBN7.2", 1),
    ("xenTro10", "GCF_000004195.4_UCB_Xtro_10.0", 2),  # two GCF at the maximum: dropped upstream
    ("xenTro10", "GCF_000004195.3_Xenopus_tropicalis_v9.0", 2),
    ("unknownAsm", "GCF_000000001.1_Fake", 9),  # not in the catalogue: ignored
]

CATALOGUE_UCSC = {
    name: {"name": name, "assembly_accession": None, "sourceName": None, "annotations": []}
    for name in ["hg38", "mm10", "ce11", "dm6", "danRer11", "rn7", "xenTro10", "hg19"]
}


def _fake_query(rows):
    def query(command, database=None):
        return list(rows)

    return query


def test_public_functions_take_keyword_only_verbose():
    for function in PUBLIC:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_assembly_accessions_match_upstream_add_accessions1(monkeypatch):
    genomes = copy.deepcopy(CATALOGUE_UCSC)
    monkeypatch.setattr(ucsc, "query_ucsc", _fake_query(EQUIVALENCES))
    upstream = ucsc.add_accessions1(copy.deepcopy(genomes))
    ours = gv.genomepy_ucsc_assembly_accessions(genomes, EQUIVALENCES)
    assert set(ours) == set(genomes)
    for name in genomes:
        assert ours[name].get("assembly_accession") == upstream[name].get("assembly_accession"), (
            name
        )
    assert ours["hg38"]["assembly_accession"] == "GCA_000001405.28"
    assert ours["ce11"]["assembly_accession"] == "GCA_000002985.3"
    assert ours["xenTro10"]["assembly_accession"] is None


def test_assembly_accessions_nan_case_is_a_documented_difference(monkeypatch):
    rows = [("hg19", "no_accession_in_this_name", 1)]
    genomes = copy.deepcopy(CATALOGUE_UCSC)
    monkeypatch.setattr(ucsc, "query_ucsc", _fake_query(rows))
    upstream = ucsc.add_accessions1(copy.deepcopy(genomes))
    ours = gv.genomepy_ucsc_assembly_accessions(genomes, rows)
    assert isinstance(upstream["hg19"]["assembly_accession"], float)
    assert math.isnan(upstream["hg19"]["assembly_accession"])
    assert ours["hg19"]["assembly_accession"] is None


def test_assembly_accessions_empty_rows_return_unchanged_copy():
    genomes = copy.deepcopy(CATALOGUE_UCSC)
    ours = gv.genomepy_ucsc_assembly_accessions(genomes, [])
    assert ours == genomes
    assert ours is not genomes


def test_sourcename_accessions_match_upstream_add_accessions2():
    genomes = {
        "hg38": {"sourceName": "GRCh38 GCA_000001405.28 Dec 2013", "annotations": []},
        "mm10": {"sourceName": "GRCm38 (GCF_000001635.26)", "annotations": []},
        "rn7": {"sourceName": "rn7 without accession", "annotations": []},
        "dm6": {"sourceName": "Release 6 GCF_000001215.4 and GCA_000001215.9", "annotations": []},
    }
    upstream = ucsc.add_accessions2(copy.deepcopy(genomes))
    ours = gv.genomepy_ucsc_sourcename_accessions(genomes)
    for name in genomes:
        assert ours[name].get("assembly_accession") == upstream[name].get("assembly_accession"), (
            name
        )
    assert ours["dm6"]["assembly_accession"] == "GCF_000001215.4"
    assert "assembly_accession" not in ours["rn7"]


def test_sourcename_accessions_skip_none_and_reject_missing_field():
    ours = gv.genomepy_ucsc_sourcename_accessions({"hg38": {"sourceName": None}})
    assert "assembly_accession" not in ours["hg38"]
    with pytest.raises(ConfigurationError):
        gv.genomepy_ucsc_sourcename_accessions({"hg38": {"annotations": []}})


def test_annotation_links_match_upstream_add_annotation_links(monkeypatch):
    rows = [
        ("hg38", "ensGene"),
        ("hg38", "refGene"),
        ("mm10", "knownGene"),
        ("notInCatalogue", "ncbiRefSeq"),
        ("hg38", "ncbiRefSeq"),
    ]
    genomes = {
        "hg38": {"annotations": []},
        "mm10": {"annotations": ["ensGene"]},
        "rn7": {"annotations": []},
    }
    monkeypatch.setattr(ucsc, "query_ucsc", _fake_query(rows))
    upstream = ucsc.add_annotation_links(copy.deepcopy(genomes))
    ours = gv.genomepy_ucsc_annotation_links(genomes, rows)
    for name in genomes:
        assert ours[name]["annotations"] == upstream[name]["annotations"], name
    assert ours["hg38"]["annotations"] == ["ensGene", "refGene", "ncbiRefSeq"]
    assert ours["mm10"]["annotations"] == ["ensGene", "knownGene"]


def test_annotation_links_need_annotations_list():
    with pytest.raises(ConfigurationError):
        gv.genomepy_ucsc_annotation_links({"hg38": {}}, [("hg38", "ensGene")])


def test_ensembl_add_grch37_matches_upstream():
    genomes = {"GRCh38": {"name": "Homo_sapiens", "assembly_name": "GRCh38"}}
    upstream = ensembl.add_grch37(copy.deepcopy(genomes))
    ours = gv.genomepy_ensembl_add_grch37(genomes)
    assert ours["GRCh37"] == upstream["GRCh37"]
    assert ours["GRCh38"] == genomes["GRCh38"]
    assert set(ours) == {"GRCh38", "GRCh37"}


def test_gencode_add_grch37_matches_upstream():
    ftp = "ftp://ftp.ebi.ac.uk/pub/databases/gencode"
    catalogue = {
        "GRCh38": {
            "annotations": [
                f"{ftp}/Gencode_human/release_41/GRCh38_mapping/gencode.v41.annotation.gtf.gz"
            ],
            "taxonomy_id": 9606,
            "species": "Homo sapiens",
            "text_search": "human",
        },
        "GRCh37": {"annotations": ["old"]},
    }
    upstream = gencode.add_grch37(copy.deepcopy(catalogue), ftp)
    ours = gv.genomepy_gencode_add_grch37(catalogue, ftp)
    assert ours["GRCh37"] == upstream["GRCh37"]
    assert ours["GRCh37"]["annotations"] == [
        f"{ftp}/Gencode_human/release_41/GRCh37_mapping/gencode.v41lift37.annotation.gtf.gz"
    ]


def test_gencode_add_grch37_rejects_catalogue_without_human_assembly():
    with pytest.raises(ConfigurationError):
        gv.genomepy_gencode_add_grch37({"mm10": {"annotations": []}}, "ftp://x")
    with pytest.raises(ConfigurationError):
        gv.genomepy_gencode_add_grch37(
            {"GRCh38": {"annotations": ["no-release-here.gtf"]}}, "ftp://x"
        )
    with pytest.raises(ConfigurationError):
        gv.genomepy_gencode_add_grch37({"GRCh38": {"annotations": ["a/release_41/b"]}}, " ")


def test_inputs_are_not_modified(monkeypatch):
    genomes = copy.deepcopy(CATALOGUE_UCSC)
    before = copy.deepcopy(genomes)
    gv.genomepy_ucsc_assembly_accessions(genomes, EQUIVALENCES)
    gv.genomepy_ucsc_annotation_links(genomes, [("hg38", "ensGene")])
    assert genomes == before


def test_malformed_rows_are_rejected():
    with pytest.raises(ConfigurationError):
        gv.genomepy_ucsc_assembly_accessions(CATALOGUE_UCSC, [("hg38", "GCA_000001405.28_x")])
    with pytest.raises(ConfigurationError):
        gv.genomepy_ucsc_annotation_links(CATALOGUE_UCSC, "hg38 ensGene")
