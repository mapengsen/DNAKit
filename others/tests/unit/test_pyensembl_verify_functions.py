"""Tests for pyensembl_verify_functions (the pyensembl install action).

The local-GTF install is checked against a tiny GTF written into tmp_path, with pyensembl's cache in
tmp_path; the gene and transcript IDs read back must be the ones in the GTF. The Ensembl-release
install reaches the network and only runs when DNAKIT_LIVE_TESTS=1. Ordinary tests block sockets.
"""

from __future__ import annotations

import importlib.util
import os
import socket

import pytest

from dnakit.annotation import pyensembl_verify_functions as module
from dnakit.exceptions import ConfigurationError

requires_pyensembl = pytest.mark.skipif(
    importlib.util.find_spec("pyensembl") is None, reason="pyensembl is not installed"
)
live_only = pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="set DNAKIT_LIVE_TESTS=1 for the Ensembl download",
)

GTF = """chr1\ttest\tgene\t100\t900\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_biotype "protein_coding";
chr1\ttest\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_name "GENE1-201"; transcript_biotype "protein_coding";
chr1\ttest\texon\t100\t300\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_name "GENE1-201"; transcript_biotype "protein_coding";
chr1\ttest\texon\t500\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_name "GENE1-201"; transcript_biotype "protein_coding";
chr2\ttest\tgene\t1000\t1500\t.\t-\t.\tgene_id "G2"; gene_name "GENE2"; gene_biotype "protein_coding";
chr2\ttest\ttranscript\t1000\t1500\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2"; gene_biotype "protein_coding"; transcript_name "GENE2-201"; transcript_biotype "protein_coding";
chr2\ttest\texon\t1000\t1500\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; exon_number "1"; gene_name "GENE2"; gene_biotype "protein_coding"; transcript_name "GENE2-201"; transcript_biotype "protein_coding";
"""


_ORIGINAL_GETADDRINFO = socket.getaddrinfo
_ORIGINAL_CONNECT = socket.socket.connect


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    # Name resolution and outbound connections are blocked; local socket pairs (used by asyncio in zarr) stay.
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket, "getaddrinfo", _blocked)
    monkeypatch.setattr(socket.socket, "connect", _blocked)


@requires_pyensembl
def test_local_gtf_install_indexes_the_annotation(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    cache = tmp_path / "cache"
    cache.mkdir()
    result = module.pyensembl_install_local_annotation(
        "TestRef", "testannot", str(gtf), cache_directory_path=str(cache), verbose=False
    )
    assert result["status"] == "installed"
    assert result["genes"] == 2
    assert result["transcripts"] == 2
    assert result["gtf"] == "genes.gtf"

    import pyensembl

    genome = pyensembl.Genome(
        reference_name="TestRef",
        annotation_name="testannot",
        gtf_path_or_url=str(gtf),
        cache_directory_path=str(cache),
    )
    assert sorted(genome.gene_ids()) == ["G1", "G2"]
    assert sorted(genome.transcript_ids()) == ["T1", "T2"]


@requires_pyensembl
def test_missing_gtf_is_rejected_before_any_work(tmp_path):
    with pytest.raises(ConfigurationError):
        module.pyensembl_install_local_annotation(
            "TestRef", "testannot", str(tmp_path / "absent.gtf"), verbose=False
        )


@requires_pyensembl
def test_bad_annotation_version_is_rejected(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    with pytest.raises(ConfigurationError):
        module.pyensembl_install_local_annotation(
            "TestRef", "testannot", str(gtf), annotation_version=0, verbose=False
        )


@requires_pyensembl
@live_only
def test_ensembl_release_install_downloads_and_indexes(tmp_path, monkeypatch):
    monkeypatch.setenv("PYENSEMBL_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(socket, "getaddrinfo", _ORIGINAL_GETADDRINFO)
    monkeypatch.setattr(socket.socket, "connect", _ORIGINAL_CONNECT)
    result = module.pyensembl_install_ensembl_release(110, "caenorhabditis_elegans", verbose=False)
    assert result["status"] == "installed"
    assert result["genes"] > 0
    assert result["transcripts"] >= result["genes"]
