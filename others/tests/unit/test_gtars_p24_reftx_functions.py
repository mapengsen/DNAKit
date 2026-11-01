"""Tests for dnakit.annotation.gtars_p24_reftx_functions (transcript stores built in tmp_path; no
network)."""

from __future__ import annotations

import inspect
import json

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("gtars")
pytest.importorskip("gtars.reftx")

from dnakit.annotation import gtars_p24_reftx_functions as m

CHROM = "ixSNGKM4iP6C8F3mkCr_2N0I8RB4mD5U"


def _tx(accession, gene, strand, mane):
    return {
        "accession": accession,
        "gene": gene,
        "chrom": CHROM,
        "strand": strand,
        "cds_start": 150,
        "cds_end": 950,
        "exons": [{"start": 100, "end": 400}, {"start": 600, "end": 1000}],
        "mane": mane,
    }


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    root = tmp_path_factory.mktemp("p24_reftx")
    path = root / "tx.reftx"
    m.transcript_store_build(
        path,
        transcripts=[
            _tx("NM_1.1", "G1", "+", {"select": True, "plus_clinical": False}),
            _tx("NM_2.1", "G2", "-", None),
        ],
        verbose=False,
    )
    return path


def test_public_functions_are_keyword_only_verbose():
    for name in m.__all__:
        parameters = inspect.signature(getattr(m, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_lookup_returns_record_strand_and_chromosome(store):
    found = m.transcript_record_lookup(store, "NM_2.1", verbose=False)
    assert found["found"] and found["strand"] == "-"
    assert found["chromosome_accession"].endswith(CHROM)
    assert found["transcript"]["exons"][1] == {"start": 600, "end": 1000}
    assert m.transcript_record_lookup(store, "NM_404.1", verbose=False)["found"] is False


def test_coordinates_plus_minus_and_gene_mode(store):
    plus = m.transcript_coordinate_to_genomic(store, 10, accession="NM_1.1", verbose=False)
    assert plus["genomic_pos"] == 159 and plus["strand"] == "+"
    minus = m.transcript_coordinate_to_genomic(
        store, 5, accession="NM_2.1", coordinate_type="n", verbose=False
    )
    assert minus["genomic_pos"] == 995 and minus["strand"] == "-"
    by_gene = m.transcript_coordinate_to_genomic(store, 10, gene="G1", verbose=False)
    assert by_gene["genomic_pos"] == plus["genomic_pos"] and by_gene["accession"] == "NM_1.1"
    shifted = m.transcript_coordinate_to_genomic(
        store, 10, accession="NM_1.1", datum=1, verbose=False
    )
    assert shifted["genomic_pos"] == 959


def test_coordinate_errors(store):
    with pytest.raises(ConfigurationError):
        m.transcript_coordinate_to_genomic(store, 10, accession="NOPE.1", verbose=False)
    with pytest.raises(ConfigurationError):
        m.transcript_coordinate_to_genomic(store, 10, verbose=False)
    with pytest.raises(ConfigurationError):
        m.transcript_coordinate_to_genomic(store, 10, gene="G1", coordinate_type="n", verbose=False)
    with pytest.raises(ConfigurationError):
        m.transcript_coordinate_to_genomic(
            store, 10, accession="NM_1.1", coordinate_type="x", verbose=False
        )


def test_build_from_cdot_and_mane(tmp_path):
    cdot = tmp_path / "c.json"
    cdot.write_text(
        json.dumps(
            {
                "transcripts": {
                    "NM_9.1": {
                        "id": "NM_9.1",
                        "gene_name": "GX",
                        "contig": "NC_000001.11",
                        "strand": 1,
                        "cds_start": 150,
                        "cds_end": 950,
                        "exons": [[100, 400], [600, 1000]],
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    mane = tmp_path / "mane.tsv"
    mane.write_text(
        "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\t"
        "MANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n"
        "1\tENSG1\tHGNC:1\tGX\tgene x\tNM_9.1\tNP_9.1\tENST1\tENSP1\tMANE Select\t"
        "NC_000001.11\t100\t1000\t+\n",
        encoding="utf-8",
    )
    out = tmp_path / "built.reftx"
    result = m.transcript_store_build(
        out,
        cdot_files=[cdot],
        mane_summary=mane,
        chrom_mappings={"NC_000001.11": CHROM},
        verbose=False,
    )
    assert result["cdot_transcripts_loaded"] == 1 and result["transcript_count"] == 1
    assert m.transcript_record_lookup(out, "NM_9.1", verbose=False)["transcript"]["gene"] == "GX"
    with pytest.raises(ConfigurationError):
        m.transcript_store_build(out, cdot_files=[cdot], verbose=False)


def test_build_input_checks(tmp_path):
    with pytest.raises(ConfigurationError):
        m.transcript_store_build(tmp_path / "e.reftx", verbose=False)
    with pytest.raises(ConfigurationError):
        m.transcript_store_build(
            tmp_path / "e.reftx", transcripts=[{"accession": "x"}], verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.transcript_record_lookup(tmp_path / "missing.reftx", "NM_1.1", verbose=False)
