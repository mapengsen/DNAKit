"""Tests for dnakit.annotation.gtars_gap_functions (HGVS to VRS allele identifier on local stores).

The fixture builds a one-transcript transcript store and a refget store from a random 400 bp FASTA in
tmp_path; nothing is downloaded.
"""

from __future__ import annotations

import inspect
import random
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("gtars")
reftx = pytest.importorskip("gtars.reftx")
refget = pytest.importorskip("gtars.refget")

from dnakit.annotation import gtars_gap_functions as facades  # noqa: E402

CONTIG = "NC_TEST.1"


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    root = tmp_path_factory.mktemp("gtars_stores")
    random.seed(5)
    sequence = "".join(random.choice("ACGT") for _ in range(400))
    fasta = root / "genome.fa"
    fasta.write_text(f">{CONTIG}\n{sequence}\n", encoding="utf-8")
    digest = refget.digest_fasta(str(fasta)).digest
    chrom = refget.sha512t24u_digest(sequence)
    builder = reftx.TxStoreBuilder()
    builder.add_transcript(
        {
            "accession": "NM_TEST.1",
            "gene": "TESTG",
            "chrom": chrom,
            "strand": "+",
            "cds_start": 50,
            "cds_end": 150,
            "exons": [{"start": 0, "end": 200}],
            "mane": None,
        }
    )
    builder.add_chrom_mapping(CONTIG, chrom)
    transcript_store = root / "tx.reftx"
    builder.build(str(transcript_store))
    store_dir = root / "refget_store"
    store = refget.RefgetStore.on_disk(str(store_dir))
    store.add_sequence_collection_from_fasta(str(fasta))
    return {"transcript_store": str(transcript_store), "refget_store": str(store_dir), "digest": digest}


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_deletion_gets_a_vrs_allele_identifier(stores):
    result = facades.hgvs_to_vrs_allele_id(
        "NM_TEST.1:c.10del", stores["transcript_store"], stores["refget_store"], stores["digest"], verbose=False
    )
    assert result["vrs_allele_id"].startswith("ga4gh:VA.")
    assert result["backend"] == "gtars"
    assert result["version"] != "unknown"


def test_identifier_is_deterministic_and_distinguishes_variants(stores):
    args = (stores["transcript_store"], stores["refget_store"], stores["digest"])
    first = facades.hgvs_to_vrs_allele_id("NM_TEST.1:c.20_21insT", *args, verbose=False)["vrs_allele_id"]
    again = facades.hgvs_to_vrs_allele_id("NM_TEST.1:c.20_21insT", *args, verbose=False)["vrs_allele_id"]
    other = facades.hgvs_to_vrs_allele_id("NM_TEST.1:c.10del", *args, verbose=False)["vrs_allele_id"]
    assert first == again
    assert first != other


def test_reference_mismatch_is_a_configuration_error(stores):
    with pytest.raises(ConfigurationError):
        facades.hgvs_to_vrs_allele_id(
            "NM_TEST.1:c.5A>G", stores["transcript_store"], stores["refget_store"], stores["digest"], verbose=False
        )


def test_inputs_are_checked(stores, tmp_path):
    with pytest.raises(ConfigurationError):
        facades.hgvs_to_vrs_allele_id(
            "NM_TEST.1:c.10del", str(tmp_path / "none.reftx"), stores["refget_store"], stores["digest"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.hgvs_to_vrs_allele_id(
            "NM_TEST.1:c.10del", stores["transcript_store"], str(tmp_path / "none"), stores["digest"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.hgvs_to_vrs_allele_id("", stores["transcript_store"], stores["refget_store"], stores["digest"], verbose=False)
    with pytest.raises(ConfigurationError):
        facades.hgvs_to_vrs_allele_id("NM_TEST.1:c.10del", stores["transcript_store"], stores["refget_store"], " ", verbose=False)
