"""Tests for dnakit.annotation.genomicranges_gap_functions (UCSC GTF URL and table reading).

The network path is only exercised when DNAKIT_LIVE_TESTS=1. Other tests use a local GTF file.
"""

from __future__ import annotations

import gzip
import inspect
import os

import pytest

from dnakit.exceptions import ConfigurationError

pytest.importorskip("genomicranges")
pytest.importorskip("pandas")

from dnakit.annotation import genomicranges_gap_functions as facades

GTF_LINES = (
    'chr1\tdemo\tgene\t100\t500\t.\t+\t.\tgene_id "g1"; gene_name "G1";\n'
    'chr1\tdemo\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr2\tdemo\texon\t50\t90\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n'
)


def test_public_functions_are_keyword_only_verbose():
    for name in facades.__all__:
        parameters = inspect.signature(getattr(facades, name)).parameters
        assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameters["verbose"].default is True, name


def test_ucsc_url_for_each_annotation():
    assert facades.ucsc_gtf_url("hg38", verbose=False) == (
        "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
    )
    assert facades.ucsc_gtf_url("mm10", "ensGene", verbose=False).endswith("/mm10.ensGene.gtf.gz")


def test_ucsc_url_rejects_unknown_annotation_and_bad_genome():
    with pytest.raises(ConfigurationError):
        facades.ucsc_gtf_url("hg38", "notATable", verbose=False)
    with pytest.raises(ConfigurationError):
        facades.ucsc_gtf_url("../hg38", verbose=False)


def test_read_local_gtf_writes_table_and_summary(tmp_path):
    source = tmp_path / "genes.gtf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(GTF_LINES)
    out = tmp_path / "genes.tsv"
    result = facades.read_ucsc_gtf("hg38", str(out), local_gtf=str(source), verbose=False)
    assert result["row_count"] == 3
    assert result["rows_per_chromosome"] == {"chr1": 2, "chr2": 1}
    assert result["source"] == str(source.resolve())
    assert out.is_file()
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert header.split("\t")[:3] == ["seqnames", "source", "feature"]


def test_read_local_plain_gtf(tmp_path):
    source = tmp_path / "plain.gtf"
    source.write_text(GTF_LINES, encoding="utf-8")
    result = facades.read_ucsc_gtf(
        "hg38", str(tmp_path / "plain.tsv"), local_gtf=str(source), verbose=False
    )
    assert result["row_count"] == 3


def test_read_refuses_existing_output_and_missing_local_file(tmp_path):
    existing = tmp_path / "exists.tsv"
    existing.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        facades.read_ucsc_gtf(
            "hg38", str(existing), local_gtf=str(tmp_path / "none.gtf"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        facades.read_ucsc_gtf(
            "hg38", str(tmp_path / "new.tsv"), local_gtf=str(tmp_path / "none.gtf"), verbose=False
        )


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="network test; set DNAKIT_LIVE_TESTS=1"
)
def test_live_ucsc_download_summary(tmp_path):
    result = facades.read_ucsc_gtf("hg38", str(tmp_path / "refgene.tsv"), verbose=False)
    assert result["row_count"] > 0
    assert result["source"].endswith("hg38.refGene.gtf.gz")
