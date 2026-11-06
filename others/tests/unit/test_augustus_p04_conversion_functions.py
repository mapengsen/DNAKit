"""Tests for dnakit.annotation.augustus_p04_conversion_functions (AUGUSTUS format and sequence helpers)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.annotation import augustus_p04_conversion_functions as c
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
GENOME = DATA / "genome.fa"
GFF3 = DATA / "genes.gff3"
requires_augustus = pytest.mark.skipif(
    shutil.which("gtf2gff.pl") is None or shutil.which("perl") is None,
    reason="AUGUSTUS (conda install -c bioconda augustus) is not on PATH",
)
GTF = (
    "chr1\tAUGUSTUS\tgene\t12656\t14013\t0.04\t+\t.\tg50\n"
    "chr1\tAUGUSTUS\ttranscript\t12656\t14013\t0.04\t+\t.\tg50.t1\n"
    'chr1\tAUGUSTUS\tstart_codon\t12868\t12870\t.\t+\t0\ttranscript_id "g50.t1"; gene_id "g50";\n'
    'chr1\tAUGUSTUS\tCDS\t12868\t12993\t0.8\t+\t0\ttranscript_id "g50.t1"; gene_id "g50";\n'
    'chr1\tAUGUSTUS\tCDS\t13249\t13479\t1\t+\t0\ttranscript_id "g50.t1"; gene_id "g50";\n'
    'chr1\tAUGUSTUS\tstop_codon\t13477\t13479\t.\t+\t0\ttranscript_id "g50.t1"; gene_id "g50";\n'
)


def test_missing_program_and_validation(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        c.augustus_gtf_to_gff(
            str(GFF3), str(tmp_path / "o"), augustus_bin_dir=str(tmp_path), verbose=False
        )
    with pytest.raises(ConfigurationError):
        c.augustus_filter_genbank_genes(
            str(GFF3), str(GFF3), str(tmp_path / "o"), mode="other", verbose=False
        )
    with pytest.raises(ConfigurationError):
        c.augustus_extract_cds_sequences(
            str(GENOME),
            str(GFF3),
            str(tmp_path / "x"),
            method="stringtie2fa",
            annotation_format="gff3",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        c.augustus_gff_to_genbank(str(GFF3), str(GENOME), -1, str(tmp_path / "o"), verbose=False)


@requires_augustus
def test_gtf_conversions(tmp_path: Path) -> None:
    gtf = tmp_path / "g.gtf"
    gtf.write_text(GTF)
    gff3 = c.augustus_gtf_to_gff(
        str(gtf), str(tmp_path / "g.gff3"), gff3=True, print_exons=True, verbose=False
    )
    text = Path(gff3["outputs"]["output"]).read_text()
    assert "ID=g50;" in text and "\texon\t" in text
    named = c.augustus_add_names_to_gff3(
        gff3["outputs"]["output"],
        str(tmp_path / "named.gff3"),
        only_gene_structure=True,
        verbose=False,
    )
    assert "Name=g50.t1" in "\n".join(named["preview"])
    fixed = c.augustus_fix_joingenes_gtf(str(gtf), str(tmp_path / "fixed.gtf"), verbose=False)
    assert fixed["status"] == "completed"
    proteins = c.augustus_extract_cds_sequences(
        str(GENOME), str(gtf), str(tmp_path / "prot"), method="gtf2aa", verbose=False
    )
    assert proteins["records"] == {"protein": 1}


@requires_augustus
def test_genbank_round_trip(tmp_path: Path) -> None:
    gtf = tmp_path / "g.gtf"
    gtf.write_text(GTF)
    genbank = c.augustus_gff_to_genbank(
        str(gtf), str(GENOME), 300, str(tmp_path / "t.gb"), verbose=False
    )
    assert genbank["loci"] == 1
    loci = c.augustus_genbank_loci(
        genbank["outputs"]["output"], str(tmp_path / "loci.tsv"), verbose=False
    )
    assert loci["lines"] == 1
    back = c.augustus_genbank_to_gff(
        genbank["outputs"]["output"], str(tmp_path / "back.gff"), verbose=False
    )
    assert back["lines"] == 2  # two CDS parts
    kept = c.augustus_filter_genbank_genes(
        genbank["outputs"]["output"],
        loci["outputs"]["output"],
        str(tmp_path / "kept.gb"),
        mode="include_loci",
        verbose=False,
    )
    assert kept["status"] == "completed"
    split = c.augustus_random_split(
        genbank["outputs"]["output"],
        1,
        str(tmp_path / "train.gb"),
        str(tmp_path / "test.gb"),
        verbose=False,
    )
    assert split["train_loci"] + split["test_loci"] == 1


@requires_augustus
def test_sequence_and_masking_helpers(tmp_path: Path) -> None:
    summary = c.augustus_summarize_base_content(
        str(GENOME), str(tmp_path / "acgt.txt"), verbose=False
    )
    assert any("total 200000bp" in line for line in summary["preview"])
    mrna = c.augustus_extract_mrna_sequences(
        str(GFF3), str(GENOME), str(tmp_path / "mrna.fa"), verbose=False
    )
    assert mrna["records"] >= 1
    merged = c.augustus_merge_masking(
        str(GENOME), str(GENOME), str(tmp_path / "merged.fa"), verbose=False
    )
    assert merged["records"] == 2
    compared = c.augustus_compare_masking(
        str(GENOME), str(GENOME), str(tmp_path / "cmp.txt"), verbose=False
    )
    assert compared["preview"][0].startswith("Masked in both files")
    checked = c.augustus_check_gff_names_in_fasta(
        str(GFF3), str(GENOME), str(tmp_path / "chk.txt"), verbose=False
    )
    assert (
        checked["names_missing"] is True and checked["wrong_column_count"] is True
    )  # '##gff-version' line, 'demo' names
    clean = tmp_path / "clean.gff"
    clean.write_text("chr1\t.\tgene\t1\t2\t.\t+\t.\tID=a\n")
    ok = c.augustus_check_gff_names_in_fasta(
        str(clean), str(GENOME), str(tmp_path / "chk2.txt"), verbose=False
    )
    assert ok["names_missing"] is False and ok["wrong_column_count"] is False
    flanking = c.augustus_compute_flanking_region(
        str(tmp_path / "g.gtf") if (tmp_path / "g.gtf").exists() else _write(tmp_path),
        str(tmp_path / "flank.txt"),
        verbose=False,
    )
    assert flanking["numbers"][1] == 1.0


def _write(tmp_path: Path) -> str:
    path = tmp_path / "g.gtf"
    path.write_text(GTF)
    return str(path)


@requires_augustus
def test_small_format_converters(tmp_path: Path) -> None:
    cegma = tmp_path / "cegma.gff"
    cegma.write_text(
        "CHR\tcegma\tFirst\t100\t200\t5.0\t+\t0\tK1.1\nCHR\tcegma\tExon\t100\t200\t5.0\t+\t.\tK1.1\n"
    )
    converted = c.augustus_cegma_to_gff(str(cegma), str(tmp_path / "cds.gtf"), verbose=False)
    assert converted["preview"] == ['CHR\tcegma\tCDS\t100\t200\t5.0\t+\t0\ttranscript_id "g1.K1.1"']
    genepred = tmp_path / "gp.txt"
    genepred.write_text("g1\tchr1\t+\t100\t900\t0\t1\t2\t100,500\t300,900\n")
    big = c.augustus_convert_genepred(str(genepred), str(tmp_path / "big.gp"), verbose=False)
    assert big["lines"] == 1
    dependencies = tmp_path / "dep.txt"
    dependencies.write_text("q1 g1 g2\nq2 g1 g2\n")
    clusters = c.augustus_cluster_gene_dependencies(
        str(dependencies), str(tmp_path / "clusters.txt"), threshold=1, verbose=False
    )
    assert clusters["status"] == "completed"
