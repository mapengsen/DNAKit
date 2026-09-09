"""Optional real-binary integration: skipped unless the named programs are on PATH."""

from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.bioinformatics import alignment, cloning, ngs, optimization, reads


def _programs(*names: str) -> None:
    if missing := [name for name in names if shutil.which(name) is None]:
        pytest.skip(f"Optional real tools not installed: {missing}")


def _sequence(length: int) -> str:
    rng = random.Random(12)
    return "".join(rng.choices("ACGT", k=length))


def test_real_mafft_and_nucleotide_hmm_search(tmp_path: Path) -> None:
    _programs("mafft", "hmmbuild", "nhmmer")
    pytest.importorskip("Bio")
    sequence = _sequence(120)
    aligned = alignment.align_sequences(
        {"first": sequence, "second": sequence[:50] + "A" + sequence[50:]}, tmp_path / "msa"
    )
    assert aligned.ids == ("first", "second")
    assert aligned.sequences[0].replace("-", "") == sequence
    stockholm = alignment.write_alignment(aligned, tmp_path / "input.sto", format="stockholm")
    profile = alignment.build_profile(stockholm, tmp_path / "profile")
    fasta = tmp_path / "target.fa"
    fasta.write_text(">target\n" + sequence + "\n")
    search = alignment.search_profile(
        profile.artifacts["profile.hmm"], fasta, tmp_path / "hmm-search"
    )
    hits = tuple(alignment.parse_nhmmer_table(search.artifacts["hits.tbl"]))
    assert any(hit["target"] == "target" and hit["strand"] == "+" for hit in hits)


def test_real_blast_database_search_and_constraint(tmp_path: Path) -> None:
    _programs("blastn", "makeblastdb")
    pytest.importorskip("Bio")
    sequence = _sequence(200)
    fasta = tmp_path / "target.fa"
    fasta.write_text(">target\n" + sequence + "\n")
    database = alignment.build_blast_database(fasta, tmp_path / "database")
    prefix = database.directory / "database"
    result = alignment.local_blast(fasta, prefix, tmp_path / "search")
    queries = tuple(alignment.parse_search_results(result.artifacts["results.xml"]))
    assert queries[0][0].id == "target"
    rule = optimization.blast_match_specification(prefix, minimum_length=30)
    assert rule.evaluate(sequence).violation > 0
    assert rule.evaluate("A" * 80).violation == 0


def test_real_reads_to_bam_and_variants(tmp_path: Path) -> None:
    _programs("bwa", "samtools", "bcftools")
    pytest.importorskip("pysam")
    pytest.importorskip("numpy")
    reference = _sequence(400)
    alternate = next(base for base in "ACGT" if base != reference[150])
    mutant = reference[:150] + alternate + reference[151:]
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + reference + "\n")
    reads.index_reference(fasta)
    fastq = tmp_path / "reads.fq"
    reads.write_fastx(
        [reads.Read(f"read{i}", mutant[100 + i : 200 + i], (40,) * 100) for i in range(12)],
        fastq,
    )
    result = ngs.dna_seq_workflow(fasta, fastq, tmp_path / "workflow", ploidy=1, qc=False)
    summary = reads.alignment_summary(result["bam"].artifacts["aligned.bam"])
    assert summary["mapped"] == 12
    pysam = pytest.importorskip("pysam")

    with pysam.VariantFile(str(result["variants"].artifacts["variants.vcf.gz"])) as handle:
        variants = list(handle)
    assert any(record.pos == 151 and record.alts == (alternate,) for record in variants)


def test_real_pydna_overlap_assembly() -> None:
    pytest.importorskip("pydna")
    sequence = _sequence(120)
    fragments = [
        cloning.Duplex.from_sequence(sequence[:80]),
        cloning.Duplex.from_sequence(sequence[50:]),
    ]
    products = cloning.assemble_duplexes(fragments, minimum_overlap=25)
    assert any(product.watson == sequence for product in products)
