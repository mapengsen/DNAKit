"""miniprot, Liftoff and Flye facades (each skips when its program is missing)."""

from __future__ import annotations

import random
from pathlib import Path

import pytest
from _p16_data import genome, need, noisy_reads, rand_dna, write_fasta

import dnakit.annotation.p16_assembly_annotation_cli_functions as an
from dnakit.exceptions import ConfigurationError


def test_miniprot_index_and_alignment(tmp_path: Path) -> None:
    program = need("miniprot")
    bio_seq = pytest.importorskip("Bio.Seq")
    contigs = genome(tmp_path)
    protein = str(bio_seq.Seq(contigs["chrA"][3000:3600]).translate()).replace("*", "A")
    proteins = write_fasta(tmp_path / "prot.faa", {"prot1": protein})
    index = an.build_miniprot_index(
        str(tmp_path / "genome.fa"), str(tmp_path / "g.mpi"), miniprot_path=program, verbose=False
    )
    paf = an.align_proteins_to_genome_miniprot(
        index["index"], proteins, str(tmp_path / "a.paf"), miniprot_path=program, verbose=False
    )
    fields = paf["head"][0].split("\t")
    assert fields[0] == "prot1" and fields[5] == "chrA" and abs(int(fields[7]) - 3000) < 20
    gff = an.align_proteins_to_genome_miniprot(
        str(tmp_path / "genome.fa"),
        proteins,
        str(tmp_path / "a.gff"),
        output_format="gff",
        translate=True,
        miniprot_path=program,
        verbose=False,
    )
    assert any("\tmRNA\t" in line for line in gff["head"]) or any(
        "##PAF" in line or "\tCDS\t" in line for line in gff["head"]
    )
    with pytest.raises(ConfigurationError):
        an.align_proteins_to_genome_miniprot(
            str(tmp_path / "genome.fa"),
            proteins,
            str(tmp_path / "b"),
            output_format="sam",
            verbose=False,
        )


def test_liftoff_lifts_genes_to_a_modified_target(tmp_path: Path) -> None:
    program = need("liftoff")
    need("minimap2")
    contigs = genome(tmp_path)
    rng = random.Random(8)
    reference = contigs["chrA"]
    target = reference[:500] + rand_dna(rng, 120) + reference[500:]
    write_fasta(tmp_path / "target.fa", {"chrA": target})
    (tmp_path / "ref.gff3").write_text(
        "##gff-version 3\nchrA\ttest\tgene\t1000\t2500\t.\t+\t.\tID=gene1\nchrA\ttest\tmRNA\t1000\t2500\t.\t+\t.\tID=m1;Parent=gene1\n"
        "chrA\ttest\texon\t1000\t2500\t.\t+\t.\tID=e1;Parent=m1\nchrA\ttest\tCDS\t1000\t2500\t.\t+\t0\tID=c1;Parent=m1\n"
    )
    result = an.lift_over_annotation_liftoff(
        str(tmp_path / "target.fa"),
        str(tmp_path / "genome.fa"),
        str(tmp_path / "ref.gff3"),
        str(tmp_path / "lift.gff3"),
        extra_copies=True,
        liftoff_path=program,
        verbose=False,
    )
    assert result["lifted_genes"] == 1 and result["unmapped_features"] == 0
    gene = next(
        line for line in (tmp_path / "lift.gff3").read_text().splitlines() if "\tgene\t" in line
    ).split("\t")
    assert int(gene[3]) == 1120 and int(gene[4]) == 2620  # shifted by the 120 bp insertion


def test_flye_assembles_and_polishes(tmp_path: Path) -> None:
    program = need("flye")
    contigs = genome(tmp_path)
    reads = noisy_reads({"chrA": contigs["chrA"]}, 40, tmp_path / "reads.fq")
    result = an.assemble_long_reads_flye(
        reads,
        str(tmp_path / "flye"),
        read_type="nano-raw",
        genome_size="20k",
        iterations=0,
        flye_path=program,
        verbose=False,
    )
    assert result["contigs"] and max(c["length"] for c in result["contigs"]) > 15000
    polished = an.polish_assembly_flye(
        str(tmp_path / "flye" / "assembly.fasta"),
        reads,
        str(tmp_path / "pol"),
        iterations=1,
        flye_path=program,
        verbose=False,
    )
    assert Path(polished["polished"]).stat().st_size > 15000
