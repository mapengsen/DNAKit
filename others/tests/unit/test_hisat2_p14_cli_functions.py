"""Tests for the HISAT2 facades; program-dependent tests skip when HISAT2 is not installed."""

from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from dnakit.alignment import hisat2_p14_cli_functions as hs
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

HAS = shutil.which("hisat2") is not None and shutil.which("hisat2-build") is not None
needs_hisat2 = pytest.mark.skipif(not HAS, reason="HISAT2 is not installed")

GTF = (
    'chr1\tdemo\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tdemo\texon\t601\t800\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tdemo\texon\t1001\t1200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)


@pytest.fixture()
def genome(tmp_path):
    rng = random.Random(8)
    sequence = "".join(rng.choice("ACGT") for _ in range(3000))
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "\n".join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + "\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GTF)
    return fasta, gtf, sequence


def test_missing_program_raises_install_hint(tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(GTF)
    with pytest.raises(BackendUnavailableError):
        hs.hisat2_extract_exons(str(gtf), str(tmp_path / "e.txt"), hisat2_path=str(tmp_path))


def test_validation(tmp_path, genome):
    fasta, gtf, _ = genome
    with pytest.raises(ConfigurationError):
        hs.hisat2_align_reads(str(tmp_path / "noindex"), str(tmp_path / "o.sam"), unpaired_reads=[str(fasta)])
    with pytest.raises(ConfigurationError):
        hs.hisat2_extract_snps_haplotypes_vcf(str(fasta), [], str(tmp_path / "x"))
    with pytest.raises(ConfigurationError):
        hs.hisat2_extract_snps_haplotypes_vcf(str(fasta), [str(gtf)], str(tmp_path / "x"), reference_type="bad")


@needs_hisat2
def test_annotation_helpers_and_index_inspection(tmp_path, genome):
    fasta, gtf, _ = genome
    splice = hs.hisat2_extract_splice_sites(str(gtf), str(tmp_path / "ss.txt"), verbose=False)
    assert splice["rows"] == [["chr1", "399", "600", "+"], ["chr1", "799", "1000", "+"]]
    exons = hs.hisat2_extract_exons(str(gtf), str(tmp_path / "ex.txt"), verbose=False)
    assert exons["rows"][0] == ["chr1", "200", "399", "+"]
    subprocess.run(["hisat2-build", "--ss", str(tmp_path / "ss.txt"), "--exon", str(tmp_path / "ex.txt"), str(fasta), str(tmp_path / "idx")], check=True, capture_output=True)
    sites = hs.hisat2_inspect_index(str(tmp_path / "idx"), str(tmp_path / "ins_ss.txt"), view="splice_sites", verbose=False)
    assert sites["row_count"] == 2
    assert hs.hisat2_inspect_index(str(tmp_path / "idx"), str(tmp_path / "ins_ex.txt"), view="exons", verbose=False)["row_count"] == 3
    names = hs.hisat2_inspect_index(str(tmp_path / "idx"), str(tmp_path / "names.txt"), view="names", verbose=False)
    assert (tmp_path / "names.txt").read_text().split() == ["chr1"]
    assert names["status"] == "completed"


@needs_hisat2
def test_align_spliced_reads_and_read_statistics(tmp_path, genome):
    fasta, gtf, sequence = genome
    subprocess.run(["hisat2-build", str(fasta), str(tmp_path / "idx")], check=True, capture_output=True)
    spliced = sequence[350:400] + sequence[600:650]
    plain = sequence[1500:1600]
    fastq = tmp_path / "reads.fq"
    fastq.write_text(f"@spliced\n{spliced}\n+\n{'I' * 100}\n@plain\n{plain}\n+\n{'I' * 100}\n")
    result = hs.hisat2_align_reads(
        str(tmp_path / "idx"), str(tmp_path / "a.sam"), unpaired_reads=[str(fastq)], summary_file=str(tmp_path / "sum.txt"),
        novel_splice_sites_out=str(tmp_path / "novel.txt"), verbose=False,
    )
    assert result["status"] == "completed"
    rows = {line.split("\t")[0]: line.split("\t") for line in (tmp_path / "a.sam").read_text().splitlines() if not line.startswith("@")}
    assert rows["plain"][3] == "1501" and rows["plain"][5] == "100M"
    assert rows["spliced"][5] == "50M200N50M"
    stats = hs.hisat2_read_statistics(str(fastq), verbose=False)
    assert (stats["reads_examined"], stats["minimum_length"], stats["maximum_length"]) == (2, 100, 100)


@needs_hisat2
def test_extract_variants_and_simulate_dna_reads(tmp_path, genome):
    fasta, gtf, sequence = genome
    vcf = tmp_path / "v.vcf"
    rows = []
    for position in (200, 500, 900):
        reference = sequence[position - 1]
        alternative = {"A": "G", "C": "T", "G": "A", "T": "C"}[reference]
        rows.append(f"chr1\t{position}\trs{position}\t{reference}\t{alternative}\t.\t.\t.\n")
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "".join(rows))
    extracted = hs.hisat2_extract_snps_haplotypes_vcf(str(fasta), [str(vcf)], str(tmp_path / "snps"), verbose=False)
    assert extracted["status"] == "completed"
    assert len((tmp_path / "snps.snp").read_text().splitlines()) == 3
    # the DNA simulator keeps only chromosomes named 1-22, X, Y, so rename the chromosome
    genome_numeric = tmp_path / "g_num.fa"
    genome_numeric.write_text(fasta.read_text().replace(">chr1", ">1"))
    gtf_numeric = tmp_path / "g_num.gtf"
    gtf_numeric.write_text(GTF.replace("chr1", "1"))
    snp_numeric = tmp_path / "snps_num.snp"
    snp_numeric.write_text((tmp_path / "snps.snp").read_text().replace("chr1", "1"))
    simulated = hs.hisat2_simulate_reads(str(genome_numeric), str(gtf_numeric), str(snp_numeric), str(tmp_path / "sim"), dna=True, single_end=True, num_fragments=50, read_length=40, random_seed=1, verbose=False)
    assert simulated["status"] == "completed" and simulated["read_type"] == "dna"
    assert any(path.endswith("sim_1.fa") or path.endswith("sim.sam") for path in simulated["outputs"])


@pytest.mark.skipif(shutil.which("hisat2-repeat") is None, reason="hisat2-repeat is not installed (not in the bioconda package)")
def test_build_repeat_index(tmp_path):
    rng = random.Random(1)
    repeat = "".join(rng.choice("ACGT") for _ in range(300))
    pieces = []
    for _ in range(8):
        pieces.append("".join(rng.choice("ACGT") for _ in range(500)))
        pieces.append(repeat)
    sequence = "".join(pieces)
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">r1\n" + "\n".join(sequence[i:i + 70] for i in range(0, len(sequence), 70)) + "\n")
    result = hs.hisat2_build_repeat_index(str(fasta), str(tmp_path / "rep"), seed_length=30, verbose=False)
    assert result["status"] == "completed"
    assert any(path.endswith(".rep.fa") for path in result["outputs"])


def test_stand_in_programs_receive_the_expected_arguments(tmp_path):
    """hisat2-repeat and the simulator are called with the documented option names (stand-in executables)."""
    folder = tmp_path / "bin"
    folder.mkdir()
    for name in ("hisat2-repeat", "hisat2_simulate_reads.py"):
        script = folder / name
        script.write_text('#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a" >> "$0.log"; done\nfor last; do :; done\n: > "$last.rep.fa"\n: > "$last.sam"\n')
        script.chmod(0o755)
    reference = tmp_path / "r.fa"
    reference.write_text(">r\nACGT\n")
    result = hs.hisat2_build_repeat_index(str(reference), str(tmp_path / "out"), seed_length=40, repeat_count=3, forward_only=True, hisat2_path=str(folder), verbose=False)
    arguments = (folder / "hisat2-repeat.log").read_text().split("\n")
    arguments = arguments[arguments.index("--version") + 1:] if "--version" in arguments else arguments
    assert result["status"] == "completed" and arguments[:5] == ["--seed-length", "40", "--repeat-count", "3", "--forward-only"]
    gtf = tmp_path / "g.gtf"
    gtf.write_text(GTF)
    snp = tmp_path / "s.snp"
    snp.write_text("")
    simulated = hs.hisat2_simulate_reads(str(reference), str(gtf), str(snp), str(tmp_path / "sim"), dna=True, num_fragments=5, hisat2_path=str(folder), verbose=False)
    assert simulated["status"] == "completed"
    lines = (folder / "hisat2_simulate_reads.py.log").read_text().split("\n")
    lines = lines[lines.index("--version") + 1:] if "--version" in lines else lines
    assert lines[:3] == ["--dna", "--num-fragment", "5"]
