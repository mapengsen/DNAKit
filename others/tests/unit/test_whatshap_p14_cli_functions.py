"""Tests for the WhatsHap facades; program-dependent tests skip when WhatsHap or pysam is missing."""

from __future__ import annotations

import random
import shutil
import subprocess

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import whatshap_p14_cli_functions as wh

pysam = pytest.importorskip("pysam")
needs_whatshap = pytest.mark.skipif(
    shutil.which("whatshap") is None, reason="WhatsHap is not installed"
)


@pytest.fixture(scope="module")
def phasing_data(tmp_path_factory):
    """A 4 kb reference, 25 heterozygous SNVs, and 500 error-free reads drawn from two haplotypes."""
    folder = tmp_path_factory.mktemp("whatshap")
    rng = random.Random(5)
    length = 4000
    reference = "".join(rng.choice("ACGT") for _ in range(length))
    (folder / "ref.fa").write_text(
        ">chr1\n" + "\n".join(reference[i : i + 60] for i in range(0, length, 60)) + "\n"
    )
    pysam.faidx(str(folder / "ref.fa"))
    sites = list(range(100, length - 200, 150))
    flip = {"A": "G", "C": "T", "G": "A", "T": "C"}
    second = list(reference)
    for position in sites:
        second[position] = flip[reference[position]]
    haplotypes = [reference, "".join(second)]
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": length}],
        "RG": [{"ID": "rg1", "SM": "S1"}],
    }
    with pysam.AlignmentFile(str(folder / "raw.bam"), "wb", header=header) as out:
        for index in range(500):
            start = rng.randrange(0, length - 300)
            read = pysam.AlignedSegment(out.header)
            read.query_name = f"r{index}"
            read.query_sequence = haplotypes[index % 2][start : start + 300]
            read.flag = 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = [(0, 300)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 300)
            read.set_tag("RG", "rg1")
            out.write(read)
    pysam.sort("-o", str(folder / "reads.bam"), str(folder / "raw.bam"))
    pysam.index(str(folder / "reads.bam"))
    with (folder / "vars.vcf").open("w") as handle:
        handle.write(
            f'##fileformat=VCFv4.2\n##contig=<ID=chr1,length={length}>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n'
        )
        for position in sites:
            handle.write(
                f"chr1\t{position + 1}\t.\t{reference[position]}\t{flip[reference[position]]}\t50\tPASS\t.\tGT\t0/1\n"
            )
    return folder, sites


@pytest.fixture(scope="module")
def phased(phasing_data, tmp_path_factory):
    if shutil.which("whatshap") is None:
        pytest.skip("WhatsHap is not installed")
    folder, _ = phasing_data
    out = tmp_path_factory.mktemp("phased")
    subprocess.run(
        [
            "whatshap",
            "phase",
            "-o",
            str(out / "ph.vcf"),
            "--reference",
            str(folder / "ref.fa"),
            str(folder / "vars.vcf"),
            str(folder / "reads.bam"),
        ],
        check=True,
        capture_output=True,
    )
    return out / "ph.vcf"


def test_missing_program_and_validation(tmp_path):
    vcf = tmp_path / "a.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    with pytest.raises(BackendUnavailableError):
        wh.whatshap_unphase_vcf(
            str(vcf), str(tmp_path / "o.vcf"), whatshap_path=str(tmp_path / "nope")
        )
    with pytest.raises(ConfigurationError):
        wh.whatshap_compare_phasing([str(vcf)], tsv_pairwise=str(tmp_path / "t"))
    with pytest.raises(ConfigurationError):
        wh.whatshap_compare_phasing([str(vcf), str(vcf)])
    with pytest.raises(ConfigurationError):
        wh.whatshap_split_reads(str(vcf), str(vcf))
    with pytest.raises(ConfigurationError):
        wh.whatshap_genotype_variants(
            str(vcf), [str(vcf)], str(tmp_path / "o.vcf"), recombination_rate=1e-8
        )
    with pytest.raises(ConfigurationError):
        wh.whatshap_polyphase(str(vcf), [str(vcf)], str(tmp_path / "o.vcf"), ploidy=1)


@needs_whatshap
def test_stats_unphase_compare(phased, tmp_path):
    stats = wh.whatshap_phasing_statistics(
        str(phased),
        str(tmp_path / "st.tsv"),
        gtf=str(tmp_path / "st.gtf"),
        block_list=str(tmp_path / "bl.txt"),
        verbose=False,
    )
    row = stats["statistics"][0]
    assert (row["variants"], row["phased"], row["blocks"]) == ("25", "25", "1")
    unphased = wh.whatshap_unphase_vcf(str(phased), str(tmp_path / "un.vcf"), verbose=False)
    assert unphased["status"] == "completed"
    genotypes = [
        line.split("\t")[9].split(":")[0]
        for line in (tmp_path / "un.vcf").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert set(genotypes) == {"0/1"}
    compared = wh.whatshap_compare_phasing(
        [str(phased), str(phased)],
        tsv_pairwise=str(tmp_path / "cmp.tsv"),
        names=["a", "b"],
        verbose=False,
    )
    assert compared["status"] == "completed" and compared["pairwise"][0]["all_switches"] == "0"


@needs_whatshap
def test_haplotag_split_haplotagphase_and_genotype(phasing_data, phased, tmp_path):
    folder, sites = phasing_data
    compressed = tmp_path / "ph.vcf"
    shutil.copy(phased, compressed)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    tagged = wh.whatshap_haplotag_reads(
        str(compressed) + ".gz",
        str(folder / "reads.bam"),
        str(tmp_path / "tag.bam"),
        reference=str(folder / "ref.fa"),
        haplotag_list=str(tmp_path / "tags.tsv"),
        verbose=False,
    )
    assert tagged["status"] == "completed"
    haplotypes = {}
    for line in (tmp_path / "tags.tsv").read_text().splitlines()[1:]:
        fields = line.split("\t")
        if fields[1] != "none":
            haplotypes[fields[0]] = fields[1]
    assert len(haplotypes) > 400
    # reads alternate between the two haplotypes: even and odd read numbers must carry different tags
    assert len({haplotypes[name] for name in haplotypes if int(name[1:]) % 2 == 0}) == 1
    assert {haplotypes[name] for name in haplotypes if int(name[1:]) % 2 == 0} != {
        haplotypes[name] for name in haplotypes if int(name[1:]) % 2 == 1
    }
    split = wh.whatshap_split_reads(
        str(folder / "reads.bam"),
        str(tmp_path / "tags.tsv"),
        output_h1=str(tmp_path / "h1.bam"),
        output_h2=str(tmp_path / "h2.bam"),
        verbose=False,
    )
    assert split["status"] == "completed"
    assert sum(1 for _ in pysam.AlignmentFile(str(tmp_path / "h1.bam"))) + sum(
        1 for _ in pysam.AlignmentFile(str(tmp_path / "h2.bam"))
    ) == len(haplotypes)
    pysam.index(str(tmp_path / "tag.bam"))
    unphased = wh.whatshap_unphase_vcf(str(phased), str(tmp_path / "un.vcf"), verbose=False)
    rephased = wh.whatshap_phase_from_haplotags(
        str(tmp_path / "un.vcf"),
        str(tmp_path / "tag.bam"),
        str(tmp_path / "hp.vcf"),
        reference=str(folder / "ref.fa"),
        verbose=False,
    )
    assert unphased["status"] == rephased["status"] == "completed"
    assert sum(
        "|" in line.split("\t")[9]
        for line in (tmp_path / "hp.vcf").read_text().splitlines()
        if not line.startswith("#")
    ) == len(sites)
    genotyped = wh.whatshap_genotype_variants(
        str(folder / "vars.vcf"),
        [str(folder / "reads.bam")],
        str(tmp_path / "gt.vcf"),
        reference=str(folder / "ref.fa"),
        verbose=False,
    )
    assert genotyped["status"] == "completed"
    calls = [
        line.split("\t")[9].split(":")[0]
        for line in (tmp_path / "gt.vcf").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert set(calls) <= {"0/1", "1/0"}


@needs_whatshap
def test_snv_candidates_learn_and_polyphase(phasing_data, phased, tmp_path):
    folder, sites = phasing_data
    found = wh.whatshap_find_snv_candidates(
        str(folder / "ref.fa"),
        str(folder / "reads.bam"),
        str(tmp_path / "snv.vcf"),
        technology="illumina",
        verbose=False,
    )
    assert found["status"] == "completed"
    positions = {
        int(line.split("\t")[1]) - 1
        for line in (tmp_path / "snv.vcf").read_text().splitlines()
        if not line.startswith("#")
    }
    assert set(sites) <= positions
    learned = wh.whatshap_learn_error_profile(
        str(folder / "reads.bam"),
        str(phased),
        str(folder / "ref.fa"),
        str(tmp_path / "learn.txt"),
        verbose=False,
    )
    assert learned["status"] == "completed"
    poly = wh.whatshap_polyphase(
        str(folder / "vars.vcf"),
        [str(folder / "reads.bam")],
        str(tmp_path / "poly.vcf"),
        ploidy=2,
        reference=str(folder / "ref.fa"),
        verbose=False,
    )
    assert poly["status"] == "completed"


def test_pedigree_and_progeny_programs_receive_the_expected_arguments(tmp_path):
    """Pedigree commands need family data; here a stand-in whatshap records the argument vector."""
    script = tmp_path / "whatshap"
    script.write_text(
        '#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a" >> "$0.log"; done\nwhile [ $# -gt 0 ]; do case "$1" in -o|--output) : > "$2";; esac; shift; done\n'
    )
    script.chmod(0o755)
    vcf = tmp_path / "family.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    bam = tmp_path / "a.bam"
    bam.write_bytes(b"")
    ped = tmp_path / "family.ped"
    ped.write_text("F1 child father mother 0 0\n")
    phased = wh.whatshap_phase_with_pedigree(
        str(vcf),
        [str(bam)],
        str(ped),
        str(tmp_path / "o.vcf"),
        recombination_rate=1e-8,
        no_genetic_haplotyping=True,
        whatshap_path=str(script),
        verbose=False,
    )
    assert phased["status"] == "completed"
    lines = (tmp_path / "whatshap.log").read_text().split("\n")
    lines = lines[lines.index("phase") :]
    assert (
        lines[0] == "phase"
        and "--ped" in lines
        and "--recombrate" in lines
        and "--no-genetic-haplotyping" in lines
    )
    (tmp_path / "whatshap.log").unlink()
    genotyped = wh.whatshap_genotype_variants(
        str(vcf),
        [str(bam)],
        str(tmp_path / "g.vcf"),
        pedigree=str(ped),
        use_ped_samples=True,
        whatshap_path=str(script),
        verbose=False,
    )
    assert genotyped["status"] == "completed" and "--ped" in (tmp_path / "whatshap.log").read_text()
    (tmp_path / "whatshap.log").unlink()
    progeny = wh.whatshap_polyphase_with_progeny(
        str(vcf),
        str(ped),
        str(tmp_path / "p.vcf"),
        ploidy=4,
        scoring_window=10,
        whatshap_path=str(script),
        verbose=False,
    )
    text = (tmp_path / "whatshap.log").read_text()
    assert (
        progeny["status"] == "completed"
        and "polyphasegenetic" in text
        and "--scoring-window" in text
    )
