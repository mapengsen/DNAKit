"""Tests for the SPAdes toolkit facades; program-dependent tests skip when SPAdes is not installed."""

from __future__ import annotations

import random
import shutil

import pytest

from dnakit.assembly import spades_p14_cli_functions as sp
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

HAS = shutil.which("spades.py") is not None and shutil.which("spades-gbuilder") is not None
needs_spades = pytest.mark.skipif(not HAS, reason="SPAdes is not installed")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    folder = tmp_path_factory.mktemp("spades")
    rng = random.Random(7)
    repeat = "".join(rng.choice("ACGT") for _ in range(400))
    genome = (
        "".join(rng.choice("ACGT") for _ in range(6000))
        + repeat
        + "".join(rng.choice("ACGT") for _ in range(6000))
        + repeat
        + "".join(rng.choice("ACGT") for _ in range(6000))
    )
    (folder / "genome.fa").write_text(">g\n" + genome + "\n")
    first, second = [], []
    for index in range(5000):
        start = rng.randrange(0, len(genome) - 300)
        fragment = genome[start : start + 300]
        a, b = fragment[:100], fragment[::-1].translate(COMPLEMENT)[:100]
        if rng.random() < 0.3:
            position = rng.randrange(100)
            a = a[:position] + rng.choice("ACGT") + a[position + 1 :]
        first.append(f"@p{index}/1\n{a}\n+\n{'I' * 100}\n")
        second.append(f"@p{index}/2\n{b}\n+\n{'I' * 100}\n")
    (folder / "r1.fq").write_text("".join(first))
    (folder / "r2.fq").write_text("".join(second))
    long_reads = "".join(
        f">l{i}\n{genome[start : start + 2000]}\n"
        for i, start in enumerate((100, 5000, 12000, 16000))
    )
    (folder / "long.fa").write_text(long_reads)
    (folder / "q.fa").write_text(">q\n" + genome[3000:3200] + "\n")
    return folder, genome


def test_missing_program_and_validation(tmp_path):
    reads = tmp_path / "r.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    with pytest.raises(BackendUnavailableError):
        sp.spades_estimate_kmer_count(21, unpaired_reads=[str(reads)], spades_path=str(tmp_path))
    with pytest.raises(ConfigurationError):
        sp.spades_estimate_kmer_count(21)
    with pytest.raises(ConfigurationError):
        sp.spades_estimate_kmer_count(21, reads_1=[str(reads)], reads_2=[])
    with pytest.raises(ConfigurationError):
        sp.spades_build_assembly_graph(
            str(tmp_path / "g"),
            kmer=21,
            fasta=str(reads),
            reads_1=[str(reads)],
            reads_2=[str(reads)],
        )
    with pytest.raises(ConfigurationError):
        sp.spades_build_assembly_graph(
            str(tmp_path / "g"), kmer=21, fasta=str(reads), graph_format="bad"
        )


@needs_spades
def test_kmer_estimate_filter_and_graph_pipeline(data, tmp_path):
    folder, _genome = data
    r1, r2 = [str(folder / "r1.fq")], [str(folder / "r2.fq")]
    estimate = sp.spades_estimate_kmer_count(21, reads_1=r1, reads_2=r2, verbose=False)
    assert estimate["status"] == "completed" and 17_000 < estimate["distinct_kmers"] < 80_000
    filtered = sp.spades_filter_reads_by_kmer_coverage(
        str(tmp_path / "filtered"),
        kmer=21,
        coverage_threshold=3,
        reads_1=r1,
        reads_2=r2,
        verbose=False,
    )
    assert filtered["status"] == "completed" and any(
        path.endswith("dataset.yaml") for path in filtered["outputs"]
    )
    graph = sp.spades_build_assembly_graph(
        str(tmp_path / "genome.gfa"),
        kmer=21,
        fasta=str(folder / "genome.fa"),
        graph_format="gfa",
        verbose=False,
    )
    assert graph["status"] == "completed"
    text = (tmp_path / "genome.gfa").read_text()
    assert text.count("\nS\t") + text.startswith("S\t") >= 3
    from_reads = sp.spades_build_assembly_graph(
        str(tmp_path / "reads.gfa"),
        kmer=21,
        reads_1=r1,
        reads_2=r2,
        graph_format="gfa",
        infer_coverage=True,
        verbose=False,
    )
    assert from_reads["status"] == "completed"
    simplified = sp.spades_simplify_assembly_graph(
        str(tmp_path / "genome.gfa"),
        str(tmp_path / "simp"),
        kmer=21,
        read_length=100,
        verbose=False,
    )
    assert simplified["status"] == "completed" and (tmp_path / "simp.gfa").is_file()
    split = sp.spades_split_assembly_graph(
        str(tmp_path / "genome.gfa"), str(tmp_path / "components"), verbose=False
    )
    assert split["components"] == 1
    mapped = sp.spades_map_reads_to_graph(
        str(tmp_path / "genome.gfa"),
        str(tmp_path / "map.txt"),
        kmer=21,
        reads_1=r1,
        reads_2=r2,
        verbose=False,
    )
    assert mapped["status"] == "completed" and "weight" in (tmp_path / "map.txt").read_text()
    aligned = sp.spaligner_align_long_reads(
        str(folder / "long.fa"),
        str(tmp_path / "genome.gfa"),
        str(tmp_path / "spaligner"),
        datatype="nanopore",
        kmer=21,
        verbose=False,
    )
    assert aligned["status"] == "completed"
    assert len((tmp_path / "spaligner" / "alignment.tsv").read_text().splitlines()) == 4
    pathracer = sp.pathracer_search_graph(
        str(folder / "q.fa"),
        str(tmp_path / "genome.gfa"),
        str(tmp_path / "pathracer"),
        query_type="nt",
        threads=2,
        verbose=False,
    )
    assert pathracer["status"] == "completed" and any(
        path.endswith(".seqs.fa") for path in pathracer["outputs"]
    )


@needs_spades
def test_error_correction(data, tmp_path):
    folder, _ = data
    result = sp.spades_correct_reads(
        str(tmp_path / "ec"),
        reads_1=[str(folder / "r1.fq")],
        reads_2=[str(folder / "r2.fq")],
        threads=2,
        phred_offset=33,
        verbose=False,
    )
    assert result["status"] == "completed" and result["corrector"] == "bayeshammer"
    assert any(path.endswith(".cor.fastq.gz") for path in result["corrected_files"])


def test_stand_in_programs_receive_the_expected_arguments(tmp_path):
    """Programs that need special data (SPlitteR, BinSPreader, PathRacer seq-fs) are checked with stand-ins."""
    folder = tmp_path / "bin"
    folder.mkdir()
    for name in ("splitter", "binspreader", "pathracer-seq-fs"):
        script = folder / name
        script.write_text(
            '#!/bin/sh\nfor a in "$@"; do printf "%s\\n" "$a" >> "$0.log"; done\nfor a in "$@"; do case "$a" in *.tsv) [ -e "$a" ] || : > "$a";; esac; done\n'
        )
        script.chmod(0o755)
    graph = tmp_path / "g.gfa"
    graph.write_text("S\t1\tACGT\n")
    yaml = tmp_path / "slr.yaml"
    yaml.write_text("[]\n")
    result = sp.splitter_resolve_repeats(
        str(graph),
        str(yaml),
        str(tmp_path / "split_out"),
        resolution_mode="meta",
        graph_type="mdbg",
        score_threshold=2,
        threads=3,
        spades_path=str(folder),
        verbose=False,
    )
    assert result["status"] == "completed"
    lines = (folder / "splitter.log").read_text().split("\n")
    lines = lines[lines.index("--no-such-flag") + 1 :] if "--no-such-flag" in lines else lines
    assert (
        "-Mmeta" in lines
        and "-Gmdbg" in lines
        and lines[lines.index("--score") + 1] == "2.0"
        and lines[lines.index("-t") + 1] == "3"
    )
    binning = tmp_path / "bins.tsv"
    binning.write_text("c1\tb1\n")
    refined = sp.binspreader_refine_bins(
        str(graph),
        str(binning),
        str(tmp_path / "out.tsv"),
        assignment_strategy="mle",
        refiner="prop",
        multiple_assignment=True,
        spades_path=str(folder),
        verbose=False,
    )
    assert refined["status"] == "completed"
    arguments = (folder / "binspreader.log").read_text().split("\n")
    assert "-Smle" in arguments and "-Rprop" in arguments and "-m" in arguments
    with pytest.raises(ConfigurationError):
        sp.binspreader_refine_bins(
            str(graph),
            str(binning),
            str(tmp_path / "o2.tsv"),
            split_reads=True,
            spades_path=str(folder),
        )
    hmm = tmp_path / "x.hmm"
    hmm.write_text("HMMER3/f\n")
    sequences = tmp_path / "s.fa"
    sequences.write_text(">s\nACGT\n")
    fs = sp.pathracer_align_hmm_to_sequences(
        str(hmm),
        str(sequences),
        str(tmp_path / "fs_out"),
        max_frameshifts=4,
        exhaustive=True,
        spades_path=str(folder),
        verbose=False,
    )
    fs_lines = (folder / "pathracer-seq-fs.log").read_text().split("\n")
    assert fs["status"] == "completed" and "--max-fs" in fs_lines and "--exhaustive" in fs_lines
