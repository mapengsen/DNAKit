"""Tests for dnakit.microbiome.khmer_p27_label_read_functions (labelled graphs, HLL merge, aligner, read files).

Offline: inputs are small random sequences written into tmp_path. Skipped when khmer is missing.
"""

from __future__ import annotations

import inspect
import random
from pathlib import Path

import pytest

khmer = pytest.importorskip("khmer")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.microbiome import khmer_p27_label_read_functions as m  # noqa: E402

V = {"verbose": False}
K = 15


def _dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture()
def contigs() -> list[str]:
    """Two contigs sharing a 40-base stem, then diverging: one branch point."""
    stem = _dna(40, 1)
    return [stem + _dna(40, 2), stem + _dna(40, 3)]


@pytest.fixture()
def labelled_reads(tmp_path: Path, contigs: list[str]) -> Path:
    path = tmp_path / "lab.fa"
    path.write_text(f">c1\n{contigs[0]}\n>c2\n{contigs[1]}\n")
    return path


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(m.__all__) == 11
    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_label_reads_files_strings_save_and_load(
    labelled_reads: Path, contigs: list[str], tmp_path: Path
) -> None:
    from_file = m.khmer_label_reads(reads_path=labelled_reads, ksize=K, **V)
    assert from_file["n_labels"] == 2
    assert {tuple(v) for v in from_file["tag_labels"].values()} <= {(0,), (1,), (0, 1)}
    strings = m.khmer_label_reads(
        sequences=contigs, labels=[3, 4], ksize=K, save_labels_path=tmp_path / "lt.bin", **V
    )
    assert strings["n_labels"] == 2
    assert {label for v in strings["tag_labels"].values() for label in v} == {3, 4}
    loaded = m.khmer_label_reads(load_labels_path=tmp_path / "lt.bin", ksize=K, **V)
    assert loaded["n_labels"] == 2
    with pytest.raises(ConfigurationError):
        m.khmer_label_reads(ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_label_reads(sequences=contigs, labels=[1], ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_label_reads(reads_path=labelled_reads, labels=[1, 2], ksize=K, **V)


def test_sweep_neighborhood_labels_and_tags_and_tag_labels(
    labelled_reads: Path, contigs: list[str]
) -> None:
    labels = m.khmer_sweep_neighborhood(
        contigs[0][:30], reads_path=labelled_reads, sweep_range=40, ksize=K, **V
    )
    assert labels["labels"] == [0, 1]
    tags = m.khmer_sweep_neighborhood(
        contigs[0][:30], output="tags", reads_path=labelled_reads, sweep_range=40, ksize=K, **V
    )
    assert tags["tags"]
    tag = tags["tags"][0]
    found = m.khmer_tag_labels(tag, reads_path=labelled_reads, ksize=K, **V)
    assert found["labels"] in ([0], [1], [0, 1])
    assert m.khmer_tag_labels("A" * K, reads_path=labelled_reads, ksize=K, **V)["labels"] == []
    with pytest.raises(ConfigurationError):
        m.khmer_sweep_neighborhood(contigs[0], output="x", reads_path=labelled_reads, ksize=K, **V)


def test_label_partitioned_reads(tmp_path: Path) -> None:
    part = tmp_path / "p.part"
    part.write_text(f">a\t5\n{_dna(60, 1)}\n>b\t9\n{_dna(60, 2)}\n")
    result = m.khmer_label_partitioned_reads(part, ksize=K, **V)
    assert result["n_labels"] == 2
    assert {label for v in result["tag_labels"].values() for label in v} == {5, 9}


def test_label_across_high_degree_nodes_and_assembly(contigs: list[str]) -> None:
    result = m.khmer_label_across_high_degree_nodes(contigs, ksize=K, **V)
    assert result["n_high_degree"] == 1
    assert next(iter(result["high_degree_labels"].values())) == [1, 2]
    assembled = m.khmer_assemble_labeled_path(contigs[0][:K], contigs, ksize=K, **V)
    assert assembled["n_paths"] >= 1
    assert all(len(path) >= K for path in assembled["paths"])
    with pytest.raises(ConfigurationError):
        m.khmer_label_across_high_degree_nodes([], ksize=K, **V)


def test_merge_hll_counters_estimates() -> None:
    first, second = _dna(300, 1), _dna(300, 2)
    result = m.khmer_merge_hll_counters(
        first_sequences=[first], second_sequences=[second], ksize=K, **V
    )
    expected = 2 * (300 - K + 1)
    assert result["first_estimate"] == pytest.approx(300 - K + 1, rel=0.1)
    assert result["merged_estimate"] == pytest.approx(expected, rel=0.1)
    same = m.khmer_merge_hll_counters(
        first_sequences=[first], second_sequences=[first], ksize=K, **V
    )
    assert same["merged_estimate"] == same["first_estimate"]
    with pytest.raises(ConfigurationError):
        m.khmer_merge_hll_counters(first_sequences=[first], ksize=K, **V)


@pytest.mark.parametrize("forward_only", [False, True])
def test_align_read_matches_khmer_and_corrects_an_error(forward_only: bool) -> None:
    genome = _dna(100, 8)
    read = genome[10:50]
    mutated = read[:20] + ("A" if read[20] != "A" else "C") + read[21:]
    result = m.khmer_align_read(
        mutated, sequences=[genome] * 5, ksize=K, forward_only=forward_only, **V
    )
    assert result["graph_alignment"] == read
    graph = khmer.Countgraph(K, 1_000_000, 4)
    for _ in range(5):
        graph.consume(genome)
    aligner = khmer.ReadAligner(graph, 2, 1.0)
    outcome = aligner.align_forward(mutated) if forward_only else aligner.align(mutated)
    assert result["score"] == pytest.approx(outcome[0])
    with pytest.raises(ConfigurationError):
        m.khmer_align_read(read, ksize=K, **V)


def test_extract_and_split_paired_reads(tmp_path: Path) -> None:
    mixed = tmp_path / "mixed.fa"
    mixed.write_text(
        ">p1/1\nACGTACGTAC\n>p1/2\nGGGGCCCCAA\n>s1\nTTTTTTTTTT\n>p2/1\nAAAACCCCGG\n>p2/2\nCCCCAAAATT\n"
    )
    result = m.khmer_extract_paired_reads(mixed, tmp_path / "pe.fa", tmp_path / "se.fa", **V)
    assert (result["n_pairs"], result["n_singletons"]) == (2, 1)
    assert (tmp_path / "pe.fa").read_text().count(">") == 4
    assert (tmp_path / "se.fa").read_text() == ">s1\nTTTTTTTTTT\n"
    inter = tmp_path / "inter.fa"
    inter.write_text(">p1/1\nACGTACGTAC\n>p1/2\nGGGGCCCCAA\n>p2/1\nAAAACCCCGG\n>p2/2\nCCCCAAAATT\n")
    split = m.khmer_split_paired_reads(inter, tmp_path / "l.fa", tmp_path / "r.fa", **V)
    assert split["n_pairs"] == 2
    assert (tmp_path / "l.fa").read_text() == ">p1/1\nACGTACGTAC\n>p2/1\nAAAACCCCGG\n"
    assert (tmp_path / "r.fa").read_text() == ">p1/2\nGGGGCCCCAA\n>p2/2\nCCCCAAAATT\n"


def test_split_paired_reads_orphans_and_errors(tmp_path: Path) -> None:
    mixed = tmp_path / "mixed.fa"
    mixed.write_text(">p1/1\nACGTACGTAC\n>p1/2\nGGGGCCCCAA\n>s1\nTTTTTTTTTT\n")
    with pytest.raises(ConfigurationError):
        m.khmer_split_paired_reads(mixed, tmp_path / "l.fa", tmp_path / "r.fa", **V)
    assert not (tmp_path / "l.fa").exists()
    result = m.khmer_split_paired_reads(
        mixed, tmp_path / "l.fa", tmp_path / "r.fa", orphan_output_path=tmp_path / "o.fa", **V
    )
    assert result["n_orphans"] == 1 and (tmp_path / "o.fa").read_text() == ">s1\nTTTTTTTTTT\n"
    singles = tmp_path / "singles.fa"
    singles.write_text(">a\nACGT\n>b\nACGT\n")
    with pytest.raises(ConfigurationError):
        m.khmer_extract_paired_reads(singles, tmp_path / "pe2.fa", tmp_path / "se2.fa", **V)
    assert not (tmp_path / "pe2.fa").exists()
    with pytest.raises(ConfigurationError):
        m.khmer_split_paired_reads(mixed, tmp_path / "x.fa", tmp_path / "x.fa", **V)


def test_read_statistics_fasta_fastq_and_empty(tmp_path: Path) -> None:
    fasta = tmp_path / "a.fa"
    fasta.write_text(">r1\nACGTACGT\n>r2\nACGT\n")
    fastq = tmp_path / "b.fq"
    fastq.write_text("@q1\nACGTAC\n+\nIIIIII\n")
    empty = tmp_path / "e.fa"
    empty.write_text("")
    result = m.khmer_read_statistics([fasta, fastq, empty], output_csv=tmp_path / "s.csv", **V)
    assert [(r["bp"], r["seqs"]) for r in result["files"]] == [(12, 2), (6, 1), (0, 0)]
    assert result["total"] == {"bp": 18, "seqs": 3, "avg_len": 6.0}
    lines = (tmp_path / "s.csv").read_text().splitlines()
    assert lines[0] == "bp,seqs,avg_len,filename"
    assert lines[1].startswith("12,2,6.0,") and len(lines) == 3
    with pytest.raises(ConfigurationError):
        m.khmer_read_statistics([], **V)
