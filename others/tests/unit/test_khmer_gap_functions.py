"""Tests for dnakit.microbiome.khmer_gap_functions (khmer table sizing, countgraph and nodegraph building).

Offline: inputs are small reads written into tmp_path, or the shared reads.fastq. Skipped when khmer is missing.
"""

from __future__ import annotations

import inspect
from collections import Counter
from pathlib import Path

import pytest

khmer = pytest.importorskip("khmer")

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import khmer_gap_functions as facade

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
READS_FASTQ = DATA / "reads.fastq"
FASTA = ">r1\nACGTACGTAC\n>r2\nACGTTTTTNN\n>r3\nAC\n>r4\nacgttgca\n"
K = 4


@pytest.fixture()
def reads(tmp_path: Path) -> Path:
    path = tmp_path / "reads.fa"
    path.write_text(FASTA, encoding="utf-8")
    return path


def _canonical(kmer: str) -> str:
    reverse = kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(kmer, reverse)


def _expected_counts(sequences: list[str], k: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for sequence in sequences:
        if len(sequence) < k or not set(sequence) <= set("ACGT"):
            continue
        for index in range(len(sequence) - k + 1):
            counts[_canonical(sequence[index : index + k])] += 1
    return counts


def test_public_functions_take_keyword_only_verbose() -> None:
    for name in facade.__all__:
        parameter = inspect.signature(getattr(facade, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, name
        assert parameter.default is True, name


def test_countgraph_table_size_from_memory_budget() -> None:
    result = facade.khmer_table_size_from_memory(
        "countgraph", max_memory_bytes=1_000_000, verbose=False
    )
    assert result["buckets_per_byte"] == 1
    assert result["table_size"] == 250_000
    assert result["estimated_memory_bytes"] == pytest.approx(1_000_000)


def test_nodegraph_table_size_from_memory_budget() -> None:
    result = facade.khmer_table_size_from_memory(
        "nodegraph", max_memory_bytes=1_000_000, verbose=False
    )
    assert result["buckets_per_byte"] == 8
    assert result["table_size"] == 2_000_000
    assert result["estimated_memory_bytes"] == pytest.approx(1_000_000)


def test_table_size_from_explicit_tablesize() -> None:
    result = facade.khmer_table_size_from_memory("countgraph", max_tablesize=5000, verbose=False)
    assert result["table_size"] == 5000
    assert result["estimated_memory_bytes"] == pytest.approx(20_000)


def test_table_size_needs_exactly_one_input() -> None:
    with pytest.raises(ConfigurationError):
        facade.khmer_table_size_from_memory("countgraph", verbose=False)
    with pytest.raises(ConfigurationError):
        facade.khmer_table_size_from_memory(
            "countgraph", max_memory_bytes=1000, max_tablesize=10, verbose=False
        )


def test_table_size_rejects_unknown_graph_type() -> None:
    with pytest.raises(ConfigurationError):
        facade.khmer_table_size_from_memory("bloom", max_tablesize=10, verbose=False)


def test_countgraph_counts_match_canonical_kmer_counts(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "cg.ct"
    result = facade.khmer_build_countgraph(
        str(reads), str(out), ksize=K, table_size=100_000, verbose=False
    )
    assert result["kind"] == "countgraph"
    assert result["reads"] == 4
    assert result["skipped_non_acgt_reads"] == 1  # r2 has N
    assert result["skipped_short_reads"] == 1  # r3 is shorter than k
    assert result["used_reads"] == 2
    table = khmer.load_countgraph(str(out))
    expected = _expected_counts(["ACGTACGTAC", "ACGTTGCA"], K)
    for kmer, count in expected.items():
        assert table.get(kmer) == count, kmer
    assert table.get("AAAA") == 0
    assert result["kmers_consumed"] == sum(expected.values())


def test_small_count_variant_is_reported(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "small.ct"
    result = facade.khmer_build_countgraph(
        str(reads), str(out), ksize=K, table_size=100_000, small_count=True, verbose=False
    )
    assert result["kind"] == "smallcountgraph"
    assert out.is_file()


def test_nodegraph_records_presence_of_used_kmers(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "ng.pt"
    result = facade.khmer_build_nodegraph(
        str(reads), str(out), ksize=K, table_size=100_000, verbose=False
    )
    assert result["kind"] == "nodegraph"
    table = khmer.load_nodegraph(str(out))
    for kmer in _expected_counts(["ACGTACGTAC", "ACGTTGCA"], K):
        assert table.get(kmer) >= 1, kmer
    assert table.get("GGGG") == 0  # absent k-mer (canonical CCCC) is not recorded


def test_build_on_shared_fastq_counts_every_read(tmp_path: Path) -> None:
    if not READS_FASTQ.is_file():
        pytest.skip("shared reads.fastq is not present")
    records = READS_FASTQ.read_text().count("\n") // 4
    out = tmp_path / "shared.pt"
    result = facade.khmer_build_nodegraph(
        str(READS_FASTQ), str(out), ksize=21, table_size=1_000_000, verbose=False
    )
    assert result["reads"] == records
    assert (
        result["used_reads"] + result["skipped_non_acgt_reads"] + result["skipped_short_reads"]
        == records
    )


def test_read_limit_is_enforced(reads: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.khmer_build_countgraph(
            str(reads), str(tmp_path / "x.ct"), ksize=K, table_size=1000, max_reads=1, verbose=False
        )


def test_build_refuses_existing_output(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "exists.ct"
    out.write_text("x")
    with pytest.raises(ConfigurationError):
        facade.khmer_build_countgraph(str(reads), str(out), ksize=K, table_size=1000, verbose=False)


def test_build_rejects_invalid_ksize(reads: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.khmer_build_nodegraph(str(reads), str(tmp_path / "x.pt"), ksize=33, verbose=False)
