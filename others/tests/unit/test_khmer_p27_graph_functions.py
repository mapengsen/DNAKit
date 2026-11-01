"""Tests for dnakit.microbiome.khmer_p27_graph_functions (khmer k-mer queries on countgraphs and nodegraphs).

Offline: inputs are small random sequences written into tmp_path. Skipped when khmer is missing.
"""

from __future__ import annotations

import inspect
import random
from pathlib import Path

import pytest

khmer = pytest.importorskip("khmer")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.microbiome import khmer_p27_graph_functions as m  # noqa: E402

V = {"verbose": False}
K = 15


def _dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture()
def reads() -> list[str]:
    a = _dna(120, 1)
    return [a, _dna(120, 2), a[:60] + _dna(60, 3)]


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(m.__all__) == 10
    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_kmer_hashes_match_khmer_and_set_mode() -> None:
    seq = _dna(40, 5)
    graph = khmer.Countgraph(K, 1000, 1)
    result = m.khmer_kmer_hashes(seq, ksize=K, **V)
    assert result["hashes"] == list(graph.get_kmer_hashes(seq))
    assert result["n_hashes"] == 40 - K + 1
    as_set = m.khmer_kmer_hashes(seq + seq[:20], ksize=K, as_set=True, **V)
    assert as_set["hashes"] == sorted(set(as_set["hashes"]))
    assert m.khmer_kmer_hashes(seq.lower(), ksize=K, **V)["hashes"] == result["hashes"]


def test_kmer_hashes_rejects_short_and_non_dna() -> None:
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_hashes("ACGT", ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_hashes("N" * 30, ksize=K, **V)


def test_sequence_abundance_statistics_and_csv(reads: list[str], tmp_path: Path) -> None:
    csv_path = tmp_path / "median.csv"
    result = m.khmer_sequence_abundance(
        sequences_to_query=[reads[0]],
        sequences=reads,
        ksize=K,
        median_threshold=2,
        include_counts=True,
        output_csv=csv_path,
        **V,
    )
    row = result["results"][0]
    graph = khmer.Countgraph(K, 1_000_000, 4)
    for read in reads:
        graph.consume(read)
    assert row["counts"] == list(graph.get_kmer_counts(reads[0]))
    assert (row["min"], row["max"]) == (1, 2)
    assert row["median"] == graph.get_median_count(reads[0])[0]
    assert row["median_at_least"] is (graph.get_median_count(reads[0])[0] >= 2)
    assert csv_path.read_text().splitlines()[0] == "name,median,average,stddev,seqlen"
    with pytest.raises(ConfigurationError):
        m.khmer_sequence_abundance(sequences=reads, ksize=K, **V)


def test_sequence_abundance_from_saved_countgraph_and_query_file(
    reads: list[str], tmp_path: Path
) -> None:
    graph = khmer.Countgraph(K, 100_000, 4)
    for read in reads:
        graph.consume(read)
    saved = tmp_path / "cg"
    graph.save(str(saved))
    query = tmp_path / "q.fa"
    query.write_text(f">q1\n{reads[1]}\n>short\nACGT\n>q2\n{reads[0][:30]}N{reads[0][31:60]}\n")
    result = m.khmer_sequence_abundance(query_path=query, graph_path=saved, **V)
    assert result["n_queries"] == 2
    assert result["skipped_shorter_than_k"] == 1
    assert [r["name"] for r in result["results"]] == ["q1", "q2"]


def test_trim_by_abundance_modes(reads: list[str]) -> None:
    low = m.khmer_trim_by_abundance(
        reads[0], 2, mode="trim_on_abundance", sequences=reads, ksize=K, **V
    )
    assert low["trim_position"] == 60
    assert low["trimmed_sequence"] == reads[0][: low["trim_position"]]
    high = m.khmer_trim_by_abundance(reads[0], 1, sequences=reads, ksize=K, **V)
    assert high["trim_position"] == 0
    with pytest.raises(ConfigurationError):
        m.khmer_trim_by_abundance(reads[0], 1, mode="other", sequences=reads, ksize=K, **V)


def test_spectral_error_positions_flag_a_substitution() -> None:
    truth = _dna(100, 8)
    error = truth[:50] + ("A" if truth[50] != "A" else "C") + truth[51:]
    result = m.khmer_spectral_error_positions(error, sequences=[truth] * 5 + [error], ksize=K, **V)
    graph = khmer.Countgraph(K, 1_000_000, 4)
    for text in [truth] * 5 + [error]:
        graph.consume(text)
    assert result["positions"] == [int(p) for p in graph.find_spectral_error_positions(error, 1)]
    assert result["positions"]


def test_merge_nodegraphs_is_union(tmp_path: Path) -> None:
    one = khmer.Nodegraph(K, 100_000, 4)
    two = khmer.Nodegraph(K, 100_000, 4)
    seqs = (_dna(60, 1), _dna(60, 2))
    one.consume(seqs[0])
    two.consume(seqs[1])
    one.save(str(tmp_path / "a"))
    two.save(str(tmp_path / "b"))
    result = m.khmer_merge_nodegraphs(tmp_path / "a", tmp_path / "b", tmp_path / "u", **V)
    merged = khmer.load_nodegraph(str(tmp_path / "u"))
    assert merged.get(seqs[0][:K]) and merged.get(seqs[1][:K])
    assert (
        result["unique_kmers_merged"]
        == result["unique_kmers_first"] + result["unique_kmers_second"]
    )
    with pytest.raises(ConfigurationError):
        m.khmer_merge_nodegraphs(tmp_path / "a", tmp_path / "b", tmp_path / "u", **V)


def test_graph_summary_for_graph_and_table(tmp_path: Path) -> None:
    graph = khmer.Countgraph(K, 100_000, 4)
    graph.consume(_dna(60, 1))
    graph.save(str(tmp_path / "cg"))
    result = m.khmer_graph_summary(tmp_path / "cg", **V)
    assert result["ksize"] == K
    assert abs(result["estimated_unique_kmers"] - 46) <= 1
    table = khmer.Counttable(K, 100_000, 4)
    table.consume(_dna(60, 1))
    table.save(str(tmp_path / "ct"))
    assert m.khmer_graph_summary(tmp_path / "ct", graph_type="counttable", **V)["ksize"] == K
    with pytest.raises(ConfigurationError):
        m.khmer_graph_summary(tmp_path / "cg", graph_type="nodegraph", **V)


def test_neighbors_degree_and_component(reads: list[str]) -> None:
    graph = khmer.Nodegraph(K, 1_000_000, 4)
    for read in reads:
        graph.consume(read)
    start = reads[1][:K]
    result = m.khmer_kmer_neighbors(start, sequences=reads, ksize=K, **V)
    assert result["degree"] == graph.kmer_degree(start) == 1
    assert [graph.hash(n) for n in result["neighbors"]] == list(graph.neighbors(start))
    size = m.khmer_connected_component_size(
        start, sequences=reads, ksize=K, count_within_radius=5, **V
    )
    assert size["component_size"] == graph.calc_connected_graph_size(start, 0)
    assert size["kmers_within_radius"] == graph.count_kmers_within_radius(start, 5, 0)
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_neighbors(start[:-1], sequences=reads, ksize=K, **V)


def test_high_degree_kmers_and_linear_path(reads: list[str]) -> None:
    branch = m.khmer_high_degree_kmers(reads[2], sequences=reads, ksize=K, **V)
    assert branch["n_high_degree"] == 1
    isolated = reads[1][:K]
    path = m.khmer_linear_path(isolated, sequences=reads, ksize=K, **V)
    assert path["path_sequence"] == reads[1]
    kmers = m.khmer_linear_path(isolated, output="kmers", sequences=reads, ksize=K, **V)
    assert kmers["path_size"] == len(reads[1]) - K + 1
    assert kmers["branch_kmers"] == []
    stopped = m.khmer_linear_path(
        reads[0][:K], output="kmers", stop_sequence=reads[2], sequences=reads, ksize=K, **V
    )
    assert stopped["branch_kmers"]


def test_graph_input_must_be_exactly_one_source(reads: list[str], tmp_path: Path) -> None:
    graph = khmer.Nodegraph(K, 1000, 1)
    graph.save(str(tmp_path / "ng"))
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_neighbors(reads[0][:K], ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_neighbors(reads[0][:K], graph_path=tmp_path / "ng", sequences=reads, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_kmer_neighbors(reads[0][:K], graph_path=tmp_path / "missing", **V)
