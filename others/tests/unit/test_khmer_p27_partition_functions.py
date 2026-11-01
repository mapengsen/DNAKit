"""Tests for dnakit.microbiome.khmer_p27_partition_functions (tags, partitions, stop tags).

Offline: reads are small random sequences written into tmp_path. Skipped when khmer is missing.
"""

from __future__ import annotations

import inspect
import random
from pathlib import Path

import pytest

khmer = pytest.importorskip("khmer")

from dnakit.exceptions import ConfigurationError  # noqa: E402
from dnakit.microbiome import khmer_p27_partition_functions as m  # noqa: E402

V = {"verbose": False}
K = 15


def _dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture()
def reads(tmp_path: Path) -> Path:
    """Three reads: a and c share 60 bases (one partition), b is separate."""
    a = _dna(120, 1)
    path = tmp_path / "t.fa"
    path.write_text(f">a\n{a}\n>b\n{_dna(120, 2)}\n>c\n{a[:60]}{_dna(60, 3)}\n")
    return path


def _seqs(path: Path) -> list[str]:
    return path.read_text().split("\n")[1::2][:3]


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(m.__all__) == 10
    for name in m.__all__:
        parameter = inspect.signature(getattr(m, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_tag_reads_save_load_text_and_extra_tags(reads: Path, tmp_path: Path) -> None:
    result = m.khmer_tag_reads(
        reads_path=reads,
        ksize=K,
        extra_tag_kmers=["A" * K],
        save_tagset_path=tmp_path / "ts.bin",
        tagset_text_path=tmp_path / "ts.txt",
        **V,
    )
    graph = khmer.Nodegraph(K, 1_000_000, 4)
    graph.consume_seqfile_and_tag(str(reads))
    assert result["n_tags"] == graph.n_tags() + 1
    assert "A" * K in result["tags"]
    assert sorted((tmp_path / "ts.txt").read_text().split()) == result["tags"]
    loaded = m.khmer_tag_reads(load_tagset_path=tmp_path / "ts.bin", ksize=K, **V)
    assert loaded["n_tags"] == result["n_tags"]
    with pytest.raises(ConfigurationError):
        m.khmer_tag_reads(ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_tag_reads(reads_path=reads, ksize=K, extra_tag_kmers=["ACGT"], **V)


def test_tag_positions_reachable_tags_and_partition_id(reads: Path) -> None:
    a = _seqs(reads)[0]
    positions = m.khmer_tag_positions(a, reads_path=reads, ksize=K, **V)["tags"]
    graph = khmer.Nodegraph(K, 1_000_000, 4)
    graph.consume_seqfile_and_tag(str(reads))
    assert [(t["position"], graph.hash(t["kmer"])) for t in positions] == [
        (int(p), int(h)) for p, h in graph.get_tags_and_positions(a)
    ]
    tag = positions[0]["kmer"]
    tags = m.khmer_reachable_tags(tag, reads_path=reads, ksize=K, **V)
    assert tags["partition_id"] is None
    assert sorted(graph.hash(t) for t in tags["tags"]) == sorted(
        int(h) for h in graph.find_all_tags_list(tag)
    )
    assigned = m.khmer_reachable_tags(tag, reads_path=reads, ksize=K, assign_partition=True, **V)
    assert assigned["partition_id"] == graph.assign_partition_id(graph.find_all_tags(tag))


def test_unique_paths_empty_graph_and_known_graph(reads: Path) -> None:
    a = _seqs(reads)[0]
    assert m.khmer_unique_paths(a, ksize=K, **V)["paths"] == [a]
    assert m.khmer_unique_paths(a, reads_path=reads, ksize=K, **V)["paths"] == []
    with pytest.raises(ConfigurationError):
        m.khmer_unique_paths(a, min_length=K - 1, ksize=K, **V)


def test_partition_reads_annotates_like_khmer(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    result = m.khmer_partition_reads(reads, output_dir=out, ksize=K, subset_size=3, **V)
    assert result["partition_count"] == {"n_partitions": 2, "n_unassigned_tags": 0}
    assert result["n_subsets"] == 3
    graph = khmer.Nodegraph(K, 1_000_000, 4)
    graph.consume_seqfile_and_tag(str(reads))
    sub = graph.do_subset_partition(0, 0, False, False)
    graph.merge_subset(sub)
    graph.output_partitions(str(reads), str(tmp_path / "ref.part"))
    names = [line for line in (out / "t.fa.part").read_text().splitlines() if line.startswith(">")]
    ref = [
        line for line in (tmp_path / "ref.part").read_text().splitlines() if line.startswith(">")
    ]
    assert [n.split("\t")[0] for n in names] == [">a", ">b", ">c"]
    assert names[0] == names[2].replace(">c", ">a") and names[0] != names[1].replace(">b", ">a")
    assert len({n.split("\t")[1] for n in names}) == 2 == len({n.split("\t")[1] for n in ref})


def test_partition_reads_saved_maps_query_and_options(reads: Path, tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    maps = tmp_path / "maps"
    maps.mkdir()
    tag = m.khmer_tag_reads(reads_path=reads, ksize=K, **V)["tags"][1]
    first = m.khmer_partition_reads(
        reads,
        output_dir=out,
        ksize=K,
        subset_size=3,
        save_partitionmap_dir=maps,
        query_kmers=[tag],
        **V,
    )
    assert len(first["saved_partitionmaps"]) == 3
    assert len(first["per_subset"]) == 3
    out2 = tmp_path / "out2"
    out2.mkdir()
    again = m.khmer_partition_reads(
        reads, output_dir=out2, ksize=K, partitionmap_paths=first["saved_partitionmaps"], **V
    )
    assert (out / "t.fa.part").read_text() == (out2 / "t.fa.part").read_text()
    assert again["partition_count"] is None
    abundance = m.khmer_partition_reads(reads, annotate=False, ksize=K, abundance_range=[1, 1], **V)
    assert abundance["partition_count"]["n_partitions"] >= 1
    with pytest.raises(ConfigurationError):
        m.khmer_partition_reads(reads, ksize=K, **V)
    with pytest.raises(ConfigurationError):
        m.khmer_partition_reads(reads, output_dir=out, ksize=K, **V)


def test_join_partitions_and_load_partitioned_reads(reads: Path, tmp_path: Path) -> None:
    a, b, _ = _seqs(reads)
    tags_a = m.khmer_tag_positions(a, reads_path=reads, ksize=K, **V)["tags"]
    tags_b = m.khmer_tag_positions(b, reads_path=reads, ksize=K, **V)["tags"]
    joined = m.khmer_join_partitions(reads, tags_a[0]["kmer"], tags_b[0]["kmer"], ksize=K, **V)
    before, after = joined["partition_ids_before"], joined["partition_ids_after"]
    assert before[0] != before[1]
    assert after[0] == after[1]
    assert joined["partition_count_after"][0] == joined["partition_count_before"][0] - 1
    with pytest.raises(ConfigurationError):
        m.khmer_join_partitions(reads, "A" * K, tags_b[0]["kmer"], ksize=K, **V)
    out = tmp_path / "out"
    out.mkdir()
    m.khmer_partition_reads(reads, output_dir=out, ksize=K, **V)
    loaded = m.khmer_load_partitioned_reads(
        out / "t.fa.part", ksize=K, query_kmers=[tags_a[0]["kmer"]], **V
    )
    assert loaded["n_reads"] == 3 and loaded["n_partitions"] == 2
    assert loaded["partition_ids"][tags_a[0]["kmer"]] == 0
    assert len(loaded["tag_partition_ids"]) == 2
    assert set(loaded["tag_partition_ids"].values()) == {
        v for v in loaded["tag_partition_ids"].values() if v > 0
    }


def _knot_reads(tmp_path: Path) -> Path:
    repeat = _dna(100, 99)
    path = tmp_path / "knot.fa"
    path.write_text(
        "".join(f">k{i}\n{_dna(120, 100 + i)}{repeat}{_dna(120, 300 + i)}\n" for i in range(40))
    )
    return path


def test_find_knots_stages_and_stop_tag_files(tmp_path: Path) -> None:
    knots = _knot_reads(tmp_path)
    initial = m.khmer_find_knots(
        knots,
        stage="initial",
        ksize=21,
        save_stoptags_path=tmp_path / "st.bin",
        stoptags_text_path=tmp_path / "st.txt",
        **V,
    )
    assert initial["n_stop_tags"] > 0
    assert sorted((tmp_path / "st.txt").read_text().split()) == initial["stop_tags"]
    ref = khmer.Nodegraph(21, 1, 1)
    ref.load_stop_tags(str(tmp_path / "st.bin"))
    assert sorted(ref.get_stop_tags()) == initial["stop_tags"]
    full = m.khmer_find_knots(knots, stage="knots", ksize=21, **V)
    assert (
        full["excursion_kmer_count_threshold"] == 2
        and initial["excursion_kmer_count_threshold"] == 5
    )
    extra = m.khmer_find_knots(
        knots,
        stage="initial",
        ksize=21,
        load_stoptags_path=tmp_path / "st.bin",
        extra_stop_kmers=["A" * 21],
        **V,
    )
    assert "A" * 21 in extra["stop_tags"]
    with pytest.raises(ConfigurationError):
        m.khmer_find_knots(knots, stage="middle", ksize=21, **V)


def test_filter_stoptags_trims_and_drops(tmp_path: Path) -> None:
    seq = _dna(80, 4)
    reads = tmp_path / "r.fa"
    reads.write_text(f">r1\n{seq}\n>r2\n{_dna(30, 5)}N{_dna(30, 6)}\n>r3\n{_dna(60, 7)}\n")
    stop = seq[40 : 40 + K]
    result = m.khmer_filter_stoptags(reads, tmp_path / "f.fa", stop_kmers=[stop], ksize=K, **V)
    assert result["n_dropped_with_n"] == 1
    assert result["n_reads"] == 3 and result["n_kept"] == 2 and result["n_trimmed"] == 1
    lines = (tmp_path / "f.fa").read_text().splitlines()
    assert lines[0] == ">r1" and seq.startswith(lines[1]) and len(lines[1]) < len(seq)
    with pytest.raises(ConfigurationError):
        m.khmer_filter_stoptags(reads, tmp_path / "g.fa", ksize=K, **V)


def test_extract_partitions_distribution_groups_and_unassigned(tmp_path: Path) -> None:
    part = tmp_path / "x.part"
    part.write_text(
        "".join(f">r{i}\t{pid}\nACGTACGTAC\n" for i, pid in enumerate([1, 1, 1, 2, 2, 3, 0]))
    )
    result = m.khmer_extract_partitions(
        part, str(tmp_path / "ep"), min_partition_size=1, write_unassigned=True, **V
    )
    assert result["n_partitions"] == 3 and result["n_unassigned_reads"] == 1
    dist = (tmp_path / "ep.dist").read_text().splitlines()
    assert dist == ["1 1 1 1", "2 1 2 3", "3 1 3 6"] or dist[0].startswith("1 ")
    assert len(result["group_files"]) == 1
    group = Path(result["group_files"][0]).read_text()
    assert group.count(">") == 5
    assert Path(result["unassigned_path"]).read_text().count(">") == 1
    only_dist = m.khmer_extract_partitions(part, str(tmp_path / "nd"), write_groups=False, **V)
    assert only_dist["group_files"] == []
    with pytest.raises(ConfigurationError):
        m.khmer_extract_partitions(part, str(tmp_path / "ep"), **V)
    bad = tmp_path / "bad.part"
    bad.write_text(">r1\nACGTACGTAC\n")
    with pytest.raises(ConfigurationError):
        m.khmer_extract_partitions(bad, str(tmp_path / "bp"), **V)
