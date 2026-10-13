"""Tests for the G15 sourmash and khmer microbiome wrappers.

Every wrapped function is checked against a hand calculation or a direct call to the upstream
package. Inputs are built inside pytest's tmp_path; nothing is downloaded.
"""

from __future__ import annotations

import inspect
import math
import random
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import khmer_functions as kf
from dnakit.microbiome import sourmash_functions as sf


def _rc(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _random_dna(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _sketch(identifier: str, hashes: list[int], *, ksize: int = 3, scaled: int = 1) -> dict:
    return {"id": identifier, "ksize": ksize, "scaled": scaled, "hashes": hashes}


def _minhash(sourmash: object, hashes: list[int], *, ksize: int = 3, scaled: int = 1) -> object:
    minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)  # type: ignore[attr-defined]
    minhash.add_many(hashes)
    return minhash


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------- sourmash_functions


def test_public_signatures_are_keyword_only_verbose() -> None:
    for module in (sf, kf):
        for name in module.__all__:
            function = getattr(module, name)
            parameter = inspect.signature(function).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is True


def test_compare_minhash_matrix_jaccard_matches_hand_and_upstream() -> None:
    sourmash = pytest.importorskip("sourmash")
    left = _sketch("left", [1, 2, 3, 4])
    right = _sketch("right", [3, 4, 5, 6, 7, 8])
    result = sf.compare_minhash_matrix([left, right], verbose=False)
    # Shared {3, 4} and union of 8 hashes give Jaccard 2 / 8.
    assert result["ids"] == ["left", "right"]
    assert result["matrix"][0][0] == pytest.approx(1.0)
    assert result["matrix"][0][1] == pytest.approx(0.25)
    assert result["matrix"][1][0] == pytest.approx(0.25)
    signatures = [
        sourmash.signature.SourmashSignature(_minhash(sourmash, [1, 2, 3, 4]), name="left"),
        sourmash.signature.SourmashSignature(_minhash(sourmash, [3, 4, 5, 6, 7, 8]), name="right"),
    ]
    direct = sourmash.compare.compare_serial(signatures, True)
    flat_result = [value for row in result["matrix"] for value in row]
    flat_direct = [float(value) for row in direct.tolist() for value in row]
    assert flat_result == pytest.approx(flat_direct)


def test_compare_minhash_matrix_containment_is_directional() -> None:
    sourmash = pytest.importorskip("sourmash")
    left = _sketch("left", [1, 2, 3, 4])
    right = _sketch("right", [3, 4, 5, 6, 7, 8])
    result = sf.compare_minhash_matrix([left, right], metric="containment", verbose=False)
    # 2 of the 4 left hashes are in right; 2 of the 6 right hashes are in left.
    assert result["matrix"][0][1] == pytest.approx(0.5)
    assert result["matrix"][1][0] == pytest.approx(1 / 3)
    left_mh = _minhash(sourmash, [1, 2, 3, 4])
    right_mh = _minhash(sourmash, [3, 4, 5, 6, 7, 8])
    assert result["matrix"][0][1] == pytest.approx(left_mh.contained_by(right_mh))
    assert result["matrix"][1][0] == pytest.approx(right_mh.contained_by(left_mh))
    maximum = sf.compare_minhash_matrix([left, right], metric="max_containment", verbose=False)
    assert maximum["matrix"][0][1] == pytest.approx(0.5)
    assert maximum["matrix"][1][0] == pytest.approx(0.5)


def test_compare_minhash_matrix_validates_input() -> None:
    pytest.importorskip("sourmash")
    left = _sketch("left", [1, 2])
    with pytest.raises(ConfigurationError, match="metric must be one of"):
        sf.compare_minhash_matrix([left, left], metric="euclidean", verbose=False)
    with pytest.raises(ConfigurationError, match="at least two"):
        sf.compare_minhash_matrix([left], verbose=False)
    with pytest.raises(ConfigurationError, match="share ksize and scaled"):
        sf.compare_minhash_matrix([left, _sketch("b", [1], ksize=5)], verbose=False)
    with pytest.raises(ConfigurationError, match="sketch mapping"):
        sf.compare_minhash_matrix([left, "not a sketch"], verbose=False)  # type: ignore[list-item]


def test_estimate_mutation_rate_jaccard_matches_upstream() -> None:
    pytest.importorskip("sourmash")
    from sourmash.minhash import jaccard_to_distance

    result = sf.estimate_mutation_rate(
        0.5, measure="jaccard", ksize=21, scaled=100, n_unique_kmers=1000, verbose=False
    )
    direct = jaccard_to_distance(0.5, 21, 100, n_unique_kmers=1000)
    assert result["distance"] == pytest.approx(float(direct.dist))
    assert result["ani"] == pytest.approx(1.0 - float(direct.dist))
    assert result["p_nothing_in_common"] == pytest.approx(float(direct.p_nothing_in_common))
    assert result["distance_low"] is None
    assert 0.0 < result["distance"] < 1.0


def test_estimate_mutation_rate_containment_interval_matches_upstream() -> None:
    pytest.importorskip("sourmash")
    from sourmash.minhash import containment_to_distance

    result = sf.estimate_mutation_rate(
        0.5,
        measure="containment",
        ksize=21,
        scaled=100,
        sequence_len_bp=100_000,
        verbose=False,
    )
    direct = containment_to_distance(
        0.5, 21, 100, sequence_len_bp=100_000, confidence=0.95, estimate_ci=True
    )
    assert result["distance"] == pytest.approx(float(direct.dist))
    assert result["distance_low"] == pytest.approx(float(direct.dist_low))
    assert result["distance_high"] == pytest.approx(float(direct.dist_high))
    assert result["distance_low"] <= result["distance"] <= result["distance_high"]


def test_estimate_mutation_rate_validates_input() -> None:
    pytest.importorskip("sourmash")
    with pytest.raises(ConfigurationError, match="strictly between 0 and 1"):
        sf.estimate_mutation_rate(1.0, n_unique_kmers=100, verbose=False)
    with pytest.raises(ConfigurationError, match="sequence_len_bp or n_unique_kmers"):
        sf.estimate_mutation_rate(0.5, verbose=False)
    with pytest.raises(ConfigurationError, match="measure must be one of"):
        sf.estimate_mutation_rate(0.5, measure="mash", n_unique_kmers=100, verbose=False)


def test_set_size_error_probability_matches_binomial_sum() -> None:
    pytest.importorskip("sourmash")
    n, scaled, tolerance = 200, 10, 0.25
    # Sketch size is Binomial(n, 1 / scaled). The estimate is within tolerance when
    # |k * scaled - n| <= tolerance * n, and the boundary counts as inside.
    inside = [k for k in range(n + 1) if abs(k * scaled - n) <= tolerance * n]
    expected = sum(
        math.comb(n, k) * (1 / scaled) ** k * (1 - 1 / scaled) ** (n - k) for k in inside
    )
    result = sf.set_size_error_probability(n, scaled, relative_error=tolerance, verbose=False)
    assert result["probability_within_tolerance"] == pytest.approx(expected, rel=1e-9)
    assert result["expected_sketch_size"] == pytest.approx(20.0)
    with pytest.raises(ConfigurationError, match="relative_error must be strictly"):
        sf.set_size_error_probability(n, scaled, relative_error=0.0, verbose=False)


def test_hash_canonical_kmer_matches_scaled_one_sketch() -> None:
    sourmash = pytest.importorskip("sourmash")
    from sourmash.minhash import hash_murmur

    sequence = "ACGTTGCAAG"
    k = 5
    sketch = sourmash.MinHash(n=0, ksize=k, scaled=1)
    sketch.add_sequence(sequence, force=True)
    expected = {
        sf.hash_canonical_kmer(sequence[i : i + k], verbose=False)["hash"]
        for i in range(len(sequence) - k + 1)
    }
    assert {int(value) for value in sketch.hashes} == expected
    reverse = sf.hash_canonical_kmer("TACGT", verbose=False)
    lower = sf.hash_canonical_kmer("acgta", verbose=False)
    assert reverse["canonical"] == "ACGTA"
    assert reverse["strand"] == "reverse_complement"
    assert lower["strand"] == "forward"
    assert reverse["hash"] == lower["hash"] == hash_murmur("ACGTA", 42)
    assert sf.hash_canonical_kmer("ACGT", verbose=False)["strand"] == "self_reverse_complement"
    with pytest.raises(ConfigurationError, match="only A, C, G and T"):
        sf.hash_canonical_kmer("ACGN", verbose=False)
    with pytest.raises(ConfigurationError, match="at most 100"):
        sf.hash_canonical_kmer("A" * 101, verbose=False)


def test_downsample_minhash_sketch_keeps_hashes_below_cut() -> None:
    sourmash = pytest.importorskip("sourmash")
    # At scaled=4 the cut is 2**62. The chosen hashes are clear of the boundary.
    hashes = [10, 2**40, 2**62 - 2**40, 2**63, 2**64 - 1]
    sketch = _sketch("genome", hashes, scaled=1)
    result = sf.downsample_minhash_sketch(sketch, scaled=4, verbose=False)
    assert result["hashes"] == [10, 2**40, 2**62 - 2**40]
    assert (result["input_hash_count"], result["hash_count"], result["scaled"]) == (5, 3, 4)
    direct = sourmash.minhash.flatten_and_downsample_scaled(_minhash(sourmash, hashes), 4)
    assert result["hashes"] == sorted(int(value) for value in direct.hashes)
    with pytest.raises(ConfigurationError, match="cannot be upsampled"):
        sf.downsample_minhash_sketch(_sketch("fine", [1], scaled=4), scaled=1, verbose=False)


def test_intersect_minhash_sketches_hand_and_coarser_scale() -> None:
    sourmash = pytest.importorskip("sourmash")
    result = sf.intersect_minhash_sketches(
        _sketch("left", [1, 2, 3, 4]), _sketch("right", [3, 4, 5, 6]), verbose=False
    )
    assert result["hashes"] == [3, 4]
    assert (result["left_id"], result["right_id"], result["scaled"]) == ("left", "right", 1)
    coarse = sf.intersect_minhash_sketches(
        _sketch("left", [1, 2, 3, 4], scaled=1),
        _sketch("right", [3, 4, 5], scaled=2),
        verbose=False,
    )
    direct = sourmash.minhash.flatten_and_intersect_scaled(
        _minhash(sourmash, [1, 2, 3, 4]), _minhash(sourmash, [3, 4, 5], scaled=2)
    )
    assert coarse["scaled"] == 2
    assert coarse["hashes"] == sorted(int(value) for value in direct.hashes) == [3, 4]
    with pytest.raises(ConfigurationError, match="share ksize"):
        sf.intersect_minhash_sketches(
            _sketch("left", [1]), _sketch("right", [1], ksize=5), verbose=False
        )


def test_search_minhash_database_threshold_and_order() -> None:
    pytest.importorskip("sourmash")
    query = _sketch("query", [1, 2, 3, 4, 5, 6])
    database = [
        _sketch("small", [1, 2, 3]),
        _sketch("big", [1, 2, 3, 4, 5, 6, 7]),
        _sketch("none", [9, 10]),
    ]
    result = sf.search_minhash_database(query, database, threshold_bp=3, verbose=False)
    assert [match["id"] for match in result["matches"]] == ["big", "small"]
    assert [match["intersect_bp"] for match in result["matches"]] == [6, 3]
    small = result["matches"][1]
    # 3 shared hashes out of 6 query hashes and out of 3 database hashes.
    assert small["fraction_of_query_in_match"] == pytest.approx(0.5)
    assert small["fraction_of_match_in_query"] == pytest.approx(1.0)
    assert result["database_size"] == 3
    # The query has 6 bp at scale 1, so a 7 bp threshold is unattainable and sourmash refuses it.
    with pytest.raises(ConfigurationError, match="unattainable"):
        sf.search_minhash_database(query, database, threshold_bp=7, verbose=False)
    none = sf.search_minhash_database(query, database, threshold_bp=4, verbose=False)
    assert [match["id"] for match in none["matches"]] == ["big"]


def test_lowest_common_lineage_hand_and_upstream() -> None:
    pytest.importorskip("sourmash")
    from sourmash.lca import lca_utils

    lineages = [
        "d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales",
        "d__Bacteria;p__Firmicutes;c__Clostridia",
    ]
    result = sf.lowest_common_lineage(lineages, verbose=False)
    assert result["lineage_string"] == "d__Bacteria;p__Firmicutes"
    assert result["lowest_rank"] == "phylum"
    assert result["distinct_children_below_lca"] == 2
    tree = lca_utils.build_tree([lca_utils.make_lineage(text) for text in lineages])
    upstream, children = lca_utils.find_lca(tree)
    assert [item.name for item in upstream] == ["d__Bacteria", "p__Firmicutes"]
    assert children == result["distinct_children_below_lca"]
    single = sf.lowest_common_lineage([lineages[0], lineages[0]], verbose=False)
    assert single["lowest_rank"] == "order"
    assert single["distinct_children_below_lca"] == 0
    with pytest.raises(ConfigurationError, match="';' as the separator"):
        sf.lowest_common_lineage(["d__Bacteria,p__Firmicutes"], verbose=False)
    with pytest.raises(ConfigurationError, match="at least one"):
        sf.lowest_common_lineage([], verbose=False)
    with pytest.raises(ConfigurationError, match="more than 8 ranks"):
        sf.lowest_common_lineage([";".join(["x__a"] * 9)], verbose=False)


# ------------------------------------------------------------------- khmer_functions


def test_hash_kmer_murmur3_matches_khmer_and_is_strand_symmetric() -> None:
    khmer = pytest.importorskip("khmer")
    result = kf.hash_kmer_murmur3("acgta", verbose=False)
    assert result["kmer"] == "ACGTA"
    assert result["hash"] == khmer.hash_murmur3("ACGTA")
    # khmer hashes lowercase input as different bases; the wrapper uppercases it first.
    assert khmer.hash_murmur3("acgta") != result["hash"]
    assert kf.hash_kmer_murmur3("TACGT", verbose=False)["hash"] == result["hash"]
    raw = kf.hash_kmer_murmur3("ACGTA", canonical=False, verbose=False)
    assert raw["hash"] == khmer.hash_no_rc_murmur3("ACGTA")
    rng = random.Random(5)
    for _ in range(40):
        length = rng.randint(1, 32)
        sequence = _random_dna(rng.randint(0, 10**6), length)
        forward = kf.hash_kmer_murmur3(sequence, verbose=False)["hash"]
        assert forward == kf.hash_kmer_murmur3(_rc(sequence), verbose=False)["hash"]
    with pytest.raises(ConfigurationError, match="only A, C, G and T"):
        kf.hash_kmer_murmur3("ACGU", verbose=False)
    with pytest.raises(ConfigurationError, match="at most 32"):
        kf.hash_kmer_murmur3("A" * 33, verbose=False)


def test_encode_kmer_2bit_hand_values_and_upstream() -> None:
    khmer = pytest.importorskip("khmer")
    # khmer order is A=0, T=1, C=2, G=3, first base in the highest bits.
    assert kf.encode_kmer_2bit("ACGT", verbose=False)["code"] == 45  # 00 10 11 01
    assert kf.encode_kmer_2bit("AAAA", verbose=False)["code"] == 0
    assert kf.encode_kmer_2bit("TTTT", canonical=False, verbose=False)["code"] == 85  # 01 01 01 01
    assert kf.encode_kmer_2bit("CCCC", canonical=False, verbose=False)["code"] == 170  # 10 10 10 10
    assert kf.encode_kmer_2bit("GGGG", canonical=False, verbose=False)["code"] == 255  # 11 11 11 11
    # Canonical: TTTT (85) pairs with AAAA (0), and GGGG (255) pairs with CCCC (170).
    assert kf.encode_kmer_2bit("TTTT", verbose=False)["code"] == 0
    assert kf.encode_kmer_2bit("GGGG", verbose=False)["code"] == 170
    assert kf.encode_kmer_2bit("GGGGA", canonical=False, verbose=False)["code"] == 1020
    # The canonical code is the smaller of GGGGA (1020) and its reverse complement TCCCC.
    assert kf.encode_kmer_2bit("GGGGA", verbose=False)["code"] == 426  # 01 10 10 10 10
    rng = random.Random(9)
    for _ in range(30):
        k = rng.randint(1, 32)
        sequence = _random_dna(rng.randint(0, 10**6), k)
        code = kf.encode_kmer_2bit(sequence, canonical=False, verbose=False)["code"]
        assert code == khmer.forward_hash_no_rc(sequence, k)
        assert kf.encode_kmer_2bit(sequence, verbose=False)["code"] == khmer.forward_hash(
            sequence, k
        )


def test_decode_kmer_2bit_round_trip_and_bounds() -> None:
    pytest.importorskip("khmer")
    assert kf.decode_kmer_2bit(45, ksize=4, verbose=False)["kmer"] == "ACGT"
    rng = random.Random(21)
    for _ in range(40):
        k = rng.randint(1, 32)
        sequence = _random_dna(rng.randint(0, 10**6), k)
        plain = kf.encode_kmer_2bit(sequence, canonical=False, verbose=False)["code"]
        assert kf.decode_kmer_2bit(plain, ksize=k, verbose=False)["kmer"] == sequence
        canonical = kf.encode_kmer_2bit(sequence, verbose=False)["code"]
        assert kf.decode_kmer_2bit(canonical, ksize=k, verbose=False)["kmer"] in {
            sequence,
            _rc(sequence),
        }
    with pytest.raises(ConfigurationError, match="value must be an integer"):
        kf.decode_kmer_2bit(4**4, ksize=4, verbose=False)
    with pytest.raises(ConfigurationError, match="ksize must be an integer"):
        kf.decode_kmer_2bit(1, ksize=0, verbose=False)


def _trim_fixture(tmp_path: Path, *, lowercase: bool = False) -> Path:
    ref = "ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG"  # 30 bases
    novel = "CAGTCGATGCATGCAA"  # 16 bases, all 9-mers distinct and seen once
    assert len({novel[i : i + 9] for i in range(len(novel) - 8)}) == len(novel) - 8
    bad = ref + "GGGGGGGGGGGGCCCCCCCCCCCC"  # novel 3' tail after the reference
    records = [
        ("ref1", ref),
        ("ref2", ref),
        ("ref3", ref),
        ("bad", bad),
        ("novel", novel),
        ("nonacgt", "ACGTNNNNACGTACGTACGT"),
        ("short", "ACGT"),
    ]
    text = "".join(
        f">{name}\n{sequence.lower() if lowercase else sequence}\n" for name, sequence in records
    )
    return _write(tmp_path / "reads.fa", text)


def test_trim_low_abundance_reads_trims_tail_and_drops_novel(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    source = _trim_fixture(tmp_path)
    output = tmp_path / "trimmed.fa"
    result = kf.trim_low_abundance_reads(source, output, ksize=9, cutoff=2, verbose=False)
    # The 'bad' read keeps 30 bases: its first k-mer with count below 2 starts at base 22 and
    # covers the first tail base, so the cut is at base 30. 'novel' has count-1 k-mers only and
    # is trimmed to nothing, which leaves it shorter than k.
    assert result["reads"] == 7
    assert result["used_reads"] == 5
    assert result["skipped_non_acgt_reads"] == 1
    assert result["skipped_short_reads"] == 1
    assert result["kept_untrimmed_reads"] == 3
    assert result["trimmed_reads"] == 1
    assert result["dropped_short_reads"] == 1
    assert result["kept_reads"] == 4
    ref = "ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG"
    assert output.read_text(encoding="utf-8") == (
        f">ref1\n{ref}\n>ref2\n{ref}\n>ref3\n{ref}\n>bad\n{ref}\n"
    )


def test_trim_low_abundance_reads_lowercase_input_matches_uppercase(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    upper_dir = tmp_path / "upper"
    lower_dir = tmp_path / "lower"
    upper_dir.mkdir()
    lower_dir.mkdir()
    upper = kf.trim_low_abundance_reads(
        _trim_fixture(upper_dir), upper_dir / "out.fa", ksize=9, cutoff=2, verbose=False
    )
    lower = kf.trim_low_abundance_reads(
        _trim_fixture(lower_dir, lowercase=True),
        lower_dir / "out.fa",
        ksize=9,
        cutoff=2,
        verbose=False,
    )
    ignored = {"source", "output"}
    assert {k: v for k, v in upper.items() if k not in ignored} == {
        k: v for k, v in lower.items() if k not in ignored
    }
    assert (upper_dir / "out.fa").read_text(encoding="utf-8") == (lower_dir / "out.fa").read_text(
        encoding="utf-8"
    )


def test_trim_low_abundance_reads_fastq_keeps_quality(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    ref = "ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG"
    bad = ref + "GGGGGGGGGGGGCCCCCCCCCCCC"
    fastq = "".join(
        f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"
        for name, sequence in [("r1", ref), ("r2", ref), ("r3", ref), ("bad", bad)]
    )
    source = _write(tmp_path / "reads.fq", fastq)
    output = tmp_path / "trimmed.fastq"
    result = kf.trim_low_abundance_reads(source, output, ksize=9, cutoff=2, verbose=False)
    assert result["trimmed_reads"] == 1
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[-4:] == ["@bad", ref, "+", "I" * 30]
    fasta_dir = tmp_path / "fasta_only"
    fasta_dir.mkdir()
    with pytest.raises(ConfigurationError, match="FASTQ output needs FASTQ input"):
        kf.trim_low_abundance_reads(
            _trim_fixture(fasta_dir), fasta_dir / "out.fq", ksize=9, cutoff=2, verbose=False
        )
    assert not (fasta_dir / "out.fq").exists()


def test_trim_low_abundance_reads_validates_arguments(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    source = _trim_fixture(tmp_path)
    with pytest.raises(ConfigurationError, match="ksize must be an integer"):
        kf.trim_low_abundance_reads(source, tmp_path / "x.fa", ksize=33, verbose=False)
    with pytest.raises(ConfigurationError, match="cutoff must be an integer"):
        kf.trim_low_abundance_reads(source, tmp_path / "x.fa", cutoff=0, verbose=False)
    with pytest.raises(ConfigurationError, match="must end with"):
        kf.trim_low_abundance_reads(source, tmp_path / "x.txt", verbose=False)
    existing = _write(tmp_path / "existing.fa", "keep me\n")
    with pytest.raises(ConfigurationError, match="will not be overwritten"):
        kf.trim_low_abundance_reads(source, existing, ksize=9, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep me\n"
    with pytest.raises(ConfigurationError, match="max_reads=2"):
        kf.trim_low_abundance_reads(source, tmp_path / "y.fa", ksize=9, max_reads=2, verbose=False)


def test_inspect_khmer_graph_matches_upstream(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    countgraph = khmer.Countgraph(7, 1009, 3)
    countgraph.consume("ACGTTGCAAGCTAGCTTTG")
    count_path = tmp_path / "counts.cg"
    countgraph.save(str(count_path))
    nodegraph = khmer.Nodegraph(7, 1009, 3)
    nodegraph.consume("ACGTTGCAAGCTAGCTTTGCAAGG")
    node_path = tmp_path / "nodes.ng"
    nodegraph.save(str(node_path))

    counts = kf.inspect_khmer_graph(count_path, kind="countgraph", verbose=False)
    info = khmer.extract_countgraph_info(str(count_path))
    assert counts["ksize"] == info.ksize == 7
    assert counts["table_size"] == info.table_size
    assert counts["n_tables"] == info.n_tables == 3
    assert counts["n_occupied"] == countgraph.n_occupied() == info.n_occupied
    assert counts["table_type_flag"] == 1

    nodes = kf.inspect_khmer_graph(node_path, kind="nodegraph", verbose=False)
    header = khmer.extract_nodegraph_info(str(node_path))
    assert nodes["table_size"] == header[1]
    assert nodes["n_tables"] == header[2]
    assert nodes["n_occupied"] == nodegraph.n_occupied()
    assert nodes["table_type_flag"] == 2

    with pytest.raises(ConfigurationError, match="is not a khmer countgraph"):
        kf.inspect_khmer_graph(node_path, kind="countgraph", verbose=False)
    with pytest.raises(ConfigurationError, match="is not a khmer nodegraph"):
        kf.inspect_khmer_graph(count_path, kind="nodegraph", verbose=False)
    with pytest.raises(ConfigurationError, match="kind must be one of"):
        kf.inspect_khmer_graph(count_path, kind="bloom", verbose=False)


def test_query_saved_countgraph_matches_direct_lookup(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    countgraph = khmer.Countgraph(5, 1009, 2)
    sequence = "ACGTTGCAAGCTAG"
    countgraph.consume(sequence)
    countgraph.consume(sequence)
    path = tmp_path / "saved.cg"
    countgraph.save(str(path))
    result = kf.query_saved_countgraph(path, ["acgtt", "AACGT", "GGGGG"], verbose=False)
    # ACGTT occurs once in the sequence and AACGT is its reverse complement; each was consumed
    # twice, so both report 2. GGGGG is absent.
    assert result["counts"] == {"AACGT": 2, "ACGTT": 2, "GGGGG": 0}
    loaded = khmer.load_countgraph(str(path))
    assert result["counts"]["ACGTT"] == loaded.get("ACGTT")
    with pytest.raises(ConfigurationError, match="table ksize is 5"):
        kf.query_saved_countgraph(path, ["ACGT"], verbose=False)
    with pytest.raises(ConfigurationError, match="not a khmer countgraph"):
        kf.query_saved_countgraph(_write(tmp_path / "x.fa", ">a\nACGT\n"), ["ACGTA"], verbose=False)


def test_estimate_kmer_table_size_matches_khmer_args() -> None:
    pytest.importorskip("khmer")
    khmer_args = pytest.importorskip("khmer.khmer_args")
    result = kf.estimate_kmer_table_size(1_000_000, fp_rate=0.01, verbose=False)
    direct = khmer_args.optimal_size(1_000_000, fp_rate=0.01)
    assert result["num_htables"] == direct.num_htables
    assert result["htable_size"] == direct.htable_size
    assert result["mem_use"] == direct.mem_use
    assert result["expected_fp_rate"] == pytest.approx(0.01, rel=1e-4)
    capped = kf.estimate_kmer_table_size(1_000_000, mem_cap=10_000_000, verbose=False)
    direct_cap = khmer_args.optimal_size(1_000_000, mem_cap=10_000_000)
    assert capped["htable_size"] == direct_cap.htable_size
    assert capped["expected_fp_rate"] == pytest.approx(direct_cap.fp_rate)
    with pytest.raises(ConfigurationError, match="exactly one of"):
        kf.estimate_kmer_table_size(1_000, mem_cap=100, fp_rate=0.1, verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one of"):
        kf.estimate_kmer_table_size(1_000, verbose=False)
    with pytest.raises(ConfigurationError, match="fp_rate must be a number"):
        kf.estimate_kmer_table_size(1_000, fp_rate=1.5, verbose=False)


def test_expected_kmer_table_collisions_matches_formula_and_never_exits(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    ref = "ACGTTGCAAGCTAGCTTTGCAAGGCCTAGG"
    source = _write(tmp_path / "ref.fa", f">a\n{ref}\n>b\n{ref}\n")
    result = kf.expected_kmer_table_collisions(
        source, ksize=9, table_size=10007, n_tables=2, verbose=False
    )
    graph = khmer.Nodegraph(9, 10007, 2)
    graph.consume(ref)
    sizes = graph.hashsizes()
    expected = (graph.n_occupied() / min(sizes)) ** len(sizes)
    assert result["hashsizes"] == [int(size) for size in sizes]
    assert result["occupied_bins"] == graph.n_occupied()
    assert result["expected_false_positive_rate"] == pytest.approx(expected)
    assert result["within_threshold"] is True
    strict = kf.expected_kmer_table_collisions(
        source, ksize=9, table_size=10007, n_tables=2, max_false_positive=1e-12, verbose=False
    )
    assert strict["within_threshold"] is False
    with pytest.raises(ConfigurationError, match="max_false_positive"):
        kf.expected_kmer_table_collisions(source, ksize=9, max_false_positive=0.0, verbose=False)
