"""Contract tests for optional sourmash sketch and khmer k-mer backends on microbial DNA."""

from __future__ import annotations

import random
import sys
from collections import Counter
from pathlib import Path
from typing import ClassVar

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.microbiome.backends import (
    compare_minhash_sketches,
    estimate_unique_kmers,
    kmer_abundance_histogram,
    normalize_reads_by_median,
    query_kmer_counts,
    read_sourmash_signatures,
    sketch_fasta_minhash,
    write_sourmash_signatures,
)


def _random_dna(seed: int, length: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _canonical_kmer_counts(sequence: str, k: int) -> Counter:
    """Exact canonical k-mer occurrence counts, used as the reference for khmer."""

    counts: Counter = Counter()
    for index in range(len(sequence) - k + 1):
        kmer = sequence[index : index + k]
        counts[min(kmer, _reverse_complement(kmer))] += 1
    return counts


def _write(tmp_path: Path, text: str, name: str = "reads.fa") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _sketch(hashes: list[int], *, ksize: int = 3, scaled: int = 1, identifier: str = "s") -> dict:
    return {"id": identifier, "ksize": ksize, "scaled": scaled, "hashes": hashes}


def test_sketch_fasta_minhash_matches_direct_sourmash(tmp_path: Path) -> None:
    sourmash = pytest.importorskip("sourmash")
    genome = _random_dna(7, 5_000)
    path = _write(tmp_path, f">genome_1 description\n{genome[:2_500]}\n{genome[2_500:]}\n")

    result = sketch_fasta_minhash(path, ksize=21, scaled=100, verbose=False)

    direct = sourmash.MinHash(n=0, ksize=21, scaled=100)
    direct.add_sequence(genome, force=True)
    sketch = result["sketches"][0]
    assert result["records"] == 1
    assert sketch["id"] == "genome_1"
    assert sketch["length"] == 5_000
    assert sketch["hashes"] == sorted(int(value) for value in direct.hashes)
    assert sketch["hash_count"] == len(direct)


def test_sketch_counts_canonical_kmers_by_hand(tmp_path: Path) -> None:
    pytest.importorskip("sourmash")
    # ACGTACGT has 3-mers ACG, CGT, GTA, TAC, ACG, CGT. Reverse complements pair
    # ACG with CGT and GTA with TAC, so two canonical k-mers remain.
    # In ACGNNACG only ACG is free of N, and it occurs twice: one canonical k-mer.
    # A 2-base record is shorter than k and yields no hashes.
    path = _write(tmp_path, ">pair\nACGTACGT\n>lower\nacgtacgt\n>with_n\nACGNNACG\n>short\nAC\n")

    sketches = sketch_fasta_minhash(path, ksize=3, scaled=1, verbose=False)["sketches"]

    counts = {sketch["id"]: sketch["hash_count"] for sketch in sketches}
    assert counts == {"pair": 2, "lower": 2, "with_n": 1, "short": 0}


def test_sketch_fasta_minhash_validates_limits_and_input(tmp_path: Path) -> None:
    pytest.importorskip("sourmash")
    path = _write(tmp_path, ">a\nACGTACGT\n>b\nACGTACGT\n")

    with pytest.raises(ConfigurationError, match="ksize must be an integer"):
        sketch_fasta_minhash(path, ksize=True, verbose=False)
    with pytest.raises(ConfigurationError, match="ksize must be an integer"):
        sketch_fasta_minhash(path, ksize=101, verbose=False)
    with pytest.raises(ConfigurationError, match="scaled must be an integer"):
        sketch_fasta_minhash(path, scaled=0, verbose=False)
    with pytest.raises(ConfigurationError, match="max_records=1"):
        sketch_fasta_minhash(path, max_records=1, verbose=False)
    with pytest.raises(ConfigurationError, match="max_bases=10"):
        sketch_fasta_minhash(path, max_bases=10, verbose=False)
    bad_inputs = {
        "no_header.fa": "ACGT\n",
        "digits.fa": ">a\nAC1T\n",
        "nul.fa": ">a\nAC\x00T\n",
        "no_name.fa": ">\nACGT\n",
        "empty.fa": "\n",
    }
    for name, text in bad_inputs.items():
        with pytest.raises(ConfigurationError):
            sketch_fasta_minhash(_write(tmp_path, text, name), verbose=False)


def test_compare_minhash_sketches_matches_hand_calculation() -> None:
    pytest.importorskip("sourmash")
    # Left has {1, 2, 3, 4}; right has {3, 4, 5, 6, 7, 8}. Two hashes are shared.
    # Jaccard = 2 / 8. Containment is directional: 2 / 4 for left, 2 / 6 for right.
    result = compare_minhash_sketches(
        _sketch([1, 2, 3, 4], identifier="left"),
        _sketch([3, 4, 5, 6, 7, 8], identifier="right"),
        verbose=False,
    )
    assert result["shared_hash_count"] == 2
    assert result["jaccard"] == pytest.approx(0.25)
    assert result["left_contained_by_right"] == pytest.approx(0.5)
    assert result["right_contained_by_left"] == pytest.approx(1 / 3)
    assert (result["left_id"], result["right_id"]) == ("left", "right")


def test_compare_minhash_sketches_matches_sourmash_on_fasta_sketches(tmp_path: Path) -> None:
    sourmash = pytest.importorskip("sourmash")
    left_seq = _random_dna(11, 4_000)
    right_seq = left_seq[:2_000] + _random_dna(12, 2_000)
    path = _write(tmp_path, f">left\n{left_seq}\n>right\n{right_seq}\n")
    left, right = sketch_fasta_minhash(path, ksize=21, scaled=100, verbose=False)["sketches"]

    result = compare_minhash_sketches(left, right, verbose=False)

    left_direct = sourmash.MinHash(n=0, ksize=21, scaled=100)
    left_direct.add_sequence(left_seq, force=True)
    right_direct = sourmash.MinHash(n=0, ksize=21, scaled=100)
    right_direct.add_sequence(right_seq, force=True)
    assert result["jaccard"] == pytest.approx(
        left_direct.similarity(right_direct, ignore_abundance=True)
    )
    assert result["left_contained_by_right"] == pytest.approx(
        left_direct.contained_by(right_direct)
    )
    assert result["right_contained_by_left"] == pytest.approx(
        right_direct.contained_by(left_direct)
    )
    assert result["shared_hash_count"] == len(set(left["hashes"]) & set(right["hashes"]))


def test_compare_minhash_sketches_rejects_incompatible_sketches() -> None:
    pytest.importorskip("sourmash")
    base = _sketch([1, 2, 3, 4])
    with pytest.raises(ConfigurationError, match="share ksize and scaled"):
        compare_minhash_sketches(base, _sketch([1, 2], ksize=5), verbose=False)
    with pytest.raises(ConfigurationError, match="share ksize and scaled"):
        compare_minhash_sketches(base, _sketch([1, 2], scaled=10), verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty list"):
        compare_minhash_sketches(base, _sketch([]), verbose=False)
    with pytest.raises(ConfigurationError, match="unsigned 64-bit"):
        compare_minhash_sketches(base, _sketch([True, 2]), verbose=False)
    with pytest.raises(ConfigurationError, match="unsigned 64-bit"):
        compare_minhash_sketches(base, _sketch([2**64]), verbose=False)
    with pytest.raises(ConfigurationError, match="scaled range"):
        compare_minhash_sketches(
            _sketch([2**63], scaled=100), _sketch([1, 2], scaled=100), verbose=False
        )
    with pytest.raises(ConfigurationError, match="must be a sketch mapping"):
        compare_minhash_sketches(base, [1, 2, 3], verbose=False)


def test_write_and_read_sourmash_signatures_round_trip(tmp_path: Path) -> None:
    sourmash = pytest.importorskip("sourmash")
    path = _write(tmp_path, f">a\n{_random_dna(21, 3_000)}\n>b\n{_random_dna(22, 3_000)}\n")
    sketches = sketch_fasta_minhash(path, ksize=21, scaled=100, verbose=False)["sketches"]
    target = tmp_path / "sketches.sig.gz"

    written = write_sourmash_signatures(sketches, target, verbose=False)
    assert written["signatures"] == 2

    read = read_sourmash_signatures(target, verbose=False)
    assert [sketch["id"] for sketch in read["sketches"]] == ["a", "b"]
    assert [sketch["hashes"] for sketch in read["sketches"]] == [s["hashes"] for s in sketches]
    assert all(sketch["length"] is None and sketch["md5sum"] for sketch in read["sketches"])
    upstream = list(sourmash.load_file_as_signatures(str(target)))
    assert [sorted(int(value) for value in sig.minhash.hashes) for sig in upstream] == [
        s["hashes"] for s in sketches
    ]

    with pytest.raises(ConfigurationError, match="already exists"):
        write_sourmash_signatures(sketches, target, verbose=False)
    with pytest.raises(ConfigurationError, match="must end with"):
        write_sourmash_signatures(sketches, tmp_path / "sketches.json", verbose=False)
    with pytest.raises(ConfigurationError, match="max_signatures=1"):
        read_sourmash_signatures(target, max_signatures=1, verbose=False)


def test_read_sourmash_signatures_rejects_non_dna_signatures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sourmash = pytest.importorskip("sourmash")

    class _Minhash:
        moltype = "protein"
        scaled = 1
        ksize = 9
        hashes: ClassVar[dict[int, int]] = {5: 0}

    class _Signature:
        name = "prot"
        minhash = _Minhash()

        def md5sum(self) -> str:
            return "0" * 32

    # Upstream reload of short protein signatures is unreliable in this environment, so the
    # moltype rule is checked with a stand-in loader that returns one protein-typed signature.
    monkeypatch.setattr(
        sourmash, "load_file_as_signatures", lambda *_args, **_kwargs: [_Signature()]
    )
    path = _write(tmp_path, ">placeholder\nACGT\n", "protein.sig")
    with pytest.raises(ConfigurationError, match="only DNA signatures"):
        read_sourmash_signatures(path, verbose=False)


def test_read_sourmash_signatures_rejects_malformed_files(tmp_path: Path) -> None:
    pytest.importorskip("sourmash")
    (tmp_path / "garbage.sig").write_text("not json", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="could not read"):
        read_sourmash_signatures(tmp_path / "garbage.sig", verbose=False)
    fasta = _write(tmp_path, ">a\nACGT\n")
    with pytest.raises(ConfigurationError, match="could not read"):
        read_sourmash_signatures(fasta, verbose=False)


def test_missing_sourmash_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path, ">a\nACGTACGT\n")
    monkeypatch.setitem(sys.modules, "sourmash", None)
    with pytest.raises(BackendUnavailableError, match="bio-sourmash"):
        sketch_fasta_minhash(path, verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-sourmash"):
        compare_minhash_sketches(_sketch([1, 2]), _sketch([2, 3]), verbose=False)


def _fasta(tmp_path: Path, records: list[tuple[str, str]], name: str = "reads.fa") -> Path:
    text = "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records)
    return _write(tmp_path, text, name)


def test_query_kmer_counts_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    # ACGTACGT has canonical 3-mer occurrences: ACG class 4 (ACG, CGT x2 each), GTA class 2
    # (GTA, TAC). The lowercase copy adds the same counts; ACGNNACG and AC are skipped.
    path = _fasta(
        tmp_path,
        [("p", "ACGTACGT"), ("n", "ACGNNACG"), ("s", "AC"), ("l", "acgtacgt")],
    )

    result = query_kmer_counts(path, ["ACG", "CGT", "GTA", "TAC", "AAA"], verbose=False)

    assert result["counts"] == {"ACG": 8, "CGT": 8, "GTA": 4, "TAC": 4, "AAA": 0}
    assert (result["reads"], result["used_reads"]) == (4, 2)
    assert (result["skipped_non_acgt_reads"], result["skipped_short_reads"]) == (1, 1)
    assert result["kmer_positions"] == 12


def test_query_kmer_counts_matches_exact_canonical_reference(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    sequence = _random_dna(3, 2_000)
    path = _fasta(tmp_path, [("chr", sequence)])
    queries = [sequence[start : start + 11] for start in (0, 100, 500, 1_000, 1_500)]

    result = query_kmer_counts(path, queries, verbose=False)

    exact = _canonical_kmer_counts(sequence, 11)
    expected = {query: exact[min(query, _reverse_complement(query))] for query in queries}
    assert result["counts"] == expected
    assert result["ksize"] == 11


def test_query_kmer_counts_and_estimates_validate_inputs(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    path = _fasta(tmp_path, [("a", "ACGTACGT"), ("b", "ACGTTTGA")])
    with pytest.raises(ConfigurationError, match="k-mer length must be in"):
        query_kmer_counts(path, ["A" * 33], verbose=False)
    with pytest.raises(ConfigurationError, match="same length"):
        query_kmer_counts(path, ["ACG", "ACGT"], verbose=False)
    with pytest.raises(ConfigurationError, match="uppercase A, C, G and T"):
        query_kmer_counts(path, ["acg"], verbose=False)
    with pytest.raises(ConfigurationError, match="uppercase A, C, G and T"):
        query_kmer_counts(path, ["ACN"], verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty list"):
        query_kmer_counts(path, [], verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty list"):
        query_kmer_counts(path, "ACG", verbose=False)
    with pytest.raises(ConfigurationError, match="table_size must be an integer"):
        query_kmer_counts(path, ["ACG"], table_size=True, verbose=False)
    with pytest.raises(ConfigurationError, match="table_size must be an integer"):
        query_kmer_counts(path, ["ACG"], table_size=999, verbose=False)
    with pytest.raises(ConfigurationError, match="n_tables must be an integer"):
        query_kmer_counts(path, ["ACG"], n_tables=0, verbose=False)
    with pytest.raises(ConfigurationError, match="max_reads=1"):
        query_kmer_counts(path, ["ACG"], max_reads=1, verbose=False)
    with pytest.raises(ConfigurationError, match="max_bases=5"):
        query_kmer_counts(path, ["ACG"], max_bases=5, verbose=False)
    with pytest.raises(ConfigurationError, match="ksize must be an integer in"):
        estimate_unique_kmers(path, ksize=33, verbose=False)
    bad = _write(tmp_path, "this is not a read file\n", "bad.txt")
    with pytest.raises(ConfigurationError, match="could not parse"):
        query_kmer_counts(bad, ["ACG"], verbose=False)


def test_estimate_unique_kmers_matches_hand_count_and_upstream(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    small = _fasta(tmp_path, [("chr", "ACGTACGT")], "small.fa")
    # Canonical 3-mers of ACGTACGT are {ACG, GTA}: two distinct k-mers.
    assert estimate_unique_kmers(small, ksize=3, table_size=1_000, verbose=False)[
        "unique_kmers_estimate"
    ] == pytest.approx(2, abs=0.5)

    sequence = _random_dna(5, 20_000)
    path = _fasta(tmp_path, [("chr", sequence)], "large.fa")
    result = estimate_unique_kmers(path, ksize=21, table_size=1_000_000, verbose=False)

    direct = khmer.Nodegraph(21, 1_000_000, 4)
    direct.consume(sequence)
    exact = len(_canonical_kmer_counts(sequence, 21))
    assert result["unique_kmers_estimate"] == direct.n_unique_kmers()
    assert result["unique_kmers_estimate"] == pytest.approx(exact, rel=0.01)
    assert result["kmer_positions"] == 20_000 - 21 + 1


def test_kmer_abundance_histogram_matches_exact_counts(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    # The same three reads as the normalization example: 36 k-mer positions in total.
    path = _fasta(
        tmp_path,
        [("r1", "ACGTTGCAAGCTAGCT"), ("r2", "ACGTTGCAAGCTAGCT"), ("r3", "TTTTGGGGCCCCAAAA")],
    )

    result = kmer_abundance_histogram(path, ksize=5, verbose=False)

    # Twelve canonical 5-mers occur twice and three occur four times: 2 * 12 + 4 * 3 = 36.
    assert result["histogram"] == [{"count": 2, "kmers": 12}, {"count": 4, "kmers": 3}]
    assert result["distinct_kmers"] == 15
    reads = ["ACGTTGCAAGCTAGCT", "ACGTTGCAAGCTAGCT", "TTTTGGGGCCCCAAAA"]
    exact: Counter = Counter()
    for read in reads:
        exact.update(_canonical_kmer_counts(read, 5))
    assert result["distinct_kmers"] == len(exact)
    expected = Counter(exact.values())
    assert result["histogram"] == [
        {"count": count, "kmers": kmers} for count, kmers in sorted(expected.items())
    ]


def test_kmer_abundance_histogram_rejects_unsupported_reads(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    lowercase = _fasta(tmp_path, [("a", "acgtacgt")], "lower.fa")
    with pytest.raises(ConfigurationError, match="uppercase A, C, G and T"):
        kmer_abundance_histogram(lowercase, ksize=3, verbose=False)
    with_n = _fasta(tmp_path, [("n", "ACGNNACGT")], "n.fa")
    with pytest.raises(ConfigurationError, match="uppercase A, C, G and T"):
        kmer_abundance_histogram(with_n, ksize=3, verbose=False)
    short = _fasta(tmp_path, [("s", "AC")], "short.fa")
    with pytest.raises(ConfigurationError, match="shorter than ksize"):
        kmer_abundance_histogram(short, ksize=3, verbose=False)


def test_normalize_reads_by_median_follows_hand_example(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    # Median counts before each read: r1 = 0 (kept), r2 = 2 (its 5-mers repeat inside the
    # read, so they were already counted once), r3 = 0 (kept).
    path = _fasta(
        tmp_path,
        [
            ("r1", "ACGTTGCAAGCTAGCT"),
            ("r2", "ACGTTGCAAGCTAGCT"),
            ("n", "ACGNNACGTAG"),
            ("s", "AC"),
            ("r3", "TTTTGGGGCCCCAAAA"),
        ],
    )
    target = tmp_path / "normalized.fa"

    result = normalize_reads_by_median(path, target, ksize=5, cutoff=2, verbose=False)

    assert target.read_text(encoding="utf-8") == (">r1\nACGTTGCAAGCTAGCT\n>r3\nTTTTGGGGCCCCAAAA\n")
    assert result["kept_reads"] == 2
    assert result["discarded_high_median_reads"] == 1
    assert (result["skipped_non_acgt_reads"], result["skipped_short_reads"]) == (1, 1)

    wide = normalize_reads_by_median(path, tmp_path / "wide.fa", ksize=5, cutoff=3, verbose=False)
    assert (wide["kept_reads"], wide["discarded_high_median_reads"]) == (3, 0)


def test_normalize_reads_by_median_validates_output(tmp_path: Path) -> None:
    pytest.importorskip("khmer")
    path = _fasta(tmp_path, [("r1", "ACGTTGCAAGCTAGCT")])
    target = tmp_path / "normalized.fa"
    normalize_reads_by_median(path, target, ksize=5, cutoff=2, verbose=False)
    with pytest.raises(ConfigurationError, match="already exists"):
        normalize_reads_by_median(path, target, ksize=5, cutoff=2, verbose=False)
    with pytest.raises(ConfigurationError, match="must end with"):
        normalize_reads_by_median(path, tmp_path / "out.txt", ksize=5, cutoff=2, verbose=False)
    with pytest.raises(ConfigurationError, match="cutoff must be an integer"):
        normalize_reads_by_median(path, tmp_path / "x.fa", cutoff=0, verbose=False)


def test_missing_khmer_raises_install_hint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = _fasta(tmp_path, [("a", "ACGTACGT")])
    monkeypatch.setitem(sys.modules, "khmer", None)
    with pytest.raises(BackendUnavailableError, match="bio-khmer"):
        query_kmer_counts(path, ["ACG"], verbose=False)
    with pytest.raises(BackendUnavailableError, match="bio-khmer"):
        normalize_reads_by_median(path, tmp_path / "out.fa", verbose=False)
