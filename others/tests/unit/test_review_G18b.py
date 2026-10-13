"""Tests for the G18b review wrappers: sourmash mutation model and k-mer LCA, khmer reads, primer3.

Every wrapped function is checked against a hand calculation or a direct call to the upstream
package. Inputs are written inside pytest's tmp_path; nothing is downloaded and no network is used,
so there is no live test for this group.
"""

from __future__ import annotations

import inspect
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.microbiome import khmer_review_functions as kr
from dnakit.microbiome import sourmash_review_functions as sr
from dnakit.molecular_lab import primer3_review_functions as pr

PAIRED_FASTQ = """@r1/1
ACGTACGT
+
IIIIIIII
@r1/2
TTGGCCAA
+
IIIIIIII
@r2/1
GGGGAAAA
+
IIIIIIII
@r3/1
CCCCTTTT
+
IIIIIIII
@r3/2
AAAACCCC
+
IIIIIIII
"""

PAIRED_FASTA = ">f1/1\nACGTACGT\n>f1/2\nTTGGCCAA\n>f2\nGGGGAAAA\n"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _hand_variance(length: float, k: int, rate: float) -> float:
    """Closed-form variance of the mutated k-mer count, written out term by term."""

    q = 1 - (1 - rate) ** k
    term_1 = length * (1 - q) * (q * (2 * k + 2 / rate - 1) - 2 * k)
    term_2 = k * (k - 1) * (1 - q) ** 2
    term_3 = (2 * (1 - q) / rate**2) * ((1 + (k - 1) * (1 - q)) * rate - q)
    return term_1 + term_2 + term_3


# ---------------------------------------------------------------- public surface


def test_public_signatures_are_keyword_only_verbose() -> None:
    for module in (sr, kr, pr):
        for name in module.__all__:
            parameter = inspect.signature(getattr(module, name)).parameters["verbose"]
            assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
            assert parameter.default is True


def test_every_output_is_json_serializable_with_provenance(tmp_path: Path) -> None:
    pytest.importorskip("sourmash")
    pytest.importorskip("khmer")
    pytest.importorskip("primer3")
    reads = _write(tmp_path / "pairs.fq", PAIRED_FASTQ)
    results = [
        sr.kmer_mutation_probability(0.1, 2, verbose=False),
        sr.expected_mutated_kmers(2, 0.1, n_unique_kmers=100, verbose=False),
        sr.variance_mutated_kmers(2, 0.1, sequence_len_bp=11, verbose=False),
        sr.probability_no_shared_kmers(0.1, 2, 10, n_unique_kmers=100, verbose=False),
        sr.bp_threshold_to_hashes(5000, 1000, 10, verbose=False),
        sr.count_kmer_lca_assignments({1: ["d__Bacteria"]}, verbose=False),
        kr.clean_khmer_read_sequences({"r1": "ACGT"}, verbose=False),
        kr.read_paired_fastx(reads, verbose=False),
        pr.validate_acgt_sequence("ACGT", verbose=False),
    ]
    for result in results:
        assert result["status"] == "completed"
        assert result["operation"]
        assert result["provenance"]
        json.dumps(result)


# ---------------------------------------------------------------- sourmash mutation model


def test_kmer_mutation_probability_matches_hand_and_upstream() -> None:
    distance_utils = pytest.importorskip("sourmash.distance_utils")
    result = sr.kmer_mutation_probability(0.1, 2, verbose=False)
    # 1 - (1 - 0.1) ** 2 = 0.19
    assert result["kmer_mutation_probability"] == pytest.approx(0.19, abs=1e-12)
    assert result["kmer_mutation_probability"] == pytest.approx(
        distance_utils.r1_to_q(2, 0.1), abs=1e-15
    )
    assert sr.kmer_mutation_probability(0.3, 1, verbose=False)["kmer_mutation_probability"] == (
        pytest.approx(0.3)
    )
    assert sr.kmer_mutation_probability(0.0, 5, verbose=False)["kmer_mutation_probability"] == 0.0
    assert sr.kmer_mutation_probability(1.0, 5, verbose=False)["kmer_mutation_probability"] == 1.0


def test_expected_mutated_kmers_matches_hand_and_upstream() -> None:
    distance_utils = pytest.importorskip("sourmash.distance_utils")
    # 100 distinct k-mers, each mutated with probability 0.19
    result = sr.expected_mutated_kmers(2, 0.1, n_unique_kmers=100, verbose=False)
    assert result["expected_mutated_kmers"] == pytest.approx(19.0, abs=1e-9)
    assert result["expected_unmutated_kmers"] == pytest.approx(81.0, abs=1e-9)
    assert result["expected_mutated_kmers"] == pytest.approx(
        distance_utils.exp_n_mutated(100, 2, 0.1), abs=1e-12
    )
    # sequence_len_bp = 101 with k = 2 gives 101 - 1 = 100 distinct k-mers
    by_length = sr.expected_mutated_kmers(2, 0.1, sequence_len_bp=101, verbose=False)
    assert by_length["n_unique_kmers"] == 100
    assert by_length["expected_mutated_kmers"] == pytest.approx(19.0, abs=1e-9)


def test_expected_mutated_kmers_rejects_ambiguous_or_short_input() -> None:
    with pytest.raises(ConfigurationError, match="exactly one"):
        sr.expected_mutated_kmers(2, 0.1, sequence_len_bp=50, n_unique_kmers=49, verbose=False)
    with pytest.raises(ConfigurationError, match="exactly one"):
        sr.expected_mutated_kmers(2, 0.1, verbose=False)
    with pytest.raises(ConfigurationError, match="sequence_len_bp"):
        sr.expected_mutated_kmers(5, 0.1, sequence_len_bp=4, verbose=False)
    with pytest.raises(ConfigurationError, match="mutation_rate"):
        sr.expected_mutated_kmers(2, 1.5, n_unique_kmers=10, verbose=False)


def test_variance_mutated_kmers_matches_hand_formula_and_upstream() -> None:
    distance_utils = pytest.importorskip("sourmash.distance_utils")
    result = sr.variance_mutated_kmers(2, 0.1, sequence_len_bp=11, verbose=False)
    expected = _hand_variance(10, 2, 0.1)  # 2.8512 by hand
    assert expected == pytest.approx(2.8512, abs=1e-9)
    assert result["variance_mutated_kmers"] == pytest.approx(expected, abs=1e-9)
    assert result["variance_mutated_kmers"] == pytest.approx(
        distance_utils.var_n_mutated(10, 2, 0.1), abs=1e-12
    )
    assert result["standard_deviation_mutated_kmers"] == pytest.approx(math.sqrt(expected))
    assert (
        sr.variance_mutated_kmers(2, 0.0, n_unique_kmers=10, verbose=False)[
            "variance_mutated_kmers"
        ]
        == 0.0
    )


def test_probability_no_shared_kmers_matches_hand_formula() -> None:
    distance_utils = pytest.importorskip("sourmash.distance_utils")
    result = sr.probability_no_shared_kmers(0.1, 2, 10, n_unique_kmers=100, verbose=False)
    # exp((100 - 19) * ln(1 - 1/10))
    hand = math.exp((100 - 19) * math.log(0.9))
    assert result["probability_nothing_in_common"] == pytest.approx(hand, rel=1e-12)
    assert result["probability_nothing_in_common"] == pytest.approx(
        distance_utils.get_exp_probability_nothing_common(0.1, 2, 10, n_unique_kmers=100),
        rel=1e-12,
    )
    assert (
        sr.probability_no_shared_kmers(0.0, 2, 10, n_unique_kmers=100, verbose=False)[
            "probability_nothing_in_common"
        ]
        == 0.0
    )
    assert (
        sr.probability_no_shared_kmers(1.0, 2, 10, n_unique_kmers=100, verbose=False)[
            "probability_nothing_in_common"
        ]
        == 1.0
    )
    with pytest.raises(ConfigurationError, match="scaled"):
        sr.probability_no_shared_kmers(0.1, 2, 1, n_unique_kmers=100, verbose=False)


def test_bp_threshold_to_hashes_matches_hand_and_upstream() -> None:
    search = pytest.importorskip("sourmash.search")
    result = sr.bp_threshold_to_hashes(5000, 1000, 10, verbose=False)
    # 5000 bp / scaled 1000 = 5 hashes; 5 / 10 hashes of the query = 0.5
    assert result["threshold_hashes"] == pytest.approx(5.0)
    assert result["containment_threshold"] == pytest.approx(0.5)
    upstream = search.calc_threshold_from_bp(5000, 1000, 10)
    assert (result["containment_threshold"], result["threshold_hashes"]) == pytest.approx(upstream)
    zero = sr.bp_threshold_to_hashes(0, 1000, 10, verbose=False)
    assert (zero["threshold_hashes"], zero["containment_threshold"]) == (0.0, 0.0)
    with pytest.raises(ConfigurationError, match="could not convert"):
        sr.bp_threshold_to_hashes(50_000, 1000, 10, verbose=False)


# ---------------------------------------------------------------- sourmash k-mer LCA


def test_count_kmer_lca_assignments_hand_example_and_order() -> None:
    assignments = {
        1: ["d__Bacteria;p__Firmicutes", "d__Bacteria;p__Bacteroidetes"],
        2: ["d__Bacteria;p__Firmicutes;c__Bacilli"],
        3: ["d__Archaea", "d__Bacteria"],
    }
    result = sr.count_kmer_lca_assignments(assignments, verbose=False)
    assert result["hash_count"] == 3
    assert result["weighted"] is False
    assert result["total_count"] == 3
    # hash 1 -> d__Bacteria, hash 2 -> its full lineage, hash 3 -> root (empty lineage)
    assert [(item["lineage"], item["lowest_rank"], item["count"]) for item in result["counts"]] == [
        ("", None, 1),
        ("d__Bacteria", "superkingdom", 1),
        ("d__Bacteria;p__Firmicutes;c__Bacilli", "class", 1),
    ]
    shuffled = dict(reversed(list(assignments.items())))
    assert sr.count_kmer_lca_assignments(shuffled, verbose=False) == result


def test_count_kmer_lca_assignments_weights_and_upstream() -> None:
    lca_utils = pytest.importorskip("sourmash.lca.lca_utils")
    assignments = {
        1: ["d__Bacteria;p__Firmicutes", "d__Bacteria;p__Bacteroidetes"],
        2: ["d__Bacteria;p__Firmicutes;c__Bacilli"],
        3: ["d__Archaea", "d__Bacteria"],
    }
    result = sr.count_kmer_lca_assignments(
        assignments, hashval_counts={1: 5, 2: 2, 3: 1}, verbose=False
    )
    assert result["weighted"] is True
    assert result["total_count"] == 8
    assert [(item["lineage"], item["count"]) for item in result["counts"]] == [
        ("d__Bacteria", 5),
        ("d__Bacteria;p__Firmicutes;c__Bacilli", 2),
        ("", 1),
    ]
    parsed = {
        key: [lca_utils.make_lineage(text) for text in texts] for key, texts in assignments.items()
    }
    upstream = lca_utils.count_lca_for_assignments(parsed, {1: 5, 2: 2, 3: 1})
    by_name = {
        ";".join(item.name for item in lineage): count for lineage, count in upstream.items()
    }
    assert by_name == {item["lineage"]: item["count"] for item in result["counts"]}


def test_count_kmer_lca_assignments_rejects_bad_input() -> None:
    with pytest.raises(ConfigurationError, match="non-empty mapping"):
        sr.count_kmer_lca_assignments({}, verbose=False)
    with pytest.raises(ConfigurationError, match="no lineage"):
        sr.count_kmer_lca_assignments({1: []}, verbose=False)
    with pytest.raises(ConfigurationError, match="','"):
        sr.count_kmer_lca_assignments({1: ["d__Bacteria,p__Firmicutes"]}, verbose=False)
    with pytest.raises(ConfigurationError, match="integer"):
        sr.count_kmer_lca_assignments({2**64: ["d__Bacteria"]}, verbose=False)
    with pytest.raises(ConfigurationError, match="exactly the hashes"):
        sr.count_kmer_lca_assignments({1: ["d__Bacteria"]}, hashval_counts={2: 1}, verbose=False)
    with pytest.raises(ConfigurationError, match="integer"):
        sr.count_kmer_lca_assignments({1: ["d__Bacteria"]}, hashval_counts={1: 0}, verbose=False)


# ---------------------------------------------------------------- khmer reads


def test_clean_khmer_read_sequences_matches_upstream_and_counts_n() -> None:
    utils = pytest.importorskip("khmer.utils")
    result = kr.clean_khmer_read_sequences({"r2": "ACGT", "r1": "acgNt"}, verbose=False)
    assert [item["id"] for item in result["reads"]] == ["r1", "r2"]
    first = result["reads"][0]
    assert first["sequence"] == "ACGNT"
    assert first["cleaned_sequence"] == "ACGAT"
    assert first["replaced_n_count"] == 1
    assert result["reads_with_n"] == 1
    assert result["replaced_n_total"] == 1
    record = SimpleNamespace(name="r1", sequence="ACGNT")
    upstream = next(iter(utils.clean_input_reads([record]))).cleaned_seq
    assert upstream == first["cleaned_sequence"]


def test_clean_khmer_read_sequences_rejects_other_symbols() -> None:
    with pytest.raises(ConfigurationError, match="only A, C, G, T and N"):
        kr.clean_khmer_read_sequences({"r1": "ACGR"}, verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty mapping"):
        kr.clean_khmer_read_sequences({}, verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty DNA string"):
        kr.clean_khmer_read_sequences({"r1": ""}, verbose=False)


def test_read_paired_fastx_pairs_and_singletons(tmp_path: Path) -> None:
    khmer = pytest.importorskip("khmer")
    utils = pytest.importorskip("khmer.utils")
    path = _write(tmp_path / "pairs.fq", PAIRED_FASTQ)
    result = kr.read_paired_fastx(path, verbose=False)
    assert result["unit_count"] == 3
    assert result["pair_count"] == 2
    assert result["singleton_count"] == 1
    kinds = [(unit["index"], unit["is_pair"]) for unit in result["units"]]
    assert kinds == [(0, True), (2, False), (3, True)]
    first = result["units"][0]
    assert first["read1"] == {"name": "r1/1", "sequence": "ACGTACGT"}
    assert first["read2"] == {"name": "r1/2", "sequence": "TTGGCCAA"}
    assert result["units"][1]["read2"] is None
    upstream = [
        (n, is_pair)
        for n, is_pair, _r1, _r2 in utils.broken_paired_reader(khmer.ReadParser(str(path)))
    ]
    assert upstream == kinds


def test_read_paired_fastx_options_and_fasta(tmp_path: Path) -> None:
    fastq = _write(tmp_path / "pairs.fq", PAIRED_FASTQ)
    fasta = _write(tmp_path / "pairs.fa", PAIRED_FASTA)
    assert kr.read_paired_fastx(fastq, min_length=9, verbose=False)["unit_count"] == 0
    assert kr.read_paired_fastx(fastq, min_length=8, verbose=False)["unit_count"] == 3
    singles = kr.read_paired_fastx(fastq, force_single=True, verbose=False)
    assert (singles["pair_count"], singles["singleton_count"]) == (0, 5)
    fasta_result = kr.read_paired_fastx(fasta, verbose=False)
    assert [unit["is_pair"] for unit in fasta_result["units"]] == [True, False]
    assert fasta_result["units"][0]["read1"]["sequence"] == "ACGTACGT"


def test_read_paired_fastx_errors(tmp_path: Path) -> None:
    fastq = _write(tmp_path / "pairs.fq", PAIRED_FASTQ)
    with pytest.raises(ConfigurationError, match="Unpaired reads"):
        kr.read_paired_fastx(fastq, require_paired=True, verbose=False)
    with pytest.raises(ConfigurationError, match="cannot both"):
        kr.read_paired_fastx(fastq, require_paired=True, force_single=True, verbose=False)
    with pytest.raises(ConfigurationError, match="max_reads"):
        kr.read_paired_fastx(fastq, max_reads=3, verbose=False)
    with pytest.raises(ConfigurationError, match="does not exist"):
        kr.read_paired_fastx(tmp_path / "missing.fq", verbose=False)


# ---------------------------------------------------------------- primer3


def test_validate_acgt_sequence_matches_upstream() -> None:
    primer3 = pytest.importorskip("primer3")
    result = pr.validate_acgt_sequence("acgt", verbose=False)
    assert result["sequence"] == "ACGT"
    assert result["length"] == 4
    assert result["sequence"] == str(primer3.p3helpers.ensure_acgt_uppercase("acgt"))
    assert result["provenance"]["package"] == "primer3-py"


@pytest.mark.parametrize("bad", ["ACGTN", "ACGU", "AC GT", "ACGR", "AC-GT"])
def test_validate_acgt_sequence_rejects_non_acgt(bad: str) -> None:
    pytest.importorskip("primer3")
    with pytest.raises(ConfigurationError, match="only A, C, G and T"):
        pr.validate_acgt_sequence(bad, verbose=False)


def test_validate_acgt_sequence_rejects_empty_and_non_string() -> None:
    with pytest.raises(ConfigurationError, match="non-empty string"):
        pr.validate_acgt_sequence("", verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty string"):
        pr.validate_acgt_sequence(5, verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="non-empty string"):
        pr.validate_acgt_sequence("A" * 100_001, verbose=False)
