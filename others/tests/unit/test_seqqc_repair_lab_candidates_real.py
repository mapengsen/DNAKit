"""Hand-verified tests for the native sequencing_qc, replication_repair and molecular_lab candidates."""

from __future__ import annotations

import itertools
import math
import struct
from pathlib import Path

import pytest

from dnakit import molecular_lab as ml
from dnakit import replication_repair as rr
from dnakit import sequencing_qc as sq
from dnakit.exceptions import ConfigurationError, InputFormatError


def run(fn, value, **options):
    """Call a candidate function, assert an 'ok' envelope and return its result."""

    envelope = fn(value, options=options or None)
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


def fq(*reads: tuple[str, str, str], comment: str = "") -> str:
    return "".join(f"@{n}{(' ' + comment) if comment else ''}\n{s}\n+\n{q}\n" for n, s, q in reads)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# Q40 = 'I', Q20 = '5', Q0 = '!'
FQ_A = fq(
    ("r1", "ACGTACGTAC", "IIIIIIIIII"),
    ("r2", "ACGTACGTAC", "IIIIIIIIII"),
    ("r3", "TTTTTGGGGG", "5555555555"),
    ("r4", "NNACGTACGT", "!!IIIIIIII"),
)


# =========================================================================== sequencing_qc: yield / length / quality


def test_calculate_yield_counts_bases_gc_and_q_fractions():
    r = run(sq.calculate_yield, FQ_A)
    assert (r["reads"], r["bases"], r["n_bases"]) == (4, 40, 2)
    assert r["mean_read_length"] == 10.0
    assert r["gc_fraction"] == 0.5  # 19 G/C among 38 A/C/G/T
    assert r["q20_bases"] == 38 and r["q30_bases"] == 28
    assert r["q20_fraction"] == 0.95 and r["q30_fraction"] == 0.7


def test_calculate_yield_coverage_and_pairs():
    r = run(sq.calculate_yield, {"read1": FQ_A, "read2": FQ_A}, genome_size=20)
    assert r["paired"] and r["read_pairs"] == 4 and r["total_bases"] == 80
    assert r["expected_coverage"] == 4.0


def test_calculate_read_length_n50_and_l50():
    reads = ["A" * n for n in (100, 200, 300, 400, 500)]
    r = run(sq.calculate_read_length, reads)
    assert (r["min"], r["max"], r["mean"], r["median"]) == (100, 500, 300.0, 300.0)
    assert r["n50"] == 400  # 500 + 400 = 900 >= 750
    assert r["n90"] == 200  # 500 + 400 + 300 + 200 = 1400 >= 1350
    assert r["l50"] == 2
    assert r["stdev"] == pytest.approx(158.113883)


def test_calculate_quality_per_position_percentiles():
    r = run(sq.calculate_quality, FQ_A)
    first = r["per_position"][0]  # qualities 40, 40, 20, 0
    assert first["mean"] == 25.0
    assert first["median"] == 30.0  # sorted 0,20,40,40 -> position 1.5 -> 20 + 0.5*20
    assert first["q1"] == 15.0 and first["q3"] == 40.0
    assert (first["min"], first["max"]) == (0, 40)
    per_read = r["per_read_mean_quality"]  # read means 40, 40, 20, 32
    assert per_read["distribution"] == {"20": 1, "32": 1, "40": 2}
    assert per_read["mode"] == 40 and per_read["mean"] == 33.0


def test_calculate_qscore_summary_both_read_q_methods():
    r = run(sq.calculate_qscore_summary, FQ_A)
    assert r["mean_q"] == 33.0 and r["median_q"] == 40.0
    assert r["fraction_bases_at_least"] == {"Q10": 0.95, "Q20": 0.95, "Q30": 0.7, "Q40": 0.7}
    # read r4 has mean error probability (2*1 + 8*1e-4)/10 = 0.20008 -> Q 6.988 < 7
    assert r["reads_with_mean_q_at_least"]["Q7"] == 3
    assert r["reads_with_mean_q_at_least"]["Q30"] == 2
    arithmetic = run(sq.calculate_qscore_summary, FQ_A, read_q_method="arithmetic")
    assert arithmetic["reads_with_mean_q_at_least"]["Q7"] == 4  # r4 arithmetic mean is 32
    assert arithmetic["mean_read_q"] == 33.0


def test_calculate_basecalling_qc_n_fraction_errors_and_pass_filter():
    r = run(sq.calculate_basecalling_qc, FQ_A)
    assert r["no_call_fraction"] == 0.05 and r["reads_with_no_calls"] == 1
    assert r["q30_fraction"] == 0.7
    assert r["expected_errors_per_read"] == pytest.approx((0.001 + 0.001 + 0.1 + 2.0008) / 4)
    assert r["reads_below_min_mean_q"] == 0
    assert (
        run(sq.calculate_basecalling_qc, FQ_A, min_mean_quality=25)["reads_below_min_mean_q"] == 1
    )
    assert r["pass_filter_fraction"] is None
    flagged = fq(("a", "ACGT", "IIII"), comment="1:N:0:AC") + fq(
        ("b", "ACGT", "IIII"), comment="1:Y:0:AC"
    )
    pf = run(sq.calculate_basecalling_qc, flagged)
    assert pf["pass_filter_fraction"] == 0.5 and pf["filtered_reads"] == 1


def test_calculate_error_profile_is_the_mean_phred_error_probability():
    r = run(sq.calculate_error_profile, FQ_A)
    assert r["per_position"][0]["expected_error_rate"] == pytest.approx(
        (1e-4 + 1e-4 + 1e-2 + 1.0) / 4
    )
    assert r["per_position"][2]["expected_error_rate"] == pytest.approx((3e-4 + 1e-2) / 4)
    assert r["overall_expected_error_rate"] == pytest.approx(2.1028 / 40)
    assert r["per_position"][0]["n_fraction"] == 0.25


def test_calculate_gc_bias_histogram_and_reference_shift():
    r = run(sq.calculate_gc_bias, ["GGGG", "AAAA", "GCAT", "GGAT"], reference_gc=0.4)
    assert r["mean_gc_percent"] == 50.0 and r["mode_gc_percent"] == 50
    assert r["histogram"] == {"0": 1, "50": 2, "100": 1}
    assert r["stdev_gc_percent"] == pytest.approx(math.sqrt(5000 / 3))
    assert r["gc_percent_by_position"][0] == 75.0  # G, A, G, G
    assert r["gc_shift_vs_reference_percent"] == pytest.approx(10.0)
    assert r["deviation_from_normal_percent"] is not None


# =========================================================================== adapters, duplication, complexity


def test_calculate_adapter_content_cumulative_percentage():
    reads = ["AAAAAGATCGGAAGAGCTT", "A" * 19]  # AGATCGGAAGAG starts at index 4 -> position 5
    r = run(sq.calculate_adapter_content, reads)
    uni = r["adapters"]["Illumina Universal Adapter"]
    assert uni["reads_with_adapter"] == 1 and uni["fraction_reads_with_adapter"] == 0.5
    assert uni["cumulative_percent_by_position"][3] == 0.0
    assert uni["cumulative_percent_by_position"][4] == 50.0
    assert r["reads_with_any_adapter"] == 1


def test_calculate_adapter_content_mismatch_tolerance_and_custom_adapters():
    reads = ["CCCCAGATCGGTAGAGCCCC"]  # one mismatch (A->T) inside the 12-mer
    assert (
        run(sq.calculate_adapter_content, reads, adapters={"x": "AGATCGGAAGAG"})[
            "reads_with_any_adapter"
        ]
        == 0
    )
    r = run(sq.calculate_adapter_content, reads, adapters={"x": "AGATCGGAAGAG"}, max_mismatches=1)
    assert r["adapters"]["x"]["cumulative_percent_by_position"][4] == 100.0


def test_calculate_duplicate_rate_variants():
    assert run(sq.calculate_duplicate_rate, FQ_A)["duplicate_rate"] == 0.25
    reads = ["ACGTACGTAC", "GTACGTACGT", "ACGTACGTAC"]  # 2nd is the reverse complement of the 1st
    assert run(sq.calculate_duplicate_rate, reads)["duplicate_rate"] == pytest.approx(1 / 3)
    assert run(sq.calculate_duplicate_rate, reads, reverse_complement_aware=True)[
        "duplicate_rate"
    ] == pytest.approx(2 / 3)
    assert (
        run(sq.calculate_duplicate_rate, ["AAAAC", "AAAAG"], prefix_length=4)["duplicate_rate"]
        == 0.5
    )
    pair = run(sq.calculate_duplicate_rate, {"read1": ["AAAA", "AAAA"], "read2": ["CCCC", "GGGG"]})
    assert pair["paired"] and pair["duplicate_rate"] == 0.0


def test_calculate_duplication_levels():
    reads = ["AAAA"] * 3 + ["CCCC", "GGGG"]
    r = run(sq.calculate_duplication, reads)
    levels = {row["level"]: row for row in r["levels"]}
    assert r["percent_remaining_after_deduplication"] == 60.0
    assert levels["3"]["percent_of_total_reads"] == 60.0
    assert levels["3"]["percent_of_deduplicated"] == pytest.approx(100 / 3, abs=1e-5)
    assert levels["1"]["percent_of_total_reads"] == 40.0


def test_estimate_library_complexity_lander_waterman():
    r = run(sq.estimate_library_complexity, {"total_reads": 1_000_000, "unique_reads": 632_121})
    x = r["estimated_library_size"]
    # N = X -> unique = X (1 - e^-1) = 0.632121 X, so X ~ 1e6
    assert x == pytest.approx(1_000_000, rel=2e-6)
    assert 632_121 / x == pytest.approx(1 - math.exp(-1_000_000 / x), rel=1e-9)
    assert r["predicted_unique_reads"]["2x"] == pytest.approx(x * (1 - math.exp(-2e6 / x)))
    assert r["predicted_unique_reads"]["2x"] == pytest.approx(
        1_000_000 * (1 - math.exp(-2)), rel=1e-5
    )
    from_reads = run(sq.estimate_library_complexity, ["AAAA", "AAAA", "CCCC", "GGGG"])
    assert from_reads["total_reads"] == 4 and from_reads["unique_reads"] == 3
    assert 3 / from_reads["estimated_library_size"] == pytest.approx(
        1 - math.exp(-4 / from_reads["estimated_library_size"])
    )
    all_unique = fq(("a", "ACGT", "IIII"), ("b", "TTGA", "IIII"))
    assert sq.estimate_library_complexity(all_unique)["result"]["estimated_library_size"] is None


# =========================================================================== barcodes, demultiplexing, indexes


def test_calculate_barcode_quality_whitelist_and_entropy():
    reads = fq(
        ("a", "ACGTTTTT", "IIII####"),
        ("b", "ACGATTTT", "IIII####"),
        ("c", "TTTTAAAA", "IIII####"),
        ("d", "ACGTCCCC", "IIII####"),
    )
    r = run(sq.calculate_barcode_quality, reads, barcode_length=4, whitelist=["ACGT", "TTTT"])
    assert r["distinct_barcodes"] == 3
    assert r["barcode_entropy_bits"] == pytest.approx(1.5)  # counts 2, 1, 1
    assert r["mean_quality"] == 40.0 and r["q30_fraction"] == 1.0
    wl = r["whitelist"]
    assert (wl["exact_match_fraction"], wl["corrected_fraction"], wl["invalid_fraction"]) == (
        0.75,
        0.25,
        0.0,
    )
    r0 = run(
        sq.calculate_barcode_quality,
        reads,
        barcode_length=4,
        whitelist=["ACGT", "TTTT"],
        max_mismatches=0,
    )
    assert r0["whitelist"]["invalid_fraction"] == 0.25


def test_assign_sample_barcodes_mismatch_tolerance_and_ties():
    reads = [("1", "ACGTAAAA"), ("2", "ACGAAAAA"), ("3", "TTTTAAAA"), ("4", "GGGGAAAA")]
    r = run(
        sq.assign_sample_barcodes, reads, samples={"S1": "ACGT", "S2": "TTTT"}, barcode_length=4
    )
    assert [a["sample"] for a in r["assignments"]] == ["S1", "S1", "S2", None]
    assert [a["distance"] for a in r["assignments"]] == [0, 1, 0, 3]
    assert r["reads_per_sample"] == {"S1": 2, "S2": 1} and r["unassigned"] == 1
    assert r["min_distance_between_sample_barcodes"] == 3
    tie = run(
        sq.assign_sample_barcodes,
        [("t", "AAATCCCC")],
        samples={"A": "AAAA", "B": "AATT"},
        barcode_length=4,
    )
    assert tie["ambiguous"] == 1 and tie["assignments"][0]["status"] == "ambiguous"
    strict = run(
        sq.assign_sample_barcodes,
        reads,
        samples={"S1": "ACGT", "S2": "TTTT"},
        barcode_length=4,
        max_mismatches=0,
    )
    assert strict["reads_per_sample"] == {"S1": 1, "S2": 1}


def test_assign_sample_barcodes_from_header_index():
    reads = fq(("a", "AAAA", "IIII"), comment="1:N:0:ACGT+TTGG") + fq(
        ("b", "AAAA", "IIII"), comment="1:N:0:TTGG+ACGT"
    )
    r = run(
        sq.assign_sample_barcodes,
        reads,
        samples={"S1": "ACGT+TTGG", "S2": "TTGG+ACGT"},
        barcode_location="header",
    )
    assert r["reads_per_sample"] == {"S1": 1, "S2": 1}


def test_demultiplex_barcoded_reads_trims_and_writes(tmp_path):
    reads = fq(
        ("a", "ACGTAAAA", "IIIIJJJJ"), ("b", "TTTTCCCC", "IIIIJJJJ"), ("c", "GGGGCCCC", "IIIIJJJJ")
    )
    r = run(
        sq.demultiplex_barcoded_reads,
        reads,
        samples={"S1": "ACGT", "S2": "TTTT"},
        barcode_length=4,
        output_directory=str(tmp_path),
    )
    assert r["samples"]["S1"]["records"][0]["sequence"] == "AAAA"
    assert r["samples"]["S1"]["records"][0]["quality"] == "JJJJ"
    assert r["samples"]["unassigned"]["records"][0]["sequence"] == "GGGGCCCC"  # not trimmed
    assert (tmp_path / "S1.fastq").read_text() == "@a\nAAAA\n+\nJJJJ\n"
    with pytest.raises(ConfigurationError):  # refuses to overwrite
        sq.demultiplex_barcoded_reads(
            reads,
            options={
                "samples": {"S1": "ACGT"},
                "barcode_length": 4,
                "output_directory": str(tmp_path),
            },
        )


def test_assign_sample_barcodes_read_end_location():
    r = run(
        sq.assign_sample_barcodes,
        [("a", "CCCCACGT"), ("b", "CCCCTTTA")],
        samples={"S1": "ACGT", "S2": "TTTT"},
        barcode_length=4,
        barcode_location="read_end",
    )
    assert [a["sample"] for a in r["assignments"]] == ["S1", "S2"] and r["assignments"][1][
        "distance"
    ] == 1


def test_calculate_index_hopping_rate():
    headers = (
        [("s1", "A", "I", "1:N:0:AAAA+CCCC")] * 6
        + [("s2", "A", "I", "1:N:0:GGGG+TTTT")] * 2
        + [("h", "A", "I", "1:N:0:AAAA+TTTT")]
        + [("u", "A", "I", "1:N:0:CCCC+AAAA")]
    )
    text = "".join(f"@{n} {c}\n{s}\n+\n{q}\n" for n, s, q, c in headers)
    samples = {"S1": "AAAA+CCCC", "S2": "GGGG+TTTT"}
    for fn in (sq.calculate_index_hopping, sq.calculate_index_hopping_tenx):
        r = run(fn, text, samples=samples)
        assert (
            r["expected_combination_reads"],
            r["index_hopped_reads"],
            r["unknown_index_reads"],
        ) == (8, 1, 1)
        assert r["index_hopping_rate"] == pytest.approx(1 / 9, abs=1e-6)
        assert r["unknown_index_fraction"] == 0.1


# =========================================================================== insert size


S30 = "GATTACAGCTAGGCTTACGGATCCATGCAA"


def test_calculate_insert_size_from_read_pair_overlap():
    pair = {"read1": [S30[:20]], "read2": [revcomp(S30[10:])]}  # fragment length 30, overlap 10
    r = run(sq.calculate_insert_size, pair, min_overlap=10, max_mismatch_rate=0.0)
    assert r["count"] == 1 and r["median"] == 30.0
    no_overlap = run(
        sq.calculate_insert_size,
        {
            "read1": [S30[:20], "ACGTACGTAC" + "GGGGGGGGGG"],
            "read2": [revcomp(S30[10:]), "TTTTTTTTTTTTTTTTTTTT"],
        },
        min_overlap=10,
        max_mismatch_rate=0.0,
    )
    assert no_overlap["pairs_without_overlap"] == 1 and no_overlap["count"] == 1


def test_calculate_insert_size_from_sam_and_fragments_and_numbers():
    sam = "\n".join(
        [
            "@HD\tVN:1.6",
            "p1\t99\tchr1\t100\t60\t50M\t=\t250\t200\t*\t*",
            "p1\t147\tchr1\t250\t60\t50M\t=\t100\t-200\t*\t*",
            "p2\t99\tchr1\t100\t60\t50M\t=\t350\t300\t*\t*",
            "p3\t99\tchr1\t100\t60\t50M\t=\t450\t400\t*\t*",
            "p4\t97\tchr1\t100\t60\t50M\t=\t450\t999\t*\t*",  # not a proper pair -> skipped
        ]
    )
    r = run(sq.calculate_insert_size, sam)
    assert (r["count"], r["min"], r["max"], r["mean"], r["median"], r["stdev"], r["mode"]) == (
        3,
        200,
        400,
        300.0,
        300.0,
        100.0,
        200,
    )
    assert r["non_proper_pairs_skipped"] == 1
    fragments = "chr1\t100\t250\tAAAC\t2\nchr1\t100\t300\tAAAC\t1\n"
    for fn in (sq.calculate_insert_size, sq.calculate_insert_size_tenx):
        f = run(fn, fragments)
        assert (
            f["count"] == 3
            and f["median"] == 150.0
            and f["mean"] == pytest.approx(500 / 3, abs=1e-5)
        )
    assert run(sq.calculate_insert_size, [100, 100, 200, 300])["median"] == 150.0


def test_calculate_insert_size_distribution_histogram():
    r = run(sq.calculate_insert_size_distribution, [100, 100, 200, 300], bin_width=100)
    assert [(h["start"], h["count"], h["cumulative_fraction"]) for h in r["histogram"]] == [
        (100, 2, 0.5),
        (200, 1, 0.75),
        (300, 1, 1.0),
    ]


# =========================================================================== contamination, platform, errors


REF = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"


def test_calculate_contamination_kmer_screen():
    other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
    reads = [REF[:30], REF[20:50], revcomp(REF[10:40]), other[:30]]
    r = run(
        sq.calculate_contamination,
        reads,
        references={"host": REF, "bug": other},
        kmer_size=12,
        expected_reference="host",
    )
    assert r["hits"]["host"]["reads"] == 3 and r["hits"]["bug"]["reads"] == 1
    assert r["contamination_fraction"] == 0.25
    assert r["reads_without_hit"] == 0


def test_detect_cross_sample_contamination_containment():
    other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
    samples = {"A": [REF[:40]], "B": [REF[20:60]], "C": [other[:40]]}
    r = run(sq.detect_cross_sample_contamination, samples, kmer_size=12, threshold=0.2)
    # A has 29 12-mers; A and B share REF[20:40] -> 9 12-mers
    assert r["containment"]["A"]["B"] == pytest.approx(9 / 29, abs=1e-5)
    assert r["containment"]["A"]["C"] == 0.0
    assert {(f["sample"], f["found_in"]) for f in r["flagged_pairs"]} == {("A", "B"), ("B", "A")}


def test_estimate_read_error_model_against_reference():
    def mutate(base: str) -> str:
        return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]

    r1 = REF[:30]
    r1 = r1[:10] + mutate(r1[10]) + r1[11:]
    r2 = REF[20:50]
    seg = REF[40:70]
    original = revcomp(seg)
    flipped = (
        mutate(original[0]) + original[1:]
    )  # mismatch at read position 0 (reverse-strand read)
    out = run(sq.estimate_read_error_model, [r1, r2, flipped], reference=REF)["model"]
    assert out["reads_aligned"] == 3 and out["bases_compared"] == 90 and out["mismatches"] == 2
    assert out["overall_mismatch_rate"] == pytest.approx(2 / 90, abs=1e-6)
    by_pos = {row["position"]: row for row in out["per_position_mismatch_rate"]}
    assert by_pos[11]["mismatch_rate"] == pytest.approx(1 / 3, abs=1e-6)
    assert by_pos[1]["mismatch_rate"] == pytest.approx(1 / 3, abs=1e-6)
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverse_read_base = complement[
        mutate(original[0])
    ]  # the changed base seen on the reference strand
    assert out["substitution_counts"] == {
        f"{REF[10]}>{mutate(REF[10])}": 1,
        f"{REF[69]}>{reverse_read_base}": 1,
    }
    phred_only = run(sq.estimate_read_error_model, FQ_A)["model"]
    assert phred_only["source"] == "phred_quality"


def test_detect_platform_from_read_names():
    cases = {
        "@A00123:45:HXYZ123DSXX:1:1101:1000:2000 1:N:0:ACGT": "illumina",
        "@3a8e0d4c-1b2e-4f1c-9d7a-2c5b8e9f0a11 runid=abc ch=12": "oxford_nanopore",
        "@m64011_190830_220126/1/ccs": "pacbio",
        "@V300012345L1C001R0010000001/1": "mgi_dnbseq",
        "@ABCDE:12345:67890": "ion_torrent",
        "@SRR1234567.1 length=36": "sra_renamed",
    }
    for header, platform in cases.items():
        text = f"{header}\nACGT\n+\nIIII\n"
        r = run(sq.detect_platform, text)
        assert r["platform"] == platform, header
    illumina = run(
        sq.detect_platform, "@A00123:45:HXYZ123DSXX:1:1101:1000:2000 1:N:0:ACGT\nACGT\n+\nIIII\n"
    )
    assert illumina["evidence"] == {
        "instrument_id": "A00123",
        "run_number": "45",
        "flowcell": "HXYZ123DSXX",
        "lane": "1",
    }


# =========================================================================== UMIs


def test_correct_umi_errors_methods():
    counts = {"AAAAAA": 100, "AAAAAT": 1, "CCCCCC": 50, "CCCCCG": 30}
    directional = run(sq.correct_umi_errors, counts)
    # 100 >= 2*1-1 merges AAAAAT; 50 >= 2*30-1 = 59 is false, so CCCCCG stays separate
    assert directional["distinct_umis_after"] == 3 and directional["corrections"] == {
        "AAAAAT": "AAAAAA"
    }
    assert run(sq.correct_umi_errors, counts, method="cluster")["distinct_umis_after"] == 2
    adjacency = run(sq.correct_umi_errors, counts, method="adjacency")
    assert adjacency["distinct_umis_after"] == 2 and adjacency["corrections"]["CCCCCG"] == "CCCCCC"
    assert run(sq.correct_umi_errors, ["AAAA", "AAAA", "AAAA", "AAAT"])["corrections"] == {
        "AAAT": "AAAA"
    }
    assert run(sq.correct_umi_errors, counts, threshold=0)["distinct_umis_after"] == 4
    # chain A(100) - B(60) - C(30): 100 >= 2*60-1 is false (A stays alone) but 60 >= 2*30-1 merges C into B; C is 2 away from A
    chain = run(sq.correct_umi_errors, {"AAAA": 100, "AAAT": 60, "AATT": 30})
    assert chain["distinct_umis_after"] == 2 and chain["corrections"] == {"AATT": "AAAT"}
    assert (
        run(sq.correct_umi_errors, {"AAAA": 100, "AAAT": 60, "AATT": 30}, method="cluster")[
            "distinct_umis_after"
        ]
        == 1
    )


def test_collapse_umi_reads_consensus_and_grouping():
    reads = fq(
        ("r1_AAAA", "ACGTACGT", "IIIIIIII"),
        ("r2_AAAA", "ACGTACGT", "IIIIIIII"),
        ("r3_AAAA", "ACGAACGT", "III5IIII"),
        ("r4_CCCC", "ACGTACGT", "IIIIIIII"),
    )
    r = run(sq.collapse_umi_reads, reads, group_by="umi")
    assert r["reads_in"] == 4 and r["molecules_out"] == 2 and r["collapse_ratio"] == 2.0
    group = next(g for g in r["groups"] if g["umi"] == "AAAA")
    assert group["reads"] == 3 and group["consensus_sequence"] == "ACGTACGT"
    # position 4: T (40+40) against A (20): 80 - 20 = 60 -> chr(33 + 60); elsewhere 120 -> capped at 93
    assert group["consensus_quality"] == "~~~]~~~~"
    split = run(sq.collapse_umi_reads, reads, group_by="umi_and_start", start_length=4)
    assert split["molecules_out"] == 3  # r3 starts with ACGA


def test_extract_umi_sequences_length_and_pattern():
    r = run(sq.extract_umi_sequences, fq(("r1", "ACGTAAAATTTT", "ABCDEFGHIJKL")), umi_length=4)
    rec = r["records"][0]
    assert (rec["id"], rec["umi"], rec["sequence"], rec["quality"]) == (
        "r1_ACGT",
        "ACGT",
        "AAAATTTT",
        "EFGHIJKL",
    )
    p = run(sq.extract_umi_sequences, ["ACGTACGGGG"], pattern="CCNNXX")
    assert (
        p["records"][0]["cell_barcode"] == "AC"
        and p["records"][0]["umi"] == "GT"
        and p["records"][0]["sequence"] == "GGGG"
    )
    three = run(sq.extract_umi_sequences, ["AAAAAAGT"], umi_length=2, umi_location="3prime")
    assert three["records"][0]["umi"] == "GT" and three["records"][0]["sequence"] == "AAAAAA"


# =========================================================================== reports and plots


def test_generate_qc_report_flags_follow_fastqc_thresholds():
    r = run(sq.generate_qc_report, FQ_A)
    flags = r["sections"]["reads"]["flags"]
    assert flags["per_base_n_content"] == "fail"  # 25% N at positions 1-2 (> 20%)
    assert flags["sequence_duplication"] == "warn"  # 25% non-unique (> 20%, <= 50%)
    assert flags["per_base_quality"] == "pass" and flags["adapter_content"] == "pass"
    assert "per_base_n_content: fail" in r["report_text"]


def test_generate_qc_report_writes_only_on_request(tmp_path):
    target = tmp_path / "report.txt"
    r = run(sq.generate_qc_report, FQ_A, output_path=str(target))
    assert target.read_text() == r["report_text"]
    with pytest.raises(ConfigurationError):
        sq.generate_qc_report(FQ_A, options={"output_path": str(target)})
    assert run(sq.generate_qc_report, FQ_A, output_path=str(target), overwrite=True)[
        "output_path"
    ] == str(target)


def test_plot_functions_return_svg(tmp_path):
    q = run(sq.plot_per_base_quality, FQ_A)
    assert q["svg"].startswith("<svg") and q["svg"].count("<rect") >= 10 and q["positions"] == 10
    s = run(sq.plot_qc_summary, FQ_A)
    assert s["metrics"]["bases >= Q30"] == 70.0 and s["metrics"]["duplicated reads"] == 25.0
    assert "<svg" in s["svg"] and "70.00%" in s["svg"]
    out = tmp_path / "q.svg"
    assert run(sq.plot_per_base_quality, FQ_A, output_path=str(out))["output_path"] == str(out)
    assert out.read_text().startswith("<svg")


# =========================================================================== 10x variants

BC = "ACGTACGTACGTACGT"
R1_A = fq(
    ("p1", BC + "AAAAAAAAAAAA", "I" * 28),
    ("p2", BC + "AAAAAAAAAAAA", "I" * 28),
    ("p3", "TTTTTTTTTTTTTTTT" + "CCCCCCCCCCCC", "I" * 16 + "5" * 12),
)
R2_A = fq(
    ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
    ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
    ("p3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", "5" * 30),
)
TENX = {"read1": R1_A, "read2": R2_A}


def test_tenx_yield_regions_and_valid_barcodes():
    r = run(sq.calculate_yield_tenx, TENX, whitelist=[BC])
    assert r["chemistry"] == "v3" and (r["barcode_length"], r["umi_length"]) == (16, 12)
    assert (
        r["regions"]["barcode"]["bases"] == 48
        and r["regions"]["umi"]["bases"] == 36
        and r["regions"]["rna_read"]["bases"] == 90
    )
    assert r["distinct_barcodes"] == 2 and r["median_reads_per_barcode"] == 1.5
    assert r["valid_barcode_fraction"] == pytest.approx(2 / 3, abs=1e-6)


def test_tenx_quality_basecalling_and_error_profile():
    q = run(sq.calculate_quality_tenx, TENX)
    assert q["regions"]["umi"]["per_position"][0]["mean"] == pytest.approx((40 + 40 + 20) / 3)
    b = run(sq.calculate_basecalling_qc_tenx, TENX)
    assert b["regions"]["barcode"]["q30_fraction"] == 1.0
    assert b["regions"]["umi"]["q30_fraction"] == pytest.approx(
        24 / 36
    )  # p3's 12 UMI bases are Q20
    assert b["regions"]["rna_read"]["q30_fraction"] == pytest.approx(60 / 90)
    e = run(sq.calculate_error_profile_tenx, TENX)
    assert e["regions"]["rna_read"]["per_position"][0]["expected_error_rate"] == pytest.approx(
        (1e-4 + 1e-4 + 1e-2) / 3
    )


def test_tenx_barcode_quality_umis_and_read_length():
    r = run(sq.calculate_barcode_quality_tenx, TENX, whitelist=[BC])
    assert r["cell_barcode"]["whitelist"]["exact_match_fraction"] == pytest.approx(2 / 3, abs=1e-6)
    assert r["umi"]["homopolymer_umi_reads"] == 3 and r["umi"]["distinct_umis"] == 2
    length = run(sq.calculate_read_length_tenx, TENX)
    assert length["read1"]["median"] == 28.0 and length["read2"]["median"] == 30.0


def test_tenx_duplication_adapter_gc_contamination():
    d = run(sq.calculate_duplication_tenx, TENX)
    assert d["duplicate_rate"] == pytest.approx(1 / 3, abs=1e-6)
    a = run(sq.calculate_adapter_content_tenx, TENX)
    assert a["adapters"]["Illumina Universal Adapter"]["reads_with_adapter"] == 2
    assert a["adapters"]["Illumina Universal Adapter"]["cumulative_percent_by_position"][
        8
    ] == pytest.approx(200 / 3, abs=1e-5)  # hit at index 8
    g = run(sq.calculate_gc_bias_tenx, TENX)
    assert g["histogram"]["100"] == 1  # p3 is all C
    host = "GATTACAGATCGGAAGAGCGGTTCAGCAGGTTACCAGTCAGGCCTAAGGTCCATTGGCC"
    c = run(sq.calculate_contamination_tenx, TENX, references={"host": host}, kmer_size=12)
    assert c["hits"]["host"]["reads"] == 2 and c["reads_without_hit"] == 1


def test_tenx_platform_report_and_plot():
    named = {
        "read1": fq(("A00123:45:HXYZ123DSXX:1:1101:1000:2000", BC + "AAAAAAAAAAAA", "I" * 28)),
        "read2": fq(("A00123:45:HXYZ123DSXX:1:1101:1000:2000", "ACGT" * 5, "I" * 20)),
    }
    p = run(sq.detect_platform_tenx, named)
    assert (
        p["platform"] == "illumina"
        and p["tenx"]["read1_modal_length"] == 28
        and p["tenx"]["layout"].startswith("16 bp barcode + 12 bp UMI")
    )
    r = run(sq.generate_qc_report_tenx, TENX)
    assert r["yield"]["distinct_barcodes"] == 2 and "Q30 fraction" in r["report_text"]
    plot = run(sq.plot_qc_summary_tenx, TENX)
    assert plot["metrics"]["Q30 bases in barcode"] == 100.0 and plot["svg"].startswith("<svg")
    v2 = run(
        sq.calculate_yield_tenx,
        {"read1": fq(("a", "A" * 26, "I" * 26)), "read2": fq(("a", "ACGT", "IIII"))},
        chemistry="v2",
    )
    assert v2["regions"]["umi"]["bases"] == 10


# =========================================================================== sequencing_qc: leftovers and errors


@pytest.mark.parametrize(
    "name",
    [
        "basecall_nanopore_reads",
        "calculate_phasing_qc",
        "calculate_phasing_qc_tenx",
        "correct_long_reads",
    ],
)
def test_sequencing_qc_leftovers_are_not_implemented(name):
    envelope = getattr(sq, name)("ACGT", options=None)
    assert envelope["status"] == "not_implemented"
    assert envelope["required_dependency"] and envelope["message"]
    assert "result" not in envelope


@pytest.mark.parametrize(
    ("fn", "value", "options", "error"),
    [
        (sq.calculate_yield, "", None, InputFormatError),
        (
            sq.calculate_yield,
            "@r1\nACGT\n+\nIII\n",
            None,
            InputFormatError,
        ),  # quality shorter than sequence
        (sq.calculate_yield, "@r1\nACGT\n+\n", None, InputFormatError),  # incomplete record
        (sq.calculate_yield, 12345, None, InputFormatError),
        (sq.calculate_quality, ["ACGT"], None, InputFormatError),  # no qualities
        (sq.calculate_yield, ["ACGT"], {"phred_offset": 40}, ConfigurationError),
        (sq.calculate_yield, ["ACGT"], "not-a-mapping", ConfigurationError),
        (
            sq.calculate_yield,
            {"read1": ["ACGT"], "read2": ["ACGT", "ACGT"]},
            None,
            InputFormatError,
        ),
        (
            sq.calculate_contamination,
            ["ACGTACGTACGTACGTACGTACGT"],
            None,
            ConfigurationError,
        ),  # references required
        (sq.assign_sample_barcodes, ["ACGT"], None, ConfigurationError),  # samples required
        (sq.calculate_yield_tenx, ["ACGT"], None, InputFormatError),  # needs a pair
        (
            sq.calculate_yield_tenx,
            {"read1": ["ACGT"], "read2": ["ACGT"]},
            None,
            InputFormatError,
        ),  # read 1 too short
        (sq.detect_cross_sample_contamination, ["ACGT"], None, InputFormatError),
        (sq.calculate_insert_size, "this is not sam", None, InputFormatError),
        (
            sq.estimate_library_complexity,
            {"total_reads": 5, "unique_reads": 9},
            None,
            ConfigurationError,
        ),
        (sq.calculate_adapter_content, ["ACGT"], {"adapters": {"x": "ACGU!"}}, ConfigurationError),
        (
            sq.extract_umi_sequences,
            ["ACGT"],
            {"umi_length": 4},
            InputFormatError,
        ),  # read not longer than the UMI
        (sq.plot_per_base_quality, ["ACGT"], None, InputFormatError),
    ],
)
def test_sequencing_qc_invalid_input_raises(fn, value, options, error):
    with pytest.raises(error):
        fn(value, options=options)


def test_sequencing_qc_reads_a_fastq_file(tmp_path):
    path = tmp_path / "x.fastq"
    path.write_text(FQ_A)
    assert run(sq.calculate_yield, str(path))["reads"] == 4
    import gzip

    gz = tmp_path / "x.fastq.gz"
    gz.write_bytes(gzip.compress(FQ_A.encode()))
    assert run(sq.calculate_yield, gz)["bases"] == 40


# =========================================================================== replication_repair


def two_replichore(n: int = 1000) -> str:
    return "GGGC" * (n // 8) + "CCCG" * (n // 8)


def test_replication_calculate_gc_skew_and_cumulative():
    r = run(rr.calculate, "GGGCCCCG", window=4)
    assert [w["skew"] for w in r["window_table"]] == [0.5, -0.5]
    assert [w["cumulative_skew"] for w in r["window_table"]] == [0.5, 0.0]
    assert r["origin_window"] == [4, 8] and r["terminus_window"] == [0, 4]
    assert r["overall_skew"] == 0.0
    at = run(rr.calculate, "AAATTTTA", window=4, skew_type="AT")
    assert [w["skew"] for w in at["window_table"]] == [0.5, -0.5]


def test_replication_calculate_origin_and_terminus_from_cumulative_skew():
    r = run(rr.calculate, two_replichore(), window=100)
    # windows 0-4 have skew +0.5 (cumulative rises to 2.5), windows 5-9 have -0.5 (back to 0)
    assert r["terminus_window"] == [400, 500] and r["terminus_estimate"] == 450
    assert r["origin_window"] == [900, 1000] and r["origin_estimate"] == 950
    assert r["cumulative_skew_range"] == 2.5
    assert "heuristic" in r and "estimate" in r["heuristic"]
    assert run(rr.calculate, ">chr\nGGGC\nCCCG\n", window=4)["windows"] == 2


def test_replication_predict_dnaa_box_cluster():
    filler = "CGGC" * 150
    cluster = "TTATCCACAGG" * 5
    seq = filler[:600] + cluster + filler[600:]
    r = run(rr.predict, seq, window=100)
    best = r["candidate_regions"][0]
    assert best["dnaa_boxes"] == 5 and best["start"] == 600
    assert r["boxes_found"] == 5
    near = filler[:300] + "TTATCCACG" + filler[300:]  # one mismatch only
    assert run(rr.predict, near)["boxes_found"] == 0
    assert run(rr.predict, near, max_mismatches=1)["boxes_found"] == 1
    rc = "CGGC" * 20 + "TGTGGATAA" + "CGGC" * 20  # reverse-strand box
    hit = run(rr.predict, rc)["candidate_regions"][0]["box_positions"]
    assert hit == [(80, "-")]


def test_replication_detect_gatc_and_microsatellites():
    seq = "AAGATCTT" + "A" * 10 + "G" + "CA" * 6 + "GATC" + "TTT"
    r = run(rr.detect, seq)
    assert r["gatc_sites"]["positions"] == [2, seq.rindex("GATC")]
    tracts = {(t["unit"], t["start"], t["copies"]) for t in r["microsatellite_tracts"]["tracts"]}
    assert tracts == {("A", 8, 10), ("CA", 19, 6)}
    assert run(rr.detect, seq, min_copies={1: 12})["microsatellite_tracts"]["count"] == 1


def test_replication_annotate_leading_and_lagging_strand():
    feats = [
        {"id": "g1", "start": 100, "end": 200, "strand": "+"},  # right replichore, + -> leading
        {"id": "g2", "start": 300, "end": 400, "strand": "-"},  # right replichore, - -> lagging
        {"id": "g3", "start": 700, "end": 800, "strand": "-"},  # left replichore, - -> leading
        {"id": "g4", "start": 900, "end": 950, "strand": "+"},  # left replichore, + -> lagging
    ]
    r = run(rr.annotate, {"length": 1000, "origin": 0, "terminus": 500, "features": feats})
    assert [f["leading_strand"] for f in r["features"]] == [True, False, True, False]
    assert [f["replichore"] for f in r["features"]] == ["right", "right", "left", "left"]
    assert r["leading_strand_fraction"] == 0.5
    wrapped = run(
        rr.annotate,
        {
            "length": 1000,
            "features": [
                {"id": "a", "start": 940, "end": 960, "strand": "+"},
                {"id": "b", "start": 490, "end": 510, "strand": "+"},
            ],
        },
        origin=800,
        terminus=300,
    )
    assert [f["replichore"] for f in wrapped["features"]] == [
        "right",
        "left",
    ]  # a: offset 150 < 500; b: offset 700
    assert [f["leading_strand"] for f in wrapped["features"]] == [True, False]
    est = run(
        rr.annotate,
        {
            "sequence": two_replichore(),
            "features": [{"id": "g", "start": 10, "end": 20, "strand": "+"}],
        },
    )
    # window 250: cumulative skew 0.5, 1.0, 0.5, 0.0 -> terminus (max) in [250, 500), origin (min) in [750, 1000)
    assert "heuristic" in est["origin_source"] and (est["origin"], est["terminus"]) == (875, 375)


def test_replication_compare_skew_profiles():
    seq = two_replichore()
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))
    same = run(rr.compare, {"a": seq, "b": seq}, window=50)
    assert same["cumulative_skew_correlation"] == pytest.approx(1.0)
    assert same["origin_relative_position_difference"] == 0.0
    flipped = run(rr.compare, [seq, comp], window=50)
    assert flipped["cumulative_skew_correlation"] == pytest.approx(-1.0)


def test_replication_map_timing_log2_ratio():
    r = run(rr.map, {"s_phase": [20, 10, 5, 10], "g1": [10, 10, 10, 10]}, window_size=1000)
    assert [w["log2_s_over_g1"] for w in r["window_table"]] == [1.0, 0.0, -1.0, 0.0]
    assert [w["timing"] for w in r["window_table"]] == [
        "early",
        "intermediate",
        "late",
        "intermediate",
    ]
    assert r["early_fraction"] == 0.25 and r["late_fraction"] == 0.25
    assert r["window_table"][2]["start"] == 2000
    peaks = run(rr.map, {"s_phase": [10, 10, 40, 10, 10], "g1": [10] * 5})
    assert [p["window"] for p in peaks["early_peaks"]] == [2]
    smooth = run(rr.map, {"s_phase": [10, 10, 40, 10, 10], "g1": [10] * 5}, smoothing_windows=3)
    assert smooth["early_peaks"] == []  # the running median removes the single-window spike
    with pytest.raises(InputFormatError):
        rr.map({"s_phase": [0, 1, 1], "g1": [1, 1, 1]})
    assert run(rr.map, {"s_phase": [0, 1, 1], "g1": [1, 1, 1]}, pseudocount=0.5)["windows"] == 3


def test_replication_simulate_timing_and_termination():
    one = run(rr.simulate, 1000, origins=[500], fork_speed=10, window=1)
    assert (
        one["s_phase_duration"] == 50.0 and one["termination_sites"] == []
    )  # farthest base (position 0) is 500 bp away
    two = run(rr.simulate, {"length": 1000}, origins=[250, 750], fork_speed=10, window=1)
    assert two["s_phase_duration"] == 25.0
    assert two["termination_sites"] == [{"position": 500.0, "time": 25.0}]
    late = run(
        rr.simulate,
        1000,
        origins=[{"position": 250, "firing_time": 0}, {"position": 750, "firing_time": 10}],
        fork_speed=10,
        window=1,
    )
    assert late["termination_sites"] == [{"position": 550.0, "time": 30.0}]
    assert late["s_phase_duration"] == pytest.approx(34.9)  # position 999: 10 + 249/10
    circular = run(rr.simulate, 1000, origins=[0], fork_speed=10, window=1, circular=True)
    assert circular["s_phase_duration"] == 50.0  # farthest point is the antipode, 500 bp away
    assert (
        run(rr.simulate, 1000, n_origins=3, seed=7)["origins"]
        == run(rr.simulate, 1000, n_origins=3, seed=7)["origins"]
    )
    assert run(rr.simulate, "ACGT" * 100, origins=[0], fork_speed=100)["length"] == 400
    curve = one["replication_curve"]
    assert curve[0]["fraction_replicated"] > 0 and curve[-1]["fraction_replicated"] == 1.0


def test_replication_build_recombination_map():
    haldane = run(
        rr.build_recombination_map,
        [
            {"id": "m1", "position": 0},
            {"id": "m2", "position": 1_000_000, "recombination_fraction": 0.1},
        ],
    )
    assert haldane["map"][1]["delta_cm"] == pytest.approx(
        100 * -0.5 * math.log(0.8), abs=1e-5
    )  # 11.157 cM
    assert haldane["map"][1]["rate_cm_per_mb"] == pytest.approx(11.157178, abs=1e-5)
    kosambi = run(
        rr.build_recombination_map,
        [{"position": 0}, {"position": 1_000_000, "recombination_fraction": 0.1}],
        map_function="kosambi",
    )
    assert kosambi["map"][1]["delta_cm"] == pytest.approx(
        100 * 0.25 * math.log(1.2 / 0.8), abs=1e-5
    )  # 10.137 cM
    counts = run(
        rr.build_recombination_map,
        [{"position": 0}, {"position": 2_000_000, "recombinants": 5, "meioses": 50}],
    )
    assert counts["map"][1]["recombination_fraction"] == 0.1 and counts["map"][1][
        "rate_cm_per_mb"
    ] == pytest.approx(11.157178 / 2, abs=1e-5)
    given = run(
        rr.build_recombination_map,
        {
            "markers": [
                {"position": 0, "cm": 0},
                {"position": 1_000_000, "cm": 1.5},
                {"position": 3_000_000, "cm": 4.5},
            ]
        },
    )
    assert [m.get("rate_cm_per_mb") for m in given["map"][1:]] == [1.5, 1.5]
    assert given["total_cm"] == 4.5 and given["mean_rate_cm_per_mb"] == pytest.approx(1.5)


X25 = "ACGTAGCTTGCAAGTCCTGAATCGA"


@pytest.mark.parametrize(
    ("reference", "edited", "extra", "pathway", "microhomology"),
    [
        (
            "GGATGCATGCCGG",
            "GGATGCCGG",
            {},
            "MMEJ",
            4,
        ),  # TGCA deleted; prefix 6 + suffix 7 - 9 = 4 bp slide
        ("GGATCCAGTCAATG", "GGATTCAATG", {}, "NHEJ", 0),  # CCAG deleted, no junction homology
        ("GGATCCAG", "GGATACCAG", {}, "NHEJ", None),  # 1 bp insertion
        ("GG" + X25 + "TTTT" + X25 + "CC", "GG" + X25 + "CC", {}, "SSA", 25),  # 25 bp direct repeat
        ("GGATCCAG", "GGATAAAG", {"donor": "GGATAAAG"}, "HDR", None),
        ("GGATCCAG", "GGATCCAG", {}, "no_edit", None),
        ("GGATCCAG", "GGATGCAG", {}, "undetermined", None),  # substitution only
        ("GGATCCAGTT", "GGATAAAAAGTT", {}, "undetermined", None),  # complex indel
    ],
)
def test_replication_classify_repair_pathway(reference, edited, extra, pathway, microhomology):
    assert len(X25) == 25
    r = run(rr.classify_repair_pathway, {"reference": reference, "edited": edited, **extra})
    assert r["pathway"] == pathway
    assert r.get("microhomology_length") == microhomology or microhomology is None
    assert "rule" in r


def test_replication_classify_repair_pathway_options():
    pair = {"reference": "GGATGCATGCCGG", "edited": "GGATGCCGG"}
    assert run(rr.classify_repair_pathway, pair, min_microhomology=5)["pathway"] == "NHEJ"
    long_ins = {"reference": "GGATCC", "edited": "GGAT" + "A" * 30 + "CC"}
    assert run(rr.classify_repair_pathway, long_ins)["pathway"] == "undetermined"
    assert run(rr.classify_repair_pathway, long_ins, max_nhej_insertion=40)["pathway"] == "NHEJ"


def test_replication_plot_marks_origin_and_terminus(tmp_path):
    out = tmp_path / "skew.svg"
    r = run(rr.plot, two_replichore(), window=100, output_path=str(out))
    assert r["svg"].startswith("<svg") and "<polyline" in r["svg"]
    assert (r["origin_estimate"], r["terminus_estimate"]) == (950, 450)
    assert out.read_text() == r["svg"]
    with pytest.raises(ConfigurationError):
        rr.plot(two_replichore(), options={"window": 100, "output_path": str(out)})


@pytest.mark.parametrize(
    ("fn", "value", "options", "error"),
    [
        (rr.calculate, "", None, InputFormatError),
        (rr.calculate, "ACGT!", None, InputFormatError),
        (rr.calculate, "ACGT", {"window": 100}, ConfigurationError),
        (rr.calculate, "ACGTACGT", {"skew_type": "XY"}, ConfigurationError),
        (rr.calculate, 42, None, InputFormatError),
        (rr.annotate, {"length": 100, "features": []}, None, InputFormatError),
        (
            rr.annotate,
            {
                "length": 100,
                "features": [{"id": "x", "start": 5, "end": 500, "strand": "+"}],
                "origin": 0,
                "terminus": 50,
            },
            None,
            InputFormatError,
        ),
        (
            rr.annotate,
            {"features": [{"id": "x", "start": 5, "end": 9, "strand": "+"}]},
            None,
            InputFormatError,
        ),
        (rr.compare, ["ACGT"], None, InputFormatError),
        (rr.map, {"s_phase": [1, 2], "g1": [1]}, None, InputFormatError),
        (
            rr.map,
            {"s_phase": [1, 2, 3], "g1": [1, 1, 1]},
            {"smoothing_windows": 2},
            ConfigurationError,
        ),
        (rr.simulate, 1000, {"origins": [5000]}, ConfigurationError),
        (rr.simulate, "x", None, InputFormatError),
        (rr.build_recombination_map, [{"position": 5}], None, InputFormatError),
        (
            rr.build_recombination_map,
            [{"position": 5}, {"position": 3, "cm": 1}],
            None,
            InputFormatError,
        ),
        (
            rr.build_recombination_map,
            [{"position": 0}, {"position": 9, "recombination_fraction": 0.7}],
            None,
            InputFormatError,
        ),
        (rr.classify_repair_pathway, {"reference": "ACGT"}, None, InputFormatError),
        (rr.detect, "ACGT", {"min_copies": {1: 1}}, ConfigurationError),
        (rr.predict, "ACGT", {"window": 5}, ConfigurationError),
        (rr.plot, "ACGTACGT", {"window": 4, "output_path": "x.png"}, ConfigurationError),
    ],
)
def test_replication_repair_invalid_input_raises(fn, value, options, error):
    with pytest.raises(error):
        fn(value, options=options)


# =========================================================================== molecular_lab


def nn_tm(seq: str, sodium: float = 0.05, strand: float = 250e-9) -> float:
    """Independent SantaLucia (1998) unified nearest-neighbour Tm for a non-self-complementary duplex."""

    table = {
        "AA": (-7.9, -22.2),
        "TT": (-7.9, -22.2),
        "AT": (-7.2, -20.4),
        "TA": (-7.2, -21.3),
        "CA": (-8.5, -22.7),
        "TG": (-8.5, -22.7),
        "GT": (-8.4, -22.4),
        "AC": (-8.4, -22.4),
        "CT": (-7.8, -21.0),
        "AG": (-7.8, -21.0),
        "GA": (-8.2, -22.2),
        "TC": (-8.2, -22.2),
        "CG": (-10.6, -27.2),
        "GC": (-9.8, -24.4),
        "GG": (-8.0, -19.9),
        "CC": (-8.0, -19.9),
    }
    dh = ds = 0.0
    for a, b in itertools.pairwise(seq):
        h, s = table[a + b]
        dh += h
        ds += s
    for end in (seq[0], seq[-1]):
        dh, ds = (dh + 0.1, ds - 2.8) if end in "GC" else (dh + 2.3, ds + 4.1)
    ds += 0.368 * (len(seq) - 1) * math.log(sodium)
    return dh * 1000 / (ds + 1.987 * math.log(strand / 4)) - 273.15


PRIMER = "CGTTGACCGTAGCTAGGCTA"


def test_lab_calculate_dilution_solves_each_unknown():
    r = run(ml.calculate, {"c1": 10, "c2": 1, "v2": 100}, calculation="dilution")
    assert (r["solved"], r["v1"], r["diluent_volume"], r["dilution_factor"]) == (
        "v1",
        10.0,
        90.0,
        10.0,
    )
    assert run(ml.calculate, {"c1": 10, "v1": 5, "v2": 50}, calculation="dilution")["c2"] == 1.0
    assert run(ml.calculate, {"v1": 5, "c2": 1, "v2": 50}, calculation="dilution")["c1"] == 10.0
    assert run(ml.calculate, {"c1": 10, "v1": 5, "c2": 2}, calculation="dilution")["v2"] == 25.0


def test_lab_calculate_molarity_and_copy_number():
    r = run(
        ml.calculate, {"length": 1000, "mass_concentration_ng_per_ul": 100}, calculation="molarity"
    )
    # 100 ng/uL = 0.1 g/L; MW of 1000 bp = 660000 g/mol -> 1.51515e-7 M
    assert r["molar_concentration_nm"] == pytest.approx(151.515, abs=1e-3) and r[
        "pmol_per_ul"
    ] == pytest.approx(0.151515, abs=1e-5)
    back = run(
        ml.calculate,
        {"length": 1000, "molar_concentration_nm": 151.515151515},
        calculation="molarity",
    )
    assert back["mass_concentration_ng_per_ul"] == pytest.approx(100.0, abs=1e-4)
    ss = run(
        ml.calculate,
        {"length": 100, "strand": "ss", "mass_concentration_ng_per_ul": 33},
        calculation="molarity",
    )
    assert ss["molecular_weight"] == 33000.0 and ss["molar_concentration_nm"] == pytest.approx(
        1000.0
    )
    exact = run(
        ml.calculate,
        {"sequence": "ACGT", "strand": "ss", "mass_concentration_ng_per_ul": 1173.84},
        calculation="molarity",
    )
    assert exact["molecular_weight"] == pytest.approx(1173.84) and exact[
        "molar_concentration_nm"
    ] == pytest.approx(1e6)
    c = run(ml.calculate, {"length": 1000, "mass_ng": 1}, calculation="copy_number")
    assert c["copies"] == pytest.approx(6.02214076e23 * 1e-9 / 660000, rel=1e-8)  # 9.12e8
    assert c["copies"] == pytest.approx(9.1244e8, rel=1e-4)


def test_lab_calculate_ligation_and_a260_and_pcr_yield():
    lig = run(
        ml.calculate,
        {
            "vector_bp": 3000,
            "insert_bp": 1000,
            "vector_ng": 50,
            "vector_conc_ng_per_ul": 25,
            "insert_conc_ng_per_ul": 10,
        },
        calculation="ligation",
    )
    assert (
        lig["insert_ng"] == 50.0
        and lig["vector_volume_ul"] == 2.0
        and lig["insert_volume_ul"] == 5.0
    )
    assert lig["insert_pmol"] / lig["vector_pmol"] == pytest.approx(3.0, rel=1e-4)
    assert (
        run(ml.calculate, {"a260": 0.5, "dilution_factor": 10}, calculation="a260")[
            "mass_concentration_ng_per_ul"
        ]
        == 250.0
    )
    assert (
        run(ml.calculate, {"a260": 0.5, "nucleic_acid": "RNA"}, calculation="a260")[
            "mass_concentration_ng_per_ul"
        ]
        == 20.0
    )
    oligo = run(ml.calculate, {"a260": 0.403, "sequence": "ACGT"}, calculation="a260")
    assert oligo["extinction_coefficient"] == 40300.0  # (AC + CG + GT) - (C + G) = 59200 - 18900
    assert oligo["molar_concentration_um"] == pytest.approx(10.0) and oligo[
        "mass_concentration_ng_per_ul"
    ] == pytest.approx(11.7384)
    assert (
        run(ml.calculate, {"template_copies": 1000, "cycles": 10}, calculation="pcr_yield")[
            "product_copies"
        ]
        == 1_024_000
    )
    assert run(
        ml.calculate,
        {"template_copies": 1000, "cycles": 10, "efficiency": 0.9},
        calculation="pcr_yield",
    )["product_copies"] == pytest.approx(613106.6, rel=1e-5)


def test_lab_calculate_tm_methods_and_gel_size():
    assert (
        run(ml.calculate, {"sequence": "ATGC"}, calculation="tm", method="wallace")["tm_celsius"]
        == 12.0
    )
    assert run(ml.calculate, {"sequence": "ATGCATGCATGCATGC"}, calculation="tm", method="basic")[
        "tm_celsius"
    ] == pytest.approx(43.375)
    nn = run(ml.calculate, {"sequence": PRIMER}, calculation="tm")
    assert nn["tm_celsius"] == pytest.approx(nn_tm(PRIMER), abs=0.02)
    salty = run(ml.calculate, {"sequence": PRIMER}, calculation="tm", sodium_mm=150)
    assert (
        salty["tm_celsius"] == pytest.approx(nn_tm(PRIMER, sodium=0.15), abs=0.02)
        and salty["tm_celsius"] > nn["tm_celsius"]
    )
    gel = run(
        ml.calculate,
        {
            "ladder": [{"bp": 1000, "distance": 10}, {"bp": 100, "distance": 20}],
            "distances": [15, 10],
        },
        calculation="gel_size",
    )
    assert gel["slope_log10bp_per_distance"] == -0.1 and gel["r_squared"] == 1.0
    assert gel["estimated_bp"] == [pytest.approx(316.228, abs=1e-3), pytest.approx(1000.0)]


def test_lab_analyze_oligo_properties():
    r = run(ml.analyze, "ACGT")["sequences"]["sequence"]
    assert r["composition"] == {"A": 1, "C": 1, "G": 1, "T": 1} and r["gc_fraction"] == 0.5
    assert r["reverse_complement"] == "ACGT" and r["is_self_complementary"] is True
    assert r["molecular_weight_da"] == pytest.approx(1173.84)
    assert r["extinction_coefficient_260nm"] == 40300.0 and r["nmol_per_od260"] == pytest.approx(
        1e6 / 40300, rel=1e-5
    )
    p = run(ml.analyze, {"p": PRIMER}, tm_method="nearest_neighbor")["sequences"]["p"]
    assert (
        p["tm_celsius"] == pytest.approx(nn_tm(PRIMER), abs=0.02) and p["longest_homopolymer"] == 2
    )
    short = run(ml.analyze, "ACGT", tm_method="basic")["sequences"]["sequence"]
    assert short["tm_celsius"] is None and "14 nt" in short["tm_note"]


def test_lab_validate_primer_rules():
    good = run(ml.validate, PRIMER, min_tm=50)
    checks = {c["check"]: c["passed"] for c in good["primers"]["sequence"]["checks"]}
    assert checks == {
        "length": True,
        "gc_content": True,
        "three_prime_gc_clamp": False,
        "homopolymer_runs": True,
        "dinucleotide_repeats": True,
        "melting_temperature": True,
    }
    bad = run(ml.validate, "AAAAAAAAAAAAAAAAAAAA")["primers"]["sequence"]
    failed = {c["check"] for c in bad["checks"] if not c["passed"]}
    assert (
        not bad["passed"]
        and {"gc_content", "three_prime_gc_clamp", "homopolymer_runs", "melting_temperature"}
        <= failed
    )
    pair = run(ml.validate, {"forward": "GCGCGCGCGCGCGCGCGCGC", "reverse": "ATATATATATATATATATGC"})
    assert pair["kind"] == "primer_pair" and pair["pair_checks"][0]["passed"] is False
    short = run(ml.validate, "ACGTAC")["primers"]["sequence"]["checks"][0]
    assert short["passed"] is False and "6 nt" in short["detail"]


def test_lab_optimize_primer_length_to_target_tm():
    template = "CGTTGACCGTAGCTAGGCTAATCGGATTCGCAGGT"
    r = run(ml.optimize, template, target_tm=60, min_length=18, max_length=30)
    expected = {n: nn_tm(template[:n]) for n in range(18, 31)}
    for row in r["candidates"]:
        assert row["tm_celsius"] == pytest.approx(expected[row["length"]], abs=0.02)
    assert r["best"]["length"] == min(expected, key=lambda n: (abs(expected[n] - 60), n))
    rev = run(ml.optimize, template, direction="reverse", min_length=18, max_length=20)
    assert rev["best"]["primer"] == revcomp(template[-rev["best"]["length"] :])
    clamp = run(ml.optimize, template, require_gc_clamp=True, min_length=18, max_length=30)
    assert clamp["best"]["gc_clamp"] is True
    with pytest.raises(InputFormatError):
        ml.optimize("ACGT", options=None)


def test_lab_plan_pcr_master_mix():
    r = run(ml.plan, {}, protocol="pcr", n_reactions=10)
    rows = {c["component"]: c for c in r["components"]}
    assert (
        rows["reaction buffer"]["volume_per_reaction_ul"] == 5.0
        and rows["reaction buffer"]["master_mix_volume_ul"] == 55.0
    )
    assert rows["dNTP mix"]["master_mix_volume_ul"] == 5.5
    assert (
        rows["forward primer"]["volume_per_reaction_ul"] == 1.25
        and rows["forward primer"]["master_mix_volume_ul"] == 13.75
    )
    assert rows["DNA polymerase"]["master_mix_volume_ul"] == 1.375
    assert rows["water"]["volume_per_reaction_ul"] == 15.875 and rows["water"][
        "master_mix_volume_ul"
    ] == pytest.approx(174.625)
    assert (
        r["master_mix_total_ul"] == pytest.approx(24 * 11)
        and r["master_mix_aliquot_per_tube_ul"] == 24.0
    )
    with pytest.raises(ConfigurationError):
        ml.plan(
            {}, options={"protocol": "pcr", "template_volume_ul": 30}
        )  # template alone exceeds the 25 uL reaction


def test_lab_plan_serial_dilution_and_ligation():
    r = run(ml.plan, {"start_concentration": 100, "steps": 3}, protocol="serial_dilution")
    assert [t["concentration"] for t in r["tubes"]] == [100, 10, 1, 0.1]
    assert (
        r["tubes"][1]["transfer_from_previous_ul"] == 10.0
        and r["tubes"][1]["diluent_ul"] == 90.0
        and r["total_diluent_ul"] == 270.0
    )
    lig = run(
        ml.plan,
        {
            "vector_bp": 3000,
            "insert_bp": 1000,
            "vector_ng": 50,
            "vector_conc_ng_per_ul": 25,
            "insert_conc_ng_per_ul": 10,
        },
        protocol="ligation",
    )
    assert (
        lig["vector_volume_ul"],
        lig["insert_volume_ul"],
        lig["ligase_buffer_10x_ul"],
        lig["ligase_ul"],
        lig["water_ul"],
    ) == (2.0, 5.0, 2.0, 1.0, 10.0)


def test_lab_parse_quantities_and_wells():
    r = run(ml.parse, ["50 ng/uL", "1.5 mM", "200 µL", "10 nmol", "5X", "2 mg/mL"])["quantities"]
    assert (r[0]["dimension"], r[0]["si_value"], r[0]["si_unit"]) == (
        "mass_concentration",
        0.05,
        "g/L",
    )
    assert (r[1]["dimension"], r[1]["si_value"]) == ("molar_concentration", 0.0015)
    assert (r[2]["dimension"], r[2]["si_value"]) == ("volume", 2e-4)
    assert (r[3]["dimension"], r[3]["si_value"]) == ("amount", 1e-8)
    assert r[4]["dimension"] == "fold_concentration" and r[5]["si_value"] == 2.0
    assert run(ml.parse, "A1:B3", kind="wells")["wells"] == ["A1", "A2", "A3", "B1", "B2", "B3"]
    assert run(ml.parse, "A1:B2, H12", kind="wells", order="column")["wells"] == [
        "A1",
        "B1",
        "A2",
        "B2",
        "H12",
    ]
    assert run(ml.parse, "P24", kind="wells", plate=384)["wells"] == ["P24"]


def test_lab_export_formats(tmp_path):
    records = [
        {"id": "p1", "sequence": "ACGT", "conc": 5},
        {"id": "p2", "sequence": "GG", "conc": None},
    ]
    assert run(ml.export, records)["text"] == "id,sequence,conc\np1,ACGT,5\np2,GG,\n"
    assert (
        run(ml.export, records, format="tsv", columns=["sequence", "id"])["text"]
        == "sequence\tid\nACGT\tp1\nGG\tp2\n"
    )
    assert run(ml.export, records, format="fasta")["text"] == ">p1\nACGT\n>p2\nGG\n"
    assert (
        run(ml.export, {"a": "ACGTACGT"}, format="fasta", line_width=10)["text"] == ">a\nACGTACGT\n"
    )
    import json

    assert json.loads(run(ml.export, records, format="json", columns=["id"])["text"]) == [
        {"id": "p1"},
        {"id": "p2"},
    ]
    target = tmp_path / "out.csv"
    out = run(ml.export, records, output_path=str(target))
    assert target.read_text() == out["text"]
    with pytest.raises(ConfigurationError):
        ml.export(records, options={"output_path": str(target)})
    assert (
        run(ml.export, [{"id": "q,1", "note": 'say "hi"'}])["text"]
        == 'id,note\n"q,1","say ""hi"""\n'
    )


def test_lab_plot_virtual_gel(tmp_path):
    r = run(ml.plot, {"ladder": [1000, 500, 100], "PCR": [650]})
    svg = r["svg"]
    import re

    def y_of(bp: int) -> float:
        m = re.search(rf'<rect x="[\d.]+" y="([\d.]+)"[^>]*><title>{bp} bp</title>', svg)
        assert m, bp
        return float(m.group(1))

    assert y_of(1000) < y_of(650) < y_of(500) < y_of(100)  # larger fragments migrate less
    # semilog scale: equal spacing for each 10-fold step (min_bp=100, max_bp=10000)
    decades = run(ml.plot, {"x": [10000, 1000, 100]})["svg"]
    ys = [
        float(
            re.search(rf'<rect x="[\d.]+" y="([\d.]+)"[^>]*><title>{bp} bp</title>', decades).group(
                1
            )
        )
        for bp in (10000, 1000, 100)
    ]
    assert ys[1] - ys[0] == pytest.approx(ys[2] - ys[1], abs=0.2) and ys[0] < ys[1] < ys[2]
    assert r["lanes"] == 2 and svg.startswith("<svg")
    out = tmp_path / "gel.svg"
    assert run(ml.plot, [[100, 200]], output_path=str(out))["output_path"] == str(out)


def _gauss(n: int, centres: dict[int, float], sigma: float = 2.0) -> list[float]:
    return [
        sum(h * math.exp(-((i - c) ** 2) / (2 * sigma**2)) for c, h in centres.items())
        for i in range(n)
    ]


def test_lab_call_bases_from_traces_with_mixed_peak():
    n = 100
    traces = {
        "A": _gauss(n, {10: 100, 90: 100, 50: 60}),
        "C": _gauss(n, {30: 100}),
        "G": _gauss(n, {50: 100}),
        "T": _gauss(n, {70: 100}),
    }
    given = run(
        ml.call_bases_from_sanger, {"traces": traces, "peak_locations": [10, 30, 50, 70, 90]}
    )
    assert given["sequence"] == "ACRTA"  # position 50: G 100 with A 60 (ratio 0.6 >= 0.5) -> R
    assert (
        given["mixed_positions"][0]["bases"] == "GA"
        and given["secondary_to_primary_peak_ratio"][2] == 0.6
    )
    detected = run(ml.call_bases_from_sanger, {"traces": traces})
    assert detected["sequence"] == "ACRTA" and detected["peaks_detected_automatically"] is True
    pure = run(ml.call_bases_from_sanger, {"traces": traces}, mixed_peak_ratio=0.7)
    assert pure["sequence"] == "ACGTA"


def write_abif(
    path: Path, order: str, traces: dict[str, list[float]], peaks: list[int], called: str
) -> None:
    """Minimal ABIF writer: header, data blobs, then a directory of (DATA9-12, FWO_1, PLOC2, PBAS2)."""

    entries = []
    blob = b""
    base = 34

    def add(name: bytes, number: int, etype: int, esize: int, count: int, data: bytes) -> None:
        nonlocal blob
        if len(data) <= 4:
            entries.append((name, number, etype, esize, count, len(data), data.ljust(4, b"\0")))
        else:
            entries.append(
                (name, number, etype, esize, count, len(data), struct.pack(">I", base + len(blob)))
            )
            blob += data

    for k, b in enumerate(order):
        series = [round(x) for x in traces[b]]
        add(b"DATA", 9 + k, 4, 2, len(series), struct.pack(f">{len(series)}h", *series))
    add(b"FWO_", 1, 2, 1, 4, order.encode())
    add(b"PLOC", 2, 4, 2, len(peaks), struct.pack(f">{len(peaks)}h", *peaks))
    add(b"PBAS", 2, 2, 1, len(called), called.encode())
    directory = b"".join(
        struct.pack(">4sI2H3I", n, num, t, s, c, d, 0)[:20] + off + struct.pack(">I", 0)
        for n, num, t, s, c, d, off in entries
    )
    dir_offset = base + len(blob)
    header = (
        b"ABIF"
        + struct.pack(">H", 101)
        + struct.pack(">4sI2H4I", b"tdir", 1, 1023, 28, len(entries), len(directory), dir_offset, 0)
    )
    path.write_bytes(header + blob + directory)


def test_lab_call_bases_from_ab1_file(tmp_path):
    n = 100
    by_base = {
        "A": _gauss(n, {10: 1000, 90: 1000}),
        "C": _gauss(n, {30: 1000}),
        "G": _gauss(n, {50: 1000}),
        "T": _gauss(n, {70: 1000}),
    }
    path = tmp_path / "x.ab1"
    write_abif(path, "GATC", by_base, [10, 30, 50, 70, 90], "ACGTA")
    r = run(ml.call_bases_from_sanger, str(path))
    assert r["sequence"] == "ACGTA" and r["channel_order"] == "GATC" and r["source"] == "ABIF file"
    assert r["agreement_with_instrument_basecalls"] == 1.0
    bad = tmp_path / "bad.ab1"
    bad.write_bytes(b"NOPE" + b"\0" * 40)
    with pytest.raises(InputFormatError):
        ml.call_bases_from_sanger(str(bad))


SANGER = "ATGGCTAGCTAGGATCCGATTACAGGCTTAACGGATCGCTAGGCATCGATTCGGCAAGTCTTGACCATGCAGG"


def test_lab_assemble_sanger_reads_both_strands():
    reads = {"f1": SANGER[:30], "f2": SANGER[20:50], "r3": revcomp(SANGER[40:])}
    r = run(ml.assemble_sanger_reads, reads, min_overlap=10)
    assert len(r["contigs"]) == 1 and r["singlets"] == []
    contig = r["contigs"][0]
    assert contig["sequence"] in (SANGER, revcomp(SANGER)) and contig["length"] == len(SANGER)
    by_id = {x["id"]: x for x in contig["reads"]}
    assert set(by_id) == {"f1", "f2", "r3"}
    assert (
        by_id["f1"]["strand"] != by_id["r3"]["strand"]
    )  # r3 is the reverse complement of the other two
    assert by_id["f1"]["strand"] == by_id["f2"]["strand"]
    no_rc = run(ml.assemble_sanger_reads, reads, min_overlap=10, both_strands=False)
    assert len(no_rc["contigs"]) == 1 and [s["id"] for s in no_rc["singlets"]] == ["r3"]


def test_lab_assemble_sanger_reads_identity_threshold_and_containment():
    a, b = SANGER[:40], SANGER[20:60]
    mutated = b[:5] + ("C" if b[5] != "C" else "G") + b[6:]  # 1 mismatch in the 20 bp overlap
    ok = run(ml.assemble_sanger_reads, [a, mutated], min_overlap=10, min_identity=0.95)
    assert len(ok["contigs"]) == 1 and ok["contigs"][0]["length"] == 60
    assert ok["contigs"][0]["sequence"][:40] == a  # the earlier read's base is kept in the overlap
    strict = run(ml.assemble_sanger_reads, [a, mutated], min_overlap=10, min_identity=0.99)
    assert strict["contigs"] == [] and len(strict["singlets"]) == 2
    contained = run(ml.assemble_sanger_reads, [SANGER[:50], SANGER[10:30]], min_overlap=10)
    assert contained["contigs"][0]["sequence"] == SANGER[:50]


@pytest.mark.parametrize(
    ("fn", "value", "options", "error"),
    [
        (ml.calculate, {}, None, ConfigurationError),  # calculation is required
        (ml.calculate, {"c1": 1, "v1": 1}, {"calculation": "dilution"}, ConfigurationError),
        (
            ml.calculate,
            {"c1": 1, "v1": 1, "c2": 5},
            {"calculation": "dilution"},
            ConfigurationError,
        ),  # c2 > c1
        (ml.calculate, {"length": 10}, {"calculation": "molarity"}, ConfigurationError),
        (
            ml.calculate,
            {"sequence": "ACGTAC"},
            {"calculation": "tm", "method": "basic"},
            ConfigurationError,
        ),
        (ml.calculate, {"sequence": "ACGU!"}, {"calculation": "tm"}, InputFormatError),
        (
            ml.calculate,
            {"template_copies": 10, "cycles": 5, "efficiency": 1.5},
            {"calculation": "pcr_yield"},
            ConfigurationError,
        ),
        (
            ml.calculate,
            {"ladder": [{"bp": 100, "distance": 5}]},
            {"calculation": "gel_size"},
            ConfigurationError,
        ),
        (ml.analyze, "", None, InputFormatError),
        (ml.analyze, "ACGTN", None, InputFormatError),
        (ml.analyze, 42, None, InputFormatError),
        (ml.validate, "ACGT!", None, InputFormatError),
        (ml.parse, "fifty", None, InputFormatError),
        (ml.parse, "5 parsecs", None, InputFormatError),
        (ml.parse, "I1", {"kind": "wells"}, InputFormatError),
        (ml.parse, "A1", {"kind": "wells", "plate": 12}, ConfigurationError),
        (ml.plan, {}, {"protocol": "gel"}, ConfigurationError),
        (
            ml.plan,
            {"start_concentration": 1, "steps": 2, "dilution_factor": 1},
            {"protocol": "serial_dilution"},
            ConfigurationError,
        ),
        (ml.export, 42, None, InputFormatError),
        (ml.export, [{"id": "a"}], {"format": "fasta"}, InputFormatError),
        (ml.export, [{"id": "a"}], {"format": "xml"}, ConfigurationError),
        (ml.plot, {"x": [-5]}, None, InputFormatError),
        (ml.optimize, "ACGTACGTACGTACGTACGT", {"direction": "up"}, ConfigurationError),
        (ml.call_bases_from_sanger, {"traces": {"A": [1], "C": [1]}}, None, InputFormatError),
        (
            ml.call_bases_from_sanger,
            {"traces": {b: [0.0, 0.0, 0.0] for b in "ACGT"}},
            None,
            InputFormatError,
        ),
        (ml.assemble_sanger_reads, ["ACGT"], {"min_identity": 0.1}, ConfigurationError),
        (ml.assemble_sanger_reads, [], None, InputFormatError),
    ],
)
def test_molecular_lab_invalid_input_raises(fn, value, options, error):
    with pytest.raises(error):
        fn(value, options=options)


# =========================================================================== envelope and registry


@pytest.mark.parametrize(
    ("fn", "operation", "module"),
    [
        (sq.calculate_yield, "calculate_yield", "sequencing_qc"),
        (rr.calculate, "calculate", "replication_repair"),
        (ml.parse, "parse", "molecular_lab"),
    ],
)
def test_envelope_shape(fn, operation, module):
    value = {"sequencing_qc": FQ_A, "replication_repair": "GGGCCCCG", "molecular_lab": "10 mM"}[
        module
    ]
    envelope = fn(value)
    assert (
        envelope["status"] == "ok"
        and envelope["module"] == module
        and envelope["operation"] == operation
    )
    assert envelope["provenance"]["implementation"] == "native" and envelope["result"]
    assert "input_summary" in envelope and isinstance(envelope["warnings"], list)
