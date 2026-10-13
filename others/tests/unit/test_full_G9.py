"""Tests for the cutadapt and umi-tools wrappers in dnakit.sequencing_qc (group G9).

Every wrapper is checked against a hand calculation or a direct call to the upstream function.
Missing optional packages are simulated with sys.modules entries set to None.
"""

from __future__ import annotations

import inspect
import itertools
import json
import os
import random
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.sequencing_qc import cutadapt_functions as cut
from dnakit.sequencing_qc import umi_tools_functions as umi

SRC = Path(__file__).resolve().parents[3] / "src"

CUTADAPT_NAMES = [
    "adapter_kmer_positions",
    "adapter_kmer_table",
    "edit_neighbours",
    "expand_adapter_braces",
    "hamming_neighbours",
    "nextseq_trim_bounds",
    "parse_adapter_spec",
    "poly_a_tail_start",
    "read_adapter_fasta",
    "sequence_kmer_chunks",
]
UMI_NAMES = [
    "average_umi_distance",
    "barcode_indel_matches",
    "call_cell_whitelist",
    "cigar_splice_offset",
    "count_umis_per_gene",
    "detect_bam_tags",
    "extract_barcode_positions",
    "find_umi_neighbours",
    "mask_low_quality_umi",
    "parse_umi_read_name",
    "read_barcode_whitelist",
    "tag_read_name_with_umi",
    "umi_connected_component",
]


def _hide_package(monkeypatch: pytest.MonkeyPatch, package: str) -> None:
    """Make a package and every already-imported submodule unimportable for one test."""

    for name in list(sys.modules):
        if name == package or name.startswith(package + "."):
            monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, package, None)


def test_public_names_match_all_and_take_native_verbose() -> None:
    for module, names in ((cut, CUTADAPT_NAMES), (umi, UMI_NAMES)):
        assert sorted(module.__all__) == names
        for name in names:
            parameters = inspect.signature(getattr(module, name)).parameters
            assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY
            assert parameters["verbose"].default is True


# cutadapt -----------------------------------------------------------------------------------


def test_hamming_neighbours_exact_sphere_and_environment_match_hand_counts() -> None:
    pytest.importorskip("cutadapt.align")
    from cutadapt import align

    exact = cut.hamming_neighbours("AC", 1, exact=True, verbose=False)
    # Two positions times three other bases: AA AG AT CC GC TC.
    assert [row["sequence"] for row in exact["neighbours"]] == ["AA", "AG", "AT", "CC", "GC", "TC"]
    assert exact["count"] == 6
    assert sorted(align.hamming_sphere("AC", 1)) == [row["sequence"] for row in exact["neighbours"]]

    environment = cut.hamming_neighbours("AC", 1, verbose=False)
    assert environment["count"] == 7  # the query itself plus the six exact neighbours
    assert environment["neighbours"][0] == {"sequence": "AC", "distance": 0, "matches": 2}


def test_hamming_neighbours_rejects_bad_input_and_respects_limit() -> None:
    pytest.importorskip("cutadapt.align")
    with pytest.raises(ConfigurationError):
        cut.hamming_neighbours("ACGN", 1, verbose=False)
    with pytest.raises(ConfigurationError):
        cut.hamming_neighbours("ACGT", True, verbose=False)
    with pytest.raises(ConfigurationError, match="max_neighbours"):
        cut.hamming_neighbours("ACGTACGT", 2, exact=True, max_neighbours=10, verbose=False)


def test_edit_neighbours_matches_brute_force_with_upstream_edit_distance() -> None:
    pytest.importorskip("cutadapt.align")
    from cutadapt.align import edit_distance

    result = cut.edit_neighbours("AC", 1, verbose=False)
    # Hand count for distance <= 1: 1 (itself) + 6 substitutions + 10 insertions + 2 deletions.
    assert result["count"] == 19
    brute = {
        "".join(letters)
        for length in (1, 2, 3)
        for letters in itertools.product("ACGT", repeat=length)
        if edit_distance("".join(letters), "AC") <= 1
    }
    assert {row["sequence"] for row in result["neighbours"]} == brute
    assert all(
        row["distance"] == edit_distance(row["sequence"], "AC") for row in result["neighbours"]
    )


def test_adapter_kmer_positions_matches_direct_upstream_call() -> None:
    heuristic = pytest.importorskip("cutadapt.kmer_heuristic")

    result = cut.adapter_kmer_positions("ACGTACGT", min_overlap=3, verbose=False)
    raw = heuristic.create_positions_and_kmers("ACGTACGT", 3, 0.1, True, False, True)
    expected = sorted((int(start), None, sorted(kmers)) for start, _stop, kmers in raw)
    observed = [(entry["start"], entry["stop"], entry["kmers"]) for entry in result["entries"]]
    assert observed == expected
    assert [entry["start"] for entry in result["entries"]] == [-8, -4, -3, 0]
    assert result["kmer_count"] == 4


def test_adapter_kmer_table_parses_the_upstream_table() -> None:
    heuristic = pytest.importorskip("cutadapt.kmer_heuristic")

    result = cut.adapter_kmer_table("ACGTACGT", verbose=False)
    raw = heuristic.create_positions_and_kmers("ACGTACGT", 3, 0.1, True, False, True)
    table = heuristic.kmer_probability_analysis(raw, 150)
    assert result["table_text"] == table
    assert result["row_count"] == 4
    assert result["rows"][0] == {
        "kmer": "ACG",
        "start": -3,
        "stop": None,
        "considered_sites": 1,
        "hit_chance_percent": 1.56,
    }
    assert result["profile_hit_chance_percent"] == pytest.approx(2.54)


def test_sequence_kmer_chunks_hand_calculation() -> None:
    pytest.importorskip("cutadapt.kmer_heuristic")
    # 10 bases in 3 chunks: sizes 4, 3, 3 -> ACGT, ACG, TAC.
    result = cut.sequence_kmer_chunks("ACGTACGTAC", 3, verbose=False)
    assert result["kmers"] == ["ACG", "ACGT", "TAC"]
    assert result["count"] == 3


def test_poly_a_tail_start_forward_and_reverse_complement() -> None:
    qualtrim = pytest.importorskip("cutadapt.qualtrim")

    tail = cut.poly_a_tail_start("ACGTAAAAAAAAA", verbose=False)
    assert (tail["index"], tail["length"]) == (4, 9)
    assert tail["index"] == qualtrim.poly_a_trim_index("ACGTAAAAAAAAA")

    short = cut.poly_a_tail_start("ACGTAA", verbose=False)  # tails under three bases are ignored
    assert (short["index"], short["length"]) == (6, 0)

    head = cut.poly_a_tail_start("TTTTTACGT", reverse_complement=True, verbose=False)
    assert (head["index"], head["length"]) == (5, 5)
    assert head["index"] == qualtrim.poly_a_trim_index("TTTTTACGT", revcomp=True)


def test_nextseq_trim_bounds_trims_the_g_run_and_matches_upstream() -> None:
    qualtrim = pytest.importorskip("cutadapt.qualtrim")
    dnaio = pytest.importorskip("dnaio")

    # The six trailing G bases score as cutoff - 1 (+1 each); the T at index 3 scores -20.
    result = cut.nextseq_trim_bounds("ACGTGGGGGG", "I" * 10, cutoff=20, verbose=False)
    assert result["index"] == 4
    assert result["retained_length"] == 4
    record = dnaio.SequenceRecord("read", "ACGTGGGGGG", "I" * 10)
    assert result["index"] == qualtrim.nextseq_trim_index(record, 20, 33)


def test_expand_adapter_braces_hand_examples_and_errors() -> None:
    pytest.importorskip("cutadapt.parser")
    assert cut.expand_adapter_braces("TGA{5}CT", verbose=False)["expanded"] == "TGAAAAACT"
    assert cut.expand_adapter_braces("AC{3}G{2}", verbose=False)["expanded"] == "ACCCGG"
    with pytest.raises(ConfigurationError):
        cut.expand_adapter_braces("{5}", verbose=False)
    with pytest.raises(ConfigurationError, match="max_length"):
        cut.expand_adapter_braces("A{1000}", max_length=100, verbose=False)


def test_read_adapter_fasta_reads_local_file_in_order(tmp_path: Path) -> None:
    parser = pytest.importorskip("cutadapt.parser")
    fasta = tmp_path / "adapters.fa"
    fasta.write_text(">ad1 description\nACGTACGT\n>ad2\nTTGGAACC\n", encoding="utf-8")

    result = cut.read_adapter_fasta(fasta, verbose=False)
    assert result["records"] == [
        {"name": "ad1", "sequence": "ACGTACGT"},
        {"name": "ad2", "sequence": "TTGGAACC"},
    ]
    assert result["records"] == [
        {"name": name, "sequence": sequence}
        for name, sequence in parser.read_adapters_fasta(str(fasta))
    ]

    bad = tmp_path / "not_fasta.txt"
    bad.write_text("ACGT\nTTTT\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        cut.read_adapter_fasta(bad, verbose=False)
    with pytest.raises(ConfigurationError, match="max_records"):
        cut.read_adapter_fasta(fasta, max_records=1, verbose=False)


def test_parse_adapter_spec_regular_linked_and_file_notation(tmp_path: Path) -> None:
    parser = pytest.importorskip("cutadapt.parser")

    front = cut.parse_adapter_spec("ACGTN", adapter_type="front", verbose=False)
    assert front["adapters"][0] == {
        "name": "1",
        "adapter_class": "FrontAdapter",
        "description": "regular 5'",
        "linked": False,
        "sequence": "ACGTN",
        "max_error_rate": 0.1,
        "min_overlap": 3,
        "indels": True,
        "read_wildcards": False,
        "adapter_wildcards": True,
    }

    linked = cut.parse_adapter_spec("ACGT...TTTT", adapter_type="back", verbose=False)
    row = linked["adapters"][0]
    assert (row["linked"], row["front_sequence"], row["back_sequence"]) == (True, "ACGT", "TTTT")

    fasta = tmp_path / "adapters.fa"
    fasta.write_text(">ad1\nACGTACGT\n>ad2\nTTGGAACC\n", encoding="utf-8")
    from_file = cut.parse_adapter_spec(
        f"file:{fasta};max_errors=0.2", adapter_type="back", verbose=False
    )
    assert [row["name"] for row in from_file["adapters"]] == ["ad1", "ad2"]
    assert from_file["adapters"][0]["max_error_rate"] == pytest.approx(0.2)
    direct = parser.make_adapter("ACGTACGT", "back", {"max_errors": 0.2, "min_overlap": 3})
    assert direct.max_error_rate == pytest.approx(0.2)

    with pytest.raises(ConfigurationError):
        cut.parse_adapter_spec("ACGT...TTTT", adapter_type="anywhere", verbose=False)


def test_cutadapt_missing_extra_names_bio_cutadapt(monkeypatch: pytest.MonkeyPatch) -> None:
    _hide_package(monkeypatch, "cutadapt")
    with pytest.raises(BackendUnavailableError, match="bio-cutadapt"):
        cut.hamming_neighbours("ACGT", 1, verbose=False)


# umi-tools ----------------------------------------------------------------------------------


def test_average_umi_distance_matches_hand_and_upstream() -> None:
    umi_methods = pytest.importorskip("umi_tools.umi_methods")
    result = umi.average_umi_distance(["AAAA", "AAAT", "CCCC"], verbose=False)
    # Pairs: AAAA-AAAT 1, AAAA-CCCC 4, AAAT-CCCC 4; mean 3.0.
    assert result["average_distance"] == pytest.approx(3.0)
    assert result["pair_count"] == 3
    upstream = umi_methods.get_average_umi_distance([b"AAAA", b"AAAT", b"CCCC"])
    assert result["average_distance"] == pytest.approx(float(upstream))
    single = umi.average_umi_distance(["AAAA"], verbose=False)
    assert single["average_distance"] is None and single["pair_count"] == 0


def test_find_umi_neighbours_matches_brute_force_on_random_umis() -> None:
    umi_methods = pytest.importorskip("umi_tools.umi_methods")
    rng = random.Random(7)
    umis = sorted({"".join(rng.choice("ACGT") for _ in range(5)) for _ in range(60)})
    for max_distance in (1, 2):
        result = umi.find_umi_neighbours(umis, max_distance=max_distance, verbose=False)
        brute = [
            (first, second, int(umi_methods.edit_distance(first.encode(), second.encode())))
            for first, second in itertools.combinations(umis, 2)
        ]
        expected = [(a, b, d) for a, b, d in brute if d <= max_distance]
        observed = [(row["umi_a"], row["umi_b"], row["distance"]) for row in result["pairs"]]
        assert observed == expected
        assert len(expected) > 0


def test_find_umi_neighbours_small_hand_example() -> None:
    pytest.importorskip("umi_tools.umi_methods")
    result = umi.find_umi_neighbours(["AAAA", "AAAT", "CCCC", "AAAG", "AAAA"], verbose=False)
    assert result["unique_umis"] == 4
    assert [(row["umi_a"], row["umi_b"]) for row in result["pairs"]] == [
        ("AAAA", "AAAG"),
        ("AAAA", "AAAT"),
        ("AAAG", "AAAT"),
    ]


def test_umi_connected_component_matches_upstream_traversal() -> None:
    network = pytest.importorskip("umi_tools.network")
    adjacency = {"A": ["B"], "B": ["A", "C"], "D": ["E"]}
    result = umi.umi_connected_component(adjacency, "A", verbose=False)
    assert result["component"] == ["A", "B", "C"]
    assert result["size"] == 3
    assert result["component"] == sorted(
        network.breadth_first_search("A", {"A": {"B"}, "B": {"A", "C"}, "C": {"B"}, "D": {"E"}})
    )
    assert umi.umi_connected_component(adjacency, "D", verbose=False)["component"] == ["D", "E"]


def test_mask_low_quality_umi_uses_threshold_inclusively() -> None:
    pytest.importorskip("umi_tools.extract_methods")
    # '!' is Q0 and '+' is Q10; a base is masked only when its quality is below the threshold.
    result = umi.mask_low_quality_umi("ACGT", "II!I", threshold=10, verbose=False)
    assert result["masked_umi"] == "ACNT"
    assert result["masked_positions"] == [2]
    boundary = umi.mask_low_quality_umi("AC", "+!", threshold=10, verbose=False)
    assert boundary["masked_umi"] == "AN"
    assert boundary["masked_positions"] == [1]


def test_extract_barcode_positions_hand_example_and_overlap_check() -> None:
    pytest.importorskip("umi_tools.extract_methods")
    result = umi.extract_barcode_positions(
        "ACGTACGT",
        "IIIIHHHH",
        umi_positions=[0, 1],
        cell_positions=[2],
        discard_positions=[7],
        verbose=False,
    )
    assert result["read_sequence"] == "TACG"
    assert result["read_qualities"] == "IHHH"
    assert (result["umi_sequence"], result["umi_qualities"]) == ("AC", "II")
    assert (result["cell_sequence"], result["cell_qualities"]) == ("G", "I")

    kept = umi.extract_barcode_positions(
        "ACGTACGT",
        "IIIIHHHH",
        umi_positions=[0, 1],
        cell_positions=[2],
        discard_positions=[7],
        retain_umi=True,
        verbose=False,
    )
    # With retain_umi upstream keeps the UMI bases and also copies the qualities of every
    # retained base into umi_qualities; the wrapper reports that value unchanged.
    assert kept["read_sequence"] == "ACTACG"
    assert kept["umi_qualities"] == "IIIHHH"

    with pytest.raises(ConfigurationError, match="overlap"):
        umi.extract_barcode_positions(
            "ACGT", "IIII", umi_positions=[0], cell_positions=[0], verbose=False
        )


def test_tag_read_name_with_umi_appends_before_first_space() -> None:
    pytest.importorskip("umi_tools.extract_methods")
    tagged = umi.tag_read_name_with_umi(
        "read1 extra words", "ACGT", cell_barcode="TTGG", verbose=False
    )
    assert tagged["tagged_name"] == "read1_TTGG_ACGT extra words"
    plain = umi.tag_read_name_with_umi("read1", "ACGT", verbose=False)
    assert plain["tagged_name"] == "read1_ACGT"
    assert plain["cell_barcode"] is None


def test_cigar_splice_offset_matches_hand_and_upstream() -> None:
    sam = pytest.importorskip("umi_tools.sam_methods")
    assert umi.cigar_splice_offset("10M100N20M", verbose=False)["splice_offset"] == 10
    assert umi.cigar_splice_offset("5S10M100N20M", verbose=False)["splice_offset"] == 15
    assert umi.cigar_splice_offset("30M", verbose=False)["splice_offset"] is None
    reverse = umi.cigar_splice_offset("10M100N20M", reverse=True, verbose=False)
    assert reverse["splice_offset"] == 20
    upstream = sam.find_splice([(0, 10), (3, 100), (0, 20)])
    assert upstream == 10
    with pytest.raises(ConfigurationError):
        umi.cigar_splice_offset("10Q", verbose=False)


def test_parse_umi_read_name_cell_and_umi_modes() -> None:
    pytest.importorskip("umi_tools.sam_methods")
    both = umi.parse_umi_read_name("read1_CELL1_UMI1", cell_barcode=True, verbose=False)
    assert (both["cell_barcode"], both["umi"]) == ("CELL1", "UMI1")
    only = umi.parse_umi_read_name("read1_UMI1", verbose=False)
    assert (only["cell_barcode"], only["umi"]) == (None, "UMI1")
    with pytest.raises(ConfigurationError):
        umi.parse_umi_read_name("UMI1", cell_barcode=True, verbose=False)


def test_barcode_indel_matches_reports_ambiguity() -> None:
    pytest.importorskip("umi_tools.whitelist_methods")
    ambiguous = umi.barcode_indel_matches("ACGT", ["ACGA", "ACGTT", "TTTT"], verbose=False)
    assert ambiguous["matches"] == ["ACGA", "ACGTT"]
    assert ambiguous["ambiguous"] is True
    single = umi.barcode_indel_matches("ACGT", ["ACGA", "TTTT"], verbose=False)
    assert single["matches"] == ["ACGA"] and single["ambiguous"] is False
    exact = umi.barcode_indel_matches("ACGT", ["ACGT", "TTTT"], verbose=False)
    assert exact["exact_in_whitelist"] is True and exact["matches"] == []


def test_read_barcode_whitelist_modes_and_error_maps(tmp_path: Path) -> None:
    whitelist_methods = pytest.importorskip("umi_tools.whitelist_methods")
    tsv = tmp_path / "whitelist.tsv"
    tsv.write_text("# comment\nAAAA\tCCCC,CCCG\nGGGG\tTTTT\n", encoding="utf-8")

    plain = umi.read_barcode_whitelist(tsv, verbose=False)
    assert plain["barcodes"] == ["AAAA", "GGGG"]
    assert plain["error_map"] is None

    from_file = umi.read_barcode_whitelist(tsv, error_correction="file", verbose=False)
    assert from_file["error_map"] == {"CCCC": "AAAA", "CCCG": "AAAA", "TTTT": "GGGG"}

    derived = umi.read_barcode_whitelist(tsv, error_correction="derived", verbose=False)
    # Two barcodes x four positions x four replacements (three bases plus N) = 32 entries.
    assert len(derived["error_map"]) == 32
    assert derived["error_map"]["CAAA"] == "AAAA"
    assert derived["error_map"]["TGGG"] == "GGGG"
    _whitelist, upstream_map = whitelist_methods.getUserDefinedBarcodes(
        str(tsv), None, deriveErrorCorrection=True, threshold=1
    )
    assert derived["error_map"] == {str(k): str(v) for k, v in sorted(upstream_map.items())}

    blank = tmp_path / "blank.tsv"
    blank.write_text("AAAA\n\nGGGG\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        umi.read_barcode_whitelist(blank, verbose=False)


def test_call_cell_whitelist_matches_direct_upstream_calls() -> None:
    whitelist_methods = pytest.importorskip("umi_tools.whitelist_methods")
    counts = {"AAAA": 100, "AAAT": 90, "CCCC": 10, "GGGG": 9, "TTTT": 1, "ACGT": 50}
    upstream = Counter(counts)

    distance = umi.call_cell_whitelist(counts, knee_method="distance", verbose=False)
    direct, _map = whitelist_methods.getCellWhitelist(upstream, knee_method="distance")
    assert distance["whitelist"] == sorted(direct)
    assert distance["cell_count"] == 3

    density = umi.call_cell_whitelist(counts, knee_method="density", verbose=False)
    direct_density, _none = whitelist_methods.getCellWhitelist(upstream, knee_method="density")
    assert density["whitelist"] == sorted(direct_density)

    with pytest.raises(ConfigurationError, match="expect_cells"):
        umi.call_cell_whitelist(counts, knee_method="distance", expect_cells=3, verbose=False)


def test_detect_bam_tags_reports_presence_of_nh_x0_xt(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    pytest.importorskip("umi_tools.dedup")
    sam = tmp_path / "reads.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\n"
        "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:1\n"
        "r2\t0\tchr1\t9\t60\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:1\n",
        encoding="utf-8",
    )
    result = umi.detect_bam_tags(sam, verbose=False)
    assert result["tags"] == {"NH": True, "X0": False, "XT": False}
    with pytest.raises(ConfigurationError):
        umi.detect_bam_tags(tmp_path / "missing.txt", verbose=False)


def test_count_umis_per_gene_cell_and_umi_modes(tmp_path: Path) -> None:
    pytest.importorskip("umi_tools.sam_methods")
    table = tmp_path / "assignments.tsv"
    table.write_text(
        "r1_CELL1_UMI1\tGENE_A\n"
        "r2_CELL1_UMI1\tGENE_A\n"
        "r3_CELL1_UMI2\tGENE_A\n"
        "r4_CELL2_UMI1\tGENE_B\n"
        "r5_CELL1_UMI1\tGENE_B\n",
        encoding="utf-8",
    )
    cells = umi.count_umis_per_gene(table, cell_barcode=True, verbose=False)
    assert cells["genes"] == [
        {"gene": "GENE_A", "reads": 3, "unique_umis": 2},
        {"gene": "GENE_B", "reads": 2, "unique_umis": 2},
    ]
    umis = umi.count_umis_per_gene(table, verbose=False)
    assert umis["genes"][1] == {"gene": "GENE_B", "reads": 2, "unique_umis": 1}


def test_count_umis_per_gene_merges_unsorted_groups_and_rejects_bad_rows(tmp_path: Path) -> None:
    pytest.importorskip("umi_tools.sam_methods")
    unsorted = tmp_path / "unsorted.tsv"
    unsorted.write_text("r1_UMI1\tGENE_A\nr2_UMI2\tGENE_B\nr3_UMI3\tGENE_A\n", encoding="utf-8")
    result = umi.count_umis_per_gene(unsorted, verbose=False)
    assert result["genes"] == [
        {"gene": "GENE_A", "reads": 2, "unique_umis": 2},
        {"gene": "GENE_B", "reads": 1, "unique_umis": 1},
    ]
    bad = tmp_path / "bad.tsv"
    bad.write_text("r1_UMI1\tGENE_A\textra\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="two non-empty"):
        umi.count_umis_per_gene(bad, verbose=False)


def test_umi_tools_missing_extra_names_bio_umi_tools(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("pysam")
    _hide_package(monkeypatch, "umi_tools")
    with pytest.raises(BackendUnavailableError, match="bio-umi-tools"):
        umi.average_umi_distance(["AAAA", "AAAT"], verbose=False)


# Determinism across hash seeds -------------------------------------------------------------

_DETERMINISM_SCRIPT = """
import json
from pathlib import Path
from dnakit.sequencing_qc import cutadapt_functions as cut
from dnakit.sequencing_qc import umi_tools_functions as umi
work = Path(__import__("sys").argv[1])
table = work / "assign.tsv"
table.write_text("r1_C1_U1\\tG_A\\nr2_C2_U2\\tG_B\\nr3_C1_U3\\tG_A\\n", encoding="utf-8")
counts = {"AAAA": 100, "AAAT": 90, "CCCC": 10, "GGGG": 9, "TTTT": 1, "ACGT": 50}
out = {
    "kmers": cut.adapter_kmer_positions("ACGTACGTTG", min_overlap=3, verbose=False),
    "table": cut.adapter_kmer_table("ACGTACGTTG", verbose=False)["rows"],
    "chunks": cut.sequence_kmer_chunks("ACGTACGTACG", 4, verbose=False),
    "neighbours": umi.find_umi_neighbours(
        ["AAAAT", "AAATT", "CCCCA", "CCCAA", "GGGGG"], max_distance=2, verbose=False
    ),
    "component": umi.umi_connected_component(
        {"A": ["B", "C"], "B": ["A"], "C": ["A", "D"]}, "B", verbose=False
    ),
    "genes": umi.count_umis_per_gene(table, cell_barcode=True, verbose=False),
    "derived": umi.read_barcode_whitelist(
        work / "wl.tsv", error_correction="derived", verbose=False
    ),
    "density": umi.call_cell_whitelist(counts, knee_method="density", verbose=False),
}
print(json.dumps(out, sort_keys=True))
"""


def test_grouped_results_are_identical_across_hash_seeds(tmp_path: Path) -> None:
    pytest.importorskip("cutadapt.kmer_heuristic")
    pytest.importorskip("umi_tools.whitelist_methods")
    (tmp_path / "wl.tsv").write_text("AAAA\tCCCC\nGGGG\tTTTT\n", encoding="utf-8")
    outputs = []
    for seed in ("1", "2"):
        env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(SRC)}
        done = subprocess.run(
            [sys.executable, "-c", _DETERMINISM_SCRIPT, str(tmp_path)],
            capture_output=True,
            text=True,
            env=env,
            check=True,
            timeout=120,
        )
        outputs.append(json.loads(done.stdout))
    assert outputs[0] == outputs[1]
