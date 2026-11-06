"""Regression tests for the epigenome, genome3d, design_crispr (primers, data storage) and regulation fixes.

Expected numbers are derived by hand; nearest-neighbour Tm values are compared with ``dnakit.thermodynamics`` directly.
"""

from __future__ import annotations

import pytest

from dnakit import design_crispr as crispr
from dnakit import epigenome as epi
from dnakit import genome3d as g3d
from dnakit.core import DNASequence
from dnakit.design_crispr import _impl_local as crispr_impl
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.thermodynamics import ThermodynamicConditions, melting_temperature


def ok(envelope: dict) -> object:
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# ----------------------------------------------------------------------------- epigenome interval_overlap
def test_interval_overlap_counts_each_input_interval_once() -> None:
    a = "chr1\t100\t300\nchr1\t400\t450\nchr1\t500\t600\nchr2\t0\t10"
    b = "chr1\t120\t150\nchr1\t200\t250\nchr1\t450\t500\nchr1\t590\t700"
    result = ok(epi.calculate(a, options={"entity": "interval_overlap", "b": b}))
    # chr1:100-300 spans two b intervals (counted once), 400-450 and 450-500 are book-ended (no shared base),
    # 500-600 overlaps 590-700, chr2 has no b interval
    assert result["a_intervals_overlapping_b"] == 2
    assert result["overlap_bp"] == 30 + 50 + 10  # 120-150, 200-250 and 590-600
    assert result["a_bp"] == 200 + 50 + 100 + 10
    assert result["b_bp"] == 30 + 50 + 50 + 110


# ----------------------------------------------------------------------------- genome3d TAD edge domains
def _domains(n: int, k: int) -> list[list[float]]:
    size = n // k
    return [
        [10.0 / (1 + abs(i - j)) if i // size == j // size else 1.0 for j in range(n)]
        for i in range(n)
    ]


def test_tads_include_edge_domains_by_default() -> None:
    result = ok(g3d.annotate(_domains(15, 3), options={"window": 2}))
    assert result["boundary_bins"] == [5, 10]
    assert result["tads"] == [
        {"start_bin": 0, "end_bin": 5, "edge": True},
        {"start_bin": 5, "end_bin": 10, "edge": False},
        {"start_bin": 10, "end_bin": 15, "edge": True},
    ]


def test_tads_single_boundary_gives_two_edge_domains() -> None:
    result = ok(g3d.annotate(_domains(16, 2), options={"window": 2}))
    assert result["boundary_bins"] == [8]
    assert result["tads"] == [
        {"start_bin": 0, "end_bin": 8, "edge": True},
        {"start_bin": 8, "end_bin": 16, "edge": True},
    ]
    legacy = ok(g3d.annotate(_domains(16, 2), options={"window": 2, "include_edge_domains": False}))
    assert legacy["tads"] == []


def test_tads_without_boundary_is_one_edge_domain() -> None:
    flat = [[1.0] * 12 for _ in range(12)]
    result = ok(g3d.annotate(flat, options={"window": 2}))
    assert result["boundary_bins"] == []
    assert result["tads"] == [{"start_bin": 0, "end_bin": 12, "edge": True}]
    assert (
        ok(g3d.annotate(flat, options={"window": 2, "include_edge_domains": False}))["tads"] == []
    )


def test_tads_interior_only_option_and_validation() -> None:
    result = ok(g3d.annotate(_domains(15, 3), options={"window": 2, "include_edge_domains": False}))
    assert result["tads"] == [{"start_bin": 5, "end_bin": 10, "edge": False}]
    with pytest.raises(ConfigurationError):
        g3d.annotate(_domains(15, 3), options={"include_edge_domains": "yes"})


# ----------------------------------------------------------------------------- calculate_primer_tm
def test_wallace_values_by_hand() -> None:
    # ATGCGTAC: A,T,T,A = 4 weak -> 8; G,C,G,C = 4 strong -> 16; total 24
    row = ok(crispr.calculate_primer_tm("ATGCGTAC"))["primers"][0]
    assert (row["tm_celsius"], row["method"], row["length"], row["gc_fraction"]) == (
        24.0,
        "wallace",
        8,
        0.5,
    )
    assert (
        ok(crispr.calculate_primer_tm("AAAAAAAAAAAAA"))["primers"][0]["tm_celsius"] == 26.0
    )  # 13 nt: still Wallace
    explicit = ok(crispr.calculate_primer_tm("a t g c", options={"method": "wallace"}))["primers"][
        0
    ]
    assert explicit["sequence"] == "ATGC" and explicit["tm_celsius"] == 12.0


def test_auto_switches_to_nearest_neighbour_above_13_nt() -> None:
    seq = "ACGTACGTACGTAC"  # 14 nt
    row = ok(crispr.calculate_primer_tm(seq))["primers"][0]
    expected = melting_temperature(
        DNASequence(seq), method="nearest_neighbor", verbose=False
    ).tm_celsius
    assert row["method"] == "nearest_neighbor"
    assert row["tm_celsius"] == pytest.approx(expected, abs=0.006)
    assert row["tm_celsius"] != 2 * 7 + 4 * 7


def test_nearest_neighbour_conditions_match_thermodynamics() -> None:
    seq = "ACGTACGTACGTACGTACGT"
    options = {
        "method": "nearest_neighbor",
        "sodium_molar": 0.1,
        "potassium_molar": 0.02,
        "strand_concentration_molar": 5e-7,
    }
    result = ok(crispr.calculate_primer_tm(seq, options=options))
    conditions = ThermodynamicConditions(
        sodium_molar=0.1, potassium_molar=0.02, strand_concentration_molar=5e-7
    )
    expected = melting_temperature(
        DNASequence(seq), method="nearest_neighbor", conditions=conditions, verbose=False
    ).tm_celsius
    assert result["primers"][0]["tm_celsius"] == pytest.approx(expected, abs=0.006)
    assert result["conditions"] == {
        "sodium_molar": 0.1,
        "potassium_molar": 0.02,
        "strand_concentration_molar": 5e-7,
    }
    default = ok(crispr.calculate_primer_tm(seq, options={"method": "nearest_neighbor"}))[
        "primers"
    ][0]["tm_celsius"]
    assert default < result["primers"][0]["tm_celsius"]  # more salt, more Tm


def test_primer_tm_input_forms_and_errors() -> None:
    rows = ok(crispr.calculate_primer_tm({"fwd": "ATGC", "rev": "GGGG"}))["primers"]
    assert [(r["id"], r["tm_celsius"]) for r in rows] == [("fwd", 12.0), ("rev", 16.0)]
    rows = ok(crispr.calculate_primer_tm([{"id": "p", "sequence": "AATT"}, "CCCC"]))["primers"]
    assert [(r["id"], r["tm_celsius"]) for r in rows] == [("p", 8.0), ("primer2", 16.0)]
    assert (
        ok(crispr.calculate_primer_tm({"id": "one", "sequence": "ATAT"}))["primers"][0]["id"]
        == "one"
    )
    for bad in ("ACGN", "ACGR", "", "ACGU"):
        with pytest.raises(InputFormatError):
            crispr.calculate_primer_tm(bad)
    with pytest.raises(ConfigurationError):
        crispr.calculate_primer_tm("A" * 70, options={"method": "nearest_neighbor"})
    with pytest.raises(ConfigurationError):
        crispr.calculate_primer_tm("ACGT", options={"method": "magic"})
    assert (
        ok(crispr.calculate_primer_tm("ACGU", options={"rna": True}))["primers"][0]["sequence"]
        == "ACGT"
    )


# ----------------------------------------------------------------------------- evaluate_primer_pair
def test_evaluate_primer_pair_values_by_hand() -> None:
    result = ok(
        crispr.evaluate_primer_pair(
            {"forward": "ATGCATGCAT", "reverse": "GGGGGGGGGG"},
            options={"method": "wallace", "min_length": 10},
        )
    )
    fwd, rev = result["forward"], result["reverse"]
    assert (fwd["tm_celsius"], rev["tm_celsius"]) == (28.0, 40.0)  # 2*6+4*4 and 4*10
    assert result["delta_tm_celsius"] == 12.0
    assert (fwd["length"], fwd["gc_fraction"], rev["gc_fraction"]) == (10, 0.4, 1.0)
    assert (fwd["gc_clamp"], rev["gc_clamp"]) == (2, 5)  # TGCAT -> G, C ; GGGGG
    assert fwd["longest_self_complementary"] == 10  # ATGCATGCAT is its own reverse complement
    assert rev["longest_self_complementary"] == 0
    assert (
        result["longest_cross_complementary"] == 1
    )  # only a single C of CCCCCCCCCC matches ATGCATGCAT
    flagged = {(f["flag"], f["primer"]) for f in result["flags"]}
    assert flagged == {
        ("tm_mismatch", None),
        ("gc_out_of_range", "reverse"),
        ("gc_clamp_out_of_range", "reverse"),
        ("self_complementarity", "forward"),
    }
    assert result["pass"] is False


def test_evaluate_primer_pair_does_not_mix_in_key_names_and_accepts_lists() -> None:
    opts = {"method": "wallace", "min_length": 1}
    by_map = ok(crispr.evaluate_primer_pair({"forward": "ATGC", "reverse": "GGCC"}, options=opts))
    by_list = ok(crispr.evaluate_primer_pair(["ATGC", "GGCC"], options=opts))
    assert by_map == by_list
    assert by_map["forward"]["length"] == 4 and by_map["forward"]["tm_celsius"] == 12.0
    assert by_map["reverse"]["tm_celsius"] == 16.0 and by_map["delta_tm_celsius"] == 4.0
    for bad in (
        "ATGC",
        ["ATGC"],
        {"forward": "ATGC"},
        {"forward": "ATGC", "reverse": "ATGN"},
        {"fwd": "A", "rev": "C"},
    ):
        with pytest.raises(InputFormatError):
            crispr.evaluate_primer_pair(bad)


def test_evaluate_primer_pair_cross_complementarity_and_thresholds() -> None:
    options = {
        "method": "wallace",
        "min_length": 1,
        "min_gc": 0.0,
        "max_gc": 1.0,
        "min_gc_clamp": 0,
        "max_gc_clamp": 5,
        "max_delta_tm": 100.0,
    }
    result = ok(crispr.evaluate_primer_pair(["AAAACCCC", "GGGGTTTT"], options=options))
    assert result["longest_cross_complementary"] == 8  # revcomp(GGGGTTTT) == AAAACCCC
    assert [f["flag"] for f in result["flags"]] == ["cross_complementarity"]
    relaxed = ok(
        crispr.evaluate_primer_pair(
            ["AAAACCCC", "GGGGTTTT"], options={**options, "max_cross_complementary": 8}
        )
    )
    assert relaxed["flags"] == [] and relaxed["pass"] is True


def test_evaluate_primer_pair_nearest_neighbour_default() -> None:
    fwd, rev = "GACGTTAGCATCAGGCTAAC", "CCATGTTCAGCTAGGATCAT"
    result = ok(crispr.evaluate_primer_pair([fwd, rev]))
    for role, seq in (("forward", fwd), ("reverse", rev)):
        expected = melting_temperature(
            DNASequence(seq), method="nearest_neighbor", verbose=False
        ).tm_celsius
        assert result[role]["method"] == "nearest_neighbor"
        assert result[role]["tm_celsius"] == pytest.approx(expected, abs=0.006)
    assert result["delta_tm_celsius"] == pytest.approx(
        abs(result["forward"]["tm_celsius"] - result["reverse"]["tm_celsius"]), abs=0.006
    )
    assert result["forward"]["gc_clamp"] == 2  # last five bases "CTAAC": C, C
    assert result["thresholds"]["max_delta_tm"] == 5.0


# ----------------------------------------------------------------------------- encode_data_into_dna
def test_bytes_round_trip_all_256_values(tmp_path) -> None:
    data = bytes(range(256))
    dna = ok(crispr.encode_data_into_dna(data))
    assert len(dna) == 4 * 256 and set(dna) <= set("ACGT")
    assert dna[:4] == "AAAA" and dna[4:8] == "AAAC" and dna[-4:] == "TTTT"  # 0x00, 0x01 and 0xff
    assert crispr_impl.dna_to_bytes(dna) == data
    assert ok(crispr.encode_data_into_dna(bytearray(data))) == dna
    path = tmp_path / "blob.bin"
    path.write_bytes(data)
    assert ok(crispr.encode_data_into_dna(path)) == dna
    assert ok(crispr.encode_data_into_dna(str(path), options={"source": "path"})) == dna
    assert (
        ok(crispr.encode_data_into_dna(str(path))) != dna
    )  # a plain string is text, not a path unless asked


def test_text_encoding_unchanged_and_matches_runtime_decoder() -> None:
    from dnakit import _candidate_runtime as runtime

    assert ok(crispr.encode_data_into_dna("hi")) == "CGGACGGC"
    assert ok(crispr.encode_data_into_dna(b"hi")) == "CGGACGGC"
    dna = ok(crispr.encode_data_into_dna("Hi! é"))
    assert crispr_impl.dna_to_bytes(dna).decode("utf-8") == "Hi! é"
    assert (
        runtime._local_result("design_crispr", "decode_dna_storage", dna, {}) == "Hi! é"
    )  # same A/C/G/T = 0/1/2/3 code
    assert ok(crispr.encode_data_into_dna({"b": 1, "a": 2})) == crispr_impl.bytes_to_dna(
        b'{"a": 2, "b": 1}'
    )


def test_encode_errors(tmp_path) -> None:
    with pytest.raises(InputFormatError):
        crispr.encode_data_into_dna(tmp_path / "missing.bin")
    with pytest.raises(InputFormatError):
        crispr.encode_data_into_dna("not a file", options={"source": "path"})
    with pytest.raises(InputFormatError):
        crispr.encode_data_into_dna(12)
    with pytest.raises(ConfigurationError):
        crispr.encode_data_into_dna(b"\0" * 1_000_001)
    with pytest.raises(InputFormatError):
        crispr_impl.dna_to_bytes("ACG")


# ----------------------------------------------------------------------------- regulation docstring
def test_in_silico_mutagenesis_docstring_describes_nested_positions() -> None:
    from dnakit import regulation as reg

    doc = reg.perform_in_silico_mutagenesis.__doc__ or ""
    assert "positions" in doc and "effects" in doc
    gata = {"A": [0, 10, 0, 10], "C": [0, 0, 0, 0], "G": [10, 0, 0, 0], "T": [0, 0, 10, 0]}
    result = ok(
        reg.perform_in_silico_mutagenesis(
            "AGATAAC", options={"motifs": {"GATA": gata}, "pseudocount": 0.5}
        )
    )
    motif = result["motifs"][0]
    assert "effects" not in motif and set(motif["positions"][0]) == {"position", "ref", "effects"}
