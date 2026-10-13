"""Contract tests for the G2 wrappers: parasail, edlib and pyhmmer alignment helpers, plus mappy."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.alignment import edlib_functions, mappy_functions, parasail_functions, pyhmmer_functions
from dnakit.alignment.edlib_functions import edlib_nice_alignment
from dnakit.alignment.parasail_functions import (
    parasail_nw_banded_score,
    parasail_sg_db_align,
    parasail_sg_de_align,
    parasail_sg_dx_align,
    parasail_sg_qb_align,
    parasail_sg_qb_db_align,
    parasail_sg_qb_de_align,
    parasail_sg_qe_align,
    parasail_sg_qe_db_align,
    parasail_sg_qe_de_align,
    parasail_sg_qx_align,
    parasail_ssw_align,
)
from dnakit.alignment.pyhmmer_functions import pyhmmer_align_to_hmm
from dnakit.bioinformatics.reads import read_fastx
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.sequence import reverse_complement_dna

# Parasail end-gap modes: the sequence ends that may stay unaligned.
SG_MODES = {
    "qb": parasail_sg_qb_align,
    "qe": parasail_sg_qe_align,
    "qx": parasail_sg_qx_align,
    "db": parasail_sg_db_align,
    "de": parasail_sg_de_align,
    "dx": parasail_sg_dx_align,
    "qb_db": parasail_sg_qb_db_align,
    "qb_de": parasail_sg_qb_de_align,
    "qe_db": parasail_sg_qe_db_align,
    "qe_de": parasail_sg_qe_de_align,
}
FREE_ENDS = {
    "qb": {"query_begin"},
    "qe": {"query_end"},
    "qx": {"query_begin", "query_end"},
    "db": {"database_begin"},
    "de": {"database_end"},
    "dx": {"database_begin", "database_end"},
    "qb_db": {"query_begin", "database_begin"},
    "qb_de": {"query_begin", "database_end"},
    "qe_db": {"query_end", "database_begin"},
    "qe_de": {"query_end", "database_end"},
}
# Core 8-mer flanked by six extra bases on both sides (a 6-base gap costs 10 + 5 * 1 = 15).
CORE = "ACGTTGCA"
FLANK = "GGGGGG"
INSIDE_TARGET = (CORE, FLANK + CORE + FLANK)  # query inside a longer target
OVERHANG_QUERY = (FLANK + CORE + FLANK, CORE)  # query longer than the target


def _matrix(parasail, match=5, mismatch=-4):
    return parasail.matrix_create("ACGTN", match, mismatch)


def _signature_contract(function):
    parameters = inspect.signature(function).parameters
    assert parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["verbose"].default is True


@pytest.mark.parametrize("mode", sorted(SG_MODES))
def test_parasail_sg_modes_match_upstream_and_hand_calculation(mode):
    parasail = pytest.importorskip("parasail")
    wrapper = SG_MODES[mode]
    upstream = getattr(parasail, f"sg_{mode}_trace_scan_sat")

    # The query sits inside the target: each unaligned target flank costs 15 unless its end is free.
    query, target = INSIDE_TARGET
    expected_inside = 40 - 15 * (2 - len(FREE_ENDS[mode] & {"database_begin", "database_end"}))
    result = wrapper(query, target)
    assert (
        result["score"]
        == expected_inside
        == upstream(query, target, 10, 1, _matrix(parasail)).score
    )
    assert set(result["free_end_gaps"]) == FREE_ENDS[mode]
    assert result["operation"] == f"parasail_sg_{mode}_align"

    # The target is shorter: each unaligned query flank costs 15 unless its end is free.
    query, target = OVERHANG_QUERY
    expected_overhang = 40 - 15 * (2 - len(FREE_ENDS[mode] & {"query_begin", "query_end"}))
    result = wrapper(query, target)
    assert result["score"] == expected_overhang
    assert result["score"] == upstream(query, target, 10, 1, _matrix(parasail)).score


def test_parasail_qb_de_overlap_hand_calculation():
    parasail = pytest.importorskip("parasail")
    query, target = "TTTTACGTACGGA", "ACGTACGGATTTT"
    result = parasail_sg_qb_de_align(query, target)
    # The query suffix ACGTACGGA (9 bases) overlaps the target prefix; both end gaps are free.
    assert result["score"] == 45
    assert (result["query_start"], result["query_end"]) == (4, 13)
    assert (result["target_start"], result["target_end"]) == (0, 9)
    # Columns and identity cover the aligned core; the CIGAR also lists the free end gaps.
    assert result["matches"] == 9 and result["alignment_columns"] == 9
    assert result["identity"] == 1.0
    assert result["cigar"] == "4I9=4D"
    assert result["query_aligned"] == "TTTTACGTACGGA----"
    assert result["target_aligned"] == "----ACGTACGGATTTT"
    upstream = parasail.sg_qb_de_trace_scan_sat(query, target, 10, 1, _matrix(parasail))
    assert upstream.score == result["score"]


def test_parasail_both_begins_or_both_ends_cannot_be_skipped_at_once():
    pytest.importorskip("parasail")
    # Both prefixes skipped would score 45; the DP starts on one boundary, so only one is skipped.
    assert parasail_sg_qb_db_align("TTTTACGTACGGA", "CCACGTACGGA")["score"] == 37
    assert parasail_sg_qb_db_align("TTTTACGTACGGA", "ACGTACGGA")["score"] == 45
    assert parasail_sg_qe_de_align("ACGTACGGATTTT", "ACGTACGGAGGGG")["score"] == 32
    assert parasail_sg_qe_de_align("ACGTACGGATTTT", "ACGTACGGA")["score"] == 45


def test_parasail_variants_agree_on_small_inputs():
    pytest.importorskip("parasail")
    query, target = "TTTTACGTACGGA", "ACGTACGGATTTT"
    reference = parasail_sg_qb_de_align(query, target)
    for method in ("scan", "striped", "diag"):
        for width in (8, 16, 32, 64, "sat"):
            result = parasail_sg_qb_de_align(query, target, method=method, width=width)
            assert result["score"] == reference["score"]
            assert result["cigar"] == reference["cigar"]
            assert result["method"] == method and result["width"] == width


def test_parasail_narrow_width_saturation_raises():
    pytest.importorskip("parasail")
    perfect = "ACGT" * 8  # 32 matches score 160, above the 8-bit lane limit of 127
    with pytest.raises(BackendExecutionError, match="saturated"):
        parasail_sg_qb_align(perfect, perfect, width=8)
    assert parasail_sg_qb_align(perfect, perfect, width=16)["score"] == 160


def test_parasail_sg_validation_errors():
    pytest.importorskip("parasail")
    with pytest.raises(ConfigurationError, match="method"):
        parasail_sg_qb_align("ACGT", "ACGT", method="bogus")
    with pytest.raises(ConfigurationError, match="width"):
        parasail_sg_qb_align("ACGT", "ACGT", width=7)
    with pytest.raises(ConfigurationError, match="width"):
        parasail_sg_qb_align("ACGT", "ACGT", width=True)
    with pytest.raises(ConfigurationError, match="unsupported"):
        parasail_sg_qb_align("ACGTR", "ACGT")
    with pytest.raises(ConfigurationError, match="max_cells"):
        parasail_sg_qb_align("ACGTACGT", "ACGTACGT", max_cells=10)
    with pytest.raises(ConfigurationError, match="match"):
        parasail_sg_qb_align("ACGT", "ACGT", match=5000)


def test_parasail_nw_banded_hand_calculation_and_band_restriction():
    parasail = pytest.importorskip("parasail")
    same = "ACGTTGCAAGTC"
    banded = parasail_nw_banded_score(same, same, band=1)
    assert banded["score"] == 5 * len(same) == 60
    assert banded["score"] == parasail.nw_banded(same, same, 10, 1, 1, _matrix(parasail)).score

    # Equal lengths with a 6-base shift: the optimal path leaves a narrow band.
    query = "ACGTTGCAAGTCCATGGATCCAGTTACGGTAACGTTAGCCA"
    target = query[:12] + "TTTTTT" + query[12:30] + query[36:]
    full = parasail.nw(query, target, 10, 1, _matrix(parasail)).score
    narrow = parasail_nw_banded_score(query, target, band=1)["score"]
    assert full == 145 and narrow == 81
    assert parasail_nw_banded_score(query, target, band=6)["score"] == full
    with pytest.raises(ConfigurationError, match="band"):
        parasail_nw_banded_score(query, target, band=0)


def test_parasail_ssw_hand_calculation():
    parasail = pytest.importorskip("parasail")
    result = parasail_ssw_align("ACGTTGCAAGT", "GGGGACGTTGCAAGTGGGG")
    assert result["found"] is True
    assert result["score"] == 11 * 5 == 55
    assert result["cigar"] == "11="
    assert (result["query_start"], result["query_end"]) == (0, 11)
    assert (result["target_start"], result["target_end"]) == (4, 15)
    assert result["identity"] == 1.0
    raw = parasail.ssw("ACGTTGCAAGT", "GGGGACGTTGCAAGTGGGG", 10, 1, _matrix(parasail))
    assert raw.score1 == result["score"]


def test_parasail_ssw_deletion_and_no_hit():
    pytest.importorskip("parasail")
    result = parasail_ssw_align("ACGTTGCAAGTCCATG", "ACGTTGCAAGGGTCCATG")
    # 16 matches (80) minus a 2-base gap (10 + 1 = 11) gives 69.
    assert result["score"] == 69
    assert result["cigar"] == "9=2D7="
    assert result["query_aligned"] == "ACGTTGCAA--GTCCATG"
    assert result["target_aligned"] == "ACGTTGCAAGGGTCCATG"
    assert result["deletions"] == 2 and result["mismatches"] == 0
    miss = parasail_ssw_align("AAAA", "CCCC")
    assert miss["found"] is False and miss["score"] == 0 and miss["query_start"] is None


def test_parasail_ssw_16_bit_limit_raises():
    pytest.importorskip("parasail")
    long_read = "ACGT" * 900  # 3600 matches at score 10 is 36000, above the 16-bit SSW limit
    with pytest.raises(BackendExecutionError, match="16-bit"):
        parasail_ssw_align(long_read, long_read, match=10)


@pytest.mark.parametrize(
    "wrapper",
    [
        parasail_nw_banded_score,
        parasail_ssw_align,
        *SG_MODES.values(),
    ],
)
def test_parasail_wrappers_have_keyword_verbose(wrapper):
    _signature_contract(wrapper)


def test_edlib_nice_alignment_matches_upstream_rows():
    edlib = pytest.importorskip("edlib")
    query, target = "ACGTTGCA", "GGACGTAGCAGG"
    result = edlib_nice_alignment(query, target)
    assert result["found"] is True
    assert result["edit_distance"] == 5
    assert result["cigar"] == "2D4=1X3=2D"
    assert result["query_aligned"] == "--ACGTTGCA--"
    assert result["matched_aligned"] == "--||||.|||--"
    assert result["target_aligned"] == "GGACGTAGCAGG"
    upstream = edlib.getNiceAlignment(
        edlib.align(query, target, mode="NW", task="path"), query, target
    )
    assert upstream == {
        "query_aligned": result["query_aligned"],
        "matched_aligned": result["matched_aligned"],
        "target_aligned": result["target_aligned"],
    }


def test_edlib_nice_alignment_half_global_span_and_gap_symbol():
    pytest.importorskip("edlib")
    result = edlib_nice_alignment("ACGTTGCA", "GGACGTAGCAGG", mode="HW", gap_symbol="*")
    assert result["edit_distance"] == 1
    assert (result["target_start"], result["target_end"]) == (2, 10)
    assert result["query_aligned"] == "ACGTTGCA"
    assert result["matched_aligned"] == "||||.|||"
    assert result["target_aligned"] == "ACGTAGCA"
    gapped = edlib_nice_alignment("ACGT", "ACTT", gap_symbol="*")
    assert gapped["matched_aligned"] == "||.|" and "*" not in gapped["query_aligned"]


def test_edlib_nice_alignment_limits_and_validation():
    pytest.importorskip("edlib")
    missed = edlib_nice_alignment("ACGT", "TTTT", k=0)
    assert missed["found"] is False and missed["query_aligned"] is None
    with pytest.raises(ConfigurationError, match="gap_symbol"):
        edlib_nice_alignment("ACGT", "ACGT", gap_symbol="A")
    with pytest.raises(ConfigurationError, match="mode"):
        edlib_nice_alignment("ACGT", "ACGT", mode="SHW")
    with pytest.raises(ConfigurationError, match="max_length"):
        edlib_nice_alignment("ACGTACGT", "ACGT", max_length=4)


def _dna_hmm(tmp_path: Path, alphabet_name: str = "dna") -> Path:
    pyhmmer = pytest.importorskip("pyhmmer")
    import io

    from pyhmmer.easel import Alphabet

    alphabet = Alphabet.dna() if alphabet_name == "dna" else Alphabet.amino()
    text = (
        ">a\nACGTACGTAC\n>b\nACGTTCGTAC\n>c\nACGAACGTAC\n"
        if alphabet_name == "dna"
        else ">a\nMKVLAAGIVG\n>b\nMKVLSAGIVG\n>c\nMKVIAAGIVG\n"
    )
    msa = pyhmmer.easel.MSAFile(io.BytesIO(text.encode()), digital=True, alphabet=alphabet).read()
    msa.name = b"probe"
    hmm, _profile, _optimized = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    path = tmp_path / f"{alphabet_name}.hmm"
    with path.open("wb") as handle:
        hmm.write(handle, binary=False)
    return path


def test_pyhmmer_align_to_hmm_matches_upstream_and_consensus(tmp_path):
    pyhmmer = pytest.importorskip("pyhmmer")
    profile = _dna_hmm(tmp_path)
    sequences = {"s1": "ACGTACGTAC", "s2": "ACGTTCGTAC", "s3": "GGACGTACGTACGG"}
    result = pyhmmer_align_to_hmm(profile, sequences, trim=True)
    assert result["profile_length"] == 10 and result["consensus"] == "acgtacgtac"
    assert result["sequence_count"] == 3 and result["alignment_columns"] == 10
    assert result["rows"][0] == {"name": "s1", "aligned": "ACGTACGTAC"}
    assert result["rows"][1] == {"name": "s2", "aligned": "ACGTTCGTAC"}

    with pyhmmer.plan7.HMMFile(str(profile)) as handle:
        hmm = next(iter(handle))
    digital = [
        pyhmmer.easel.TextSequence(name=name, sequence=text).digitize(hmm.alphabet)
        for name, text in sequences.items()
    ]
    upstream = pyhmmer.hmmalign(hmm, digital, trim=True)
    assert [str(row) for row in upstream.alignment] == [row["aligned"] for row in result["rows"]]


def test_pyhmmer_align_to_hmm_untrimmed_and_validation(tmp_path):
    pytest.importorskip("pyhmmer")
    profile = _dna_hmm(tmp_path)
    untrimmed = pyhmmer_align_to_hmm(profile, {"s1": "GGACGTACGTACGG"})
    assert untrimmed["rows"][0]["aligned"] == "ggACGTACGTACgg"
    assert untrimmed["alignment_columns"] == 14
    with pytest.raises(ConfigurationError, match="max_sequences"):
        pyhmmer_align_to_hmm(profile, {"a": "ACGT", "b": "ACGT"}, max_sequences=1)
    with pytest.raises(ConfigurationError, match="trim"):
        pyhmmer_align_to_hmm(profile, {"a": "ACGT"}, trim="yes")


def test_pyhmmer_align_to_hmm_rejects_protein_profile(tmp_path):
    pytest.importorskip("pyhmmer")
    protein = _dna_hmm(tmp_path, alphabet_name="amino")
    with pytest.raises(ConfigurationError, match="DNA profile"):
        pyhmmer_align_to_hmm(protein, {"a": "ACGT"})


def test_mappy_duplicates_match_dnakit_functions(tmp_path):
    mappy = pytest.importorskip("mappy")
    assert mappy_functions.__all__ == ()
    assert mappy.revcomp("ACGTNACGTT") == reverse_complement_dna("ACGTNACGTT") == "AACGTNACGT"

    fasta = tmp_path / "reads.fa"
    fasta.write_text(">r1 desc\nACGTNACGT\n>r2\nGGCC\n", encoding="ascii")
    upstream = [(name, seq) for name, seq, *_rest in mappy.fastx_read(str(fasta))]
    dnakit_reads = [(read.id, read.sequence) for read in read_fastx(fasta, format="fasta")]
    assert upstream == dnakit_reads == [("r1", "ACGTNACGT"), ("r2", "GGCC")]


def test_wrapper_modules_export_expected_names():
    assert parasail_functions.__all__ == sorted(
        [
            "parasail_nw_banded_score",
            "parasail_ssw_align",
            *(f"parasail_sg_{mode}_align" for mode in SG_MODES),
        ]
    )
    assert edlib_functions.__all__ == ["edlib_nice_alignment"]
    assert pyhmmer_functions.__all__ == ["pyhmmer_align_to_hmm"]
