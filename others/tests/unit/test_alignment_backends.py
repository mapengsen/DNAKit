"""Contract tests for optional minimap2, edlib, parasail and pyhmmer alignment adapters."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.alignment.backends import (
    align_dna_edlib,
    align_dna_parasail,
    build_dna_hmm,
    build_minimap2_index,
    calculate_edit_distance_matrix,
    fetch_minimap2_reference_region,
    locate_approximate_pattern,
    map_query_sequences_paf,
    read_hmm_profile,
    score_dna_parasail,
    score_dna_parasail_matrix,
    search_dna_hmm,
)
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

CHR1 = (
    "ACGTTGCAAGGCTTAGCCATGGATCCAGTTACGGTAACGTAGCTAGGCTAACGTTGCAAGCTTAGGCCATGCAGTACGATCGATGCTAGC"
    "TAGCAGGTTTAACCGGTTAACGGCAATTGCAGCGTTAGGC"
)
CHR2 = (
    "GGATCCTTAGGCAACCGGTTAAGCTTGCAGTCAGTCAGGATCGGATCCAAGCTAGCTTGCAAGGCTTTAGCAGCATCGAATCGCGATCGTAG"
    "CGTACGTCAGTAGCTGACTGATCGTAGCTAGCTAG"
)
MSA_TEXT = """>seqA
ACGTACGTTAGCAAGGCTTA
>seqB
ACGTACGATAGCAAGGCTTA
>seqC
ACGTTCGTTAGCAAGG-TTA
>seqD
ACGTACGTTAGCTAGGCTTA
"""
TARGET_TEXT = """>t_fwd
GGGGACGTACGTTAGCAAGGCTTAGGGG
>t_rev
CCCCTAAGCCTTGCTAACGTACGTCCCC
>t_none
TTTTTTTTTTTTTTTTTTTTTTTTTTTT
>t_mut
CCCCACGTACGATAGCAAGGCTTACCCC
"""


def _revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


@pytest.fixture
def reference(tmp_path: Path) -> Path:
    path = tmp_path / "reference.fa"
    path.write_text(f">chr1\n{CHR1}\n>chr2\n{CHR2}\n", encoding="utf-8")
    return path


def _records_by_name(result: dict) -> dict[str, dict]:
    return {record["qname"]: record for record in result["records"]}


def test_map_query_sequences_paf_reports_exact_and_reverse_hits(reference: Path) -> None:
    mappy = pytest.importorskip("mappy")
    query = CHR1[10:90]
    result = map_query_sequences_paf(
        reference, {"exact": query, "reverse": _revcomp(query)}, verbose=False
    )
    records = _records_by_name(result)
    # The 80-base query is chr1 positions 11-90 (1-based), i.e. zero-based [10, 90).
    assert records["exact"]["tname"] == "chr1"
    assert (records["exact"]["tstart"], records["exact"]["tend"]) == (10, 90)
    assert (records["exact"]["qstart"], records["exact"]["qend"]) == (0, 80)
    assert records["exact"]["strand"] == "+"
    assert records["exact"]["nmatch"] == 80 and records["exact"]["alen"] == 80
    assert records["exact"]["edit_distance"] == 0 and records["exact"]["cigar"] == "80M"
    assert records["reverse"]["strand"] == "-"
    assert (records["reverse"]["tstart"], records["reverse"]["tend"]) == (10, 90)

    aligner = mappy.Aligner(str(reference), preset="map-ont", n_threads=1)
    upstream = next(iter(aligner.map(query)))
    assert records["exact"]["tstart"] == upstream.r_st
    assert records["exact"]["tend"] == upstream.r_en
    assert records["exact"]["mapq"] == upstream.mapq


def test_map_query_sequences_paf_counts_one_substitution(reference: Path) -> None:
    pytest.importorskip("mappy")
    query = CHR1[10:90]
    position = 30
    substituted = "A" if query[position] != "A" else "C"
    mutated = query[:position] + substituted + query[position + 1 :]
    record = map_query_sequences_paf(reference, {"mut": mutated}, verbose=False)["records"][0]
    assert record["edit_distance"] == 1
    assert record["nmatch"] == 79
    assert record["alen"] == 80
    assert record["cigar"] == "80M"


def test_map_query_sequences_paf_rejects_bad_input(reference: Path) -> None:
    pytest.importorskip("mappy")
    with pytest.raises(ConfigurationError, match="preset must be one of"):
        map_query_sequences_paf(reference, {"q": CHR1[:40]}, preset="bogus", verbose=False)
    with pytest.raises(ConfigurationError, match="max_queries"):
        map_query_sequences_paf(
            reference, {"a": CHR1[:40], "b": CHR1[:40]}, max_queries=1, verbose=False
        )
    with pytest.raises(ConfigurationError, match="without whitespace"):
        map_query_sequences_paf(reference, {"bad name": CHR1[:40]}, verbose=False)
    with pytest.raises(ConfigurationError, match="non-empty"):
        map_query_sequences_paf(reference, {"empty": ""}, verbose=False)
    with pytest.raises(ConfigurationError, match="unsupported"):
        map_query_sequences_paf(reference, {"iupac": "ACGRT"}, verbose=False)
    with pytest.raises(ConfigurationError, match="max_queries"):
        map_query_sequences_paf(reference, {"q": CHR1[:40]}, max_queries=True, verbose=False)


def test_build_minimap2_index_writes_reloadable_index(reference: Path, tmp_path: Path) -> None:
    pytest.importorskip("mappy")
    output = tmp_path / "reference.mmi"
    result = build_minimap2_index(reference, output, verbose=False)
    assert result["contigs"] == ["chr1", "chr2"]
    assert result["contig_count"] == 2
    assert output.is_file() and result["index_bytes"] == output.stat().st_size

    with pytest.raises(ConfigurationError, match="already exists"):
        build_minimap2_index(reference, output, verbose=False)
    rebuilt = build_minimap2_index(reference, output, overwrite=True, verbose=False)
    assert rebuilt["contig_count"] == 2

    record = map_query_sequences_paf(output, {"exact": CHR1[10:90]}, verbose=False)["records"][0]
    assert (record["tstart"], record["tend"]) == (10, 90)


def test_fetch_minimap2_reference_region_slices_one_contig(reference: Path) -> None:
    pytest.importorskip("mappy")
    result = fetch_minimap2_reference_region(reference, "chr2", 0, 12, verbose=False)
    assert result["sequence"] == CHR2[:12]
    assert result["length"] == 12
    middle = fetch_minimap2_reference_region(reference, "chr1", 10, 90, verbose=False)
    assert middle["sequence"] == CHR1[10:90]
    with pytest.raises(ConfigurationError, match="extends past the contig"):
        fetch_minimap2_reference_region(reference, "chr2", 120, 200, verbose=False)
    with pytest.raises(ConfigurationError, match="not in the reference"):
        fetch_minimap2_reference_region(reference, "chrX", 0, 10, verbose=False)
    with pytest.raises(ConfigurationError, match="greater than start"):
        fetch_minimap2_reference_region(reference, "chr1", 10, 10, verbose=False)


def test_missing_mappy_raises_install_hint(
    reference: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(sys.modules, "mappy", None)
    with pytest.raises(BackendUnavailableError, match="bio-mappy"):
        map_query_sequences_paf(reference, {"q": CHR1[:40]}, verbose=False)


def test_align_dna_edlib_matches_hand_calculation() -> None:
    edlib = pytest.importorskip("edlib")
    result = align_dna_edlib("ACGTACGT", "ACGTTACGT", verbose=False)
    # One extra T in the target: a single target-only base, so one 'D' operation.
    assert result["found"] is True
    assert result["edit_distance"] == 1
    assert result["cigar"] == "4=1D4="
    assert (result["matches"], result["deletions"]) == (8, 1)
    assert (result["target_start"], result["target_end"]) == (0, 9)
    assert result["identity"] == pytest.approx(8 / 9)

    upstream = edlib.align("ACGTACGT", "ACGTTACGT", mode="NW", task="path")
    assert result["edit_distance"] == upstream["editDistance"]
    assert result["cigar"] == upstream["cigar"]


def test_align_dna_edlib_iupac_equalities_are_optional() -> None:
    pytest.importorskip("edlib")
    literal = align_dna_edlib("ACRT", "ACGT", verbose=False)
    assert literal["edit_distance"] == 1 and literal["mismatches"] == 1
    ambiguous = align_dna_edlib("ACRT", "ACGT", iupac=True, verbose=False)
    assert ambiguous["edit_distance"] == 0 and ambiguous["cigar"] == "4="
    reversed_roles = align_dna_edlib("ACGT", "ACRT", iupac=True, verbose=False)
    assert reversed_roles["edit_distance"] == 0


def test_align_dna_edlib_infix_mode_reports_target_span() -> None:
    pytest.importorskip("edlib")
    result = align_dna_edlib("ACGT", "TTACGTTT", mode="HW", verbose=False)
    assert result["edit_distance"] == 0
    assert (result["target_start"], result["target_end"]) == (2, 6)
    assert (result["query_start"], result["query_end"]) == (0, 4)


def test_align_dna_edlib_reports_distance_beyond_k_as_not_found() -> None:
    pytest.importorskip("edlib")
    result = align_dna_edlib("AAAA", "CCCC", k=1, verbose=False)
    assert result["found"] is False
    assert result["edit_distance"] is None and result["cigar"] is None


def test_align_dna_edlib_rejects_bad_arguments() -> None:
    with pytest.raises(ConfigurationError, match="mode must be one of"):
        align_dna_edlib("ACGT", "ACGT", mode="SHW", verbose=False)
    with pytest.raises(ConfigurationError, match="k must be an integer"):
        align_dna_edlib("ACGT", "ACGT", k=True, verbose=False)
    with pytest.raises(ConfigurationError, match="iupac must be a boolean"):
        align_dna_edlib("ACGT", "ACGT", iupac=1, verbose=False)


def test_locate_approximate_pattern_lists_lowest_distance_spans() -> None:
    edlib = pytest.importorskip("edlib")
    result = locate_approximate_pattern("ACGT", "TTACGTTTACGTT", max_edits=0, verbose=False)
    assert result["found"] is True and result["edit_distance"] == 0
    assert result["locations"] == [{"start": 2, "end": 6}, {"start": 8, "end": 12}]

    upstream = edlib.align("ACGT", "TTACGTTTACGTT", mode="HW", task="locations", k=0)
    assert [(loc["start"], loc["end"] - 1) for loc in result["locations"]] == upstream["locations"]

    missing = locate_approximate_pattern("GGGG", "TTACGTTTACGTT", max_edits=0, verbose=False)
    assert missing["found"] is False and missing["locations"] == []


def test_calculate_edit_distance_matrix_is_symmetric_and_matches_upstream() -> None:
    edlib = pytest.importorskip("edlib")
    sequences = {"a": "ACGT", "b": "ACGT", "c": "ACGA", "d": "TTTT"}
    result = calculate_edit_distance_matrix(sequences, verbose=False)
    assert result["names"] == ["a", "b", "c", "d"]
    assert result["matrix"] == [[0, 0, 1, 3], [0, 0, 1, 3], [1, 1, 0, 4], [3, 3, 4, 0]]
    assert result["pair_count"] == 6
    for i, (_, left) in enumerate(sequences.items()):
        for j, (_, right) in enumerate(sequences.items()):
            upstream = edlib.align(left, right, mode="NW", task="distance")["editDistance"]
            assert result["matrix"][i][j] == upstream


def test_missing_edlib_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "edlib", None)
    with pytest.raises(BackendUnavailableError, match="bio-edlib"):
        align_dna_edlib("ACGT", "ACGT", verbose=False)


def test_align_dna_parasail_global_matches_hand_score() -> None:
    parasail = pytest.importorskip("parasail")
    identical = align_dna_parasail("ACGT", "ACGT", match=5, mismatch=-4, verbose=False)
    assert identical["score"] == 20 and identical["cigar"] == "4="
    assert identical["identity"] == pytest.approx(1.0)

    # Match 2, mismatch -1, gap_open 5: 8 matches (16) minus one gap base (5) = 11.
    gapped = align_dna_parasail(
        "ACGTACGT", "ACGTTACGT", match=2, mismatch=-1, gap_open=5, gap_extend=1, verbose=False
    )
    assert gapped["score"] == 11
    assert gapped["cigar"] == "3=1D5="
    assert gapped["query_aligned"] == "ACG-TACGT" and gapped["deletions"] == 1

    upstream = parasail.nw_trace_scan_sat(
        "ACGTACGT", "ACGTTACGT", 5, 1, parasail.matrix_create("ACGTN", 2, -1)
    )
    assert gapped["score"] == upstream.score


def test_parasail_gap_cost_is_open_plus_extensions() -> None:
    pytest.importorskip("parasail")
    # Two-base gap: open + 1 * extend. Score 16 - (5 + 1) = 10, then 16 - (10 + 2) = 4.
    short = align_dna_parasail(
        "ACGTTTACGT", "ACGTACGT", match=2, mismatch=-1, gap_open=5, gap_extend=1, verbose=False
    )
    assert short["score"] == 10 and short["cigar"] == "3=2I5="
    long = align_dna_parasail(
        "ACGTTTACGT", "ACGTACGT", match=2, mismatch=-1, gap_open=10, gap_extend=2, verbose=False
    )
    assert long["score"] == 4


def test_parasail_local_region_matches_hand_calculation() -> None:
    pytest.importorskip("parasail")
    result = align_dna_parasail(
        "TTTGATTACATTT",
        "GGGGATTACAGGG",
        mode="local",
        match=2,
        mismatch=-1,
        gap_open=5,
        gap_extend=1,
        verbose=False,
    )
    # GATTACA is query[3:10] and target[3:10]; seven matches at score 2 each.
    assert result["score"] == 14 and result["cigar"] == "7="
    assert (result["query_start"], result["query_end"]) == (3, 10)
    assert (result["target_start"], result["target_end"]) == (3, 10)


def test_parasail_semiglobal_trims_free_end_gaps() -> None:
    parasail = pytest.importorskip("parasail")
    result = align_dna_parasail(
        "GATTACA",
        "TTTGATTACATTT",
        mode="semiglobal",
        match=2,
        mismatch=-1,
        gap_open=5,
        gap_extend=1,
        verbose=False,
    )
    assert result["score"] == 14
    assert (result["query_start"], result["query_end"]) == (0, 7)
    assert (result["target_start"], result["target_end"]) == (3, 10)
    assert result["identity"] == pytest.approx(1.0)
    upstream = parasail.sg_trace_scan_sat(
        "GATTACA", "TTTGATTACATTT", 5, 1, parasail.matrix_create("ACGTN", 2, -1)
    )
    assert result["score"] == upstream.score


def test_score_dna_parasail_end_cells_match_alignment_spans() -> None:
    pytest.importorskip("parasail")
    options = {"match": 2, "mismatch": -1, "gap_open": 5, "gap_extend": 1, "verbose": False}
    for mode, query, target in [
        ("local", "TTTGATTACATTT", "GGGGATTACAGGG"),
        ("semiglobal", "GATTACA", "TTTGATTACATTT"),
    ]:
        aligned = align_dna_parasail(query, target, mode=mode, **options)
        scored = score_dna_parasail(query, target, mode=mode, **options)
        assert scored["score"] == aligned["score"]
        assert scored["query_end"] == aligned["query_end"]
        assert scored["target_end"] == aligned["target_end"]
    full = score_dna_parasail("ACGT", "ACGTTT", **options)
    assert (full["query_end"], full["target_end"]) == (4, 6)


def test_local_alignment_without_similarity_reports_no_region() -> None:
    pytest.importorskip("parasail")
    scored = score_dna_parasail("AAAA", "CCCC", mode="local", verbose=False)
    assert scored["score"] == 0 and scored["query_end"] is None
    aligned = align_dna_parasail("AAAA", "CCCC", mode="local", verbose=False)
    assert aligned["found"] is False and aligned["identity"] is None


def test_score_dna_parasail_matrix_is_symmetric_and_matches_pairwise_calls() -> None:
    pytest.importorskip("parasail")
    sequences = {"x": "ACGTACGT", "y": "ACGTTACGT", "z": "TTTTGGGG"}
    options = {"mode": "local", "match": 2, "mismatch": -1, "gap_open": 5, "gap_extend": 1}
    result = score_dna_parasail_matrix(sequences, verbose=False, **options)
    scores = result["scores"]
    assert scores[0][1] == scores[1][0]
    names = list(sequences)
    for i, left in enumerate(names):
        for j, right in enumerate(names):
            pair = score_dna_parasail(sequences[left], sequences[right], verbose=False, **options)
            assert scores[i][j] == pair["score"]
    assert result["pair_count"] == 6


def test_parasail_dnafull_matrix_uses_builtin_scores() -> None:
    parasail = pytest.importorskip("parasail")
    result = align_dna_parasail("ACGT", "ACGT", matrix="dnafull", verbose=False)
    assert result["score"] == 20 and result["matrix"] == "dnafull"
    upstream = parasail.nw_scan_sat("ACGT", "ACGT", 10, 1, parasail.dnafull)
    assert result["score"] == upstream.score


def test_parasail_rejects_iupac_symbols_and_oversized_pairs() -> None:
    pytest.importorskip("parasail")
    with pytest.raises(ConfigurationError, match="unsupported symbols: R"):
        align_dna_parasail("ACGR", "ACGT", verbose=False)
    with pytest.raises(ConfigurationError, match="max_cells"):
        align_dna_parasail("ACGT", "ACGT", max_cells=10, verbose=False)
    with pytest.raises(ConfigurationError, match="mode must be one of"):
        align_dna_parasail("ACGT", "ACGT", mode="glocal", verbose=False)


def test_missing_parasail_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "parasail", None)
    with pytest.raises(BackendUnavailableError, match="bio-parasail"):
        align_dna_parasail("ACGT", "ACGT", verbose=False)


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_build_dna_hmm_matches_direct_pyhmmer(tmp_path: Path) -> None:
    pyhmmer = pytest.importorskip("pyhmmer")
    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    output = tmp_path / "msa.hmm"
    result = build_dna_hmm(alignment, output_path=output, verbose=False)
    assert result["name"] == "msa"
    assert result["length"] == 20  # four aligned rows of twenty columns
    assert result["nseq"] == 4
    assert output.is_file()

    alphabet = pyhmmer.easel.Alphabet.dna()
    with pyhmmer.easel.MSAFile(str(alignment), digital=True, alphabet=alphabet) as handle:
        msa = handle.read()
    msa.name = b"msa"
    upstream, _profile, _optimized = pyhmmer.plan7.Builder(alphabet).build_msa(
        msa, pyhmmer.plan7.Background(alphabet)
    )
    assert result["consensus"] == upstream.consensus
    assert result["length"] == upstream.M

    with pytest.raises(ConfigurationError, match="already exists"):
        build_dna_hmm(alignment, output_path=output, verbose=False)
    rebuilt = build_dna_hmm(alignment, output_path=output, overwrite=True, verbose=False)
    assert rebuilt["length"] == 20


def test_build_dna_hmm_rejects_unaligned_input_and_bad_names(tmp_path: Path) -> None:
    pytest.importorskip("pyhmmer")
    unaligned = _write(tmp_path, "unaligned.fa", ">a\nACGTAC\n>b\nACGT\n")
    with pytest.raises(ConfigurationError, match="could not read"):
        build_dna_hmm(unaligned, verbose=False)
    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    with pytest.raises(ConfigurationError, match="HMM name must"):
        build_dna_hmm(alignment, name="bad name", verbose=False)


def test_search_dna_hmm_finds_forward_reverse_and_mutated_targets(tmp_path: Path) -> None:
    pyhmmer = pytest.importorskip("pyhmmer")
    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    profile = tmp_path / "msa.hmm"
    build_dna_hmm(alignment, output_path=profile, verbose=False)
    targets = _write(tmp_path, "targets.fa", TARGET_TEXT)

    result = search_dna_hmm(profile, targets, verbose=False)
    assert result["sequences_searched"] == 4
    by_target = {row["target"]: row for row in result["hits"]}
    assert set(by_target) == {"t_fwd", "t_rev", "t_mut"}
    assert (by_target["t_fwd"]["start"], by_target["t_fwd"]["end"]) == (5, 24)
    assert by_target["t_fwd"]["strand"] == "+"
    assert (by_target["t_rev"]["start"], by_target["t_rev"]["end"]) == (5, 24)
    assert by_target["t_rev"]["strand"] == "-"
    assert by_target["t_mut"]["strand"] == "+"

    forward = search_dna_hmm(profile, targets, strand="forward", verbose=False)
    assert sorted(row["target"] for row in forward["hits"]) == ["t_fwd", "t_mut"]
    reverse = search_dna_hmm(profile, targets, strand="reverse", verbose=False)
    assert [row["target"] for row in reverse["hits"]] == ["t_rev"]

    handle = pyhmmer.plan7.HMMFile(str(profile))
    with handle:
        hmm = next(iter(handle))
    seqs = list(pyhmmer.easel.SequenceFile(str(targets), digital=True, alphabet=hmm.alphabet))
    upstream = next(iter(pyhmmer.nhmmer(hmm, seqs, cpus=1, E=10.0)))
    upstream_fwd = next(hit for hit in upstream if hit.name == "t_fwd")
    alignment_fwd = upstream_fwd.domains[0].alignment
    assert (alignment_fwd.target_from, alignment_fwd.target_to) == (5, 24)


def test_search_dna_hmm_rejects_protein_profile_and_limits(tmp_path: Path) -> None:
    pyhmmer = pytest.importorskip("pyhmmer")
    amino = pyhmmer.easel.Alphabet.amino()
    protein = _write(tmp_path, "prot.fa", ">p1\nMKVLA\n>p2\nMKVLG\n")
    with pyhmmer.easel.MSAFile(str(protein), digital=True, alphabet=amino) as handle:
        msa = handle.read()
    msa.name = b"prot"
    hmm, _profile, _optimized = pyhmmer.plan7.Builder(amino).build_msa(
        msa, pyhmmer.plan7.Background(amino)
    )
    protein_hmm = tmp_path / "prot.hmm"
    with protein_hmm.open("xb") as handle:
        hmm.write(handle, binary=False)
    targets = _write(tmp_path, "targets.fa", TARGET_TEXT)
    with pytest.raises(ConfigurationError, match="needs a DNA profile"):
        search_dna_hmm(protein_hmm, targets, verbose=False)

    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    profile = tmp_path / "msa.hmm"
    build_dna_hmm(alignment, output_path=profile, verbose=False)
    with pytest.raises(ConfigurationError, match="max_hits"):
        search_dna_hmm(profile, targets, max_hits=1, verbose=False)
    with pytest.raises(ConfigurationError, match="max_sequences"):
        search_dna_hmm(profile, targets, max_sequences=2, verbose=False)
    with pytest.raises(ConfigurationError, match="evalue must be a number"):
        search_dna_hmm(profile, targets, evalue=0.0, verbose=False)


def test_read_hmm_profile_reports_each_profile(tmp_path: Path) -> None:
    pytest.importorskip("pyhmmer")
    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    first = tmp_path / "first.hmm"
    second = tmp_path / "second.hmm"
    build_dna_hmm(alignment, name="first", output_path=first, verbose=False)
    build_dna_hmm(alignment, name="second", output_path=second, verbose=False)
    combined = tmp_path / "combined.hmm"
    combined.write_text(first.read_text(encoding="utf-8") + second.read_text(encoding="utf-8"))

    result = read_hmm_profile(combined, verbose=False)
    assert result["profile_count"] == 2
    assert [profile["name"] for profile in result["profiles"]] == ["first", "second"]
    assert result["profiles"][0]["length"] == 20
    assert result["profiles"][0]["alphabet"] == "DNA"
    assert result["profiles"][0]["nseq"] == 4

    with pytest.raises(ConfigurationError, match="max_profiles"):
        read_hmm_profile(combined, max_profiles=1, verbose=False)


def test_missing_pyhmmer_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    alignment = _write(tmp_path, "msa.fa", MSA_TEXT)
    monkeypatch.setitem(sys.modules, "pyhmmer", None)
    with pytest.raises(BackendUnavailableError, match="bio-pyhmmer"):
        build_dna_hmm(alignment, verbose=False)
