"""Tests for dnakit.alignment.biopython_p18_alignment_functions.

Alignments, search reports and MAF files come from _biopython_p18_data.py or are generated in
tmp_path. Nothing uses the network.
"""

from __future__ import annotations

import pytest

pytest.importorskip("Bio")

import _biopython_p18_data as data
from Bio import Align, AlignIO, SearchIO

from dnakit.alignment import biopython_p18_alignment_functions as m
from dnakit.exceptions import ConfigurationError

QUIET = {"verbose": False}
ROWS = {"r1": "ACGT-ACGTA", "r2": "AC-TTACGAA", "r3": "ACGTTACG-A"}


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# ------------------------------------------------------------------ alignment files


@pytest.mark.parametrize(
    "fmt", ["a2m", "clustal", "fasta", "maf", "mauve", "nexus", "phylip", "stockholm"]
)
def test_write_aligned_sequences_is_readable_by_biopython(tmp_path, fmt):
    rows = {"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA", "s3": "ACGTTACG-A"}
    out = tmp_path / f"x.{fmt}"
    result = m.write_aligned_sequences(rows, out, file_format=fmt, **QUIET)
    assert result["sequence_count"] == 3 and result["aligned_columns"] == 10
    parsed = list(Align.parse(str(out), fmt))
    assert len(parsed) == 1
    got = {str(seq.id): str(parsed[0][i]).upper() for i, seq in enumerate(parsed[0].sequences)}
    if fmt == "mauve":
        got = {name: row for name, row in zip(rows, got.values(), strict=True)}
    assert got == rows


def test_write_aligned_sequences_errors(tmp_path):
    out = tmp_path / "x.fasta"
    m.write_aligned_sequences(ROWS, out, file_format="fasta", **QUIET)
    with pytest.raises(ConfigurationError, match="already exists"):
        m.write_aligned_sequences(ROWS, out, file_format="fasta", **QUIET)
    with pytest.raises(ConfigurationError, match="same aligned length"):
        m.write_aligned_sequences(
            {"a": "ACGT", "b": "ACG"}, tmp_path / "y.fasta", file_format="fasta", **QUIET
        )
    with pytest.raises(ConfigurationError):
        m.write_aligned_sequences(ROWS, tmp_path / "z.msf", file_format="msf", **QUIET)
    with pytest.raises(ConfigurationError):
        m.write_aligned_sequences(
            {"a b": "ACGT", "c": "ACGT"}, tmp_path / "w.fasta", file_format="fasta", **QUIET
        )


def test_read_alignment_records_text_formats(tmp_path):
    needle = _write(tmp_path, "needle.txt", data.EMBOSS_NEEDLE)
    emboss = m.read_alignment_records(needle, file_format="emboss", **QUIET)
    assert emboss["alignment_count"] == 5
    assert emboss["alignments"][0]["sequence_ids"][0] == "ref_rec"
    msf = m.read_alignment_records(_write(tmp_path, "w.msf", data.MSF), file_format="msf", **QUIET)
    assert msf["alignments"][0]["aligned_columns"] == 99
    vulgar = m.read_alignment_records(
        _write(tmp_path, "v.exn", data.EXONERATE_VULGAR), file_format="exonerate", **QUIET
    )
    assert vulgar["alignment_count"] == 3 and vulgar["alignments"][0]["score"] == 6150.0
    chain = m.read_alignment_records(
        _write(tmp_path, "c.chain", data.CHAIN), file_format="chain", **QUIET
    )
    assert chain["alignment_count"] == 19


def test_read_alignment_records_after_writing(tmp_path):
    for fmt in ("a2m", "stockholm", "nexus", "mauve"):
        out = tmp_path / f"w.{fmt}"
        m.write_aligned_sequences(
            {"s1": "ACGT-ACGTA", "s2": "AC-TTACGAA"}, out, file_format=fmt, **QUIET
        )
        result = m.read_alignment_records(out, file_format=fmt, **QUIET)
        assert result["alignment_count"] == 1
        rows = result["alignments"][0]["aligned_rows"]
        assert sorted(row.upper() for row in rows.values()) == sorted(["ACGT-ACGTA", "AC-TTACGAA"])


def test_read_alignment_records_limits_and_errors(tmp_path):
    chain = _write(tmp_path, "c.chain", data.CHAIN)
    with pytest.raises(ConfigurationError, match="max_alignments"):
        m.read_alignment_records(chain, file_format="chain", max_alignments=2, **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_alignment_records(chain, file_format="stockholm", **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_alignment_records(chain, file_format="psl", **QUIET)


def test_convert_alignment_file_coordinate_formats(tmp_path):
    psl = _write(tmp_path, "x.psl", data.PSL)
    chain = tmp_path / "x.chain"
    result = m.convert_alignment_file(
        psl, chain, source_format="psl", target_format="chain", **QUIET
    )
    assert result["alignment_count"] == 19
    assert len(list(Align.parse(str(chain), "chain"))) == 19
    big_bed = tmp_path / "x.bb"
    m.convert_alignment_file(psl, big_bed, source_format="psl", target_format="bigbed", **QUIET)
    assert len(list(Align.parse(str(big_bed), "bigbed"))) == 19
    big_psl = tmp_path / "x.bigpsl"
    m.convert_alignment_file(psl, big_psl, source_format="psl", target_format="bigpsl", **QUIET)
    assert len(list(Align.parse(str(big_psl), "bigpsl"))) == 19
    back = tmp_path / "back.psl"
    m.convert_alignment_file(big_psl, back, source_format="bigpsl", target_format="psl", **QUIET)
    assert len(list(Align.parse(str(back), "psl"))) == 19
    for variant in ("vulgar", "cigar"):
        out = tmp_path / f"x.{variant}"
        m.convert_alignment_file(
            chain,
            out,
            source_format="chain",
            target_format="exonerate",
            exonerate_variant=variant,
            **QUIET,
        )
        assert out.read_text(encoding="utf-8").count(f"{variant}:") == 19


def test_convert_alignment_file_sequence_formats(tmp_path):
    maf = _write(tmp_path, "x.maf", data.MAF)
    for target in ("clustal", "phylip", "stockholm", "mauve"):
        out = tmp_path / f"x.{target}"
        result = m.convert_alignment_file(
            maf, out, source_format="maf", target_format=target, **QUIET
        )
        assert result["alignment_count"] == 3
        assert out.stat().st_size > 0
    big_maf = tmp_path / "x.bigmaf"
    m.convert_alignment_file(maf, big_maf, source_format="maf", target_format="bigmaf", **QUIET)
    assert len(list(Align.parse(str(big_maf), "bigmaf"))) == 3
    with pytest.raises(ConfigurationError, match="single alignment"):
        m.convert_alignment_file(
            maf, tmp_path / "x.a2m", source_format="maf", target_format="a2m", **QUIET
        )
    assert not (tmp_path / "x.a2m").exists()
    with pytest.raises(ConfigurationError):
        m.convert_alignment_file(
            maf, tmp_path / "x.psl", source_format="maf", target_format="psl", **QUIET
        )


# ------------------------------------------------------------------ MAF index


def test_search_and_splice_maf_index(tmp_path):
    maf = _write(tmp_path, "x.maf", data.MAF)
    found = m.search_maf_index(maf, "hg16.chr7", [27578828], [27578850], **QUIET)
    assert found["indexed_blocks"] == 3 and found["block_count"] == 1
    assert found["blocks"][0]["row_count"] == 5
    index = tmp_path / "x.mafindex"
    first = m.search_maf_index(maf, "hg16.chr7", [27699739], [27699745], index_path=index, **QUIET)
    assert index.exists() and first["block_count"] == 1
    again = m.search_maf_index(maf, "hg16.chr7", [27699739], [27699745], index_path=index, **QUIET)
    assert again["blocks"] == first["blocks"]
    spliced = m.splice_maf_index(
        maf, "hg16.chr7", [27578828, 27699739], [27578838, 27699745], **QUIET
    )
    assert spliced["alignment_length"] >= 16
    assert "hg16.chr7" in spliced["rows"] and "panTro1.chr6" in spliced["rows"]
    assert len({len(row) for row in spliced["rows"].values()}) == 1
    reverse = m.splice_maf_index(
        maf, "hg16.chr7", [27578828, 27699739], [27578838, 27699745], strand=-1, **QUIET
    )
    assert reverse["rows"]["hg16.chr7"] != spliced["rows"]["hg16.chr7"]
    with pytest.raises(ConfigurationError):
        m.splice_maf_index(maf, "hg16.chr7", [10], [5], **QUIET)
    with pytest.raises(ConfigurationError):
        m.splice_maf_index(maf, "hg16.chr7", [1], [2], strand=0, **QUIET)
    with pytest.raises(ConfigurationError):
        m.search_maf_index(maf, "no.such.seq", [1], [2], **QUIET)


# ------------------------------------------------------------------ search reports


def test_read_search_results_hmmer_and_infernal(tmp_path):
    cases = [
        (data.HMMER3_TAB, "hmmer3-tab", 4, 11),
        (data.HMMSCAN_DOMTAB, "hmmscan3-domtab", 4, 11),
        (data.HMMER2_TEXT, "hmmer2-text", 1, 1),
        (data.HMMER3_TEXT, "hmmer3-text", 1, 10),
        (data.INFERNAL_TAB, "infernal-tab", 3, 4),
    ]
    for text, fmt, queries, hits in cases:
        path = _write(tmp_path, f"{fmt}.out", text)
        result = m.read_search_results(path, file_format=fmt, **QUIET)
        assert (result["query_count"], result["hit_count"]) == (queries, hits), fmt
        parsed = list(SearchIO.parse(str(path), fmt))
        assert [q.id for q in parsed] == [q["id"] for q in result["queries"]]
        first_hit = result["queries"][0]["hits"][0]
        assert first_hit["id"] == parsed[0][0].id
        first_hsp = first_hit["hsps"][0]
        assert first_hsp["bitscore"] == pytest.approx(parsed[0][0][0].bitscore)
        assert first_hsp["evalue"] == pytest.approx(parsed[0][0][0].evalue)


def test_read_search_results_hsp_details_and_limits(tmp_path):
    path = _write(tmp_path, "d.out", data.HMMSCAN_DOMTAB)
    result = m.read_search_results(path, file_format="hmmscan3-domtab", **QUIET)
    hsp = result["queries"][0]["hits"][0]["hsps"][0]
    assert hsp["query_range"][0] < hsp["query_range"][1]
    with pytest.raises(ConfigurationError, match="max_queries"):
        m.read_search_results(path, file_format="hmmscan3-domtab", max_queries=1, **QUIET)
    wrong = m.read_search_results(path, file_format="hmmer3-text", **QUIET)
    assert wrong["query_count"] == 0
    with pytest.raises(ConfigurationError):
        m.read_search_results(tmp_path / "missing.out", file_format="hmmer3-tab", **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_search_results(path, file_format="blast-tab-typo", **QUIET)


def test_index_search_results_memory_and_sqlite(tmp_path):
    path = _write(tmp_path, "t.out", data.HMMER3_TEXT)
    memory = m.index_search_results(
        path, file_format="hmmer3-text", query_ids=["infile_sto"], **QUIET
    )
    assert memory["query_ids"] == ["infile_sto"] and memory["queries"][0]["hit_count"] == 10
    database = tmp_path / "idx.sqlite"
    stored = m.index_search_results(path, file_format="hmmer3-text", index_path=database, **QUIET)
    assert database.exists() and stored["query_ids"] == ["infile_sto"]
    reused = m.index_search_results(
        path, file_format="hmmer3-text", index_path=database, query_ids=["infile_sto"], **QUIET
    )
    assert reused["queries"] == memory["queries"]
    with pytest.raises(ConfigurationError, match="not found"):
        m.index_search_results(path, file_format="hmmer3-text", query_ids=["nope"], **QUIET)


def test_index_search_results_infernal_and_refused_tabular_formats(tmp_path):
    path = _write(tmp_path, "cm.tbl", data.INFERNAL_TAB)
    result = m.index_search_results(path, file_format="infernal-tab", **QUIET)
    assert result["query_count"] == 3 and result["queries"] == []
    tabular = _write(tmp_path, "t.tab", data.HMMER3_TAB)
    for fmt in ("hmmer3-tab", "hmmscan3-domtab"):
        with pytest.raises(ConfigurationError, match="cannot index"):
            m.index_search_results(tabular, file_format=fmt, **QUIET)


def test_index_search_results_blast_tab_with_comments(tmp_path):
    fields = (
        "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, "
        "q. start, q. end, s. start, s. end, evalue, bit score\n"
    )
    text = (
        f"# BLASTN 2.2.28+\n# Query: q1\n# Database: db\n{fields}# 1 hits found\n"
        "q1\ts1\t100.00\t10\t0\t0\t1\t10\t1\t10\t1e-05\t20.0\n"
        f"# BLASTN 2.2.28+\n# Query: q2\n# Database: db\n{fields}# 1 hits found\n"
        "q2\ts2\t90.00\t10\t1\t0\t1\t10\t5\t14\t1e-03\t15.0\n# BLAST processed 2 queries\n"
    )
    path = _write(tmp_path, "b.txt", text)
    result = m.index_search_results(path, file_format="blast-tab", comments=True, **QUIET)
    assert result["query_ids"] == ["q1", "q2"]
    read = m.read_search_results(path, file_format="blast-tab", comments=True, **QUIET)
    assert read["hit_count"] == 2


def test_filter_search_results(tmp_path):
    path = _write(tmp_path, "d.out", data.HMMSCAN_DOMTAB)
    everything = m.read_search_results(path, file_format="hmmscan3-domtab", **QUIET)
    strict = m.filter_search_results(path, file_format="hmmscan3-domtab", max_evalue=1e-5, **QUIET)
    assert strict["hits_before"] == everything["hit_count"]
    expected = sum(
        1
        for q in SearchIO.parse(str(path), "hmmscan3-domtab")
        for hit in q
        if any(hsp.evalue <= 1e-5 for hsp in hit.hsps)
    )
    assert strict["hits_after"] == expected < strict["hits_before"]
    scored = m.filter_search_results(path, file_format="hmmscan3-domtab", min_bitscore=100, **QUIET)
    assert all(
        hsp["bitscore"] >= 100
        for q in scored["queries"]
        for hit in q["hits"]
        for hsp in hit["hsps"]
    )
    named = m.filter_search_results(
        path, file_format="hmmscan3-domtab", hit_ids=["Globin"], **QUIET
    )
    assert named["hits_after"] == 1
    multi = m.filter_search_results(
        path, file_format="hmmscan3-domtab", min_hsps_per_hit=2, **QUIET
    )
    expected_multi = sum(
        1 for q in SearchIO.parse(str(path), "hmmscan3-domtab") for hit in q if len(hit.hsps) >= 2
    )
    assert multi["hits_after"] == expected_multi >= 2
    with pytest.raises(ConfigurationError, match="at least one"):
        m.filter_search_results(path, file_format="hmmscan3-domtab", **QUIET)


def test_sort_search_results(tmp_path):
    path = _write(tmp_path, "t.out", data.HMMER3_TEXT)
    by_score = m.sort_search_results(
        path, file_format="hmmer3-text", sort_by="bitscore", reverse=True, **QUIET
    )
    scores = [hit["bitscore"] for hit in by_score["queries"][0]["hits"]]
    assert scores == sorted(scores, reverse=True)
    by_id = m.sort_search_results(path, file_format="hmmer3-text", sort_by="id", **QUIET)
    ids = [hit["id"] for hit in by_id["queries"][0]["hits"]]
    assert ids == sorted(ids)
    by_evalue = m.sort_search_results(path, file_format="hmmer3-text", sort_by="evalue", **QUIET)
    values = [hit["evalue"] for hit in by_evalue["queries"][0]["hits"]]
    assert values == sorted(values)
    with pytest.raises(ConfigurationError):
        m.sort_search_results(path, file_format="hmmer3-text", sort_by="color", **QUIET)


# ------------------------------------------------------------------ pairwise aligners


def test_fogsaa_matches_dynamic_programming_score():
    scores = {
        "match_score": 2.0,
        "mismatch_score": -1.0,
        "open_gap_score": -2.0,
        "extend_gap_score": -1.0,
    }
    left, right = "ACGTACGTTGCA", "ACGTCGTTTGCA"
    result = m.fogsaa_global_alignment(left, right, **scores, **QUIET)
    reference = Align.PairwiseAligner(mode="global", **scores).score(left, right)
    assert result["score"] == reference == 18.0
    rows = result["alignments"][0]["aligned_rows"]
    assert (
        rows["sequence_a"].replace("-", "") == left and rows["sequence_b"].replace("-", "") == right
    )
    assert m.fogsaa_global_alignment("ACGT", "ACGT", **QUIET)["score"] == 4.0
    with pytest.raises(ConfigurationError):
        m.fogsaa_global_alignment("", "ACGT", **QUIET)
    with pytest.raises(ConfigurationError):
        m.fogsaa_global_alignment("ACGT", "ACGT", max_alignments=0, **QUIET)


def test_align_nucleotide_to_protein_plain_and_frameshift():
    plain = m.align_nucleotide_to_protein("MAIVMG", "ATGGCCATTGTAATGGGC", **QUIET)
    assert plain["score"] == 6.0
    assert plain["alignments"][0]["coordinates"] == [[0, 6], [0, 18]]
    rna = m.align_nucleotide_to_protein("MAIVMG", "AUGGCCAUUGUAAUGGGC", **QUIET)
    assert rna["score"] == 6.0
    shifted = m.align_nucleotide_to_protein(
        "MLGSNESRVCPKLSPS", "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG", **QUIET
    )
    assert shifted["score"] == 13.0
    with pytest.raises(ConfigurationError):
        m.align_nucleotide_to_protein("MAIVMG", "ATGGCC", codon_table=99, **QUIET)


# ------------------------------------------------------------------ alignment analysis


def test_sort_alignment_rows():
    assert m.sort_alignment_rows(ROWS, sort_by="id", reverse=True, **QUIET)["order"] == [
        "r3",
        "r2",
        "r1",
    ]
    by_gc = m.sort_alignment_rows(ROWS, sort_by="gc_fraction", reverse=True, **QUIET)
    assert by_gc["order"] == ["r1", "r3", "r2"]
    by_gaps = m.sort_alignment_rows(ROWS, sort_by="gap_count", reverse=True, **QUIET)
    assert set(by_gaps["order"]) == set(ROWS)
    sorted_rows = m.sort_alignment_rows(ROWS, sort_by="sequence", **QUIET)["aligned_sequences"]
    assert list(sorted_rows) == sorted(ROWS, key=lambda k: ROWS[k].replace("-", ""))
    with pytest.raises(ConfigurationError):
        m.sort_alignment_rows({"a": "ACGT"}, **QUIET)


def test_alignment_substitution_counts_and_frequencies():
    counts = m.alignment_substitution_counts(ROWS, **QUIET)
    assert counts["alphabet"] == "ACGT"
    assert sum(sum(row) for row in counts["matrix"]) == 24.0
    frequencies = m.alignment_column_frequencies(ROWS, **QUIET)
    assert frequencies["frequencies"]["-"] == [0, 0, 1, 0, 1, 0, 0, 0, 1, 0]
    assert frequencies["column_count"] == 10
    for column in range(10):
        assert sum(values[column] for values in frequencies["frequencies"].values()) == 3


def test_alignment_identity_counts():
    result = m.alignment_identity_counts(ROWS, **QUIET)["counts"]
    assert result["identities"] == 23 and result["mismatches"] == 1 and result["gaps"] == 6
    assert result["aligned"] == 24 and result["internal_gaps"] == 6
    pair = m.alignment_identity_counts({"a": "ACGT-ACGTA", "b": "AC-TTACGAA"}, **QUIET)["counts"]
    assert pair["identities"] == 7 and pair["gaps"] == 2
    scored = m.alignment_identity_counts(
        {"a": "ACGT-ACGTA", "b": "AC-TTACGAA"}, scoring_matrix="NUC.4.4", **QUIET
    )["counts"]
    assert scored["positives"] == 7 and "substitution_score" in scored
    wild = m.alignment_identity_counts({"a": "ACGTN", "b": "ACGAA"}, wildcard="N", **QUIET)[
        "counts"
    ]
    assert wild["identities"] == 3
    with pytest.raises(ConfigurationError):
        m.alignment_identity_counts(ROWS, wildcard="N", scoring_matrix="NUC.4.4", **QUIET)
    with pytest.raises(ConfigurationError):
        m.alignment_identity_counts(ROWS, scoring_matrix="NOPE", **QUIET)


def test_alignment_coordinate_maps():
    result = m.alignment_coordinate_maps(ROWS, **QUIET)
    assert result["indices"]["r1"] == [0, 1, 2, 3, -1, 4, 5, 6, 7, 8]
    assert result["inverse_indices"]["r2"] == [0, 1, 3, 4, 5, 6, 7, 8, 9]
    assert result["aligned_blocks"] is None
    pair = m.alignment_coordinate_maps({"a": "AC-GT", "b": "ACAG-"}, **QUIET)
    assert pair["aligned_blocks"] == [[[0, 2], [2, 3]], [[0, 2], [3, 4]]]


def test_parse_printed_alignment():
    result = m.parse_printed_alignment(["TAGGCATACGTG", "AACG--TACGT-", "-ACGCATACTTG"], **QUIET)
    assert result["sequences"] == ["TAGGCATACGTG", "AACGTACGT", "ACGCATACTTG"]
    assert result["coordinates"] == [
        [0, 1, 4, 6, 11, 12],
        [0, 1, 4, 4, 9, 9],
        [0, 0, 3, 5, 10, 11],
    ]
    with pytest.raises(ConfigurationError, match="same number"):
        m.parse_printed_alignment(["ACGT", "AC"], **QUIET)


def test_merge_alignments_by_reference():
    result = m.merge_alignments_by_reference(
        {"q1": {"target": "ACGT", "query": "AC-T"}, "q2": {"target": "ACGT-", "query": "A-GTT"}},
        **QUIET,
    )
    assert result["aligned_sequences"] == {"reference": "ACGT-", "q1": "AC-T-", "q2": "A-GTT"}
    with pytest.raises(ConfigurationError, match="cannot be merged"):
        m.merge_alignments_by_reference(
            {"q1": {"target": "ACGT", "query": "ACGT"}, "q2": {"target": "TTTT", "query": "TTTT"}},
            **QUIET,
        )


def test_map_alignment_through_matches_biopython():
    aligner = Align.PairwiseAligner(mode="local", open_gap_score=-1, extend_gap_score=0)
    chromosome = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
    transcript = "CCCCCCCGGGGGG"
    alignment1 = aligner.align(chromosome, transcript)[0]
    alignment2 = aligner.align(transcript, "CCCCGGGG")[0]
    mapped = alignment1.map(alignment2)
    first = {"target": str(alignment1[0]), "query": str(alignment1[1])}
    second = {"target": transcript, "query": "---CCCCGGGG--"}
    result = m.map_alignment_through(first, [second], **QUIET)
    assert result["method"] == "map"
    assert result["aligned_rows"] == {"target": str(mapped[0]), "query": str(mapped[1])}
    many = m.map_alignment_through(first, [second, second], **QUIET)
    assert many["method"] == "mapall"
    with pytest.raises(ConfigurationError):
        m.map_alignment_through(first, [], **QUIET)


# ------------------------------------------------------------------ Nexus matrix operations

DNA = {"t1": "ACGT-ACG", "t2": "ACGT-ATG", "t3": "ACGA-ACG"}


def test_find_alignment_columns():
    assert m.find_alignment_columns(DNA, **QUIET)["columns"] == [0, 1, 2, 4, 5, 7]
    assert m.find_alignment_columns(DNA, kind="gap_only", **QUIET)["columns"] == [4]
    without_t3 = m.find_alignment_columns(DNA, delete_taxa=["t3"], **QUIET)["columns"]
    assert 3 in without_t3
    excluded = m.find_alignment_columns(DNA, exclude_sites=[0], **QUIET)["columns"]
    assert 0 not in excluded
    with pytest.raises(ConfigurationError):
        m.find_alignment_columns(DNA, delete_taxa=["zz"], **QUIET)
    with pytest.raises(ConfigurationError):
        m.find_alignment_columns({"a": "ACGU1", "b": "ACGUA"}, **QUIET)


def test_site_character_states():
    result = m.site_character_states(DNA, sites=[3, 4, 6], **QUIET)["sites"]
    assert result == [
        {"site": 3, "states": ["A", "T"]},
        {"site": 4, "states": ["?"]},
        {"site": 6, "states": ["C", "T"]},
    ]
    wide = m.site_character_states(DNA, sites=[4], narrow=False, **QUIET)["sites"][0]["states"]
    assert wide == ["-", "A", "C", "G", "T"] or len(wide) >= 4
    assert len(m.site_character_states(DNA, **QUIET)["sites"]) == 8
    with pytest.raises(ConfigurationError):
        m.site_character_states(DNA, sites=[99], **QUIET)


def test_subset_alignment_and_insert_gaps():
    subset = m.subset_alignment(DNA, delete_taxa=["t2"], exclude_sites=[0, 1], **QUIET)
    assert subset["aligned_sequences"] == {"t1": "GT-ACG", "t3": "GA-ACG"}
    assert subset["row_count"] == 2 and subset["column_count"] == 6
    with pytest.raises(ConfigurationError, match="Give"):
        m.subset_alignment(DNA, **QUIET)
    with pytest.raises(ConfigurationError, match="remain"):
        m.subset_alignment(DNA, delete_taxa=["t1", "t2", "t3"], **QUIET)
    inserted = m.insert_alignment_gaps(DNA, 2, count=2, **QUIET)
    assert inserted["column_count"] == 10
    assert inserted["aligned_sequences"]["t1"] == "AC--GT-ACG"
    with pytest.raises(ConfigurationError):
        m.insert_alignment_gaps(DNA, 99, **QUIET)


def test_roundtrip_with_alignio(tmp_path):
    out = tmp_path / "x.clustal"
    m.write_aligned_sequences(ROWS, out, file_format="clustal", **QUIET)
    block = AlignIO.read(str(out), "clustal")
    assert {r.id: str(r.seq) for r in block} == ROWS


def test_public_functions_have_verbose_keyword_and_docstrings():
    import inspect

    for name in m.__all__:
        function = getattr(m, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
