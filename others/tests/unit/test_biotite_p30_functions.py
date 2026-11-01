"""Tests for the Biotite sequence/alignment and nucleic-acid structure wrappers (pack P30).

Results are compared with direct Biotite calls on the same inputs. Nothing uses the network.
The two small PDB files in others/tests/data/biotite_p30 are model 1 of PDB 2KOC (RNA hairpin with
hydrogens) and PDB 1BNA (Dickerson DNA dodecamer).
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("biotite")
np = pytest.importorskip("numpy")

import biotite.sequence as bseq
import biotite.sequence.align as balign
import biotite.structure as bstruct

from dnakit.alignment import biotite_p30_functions as m
from dnakit.exceptions import ConfigurationError
from dnakit.secondary_structure import biotite_p30_structure_functions as s

QUIET = {"verbose": False}
DATA = Path(__file__).resolve().parents[1] / "data" / "biotite_p30"
ALIGNED = ["ACGTACGT", "AC-TACGT", "ACGTTCGT"]
SEQ1, SEQ2 = "ACGTACGTTAGC", "TTACGTACGTAAGC"


def _alignment(rows):
    return balign.Alignment.from_strings(rows, bseq.NucleotideSequence)


def test_start_codon_matches_biotite_and_table_11():
    result = m.codon_table_is_start_codon(["ATG", "GTG", "AAA"], table=11, **QUIET)
    table = bseq.CodonTable.load(11)
    assert result["is_start"] == [True, True, False]
    assert result["is_start"][0] == bool(table.is_start_codon(np.array([[0, 3, 2]]))[0])
    assert m.codon_table_is_start_codon(["GTG"], **QUIET)["is_start"] == [False]
    with pytest.raises(ConfigurationError):
        m.codon_table_is_start_codon(["AT"], **QUIET)


def test_profile_counts_probability_score():
    profile = bseq.SequenceProfile.from_alignment(_alignment(ALIGNED))
    counts = m.sequence_profile_from_alignment(ALIGNED, **QUIET)
    assert counts["counts"] == profile.symbols.tolist()
    assert counts["gap_counts"] == profile.gaps.tolist()
    seq = bseq.NucleotideSequence("ACGTACGT")
    assert m.sequence_profile_probability(ALIGNED, "ACGTACGT", **QUIET)["probability"] == pytest.approx(
        profile.sequence_probability(seq)
    )
    assert m.sequence_profile_score(ALIGNED, "ACGTACGT", pseudocount=1, **QUIET)["score"] == pytest.approx(
        profile.sequence_score(seq, pseudocount=1)
    )
    with pytest.raises(ConfigurationError):
        m.sequence_profile_probability(ALIGNED, "ACG", **QUIET)


def test_alignment_from_strings_trace():
    result = m.alignment_from_strings(ALIGNED, **QUIET)
    assert result["trace"][2] == [2, -1, 2]
    assert result["sequences"][1] == "ACTACGT"
    with pytest.raises(ConfigurationError):
        m.alignment_from_strings(["AC", "ACG"], **QUIET)


def test_kmer_index_and_match():
    refs = ["ACGTACGTTAGC", "TTACGTA"]
    table = balign.KmerTable.from_sequences(
        4, [bseq.NucleotideSequence(x) for x in refs]
    )
    code = int(table.kmer_alphabet.encode("ACGT"))
    index = m.kmer_index(refs, k=4, **QUIET)
    row = next(item for item in index["kmers"] if item["kmer"] == "ACGT")
    assert row["positions"] == table[code].tolist()
    bucket = m.kmer_index(refs, k=4, bucketed=True, **QUIET)
    assert bucket["table"] == "bucket"
    matches = m.kmer_match(refs, "CGTACG", k=4, **QUIET)
    expected = table.match(bseq.NucleotideSequence("CGTACG")).tolist()
    assert matches["matches"] == expected
    similar = m.kmer_match(refs, "CGTACG", k=4, similarity_threshold=3, **QUIET)
    assert similar["match_count"] >= matches["match_count"]


def test_syncmer_mincode_and_similar_kmers():
    text = "ACGTACGTTAGCATCG"
    seq = bseq.NucleotideSequence(text)
    selector = balign.CachedSyncmerSelector(seq.alphabet, 4, 2)
    positions, _ = selector.select(seq)
    result = m.syncmer_selection(text, k=4, s=2, **QUIET)
    assert result["positions"] == positions.tolist()
    assert all(len(kmer) == 4 for kmer in result["kmers"])
    for position, kmer in zip(result["positions"], result["kmers"], strict=True):
        assert text[position : position + 4] == kmer
    mincode = m.mincode_selection(text, k=4, compression=2.0, **QUIET)
    assert 0 < len(mincode["positions"]) < len(text) - 3
    similar = m.similar_kmers("ACGT", threshold=3, **QUIET)
    assert "ACGT" in similar["similar_kmers"] and "ACGA" in similar["similar_kmers"]
    assert similar["count"] == 13
    with pytest.raises(ConfigurationError):
        m.syncmer_selection(text, k=4, s=4, **QUIET)


def test_evalue_functions():
    params = m.alignment_evalue_parameters(sample_length=100, sample_size=100, seed=0, **QUIET)
    again = m.alignment_evalue_parameters(sample_length=100, sample_size=100, seed=0, **QUIET)
    assert params["lambda"] == again["lambda"] > 0
    estimator = balign.EValueEstimator(1.0, 0.1)
    result = m.alignment_log_evalue(30, 100, 100, lam=1.0, k=0.1, **QUIET)
    assert result["log_evalue"] == pytest.approx(float(estimator.log_evalue(30, 100, 100)))
    with pytest.raises(ConfigurationError):
        m.alignment_log_evalue(30, 100, 100, lam=1.0, **QUIET)


def test_seed_banded_local_gapped_ungapped():
    s1, s2 = bseq.NucleotideSequence(SEQ1), bseq.NucleotideSequence(SEQ2)
    ext = balign.SeedExtension((1, -1), 5).align(s1, s2, (4, 6))
    result = m.seed_extension_ungapped(SEQ1, SEQ2, seed=[4, 6], threshold=5, **QUIET)
    assert result["aligned"] == [str(x) for x in ext.get_gapped_sequences()]
    assert result["score"] == ext.score
    banded = balign.align_banded(s1, s2, (-3, 3), local=True)
    got = m.align_banded_pairwise(SEQ1, SEQ2, band=[-3, 3], local=True, **QUIET)
    assert got["alignments"][0]["aligned"] == [str(x) for x in banded[0].get_gapped_sequences()]
    local = balign.align_local_gapped(s1, s2, (4, 6), 5)
    got = m.align_local_gapped_from_seed(SEQ1, SEQ2, seed=[4, 6], threshold=5, **QUIET)
    assert got["alignments"][0]["score"] == local[0].score
    ungapped = m.align_ungapped_pairwise("ACGT", "ACGA", **QUIET)
    assert ungapped["score"] == 2
    with pytest.raises(ConfigurationError):
        m.align_ungapped_pairwise("ACGT", "ACG", **QUIET)
    with pytest.raises(ConfigurationError):
        m.align_banded_pairwise(SEQ1, SEQ2, band=[3, -3], **QUIET)


def test_multiple_alignment_keeps_input_order():
    seqs = ["ACGTACGTTAGCATCG", "TTACGTACGTAAGC", "ACGTACGTTAGCATCG"]
    result = m.align_multiple_progressive(seqs, **QUIET)
    assert [x.replace("-", "") for x in result["aligned"]] == seqs
    assert len({len(x) for x in result["aligned"]}) == 1


def test_alignment_edit_and_identity_functions():
    trimmed = m.alignment_remove_terminal_gaps(["--ACGT", "TTACG-"], **QUIET)
    assert trimmed["aligned"] == ["ACG", "ACG"]
    assert (trimmed["first_kept_column"], trimmed["last_kept_column_exclusive"]) == (2, 5)
    assert m.alignment_remove_gaps(ALIGNED, **QUIET)["aligned"][0] == "ACTACGT"
    alignment = _alignment(ALIGNED)
    pair = m.alignment_pairwise_identity(ALIGNED, mode="all", **QUIET)
    assert pair["identity"] == balign.get_pairwise_sequence_identity(alignment, mode="all").tolist()
    assert m.alignment_sequence_identity(ALIGNED, **QUIET)["identity"] == pytest.approx(
        balign.get_sequence_identity(alignment)
    )
    assert m.alignment_score(ALIGNED[:2], **QUIET)["score"] == balign.score(_alignment(ALIGNED[:2]))
    assert m.alignment_score(ALIGNED[:2], terminal_penalty=False, gap_penalty=-5, **QUIET)["score"] == balign.score(
        _alignment(ALIGNED[:2]), (1, -1), -5, False
    )
    with pytest.raises(ConfigurationError):
        m.alignment_pairwise_identity(ALIGNED, mode="bad", **QUIET)


@pytest.mark.parametrize("coloring", ["type", "similarity"])
def test_plot_alignment_modes(tmp_path, coloring):
    pytest.importorskip("matplotlib")
    result = m.plot_sequence_alignment(ALIGNED, str(tmp_path / f"{coloring}.png"), coloring=coloring, **QUIET)
    assert Path(result["image"]).stat().st_size > 1000
    assert result["columns"] == 8


def test_plot_alignment_values_and_errors(tmp_path):
    pytest.importorskip("matplotlib")
    values = [[0.5] * 3 for _ in range(8)]
    result = m.plot_sequence_alignment(ALIGNED, str(tmp_path / "v.svg"), coloring="values", values=values, **QUIET)
    assert Path(result["image"]).is_file()
    with pytest.raises(ConfigurationError):
        m.plot_sequence_alignment(ALIGNED, str(tmp_path / "w.png"), coloring="values", **QUIET)
    with pytest.raises(ConfigurationError):
        m.plot_sequence_alignment(ALIGNED, str(tmp_path / "x.txt"), **QUIET)


def test_write_a3m(tmp_path):
    out = tmp_path / "a.a3m"
    result = m.write_a3m_alignment(ALIGNED, str(out), seq_type="dna", **QUIET)
    text = out.read_text()
    assert text.startswith(">query\nACGTACGT\n>target_1\nAC-TACGT")
    assert result["targets"] == 2
    with pytest.raises(ConfigurationError):
        m.write_a3m_alignment(ALIGNED, str(out), seq_type="dna", **QUIET)


def test_structure_base_pairs_edges_stacking():
    rna = str(DATA / "rna_hairpin.pdb")
    pairs = s.nucleic_acid_base_pairs(rna, **QUIET)
    assert pairs["pair_count"] == 5
    assert pairs["base_pairs"][0]["base_1"]["res_id"] == 1
    assert pairs["base_pairs"][0]["base_2"]["res_id"] == 14
    detail = s.nucleic_acid_base_pairs(rna, annotate_edges=True, annotate_glycosidic_bond=True, **QUIET)
    assert detail["base_pairs"][0]["edge_1"] == "watson_crick"
    assert detail["base_pairs"][0]["glycosidic_bond"] == "cis"
    stacking = s.nucleic_acid_base_stacking(rna, **QUIET)
    assert stacking["stack_count"] == 6
    dna = str(DATA / "dna_dickerson.pdb")
    assert s.nucleic_acid_base_pairs(dna, **QUIET)["pair_count"] == 12
    with pytest.raises(ConfigurationError):
        s.nucleic_acid_base_pairs(dna, annotate_edges=True, **QUIET)
    with pytest.raises(ConfigurationError):
        s.nucleic_acid_base_pairs(str(DATA / "missing.pdb"), **QUIET)


def test_structure_dot_bracket_and_pseudoknots():
    rna = str(DATA / "rna_hairpin.pdb")
    assert s.nucleic_acid_dot_bracket(rna, **QUIET)["dot_bracket"] == ["(((((....)))))"]
    duplex = s.nucleic_acid_dot_bracket(str(DATA / "dna_dickerson.pdb"), **QUIET)
    assert duplex["dot_bracket"] == ["((((((((((((" + "))))))))))))"]
    pairs = [[0, 10], [1, 9], [2, 8], [4, 14], [5, 13]]
    assert s.base_pairs_to_dot_bracket(pairs, 15, **QUIET)["dot_bracket"] == ["(((.[[..)))..]]"]
    assert s.pseudoknot_order(pairs, **QUIET)["orders"] == [[0, 0, 0, 1, 1]]
    assert bstruct.pseudoknots(np.array(pairs)).tolist() == [[0, 0, 0, 1, 1]]
    with pytest.raises(ConfigurationError):
        s.base_pairs_to_dot_bracket([[0, 1], [1, 2]], 5, **QUIET)
    with pytest.raises(ConfigurationError):
        s.base_pairs_to_dot_bracket([[0, 9]], 5, **QUIET)


def test_plot_secondary_structure_needs_rnaplot(tmp_path):
    import shutil

    if shutil.which("RNAplot") is None:
        with pytest.raises(Exception, match="RNAplot"):
            s.plot_nucleotide_secondary_structure("GGGAAACCC", "(((...)))", str(tmp_path / "p.png"), **QUIET)
    else:
        result = s.plot_nucleotide_secondary_structure(
            "GGGAAACCC", "(((...)))", str(tmp_path / "p.png"), **QUIET
        )
        assert Path(result["image"]).stat().st_size > 1000
        assert result["base_pairs"] == [[0, 8], [1, 7], [2, 6]]
    with pytest.raises(ConfigurationError):
        s.plot_nucleotide_secondary_structure("GGG", "((", str(tmp_path / "q.png"), **QUIET)
