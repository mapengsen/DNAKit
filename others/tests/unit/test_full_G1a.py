"""Unit tests for the Biopython and scikit-bio wrappers in ``dnakit.sequence``.

Each wrapped function is checked against a hand calculation or a direct call to
the upstream package. Files are written into ``tmp_path``. No network access is
used.
"""

from __future__ import annotations

import math
import random

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.sequence.Bio_functions import (
    ambiguous_codon_amino_acids,
    bootstrap_consensus_tree,
    build_codon_alignment,
    convert_quality_score,
    expand_ambiguous_codons,
    gc_content_by_codon_position,
    gc_skew_windows,
    load_substitution_matrix,
    local_composition_complexity,
    mcdonald_kreitman_test,
    pairwise_dn_ds_matrix,
    parse_kegg_gene_entries,
    parse_primer3_output,
    parse_primersearch_output,
    read_ace_assembly,
    read_gaf_annotations,
    read_genepop_genotypes,
    read_phd_sequence,
    read_psl_alignments,
    sequence_checksum,
)
from dnakit.sequence.skbio_functions import (
    alignment_distance_matrix,
    alignment_score,
    multiple_nucleotide_alignment,
    nucleotide_distance,
    nucleotide_transition_matrix,
    sniff_sequence_format,
)

S1 = "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCA"
S2 = "ATGGCAGCCAAATTCGGGCTGAAAGATCGCAAGCCTTTCGCA"
S3 = "ATGGCTGCCAAATTTGGGCTGAAAGATCGTAAGCCTTTCGCT"


def _write(tmp_path, name: str, text: str):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --- Biopython ---------------------------------------------------------------


def test_gc_content_by_codon_position_hand_counts_and_upstream():
    pytest.importorskip("Bio.SeqUtils")
    from Bio.SeqUtils import GC123

    result = gc_content_by_codon_position("ATGCGTACGTTAGC", verbose=False)
    # 7 of 14 bases are G or C. Codon position 1 has 2 of 5 (40%), position 2 has 3 of 5
    # (60%) and position 3 has 2 of 4 (50%).
    assert result["gc_percent"] == {"total": 50.0, "first": 40.0, "second": 60.0, "third": 50.0}
    assert tuple(result["gc_percent"].values()) == tuple(float(v) for v in GC123("ATGCGTACGTTAGC"))


def test_gc_skew_windows_hand_calculation():
    pytest.importorskip("Bio.SeqUtils")
    result = gc_skew_windows("GGGGCCCCAAAATTTTGGGG", window=5, verbose=False)
    assert [item["skew"] for item in result["windows"]] == [0.6, -1.0, 0.0, 1.0]
    assert [(item["start"], item["end"]) for item in result["windows"]] == [
        (1, 5),
        (6, 10),
        (11, 15),
        (16, 20),
    ]


def test_local_composition_complexity_uniform_window_is_one():
    pytest.importorskip("Bio.SeqUtils.lcc")
    from Bio.SeqUtils.lcc import lcc_mult

    result = local_composition_complexity("ACGTACGT", window=4, verbose=False)
    # Every window of 4 contains each base once, so the complexity is maximal (1.0).
    assert result["values"] == pytest.approx([1.0] * 5)
    assert result["values"] == [float(v) for v in lcc_mult("ACGTACGT", 4)]
    assert result["value_count"] == 5


def test_ambiguous_codon_amino_acids_and_stop_only_case():
    pytest.importorskip("Bio.Data.CodonTable")
    gcn = ambiguous_codon_amino_acids("GCN", verbose=False)
    assert gcn["amino_acids"] == ["A"]
    assert gcn["stop_only"] is False
    tty = ambiguous_codon_amino_acids("TTY", verbose=False)
    assert tty["amino_acids"] == ["F"]
    stop = ambiguous_codon_amino_acids("TAR", verbose=False)
    assert stop["amino_acids"] == []
    assert stop["stop_only"] is True
    with pytest.raises(ConfigurationError):
        ambiguous_codon_amino_acids("AT", verbose=False)


def test_expand_ambiguous_codons_hand_rule():
    pytest.importorskip("Bio.Data.CodonTable")
    result = expand_ambiguous_codons(["TGA", "TAA", "TAG"], verbose=False)
    # TAR covers TAA and TAG, and TRA covers TAA and TGA: both are fully present.
    # TRR would also need TGG, which is absent, so it is not added.
    assert result["added_codons"] == ["TAR", "TRA"]
    assert result["expanded_codons"] == ["TAA", "TAG", "TAR", "TGA", "TRA"]


def test_sequence_checksum_matches_upstream_for_each_method():
    pytest.importorskip("Bio.SeqUtils.CheckSum")
    from Bio.SeqUtils import CheckSum

    expected = {
        "seguid": CheckSum.seguid("ATGCGT"),
        "crc32": CheckSum.crc32("ATGCGT"),
        "crc64": CheckSum.crc64("ATGCGT"),
        "gcg": CheckSum.gcg("ATGCGT"),
    }
    for method, value in expected.items():
        assert sequence_checksum("atg cgt", method=method, verbose=False)["value"] == value
    with pytest.raises(ConfigurationError):
        sequence_checksum("ATGCGT", method="md5", verbose=False)


def test_mcdonald_kreitman_test_hand_table_and_upstream():
    pytest.importorskip("Bio.Align.analysis")
    seqs = {
        "a1": "GCTGAAAAATTT",
        "a2": "GCCGAAAAATTA",
        "b1": "GCTGAGCAATTT",
        "b2": "GCTGAGCAATTT",
    }
    species = {"a1": "A", "a2": "A", "b1": "B", "b2": "B"}
    result = mcdonald_kreitman_test(seqs, species, verbose=False)
    # Each of the four cells holds one change (syn/nonsyn x fixed/polymorphic), so
    # the G statistic is 0 and the p-value is 1.
    assert result["p_value"] == pytest.approx(1.0)
    assert result["species"] == ["A", "B"]
    assert result["codon_count"] == 4


def test_pairwise_dn_ds_matrix_matches_direct_pair_and_zero_dn():
    pytest.importorskip("Bio.Align.analysis")
    seqs = {"s1": S1, "s2": S2, "s3": S3}
    result = pairwise_dn_ds_matrix(seqs, verbose=False)
    # The s1-s2 differences are all synonymous (GCT/GCA, TTT/TTC, CGT/CGC), so dN is 0.
    assert result["ids"] == ["s1", "s2", "s3"]
    assert result["dn"][0][1] == pytest.approx(0.0)
    assert result["saturated_pair_count"] == 0
    assert result["ds"][1][0] == result["ds"][0][1]
    assert all(result["ds"][i][i] == 0.0 for i in range(3))

    import numpy as np
    from Bio.Align import Alignment
    from Bio.Align.analysis import calculate_dn_ds
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    pair = Alignment(
        [SeqRecord(Seq(S1), id="s1"), SeqRecord(Seq(S2), id="s2")],
        np.array([[0, len(S1)], [0, len(S1)]]),
    )
    _, direct_ds = calculate_dn_ds(pair)
    assert result["ds"][0][1] == pytest.approx(direct_ds)


def test_build_codon_alignment_threads_codons_through_gapped_protein():
    pytest.importorskip("Bio.codonalign")
    protein = {"s1": "MAA-KF", "s2": "MAGKF-", "s3": "MAAKY-"}
    nucleotides = {
        "s1": "ATGGCTGCCAAATTT",
        "s2": "ATGGCAGGCAAATTT",
        "s3": "ATGGCTGCCAAATAT",
    }
    result = build_codon_alignment(protein, nucleotides, verbose=False)
    assert result["alignment_length"] == 18
    # Protein gaps become three-base codon gaps; the ungapped bases stay in order.
    assert result["sequences"]["s1"] == "ATGGCTGCC---AAATTT"
    assert result["sequences"]["s2"] == "ATGGCAGGCAAATTT---"
    assert result["sequences"]["s3"] == "ATGGCTGCCAAATAT---"
    assert {k: v.replace("-", "") for k, v in result["sequences"].items()} == nucleotides


def test_parse_primersearch_output_reads_amplicon(tmp_path):
    pytest.importorskip("Bio.Emboss.PrimerSearch")
    path = _write(
        tmp_path,
        "primersearch.txt",
        "Primer name P1\nAmplimer 1\n\tSequence: ATGCGATTTGCAAGCTTAGC\n"
        "\tAmplimer length: 42 bp\n\tForward primer: 0 mismatches\n",
    )
    result = parse_primersearch_output(path, verbose=False)
    assert result["primer_count"] == 1
    amplicon = result["primers"]["P1"][0]
    assert amplicon["sequence"] == "ATGCGATTTGCAAGCTTAGC"
    assert amplicon["length"] == 42
    assert amplicon["details"] == ["Forward primer: 0 mismatches"]


def test_read_ace_assembly_reads_contig_and_placements(tmp_path):
    pytest.importorskip("Bio.Sequencing.Ace")
    path = _write(
        tmp_path,
        "contig.ace",
        "AS 1 2\n\nCO Contig1 12 2 4 U\nACGTACGTACGT\n\nBQ\n"
        + " ".join(["30"] * 12)
        + "\n\nAF read1 U 1\nAF read2 U 2\n\n"
        "RD read1 12 0 0\nACGTACGTACGT\n\nQA 1 12 1 12\n"
        "DS CHROMAT_FILE read1 PHD_FILE read1.phd.1 TIME Thu Jan 1 00:00:00 2024\n\n"
        "RD read2 12 0 0\nACGTACGTACGT\n\nQA 1 12 1 12\n"
        "DS CHROMAT_FILE read2 PHD_FILE read2.phd.1 TIME Thu Jan 1 00:00:00 2024\n\n",
    )
    result = read_ace_assembly(path, verbose=False)
    assert result["contig_count"] == 1
    assert result["read_count"] == 2
    contig = result["contigs"][0]
    assert contig["name"] == "Contig1"
    assert contig["sequence"] == "ACGTACGTACGT"
    assert contig["placements"] == [
        {"read": "read1", "strand": "U", "padded_start": 1},
        {"read": "read2", "strand": "U", "padded_start": 2},
    ]


def test_read_phd_sequence_bases_qualities_and_peaks(tmp_path):
    pytest.importorskip("Bio.Sequencing.Phd")
    path = _write(
        tmp_path,
        "read.phd",
        "BEGIN_SEQUENCE read1\nBEGIN_COMMENT\nCHROMAT_FILE: read1.ab1\n"
        "PHRED_VERSION: 0.980528.c\nEND_COMMENT\nBEGIN_DNA\n"
        "A 10 1\nC 20 4\nG 30 2\nT 40 3\nEND_DNA\nEND_SEQUENCE\n",
    )
    result = read_phd_sequence(path, verbose=False)
    assert result["file_name"] == "read1"
    assert result["sequence"] == "ACGT"
    assert result["qualities"] == [10, 20, 30, 40]
    assert result["peak_locations"] == [1, 4, 2, 3]
    assert result["comments"]["chromat_file"] == "read1.ab1"


def test_read_genepop_genotypes_missing_data_as_null(tmp_path):
    pytest.importorskip("Bio.PopGen.GenePop")
    path = _write(
        tmp_path,
        "pop.genepop",
        "Example\nLocus1\nLocus2\nPop\nind1 , 0101 0203\nind2 , 0102 0000\n",
    )
    result = read_genepop_genotypes(path, verbose=False)
    assert result["loci"] == ["Locus1", "Locus2"]
    assert result["population_count"] == 1
    assert result["individual_count"] == 2
    individuals = result["populations"][0]["individuals"]
    assert individuals[0] == {"name": "ind1", "genotypes": [[1, 1], [2, 3]]}
    assert individuals[1]["genotypes"] == [[1, 2], [None, None]]


def test_read_gaf_annotations_counts_terms_and_aspects(tmp_path):
    pytest.importorskip("Bio.UniProt.GOA")
    first = ["UniProtKB", "P12345", "GENEA", "", "GO:0003700", "PMID:1234567", "IDA"]
    first += ["", "F", "TF A", "GENEA|TFA", "protein", "taxon:9606", "20240101", "UniProt"]
    second = ["UniProtKB", "P67890", "GENEB", "involved_in", "GO:0006355", "PMID:7654321"]
    second += ["IEA", "", "P", "Protein B", "GENEB", "protein", "taxon:9606", "20240102"]
    second += ["UniProt"]
    text = "!gaf-version: 2.2\n" + "\t".join(first) + "\n" + "\t".join(second) + "\n"
    path = _write(tmp_path, "annotations.gaf", text)
    result = read_gaf_annotations(path, verbose=False)
    assert result["record_count"] == 2
    assert result["go_term_counts"] == {"GO:0003700": 1, "GO:0006355": 1}
    assert result["aspect_counts"] == {"F": 1, "P": 1}
    assert result["records"][0]["synonyms"] == ["GENEA", "TFA"]
    assert result["records"][1]["qualifiers"] == ["involved_in"]


def test_parse_kegg_gene_entries_fields(tmp_path):
    pytest.importorskip("Bio.KEGG.Gene")
    path = _write(
        tmp_path,
        "gene.txt",
        "ENTRY       hsa:7157                    CDS       T01001\n"
        "SYMBOL      TP53\n"
        "NAME        (RefSeq) tumor protein p53\n"
        "ORTHOLOGY   K04442  tumor protein p53\n"
        "ORGANISM    hsa  Homo sapiens (human)\n"
        "POSITION    17:7661779..7687538\n"
        "DBLINKS     NCBI-GeneID: 7157\n"
        "///\n",
    )
    result = parse_kegg_gene_entries(path, verbose=False)
    entry = result["entries"][0]
    assert result["entry_count"] == 1
    assert entry["entry"] == "hsa:7157"
    assert entry["orthology"] == [{"id": "K04442", "name": "tumor protein p53"}]
    assert entry["organism"] == {"code": "hsa", "name": "Homo sapiens (human)"}
    assert entry["position"] == "17:7661779..7687538"
    assert entry["dblinks"] == [{"database": "NCBI-GeneID", "identifiers": ["7157"]}]


def test_load_substitution_matrix_nuc44_values_and_upstream():
    pytest.importorskip("Bio.Align.substitution_matrices")
    from Bio.Align import substitution_matrices

    result = load_substitution_matrix("NUC.4.4", verbose=False)
    index = {symbol: position for position, symbol in enumerate(result["alphabet"])}
    matrix = result["matrix"]
    # NUC.4.4 scores a match as +5 and a mismatch as -4.
    assert matrix[index["A"]][index["A"]] == 5.0
    assert matrix[index["A"]][index["C"]] == -4.0
    assert result["symmetric"] is True
    upstream = substitution_matrices.load("NUC.4.4")
    assert matrix[index["G"]][index["T"]] == float(upstream[("G", "T")])
    with pytest.raises(ConfigurationError):
        load_substitution_matrix("NOT-A-MATRIX", verbose=False)


def test_convert_quality_score_closed_forms():
    pytest.importorskip("Bio.SeqIO.QualityIO")
    forward = convert_quality_score(10.0, verbose=False)
    # Phred = 10 log10(10^(Q/10) + 1), so Q = 10 gives 10 log10(11).
    assert forward["output"] == pytest.approx(10 * math.log10(11))
    back = convert_quality_score(10.0, direction="phred_to_solexa", verbose=False)
    # Solexa = 10 log10(10^(Q/10) - 1), so Q = 10 gives 10 log10(9).
    assert back["output"] == pytest.approx(10 * math.log10(9))
    with pytest.raises(ConfigurationError):
        convert_quality_score(-1.0, direction="phred_to_solexa", verbose=False)


def test_read_psl_alignments_block_coordinates(tmp_path):
    pytest.importorskip("Bio.Align")
    path = _write(
        tmp_path,
        "alignments.psl",
        "psLayout version 3\n\n"
        "match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\t"
        "Q name\tQ size\tQ start\tQ end\t"
        "T name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\n"
        "     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \t"
        "\t\t\t\t\t\t\t\t\n"
        "---------------------------------------------------------------------------------------------------------------\n"
        "20\t0\t0\t0\t0\t0\t0\t0\t+\tread1\t20\t0\t20\tchr1\t1000\t100\t120\t1\t20,\t0,\t100,\n",
    )
    result = read_psl_alignments(path, verbose=False)
    assert result["alignment_count"] == 1
    alignment = result["alignments"][0]
    assert (alignment["target_id"], alignment["query_id"]) == ("chr1", "read1")
    assert (alignment["target_start"], alignment["target_end"]) == (100, 120)
    assert (alignment["query_start"], alignment["query_end"]) == (0, 20)


def test_bootstrap_consensus_tree_is_seed_deterministic_and_matches_upstream():
    pytest.importorskip("Bio.Phylo.Consensus")
    from Bio.Align import MultipleSeqAlignment
    from Bio.Phylo.Consensus import bootstrap_consensus, majority_consensus
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rows = {
        "s1": "ACGTACGTAAGGCCTTACGT",
        "s2": "ACGTACCTAAGGCATTACGA",
        "s3": "ACGTTCGTAAGGCCTAACGT",
        "s4": "ACGTTCGAAAGGCCTAACGA",
    }
    result = bootstrap_consensus_tree(rows, replicates=5, seed=1, verbose=False)
    again = bootstrap_consensus_tree(rows, replicates=5, seed=1, verbose=False)
    assert result["newick"] == again["newick"]
    assert result["newick"].endswith(";")
    for identifier in rows:
        assert identifier in result["newick"]

    alignment = MultipleSeqAlignment([SeqRecord(Seq(v), id=k) for k, v in rows.items()])
    constructor = DistanceTreeConstructor(DistanceCalculator("identity"), "nj")
    random.seed(1)
    direct = bootstrap_consensus(alignment, 5, constructor, majority_consensus)
    assert result["newick"] == direct.format("newick").strip()


def test_parse_primer3_output_pair_and_positions(tmp_path):
    pytest.importorskip("Bio.Emboss.Primer3")
    path = _write(
        tmp_path,
        "primer3.txt",
        "# PRIMER3 test\n"
        "#                      Start  Len   Tm     GC%   Sequence\n"
        "     PRODUCT SIZE: 60\n"
        "     FORWARD PRIMER    1  20  60.00  50.00  ACGTACGTACGTACGTACGT\n"
        "     REVERSE PRIMER   60  20  60.00  50.00  TGCATGCATGCATGCATGCA\n",
    )
    result = parse_primer3_output(path, verbose=False)
    pair = result["primer_pairs"][0]
    assert result["pair_count"] == 1
    assert pair["product_size"] == 60
    assert pair["forward"]["sequence"] == "ACGTACGTACGTACGTACGT"
    # Positions and melting temperatures are copied from the report as written.
    assert pair["reverse"]["start"] == 60
    assert pair["reverse"]["tm"] == pytest.approx(60.0)
    assert pair["internal"] is None


# --- scikit-bio ----------------------------------------------------------------


def test_nucleotide_distance_p_distance_and_jc69_closed_form():
    pytest.importorskip("skbio")
    left = "ACGTACGTAAGGCCTTACGT"
    right = "ACGTACCTAAGGCATTACGA"
    p_distance = nucleotide_distance(left, right, model="pdist", verbose=False)
    # Three of 20 sites differ.
    assert p_distance["distance"] == pytest.approx(0.15)
    jc = nucleotide_distance(left, right, model="jc69", verbose=False)
    # JC69: d = -3/4 ln(1 - 4/3 p).
    assert jc["distance"] == pytest.approx(-0.75 * math.log(1 - 4 / 3 * 0.15))


def test_nucleotide_distance_saturated_is_undefined_not_a_number():
    pytest.importorskip("skbio")
    result = nucleotide_distance("AAAAAAAAAA", "CCCCCCCCCC", model="jc69", verbose=False)
    assert result["status"] == "undefined"
    assert result["distance"] is None
    with pytest.raises(ConfigurationError):
        nucleotide_distance("ACGT", "ACG", verbose=False)


def test_nucleotide_transition_matrix_jc69_closed_form():
    pytest.importorskip("skbio")
    from skbio.sequence import transition

    result = nucleotide_transition_matrix(0.2, verbose=False)
    same = 0.25 + 0.75 * math.exp(-4 * 0.2 / 3)
    other = 0.25 - 0.25 * math.exp(-4 * 0.2 / 3)
    assert result["matrix"][0][0] == pytest.approx(same)
    assert result["matrix"][0][1] == pytest.approx(other)
    direct = transition.jc69(0.2, dtype="float64").data.tolist()
    flat_result = [value for row in result["matrix"] for value in row]
    flat_direct = [value for row in direct for value in row]
    assert flat_result == pytest.approx(flat_direct)
    assert all(sum(row) == pytest.approx(1.0) for row in result["matrix"])


def test_alignment_distance_matrix_hamming_hand_counts():
    pytest.importorskip("skbio")
    rows = {"x": "ACGTACGTAA", "y": "ACGTTCGTAC", "z": "ACGAACGTAA"}
    result = alignment_distance_matrix(rows, metric="hamming", verbose=False)
    # x-y differ at sites 5 and 10 (2/10); x-z at site 4 (1/10); y-z at sites 4, 5, 10 (3/10).
    assert result["ids"] == ["x", "y", "z"]
    assert result["matrix"][0][1] == pytest.approx(0.2)
    assert result["matrix"][0][2] == pytest.approx(0.1)
    assert result["matrix"][1][2] == pytest.approx(0.3)


def test_multiple_nucleotide_alignment_recovers_inputs_and_matches_upstream():
    pytest.importorskip("skbio")
    import skbio.alignment as align
    from skbio import DNA

    rows = {"s1": "ACGTACGTAAGGCC", "s2": "ACGTACCTAAGGCA", "s3": "ACGTTCGTAAGGCCTA"}
    result = multiple_nucleotide_alignment(rows, verbose=False)
    for identifier, original in rows.items():
        assert result["alignment"][identifier].replace("-", "") == original
    lengths = {len(text) for text in result["alignment"].values()}
    assert lengths == {result["alignment_length"]}
    objects = [DNA(text) for text in rows.values()]
    upstream = align.multi_align_nucl(objects, sub_score=(1.0, -1.0), gap_cost=2.0, free_ends=True)
    assert list(result["alignment"].values()) == upstream.path.to_aligned(objects)


def test_alignment_score_hand_values():
    pytest.importorskip("skbio")
    identical = alignment_score({"x": "ACGTA", "y": "ACGTA"}, verbose=False)
    assert identical["score"] == pytest.approx(5.0)
    gapped = alignment_score({"x": "AC-TA", "y": "ACGTA"}, verbose=False)
    # Four matches score +4; one gap costs 2, giving 2.
    assert gapped["score"] == pytest.approx(2.0)


def test_sniff_sequence_format_fasta_and_fastq(tmp_path):
    pytest.importorskip("skbio")
    fasta = _write(tmp_path, "reads.fa", ">a\nACGT\n>b\nACGA\n")
    fastq = _write(tmp_path, "reads.fq", "@r1\nACGT\n+\nIIII\n")
    assert sniff_sequence_format(fasta, verbose=False)["format"] == "fasta"
    assert sniff_sequence_format(fastq, verbose=False)["format"] == "fastq"
    unknown = _write(tmp_path, "notes.txt", "not a sequence file at all\n")
    with pytest.raises(ConfigurationError):
        sniff_sequence_format(unknown, verbose=False)
