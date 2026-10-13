"""Review group G18o: Biopython, scikit-bio, cutadapt and umi_tools wrappers.

Each wrapped function is checked against the real upstream call or a hand calculation. Tests
build their own small inputs in ``tmp_path``; nothing needs the network.
"""

from __future__ import annotations

import inspect
import io
from pathlib import Path

import pytest

from dnakit.alignment import mappy_review_functions as mappy_review
from dnakit.alignment import parasail_review_functions as parasail_review
from dnakit.alignment import pyhmmer_review_functions as pyhmmer_review
from dnakit.exceptions import ConfigurationError
from dnakit.sequence import Bio_review_functions as bio_review
from dnakit.sequence import pyfastx_review_functions as pyfastx_review
from dnakit.sequence import skbio_review_functions as skbio_review
from dnakit.sequencing_qc import cutadapt_review_functions as cutadapt_review
from dnakit.sequencing_qc import umi_tools_review_functions as umi_review

WRAPPERS = [
    bio_review.adam_consensus_tree,
    bio_review.bootstrap_branch_support,
    bio_review.codon_back_table,
    bio_review.complement_rna_sequence,
    bio_review.convert_search_results,
    bio_review.melting_temperature_gc,
    bio_review.read_fasta_qual_pairs,
    bio_review.read_mast_report,
    bio_review.reverse_complement_rna_sequence,
    bio_review.tree_to_graph_edges,
    bio_review.write_gaf_file,
    bio_review.write_motifs,
    bio_review.write_pairwise_alignment_file,
    bio_review.write_qual_file,
    skbio_review.build_bme_tree,
    skbio_review.build_gme_tree,
    skbio_review.improve_tree_nni,
    skbio_review.tree_path_distances,
    skbio_review.weighted_robinson_foulds_distances,
    cutadapt_review.adapter_untrimmed_interval,
    umi_review.barcode_from_bam_tag,
    umi_review.metacontig_contig_map,
    umi_review.parse_prefixed_umi_name,
    umi_review.resolve_above_knee_barcodes,
    umi_review.soft_clip_adjusted_position,
    umi_review.tag_reads_by_metacontig,
]

SAM_TEXT = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:tx1\tLN:100\n"
    "@SQ\tSN:tx2\tLN:100\n"
    "@SQ\tSN:tx3\tLN:100\n"
    "r1\t0\ttx1\t1\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRX:Z:ACGTAC-1\n"
    "r2\t0\ttx2\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\tRX:Z:GGGGTT-1\n"
    "r3\t16\ttx2\t10\t60\t3S5M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
)
GENE_MAP = "GENE1\ttx1\nGENE1\ttx2\nGENE2\ttx3\n"
MAST_XML = """<?xml version="1.0" encoding="UTF-8"?>
<mast version="5.5.0">
<sequence_dbs><sequence_db source="seqs.fa"/></sequence_dbs>
<alphabet name="DNA"/>
<settings strand_handling="combine"/>
<motifs>
<motif id="MOTIF1" alt="m1" length="4"/>
<motif id="MOTIF2" alt="m2" length="3"/>
</motifs>
<sequences>
<sequence name="s1" length="20">
<seg start="1"><hit idx="0" pos="3" rc="n"/></seg>
</sequence>
<sequence name="s2" length="20"/>
</sequences>
</mast>
"""


def _bam_from_sam(tmp_path: Path) -> tuple[Path, Path]:
    import pysam

    sam = tmp_path / "reads.sam"
    sam.write_text(SAM_TEXT, encoding="utf-8")
    bam = tmp_path / "reads.bam"
    with (
        pysam.AlignmentFile(str(sam), "r") as source,
        pysam.AlignmentFile(str(bam), "wb", template=source) as target,
    ):
        for record in source:
            target.write(record)
    pysam.index(str(bam))
    return sam, bam


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_every_wrapper_has_native_keyword_only_verbose() -> None:
    for function in WRAPPERS:
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY, function.__name__
        assert parameter.default is True, function.__name__


def test_modules_export_only_new_names_and_no_collisions() -> None:
    import inspect

    import dnakit.alignment as alignment
    import dnakit.sequence as sequence
    import dnakit.sequencing_qc as sequencing_qc

    by_package = (
        (sequence, (bio_review, skbio_review, pyfastx_review)),
        (alignment, (parasail_review, mappy_review, pyhmmer_review)),
        (sequencing_qc, (cutadapt_review, umi_review)),
    )
    exported = [name for _, modules in by_package for module in modules for name in module.__all__]
    assert len(exported) == len(set(exported))
    # Each area package must expose exactly the review function under its name (nothing shadowed).
    for package, modules in by_package:
        for module in modules:
            for name in module.__all__:
                assert name in package.__all__
                assert inspect.unwrap(getattr(package, name)) is inspect.unwrap(
                    getattr(module, name)
                )
    assert parasail_review.__all__ == () and mappy_review.__all__ == ()
    assert pyhmmer_review.__all__ == () and pyfastx_review.__all__ == ()


def test_pairwise_alignment_file_matches_hand_scores(tmp_path: Path) -> None:
    result = bio_review.write_pairwise_alignment_file(
        "ACGTACGT", "ACGAACGT", tmp_path / "pair.fa", verbose=False
    )
    assert result["score"] == 13.0  # seven matches at 2 and one mismatch at -1
    assert result["aligned_query"] == "ACGTACGT" and result["aligned_target"] == "ACGAACGT"
    assert (tmp_path / "pair.fa").read_text().startswith(">query")
    with pytest.raises(ConfigurationError):
        bio_review.write_pairwise_alignment_file(
            "ACGT", "ACGT", tmp_path / "pair.fa", verbose=False
        )


def test_convert_search_results_blast_tab_round_trip(tmp_path: Path) -> None:
    source = _write(
        tmp_path / "hits.tsv",
        "q1\tt1\t100.00\t8\t0\t0\t1\t8\t1\t8\t1e-05\t16.0\n"
        "q2\tt2\t87.50\t8\t1\t0\t1\t8\t1\t8\t2e-03\t12.5\n",
    )
    result = bio_review.convert_search_results(
        source, "blast-tab", tmp_path / "out.tsv", "blast-tab", verbose=False
    )
    assert result["query_count"] == 2
    assert (tmp_path / "out.tsv").read_text().splitlines()[1].startswith("q2\tt2\t87.50")
    with pytest.raises(ConfigurationError):
        bio_review.convert_search_results(
            source, "blast-tab", tmp_path / "x.xml", "nope", verbose=False
        )
    with pytest.raises(ConfigurationError):
        bio_review.convert_search_results(
            source, "blast-tab", tmp_path / "y.xml", "blast-xml", verbose=False
        )


def test_codon_back_table_standard_code() -> None:
    from Bio.Data import CodonTable

    result = bio_review.codon_back_table(1, verbose=False)
    back = result["back_table"]
    assert result["codon_count"] == 61 and result["table_name"] == "Standard"
    assert back["M"] == "ATG" and back["*"] == "TAA"
    assert back["L"] == "TTG"  # last of CTA, CTC, CTG, CTT, TTA, TTG in sorted order
    assert len(back) == 21
    with pytest.raises(ConfigurationError):
        bio_review.codon_back_table(1, stop_codon="ATG", verbose=False)
    assert CodonTable.unambiguous_dna_by_id[1].forward_table["ATG"] == "M"


def test_adam_consensus_matches_upstream_for_identical_trees() -> None:
    from Bio import Phylo
    from Bio.Phylo.Consensus import adam_consensus

    newick = "((A:1,B:1):1,(C:1,D:1):1);"
    result = bio_review.adam_consensus_tree([newick, newick], verbose=False)
    upstream = adam_consensus(
        [Phylo.read(io.StringIO(newick), "newick"), Phylo.read(io.StringIO(newick), "newick")]
    )
    assert result["newick"] == upstream.format("newick").strip()
    assert result["tree_count"] == 2


def test_bootstrap_support_hand_calculation() -> None:
    target = "((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.3);"
    replicates = [target, "((A:0.1,C:0.2):0.3,(B:0.1,D:0.2):0.3);", target]
    result = bio_review.bootstrap_branch_support(target, replicates, verbose=False)
    # (A,B) and (C,D) occur in two of three replicates: 2/3 = 66.67 percent.
    assert result["newick"] == "((A:0.1,B:0.2)66.67:0.3,(C:0.1,D:0.2)66.67:0.3)100.00:0;"
    assert result["replicate_count"] == 3


def test_tree_to_graph_edges_hand_listing() -> None:
    result = bio_review.tree_to_graph_edges("((A:0.1,B:0.2):0.3,C:0.5);", verbose=False)
    assert result["node_count"] == 5 and result["edge_count"] == 4 and result["tip_count"] == 3
    edges = {(edge["source"], edge["target"]): edge["branch_length"] for edge in result["edges"]}
    assert edges == {(0, 1): 0.3, (0, 4): 0.5, (1, 2): 0.1, (1, 3): 0.2}
    names = {node["id"]: node["name"] for node in result["nodes"]}
    assert names[2] == "A" and names[4] == "C" and names[0] is None


def test_rna_complement_and_reverse_complement_hand_values() -> None:
    assert bio_review.complement_rna_sequence("ACGU", verbose=False)["sequence"] == "UGCA"
    assert bio_review.reverse_complement_rna_sequence("AAGC", verbose=False)["sequence"] == "GCUU"
    with pytest.raises(ConfigurationError):
        bio_review.complement_rna_sequence("ACGT", verbose=False)


def test_fasta_qual_pairs_round_trip(tmp_path: Path) -> None:
    fasta = _write(tmp_path / "f.fa", ">r1 desc\nACGT\n>r2\nGGCC\n")
    qual = _write(tmp_path / "q.qual", ">r1 desc\n30 20 10 5\n>r2\n40 40 40 40\n")
    result = bio_review.read_fasta_qual_pairs(fasta, qual, verbose=False)
    assert result["record_count"] == 2
    assert result["records"][0]["qualities"] == [30, 20, 10, 5]
    written = bio_review.write_qual_file(
        [{"id": "r1", "sequence": "ACGT", "qualities": [30, 20, 10, 5]}],
        tmp_path / "out.qual",
        verbose=False,
    )
    assert written["record_count"] == 1
    assert (tmp_path / "out.qual").read_text().splitlines()[1] == "30 20 10 5"
    with pytest.raises(ConfigurationError):
        bio_review.write_qual_file(
            [{"id": "r1", "sequence": "ACGT", "qualities": [30]}],
            tmp_path / "bad.qual",
            verbose=False,
        )


def test_melting_temperature_gc_matches_upstream() -> None:
    from Bio.SeqUtils import MeltingTemp

    sequence = "ACGTACGTACGTACGTAAGC"
    result = bio_review.melting_temperature_gc(sequence, verbose=False)
    expected = MeltingTemp.Tm_GC(sequence, check=True, strict=True, valueset=7, Na=50)
    assert result["tm_celsius"] == pytest.approx(expected)
    with pytest.raises(ConfigurationError):
        bio_review.melting_temperature_gc("ACGTN", verbose=False)
    with pytest.raises(ConfigurationError):
        bio_review.melting_temperature_gc(sequence, valueset=9, verbose=False)


def test_gaf_file_has_seventeen_fields(tmp_path: Path) -> None:
    record = {
        "DB": "UniProtKB",
        "DB_Object_ID": "P12345",
        "DB_Object_Symbol": "GENE1",
        "GO_ID": "GO:0005634",
        "Evidence": "IDA",
        "Aspect": "C",
        "Taxon_ID": ["taxon:9606"],
    }
    result = bio_review.write_gaf_file([record], tmp_path / "x.gaf", verbose=False)
    assert result["record_count"] == 1
    fields = (tmp_path / "x.gaf").read_text().rstrip("\n").split("\t")
    assert len(fields) == 17
    assert fields[0] == "UniProtKB" and fields[4] == "GO:0005634" and fields[12] == "taxon:9606"
    with pytest.raises(ConfigurationError):
        bio_review.write_gaf_file([{"DB": "a\tb"}], tmp_path / "bad.gaf", verbose=False)


def test_write_motifs_consensus_and_pfm_limit() -> None:
    result = bio_review.write_motifs({"M1": ["ACGT", "ACGA", "ACGT"]}, "jaspar", verbose=False)
    assert result["motifs"] == [{"name": "M1", "length": 4, "consensus": "ACGT"}]
    assert result["text"].startswith(">M1 M1\nA [  3.00   0.00   0.00   1.00]")
    with pytest.raises(ConfigurationError):
        bio_review.write_motifs({"M1": ["ACGT"], "M2": ["AAAA"]}, "pfm", verbose=False)
    with pytest.raises(ConfigurationError):
        bio_review.write_motifs({"M1": ["ACGT", "ACG"]}, "jaspar", verbose=False)


def test_mast_report_hand_diagram(tmp_path: Path) -> None:
    report = _write(tmp_path / "mast.xml", MAST_XML)
    result = bio_review.read_mast_report(report, verbose=False)
    assert result["alphabet"] == "DNA" and result["strand_handling"] == "combine"
    assert [motif["name"] for motif in result["motifs"]] == ["1", "2"]
    diagrams = {item["name"]: item["diagram"] for item in result["sequences"]}
    # hit at position 3 of motif 1 (length 4): gap 2, forward hit, tail 20 - 3 - 4 + 1 = 14
    assert diagrams == {"s1": "2-[+1]-14", "s2": "20"}
    with pytest.raises(ConfigurationError):
        bio_review.read_mast_report(_write(tmp_path / "bad.xml", "<not-mast/>"), verbose=False)


def _rounded(tree: object) -> str:
    """Upstream Newick with branch lengths rounded as the wrapper does."""

    for node in tree.postorder(include_self=True):  # type: ignore[attr-defined]
        if node.length is not None:
            node.length = round(float(node.length), 10)
    return str(tree).strip()


def test_skbio_tree_builders_match_upstream_on_hand_p_distances() -> None:
    from skbio import DistanceMatrix
    from skbio.tree import bme, gme

    sequences = {
        "A": "ACGTACGTAC",
        "B": "ACGTACGAAC",
        "C": "TCGTACGTAC",
        "D": "TCGTTCGTAG",
    }
    # Hand p-distances (mismatches / 10): AB 0.1, AC 0.1, AD 0.3, BC 0.2, BD 0.4, CD 0.2.
    matrix = DistanceMatrix(
        [
            [0.0, 0.1, 0.1, 0.3],
            [0.1, 0.0, 0.2, 0.4],
            [0.1, 0.2, 0.0, 0.2],
            [0.3, 0.4, 0.2, 0.0],
        ],
        ["A", "B", "C", "D"],
    )
    assert skbio_review.build_bme_tree(sequences, verbose=False)["newick"] == _rounded(bme(matrix))
    assert skbio_review.build_gme_tree(sequences, verbose=False)["newick"] == _rounded(gme(matrix))
    with pytest.raises(ConfigurationError):
        skbio_review.build_bme_tree({"A": "ACGT", "B": "ACG"}, verbose=False)


def test_skbio_nni_rejects_mismatched_tips() -> None:
    sequences = {"A": "ACGT", "B": "ACGA", "C": "AGGA", "D": "TGGA"}
    with pytest.raises(ConfigurationError):
        skbio_review.improve_tree_nni("((A:1,C:1):1,(B:1,X:1):1);", sequences, verbose=False)
    with pytest.raises(ConfigurationError):
        skbio_review.improve_tree_nni(
            "((A:1,C:1):1,B:1);", {"A": "ACGT", "B": "ACGA", "C": "AGGA"}, verbose=False
        )
    result = skbio_review.improve_tree_nni("((A:1,C:1):1,(B:1,D:1):1);", sequences, verbose=False)
    assert result["tip_count"] == 4 and result["newick"].endswith(";")


def test_tree_path_and_wrf_distances_match_upstream() -> None:
    from skbio import TreeNode
    from skbio.tree import path_dists, wrf_dists

    newicks = ["((A:1,B:1):1,(C:1,D:1):1);", "((A:1,C:1):1,(B:1,D:1):1);"]
    path = skbio_review.tree_path_distances(newicks, ids=["t1", "t2"], verbose=False)
    # Pairs AB, CD, AC, BD differ by 2 each: sqrt(4 * 2**2) = 4.
    assert path["matrix"] == [[0.0, 4.0], [4.0, 0.0]] and path["ids"] == ["t1", "t2"]
    trees = [TreeNode.read(io.StringIO(text)) for text in newicks]
    upstream = wrf_dists(trees, ids=["t1", "t2"])
    wrf = skbio_review.weighted_robinson_foulds_distances(newicks, ids=["t1", "t2"], verbose=False)
    assert wrf["matrix"] == upstream.data.tolist()
    assert path_dists(trees, ids=["t1", "t2"]).data[0, 1] == 4.0
    with pytest.raises(ConfigurationError):
        skbio_review.tree_path_distances(newicks, ids=["same", "same"], verbose=False)


def test_adapter_untrimmed_interval_hand_values() -> None:
    result = cutadapt_review.adapter_untrimmed_interval(
        "ACGTACGTTTAGATCGGAA", ["AGATCGG"], verbose=False
    )
    assert (result["start"], result["stop"]) == (0, 10)
    assert result["untrimmed_sequence"] == "ACGTACGTTT"
    assert result["adapter_matches"] == [
        {"adapter": "AGATCGG", "start": 10, "stop": 17, "errors": 0}
    ]
    five = cutadapt_review.adapter_untrimmed_interval(
        "AGATCGGACGTACGT", ["AGATCGG"], position="5prime", verbose=False
    )
    assert (five["start"], five["stop"]) == (7, 15)
    unmatched = cutadapt_review.adapter_untrimmed_interval("ACGTACGT", ["GGGGGGGG"], verbose=False)
    assert unmatched["length"] == 8 and unmatched["matched_count"] == 0
    with pytest.raises(ConfigurationError):
        cutadapt_review.adapter_untrimmed_interval("ACGT", ["AG"], position="middle", verbose=False)


def test_barcode_from_bam_tag_hand_values() -> None:
    result = umi_review.barcode_from_bam_tag(
        {"RX": "ACGTAC-1", "CB": "AAACCT-1"}, umi_tag_split="-", verbose=False
    )
    assert result["umi"] == "ACGTAC" and result["cell_barcode"] is None
    cell = umi_review.barcode_from_bam_tag(
        {"RX": "ACGTAC-1", "CB": "AAACCT-1"}, cell_barcode=True, cell_tag="CB", verbose=False
    )
    assert cell["umi"] == "ACGTAC-1" and cell["cell_barcode"] == "AAACCT"
    with pytest.raises(ConfigurationError):
        umi_review.barcode_from_bam_tag({"CB": "AAA"}, verbose=False)


def test_metacontig_map_and_read_counts(tmp_path: Path) -> None:
    sam, bam = _bam_from_sam(tmp_path)
    gene_map = _write(tmp_path / "g2t.tsv", GENE_MAP)
    mapped = umi_review.metacontig_contig_map(sam, gene_map, verbose=False)
    assert mapped["metacontigs"] == {"GENE1": ["tx1", "tx2"], "GENE2": ["tx3"]}
    counted = umi_review.tag_reads_by_metacontig(bam, gene_map, verbose=False)
    assert counted["reads_per_metacontig"] == {"GENE1": 3, "GENE2": 0}
    assert counted["read_count"] == 3
    unindexed = _write(tmp_path / "plain.sam", SAM_TEXT)
    with pytest.raises(ConfigurationError):
        umi_review.tag_reads_by_metacontig(unindexed, gene_map, verbose=False)


def test_soft_clip_position_hand_values() -> None:
    forward = umi_review.soft_clip_adjusted_position("3S5M2N2M", reference_start=10, verbose=False)
    assert forward["position"] == 7 and forward["start"] == 7
    assert forward["is_spliced"] is True and forward["splice_offset"] == 8  # 3S + 5M
    reverse = umi_review.soft_clip_adjusted_position(
        "5M2S", reference_start=10, is_reverse=True, verbose=False
    )
    assert reverse["start"] == 10 and reverse["position"] == 17  # end 15 plus clip 2
    assert reverse["is_spliced"] is False and reverse["splice_offset"] is None
    with pytest.raises(ConfigurationError):
        umi_review.soft_clip_adjusted_position("5Q", verbose=False)


def test_prefixed_umi_name_and_above_knee_resolution() -> None:
    parsed = umi_review.parse_prefixed_umi_name(
        "HWI:1:UMI_ACGT:CELL_TTAA", cell_barcode=True, verbose=False
    )
    assert parsed["umi"] == "ACGT" and parsed["cell_barcode"] == "TTAA"
    with pytest.raises(ConfigurationError):
        umi_review.parse_prefixed_umi_name("HWI:1:ACGT", verbose=False)
    discard = umi_review.resolve_above_knee_barcodes(
        {"AAAA": 10, "AAAT": 5}, ["AAAA", "AAAT"], {}, verbose=False
    )
    assert discard["kept"] == ["AAAA"] and discard["discarded"] == ["AAAT"]
    correct = umi_review.resolve_above_knee_barcodes(
        {"AAAA": 10, "AAAT": 5},
        ["AAAA", "AAAT"],
        {"AAAA": ["CCCC"]},
        resolution_method="correct",
        verbose=False,
    )
    assert correct["true_to_false_map"] == {"AAAA": ["AAAT", "CCCC"]}


def test_no_wrapper_writes_over_an_existing_file(tmp_path: Path) -> None:
    existing = _write(tmp_path / "keep.fa", "keep")
    with pytest.raises(ConfigurationError):
        bio_review.write_pairwise_alignment_file("ACGT", "ACGT", existing, verbose=False)
    assert existing.read_text() == "keep"
