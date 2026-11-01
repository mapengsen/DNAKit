"""Tests for dnakit.sequence.skbio_p20_io_functions (scikit-bio file I/O, sequences, alignments)."""

from __future__ import annotations

import pytest

pytest.importorskip("skbio")

from dnakit.exceptions import ConfigurationError
from dnakit.sequence import skbio_p20_io_functions as io_fn

NODES_ROW = (
    "\t|\t".join(["{id}", "{parent}", "{rank}", "", "8", "0", "1", "0", "0", "0", "0", "0", ""])
    + "\t|\n"
)


def test_distance_matrix_round_trip_lsmat_and_binary(tmp_path):
    matrix = [[0, 1, 2], [1, 0, 1.5], [2, 1.5, 0]]
    for fmt, name in (("lsmat", "d.lsmat"), ("binary_dm", "d.h5")):
        out = io_fn.write_skbio_distance_matrix(
            matrix, str(tmp_path / name), ids=["a", "b", "c"], file_format=fmt, verbose=False
        )
        assert out["size"] == 3
        back = io_fn.read_skbio_file(str(tmp_path / name), fmt, verbose=False)
        assert back["ids"] == ["a", "b", "c"]
        assert back["matrix"] == [[0.0, 1.0, 2.0], [1.0, 0.0, 1.5], [2.0, 1.5, 0.0]]


def test_phylip_distance_matrix_read(tmp_path):
    path = tmp_path / "d.phy"
    path.write_text("2\na 0 1\nb 1 0\n")
    assert io_fn.read_skbio_file(str(path), "phylip_dm", verbose=False)["ids"] == ["a", "b"]


def test_read_alignments_tree_and_qseq(tmp_path):
    (tmp_path / "c.aln").write_text("CLUSTAL W (1.83)\n\ns1 ACGT-A\ns2 ACGTTA\n\n")
    (tmp_path / "p.phy").write_text(" 2 6\ns1        ACGT-A\ns2        ACGTTA\n")
    (tmp_path / "t.nwk").write_text("((a:1,b:2)90:1,c:3);")
    (tmp_path / "q.txt").write_text("M1\t1\t1\t1\t1\t1\t1\t1\tACG.\tIIII\t1\n")
    for name, fmt in (("c.aln", "clustal"), ("p.phy", "phylip")):
        result = io_fn.read_skbio_file(str(tmp_path / name), fmt, verbose=False)
        assert result["alignment"] == {"s1": "ACGT-A", "s2": "ACGTTA"}
        assert result["alignment_length"] == 6
    tree = io_fn.read_skbio_file(str(tmp_path / "t.nwk"), "newick", verbose=False)
    assert tree["tip_count"] == 3
    reads = io_fn.read_skbio_file(str(tmp_path / "q.txt"), "qseq", verbose=False)
    assert reads["records"][0]["sequence"] == "ACG." and reads["records"][0]["quality"] == [40] * 4


def test_ordination_metadata_and_biom_round_trip(tmp_path):
    io_fn.write_skbio_ordination(
        [[0.1, 0.2], [0.3, -0.2], [-0.4, 0.0]],
        [2.0, 1.0],
        str(tmp_path / "o.ord"),
        sample_ids=["a", "b", "c"],
        proportion_explained=[0.6, 0.4],
        verbose=False,
    )
    ordination = io_fn.read_skbio_file(str(tmp_path / "o.ord"), "ordination", verbose=False)[
        "ordination"
    ]
    assert ordination["samples"]["index"] == ["a", "b", "c"]
    io_fn.write_skbio_sample_metadata(
        {"g": ["x", "y"], "n": [1, 2]}, ["a", "b"], str(tmp_path / "m.tsv"), verbose=False
    )
    meta = io_fn.read_skbio_file(str(tmp_path / "m.tsv"), "sample_metadata", verbose=False)
    assert meta["column_types"] == {"g": "categorical", "n": "numeric"} and meta["ids"] == [
        "a",
        "b",
    ]
    pytest.importorskip("biom")
    io_fn.write_skbio_biom_table([[1, 2], [3, 4]], str(tmp_path / "t.biom"), verbose=False)
    table = io_fn.read_skbio_file(str(tmp_path / "t.biom"), "biom", verbose=False)
    assert table["counts_features_by_samples"] == [[1.0, 3.0], [2.0, 4.0]]


def test_taxdump_read_and_errors(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text(
        NODES_ROW.format(id=1, parent=1, rank="no rank")
        + NODES_ROW.format(id=2, parent=1, rank="superkingdom")
    )
    result = io_fn.read_skbio_file(str(nodes), "taxdump", verbose=False)
    assert result["total"] == 2 and result["data"][1][1] == "superkingdom"
    with pytest.raises(ConfigurationError):
        io_fn.read_skbio_file(str(tmp_path / "missing"), "newick", verbose=False)
    with pytest.raises(ConfigurationError):
        io_fn.read_skbio_file(str(nodes), "nope", verbose=False)
    with pytest.raises(ConfigurationError):
        io_fn.read_skbio_file(str(nodes), "newick", verbose=False)


def test_sequence_queries():
    assert (
        io_fn.lowercase_masked_sequence("ACGTACGT", [[1, 3], [6, 8]], verbose=False)["sequence"]
        == "AcgTACgt"
    )
    assert io_fn.find_subsequence("ACGT", "GT", verbose=False)["index"] == 2
    assert io_fn.find_subsequence("ACGT", "TT", verbose=False)["found"] is False
    segments = io_fn.split_contiguous_segments(
        "ACG-TNA", included_characters="-N", invert=True, verbose=False
    )
    assert [(s["start"], s["end"], s["sequence"]) for s in segments["segments"]] == [
        (0, 3, "ACG"),
        (4, 5, "T"),
        (6, 7, "A"),
    ]
    short = io_fn.split_contiguous_segments(
        "ACG-TNA", included_intervals=[[0, 3], [6, 7]], min_length=2, verbose=False
    )
    assert [s["sequence"] for s in short["segments"]] == ["ACG"]
    assert io_fn.locate_sequence_positions("ACG-TNA", "gaps", verbose=False)["positions"] == [3]
    assert io_fn.locate_sequence_positions("ACG-TNA", "degenerates", verbose=False)[
        "positions"
    ] == [5]
    assert io_fn.locate_sequence_positions("ACG-TNA", "definites", verbose=False)["runs"] == [
        [0, 3],
        [4, 5],
        [6, 7],
    ]
    assert (
        io_fn.locate_sequence_positions("NNN", "definites", verbose=False)["has_definites"] is False
    )
    motifs = io_fn.find_nucleotide_motifs(
        "ACGGGGAGGCGGAG", "purine-run", min_length=2, verbose=False
    )
    assert [(m["start"], m["end"]) for m in motifs["motifs"]] == [(2, 9), (10, 14)]
    assert (
        io_fn.is_reverse_complement_pair("ACGT", "ACGT", verbose=False)["is_reverse_complement"]
        is True
    )
    assert (
        io_fn.is_reverse_complement_pair("AACG", "ACGT", verbose=False)["is_reverse_complement"]
        is False
    )
    with pytest.raises(ConfigurationError):
        io_fn.find_subsequence("ACGU", "A", verbose=False)
    with pytest.raises(ConfigurationError):
        io_fn.split_contiguous_segments("ACGT", verbose=False)


def test_msa_functions_and_alignment_paths():
    alignment = {"a": "AC-T", "b": "ACGT", "c": "ACGA"}
    assert io_fn.msa_consensus_sequence(alignment, verbose=False)["consensus"] == "ACGT"
    conservation = io_fn.msa_conservation(alignment, verbose=False)["conservation"]
    assert conservation[:3] == [1.0, 1.0, None] and conservation[3] == pytest.approx(
        0.5408520829727552
    )
    assert io_fn.msa_gap_frequencies(alignment, verbose=False)["frequencies"] == [0, 0, 1, 0]
    per_sequence = io_fn.msa_gap_frequencies(
        alignment, axis="position", relative=True, verbose=False
    )
    assert per_sequence["frequencies"] == [0.25, 0.0, 0.0] and per_sequence["names"] == [
        "a",
        "b",
        "c",
    ]
    path = io_fn.parse_cigar_alignment_path("3M1D2M", verbose=False)
    assert (
        path["lengths"] == [3, 1, 2] and path["alignment_length"] == 6 and path["cigar"] == "3M1D2M"
    )
    gapped = io_fn.apply_alignment_path(
        {"x": "ACGTAC", "y": "ACGAC"}, cigar="3M1D2M", verbose=False
    )
    assert gapped["alignment"] == {"x": "ACGTAC", "y": "ACG-AC"}
    manual = io_fn.apply_alignment_path(
        {"x": "ACGTAC", "y": "ACGAC"},
        lengths=[3, 1, 2],
        states=[[0, 0], [0, 1], [0, 0]],
        starts=[0, 0],
        verbose=False,
    )
    assert manual["alignment"] == gapped["alignment"]
    with pytest.raises(ConfigurationError):
        io_fn.msa_consensus_sequence({"a": "AC", "b": "ACG"}, verbose=False)
    with pytest.raises(ConfigurationError):
        io_fn.apply_alignment_path({"x": "AC", "y": "AC"}, verbose=False)
