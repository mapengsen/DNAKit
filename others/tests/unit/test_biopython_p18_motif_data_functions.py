"""Tests for dnakit.regulation.biopython_p18_motif_data_functions (motifs, GEO, CEL, Cluster)."""

from __future__ import annotations

import pytest

pytest.importorskip("Bio")

import _biopython_p18_data as data
from Bio import Cluster, motifs

from dnakit.exceptions import ConfigurationError
from dnakit.regulation import biopython_p18_motif_data_functions as m

QUIET = {"verbose": False}

GEO = (
    "^PLATFORM = GPL1\n!Platform_title = demo platform\n#ID = probe id\n#GB_ACC = accession\n"
    "!platform_table_begin\nID\tGB_ACC\np1\tAB000001\np2\tAB000002\np3\tAB000003\n"
    "!platform_table_end\n"
    "^SAMPLE = GSM1\n!Sample_title = first sample\n!Sample_platform_id = GPL1\n"
    "#ID_REF = probe id\n#VALUE = signal\n"
    "!sample_table_begin\nID_REF\tVALUE\np1\t1.5\np2\t2.5\n!sample_table_end\n"
)

CLUSTER_TABLE = "YORF\tNAME\tGWEIGHT\ts1\ts2\ts3\n" + "".join(
    f"g{i}\tname{i}\t1\t{i % 3}\t{(i * 2) % 5}\t{(i * 3) % 7}\n" for i in range(12)
)


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# ------------------------------------------------------------------ motifs


@pytest.mark.parametrize(
    ("text", "fmt", "count", "first_length"),
    [
        (data.ALIGNACE, "alignace", 16, 12),
        (data.MEME_XML, "meme", 2, 13),
        (data.CLUSTERBUSTER, "clusterbuster", 3, 6),
        (data.PFM_FOUR_COLUMNS, "pfm-four-columns", 8, 8),
        (data.PFM_FOUR_ROWS, "pfm-four-rows", 9, 6),
        (data.XMS, "xms", 1, 14),
    ],
)
def test_read_motif_collection_matches_biopython(tmp_path, text, fmt, count, first_length):
    path = _write(tmp_path, f"m.{fmt}", text)
    result = m.read_motif_collection(path, file_format=fmt, **QUIET)
    expected = list(motifs.parse(path.open(encoding="utf-8"), fmt))
    assert result["motif_count"] == len(expected) == count
    first = result["motifs"][0]
    assert first["length"] == first_length
    assert first["consensus"] == str(expected[0].consensus)
    assert first["gc_content"] == pytest.approx(expected[0].counts.gc_content)
    for letter in expected[0].alphabet:
        assert first["counts"][letter] == [float(v) for v in expected[0].counts[letter]]


def test_read_motif_collection_meme_and_alignace_details(tmp_path):
    meme = _write(tmp_path, "m.xml", data.MEME_XML)
    result = m.read_motif_collection(meme, file_format="meme", include_instances=True, **QUIET)
    first = result["motifs"][0]
    assert first["evalue"] == pytest.approx(0.19) and first["num_occurrences"] == 7
    assert first["instance_count"] == 7 == len(first["instances"])
    assert first["instances"][0]["sequence_name"] == "INO1"
    ace = _write(tmp_path, "a.out", data.ALIGNACE)
    first_ace = m.read_motif_collection(ace, file_format="alignace", **QUIET)["motifs"][0]
    assert first_ace["score"] == pytest.approx(57.9079) and len(first_ace["mask"]) == 12
    assert "instances" not in first_ace


def test_read_motif_collection_errors(tmp_path):
    path = _write(tmp_path, "m.xml", data.MEME_XML)
    with pytest.raises(ConfigurationError, match="max_motifs"):
        m.read_motif_collection(path, file_format="meme", max_motifs=1, **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_motif_collection(path, file_format="xms", **QUIET)
    with pytest.raises(ConfigurationError):
        m.read_motif_collection(path, file_format="jaspar", **QUIET)


def test_motif_gc_content():
    result = m.motif_gc_content(
        {
            "sites": ["ACGT", "ACGA"],
            "counts": {"A": [1, 0], "C": [0, 3], "G": [0, 0], "T": [2, 0]},
        },
        **QUIET,
    )
    by_name = {item["name"]: item for item in result["motifs"]}
    assert (
        by_name["sites"]["gc_content"] == pytest.approx(0.5)
        and by_name["sites"]["input"] == "sites"
    )
    assert by_name["counts"]["gc_content"] == pytest.approx(0.5)
    reference = motifs.create(["ACGT", "ACGA"]).counts.gc_content
    assert by_name["sites"]["gc_content"] == pytest.approx(reference)
    all_gc = m.motif_gc_content({"gc": ["GGCC", "CCGG"]}, **QUIET)["motifs"][0]["gc_content"]
    assert all_gc == 1.0
    with pytest.raises(ConfigurationError):
        m.motif_gc_content({"bad": ["ACGX"]}, **QUIET)
    with pytest.raises(ConfigurationError):
        m.motif_gc_content({"bad": ["ACG", "AC"]}, **QUIET)
    with pytest.raises(ConfigurationError):
        m.motif_gc_content({"bad": {"A": [1], "C": [1]}}, **QUIET)


# ------------------------------------------------------------------ GEO, CEL, UniGene


def test_read_geo_soft(tmp_path):
    path = _write(tmp_path, "g.soft", GEO)
    result = m.read_geo_soft(path, max_table_rows=1, **QUIET)
    assert [r["entity_type"] for r in result["records"]] == ["PLATFORM", "SAMPLE"]
    platform, sample = result["records"]
    assert platform["table_header"] == ["ID", "GB_ACC"]
    assert platform["table_row_count"] == 3 and len(platform["table_rows"]) == 1
    assert sample["entity_id"] == "GSM1" and sample["table_rows"] == [["p1", "1.5"]]
    assert sample["table_row_count"] == 2
    assert sample["entity_attributes"]["Sample_title"] == "first sample"
    assert sample["col_defs"]["VALUE"] == "signal"
    with pytest.raises(ConfigurationError, match="max_records"):
        m.read_geo_soft(path, max_records=1, **QUIET)


def test_read_affymetrix_cel_versions(tmp_path):
    v3 = _write(tmp_path, "v3.CEL", data.CEL_V3)
    v4 = tmp_path / "v4.CEL"
    v4.write_bytes(data.CEL_V4)
    first = m.read_affymetrix_cel(v3, include_values=True, **QUIET)
    second = m.read_affymetrix_cel(v4, **QUIET)
    assert (first["version"], second["version"]) == (3, 4)
    assert first["rows"] == first["columns"] == 5 and first["cell_count"] == 25
    assert len(first["intensities"]) == 5 and len(first["intensities"][0]) == 5
    assert first["intensity_mean"] == pytest.approx(sum(sum(r) for r in first["intensities"]) / 25)
    assert "intensities" not in second and second["algorithm"] == "Percentile"
    with pytest.raises(ConfigurationError, match="max_cells"):
        m.read_affymetrix_cel(v3, include_values=True, max_cells=10, **QUIET)
    junk = _write(tmp_path, "junk.CEL", "not a cel file\n")
    with pytest.raises(ConfigurationError):
        m.read_affymetrix_cel(junk, **QUIET)


def test_read_unigene(tmp_path):
    path = _write(tmp_path, "u.data", data.UNIGENE)
    result = m.read_unigene(path, **QUIET)
    assert result["record_count"] == 4
    first = result["records"][0]
    assert first["ID"] == "Mdm.1" and first["title"].startswith("Similar to cytochrome")
    assert len(first["sequence"]) == 8 and first["sequence"][0]["acc"] == "XM_001377126.1"
    assert first["sequence"][5]["seqtype"] == "EST"
    with pytest.raises(ConfigurationError, match="max_records"):
        m.read_unigene(path, max_records=2, **QUIET)


# ------------------------------------------------------------------ Cluster


def test_read_cluster_expression_data_matches_biopython(tmp_path):
    path = _write(tmp_path, "e.txt", CLUSTER_TABLE)
    result = m.read_cluster_expression_data(path, include_data=True, **QUIET)
    record = Cluster.read(path.open(encoding="utf-8"))
    assert (result["gene_count"], result["experiment_count"]) == record.data.shape == (12, 3)
    assert result["gene_ids"][:2] == ["g0", "g1"] and result["gene_names"][1] == "name1"
    assert result["experiment_ids"] == ["s1", "s2", "s3"]
    assert result["data"] == record.data.tolist()
    assert result["gene_weights"] is not None
    brief = m.read_cluster_expression_data(path, **QUIET)
    assert "data" not in brief
    with pytest.raises(ConfigurationError):
        m.read_cluster_expression_data(_write(tmp_path, "bad.txt", "x\ty\nnot\tnumbers\n"), **QUIET)


def test_cluster_expression_data_hierarchical(tmp_path):
    path = _write(tmp_path, "e.txt", CLUSTER_TABLE)
    prefix = tmp_path / "hier"
    result = m.cluster_expression_data(path, prefix, cluster_samples=True, **QUIET)
    names = sorted(p.rsplit(".", 1)[1] for p in result["output_files"])
    assert names == ["atr", "cdt", "gtr"]
    assert result["gene_tree_nodes"] == 11 and result["sample_tree_nodes"] == 2
    gtr = (tmp_path / "hier.gtr").read_text(encoding="utf-8").splitlines()
    assert len(gtr) == 11 and gtr[0].startswith("NODE1X\t")
    record = Cluster.read(path.open(encoding="utf-8"))
    tree = record.treecluster(method="m", dist="e")
    # TreeView stores the similarity 1 - distance in the GTR file
    assert tree[0].distance == pytest.approx(1.0 - float(gtr[0].split("\t")[3]))
    with pytest.raises(ConfigurationError, match="already exist"):
        m.cluster_expression_data(path, prefix, **QUIET)


def test_cluster_expression_data_kmeans_and_errors(tmp_path):
    path = _write(tmp_path, "e.txt", CLUSTER_TABLE)
    initial = [i % 3 for i in range(12)]
    result = m.cluster_expression_data(
        path,
        tmp_path / "km",
        algorithm="kmeans",
        cluster_count=3,
        initial_clusters=initial,
        **QUIET,
    )
    assert sorted(p.rsplit(".", 1)[1] for p in result["output_files"]) == ["cdt", "kgg"]
    assert len(result["cluster_labels"]) == 12 and set(result["cluster_labels"]) <= {0, 1, 2}
    record = Cluster.read(path.open(encoding="utf-8"))
    labels, error, _ = record.kcluster(nclusters=3, initialid=initial)
    assert list(result["cluster_labels"]) == [int(v) for v in labels]
    assert result["within_cluster_error"] == pytest.approx(error)
    with pytest.raises(ConfigurationError, match="every cluster number"):
        m.cluster_expression_data(
            path,
            tmp_path / "k2",
            algorithm="kmeans",
            cluster_count=3,
            initial_clusters=[0] * 12,
            **QUIET,
        )
    with pytest.raises(ConfigurationError):
        m.cluster_expression_data(
            path, tmp_path / "k3", algorithm="kmeans", cluster_count=99, **QUIET
        )
    with pytest.raises(ConfigurationError):
        m.cluster_expression_data(path, tmp_path / "k4", distance="manhattan", **QUIET)
    assert not list(tmp_path.glob("k[234]*"))


def test_public_functions_have_verbose_keyword_and_docstrings():
    import inspect

    for name in m.__all__:
        function = getattr(m, name)
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True
        assert (function.__doc__ or "").strip()
