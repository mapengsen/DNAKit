"""Tests for dnakit.annotation.gffutils_p37_functions (look-ups, derived features, database edits, GFF writing)."""

from __future__ import annotations

import pytest

gffutils = pytest.importorskip("gffutils")

from dnakit.annotation import gffutils_p37_functions as g
from dnakit.exceptions import ConfigurationError

GFF = """##gff-version 3
chr1\tdemo\tgene\t101\t900\t.\t+\t.\tID=geneA;Name=alpha
chr1\tdemo\tmRNA\t101\t900\t.\t+\t.\tID=tA1;Parent=geneA
chr1\tdemo\texon\t101\t200\t.\t+\t.\tID=eA1;Parent=tA1
chr1\tdemo\texon\t301\t400\t.\t+\t.\tID=eA2;Parent=tA1
chr1\tdemo\texon\t601\t900\t.\t+\t.\tID=eA3;Parent=tA1
chr1\tdemo\tCDS\t151\t200\t.\t+\t0\tID=cA1;Parent=tA1
chr1\tdemo\tgene\t1201\t1500\t.\t-\t.\tID=geneB;Name=beta
chr1\tdemo\tmRNA\t1201\t1500\t.\t-\t.\tID=tB1;Parent=geneB
chr1\tdemo\texon\t1201\t1300\t.\t-\t.\tID=eB1;Parent=tB1
chr1\tdemo\texon\t1401\t1500\t.\t-\t.\tID=eB2;Parent=tB1
"""


@pytest.fixture()
def gff(tmp_path):
    path = tmp_path / "demo.gff3"
    path.write_text(GFF, encoding="utf-8")
    return path


def test_fetch_and_list_match_gffutils(gff):
    fetched = g.fetch_gffutils_feature(gff, "tA1", verbose=False)
    db = gffutils.create_db(str(gff), ":memory:", verbose=False)
    assert fetched["gff_line"] == str(db["tA1"])
    assert fetched["feature"]["start"] == 101
    assert g.list_gffutils_features(gff, featuretype="gene", verbose=False)["count"] == 2
    assert g.list_gffutils_features(gff, verbose=False)["count"] == sum(
        1 for _ in db.all_features()
    )
    assert g.list_gffutils_features(gff, seqid="chrX", verbose=False)["count"] == 0
    with pytest.raises(ConfigurationError):
        g.fetch_gffutils_feature(gff, "missing", verbose=False)


def test_parent_children_groups_and_parents(gff):
    groups = g.group_gffutils_parent_children(gff, verbose=False)["groups"]
    assert [group["parent"]["id"] for group in groups] == ["geneA", "geneB"]
    db = gffutils.create_db(str(gff), ":memory:", verbose=False)
    direct = next(iter(db.iter_by_parent_childs(featuretype="gene")))
    # The first element is the parent; only the remaining five records are children.
    assert [child["id"] for child in groups[0]["children"]] == [child.id for child in direct[1:]]
    parents = g.fetch_gffutils_parents(gff, "eA1", verbose=False)["parents"]
    assert {p["id"] for p in parents} == {"tA1", "geneA"}
    assert [
        p["id"] for p in g.fetch_gffutils_parents(gff, "eA1", level=1, verbose=False)["parents"]
    ] == ["tA1"]


def test_interfeatures_introns_splice_sites_and_children_bp(gff):
    gap = g.get_gffutils_interfeatures(gff, verbose=False)["interfeatures"]
    assert (gap[0]["start"], gap[0]["end"]) == (901, 1200)
    introns = g.create_gffutils_introns(gff, verbose=False)["features"]
    assert [(i["start"], i["end"]) for i in introns] == [(201, 300), (401, 600), (1301, 1400)]
    sites = g.create_gffutils_splice_sites(gff, verbose=False)
    assert sites["count"] == 6
    assert g.gffutils_children_bp(gff, "geneA", verbose=False)["base_pairs"] == 500
    assert (
        g.gffutils_children_bp(gff, "geneA", child_featuretype="CDS", verbose=False)["base_pairs"]
        == 50
    )


def test_derived_features_are_written_only_to_new_files(gff, tmp_path):
    out = tmp_path / "introns.gff3"
    g.create_gffutils_introns(gff, output_path=out, verbose=False)
    assert out.read_text().startswith("##gff-version 3") and out.read_text().count("intron") == 3
    with pytest.raises(ConfigurationError):
        g.create_gffutils_introns(gff, output_path=out, verbose=False)


def test_database_edits_write_new_databases(gff, tmp_path):
    deleted = g.delete_gffutils_features(gff, tmp_path / "d.db", featuretype="exon", verbose=False)
    assert (deleted["features_before"], deleted["features_after"], deleted["deleted"]) == (10, 5, 5)
    assert (
        g.delete_gffutils_features(gff, tmp_path / "d2.db", feature_ids=["cA1"], verbose=False)[
            "features_after"
        ]
        == 9
    )
    with pytest.raises(ConfigurationError):
        g.delete_gffutils_features(gff, tmp_path / "d3.db", verbose=False)
    updated = g.update_gffutils_features(
        gff, ["chr1\tx\tgene\t5\t50\t.\t+\t.\tID=geneC"], tmp_path / "u.db", verbose=False
    )
    assert (updated["features_before"], updated["features_after"]) == (10, 11)
    assert gffutils.FeatureDB(str(tmp_path / "u.db"))["geneC"].end == 50
    related = g.add_gffutils_relation(gff, "geneB", "tA1", 1, tmp_path / "r.db", verbose=False)
    assert related["children_after"] == related["children_before"] + 1
    assert gff.read_text() == GFF  # the input is untouched


def test_feature_sequence_follows_strand(gff, tmp_path):
    fasta = tmp_path / "genome.fa"
    sequence = "".join("ACGT"[i * 7 % 4] for i in range(1600))
    fasta.write_text(f">chr1\n{sequence}\n")
    plus = g.extract_gffutils_feature_sequence(gff, fasta, "eA1", verbose=False)
    assert plus["sequence"] == sequence[100:200]
    minus = g.extract_gffutils_feature_sequence(gff, fasta, "eB1", verbose=False)
    complement = str.maketrans("ACGT", "TGCA")
    assert minus["sequence"] == sequence[1200:1300].translate(complement)[::-1]
    raw = g.extract_gffutils_feature_sequence(gff, fasta, "eB1", use_strand=False, verbose=False)
    assert raw["sequence"] == sequence[1200:1300]


def test_write_gff_and_canonical_order(gff, tmp_path):
    plain = g.write_gffutils_gff(gff, tmp_path / "all.gff3", verbose=False)
    assert plain["records_written"] == 10
    only = g.write_gffutils_gff(gff, tmp_path / "genes.gff3", featuretype="gene", verbose=False)
    assert only["records_written"] == 2
    canonical = g.write_gffutils_gff(
        gff, tmp_path / "canon.gff3", canonical_gene_order=True, verbose=False
    )
    lines = [
        x.split("\t")[2]
        for x in (tmp_path / "canon.gff3").read_text().splitlines()
        if x and not x.startswith("#")
    ]
    assert lines[:2] == ["gene", "mRNA"] and canonical["records_written"] == 10


def test_remove_duplicates_merges_repeated_ids(tmp_path):
    source = tmp_path / "dup.gff3"
    line = "chr1\tx\tgene\t1\t100\t.\t+\t.\tID=g1;Name=a\n"
    source.write_text(
        "##gff-version 3\n" + line + line + "chr1\tx\tgene\t200\t300\t.\t+\t.\tID=g2\n"
    )
    result = g.remove_gffutils_duplicates(source, tmp_path / "out.gff3", verbose=False)
    assert (result["records_before"], result["records_after"], result["removed"]) == (3, 2, 1)


def test_plot_gene_model_saves_image(gff, tmp_path):
    pytest.importorskip("pybedtools")
    pytest.importorskip("matplotlib")
    result = g.plot_gffutils_gene_model(
        gff, "geneA", tmp_path / "gene.png", utrs=["UTR"], verbose=False
    )
    assert (tmp_path / "gene.png").stat().st_size > 1000
    assert result["transcripts"][0]["cds_count"] == 1
