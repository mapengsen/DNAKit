from __future__ import annotations

import random

import pytest

pyensembl = pytest.importorskip("pyensembl")

from dnakit.annotation import pyensembl_p36_functions as m
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_RNG = random.Random(3)
CHR1 = "".join(_RNG.choice("ACGT") for _ in range(2000))
CHR2 = "".join(_RNG.choice("ACGT") for _ in range(1000))


def _attrs(**kw):
    return " ".join(f'{k} "{v}";' for k, v in kw.items())


def _gtf_text():
    rows = []

    def add(contig, feature, start, end, strand, **kw):
        rows.append(f"{contig}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{_attrs(**kw)}")

    gene = dict(gene_id="GA", gene_name="AAA", gene_biotype="protein_coding", gene_version="1")
    add("1", "gene", 100, 600, "+", **gene)
    for tid, tname, exons, coding in (
        ("TA1", "AAA-201", [(100, 200), (300, 400), (500, 600)], True),
        ("TA2", "AAA-202", [(100, 200), (500, 600)], False),
    ):
        tx = dict(gene, transcript_id=tid, transcript_name=tname, transcript_version="1",
                  transcript_biotype="protein_coding" if coding else "retained_intron")
        add("1", "transcript", 100, 600, "+", **tx)
        for number, (s, e) in enumerate(exons, 1):
            add("1", "exon", s, e, "+", **dict(tx, exon_number=number, exon_id=f"E{tid}{number}"))
        if coding:
            for number, (s, e) in enumerate([(150, 200), (300, 400), (500, 550)], 1):
                add("1", "CDS", s, e, "+", **dict(tx, exon_number=number, protein_id="PA1", protein_version="1"))
    gene = dict(gene_id="GB", gene_name="BBB", gene_biotype="lncRNA", gene_version="1")
    add("1", "gene", 1000, 1300, "-", **gene)
    tx = dict(gene, transcript_id="TB1", transcript_name="BBB-201", transcript_biotype="lncRNA", transcript_version="1")
    add("1", "transcript", 1000, 1300, "-", **tx)
    add("1", "exon", 1000, 1100, "-", **dict(tx, exon_number=2, exon_id="EB2"))
    add("1", "exon", 1200, 1300, "-", **dict(tx, exon_number=1, exon_id="EB1"))
    gene = dict(gene_id="GC", gene_name="CCC", gene_biotype="protein_coding", gene_version="1")
    add("2", "gene", 50, 300, "+", **gene)
    tx = dict(gene, transcript_id="TC1", transcript_name="CCC-201", transcript_biotype="protein_coding",
              transcript_version="1")
    add("2", "transcript", 50, 300, "+", **tx)
    add("2", "exon", 50, 300, "+", **dict(tx, exon_number=1, exon_id="EC1"))
    add("2", "CDS", 60, 299, "+", **dict(tx, exon_number=1, protein_id="PC1", protein_version="1"))
    return "\n".join(rows) + "\n"


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    folder = tmp_path_factory.mktemp("p36ens")
    (folder / "g.gtf").write_text(_gtf_text())
    (folder / "g.fa").write_text(f">1\n{CHR1}\n>2\n{CHR2}\n")
    (folder / "tx.fa").write_text(">TA1 cdna\nACGTACGTAC\nGGGG\n>TB1\nTTTTAAAA\n")
    (folder / "hgnc.tsv").write_text(
        "symbol\talias_symbol\tprev_symbol\tensembl_gene_id\tstatus\n"
        "BRCA1\tRNF53|BRCC1\tFANCS\tENSG00000012048\tApproved\n"
        "GONE\t\t\t\tApproved\n"
        "OLD\t\t\tENSG00000000001\tWithdrawn\n"
    )
    return {"gtf": folder / "g.gtf", "fa": folder / "g.fa", "tx": folder / "tx.fa",
            "hgnc": folder / "hgnc.tsv", "cache": folder / "cache"}


def _genome(files):
    genome = pyensembl.Genome(
        reference_name="check", annotation_name="check", gtf_path_or_url=str(files["gtf"]),
        cache_directory_path=str(files["cache"]),
    )
    genome.index()
    return genome


def test_hgnc_aliases_read_offline(files):
    result = m.load_hgnc_gene_aliases(files["hgnc"], names=["RNF53", "FANCS", "none"], verbose=False)
    assert result["aliases"] == {"RNF53": ["ENSG00000012048"], "FANCS": ["ENSG00000012048"]}
    assert result["missing"] == ["none"] and result["total_names"] == 4  # Withdrawn and ID-less rows dropped
    everything = m.load_hgnc_gene_aliases(files["hgnc"], max_names=2, verbose=False)
    assert len(everything["aliases"]) == 2
    with pytest.raises(ConfigurationError):
        m.load_hgnc_gene_aliases(verbose=False)
    with pytest.raises(ConfigurationError):
        m.load_hgnc_gene_aliases(files["hgnc"], download=True, verbose=False)


def test_genomic_sequence_strands_and_masks(files):
    plus = m.fetch_genomic_sequence(files["fa"], "1", 100, 110, verbose=False)
    minus = m.fetch_genomic_sequence(files["fa"], "1", 100, 110, strand="-", verbose=False)
    assert plus["sequence"] == CHR1[99:110]
    comp = str.maketrans("ACGT", "TGCA")
    assert minus["sequence"] == CHR1[99:110].translate(comp)[::-1] and plus["length"] == 11
    assert m.fetch_genomic_sequence(files["fa"], "2", 1, 5, verbose=False)["sequence"] == CHR2[:5]
    with pytest.raises(BackendExecutionError):
        m.fetch_genomic_sequence(files["fa"], "1", 1990, 2100, verbose=False)
    with pytest.raises(BackendExecutionError):
        m.fetch_genomic_sequence(files["fa"], "9", 1, 5, verbose=False)
    with pytest.raises(ConfigurationError):
        m.fetch_genomic_sequence(files["fa"], "1", 10, 5, verbose=False)


def test_transcript_sequences_from_fasta(files):
    result = m.fetch_transcript_sequence(files["tx"], ["TA1", "TB1", "TZ"], verbose=False)
    assert result["sequences"] == {"TA1": "ACGTACGTACGGGG", "TB1": "TTTTAAAA", "TZ": None}
    assert result["missing"] == ["TZ"] and result["lengths"]["TA1"] == 14
    with pytest.raises(ConfigurationError):
        m.fetch_transcript_sequence(files["tx"], [], verbose=False)


def test_features_at_locus_match_pyensembl(files):
    genome = _genome(files)
    cache = files["cache"]
    tx = m.query_annotation_at_locus(files["gtf"], "1", 150, 160, cache_dir=cache, verbose=False)
    assert [r["transcript_id"] for r in tx["results"]] == [t.transcript_id for t in genome.transcripts_at_locus("1", 150, 160)]
    ids = m.query_annotation_at_locus(files["gtf"], "1", 150, 160, output="ids", cache_dir=cache, verbose=False)
    assert ids["results"] == ["TA1", "TA2"]
    names = m.query_annotation_at_locus(files["gtf"], "1", 150, output="names", cache_dir=cache, verbose=False)
    assert names["results"] == ["AAA-201", "AAA-202"] and names["locus"]["end"] == 150
    exons = m.query_annotation_at_locus(files["gtf"], "1", 1050, 1060, feature="exon", cache_dir=cache, verbose=False)
    assert [e["exon_id"] for e in exons["results"]] == ["EB2"]
    proteins = m.query_annotation_at_locus(files["gtf"], "1", 160, 170, feature="protein", output="ids",
                                           cache_dir=cache, verbose=False)
    assert proteins["results"] == ["PA1"]
    assert m.query_annotation_at_locus(files["gtf"], "1", 700, 800, cache_dir=cache, verbose=False)["count"] == 0
    stranded = m.query_annotation_at_locus(files["gtf"], "1", 150, 160, strand="-", cache_dir=cache, verbose=False)
    assert stranded["count"] == 0
    capped = m.query_annotation_at_locus(files["gtf"], "1", 150, 160, max_results=1, cache_dir=cache, verbose=False)
    assert capped["truncated"] and capped["count"] == 2 and len(capped["results"]) == 1
    for bad in ({"feature": "exon", "output": "names"}, {"feature": "protein", "output": "records"}):
        with pytest.raises(ConfigurationError):
            m.query_annotation_at_locus(files["gtf"], "1", 1, cache_dir=cache, verbose=False, **bad)
    with pytest.raises(ConfigurationError):
        m.query_annotation_at_locus(files["gtf"], "7", 1, cache_dir=cache, verbose=False)


def test_merged_intervals_and_nearest_transcript(files):
    cache = files["cache"]
    merged = m.merge_gene_intervals(files["gtf"], "1", cache_dir=cache, verbose=False)
    assert merged["intervals"] == [[100, 600], [1000, 1300]] and merged["covered_bases"] == 802
    assert m.merge_gene_intervals(files["gtf"], "1", strand="-", cache_dir=cache, verbose=False)["intervals"] == [[1000, 1300]]
    near = m.find_nearest_transcript(files["gtf"], "1", 800, cache_dir=cache, verbose=False)
    assert near["distance"] == 200 and near["transcript"]["transcript_id"] == "TA1"
    inside = m.find_nearest_transcript(files["gtf"], "1", 1100, cache_dir=cache, verbose=False)
    assert inside["distance"] == 0 and inside["transcript"]["transcript_id"] == "TB1"
    stranded = m.find_nearest_transcript(files["gtf"], "1", 800, strand="-", cache_dir=cache, verbose=False)
    assert stranded["transcript"]["transcript_id"] == "TB1" and stranded["distance"] == 200


def test_get_annotation_by_id(files):
    cache = files["cache"]
    gene = m.get_annotation_by_id(files["gtf"], "gene", "GA", cache_dir=cache, verbose=False)["record"]
    assert gene["gene_name"] == "AAA" and gene["transcript_ids"] == ["TA1", "TA2"] and gene["strand"] == "+"
    tx = m.get_annotation_by_id(files["gtf"], "transcript", "TA1", cache_dir=cache, verbose=False)["record"]
    assert tx["transcript_name"] == "AAA-201" and tx["exon_ids"] == ["ETA11", "ETA12", "ETA13"]
    exon = m.get_annotation_by_id(files["gtf"], "exon", "EB1", cache_dir=cache, verbose=False)["record"]
    assert (exon["start"], exon["end"], exon["strand"], exon["gene_id"]) == (1200, 1300, "-", "GB")
    with pytest.raises(BackendExecutionError):
        m.get_annotation_by_id(files["gtf"], "gene", "NOPE", cache_dir=cache, verbose=False)
    with pytest.raises(ConfigurationError):
        m.get_annotation_by_id(files["gtf"], "protein", "PA1", cache_dir=cache, verbose=False)


def test_contigs_and_feature_lists_match_pyensembl(files):
    genome = _genome(files)
    cache = files["cache"]
    assert m.list_annotation_contigs(files["gtf"], cache_dir=cache, verbose=False)["contigs"] == genome.contigs()
    for feature, method in (("gene", "gene_ids"), ("transcript", "transcript_ids"), ("exon", "exon_ids"), ("protein", "protein_ids")):
        got = m.list_annotation_features(files["gtf"], feature, cache_dir=cache, verbose=False)
        assert got["results"] == getattr(genome, method)() and got["count"] == len(got["results"])
    plus = m.list_annotation_features(files["gtf"], "transcript", contig="1", strand="+", cache_dir=cache, verbose=False)
    assert plus["results"] == ["TA1", "TA2"]
    lnc = m.list_annotation_features(files["gtf"], "gene", biotype="lncRNA", output="records", cache_dir=cache, verbose=False)
    assert [g["gene_id"] for g in lnc["results"]] == ["GB"] and lnc["results"][0]["biotype"] == "lncRNA"
    names = m.list_annotation_features(files["gtf"], "transcript", output="names", cache_dir=cache, verbose=False)
    assert names["results"] == ["AAA-201", "AAA-202", "BBB-201", "CCC-201"]
    assert m.list_annotation_features(files["gtf"], "exon", contig="2", cache_dir=cache, verbose=False)["results"] == ["EC1"]
    capped = m.list_annotation_features(files["gtf"], "exon", max_results=3, cache_dir=cache, verbose=False)
    assert capped["truncated"] and len(capped["results"]) == 3
    for bad in ({"feature": "protein", "output": "names"}, {"feature": "exon", "biotype": "x"}, {"feature": "protein", "output": "records"}):
        with pytest.raises(ConfigurationError):
            m.list_annotation_features(files["gtf"], cache_dir=cache, verbose=False, **bad)


@pytest.mark.parametrize(
    ("source", "identifier", "target", "expected"),
    [
        ("protein_id", "PA1", "gene", ["GA"]),
        ("protein_id", "PA1", "transcript", ["TA1"]),
        ("transcript_id", "TA1", "gene", ["GA"]),
        ("transcript_name", "AAA-202", "gene", ["GA"]),
        ("transcript_id", "TA1", "exon", ["ETA11", "ETA12", "ETA13"]),
        ("transcript_name", "BBB-201", "exon", ["EB1", "EB2"]),
        ("transcript_name", "AAA-201", "transcript", ["TA1"]),
        ("exon_id", "EB1", "gene", ["GB"]),
        ("exon_id", "ETA11", "transcript", ["TA1"]),
        ("gene_id", "GA", "transcript", ["TA1", "TA2"]),
        ("gene_name", "AAA", "transcript", ["TA1", "TA2"]),
        ("gene_id", "GA", "exon", ["ETA11", "ETA12", "ETA13", "ETA21", "ETA22"]),
        ("gene_name", "CCC", "exon", ["EC1"]),
    ],
)
def test_related_features(files, source, identifier, target, expected):
    result = m.related_annotation_features(files["gtf"], source, identifier, target, cache_dir=files["cache"], verbose=False)
    key = f"{target}_id"
    assert sorted(r[key] for r in result["results"]) == sorted(expected)
    assert result["count"] == len(expected)


def test_related_features_rejects_unsupported_pairs(files):
    cache = files["cache"]
    for source, target in (("gene_id", "gene"), ("protein_id", "exon"), ("exon_id", "exon")):
        with pytest.raises(ConfigurationError):
            m.related_annotation_features(files["gtf"], source, "x", target, cache_dir=cache, verbose=False)
    with pytest.raises(BackendExecutionError):
        m.related_annotation_features(files["gtf"], "protein_id", "NOPE", "gene", cache_dir=cache, verbose=False)


def test_transcript_offset_equals_pyensembl_and_intron_is_reported(files):
    genome = _genome(files)
    transcript = genome.transcript_by_id("TA1")
    for position in (100, 160, 350, 600):
        got = m.map_position_to_transcript_offset(files["gtf"], "TA1", position, cache_dir=files["cache"], verbose=False)
        assert got["offset"] == transcript.spliced_offset(position)
    minus = m.map_position_to_transcript_offset(files["gtf"], "TB1", 1250, cache_dir=files["cache"], verbose=False)
    assert minus["offset"] == genome.transcript_by_id("TB1").spliced_offset(1250) == 50
    intron = m.map_position_to_transcript_offset(files["gtf"], "TA1", 250, cache_dir=files["cache"], verbose=False)
    assert intron["offset"] is None and "250" in intron["reason"]


def test_locus_offset_and_comparison():
    assert m.genomic_locus_offset("1", 10, 20, "+", position=13, verbose=False)["offset"] == 3
    assert m.genomic_locus_offset("1", 10, 20, "-", position=13, verbose=False)["offset"] == 7
    assert m.genomic_locus_offset("1", 10, 20, "+", range_start=12, range_end=14, verbose=False)["offset"] == [2, 4]
    with pytest.raises(BackendExecutionError):
        m.genomic_locus_offset("1", 10, 20, "+", position=99, verbose=False)
    with pytest.raises(ConfigurationError):
        m.genomic_locus_offset("1", 10, 20, "+", verbose=False)
    with pytest.raises(ConfigurationError):
        m.genomic_locus_offset("1", 10, 20, "+", position=12, range_start=11, range_end=13, verbose=False)
    a = {"contig": "1", "start": 10, "end": 20, "strand": "+"}
    overlap = m.compare_genomic_loci(a, {"contig": "1", "start": 15, "end": 40, "strand": "+"}, verbose=False)
    assert (overlap["distance"], overlap["overlap_length"]) == (0, 6)
    apart = m.compare_genomic_loci(a, {"contig": "1", "start": 30, "end": 40, "strand": "+"}, verbose=False)
    assert (apart["distance"], apart["overlap_length"]) == (10, 0)
    other = m.compare_genomic_loci(a, {"contig": "2", "start": 10, "end": 20}, verbose=False)
    assert other["distance"] is None and other["overlap_length"] == 0
    opposite = {"contig": "1", "start": 15, "end": 40, "strand": "-"}
    assert m.compare_genomic_loci(a, opposite, verbose=False)["overlap_length"] == 0
    assert m.compare_genomic_loci(a, opposite, ignore_strand=True, verbose=False)["overlap_length"] == 6
    with pytest.raises(ConfigurationError):
        m.compare_genomic_loci(a, {"contig": "1", "start": 5, "end": 1}, verbose=False)
