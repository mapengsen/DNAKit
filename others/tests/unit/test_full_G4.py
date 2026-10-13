"""Unit tests for the G4 wrappers: polars-bio, gtars, ncls, and genomicranges annotation functions.

Each wrapped function is checked against a hand calculation or a direct call to the upstream
package. All inputs are written into pytest's tmp_path; no network access is used.
"""

from __future__ import annotations

import pytest

from dnakit.annotation import genomicranges_functions as gr_api
from dnakit.annotation import gtars_functions as gt_api
from dnakit.annotation import polars_bio_functions as pb_api
from dnakit.exceptions import ConfigurationError

A_BED = "chr1\t100\t200\ta\nchr1\t150\t300\tb\nchr1\t400\t500\tc\nchr2\t0\t50\td\n"
B_BED = "chr1\t180\t220\tx\nchr1\t450\t460\ty\nchr2\t10\t20\tz\n"
GTF = "\n".join(
    [
        'chr1\tensembl\tgene\t11\t100\t.\t+\t.\tgene_id "G1"; gene_name "ABC";',
        'chr1\tensembl\texon\t11\t50\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        'chr1\tensembl\tgene\t200\t300\t.\t-\t.\tgene_id "G2"; gene_name "DEF";',
        "",
    ]
)
GFF3 = "\n".join(
    [
        "##gff-version 3",
        "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=GENE1",
        "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1",
        "",
    ]
)
VCF = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=1000>",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
        "chr1\t10\trs1\tA\tG\t50\tPASS\tDP=12\tGT\t0/1\t1/1",
        "chr1\t25\t.\tC\tT\t40\tPASS\tDP=8\tGT\t0/0\t0/1",
        "",
    ]
)
FASTA = ">chr1 first\nACGTACGTAC\nGGTT\n>chr2\nTTTTAAAACCCC\n"
FASTQ = "@r1\nACGTN\n+\nIIII#\n@r2\nGGCCAA\n+\nIIIIII\n"
SAM = "\n".join(
    [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        "r1\t0\tchr1\t11\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0",
        "r2\t0\tchr1\t13\t60\t4M\t*\t0\t0\tGTAC\tIIII\tNM:i:1",
        "",
    ]
)
PAIRS = "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
PAIRS += "r1\tchr1\t100\tchr1\t500\t+\t-\nr2\tchr1\t200\tchr2\t40\t+\t+\n"
STO = "# STOCKHOLM 1.0\n#=GF ID test_family\nseq1  ACGT-ACGT\nseq2  ACGTTACGA\n//\n"
GENES = "\n".join(
    [
        'chr1\tens\tgene\t1001\t2000\t.\t+\t.\tgene_id "G1"; gene_name "ABC"; '
        'gene_type "protein_coding";',
        'chr1\tens\texon\t1001\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_type "protein_coding";',
        'chr1\tens\texon\t1801\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'gene_type "protein_coding";',
        'chr1\tens\tgene\t5001\t6000\t.\t-\t.\tgene_id "G2"; gene_name "DEF"; '
        'gene_type "protein_coding";',
        'chr1\tens\texon\t5001\t6000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; '
        'gene_type "protein_coding";',
        "",
    ]
)
REGIONS = "chr1\t950\t1010\tp1\nchr1\t1500\t1600\tp2\nchr1\t5950\t6100\tp3\nchr1\t9000\t9100\tp4\n"
GENOME = ">chr1\nACGTACGTAC\nGGTT\n>chr2\nTTTTAAAACCCC\n"


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _genome_fasta(tmp_path):
    return _write(tmp_path, "genome.fa", GENOME)


# polars-bio interval arithmetic.


def test_overlap_intervals_matches_hand_pairs_and_upstream(tmp_path):
    pb = pytest.importorskip("polars_bio")
    left = _write(tmp_path, "a.bed", A_BED)
    right = _write(tmp_path, "b.bed", B_BED)
    result = pb_api.overlap_intervals(left, right, verbose=False)
    pairs = {
        (row["start_1"], row["end_1"], row["start_2"], row["end_2"]) for row in result["pairs"]
    }
    assert pairs == {
        (100, 200, 180, 220),
        (150, 300, 180, 220),
        (400, 500, 450, 460),
        (0, 50, 10, 20),
    }
    upstream = pb.overlap(
        pb.read_bed(str(left), use_zero_based=True),
        pb.read_bed(str(right), use_zero_based=True),
        output_type="polars.DataFrame",
    )
    assert result["pair_count"] == upstream.height == 4


def test_count_interval_overlaps_counts_each_source_interval(tmp_path):
    pytest.importorskip("polars_bio")
    left = _write(tmp_path, "a.bed", A_BED)
    right = _write(tmp_path, "b.bed", B_BED)
    result = pb_api.count_interval_overlaps(left, right, verbose=False)
    counts = {(row["chrom"], row["start"]): row["count"] for row in result["intervals"]}
    assert counts == {("chr1", 100): 1, ("chr1", 150): 1, ("chr1", 400): 1, ("chr2", 0): 1}
    assert result["total_overlaps"] == 4


def test_cluster_intervals_numbers_clusters_in_coordinate_order(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "a.bed", A_BED)
    result = pb_api.cluster_intervals(path, verbose=False)
    assert result["cluster_count"] == 3
    by_name = {row["name"]: row for row in result["intervals"]}
    assert by_name["a"]["cluster"] == by_name["b"]["cluster"] == 0
    assert (by_name["a"]["cluster_start"], by_name["a"]["cluster_end"]) == (100, 300)
    assert by_name["c"]["cluster"] == 1
    assert by_name["d"]["cluster"] == 2


def test_complement_intervals_returns_gaps_inside_chromosomes(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "a.bed", A_BED)
    result = pb_api.complement_intervals(
        path, chrom_sizes={"chr1": 1000, "chr2": 100}, verbose=False
    )
    gaps = [(row["chrom"], row["start"], row["end"]) for row in result["gaps"]]
    assert gaps == [
        ("chr1", 0, 100),
        ("chr1", 300, 400),
        ("chr1", 500, 1000),
        ("chr2", 50, 100),
    ]


def test_subtract_intervals_removes_covered_bases(tmp_path):
    pytest.importorskip("polars_bio")
    left = _write(tmp_path, "a.bed", A_BED)
    right = _write(tmp_path, "b.bed", B_BED)
    result = pb_api.subtract_intervals(left, right, verbose=False)
    pieces = {(row["chrom"], row["start"], row["end"]) for row in result["pieces"]}
    assert ("chr1", 100, 180) in pieces
    assert ("chr1", 150, 180) in pieces
    assert ("chr1", 220, 300) in pieces
    assert ("chr1", 400, 450) in pieces
    assert ("chr1", 460, 500) in pieces
    assert ("chr2", 0, 10) in pieces
    assert ("chr2", 20, 50) in pieces
    assert result["piece_count"] == 7


# polars-bio gene and variant files.


def test_read_gene_annotation_gtf_and_gff(tmp_path):
    pb = pytest.importorskip("polars_bio")
    gtf = _write(tmp_path, "genes.gtf", GTF)
    result = pb_api.read_gene_annotation(gtf, attributes=["gene_id", "gene_name"], verbose=False)
    genes = {(row["type"], row["gene_id"], row["start"], row["end"]) for row in result["features"]}
    assert ("gene", "G1", 10, 100) in genes
    assert ("exon", "G1", 10, 50) in genes
    assert ("gene", "G2", 199, 300) in genes
    upstream = pb.read_gtf(str(gtf), attr_fields=["gene_id", "gene_name"], use_zero_based=True)
    assert result["feature_count"] == upstream.height == 3
    gff = _write(tmp_path, "genes.gff3", GFF3)
    gff_result = pb_api.read_gene_annotation(gff, file_format="gff", verbose=False)
    assert gff_result["feature_count"] == 2
    assert {row["start"] for row in gff_result["features"]} == {99}


def test_read_variant_calls_vcf_and_bcf(tmp_path):
    pysam = pytest.importorskip("pysam")
    vcf = _write(tmp_path, "v.vcf", VCF)
    result = pb_api.read_variant_calls(vcf, info_fields=["DP"], format_fields=["GT"], verbose=False)
    assert result["variant_count"] == 2
    first = result["variants"][0]
    assert (first["chrom"], first["start"], first["end"], first["id"]) == ("chr1", 9, 10, "rs1")
    assert first["DP"] == 12
    assert first["genotypes"] == {"GT": ["0/1", "1/1"]}
    bcf_path = tmp_path / "v.bcf"
    source = pysam.VariantFile(str(vcf))
    target = pysam.VariantFile(str(bcf_path), "wb", header=source.header)
    for record in source:
        target.write(record)
    target.close()
    source.close()
    bcf = pb_api.read_variant_calls(
        bcf_path, file_format="bcf", info_fields=["DP"], format_fields=["GT"], verbose=False
    )
    assert [(row["start"], row["end"]) for row in bcf["variants"]] == [(9, 10), (24, 25)]


def test_write_vcf_variants_round_trip_and_no_overwrite(tmp_path):
    pytest.importorskip("polars_bio")
    source = _write(tmp_path, "v.vcf", VCF)
    target = tmp_path / "copy.vcf"
    result = pb_api.write_vcf_variants(
        source, target, info_fields=["DP"], format_fields=["GT"], verbose=False
    )
    assert result["records_written"] == 2
    again = pb_api.read_variant_calls(
        target, info_fields=["DP"], format_fields=["GT"], verbose=False
    )
    original = pb_api.read_variant_calls(
        source, info_fields=["DP"], format_fields=["GT"], verbose=False
    )
    assert [(row["start"], row["ref"], row["alt"]) for row in again["variants"]] == [
        (row["start"], row["ref"], row["alt"]) for row in original["variants"]
    ]
    with pytest.raises(ConfigurationError, match="already exists"):
        pb_api.write_vcf_variants(source, target, verbose=False)


# polars-bio sequence files.


def test_read_fasta_sequences_lengths_and_sequences(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "s.fa", FASTA)
    result = pb_api.read_fasta_sequences(path, verbose=False)
    sequences = {row["name"]: (row["sequence"], row["length"]) for row in result["sequences"]}
    assert sequences == {"chr1": ("ACGTACGTACGGTT", 14), "chr2": ("TTTTAAAACCCC", 12)}


def test_read_fastq_reads_keeps_quality_strings(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "r.fq", FASTQ)
    result = pb_api.read_fastq_reads(path, verbose=False)
    reads = {row["name"]: row for row in result["reads"]}
    assert reads["r1"]["sequence"] == "ACGTN"
    assert reads["r1"]["quality_scores"] == "IIII#"
    assert reads["r2"]["length"] == 6
    assert result["read_count"] == 2


def test_write_fasta_sequences_writes_uppercase_and_rejects_bad_letters(tmp_path):
    pb = pytest.importorskip("polars_bio")
    target = tmp_path / "out.fa"
    result = pb_api.write_fasta_sequences({"x1": "acgtn", "x2": "GGCC"}, target, verbose=False)
    assert result["records_written"] == 2
    back = pb.read_fasta(str(target)).to_dicts()
    assert {row["name"]: row["sequence"] for row in back} == {"x1": "ACGTN", "x2": "GGCC"}
    with pytest.raises(ConfigurationError):
        pb_api.write_fasta_sequences({"bad": "ACGX1"}, tmp_path / "bad.fa", verbose=False)


# polars-bio alignments, pileups, quality, pairs, matrices, and alignments of sequences.


def test_read_alignment_records_sam_and_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    sam = _write(tmp_path, "a.sam", SAM)
    result = pb_api.read_alignment_records(sam, verbose=False)
    records = {row["name"]: row for row in result["alignments"]}
    assert (records["r1"]["start"], records["r1"]["end"], records["r1"]["cigar"]) == (10, 15, "5M")
    assert (records["r2"]["start"], records["r2"]["end"]) == (12, 16)
    bam_path = tmp_path / "a.bam"
    source = pysam.AlignmentFile(str(sam))
    target = pysam.AlignmentFile(str(bam_path), "wb", header=source.header)
    for read in source:
        target.write(read)
    target.close()
    source.close()
    bam = pb_api.read_alignment_records(bam_path, file_format="bam", verbose=False)
    assert [(row["start"], row["end"]) for row in bam["alignments"]] == [(10, 15), (12, 16)]


def test_alignment_depth_runs_match_hand_calculation(tmp_path):
    pytest.importorskip("polars_bio")
    sam = _write(tmp_path, "a.sam", SAM)
    result = pb_api.alignment_depth(sam, verbose=False)
    runs = [(row["start"], row["end"], row["depth"]) for row in result["runs"]]
    assert runs == [(10, 12, 1), (12, 15, 2), (15, 16, 1)]


def test_fastqc_summary_basic_statistics(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "r.fq", FASTQ)
    result = pb_api.fastqc_summary(path, verbose=False)
    stats = result["basic_statistics"]
    assert stats["n_seq"] == 2
    assert stats["total_bases"] == 11
    assert (stats["min_len"], stats["max_len"]) == (5, 6)
    # Six G or C among the ten A, C, G, and T bases; the N base is excluded.
    assert stats["gc_pct"] == pytest.approx(60.0)
    assert "basic_stats" in result["module_status"]


def test_read_contact_pairs_zero_based_positions(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "c.pairs", PAIRS)
    result = pb_api.read_contact_pairs(path, verbose=False)
    pairs = {row["readID"]: (row["pos1"], row["pos2"]) for row in result["pairs"]}
    assert pairs == {"r1": (99, 499), "r2": (199, 39)}


def test_read_bed_intervals_keeps_names(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "a.bed", A_BED)
    result = pb_api.read_bed_intervals(path, verbose=False)
    assert result["interval_count"] == 4
    assert {row["name"]: (row["start"], row["end"]) for row in result["intervals"]}["b"] == (
        150,
        300,
    )


def test_read_cooler_matrix_totals(tmp_path):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    bins = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 100], "end": [100, 200]})
    pixels = pd.DataFrame({"bin1_id": [0, 0, 1], "bin2_id": [0, 1, 1], "count": [3, 2, 5]})
    path = tmp_path / "t.cool"
    cooler.create_cooler(str(path), bins, pixels, ordered=False)
    result = pb_api.read_cooler_matrix(path, verbose=False)
    assert result["contact_count"] == 3
    assert result["total_counts"] == 10
    assert (result["contacts"][0]["start1"], result["contacts"][0]["end1"]) == (0, 100)


def test_read_stockholm_alignments_sequences(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "m.sto", STO)
    result = pb_api.read_stockholm_alignments(path, verbose=False)
    sequences = {row["name"]: row["sequence"] for row in result["alignments"]}
    assert sequences == {"seq1": "ACGT-ACGT", "seq2": "ACGTTACGA"}
    assert result["alignment_columns"] == 9


def test_polars_bio_row_limit_is_enforced(tmp_path):
    pytest.importorskip("polars_bio")
    path = _write(tmp_path, "a.bed", A_BED)
    with pytest.raises(ConfigurationError, match="max_rows"):
        pb_api.read_bed_intervals(path, max_rows=2, verbose=False)


# gtars genome sequences and region sets.


def test_region_gc_content_matches_hand_fractions(tmp_path):
    pytest.importorskip("gtars")
    regions = _write(tmp_path, "r.bed", "chr2\t0\t4\nchr2\t8\t12\nchr1\t10\t14\n")
    result = gt_api.region_gc_content(regions, _genome_fasta(tmp_path), verbose=False)
    gc = {(row["chrom"], row["start"]): row["gc_fraction"] for row in result["regions"]}
    assert gc == {("chr1", 10): 0.5, ("chr2", 0): 0.0, ("chr2", 8): 1.0}
    assert result["mean_gc_fraction"] == pytest.approx(0.5)


def test_region_dinucleotide_frequencies_counts(tmp_path):
    gtars = pytest.importorskip("gtars")
    regions = _write(tmp_path, "r.bed", "chr1\t0\t10\nchr1\t2\t6\n")
    fasta = _genome_fasta(tmp_path)
    result = gt_api.region_dinucleotide_frequencies(regions, fasta, raw_counts=True, verbose=False)
    first = {row["label"]: row["frequencies"] for row in result["regions"]}["chr1_0_10"]
    # ACGTACGTAC has AC three times, CG twice, GT twice, and TA twice.
    assert (first["AC"], first["CG"], first["GT"], first["TA"]) == (3.0, 2.0, 2.0, 2.0)
    assert sum(first.values()) == 9.0
    assert result["dinucleotides"][1] == "AC"
    models = gtars.models
    direct = gtars.genomic_distributions.calc_dinucl_freq(
        models.RegionSet(str(regions)), models.GenomeAssembly(str(fasta)), raw_counts=True
    )
    assert direct["frequencies"][0] == [first[name] for name in result["dinucleotides"]]


def test_partition_region_counts_observed_and_expected(tmp_path):
    pytest.importorskip("gtars")
    regions = _write(tmp_path, "r.bed", REGIONS)
    genes = _write(tmp_path, "genes.gtf", GENES)
    observed = gt_api.partition_region_counts(regions, genes, verbose=False)
    counts = {row["partition"]: row["observed"] for row in observed["partitions"]}
    # p1 and p3 overlap the core promoters of G1 and G2, p2 lies in an intron, p4 is intergenic.
    assert counts == {
        "promoterCore": 2,
        "promoterProx": 0,
        "exon": 0,
        "intron": 1,
        "intergenic": 1,
    }
    assert observed["total"] == 4
    expected = gt_api.partition_region_counts(
        regions, genes, chrom_sizes={"chr1": 10000}, verbose=False
    )
    assert expected["mode"] == "expected"
    assert sum(row["observed"] for row in expected["partitions"]) == 4
    assert all(row["p_value"] is None or 0 <= row["p_value"] <= 1 for row in expected["partitions"])


def test_signal_matrix_summary_uses_overlapping_labels(tmp_path):
    pytest.importorskip("gtars")
    regions = _write(tmp_path, "u.bed", "chr1\t0\t10\nchr1\t20\t30\n")
    signal = _write(tmp_path, "sig.tsv", "region\tcondA\tcondB\nchr1_0_100\t1.0\t2.0\n")
    result = gt_api.signal_matrix_summary(regions, signal, verbose=False)
    assert result["condition_names"] == ["condA", "condB"]
    values = {row["label"]: row["values"] for row in result["regions"]}
    assert values == {
        "chr1_0_10": {"condA": 1.0, "condB": 2.0},
        "chr1_20_30": {"condA": 1.0, "condB": 2.0},
    }


def test_consensus_regions_counts_files_per_region(tmp_path):
    pytest.importorskip("gtars")
    first = _write(tmp_path, "u.bed", "chr1\t0\t10\nchr1\t20\t30\n")
    second = _write(tmp_path, "u2.bed", "chr1\t0\t10\nchr1\t20\t30\n")
    third = _write(tmp_path, "r.bed", "chr1\t20\t30\nchr2\t5\t9\n")
    result = gt_api.consensus_regions([first, second, third], verbose=False)
    counts = {(row["chrom"], row["start"], row["end"]): row["count"] for row in result["regions"]}
    assert counts == {("chr1", 0, 10): 2, ("chr1", 20, 30): 3, ("chr2", 5, 9): 1}
    with pytest.raises(ConfigurationError, match="At least two"):
        gt_api.consensus_regions([first], verbose=False)


def test_lola_enrichment_matches_direct_upstream(tmp_path):
    gtars = pytest.importorskip("gtars")
    user = _write(tmp_path, "user.bed", "chr1\t0\t10\nchr1\t20\t30\n")
    universe = _write(tmp_path, "universe.bed", "chr1\t0\t100\nchr1\t100\t200\nchr2\t0\t50\n")
    db_a = _write(tmp_path, "dbA.bed", "chr1\t0\t10\n")
    db_b = _write(tmp_path, "dbB.bed", "chr1\t20\t30\nchr2\t0\t50\n")
    result = gt_api.lola_enrichment(user, universe, [db_a, db_b], verbose=False)
    assert result["database_count"] == 2
    lola = gtars.lola
    direct = lola.run_lola(
        [[("chr1", 0, 10), ("chr1", 20, 30)]],
        [("chr1", 0, 100), ("chr1", 100, 200), ("chr2", 0, 50)],
        lola.RegionDB.from_bed_files([str(db_a), str(db_b)], ["dbA.bed", "dbB.bed"]),
    )
    got = [row["p_value_log10"] for row in result["results"]]
    assert got == pytest.approx([float(value) for value in direct["pValueLog"]])
    assert [row["db_name"] for row in result["results"]] == ["dbA.bed", "dbB.bed"]


def test_tokenize_fragments_maps_records_to_universe(tmp_path):
    pytest.importorskip("gtars")
    fragments = _write(
        tmp_path, "frags.tsv", "chr1\t10\t20\tAAA\t1\nchr1\t150\t160\tAAA\t2\nchr2\t5\t9\tBBB\t3\n"
    )
    universe = _write(tmp_path, "univ.bed", "chr1\t0\t100\nchr1\t100\t200\nchr2\t0\t50\n")
    result = gt_api.tokenize_fragments(fragments, universe, verbose=False)
    tokens = {row["barcode"]: row["token_ids"] for row in result["barcodes"]}
    assert tokens == {"AAA": [0, 1], "BBB": [2]}
    with pytest.raises(ConfigurationError, match="BED columns"):
        bad = _write(tmp_path, "bad.bed", "chr1\t0\t100\tname\n")
        gt_api.tokenize_fragments(fragments, bad, verbose=False)


# gtars sequences, digests, and identifiers.


def test_digest_fasta_sequences_matches_direct_digests(tmp_path):
    gtars = pytest.importorskip("gtars")
    fasta = _write(tmp_path, "s.fa", FASTA)
    result = gt_api.digest_fasta_sequences(fasta, verbose=False)
    by_name = {row["name"]: row for row in result["sequences"]}
    assert by_name["chr1"]["length"] == 14
    assert by_name["chr1"]["sha512t24u"] == gtars.refget.sha512t24u_digest("ACGTACGTACGGTT")
    assert by_name["chr2"]["md5"] == gtars.refget.md5_digest("TTTTAAAACCCC")
    assert result["collection_digest"] == str(gtars.refget.digest_fasta(str(fasta)).digest)


def test_fasta_index_records_lengths(tmp_path):
    gtars = pytest.importorskip("gtars")
    fasta = _write(tmp_path, "s.fa", FASTA)
    result = gt_api.fasta_index_records(fasta, verbose=False)
    lengths = {row["name"]: row["length"] for row in result["records"]}
    assert lengths == {"chr1": 14, "chr2": 12}
    direct = {item.name: item.length for item in gtars.refget.compute_fai(str(fasta))}
    assert lengths == direct


def test_sequence_digest_algorithms_and_case(tmp_path):
    gtars = pytest.importorskip("gtars")
    sha = gt_api.sequence_digest("acgt", verbose=False)
    assert sha["digest"] == gtars.refget.sha512t24u_digest("ACGT")
    assert sha["length"] == 4
    md5 = gt_api.sequence_digest("ACGT", algorithm="md5", verbose=False)
    assert md5["digest"] == gtars.refget.md5_digest("ACGT")
    with pytest.raises(ConfigurationError):
        gt_api.sequence_digest("ACGT", algorithm="sha1", verbose=False)


def test_normalize_variant_allele_snv_and_deletion(tmp_path):
    gtars = pytest.importorskip("gtars")
    snv = gt_api.normalize_variant_allele("ACGTACGTAC", 2, "G", "T", verbose=False)
    assert (snv["start"], snv["end"], snv["allele"]) == (2, 3, "T")
    deletion = gt_api.normalize_variant_allele("ACGTTTTA", 2, "GT", "G", verbose=False)
    direct = gtars.vrs.normalize_allele("ACGTTTTA", 2, "GT", "G")
    assert (deletion["start"], deletion["end"], deletion["allele"]) == (
        direct["start"],
        direct["end"],
        direct["allele"],
    )
    with pytest.raises(ConfigurationError, match="does not match"):
        gt_api.normalize_variant_allele("ACGTACGTAC", 2, "A", "T", verbose=False)


def test_vrs_allele_id_matches_direct_upstream(tmp_path):
    gtars = pytest.importorskip("gtars")
    digest = gtars.refget.sha512t24u_digest("ACGTACGTACGGTT")
    prefixed = gt_api.vrs_allele_id(digest, 1, 2, "g", verbose=False)
    assert prefixed["identifier"] == gtars.vrs.vrs_id(digest, 1, 2, "G")
    assert prefixed["identifier"].startswith("ga4gh:VA.")
    bare = gt_api.vrs_allele_id(digest, 1, 2, "G", prefixed=False, verbose=False)
    assert bare["identifier"] == gtars.vrs.vrs_digest(digest, 1, 2, "G")


def test_parse_hgvs_variant_components(tmp_path):
    pytest.importorskip("gtars")
    result = gt_api.parse_hgvs_variant("NC_000001.11:g.12345A>G", verbose=False)
    assert result["accession"] == "NC_000001.11"
    assert result["reference_type"] == "g"
    edit = result["variant"]["pos_edit"]["edit"]
    assert (edit["kind"], edit["ref"], edit["alt"]) == ("substitution", "A", "G")
    assert result["variant"]["pos_edit"]["start"]["position"]["base"] == 12345
    with pytest.raises(ConfigurationError):
        gt_api.parse_hgvs_variant("not an hgvs string", verbose=False)


# genomicranges GTF ranges, windows, and tiling.


def test_read_gtf_ranges_restores_gtf_ends(tmp_path):
    pytest.importorskip("genomicranges")
    path = _write(tmp_path, "genes.gtf", GTF)
    result = gr_api.read_gtf_ranges(path, verbose=False)
    spans = {
        (row["feature"], row["gene_id"], row["start"], row["end"], row["strand"])
        for row in result["features"]
    }
    assert ("gene", "G1", 10, 100, "+") in spans
    assert ("exon", "G1", 10, 50, "+") in spans
    assert ("gene", "G2", 199, 300, "-") in spans
    assert result["feature_count"] == 3


def test_upstream_downstream_regions_tss_and_tes(tmp_path):
    pytest.importorskip("genomicranges")
    path = _write(tmp_path, "r.bed", REGIONS)
    tss = gr_api.upstream_downstream_regions(path, upstream=100, downstream=50, verbose=False)
    windows = {(row["source_start"]): (row["start"], row["end"]) for row in tss["windows"]}
    # Start 950 (0-based) is the TSS on the default strand: [TSS-100, TSS+50).
    assert windows[950] == (850, 1000)
    assert windows[9000] == (8900, 9050)
    tes = gr_api.upstream_downstream_regions(
        path, upstream=100, downstream=50, site="TES", verbose=False
    )
    tes_windows = {row["source_start"]: (row["start"], row["end"]) for row in tes["windows"]}
    # End 1010 (exclusive) is the last base at 1-based 1010, so the window is 1-based 910..1059.
    assert tes_windows[950] == (909, 1059)


def test_upstream_downstream_regions_minus_strand(tmp_path):
    pytest.importorskip("genomicranges")
    path = _write(tmp_path, "minus.bed", "chr1\t500\t600\tm\t0\t-\n")
    result = gr_api.upstream_downstream_regions(path, upstream=100, downstream=50, verbose=False)
    window = result["windows"][0]
    # On the minus strand the TSS is the exclusive end 600: [TSS-50, TSS+100).
    assert (window["start"], window["end"], window["strand"]) == (550, 700, "-")


def test_split_interval_bins_covers_range_once():
    pytest.importorskip("genomicranges")
    result = gr_api.split_interval_bins(0, 100, 30, verbose=False)
    bins = [(row["start"], row["end"]) for row in result["bins"]]
    assert bins == [(0, 30), (30, 60), (60, 90), (90, 100)]
    with pytest.raises(ConfigurationError):
        gr_api.split_interval_bins(10, 10, 5, verbose=False)


def test_slide_interval_windows_moves_by_step():
    pytest.importorskip("genomicranges")
    result = gr_api.slide_interval_windows(0, 100, 30, 10, verbose=False)
    starts = [row["start"] for row in result["windows"]]
    assert starts == [0, 10, 20, 30, 40, 50, 60, 70]
    assert all(row["end"] - row["start"] == 30 for row in result["windows"])
    short = gr_api.slide_interval_windows(0, 20, 30, 10, verbose=False)
    assert short["windows"] == [{"index": 0, "start": 0, "end": 20}]
