"""Unit tests for the G18a bioframe review wrappers.

Each wrapped function is checked against a hand calculation or a direct upstream call. Fixtures are
written inside tmp_path. Optional packages skip their tests through pytest.importorskip; the
bigBed test uses the test.bigBed file shipped with the pyBigWigTest package when it is present.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from dnakit.annotation import bioframe_review_functions as review
from dnakit.exceptions import ConfigurationError

BED6_A = (
    "chr1\t0\t100\ta\t5\t+\nchr1\t50\t150\tb\t3\t+\n"
    "chr1\t300\t400\tc\t1\t-\nchr2\t10\t60\td\t2\t+\n"
)
VIEW = "chr1\t0\t200\tv1\nchr2\t0\t100\tv2\n"
BED4_B = "chr1\t80\t120\tx\nchr2\t0\t60\ty\n"
TILE = "chr1\t0\t100\nchr1\t100\t200\n"
FASTA_CHR1 = "ACGTNNGGCCGATTACAGAATTCGC"
FASTA_CHR2 = "GGGGCCCCAAAATTTTGCGCGCAT"
FASTA = f">chr1\n{FASTA_CHR1}\n>chr2\n{FASTA_CHR2}\n"
SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:1000\n"
    "r1\t0\tchr1\t11\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tNM:i:0\n"
    "r2\t16\tchr1\t101\t30\t4M2D4M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
)
MRNA = (
    "matches\tmisMatches\trepMatches\tnCount\tqNumInsert\tqBaseInsert\ttNumInsert\t"
    "tBaseInsert\tstrand\tqName\tqSize\tqStart\tqEnd\ttName\ttSize\ttStart\ttEnd\t"
    "blockCount\tblockSizes\tqStarts\ttStarts\n"
    "100\t0\t0\t0\t0\t0\t0\t0\t+\tNM1\t500\t0\t100\tchr1\t1000\t20\t120\t1\t100,\t0,\t20,\n"
)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _spans(rows, *keys):
    return [tuple(row[key] for key in keys) for row in rows]


@pytest.fixture
def bioframe():
    return pytest.importorskip("bioframe")


def test_every_public_function_takes_keyword_only_verbose():
    for name in review.__all__:
        parameter = inspect.signature(getattr(review, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_validate_bioframe_bed(tmp_path, bioframe):
    good = _write(tmp_path / "good.bed", BED6_A)
    result = review.validate_bioframe_bed(good, verbose=False)
    assert result["passed"] is True
    assert result["record_count"] == 4
    assert result["provenance"]["package"].startswith("bioframe")

    reversed_ = _write(tmp_path / "rev.bed", "chr1\t20\t10\n")
    assert review.validate_bioframe_bed(reversed_, verbose=False)["passed"] is False

    two_columns = _write(tmp_path / "two.bed", "chr1\t10\n")
    with pytest.raises(ConfigurationError):
        review.validate_bioframe_bed(two_columns, verbose=False)


def test_check_bioframe_cataloged(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", VIEW)
    named_a = _write(tmp_path / "a.bed", BED6_A)
    assert review.check_bioframe_cataloged(named_a, view, verbose=False)["passed"] is False
    labelled = _write(tmp_path / "lab.bed", "chr1\t0\t10\tv1\nchr2\t0\t5\tv2\n")
    assert review.check_bioframe_cataloged(labelled, view, verbose=False)["passed"] is True


def test_check_bioframe_contained(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", VIEW)
    inside = _write(tmp_path / "in.bed", "chr1\t10\t20\nchr2\t0\t50\n")
    assert review.check_bioframe_contained(inside, view, verbose=False)["passed"] is True
    outside = _write(tmp_path / "out.bed", BED6_A)
    assert review.check_bioframe_contained(outside, view, verbose=False)["passed"] is False


def test_check_bioframe_covering_and_tiling(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", "chr1\t0\t200\tv\n")
    tiles = _write(tmp_path / "tiles.bed", TILE)
    assert review.check_bioframe_covering(tiles, view, verbose=False)["passed"] is True
    assert review.check_bioframe_tiling(tiles, view, verbose=False)["passed"] is True
    half = _write(tmp_path / "half.bed", "chr1\t0\t100\n")
    assert review.check_bioframe_covering(half, view, verbose=False)["passed"] is False
    assert review.check_bioframe_tiling(half, view, verbose=False)["passed"] is False


def test_check_bioframe_overlapping_and_sorted(tmp_path, bioframe):
    overlap = _write(tmp_path / "ov.bed", "chr1\t0\t10\nchr1\t5\t20\n")
    touch = _write(tmp_path / "touch.bed", "chr1\t0\t10\nchr1\t10\t20\n")
    assert review.check_bioframe_overlapping(overlap, verbose=False)["passed"] is True
    assert review.check_bioframe_overlapping(touch, verbose=False)["passed"] is False
    sorted_file = _write(tmp_path / "sorted.bed", "chr1\t0\t10\nchr1\t5\t20\nchr2\t0\t1\n")
    unsorted = _write(tmp_path / "unsorted.bed", "chr1\t5\t20\nchr1\t0\t10\n")
    assert review.check_bioframe_sorted(sorted_file, verbose=False)["passed"] is True
    assert review.check_bioframe_sorted(unsorted, verbose=False)["passed"] is False


def test_check_bioframe_viewframe(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", VIEW)
    bed3 = _write(tmp_path / "bed3.bed", "chr1\t0\t100\n")
    assert review.check_bioframe_viewframe(view, verbose=False)["passed"] is True
    assert review.check_bioframe_viewframe(bed3, verbose=False)["passed"] is False


def test_resolve_bioframe_region(bioframe):
    span = review.resolve_bioframe_region("chr1:100-200", verbose=False)
    assert (span["chrom"], span["start"], span["end"], span["length"]) == ("chr1", 100, 200, 100)
    commas = review.resolve_bioframe_region("chr1:1,000-2,000", verbose=False)
    assert (commas["start"], commas["end"]) == (1000, 2000)
    whole = review.resolve_bioframe_region("chr1", chromsizes={"chr1": 1000}, verbose=False)
    assert (whole["start"], whole["end"]) == (0, 1000)
    with pytest.raises(ConfigurationError):
        review.resolve_bioframe_region("chr1", verbose=False)
    with pytest.raises(ConfigurationError):
        review.resolve_bioframe_region("chr1:5-2000", chromsizes={"chr1": 1000}, verbose=False)


def test_format_bioframe_region(bioframe):
    result = review.format_bioframe_region("chr1", 0, 100, verbose=False)
    assert result["region"] == "chr1:0-100"
    with pytest.raises(ConfigurationError):
        review.format_bioframe_region("chr1", 100, 100, verbose=False)
    with pytest.raises(ConfigurationError):
        review.format_bioframe_region("chr1", -1, 100, verbose=False)


def test_sanitize_bioframe_intervals(tmp_path, bioframe):
    bad = _write(tmp_path / "bad.bed", "chr1\t20\t10\tn\nchr1\t0\t5\tm\n")
    with pytest.raises(ConfigurationError):
        review.sanitize_bioframe_intervals(bad, verbose=False)
    flipped = review.sanitize_bioframe_intervals(bad, start_exceed_end_action="flip", verbose=False)
    assert _spans(flipped["intervals"], "start", "end") == [(10, 20), (0, 5)]
    dropped = review.sanitize_bioframe_intervals(bad, start_exceed_end_action="drop", verbose=False)
    assert dropped["interval_count"] == 1
    with pytest.raises(ConfigurationError):
        review.sanitize_bioframe_intervals(bad, start_exceed_end_action="swap", verbose=False)


def test_sort_bioframe_intervals(tmp_path, bioframe):
    path = _write(tmp_path / "s.bed", "chr10\t5\t9\nchr2\t1\t3\nchr1\t7\t8\nchr2\t0\t1\n")
    result = review.sort_bioframe_intervals(path, verbose=False)
    assert _spans(result["intervals"], "chrom", "start") == [
        ("chr1", 7),
        ("chr10", 5),
        ("chr2", 0),
        ("chr2", 1),
    ]


def test_setdiff_bioframe_intervals(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    b = _write(tmp_path / "b.bed", BED4_B)
    result = review.setdiff_bioframe_intervals(a, b, verbose=False)
    assert _spans(result["intervals"], "chrom", "start", "end") == [("chr1", 300, 400)]


def test_trim_bioframe_intervals(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", VIEW)
    kept = _write(tmp_path / "a.bed", BED6_A)
    result = review.trim_bioframe_intervals(kept, view, verbose=False)
    assert result["interval_count"] == 3
    assert result["dropped_count"] == 1
    clip = _write(tmp_path / "clip.bed", "chr1\t150\t250\tz\nchr1\t300\t400\tc\n")
    clipped = review.trim_bioframe_intervals(clip, view, verbose=False)
    assert _spans(clipped["intervals"], "start", "end") == [(150, 200)]
    assert clipped["dropped_count"] == 1


def test_complement_bioframe_intervals(tmp_path, bioframe):
    view = _write(tmp_path / "view.bed", VIEW)
    b = _write(tmp_path / "b.bed", BED4_B)
    result = review.complement_bioframe_intervals(b, view, verbose=False)
    gaps = sorted(_spans(result["gaps"], "chrom", "start", "end"))
    assert gaps == [("chr1", 0, 80), ("chr1", 120, 200), ("chr2", 60, 100)]


def test_count_bioframe_overlaps(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    b = _write(tmp_path / "b.bed", BED4_B)
    result = review.count_bioframe_overlaps(a, b, verbose=False)
    assert [row["count"] for row in result["intervals"]] == [1, 1, 0, 1]


def test_select_bioframe_region(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    middle = review.select_bioframe_region(a, "chr1:90-120", verbose=False)
    assert [row["name"] for row in middle["intervals"]] == ["a", "b"]
    chrom = review.select_bioframe_region(a, "chr2", verbose=False)
    assert [row["name"] for row in chrom["intervals"]] == ["d"]


def test_assign_bioframe_views(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    view = _write(tmp_path / "view.bed", VIEW)
    result = review.assign_bioframe_views(a, view, verbose=False)
    assigned = [row["view_region"] for row in result["intervals"]]
    assert assigned == ["v1", "v1", None, "v2"]
    dropped = review.assign_bioframe_views(a, view, drop_unassigned=True, verbose=False)
    assert dropped["interval_count"] == 3


def test_mark_and_merge_bioframe_runs(tmp_path, bioframe):
    runs = _write(tmp_path / "runs.bed", "chr1\t0\t100\tx\nchr1\t100\t200\tx\nchr1\t200\t300\ty\n")
    marked = review.mark_bioframe_runs(runs, "name", verbose=False)
    assert [row["run"] for row in marked["intervals"]] == [0, 0, 1]
    merged = review.merge_bioframe_runs(runs, "name", verbose=False)
    assert _spans(merged["intervals"], "start", "end") == [(0, 200), (200, 300)]
    overlapping = _write(tmp_path / "ov.bed", BED6_A)
    with pytest.raises(ConfigurationError):
        review.mark_bioframe_runs(overlapping, "name", verbose=False)
    with pytest.raises(ConfigurationError):
        review.merge_bioframe_runs(runs, "score", verbose=False)


def test_pair_bioframe_intervals_by_distance(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    # Midpoints on chr1: a=50, b=100, c=350, so the pair distances are 50, 250 and 300.
    near = review.pair_bioframe_intervals_by_distance(a, min_sep=0, max_sep=280, verbose=False)
    assert near["pair_count"] == 2
    far = review.pair_bioframe_intervals_by_distance(a, min_sep=260, max_sep=400, verbose=False)
    assert far["pair_count"] == 1
    assert (far["pairs"][0]["start_1"], far["pairs"][0]["start_2"]) == (0, 300)
    with pytest.raises(ConfigurationError):
        review.pair_bioframe_intervals_by_distance(a, min_sep=5, max_sep=5, verbose=False)


def test_sequence_gc_fraction():
    assert review.sequence_gc_fraction("GGCCAATT", verbose=False)["gc_fraction"] == 0.5
    assert review.sequence_gc_fraction("GGNN", verbose=False)["gc_fraction"] == 1.0
    unmapped = review.sequence_gc_fraction("GGNN", mapped_only=False, verbose=False)
    assert unmapped["gc_fraction"] == 0.5
    assert review.sequence_gc_fraction("NNNN", verbose=False)["gc_fraction"] is None
    with pytest.raises(ConfigurationError):
        review.sequence_gc_fraction("ACGZ", verbose=False)


def test_frac_gc_and_frac_mapped(tmp_path, bioframe):
    fasta = _write(tmp_path / "g.fa", FASTA)
    regions = _write(tmp_path / "small.bed", "chr1\t0\t10\tr1\nchr2\t0\t12\tr3\n")
    gc = review.frac_gc_bioframe_intervals(regions, fasta, verbose=False)
    # chr1[0:10] = ACGTNNGGCC has six G or C among eight non-N bases.
    assert gc["intervals"][0]["frac_gc"] == pytest.approx(6 / 8)
    # chr2[0:12] = GGGGCCCCAAAA has eight G or C among twelve bases.
    assert gc["intervals"][1]["frac_gc"] == pytest.approx(8 / 12)
    unmapped = review.frac_gc_bioframe_intervals(regions, fasta, mapped_only=False, verbose=False)
    assert unmapped["intervals"][0]["frac_gc"] == pytest.approx(6 / 10)

    mapped = review.frac_mapped_bioframe_intervals(regions, fasta, verbose=False)
    assert mapped["intervals"][0]["frac_mapped"] == pytest.approx(8 / 10)
    outside = _write(tmp_path / "out.bed", "chr1\t30\t40\tq\n")
    assert review.frac_mapped_bioframe_intervals(outside, fasta, verbose=False)["intervals"][0][
        "frac_mapped"
    ] == pytest.approx(0.0)


def test_digest_bioframe_fasta(tmp_path, bioframe):
    pytest.importorskip("Bio.Restriction")
    fasta = _write(tmp_path / "g.fa", FASTA)
    result = review.digest_bioframe_fasta(fasta, "EcoRI", verbose=False)
    site = FASTA_CHR1.find("GAATTC")
    chr1 = [row for row in result["fragments"] if row["chrom"] == "chr1"]
    assert _spans(chr1, "start", "end") == [(0, site + 1), (site + 1, len(FASTA_CHR1))]
    # Biopython reports the 1-based base after the cut; the 0-based boundary is one less.
    biopython = pytest.importorskip("Bio.Restriction")
    from Bio.Seq import Seq

    assert biopython.EcoRI.search(Seq(FASTA_CHR1)) == [site + 2]
    chr2 = [row for row in result["fragments"] if row["chrom"] == "chr2"]
    assert _spans(chr2, "start", "end") == [(0, len(FASTA_CHR2))]
    with pytest.raises(ConfigurationError):
        review.digest_bioframe_fasta(fasta, "NoSuchEnzyme", verbose=False)


def test_summarize_bioframe_fasta(tmp_path, bioframe):
    fasta = _write(tmp_path / "g.fa", FASTA)
    result = review.summarize_bioframe_fasta(fasta, verbose=False)
    assert [(row["name"], row["length"]) for row in result["sequences"]] == [
        ("chr1", len(FASTA_CHR1)),
        ("chr2", len(FASTA_CHR2)),
    ]
    assert result["total_length"] == len(FASTA_CHR1) + len(FASTA_CHR2)


def test_list_and_describe_assemblies(bioframe):
    listing = review.list_bioframe_assemblies(verbose=False)
    assert listing["assembly_count"] > 10
    builds = {(row["provider"], row["provider_build"]) for row in listing["assemblies"]}
    assert ("ucsc", "hg38") in builds
    hg38 = review.describe_bioframe_assembly("ucsc.hg38", verbose=False)
    assert hg38["organism"] == "homo sapiens"
    chr1 = next(row for row in hg38["sequences"] if row["name"] == "chr1")
    assert chr1["length"] == 248956422
    with pytest.raises(ConfigurationError):
        review.describe_bioframe_assembly("no-such-build", verbose=False)
    with pytest.raises(ConfigurationError):
        review.describe_bioframe_assembly("ucsc.hg38", roles="assembled", verbose=False)


def test_make_bioframe_chromarms():
    arms = review.make_bioframe_chromarms(
        {"chr1": 1000, "chr2": 500}, {"chr1": 400, "chr2": 200}, verbose=False
    )
    assert _spans(arms["arms"], "chrom", "start", "end", "name") == [
        ("chr1", 0, 400, "chr1_p"),
        ("chr1", 400, 1000, "chr1_q"),
        ("chr2", 0, 200, "chr2_p"),
        ("chr2", 200, 500, "chr2_q"),
    ]
    with pytest.raises(ConfigurationError):
        review.make_bioframe_chromarms({"chr1": 1000, "chr2": 500}, {"chr1": 400}, verbose=False)
    with pytest.raises(ConfigurationError):
        review.make_bioframe_chromarms({"chr1": 1000}, {"chr1": 1000}, verbose=False)


def test_read_bioframe_chromsizes(tmp_path, bioframe):
    path = _write(tmp_path / "sizes.txt", "chr1\t1000\nchr2\t500\nchrUn_x\t10\n")
    filtered = review.read_bioframe_chromsizes(path, verbose=False)
    assert [row["chrom"] for row in filtered["chromsizes"]] == ["chr1", "chr2"]
    everything = review.read_bioframe_chromsizes(path, filter_chroms=False, verbose=False)
    assert everything["chromosome_count"] == 3
    assert everything["chromsizes"][0]["length"] == 1000


def test_read_bioframe_tabix(tmp_path, bioframe):
    pysam = pytest.importorskip("pysam")
    plain = tmp_path / "tab.bed"
    plain.write_text("chr1\t10\t20\tn\nchr1\t120\t130\tm\nchr2\t90\t99\tk\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        review.read_bioframe_tabix(plain, verbose=False)
    packed = tmp_path / "tab.bed.gz"
    pysam.tabix_compress(str(plain), str(packed), force=True)
    pysam.tabix_index(str(packed), preset="bed", force=True)
    whole = review.read_bioframe_tabix(packed, verbose=False)
    assert _spans(whole["intervals"], "chrom", "start", "end") == [
        ("chr1", 10, 20),
        ("chr1", 120, 130),
        ("chr2", 90, 99),
    ]
    region = review.read_bioframe_tabix(packed, chrom="chr1", start=0, end=100, verbose=False)
    assert [row["name"] for row in region["intervals"]] == ["n"]
    with pytest.raises(ConfigurationError):
        review.read_bioframe_tabix(packed, start=0, end=100, verbose=False)


def test_read_bioframe_alignments(tmp_path, bioframe):
    pytest.importorskip("pysam")
    sam = _write(tmp_path / "x.sam", SAM)
    result = review.read_bioframe_alignments(sam, verbose=False)
    assert result["alignment_count"] == 2
    assert [row["QNAME"] for row in result["alignments"]] == ["r1", "r2"]
    assert [row["POS"] for row in result["alignments"]] == [10, 100]
    assert result["alignments"][0]["CIGAR"] == "8M"
    with pytest.raises(ConfigurationError):
        review.read_bioframe_alignments(sam, chrom="chr1", start=0, end=50, verbose=False)
    with pytest.raises(ConfigurationError):
        review.read_bioframe_alignments(sam, start=0, end=50, verbose=False)


def test_read_bioframe_bigbed(bioframe):
    pytest.importorskip("pyBigWig")
    fixture_package = pytest.importorskip("pyBigWigTest")
    path = Path(fixture_package.__file__).parent / "test.bigBed"
    if not path.is_file():
        pytest.skip("pyBigWigTest does not ship test.bigBed")
    chroms = sorted(bioframe.io.fileops.pyBigWig.open(str(path)).chroms())
    result = review.read_bioframe_bigbed(path, chroms[0], verbose=False)
    direct = bioframe.read_bigbed(str(path), chroms[0], engine="pybigwig")
    assert result["interval_count"] == len(direct) > 0
    assert all(row["start"] < row["end"] for row in result["intervals"])
    with pytest.raises(ConfigurationError):
        review.read_bioframe_bigbed(path, "chrZZZ_not_here", verbose=False)


def test_read_bioframe_table(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    auto = review.read_bioframe_table(a, verbose=False)
    assert auto["schema"] == "bed6"
    assert auto["row_count"] == 4
    with pytest.raises(ConfigurationError):
        review.read_bioframe_table(a, schema="bed4", verbose=False)
    commented = _write(tmp_path / "c.bed", "# header\nchr1\t1\t2\n")
    assert review.read_bioframe_table(commented, verbose=False)["row_count"] == 1
    with pytest.raises(ConfigurationError):
        review.read_bioframe_table(a, max_rows=2, verbose=False)


def test_write_bioframe_bed(tmp_path, bioframe):
    source = _write(tmp_path / "a.bed", BED6_A)
    target = tmp_path / "out.bed"
    result = review.write_bioframe_bed(source, target, verbose=False)
    assert result["record_count"] == 4
    assert target.read_text(encoding="utf-8") == BED6_A
    with pytest.raises(ConfigurationError):
        review.write_bioframe_bed(source, target, verbose=False)
    unsorted = _write(tmp_path / "u.bed", "chr1\t5\t10\nchr1\t0\t3\n")
    rejected = tmp_path / "rejected.bed"
    with pytest.raises(ConfigurationError):
        review.write_bioframe_bed(unsorted, rejected, require_sorted=True, verbose=False)
    assert not rejected.exists()
    with pytest.raises(ConfigurationError):
        review.write_bioframe_bed(source, tmp_path / "bad.bed", schema="bed7", verbose=False)


def test_read_pairix_is_not_wrapped_and_names_stay_internal():
    assert not hasattr(review, "read_bioframe_pairix")
    assert "parse_bioframe_region" not in review.__all__


def test_public_results_are_deterministic(tmp_path, bioframe):
    a = _write(tmp_path / "a.bed", BED6_A)
    first = review.pair_bioframe_intervals_by_distance(a, min_sep=0, max_sep=400, verbose=False)
    second = review.pair_bioframe_intervals_by_distance(a, min_sep=0, max_sep=400, verbose=False)
    assert first == second
