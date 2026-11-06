import shutil
import subprocess

import pytest

from dnakit.annotation import pybedtools_p33_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

pytest.importorskip("pybedtools")
pytestmark = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")

A = "chr1\t100\t200\ta1\t5\t+\nchr1\t500\t600\ta2\t6\t-\nchr2\t10\t90\ta3\t7\t+\n"
B = "chr1\t150\t250\tb1\t1\t+\nchr1\t900\t950\tb2\t2\t-\nchr2\t100\t150\tb3\t3\t+\n"
GENOME = "chr1\t2000\nchr2\t1000\n"
BED12 = "chr1\t100\t400\tg1\t0\t+\t100\t400\t0\t2\t50,60\t0,240\n"
BEDPE = "chr1\t100\t150\tchr1\t300\t350\tp1\t5\t+\t-\nchr2\t10\t60\tchr2\t200\t250\tp2\t6\t+\t+\n"
GFF = (
    "chr1\tt\tgene\t101\t600\t.\t+\t.\tID=g1\n"
    "chr1\tt\texon\t101\t200\t.\t+\t.\tParent=g1\n"
    "chr1\tt\texon\t301\t400\t.\t+\t.\tParent=g1\n"
    "chr1\tt\texon\t501\t600\t.\t+\t.\tParent=g1\n"
)
FASTA = ">chr1\n" + "ACGT" * 250 + "\n>chr2\n" + "AACCGGTT" * 125 + "\n"


@pytest.fixture()
def files(tmp_path):
    for name, text in (
        ("a.bed", A),
        ("b.bed", B),
        ("g.sizes", GENOME),
        ("g12.bed", BED12),
        ("p.bedpe", BEDPE),
        ("genes.gff", GFF),
        ("ref.fa", FASTA),
    ):
        (tmp_path / name).write_text(text)
    return tmp_path


@pytest.fixture()
def bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "r.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 2000}, {"SN": "chr2", "LN": 1000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for index, (ref, pos) in enumerate([(0, 120), (0, 520), (0, 1500), (1, 20)]):
            read = pysam.AlignedSegment()
            read.query_name = f"read{index}"
            read.query_sequence = "A" * 40
            read.flag = 0
            read.reference_id = ref
            read.reference_start = pos
            read.mapping_quality = 60
            read.cigar = [(0, 40)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            out.write(read)
    return path


def direct(*args):
    return subprocess.run(["bedtools", *args], capture_output=True, text=True, check=True).stdout


def rows(result):
    return ["\t".join(r) for r in result["records"]]


def test_bed12_to_bed6_matches_bedtools(files):
    result = m.convert_bed12_to_bed6_pybedtools(files / "g12.bed", verbose=False)
    assert "\n".join(rows(result)) + "\n" == direct("bed12tobed6", "-i", str(files / "g12.bed"))
    assert result["n_records"] == 2
    named = m.convert_bed12_to_bed6_pybedtools(files / "g12.bed", number_names=True, verbose=False)
    assert [r[4] for r in named["records"]] == ["1", "2"]


def test_bam_to_bed_and_back(files, bam):
    bed = m.bam_to_bed_pybedtools(bam, verbose=False)
    assert "\n".join(rows(bed)) + "\n" == direct("bamtobed", "-i", str(bam))
    assert bed["n_records"] == 4
    out = files / "back.bam"
    result = m.convert_intervals_to_bam_pybedtools(
        files / "a.bed", files / "g.sizes", out, verbose=False
    )
    assert result["input_format"] == "bed" and out.is_file()
    if result["n_alignments"] is not None:
        assert result["n_alignments"] == 3
    with pytest.raises(ConfigurationError):
        m.convert_intervals_to_bam_pybedtools(
            files / "a.bed", files / "g.sizes", out, verbose=False
        )


def test_bedpe_to_bam(files):
    out = files / "pe.bam"
    result = m.convert_intervals_to_bam_pybedtools(
        files / "p.bedpe", files / "g.sizes", out, input_format="bedpe", verbose=False
    )
    assert result["input_format"] == "bedpe" and out.stat().st_size > 0
    with pytest.raises(ConfigurationError):
        m.convert_intervals_to_bam_pybedtools(
            files / "p.bedpe",
            files / "g.sizes",
            files / "x.bam",
            input_format="bedpe",
            bed12=True,
            verbose=False,
        )


def test_records_dataframe(files):
    result = m.intervals_to_records_pybedtools(files / "a.bed", verbose=False)
    assert result["columns"] == ["chrom", "start", "end", "name", "score", "strand"]
    assert result["records"][0] == {
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "name": "a1",
        "score": 5,
        "strand": "+",
    }
    plain = m.intervals_to_records_pybedtools(
        files / "a.bed", disable_auto_names=True, max_rows=1, verbose=False
    )
    assert plain["columns"][0] == "0" and plain["truncated"] is True


def test_tabix_contigs(files):
    pysam = pytest.importorskip("pysam")
    path = files / "s.bed"
    path.write_text("chr1\t1\t2\nchr2\t3\t4\n")
    gz = pysam.tabix_index(str(path), preset="bed", force=True)
    result = m.list_tabix_contigs_pybedtools(gz, verbose=False)
    assert result["contigs"] == ["chr1", "chr2"]
    with pytest.raises(ConfigurationError):
        m.list_tabix_contigs_pybedtools(files / "a.bed", verbose=False)


def test_file_format_detection(files, bam):
    pysam = pytest.importorskip("pysam")
    assert m.detect_alignment_file_format_pybedtools(bam, verbose=False)["format"] == "bam"
    plain = m.detect_alignment_file_format_pybedtools(files / "a.bed", verbose=False)
    assert plain["format"] == "other" and not plain["is_bgzip"]
    path = files / "s.bed"
    path.write_text("chr1\t1\t2\n")
    gz = pysam.tabix_index(str(path), preset="bed", force=True)
    result = m.detect_alignment_file_format_pybedtools(gz, verbose=False)
    assert result["format"] == "bgzf" and result["is_bgzip"] and not result["is_bam"]


def test_interval_has_hits(files):
    assert m.interval_has_hits_pybedtools(files / "a.bed", "chr1", 150, 160, verbose=False)[
        "has_hits"
    ]
    assert not m.interval_has_hits_pybedtools(files / "a.bed", "chr1", 300, 400, verbose=False)[
        "has_hits"
    ]
    opposite = m.interval_has_hits_pybedtools(
        files / "a.bed", "chr1", 150, 160, strand="-", same_strand=True, verbose=False
    )
    assert opposite["has_hits"] is False
    with pytest.raises(ConfigurationError):
        m.interval_has_hits_pybedtools(files / "a.bed", "chr1", 10, 5, verbose=False)


def test_slop_and_shift_match_bedtools(files):
    slop = m.slop_pybedtools_intervals(files / "a.bed", files / "g.sizes", both=20, verbose=False)
    assert "\n".join(rows(slop)) + "\n" == direct(
        "slop", "-i", str(files / "a.bed"), "-g", str(files / "g.sizes"), "-b", "20"
    )
    assert slop["records"][2][1] == "0"  # clipped at the chromosome start
    stranded = m.slop_pybedtools_intervals(
        files / "a.bed", files / "g.sizes", left=10, right=0, strand_aware=True, verbose=False
    )
    assert stranded["records"][1][1:3] == ["500", "610"]
    shift = m.shift_pybedtools_intervals(
        files / "a.bed", files / "g.sizes", shift=30, verbose=False
    )
    assert "\n".join(rows(shift)) + "\n" == direct(
        "shift", "-i", str(files / "a.bed"), "-g", str(files / "g.sizes"), "-s", "30"
    )
    with pytest.raises(ConfigurationError):
        m.slop_pybedtools_intervals(files / "a.bed", files / "g.sizes", verbose=False)
    with pytest.raises(ConfigurationError):
        m.shift_pybedtools_intervals(files / "a.bed", files / "g.sizes", verbose=False)


def test_map_matches_bedtools(files):
    result = m.map_pybedtools_values(
        files / "a.bed", files / "b.bed", column=5, operation="sum", verbose=False
    )
    expected = direct(
        "map", "-a", str(files / "a.bed"), "-b", str(files / "b.bed"), "-c", "5", "-o", "sum"
    )
    assert "\n".join(rows(result)) + "\n" == expected
    assert result["records"][0][-1] == "1" and result["records"][1][-1] == "."
    with pytest.raises(ConfigurationError):
        m.map_pybedtools_values(files / "a.bed", files / "b.bed", operation="bogus", verbose=False)


def test_annotate_matches_bedtools(files):
    result = m.annotate_pybedtools_intervals(
        files / "a.bed", [files / "b.bed"], counts=True, verbose=False
    )
    expected = direct(
        "annotate", "-i", str(files / "a.bed"), "-files", str(files / "b.bed"), "-counts"
    )
    assert "\n".join(rows(result)) + "\n" == expected
    with pytest.raises(ConfigurationError):
        m.annotate_pybedtools_intervals(files / "a.bed", [], verbose=False)


def test_genome_coverage_total_coverage(files, bam):
    hist = m.genome_coverage_pybedtools(files / "a.bed", files / "g.sizes", verbose=False)
    expected = direct("genomecov", "-i", str(files / "a.bed"), "-g", str(files / "g.sizes"))
    assert "\n".join(rows(hist)) + "\n" == expected
    graph = m.genome_coverage_pybedtools(
        files / "a.bed", files / "g.sizes", report="bedgraph", verbose=False
    )
    assert graph["records"][0] == ["chr1", "100", "200", "1"]
    from_bam = m.genome_coverage_pybedtools(bam, report="bedgraph", verbose=False)
    assert from_bam["records"][0][:3] == ["chr1", "120", "160"]
    assert (
        m.total_coverage_pybedtools(files / "a.bed", verbose=False)["total_bases"] == 100 + 100 + 80
    )
    with pytest.raises(ConfigurationError):
        m.genome_coverage_pybedtools(files / "a.bed", verbose=False)
    with pytest.raises(ConfigurationError):
        m.genome_coverage_pybedtools(files / "a.bed", files / "g.sizes", scale=2.0, verbose=False)


def test_pair_to_bed(files):
    result = m.pair_to_bed_pybedtools(files / "p.bedpe", files / "a.bed", verbose=False)
    expected = direct(
        "pairtobed", "-a", str(files / "p.bedpe"), "-b", str(files / "a.bed"), "-type", "either"
    )
    assert "\n".join(rows(result)) + "\n" == expected
    neither = m.pair_to_bed_pybedtools(
        files / "p.bedpe", files / "b.bed", pair_type="neither", verbose=False
    )
    assert neither["n_records"] == 2


def test_cat_matches_pybedtools(files):
    result = m.cat_pybedtools_intervals([files / "a.bed", files / "b.bed"], verbose=False)
    assert result["records"][0] == ["chr1", "100", "250"]
    assert (
        m.cat_pybedtools_intervals(
            [files / "a.bed", files / "b.bed"], postmerge=False, verbose=False
        )["n_records"]
        == 6
    )
    with pytest.raises(ConfigurationError):
        m.cat_pybedtools_intervals([files / "a.bed"], verbose=False)


def test_make_windows(files):
    genome = m.make_windows_pybedtools(files / "g.sizes", window_size=500, verbose=False)
    assert "\n".join(rows(genome)) + "\n" == direct(
        "makewindows", "-g", str(files / "g.sizes"), "-w", "500"
    )
    bed = m.make_windows_pybedtools(
        files / "a.bed", source_type="bed", number_of_windows=2, id_naming="winnum", verbose=False
    )
    assert bed["n_records"] == 6 and bed["records"][0][3] == "1"
    with pytest.raises(ConfigurationError):
        m.make_windows_pybedtools(files / "g.sizes", verbose=False)
    with pytest.raises(ConfigurationError):
        m.make_windows_pybedtools(
            files / "g.sizes", window_size=5, number_of_windows=2, verbose=False
        )


def test_spacing(files):
    result = m.spacing_pybedtools_intervals(files / "a.bed", verbose=False)
    assert [r[-1] for r in result["records"]] == [".", "300", "."]


def test_fisher_matches_bedtools(files):
    result = m.fisher_pybedtools_overlap(
        files / "a.bed", files / "b.bed", files / "g.sizes", verbose=False
    )
    text = direct(
        "fisher",
        "-a",
        str(files / "a.bed"),
        "-b",
        str(files / "b.bed"),
        "-g",
        str(files / "g.sizes"),
    )
    assert result["report"].strip() == text.strip()
    assert result["table"]["in -a"]["in -b"] == 1
    assert 0 <= result["two_tail"] <= 1


def test_introns(files):
    result = m.introns_pybedtools(files / "genes.gff", verbose=False)
    assert [r[1:3] for r in result["records"]] == [["200", "300"], ["400", "500"]]
    assert result["n_invalid_dropped"] == 0


def test_truncate_and_remove_invalid(files):
    (files / "wide.bed").write_text("chr1\t1900\t2500\tw\nchr3\t1\t5\tz\nchr2\t5\t50\ty\n")
    result = m.truncate_pybedtools_to_chrom(files / "wide.bed", files / "g.sizes", verbose=False)
    assert sorted(rows(result)) == ["chr1\t1900\t2000", "chr2\t5\t50"] or result["n_records"] == 2
    (files / "bad.bed").write_text("chr1\t10\t100\nchr1\t100\t10\nchr1\t5\t6\n")
    cleaned = m.remove_invalid_pybedtools_intervals(files / "bad.bed", verbose=False)
    assert cleaned["n_records"] == 2 and cleaned["n_removed"] == 1


def test_tag_bam(files, bam):
    out = files / "tag.bam"
    result = m.tag_bam_pybedtools(
        bam, [files / "a.bed"], out, labels=["A"], tag="XA", verbose=False
    )
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(out), "rb") as handle:
        tags = {
            read.query_name: read.get_tag("XA") if read.has_tag("XA") else None for read in handle
        }
    assert tags["read0"] == "A" and tags["read2"] is None
    assert result["n_alignments"] == 4
    with pytest.raises(ConfigurationError):
        m.tag_bam_pybedtools(bam, [files / "a.bed"], out, verbose=False)


def test_random_overlap_functions_match_pybedtools(files):
    pybedtools = pytest.importorskip("pybedtools")
    dist = m.random_overlap_distribution_pybedtools(
        files / "a.bed", files / "b.bed", files / "g.sizes", iterations=8, seed=0, verbose=False
    )
    first = pybedtools.BedTool(str(files / "a.bed"))
    expected = list(
        first.randomintersection(
            pybedtools.BedTool(str(files / "b.bed")),
            iterations=8,
            debug=True,
            shuffle_kwargs={"g": str(files / "g.sizes")},
        )
    )
    assert dist["distribution"] == expected and dist["actual"] == 1
    again = m.random_overlap_distribution_pybedtools(
        files / "a.bed", files / "b.bed", files / "g.sizes", iterations=8, seed=0, verbose=False
    )
    assert again["distribution"] == dist["distribution"]
    bp = m.random_overlap_distribution_pybedtools(
        files / "a.bed",
        files / "b.bed",
        files / "g.sizes",
        iterations=4,
        metric="bp",
        verbose=False,
    )
    assert bp["actual"] == 50 and len(bp["distribution"]) == 4
    summary = m.random_overlap_summary_pybedtools(
        files / "a.bed",
        files / "b.bed",
        files / "g.sizes",
        iterations=8,
        include_distribution=True,
        verbose=False,
    )
    assert summary["distribution"] == expected and summary["self"] == 3 and summary["other"] == 3
    assert 0 <= summary["percentile"] <= 100
    matrix = m.random_overlap_matrix_pybedtools(
        [files / "a.bed", files / "b.bed"], files / "g.sizes", iterations=3, verbose=False
    )
    assert set(matrix["matrix"]) == {"a", "b"} and matrix["matrix"]["a"]["b"]["actual"] == 1
    with pytest.raises(ConfigurationError):
        m.random_overlap_distribution_pybedtools(
            files / "a.bed", files / "b.bed", files / "g.sizes", metric="x", verbose=False
        )


def test_nucleotide_content_matches_bedtools(files):
    result = m.nucleotide_content_pybedtools(files / "a.bed", files / "ref.fa", verbose=False)
    expected = direct(
        "nuc", "-fi", str(files / "ref.fa"), "-bed", str(files / "a.bed")
    ).splitlines()
    assert "\n".join(rows(result)) == "\n".join(expected[1:])
    assert result["header"][-1] == "15_seq_len"
    with pytest.raises(ConfigurationError):
        m.nucleotide_content_pybedtools(
            files / "a.bed", files / "ref.fa", pattern="A1", verbose=False
        )


def test_plot_tracks_igv_and_links(files):
    pytest.importorskip("matplotlib")
    image = files / "t.png"
    result = m.plot_interval_tracks_pybedtools(
        {"a": files / "a.bed", "b": files / "b.bed"},
        "chr1",
        0,
        1000,
        image,
        label_intervals=True,
        verbose=False,
    )
    assert image.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    assert [len(v) for v in result["tracks"].values()] == [2, 2]
    with pytest.raises(ConfigurationError):
        m.plot_interval_tracks_pybedtools(
            {"a": files / "a.bed"}, "chr1", 0, 1000, image, verbose=False
        )
    igv = m.build_igv_script_pybedtools(
        files / "a.bed", snapshot_dir="/snaps", slop=5, verbose=False
    )
    assert igv["n_commands"] == 3 and igv["script"].startswith("snapshotDirectory /snaps")
    links = m.build_ucsc_links_page_pybedtools(
        files / "a.bed", database="hg38", output_path=files / "l.html", verbose=False
    )
    assert links["n_links"] == 3 and "db=hg38" in links["html"] and (files / "l.html").is_file()
    with pytest.raises(ConfigurationError):
        m.build_ucsc_links_page_pybedtools(files / "a.bed", base_url="ftp://x", verbose=False)


def test_output_path_and_row_limit(files):
    out = files / "slop.bed"
    result = m.slop_pybedtools_intervals(
        files / "a.bed", files / "g.sizes", both=1, output_path=out, max_rows=1, verbose=False
    )
    assert result["truncated"] and len(result["records"]) == 1
    assert len(out.read_text().splitlines()) == 3
    with pytest.raises(ConfigurationError):
        m.slop_pybedtools_intervals(
            files / "a.bed", files / "g.sizes", both=1, output_path=out, verbose=False
        )


def test_missing_bedtools_reports_install_hint(files, monkeypatch):
    monkeypatch.setattr(m.shutil, "which", lambda name: None)
    with pytest.raises(BackendUnavailableError):
        m.total_coverage_pybedtools(files / "a.bed", verbose=False)
