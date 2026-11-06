import types

import pytest

from dnakit.annotation import pyranges1_p33_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

try:
    import pyranges1 as pr

    NEW = hasattr(pr, "read_parquet") and hasattr(pr.PyRanges, "liftover_ranges")
except ImportError:  # pragma: no cover - optional dependency
    pr = None
    NEW = False

needs = pytest.mark.skipif(not NEW, reason="pyranges1 >= 1.5.0 not installed")


def _attrs(gene, transcript, extra=""):
    return f'gene_id "{gene}"; transcript_id "{transcript}";{extra}'


GTF_ROWS = [
    ("gene", 101, 1000, ".", 'gene_id "g1";'),
    ("transcript", 101, 1000, ".", _attrs("g1", "t1", ' tag "Ensembl_canonical"; tag "basic";')),
    ("exon", 101, 300, ".", _attrs("g1", "t1")),
    ("exon", 501, 1000, ".", _attrs("g1", "t1")),
    ("CDS", 201, 300, "0", _attrs("g1", "t1")),
    ("CDS", 501, 800, "0", _attrs("g1", "t1")),
    ("UTR", 101, 200, ".", _attrs("g1", "t1")),
    ("UTR", 801, 1000, ".", _attrs("g1", "t1")),
    ("transcript", 101, 600, ".", _attrs("g1", "t2", ' tag "basic";')),
    ("exon", 101, 600, ".", _attrs("g1", "t2")),
]


@pytest.fixture()
def files(tmp_path):
    def write(name, text):
        (tmp_path / name).write_text(text)

    write(
        "iv.bed",
        "chr1\t100\t200\ta\t1\t+\nchr1\t150\t300\tb\t2\t+\nchr1\t500\t600\tc\t3\t-\nchr2\t10\t50\td\t4\t+\n",
    )
    write("other.bed", "chr1\t180\t250\tx\t1\t+\nchr1\t550\t700\ty\t2\t-\n")
    lines = [
        "\t".join(["chr1", "src", f, str(s), str(e), ".", "+", fr, at])
        for f, s, e, fr, at in GTF_ROWS
    ]
    write("ann.gtf", "\n".join(lines) + "\n")
    exons = [line for line in lines if "\texon\t" in line and '"t1"' in line]
    write("exons.gtf", "\n".join(exons) + "\n")
    write("q.bed", "t1\t0\t50\tq1\t0\t+\nt1\t190\t210\tq2\t0\t+\n")
    write("g.bed", "chr1\t110\t130\tq3\t0\t+\n")
    write(
        "np.narrowPeak",
        "chr1\t100\t200\tp1\t500\t.\t5.0\t4.0\t3.0\t50\nchr1\t300\t400\tp2\t600\t.\t6.0\t5.0\t4.0\t40\n",
    )
    write("bp.broadPeak", "chr1\t100\t200\tp1\t500\t.\t5.0\t4.0\t3.0\n")
    write(
        "gp.gappedPeak", "chr1\t100\t200\tp1\t500\t.\t100\t200\t0\t2\t20,20\t0,80\t5.0\t4.0\t3.0\n"
    )
    write("t.wig", "fixedStep chrom=chr1 start=101 step=10 span=10\n1.5\n2.5\n")
    write("c.chain", "chain 1000 chr1 1000 + 100 900 chrA 1500 + 200 1000 1\n400 100 100\n300\n\n")
    write("l.bed", "chr1\t150\t250\tm\t0\t+\nchr1\t450\t480\tn\t0\t+\nchr1\t5\t20\tz\t0\t+\n")
    write("sizes.txt", "chr1\t1000\nchr2\t500\n")
    write(
        "dups.bed",
        "chr1\t10\t20\ta\t1\t+\nchr1\t10\t20\tb\t2\t+\nchr1\t10\t20\tc\t3\t-\nchr1\t30\t40\td\t4\t+\n",
    )
    write("dot.bed", "chr1\t10\t20\ta\t1\t.\nchr1\t30\t40\tb\t2\t+\n")
    write("sig.bedgraph", "chr1\t100\t200\t2.5\nchr1\t300\t400\t4.0\n")
    return tmp_path


def test_old_release_reports_upgrade_hint(monkeypatch):
    fake = types.SimpleNamespace(PyRanges=types.SimpleNamespace())
    monkeypatch.setattr(m, "require", lambda module, extra: fake)
    with pytest.raises(BackendUnavailableError):
        m._pyranges()


@needs
def test_readers(files):
    wig = m.read_pyranges1_intervals(files / "t.wig", "wig", verbose=False)
    assert wig["rows"] == [
        {"Chromosome": "chr1", "Start": 100, "End": 110, "Value": 1.5},
        {"Chromosome": "chr1", "Start": 110, "End": 120, "Value": 2.5},
    ]
    broad = m.read_pyranges1_intervals(files / "bp.broadPeak", "broadPeak", verbose=False)
    assert broad["rows"][0]["SignalValue"] == 5.0 and broad["rows"][0]["QValue"] == 3.0
    gapped = m.read_pyranges1_intervals(files / "gp.gappedPeak", "gappedPeak", verbose=False)
    assert gapped["rows"][0]["BlockSizes"] == "20,20"
    chain = m.read_pyranges1_intervals(files / "c.chain", "chain", verbose=False)
    assert [(r["Start"], r["End"], r["QueryStart"]) for r in chain["rows"]] == [
        (100, 500, 200),
        (600, 900, 700),
    ]
    import pandas as pd

    pd.DataFrame({"Chromosome": ["chr1"], "Start": [1], "End": [5]}).to_parquet(files / "x.parquet")
    table = m.read_pyranges1_intervals(files / "x.parquet", "parquet", verbose=False)
    assert table["rows"] == [{"Chromosome": "chr1", "Start": 1, "End": 5}]
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_intervals(files / "t.wig", "bogus", verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_intervals(files / "iv.bed", "parquet", verbose=False)


@needs
def test_gene_annotation_helpers(files):
    tags = m.has_pyranges1_tag(files / "ann.gtf", "Ensembl_canonical", verbose=False)
    assert tags["n_tagged"] == 1
    flagged = [r["transcript_id"] for r in tags["rows"] if r["has_tag"]]
    assert flagged == ["t1"]
    selected = m.select_pyranges1_transcripts(files / "ann.gtf", verbose=False)
    assert {r["transcript_id"] for r in selected["rows"]} == {None, "t1"}
    longest = m.select_pyranges1_transcripts(files / "ann.gtf", how="longest", verbose=False)
    assert longest["n_rows"] == selected["n_rows"] == 8
    out = files / "sel.tsv"
    m.select_pyranges1_transcripts(files / "ann.gtf", output_path=out, verbose=False)
    assert len(out.read_text().splitlines()) == 9
    utr = m.label_pyranges1_utrs(files / "ann.gtf", verbose=False)
    assert (
        utr["feature_counts"]["five_prime_utr"] == 1
        and utr["feature_counts"]["three_prime_utr"] == 1
    )
    assert "UTR" not in utr["feature_counts"]


@needs
def test_interval_arithmetic(files):
    outer = m.outer_pyranges1_ranges(files / "iv.bed", verbose=False)
    assert [(r["Chromosome"], r["Start"], r["End"]) for r in outer["rows"]] == [
        ("chr1", 100, 300),
        ("chr1", 500, 600),
        ("chr2", 10, 50),
    ]
    resized = m.resize_pyranges1_ranges(files / "iv.bed", 10, fix="start", verbose=False)
    assert [(r["Start"], r["End"]) for r in resized["rows"]] == [
        (100, 110),
        (150, 160),
        (590, 600),
        (10, 20),
    ]
    five = m.end_positions_pyranges1(files / "iv.bed", verbose=False)
    three = m.end_positions_pyranges1(files / "iv.bed", end="three", verbose=False)
    assert [(r["Start"], r["End"]) for r in five["rows"]] == [
        (100, 101),
        (150, 151),
        (599, 600),
        (10, 11),
    ]
    assert [(r["Start"], r["End"]) for r in three["rows"]] == [
        (199, 200),
        (299, 300),
        (500, 501),
        (49, 50),
    ]
    lengths = m.interval_lengths_pyranges1(files / "iv.bed", verbose=False)
    assert lengths["lengths"] == [100, 150, 100, 40] and lengths["sum"] == 390
    up = m.flank_pyranges1_ranges(files / "iv.bed", 20, verbose=False)
    down = m.flank_pyranges1_ranges(files / "iv.bed", 20, side="downstream", gap=5, verbose=False)
    assert (up["rows"][0]["Start"], up["rows"][0]["End"]) == (80, 100)
    assert (up["rows"][2]["Start"], up["rows"][2]["End"]) == (
        600,
        620,
    )  # minus strand: 5' side is right
    assert (down["rows"][0]["Start"], down["rows"][0]["End"]) == (205, 225)
    windows = m.window_pyranges1_ranges(files / "iv.bed", 60, add_window_id=True, verbose=False)
    assert windows["n_rows"] == 2 + 3 + 2 + 1
    disjoint = m.max_disjoint_pyranges1_overlaps(files / "iv.bed", verbose=False)
    assert disjoint["n_rows"] == 3
    sliced = m.slice_pyranges1_ranges(files / "iv.bed", start=0, end=10, verbose=False)
    assert next((r["Start"], r["End"]) for r in sliced["rows"]) == (100, 110)
    summary = m.summarize_pyranges1_intervals(files / "iv.bed", verbose=False)
    assert summary["index"][0] == "count"
    with pytest.raises(ConfigurationError):
        m.resize_pyranges1_ranges(files / "iv.bed", 10, fix="middle", verbose=False)
    with pytest.raises(ConfigurationError):
        m.flank_pyranges1_ranges(files / "iv.bed", 0, verbose=False)


@needs
def test_spliced_slice_and_cumsum(files):
    sliced = m.slice_pyranges1_ranges(
        files / "exons.gtf",
        start=0,
        end=250,
        input_format="gtf",
        group_by="transcript_id",
        verbose=False,
    )
    assert [(r["Start"], r["End"]) for r in sliced["rows"]] == [(100, 300), (500, 550)]
    cumulative = m.group_cumsum_pyranges1(
        files / "exons.gtf", input_format="gtf", group_by="transcript_id", verbose=False
    )
    assert [(r["Start"], r["End"]) for r in cumulative["rows"]] == [(0, 200), (200, 700)]


@needs
def test_join_metrics_and_combination(files):
    joined = m.join_pyranges1_overlaps(files / "iv.bed", files / "other.bed", verbose=False)
    assert [(r["Name"], r["Name_b"]) for r in joined["rows"]] == [
        ("a", "x"),
        ("b", "x"),
        ("c", "y"),
    ]
    left = m.join_pyranges1_overlaps(
        files / "iv.bed", files / "other.bed", join_type="left", overlap_column="ov", verbose=False
    )
    assert left["n_rows"] == 4 and "ov" in left["columns"]
    metrics = m.interval_metrics_pyranges1(
        files / "iv.bed", files / "other.bed", metrics=["fraction", "jaccard"], verbose=False
    )
    assert metrics["n_rows"] == 3 and {"fraction", "jaccard"} <= set(metrics["columns"])
    first = metrics["rows"][0]
    assert first["fraction"] == pytest.approx(20 / 100)  # overlap 180-200 over the length of a
    combined = m.combine_pyranges1_interval_columns(
        files / "iv.bed", files / "other.bed", verbose=False
    )
    assert next((r["Start"], r["End"]) for r in combined["rows"]) == (180, 200)
    union = m.combine_pyranges1_interval_columns(
        files / "iv.bed", files / "other.bed", function="union", verbose=False
    )
    assert next((r["Start"], r["End"]) for r in union["rows"]) == (100, 250)


@needs
def test_coordinate_mapping_both_directions(files):
    to_global = m.map_pyranges1_coordinates(
        files / "q.bed", files / "exons.gtf", "transcript_id", verbose=False
    )
    spans = [(r["Chromosome"], r["Start"], r["End"], r["Name"]) for r in to_global["rows"]]
    assert spans == [("chr1", 100, 150, "q1"), ("chr1", 290, 300, "q2"), ("chr1", 500, 510, "q2")]
    to_local = m.map_pyranges1_coordinates(
        files / "g.bed", files / "exons.gtf", "transcript_id", direction="to_local", verbose=False
    )
    assert [(r["Chromosome"], r["Start"], r["End"]) for r in to_local["rows"]] == [("t1", 10, 30)]
    with pytest.raises(ConfigurationError):
        m.map_pyranges1_coordinates(
            files / "g.bed",
            files / "exons.gtf",
            "transcript_id",
            direction="sideways",
            verbose=False,
        )


@needs
def test_liftover_matches_chain_offsets(files):
    result = m.liftover_pyranges1_ranges(
        files / "l.bed", files / "c.chain", return_unmapped=True, verbose=False
    )
    lifted = {r["Name"]: (r["Chromosome"], r["Start"], r["End"]) for r in result["rows"]}
    assert lifted == {"m": ("chrA", 250, 350), "n": ("chrA", 550, 580)}
    assert result["n_unmapped"] == 1 and result["unmapped_rows"][0]["Name"] == "z"
    with pytest.raises(ConfigurationError):
        m.liftover_pyranges1_ranges(
            files / "l.bed", files / "c.chain", min_match=2.0, verbose=False
        )


@needs
def test_strand_and_duplicate_repair(files):
    fixed = m.make_pyranges1_strand_valid(files / "dot.bed", verbose=False)
    assert fixed["n_changed"] == 1 and {r["Strand"] for r in fixed["rows"]} == {"+"}
    stranded = m.drop_pyranges1_duplicate_positions(files / "dups.bed", verbose=False)
    assert stranded["n_dropped"] == 1 and [r["Name"] for r in stranded["rows"]] == ["a", "c", "d"]
    unstranded = m.drop_pyranges1_duplicate_positions(
        files / "dups.bed", use_strand="no", keep="last", verbose=False
    )
    assert [r["Name"] for r in unstranded["rows"]] == ["c", "d"]
    none = m.drop_pyranges1_duplicate_positions(
        files / "dups.bed", keep="none", use_strand="no", verbose=False
    )
    assert [r["Name"] for r in none["rows"]] == ["d"]


@needs
def test_writers_round_trip(files):
    for source, kind, name in (
        ("np.narrowPeak", "narrowPeak", "out.narrowPeak"),
        ("bp.broadPeak", "broadPeak", "out.broadPeak"),
        ("gp.gappedPeak", "gappedPeak", "out.gappedPeak"),
    ):
        result = m.write_pyranges1_intervals(
            files / source, files / name, kind, input_format=kind, verbose=False
        )
        assert result["output_format"] == kind
        assert (files / name).read_text() == (files / source).read_text()
    table = m.write_pyranges1_intervals(
        files / "iv.bed", files / "out.parquet", "parquet", verbose=False
    )
    back = m.read_pyranges1_intervals(files / "out.parquet", "parquet", verbose=False)
    assert table["n_rows"] == back["n_rows"] == 4
    graph = m.write_pyranges1_intervals(
        files / "sig.bedgraph",
        files / "out.bedgraph",
        "bedgraph",
        input_format="bedgraph",
        value_column="Value",
        verbose=False,
    )
    assert graph["n_rows"] == 2 and (files / "out.bedgraph").read_text().count("\n") >= 2
    with pytest.raises(ConfigurationError):
        m.write_pyranges1_intervals(
            files / "iv.bed", files / "out.parquet", "parquet", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.write_pyranges1_intervals(files / "iv.bed", files / "o.x", "wig", verbose=False)


@needs
def test_pairs_writer(files):
    (files / "in.pairs").write_text(
        "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n"
        "r1\tchr1\t100\tchr2\t200\t+\t-\n"
    )
    result = m.write_pyranges1_intervals(
        files / "in.pairs", files / "out.pairs", "pairs", input_format="pairs", verbose=False
    )
    assert result["n_rows"] == 1
    assert "r1\tchr1\t100\tchr2\t200\t+\t-" in (files / "out.pairs").read_text()


@needs
def test_bigwig_round_trip_and_stats(files):
    pytest.importorskip("pyBigWig")
    pytest.importorskip("runlength")
    pytest.importorskip("pyarrow")
    out = files / "cov.bw"
    result = m.coverage_to_bigwig_pyranges1(
        files / "iv.bed", out, chromsizes_path=files / "sizes.txt", rpm=False, verbose=False
    )
    assert result["n_rows"] == 4 and out.stat().st_size > 0
    stats = m.bigwig_stats_pyranges1(files / "iv.bed", out, stats=["mean", "max"], verbose=False)
    assert [r["mean"] for r in stats["rows"]][3] == pytest.approx(1.0)
    assert stats["rows"][0]["max"] == 2.0  # a and b overlap at 150-200
    binned = m.bigwig_stats_pyranges1(files / "iv.bed", out, bins=2, verbose=False)
    assert "mean_0" in binned["columns"] and "mean_1" in binned["columns"]
    with pytest.raises(ConfigurationError):
        m.bigwig_stats_pyranges1(files / "iv.bed", out, stats="median", verbose=False)


@needs
def test_coverage_runs(files):
    pytest.importorskip("runlength")
    runs = m.coverage_runs_pyranges1(files / "iv.bed", verbose=False)
    first = [r for r in runs["rows"] if r["Chromosome"] == "chr1" and r["Strand"] == "+"]
    assert [(r["start"], r["end"], r["value"]) for r in first] == [
        (0, 100, 0),
        (100, 150, 1),
        (150, 200, 2),
        (200, 300, 1),
    ]
    plain = m.coverage_runs_pyranges1(files / "iv.bed", use_strand="no", verbose=False)
    assert plain["columns"][0] == "Chromosome" and "Strand" not in plain["columns"]


@needs
def test_row_limit_and_output_file(files):
    out = files / "w.tsv"
    result = m.window_pyranges1_ranges(
        files / "iv.bed", 50, output_path=out, max_rows=2, verbose=False
    )
    assert result["truncated"] and len(result["rows"]) == 2
    assert len(out.read_text().splitlines()) == result["n_rows"] + 1
    with pytest.raises(ConfigurationError):
        m.window_pyranges1_ranges(files / "iv.bed", 50, output_path=out, verbose=False)
