import shutil
import subprocess

import pytest

from dnakit.annotation import bedtools_p13_cli_functions as b
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs = pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools not installed")

A = "chr1\t100\t200\ta1\t5\t+\nchr1\t500\t600\ta2\t6\t-\nchr2\t10\t90\ta3\t7\t+\n"
B = "chr1\t150\t250\tb1\t1\t+\nchr1\t900\t950\tb2\t2\t-\nchr2\t100\t150\tb3\t3\t+\n"
GENOME = "chr1\t2000\nchr2\t1000\n"


@pytest.fixture()
def files(tmp_path):
    for name, text in (("a.bed", A), ("b.bed", B), ("g.sizes", GENOME)):
        (tmp_path / name).write_text(text)
    return tmp_path


def direct(*args):
    return subprocess.run(["bedtools", *args], capture_output=True, text=True, check=True).stdout


def test_missing_program_and_validation(files):
    with pytest.raises(BackendUnavailableError):
        b.bedtools_window(
            str(files / "a.bed"),
            str(files / "b.bed"),
            str(files / "o"),
            bedtools_path=str(files / "nope"),
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        b.bedtools_window(
            str(files / "a.bed"),
            str(files / "b.bed"),
            str(files / "o"),
            window=5,
            left=1,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        b.bedtools_bam_to_bedpe(str(files / "a.bed"), str(files / "a.bed"), verbose=False)


@needs
def test_window_count_matches_bedtools(files):
    out = files / "w.bed"
    result = b.bedtools_window(
        str(files / "a.bed"),
        str(files / "b.bed"),
        str(out),
        window=100,
        report="count",
        verbose=False,
    )
    assert result["status"] == "completed" and result["version"]
    assert out.read_text() == direct(
        "window", "-a", str(files / "a.bed"), "-b", str(files / "b.bed"), "-w", "100", "-c"
    )
    with pytest.raises(ConfigurationError):  # never overwrites
        b.bedtools_window(str(files / "a.bed"), str(files / "b.bed"), str(out), verbose=False)


@needs
def test_random_shuffle_sample_are_seeded(files):
    r = b.bedtools_random(
        str(files / "g.sizes"), str(files / "r.bed"), count=7, length=50, seed=3, verbose=False
    )
    assert r["record_count"] == 7
    assert (files / "r.bed").read_text() == direct(
        "random", "-g", str(files / "g.sizes"), "-l", "50", "-n", "7", "-seed", "3"
    )
    s = b.bedtools_shuffle(
        str(files / "a.bed"), str(files / "g.sizes"), str(files / "s.bed"), seed=2, verbose=False
    )
    assert s["record_count"] == 3
    m = b.bedtools_sample(
        str(files / "a.bed"), str(files / "m.bed"), count=2, seed=2, verbose=False
    )
    assert m["record_count"] == 2


@needs
def test_multiinter_unionbedg_pairtopair_reldist_coverage_split(files):
    res = b.bedtools_multiinter(
        [str(files / "a.bed"), str(files / "b.bed")],
        str(files / "mi.txt"),
        names=["a", "b"],
        header=True,
        verbose=False,
    )
    assert res["status"] == "completed" and (files / "mi.txt").read_text().startswith("chrom")
    (files / "x.bg").write_text("chr1\t0\t100\t1\nchr1\t200\t300\t2\n")
    (files / "y.bg").write_text("chr1\t50\t250\t4\n")
    u = b.bedtools_unionbedg(
        [str(files / "x.bg"), str(files / "y.bg")], str(files / "u.txt"), verbose=False
    )
    assert u["record_count"] >= 3
    (files / "a.bedpe").write_text("chr1\t10\t20\tchr1\t100\t110\tp1\t1\t+\t-\n")
    pp = b.bedtools_pairtopair(
        str(files / "a.bedpe"), str(files / "a.bedpe"), str(files / "pp.txt"), verbose=False
    )
    assert pp["record_count"] == 1
    rd = b.bedtools_reldist(
        str(files / "a.bed"), str(files / "b.bed"), str(files / "rd.txt"), verbose=False
    )
    assert rd["summary"] and sum(row["fraction"] for row in rd["summary"]) == pytest.approx(
        1.0, abs=0.05
    )
    det = b.bedtools_reldist(
        str(files / "a.bed"),
        str(files / "b.bed"),
        str(files / "rdd.txt"),
        detail=True,
        verbose=False,
    )
    assert det["summary"] == [] and det["record_count"] >= 1
    h = b.bedtools_coverage_profile(
        str(files / "a.bed"), str(files / "b.bed"), str(files / "h.txt"), verbose=False
    )
    d = b.bedtools_coverage_profile(
        str(files / "a.bed"),
        str(files / "b.bed"),
        str(files / "d.txt"),
        report="depth",
        verbose=False,
    )
    assert d["record_count"] == 100 + 100 + 80 and h["report"] == "histogram"
    sp = b.bedtools_split(
        str(files / "a.bed"), str(files / "split"), 2, algorithm="simple", verbose=False
    )
    assert sum(f["record_count"] for f in sp["files"]) == 3 and len(sp["outputs"]) == 2


@needs
@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
def test_multicov_and_bamtobed_bedpe(files):
    sam = files / "p.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:2000\n"
        "p1\t99\tchr1\t101\t60\t50M\t=\t301\t250\t" + "A" * 50 + "\t" + "I" * 50 + "\n"
        "p1\t147\tchr1\t301\t60\t50M\t=\t101\t-250\t" + "A" * 50 + "\t" + "I" * 50 + "\n"
    )
    bam = files / "p.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    mc = b.bedtools_multicov([str(bam)], str(files / "a.bed"), str(files / "mc.txt"), verbose=False)
    assert mc["record_count"] == 3
    q = files / "q.bam"
    subprocess.run(["samtools", "sort", "-n", "-o", str(q), str(sam)], check=True)
    bp = b.bedtools_bam_to_bedpe(str(q), str(files / "bp.txt"), verbose=False)
    assert bp["record_count"] == 1
