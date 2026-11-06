"""Tests for dnakit.sequence.gimmemotifs_gap_functions (GimmeMotifs wrappers, Fig. 2 coverage gap).

Offline by design: inputs are the shared consistency data and hand-written PFM files, the gimme
command-line program is used only when installed, and the JASPAR request is mocked. Set
DNAKIT_LIVE_TESTS=1 to run the JASPAR lookup against the real service.
"""

from __future__ import annotations

import inspect
import os
from pathlib import Path

import pytest

pytest.importorskip("gimmemotifs")

from dnakit.exceptions import (
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.sequence import gimmemotifs_gap_functions as gap

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)

PFM_TEXT = """>MA0001.1 AR1
0.2 0.2 0.2 0.4
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
>MA0002.1 GATA
0.1 0.1 0.7 0.1
0.1 0.7 0.1 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
0.4 0.2 0.2 0.2
0.2 0.4 0.2 0.2
"""

DB_PFM_TEXT = """>DB1 long
0.25 0.25 0.25 0.25
0.1 0.1 0.7 0.1
0.8 0.1 0.05 0.05
0.05 0.8 0.1 0.05
0.6 0.2 0.1 0.1
0.1 0.7 0.1 0.1
0.25 0.25 0.25 0.25
0.25 0.25 0.25 0.25
"""


def _need(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"shared test data not present: {path}")
    return path


@pytest.fixture()
def pfm(tmp_path: Path) -> Path:
    path = tmp_path / "motifs.pfm"
    path.write_text(PFM_TEXT, encoding="utf-8")
    return path


@pytest.fixture()
def db_pfm(tmp_path: Path) -> Path:
    path = tmp_path / "database.pfm"
    path.write_text(DB_PFM_TEXT, encoding="utf-8")
    return path


@pytest.fixture()
def gimme_available() -> None:
    try:
        gap._gimme_program(None)
    except BackendUnavailableError:
        pytest.skip("gimme command-line program is not installed")


def _first_reference(bam: Path) -> str:
    pysam = pytest.importorskip("pysam")
    with pysam.AlignmentFile(str(bam)) as handle:
        return handle.references[0]


def test_public_functions_are_verbose_keyword_only() -> None:
    assert len(gap.__all__) == 21
    for name in gap.__all__:
        signature = inspect.signature(getattr(gap, name))
        assert "verbose" in signature.parameters, name
        assert signature.parameters["verbose"].kind is inspect.Parameter.KEYWORD_ONLY, name
        assert signature.parameters["verbose"].default is True, name


def test_missing_gimme_program_gives_install_hint(tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        gap.gimmemotifs_pfm_scan(
            str(_need(DATA / "genes.fa")),
            pfm_file=str(_write(tmp_path, "m.pfm")),
            output_dir=str(tmp_path / "scan"),
            cutoff=0.8,
            gimme_path=str(tmp_path / "no-such-gimme"),
        )
    assert "gimme" in str(info.value)


def _write(tmp_path: Path, name: str) -> Path:
    path = tmp_path / name
    path.write_text(PFM_TEXT, encoding="utf-8")
    return path


@pytest.mark.parametrize("bad_cutoff", [0, -0.1, 1.5, "0.8"])
def test_pfm_scan_validates_cutoff_before_running(tmp_path: Path, pfm: Path, bad_cutoff) -> None:
    with pytest.raises(ConfigurationError):
        gap.gimmemotifs_pfm_scan(
            str(_need(DATA / "genes.fa")),
            pfm_file=str(pfm),
            output_dir=str(tmp_path / "scan"),
            cutoff=bad_cutoff,
        )


def test_pfm_scan_rejects_unknown_format(tmp_path: Path, pfm: Path) -> None:
    with pytest.raises(ConfigurationError):
        gap.gimmemotifs_pfm_scan(
            str(_need(DATA / "genes.fa")),
            pfm_file=str(pfm),
            output_dir=str(tmp_path / "scan"),
            cutoff=0.8,
            output_format="fasta",
        )


def test_pfm_scan_table_with_gimme(tmp_path: Path, pfm: Path, gimme_available) -> None:
    result = gap.gimmemotifs_pfm_scan(
        str(_need(DATA / "genes.fa")),
        pfm_file=str(pfm),
        output_dir=str(tmp_path / "scan"),
        cutoff=0.8,
        threads=1,
    )
    assert result["status"] == "complete"
    scan = Path(result["outputs"]["scan"])
    assert scan.is_file()
    assert result["row_count"] == len(result["rows"]) > 0
    assert "sequence" in result["rows"][0]
    assert "MA0001.1 AR1" in result["rows"][0]


def test_match_motifs_table_with_gimme(
    tmp_path: Path, pfm: Path, db_pfm: Path, gimme_available
) -> None:
    result = gap.gimmemotifs_match_motifs(
        str(pfm), database_pfm=str(db_pfm), output_path=str(tmp_path / "match.tsv")
    )
    assert Path(result["outputs"]["table"]).is_file()
    assert result["row_count"] == 2
    assert {row["Motif"] for row in result["rows"]} == {"MA0001.1 AR1", "MA0002.1 GATA"}


def test_match_motifs_requires_database(tmp_path: Path, pfm: Path) -> None:
    with pytest.raises(TypeError):
        gap.gimmemotifs_match_motifs(str(pfm), output_path=str(tmp_path / "m.tsv"))  # type: ignore[call-arg]


def test_logo_images_with_gimme(tmp_path: Path, pfm: Path, gimme_available) -> None:
    result = gap.gimmemotifs_logo_images(
        str(pfm), output_dir=str(tmp_path / "logos"), motif_ids=["MA0001.1 AR1"], kind="information"
    )
    images = [Path(value) for value in result["outputs"].values()]
    assert images and all(path.suffix == ".png" and path.is_file() for path in images)
    assert all(str(path).startswith(str(tmp_path / "logos")) for path in images)


def test_motif_location_plots_with_gimme(tmp_path: Path, pfm: Path, gimme_available) -> None:
    result = gap.gimmemotifs_motif_location_plots(
        str(_need(DATA / "genes.fa")),
        pfm_file=str(pfm),
        output_dir=str(tmp_path / "loc"),
        size=200,
        motif_ids=["MA0002.1 GATA"],
        cutoff=0.7,
    )
    assert any(path.endswith("_histogram.svg") for path in result["outputs"].values())


def test_cluster_pfm_with_gimme(tmp_path: Path, pfm: Path, gimme_available) -> None:
    result = gap.gimmemotifs_cluster_pfm(
        str(pfm), output_dir=str(tmp_path / "clusters"), threshold=0.7, threads=1
    )
    assert Path(result["outputs"]["clustered_pfm"]).is_file()
    assert Path(result["outputs"]["report"]).is_file()


def test_combine_peak_summits_merges_nearby_summits(tmp_path: Path) -> None:
    summits = tmp_path / "summits.narrowPeak"
    summits.write_text(
        "chr1\t1000\t1100\tp1\t50\t.\t5.0\t3.0\t2.0\t50\n"
        "chr1\t1050\t1150\tp2\t80\t.\t6.0\t4.0\t2.0\t60\n"
        "chr1\t5000\t5100\tp3\t30\t.\t4.0\t2.0\t1.0\t20\n",
        encoding="utf-8",
    )
    output = tmp_path / "combined.tsv"
    result = gap.gimmemotifs_combine_peak_summits(
        [str(summits)],
        genome_fasta=str(_need(DATA / "genome.fa")),
        window=200,
        output_path=str(output),
    )
    assert result["row_count"] >= 1
    assert output.is_file()
    assert result["columns"]


def test_coverage_table_counts_reads_in_regions(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    bam = _need(DATA / "reads.bam")
    if not (DATA / "reads.bam.bai").exists():
        pytest.skip("BAM index not present")
    chrom = _first_reference(bam)
    summits = tmp_path / "summits.narrowPeak"
    summits.write_text(
        f"{chrom}\t1000\t1100\tp1\t50\t.\t5.0\t3.0\t2.0\t50\n"
        f"{chrom}\t3000\t3100\tp2\t60\t.\t5.0\t3.0\t2.0\t50\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"
    result = gap.gimmemotifs_coverage_table(
        str(summits), [str(bam)], output_path=str(output), log_transform=False, threads=1
    )
    assert result["row_count"] > 0
    assert output.is_file()


def test_coverage_table_requires_bam_index(tmp_path: Path) -> None:
    _need(DATA / "reads.bam")
    bam = tmp_path / "unindexed.bam"
    bam.write_bytes((DATA / "reads.bam").read_bytes())
    summits = tmp_path / "summits.narrowPeak"
    summits.write_text("chr1\t1000\t1100\tp1\t50\t.\t5.0\t3.0\t2.0\t50\n", encoding="utf-8")
    with pytest.raises(ConfigurationError, match="index"):
        gap.gimmemotifs_coverage_table(str(summits), [str(bam)], threads=1)


def test_histogram_plot_writes_svg(tmp_path: Path) -> None:
    result = gap.gimmemotifs_histogram_plot(
        [1.0, 2.0, 2.5, 3.0, 9.0],
        str(tmp_path / "hist"),
        x_range=(0.0, 10.0),
        bins=5,
        title="positions",
        x_label="Position",
    )
    image = Path(result["outputs"]["image"])
    assert image.name == "hist.svg" and image.is_file()
    assert result["count"] == 5 and result["maximum"] == 9.0


def test_histogram_plot_rejects_empty_values(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        gap.gimmemotifs_histogram_plot([], str(tmp_path / "empty"))


def test_roc_curve_plot_single_and_multiple(tmp_path: Path) -> None:
    single = gap.gimmemotifs_roc_curve_plot(
        str(tmp_path / "roc_single"), plot_x=[0.0, 0.2, 1.0], plot_y=[0.0, 0.6, 1.0]
    )
    assert Path(single["outputs"]["image"]).name == "roc_single.png"
    assert Path(single["outputs"]["image"]).is_file()
    multi = gap.gimmemotifs_roc_curve_plot(
        str(tmp_path / "roc_multi.png"),
        plot_x=[[0.0, 0.5, 1.0], [0.0, 0.1, 1.0]],
        plot_y=[[0.0, 0.8, 1.0], [0.0, 0.3, 1.0]],
        curve_ids=["a", "b"],
    )
    assert multi["curves"] == 2
    assert Path(multi["outputs"]["image"]).is_file()


def test_match_plot_writes_png(tmp_path: Path, pfm: Path, db_pfm: Path) -> None:
    result = gap.gimmemotifs_match_plot(
        [
            {
                "query_pfm": str(pfm),
                "query_id": "MA0001.1 AR1",
                "database_pfm": str(db_pfm),
                "database_id": "DB1 long",
                "p_value": 0.001,
            }
        ],
        str(tmp_path / "match.png"),
    )
    assert Path(result["outputs"]["image"]).is_file()
    assert result["pairs"] == 1


def test_motif_localization_reports_p_value(tmp_path: Path, pfm: Path) -> None:
    result = gap.gimmemotifs_motif_localization(
        str(_need(DATA / "genes.fa")),
        pfm_file=str(pfm),
        motif_id="MA0001.1 AR1",
        output_path=str(tmp_path / "loc"),
        size=200,
        cutoff=0.7,
    )
    assert result["motif_id"] == "MA0001.1 AR1"
    assert 0.0 <= result["p_value"] <= 1.0
    assert Path(result["outputs"]["image"]).suffix == ".svg"


def test_motif_localization_rejects_unknown_motif(tmp_path: Path, pfm: Path) -> None:
    with pytest.raises(ConfigurationError, match="not in"):
        gap.gimmemotifs_motif_localization(
            str(_need(DATA / "genes.fa")),
            pfm_file=str(pfm),
            motif_id="missing",
            output_path=str(tmp_path / "loc"),
            size=200,
        )


def test_cluster_motif_report_writes_html(tmp_path: Path, pfm: Path) -> None:
    result = gap.gimmemotifs_cluster_motif_report(
        str(pfm), output_dir=str(tmp_path / "report"), threshold=0.7, title="test"
    )
    assert Path(result["outputs"]["clustered_pfm"]).is_file()
    assert Path(result["outputs"]["report"]).is_file()
    assert result["cluster_count"] == len(result["clusters"]) == 2


def test_motif_image_series_writes_logos(tmp_path: Path, pfm: Path) -> None:
    result = gap.gimmemotifs_motif_image_series(
        ["MA0001.1 AR1", "MA0002.1 GATA"], pfm_file=str(pfm), output_dir=str(tmp_path / "series")
    )
    assert result["image_count"] == 2
    for image in result["images"].values():
        assert Path(image).is_file()


def test_create_background_random_from_input(tmp_path: Path) -> None:
    result = gap.gimmemotifs_create_background(
        str(tmp_path / "random.fa"),
        background_type="random",
        input_file=str(_need(DATA / "genes.fa")),
        number=50,
    )
    text = Path(result["outputs"]["background"]).read_text()
    assert text.count(">") == 50


def test_create_background_genomic_from_genome_fasta(tmp_path: Path) -> None:
    result = gap.gimmemotifs_create_background(
        str(tmp_path / "genomic.fa"),
        background_type="genomic",
        genome=str(_need(DATA / "genome.fa")),
        size=200,
        number=20,
    )
    assert Path(result["outputs"]["background"]).read_text().count(">") == 20


@pytest.mark.parametrize(
    "kwargs",
    [
        {"background_type": "promoter", "genome": "hg38"},
        {"background_type": "markov"},
        {"background_type": "random"},
        {"background_type": "random", "input_file": "x", "fmt": "bed"},
        {"background_type": "genomic", "genome": None},
        {"background_type": "genomic", "genome": "hg38", "size": 200},
        {"background_type": "gc", "genome": "/nonexistent/genome.fa"},
    ],
)
def test_create_background_rejects_invalid_inputs(tmp_path: Path, kwargs) -> None:
    with pytest.raises(ConfigurationError):
        gap.gimmemotifs_create_background(str(tmp_path / "bg.fa"), **kwargs)


def test_gc_bin_index_and_bed_from_genome_fasta(tmp_path: Path) -> None:
    genome = str(_need(DATA / "genome.fa"))
    index = gap.gimmemotifs_gc_bin_index(genome, str(tmp_path / "gc.index"), min_bin_size=10)
    assert Path(index["outputs"]["index"]).exists()
    bed = gap.gimmemotifs_gc_bin_bed(
        str(tmp_path / "gcbin.bed"), genome=genome, number=20, random_state=1, min_bin_size=10
    )
    lines = Path(bed["outputs"]["bed"]).read_text().splitlines()
    assert len(lines) >= 1
    assert all(len(line.split("\t")) >= 3 for line in lines)


def test_matched_gc_bed_follows_input_regions(tmp_path: Path) -> None:
    result = gap.gimmemotifs_matched_gc_bed(
        str(tmp_path / "matched.bed"),
        match_bed=str(_need(DATA / "regions_a.bed")),
        genome=str(_need(DATA / "genome.fa")),
        number=10,
        min_bin_size=10,
    )
    assert Path(result["outputs"]["bed"]).is_file()


def test_jaspar_lookup_uses_the_gimmemotifs_request(monkeypatch: pytest.MonkeyPatch) -> None:
    import gimmemotifs.utils as utils

    calls: list[str] = []

    class _Response:
        ok = True

        def json(self) -> dict:
            return {"matrix_id": "MA0001", "name": "AR1"}

    def fake_get(url: str, *args, **kwargs):
        calls.append(url)
        return _Response()

    monkeypatch.setattr(utils.requests, "get", fake_get)
    result = gap.gimmemotifs_jaspar_motif_info("MA0001.1")
    assert calls and calls[0].endswith("/MA0001.1?format=json")
    assert result["info"]["name"] == "AR1"


@pytest.mark.skipif(os.environ.get("DNAKIT_LIVE_TESTS") != "1", reason="live JASPAR request")
def test_live_jaspar_lookup() -> None:
    result = gap.gimmemotifs_jaspar_motif_info("MA0001.1")
    assert result["info"]["matrix_id"] == "MA0001"


def test_denovo_functions_need_an_installed_tool(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    import gimmemotifs.tools as tools

    class _NotInstalled:
        def is_installed(self) -> bool:
            return False

    monkeypatch.setattr(tools, "get_tool", lambda name: _NotInstalled())
    fasta = str(_need(DATA / "genes.fa"))
    with pytest.raises(BackendUnavailableError, match="de novo"):
        gap.gimmemotifs_denovo_motifs(fasta, output_dir=str(tmp_path / "denovo"))
    with pytest.raises(BackendUnavailableError):
        gap.gimmemotifs_denovo_predict(fasta, fasta, str(tmp_path / "motifs.pfm"))
    with pytest.raises(BackendUnavailableError):
        gap.gimmemotifs_parallel_denovo_predict(fasta, str(tmp_path / "pp.pfm"))
    assert not (tmp_path / "denovo").exists()
