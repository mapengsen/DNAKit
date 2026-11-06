"""Tests for dnakit.epigenome.macs3_gap_functions (MACS3 subcommands as typed facades).

The shared BAM has forward-strand reads only, so MACS3 cannot build a fragment model from it; the model tests use a
strand-symmetric BAM written into tmp_path with a known spacing between the plus and minus strand reads.
"""

from __future__ import annotations

import random
import shutil
import sys
from pathlib import Path

import pytest

from dnakit.epigenome import macs3_gap_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
READS = DATA / "reads.bam"
_READ_LENGTH = 50
_SPACING = 200


def _macs3_installed() -> bool:
    return shutil.which("macs3") is not None or (Path(sys.executable).parent / "macs3").is_file()


requires_macs3 = pytest.mark.skipif(
    not _macs3_installed(), reason="macs3 (python -m pip install macs3) is not installed"
)
requires_shared_bam = pytest.mark.skipif(
    not READS.exists(), reason="shared consistency/data/reads.bam is missing"
)


def _strand_pair_bam(path: Path) -> Path:
    """A sorted, indexed BAM with plus reads and minus reads 200 bp apart in 40 clusters per 1 kb, plus background."""
    pysam = pytest.importorskip("pysam")
    rng = random.Random(5)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 200000}]}
    records = []

    def add(position: int, reverse: bool) -> None:
        segment = pysam.AlignedSegment()
        segment.query_name = f"read{len(records)}"
        segment.query_sequence = "A" * _READ_LENGTH
        segment.query_qualities = pysam.qualitystring_to_array("I" * _READ_LENGTH)
        segment.reference_id = 0
        segment.reference_start = position
        segment.cigarstring = f"{_READ_LENGTH}M"
        segment.mapping_quality = 60
        segment.flag = 16 if reverse else 0
        records.append(segment)

    for cluster in range(0, 190_000, 1000):
        start = cluster + rng.randint(0, 300)
        for _ in range(40):
            add(start + rng.randint(0, 150), False)
            add(start + _SPACING - _READ_LENGTH + rng.randint(0, 150), True)
    for _ in range(300):
        add(rng.randint(0, 199_000), rng.random() < 0.5)
    records.sort(key=lambda record: record.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for record in records:
            handle.write(record)
    pysam.index(str(path))
    return path


@pytest.fixture
def pair_bam(tmp_path: Path) -> Path:
    return _strand_pair_bam(tmp_path / "pairs.bam")


# ---------------------------------------------------------------------------
# callpeak
# ---------------------------------------------------------------------------


@requires_macs3
@requires_shared_bam
def test_callpeak_on_shared_bam_writes_peaks_summits_and_tracks(tmp_path: Path) -> None:
    out = tmp_path / "callpeak"
    result = m.call_chip_peaks(
        [str(READS)],
        str(out),
        file_format="BAM",
        genome_size=200000,
        nomodel=True,
        extsize=200,
        name="shared",
        save_bedgraph=True,
        call_summits=True,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert "abs_summit" in (out / "shared_peaks.xls").read_text()
    assert (out / "shared_peaks.narrowPeak").is_file()
    assert (out / "shared_summits.bed").is_file()
    assert (out / "shared_treat_pileup.bdg").read_text().startswith("chr1\t")
    assert (out / "shared_control_lambda.bdg").is_file()
    assert result["model_script"] is None
    assert result["version"] is not None


@requires_macs3
def test_callpeak_finds_the_enriched_clusters(tmp_path: Path, pair_bam: Path) -> None:
    result = m.call_chip_peaks(
        [str(pair_bam)],
        str(tmp_path / "peaks"),
        file_format="BAM",
        genome_size=2_000_000,
        nomodel=True,
        extsize=_SPACING,
        name="planted",
    )
    assert result["status"] == "completed", result["log_tail"]
    peaks = [
        line
        for line in (tmp_path / "peaks" / "planted_peaks.narrowPeak").read_text().splitlines()
        if line
    ]
    assert len(peaks) >= 1
    assert all(line.startswith("chr1\t") for line in peaks)


@requires_macs3
def test_callpeak_builds_the_model_script_when_strand_pairs_exist(
    tmp_path: Path, pair_bam: Path
) -> None:
    result = m.call_chip_peaks(
        [str(pair_bam)],
        str(tmp_path / "modeled"),
        file_format="BAM",
        genome_size=2_000_000,
        name="modeled",
    )
    assert result["status"] == "completed", result["log_tail"]
    script = Path(result["model_script"])
    assert script.name == "modeled_model.r"
    assert "altd  <- c(" in script.read_text()


def test_callpeak_rejects_qvalue_and_pvalue_together(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        m.call_chip_peaks([str(READS)], str(tmp_path / "x"), qvalue=0.05, pvalue=0.01)


def test_callpeak_refuses_an_existing_output_directory(tmp_path: Path) -> None:
    existing = tmp_path / "exists"
    existing.mkdir()
    with pytest.raises(ConfigurationError):
        m.call_chip_peaks([str(READS)], str(existing), file_format="BAM", nomodel=True)


def test_callpeak_validates_genome_size_and_mfold(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        m.call_chip_peaks([str(READS)], str(tmp_path / "a"), genome_size="human")
    with pytest.raises(ConfigurationError):
        m.call_chip_peaks([str(READS)], str(tmp_path / "b"), mfold=(50, 5))
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()


# ---------------------------------------------------------------------------
# predictd and the fragment-model plot
# ---------------------------------------------------------------------------


@requires_macs3
def test_predictd_recovers_the_strand_spacing_and_writes_the_r_script(
    tmp_path: Path, pair_bam: Path
) -> None:
    out = tmp_path / "predictd"
    result = m.predict_fragment_size(
        [str(pair_bam)],
        str(out),
        file_format="BAM",
        genome_size=2_000_000,
        script_name="fragment.R",
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["fragment_length_source"] == "single_end_model"
    assert abs(result["fragment_length"] - _SPACING) <= 10
    script = out / "fragment.R"
    assert result["model_script"] == str(script)
    assert "xcorr <- c(" in script.read_text()


@requires_macs3
@requires_shared_bam
def test_predictd_on_the_forward_only_shared_bam_reports_no_model(tmp_path: Path) -> None:
    result = m.predict_fragment_size(
        [str(READS)], str(tmp_path / "shared"), file_format="BAM", genome_size=200000
    )
    assert result["status"] == "failed"
    assert result["model_script"] is None
    assert result["fragment_length"] is None
    assert "Can't find enough pairs" in result["log_tail"]


@requires_macs3
def test_predictd_paired_end_reports_the_average_insertion_length(tmp_path: Path) -> None:
    bedpe = tmp_path / "fragments.bedpe"
    bedpe.write_text("".join(f"chr1\t{start}\t{start + 180}\n" for start in range(1000, 9000, 100)))
    result = m.predict_fragment_size([str(bedpe)], str(tmp_path / "pe"), file_format="BEDPE")
    assert result["status"] == "completed", result["log_tail"]
    assert result["fragment_length_source"] == "paired_end_average"
    assert result["fragment_length"] == 180
    assert result["model_script"] is None


@requires_macs3
def test_plot_fragment_model_returns_the_model_data(tmp_path: Path, pair_bam: Path) -> None:
    result = m.plot_fragment_model(
        [str(pair_bam)],
        str(tmp_path / "plot"),
        file_format="BAM",
        genome_size=2_000_000,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert Path(result["script"]).name == "predictd_model.R"
    data = result["data"]
    assert (
        len(data["forward_percent"])
        == len(data["reverse_percent"])
        == len(data["distance_from_middle"])
    )
    assert sum(data["forward_percent"]) == pytest.approx(100.0, rel=1e-6)
    assert sum(data["reverse_percent"]) == pytest.approx(100.0, rel=1e-6)
    assert abs(result["fragment_length"] - _SPACING) <= 10
    if shutil.which("Rscript") is None:
        assert result["image"] is None
        assert "Rscript is not on PATH" in result["image_note"]
    else:
        assert result["image"] is not None and Path(result["image"]).is_file()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not on PATH")
@requires_macs3
def test_plot_fragment_model_renders_the_pdf_with_rscript(tmp_path: Path, pair_bam: Path) -> None:
    result = m.plot_fragment_model(
        [str(pair_bam)], str(tmp_path / "pdf"), file_format="BAM", genome_size=2_000_000
    )
    assert result["status"] == "completed", result["log_tail"]
    assert result["image"] is not None
    assert Path(result["image"]).read_bytes()[:4] == b"%PDF"


# ---------------------------------------------------------------------------
# program missing
# ---------------------------------------------------------------------------


def test_missing_macs3_raises_backend_error_with_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty_path = tmp_path / "empty-path"
    empty_path.mkdir()
    monkeypatch.setenv("PATH", str(empty_path))
    target = tmp_path / "predictd"
    with pytest.raises(BackendUnavailableError) as raised:
        m.predict_fragment_size(
            [str(READS)],
            str(target),
            file_format="BAM",
            macs3_path=str(tmp_path / "no-such-macs3"),
        )
    assert "macs3" in str(raised.value)
    assert raised.value.hint is not None and "pip install macs3" in raised.value.hint
    assert not target.exists(), "no output directory may be created when the program is missing"


# ---------------------------------------------------------------------------
# bedGraph and track operations
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def shared_tracks(tmp_path_factory: pytest.TempPathFactory) -> dict[str, str]:
    """Treatment pileup and control lambda tracks from the shared reads, and the ppois and qpois score tracks."""
    base = tmp_path_factory.mktemp("tracks")
    called = m.call_chip_peaks(
        [str(READS)],
        str(base / "callpeak"),
        file_format="BAM",
        genome_size=200_000,
        nomodel=True,
        extsize=200,
        name="s",
        save_bedgraph=True,
    )
    assert called["status"] == "completed", called["log_tail"]
    scored = m.compare_bedgraph_tracks(
        str(base / "callpeak" / "s_treat_pileup.bdg"),
        str(base / "callpeak" / "s_control_lambda.bdg"),
        str(base / "scores"),
        name="cmp",
        methods=("ppois", "qpois"),
    )
    assert scored["status"] == "completed", scored["log_tail"]
    return {
        "treat": str(base / "callpeak" / "s_treat_pileup.bdg"),
        "control": str(base / "callpeak" / "s_control_lambda.bdg"),
        "ppois": str(base / "scores" / "cmp_ppois.bdg"),
        "qpois": str(base / "scores" / "cmp_qpois.bdg"),
    }


def _rows(path: Path) -> list[list[str]]:
    return [
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("track")
    ]


@requires_macs3
@requires_shared_bam
def test_bdgcmp_writes_one_score_track_per_method(shared_tracks: dict[str, str]) -> None:
    for key in ("ppois", "qpois"):
        rows = _rows(Path(shared_tracks[key]))
        assert rows and all(len(row) == 4 for row in rows)


@requires_macs3
def test_bdgopt_multiply_scales_scores_and_keeps_intervals(
    tmp_path: Path, shared_tracks: dict[str, str]
) -> None:
    result = m.operate_bedgraph_scores(
        shared_tracks["ppois"], str(tmp_path / "half.bdg"), method="multiply", extra_params=[0.5]
    )
    assert result["status"] == "completed", result["log_tail"]
    before = _rows(Path(shared_tracks["ppois"]))
    after = _rows(tmp_path / "half.bdg")
    assert len(after) == len(before)
    for scaled, original in zip(after, before, strict=True):
        assert scaled[:3] == original[:3]
        assert float(scaled[3]) == pytest.approx(float(original[3]) * 0.5, abs=2e-5)


@requires_macs3
def test_bdgopt_p2q_converts_without_a_parameter(
    tmp_path: Path, shared_tracks: dict[str, str]
) -> None:
    result = m.operate_bedgraph_scores(
        shared_tracks["ppois"], str(tmp_path / "q.bdg"), method="p2q"
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _rows(tmp_path / "q.bdg")


@requires_macs3
def test_bdgpeakcall_and_bdgbroadcall_write_their_tables(
    tmp_path: Path, shared_tracks: dict[str, str]
) -> None:
    peaks = m.call_peaks_from_bedgraph(
        shared_tracks["ppois"],
        str(tmp_path / "peaks.narrowPeak"),
        cutoff=2.0,
        min_length=50,
        max_gap=30,
        call_summits=True,
        trackline=False,
    )
    assert peaks["status"] == "completed", peaks["log_tail"]
    assert not (tmp_path / "peaks.narrowPeak").read_text().startswith("track")
    broad = m.call_broad_peaks_from_bedgraph(
        shared_tracks["ppois"],
        str(tmp_path / "broad.gappedPeak"),
        cutoff_peak=2.0,
        cutoff_link=1.0,
        min_length=50,
    )
    assert broad["status"] == "completed", broad["log_tail"]
    assert {len(row) for row in _rows(tmp_path / "broad.gappedPeak")} <= {15}


@requires_macs3
def test_bdgdiff_writes_condition_specific_and_common_regions(
    tmp_path: Path, shared_tracks: dict[str, str]
) -> None:
    outputs = [tmp_path / name for name in ("cond1.bed", "cond2.bed", "common.bed")]
    result = m.differential_bedgraph_peaks(
        shared_tracks["treat"],
        shared_tracks["treat"],
        shared_tracks["control"],
        shared_tracks["control"],
        str(outputs[0]),
        str(outputs[1]),
        str(outputs[2]),
        cutoff=2.0,
        min_length=50,
        max_gap=100,
    )
    assert result["status"] == "completed", result["log_tail"]
    assert all(path.is_file() for path in outputs)


@requires_macs3
def test_cmbreps_combines_two_replicate_tracks(
    tmp_path: Path, shared_tracks: dict[str, str]
) -> None:
    second = m.compare_bedgraph_tracks(
        shared_tracks["treat"],
        shared_tracks["control"],
        str(tmp_path / "rep2"),
        name="rep2",
        methods=("ppois",),
    )
    assert second["status"] == "completed", second["log_tail"]
    result = m.combine_replicate_bedgraphs(
        [shared_tracks["ppois"], second["score_files"][0]],
        str(tmp_path / "combined.bdg"),
        method="fisher",
    )
    assert result["status"] == "completed", result["log_tail"]
    assert _rows(tmp_path / "combined.bdg")


def test_bedgraph_facades_validate_before_running(tmp_path: Path) -> None:
    track = tmp_path / "track.bdg"
    track.write_text("chr1\t0\t10\t1.0\n")
    with pytest.raises(ConfigurationError):
        m.operate_bedgraph_scores(str(track), str(tmp_path / "a.bdg"), method="multiply")
    with pytest.raises(ConfigurationError):
        m.operate_bedgraph_scores(
            str(track), str(tmp_path / "b.bdg"), method="p2q", extra_params=[1.0]
        )
    with pytest.raises(ConfigurationError):
        m.combine_replicate_bedgraphs([str(track)], str(tmp_path / "c.bdg"))
    with pytest.raises(ConfigurationError):
        m.compare_bedgraph_tracks(
            str(track), str(track), str(tmp_path / "scores"), methods=("bogus",)
        )
    assert not (tmp_path / "a.bdg").exists() and not (tmp_path / "scores").exists()


# ---------------------------------------------------------------------------
# alignment operations
# ---------------------------------------------------------------------------


@requires_macs3
@requires_shared_bam
def test_pileup_extends_each_read_to_the_extension_size(tmp_path: Path) -> None:
    out = tmp_path / "pileup.bdg"
    result = m.pileup_reads_to_bedgraph([str(READS)], str(out), file_format="BAM", extsize=200)
    assert result["status"] == "completed", result["log_tail"]
    rows = _rows(out)
    total = sum((int(end) - int(start)) * float(value) for _, start, end, value in rows)
    assert total == pytest.approx(5000 * 200, rel=0.01)


@requires_macs3
@requires_shared_bam
def test_filterdup_keeps_all_reads_or_drops_duplicates(tmp_path: Path) -> None:
    everything = m.filter_duplicate_reads(
        [str(READS)], str(tmp_path / "all.bed"), file_format="BAM", keep_duplicates="all"
    )
    assert everything["status"] == "completed", everything["log_tail"]
    assert len(_rows(tmp_path / "all.bed")) == 5000
    one = m.filter_duplicate_reads(
        [str(READS)], str(tmp_path / "one.bed"), file_format="BAM", keep_duplicates=1
    )
    assert one["status"] == "completed", one["log_tail"]
    assert 4000 < len(_rows(tmp_path / "one.bed")) < 5000


@requires_macs3
@requires_shared_bam
def test_randsample_by_number_and_by_percentage(tmp_path: Path) -> None:
    by_number = m.random_sample_reads(
        [str(READS)], str(tmp_path / "n.bed"), number=100, seed=7, file_format="BAM"
    )
    assert by_number["status"] == "completed", by_number["log_tail"]
    assert 90 <= len(_rows(tmp_path / "n.bed")) <= 110
    by_percent = m.random_sample_reads(
        [str(READS)], str(tmp_path / "p.bed"), percentage=10.0, file_format="BAM"
    )
    assert by_percent["status"] == "completed", by_percent["log_tail"]
    assert 400 <= len(_rows(tmp_path / "p.bed")) <= 600


def test_randsample_needs_exactly_one_amount(tmp_path: Path) -> None:
    track = tmp_path / "reads.bed"
    track.write_text("chr1\t0\t10\n")
    with pytest.raises(ConfigurationError):
        m.random_sample_reads([str(track)], str(tmp_path / "x.bed"))
    with pytest.raises(ConfigurationError):
        m.random_sample_reads([str(track)], str(tmp_path / "y.bed"), number=10, percentage=10.0)


@requires_macs3
@requires_shared_bam
def test_refinepeak_places_summits_near_the_candidate(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.bed"
    candidates.write_text("chr1\t1000\t1300\tp1\t10\t.\n")
    result = m.refine_peak_summits(
        str(candidates), [str(READS)], str(tmp_path / "refined.bed"), file_format="BAM"
    )
    assert result["status"] == "completed", result["log_tail"]
    rows = _rows(tmp_path / "refined.bed")
    assert rows and all(row[0] == "chr1" and int(row[2]) - int(row[1]) == 1 for row in rows)


def _with_md_tags(source: Path, target: Path) -> Path:
    """A copy of a BAM with MD tags computed against the shared reference (samtools calmd equivalent for M-only reads)."""
    pysam = pytest.importorskip("pysam")
    reference = pysam.FastaFile(str(DATA / "genome.fa"))

    def md_for(read: str, ref: str) -> str:
        parts, run = [], 0
        for base, expected in zip(read.upper(), ref.upper(), strict=True):
            if base == expected:
                run += 1
            else:
                parts += [str(run), expected]
                run = 0
        return "".join(parts) + str(run)

    with (
        pysam.AlignmentFile(str(source)) as reader,
        pysam.AlignmentFile(str(target), "wb", template=reader) as writer,
    ):
        for read in reader:
            if not read.is_unmapped:
                ref = reference.fetch(read.reference_name, read.reference_start, read.reference_end)
                read.set_tag("MD", md_for(read.query_sequence, ref), value_type="Z")
            writer.write(read)
    pysam.index(str(target))
    return target


@requires_macs3
@requires_shared_bam
def test_callvar_writes_a_vcf_for_reads_with_md_tags(tmp_path: Path) -> None:
    bam = _with_md_tags(READS, tmp_path / "reads.md.bam")
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t1000\t3000\tpk\t10\n")
    result = m.call_variants_in_peaks(str(peaks), str(bam), str(tmp_path / "calls.vcf"))
    assert result["status"] == "completed", result["log_tail"]
    text = (tmp_path / "calls.vcf").read_text()
    assert "##contig=<ID=chr1,length=120000" in text and "#CHROM\tPOS" in text
    assert result["warnings"] == []


@requires_macs3
@requires_shared_bam
def test_callvar_warns_when_reads_lack_md_tags(tmp_path: Path) -> None:
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t1000\t3000\tpk\t10\n")
    result = m.call_variants_in_peaks(str(peaks), str(READS), str(tmp_path / "plain.vcf"))
    assert result["status"] == "completed", result["log_tail"]
    assert result["warnings"] and "MD tags" in result["warnings"][0]


def test_callvar_validates_fermi_mode(tmp_path: Path) -> None:
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t1000\t3000\tpk\t10\n")
    bam = tmp_path / "reads.bam"
    bam.write_bytes(b"")
    with pytest.raises(ConfigurationError):
        m.call_variants_in_peaks(str(peaks), str(bam), str(tmp_path / "x.vcf"), fermi="maybe")


def _open_region_fragments(path: Path) -> Path:
    """Paired-end style fragments on chr1: 35 open regions of short fragments over a nucleosomal background."""
    rng = random.Random(3)
    rows = []
    for centre in range(4000, 4000 + 35 * 3000, 3000):
        for _ in range(400):
            start = centre + rng.randint(-150, 150)
            rows.append((start, start + rng.randint(40, 110)))
    for _ in range(9000):
        start = rng.randint(0, 119_000)
        rows.append((start, start + rng.randint(160, 320)))
    rows.sort()
    path.write_text("".join(f"chr1\t{start}\t{end}\n" for start, end in rows))
    return path


@requires_macs3
def test_hmmratac_calls_accessible_regions(tmp_path: Path) -> None:
    fragments = _open_region_fragments(tmp_path / "fragments.bedpe")
    out = tmp_path / "hmm"
    result = m.call_atac_peaks_hmmratac(
        [str(fragments)],
        str(out),
        name="toy",
        file_format="BEDPE",
        lower_cutoff=2,
        upper_cutoff=6,
        prescan_cutoff=1.2,
        save_states=True,
    )
    assert result["status"] == "completed", result["log_tail"]
    regions = _rows(out / "toy_accessible_regions.narrowPeak")
    assert regions and all(row[0] == "chr1" for row in regions)
    assert (out / "toy_cutoff_analysis.tsv").is_file()
    assert (out / "toy_model.json").is_file()
    assert (out / "toy_states.bed").is_file()


@requires_macs3
def test_hmmratac_default_cutoffs_stop_when_no_training_region_qualifies(tmp_path: Path) -> None:
    fragments = _open_region_fragments(tmp_path / "fragments.bedpe")
    result = m.call_atac_peaks_hmmratac(
        [str(fragments)], str(tmp_path / "default"), name="d", file_format="BEDPE"
    )
    assert result["status"] == "failed"
    assert "Not enough training regions" in result["log_tail"]


def test_missing_macs3_for_alignment_and_track_facades(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty_path = tmp_path / "empty-path"
    empty_path.mkdir()
    monkeypatch.setenv("PATH", str(empty_path))
    reads = tmp_path / "reads.bam"
    reads.write_bytes(b"")
    target = tmp_path / "pileup.bdg"
    with pytest.raises(BackendUnavailableError) as raised:
        m.pileup_reads_to_bedgraph(
            [str(reads)], str(target), file_format="BAM", macs3_path=str(tmp_path / "no-such-macs3")
        )
    assert raised.value.hint is not None and "pip install macs3" in raised.value.hint
    assert not target.exists()


def test_existing_output_files_are_never_overwritten(tmp_path: Path) -> None:
    reads = tmp_path / "reads.bam"
    reads.write_bytes(b"")
    existing = tmp_path / "pileup.bdg"
    existing.write_text("keep me\n")
    with pytest.raises(ConfigurationError):
        m.pileup_reads_to_bedgraph([str(reads)], str(existing), file_format="BAM")
    assert existing.read_text() == "keep me\n"
