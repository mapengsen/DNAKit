"""Tests for dnakit.annotation.pyranges1_bam_gap_functions (BAM alignments with pysam)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")
pyranges1 = pytest.importorskip("pyranges1")

from dnakit.annotation import pyranges1_bam_gap_functions as m
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)
SHARED_BAM = DATA / "reads.bam"
HEADER = {
    "HD": {"VN": "1.6", "SO": "unsorted"},
    "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 500}],
}


def _segment(header: pysam.AlignmentHeader, spec: dict) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment(header)
    seg.query_name = spec["name"]
    seg.flag = spec["flag"]
    seg.reference_id = spec.get("ref", -1)
    seg.reference_start = spec.get("pos", -1)
    seg.mapping_quality = spec.get("mapq", 60)
    if spec.get("cigar") is not None:
        seg.cigarstring = spec["cigar"]
    seg.query_sequence = spec["seq"]
    seg.query_qualities = pysam.qualitystring_to_array("I" * len(spec["seq"]))
    return seg


def _write_bam(path: Path) -> Path:
    """Synthetic BAM: one read per filter (strand, soft clips, deletion, unmapped, duplicate)."""
    reads = [
        {"name": "fwd", "flag": 0, "ref": 0, "pos": 99, "cigar": "10M", "seq": "ACGTACGTAC"},
        {
            "name": "rev",
            "flag": 16,
            "ref": 0,
            "pos": 199,
            "cigar": "5S10M2D5M",
            "seq": "ACGTA" + "C" * 10 + "GTACG",
        },
        {"name": "unmapped", "flag": 4, "seq": "ACGTACGTAC"},
        {"name": "dup", "flag": 1024, "ref": 0, "pos": 300, "cigar": "8M", "seq": "ACGTACGT"},
        {
            "name": "lowmq",
            "flag": 0,
            "ref": 1,
            "pos": 10,
            "cigar": "6M",
            "seq": "ACGTAC",
            "mapq": 3,
        },
        {"name": "chr2_fwd", "flag": 0, "ref": 1, "pos": 50, "cigar": "4M", "seq": "GGCC"},
        {"name": "secondary", "flag": 256, "ref": 0, "pos": 400, "cigar": "4M", "seq": "TTAA"},
        {"name": "placed_unmapped", "flag": 4, "ref": 0, "pos": 500, "seq": "ACGT"},
    ]
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for spec in reads:
            out.write(_segment(out.header, spec))
    return path


@pytest.fixture()
def synthetic(tmp_path: Path) -> Path:
    return _write_bam(tmp_path / "synthetic.bam")


def test_default_filters_and_interval_convention(synthetic: Path) -> None:
    result = m.read_pyranges1_bam(str(synthetic), verbose=False)
    assert result["status"] == "completed"
    assert result["columns"] == ["Chromosome", "Start", "End", "Strand", "Flag"]
    # unmapped (4), duplicate (1024) and the unplaced/unmapped reads are dropped by the default
    # filter 1540
    assert [r["Flag"] for r in result["records"]] == [0, 16, 0, 0, 256]
    assert result["records"][0] == {
        "Chromosome": "chr1",
        "Start": 99,
        "End": 109,
        "Strand": "+",
        "Flag": 0,
    }
    # 5S10M2D5M: soft clips do not count, deletions do: reference span 10 + 2 + 5 = 17
    assert result["records"][1] == {
        "Chromosome": "chr1",
        "Start": 199,
        "End": 216,
        "Strand": "-",
        "Flag": 16,
    }
    assert result["record_count"] == 5 and result["truncated"] is False
    assert result["provenance"]["backend"].startswith("pysam ")
    assert result["provenance"]["package"].startswith("pyranges1 ")


def test_mapq_required_and_filter_flags(synthetic: Path) -> None:
    by_mapq = m.read_pyranges1_bam(str(synthetic), mapq=10, verbose=False)
    assert [r["Chromosome"] for r in by_mapq["records"]] == ["chr1", "chr1", "chr2", "chr1"]
    assert by_mapq["record_count"] == 4
    reverse_only = m.read_pyranges1_bam(str(synthetic), required_flag=16, verbose=False)
    assert [r["Start"] for r in reverse_only["records"]] == [199]
    unfiltered = m.read_pyranges1_bam(str(synthetic), filter_flag=0, verbose=False)
    # duplicates are kept when filter_flag is 0; reads without a reference span are never reported
    assert unfiltered["record_count"] == 6
    assert 1024 in {r["Flag"] for r in unfiltered["records"]}
    assert all(r["End"] > r["Start"] for r in unfiltered["records"])


def test_full_columns_add_query_fields(synthetic: Path) -> None:
    result = m.read_pyranges1_bam(str(synthetic), sparse=False, verbose=False)
    assert result["columns"] == [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "Flag",
        "QueryStart",
        "QueryEnd",
        "QuerySequence",
        "Name",
        "Cigar",
        "Quality",
    ]
    rev = result["records"][1]
    assert rev["Name"] == "rev"
    assert rev["Cigar"] == "5S10M2D5M"
    assert (rev["QueryStart"], rev["QueryEnd"]) == (5, 20)
    assert rev["QuerySequence"] == "ACGTA" + "C" * 10 + "GTACG"
    assert rev["Quality"] == "I" * 20


def test_max_rows_truncates_but_keeps_the_count(synthetic: Path) -> None:
    result = m.read_pyranges1_bam(str(synthetic), max_rows=2, verbose=False)
    assert len(result["records"]) == 2
    assert result["record_count"] == 5
    assert result["truncated"] is True


def test_shared_reads_bam_is_read_completely() -> None:
    result = m.read_pyranges1_bam(str(SHARED_BAM), max_rows=3, verbose=False)
    assert result["record_count"] == 5000
    assert result["truncated"] is True
    assert all(r["End"] > r["Start"] >= 0 and r["Strand"] in {"+", "-"} for r in result["records"])


def test_input_validation(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_bam(str(tmp_path / "missing.bam"), verbose=False)
    text = tmp_path / "not_a_bam.bam"
    text.write_text("not a bam file\n", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_bam(str(text), verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_bam(str(SHARED_BAM), mapq=300, verbose=False)
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_bam(str(SHARED_BAM), sparse="yes", verbose=False)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError):
        m.read_pyranges1_bam(str(SHARED_BAM), max_rows=0, verbose=False)


def test_missing_pysam_raises_with_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pysam", None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        m.read_pyranges1_bam(str(SHARED_BAM), verbose=False)
    assert "bio-reads" in str(excinfo.value.hint or "")


def test_missing_pyranges1_raises_with_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pyranges1", None)
    with pytest.raises(BackendUnavailableError) as excinfo:
        m.read_pyranges1_bam(str(SHARED_BAM), verbose=False)
    assert "bio-pyranges1" in str(excinfo.value.hint or "")


def test_upstream_read_bam_exits_without_bamread() -> None:
    """Reproduction behind the pysam implementation: pyranges1.read_bam exits without bamread."""
    if importlib.util.find_spec("bamread") is not None:
        pytest.skip("bamread is installed; upstream read_bam would run")
    assert pyranges1.readers.read_bam is pyranges1.read_bam
    with pytest.raises(SystemExit) as excinfo:
        pyranges1.read_bam(str(SHARED_BAM))
    assert excinfo.value.code == 1


def test_ensure_pyranges_aliases_are_one_function() -> None:
    """Evidence for category (d): every ensure_pyranges alias is the same object (PyRanges(df) plus
    isinstance).
    """
    import pyranges1.cli as cli
    import pyranges1.core.empty as empty
    import pyranges1.core.example_data as example_data
    import pyranges1.core.loci_getter as loci_getter
    import pyranges1.core.multioverlap as multioverlap
    import pyranges1.core.out as out
    import pyranges1.core.pyranges_main as pyranges_main
    import pyranges1.core.random as random_module
    import pyranges1.ext.orfs as ext_orfs
    import pyranges1.ext.stats as ext_stats
    import pyranges1.methods.map_to_global as map_to_global
    import pyranges1.methods.split as split_module
    import pyranges1.methods.tile_genome as tile_genome
    import pyranges1.readers as readers
    from pyranges1.core import pyranges_helpers

    modules = [
        cli,
        empty,
        example_data,
        loci_getter,
        multioverlap,
        out,
        pyranges_helpers,
        pyranges_main,
        random_module,
        ext_orfs,
        ext_stats,
        map_to_global,
        split_module,
        tile_genome,
        readers,
    ]
    assert all(mod.ensure_pyranges is pyranges_helpers.ensure_pyranges for mod in modules)
    assert not hasattr(pyranges1, "ensure_pyranges")
