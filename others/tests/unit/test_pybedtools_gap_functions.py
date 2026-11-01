"""Tests for pybedtools_gap_functions: chromosome-size lookups, UCSC and samtools facades, BAM
coverage, Venn and heat-map helpers, feature colouring.

Chromosome-size lookups are tested offline with pybedtools' built-in registry and with fakes that
replace pybedtools or genomepy callables; the socket guard below blocks real connections. Tests that
need bedtools (installed here) or pysam skip with a reason when they are missing. Tests of the UCSC
and samtools facades use stand-in executables that record their arguments under tmp_path, so the
argument lists are checked without the real programs; tests that need the real programs skip when
they are not on PATH.
"""

from __future__ import annotations

import shutil
import socket
import sys
from pathlib import Path

import pytest

from dnakit.bioinformatics import pybedtools_gap_functions as module
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def _blocked(*args, **kwargs):
        raise RuntimeError("network disabled in unit tests")

    monkeypatch.setattr(socket.socket, "connect", _blocked)


_REQUIRES_BEDTOOLS = pytest.mark.skipif(
    shutil.which("bedtools") is None, reason="bedtools is not installed"
)
_REQUIRES_SAMTOOLS = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools is not installed"
)
_REQUIRES_RSCRIPT = pytest.mark.skipif(
    shutil.which("Rscript") is None, reason="Rscript is not installed"
)
_REQUIRES_BEDTOBIGBED = pytest.mark.skipif(
    shutil.which("bedToBigBed") is None, reason="UCSC bedToBigBed is not installed"
)
_REQUIRES_BEDGRAPHTOBIGWIG = pytest.mark.skipif(
    shutil.which("bedGraphToBigWig") is None, reason="UCSC bedGraphToBigWig is not installed"
)

# Writes its stand-in output to the last argument that does not start with '-' and records the
# argument list.
_STAND_IN_HEAD = (
    'out=""\n'
    'for a in "$@"; do case "$a" in -*) ;; *) out="$a";; esac; done\n'
    'printf "%s\\n" "$@" > "$(dirname "$0")/argv.txt"\n'
    'printf "stand-in output\\n" > "$out"\n'
)


def _write_stand_in(directory: Path, name: str, body: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / name
    script.write_text("#!/bin/sh\n" + body)
    script.chmod(0o755)
    return script


def _argv(directory: Path) -> list[str]:
    return (directory / "argv.txt").read_text().splitlines()


def _empty_path(monkeypatch, tmp_path: Path) -> None:
    empty = tmp_path / "empty_bin"
    empty.mkdir(exist_ok=True)
    monkeypatch.setenv("PATH", str(empty))


def _write_bam(path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
        for index, position in enumerate([10, 20, 20, 500]):
            read = pysam.AlignedSegment()
            read.query_name = f"r{index}"
            read.query_sequence = "ACGTACGTAC"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = position
            read.mapping_quality = 60
            read.cigarstring = "10M"
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            handle.write(read)


def _venn_example(tmp_path: Path) -> tuple[Path, Path]:
    """The two-file example from pybedtools' cleaned_intersect documentation."""
    x = tmp_path / "x.bed"
    x.write_text("chr1\t1\t100\tf1\t0\t+\nchr1\t500\t6000\tf2\t0\t+\n")
    y = tmp_path / "y.bed"
    y.write_text("chr1\t50\t100\tg1\t0\t+\nchr1\t80\t200\tg2\t0\t+\nchr9\t777\t888\tg3\t0\t+\n")
    return x, y


def test_public_names_are_exported_and_callable():
    for name in module.__all__:
        assert callable(getattr(module, name)), name
    assert len(module.__all__) == 15
    assert "truncator" not in module.__all__
    assert "binary_summary" not in module.__all__


# ---- chromosome sizes
# -------------------------------------------------------------------------------------


def test_chromsizes_from_registry_is_offline_and_writes_file(tmp_path):
    pytest.importorskip("pybedtools")
    out = tmp_path / "dm3.chrom.sizes"
    result = module.chromsizes("dm3", output_path=str(out), verbose=False)
    assert result["status"] == "completed"
    assert result["data"]["chr2L"] == 23011544
    assert result["service"].startswith("pybedtools built-in genome registry")
    assert result["retrieved_at"].endswith("+00:00")
    assert "chr2L\t23011544" in out.read_text().splitlines()


def test_chromsizes_lookup_uses_pybedtools_with_fake(monkeypatch):
    pybedtools = pytest.importorskip("pybedtools")
    seen = []

    def fake(name):
        seen.append(name)
        return {"chrS": (0, 42), "chrT": (0, 7)}

    monkeypatch.setattr(pybedtools, "chromsizes", fake)
    result = module.chromsizes("fakeGenomeX", verbose=False)
    assert seen == ["fakeGenomeX"]
    assert result["data"] == {"chrS": 42, "chrT": 7}
    assert "network" in result["service"]


def test_chromsizes_maps_lookup_failure(monkeypatch):
    pybedtools = pytest.importorskip("pybedtools")

    def failing(name):
        raise ValueError("no sizes anywhere")

    monkeypatch.setattr(pybedtools, "chromsizes", failing)
    with pytest.raises(BackendExecutionError, match="no sizes anywhere"):
        module.chromsizes("fakeGenomeX", verbose=False)


def test_get_chromsizes_from_ucsc_passes_arguments_with_fake(monkeypatch, tmp_path):
    pybedtools = pytest.importorskip("pybedtools")
    seen = {}

    def fake(genome, saveas=None, timeout=None, host_url=None, **kwargs):
        seen.update(genome=genome, saveas=saveas, timeout=timeout, host_url=host_url)
        if saveas is not None:
            Path(saveas).write_text("chr1\t10\n")
        return {"chr1": (0, 10)}

    monkeypatch.setattr(pybedtools, "get_chromsizes_from_ucsc", fake)
    out = tmp_path / "ucsc.sizes"
    result = module.get_chromsizes_from_ucsc(
        "hg38", timeout=5.0, host_url="ucsc.example.org", output_path=str(out), verbose=False
    )
    assert seen == {
        "genome": "hg38",
        "saveas": str(out),
        "timeout": 5.0,
        "host_url": "ucsc.example.org",
    }
    assert result["data"] == {"chr1": 10}
    assert result["request"]["host_url"] == "ucsc.example.org"
    assert result["output_path"] == str(out)
    assert out.is_file()


def test_get_chromsizes_from_ucsc_maps_errors_and_rejects_bad_input(monkeypatch):
    pybedtools = pytest.importorskip("pybedtools")

    def failing(genome, **kwargs):
        raise OSError("Can't find mysql at path mysql")

    monkeypatch.setattr(pybedtools, "get_chromsizes_from_ucsc", failing)
    with pytest.raises(BackendExecutionError, match="mysql"):
        module.get_chromsizes_from_ucsc("hg38", verbose=False)
    with pytest.raises(ConfigurationError):
        module.get_chromsizes_from_ucsc("hg38", host_url="bad host;rm", verbose=False)
    with pytest.raises(ConfigurationError):
        module.get_chromsizes_from_ucsc("hg38", timeout=0, verbose=False)
    with pytest.raises(ConfigurationError):
        module.chromsizes("../etc", verbose=False)


def test_get_chromsizes_from_genomepy_reads_installed_sizes(monkeypatch, tmp_path):
    genomepy = pytest.importorskip("genomepy")
    sizes_file = tmp_path / "installed.sizes"
    sizes_file.write_text("chr1\t500\nchr2\t900\n")

    class FakeGenome:
        def __init__(self, name):
            self.name = name
            self.sizes_file = str(sizes_file)

    monkeypatch.setattr(genomepy, "Genome", FakeGenome)
    result = module.get_chromsizes_from_genomepy("fakeInstalled", verbose=False)
    assert result["status"] == "completed"
    assert result["data"] == {"chr1": 500, "chr2": 900}
    assert result["service"].startswith("genomepy")


def test_get_chromsizes_from_genomepy_reports_not_found(monkeypatch):
    genomepy = pytest.importorskip("genomepy")

    class MissingGenome:
        def __init__(self, name):
            raise FileNotFoundError(name)

    monkeypatch.setattr(genomepy, "Genome", MissingGenome)
    result = module.get_chromsizes_from_genomepy("notInstalled", verbose=False)
    assert result["status"] == "not_found"
    assert result["data"] is None


def test_genomepy_missing_reports_install_hint(monkeypatch):
    monkeypatch.setitem(sys.modules, "genomepy", None)
    with pytest.raises(BackendUnavailableError) as caught:
        module.get_chromsizes_from_genomepy("hg38", verbose=False)
    assert "dnakit[bio-genomepy]" in (caught.value.hint or "")


# ---- UCSC and samtools facades: missing program
# ------------------------------------------------------------


def test_external_facades_report_missing_programs(monkeypatch, tmp_path):
    pytest.importorskip("pybedtools")
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\tn\t5\t+\n")
    bedgraph = tmp_path / "in.bedGraph"
    bedgraph.write_text("chr1\t0\t10\t1\n")
    bigwig = tmp_path / "in.bw"
    bigwig.write_bytes(b"x")
    bigbed = tmp_path / "in.bb"
    bigbed.write_bytes(b"x")
    wig = tmp_path / "in.wig"
    wig.write_text("fixedStep chrom=chr1 start=1 step=1\n1\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"x")
    _empty_path(monkeypatch, tmp_path)
    cases = [
        (
            lambda: module.bigbed(str(bed), str(sizes), str(tmp_path / "a.bb"), verbose=False),
            "ucsc-bedtobigbed",
        ),
        (
            lambda: module.bigbed_to_bed(str(bigbed), str(tmp_path / "a.bed"), verbose=False),
            "ucsc-bigbedtobed",
        ),
        (
            lambda: module.bedgraph_to_bigwig(
                str(bedgraph), str(sizes), str(tmp_path / "a.bw"), verbose=False
            ),
            "ucsc-bedgraphtobigwig",
        ),
        (
            lambda: module.bigwig_to_bedgraph(
                str(bigwig), str(tmp_path / "a.bedGraph"), verbose=False
            ),
            "ucsc-bigwigtobedgraph",
        ),
        (
            lambda: module.wig_to_bigwig(
                str(wig), str(sizes), str(tmp_path / "a.bw"), verbose=False
            ),
            "ucsc-wigtobigwig",
        ),
        (lambda: module.mapped_read_count(str(bam), verbose=False), "samtools"),
        (
            lambda: module.bam_to_bigwig(
                str(bam), str(sizes), str(tmp_path / "a.bw"), verbose=False
            ),
            "bedtools",
        ),
    ]
    for call, marker in cases:
        with pytest.raises(BackendUnavailableError) as caught:
            call()
        assert marker in (caught.value.hint or "") or marker in str(caught.value), marker


# ---- UCSC and samtools facades: argument lists with stand-in programs
# -----------------------------------


def test_bigbed_builds_fixed_arguments_with_stand_in(tmp_path):
    pytest.importorskip("pybedtools")
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\tn\t5\t+\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    autosql = tmp_path / "schema.as"
    autosql.write_text("table x\n")
    tool_dir = tmp_path / "tools"
    tool = _write_stand_in(tool_dir, "bedToBigBed", _STAND_IN_HEAD)
    out = tmp_path / "out.bb"
    result = module.bigbed(
        str(bed),
        str(sizes),
        str(out),
        unc=True,
        autosql=str(autosql),
        bedtobigbed_path=str(tool),
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["outputs"] == [str(out)]
    assert _argv(tool_dir) == [
        str(bed.resolve()),
        str(sizes.resolve()),
        str(out.resolve()),
        "-blockSize=256",
        "-itemsPerSlot=512",
        "-type=bed6",
        "-unc",
        f"-as={autosql.resolve()}",
    ]


def test_bigbed_rejects_bad_type_and_existing_output(tmp_path):
    pytest.importorskip("pybedtools")
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    existing = tmp_path / "exists.bb"
    existing.write_text("x")
    with pytest.raises(ConfigurationError):
        module.bigbed(str(bed), str(sizes), str(existing), verbose=False)
    with pytest.raises(ConfigurationError):
        module.bigbed(
            str(bed), str(sizes), str(tmp_path / "n.bb"), bed_type="bed6 -x", verbose=False
        )


def test_bigbed_to_bed_uses_equals_style_options(tmp_path):
    pytest.importorskip("pybedtools")
    source = tmp_path / "in.bb"
    source.write_bytes(b"x")
    tool_dir = tmp_path / "tools"
    tool = _write_stand_in(tool_dir, "bigBedToBed", _STAND_IN_HEAD)
    out = tmp_path / "out.bed"
    result = module.bigbed_to_bed(
        str(source),
        str(out),
        chrom="chr1",
        start=5,
        end=50,
        max_items=10,
        bigbedtobed_path=str(tool),
        verbose=False,
    )
    assert result["status"] == "completed"
    assert _argv(tool_dir) == [
        "-chrom=chr1",
        "-start=5",
        "-end=50",
        "-maxItems=10",
        str(source.resolve()),
        str(out.resolve()),
    ]
    with pytest.raises(ConfigurationError):
        module.bigbed_to_bed(str(source), str(tmp_path / "x.bed"), start=5, verbose=False)


def test_bigwig_facades_build_fixed_arguments_with_stand_ins(tmp_path):
    pytest.importorskip("pybedtools")
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    wig = tmp_path / "in.wig"
    wig.write_text("fixedStep chrom=chr1 start=1 step=1\n1\n")
    bigwig = tmp_path / "in.bw"
    bigwig.write_bytes(b"x")

    bg_dir = tmp_path / "bg_tool"
    bg_tool = _write_stand_in(bg_dir, "bedGraphToBigWig", _STAND_IN_HEAD)
    module.bedgraph_to_bigwig(
        str(graph),
        str(sizes),
        str(tmp_path / "a.bw"),
        bedgraphtobigwig_path=str(bg_tool),
        verbose=False,
    )
    assert _argv(bg_dir) == [
        str(graph.resolve()),
        str(sizes.resolve()),
        str((tmp_path / "a.bw").resolve()),
    ]

    w_dir = tmp_path / "wig_tool"
    w_tool = _write_stand_in(w_dir, "wigToBigWig", _STAND_IN_HEAD)
    module.wig_to_bigwig(
        str(wig), str(sizes), str(tmp_path / "b.bw"), wigtobigwig_path=str(w_tool), verbose=False
    )
    assert _argv(w_dir) == [
        str(wig.resolve()),
        str(sizes.resolve()),
        str((tmp_path / "b.bw").resolve()),
    ]

    b2g_dir = tmp_path / "b2g_tool"
    b2g_tool = _write_stand_in(b2g_dir, "bigWigToBedGraph", _STAND_IN_HEAD)
    module.bigwig_to_bedgraph(
        str(bigwig),
        str(tmp_path / "c.bedGraph"),
        chrom="chr1",
        udc_dir="/cache",
        bigwigtobedgraph_path=str(b2g_tool),
        verbose=False,
    )
    assert _argv(b2g_dir) == [
        "-chrom=chr1",
        "-udcDir=/cache",
        str(bigwig.resolve()),
        str((tmp_path / "c.bedGraph").resolve()),
    ]


def test_failed_program_run_is_reported_not_raised(tmp_path):
    pytest.importorskip("pybedtools")
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    tool = _write_stand_in(tmp_path / "bad", "bedGraphToBigWig", 'echo "bad input" >&2\nexit 3\n')
    result = module.bedgraph_to_bigwig(
        str(graph),
        str(sizes),
        str(tmp_path / "out.bw"),
        bedgraphtobigwig_path=str(tool),
        verbose=False,
    )
    assert result["status"] == "failed"
    assert result["return_code"] == 3
    assert result["outputs"] == []
    assert "bad input" in result["log_tail"]


def test_mapped_read_count_parses_samtools_output(tmp_path):
    pytest.importorskip("pybedtools")
    bam = tmp_path / "reads.bam"
    bam.write_bytes(b"x")
    tool_dir = tmp_path / "samtools_dir"
    tool = _write_stand_in(
        tool_dir, "samtools", 'printf "%s\\n" "$@" > "$(dirname "$0")/argv.txt"\nprintf "42\\n"\n'
    )
    result = module.mapped_read_count(str(bam), samtools_path=str(tool), verbose=False)
    assert result["status"] == "completed"
    assert result["mapped_reads"] == 42
    assert _argv(tool_dir) == ["view", "-c", "-F", "0x4", str(bam.resolve())]
    assert result["outputs"] == []


def test_mapped_read_count_reports_samtools_failure(tmp_path):
    pytest.importorskip("pybedtools")
    bam = tmp_path / "reads.bam"
    bam.write_bytes(b"x")
    tool = _write_stand_in(tmp_path / "bad", "samtools", 'echo "not a BAM" >&2\nexit 1\n')
    result = module.mapped_read_count(str(bam), samtools_path=str(tool), verbose=False)
    assert result["status"] == "failed"
    assert result["mapped_reads"] is None
    assert "not a BAM" in result["log_tail"]


@_REQUIRES_BEDTOOLS
def test_bam_to_bigwig_uses_real_bedtools_and_removes_intermediates(tmp_path):
    pytest.importorskip("pybedtools")
    bam = tmp_path / "reads.bam"
    _write_bam(bam)
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n")
    seen = tmp_path / "seen.bedgraph"
    tool_dir = tmp_path / "converter"
    tool = _write_stand_in(tool_dir, "bedGraphToBigWig", f'cp "$1" "{seen}"\n' + _STAND_IN_HEAD)
    out = tmp_path / "coverage.bw"
    result = module.bam_to_bigwig(
        str(bam), str(sizes), str(out), bedgraphtobigwig_path=str(tool), verbose=False
    )
    assert result["status"] == "completed"
    assert result["outputs"] == [str(out)]
    rows = [line.split("\t") for line in seen.read_text().splitlines() if line.strip()]
    assert [(r[0], int(r[1]), int(r[2]), float(r[3])) for r in rows][:2] == [
        ("chr1", 10, 20, 1.0),
        ("chr1", 20, 30, 2.0),
    ]
    assert not list(tmp_path.glob("coverage.bw.tmp-*"))


@_REQUIRES_BEDTOOLS
def test_bam_to_bigwig_scale_uses_samtools_count(tmp_path):
    pytest.importorskip("pybedtools")
    bam = tmp_path / "reads.bam"
    _write_bam(bam)
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t1000\n")
    seen = tmp_path / "seen_scaled.bedgraph"
    converter = _write_stand_in(
        tmp_path / "converter", "bedGraphToBigWig", f'cp "$1" "{seen}"\n' + _STAND_IN_HEAD
    )
    samtools = _write_stand_in(tmp_path / "samtools_dir", "samtools", 'printf "4\\n"\n')
    result = module.bam_to_bigwig(
        str(bam),
        str(sizes),
        str(tmp_path / "scaled.bw"),
        scale=True,
        samtools_path=str(samtools),
        bedgraphtobigwig_path=str(converter),
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["mapped_reads"] == 4
    assert result["scale_factor"] == pytest.approx(250000.0)
    first = seen.read_text().splitlines()[0].split("\t")
    assert float(first[3]) == pytest.approx(250000.0)


@_REQUIRES_SAMTOOLS
def test_mapped_read_count_real_samtools(tmp_path):
    pytest.importorskip("pysam")
    bam = tmp_path / "reads.bam"
    _write_bam(bam)
    result = module.mapped_read_count(str(bam), verbose=False)
    assert result["mapped_reads"] == 4


@_REQUIRES_BEDTOBIGBED
def test_bigbed_round_trip_with_real_ucsc_tools(tmp_path):
    pytest.importorskip("pybedtools")
    bed = tmp_path / "in.bed"
    bed.write_text("chr1\t0\t10\tfeat\t0\t+\nchr1\t20\t30\tfeat2\t0\t-\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    out = tmp_path / "out.bb"
    assert module.bigbed(str(bed), str(sizes), str(out), verbose=False)["status"] == "completed"
    assert out.read_bytes()[:4] == bytes([0xEB, 0xF2, 0x89, 0x87])
    if shutil.which("bigBedToBed") is not None:
        back = module.bigbed_to_bed(str(out), str(tmp_path / "back.bed"), verbose=False)
        assert back["status"] == "completed"


@_REQUIRES_BEDGRAPHTOBIGWIG
def test_bedgraph_to_bigwig_with_real_ucsc_tools(tmp_path):
    graph = tmp_path / "in.bedGraph"
    graph.write_text("chr1\t0\t10\t1.5\nchr1\t10\t20\t2.5\n")
    sizes = tmp_path / "sizes.txt"
    sizes.write_text("chr1\t100\n")
    out = tmp_path / "out.bw"
    result = module.bedgraph_to_bigwig(str(graph), str(sizes), str(out), verbose=False)
    assert result["status"] == "completed"
    assert out.read_bytes()[:4] == bytes([0x26, 0xFC, 0x8F, 0x88])


# ---- matrix, Venn and colour helpers
# ----------------------------------------------------------------------


def test_sort_binary_matrix_matches_pybedtools():
    pytest.importorskip("numpy")
    plotting = pytest.importorskip("pybedtools.contrib.plotting")
    import numpy as np

    matrix = [[1, 0, 0], [0, 1, 1], [1, 1, 0], [0, 0, 1], [1, 0, 1]]
    ours = module.sort_binary_matrix(matrix, verbose=False)
    theirs = plotting.sort_binary_matrix(np.array(matrix)).tolist()
    assert ours == theirs
    assert module.sort_binary_matrix([[1, 1], [1, 0], [0, 1]], verbose=False) == [0, 1, 2]
    assert module.sort_binary_matrix([], verbose=False) == []
    with pytest.raises(ConfigurationError):
        module.sort_binary_matrix([[1, 2]], verbose=False)
    with pytest.raises(ConfigurationError):
        module.sort_binary_matrix([[1, 0], [1]], verbose=False)


@_REQUIRES_BEDTOOLS
def test_cleaned_intersect_matches_documented_example(tmp_path):
    pytest.importorskip("pybedtools")
    venn = pytest.importorskip("pybedtools.contrib.venn_maker")
    import pybedtools

    x, y = _venn_example(tmp_path)
    result = module.cleaned_intersect([str(x), str(y)], verbose=False)
    assert result["sets"] == [["chr1|1|100", "chr1|500|6000"], ["chr1|1|100", "chr9|777|888"]]
    assert result["set_sizes"] == [2, 2]
    upstream = venn.cleaned_intersect([pybedtools.BedTool(str(x)), pybedtools.BedTool(str(y))])
    assert [len(part) for part in upstream] == result["set_sizes"]
    with pytest.raises(ConfigurationError):
        module.cleaned_intersect([str(x)], verbose=False)
    with pytest.raises(ConfigurationError):
        module.cleaned_intersect([str(x)] * 5, verbose=False)


@_REQUIRES_BEDTOOLS
def test_venn_maker_writes_script_and_region_counts(tmp_path, monkeypatch):
    pytest.importorskip("pybedtools")
    x, y = _venn_example(tmp_path)
    z = tmp_path / "z.bed"
    z.write_text("chr1\t90\t120\th1\t0\t+\nchr9\t800\t850\th2\t0\t+\n")
    script = tmp_path / "venn.R"
    real_which = shutil.which
    monkeypatch.setattr(
        module.shutil,
        "which",
        lambda name, *a, **k: None if name == "Rscript" else real_which(name, *a, **k),
    )
    result = module.venn_maker(
        [str(x), str(y), str(z)],
        str(script),
        names=["x", "y", "z"],
        figure_path=str(tmp_path / "venn.tiff"),
        verbose=False,
    )
    text = script.read_text()
    assert "library(VennDiagram)" in text and "venn.diagram(" in text
    assert 'category.names = c("x","y","z")' in text
    assert result["rendered"] is False
    assert "Rscript was not found" in result["render_note"]
    assert result["data"]["regions"] == {"x": 1, "x&y&z": 1, "y&z": 1}
    assert result["data"]["set_sizes"] == {"x": 2, "y": 2, "z": 2}


@_REQUIRES_BEDTOOLS
def test_venn_maker_renders_with_stand_in_rscript(tmp_path, monkeypatch):
    pytest.importorskip("pybedtools")
    x, y = _venn_example(tmp_path)
    script = tmp_path / "venn2.R"
    figure = tmp_path / "venn2.tiff"
    tool = _write_stand_in(
        tmp_path / "r_tool",
        "Rscript",
        'fig=$(sed -n \'s/.*filename="\\([^"]*\\)".*/\\1/p\' "$1")\nprintf "tiff" > "$fig"\n',
    )
    real_which = shutil.which
    monkeypatch.setattr(
        module.shutil,
        "which",
        lambda name, *a, **k: str(tool) if name == "Rscript" else real_which(name, *a, **k),
    )
    result = module.venn_maker(
        [str(x), str(y)], str(script), figure_path=str(figure), verbose=False
    )
    assert result["rendered"] is True
    assert figure.read_bytes() == b"tiff"
    assert result["data"]["names"] == ["a", "b"]


def test_venn_maker_rejects_bad_names_and_arguments(tmp_path):
    pytest.importorskip("pybedtools")
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t1\t10\n")
    other = tmp_path / "b.bed"
    other.write_text("chr1\t5\t15\n")
    with pytest.raises(ConfigurationError):
        module.venn_maker(
            [str(bed), str(other)], str(tmp_path / "s.R"), names=["a", 'b"c'], verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.venn_maker(
            [str(bed), str(other)], str(tmp_path / "s.R"), names=["a", "a"], verbose=False
        )
    with pytest.raises(ConfigurationError):
        module.venn_maker(
            [str(bed), str(other)],
            str(tmp_path / "s.R"),
            additional_args=["x=1); system('rm -rf /')"],
            verbose=False,
        )


@_REQUIRES_BEDTOOLS
def test_binary_heatmap_counts_classes_and_writes_png(tmp_path):
    pytest.importorskip("pybedtools")
    pytest.importorskip("matplotlib")
    x, y = _venn_example(tmp_path)
    out = tmp_path / "heatmap.png"
    result = module.binary_heatmap([str(x), str(y)], ["a", "b"], str(out), verbose=False)
    assert result["status"] == "completed"
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    counts = result["data"]["class_counts"]
    assert sum(counts.values()) == result["data"]["interval_count"]
    assert set(counts) <= {"a", "b", "a,b"}
    assert len(result["data"]["matrix"]) == result["data"]["interval_count"]
    assert sorted(result["data"]["row_order"]) == list(range(result["data"]["interval_count"]))


def test_add_color_matches_upstream_and_validates(tmp_path):
    pybedtools = pytest.importorskip("pybedtools")
    pytest.importorskip("matplotlib")
    featurefuncs = pytest.importorskip("pybedtools.featurefuncs")
    import matplotlib
    from matplotlib.colors import Normalize

    source = tmp_path / "scores.bed"
    source.write_text("chr1\t0\t10\ta\t0\t+\nchr1\t20\t30\tb\t500\t+\nchr1\t40\t50\tc\t1000\t-\n")
    out = tmp_path / "coloured.bed"
    result = module.add_color(str(source), str(out), colormap="viridis", verbose=False)
    assert result["status"] == "completed"
    assert result["data"]["feature_count"] == 3
    assert result["data"]["score_min"] == 0.0 and result["data"]["score_max"] == 1000.0
    upstream = pybedtools.BedTool(str(source)).each(
        featurefuncs.add_color,
        cmap=matplotlib.colormaps["viridis"],
        norm=Normalize(vmin=0.0, vmax=1000.0),
    )
    expected = [line.split("\t")[8] for line in str(upstream).splitlines()]
    assert result["data"]["rgb"] == expected
    assert result["data"]["rgb"][0] != result["data"]["rgb"][2]
    with pytest.raises(ConfigurationError):
        module.add_color(str(source), str(tmp_path / "x.bed"), colormap="not_a_map", verbose=False)
    with pytest.raises(ConfigurationError):
        module.add_color(str(source), str(out), verbose=False)
