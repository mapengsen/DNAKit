"""Contract tests for the BEDOPS command-line facades (no network; BEDOPS tests skip without bedops)."""

from __future__ import annotations

import inspect
import shutil
import subprocess
from pathlib import Path

import pytest

import dnakit.annotation.bedops_p05_cli_functions as bo
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

requires_bedops = pytest.mark.skipif(
    any(
        shutil.which(name) is None
        for name in ("bedops", "bedmap", "sort-bed", "starch", "convert2bed")
    ),
    reason="the BEDOPS executables are not on PATH",
)

# Sorted BED (sort-bed order): reference intervals, map intervals with scores, and a second chromosome.
REFERENCE = "chr1\t100\t200\tr1\nchr1\t300\t400\tr2\nchr1\t1000\t1100\tr3\nchr2\t10\t60\tr4\n"
MAP = (
    "chr1\t150\t250\tm1\t2\t+\n"
    "chr1\t180\t350\tm2\t4\t-\n"
    "chr1\t320\t330\tm3\t6\t+\n"
    "chr2\t0\t20\tm4\t10\t+\n"
)
VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t101\trs1\tA\tG\t50\tPASS\t.\n"
    "chr1\t51\trs2\tC\tT\t60\tPASS\t.\n"
)


def _write(path: Path, text: str) -> str:
    path.write_text(text, encoding="utf-8")
    return str(path)


def _rows(path: str | Path) -> list[list[str]]:
    return [line.split("\t") for line in Path(path).read_text().splitlines()]


# ---------------------------------------------------------------------------
# Contract without BEDOPS
# ---------------------------------------------------------------------------


def test_public_functions_are_keyword_only_verbose() -> None:
    assert len(bo.__all__) == 13
    for name in bo.__all__:
        parameter = inspect.signature(getattr(bo, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True
        assert getattr(bo, name).__doc__


def test_operation_vocabularies_cover_the_work_list() -> None:
    assert set(bo.SET_OPERATIONS) == {"element-of", "not-element-of", "symmdiff", "everything"}
    for statistic in (
        "cv",
        "kth",
        "mad",
        "max",
        "mean",
        "median",
        "min",
        "stdev",
        "sum",
        "tmean",
        "variance",
        "wmean",
        "bases",
        "count",
    ):
        assert statistic in bo.SCORE_STATISTICS
    for field in (
        "echo",
        "echo-map",
        "echo-map-id",
        "echo-map-id-uniq",
        "echo-map-range",
        "echo-map-score",
        "echo-map-size",
        "echo-overlap-size",
        "echo-ref-name",
        "echo-ref-row-id",
        "max-element",
        "min-element",
    ):
        assert field in bo.ECHO_FIELDS
    assert set(bo.CONVERT_INPUT_FORMATS) == {
        "bam",
        "gff",
        "gtf",
        "gvf",
        "psl",
        "rmsk",
        "sam",
        "vcf",
        "wig",
    }


def test_missing_program_is_reported_with_install_hint(tmp_path: Path) -> None:
    source = _write(tmp_path / "a.bed", REFERENCE)
    with pytest.raises(BackendUnavailableError) as caught:
        bo.bedextract_list_chromosomes(source, bedops_bin_dir=str(tmp_path), verbose=False)
    assert "bedops" in str(caught.value).lower() or "bedextract" in str(caught.value).lower()
    with pytest.raises(BackendUnavailableError):
        bo.bedextract_list_chromosomes(
            source, bedops_bin_dir=str(tmp_path / "missing"), verbose=False
        )


def test_arguments_are_validated_before_any_program_runs(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", REFERENCE)
    b = _write(tmp_path / "b.bed", MAP)
    existing = _write(tmp_path / "exists.txt", "x")
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation([a, b], str(tmp_path / "o1"), operation="intersect", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation([a], str(tmp_path / "o2"), operation="symmdiff", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation(
            [a, b], str(tmp_path / "o3"), operation="symmdiff", overlap=1, verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation([a, b], existing, verbose=False)  # outputs are never overwritten
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation(
            [a, str(tmp_path / "nope.bed")], str(tmp_path / "o4"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation([a, b], str(tmp_path / "o5"), overlap="half", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedops_set_operation([a, b], str(tmp_path / "o6"), chromosome="--ec", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedmap_score_statistic(
            a, str(tmp_path / "o7"), map_path=b, statistic="kth", verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.bedmap_score_statistic(
            a, str(tmp_path / "o8"), map_path=b, statistic="mean", kth_fraction=0.5, verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.bedmap_score_statistic(
            a,
            str(tmp_path / "o9"),
            map_path=b,
            statistic="tmean",
            trim_low=0.7,
            trim_high=0.7,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        bo.bedmap_score_statistic(
            a,
            str(tmp_path / "o10"),
            map_path=b,
            overlap_mode="fraction-ref",
            overlap_value=2,
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        bo.bedmap_echo_elements(
            a, str(tmp_path / "o11"), map_path=b, field="echo-everything", verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.convert_to_bed(a, str(tmp_path / "o12"), input_format="bed", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.convert_to_bed(
            a,
            str(tmp_path / "o13"),
            input_format="vcf",
            format_options=["zero-indexed"],
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        bo.convert_to_bed(
            a,
            str(tmp_path / "o14"),
            input_format="vcf",
            sort_output=False,
            output_format="starch",
            verbose=False,
        )
    with pytest.raises(ConfigurationError):
        bo.convert_to_bed(
            a, str(tmp_path / "o15"), input_format="vcf", attribute_key="gene_id", verbose=False
        )
    with pytest.raises(ConfigurationError):
        bo.starchstrip_chromosomes(a, str(tmp_path / "o16"), [], verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedops_sort_bed_filter([a], str(tmp_path / "o17"), mode="all", verbose=False)
    with pytest.raises(ConfigurationError):
        bo.bedops_closest_features(a, b, str(tmp_path / "o18"), delimiter="-x", verbose=False)
    assert not any(path.name.startswith("o") for path in tmp_path.iterdir())


# ---------------------------------------------------------------------------
# With BEDOPS
# ---------------------------------------------------------------------------


@requires_bedops
def test_set_operations(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", REFERENCE)
    b = _write(tmp_path / "b.bed", MAP)
    element_of = bo.bedops_set_operation(
        [a, b], str(tmp_path / "e.bed"), operation="element-of", overlap=1, verbose=False
    )
    assert element_of["status"] == "completed" and element_of["backend"] == "bedops"
    assert element_of["version"] and element_of["provenance"]["tested_version"] == "2.4.42"
    assert [row[3] for row in _rows(tmp_path / "e.bed")] == ["r1", "r2", "r4"]
    assert element_of["line_count"] == 3
    not_element = bo.bedops_set_operation(
        [a, b], str(tmp_path / "n.bed"), operation="not-element-of", overlap=1, verbose=False
    )
    assert [row[3] for row in _rows(tmp_path / "n.bed")] == ["r3"] and not_element[
        "line_count"
    ] == 1
    # By default BEDOPS requires 100% of the reference element to be covered by one element.
    bo.bedops_set_operation(
        [a, b], str(tmp_path / "e100.bed"), operation="element-of", verbose=False
    )
    assert _rows(tmp_path / "e100.bed") == []
    bo.bedops_set_operation(
        [a, b], str(tmp_path / "e50.bed"), operation="element-of", overlap="50%", verbose=False
    )
    assert [row[3] for row in _rows(tmp_path / "e50.bed")] == ["r1", "r2"]
    bo.bedops_set_operation([a, b], str(tmp_path / "s.bed"), operation="symmdiff", verbose=False)
    assert _rows(tmp_path / "s.bed")[0][:3] == ["chr1", "100", "150"]
    bo.bedops_set_operation([a, b], str(tmp_path / "u.bed"), operation="everything", verbose=False)
    assert len(_rows(tmp_path / "u.bed")) == 8


@requires_bedops
def test_chromosome_listing_and_extraction_and_sorting(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", REFERENCE)
    assert bo.bedextract_list_chromosomes(a, verbose=False)["chromosomes"] == ["chr1", "chr2"]
    result = bo.bedextract_chromosome(a, "chr2", str(tmp_path / "chr2.bed"), verbose=False)
    assert (
        _rows(tmp_path / "chr2.bed") == [["chr2", "10", "60", "r4"]] and result["line_count"] == 1
    )
    bo.bedextract_chromosome(a, "chrZ", str(tmp_path / "none.bed"), verbose=False)
    assert (tmp_path / "none.bed").read_text() == ""
    messy = _write(
        tmp_path / "messy.bed", "chr2\t5\t9\tx\nchr1\t1\t2\ty\nchr1\t1\t2\ty\nchr1\t7\t8\tz\n"
    )
    bo.bedops_sort_bed_filter([messy], str(tmp_path / "uniq.bed"), mode="unique", verbose=False)
    assert [row[3] for row in _rows(tmp_path / "uniq.bed")] == ["y", "z", "x"]
    bo.bedops_sort_bed_filter([messy], str(tmp_path / "dup.bed"), mode="duplicates", verbose=False)
    assert [row[3] for row in _rows(tmp_path / "dup.bed")] == ["y"]


@requires_bedops
def test_closest_features(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", "chr1\t500\t600\tq1\n")
    b = _write(tmp_path / "b.bed", "chr1\t100\t200\tl\nchr1\t800\t900\tr\n")
    result = bo.bedops_closest_features(
        a, b, str(tmp_path / "c.txt"), distances=True, verbose=False
    )
    line = (tmp_path / "c.txt").read_text().strip()
    assert line == "chr1\t500\t600\tq1|chr1\t100\t200\tl|-301|chr1\t800\t900\tr|201"
    assert result["line_count"] == 1
    bo.bedops_closest_features(
        a,
        b,
        str(tmp_path / "c2.txt"),
        closest=True,
        no_reference=True,
        delimiter=";",
        verbose=False,
    )
    assert (tmp_path / "c2.txt").read_text().strip() == "chr1\t800\t900\tr"


@requires_bedops
def test_bedmap_statistics_and_echo(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", REFERENCE)
    b = _write(tmp_path / "b.bed", MAP)
    values = {}
    for statistic in (
        "mean",
        "sum",
        "min",
        "max",
        "median",
        "count",
        "bases",
        "variance",
        "stdev",
        "cv",
        "wmean",
    ):
        result = bo.bedmap_score_statistic(
            a, str(tmp_path / f"{statistic}.txt"), map_path=b, statistic=statistic, verbose=False
        )
        values[statistic] = result["values"]
    assert values["count"] == [2, 2, 0, 1]
    assert values["sum"][:2] == [6, 10] and values["sum"][2] == "NAN"
    assert values["mean"][:2] == [3.0, 5.0]
    assert values["min"][:2] == [2, 4] and values["max"][:2] == [4, 6]
    assert values["bases"][:2] == [70, 60]
    kth = bo.bedmap_score_statistic(
        a, str(tmp_path / "kth.txt"), map_path=b, statistic="kth", kth_fraction=0.5, verbose=False
    )
    assert kth["command"][1:3] == ["--kth", "0.5"]
    tmean = bo.bedmap_score_statistic(
        a,
        str(tmp_path / "tm.txt"),
        map_path=b,
        statistic="tmean",
        trim_low=0.0,
        trim_high=0.0,
        verbose=False,
    )
    assert tmean["values"][:2] == [3.0, 5.0]
    mad = bo.bedmap_score_statistic(
        a, str(tmp_path / "mad.txt"), map_path=b, statistic="mad", mad_multiplier=2, verbose=False
    )
    assert mad["status"] == "completed"
    frac = bo.bedmap_score_statistic(
        a,
        str(tmp_path / "frac.txt"),
        map_path=b,
        statistic="count",
        overlap_mode="fraction-ref",
        overlap_value=0.5,
        verbose=False,
    )
    assert frac["values"] == [1, 1, 0, 0]
    echo_ids = bo.bedmap_echo_elements(
        a, str(tmp_path / "ids.txt"), map_path=b, field="echo-map-id", verbose=False
    )
    assert echo_ids["preview"] == ["m1;m2", "m2;m3", "", "m4"]
    scores = bo.bedmap_echo_elements(
        a,
        str(tmp_path / "scores.txt"),
        map_path=b,
        field="echo-map-score",
        delimiter="\t",
        verbose=False,
    )
    assert scores["preview"][0] == "2.000000;4.000000"
    names = bo.bedmap_echo_elements(
        a, str(tmp_path / "names.txt"), map_path=b, field="echo-ref-name", verbose=False
    )
    assert names["preview"][0] == "chr1:100-200"
    best = bo.bedmap_echo_elements(
        a,
        str(tmp_path / "maxel.txt"),
        map_path=b,
        field="max-element",
        skip_unmapped=True,
        verbose=False,
    )
    assert best["line_count"] == 3 and best["preview"][0].startswith("chr1\t180\t350\tm2")
    overlap_size = bo.bedmap_echo_elements(
        a, str(tmp_path / "ov.txt"), map_path=b, field="echo-overlap-size", verbose=False
    )
    assert overlap_size["preview"][0] == "50;20"


@requires_bedops
def test_starch_archive_family(tmp_path: Path) -> None:
    a = _write(tmp_path / "a.bed", REFERENCE)
    b = _write(tmp_path / "b.bed", MAP)
    first = bo.starch_compress(
        a, str(tmp_path / "a.starch"), compression="gzip", note="unit test", verbose=False
    )
    assert first["element_count"] == 4 and (tmp_path / "a.starch").is_file()
    bo.starch_compress(b, str(tmp_path / "b.starch"), verbose=False)
    back = bo.unstarch_to_bed(str(tmp_path / "a.starch"), str(tmp_path / "back.bed"), verbose=False)
    assert (tmp_path / "back.bed").read_text() == REFERENCE and back["line_count"] == 4
    bo.unstarch_to_bed(
        str(tmp_path / "a.starch"), str(tmp_path / "chr2.bed"), chromosome="chr2", verbose=False
    )
    assert (tmp_path / "chr2.bed").read_text() == "chr2\t10\t60\tr4\n"

    def summary(name: str, **kw: str) -> object:
        return bo.unstarch_summary(
            str(tmp_path / "a.starch"), statistic=name, verbose=False, **kw
        )["value"]

    assert summary("elements") == 4 and summary("elements", chromosome="chr1") == 3
    assert summary("bases") == 100 + 100 + 100 + 50 and summary("bases-uniq") == 350
    assert summary("has-duplicate") is False and summary("has-nested") is False
    assert summary("list-chromosomes") == ["chr1", "chr2"]
    assert summary("list-json")["archive"]["type"] == "starch"
    nested = _write(tmp_path / "nested.bed", "chr1\t0\t100\ta\nchr1\t10\t20\tb\nchr1\t10\t20\tb\n")
    bo.starch_compress(nested, str(tmp_path / "n.starch"), verbose=False)
    assert (
        bo.unstarch_summary(str(tmp_path / "n.starch"), statistic="has-nested", verbose=False)[
            "value"
        ]
        is True
    )
    assert (
        bo.unstarch_summary(str(tmp_path / "n.starch"), statistic="has-duplicate", verbose=False)[
            "value"
        ]
        is True
    )
    merged = bo.starchcat_archives(
        [str(tmp_path / "a.starch"), str(tmp_path / "b.starch")],
        str(tmp_path / "ab.starch"),
        verbose=False,
    )
    assert merged["element_count"] == 8
    stripped = bo.starchstrip_chromosomes(
        str(tmp_path / "ab.starch"), str(tmp_path / "only2.starch"), ["chr2"], verbose=False
    )
    assert stripped["chromosomes_kept"] == ["chr2"] and stripped["element_count"] == 2
    dropped = bo.starchstrip_chromosomes(
        str(tmp_path / "ab.starch"),
        str(tmp_path / "no1.starch"),
        ["chr1"],
        mode="exclude",
        verbose=False,
    )
    assert dropped["chromosomes_kept"] == ["chr2"]
    bo.unstarch_to_bed(str(tmp_path / "no1.starch"), str(tmp_path / "no1.bed"), verbose=False)
    assert len(_rows(tmp_path / "no1.bed")) == 2  # readable even though chr1 was the first stream
    with pytest.raises(ConfigurationError):
        bo.starchstrip_chromosomes(
            str(tmp_path / "ab.starch"),
            str(tmp_path / "none.starch"),
            ["chr1", "chr2"],
            mode="exclude",
            verbose=False,
        )


@requires_bedops
def test_convert_to_bed_and_starch(tmp_path: Path) -> None:
    vcf = _write(tmp_path / "v.vcf", VCF)
    bed = bo.convert_to_bed(vcf, str(tmp_path / "v.bed"), input_format="vcf", verbose=False)
    assert [row[:4] for row in _rows(tmp_path / "v.bed")] == [
        ["chr1", "50", "51", "rs2"],
        ["chr1", "100", "101", "rs1"],
    ]
    assert bed["input_format"] == "vcf" and bed["line_count"] == 2
    unsorted = bo.convert_to_bed(
        vcf, str(tmp_path / "u.bed"), input_format="vcf", sort_output=False, verbose=False
    )
    assert _rows(tmp_path / "u.bed")[0][3] == "rs1" and unsorted["line_count"] == 2
    starch = bo.convert_to_bed(
        vcf,
        str(tmp_path / "v.starch"),
        input_format="vcf",
        output_format="starch",
        starch_note="vcf",
        verbose=False,
    )
    assert starch["element_count"] == 2
    gtf = _write(
        tmp_path / "g.gtf",
        'chr1\tdemo\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
    )
    bo.convert_to_bed(
        gtf,
        str(tmp_path / "g.bed"),
        input_format="gtf",
        attribute_key="transcript_id",
        verbose=False,
    )
    assert _rows(tmp_path / "g.bed")[0][:4] == ["chr1", "100", "200", "t1"]
    wig = _write(tmp_path / "s.wig", "fixedStep chrom=chr1 start=101 step=10 span=10\n1.5\n2.5\n")
    bo.convert_to_bed(wig, str(tmp_path / "s.bed"), input_format="wig", verbose=False)
    assert [row[:3] for row in _rows(tmp_path / "s.bed")] == [
        ["chr1", "100", "110"],
        ["chr1", "110", "120"],
    ]
    rmsk = _write(
        tmp_path / "r.out",
        "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat\n"
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID\n\n"
        "  463  10.2  0.0  1.2  chr1        101     200  (900)  +  L1PA2         LINE/L1                 1   100  (10)    1\n",
    )
    bo.convert_to_bed(rmsk, str(tmp_path / "r.bed"), input_format="rmsk", verbose=False)
    assert _rows(tmp_path / "r.bed")[0][:4] == ["chr1", "100", "200", "L1PA2"]


@requires_bedops
@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools is not on PATH")
def test_convert_sam_to_bed(tmp_path: Path) -> None:
    sam = _write(
        tmp_path / "r.sam",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
        "r1\t0\tchr1\t11\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n",
    )
    result = bo.convert_to_bed(sam, str(tmp_path / "r.bed"), input_format="sam", verbose=False)
    assert (
        _rows(tmp_path / "r.bed")[0][:4] == ["chr1", "10", "20", "r1"] and result["line_count"] == 1
    )


def test_sam_conversion_needs_samtools(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sam = _write(tmp_path / "r.sam", "@HD\tVN:1.6\n")
    monkeypatch.setattr(bo.shutil, "which", lambda name: None)
    with pytest.raises(BackendUnavailableError):
        bo.convert_to_bed(sam, str(tmp_path / "r.bed"), input_format="sam", verbose=False)


@requires_bedops
def test_failed_run_raises_and_removes_partial_output(tmp_path: Path) -> None:
    from dnakit.exceptions import BackendExecutionError

    bad = _write(tmp_path / "bad.bed", "chr1\tnot-a-number\t10\n")
    with pytest.raises(BackendExecutionError):
        bo.bedops_set_operation(
            [bad, bad],
            str(tmp_path / "out.bed"),
            operation="element-of",
            error_check=True,
            verbose=False,
        )
    assert not (tmp_path / "out.bed").exists()
    assert subprocess.run(["bedops", "--version"], capture_output=True).returncode == 0
