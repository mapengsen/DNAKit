"""Tests for dnakit.variants.picard_vcf_intervals_p08_cli_functions (Picard VCF/interval tools).

The tests run the real Picard (3.5.0 tested) on a small synthetic dataset. They skip when Picard
(`picard` on PATH, or `PICARD_JAR` plus `java`) is not installed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

import _picard_p08_data as data

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import picard_vcf_intervals_p08_cli_functions as m

requires_picard = pytest.mark.skipif(
    not data.picard_available(),
    reason="Picard (picard on PATH, or PICARD_JAR with java) is not installed",
)


@pytest.fixture(scope="module")
def files(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return data.write_dataset(tmp_path_factory.mktemp("picard_vcf"))


def _ok(result: dict) -> dict:
    assert result["status"] == "completed", result["log_tail"][-1500:]
    assert result["backend"] == "picard"
    return result


def _data_lines(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line and not line.startswith(("#", "@"))
    ]


def _interval_list(path: Path, rows: list[tuple[str, int, int, str]]) -> Path:
    header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1200\n@SQ\tSN:chr2\tLN:800\n"
    path.write_text(header + "".join(f"{c}\t{s}\t{e}\t+\t{n}\n" for c, s, e, n in rows))
    return path


def test_missing_picard_and_overwrite_protection(files: dict[str, Path], tmp_path: Path) -> None:
    with pytest.raises(BackendUnavailableError):
        m.picard_filter_vcf(
            str(files["vcf"]),
            str(tmp_path / "o.vcf"),
            picard_path=str(tmp_path / "x"),
            verbose=False,
        )
    target = tmp_path / "o.vcf"
    target.write_text("keep")
    with pytest.raises(ConfigurationError):
        m.picard_filter_vcf(str(files["vcf"]), str(target), verbose=False)
    assert target.read_text() == "keep"


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (
            lambda f, t: m.picard_interval_list_tools(
                [str(f["sam"])], str(t / "o"), action="MERGE", verbose=False
            ),
            "action",
        ),
        (
            lambda f, t: m.picard_interval_list_tools(
                [str(f["sam"])], str(t / "o"), action="SUBTRACT", verbose=False
            ),
            "second_input_paths",
        ),
        (
            lambda f, t: m.picard_interval_list_to_bed(
                str(f["sam"]), str(t / "o"), score=2000, verbose=False
            ),
            "score",
        ),
        (
            lambda f, t: m.picard_scatter_intervals_by_ns(
                str(f["fasta"]), str(t / "o"), output_type="X", verbose=False
            ),
            "output_type",
        ),
        (
            lambda f, t: m.picard_normalize_fasta(
                str(f["fasta"]), str(t / "o"), line_length=0, verbose=False
            ),
            "line_length",
        ),
    ],
)
def test_arguments_are_validated(
    files: dict[str, Path], tmp_path: Path, call, message: str
) -> None:
    with pytest.raises(ConfigurationError) as error:
        call(files, tmp_path)
    assert message in str(error.value)


@requires_picard
def test_filter_vcf_removes_low_quality_genotypes(files: dict[str, Path], tmp_path: Path) -> None:
    plain = tmp_path / "plain.vcf"
    plain.write_text((files["vcf"].parent / "sample.vcf").read_text())
    result = _ok(m.picard_filter_vcf(str(plain), str(tmp_path / "f.vcf"), min_gq=34, verbose=False))
    assert result["records"] == 8 and result["index"] is not None
    rows = _data_lines(tmp_path / "f.vcf")
    flagged = ["LowGQ" in row[9] for row in rows]
    assert flagged == [True] * 4 + [False] * 4  # GQ is 30 + index; 34 is the cut-off
    site = _ok(m.picard_filter_vcf(str(plain), str(tmp_path / "g.vcf"), min_qd=8.5, verbose=False))
    filters = [row[6] for row in _data_lines(Path(site["outputs"]["output"]))]
    assert filters == ["LowQD"] * 4 + ["PASS"] * 4  # QD is 5 + index


@requires_picard
def test_make_sites_only_vcf_drops_the_genotype_columns(
    files: dict[str, Path], tmp_path: Path
) -> None:
    result = _ok(
        m.picard_make_sites_only_vcf(str(files["vcf"]), str(tmp_path / "sites.vcf"), verbose=False)
    )
    rows = _data_lines(Path(result["outputs"]["output"]))
    assert result["records"] == 8 and all(len(row) == 8 for row in rows)


@requires_picard
def test_vcf_sample_name_map_and_rename(files: dict[str, Path], tmp_path: Path) -> None:
    other = tmp_path / "other.vcf"
    other.write_text((files["vcf"].parent / "sample.vcf").read_text().replace("sample1", "sample2"))
    zipped = data.index_vcf(other)
    result = _ok(
        m.picard_make_vcf_sample_name_map(
            [str(files["vcf"]), str(zipped)], str(tmp_path / "map.tsv"), verbose=False
        )
    )
    assert {entry["sample"] for entry in result["entries"]} == {"sample1", "sample2"}
    assert {Path(entry["path"]).name for entry in result["entries"]} == {
        "sample.vcf.gz",
        "other.vcf.gz",
    }
    renamed = _ok(
        m.picard_rename_sample_in_vcf(
            str(files["vcf"]),
            str(tmp_path / "renamed.vcf"),
            "NEW",
            old_sample_name="sample1",
            verbose=False,
        )
    )
    header = next(
        line
        for line in Path(renamed["outputs"]["output"]).read_text().splitlines()
        if line.startswith("#CHROM")
    )
    assert header.endswith("\tNEW")
    wrong = m.picard_rename_sample_in_vcf(
        str(files["vcf"]), str(tmp_path / "bad.vcf"), "NEW", old_sample_name="other", verbose=False
    )
    assert wrong["status"] == "failed"


@requires_picard
def test_split_vcfs_separates_snps_from_indels(files: dict[str, Path], tmp_path: Path) -> None:
    mixed = tmp_path / "mixed.vcf"
    lines = (files["vcf"].parent / "sample.vcf").read_text().splitlines()
    head = [line for line in lines if line.startswith("#")]
    body = [line for line in lines if not line.startswith("#")]
    ref, pos = body[0].split("\t")[3], body[0].split("\t")[1]
    info = "DP=10;QD=5;FS=1.0\tGT:DP:GQ\t0/1:8:30"
    indel = f"chr1\t{int(pos) + 1}\t.\t{ref}\t{ref}A\t50\t.\t{info}"
    mixed.write_text("\n".join([*head, *body, indel]) + "\n")
    result = _ok(
        m.picard_split_vcfs(
            str(mixed), str(tmp_path / "snp.vcf"), str(tmp_path / "indel.vcf"), verbose=False
        )
    )
    assert result["snp_records"] == 8 and result["indel_records"] == 1
    other = _ok(
        m.picard_split_vcfs(
            str(files["vcf"]),
            str(tmp_path / "s2.vcf"),
            str(tmp_path / "i2.vcf"),
            strict=False,
            verbose=False,
        )
    )
    assert other["indel_records"] == 0


def _lift_files(files: dict[str, Path], folder: Path) -> dict[str, Path]:
    chain = folder / "a_to_b.chain"
    chain.write_text(
        "chain 1000 chr1 1200 + 0 1200 chrA 1200 + 0 1200 1\n1200\n\n"
        "chain 1000 chr2 800 + 0 800 chrB 800 + 0 800 2\n800\n\n"
    )
    target = folder / "target.fa"
    target.write_text(
        files["fasta"].read_text().replace(">chr1", ">chrA").replace(">chr2", ">chrB")
    )
    (folder / "target.dict").write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:1200\n@SQ\tSN:chrB\tLN:800\n"
    )
    (folder / "target.fa.fai").write_text(
        files["fasta"]
        .with_name("ref.fa.fai")
        .read_text()
        .replace("chr1", "chrA")
        .replace("chr2", "chrB")
    )
    return {"chain": chain, "fasta": target, "dict": folder / "target.dict"}


@requires_picard
def test_liftover_vcf_and_haplotype_map_rename_the_contigs(
    files: dict[str, Path], tmp_path: Path
) -> None:
    lift = _lift_files(files, tmp_path)
    plain = tmp_path / "plain.vcf"
    plain.write_text((files["vcf"].parent / "sample.vcf").read_text())
    result = _ok(
        m.picard_liftover_vcf(
            str(plain),
            str(tmp_path / "lifted.vcf"),
            str(tmp_path / "rej.vcf"),
            str(lift["chain"]),
            str(lift["fasta"]),
            write_original_position=True,
            verbose=False,
        )
    )
    assert result["lifted_records"] == 8 and result["rejected_records"] == 0
    rows = _data_lines(Path(result["outputs"]["output"]))
    assert {row[0] for row in rows} == {"chrA"}
    assert [row[1] for row in rows] == [row[1] for row in _data_lines(plain)]
    assert all("OriginalContig=chr1" in row[7] for row in rows)
    hap = _ok(
        m.picard_liftover_haplotype_map(
            str(files["haplotype_map"]),
            str(tmp_path / "hap.map"),
            str(lift["chain"]),
            str(lift["dict"]),
            verbose=False,
        )
    )
    assert hap["snps"] == 8
    assert {row[0] for row in _data_lines(Path(hap["outputs"]["output"]))} == {"chrA"}


@requires_picard
def test_interval_list_to_bed_converts_coordinates(tmp_path: Path) -> None:
    source = _interval_list(
        tmp_path / "a.interval_list", [("chr2", 10, 50, "c"), ("chr1", 100, 200, "a")]
    )
    result = _ok(
        m.picard_interval_list_to_bed(str(source), str(tmp_path / "a.bed"), score=7, verbose=False)
    )
    assert _data_lines(Path(result["outputs"]["output"])) == [
        ["chr1", "99", "200", "a", "7", "+"],
        ["chr2", "9", "50", "c", "7", "+"],
    ]
    unsorted = _ok(
        m.picard_interval_list_to_bed(
            str(source), str(tmp_path / "u.bed"), sort=False, verbose=False
        )
    )
    assert [row[0] for row in _data_lines(Path(unsorted["outputs"]["output"]))] == ["chr2", "chr1"]


@requires_picard
def test_interval_list_tools_actions(tmp_path: Path) -> None:
    first = _interval_list(
        tmp_path / "a.interval_list",
        [("chr1", 100, 200, "a"), ("chr1", 300, 400, "b"), ("chr2", 10, 50, "c")],
    )
    second = _interval_list(tmp_path / "b.interval_list", [("chr1", 150, 180, "x")])
    concat = _ok(
        m.picard_interval_list_tools(
            [str(first), str(second)], str(tmp_path / "c.il"), verbose=False
        )
    )
    assert concat["intervals"] == 4
    overlaps = _ok(
        m.picard_interval_list_tools(
            [str(first)],
            str(tmp_path / "o.il"),
            action="OVERLAPS",
            second_input_paths=[str(second)],
            count_value="INTERVALS",
            verbose=False,
        )
    )
    assert [row[:3] for row in _data_lines(Path(overlaps["outputs"]["output"]))] == [
        ["chr1", "100", "200"]
    ]
    assert overlaps["count"] == 1
    subtract = _ok(
        m.picard_interval_list_tools(
            [str(first)],
            str(tmp_path / "s.il"),
            action="SUBTRACT",
            second_input_paths=[str(second)],
            count_value="BASES",
            verbose=False,
        )
    )
    assert subtract["count"] == 101 + 101 + 41 - 31  # bases of a minus the 31 bases of x
    padded = _ok(
        m.picard_interval_list_tools(
            [str(second)], str(tmp_path / "p.il"), action="UNION", padding=10, verbose=False
        )
    )
    assert _data_lines(Path(padded["outputs"]["output"]))[0][1:3] == ["140", "190"]
    inverted = _ok(
        m.picard_interval_list_tools(
            [str(second)],
            str(tmp_path / "i.il"),
            action="UNION",
            invert=True,
            count_value="BASES",
            verbose=False,
        )
    )
    assert inverted["count"] == 2000 - 31


@requires_picard
def test_scatter_intervals_by_ns_alternates_n_and_acgt(
    files: dict[str, Path], tmp_path: Path
) -> None:
    sequence = "ACGT" * 25 + "N" * 20 + "TTGCA" * 20
    reference = tmp_path / "gap.fa"
    reference.write_text(f">chrN\n{sequence}\n")
    (tmp_path / "gap.fa.fai").write_text(
        f"chrN\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n"
    )
    (tmp_path / "gap.dict").write_text(f"@HD\tVN:1.6\n@SQ\tSN:chrN\tLN:{len(sequence)}\n")
    result = _ok(
        m.picard_scatter_intervals_by_ns(str(reference), str(tmp_path / "n.il"), verbose=False)
    )
    rows = _data_lines(Path(result["outputs"]["output"]))
    assert [(row[1], row[2], row[4]) for row in rows] == [
        ("1", "100", "ACGTmer"),
        ("101", "120", "Nmer"),
        ("121", "220", "ACGTmer"),
    ]
    only_n = _ok(
        m.picard_scatter_intervals_by_ns(
            str(reference), str(tmp_path / "only.il"), output_type="N", verbose=False
        )
    )
    assert only_n["intervals"] == 1
    clean = _ok(
        m.picard_scatter_intervals_by_ns(
            str(files["fasta"]), str(tmp_path / "clean.il"), verbose=False
        )
    )
    assert [row[4] for row in _data_lines(Path(clean["outputs"]["output"]))] == [
        "ACGTmer",
        "ACGTmer",
    ]


@requires_picard
def test_normalize_fasta_sets_line_length_and_keeps_the_sequence(tmp_path: Path) -> None:
    source = tmp_path / "messy.fa"
    source.write_text(">s1 desc\n" + "ACGT" * 10 + "\nAC\n>s2\nACGTAC\nGT\n")
    result = _ok(
        m.picard_normalize_fasta(
            str(source),
            str(tmp_path / "n.fa"),
            line_length=16,
            truncate_sequence_names_at_whitespace=True,
            verbose=False,
        )
    )
    lines = Path(result["outputs"]["output"]).read_text().splitlines()
    assert lines == [
        "&gt;s1".replace("&gt;", ">"),
        "ACGTACGTACGTACGT",
        "ACGTACGTACGTACGT",
        "ACGTACGTAC",
        ">s2",
        "ACGTACGT",
    ]
    assert result["records"] == 2


@requires_picard
def test_sort_gff_orders_by_contig_and_position(tmp_path: Path) -> None:
    source = tmp_path / "x.gff3"
    source.write_text(
        "##gff-version 3\n"
        "chr2\tsrc\tgene\t10\t50\t.\t+\t.\tID=g2\n"
        "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g1b\n"
        "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1a\n"
    )
    result = _ok(m.picard_sort_gff(str(source), str(tmp_path / "sorted.gff3"), verbose=False))
    ids = [row[8] for row in _data_lines(Path(result["outputs"]["output"])) if len(row) > 8]
    assert ids == ["ID=g1a", "ID=g1b", "ID=g2"] and result["features"] >= 3
    dictionary = tmp_path / "order.dict"
    dictionary.write_text("@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:800\n@SQ\tSN:chr1\tLN:1200\n")
    ordered = _ok(
        m.picard_sort_gff(
            str(source),
            str(tmp_path / "ordered.gff3"),
            sequence_dictionary=str(dictionary),
            n_records_in_memory=2,
            verbose=False,
        )
    )
    ids = [row[8] for row in _data_lines(Path(ordered["outputs"]["output"])) if len(row) > 8]
    assert ids == ["ID=g2", "ID=g1a", "ID=g1b"]
