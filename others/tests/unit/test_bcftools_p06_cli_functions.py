"""Tests for dnakit.variants.bcftools_p06_cli_functions (bcftools convert, index, norm and polysomy facades)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest
from _p06_htslib_data import DATA, body, data_vcf_gz, requires_bcftools, requires_data, write_vcf

from dnakit.exceptions import BackendUnavailableError, ConfigurationError
from dnakit.variants import bcftools_p06_cli_functions as m

pytestmark = [requires_bcftools, requires_data]


def _bcftools(*arguments: str) -> str:
    return subprocess.run(
        ["bcftools", *arguments], check=True, capture_output=True, text=True
    ).stdout


def _small_vcf(tmp_path: Path) -> Path:
    """Nine biallelic SNVs of three samples, from the shared VCF."""
    source = data_vcf_gz(tmp_path)
    target = tmp_path / "small.vcf"
    subprocess.run(
        [
            "bcftools",
            "view",
            "-r",
            "chr1:1-5000",
            "-s",
            "S00,S01,S02",
            "-v",
            "snps",
            "-o",
            str(target),
            str(source),
        ],
        check=True,
        capture_output=True,
    )
    return target


def _genotypes(path: Path) -> list[list[str]]:
    return [line.split("\t")[9:] for line in body(path)]


@pytest.mark.parametrize(
    ("converter_flag", "suffixes", "input_format"),
    [
        ("-g", ("gen.gz", "samples"), "gensample"),
        ("--hapsample", ("hap.gz", "samples"), "hapsample"),
        ("-h", ("hap.gz", "legend.gz", "samples"), "haplegendsample"),
    ],
)
def test_convert_text_formats_back_to_the_source_genotypes(
    tmp_path: Path, converter_flag: str, suffixes: tuple[str, ...], input_format: str
) -> None:
    small = _small_vcf(tmp_path)
    prefix = tmp_path / "panel"
    subprocess.run(
        ["bcftools", "convert", converter_flag, str(prefix), str(small)],
        check=True,
        capture_output=True,
    )
    files = [f"{prefix}.{suffix}" for suffix in suffixes]
    out = tmp_path / "back.vcf"
    result = m.bcftools_convert_to_vcf(input_format, files, out, verbose=False)
    assert result["status"] == "completed" and result["records"] == 9
    expected = [
        [g.split(":")[0].replace("1/0", "0/1").replace("./.", "0/0") for g in row]
        for row in _genotypes(small)
    ]
    got = [
        [g.split(":")[0].replace("1/0", "0/1").replace("|", "/") for g in row]
        for row in _genotypes(out)
    ]
    if input_format == "gensample":
        assert got == expected
    else:
        assert sorted(sorted(g.replace(".", "0").split("/")) for row in got for g in row) == sorted(
            sorted(g.split("/")) for row in expected for g in row
        )
    direct = _bcftools(
        "convert",
        "-G"
        if input_format == "gensample"
        else "--hapsample2vcf"
        if input_format == "hapsample"
        else "-H",
        ",".join(files),
    )
    assert body(out) == [line for line in direct.splitlines() if not line.startswith("#")]


def test_convert_tsv_and_gvcf(tmp_path: Path) -> None:
    reference = DATA / "genome.fa"
    if not reference.exists():
        pytest.skip("genome.fa missing")
    raw = tmp_path / "raw.tsv"
    raw.write_text("rs1\tchr1\t3\tAC\nrs2\tchr1\t5\tTT\nrs3\tchr1\t7\t--\n", encoding="utf-8")
    out = tmp_path / "tsv.vcf"
    result = m.bcftools_convert_to_vcf(
        "tsv", [raw], out, reference_fasta=reference, sample_names=["SAMP1"], verbose=False
    )
    assert result["records"] == 3
    assert result["summary"]["Rows total"] == 3 and result["summary"]["Missing GTs"] == 1
    assert [line.split("\t")[9] for line in body(out)] == ["1/2", "0/0", "./."]
    gvcf = tmp_path / "g.vcf"
    gvcf.write_text(
        '##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##contig=<ID=chr1,length=120000>\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1\t.\tT\t<*>\t.\t.\tEND=5\tGT\t0/0\nchr1\t6\t.\tT\tC,<*>\t50\t.\t.\tGT\t1/1\n",
        encoding="utf-8",
    )
    expanded = m.bcftools_convert_to_vcf(
        "gvcf", [gvcf], tmp_path / "e.vcf", reference_fasta=reference, verbose=False
    )
    assert expanded["records"] == 6
    assert [line.split("\t")[1] for line in body(tmp_path / "e.vcf")] == [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
    ]


def test_convert_validation(tmp_path: Path) -> None:
    some = tmp_path / "a.txt"
    some.write_text("x", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_to_vcf("nonsense", [some], tmp_path / "o.vcf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_to_vcf("gensample", [some], tmp_path / "o.vcf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_to_vcf("tsv", [some], tmp_path / "o.vcf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_to_vcf(
            "hapsample", [some, some], tmp_path / "o.vcf", columns="ID", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_convert_to_vcf(
            "gensample", [some, some], tmp_path / "o.vcf", haploid_to_diploid=True, verbose=False
        )


def test_index_records_per_contig_and_total(tmp_path: Path) -> None:
    source = data_vcf_gz(tmp_path)
    detailed = m.bcftools_index_records(source, verbose=False)
    direct = [line.split("\t") for line in _bcftools("index", "-s", str(source)).splitlines()]
    assert [(r["contig"], r["length"], r["records"]) for r in detailed["contigs"]] == [
        (c, int(n), int(r)) for c, n, r in direct
    ]
    assert detailed["total_records"] == 400
    total = m.bcftools_index_records(source, per_contig=False, verbose=False)
    assert total["total_records"] == 400 and "contigs" not in total
    with pytest.raises(ConfigurationError):
        m.bcftools_index_records(
            source, per_contig=False, include_empty_contigs=True, verbose=False
        )


def test_index_records_all_contigs_option_matches_bcftools(tmp_path: Path) -> None:
    source = data_vcf_gz(tmp_path)
    with_all = m.bcftools_index_records(source, include_empty_contigs=True, verbose=False)
    assert "-a" in with_all["command"]
    direct = [line.split("\t") for line in _bcftools("index", "-s", "-a", str(source)).splitlines()]
    assert [(r["contig"], r["records"]) for r in with_all["contigs"]] == [
        (c, int(r)) for c, _, r in direct
    ]


def test_norm_left_aligns_and_matches_bcftools(tmp_path: Path) -> None:
    reference = DATA / "genome.fa"
    source = data_vcf_gz(tmp_path)
    out = tmp_path / "norm.vcf"
    result = m.bcftools_norm_variants(
        source, out, reference_fasta=reference, check_ref="w", verbose=False
    )
    assert result["records"] == 400
    direct = _bcftools("norm", "-f", str(reference), "-c", "w", str(source))
    assert body(out) == [line for line in direct.splitlines() if not line.startswith("#")]


def test_norm_splits_atomizes_and_removes_duplicates(tmp_path: Path) -> None:
    text = (
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tC,G\t.\t.\t.\tGT\t1/2\n"
        "chr1\t20\t.\tAC\tGT\t.\t.\t.\tGT\t0/1\n"
        "chr1\t30\t.\tA\tT\t.\t.\t.\tGT\t0/1\n"
        "chr1\t30\t.\tA\tT\t.\t.\t.\tGT\t0/1\n"
    )
    source = write_vcf(tmp_path / "multi.vcf", text, compress=True, index=True)
    split = m.bcftools_norm_variants(
        source, tmp_path / "split.vcf", multiallelics="-any", verbose=False
    )
    assert split["records"] == 5 and split["summary"] == {}
    assert [line.split("\t")[4] for line in body(tmp_path / "split.vcf")][:2] == ["C", "G"]
    atoms = m.bcftools_norm_variants(source, tmp_path / "atoms.vcf", atomize=True, verbose=False)
    assert [line.split("\t")[1:5] for line in body(tmp_path / "atoms.vcf")][2:4] == [
        ["20", ".", "A", "G"],
        ["21", ".", "C", "T"],
    ]
    assert atoms["records"] == 6
    dedup = m.bcftools_norm_variants(
        source, tmp_path / "dedup.vcf", remove_duplicates="exact", verbose=False
    )
    assert dedup["records"] == 3
    with pytest.raises(ConfigurationError):
        m.bcftools_norm_variants(source, tmp_path / "none.vcf", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_norm_variants(source, tmp_path / "bad.vcf", multiallelics="x", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_norm_variants(
            source, tmp_path / "bad2.vcf", check_ref="e", atomize=True, verbose=False
        )


def test_polysomy_reports_a_build_without_gsl(tmp_path: Path) -> None:
    fake = tmp_path / "bcftools"
    fake.write_text(
        "#!/bin/sh\necho \"[E::main] unrecognized command 'polysomy'\" >&2\nexit 1\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    source = data_vcf_gz(tmp_path)
    with pytest.raises(BackendUnavailableError):
        m.bcftools_polysomy(source, tmp_path / "poly", bcftools_path=str(fake), verbose=False)
    assert not (tmp_path / "poly").exists()
    with pytest.raises(ConfigurationError):
        m.bcftools_polysomy(source, tmp_path / "poly2", cn_penalty=-1.0, verbose=False)


def test_missing_program_has_an_install_hint(tmp_path: Path) -> None:
    source = data_vcf_gz(tmp_path)
    with pytest.raises(BackendUnavailableError) as raised:
        m.bcftools_index_records(source, bcftools_path=str(tmp_path / "absent"), verbose=False)
    assert "bcftools" in str(raised.value)


def _polysomy_available() -> bool:
    # Version note: avoid a collection-time failure when optional bcftools is absent.
    if shutil.which("bcftools") is None:
        return False
    done = subprocess.run(["bcftools", "polysomy"], capture_output=True, text=True, check=False)
    return "unrecognized command" not in done.stderr


def _baf_vcf(path: Path) -> Path:
    """Chromosome 1 with peaks at BAF 0, 0.5, 1 (two copies) and chromosome 2 with 0, 1/3, 2/3, 1 (three copies)."""
    import random

    rng = random.Random(3)
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1,length=1000000>",
        "##contig=<ID=2,length=1000000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
        '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAF">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    for chrom, centers in (("1", [0.0, 0.5, 1.0]), ("2", [0.0, 1 / 3, 2 / 3, 1.0])):
        for index in range(2500):
            center = rng.choice(centers)
            if center == 0.0:
                value, genotype = abs(rng.gauss(0, 0.02)), "0/0"
            elif center == 1.0:
                value, genotype = 1 - abs(rng.gauss(0, 0.02)), "1/1"
            else:
                value, genotype = rng.gauss(center, 0.025), "0/1"
            lines.append(
                f"{chrom}\t{1000 + index * 100}\t.\tA\tG\t.\t.\t.\tGT:BAF\t{genotype}:{min(1.0, max(0.0, value)):.4f}"
            )
    return write_vcf(path, "\n".join(lines) + "\n", compress=True, index=True)


@pytest.mark.skipif(
    not _polysomy_available(), reason="this bcftools was built without GSL (no polysomy command)"
)
def test_polysomy_estimates_copy_numbers(tmp_path: Path) -> None:
    source = _baf_vcf(tmp_path / "baf.vcf")
    result = m.bcftools_polysomy(source, tmp_path / "poly", sample="S1", verbose=False)
    assert {row["chromosome"]: row["copy_number"] for row in result["copy_numbers"]} == {
        "1": 2,
        "2": 3,
    }
    assert sorted(result["files"]) == ["dist.dat", "dist.py"]
    with pytest.raises(ConfigurationError):
        m.bcftools_polysomy(source, tmp_path / "poly", verbose=False)
