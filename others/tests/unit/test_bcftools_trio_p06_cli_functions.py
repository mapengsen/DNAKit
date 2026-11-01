"""Tests for dnakit.variants.bcftools_trio_p06_cli_functions (bcftools trio and family plugin facades).

The toy trio has 40 sites that cycle through seven genotype patterns; patterns 4 and 6 break Mendelian
inheritance, which makes 11 inconsistent sites (indices 4, 6, 11, 13, 18, 20, 25, 27, 32, 34, 39).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from _p06_htslib_data import (
    PED_TEXT,
    body,
    requires_bcftools,
    requires_plugins,
    trio_vcf_text,
    write_vcf,
)

from dnakit.exceptions import ConfigurationError
from dnakit.variants import bcftools_trio_p06_cli_functions as m

pytestmark = [requires_bcftools, requires_plugins]

ERROR_SITES = 11


def _direct(*arguments: str) -> str:
    return subprocess.run(["bcftools", *arguments], check=True, capture_output=True, text=True).stdout


@pytest.fixture
def trio(tmp_path: Path) -> Path:
    return write_vcf(tmp_path / "trio.vcf", trio_vcf_text(), compress=True, index=True)


@pytest.fixture
def ped(tmp_path: Path) -> Path:
    path = tmp_path / "trio.ped"
    path.write_text(PED_TEXT, encoding="utf-8")
    return path


def _rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in body(path)]


def _format_values(row: list[str], key: str, sample: int = 9) -> str:
    return row[sample].split(":")[row[8].split(":").index(key)]


def test_parental_origin_of_a_duplication(tmp_path: Path) -> None:
    lines = [
        "##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\tF\tM",
    ]
    for index in range(30):
        lines.append(
            f"chr1\t{1000 + 100 * index}\t.\tA\tG\t50\tPASS\t.\tGT:AD:PL:DP\t"
            "0/1:10,20:60,0,60:30\t0/1:15,15:60,0,60:30\t0/0:30,0:0,60,255:30"
        )
    source = write_vcf(tmp_path / "cnv.vcf", "\n".join(lines) + "\n", compress=True, index=True)
    result = m.bcftools_parental_origin(
        source, tmp_path / "origin.txt", trio=["P", "F", "M"], cnv_type="dup", region="chr1:1-5000", verbose=False
    )
    assert result["predicted_origin"] == "paternal" and result["markers"] == 30
    direct = _direct("+parental-origin", "-p", "P,F,M", "-t", "dup", "-r", "chr1:1-5000", "-b", "0.01", str(source))
    assert [r for r in (tmp_path / "origin.txt").read_text(encoding="utf-8").splitlines() if not r.startswith("#")] == [
        r for r in direct.splitlines() if not r.startswith("#")
    ]
    with pytest.raises(ConfigurationError):
        m.bcftools_parental_origin(
            source, tmp_path / "bad.txt", trio=["P", "F"], cnv_type="dup", region="chr1:1-5000", verbose=False
        )
    with pytest.raises(ConfigurationError):
        m.bcftools_parental_origin(
            source, tmp_path / "bad.txt", trio=["P", "F", "M"], cnv_type="gain", region="chr1:1-5000", verbose=False
        )


def test_trio_stats_counts_mendelian_errors(tmp_path: Path, trio: Path, ped: Path) -> None:
    from_ped = m.bcftools_trio_stats(trio, tmp_path / "ped.txt", ped_path=ped, verbose=False)
    from_names = m.bcftools_trio_stats(trio, tmp_path / "names.txt", trio=["C", "F", "M"], verbose=False)
    first = from_ped["trios"][0]
    assert (first["child"], first["father"], first["mother"]) == ("C", "F", "M")
    assert first["dnm_or_mendel_errors"] == ERROR_SITES
    assert from_names["trios"][0]["dnm_or_mendel_errors"] == ERROR_SITES
    assert first["valid_gts"] == 40
    ranged = m.bcftools_trio_stats(trio, tmp_path / "range.txt", ped_path=ped, include="GQ>{10,20}", verbose=False)
    assert len(ranged["trios"]) == 2
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_stats(trio, tmp_path / "bad.txt", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_stats(trio, tmp_path / "bad.txt", ped_path=ped, trio=["C", "F", "M"], verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_stats(trio, tmp_path / "bad.txt", ped_path=ped, debug=["nothing"], verbose=False)


def test_trio_switch_rate_on_phased_genotypes(tmp_path: Path, ped: Path) -> None:
    source = write_vcf(tmp_path / "phased.vcf", trio_vcf_text(phased=True), compress=True, index=True)
    result = m.bcftools_trio_switch_rate(source, tmp_path / "switch.txt", ped_path=ped, verbose=False)
    row = result["trios"][0]
    assert (row["father"], row["mother"], row["child"]) == ("F", "M", "C")
    assert row["tested"] == 12 and row["mendelian_errors"] == 6 and row["switches"] == 0
    direct = _direct("+trio-switch-rate", str(source), "--", "-p", str(ped))
    assert [r for r in (tmp_path / "switch.txt").read_text(encoding="utf-8").splitlines() if r.startswith("TRIO")] == [
        r for r in direct.splitlines() if r.startswith("TRIO")
    ]


def test_mendelian_counts_and_vcf_modes(tmp_path: Path, trio: Path, ped: Path) -> None:
    counts = m.bcftools_mendelian_consistency(trio, tmp_path / "counts.txt", ped_path=ped, verbose=False)
    summary = {row[0]: row[1] for row in counts["report"]["rows"] if isinstance(row[1], int)}
    assert summary["sites_merr"] == ERROR_SITES and summary["sites_good"] == 40 - ERROR_SITES
    annotated = m.bcftools_mendelian_consistency(
        trio, tmp_path / "annotated.vcf", trio=["C", "F", "M"], mode="a", verbose=False
    )
    assert annotated["records"] == 40
    assert sum("MERR=1" in row[7] for row in _rows(tmp_path / "annotated.vcf")) == ERROR_SITES
    errors = m.bcftools_mendelian_consistency(trio, tmp_path / "e.vcf", ped_path=ped, mode="e", verbose=False)
    clean = m.bcftools_mendelian_consistency(trio, tmp_path / "E.vcf", ped_path=ped, mode="E", verbose=False)
    assert (errors["records"], clean["records"]) == (ERROR_SITES, 40 - ERROR_SITES)
    deleted = m.bcftools_mendelian_consistency(trio, tmp_path / "d.vcf", ped_path=ped, mode="d", verbose=False)
    assert deleted["records"] == 40
    assert sum(row[9].startswith("./.") for row in _rows(tmp_path / "d.vcf")) == ERROR_SITES
    direct = _direct("+mendelian2", str(trio), "-P", str(ped), "--rules", "GRCh37", "-m", "E")
    assert body(tmp_path / "E.vcf") == [line for line in direct.splitlines() if not line.startswith("#")]
    with pytest.raises(ConfigurationError):
        m.bcftools_mendelian_consistency(trio, tmp_path / "bad.vcf", ped_path=ped, mode="xyz", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_mendelian_consistency(trio, tmp_path / "bad.vcf", mode="a", verbose=False)


def test_trio_denovo_naive_flags_the_inconsistent_sites(tmp_path: Path, trio: Path) -> None:
    result = m.bcftools_trio_denovo(trio, tmp_path / "naive.vcf", trio=["C", "F", "M"], model="naive", verbose=False)
    flagged = [row for row in _rows(tmp_path / "naive.vcf") if "DNM" in row[8].split(":")]
    assert result["model"] == "naive" and len(flagged) == ERROR_SITES
    direct = _direct("+trio-dnm3", "-p", "C,F,M", "--use-NAIVE", "--dnm-tag", "DNM:flag", "-X", "GRCh37", str(trio))
    assert result["records"] == len([line for line in direct.splitlines() if not line.startswith("#")])


@pytest.mark.parametrize("model", ["dng", "alm", "dmm"])
def test_trio_denovo_scores_every_site(tmp_path: Path, trio: Path, model: str) -> None:
    extra = {"max_qm": -1.0} if model == "dmm" else {}
    result = m.bcftools_trio_denovo(
        trio, tmp_path / f"{model}.vcf", trio=["C", "F", "M"], model=model, **extra, verbose=False
    )
    rows = _rows(tmp_path / f"{model}.vcf")
    assert result["records"] == 40 and all("DNM" in row[8].split(":") for row in rows)
    scores = {int(row[1]): float(_format_values(row, "DNM")) for row in rows}
    assert len(scores) == 40
    if model in {"dng", "alm"}:
        # the plain de novo site (child 0/1, both parents 0/0) outscores every transmitted heterozygous site
        assert scores[1400] > max(scores[1000], scores[1200])
    direct = _direct(
        "+trio-dnm3", "-p", "C,F,M", f"--use-{model.upper()}", "--dnm-tag", "DNM:log", "-X", "GRCh37",
        "--mrate", "1e-08", *(["--max-QM", "-1.0"] if model == "dmm" else []), str(trio),
    )
    assert body(tmp_path / f"{model}.vcf") == [line for line in direct.splitlines() if not line.startswith("#")]


def test_trio_denovo_option_validation(tmp_path: Path, trio: Path) -> None:
    out = tmp_path / "x.vcf"
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_denovo(trio, out, trio=["C", "F", "M"], model="magic", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_denovo(trio, out, trio=["C", "F", "M"], model="dmm", allele_dropout=0.1, verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_denovo(trio, out, trio=["C", "F", "M"], model="alm", max_qm=30, verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_denovo(trio, out, trio=["C", "F", "M"], dnm_tag="bad tag", verbose=False)
    with pytest.raises(ConfigurationError):
        m.bcftools_trio_denovo(trio, out, verbose=False)
    assert not out.exists()
