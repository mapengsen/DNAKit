"""Shared helpers and toy data for the pack P06 samtools / bcftools facade tests."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

DATA = (
    Path(__file__).resolve().parents[3]
    / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
)

HAVE_SAMTOOLS = shutil.which("samtools") is not None
HAVE_BCFTOOLS = shutil.which("bcftools") is not None


def have_plugins() -> bool:
    """True when bcftools is on PATH and can load its plugins (``bcftools plugin -l`` lists fixref)."""
    exe = shutil.which("bcftools")
    if exe is None:
        return False
    from dnakit.bioinformatics import _htslib_p06 as h

    try:
        done = h.invoke("bcftools", None, ["plugin", "-l"], plugin=True)
    except Exception:
        return False
    return "fixref" in done.stdout


requires_samtools = pytest.mark.skipif(not HAVE_SAMTOOLS, reason="samtools is not on PATH")
requires_bcftools = pytest.mark.skipif(not HAVE_BCFTOOLS, reason="bcftools is not on PATH")
requires_plugins = pytest.mark.skipif(
    not have_plugins(), reason="bcftools plugins are not available"
)
requires_data = pytest.mark.skipif(not DATA.is_dir(), reason="shared consistency/data is missing")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-scaled quality sum">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Strand bias">
"""


def _call(genotype: str, depth: int = 30) -> str:
    first, second = genotype.replace("|", "/").split("/")
    alt = (first == "1") + (second == "1")
    ref_reads = depth * (2 - alt) // 2
    alt_reads = depth - ref_reads
    pl = {0: "0,60,255", 1: "60,0,60", 2: "255,60,0"}[alt]
    return f"{genotype}:{pl}:{ref_reads},{alt_reads}:{depth}:60:{ref_reads * 38},{alt_reads * 38}:0"


def trio_rows(count: int = 40, *, phased: bool = False) -> list[tuple[str, str, str, str]]:
    """Genotypes (child, father, mother, unrelated) cycling through 7 patterns, 3 of 7 Mendelian errors or de novo."""
    sep = "|" if phased else "/"
    patterns = [
        (f"0{sep}1", f"0{sep}1", f"0{sep}0", f"0{sep}0"),
        (f"1{sep}1", f"0{sep}1", f"0{sep}1", f"0{sep}1"),
        (f"0{sep}1", f"1{sep}1", f"0{sep}0", f"0{sep}0"),
        (f"0{sep}0", f"0{sep}0", f"0{sep}0", f"0{sep}1"),
        (f"0{sep}1", f"0{sep}0", f"0{sep}0", f"0{sep}0"),
        (f"0{sep}0", f"0{sep}1", f"0{sep}1", f"1{sep}1"),
        (f"1{sep}1", f"0{sep}0", f"0{sep}1", f"0{sep}0"),
    ]
    return [patterns[i % len(patterns)] for i in range(count)]


def trio_vcf_text(count: int = 40, *, phased: bool = False, indels: bool = False) -> str:
    """A proband C, father F, mother M and unrelated U on chr1 with GT, PL, AD, DP, GQ, QS and SP.

    With ``indels`` every third site is a 2-bp insertion or a 3-bp deletion.
    """
    lines = [HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tF\tM\tU"]
    for index, genotypes in enumerate(trio_rows(count, phased=phased)):
        ref, alt = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")][index % 4]
        if indels and index % 3 == 0:
            ref, alt = (ref, ref + "CG") if index % 2 == 0 else (ref + "CGT", ref)
        fields = "\t".join(_call(genotype) for genotype in genotypes)
        lines.append(
            f"chr1\t{1000 + 100 * index}\t.\t{ref}\t{alt}\t50\tPASS\tDP=120\tGT:PL:AD:DP:GQ:QS:SP\t{fields}"
        )
    return "\n".join(lines) + "\n"


PED_TEXT = "fam1\tC\tF\tM\t1\t0\n"


def write_vcf(path: Path, text: str, *, compress: bool = False, index: bool = False) -> Path:
    """Write VCF text; with ``compress`` bgzip it through bcftools, with ``index`` add a CSI index."""
    path.write_text(text, encoding="utf-8")
    if not compress:
        return path
    target = path.with_name(path.name + ".gz")
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(target), str(path)], check=True, capture_output=True
    )
    if index:
        subprocess.run(["bcftools", "index", str(target)], check=True, capture_output=True)
    return target


def data_vcf_gz(tmp_path: Path) -> Path:
    """The shared 24-sample variants.vcf as an indexed bgzip file in ``tmp_path``."""
    return write_vcf(
        tmp_path / "variants.vcf",
        (DATA / "variants.vcf").read_text(encoding="utf-8"),
        compress=True,
        index=True,
    )


def body(path: Path) -> list[str]:
    """Record lines of a plain-text VCF."""
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
