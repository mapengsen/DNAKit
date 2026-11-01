"""Filtlong, NanoFilt and Porechop facades; the time filter runs without a program."""

from __future__ import annotations

import random
import stat
from pathlib import Path

import pytest
from _p16_data import genome, long_reads, need, rand_dna, write_fastq

import dnakit.sequencing_qc.p16_readqc_cli_functions as rq
from dnakit.exceptions import ConfigurationError


def _count(path: str) -> int:
    return sum(1 for line in Path(path).read_text().splitlines() if line.startswith("@lr"))


def test_filtlong_with_reference_reads_trim_and_split(tmp_path: Path) -> None:
    program = need("filtlong")
    contigs = genome(tmp_path)
    plain, _ = long_reads(tmp_path, contigs)
    rng = random.Random(1)
    shorts = []
    for i in range(400):
        name = rng.choice(list(contigs))
        start = rng.randrange(len(contigs[name]) - 100)
        shorts.append((f"s{i}", contigs[name][start : start + 100]))
    short = write_fastq(tmp_path / "short.fastq", shorts)
    plain_result = rq.filter_long_reads_with_reference(plain, str(tmp_path / "a.fq"), min_length=500, filtlong_path=program, verbose=False)
    ref_result = rq.filter_long_reads_with_reference(plain, str(tmp_path / "b.fq"), illumina_reads_1=short, trim=True, min_length=500,
                                                     filtlong_path=program, verbose=False)
    assert plain_result["reads_kept"] == 40 and 0 < ref_result["reads_kept"] <= 40
    split = rq.filter_long_reads_with_reference(plain, str(tmp_path / "c.fq"), assembly=str(tmp_path / "genome.fa"), split=100,
                                                filtlong_path=program, verbose=False)
    assert split["reads_kept"] >= 1
    with pytest.raises(ConfigurationError):
        rq.filter_long_reads_with_reference(plain, str(tmp_path / "d.fq"), illumina_reads_1=short, assembly=str(tmp_path / "genome.fa"), verbose=False)


def test_nanofilt_gc_filter(tmp_path: Path) -> None:
    program = need("NanoFilt")
    rng = random.Random(2)
    records = [("gc_low", "".join(rng.choice("AT") for _ in range(200))), ("gc_mid", "ACGT" * 50), ("gc_high", "".join(rng.choice("GC") for _ in range(200)))]
    source = write_fastq(tmp_path / "r.fastq", records)
    result = rq.filter_nanopore_reads_by_gc_summary(source, str(tmp_path / "o.fq"), min_gc=0.3, max_gc=0.7, nanofilt_path=program, verbose=False)
    assert result["reads_kept"] == 1 and (tmp_path / "o.fq").read_text().startswith("@gc_mid")
    with pytest.raises(ConfigurationError):
        rq.filter_nanopore_reads_by_gc_summary(source, str(tmp_path / "p.fq"), verbose=False)


def test_start_time_filter(tmp_path: Path) -> None:
    pytest.importorskip("dateutil")
    contigs = genome(tmp_path)
    _, zipped = long_reads(tmp_path, contigs)
    result = rq.filter_nanopore_reads_by_start_time(zipped, str(tmp_path / "t.fq"), time_from="2024-01-02T00:00:00Z", time_to="2024-01-03T23:00:00Z", verbose=False)
    assert result["reads_total"] == 40 and result["reads_kept"] == 16
    kept = (tmp_path / "t.fq").read_text()
    assert "start_time=2024-01-01" not in kept and "start_time=2024-01-04" not in kept
    nothing = rq.filter_nanopore_reads_by_start_time(zipped, str(tmp_path / "n.fq"), time_from="2030-01-01T00:00:00Z", time_to="2030-01-02T00:00:00Z", verbose=False)
    assert nothing["reads_kept"] == 0


def test_porechop_arguments_with_recording_shim(tmp_path: Path) -> None:
    """Porechop is archived and its conda build is unreliable; the command assembly is checked with a recording shim."""
    shim = tmp_path / "porechop"
    shim.write_text(
        "#!/bin/sh\nprintf '%s\\n' \"$@\" >> \"$(dirname \"$0\")/argv.txt\"\n"
        "out=''; bins=''\nwhile [ $# -gt 0 ]; do case \"$1\" in -o) out=$2; shift;; -b) bins=$2; shift;; esac; shift; done\n"
        "[ -n \"$out\" ] && : > \"$out\"\n[ -n \"$bins\" ] && mkdir -p \"$bins\" && : > \"$bins/BC01.fastq\"\nexit 0\n"
    )
    shim.chmod(shim.stat().st_mode | stat.S_IEXEC)
    reads = write_fastq(tmp_path / "r.fq", [("lr0", rand_dna(random.Random(1), 50))])
    rq.trim_nanopore_adapters_discard_middle(reads, str(tmp_path / "out.fq"), porechop_path=str(shim), verbose=False)
    assert "--discard_middle" in (tmp_path / "argv.txt").read_text().split()
    result = rq.demultiplex_nanopore_barcodes(reads, str(tmp_path / "bins"), require_two_barcodes=True, porechop_path=str(shim), verbose=False)
    argv = (tmp_path / "argv.txt").read_text().split()
    assert "-b" in argv and "--require_two_barcodes" in argv and "BC01.fastq" in result["reads_per_bin"]
