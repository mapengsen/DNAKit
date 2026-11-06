from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.assembly import megahit_p15_cli_functions as mh
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs_megahit = pytest.mark.skipif(
    shutil.which("megahit_core") is None, reason="megahit_core is not installed"
)


def reads_file(path: Path, count: int = 50, length: int = 40) -> Path:
    random.seed(3)
    path.write_text(
        "".join(
            f"@r{i}\n{''.join(random.choice('ACGT') for _ in range(length))}\n+\n{'I' * length}\n"
            for i in range(count)
        )
    )
    return path


def test_missing_program_and_validation(tmp_path: Path) -> None:
    reads = reads_file(tmp_path / "r.fq")
    with pytest.raises(BackendUnavailableError):
        mh.megahit_read_statistics(
            str(reads),
            str(tmp_path / "s.txt"),
            megahit_core_path=str(tmp_path / "missing"),
            verbose=False,
        )
    contigs = tmp_path / "c.fa"
    contigs.write_text(">c\nACGT\n")
    with pytest.raises(ConfigurationError):
        mh.megahit_local_assembly(str(contigs), str(tmp_path / "o"), verbose=False)  # no reads
    with pytest.raises(ConfigurationError):
        mh.megahit_local_assembly(
            str(contigs), str(tmp_path / "o2"), paired_reads=[[str(reads)]], verbose=False
        )


@needs_megahit
def test_read_statistics(tmp_path: Path) -> None:
    result = mh.megahit_read_statistics(
        str(reads_file(tmp_path / "r.fq")), str(tmp_path / "s.txt"), verbose=False
    )
    assert (
        result["read_count"],
        result["total_bases"],
        result["longest"],
        result["shortest"],
        result["mean_length"],
    ) == (50, 2000, 40, 40, 40.0)


@needs_megahit
def test_fastg_and_local_assembly(tmp_path: Path) -> None:
    random.seed(5)
    genome = "".join(random.choice("ACGT") for _ in range(6000))
    contigs = tmp_path / "k21.contigs.fa"
    contigs.write_text(f">k21_0 flag=1 multi=10.0 len=2000\n{genome[:2000]}\n")
    fastg = mh.megahit_contigs_to_fastg(str(contigs), 21, str(tmp_path / "g.fastg"), verbose=False)
    assert fastg["node_count"] == 2 and "NODE_1_length_2000" in (tmp_path / "g.fastg").read_text()
    first, second = tmp_path / "p1.fq", tmp_path / "p2.fq"
    rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))  # noqa: E731
    with first.open("w") as f1, second.open("w") as f2:
        for i in range(3000):
            insert = random.randint(300, 400)
            start = random.randrange(0, len(genome) - insert)
            frag = genome[start : start + insert]
            f1.write(f"@r{i}/1\n{frag[:100]}\n+\n{'I' * 100}\n")
            f2.write(f"@r{i}/2\n{rc(frag)[:100]}\n+\n{'I' * 100}\n")
    split = tmp_path / "split.fa"
    split.write_text(f">a\n{genome[:2500]}\n>b\n{genome[3000:5500]}\n")
    result = mh.megahit_local_assembly(
        str(split),
        str(tmp_path / "local"),
        paired_reads=[[str(first), str(second)]],
        min_contig_len=100,
        verbose=False,
    )
    assert result["input_contig_count"] == 2 and result["local_contig_count"] >= 1
    assert result["local_contig_bases"] > 0 and Path(result["local_contig_file"]).is_file()
