import random
import shutil

import pytest

from dnakit.assembly import hifiasm_p13_cli_functions as h
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs = pytest.mark.skipif(shutil.which("hifiasm") is None, reason="hifiasm not installed")


@pytest.fixture()
def reads(tmp_path):
    rng = random.Random(1)
    genome = "".join(rng.choice("ACGT") for _ in range(60000))
    comp = str.maketrans("ACGT", "TGCA")
    path = tmp_path / "hifi.fq"
    with path.open("w") as handle:
        for i in range(180):
            length = 8000
            start = rng.randrange(0, len(genome) - length)
            read = genome[start : start + length]
            if rng.random() < 0.5:
                read = read.translate(comp)[::-1]
            handle.write(f"@r{i}\n{read}\n+\n{'I' * length}\n")
    return path


def fake_program(tmp_path):
    script = tmp_path / "hifiasm"
    script.write_text(
        '#!/bin/sh\nif [ "$1" = "--version" ]; then echo 0.25.0-r1; exit 0; fi\necho "$@" > "$HIFIASM_ARGS"\nexit 0\n'
    )
    script.chmod(0o755)
    return script


def test_validation(tmp_path, reads):
    with pytest.raises(BackendUnavailableError):
        h.hifiasm_assemble(
            [str(reads)], str(tmp_path / "a"), hifiasm_path=str(tmp_path / "none"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        h.hifiasm_assemble(
            [str(reads)], str(tmp_path / "a"), hic_reads1_paths=[str(reads)], verbose=False
        )
    with pytest.raises(ConfigurationError):
        h.hifiasm_assemble(
            [str(reads)], str(tmp_path / "a"), paternal_kmers_path=str(reads), verbose=False
        )
    with pytest.raises(ConfigurationError):
        h.hifiasm_assemble([str(reads)], str(tmp_path / "a"), purge_level=7, verbose=False)


def test_argument_mapping_with_fake_program(tmp_path, reads, monkeypatch):
    monkeypatch.setenv("HIFIASM_ARGS", str(tmp_path / "args.txt"))
    fake = fake_program(tmp_path)
    result = h.hifiasm_assemble(
        [str(reads)],
        str(tmp_path / "o"),
        ont=True,
        ultra_long_paths=[str(reads)],
        hic_reads1_paths=[str(reads)],
        hic_reads2_paths=[str(reads)],
        dual_scaffold=True,
        purge_level=2,
        write_error_corrected_reads=True,
        write_overlaps_paf=True,
        low_quality_percent=50,
        break_low_coverage=3,
        break_high_coverage=100,
        misjoin_min_size=1000,
        telomere_motif="CCCTAA",
        chimeric_read_support=1,
        min_ont_read_length=500,
        min_ultra_long_length=200,
        hifiasm_path=str(fake),
        verbose=False,
    )
    assert result["status"] == "failed"  # no assembly graph was produced by the fake program
    args = (tmp_path / "args.txt").read_text().split()
    for flag in (
        "--ont",
        "--ul",
        "--h1",
        "--h2",
        "--dual-scaf",
        "-l",
        "--write-ec",
        "--write-paf",
        "--lowQ",
        "--b-cov",
        "--h-cov",
        "--l-msjoin",
        "--telo-m",
        "--chem-c",
        "--rl-cut",
        "--ul-cut",
    ):
        assert flag in args


@needs
def test_real_assembly_and_options(tmp_path, reads):
    result = h.hifiasm_assemble(
        [str(reads)],
        str(tmp_path / "asm"),
        bloom_filter_bits=0,
        threads=2,
        write_error_corrected_reads=True,
        write_overlaps_paf=True,
        low_quality_percent=50,
        verbose=False,
    )
    assert result["status"] == "completed" and result["version"]
    assert result["error_corrected_reads"] and result["overlaps_paf"]
    assert any(item["total_bp"] > 50000 for item in result["contig_summary"])
    with pytest.raises(ConfigurationError):
        h.hifiasm_assemble([str(reads)], str(tmp_path / "asm"), bloom_filter_bits=0, verbose=False)
