from __future__ import annotations

import random
import shutil
from pathlib import Path

import pytest

from dnakit.assembly import quast_p15_cli_functions as qu
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs_quast = pytest.mark.skipif(shutil.which("quast.py") is None, reason="quast.py is not installed")


def genome_files(folder: Path) -> tuple[Path, Path]:
    random.seed(21)
    genome = "".join(random.choice("ACGT") for _ in range(20000))
    reference = folder / "ref.fa"
    reference.write_text(f">ref\n{genome}\n")
    contigs = folder / "contigs.fa"
    contigs.write_text(f">c1\n{genome[:9000]}\n>c2\n{genome[9500:]}\n")
    return reference, contigs


def test_missing_program_and_validation(tmp_path: Path) -> None:
    reference, contigs = genome_files(tmp_path)
    with pytest.raises(BackendUnavailableError):
        qu.quast_assess_assembly([str(contigs)], str(tmp_path / "o"), reference=str(reference), quast_path=str(tmp_path / "missing"), verbose=False)
    with pytest.raises(ConfigurationError):
        qu.quast_assess_assembly([str(contigs)], str(tmp_path / "o1"), upper_bound_assembly=True, verbose=False)
    with pytest.raises(ConfigurationError):
        qu.quast_predict_genes([str(contigs)], str(tmp_path / "o2"), gene_finder=None, verbose=False)
    with pytest.raises(ConfigurationError):
        qu.quast_map_reads([str(contigs)], str(tmp_path / "o3"), verbose=False)
    with pytest.raises(ConfigurationError):
        qu.metaquast_assess_metagenome([str(contigs)], str(tmp_path / "o4"), use_blast_search=False, verbose=False)


@needs_quast
def test_reference_metrics_and_kmers(tmp_path: Path) -> None:
    reference, contigs = genome_files(tmp_path)
    result = qu.quast_assess_assembly([str(contigs)], str(tmp_path / "q"), reference=str(reference), k_mer_stats=True, k_mer_size=31, verbose=False)
    assert result["report"]["Genome fraction (%)"]["contigs"] == pytest.approx(97.5)
    assert result["report"]["# misassemblies"]["contigs"] == 0
    assert result["report"]["N50"]["contigs"] == 10500 and "K-mer-based compl. (%)" in result["report"]


@needs_quast
def test_read_mapping(tmp_path: Path) -> None:
    reference, contigs = genome_files(tmp_path)
    genome = reference.read_text().split("\n")[1]
    reads = tmp_path / "reads.fq"
    with reads.open("w") as handle:
        for i in range(2000):
            start = random.randrange(0, len(genome) - 100)
            handle.write(f"@r{i}\n{genome[start:start + 100]}\n+\n{'I' * 100}\n")
    result = qu.quast_map_reads([str(contigs)], str(tmp_path / "q"), reference=str(reference), single=str(reads), verbose=False)
    assert result["report"]["Mapped (%)"]["contigs"] > 90
