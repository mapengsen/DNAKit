"""Tests for dnakit.sequencing_qc.umi_tools_gap_wrappers (umi_tools sub-command facades)."""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import pytest

pytest.importorskip("umi_tools")
pysam = pytest.importorskip("pysam")

from dnakit.exceptions import BackendUnavailableError, ConfigurationError  # noqa: E402
from dnakit.sequencing_qc import umi_tools_gap_wrappers as wrappers  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "paper/write_paper/figure/fig2/task_coverage/consistency/data"
READS_BAM = DATA / "reads.bam"
READS_FASTQ = DATA / "reads.fastq"
_VENV_BIN = Path(sys.executable).parent
_UMIS = ("ACGTAC", "ACGTAG", "TTTTGG", "GGCCAA")


@pytest.fixture(autouse=True)
def _umi_tools_on_path(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", f"{_VENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}")


@pytest.fixture
def tagged_bam(tmp_path: Path) -> Path:
    """First 600 mapped reads of the shared reads.bam, with a UMI in RX and a gene in XT (built in tmp_path).

    The shared BAM has no UMI or gene tags, so the tags are added here; the shared file itself is not written.
    """
    if not READS_BAM.is_file():
        pytest.skip("shared reads.bam is missing")
    source = pysam.AlignmentFile(str(READS_BAM), "rb")
    unsorted = tmp_path / "unsorted.bam"
    reads = []
    for index, read in enumerate(source.fetch(until_eof=True)):
        if len(reads) >= 600:
            break
        if read.is_unmapped:
            continue
        read.set_tag("RX", _UMIS[index % len(_UMIS)], value_type="Z")
        read.set_tag("XT", f"gene{read.reference_start // 2000}", value_type="Z")
        reads.append(read)
    with pysam.AlignmentFile(str(unsorted), "wb", header=source.header) as out:
        for read in reads:
            out.write(read)
    source.close()
    sorted_bam = tmp_path / "tagged.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return sorted_bam


def _cell_fastq(tmp_path: Path, name: str, cell_length: int = 8, umi_length: int = 4) -> tuple[Path, set[str]]:
    """Synthetic barcode reads: 30 real cells with 40-120 reads each, plus 300 barcodes with 1-2 reads.

    Read sequences come from the shared reads.fastq. Returns the FASTQ path and the set of real cell barcodes.
    """
    rng = random.Random(7)
    lines = READS_FASTQ.read_text().split("\n")
    sequences = [lines[i + 1] for i in range(0, len(lines) - 3, 4) if lines[i + 1]]

    def random_bases(length: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(length))

    true_cells = [random_bases(cell_length) for _ in range(30)]
    plan: list[str] = []
    for cell in true_cells:
        plan += [cell] * rng.randint(40, 120)
    for _ in range(300):
        plan += [random_bases(cell_length)] * rng.randint(1, 2)
    rng.shuffle(plan)
    records = []
    for index, cell in enumerate(plan):
        sequence = sequences[index % len(sequences)][:60]
        umi = random_bases(umi_length)
        records.append(f"@read{index:05d}\n{cell}{umi}{sequence}\n+\n{'I' * (cell_length + umi_length + len(sequence))}\n")
    path = tmp_path / name
    path.write_text("".join(records))
    return path, set(true_cells)


@pytest.fixture
def barcode_fastq(tmp_path: Path) -> tuple[Path, set[str]]:
    return _cell_fastq(tmp_path, "barcodes.fastq")


def test_count_per_gene_counts_umi_groups(tagged_bam: Path, tmp_path: Path) -> None:
    out = tmp_path / "gene_counts.tsv"
    result = wrappers.umi_count_per_gene(str(tagged_bam), str(out), gene_tag="XT", extract_umi_method="tag",
                                         umi_tag="RX", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    lines = out.read_text().splitlines()
    assert lines[0] == "gene\tcount"
    counts = {line.split("\t")[0]: int(line.split("\t")[1]) for line in lines[1:]}
    assert counts and all(1 <= value <= len(_UMIS) for value in counts.values())
    assert result["version"] == "1.1.6"


def test_count_needs_a_gene_definition(tagged_bam: Path, tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        wrappers.umi_count_per_gene(str(tagged_bam), str(tmp_path / "x.tsv"), verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.umi_count_per_gene(str(tagged_bam), str(tmp_path / "y.tsv"), gene_tag="XT", per_contig=True,
                                    verbose=False)


def test_count_table_counts_umis_per_gene(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    table.write_text("readA_ACGT\tgene1\nreadB_TTGA\tgene1\nreadC_ACGT\tgene1\nreadD_GGCA\tgene2\n")
    out = tmp_path / "counted.tsv"
    result = wrappers.umi_count_table(str(table), str(out), verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert out.read_text().splitlines() == ["gene\tcount", "gene1\t2", "gene2\t1"]


def test_dedup_removes_duplicates(tagged_bam: Path, tmp_path: Path) -> None:
    out = tmp_path / "dedup.bam"
    result = wrappers.umi_dedup_bam(str(tagged_bam), str(out), extract_umi_method="tag", umi_tag="RX", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    before = sum(1 for _ in pysam.AlignmentFile(str(tagged_bam)))
    after = sum(1 for _ in pysam.AlignmentFile(str(out)))
    assert 0 < after < before


def test_group_tags_reads_and_writes_the_table(tagged_bam: Path, tmp_path: Path) -> None:
    out = tmp_path / "grouped.bam"
    table = tmp_path / "grouped.tsv"
    result = wrappers.umi_group_bam(str(tagged_bam), str(out), group_tsv=str(table), extract_umi_method="tag",
                                    umi_tag="RX", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    assert result["outputs"] == [str(out), str(table)]
    header = table.read_text().splitlines()[0].split("\t")
    assert header[:5] == ["read_id", "contig", "position", "gene", "umi"]
    tagged = [read for read in pysam.AlignmentFile(str(out)) if read.has_tag("BX")]
    assert tagged, "reads should carry the BX group tag"


def test_extract_moves_the_barcode_into_read_names(barcode_fastq: tuple[Path, set[str]], tmp_path: Path) -> None:
    source, _cells = barcode_fastq
    out = tmp_path / "extracted.fastq"
    result = wrappers.umi_extract_to_names(str(source), str(out), bc_pattern="CCCCCCCCNNNN", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    first_header = out.read_text().splitlines()[0]
    raw_sequence = source.read_text().splitlines()[1]
    expected_header = f"@read00000_{raw_sequence[:8]}_{raw_sequence[8:12]}"
    assert first_header == expected_header


def test_extract_refuses_bad_arguments(tmp_path: Path, barcode_fastq: tuple[Path, set[str]]) -> None:
    source, _cells = barcode_fastq
    with pytest.raises(ConfigurationError):
        wrappers.umi_extract_to_names(str(source), str(tmp_path / "a.fq"), bc_pattern="CCCCCCCC", verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.umi_extract_to_names(str(source), str(tmp_path / "b.fq"), bc_pattern="CCCCNNNN",
                                      read2_out=str(tmp_path / "r2.fq"), verbose=False)


def test_whitelist_recovers_the_real_cells(barcode_fastq: tuple[Path, set[str]], tmp_path: Path) -> None:
    source, true_cells = barcode_fastq
    out = tmp_path / "whitelist.tsv"
    result = wrappers.umi_whitelist_cells(str(source), str(out), bc_pattern="CCCCCCCCNNNN", verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    barcodes = {line.split("\t")[0] for line in out.read_text().splitlines() if line}
    assert barcodes == true_cells


def test_whitelist_refuses_a_pattern_without_cell_bases(barcode_fastq: tuple[Path, set[str]], tmp_path: Path) -> None:
    source, _cells = barcode_fastq
    with pytest.raises(ConfigurationError):
        wrappers.umi_whitelist_cells(str(source), str(tmp_path / "w.tsv"), bc_pattern="NNNNNNNN", verbose=False)
    with pytest.raises(ConfigurationError):
        wrappers.umi_whitelist_cells(str(source), str(tmp_path / "w2.tsv"), bc_pattern="CCCCCCCCNNNN",
                                     knee_method="sideways", verbose=False)


def test_prepare_for_rsem_writes_a_bam(tagged_bam: Path, tmp_path: Path) -> None:
    out = tmp_path / "prepared.bam"
    result = wrappers.umi_prepare_bam_for_rsem(str(tagged_bam), str(out), verbose=False)
    assert result["status"] == "completed", result["log_tail"]
    with pysam.AlignmentFile(str(out)) as handle:
        assert handle.header is not None


def test_missing_program_gives_install_hint(tmp_path: Path, tagged_bam: Path) -> None:
    with pytest.raises(BackendUnavailableError) as info:
        wrappers.umi_count_per_gene(str(tagged_bam), str(tmp_path / "m.tsv"), gene_tag="XT",
                                    umi_tools_path=str(tmp_path / "no-umi-tools"), verbose=False)
    assert "pip install umi_tools" in info.value.hint


def test_public_names_are_exported() -> None:
    assert set(wrappers.__all__) == {
        "umi_count_per_gene", "umi_count_table", "umi_dedup_bam", "umi_extract_to_names",
        "umi_group_bam", "umi_prepare_bam_for_rsem", "umi_whitelist_cells",
    }
