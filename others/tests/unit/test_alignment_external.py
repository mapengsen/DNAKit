"""Contracts for the formal BLAST+, minimap2 and BEDTools facades."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.alignment import (
    build_blast_database,
    consensus_spoa,
    intersect_bed,
    local_blast,
    map_long_reads,
    sort_bed_bedops,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools import default_tool_registry


def _fake_program(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "import pathlib, sys\n"
        "name = pathlib.Path(sys.argv[0]).name\n"
        "if name == 'fake-minimap2':\n"
        "    print('@SQ\\tSN:chr1\\tLN:4')\n"
        "    print('read\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII')\n"
        "elif name == 'fake-bedtools':\n"
        "    print('chr1\\t0\\t4\\tchr1\\t2\\t6')\n"
        "elif name == 'fake-sort-bed':\n"
        "    print('chr1\\t0\\t4')\n"
        "elif name == 'fake-makeblastdb':\n"
        "    output = sys.argv[sys.argv.index('-out') + 1]\n"
        "    for suffix in ('.nhr', '.nin', '.nsq'):\n"
        "        pathlib.Path(output + suffix).write_text('index', encoding='ascii')\n"
        "elif name == 'fake-blastn':\n"
        "    output = sys.argv[sys.argv.index('-out') + 1]\n"
        "    pathlib.Path(output).write_text('<BlastOutput/>\\n', encoding='ascii')\n"
        "elif name == 'fake-spoa':\n"
        "    print('>consensus')\n"
        "    print('ACGT')\n"
        "print('diagnostic', file=sys.stderr)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def test_minimap2_facade_returns_bounded_json_result(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fq"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    reads.write_text("@read\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_program(tmp_path / "fake-minimap2")

    result = map_long_reads(
        reference,
        reads,
        tmp_path / "mapped",
        program=str(executable),
        verbose=False,
    )

    assert result["backend"] == "minimap2"
    artifact = result["result"]["artifacts"]["alignments.sam"]
    assert "@SQ" in Path(artifact).read_text(encoding="ascii")


def test_blast_facades_keep_database_and_search_outputs_separate(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    query = tmp_path / "query.fa"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    query.write_text(">query\nACGT\n", encoding="ascii")
    makeblastdb = _fake_program(tmp_path / "fake-makeblastdb")
    blastn = _fake_program(tmp_path / "fake-blastn")

    database = build_blast_database(
        reference,
        tmp_path / "database",
        program=str(makeblastdb),
        verbose=False,
    )
    result = local_blast(
        query,
        tmp_path / "database" / "database",
        tmp_path / "search",
        program=str(blastn),
        verbose=False,
    )

    assert database["operation"] == "build_blast_database"
    assert result["result"]["artifacts"]["results.xml"].endswith("results.xml")


def test_bedtools_facade_uses_fixed_intersection_operation(tmp_path: Path) -> None:
    first = tmp_path / "first.bed"
    second = tmp_path / "second.bed"
    first.write_text("chr1\t0\t4\n", encoding="ascii")
    second.write_text("chr1\t2\t6\n", encoding="ascii")
    executable = _fake_program(tmp_path / "fake-bedtools")

    result = intersect_bed(first, second, tmp_path / "intersection", program=str(executable))

    assert result["result"]["operation"] == "intersect -wa -wb"
    assert Path(result["result"]["artifacts"]["intersection.bed"]).read_text().startswith("chr1")


def test_spoa_facade_captures_consensus_fasta(tmp_path: Path) -> None:
    source = tmp_path / "reads.fasta"
    source.write_text(">one\nACGT\n>two\nACGA\n", encoding="ascii")
    executable = _fake_program(tmp_path / "fake-spoa")

    result = consensus_spoa(
        source,
        tmp_path / "consensus",
        algorithm=1,
        strand_ambiguous=True,
        program=str(executable),
        verbose=False,
    )

    assert result["backend"] == "spoa"
    assert ">consensus" in Path(result["result"]["artifacts"]["consensus.fasta"]).read_text()


def test_bedops_facade_captures_sorted_bed(tmp_path: Path) -> None:
    source = tmp_path / "source.bed"
    source.write_text("chr1\t0\t4\n", encoding="ascii")
    executable = _fake_program(tmp_path / "fake-sort-bed")

    result = sort_bed_bedops(source, tmp_path / "sorted", program=str(executable), verbose=False)

    assert result["backend"] == "bedops"
    assert Path(result["result"]["artifacts"]["sorted.bed"]).read_text().startswith("chr1")


def test_external_alignment_facades_are_registered_and_validate_fresh_outputs() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.alignment.align_sequences",
        "dnakit.alignment.consensus_spoa",
        "dnakit.alignment.build_blast_database",
        "dnakit.alignment.local_blast",
        "dnakit.alignment.map_long_reads",
        "dnakit.alignment.intersect_bed",
        "dnakit.alignment.sort_bed_bedops",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True


def test_msa_facade_rejects_single_sequence_before_external_execution(tmp_path: Path) -> None:
    from dnakit.alignment import align_sequences

    with pytest.raises(ConfigurationError, match="two sequences"):
        align_sequences({"only": "ACGT"}, tmp_path / "msa", verbose=False)


def test_minimap2_facade_rejects_existing_output_directory(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reads = tmp_path / "reads.fq"
    reference.write_text(">chr1\nACGT\n", encoding="ascii")
    reads.write_text("@read\nACGT\n+\nIIII\n", encoding="ascii")
    output = tmp_path / "mapped"
    output.mkdir()
    with pytest.raises(ConfigurationError, match="already exists"):
        map_long_reads(reference, reads, output, verbose=False)
