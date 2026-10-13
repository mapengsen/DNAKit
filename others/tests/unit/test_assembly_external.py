"""Contract tests for the Shasta assembly facade."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from dnakit.assembly import (
    assemble_flye,
    assemble_hifiasm,
    assemble_megahit,
    assemble_shasta,
    evaluate_merqury,
    run_busco,
)
from dnakit.exceptions import ConfigurationError
from dnakit.tools import default_tool_registry


def _fake_shasta(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "directory = Path(sys.argv[sys.argv.index('--assemblyDirectory') + 1])\n"
        "(directory / 'Assembly.fasta').write_text('>contig1\\nACGT\\n', encoding='ascii')\n"
        "(directory / 'Assembly.gfa').write_text('H\\tVN:Z:1.0\\n', encoding='ascii')\n"
        "print('shasta diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_assemble_shasta_writes_fasta_and_gfa(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = _fake_shasta(tmp_path / "fake-shasta")

    result = assemble_shasta(
        reads,
        tmp_path / "shasta-output",
        config="Nanopore-May2022",
        threads=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "shasta"
    assert "Assembly.fasta" in result["result"]["artifacts"]
    assert "Assembly.gfa" in result["result"]["artifacts"]


def test_assemble_shasta_validates_config(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    with pytest.raises(ConfigurationError, match="config"):
        assemble_shasta(reads, tmp_path / "invalid", config="bad config", verbose=False)


def test_evaluate_merqury_writes_qv_and_completeness(tmp_path: Path) -> None:
    read_database = tmp_path / "reads.meryl"
    read_database.mkdir()
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nACGT\n", encoding="ascii")
    executable = tmp_path / "fake-merqury.sh"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "prefix = Path(sys.argv[3])\n"
        "Path(str(prefix) + '.qv').write_text('QV 40\\n', encoding='ascii')\n"
        "Path(str(prefix) + '.completeness.stats').write_text('assembly\\t99\\n', encoding='ascii')\n"
        "print('merqury diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    executable.chmod(0o755)

    result = evaluate_merqury(
        read_database,
        assembly,
        tmp_path / "merqury-output",
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "merqury"
    assert "merqury.qv" in result["result"]["artifacts"]
    assert "merqury.completeness.stats" in result["result"]["artifacts"]


def test_run_busco_captures_completeness_artifacts(tmp_path: Path) -> None:
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nACGT\n", encoding="ascii")
    executable = tmp_path / "fake-busco"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "directory = Path(sys.argv[sys.argv.index('--out_path') + 1])\n"
        "run = directory / 'busco_run'\n"
        "run.mkdir()\n"
        "(run / 'short_summary.txt').write_text('C:100%\\n', encoding='ascii')\n"
        "print('busco diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    executable.chmod(0o755)

    result = run_busco(
        assembly,
        "eukaryota_odb10",
        tmp_path / "busco-output",
        threads=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "busco_run/short_summary.txt" in result["result"]["artifacts"]


def test_assemble_flye_writes_assembly_fasta(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = tmp_path / "fake-flye"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "directory = Path(sys.argv[sys.argv.index('--out-dir') + 1])\n"
        "(directory / 'assembly.fasta').write_text('>contig1\\nACGT\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    executable.chmod(0o755)
    result = assemble_flye(
        reads,
        tmp_path / "flye-output",
        genome_size="5m",
        platform="nano-hq",
        program=str(executable),
        verbose=False,
    )
    assert result["backend"] == "flye"
    assert "assembly.fasta" in result["result"]["artifacts"]


def test_assemble_megahit_writes_final_contigs(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = tmp_path / "fake-megahit"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "directory = Path(sys.argv[sys.argv.index('-o') + 1])\n"
        "directory.mkdir()\n"
        "(directory / 'final.contigs.fa').write_text('>contig1\\nACGT\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    executable.chmod(0o755)
    result = assemble_megahit(reads, tmp_path / "megahit-output", program=str(executable), verbose=False)
    assert result["backend"] == "megahit"
    assert "megahit/final.contigs.fa" in result["result"]["artifacts"]


def test_assemble_hifiasm_writes_primary_contig_graph(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="ascii")
    executable = tmp_path / "fake-hifiasm"
    executable.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "prefix = Path(sys.argv[sys.argv.index('-o') + 1])\n"
        "Path(str(prefix) + '.bp.p_ctg.gfa').write_text('H\\tVN:Z:1.0\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    executable.chmod(0o755)
    result = assemble_hifiasm(reads, tmp_path / "hifiasm-output", program=str(executable), verbose=False)
    assert result["backend"] == "hifiasm"
    assert "assembly.bp.p_ctg.gfa" in result["result"]["artifacts"]


def test_assemble_shasta_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.assembly.assemble_shasta")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True


def test_evaluate_merqury_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.assembly.evaluate_merqury")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True


def test_run_busco_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.assembly.run_busco")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True


def test_new_assembly_facades_are_registered_as_external_tools() -> None:
    registry = default_tool_registry(verbose=False)
    for name in (
        "dnakit.assembly.assemble_flye",
        "dnakit.assembly.assemble_hifiasm",
        "dnakit.assembly.assemble_megahit",
    ):
        spec = registry.resolve(name)
        assert spec.effect == "external_process"
        assert spec.agent_compatible is True
        assert spec.requires_confirmation is True
