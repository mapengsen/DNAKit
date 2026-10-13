"""Contract tests for the bounded SeqMagick facade."""

from __future__ import annotations

import sys
from pathlib import Path

from dnakit.fasta_fastq import cluster_vsearch, convert_seqmagick
from dnakit.tools import default_tool_registry


def _fake_seqmagick(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "arguments = sys.argv\n"
        "source = Path(arguments[-2])\n"
        "target = Path(arguments[-1])\n"
        "target.write_text(source.read_text(encoding='ascii'), encoding='ascii')\n"
        "print('seqmagick diagnostic', file=sys.stderr)\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def _fake_vsearch(path: Path) -> Path:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "Path(sys.argv[sys.argv.index('--centroids') + 1]).write_text('>centroid\\nACGT\\n', encoding='ascii')\n"
        "Path(sys.argv[sys.argv.index('--uc') + 1]).write_text('S\\t0\\t4\\t*\\t*\\t*\\t*\\t*\\tcentroid\\n', encoding='ascii')\n",
        encoding="ascii",
    )
    path.chmod(0o755)
    return path


def test_convert_seqmagick(tmp_path: Path) -> None:
    source = tmp_path / "source.fasta"
    source.write_text(">r1\nACGT\n", encoding="ascii")
    executable = _fake_seqmagick(tmp_path / "seqmagick")

    result = convert_seqmagick(
        source,
        tmp_path / "seqmagick-output",
        output_format="fasta",
        ungap=True,
        upper=True,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert result["backend"] == "seqmagick"
    assert "converted.fasta" in result["result"]["artifacts"]


def test_convert_seqmagick_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.fasta_fastq.convert_seqmagick")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True


def test_cluster_vsearch_writes_centroids_and_uc(tmp_path: Path) -> None:
    source = tmp_path / "source.fasta"
    source.write_text(">r1\nACGT\n>r2\nACGA\n", encoding="ascii")
    executable = _fake_vsearch(tmp_path / "vsearch")

    result = cluster_vsearch(
        source,
        tmp_path / "vsearch-output",
        identity=0.9,
        threads=2,
        program=str(executable),
        verbose=False,
    )

    assert result["status"] == "completed"
    assert "centroids.fasta" in result["result"]["artifacts"]
    assert "clusters.uc" in result["result"]["artifacts"]


def test_cluster_vsearch_is_registered_as_external_tool() -> None:
    spec = default_tool_registry(verbose=False).resolve("dnakit.fasta_fastq.cluster_vsearch")
    assert spec.effect == "external_process"
    assert spec.agent_compatible is True
    assert spec.requires_confirmation is True
