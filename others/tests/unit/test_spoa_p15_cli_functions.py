from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from dnakit.alignment import spoa_p15_cli_functions as sp
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

needs_spoa = pytest.mark.skipif(shutil.which("spoa") is None, reason="spoa is not installed")
FASTA = ">a\nACGTACGTACGTAAC\n>b\nACGTACCTACGTAAC\n>c\nACGTACGTACGAAC\n>d\nGTTACGTACGTACGT\n"


def test_missing_program_and_validation(tmp_path: Path) -> None:
    source = tmp_path / "s.fa"
    source.write_text(FASTA)
    with pytest.raises(BackendUnavailableError):
        sp.spoa_partial_order_alignment(
            str(source), str(tmp_path / "o.fa"), spoa_path=str(tmp_path / "missing"), verbose=False
        )
    with pytest.raises(ConfigurationError):
        sp.spoa_partial_order_alignment(
            str(source), str(tmp_path / "o.fa"), result="nothing", verbose=False
        )
    with pytest.raises(ConfigurationError):
        sp.spoa_partial_order_alignment(
            str(source), str(tmp_path / "o.fa"), gap_open=5, verbose=False
        )
    with pytest.raises(ConfigurationError):
        sp.spoa_partial_order_alignment(
            str(source), str(tmp_path / "o.fa"), result="msa", min_coverage=2, verbose=False
        )


@needs_spoa
def test_msa_graph_dot_and_strand(tmp_path: Path) -> None:
    source = tmp_path / "s.fa"
    source.write_text(FASTA)
    msa = sp.spoa_partial_order_alignment(
        str(source), str(tmp_path / "m.fa"), algorithm="global", verbose=False
    )
    assert len({len(r["sequence"]) for r in msa["alignment"]}) == 1 and len(msa["alignment"]) == 4
    both = sp.spoa_partial_order_alignment(
        str(source),
        str(tmp_path / "b.fa"),
        result="consensus_and_msa",
        algorithm="global",
        strand_ambiguous=True,
        dot_path=str(tmp_path / "g.dot"),
        verbose=False,
    )
    assert both["consensus"] == "ACGTACGTACGTAAC" and (tmp_path / "g.dot").read_text().startswith(
        "digraph"
    )
    assert (
        both["alignment"][3]["sequence"].replace("-", "") == "ACGTACGTACGTAAC"
    )  # read d was reverse-complemented
    graph = sp.spoa_partial_order_alignment(
        str(source),
        str(tmp_path / "g.gfa"),
        result="consensus_and_graph",
        algorithm="global",
        verbose=False,
    )
    assert graph["path_count"] == 5 and graph["node_count"] > 10
