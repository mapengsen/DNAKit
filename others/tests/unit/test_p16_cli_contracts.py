"""Contract tests of the P16 command-line facades that need no external program and no network."""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

from dnakit.exceptions import BackendUnavailableError, ConfigurationError

MODULES = {
    "dnakit.microbiome.p16_microbial_cli_functions": 23,
    "dnakit.sequencing_qc.p16_readqc_cli_functions": 5,
    "dnakit.variants.p16_variant_cli_functions": 7,
    "dnakit.variants.freebayes_p16_cli_functions": 5,
    "dnakit.annotation.p16_assembly_annotation_cli_functions": 5,
    "dnakit.epigenome.p16_methyldackel_cli_functions": 5,
    "dnakit.phylogeny.p16_fasttree_cli_functions": 2,
    "dnakit.design_crispr.p16_crispr_screen_cli_functions": 7,
    "dnakit.representations.lucaonetasks_p16_functions": 1,
}


@pytest.mark.parametrize(("name", "count"), MODULES.items())
def test_public_functions_have_keyword_only_verbose(name: str, count: int) -> None:
    module = importlib.import_module(name)
    assert len(module.__all__) == count
    for function_name in module.__all__:
        function = getattr(module, function_name)
        assert function.__doc__, function_name
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True, (
            function_name
        )


def test_missing_program_raises_with_install_hint(tmp_path: Path, monkeypatch) -> None:
    from dnakit.phylogeny import p16_fasttree_cli_functions as ft

    alignment = tmp_path / "a.fa"
    alignment.write_text(">a\nACGT\n>b\nACGA\n")
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(BackendUnavailableError) as raised:
        ft.compute_fasttree_distance_matrix(
            str(alignment), str(tmp_path / "d.phy"), fasttree_path=None, verbose=False
        )
    assert raised.value.hint and "fasttree" in raised.value.hint.lower()


def test_explicit_missing_program_path_raises(tmp_path: Path) -> None:
    from dnakit.sequencing_qc import p16_readqc_cli_functions as rq

    reads = tmp_path / "r.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    with pytest.raises(BackendUnavailableError):
        rq.filter_long_reads_with_reference(
            str(reads), str(tmp_path / "o.fq"), filtlong_path=str(tmp_path / "nope"), verbose=False
        )


def test_outputs_are_never_overwritten(tmp_path: Path) -> None:
    from dnakit.phylogeny import p16_fasttree_cli_functions as ft

    alignment = tmp_path / "a.fa"
    alignment.write_text(">a\nACGT\n>b\nACGA\n")
    existing = tmp_path / "d.phy"
    existing.write_text("x")
    with pytest.raises(ConfigurationError, match="must not exist"):
        ft.compute_fasttree_distance_matrix(str(alignment), str(existing), verbose=False)


def test_argument_validation(tmp_path: Path) -> None:
    from dnakit.microbiome import p16_microbial_cli_functions as mb
    from dnakit.sequencing_qc import p16_readqc_cli_functions as rq
    from dnakit.variants import p16_variant_cli_functions as pv

    fasta = tmp_path / "x.fa"
    fasta.write_text(">a\nACGT\n")
    with pytest.raises(ConfigurationError, match="rank"):
        mb.summarize_kaiju_table(
            [str(fasta)], str(fasta), str(fasta), str(tmp_path / "t"), rank="kingdom", verbose=False
        )
    with pytest.raises(ConfigurationError, match="step"):
        mb.run_genomad(str(fasta), str(tmp_path / "g"), step="bogus", verbose=False)
    with pytest.raises(ConfigurationError, match="database"):
        mb.run_genomad(str(fasta), str(tmp_path / "g"), step="annotate", verbose=False)
    with pytest.raises(ConfigurationError, match="reference"):
        rq.filter_long_reads_with_reference(
            str(fasta), str(tmp_path / "o"), trim=True, verbose=False
        )
    with pytest.raises(ConfigurationError, match="Only one"):
        mb.summarize_kaiju_table(
            [str(fasta)],
            str(fasta),
            str(fasta),
            str(tmp_path / "t2"),
            min_percent=1.0,
            min_reads=2,
            verbose=False,
        )
    with pytest.raises(ConfigurationError, match="colon-separated"):
        pv.compute_mosdepth_quantized_coverage(
            str(fasta), str(tmp_path / "q"), quantize="a:b", verbose=False
        )
