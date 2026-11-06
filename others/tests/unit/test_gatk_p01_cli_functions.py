"""Contract tests for the GATK command-line facades of pack P01 (no network; real runs skip without gatk)."""

from __future__ import annotations

import inspect
import shutil
from pathlib import Path

import pytest

import dnakit.alignment.gatk_p01_alignment_cli_functions as alignment
import dnakit.sequence.gatk_p01_reference_cli_functions as reference
import dnakit.variants.gatk_p01_variants_cli_functions as variants
from dnakit.bioinformatics import _gatk_p01_core as core
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

MODULES = (alignment, reference, variants)
requires_gatk = pytest.mark.skipif(
    shutil.which("gatk") is None, reason="the gatk launcher is not on PATH"
)
SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
    "r1\t0\tchr1\t10\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
)


def _all_functions():
    for module in MODULES:
        for name in module.__all__:
            yield name, getattr(module, name), module


def test_public_function_count_and_verbose_contract() -> None:
    names = [name for name, _, _ in _all_functions()]
    assert len(names) == len(set(names)) == 80
    for name, function, _ in _all_functions():
        parameter = inspect.signature(function).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is True, name


def test_every_spec_flag_is_listed_by_the_tool_help() -> None:
    for module in MODULES:
        for name, spec in module._SPECS.items():
            for param in spec.params:
                flags = (
                    [param.flag]
                    if isinstance(param.flag, str)
                    else list((param.flag or {}).items())
                )
                for entry in flags:
                    variant, flag = entry if isinstance(entry, tuple) else (None, entry)
                    tools = (
                        [spec.tools[variant]]
                        if variant
                        else [
                            tool
                            for key, tool in spec.tools.items()
                            if param.v is None or key in param.v
                        ]
                    )
                    assert all(flag in core.TOOL_OPTIONS[tool] for tool in tools), (
                        name,
                        param.name,
                    )


def test_javascript_execution_is_not_reachable() -> None:
    assert all("JAVASCRIPT_FILE" not in options for options in core.TOOL_OPTIONS.values())
    assert (
        "includeJavascript"
        not in inspect.signature(alignment.gatk_filter_sam_reads).parameters["filter"].annotation
    )


def test_missing_gatk_raises_with_install_hint(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "in.sam"
    source.write_text(SAM)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv("DNAKIT_GATK", raising=False)
    with pytest.raises(BackendUnavailableError) as raised:
        alignment.gatk_clean_sam(str(source), str(tmp_path / "out.sam"), verbose=False)
    assert raised.value.hint and "gatk" in raised.value.hint
    with pytest.raises(BackendUnavailableError):
        alignment.gatk_clean_sam(
            str(source), str(tmp_path / "out.sam"), gatk_path="/nonexistent/gatk", verbose=False
        )
    assert not (tmp_path / "out.sam").exists()


def test_validation_happens_before_execution(tmp_path: Path) -> None:
    source = tmp_path / "in.sam"
    source.write_text(SAM)
    existing = tmp_path / "exists.sam"
    existing.write_text("x")
    common = {"gatk_path": "/nonexistent/gatk", "verbose": False}
    with pytest.raises(ConfigurationError, match="not allowed"):
        alignment.gatk_clean_sam(
            str(source), str(tmp_path / "o.sam"), options={"not-an-option": 1}, **common
        )
    with pytest.raises(ConfigurationError, match="set by a parameter"):
        alignment.gatk_clean_sam(
            str(source), str(tmp_path / "o.sam"), options={"INPUT": "x"}, **common
        )
    with pytest.raises(ConfigurationError, match="will not be overwritten"):
        alignment.gatk_clean_sam(str(source), str(existing), **common)
    with pytest.raises(ConfigurationError, match="existing file"):
        alignment.gatk_clean_sam(str(tmp_path / "missing.sam"), str(tmp_path / "o.sam"), **common)
    with pytest.raises(ConfigurationError, match="one of"):
        alignment.gatk_filter_sam_reads(
            str(source), str(tmp_path / "o.sam"), "includeJavascript", **common
        )
    with pytest.raises(ConfigurationError, match="java_max_memory"):
        alignment.gatk_clean_sam(
            str(source),
            str(tmp_path / "o.sam"),
            java_max_memory="lots",
            gatk_path=__file__,
            verbose=False,
        )


def test_mode_dependent_arguments(tmp_path: Path) -> None:
    source = tmp_path / "in.bam"
    source.write_bytes(b"x")
    common = {"gatk_path": "/nonexistent/gatk", "verbose": False}
    with pytest.raises(ConfigurationError, match="required when the mode"):
        alignment.gatk_pathseq(step="score", **common)
    with pytest.raises(ConfigurationError, match="not used when the mode"):
        alignment.gatk_pathseq(
            step="filter", input_path=str(source), taxonomy_file=str(source), **common
        )
    with pytest.raises(ConfigurationError, match="at least one of"):
        alignment.gatk_pathseq(step="filter", input_path=str(source), **common)
    with pytest.raises(ConfigurationError, match="one of"):
        alignment.gatk_split_sam(str(source), str(tmp_path / "d"), by="nothing", **common)
    with pytest.raises(BackendUnavailableError):
        alignment.gatk_split_sam(str(source), str(tmp_path / "d"), by="library", **common)
    assert not (tmp_path / "d").exists()


def test_derived_index_must_not_exist(tmp_path: Path) -> None:
    bam = tmp_path / "a.bam"
    bam.write_bytes(b"x")
    (tmp_path / "a.bam.sbi").write_text("old")
    with pytest.raises(ConfigurationError, match="already exists"):
        alignment.gatk_create_bam_splitting_index(
            str(bam), gatk_path="/nonexistent/gatk", verbose=False
        )


@requires_gatk
def test_real_clean_sam_and_count_bases(tmp_path: Path) -> None:
    source = tmp_path / "in.sam"
    source.write_text(SAM)
    result = alignment.gatk_clean_sam(str(source), str(tmp_path / "out.sam"), verbose=False)
    assert result["status"] == "completed" and result["tool"] == "CleanSam"
    fasta = tmp_path / "r.fa"
    fasta.write_text(">c\nAACCGGTTN\n")
    (tmp_path / "r.fa.fai").write_text("c\t9\t3\t9\t10\n")
    (tmp_path / "r.dict").write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:9\n")
    counts = reference.gatk_count_reference_bases(str(fasta), verbose=False)
    assert {k: v for k, v in counts["base_counts"].items() if k in "ACGT"} == {
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
    }
