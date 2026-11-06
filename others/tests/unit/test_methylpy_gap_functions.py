"""Tests for dnakit.epigenome.methylpy_gap_functions (methylpy facades).

Upstream methylpy functions are replaced by fakes in these tests (no bowtie2, samtools or network); the executables
are checked by the facade, so fake executables are created in tmp_path. Real runs need cutadapt (present in the venv)
or the external program (skipped with the reason when it is missing).
"""

from __future__ import annotations

import inspect
import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("methylpy")

from dnakit.epigenome import methylpy_gap_functions as facade
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError

_LOCAL_CUTADAPT = Path(sys.executable).parent / "cutadapt"
_CUTADAPT = shutil.which("cutadapt") or (
    str(_LOCAL_CUTADAPT) if _LOCAL_CUTADAPT.is_file() else None
)
needs_cutadapt = pytest.mark.skipif(
    _CUTADAPT is None, reason="cutadapt executable is not installed"
)
needs_bgzip_tabix = pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip/tabix (htslib) is not installed",
)
needs_samtools = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools is not installed"
)


def _gsl_installed() -> bool:
    import ctypes.util

    return ctypes.util.find_library("gsl") is not None


def _touch(directory: Path, *names: str) -> list[str]:
    paths = []
    for name in names:
        path = directory / name
        path.write_text("x\n", encoding="utf-8")
        paths.append(str(path))
    return paths


def _fake_tools(directory: Path, *names: str) -> str:
    for name in names:
        tool = directory / name
        tool.write_text("#!/bin/sh\n", encoding="utf-8")
        tool.chmod(0o755)
    return str(directory) + "/"


def _recorder(calls: list[dict]):
    def fake(*args, **kwargs):
        calls.append({"args": args, "kwargs": kwargs})

    return fake


def test_public_functions_take_keyword_only_verbose() -> None:
    assert len(facade.__all__) == len(set(facade.__all__))
    for name in facade.__all__:
        verbose = inspect.signature(getattr(facade, name)).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True


def test_missing_program_raises_before_running_upstream(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import methylpy.call_mc_se as upstream

    called: list[dict] = []
    monkeypatch.setattr(upstream, "run_mapping", _recorder(called))
    refs = _touch(tmp_path, "fwd.fa", "rev.fa", "ref.fa")
    reads = _touch(tmp_path, "reads.fq")
    with pytest.raises(BackendUnavailableError) as info:
        facade.methylpy_se_run_mapping(
            "lib1",
            reads,
            "s1",
            refs[0],
            refs[1],
            refs[2],
            path_to_samtools=str(tmp_path / "missing"),
            path_to_aligner=str(tmp_path),
            path_to_cutadapt=str(tmp_path),
        )
    assert "samtools" in str(info.value)
    assert "conda install -c bioconda samtools" in (info.value.hint or "")
    assert called == []


def test_only_arguments_set_by_the_caller_are_forwarded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import methylpy.call_mc_se as upstream

    calls: list[dict] = []
    monkeypatch.setattr(upstream, "run_mapping", _recorder(calls))
    tools = _fake_tools(tmp_path, "samtools", "bowtie2", "cutadapt")
    refs = _touch(tmp_path, "fwd.fa", "rev.fa", "ref.fa")
    reads = _touch(tmp_path, "reads.fq")
    result = facade.methylpy_se_run_mapping(
        "lib1",
        reads,
        "s1",
        refs[0],
        refs[1],
        refs[2],
        path_to_samtools=tools,
        path_to_aligner=tools,
        path_to_cutadapt=tools,
        min_mapq=20,
    )
    assert result["status"] == "completed"
    assert result["backend"] == "methylpy"
    assert result["function"] == "methylpy.call_mc_se.run_mapping"
    kwargs = calls[0]["kwargs"]
    assert kwargs["current_library"] == "lib1"
    assert kwargs["library_files"] == reads
    assert kwargs["min_mapq"] == 20
    assert kwargs["path_to_samtools"] == tools
    assert "merge_by_max_mapq" not in kwargs  # upstream default kept
    assert "aligner" not in kwargs


def test_pipeline_checks_cutadapt_only_when_trimming(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import methylpy.call_mc_se as upstream

    calls: list[dict] = []
    monkeypatch.setattr(upstream, "run_methylation_pipeline", _recorder(calls))
    tools = _fake_tools(tmp_path, "samtools", "bowtie2")
    refs = _touch(tmp_path, "fwd.fa", "rev.fa", "ref.fa")
    reads = _touch(tmp_path, "reads.fq")
    with pytest.raises(BackendUnavailableError):
        facade.methylpy_se_run_methylation_pipeline(
            reads,
            "s1",
            refs[0],
            refs[1],
            refs[2],
            path_to_samtools=tools,
            path_to_aligner=tools,
            path_to_cutadapt=str(tmp_path / "no_cutadapt"),
        )
    result = facade.methylpy_se_run_methylation_pipeline(
        reads,
        "s1",
        refs[0],
        refs[1],
        refs[2],
        path_to_samtools=tools,
        path_to_aligner=tools,
        trim_reads=False,
    )
    assert result["status"] == "completed"
    assert "trim_reads" in calls[0]["kwargs"]


def test_paired_end_mapping_forwards_both_read_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import methylpy.call_mc_pe as upstream

    calls: list[dict] = []
    monkeypatch.setattr(upstream, "run_mapping_pe", _recorder(calls))
    tools = _fake_tools(tmp_path, "samtools", "bowtie2", "cutadapt")
    refs = _touch(tmp_path, "fwd.fa", "rev.fa", "ref.fa")
    r1 = _touch(tmp_path, "r1.fq")
    r2 = _touch(tmp_path, "r2.fq")
    facade.methylpy_pe_run_mapping(
        "lib1",
        r1,
        r2,
        "s1",
        refs[0],
        refs[1],
        refs[2],
        path_to_samtools=tools,
        path_to_aligner=tools,
        path_to_cutadapt=tools,
    )
    kwargs = calls[0]["kwargs"]
    assert kwargs["library_read1_files"] == r1
    assert kwargs["library_read2_files"] == r2


def test_remove_clonal_needs_picard_jar(tmp_path: Path) -> None:
    bam = _touch(tmp_path, "in.bam")[0]
    with pytest.raises(ConfigurationError):
        facade.methylpy_remove_clonal_bam(
            bam, str(tmp_path / "out.bam"), "picard", path_to_picard=None
        )
    with pytest.raises(BackendUnavailableError) as info:
        facade.methylpy_remove_clonal_bam(
            bam, str(tmp_path / "out.bam"), "picard", path_to_picard=str(tmp_path)
        )
    assert "picard.jar" in str(info.value)


def test_bgzip_and_convert_report_missing_htslib_or_samtools(tmp_path: Path) -> None:
    allc = _touch(tmp_path, "allc.tsv")[0]
    with pytest.raises(BackendUnavailableError):
        facade.methylpy_bgzip_allc_file(
            allc, path_to_bgzip=str(tmp_path / "x"), path_to_tabix=str(tmp_path / "y")
        )
    refs = _touch(tmp_path, "ref.fa")[0]
    with pytest.raises(BackendUnavailableError):
        facade.methylpy_convert_allc_to_bigwig(
            allc,
            str(tmp_path / "out.bw"),
            refs,
            path_to_samtools=str(tmp_path / "x"),
            path_to_wigToBigWig=str(tmp_path / "y"),
        )


def test_inputs_are_checked_before_running(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError):
        facade.methylpy_dmrfind([str(tmp_path / "missing.allc")], str(tmp_path / "dmr"))
    with pytest.raises(ConfigurationError):
        facade.methylpy_se_find_multi_mappers(
            str(tmp_path / "missing.sam"), str(tmp_path / "out.sam")
        )


@pytest.mark.skipif(
    _gsl_installed(), reason="GSL is installed; the missing-library error path does not apply"
)
def test_bundled_rms_binary_reports_missing_gsl(tmp_path: Path) -> None:
    files = _touch(tmp_path, "a.allc", "b.allc")
    with pytest.raises(BackendExecutionError) as info:
        facade.methylpy_run_rms_tests(files, str(tmp_path / "rms_out"), ["a", "b"])
    assert "GSL" in str(info.value)


@needs_cutadapt
def test_quality_trim_se_uses_cutadapt_without_shell(tmp_path: Path) -> None:
    source = tmp_path / "reads.fq"
    source.write_text(
        "@r1\nACGTACGTAGATCGGAAGAGCACACGTCTG\n+\n" + "I" * 30 + "\n", encoding="utf-8"
    )
    output = tmp_path / "trimmed.fq"
    result = facade.methylpy_quality_trim_se(
        [str(source)], [str(output)], min_read_len=5, path_to_cutadapt=str(Path(_CUTADAPT).parent)
    )
    assert result["status"] == "completed"
    assert output.read_text(encoding="utf-8").splitlines()[1] == "ACGTACGT"


def test_quality_trim_se_checks_inputs_and_output_count(tmp_path: Path) -> None:
    source = _touch(tmp_path, "reads.fq")[0]
    with pytest.raises(ConfigurationError):
        facade.methylpy_quality_trim_se([source], [], path_to_cutadapt=str(tmp_path))
    with pytest.raises(ConfigurationError):
        facade.methylpy_quality_trim_se(
            [source],
            [str(tmp_path / "o.fq")],
            adapter_seq="ACGTXYZ!",
            path_to_cutadapt=str(tmp_path),
        )


@needs_cutadapt
def test_quality_trim_pe_runs_upstream_with_cutadapt(tmp_path: Path) -> None:
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    r1.write_text("@p\nACGTACGTAGATCGGAAGAGCACACGTCTGAAC\n+\n" + "I" * 33 + "\n", encoding="utf-8")
    r2.write_text("@p\nGGGGCCCCAAAATTTTGGGG\n+\n" + "I" * 20 + "\n", encoding="utf-8")
    out1 = tmp_path / "o1.fq"
    out2 = tmp_path / "o2.fq"
    result = facade.methylpy_quality_trim_pe(
        str(r1),
        str(out1),
        str(r2),
        str(out2),
        path_to_cutadapt=str(Path(_CUTADAPT).parent),
        min_read_len=5,
    )
    assert result["status"] == "completed", result
    assert out1.is_file() and out2.is_file()


@needs_bgzip_tabix
def test_bgzip_allc_file_writes_indexed_archive(tmp_path: Path) -> None:
    allc = tmp_path / "allc_chr1.tsv"
    allc.write_text("chr1\t100\t+\tCGG\t1\t3\t1\nchr1\t200\t-\tCHH\t0\t2\t0\n", encoding="utf-8")
    result = facade.methylpy_bgzip_allc_file(str(allc))
    assert result["status"] == "completed"
    assert (tmp_path / "allc_chr1.tsv.gz").is_file()
    assert (tmp_path / "allc_chr1.tsv.gz.tbi").is_file()


@needs_samtools
def test_merge_bam_files_with_samtools(tmp_path: Path) -> None:
    import subprocess

    sam = (
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
        "r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
    )
    inputs = []
    for name in ("a", "b"):
        sam_path = tmp_path / f"{name}.sam"
        sam_path.write_text(sam, encoding="utf-8")
        bam_path = tmp_path / f"{name}.bam"
        subprocess.run(["samtools", "view", "-b", "-o", str(bam_path), str(sam_path)], check=True)
        inputs.append(str(bam_path))
    output = tmp_path / "merged.bam"
    result = facade.methylpy_se_merge_bam_files(inputs, str(output))
    assert result["status"] == "completed"
    assert output.is_file()
