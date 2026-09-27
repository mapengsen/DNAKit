"""CLI-adapter contract tests using real processes with tiny synthetic backends.

These executables exercise argv, output handling and error propagation only.
They do not run scientific tools and do not validate scientific correctness.
Official CLI references were reviewed separately for PLINK 2, GCTA, Modkit,
and FastANI; no network access or optional scientific software is needed here.
"""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import (
    BackendExecutionError,
    BackendTimeoutError,
    BackendUnavailableError,
    ConfigurationError,
)
from dnakit.integrations.catalog import describe_analysis, list_analyses
from dnakit.integrations.execution import plan_analysis, run_analysis


def _executable(tmp_path: Path, body: str, *, name: str = "tiny backend") -> Path:
    """Create a local Python executable; subprocess must never interpret a shell."""
    program = tmp_path / name
    program.write_text(
        f"#!{sys.executable}\n"
        "import json, sys, time\n"
        "from pathlib import Path\n"
        "args = sys.argv[1:]\n" + textwrap.dedent(body),
        encoding="utf-8",
    )
    program.chmod(0o755)
    return program


@pytest.fixture
def plink_inputs(tmp_path: Path) -> dict[str, str]:
    prefix = tmp_path / "genotype set"
    for suffix in (".bed", ".bim", ".fam"):
        Path(str(prefix) + suffix).write_bytes(b"contract fixture\n")
    return {"genotypes": str(prefix)}


@pytest.fixture
def fasta_inputs(tmp_path: Path) -> dict[str, str]:
    query = tmp_path / "query with spaces.fa"
    reference = tmp_path / "reference.fa"
    query.write_text(">query\nACGTACGT\n", encoding="utf-8")
    reference.write_text(">reference\nACGTACGT\n", encoding="utf-8")
    return {"query": str(query), "reference": str(reference)}


def _pca_backend(tmp_path: Path, *, empty: bool = False) -> Path:
    return _executable(
        tmp_path,
        f"""
        prefix = Path(args[args.index('--out') + 1])
        content = {"" if empty else "contract-output"!r}
        Path(str(prefix) + '.eigenvec').write_text(content)
        Path(str(prefix) + '.eigenval').write_text(content)
        Path('invocation.json').write_text(json.dumps(args))
        print('contract backend completed')
        """,
    )


def test_planning_is_read_only_and_preserves_spaced_argv(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    directory = tmp_path / "planned output"
    plan = plan_analysis(
        "plink2-pca",
        plink_inputs,
        str(directory),
        parameters={"components": 3},
        threads=2,
        executable_paths={"plink2": str(tmp_path / "not-installed")},
        verbose=False,
    )
    assert not directory.exists()
    assert plan["status"] == "planned"
    step = plan["steps"][0]
    assert step["executable"] is None
    assert step["arguments"] == [
        "--bfile",
        plink_inputs["genotypes"],
        "--pca",
        "3",
        "--out",
        str(directory / "pca"),
        "--threads",
        "2",
    ]


def test_missing_backend_never_creates_output_directory(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    directory = tmp_path / "missing backend output"
    with pytest.raises(BackendUnavailableError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            executable_paths={"plink2": str(tmp_path / "does-not-exist")},
            verbose=False,
        )
    assert not directory.exists()


@pytest.mark.parametrize("suffix", [".bed", ".bim", ".fam"])
def test_plink_prefix_requires_each_component_before_writes(
    tmp_path: Path, plink_inputs: dict[str, str], suffix: str
) -> None:
    Path(plink_inputs["genotypes"] + suffix).unlink()
    directory = tmp_path / "incomplete prefix"
    with pytest.raises(ConfigurationError, match="Incomplete"):
        plan_analysis("plink2-pca", plink_inputs, str(directory), verbose=False)
    assert not directory.exists()


def test_gcta_prefix_and_quantitative_covariate_contract(tmp_path: Path) -> None:
    prefix = tmp_path / "relationship matrix"
    for suffix in (".grm.bin", ".grm.N.bin", ".grm.id"):
        Path(str(prefix) + suffix).write_bytes(b"contract fixture\n")
    phenotype = tmp_path / "phenotypes.txt"
    covariates = tmp_path / "covariates.txt"
    phenotype.write_text("F I 1.5\n", encoding="utf-8")
    covariates.write_text("F I 0.2\n", encoding="utf-8")
    inputs = {"grm": str(prefix), "phenotypes": str(phenotype), "covariates": str(covariates)}
    plan = plan_analysis("gcta-reml", inputs, str(tmp_path / "reml"), verbose=False)
    arguments = plan["steps"][0]["arguments"]
    assert arguments[arguments.index("--grm") + 1] == str(prefix)
    assert arguments[arguments.index("--qcovar") + 1] == str(covariates)
    assert "--reml" in arguments
    assert plan["steps"][0]["expected_outputs"] == ["reml.hsq"]
    Path(str(prefix) + ".grm.N.bin").unlink()
    with pytest.raises(ConfigurationError, match="Incomplete"):
        plan_analysis("gcta-reml", inputs, str(tmp_path / "reml"), verbose=False)


@pytest.mark.parametrize("value", [True, 0, 101, 1.5, "2", None])
def test_integer_parameters_are_strict_and_bounded(
    tmp_path: Path, plink_inputs: dict[str, str], value: Any
) -> None:
    directory = tmp_path / "invalid parameter"
    with pytest.raises(ConfigurationError):
        plan_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            parameters={"components": value},
            verbose=False,
        )
    assert not directory.exists()


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), -0.1, 1.1, "0.2"])
def test_numeric_parameters_reject_nonfinite_or_out_of_range(
    tmp_path: Path, plink_inputs: dict[str, str], value: Any
) -> None:
    with pytest.raises(ConfigurationError):
        plan_analysis(
            "plink2-ld-prune",
            plink_inputs,
            str(tmp_path / "invalid r2"),
            parameters={"r2": value},
            verbose=False,
        )


@pytest.mark.parametrize("threads", [True, 0, 1025, 1.2, "2"])
def test_threads_validation_happens_before_creating_output(
    tmp_path: Path, plink_inputs: dict[str, str], threads: Any
) -> None:
    directory = tmp_path / "bad threads"
    with pytest.raises(ConfigurationError):
        run_analysis("plink2-pca", plink_inputs, str(directory), threads=threads, verbose=False)
    assert not directory.exists()


def test_unknown_arguments_and_unrelated_executable_are_rejected(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    directory = str(tmp_path / "output")
    with pytest.raises(ConfigurationError):
        plan_analysis("plink2-pca", plink_inputs, directory, parameters={"shell": True})
    with pytest.raises(ConfigurationError):
        plan_analysis("plink2-pca", plink_inputs, directory, executable_paths={"sh": "/bin/sh"})
    with pytest.raises(ConfigurationError):
        plan_analysis("unreviewed-command", plink_inputs, directory)
    with pytest.raises(ConfigurationError):
        plan_analysis("plink2-pca", {**plink_inputs, "extra": "file"}, directory)


def test_kraken_database_requires_its_declared_file_set(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reads.write_text(">read\nACGT\n", encoding="utf-8")
    database = tmp_path / "database"
    database.mkdir()
    (database / "hash.k2d").write_bytes(b"contract fixture")
    with pytest.raises(ConfigurationError, match="missing"):
        plan_analysis(
            "kraken2-classify",
            {"reads": str(reads), "database": str(database)},
            str(tmp_path / "output"),
            verbose=False,
        )


@pytest.mark.parametrize("kingdom", ["--help", "bac\n--help", "unknown", ""])
def test_string_parameters_do_not_accept_options_or_invalid_enums(
    tmp_path: Path, fasta_inputs: dict[str, str], kingdom: str
) -> None:
    with pytest.raises(ConfigurationError):
        plan_analysis(
            "barrnap-rrna",
            {"assembly": fasta_inputs["query"]},
            str(tmp_path / "output"),
            parameters={"kingdom": kingdom},
            verbose=False,
        )


def test_existing_directory_is_protected_before_process_start(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    directory = tmp_path / "existing"
    directory.mkdir()
    sentinel = directory / "keep.txt"
    sentinel.write_text("original", encoding="utf-8")
    backend = _pca_backend(tmp_path)
    with pytest.raises(ConfigurationError, match="already exists"):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    assert sentinel.read_text() == "original"
    assert sorted(path.name for path in directory.iterdir()) == ["keep.txt"]


def test_real_process_uses_argv_without_shell_and_writes_success_manifest(
    tmp_path: Path, plink_inputs: dict[str, str], monkeypatch: pytest.MonkeyPatch
) -> None:
    backend = _pca_backend(tmp_path)
    # If interpolated through a shell, these names would perform substitution.
    directory = tmp_path / "result $(touch injected) `touch injected2`; output"
    original_popen = subprocess.Popen
    calls: list[tuple[Any, dict[str, Any]]] = []

    def recording_popen(command: Any, **kwargs: Any) -> Any:
        calls.append((command, kwargs))
        return original_popen(command, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", recording_popen)
    result = run_analysis(
        "plink2-pca",
        plink_inputs,
        str(directory),
        parameters={"components": 2},
        executable_paths={"plink2": str(backend)},
        verbose=False,
    )
    assert len(calls) == 1 and calls[0][1]["shell"] is False
    assert isinstance(calls[0][0], (list, tuple))
    assert calls[0][0][0] == str(backend)
    arguments = json.loads((directory / "invocation.json").read_text())
    assert arguments == result["steps"][0]["arguments"]
    assert arguments[arguments.index("--bfile") + 1] == plink_inputs["genotypes"]
    assert result["status"] == "completed"
    assert {"pca.eigenvec", "pca.eigenval"} <= set(result["artifacts"])
    assert "contract backend completed" in result["log"]
    assert json.loads((directory / "dnakit-run.json").read_text())["status"] == "completed"
    assert not (directory / "injected").exists()
    assert not (directory / "injected2").exists()


def test_nonzero_exit_records_failed_manifest_and_backend_diagnostics(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    backend = _executable(tmp_path, "print('diagnostic: sample mismatch'); sys.exit(7)\n")
    directory = tmp_path / "failed run"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    manifest = json.loads((directory / "dnakit-run.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["error_type"] == "BackendExecutionError"
    assert "diagnostic: sample mismatch" in json.dumps(manifest)
    assert "artifacts" not in manifest


def test_exit_zero_without_expected_outputs_is_not_success(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    backend = _executable(tmp_path, "print('no data written')\n")
    directory = tmp_path / "missing outputs"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    assert json.loads((directory / "dnakit-run.json").read_text())["status"] == "failed"


def test_empty_pca_files_are_not_a_successful_analysis(
    tmp_path: Path, plink_inputs: dict[str, str]
) -> None:
    backend = _pca_backend(tmp_path, empty=True)
    directory = tmp_path / "empty pca"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    assert json.loads((directory / "dnakit-run.json").read_text())["status"] == "failed"


def test_empty_fastani_output_is_valid_no_reportable_pairs(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    # FastANI's official contract omits pairs below its reporting threshold.
    backend = _executable(tmp_path, "Path(args[args.index('-o') + 1]).write_text('')\n")
    result = run_analysis(
        "fastani",
        fasta_inputs,
        str(tmp_path / "no ani pairs"),
        executable_paths={"fastANI": str(backend)},
        verbose=False,
    )
    assert result["status"] == "completed"
    assert Path(result["artifacts"]["ani.tsv"]).stat().st_size == 0


def test_stdout_artifact_does_not_mix_stderr(tmp_path: Path, fasta_inputs: dict[str, str]) -> None:
    backend = _executable(
        tmp_path,
        "print('##gff-version 3'); print('diagnostic only', file=sys.stderr)\n",
    )
    result = run_analysis(
        "barrnap-rrna",
        {"assembly": fasta_inputs["query"]},
        str(tmp_path / "stdout run"),
        executable_paths={"barrnap": str(backend)},
        verbose=False,
    )
    assert Path(result["artifacts"]["rrna.gff"]).read_text() == "##gff-version 3\n"
    assert "diagnostic only" in result["log"]


def test_timeout_is_failed_and_keeps_manifest(tmp_path: Path, plink_inputs: dict[str, str]) -> None:
    backend = _executable(tmp_path, "time.sleep(5)\n")
    directory = tmp_path / "timed out"
    with pytest.raises(BackendTimeoutError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            timeout_seconds=0.05,
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    manifest = json.loads((directory / "dnakit-run.json").read_text())
    assert manifest["status"] == "failed" and manifest["error_type"] == "BackendTimeoutError"


def test_oversized_expected_artifact_is_failed(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    backend = _executable(tmp_path, "Path(args[args.index('-o') + 1]).write_text('x' * 1024)\n")
    directory = tmp_path / "too large"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "fastani",
            fasta_inputs,
            str(directory),
            max_artifact_bytes=128,
            executable_paths={"fastANI": str(backend)},
            verbose=False,
        )
    assert json.loads((directory / "dnakit-run.json").read_text())["status"] == "failed"


def test_expected_output_cannot_be_symlink_outside_workspace(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    external = tmp_path / "outside.tsv"
    external.write_text("preexisting", encoding="utf-8")
    backend = _executable(
        tmp_path, f"Path(args[args.index('-o') + 1]).symlink_to({str(external)!r})\n"
    )
    directory = tmp_path / "linked result"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "fastani",
            fasta_inputs,
            str(directory),
            executable_paths={"fastANI": str(backend)},
            verbose=False,
        )
    assert external.read_text() == "preexisting"


def test_catalog_results_do_not_allow_callers_to_mutate_templates() -> None:
    original = describe_analysis("plink2-pca", verbose=False)
    changed = describe_analysis("plink2-pca", verbose=False)
    changed["steps"][0]["arguments"][0] = "--unreviewed-flag"
    assert describe_analysis("plink2-pca", verbose=False) == original
    found = list_analyses("plink2-pca", verbose=False)
    assert len(found) == 1
    found[0]["steps"].clear()
    assert describe_analysis("plink2-pca", verbose=False) == original


def test_modkit_requires_prebuilt_reference_index_without_creating_it(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    bam = tmp_path / "aligned.bam"
    bam.write_bytes(b"contract fixture")
    inputs = {"alignment": str(bam), "reference": fasta_inputs["reference"]}
    directory = tmp_path / "modkit output"
    index = Path(fasta_inputs["reference"] + ".fai")
    with pytest.raises(ConfigurationError, match="index"):
        plan_analysis("modkit-pileup", inputs, str(directory), verbose=False)
    assert not directory.exists() and not index.exists()
    index.write_text("reference\t8\t11\t8\t9\n", encoding="utf-8")
    plan = plan_analysis("modkit-pileup", inputs, str(directory), verbose=False)
    arguments = plan["steps"][0]["arguments"]
    assert arguments[:3] == ["pileup", str(bam), str(directory / "methylation.bed")]
    reference_flag = "--reference" if "--reference" in arguments else "--ref"
    assert arguments[arguments.index(reference_flag) + 1] == fasta_inputs["reference"]


@pytest.mark.parametrize("timeout", [True, "30", None, float("nan"), float("inf"), 0, -1])
def test_timeout_validation_happens_before_creating_output(
    tmp_path: Path, plink_inputs: dict[str, str], timeout: Any
) -> None:
    directory = tmp_path / "bad timeout"
    backend = _pca_backend(tmp_path)
    with pytest.raises(ConfigurationError):
        run_analysis(
            "plink2-pca",
            plink_inputs,
            str(directory),
            timeout_seconds=timeout,
            executable_paths={"plink2": str(backend)},
            verbose=False,
        )
    assert not directory.exists()


def test_wildcard_artifacts_must_match_the_reviewed_output_prefix(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    backend = _executable(
        tmp_path,
        "prefix = args[args.index('-o') + 1]\n"
        "Path(prefix + '.bp.p_ctg.gfa').write_text('H\\tVN:Z:1.0\\nS\\tctg\\tACGT\\n')\n",
    )
    result = run_analysis(
        "hifiasm-assembly",
        {"reads": fasta_inputs["query"]},
        str(tmp_path / "assembly"),
        executable_paths={"hifiasm": str(backend)},
        verbose=False,
    )
    assert result["status"] == "completed"
    assert "assembly.bp.p_ctg.gfa" in result["artifacts"]
    wrong_backend = _executable(
        tmp_path, "Path('unrelated.gfa').write_text('H\\tVN:Z:1.0\\n')\n", name="wrong backend"
    )
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "hifiasm-assembly",
            {"reads": fasta_inputs["query"]},
            str(tmp_path / "wrong assembly"),
            executable_paths={"hifiasm": str(wrong_backend)},
            verbose=False,
        )


def test_unexpected_output_cannot_escape_through_returned_artifacts(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    external = tmp_path / "existing data.tsv"
    external.write_text("preexisting", encoding="utf-8")
    backend = _executable(
        tmp_path,
        "Path(args[args.index('-o') + 1]).write_text('query\\treference\\t99\\t1\\t1\\n')\n"
        f"Path('extra.tsv').symlink_to({str(external)!r})\n",
    )
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "fastani",
            fasta_inputs,
            str(tmp_path / "extra symlink"),
            executable_paths={"fastANI": str(backend)},
            verbose=False,
        )
    assert external.read_text() == "preexisting"


def test_unexpected_output_is_subject_to_reported_artifact_size_limit(
    tmp_path: Path, fasta_inputs: dict[str, str]
) -> None:
    backend = _executable(
        tmp_path,
        "Path(args[args.index('-o') + 1]).write_text('query\\treference\\t99\\t1\\t1\\n')\n"
        "Path('extra.bin').write_bytes(b'x' * 16384)\n",
    )
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "fastani",
            fasta_inputs,
            str(tmp_path / "extra oversized"),
            max_artifact_bytes=4096,
            executable_paths={"fastANI": str(backend)},
            verbose=False,
        )
