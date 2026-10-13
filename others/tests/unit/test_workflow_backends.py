"""Contract tests for pipeline-definition adapters (WDL, CWL, Snakemake, nf-core, Galaxy).

These functions parse and validate workflow definitions. They do not compute on DNA.
"""

from __future__ import annotations

import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

from dnakit.exceptions import BackendError, BackendUnavailableError, ConfigurationError
from dnakit.workflows import backends
from dnakit.workflows.backends import (
    check_wdl_inputs,
    list_galaxy_histories,
    list_galaxy_tools,
    list_galaxy_workflows,
    list_nfcore_processes,
    parse_wdl_document,
    plan_snakemake_jobs,
    read_cwl_workflow,
    read_snakefile,
)

WDL_TEXT = """version 1.0

task Align {
  input {
    File reads
    String sample_id
    Int threads = 4
  }
  command <<<
    echo ~{sample_id} > out.txt
  >>>
  output {
    File bam = "out.txt"
  }
  runtime {
    cpu: threads
  }
}

task Count {
  input {
    File bam
  }
  command <<<
    wc -l ~{bam} > count.txt
  >>>
  output {
    Int lines = 1
  }
}

workflow Pipeline {
  input {
    File reads_file
    String sample
  }
  call Align { input: reads = reads_file, sample_id = sample }
  call Count { input: bam = Align.bam }
  output {
    File aligned = Align.bam
  }
}
"""

CWL_WORKFLOW = """cwlVersion: v1.2
class: Workflow
inputs:
  reads: File
  sample: string
outputs:
  bam:
    type: File
    outputSource: align/out
steps:
  align:
    run: align.cwl
    in:
      reads: reads
      sample: sample
    out: [out]
  count:
    run: count.cwl
    in:
      bam: align/out
    out: [lines]
"""

CWL_ALIGN = """cwlVersion: v1.2
class: CommandLineTool
baseCommand: [cp]
inputs:
  reads:
    type: File
    inputBinding: {position: 1}
  sample:
    type: string
outputs:
  out:
    type: stdout
stdout: aligned.txt
"""

CWL_COUNT = """cwlVersion: v1.2
class: CommandLineTool
baseCommand: [wc, -l]
inputs:
  bam:
    type: File
    inputBinding: {position: 1}
outputs:
  lines:
    type: stdout
stdout: count.txt
"""

SNAKEFILE = """include: "rules/extra.smk"

rule all:
    input:
        expand("counts/{sample}.txt", sample=["A", "B"])

rule align:
    input:
        "reads/{sample}.fastq"
    output:
        "bam/{sample}.bam"
    shell:
        "cp {input} {output}"
"""

SNAKE_RULES = """rule count:
    input:
        "bam/{sample}.bam"
    output:
        "counts/{sample}.txt"
    shell:
        "wc -l {input} > {output}"
"""


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _snake_project(
    root: Path, *, reads: tuple[str, ...] = ("A", "B"), bam: tuple[str, ...] = ()
) -> Path:
    snakefile = _write(root / "Snakefile", SNAKEFILE)
    _write(root / "rules" / "extra.smk", SNAKE_RULES)
    for sample in reads:
        _write(root / "reads" / f"{sample}.fastq", "ACGT\n")
    for sample in bam:
        _write(root / "bam" / f"{sample}.bam", "bam\n")
    return snakefile


def _process(name: str, inputs: str, script: str) -> str:
    return (
        f"process {name} {{\n    input:\n    {inputs}\n    script:\n"
        f'    """\n    {script}\n    """\n}}\n'
    )


def _nfcore_checkout(root: Path) -> Path:
    _write(root / "main.nf", "nextflow.enable.dsl = 2\n")
    _write(
        root / "modules" / "local" / "count_lines.nf",
        _process("COUNT_LINES", "path reads", "wc -l ${reads}"),
    )
    _write(
        root / "modules" / "local" / "tabix" / "main.nf",
        _process("TABIX", "path vcf", "tabix ${vcf}"),
    )
    _write(
        root / "modules" / "nf-core" / "fastqc" / "main.nf",
        _process("FASTQC", "tuple val(meta), path(reads)", "fastqc ${reads}"),
    )
    _write(
        root / "modules" / "nf-core" / "samtools" / "sort" / "main.nf",
        _process("SAMTOOLS_SORT", "tuple val(meta), path(bam)", "samtools sort ${bam}"),
    )
    return root


class _FakeGalaxy:
    """Replaces requests.get so bioblend runs its real client code against canned JSON."""

    def __init__(self) -> None:
        self.routes: dict[str, tuple[int, Any]] = {}
        self.calls: list[dict[str, Any]] = []

    def __call__(self, url: str, **kwargs: Any) -> Any:
        import requests

        self.calls.append({"url": url, "headers": dict(kwargs.get("headers") or {})})
        path = url.split("/api/", 1)[1].split("?", 1)[0]
        status, body = self.routes.get(path, (404, {"err_msg": "not found"}))
        response = requests.Response()
        response.status_code = status
        response.headers["Content-Type"] = "application/json"
        response.url = url
        response._content = json.dumps(body).encode("utf-8")
        return response


@pytest.fixture
def galaxy(monkeypatch: pytest.MonkeyPatch) -> _FakeGalaxy:
    requests = pytest.importorskip("requests")
    pytest.importorskip("bioblend")
    fake = _FakeGalaxy()
    monkeypatch.setattr(requests, "get", fake)
    return fake


# Public contract -------------------------------------------------------------


def test_public_functions_are_exported_with_keyword_only_verbose() -> None:
    assert set(backends.__all__) == {
        "check_wdl_inputs",
        "list_galaxy_histories",
        "list_galaxy_tools",
        "list_galaxy_workflows",
        "list_nfcore_processes",
        "parse_wdl_document",
        "plan_snakemake_jobs",
        "read_cwl_workflow",
        "read_snakefile",
    }
    for name in backends.__all__:
        parameter = inspect.signature(getattr(backends, name)).parameters["verbose"]
        assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
        assert parameter.default is True


def test_missing_wdl_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path / "pipe.wdl", WDL_TEXT)
    monkeypatch.setitem(sys.modules, "WDL", None)
    with pytest.raises(BackendUnavailableError, match="bio-wdl"):
        parse_wdl_document(path, verbose=False)


def test_missing_cwltool_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write(tmp_path / "wf.cwl", CWL_WORKFLOW)
    monkeypatch.setitem(sys.modules, "cwltool", None)
    with pytest.raises(BackendUnavailableError, match="bio-cwl"):
        read_cwl_workflow(path, verbose=False)


def test_missing_snakemake_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    snakefile = _snake_project(tmp_path)
    monkeypatch.setitem(sys.modules, "snakemake", None)
    with pytest.raises(BackendUnavailableError, match="bio-snakemake"):
        read_snakefile(snakefile, verbose=False)


def test_missing_nf_core_package_raises_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkout = _nfcore_checkout(tmp_path / "pipe")
    monkeypatch.setitem(sys.modules, "nf_core", None)
    with pytest.raises(BackendUnavailableError, match="bio-nf-core"):
        list_nfcore_processes(checkout, verbose=False)


def test_missing_bioblend_package_raises_install_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "bioblend", None)
    with pytest.raises(BackendUnavailableError, match="bio-bioblend"):
        list_galaxy_tools("https://galaxy.example.org", verbose=False)


# miniwdl ---------------------------------------------------------------------


def test_parse_wdl_document_lists_tasks_inputs_outputs_and_calls(tmp_path: Path) -> None:
    pytest.importorskip("WDL")
    path = _write(tmp_path / "pipe.wdl", WDL_TEXT)

    result = parse_wdl_document(path, verbose=False)

    assert result["status"] == "completed"
    assert result["operation"] == "parse_wdl_document"
    assert result["wdl_version"] == "1.0"
    assert result["imports"] == []
    assert [task["name"] for task in result["tasks"]] == ["Align", "Count"]
    align, count = result["tasks"]
    # Inputs with a default (threads) are optional; the others are required.
    assert align["inputs"] == [
        {"name": "reads", "type": "File", "required": True},
        {"name": "sample_id", "type": "String", "required": True},
        {"name": "threads", "type": "Int", "required": False},
    ]
    assert align["outputs"] == [{"name": "bam", "type": "File"}]
    assert count["inputs"] == [{"name": "bam", "type": "File", "required": True}]
    workflow = result["workflow"]
    assert workflow["name"] == "Pipeline"
    assert workflow["inputs"] == [
        {"name": "Align.threads", "type": "Int", "required": False},
        {"name": "reads_file", "type": "File", "required": True},
        {"name": "sample", "type": "String", "required": True},
    ]
    assert workflow["outputs"] == [{"name": "aligned", "type": "File"}]
    assert workflow["calls"] == [
        {"name": "Align", "task": "Align"},
        {"name": "Count", "task": "Count"},
    ]
    assert result["provenance"]["package"].startswith("miniwdl ")


def test_check_wdl_inputs_reports_missing_unknown_and_invalid_values(tmp_path: Path) -> None:
    pytest.importorskip("WDL")
    path = _write(tmp_path / "pipe.wdl", WDL_TEXT)

    valid = check_wdl_inputs(path, {"reads_file": "reads.fq", "sample": "S1"}, verbose=False)
    assert valid["valid"] is True
    assert valid["missing"] == [] and valid["unknown"] == [] and valid["invalid"] == []

    bad = check_wdl_inputs(path, {"sample": 3, "bogus": "x", "Align.threads": 8}, verbose=False)
    assert bad["valid"] is False
    assert bad["missing"] == ["reads_file"]
    assert bad["unknown"] == ["bogus"]
    assert [item["name"] for item in bad["invalid"]] == ["sample"]
    assert bad["invalid"][0]["type"] == "String"
    unhashable = check_wdl_inputs(path, {"sample": {1, 2}}, verbose=False)
    assert [item["name"] for item in unhashable["invalid"]] == ["sample"]


def test_wdl_rejects_bad_syntax_missing_imports_and_oversized_files(tmp_path: Path) -> None:
    pytest.importorskip("WDL")
    broken = _write(
        tmp_path / "broken.wdl",
        "version 1.0\ntask Broken {\n  command <<<\n    echo ~{x\n  >>>\n}\n",
    )
    with pytest.raises(ConfigurationError, match=re.escape("miniwdl rejected broken.wdl")):
        parse_wdl_document(broken, verbose=False)

    importing = _write(
        tmp_path / "importing.wdl",
        'version 1.0\nimport "missing_module.wdl" as M\nworkflow W {\n  call M.Nope\n}\n',
    )
    with pytest.raises(ConfigurationError, match=re.escape("miniwdl rejected importing.wdl")):
        parse_wdl_document(importing, verbose=False)

    path = _write(tmp_path / "pipe.wdl", WDL_TEXT)
    with pytest.raises(ConfigurationError, match="max_bytes"):
        parse_wdl_document(path, max_bytes=10, verbose=False)
    with pytest.raises(ConfigurationError, match="max_bytes"):
        parse_wdl_document(path, max_bytes=True, verbose=False)


# cwltool ---------------------------------------------------------------------


def test_read_cwl_workflow_lists_steps_sources_and_outputs(tmp_path: Path) -> None:
    pytest.importorskip("cwltool")
    workflow = _write(tmp_path / "wf.cwl", CWL_WORKFLOW)
    _write(tmp_path / "align.cwl", CWL_ALIGN)
    _write(tmp_path / "count.cwl", CWL_COUNT)

    result = read_cwl_workflow(workflow, verbose=False)

    assert result["class"] == "Workflow"
    assert result["inputs"] == [
        {"name": "reads", "type": "File"},
        {"name": "sample", "type": "string"},
    ]
    steps = {step["name"]: step for step in result["steps"]}
    assert list(steps) == ["align", "count"]
    assert steps["align"]["run"] == "align.cwl"
    assert steps["align"]["run_class"] == "CommandLineTool"
    assert steps["align"]["outputs"] == ["out"]
    assert steps["align"]["inputs"] == [
        {"name": "reads", "sources": ["reads"]},
        {"name": "sample", "sources": ["sample"]},
    ]
    assert steps["count"]["inputs"] == [{"name": "bam", "sources": ["align/out"]}]
    assert result["outputs"] == [{"name": "bam", "type": "File", "sources": ["align/out"]}]
    assert result["provenance"]["package"].startswith("cwltool ")


def test_read_cwl_workflow_rejects_invalid_documents_and_tools(tmp_path: Path) -> None:
    pytest.importorskip("cwltool")
    invalid = _write(
        tmp_path / "bad.cwl",
        "cwlVersion: v1.2\nclass: Workflow\ninputs:\n  reads: Fyle\nsteps: []\n",
    )
    with pytest.raises(ConfigurationError, match=re.escape("cwltool rejected bad.cwl")):
        read_cwl_workflow(invalid, verbose=False)

    tool = _write(tmp_path / "align.cwl", CWL_ALIGN)
    with pytest.raises(ConfigurationError, match="needs class: Workflow"):
        read_cwl_workflow(tool, verbose=False)


def test_read_cwl_workflow_steps_do_not_depend_on_hash_seed(tmp_path: Path) -> None:
    pytest.importorskip("cwltool")
    _write(tmp_path / "align.cwl", CWL_ALIGN)
    _write(tmp_path / "count.cwl", CWL_COUNT)
    workflow = _write(tmp_path / "wf.cwl", CWL_WORKFLOW)
    script = (
        "import json, sys\n"
        "from dnakit.workflows.backends import read_cwl_workflow\n"
        "print(json.dumps(read_cwl_workflow(sys.argv[1], verbose=False), sort_keys=True))\n"
    )
    outputs = set()
    for seed in ("0", "1", "3"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        completed = subprocess.run(
            [sys.executable, "-c", script, str(workflow)],
            capture_output=True,
            text=True,
            env=env,
            check=True,
        )
        outputs.add(completed.stdout.strip().splitlines()[-1])
    assert len(outputs) == 1


# Snakemake -------------------------------------------------------------------


def test_read_snakefile_lists_included_rules_sorted_by_name(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    snakefile = _snake_project(tmp_path)
    before = Path.cwd()

    result = read_snakefile(snakefile, verbose=False)

    assert Path.cwd() == before
    assert result["default_target"] == "all"
    assert result["rules"] == [
        {"name": "align", "inputs": ["reads/{sample}.fastq"], "outputs": ["bam/{sample}.bam"]},
        {"name": "all", "inputs": ["counts/A.txt", "counts/B.txt"], "outputs": []},
        {"name": "count", "inputs": ["bam/{sample}.bam"], "outputs": ["counts/{sample}.txt"]},
    ]
    assert not (tmp_path / ".snakemake").exists()  # parsing writes only to a scratch directory


def test_read_snakefile_reports_syntax_errors(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    snakefile = _write(tmp_path / "Snakefile", 'rule a:\n    input: "x"\n    output "y"\n')
    with pytest.raises(ConfigurationError, match="Snakemake rejected Snakefile"):
        read_snakefile(snakefile, verbose=False)


def test_plan_snakemake_jobs_matches_hand_calculation(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    # bam/A.bam exists, so only count A needs no align. For all targets: count A,
    # align B, count B and the all rule itself make four jobs.
    snakefile = _snake_project(tmp_path, bam=("A",))

    result = plan_snakemake_jobs(snakefile, verbose=False)

    assert result["job_count"] == 4
    assert result["jobs_by_rule"] == {"align": 1, "all": 1, "count": 2}
    assert [(job["rule"], job["outputs"]) for job in result["jobs"]] == [
        ("align", ["bam/B.bam"]),
        ("all", []),
        ("count", ["counts/A.txt"]),
        ("count", ["counts/B.txt"]),
    ]
    assert result["jobs"][0]["wildcards"] == {"sample": "B"}
    assert result["jobs"][0]["inputs"] == ["reads/B.fastq"]
    assert not (tmp_path / "bam" / "B.bam").exists()  # dry run creates no outputs


def test_plan_snakemake_jobs_accepts_targets_and_reports_limits(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    snakefile = _snake_project(tmp_path, bam=("A",))

    single = plan_snakemake_jobs(snakefile, targets=["counts/A.txt"], verbose=False)
    assert single["targets"] == ["counts/A.txt"]
    assert [(job["rule"], job["outputs"]) for job in single["jobs"]] == [
        ("count", ["counts/A.txt"])
    ]

    with pytest.raises(ConfigurationError, match="max_jobs"):
        plan_snakemake_jobs(snakefile, max_jobs=1, verbose=False)
    with pytest.raises(ConfigurationError, match="No rule to produce"):
        plan_snakemake_jobs(snakefile, targets=["missing.txt"], verbose=False)


def test_plan_snakemake_jobs_fails_on_missing_input(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    snakefile = _snake_project(tmp_path, reads=("A",), bam=("A",))
    with pytest.raises(ConfigurationError, match="Missing input files for rule align"):
        plan_snakemake_jobs(snakefile, targets=["counts/B.txt"], verbose=False)


def test_plan_snakemake_jobs_is_stable_across_hash_seeds(tmp_path: Path) -> None:
    pytest.importorskip("snakemake")
    snakefile = _snake_project(tmp_path, bam=("A",))
    script = (
        "import json, sys\n"
        "from dnakit.workflows.backends import plan_snakemake_jobs\n"
        "print(json.dumps(plan_snakemake_jobs(sys.argv[1], verbose=False), sort_keys=True))\n"
    )
    outputs = set()
    for seed in ("0", "2"):
        completed = subprocess.run(
            [sys.executable, "-c", script, str(snakefile)],
            capture_output=True,
            text=True,
            env={**os.environ, "PYTHONHASHSEED": seed},
            check=True,
        )
        outputs.add(completed.stdout.strip().splitlines()[-1])
    assert len(outputs) == 1


# nf-core ---------------------------------------------------------------------


def test_list_nfcore_processes_reads_local_and_nf_core_modules(tmp_path: Path) -> None:
    pytest.importorskip("nf_core")
    checkout = _nfcore_checkout(tmp_path / "pipe")

    result = list_nfcore_processes(checkout, verbose=False)

    assert result["component_count"] == 4
    assert result["process_count"] == 4
    assert [(row["kind"], row["name"], row["process"]) for row in result["components"]] == [
        ("local", "count_lines", "COUNT_LINES"),
        ("local", "tabix", "TABIX"),
        ("nf-core", "fastqc", "FASTQC"),
        ("nf-core", "samtools/sort", "SAMTOOLS_SORT"),
    ]
    assert result["components"][1]["path"] == "modules/local/tabix/main.nf"
    assert result["components"][0]["path"] == "modules/local/count_lines.nf"
    assert result["provenance"]["package"].startswith("nf-core ")


def test_list_nfcore_processes_rejects_non_pipeline_directories(tmp_path: Path) -> None:
    pytest.importorskip("nf_core")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ConfigurationError, match="not an nf-core pipeline checkout"):
        list_nfcore_processes(empty, verbose=False)
    with pytest.raises(ConfigurationError, match="existing directory"):
        list_nfcore_processes(tmp_path / "absent", verbose=False)
    checkout = _nfcore_checkout(tmp_path / "pipe")
    with pytest.raises(ConfigurationError, match="max_components"):
        list_nfcore_processes(checkout, max_components=1, verbose=False)


# Galaxy (bioblend, mocked server) --------------------------------------------


def test_list_galaxy_tools_uses_bioblend_and_sorts_by_id(galaxy: _FakeGalaxy) -> None:
    galaxy.routes["tools"] = (
        200,
        [
            {"id": "cat1", "name": "Concatenate", "version": "1.0.0", "panel_section_name": "Text"},
            {"id": "align0", "name": "Align", "version": "2.1", "panel_section_name": None},
        ],
    )

    result = list_galaxy_tools("https://galaxy.example.org/", api_key="SECRET-KEY", verbose=False)

    assert result["server"] == "https://galaxy.example.org"
    assert result["count"] == 2
    assert [item["id"] for item in result["items"]] == ["align0", "cat1"]
    assert result["items"][0] == {
        "id": "align0",
        "name": "Align",
        "version": "2.1",
        "section": None,
    }
    assert galaxy.calls[0]["url"] == "https://galaxy.example.org/api/tools"
    assert galaxy.calls[0]["headers"]["x-api-key"] == "SECRET-KEY"
    assert "SECRET-KEY" not in json.dumps(result)


def test_list_galaxy_histories_and_workflows_parse_server_rows(galaxy: _FakeGalaxy) -> None:
    galaxy.routes["histories"] = (
        200,
        [
            {"id": "h2", "name": "b", "count": 4, "deleted": False, "published": False},
            {"id": "h1", "name": "a", "count": 3, "deleted": False, "published": True},
        ],
    )
    galaxy.routes["workflows"] = (200, [{"id": "w1", "name": "wgs", "published": False}])

    histories = list_galaxy_histories("https://galaxy.example.org", verbose=False)
    workflows = list_galaxy_workflows("https://galaxy.example.org", verbose=False)

    assert [item["id"] for item in histories["items"]] == ["h1", "h2"]
    assert histories["items"][0] == {
        "id": "h1",
        "name": "a",
        "count": 3,
        "deleted": False,
        "published": True,
    }
    assert workflows["items"] == [{"id": "w1", "name": "wgs", "published": False}]
    assert workflows["provenance"]["package"].startswith("bioblend ")


def test_galaxy_rejects_credentials_in_url_and_oversized_responses(galaxy: _FakeGalaxy) -> None:
    pytest.importorskip("bioblend")
    for url in (
        "ftp://galaxy.example.org",
        "https://user:pw@galaxy.example.org",
        "https://galaxy.example.org/?key=1",
    ):
        with pytest.raises(ConfigurationError):
            list_galaxy_tools(url, verbose=False)
    with pytest.raises(ConfigurationError, match="max_items"):
        list_galaxy_tools("https://galaxy.example.org", max_items=0, verbose=False)

    galaxy.routes["tools"] = (200, [{"id": "a"}, {"id": "b"}])
    with pytest.raises(ConfigurationError, match="max_items is 1"):
        list_galaxy_tools("https://galaxy.example.org", max_items=1, verbose=False)


def test_galaxy_http_and_network_failures_raise_backend_error(galaxy: _FakeGalaxy) -> None:
    galaxy.routes["tools"] = (500, {"err_msg": "boom", "err_code": 0})
    with pytest.raises(BackendError, match="Galaxy tools request failed"):
        list_galaxy_tools("https://galaxy.example.org", api_key="K", verbose=False)


@pytest.mark.skipif(
    os.environ.get("DNAKIT_LIVE_TESTS") != "1",
    reason="set DNAKIT_LIVE_TESTS=1 for a live Galaxy server",
)
def test_live_galaxy_tools_listing() -> None:
    url = os.environ.get("DNAKIT_GALAXY_URL", "https://usegalaxy.org")
    result = list_galaxy_tools(url, verbose=False)
    assert result["count"] >= 1
