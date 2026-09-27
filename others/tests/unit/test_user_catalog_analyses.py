"""Review-backed CLI contracts for the user's catalog, without scientific backends.

Tiny executables test actual subprocess/argv and artifact handling only. These
tests do not claim biological validation; optional real-backend smoke tests live
in test_user_catalog_analysis_backends.py.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.integrations.catalog import describe_analysis
from dnakit.integrations.execution import plan_analysis, run_analysis

# Independent expected filenames and distinctive CLI arguments reviewed against
# the original manuals/source. Wildcards are instantiated with "fixture" below.
CONTRACTS = {
    "bakta-annotation": (
        "bakta",
        ["--db", "--skip-plot"],
        [
            "bakta/annotation.gff3",
            "bakta/annotation.tsv",
            "bakta/annotation.json",
            "bakta/annotation.fna",
        ],
    ),
    "busco-proteins": (
        "busco",
        ["proteins", "--offline", "--opt-out-run-stats"],
        ["assessment/short_summary*.json", "assessment/short_summary*.txt"],
    ),
    "busco-genome": (
        "busco",
        ["genome", "--offline", "--opt-out-run-stats"],
        ["assessment/short_summary*.json", "assessment/short_summary*.txt"],
    ),
    "augustus-genes": ("augustus", ["--gff3=on", "--strand=both"], ["genes.gff3"]),
    "diamond-blastp": (
        "diamond",
        ["blastp", "--outfmt", "--max-target-seqs"],
        ["reference.dmnd", "hits.tsv"],
    ),
    "diamond-blastx": (
        "diamond",
        ["blastx", "--outfmt", "--max-target-seqs"],
        ["reference.dmnd", "hits.tsv"],
    ),
    "eggnog-proteins": (
        "eggnog-mapper",
        ["--itype", "proteins", "--data_dir"],
        ["annotation.emapper.annotations", "annotation.emapper.seed_orthologs"],
    ),
    "miniprot-align": ("miniprot", ["--gff"], ["alignment.gff3"]),
    "liftoff-annotation": ("liftoff", ["-db", "-dir", "-u"], ["lifted.gff3", "unmapped.txt"]),
    "checkm2-predict": (
        "checkm2",
        ["predict", "--database_path", "--lowmem"],
        ["checkm2/quality_report.tsv"],
    ),
    "quast-assembly": ("quast", ["--reference", "--no-icarus", "--no-html"], ["quast/report.tsv"]),
    "vsearch-dereplicate": (
        "vsearch",
        ["--derep_fulllength", "--sizeout"],
        ["dereplicated.fasta", "clusters.uc"],
    ),
    "vsearch-cluster": (
        "vsearch",
        ["--cluster_fast", "--centroids"],
        ["centroids.fasta", "clusters.uc"],
    ),
    "mmseqs-protein-search": (
        "mmseqs2",
        ["easy-search", "--search-type", "--split-memory-limit"],
        ["hits.tsv"],
    ),
    "mmseqs-protein-cluster": (
        "mmseqs2",
        ["easy-cluster", "--cov-mode", "--split-memory-limit"],
        ["clusters_cluster.tsv", "clusters_rep_seq.fasta", "clusters_all_seqs.fasta"],
    ),
    "abricate-screen": (
        "abricate",
        ["--datadir", "--minid", "--mincov", "--threads"],
        ["hits.tsv"],
    ),
    "mlst-typing": ("mlst", ["--blastdb", "--datadir", "--json"], ["mlst.tsv", "mlst.json"]),
    "genomad-mobile-elements": (
        "genomad",
        ["end-to-end", "--splits", "--cleanup"],
        ["genomad/*_summary/*_virus_summary.tsv", "genomad/*_summary/*_plasmid_summary.tsv"],
    ),
    "repeatmasker-library": ("repeatmasker", ["-lib", "-dir", "-xsmall"], ["*.out", "*.masked"]),
    "mageck-test": (
        "mageck",
        ["test", "-k", "--norm-method"],
        ["screen.gene_summary.txt", "screen.sgrna_summary.txt"],
    ),
    "mageck-count": (
        "mageck",
        ["count", "--fastq", "--sample-label"],
        ["counts.count.txt", "counts.countsummary.txt"],
    ),
    "crispresso2-amplicon": (
        "crispresso2",
        ["--amplicon_seq", "--guide_seq", "--suppress_report"],
        [
            "CRISPResso_on_amplicon/CRISPResso_quantification_of_editing_frequency.txt",
            "CRISPResso_on_amplicon/CRISPResso2_info.json",
        ],
    ),
    "drugz-unpaired": ("drugz", ["-unpaired", "--half_window_size", "-x"], ["drugz.tsv"]),
    "gemma-lmm": (
        "gemma",
        ["-bfile", "-k", "-c", "-lmm", "4"],
        ["association.assoc.txt", "association.log.txt"],
    ),
}


def _inputs(tmp_path: Path, analysis: str) -> tuple[dict[str, str], dict]:
    definition = describe_analysis(analysis, verbose=False)
    result = {}
    for name, schema in definition["inputs"].items():
        path = tmp_path / f"input {name} ; literal"
        if schema["type"] == "directory":
            path.mkdir()
            for child in schema.get("contains", []):
                (path / child).parent.mkdir(parents=True, exist_ok=True)
                (path / child).write_text("contract fixture\n")
        elif schema["type"] == "prefix":
            for suffix in schema["suffixes"]:
                Path(str(path) + suffix).write_text("contract fixture\n")
        else:
            path.write_text("contract fixture\n")
        for suffix in schema.get("required_suffixes", []):
            Path(str(path) + suffix).write_text("contract fixture\n")
        result[name] = str(path)
    parameters = {
        name: schema["example"]
        for name, schema in definition["parameters"].items()
        if "default" not in schema
    }
    return result, parameters


@pytest.mark.parametrize("analysis", CONTRACTS)
def test_reviewed_contract_subprocess_and_provenance(tmp_path: Path, analysis: str) -> None:
    inputs, parameters = _inputs(tmp_path, analysis)
    output = tmp_path / "output with spaces ; literal"
    project_id, flags, expected = CONTRACTS[analysis]
    plan = plan_analysis(analysis, inputs, str(output), parameters=parameters, verbose=False)
    assert plan["project_id"] == project_id
    assert not output.exists()
    arguments = [arg for step in plan["steps"] for arg in step["arguments"]]
    assert set(flags) <= set(arguments)
    assert [p for s in plan["steps"] for p in s["expected_outputs"]] == expected
    for value in inputs.values():
        assert any(arg == value or arg.endswith("=" + value) for arg in arguments)
    program = tmp_path / "contract executable with spaces"
    paths = [p.replace("*", "fixture") for p in expected]
    program.write_text(
        f"#!{sys.executable}\nimport json, sys\nfrom pathlib import Path\n"
        f"root = Path({str(output)!r})\n"
        f"for name in {paths!r}:\n"
        "    path = root / name\n"
        "    path.parent.mkdir(parents=True, exist_ok=True)\n"
        "    path.write_text('contract output\\n')\n"
        "with (root / 'argv.jsonl').open('a') as handle:\n"
        "    handle.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "print('contract stdout')\n"
    )
    program.chmod(0o755)
    overrides = {s["program"]: str(program) for s in plan["steps"]}
    result = run_analysis(
        analysis,
        inputs,
        str(output),
        parameters=parameters,
        executable_paths=overrides,
        timeout_seconds=10,
        verbose=False,
    )
    assert result["status"] == "completed"
    assert result["project_id"] == project_id
    manifest = json.loads((output / "dnakit-run.json").read_text())
    assert manifest["project_id"] == project_id
    assert manifest["sources"] and manifest["validation"]
    invocations = [json.loads(line) for line in (output / "argv.jsonl").read_text().splitlines()]
    assert invocations == [step["arguments"] for step in plan["steps"]]


@pytest.mark.parametrize("analysis", CONTRACTS)
def test_every_new_recipe_missing_backend_is_read_only(tmp_path: Path, analysis: str) -> None:
    inputs, parameters = _inputs(tmp_path, analysis)
    definition = describe_analysis(analysis, verbose=False)
    missing = {s["program"]: str(tmp_path / "absent") for s in definition["steps"]}
    output = tmp_path / "missing backend"
    with pytest.raises(BackendUnavailableError):
        run_analysis(
            analysis,
            inputs,
            str(output),
            parameters=parameters,
            executable_paths=missing,
            verbose=False,
        )
    assert not output.exists()


@pytest.mark.parametrize(
    "analysis",
    ["augustus-genes", "abricate-screen", "mageck-test", "drugz-unpaired", "crispresso2-amplicon"],
)
def test_required_biological_parameters_have_no_inferred_default(
    tmp_path: Path, analysis: str
) -> None:
    inputs, _ = _inputs(tmp_path, analysis)
    with pytest.raises(ConfigurationError):
        plan_analysis(analysis, inputs, str(tmp_path / "output"), verbose=False)


def test_huge_numeric_input_is_configuration_error(tmp_path: Path) -> None:
    inputs, _ = _inputs(tmp_path, "vsearch-cluster")
    with pytest.raises(ConfigurationError):
        plan_analysis(
            "vsearch-cluster",
            inputs,
            str(tmp_path / "output"),
            parameters={"identity": 10**10000},
            verbose=False,
        )


@pytest.mark.parametrize("name", ["target", "reference"])
def test_liftoff_requires_prebuilt_indices_and_annotation_database(
    tmp_path: Path, name: str
) -> None:
    inputs, _ = _inputs(tmp_path, "liftoff-annotation")
    Path(inputs[name] + ".fai").unlink()
    with pytest.raises(ConfigurationError, match=r"existing \.fai"):
        plan_analysis("liftoff-annotation", inputs, str(tmp_path / "output"), verbose=False)


@pytest.mark.parametrize("suffix", [".nhr", ".nin", ".nsq"])
def test_mlst_requires_complete_explicit_blast_database(tmp_path: Path, suffix: str) -> None:
    inputs, _ = _inputs(tmp_path, "mlst-typing")
    Path(inputs["blast_database"] + suffix).unlink()
    with pytest.raises(ConfigurationError, match="Incomplete"):
        plan_analysis("mlst-typing", inputs, str(tmp_path / "output"), verbose=False)


def test_zero_exit_without_scientific_artifact_is_failed_manifest(tmp_path: Path) -> None:
    inputs, _ = _inputs(tmp_path, "vsearch-dereplicate")
    program = tmp_path / "empty backend"
    program.write_text(f"#!{sys.executable}\nprint('no scientific outputs')\n")
    program.chmod(0o755)
    output = tmp_path / "failed output"
    with pytest.raises(BackendExecutionError):
        run_analysis(
            "vsearch-dereplicate",
            inputs,
            str(output),
            executable_paths={"vsearch": str(program)},
            verbose=False,
        )
    manifest = json.loads((output / "dnakit-run.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["project_id"] == "vsearch"


def test_existing_output_directory_remains_untouched(tmp_path: Path) -> None:
    inputs, _ = _inputs(tmp_path, "diamond-blastp")
    output = tmp_path / "existing"
    output.mkdir()
    sentinel = output / "existing result"
    sentinel.write_text("keep me")
    with pytest.raises(ConfigurationError, match="already exists"):
        plan_analysis("diamond-blastp", inputs, str(output), verbose=False)
    assert sentinel.read_text() == "keep me"
