"""User-facing discovery, resource integrity, and non-destructive skill installation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from dnakit.cli.app import app
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.skills import catalog as catalog_module
from dnakit.skills import (
    check_skill,
    cite_skill,
    get_skill,
    install_skill,
    list_skill_categories,
    list_skills,
    read_skill_file,
    search_skills,
    skill_catalog_stats,
)


def test_chinese_search_finds_the_actual_task_and_retains_duplicate_aliases() -> None:
    result = search_skills("变异注释", limit=5)
    assert "variant-annotation" in {x["id"] for x in result["skills"]}
    assert get_skill("tooluniverse-fastq-qc")["id"] == "read-qc-quality-reports"
    assert get_skill("bwa-mem2-dna-aligner")["id"] == "read-alignment-bwa-alignment"
    assert (
        get_skill("bio-read-alignment-bowtie2-alignment")["capability"]
        != get_skill("bio-read-alignment-bwa-alignment")["capability"]
    )
    assert search_skills("zyxwvnonexistenttask")["total"] == 0


@pytest.mark.parametrize(
    "query,expected",
    [
        ("整合子", "integronfinder"),
        ("质粒", "mob-suite"),
        ("基因组大模型", "genome-language-model"),
        ("古DNA", "grapevine-adna"),
    ],
)
def test_expanded_dna_tasks_are_discoverable_in_chinese(query: str, expected: str) -> None:
    assert expected in {entry["id"] for entry in search_skills(query, limit=5)["skills"]}


def test_same_named_workflows_cannot_redirect_installation_arbitrarily(tmp_path: Path) -> None:
    assert get_skill("alignment-and-mapping")["id"] == "read-alignment-bwa-alignment"
    assert get_skill("scikit-bio")["id"] == "scikit-bio"
    with pytest.raises(ConfigurationError, match="ambiguous"):
        install_skill("phylogenetics", destination=tmp_path / "unselected")
    assert not (tmp_path / "unselected").exists()
    found = {entry["id"] for entry in search_skills("variant-annotation")["skills"]}
    assert {"variant-annotation", "regulatory-variant-analysis"} <= found


def test_expansion_installs_portable_resources_and_preserves_invocation_policy(
    tmp_path: Path,
) -> None:
    result = install_skill("dnabert2", destination=tmp_path)
    root = Path(result["path"])
    document = (root / "SKILL.md").read_text()
    assert "../../docs/" not in document
    assert "references/upstream-docs/references.md" in document
    assert (root / "references/upstream-docs/references.md").is_file()
    assert (root / "licenses/pku-metagenomics/LICENSE").is_file()

    result = install_skill("encode-toolkit-setup", destination=tmp_path)
    root = Path(result["path"])
    metadata = yaml.safe_load((root / "SKILL.md").read_text().split("---", 2)[1])
    assert metadata["disable-model-invocation"] is True
    assert metadata["name"] == "encode-toolkit-setup"
    assert (
        "GNU AFFERO GENERAL PUBLIC LICENSE"
        in (root / "licenses/encode-toolkit/LICENSE").read_text()
    )
    assert get_skill("encode-toolkit-setup")["license"] == "AGPL-3.0-only"

    result = install_skill("pod5", destination=tmp_path)
    root = Path(result["path"])
    assert not (root / "bin").exists() and not (root / "bin").is_symlink()
    assert "command -v pod5" in (root / "SKILL.md").read_text()
    assert "pod5" in get_skill("pod5")["runtime"]["commands"]


def test_deepmind_dual_license_and_scope_declaration_are_available() -> None:
    identifier = "alphagenome-single-variant-analysis"
    entry = get_skill(identifier)
    assert entry["license"] == "Apache-2.0 AND CC-BY-4.0"
    paths = {resource["path"] for resource in entry["files"]}
    assert {
        "licenses/deepmind/LICENSE",
        "licenses/deepmind/CC-BY-4.0.txt",
        "licenses/deepmind/README.md",
        "licenses/deepmind/SKILL_LICENSES.md",
    } <= paths
    text = read_skill_file(identifier, resource="licenses/deepmind/README.md")["content"]
    assert "CC-BY" in text and "Apache" in text


def test_browse_filters_pagination_and_metadata_isolation() -> None:
    first = list_skills(category="variants", limit=3)
    second = list_skills(category="variants", offset=3, limit=3)
    assert not {x["id"] for x in first["skills"]} & {x["id"] for x in second["skills"]}
    result = search_skills("variant", source="tooluniverse", limit=2)
    assert all(x["source_id"] == "tooluniverse" for x in result["skills"])
    categories = list_skill_categories()
    assert sum(x["count"] for x in categories.values()) == skill_catalog_stats()["count"]
    entry = get_skill("dnakit")
    entry["runtime"]["python_modules"].append("not-a-real-module")
    assert "not-a-real-module" not in get_skill("dnakit")["runtime"]["python_modules"]


@pytest.mark.parametrize(
    "kwargs", [{"limit": 0}, {"category": "bad-category"}, {"source": "no-source"}]
)
def test_invalid_search_arguments_fail_explicitly(kwargs: dict[str, object]) -> None:
    with pytest.raises(ConfigurationError):
        search_skills("variant", **kwargs)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "path", ["../LICENSE", "/etc/passwd", "licenses/../../LICENSE", "a\\b", "./SKILL.md"]
)
def test_reader_rejects_path_traversal(path: str) -> None:
    with pytest.raises(ConfigurationError) as caught:
        read_skill_file("dnakit", resource=path)
    assert caught.value.code == "SKILL_RESOURCE_INVALID"


def test_reader_pagination_reconstructs_the_document() -> None:
    first = read_skill_file("dnakit", limit=123)
    assert first["next_offset"] == 123
    remainder = read_skill_file("dnakit", offset=123, limit=100_000)
    assert first["content"] + remainder["content"] == read_skill_file("dnakit")["content"]
    assert remainder["next_offset"] is None
    assert "LICENSE" in first["files"]
    assert first["runtime"]["python_modules"] == ["dnakit"]


def test_damaged_bundle_is_rejected_before_install_writes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "bundle"
    (root / "dnakit").mkdir(parents=True)
    (root / "dnakit" / "SKILL.md").write_text("tampered")
    original_catalog = catalog_module._catalog()
    monkeypatch.setattr(catalog_module, "_catalog", lambda: original_catalog)
    monkeypatch.setattr(catalog_module, "_locations", lambda: (root, root / "catalog.json"))
    with pytest.raises(InputFormatError):
        read_skill_file("dnakit")
    with pytest.raises((InputFormatError, OSError)):
        install_skill("dnakit", destination=tmp_path / "installation")
    assert not (tmp_path / "installation").exists()


@pytest.mark.parametrize(
    "agent,relative",
    [
        ("codex", ".agents/skills"),
        ("claude-code", ".claude/skills"),
        ("cursor", ".cursor/skills"),
    ],
)
def test_client_installations_are_complete_and_idempotent(
    tmp_path: Path,
    agent: str,
    relative: str,
) -> None:
    identifier = "genomic-coordinates"
    result = install_skill(identifier, agent=agent, project_dir=tmp_path)
    directory = tmp_path / relative / identifier
    assert result["path"] == str(directory)
    for resource in get_skill(identifier)["files"]:
        assert (
            hashlib.sha256((directory / resource["path"]).read_bytes()).hexdigest()
            == resource["sha256"]
        )
    assert (directory / "UPSTREAM.md").is_file()
    assert install_skill(identifier, agent=agent, project_dir=tmp_path)["status"] == "unchanged"


def test_global_install_uses_the_selected_home(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    result = install_skill("dnakit", agent="codex", global_scope=True)
    assert result["path"] == str(tmp_path / ".agents/skills/dnakit")


def test_installer_protects_unmanaged_and_edited_directories(tmp_path: Path) -> None:
    unmanaged = tmp_path / "manual" / "dnakit"
    unmanaged.mkdir(parents=True)
    (unmanaged / "SKILL.md").write_text("user-owned skill")
    with pytest.raises(ConfigurationError):
        install_skill("dnakit", destination=unmanaged.parent, force=True)
    assert (unmanaged / "SKILL.md").read_text() == "user-owned skill"
    result = install_skill("dnakit", destination=tmp_path / "managed")
    target = Path(result["path"])
    (target / "SKILL.md").write_text("local edit")
    with pytest.raises(ConfigurationError):
        install_skill("dnakit", destination=target.parent)
    assert (target / "SKILL.md").read_text() == "local edit"
    install_skill("dnakit", destination=target.parent, force=True)
    assert (target / "SKILL.md").read_text().startswith("---\nname: dnakit")


def test_installer_rejects_symlink_targets(tmp_path: Path) -> None:
    real = tmp_path / "real"
    real.mkdir()
    (tmp_path / "dnakit").symlink_to(real, target_is_directory=True)
    with pytest.raises(ConfigurationError):
        install_skill("dnakit", destination=tmp_path, force=True)
    assert not list(real.iterdir())


def test_install_rollback_preserves_previous_skill(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    install_skill("dnakit", destination=tmp_path)
    target = tmp_path / "dnakit"
    (target / "SKILL.md").write_text("local edit")
    original = Path.rename

    def fail_staging(path: Path, destination: Path) -> Path:
        if path.name.startswith(".dnakit-") and "backup" not in path.name:
            raise OSError("simulated rename failure")
        return original(path, destination)

    monkeypatch.setattr(Path, "rename", fail_staging)
    with pytest.raises(OSError, match="simulated"):
        install_skill("dnakit", destination=tmp_path, force=True)
    assert (target / "SKILL.md").read_text() == "local edit"


def test_dependency_probes_do_not_claim_workflow_execution(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(catalog_module.shutil, "which", lambda _: None)
    result = check_skill("bio-read-qc-quality-reports")
    assert result["status"] == "missing_dependencies"
    assert result["execution_validation"] == "not_run"
    assert any(x["name"] == "fastqc" for x in result["checks"])


def test_citations_preserve_upstream_authorship_and_revision() -> None:
    entry = get_skill("genomic-coordinates")
    citation = cite_skill("genomic-coordinates", format="json")
    assert isinstance(citation, dict)
    assert citation["version"] == entry["upstream_revision"]
    assert "K-Dense" in citation["author"]
    assert "@software{" in cite_skill("genomic-coordinates", format="bibtex")


def test_cli_search_and_install_work_in_an_arbitrary_directory(tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["skills", "search", "坐标转换", "--json", "--limit", "5"])
    assert result.exit_code == 0, result.output
    assert "genomic-coordinates" in {x["id"] for x in json.loads(result.output)["skills"]}
    result = runner.invoke(
        app, ["skills", "install", "dnakit", "--agent", "cursor", "--project", str(tmp_path)]
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / ".cursor/skills/dnakit/SKILL.md").is_file()
