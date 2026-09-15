"""Native adaptations survive re-import, preserve provenance and protect edits."""

from __future__ import annotations

import copy
import json
import runpy
import shutil
from pathlib import Path

import pytest

from dnakit.skills import get_skill, install_skill, read_skill_file

ROOT = Path(__file__).resolve().parents[3]


def test_renamed_skill_installs_by_both_aliases_and_retains_provenance(tmp_path: Path) -> None:
    skill = get_skill("sequence-retrieval")
    assert get_skill("tooluniverse-sequence-retrieval")["id"] == skill["id"]
    assert get_skill("public-db-lookup")["id"] == skill["id"]
    assert skill["upstream_name"] == "tooluniverse-sequence-retrieval"
    assert skill["source_id"] == "tooluniverse" and skill["license"] == "Apache-2.0"
    assert skill["runtime"]["python_modules"] == ["dnakit"]
    assert skill["upstream_metadata"]["disable-model-invocation"] is True
    install_skill("tooluniverse-sequence-retrieval", destination=tmp_path)
    target = tmp_path / "sequence-retrieval"
    assert (target / "examples.md").is_file()
    assert (target / "licenses/tooluniverse/LICENSE").is_file()
    assert not (tmp_path / "tooluniverse-sequence-retrieval").exists()
    assert (
        "import dnakit" in read_skill_file("sequence-retrieval", resource="examples.md")["content"]
    )


def test_adaptation_replay_is_idempotent_and_does_not_overwrite_local_edits(tmp_path: Path) -> None:
    module = runpy.run_path(str(ROOT / "others/skills/adaptations.py"))
    manifest, apply = module["manifest"], module["apply_adaptations"]
    rule_path = ROOT / "registry/adaptations.json"
    rule = next(
        rule
        for rule in json.loads(rule_path.read_text())["adaptations"]
        if rule["id"] == "sequence-retrieval"
    )
    (tmp_path / "registry").mkdir()
    (tmp_path / "registry/adaptations.json").write_text(json.dumps({"adaptations": [rule]}))
    shutil.copytree(
        ROOT / "others/skills/adaptations/sequence-retrieval",
        tmp_path / "others/skills/adaptations/sequence-retrieval",
    )
    source = tmp_path / "skills" / rule["original_id"]
    source.mkdir(parents=True)
    (source / "SKILL.md").write_text(
        "---\nname: tooluniverse-sequence-retrieval\n"
        "description: Original fixture\n---\nUpstream instructions\n"
    )
    (source / "examples.md").write_text("Original example fixture\n")
    (source / "LICENSE").write_text("Synthetic preserved notice fixture\n")
    files = manifest(source)
    original_hashes = [{"path": f["path"], "sha256": f["sha256"]} for f in files]
    record = {
        "id": rule["original_id"],
        "source_id": rule["source_id"],
        "upstream_revision": rule["upstream_revision"],
        "upstream_files": original_hashes,
        "upstream_metadata": {"name": rule["original_id"]},
        "aliases": ["public-db-lookup"],
        "modifications": [],
        "files": files,
        "validation": {"execution": "not_run"},
    }
    records, decisions = [record], [{"canonical_id": rule["original_id"]}]
    apply(records, decisions, root=tmp_path)
    assert records[0]["id"] == "sequence-retrieval"
    assert decisions[0]["canonical_id"] == "sequence-retrieval"
    assert records[0]["upstream_files"] == original_hashes
    assert records[0]["upstream_metadata"] == {"name": rule["original_id"]}
    target = tmp_path / "skills/sequence-retrieval"
    assert (target / "LICENSE").read_text() == "Synthetic preserved notice fixture\n"
    before = copy.deepcopy(records)
    apply(records, decisions, root=tmp_path)
    assert records == before
    (target / "examples.md").write_text("User edits")
    with pytest.raises(ValueError, match="Unregistered local skill edits"):
        apply(records, decisions, root=tmp_path)
    assert (target / "examples.md").read_text() == "User edits"


def test_adapted_resources_match_maintainer_templates() -> None:
    source = ROOT / "others/skills/adaptations/sequence-retrieval"
    for path in source.rglob("*"):
        if path.is_file():
            assert (
                path.read_bytes()
                == (ROOT / "skills/sequence-retrieval" / path.relative_to(source)).read_bytes()
            )


def test_incremental_replay_retains_ambiguous_names() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    assert set(catalog["ambiguous_names"]) >= {
        "phylogenetics",
        "variant-annotation",
        "variant-calling",
    }
    assert len(catalog["ambiguous_names"]["variant-annotation"]) == 2
