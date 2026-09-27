"""Regression checks for additive imports and reproducible catalog rebuilds."""

from __future__ import annotations

import importlib
import json
from pathlib import Path

import pytest


@pytest.fixture
def catalog_scripts(monkeypatch):
    directory = Path(__file__).resolve().parents[2] / "skills"
    monkeypatch.syspath_prepend(str(directory))
    return (
        importlib.import_module("build_catalog"),
        importlib.import_module("sync_genomics_expansion"),
    )


def test_baseline_cleanup_preserves_registered_extensions(tmp_path, monkeypatch, catalog_scripts):
    build, _ = catalog_scripts
    registry = tmp_path / "registry"
    registry.mkdir()
    names = [
        "old-baseline",
        "gigwa-region-scan",
        "dnakit-genomic-breeding",
        "dnakit",
        "dnakit-plan",
    ]
    (registry / "catalog.json").write_text(json.dumps({"skills": [{"id": n} for n in names]}))
    for name in names:
        folder = tmp_path / "skills" / name
        folder.mkdir(parents=True)
        (folder / "sentinel").write_text(name)
    monkeypatch.setattr(build, "ROOT", tmp_path)
    monkeypatch.setattr(build, "REGISTRY", registry)
    monkeypatch.setattr(build, "WORK", tmp_path / "work")
    assert build.import_selected([], []) == []
    assert not (tmp_path / "skills/old-baseline").exists()
    for name in names[1:]:
        assert (tmp_path / "skills" / name / "sentinel").read_text() == name


def test_source_identity_collision_is_rejected(tmp_path, monkeypatch, catalog_scripts):
    _, sync = catalog_scripts
    (tmp_path / "registry").mkdir()
    source = {
        "id": "gigwa-mcp",
        "repository": "someone/unrelated",
        "revision": "wrong",
        "content_license": "MIT",
    }
    (tmp_path / "registry/sources.json").write_text(json.dumps({"sources": [source]}))
    monkeypatch.setattr(sync, "ROOT", tmp_path)
    with pytest.raises(ValueError, match="Source identity collision"):
        sync.sync_sources({})
    assert json.loads((tmp_path / "registry/sources.json").read_text())["sources"] == [source]


def test_extension_merge_preserves_baseline_and_is_idempotent(monkeypatch, catalog_scripts):
    _, sync = catalog_scripts
    old = {"id": "existing", "capability": "existing", "upstream_metadata": {"author": "A"}}
    added = {"id": "added", "capability": "added", "runtime": {"commands": ["tool"]}}
    monkeypatch.setattr(sync, "expansion_records", lambda: [added])
    merged = sync.merge_records([old])
    assert next(item for item in merged if item["id"] == "existing") == old
    assert sync.merge_records(merged) == merged
    with pytest.raises(ValueError, match="Duplicate catalog"):
        sync.merge_records([old, {"id": "collision", "capability": "added"}])
    with pytest.raises(ValueError, match="Managed skill identity collision"):
        sync.merge_records([{**added, "source_id": "unrelated-upstream"}])


def test_full_build_replays_extensions_after_baseline_adaptations(
    tmp_path, monkeypatch, catalog_scripts
):
    build, sync = catalog_scripts
    registry = tmp_path / "registry"
    registry.mkdir()
    sources = [{"id": "legacy"}, {"id": "gigwa-mcp"}]
    (registry / "sources.json").write_text(json.dumps({"sources": sources}))
    baseline = {
        "id": "existing",
        "source_id": "legacy",
        "category": "workflows",
        "capability": "existing",
        "license": "MIT",
        "upstream_revision": "abc",
        "upstream_url": "https://example.test",
        "aliases": [],
        "keywords": [],
    }
    native = {**baseline, "id": "dnakit", "capability": "entry"}
    planning = {**baseline, "id": "dnakit-plan", "capability": "planning"}
    added = {**baseline, "id": "added", "capability": "added"}
    events = []

    def scan(filtered, fetch):
        assert filtered == [{"id": "legacy"}]
        return [], []

    def adapt(records, decisions):
        assert {row["id"] for row in records} == {"existing", "dnakit", "dnakit-plan"}
        events.append("adapt")

    def merge(records):
        assert events == ["preflight", "adapt"]
        events.append("merge")
        return [*records, added]

    monkeypatch.setattr(build, "ROOT", tmp_path)
    monkeypatch.setattr(build, "REGISTRY", registry)
    monkeypatch.setattr(build, "WORK", tmp_path / "work")
    monkeypatch.setattr(build, "scan", scan)
    monkeypatch.setattr(build, "preflight", lambda records: events.append("preflight"))
    monkeypatch.setattr(build, "import_selected", lambda sources, selected: [baseline])
    monkeypatch.setattr(build, "entry_record", lambda: native)
    monkeypatch.setattr(build, "plan_entry_record", lambda: planning)
    monkeypatch.setattr(build, "apply_adaptations", adapt)
    monkeypatch.setattr(build, "merge_records", merge)
    monkeypatch.setattr(build, "merge_decisions", lambda decisions: decisions)
    monkeypatch.setattr(build, "extend_metadata", lambda catalog: None)
    monkeypatch.setattr(build, "write_expansion_notices", lambda: None)
    monkeypatch.setattr(build, "check", lambda: None)
    migration = importlib.import_module("write_migration")
    monkeypatch.setattr(migration, "main", lambda: None)
    monkeypatch.setattr("sys.argv", ["build_catalog.py"])
    (tmp_path / "skills").mkdir()
    build.main()
    result = json.loads((registry / "catalog.json").read_text())
    assert {item["id"] for item in result["skills"]} == {
        "existing",
        "dnakit",
        "dnakit-plan",
        "added",
    }
    assert events == ["preflight", "adapt", "merge"]
    assert sync.is_managed_source({"id": "gigwa-mcp"})


def test_failed_extension_preflight_never_starts_baseline_import(
    tmp_path, monkeypatch, catalog_scripts
):
    build, _ = catalog_scripts
    registry = tmp_path / "registry"
    registry.mkdir()
    (registry / "sources.json").write_text('{"sources": []}')
    original = '{"skills": [{"id": "baseline"}]}'
    (registry / "catalog.json").write_text(original)
    invoked = []

    def reject(records):
        assert records == [{"id": "baseline"}]
        raise ValueError("Resource checksum mismatch")

    monkeypatch.setattr(build, "REGISTRY", registry)
    monkeypatch.setattr(build, "WORK", tmp_path / "work")
    monkeypatch.setattr(build, "scan", lambda sources, fetch: ([], []))
    monkeypatch.setattr(build, "preflight", reject)
    monkeypatch.setattr(build, "import_selected", lambda *args: invoked.append("import"))
    monkeypatch.setattr("sys.argv", ["build_catalog.py"])
    with pytest.raises(ValueError, match="Resource checksum mismatch"):
        build.main()
    assert invoked == []
    assert (registry / "catalog.json").read_text() == original


def test_import_check_rejects_source_license_drift(tmp_path, catalog_scripts):
    importer = importlib.import_module("import_genomics_skills")
    registry = tmp_path / "registry"
    registry.mkdir()
    sources = [dict(source) for source in importer.SOURCES]
    sources[0]["content_license"] = "unrelated-license"
    (registry / "genomics-new-skills.json").write_text(
        json.dumps({"sources": sources, "skills": []})
    )
    with pytest.raises(ValueError, match="Pinned source metadata mismatch"):
        importer.check(tmp_path)
