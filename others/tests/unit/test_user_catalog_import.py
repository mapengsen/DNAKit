"""Guard provenance freshness and additive project/source identity handling."""

from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path

import pytest


@pytest.fixture
def scripts(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "skills"))
    return importlib.import_module("user_catalog_support"), importlib.import_module(
        "build_genomics_survey"
    )


def audit_fixture(root):
    data = {
        "entries": [{"index": 1, "input_stable_id": "tool", "project_id": "tool"}],
        "relationships": [
            {
                "来源": "Tool",
                "关系类型": "uses_component",
                "对象": "Library",
                "证据入口": "https://example.org/README",
            }
        ],
    }
    evidence = {
        "input_sha256": hashlib.sha256(
            json.dumps(data, ensure_ascii=False, sort_keys=True).encode()
        ).hexdigest(),
        "pages": [{"url": "https://example.org/README", "sha256": "a" * 64}],
        "relationships": [
            {
                "index": 1,
                "source": "Tool",
                "relation_type": "uses_component",
                "target": "Library",
                "evidence_url": "https://example.org/README",
                "relationship_verified": False,
            }
        ],
    }
    (root / "registry").mkdir()
    save(root, data, evidence)
    return data, evidence


def save(root, data, evidence):
    (root / "registry/user-catalog-335.json").write_text(json.dumps(data, ensure_ascii=False))
    (root / "registry/user-catalog-evidence.json").write_text(json.dumps(evidence))


def test_changed_input_cannot_reuse_old_source_evidence(tmp_path, scripts):
    support, _ = scripts
    data, evidence = audit_fixture(tmp_path)
    support.load_user_catalog(tmp_path)
    data["entries"][0]["notes"] = "A newly supplied unsupported claim"
    save(tmp_path, data, evidence)
    with pytest.raises(ValueError, match="does not match"):
        support.load_user_catalog(tmp_path)


@pytest.mark.parametrize("collection,key", [("pages", "url"), ("relationships", "index")])
def test_duplicate_evidence_identity_is_rejected(tmp_path, scripts, collection, key):
    support, _ = scripts
    data, evidence = audit_fixture(tmp_path)
    evidence[collection].append(dict(evidence[collection][0]))
    save(tmp_path, data, evidence)
    with pytest.raises(ValueError, match=key):
        support.load_user_catalog(tmp_path)


def test_relationship_cannot_be_attached_to_another_target(tmp_path, scripts):
    support, _ = scripts
    data, evidence = audit_fixture(tmp_path)
    evidence["relationships"][0]["target"] = "Unrelated Library"
    save(tmp_path, data, evidence)
    with pytest.raises(ValueError, match="disagrees"):
        support.load_user_catalog(tmp_path)


def test_context_review_is_bound_to_retrieved_content(tmp_path, scripts):
    support, _ = scripts
    _, evidence = audit_fixture(tmp_path)
    review = {
        **evidence["relationships"][0],
        "evidence_sha256": "a" * 64,
        "relationship_verified": True,
        "review_note": "Reviewed dependency declaration",
    }
    path = tmp_path / "registry/user-catalog-context-review.json"
    path.write_text(json.dumps({"relationships": [review]}))
    _, result = support.load_user_catalog(tmp_path)
    assert result["relationships"][0]["relationship_verified"] is True
    review["evidence_sha256"] = "b" * 64
    path.write_text(json.dumps({"relationships": [review]}))
    with pytest.raises(ValueError, match="source changed"):
        support.load_user_catalog(tmp_path)


def test_missing_relationship_evidence_is_rejected(tmp_path, scripts):
    support, _ = scripts
    data, evidence = audit_fixture(tmp_path)
    evidence["relationships"] = []
    save(tmp_path, data, evidence)
    with pytest.raises(ValueError, match="Missing or extra"):
        support.load_user_catalog(tmp_path)


def test_two_pinned_snapshots_keep_all_skills_under_one_project(tmp_path, monkeypatch, scripts):
    _, build = scripts
    parts = tmp_path / "registry/genomics-survey-parts"
    parts.mkdir(parents=True)
    (parts / "fixture.json").write_text(json.dumps({"queries": [], "candidates": []}))
    catalog = {
        "sources": [
            {"id": "first", "repository": "owner/project", "revision": "1"},
            {"id": "second", "repository": "owner/project", "revision": "2"},
        ],
        "skills": [
            {"id": "one", "source_id": "first", "category": "workflows"},
            {"id": "two", "source_id": "second", "category": "workflows"},
        ],
    }
    (tmp_path / "registry/catalog.json").write_text(json.dumps(catalog))
    integration = tmp_path / "src/dnakit/integrations"
    integration.mkdir(parents=True)
    (integration / "analyses.json").write_text('{"analyses": []}')
    (integration / "mcp_servers.json").write_text('{"servers": []}')
    monkeypatch.setattr(build, "ROOT", tmp_path)
    result = build.build()["projects"]
    assert len(result) == 1
    assert result[0]["dnakit_integration"]["bundled_skills"] == ["one", "two"]
    assert {s["revision"] for s in result[0]["pinned_sources"]} == {"1", "2"}


def test_repository_rename_does_not_merge_different_versioned_projects(scripts):
    _, build = scripts
    assert build.repo_key("https://github.com/weberlab-hhu/Helixer") == build.repo_key(
        "https://github.com/usadellab/Helixer"
    )
    assert build.repo_key("https://github.com/Gaius-Augustus/BRAKER") != build.repo_key(
        "https://github.com/Gaius-Augustus/BRAKER4"
    )
