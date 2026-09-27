"""Integrity, provenance and standalone-resource checks for the 335-catalog import."""

from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.fixture
def importer(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "others/skills"))
    return importlib.import_module("import_user_catalog_skills")


@pytest.fixture
def imported_copy(tmp_path):
    manifest = ROOT / "registry/user-catalog-skills.json"
    (tmp_path / "registry").mkdir()
    shutil.copy2(manifest, tmp_path / "registry" / manifest.name)
    for item in json.loads(manifest.read_text())["skills"]:
        shutil.copytree(ROOT / "skills" / item["id"], tmp_path / "skills" / item["id"])
    return tmp_path


@pytest.mark.parametrize(
    "field,value",
    [
        ("attribution", "unrelated author"),
        ("upstream_url", "https://example.invalid/unpinned"),
        ("runtime", {"commands": ["unreviewed-command"]}),
        ("data_license", {"license": "MIT"}),
        ("original_format", "invented-native-skill"),
    ],
)
def test_pinned_record_metadata_cannot_drift(importer, imported_copy, field, value):
    path = imported_copy / importer.MANIFEST
    manifest = json.loads(path.read_text())
    manifest["skills"][0][field] = value
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="identity mismatch"):
        importer.check(imported_copy)


def test_resource_drift_and_link_are_rejected(importer, imported_copy, tmp_path):
    path = imported_copy / "skills/onekgpd-individual-queries/assets/kgpe.json"
    original = path.read_bytes()
    path.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        importer.check(imported_copy)
    outside = tmp_path / "outside.json"
    outside.write_bytes(original)
    path.unlink()
    path.symlink_to(outside)
    with pytest.raises(ValueError, match="linked resource"):
        importer.check(imported_copy)


def test_import_refuses_existing_skill_before_network(importer, tmp_path, monkeypatch):
    selected = importer.read_selection()["skills"][0]["id"]
    folder = tmp_path / "skills" / selected
    folder.mkdir(parents=True)
    sentinel = folder / "user-change.txt"
    sentinel.write_text("preserve me")
    monkeypatch.setattr("sys.argv", ["import_user_catalog_skills.py", "--root", str(tmp_path)])
    monkeypatch.setattr(importer, "fetch", lambda *args: pytest.fail("network before preflight"))
    with pytest.raises(FileExistsError, match="existing skill"):
        importer.main()
    assert sentinel.read_text() == "preserve me"
    assert not (tmp_path / importer.MANIFEST).exists()


def run_script(folder, script, *args):
    environment = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [sys.executable, "-B", str(folder / script), *args],
        cwd=folder.parent,
        env=environment,
        text=True,
        capture_output=True,
        timeout=30,
        check=True,
    )


def test_bioconductor_helper_is_standalone(imported_copy):
    folder = imported_copy / "skills/bioconductor-bridge"
    result = run_script(folder, "bioconductor_bridge.py", "--help")
    assert "--workflow" in result.stdout
    assert (folder / "resources/clawbio/common/checksums.py").is_file()


def test_galaxy_offline_catalog_is_bundled(imported_copy):
    folder = imported_copy / "skills/galaxy-bridge"
    result = run_script(folder, "galaxy_bridge.py", "--search", "bwa")
    assert "bwa" in result.stdout.lower()
    catalog = json.loads((folder / "galaxy_catalog.json").read_text())
    assert len(catalog["tools"]) > 1000
    assert len(list((folder / "galaxy_skills").rglob("*.md"))) > 100


def test_onekgpd_known_population_metadata_is_offline(imported_copy, tmp_path):
    folder = imported_copy / "skills/onekgpd-individual-queries"
    output = tmp_path / "population.json"
    run_script(
        folder,
        "scripts/onekgpd_meta.py",
        "sample-metadata",
        "--samples",
        "HG00096",
        "--output",
        str(output),
    )
    data = json.loads(output.read_text())
    text = json.dumps(data)
    assert "HG00096" in text and "GBR" in text


def test_guide_conversions_keep_original_format_and_license(importer):
    manifest = json.loads((ROOT / importer.MANIFEST).read_text())
    entries = {r["id"]: r for r in manifest["skills"]}
    assert entries["galaxy-bioblend"]["original_format"] == "markdown-guide"
    assert entries["kiro-healthomics-genomics"]["original_format"] == "kiro-markdown-and-steering"
    assert entries["kiro-healthomics-genomics"]["license"] == "MIT-0"
    ukb = entries["ukb-ppp-region-fetch"]
    assert "Aviv Madar" in ukb["attribution"]
    assert ukb["data_license"]["license"] == "CC-BY-4.0"
    galaxy = ROOT / "skills/galaxy-bioblend"
    assert (galaxy / "references/original/galaxy-bioblend.md").is_file()
    assert not any("gene_mapping.py" in f["path"] for f in entries["galaxy-bioblend"]["files"])
