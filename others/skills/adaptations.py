"""Apply reviewed native DNAKit adaptations after importing pinned originals.

Run this file to refresh only adapted skills in the existing bundle. The full
build_catalog pipeline calls the same function. Local templates live beside this
module; original hashes, metadata, license files and citations remain intact.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def manifest(directory: Path, modes: dict | None = None) -> list[dict]:
    return [
        {
            "path": str(path.relative_to(directory)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
            "mode": (modes or {}).get(str(path.relative_to(directory)), 0o644),
        }
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    ]


def apply_adaptations(records: list[dict], decisions: list[dict], *, root: Path = ROOT) -> None:
    """Replace reviewed resources, retain old IDs as aliases, and refresh hashes."""
    rules = json.loads((root / "registry/adaptations.json").read_text())["adaptations"]
    # Validate the entire input set before changing any directory. An edit in a
    # later skill must not leave earlier skills changed but absent from the catalog.
    for rule in rules:
        matches = [entry for entry in records if entry["id"] in {rule["original_id"], rule["id"]}]
        if len(matches) != 1:
            raise ValueError(f"Adaptation requires one canonical entry: {rule['id']}")
        entry = matches[0]
        directory = root / "skills" / entry["id"]
        modes = {file["path"]: file["mode"] for file in entry["files"]}
        if manifest(directory, modes) != entry["files"]:
            raise ValueError(f"Unregistered local skill edits; review before adapting: {directory}")
    for rule in rules:
        old, new = rule["original_id"], rule["id"]
        matches = [entry for entry in records if entry["id"] in {old, new}]
        if len(matches) != 1:
            raise ValueError(f"Adaptation requires exactly one canonical entry: {old} -> {new}")
        record = matches[0]
        if (
            record["source_id"] != rule["source_id"]
            or record["upstream_revision"] != rule["upstream_revision"]
        ):
            raise ValueError(f"Upstream changed; review native adaptation: {new}")
        source = root / "skills" / record["id"]
        target = root / "skills" / new
        if source != target and target.exists():
            raise ValueError(f"Refusing to replace an existing skill: {target}")
        modes = {entry["path"]: entry["mode"] for entry in record["files"]}
        if manifest(source, modes) != record["files"]:
            raise ValueError(f"Unregistered local skill edits; review before adapting: {source}")
        template = root / "others/skills/adaptations" / new
        replacements = sorted(path for path in template.rglob("*") if path.is_file())
        if not replacements:
            raise ValueError(f"Missing adaptation templates: {new}")
        current_files = {item["path"]: item["sha256"] for item in record["files"]}
        if record["id"] == new and all(
            current_files.get(str(path.relative_to(template)))
            == hashlib.sha256(path.read_bytes()).hexdigest()
            for path in replacements
        ):
            record["runtime"] = rule["runtime"]
            record["adaptation"] = rule
            continue
        updated = copy.deepcopy(record)
        with tempfile.TemporaryDirectory(prefix=".dnakit-adaptation-", dir=source.parent) as raw:
            stage = Path(raw) / new
            shutil.copytree(source, stage)
            replaced = []
            for path in replacements:
                relative = path.relative_to(template)
                if path.is_symlink() or relative.parts[0] == "licenses":
                    raise ValueError("Adaptations must retain the original license resources")
                destination = stage / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(path.read_bytes())
                replaced.append(str(relative))
            metadata = yaml.safe_load((stage / "SKILL.md").read_text().split("---", 2)[1])
            if metadata["name"] != new:
                raise ValueError(f"Adapted metadata name disagrees with registry: {new}")
            updated.update(
                id=new,
                description=metadata["description"],
                aliases=sorted((set(record["aliases"]) | {old}) - {new}),
                runtime=rule["runtime"],
                adaptation=rule,
                modifications=sorted(
                    set(record.get("modifications", [])) | (set(replaced) - {"UPSTREAM.md"})
                ),
                files=manifest(stage, modes),
            )
            if updated["files"] != record["files"]:
                updated["validation"] = {
                    "metadata": "passed",
                    "resources": "checksummed",
                    "execution": "not_run",
                }
            backup = Path(raw) / "previous"
            source.rename(backup)
            try:
                stage.rename(target)
            except BaseException:
                backup.rename(source)
                raise
        record.clear()
        record.update(updated)
        for decision in decisions:
            if decision.get("canonical_id") == old:
                decision["canonical_id"] = new
    renamed = {rule["original_id"]: rule["id"] for rule in rules}
    for record in records:
        for dependency in record.get("required_skills", []):
            if dependency.get("id") in renamed:
                dependency["id"] = renamed[dependency["id"]]


def write_indexes(catalog: dict, *, root: Path = ROOT) -> None:
    """Refresh derived public identifiers; provenance paths stay upstream paths."""
    with (root / "registry/skills.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "id",
                "category",
                "capability",
                "source_id",
                "license",
                "upstream_revision",
                "upstream_url",
            ],
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(catalog["skills"])
    path = root / "skills/upstream-manifest.json"
    value = json.loads(path.read_text())
    value["skills"] = [
        {"id": entry["id"], "source": entry["source_id"], "revision": entry["upstream_revision"]}
        for entry in catalog["skills"]
    ]
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main() -> None:
    from build_catalog import check, disambiguate_aliases
    from write_inventory import main as write_inventory
    from write_migration import main as write_migration

    catalog_path, decisions_path = ROOT / "registry/catalog.json", ROOT / "registry/decisions.json"
    catalog, decisions = (
        json.loads(catalog_path.read_text()),
        json.loads(decisions_path.read_text()),
    )
    apply_adaptations(catalog["skills"], decisions["decisions"])
    catalog["skills"].sort(key=lambda entry: entry["id"])
    catalog["ambiguous_names"] = disambiguate_aliases(catalog["skills"])
    for path, value in ((catalog_path, catalog), (decisions_path, decisions)):
        path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    write_indexes(catalog)
    write_migration()
    write_inventory()
    check()


if __name__ == "__main__":
    main()
