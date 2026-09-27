"""Merge explicitly reviewed genomics extensions while preserving baseline entries.

The upstream import manifest and native skill specification are authoritative.
This helper is also called by build_catalog.py so a full rebuild retains them.
No scientific workflows, dependency installers or network requests are executed.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from import_genomics_skills import DATE, SOURCES, metadata, progress
from import_genomics_skills import check as check_imports
from import_user_catalog_skills import INTEGRATION as USER_INTEGRATION
from import_user_catalog_skills import check as check_user_imports

ROOT = Path(__file__).resolve().parents[2]
UPSTREAM_MANIFEST = ROOT / "registry/genomics-new-skills.json"
NATIVE_MANIFEST = ROOT / "registry/genomics-native-skills.json"
USER_MANIFEST = ROOT / "registry/user-catalog-skills.json"
INTEGRATION = "genomics-expansion"
CATEGORIES = {
    "breeding": {"zh": "基因组育种与农业", "en": "Genomic breeding and agriculture"},
    "dna-storage": {"zh": "DNA 数据存储", "en": "DNA data storage"},
    "nanotechnology": {"zh": "DNA 纳米技术", "en": "DNA nanotechnology"},
}
ALIASES = {
    "基因组育种": "genomic breeding",
    "基因组选择": "genomic selection",
    "育种": "breeding",
    "亲本选择": "parental selection",
    "种质": "germplasm",
    "DNA存储": "DNA data storage",
    "dna存储": "DNA data storage",
    "数据存储": "DNA data storage",
    "DNA纳米": "DNA nanotechnology",
    "dna纳米": "DNA nanotechnology",
    "纳米技术": "DNA nanotechnology",
    "纳米": "DNA nanotechnology",
    "法医": "forensic",
    "亲子鉴定": "kinship",
}


def read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def managed_ids() -> set[str]:
    """Only explicitly registered extensions are exempt from baseline cleanup."""
    result: set[str] = set()
    for path in (UPSTREAM_MANIFEST, NATIVE_MANIFEST, USER_MANIFEST):
        if path.exists():
            result.update(item["id"] for item in read(path)["skills"])
    return result


def is_managed_source(source: dict[str, Any]) -> bool:
    ids = {item["id"] for item in SOURCES}
    if USER_MANIFEST.exists():
        ids.update(item["id"] for item in read(USER_MANIFEST)["sources"])
    return source["id"] in ids


def imported_manifests() -> list[tuple[Path, str]]:
    result = [(UPSTREAM_MANIFEST, INTEGRATION)]
    if USER_MANIFEST.exists():
        result.append((USER_MANIFEST, USER_INTEGRATION))
    return result


def file_manifest(folder: Path) -> list[dict[str, Any]]:
    if folder.is_symlink():
        raise ValueError(f"Unexpected symlink skill directory: {folder}")
    result = []
    for path in sorted(folder.rglob("*")):
        if not path.is_file():
            continue
        if (
            path.is_symlink()
            or "__pycache__" in path.parts
            or not path.resolve().is_relative_to(folder.resolve())
        ):
            raise ValueError(f"Unexpected generated/link resource: {path}")
        blob = path.read_bytes()
        result.append(
            {
                "path": str(path.relative_to(folder)),
                "sha256": hashlib.sha256(blob).hexdigest(),
                "bytes": len(blob),
                "mode": path.stat().st_mode & 0o777,
            }
        )
    return result


def native_records(*, write_provenance: bool = True) -> list[dict[str, Any]]:
    records = []
    seen = set()
    for spec in read(NATIVE_MANIFEST)["skills"]:
        skill_id = spec["id"]
        if not re.fullmatch(r"dnakit-[a-z0-9-]{1,56}", skill_id) or skill_id in seen:
            raise ValueError(f"Invalid or duplicate native skill ID: {skill_id}")
        seen.add(skill_id)
        folder = ROOT / "skills" / skill_id
        for required in ("SKILL.md", "DNAKIT.md", "LICENSE"):
            if not (folder / required).is_file():
                raise FileNotFoundError(f"Missing native skill resource: {folder / required}")
        meta, _ = metadata((folder / "SKILL.md").read_bytes())
        if meta["name"] != skill_id:
            raise ValueError(f"Native skill ID mismatch: {skill_id}")
        functions = spec["native_functions"]
        if any(not re.fullmatch(r"dnakit\.(genomics|integrations)\.[a-z_]+", f) for f in functions):
            raise ValueError(f"Unexpected native function: {skill_id}")
        notice = (
            f"# DNAKit native workflow provenance\n\nSkill: `{skill_id}`. Date: {DATE}.\n\n"
            "This workflow and its Python adapter are maintained by DNAKit contributors under "
            "the repository MIT license. It is not an imported upstream SKILL.md.\n\n"
            "The underlying scientific methods, software, data and model authors remain credited "
            "in SKILL.md and DNAKIT.md. The skill does not confer a license to external packages "
            "or establish a full scientific validation result.\n\n"
            "Public functions:\n\n" + "\n".join(f"- `{name}`" for name in functions) + "\n"
        )
        provenance = folder / "UPSTREAM.md"
        if write_provenance:
            provenance.write_text(notice, encoding="utf-8")
        elif not provenance.exists() or provenance.read_text() != notice:
            raise ValueError(f"Native provenance needs synchronization: {skill_id}")
        citations = sorted(
            set(
                re.findall(
                    r"https://doi\.org/[^\s)>,`]+",
                    (folder / "SKILL.md").read_text() + (folder / "DNAKIT.md").read_text(),
                )
            )
        )
        records.append(
            {
                "id": skill_id,
                "role": "dna_workflow",
                "description": meta["description"],
                "summary_zh": spec["summary_zh"],
                "category": spec["category"],
                "capability": skill_id,
                "source_id": "dnakit",
                "license": "MIT",
                "upstream_name": skill_id,
                "upstream_url": f"https://github.com/mapengsen/DNAKit/tree/main/skills/{skill_id}",
                "upstream_revision": "dnakit-genomics-" + DATE,
                "attribution": (
                    "DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md"
                ),
                "citation_guidance": (
                    "Cite the installed DNAKit version and the underlying methods/software "
                    "actually used; see the skill references."
                ),
                "citations": citations,
                "native_functions": functions,
                "aliases": [],
                "alternatives": [],
                "keywords": spec["keywords"],
                "runtime": {
                    "python_modules": ["dnakit", *spec["python_modules"]],
                    "commands": [],
                    "environment": [],
                    "manual_checks": [
                        "See DNAKIT.md for optional analysis-specific dependencies "
                        "and validation limits."
                    ],
                },
                "required_skills": [],
                "validation": {
                    "metadata": "passed",
                    "resources": "checksummed",
                    "execution": "not_run",
                },
                "adaptation": {"implementation": "native-python-api", "original_id": skill_id},
                "modifications": [],
                "files": file_manifest(folder),
            }
        )
    return records


def expansion_records(*, write_provenance: bool = True) -> list[dict[str, Any]]:
    check_imports(ROOT)
    imported = read(UPSTREAM_MANIFEST)["skills"]
    additions = []
    if USER_MANIFEST.exists():
        check_user_imports(ROOT)
        additions = read(USER_MANIFEST)["skills"]
    return [*imported, *native_records(write_provenance=write_provenance), *additions]


def check_record_identities(records: list[dict[str, Any]], additions: list[dict[str, Any]]) -> None:
    managed = {item["id"]: item for item in additions}
    for item in records:
        if item["id"] in managed:
            incoming = managed[item["id"]]
            for key in ("source_id", "upstream_revision", "upstream_path", "native_functions"):
                if item.get(key) != incoming.get(key):
                    raise ValueError(f"Managed skill identity collision: {item['id']}/{key}")


def preflight(records: list[dict[str, Any]]) -> None:
    """Reject missing/drifted extensions before the baseline build replaces any files."""
    check_record_identities(records, expansion_records(write_provenance=False))


def merge_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    additions = expansion_records()
    check_record_identities(records, additions)
    managed = {item["id"] for item in additions}
    result = [entry for entry in records if entry["id"] not in managed] + additions
    ids = [entry["id"] for entry in result]
    capabilities = [entry["capability"] for entry in result]
    if len(ids) != len(set(ids)) or len(capabilities) != len(set(capabilities)):
        raise ValueError("Duplicate catalog ID/capability after genomics integration")
    return sorted(result, key=lambda item: item["id"])


def extend_metadata(catalog: dict[str, Any]) -> None:
    catalog["categories"].update(CATEGORIES)
    catalog["query_aliases"].update(ALIASES)
    catalog["snapshot_date"] = max(catalog["snapshot_date"], DATE)
    catalog["extension_manifests"] = [
        "registry/genomics-new-skills.json",
        "registry/genomics-native-skills.json",
    ]
    imported = len(read(UPSTREAM_MANIFEST)["skills"])
    native = len(read(NATIVE_MANIFEST)["skills"])
    catalog["selection"][INTEGRATION] = imported
    catalog["extension_summary"] = {
        "imported": imported,
        "native": native,
        "total": imported + native,
    }
    if USER_MANIFEST.exists():
        count = len(read(USER_MANIFEST)["skills"])
        catalog["extension_manifests"].append("registry/user-catalog-skills.json")
        catalog["selection"][USER_INTEGRATION] = count
        catalog["user_catalog_summary"] = {
            "imported": count,
            "sources": len(read(USER_MANIFEST)["sources"]),
            "selection": "registry/user-catalog-skill-selection.json",
            "review": "registry/user-catalog-skills-review.json",
        }


def sync_sources(catalog: dict[str, Any]) -> None:
    current = read(ROOT / "registry/sources.json")
    by_id = {item["id"]: item for item in current["sources"]}
    for manifest, integration in imported_manifests():
        for source in read(manifest)["sources"]:
            existing = by_id.get(source["id"])
            if existing and any(
                existing.get(key) != source[key]
                for key in ("repository", "revision", "content_license")
            ):
                raise ValueError(f"Source identity collision: {source['id']}")
            by_id[source["id"]] = {**source, "catalog_integration": integration}
    current.update(snapshot_date=DATE, sources=list(by_id.values()))
    write(ROOT / "registry/sources.json", current)
    catalog["sources"] = [{"id": "dnakit", "repository": "mapengsen/DNAKit"}, *current["sources"]]


def check_sources(catalog: dict[str, Any]) -> None:
    """The discovery catalog and source registry must agree with the pinned importer."""
    expected = {
        item["id"]: item
        for manifest, _ in imported_manifests()
        for item in read(manifest)["sources"]
    }
    for label, records in (
        ("catalog", catalog["sources"]),
        ("sources", read(ROOT / "registry/sources.json")["sources"]),
    ):
        by_id = {item["id"]: item for item in records}
        if len(by_id) != len(records):
            raise ValueError(f"Duplicate source ID in {label}")
        for source_id, source in expected.items():
            actual = by_id.get(source_id)
            if actual is None:
                raise ValueError(f"Missing pinned source in {label}: {source_id}")
            for key in source:
                if actual.get(key) != source[key]:
                    raise ValueError(f"Source metadata mismatch in {label}: {source_id}/{key}")


def merge_decisions(decisions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep reviewed baseline paths and record the explicitly imported extension paths."""
    result = [
        item for item in decisions if item.get("integration") not in {INTEGRATION, USER_INTEGRATION}
    ]
    for manifest, integration in imported_manifests():
        for item in read(manifest)["skills"]:
            converted = item.get("original_format", "SKILL.md") != "SKILL.md"
            result.append(
                {
                    "source": item["source_id"],
                    "path": item["upstream_path"]
                    if converted
                    else item["upstream_path"] + "/SKILL.md",
                    "name": item["id"],
                    "status": integration,
                    "integration": integration,
                    "reason": (
                        "Transparent licensed Markdown-guide conversion; "
                        "see selected original files. "
                        if converted
                        else "DNA task with retained license; commit-pinned complete skill import. "
                    )
                    + f"See registry/{manifest.name} for portability edits and validation limits.",
                }
            )
    return result


def notices() -> None:
    manifest = read(UPSTREAM_MANIFEST)
    section = [
        "## Genomics expansion sources",
        "",
        "| Source | License | Retained notices |",
        "| --- | --- | --- |",
    ]
    for source in manifest["sources"]:
        representative = next(
            item for item in manifest["skills"] if item["source_id"] == source["id"]
        )
        folder = ROOT / "skills" / representative["id"]
        destination = ROOT / "licenses/skills" / source["id"]
        destination.mkdir(parents=True, exist_ok=True)
        for name in ("LICENSE", "NOTICE", "CITATION.cff"):
            if (folder / name).is_file():
                (destination / name).write_bytes((folder / name).read_bytes())
        url = f"https://github.com/{source['repository']}/tree/{source['revision']}"
        section.append(
            f"| [{source['repository']}]({url}) | {source['content_license']} | "
            f"[source notices](licenses/skills/{source['id']}/) and skill-root notices |"
        )
    section += [
        "",
        "These skills retain the original author notices and pinned file hashes. "
        "DNAKit's added execution notes, frontmatter normalization and portable links are "
        "recorded in the extension manifest. Scientific backends and model/data/service "
        "licenses remain separate; skill inclusion is not execution validation.",
        "",
    ]
    start = "<!-- BEGIN DNAKit genomics expansion -->"
    end = "<!-- END DNAKit genomics expansion -->"
    path = ROOT / "THIRD_PARTY_NOTICES.md"
    text = path.read_text()
    block = start + "\n" + "\n".join(section) + "\n" + end
    if start in text:
        a, b = text.index(start), text.index(end) + len(end)
        text = text[:a] + block + text[b:]
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text)
    user_notices()


def user_notices() -> None:
    """Retain individual skill licenses; do not flatten per-leaf or data terms."""
    if not USER_MANIFEST.exists():
        return
    manifest = read(USER_MANIFEST)
    lines = [
        "## User catalog skill additions",
        "",
        "| Skill / original author | License | Notices |",
        "| --- | --- | --- |",
    ]
    for item in manifest["skills"]:
        folder = ROOT / "skills" / item["id"]
        destination = ROOT / "licenses/skills" / item["source_id"] / item["id"]
        destination.mkdir(parents=True, exist_ok=True)
        for name in (
            "LICENSE",
            "LICENSE.txt",
            "SKILL_LICENSES.md",
            "THIRD-PARTY.txt",
            "UPSTREAM.md",
            "references/upstream-readme.md",
            "bundled_slices/README.md",
        ):
            if (folder / name).is_file():
                target = destination / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes((folder / name).read_bytes())
        lines.append(
            f"| [{item['id']}]({item['upstream_url']}) / {item['attribution']} | "
            f"{item['license']} | "
            f"[retained notices](licenses/skills/{item['source_id']}/{item['id']}/) |"
        )
    lines += [
        "",
        "The 335-project catalog additions retain original resource hashes and distinct "
        "software, text and bundled-data licenses. Converted Kiro/Galaxy guides are "
        "explicitly identified; remote services and scientific backends are "
        "installed/configured separately. See registry/user-catalog-skills-review.json "
        "for selection and validation limits.",
        "",
    ]
    start = "<!-- BEGIN DNAKit user catalog skills -->"
    end = "<!-- END DNAKit user catalog skills -->"
    path = ROOT / "THIRD_PARTY_NOTICES.md"
    text = path.read_text()
    block = start + "\n" + "\n".join(lines) + "\n" + end
    if start in text:
        a, b = text.index(start), text.index(end) + len(end)
        text = text[:a] + block + text[b:]
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    path.write_text(text)


def write_auxiliary(catalog: dict[str, Any]) -> None:
    import csv

    from write_migration import main as write_migration

    decisions = read(ROOT / "registry/decisions.json")
    decisions.update(
        snapshot_date=catalog["snapshot_date"], decisions=merge_decisions(decisions["decisions"])
    )
    write(ROOT / "registry/decisions.json", decisions)

    with (ROOT / "registry/skills.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "id",
            "category",
            "capability",
            "source_id",
            "license",
            "upstream_revision",
            "upstream_url",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(catalog["skills"])
    write(
        ROOT / "skills/upstream-manifest.json",
        {
            "schema_version": 2,
            "snapshot_date": catalog["snapshot_date"],
            "catalog": "../registry/catalog.json",
            "decisions": "../registry/decisions.json",
            "skills": [
                {
                    "id": item["id"],
                    "source": item["source_id"],
                    "revision": item["upstream_revision"],
                }
                for item in catalog["skills"]
            ],
        },
    )
    write_migration()
    notices()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Verify the merged catalog without writing."
    )
    args = parser.parse_args()
    from build_catalog import check as check_catalog
    from build_catalog import disambiguate_aliases

    catalog = read(ROOT / "registry/catalog.json")
    if args.check:
        check_sources(catalog)
        expected = {item["id"]: item for item in expansion_records(write_provenance=False)}
        found = {item["id"]: item for item in catalog["skills"]}
        for skill_id, item in expected.items():
            if skill_id not in found or found[skill_id]["files"] != item["files"]:
                raise ValueError(f"Catalog extension needs synchronization: {skill_id}")
            for key, value in item.items():
                if (
                    key not in {"aliases", "keywords", "validation"}
                    and found[skill_id].get(key) != value
                ):
                    raise ValueError(f"Catalog extension metadata mismatch: {skill_id}/{key}")
        check_catalog()
        return 0
    baseline = {
        item["id"]: copy.deepcopy(item)
        for item in catalog["skills"]
        if item["id"] not in managed_ids()
    }
    catalog["skills"] = merge_records(catalog["skills"])
    progress(1, 3, "merged imported and native skill records")
    sync_sources(catalog)
    extend_metadata(catalog)
    catalog["ambiguous_names"] = disambiguate_aliases(catalog["skills"])
    # Disambiguation may adjust aliases/keywords; scientific baseline fields must survive.
    for item in catalog["skills"]:
        if item["id"] in baseline:
            for key in ("files", "upstream_metadata", "adaptation", "runtime", "upstream_revision"):
                if item.get(key) != baseline[item["id"]].get(key):
                    raise ValueError(f"Baseline scientific metadata changed: {item['id']}/{key}")
    write(ROOT / "registry/catalog.json", catalog)
    progress(2, 3, "synchronized sources and categories")
    write_auxiliary(catalog)
    check_catalog()
    progress(3, 3, "validated catalog, licenses and migration")
    print(
        json.dumps(
            {
                "skills": len(catalog["skills"]),
                "extensions": len(managed_ids()),
                "sources": dict(Counter(item["source_id"] for item in catalog["skills"])),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
