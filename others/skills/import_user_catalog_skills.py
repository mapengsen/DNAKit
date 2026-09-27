"""Import reviewed skill increments from the user-supplied 335-project catalog.

The explicit selection pins revisions, archive/file hashes, license scope and
portability edits. --check is offline and read-only. Import never overwrites an
existing skill, runs scientific code, installs dependencies or contacts a service
other than the fixed upstream GitHub download URLs.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import posixpath
import re
import sys
import tarfile
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any

import yaml
from import_genomics_skills import digest, metadata, progress

ROOT = Path(__file__).resolve().parents[2]
SELECTION = ROOT / "registry/user-catalog-skill-selection.json"
MANIFEST = Path("registry/user-catalog-skills.json")
INTEGRATION = "user-catalog-expansion"
FileMap = dict[str, tuple[bytes, int]]


def read_selection() -> dict[str, Any]:
    data = json.loads(SELECTION.read_text())
    source_ids = [s["id"] for s in data["sources"]]
    ids = [s["id"] for s in data["skills"]]
    if len(source_ids) != len(set(source_ids)) or len(ids) != len(set(ids)):
        raise ValueError("Duplicate selected skill/source identity")
    for source in data["sources"]:
        if not re.fullmatch(r"[0-9a-f]{40}", source["revision"]):
            raise ValueError("Source revision must be a complete Git commit")
    for spec in data["skills"]:
        if not re.fullmatch(r"[a-z0-9-]{1,63}", spec["id"]):
            raise ValueError("Invalid selected skill ID")
        if spec["source_id"] not in source_ids:
            raise ValueError("Unknown selected source")
    return data


def safe_path(value: str) -> str:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or ".." in path.parts or "\\" in value:
        raise ValueError(f"Unsafe relative resource path: {value}")
    if str(path) != value or value == ".":
        raise ValueError(f"Noncanonical resource path: {value}")
    return value


def fetch(url: str, limit: int = 120_000_000) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "DNAKit-skills-import"})
    with urllib.request.urlopen(request, timeout=90) as response:
        blob = response.read(limit + 1)
    if len(blob) > limit:
        raise ValueError(f"Source exceeds reviewed download limit: {url}")
    return blob


def source_files(source: dict[str, Any], cache: Path) -> FileMap:
    """Read complete pinned archives, or an explicit complete leaf file list."""
    repo, revision = source["repository"], source["revision"]
    if source.get("acquisition") == "pinned-raw-files":
        result = {}
        for item in source["raw_files"]:
            name = safe_path(item["path"])
            path = cache / source["raw_cache"] / name
            blob = (
                path.read_bytes()
                if path.exists()
                else fetch(
                    f"https://raw.githubusercontent.com/{repo}/{revision}/{name}", 10_000_000
                )
            )
            git_hash = hashlib.sha1(b"blob " + str(len(blob)).encode() + b"\0" + blob).hexdigest()
            if digest(blob) != item["sha256"] or git_hash != item["git_blob_sha1"]:
                raise ValueError(f"Pinned raw-file hash mismatch: {name}")
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as handle:
                    handle.write(blob)
            if name in result:
                raise ValueError(f"Duplicate upstream resource: {name}")
            result[name] = blob, item["mode"]
        return result
    archive = cache / (repo.replace("/", "__") + ".tar.gz")
    blob = (
        archive.read_bytes()
        if archive.exists()
        else fetch(f"https://codeload.github.com/{repo}/tar.gz/{revision}")
    )
    if digest(blob) != source["archive_sha256"]:
        raise ValueError(f"Pinned archive hash mismatch: {repo}")
    if not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        with archive.open("xb") as handle:
            handle.write(blob)
    result = {}
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if member.isdir():
                continue
            parts = PurePosixPath(member.name).parts
            if not parts or not parts[0].endswith(revision):
                raise ValueError(f"Unpinned archive member: {member.name}")
            name = safe_path(str(PurePosixPath(*parts[1:])))
            # Selected resources are checked for completeness below; links are never extracted.
            if not member.isfile():
                for spec in read_selection()["skills"]:
                    if spec["source_id"] != source["id"]:
                        continue
                    required = set(spec.get("include_paths", [])) | set(
                        spec.get("extra_files", {}).values()
                    )
                    in_leaf = "include_paths" not in spec and name.startswith(
                        spec["upstream_path"] + "/"
                    )
                    if in_leaf or name in required or name == spec["license_path"]:
                        raise ValueError(f"Selected source contains a nonregular resource: {name}")
                continue
            if name in result or member.size > 30_000_000:
                raise ValueError(f"Duplicate or oversized upstream file: {name}")
            handle = tar.extractfile(member)
            assert handle is not None
            result[name] = handle.read(), 0o755 if member.mode & 0o111 else 0o644
    return result


def make_skill(
    spec: dict[str, Any], source: dict[str, Any], upstream: FileMap
) -> tuple[dict, FileMap]:
    skill_id, folder = spec["id"], spec["upstream_path"]
    files: FileMap = {}
    paths = {}

    def retain(destination: str, original: str) -> None:
        safe_path(destination)
        safe_path(original)
        if destination in files and paths[destination] != original:
            raise ValueError(f"Conflicting selected resource: {skill_id}/{destination}")
        files[destination] = upstream[original]
        paths[destination] = original

    if "include_paths" in spec:
        for name in spec["include_paths"]:
            retain("references/original/" + name, name)
    else:
        for name in sorted(upstream):
            if name.startswith(folder + "/"):
                retain(name[len(folder) + 1 :], name)
        if "SKILL.md" not in files:
            raise ValueError(f"Missing upstream SKILL.md: {skill_id}")
    for destination, original in spec.get("extra_files", {}).items():
        retain(destination, original)
    # An individual skill license takes precedence over a repository license.
    retain("LICENSE", spec["license_path"])
    if "README.md" in upstream:
        retain("references/upstream-readme.md", "README.md")
    originals = [
        {"path": name, "upstream_path": paths[name], "sha256": digest(blob), "bytes": len(blob)}
        for name, (blob, _) in sorted(files.items())
    ]
    if "include_paths" in spec:
        meta = {"name": skill_id, "description": spec["description"]}
        body = (
            f"\n# {spec['generated_title']}\n\n"
            "This entrypoint adapts the following original licensed Markdown guides. "
            "They were not originally packaged as SKILL.md. Follow DNAKIT.md for actual "
            "execution paths, then read the appropriate complete guide.\n\n"
            + "\n".join(f"- [{p}]({p})" for p in files if p.startswith("references/original/"))
            + "\n"
        ).encode()
    else:
        meta, end = metadata(files["SKILL.md"][0])
        body = files["SKILL.md"][0][end:]
    local_meta = dict(meta)
    local_meta["name"] = skill_id
    local_meta["metadata"] = dict(meta.get("metadata", {}))
    allowed = {"name", "description", "license", "metadata", "allowed-tools"}
    for key in list(local_meta):
        if key not in allowed:
            local_meta["metadata"]["upstream-" + key] = local_meta.pop(key)
    description = str(local_meta["description"]).strip().replace("<", "(").replace(">", ")")
    if len(description) > 1024:
        raise ValueError(f"Review a concise description for {skill_id}")
    local_meta["description"] = description
    intro = (
        b"\n## DNAKit execution\n\nRead [DNAKIT.md](DNAKIT.md) "
        b"for runtime, resources and validation limits.\n"
    )
    files["SKILL.md"] = (
        b"---\n"
        + yaml.safe_dump(local_meta, sort_keys=False, allow_unicode=True).encode()
        + b"---\n"
        + intro
        + body,
        0o644,
    )
    modifications = ["SKILL.md"]
    for edit in spec.get("replacements", []):
        blob, mode = files[edit["path"]]
        old, new = edit["old"].encode(), edit["new"].encode()
        if blob.count(old) != edit["count"]:
            raise ValueError(f"Portability edit no longer matches: {skill_id}/{edit['path']}")
        files[edit["path"]] = blob.replace(old, new), mode
        modifications.append(edit["path"])
    # Keep complete local references; point out-of-folder references at the pinned commit.
    for name, original_path in paths.items():
        if not name.endswith(".md"):
            continue
        blob, mode = files[name]

        def link(
            match: re.Match[str], name: str = name, original_path: str = original_path
        ) -> str:
            target = match[1]
            if target.startswith(("https://", "http://", "mailto:", "#")):
                return match[0]
            path, sep, fragment = target.partition("#")
            local = posixpath.normpath(posixpath.join(posixpath.dirname(name), path))
            if (
                local in {"DNAKIT.md", "UPSTREAM.md"}
                or local in files
                or any(key.startswith(local.rstrip("/") + "/") for key in files)
            ):
                return match[0]
            remote = posixpath.normpath(posixpath.join(posixpath.dirname(original_path), path))
            return (
                "]("
                + f"https://github.com/{source['repository']}/blob/{source['revision']}/{remote}"
                + (sep + fragment if sep else "")
                + ")"
            )

        revised = re.sub(r"\]\(([^\s)]+)\)", link, blob.decode())
        if name == "SKILL.md" and spec["source_id"] == "user335-deepmind":
            revised = revised.replace(
                "Read the `uv` skill", "Read [the bundled `uv` skill](references/uv/SKILL.md)"
            )
        if revised.encode() != blob:
            files[name] = revised.encode(), mode
            modifications.append(name)
    license_name = spec.get("license", source["content_license"])
    attribution = spec.get("attribution", source["attribution"])
    files["DNAKIT.md"] = (
        (
            f"# DNAKit execution: {skill_id}\n\n{spec['execution_note']}\n\n"
            f"Install this complete folder with `dnakit skills install {skill_id} "
            "--destination ./analysis-skills`. "
            f"Inspect optional dependencies with `dnakit skills doctor {skill_id}`. "
            "Resolve script/resource paths against the installed folder.\n\n"
            "Import validation covers metadata, file integrity and Python syntax. "
            "It does not establish scientific accuracy, service availability "
            "or a successful pipeline run. "
            "Cite the upstream authors and the actual data/methods used; see UPSTREAM.md.\n"
        ).encode(),
        0o644,
    )
    files["UPSTREAM.md"] = (
        (
            f"# Upstream provenance\n\n- Repository: https://github.com/{source['repository']}\n"
            f"- Commit: `{source['revision']}`\n- Original path: `{folder or '.'}`\n"
            "- Original format: "
            f"{spec.get('original_format', 'SKILL.md with complete leaf resources')}\n"
            f"- Attribution: {attribution}\n- License: {license_name}; "
            "see LICENSE and retained source notices.\n"
            f"- Archive SHA-256: `{source.get('archive_sha256')}` "
            "(raw-file sources instead pin every Git blob and SHA-256).\n\n"
            "DNAKit adds this file and DNAKIT.md; normalizes frontmatter "
            "and portable documentation links; "
            "records every original/distribution hash and explicit portability edit in "
            "registry/user-catalog-skills.json and registry/user-catalog-skill-selection.json. "
            "Scientific algorithms are unchanged. Only explicitly selected original Markdown files "
            "are included for converted guides. Software, text, bundled data and "
            "remote services may have distinct terms.\n"
        ).encode(),
        0o644,
    )
    for name, (blob, _) in files.items():
        if name.endswith(".py"):
            ast.parse(blob.decode(), filename=skill_id + "/" + name)
    record = {
        "id": skill_id,
        "role": "dna_workflow",
        "description": description,
        "summary_zh": spec["summary_zh"],
        "category": spec["category"],
        "capability": skill_id,
        "source_id": source["id"],
        "license": license_name,
        "license_basis": f"Pinned {spec['license_path']} retained; "
        "any data/text scope is recorded separately.",
        "upstream_metadata": meta,
        "upstream_license_field": meta.get("license"),
        "upstream_name": meta["name"],
        "upstream_path": folder,
        "upstream_url": f"https://github.com/{source['repository']}/tree/{source['revision']}"
        + ("/" + folder if folder else ""),
        "upstream_revision": source["revision"],
        "attribution": attribution,
        "citation_guidance": "Cite original skill authors, pinned source "
        "and methods/data actually used; see UPSTREAM.md.",
        "aliases": [],
        "alternatives": [],
        "keywords": spec["keywords"],
        "runtime": {
            "python_modules": spec["python_modules"],
            "commands": spec["commands"],
            "environment": spec["environment"],
            "manual_checks": [spec["execution_note"]],
        },
        "required_skills": [],
        "validation": {
            "metadata": "passed",
            "resources": "checksummed",
            "python_syntax": "passed",
            "execution": "not_run",
        },
        "original_format": spec.get("original_format", "SKILL.md"),
        "modifications": sorted(set(modifications)),
        "upstream_files": originals,
        "files": [
            {"path": name, "sha256": digest(blob), "bytes": len(blob), "mode": mode}
            for name, (blob, mode) in sorted(files.items())
        ],
    }
    if "data_license" in spec:
        record["data_license"] = spec["data_license"]
    return record, files


def check(root: Path) -> None:
    selection = read_selection()
    manifest = json.loads((root / MANIFEST).read_text())
    if manifest.get("selection_sha256") != digest(SELECTION.read_bytes()):
        raise ValueError("Selection metadata drift; review and reproduce the import")
    if manifest["sources"] != selection["sources"]:
        raise ValueError("Pinned source metadata mismatch")
    specs = {s["id"]: s for s in selection["skills"]}
    records = manifest["skills"]
    if len(records) != len(specs) or {s["id"] for s in records} != set(specs):
        raise ValueError("Selected skill set mismatch")
    sources = {s["id"]: s for s in selection["sources"]}
    for n, record in enumerate(records, 1):
        spec = specs[record["id"]]
        source = sources[spec["source_id"]]
        for key, expected in {
            "source_id": source["id"],
            "upstream_revision": source["revision"],
            "upstream_path": spec["upstream_path"],
            "license": spec.get("license", source["content_license"]),
            "category": spec["category"],
            "summary_zh": spec["summary_zh"],
            "attribution": spec.get("attribution", source["attribution"]),
            "upstream_url": f"https://github.com/{source['repository']}/tree/{source['revision']}"
            + ("/" + spec["upstream_path"] if spec["upstream_path"] else ""),
            "original_format": spec.get("original_format", "SKILL.md"),
            "data_license": spec.get("data_license"),
            "runtime": {
                "python_modules": spec["python_modules"],
                "commands": spec["commands"],
                "environment": spec["environment"],
                "manual_checks": [spec["execution_note"]],
            },
        }.items():
            if record.get(key) != expected:
                raise ValueError(f"Selected skill identity mismatch: {record['id']}/{key}")
        folder = root / "skills" / record["id"]
        if folder.is_symlink() or not folder.is_dir():
            raise ValueError(f"Unsafe/missing skill folder: {folder}")
        actual = set()
        for path in folder.rglob("*"):
            if path.is_symlink():
                raise ValueError(f"Unsafe linked resource: {path}")
            if path.is_file():
                actual.add(str(path.relative_to(folder)))
        expected = {safe_path(item["path"]) for item in record["files"]}
        if len(expected) != len(record["files"]) or actual != expected:
            raise ValueError(f"Unexpected/missing resource set: {record['id']}")
        for item in record["files"]:
            path = folder / item["path"]
            blob = path.read_bytes()
            if digest(blob) != item["sha256"] or len(blob) != item["bytes"]:
                raise ValueError(f"Resource checksum mismatch: {path}")
            if path.stat().st_mode & 0o777 != item["mode"]:
                raise ValueError(f"Resource mode mismatch: {path}")
            if path.suffix == ".py":
                ast.parse(blob.decode(), filename=str(path))
        meta, _ = metadata((folder / "SKILL.md").read_bytes())
        if meta["name"] != record["id"] or meta["description"] != record["description"]:
            raise ValueError(f"Metadata mismatch: {record['id']}")
        progress(n, len(records), "checked " + record["id"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--cache", type=Path, default=ROOT / "temp/user-catalog-skills-pinned")
    parser.add_argument("--check", action="store_true", help="Read-only offline validation")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.check:
        check(root)
        return 0
    selection = read_selection()
    target_manifest = root / MANIFEST
    if target_manifest.exists() or target_manifest.is_symlink():
        raise FileExistsError(f"Refusing to overwrite {target_manifest}; use --check")
    for spec in selection["skills"]:
        target = root / "skills" / spec["id"]
        if target.exists() or target.is_symlink():
            raise FileExistsError(f"Refusing to overwrite existing skill: {target}")
    if (root / "skills").is_symlink() or (root / "registry").is_symlink():
        raise ValueError("Refusing symlink output parent")
    plans = []
    for n, source in enumerate(selection["sources"], 1):
        upstream = source_files(source, args.cache)
        for spec in selection["skills"]:
            if spec["source_id"] == source["id"]:
                plans.append(make_skill(spec, source, upstream))
        progress(n, len(selection["sources"]), "verified " + source["repository"])
    # All sources, syntax, destinations and edits are checked before any skill is written.
    for n, (record, files) in enumerate(plans, 1):
        folder = root / "skills" / record["id"]
        folder.mkdir(parents=True, exist_ok=False)
        for name, (blob, mode) in sorted(files.items()):
            path = folder / safe_path(name)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as handle:
                handle.write(blob)
            path.chmod(mode)
        progress(n, len(plans), "imported " + record["id"])
    manifest = {
        "schema_version": 1,
        "snapshot_date": selection["snapshot_date"],
        "selection_sha256": digest(SELECTION.read_bytes()),
        "sources": selection["sources"],
        "skills": sorted((r for r, _ in plans), key=lambda r: r["id"]),
    }
    target_manifest.parent.mkdir(parents=True, exist_ok=True)
    with target_manifest.open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    check(root)
    print(f"Imported {len(plans)} reviewed skill increments")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
