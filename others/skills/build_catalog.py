"""Reproduce the reviewed snapshot; never execute upstream workflows.

Use --fetch to obtain missing pinned Git revisions, --plan to review selection in
temp, or --check to validate the existing bundle offline. Source revisions live in
registry/sources.json. Review scope, equivalence and licenses before updating them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import posixpath
import re
import shutil
import subprocess
import tarfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

import yaml
from adaptations import apply_adaptations
from curation_rules import (
    BIO_SUPPORT,
    CATEGORIES,
    EQUIVALENT,
    PREFERRED,
    PRIORITY,
    QUERY_ALIASES,
    classify,
)
from sync_genomics_expansion import (
    extend_metadata,
    is_managed_source,
    managed_ids,
    merge_decisions,
    merge_records,
    preflight,
)
from sync_genomics_expansion import (
    notices as write_expansion_notices,
)

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = ROOT / "registry"
WORK = ROOT / "temp" / "dnakit-skills-work"
DATE = "2026-09-12"
PLAN_DATE = "2026-09-25"
OWNED_SKILLS = {"dnakit", "dnakit-plan"}
ENV = {**os.environ, "GIT_LFS_SKIP_SMUDGE": "1", "GIT_TERMINAL_PROMPT": "0"}
EQUIVALENCE = {name: group for group, names in EQUIVALENT.items() for name in names}


def expansion_review():
    path = REGISTRY / "expansion.json"
    if not path.exists():
        return {}
    rows = json.loads(path.read_text())["candidates"]
    result = {(row["source"], row["path"]): row for row in rows}
    if len(result) != len(rows):
        raise ValueError("Duplicate source/path in the expansion review")
    return result


def reviewed_metadata(data, review):
    """Apply only explicitly reviewed portability edits, keeping original attribution."""
    original = data
    reasons = []
    if review.get("strip_leading_whitespace"):
        data = data.lstrip()
        if data != original:
            reasons.append("removed whitespace before YAML frontmatter")
    original_name = None
    for replacement in review.get("text_replacements", []):
        old, new = replacement["old"].encode(), replacement["new"].encode()
        if data.count(old) != replacement.get("count", 1):
            raise ValueError("Reviewed resource-path edit no longer matches the pinned source")
        data = data.replace(old, new)
        reasons.append(replacement["reason"])
    if review.get("name_override"):
        text = data.decode("utf-8")
        header = re.match(r"\A---\r?\n(.*?)\r?\n---", text, re.S)
        if not header:
            raise ValueError("Cannot rename a skill without YAML frontmatter")
        original_name = yaml.safe_load(header[1]).get("name")
        renamed = re.sub(
            r"^name:.*$",
            "name: " + review["name_override"],
            header[1],
            count=1,
            flags=re.M,
        )
        data = (text[: header.start(1)] + renamed + text[header.end(1) :]).encode()
        reasons.append("renamed the skill identifier for portable discovery")
    metadata, normalized, repaired = parse_skill(data)
    if repaired:
        reasons.append("quoted the YAML description; workflow body unchanged")
    if data != original:
        closing = normalized.find(b"\n---", 4) + len(b"\n---")
        normalized = (
            normalized[:closing]
            + b"\n\n<!-- DNAKit metadata portability edit, 2026-09-12; see UPSTREAM.md. -->"
            + normalized[closing:]
        )
    return metadata, normalized, normalized != original, original_name, reasons


def git(source, *args, **kwargs):
    return subprocess.run(
        ["git", "-C", str(ROOT / source["cache_path"]), *args],
        check=True,
        stdout=subprocess.PIPE,
        env=ENV,
        **kwargs,
    ).stdout


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8"
    )


def parse_skill(data):
    """Repair invalid single-line YAML descriptions only; preserve the workflow body."""
    text = data.decode("utf-8").replace("\r\n", "\n")
    match = re.match(r"\A---\n(.*?)\n---(?:\n|$)", text, re.S)
    if not match:
        raise ValueError("Missing YAML frontmatter")
    header, modified = match[1], False
    try:
        metadata = yaml.safe_load(header)
    except yaml.YAMLError:
        header = re.sub(
            r"^description: (.+)$",
            lambda m: "description: " + json.dumps(m[1], ensure_ascii=False),
            header,
            flags=re.M,
        )
        metadata = yaml.safe_load(header)
        modified = True
    if not isinstance(metadata, dict) or not isinstance(metadata.get("description"), str):
        raise ValueError("Invalid description")
    if not re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", str(metadata.get("name", ""))):
        raise ValueError("Invalid skill name")
    if modified:
        text = (
            "---\n"
            + header
            + "\n---\n"
            + (
                "\n<!-- DNAKit modification, 2026-09-12: quote the YAML "
                "description for parser compatibility. "
                "Workflow body unchanged; see UPSTREAM.md. -->\n"
            )
            + text[match.end() :]
        )
    return metadata, text.encode("utf-8") if modified else data, modified


def license_decision(source, name, metadata):
    field, sid = str(metadata.get("license", "")), source["id"]
    if sid == "omicsclaw":
        return None, "Root Apache-2.0 and per-skill MIT declarations need their scope clarified."
    if sid == "tooluniverse" and name == "tooluniverse-gwas-study-explorer":
        return (
            None,
            "Workflow body states MIT while the repository declares Apache-2.0; review pending.",
        )
    if sid == "kdense":
        if field in {"MIT", "MIT license"}:
            return "MIT", "Individual skill MIT declaration; repository MIT notice retained."
        return (
            None,
            "Individual license needs a complete per-skill redistribution "
            "review; omitted from this core snapshot.",
        )
    if sid == "clawbio":
        if field.upper() == "PROPRIETARY":
            return None, "Explicit proprietary skill; no redistribution permission established."
        if field not in {"", "MIT", "Apache-2.0"}:
            return (
                None,
                "Individual license requires a separate redistribution review; "
                "GPL is not a blanket prohibition.",
            )
        return (
            field or "MIT",
            "ClawBio permits individual skill license overrides; model/data terms remain separate.",
        )
    if sid == "sciagent":
        return (
            "CC-BY-4.0",
            "README licenses original content CC-BY-4.0; frontmatter license "
            "describes underlying tools/data.",
        )
    if sid == "bionemo":
        return (
            "Apache-2.0 AND CC-BY-4.0",
            "Both licenses and NOTICE retained: Apache-2.0 code, CC-BY-4.0 "
            "content; individual AND declaration preserved.",
        )
    return source[
        "content_license"
    ], "Repository content license applies; underlying tools/data are not bundled."


def scan(sources, fetch):
    candidates, decisions, direct_duplicates = [], [], []
    reviews = expansion_review()
    for source in sources:
        cache = ROOT / source["cache_path"]
        if not cache.exists():
            if not fetch:
                raise RuntimeError(f"Missing cache {cache}; run with --fetch")
            cache.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                [
                    "git",
                    "clone",
                    "--filter=blob:none",
                    "--no-checkout",
                    "--depth=1",
                    "https://github.com/" + source["repository"] + ".git",
                    str(cache),
                ],
                check=True,
                env=ENV,
            )
        try:
            git(source, "cat-file", "-e", source["revision"] + "^{commit}")
        except subprocess.CalledProcessError:
            if not fetch:
                raise
            git(source, "fetch", "--depth=1", "origin", source["revision"])
        paths = (
            git(source, "ls-tree", "-r", "--name-only", source["revision"]).decode().splitlines()
        )
        definitions = [p for p in paths if p.endswith("/SKILL.md") or p == "SKILL.md"]
        if source["handling"] != "curated-import":
            decisions.extend(
                {
                    "source": source["id"],
                    "path": p,
                    "status": "reference_only",
                    "reason": source["policy"],
                }
                for p in definitions
            )
            continue
        objects = [source["revision"] + ":" + path for path in definitions]
        stream = io.BytesIO(
            git(source, "cat-file", "--batch", input=("\n".join(objects) + "\n").encode())
        )
        for path in definitions:
            size = int(stream.readline().split()[2])
            data = stream.read(size)
            stream.read(1)
            review = reviews.get((source["id"], path), {})
            if review.get("decision") in {"duplicate", "quality_review", "license_review"}:
                try:
                    original_metadata, _, _ = parse_skill(data)
                    original_name = original_metadata["name"]
                except (ValueError, yaml.YAMLError, UnicodeDecodeError):
                    original_name = PurePosixPath(path).parent.name
                decision = {
                    "source": source["id"],
                    "path": path,
                    "name": original_name,
                    "status": review["decision"],
                    "reason": review["reason"],
                }
                if review["decision"] == "duplicate":
                    decision["canonical_id"] = review["canonical_id"]
                    direct_duplicates.append(decision)
                else:
                    decisions.append(decision)
                continue
            try:
                metadata, normalized, modified, original_name, changes = reviewed_metadata(
                    data, review
                )
                name = metadata["name"]
            except (ValueError, yaml.YAMLError, UnicodeDecodeError) as exc:
                decisions.append(
                    {
                        "source": source["id"],
                        "path": path,
                        "status": "invalid_metadata",
                        "reason": str(exc).splitlines()[0],
                    }
                )
                continue
            decision = {"source": source["id"], "path": path, "name": name}
            category = review.get("category") or classify(source["id"], path, name)
            if category is None:
                decisions.append(
                    {
                        **decision,
                        "status": "out_of_scope",
                        "reason": (
                            "Outside reviewed DNA/genomics tasks, or a nested example/legacy copy."
                        ),
                    }
                )
                continue
            license_id, reason = license_decision(source, name, metadata)
            if license_id is None:
                decisions.append({**decision, "status": "license_review", "reason": reason})
                continue
            candidates.append(
                {
                    **decision,
                    "category": category,
                    "capability": EQUIVALENCE.get(name, name),
                    "license": license_id,
                    "license_basis": reason,
                    "metadata": metadata,
                    "original_name": original_name,
                    "changes": changes,
                    "review": review,
                    "normalized": normalized,
                    "modified": modified,
                }
            )
    groups = defaultdict(list)
    for item in candidates:
        groups[item["capability"]].append(item)
    selected = []
    for group, options in sorted(groups.items()):
        preferred = PREFERRED.get(group)
        options.sort(
            key=lambda x: (
                0 if preferred in (x["name"], x["path"]) else 1,
                PRIORITY.get(x["source"], 100),
                x["path"],
            )
        )
        winner, *alternatives = options
        winner["alternatives"] = [
            {key: x[key] for key in ("source", "path", "name")} for x in alternatives
        ]
        selected.append(winner)
        decisions.append(
            {
                "source": winner["source"],
                "path": winner["path"],
                "name": winner["name"],
                "status": "selected",
                "canonical_id": winner["name"],
                "capability": group,
                "reason": winner["license_basis"],
            }
        )
        for alternative in alternatives:
            decisions.append(
                {
                    **{k: alternative[k] for k in ("source", "path", "name")},
                    "status": "duplicate",
                    "canonical_id": winner["name"],
                    "capability": group,
                    "reason": (
                        "Same reviewed task/stage; retain the chosen implementation and "
                        "this provenance record."
                    ),
                }
            )
    winners = {x["name"]: x for x in selected}
    for decision in direct_duplicates:
        target = decision["canonical_id"]
        if target not in winners:
            raise ValueError(f"Reviewed duplicate refers to missing skill: {target}")
        winner = winners[target]
        winner["alternatives"].append({key: decision[key] for key in ("source", "path", "name")})
        decisions.append({**decision, "capability": winner["capability"]})
    if len({x["name"] for x in selected}) != len(selected):
        raise ValueError("Selected names collide; review equivalence rules")
    return sorted(selected, key=lambda x: x["name"]), sorted(
        decisions, key=lambda x: (x["source"], x["path"])
    )


def archive_files(source, paths, excluded=None):
    """Read a pinned Git archive; reject symlinks, traversal, and LFS placeholders."""
    data = git(
        source,
        "-c",
        "core.autocrlf=false",
        "archive",
        "--format=tar",
        source["revision"],
        "--",
        *paths,
    )
    result = {}
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(f"Unsafe archive path {member.name}")
            if member.name in (excluded or {}):
                expected = excluded[member.name]
                if not member.issym() or member.linkname != expected["symlink_target"]:
                    raise ValueError(f"Reviewed symlink no longer matches: {member.name}")
                continue
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError(f"Unsupported upstream resource {member.name}: not a regular file")
            if member.size > 10_000_000:
                raise ValueError(f"Review large resource {member.name} before bundling")
            content = archive.extractfile(member).read()
            if content.startswith(b"version https://git-lfs.github.com/spec/v1"):
                raise ValueError(f"Unresolved LFS resource {member.name}; do not ship a pointer")
            result[member.name] = (content, member.mode)
    return result


def portable_doc_links(data, path, old_root, source):
    """Keep internal documentation links local and pin links outside the copied tree."""
    text = data.decode("utf-8")
    in_fence = False
    lines = []

    def replace(match):
        target = match[2]
        if target.startswith(("#", "/")) or urlsplit(target).scheme:
            return match[0]
        pathname, marker, fragment = target.partition("#")
        resolved = posixpath.normpath(posixpath.join(posixpath.dirname(path), pathname))
        if resolved.startswith(old_root + "/"):
            return match[0]
        pinned = f"https://github.com/{source['repository']}/blob/{source['revision']}/{resolved}"
        return match[1] + pinned + (marker + fragment if marker else "") + match[3]

    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith(("```", "~~~")):
            in_fence = not in_fence
        if not in_fence:
            line = re.sub(r"(\]\()([^\s)]+)(\))", replace, line)
        lines.append(line)
    result = "".join(lines)
    if result != text:
        result += "\n<!-- DNAKit link portability edit, 2026-09-12; see UPSTREAM.md. -->\n"
    return result.encode("utf-8")


def runtime_for(item):
    metadata = item["metadata"]
    runtime = {
        "python_modules": [],
        "commands": [],
        "environment": [],
        "manual_checks": [
            "Read SKILL.md for complete dependencies, versions, reference "
            "files, API access, and model/data terms."
        ],
    }
    primary = str(metadata.get("primary_tool", ""))
    if metadata.get("tool_type") == "cli" and re.fullmatch(r"[A-Za-z0-9_.+-]+", primary):
        runtime["commands"].append(primary)
    if metadata.get("tool_type") == "python" and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.]*", primary):
        runtime["python_modules"].append({"biopython": "Bio"}.get(primary, primary))
    if item["source"] == "tooluniverse":
        runtime["python_modules"].append("tooluniverse")
        runtime["manual_checks"].append(
            "Configure ToolUniverse and the APIs used by this workflow; "
            "DNAKit has a separate tool catalog."
        )
    nested = metadata.get("metadata", {})
    if isinstance(nested, dict):
        requires = (nested.get("openclaw") or {}).get("requires") or {}
        runtime["commands"].extend(requires.get("bins") or [])
        runtime["environment"].extend(requires.get("env") or [])
        if requires.get("anyBins"):
            runtime["manual_checks"].append(
                "One of these commands is required: " + ", ".join(requires["anyBins"])
            )
    if item["name"] == "genomic-coordinates":
        runtime["python_minimum"] = "3.11"
    for key in ("python_modules", "commands", "environment"):
        runtime[key] = sorted(set(runtime[key]))
    return runtime


def file_manifest(directory, modes=None):
    return [
        {
            "path": str(p.relative_to(directory)),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "bytes": p.stat().st_size,
            "mode": (modes or {}).get(str(p.relative_to(directory)), 0o644),
        }
        for p in sorted(directory.rglob("*"))
        if p.is_file()
    ]


def import_selected(sources, selected):
    staging = WORK / "staged-skills"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    archives, licenses = {}, {}
    for source in sources:
        chosen = [x for x in selected if x["source"] == source["id"]]
        if not chosen:
            continue
        paths = (
            git(source, "ls-tree", "-r", "--name-only", source["revision"]).decode().splitlines()
        )
        notices = [
            p
            for p in paths
            if "/" not in p and p.upper().startswith(("LICENSE", "NOTICE", "COPYING", "CITATION"))
        ]
        notices.extend(source.get("extra_notices", []))
        directories = [str(PurePosixPath(x["path"]).parent) for x in chosen]
        print(f"Reading {source['id']}: {len(chosen)} complete skill folders", flush=True)
        archives[source["id"]] = archive_files(
            source,
            [*notices, *directories, *source.get("portable_docs", {})],
            {
                item["path"].removesuffix("SKILL.md") + path: rule
                for item in chosen
                for path, rule in item.get("review", {}).get("excluded_resources", {}).items()
            },
        )
        licenses[source["id"]] = {p: archives[source["id"]][p] for p in notices}
    source_by_id = {x["id"]: x for x in sources}
    path_mapping = {(x["source"], x["path"].removesuffix("/SKILL.md")): x["name"] for x in selected}
    for item in selected:
        for alt in item["alternatives"]:
            path_mapping[(alt["source"], alt["path"].removesuffix("/SKILL.md"))] = item["name"]
    records = []
    for item in selected:
        sid, name = item["source"], item["name"]
        source = source_by_id[sid]
        directory = staging / name
        directory.mkdir()
        prefix = item["path"].removesuffix("SKILL.md")
        modes, upstream_files = {}, []
        modified = []
        changes = item.get("changes", [])[:]
        for path, (data, mode) in archives[sid].items():
            if not path.startswith(prefix):
                continue
            relative = path[len(prefix) :]
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(item["normalized"] if relative == "SKILL.md" else data)
            target.chmod(mode)
            modes[relative] = mode
            upstream_files.append({"path": relative, "sha256": hashlib.sha256(data).hexdigest()})
            if relative == "SKILL.md" and item["modified"]:
                modified.append(relative)
        for old_root, new_root in source.get("portable_docs", {}).items():
            for path, (data, mode) in archives[sid].items():
                if not path.startswith(old_root + "/"):
                    continue
                relative = new_root + "/" + path.removeprefix(old_root + "/")
                target = directory / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                portable = (
                    portable_doc_links(data, path, old_root, source)
                    if path.endswith(".md")
                    else data
                )
                target.write_bytes(portable)
                if portable != data:
                    modified.append(relative)
                    changes.append("pinned external links in copied repository documentation")
                modes[relative] = mode
                upstream_files.append(
                    {
                        "path": relative,
                        "upstream_path": path,
                        "sha256": hashlib.sha256(data).hexdigest(),
                    }
                )
            # Only relocate repository-wide documentation links; analysis steps stay intact.
            skill_path = directory / "SKILL.md"
            before = skill_path.read_text(encoding="utf-8")
            after = before.replace("../../" + old_root + "/", new_root + "/")
            if after != before:
                after += (
                    "\n<!-- DNAKit documentation-link portability edit, 2026-09-12: "
                    "bundled repository docs locally; see UPSTREAM.md. -->\n"
                )
                skill_path.write_text(after, encoding="utf-8")
                modified.append("SKILL.md")
                changes.append("relocated repository documentation links and bundled their targets")
        notice_owners = [sid, "kdense", "tooluniverse"] if sid == "sciagent" else [sid]
        for owner in notice_owners:
            for path, (data, _) in licenses[owner].items():
                target = directory / "licenses" / owner / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
        if sid == "sciagent" or sid == "deepmind":
            (directory / "licenses" / sid / "CC-BY-4.0.txt").write_bytes(
                archives["bionemo"]["LICENSE-CC-BY-4.0"][0]
            )
        metadata = item["metadata"]
        nested = metadata.get("metadata", {})
        author = (
            (nested.get("skill-author") or nested.get("author"))
            if isinstance(nested, dict)
            else None
        )
        attribution = (str(author) + "; " if author else "") + source["attribution"]
        url = f"https://github.com/{source['repository']}/tree/{source['revision']}/{prefix.rstrip('/')}"
        modified = sorted(set(modified))
        required = [
            {
                "upstream_path": path,
                "id": path_mapping.get((sid, path)),
                "status": "bundled" if (sid, path) in path_mapping else "upstream_only",
            }
            for path in metadata.get("depends_on", [])
        ]
        runtime = runtime_for(item)
        for declarations in (source.get("runtime", {}), item.get("review", {}).get("runtime", {})):
            for key, values in declarations.items():
                runtime.setdefault(key, []).extend(values)
        unresolved = [x["upstream_path"] for x in required if x["status"] == "upstream_only"]
        if unresolved:
            runtime["manual_checks"].append(
                "Upstream cross-skill dependencies not bundled: " + ", ".join(unresolved)
            )
        record = {
            "id": name,
            "role": (
                item.get("review", {}).get("role")
                or (
                    "workflow_dependency"
                    if sid == "bioskills" and prefix.rstrip("/") in BIO_SUPPORT
                    else "dna_workflow"
                )
            ),
            "description": metadata["description"].strip(),
            "summary_zh": CATEGORIES[item["category"]]["zh"]
            + "\uff1a"
            + name.removeprefix("bio-").removeprefix("tooluniverse-").replace("-", " "),
            "category": item["category"],
            "capability": item["capability"],
            "source_id": sid,
            "license": item["license"],
            "license_basis": item["license_basis"],
            "upstream_license_field": metadata.get("license"),
            "upstream_metadata": metadata,
            "upstream_name": item.get("original_name") or name,
            "upstream_url": url,
            "upstream_path": prefix.rstrip("/"),
            "upstream_revision": source["revision"],
            "attribution": attribution,
            "citation_guidance": source["citation_guidance"],
            "aliases": sorted(
                {
                    x["name"]
                    for x in item["alternatives"]
                    if re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", x["name"])
                }
                - {name}
            ),
            "alternatives": item["alternatives"],
            "keywords": [item["category"], CATEGORIES[item["category"]]["en"]],
            "runtime": runtime,
            "required_skills": required,
            "validation": {
                "metadata": "passed",
                "resources": "checksummed",
                "execution": "not_run",
            },
            "modifications": modified,
            "excluded_upstream_resources": item.get("review", {}).get("excluded_resources", {}),
            "upstream_files": sorted(upstream_files, key=lambda x: x["path"]),
        }
        for rule in item.get("review", {}).get("excluded_resources", {}).values():
            changes.append(rule["reason"])
        change = "; ".join(dict.fromkeys(changes)) if changes else "none"
        # Keep the original first-round notice stable when this was its sole edit.
        if changes == ["quoted the YAML description; workflow body unchanged"]:
            change = "quoted the YAML description in SKILL.md; workflow body unchanged"
        provenance = (
            f"# Upstream provenance: {name}\n\n- Source: [{source['repository']}]({url})\n"
            f"- Revision: `{source['revision']}`\n- Original path: `{prefix.rstrip('/')}`\n"
            f"- Attribution: {attribution}\n- Content license: **{item['license']}**. "
            f"{item['license_basis']}\n"
            f"- Snapshot: {DATE}\n- DNAKit changes to upstream files: {change}.\n"
            "- DNAKit additions: this file and copies of source notices in `licenses/`.\n\n"
            "## Runtime and validation\n\nThe original workflow has not been "
            "executed during import. "
            "Metadata parsing and checksums do not establish scientific "
            "correctness. Install the upstream tools, "
            "reference data and API credentials specified in SKILL.md. Their "
            "licenses remain separate. "
            "This bundle does not include those tools, model weights or databases. "
            "Resolve original resource paths from this skill directory. "
            "Repository-wide dependencies may require the upstream runtime. "
            f"Use `dnakit skills info {name} --json` for cross-skill dependencies.\n\n"
            f"## Citation\n\n{source['citation_guidance']}\n\n"
            f"{attribution}. {name}. Revision {source['revision']}. {url}\n\n"
            "Also cite the actual tools, models, datasets and papers used. "
            "Inclusion does not imply upstream endorsement.\n"
        )
        (directory / "UPSTREAM.md").write_text(provenance, encoding="utf-8")
        record["files"] = file_manifest(directory, modes)
        records.append(record)
    old_catalog = REGISTRY / "catalog.json"
    if old_catalog.exists():
        old_ids = {x["id"] for x in json.loads(old_catalog.read_text())["skills"]} - OWNED_SKILLS
    else:
        baseline = WORK / "baseline-skills"
        old_ids = (
            {p.name for p in baseline.iterdir() if p.is_dir() and (p / "UPSTREAM.md").exists()}
            if baseline.exists()
            else set()
        ) - OWNED_SKILLS
    # Explicit extension manifests own these folders, including native workflows.
    # Their checked resources are merged after baseline adaptations are replayed.
    old_ids -= managed_ids()
    for item in selected:
        target = ROOT / "skills" / item["name"]
        if target.exists() and item["name"] not in old_ids:
            raise RuntimeError(f"Unmanaged target skill: {target}")
    for name in old_ids:
        target = ROOT / "skills" / name
        if target.exists():
            shutil.rmtree(target)
    for directory in staging.iterdir():
        shutil.move(str(directory), ROOT / "skills" / directory.name)
    for sid, files in licenses.items():
        for path, (data, _) in files.items():
            destination = ROOT / "licenses" / "skills" / sid / path
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
    return records


def entry_record():
    directory = ROOT / "skills" / "dnakit"
    (directory / "LICENSE").write_bytes((ROOT / "LICENSE").read_bytes())
    metadata, _, _ = parse_skill((directory / "SKILL.md").read_bytes())
    return {
        "id": "dnakit",
        "description": metadata["description"],
        "summary_zh": "DNAKit 统一技能检索与工作流入口",
        "category": "workflows",
        "capability": "dnakit-discovery",
        "source_id": "dnakit",
        "license": "MIT",
        "upstream_name": "DNAKit skill discovery",
        "upstream_url": "https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit",
        "upstream_revision": "catalog-" + DATE,
        "attribution": "DNAKit contributors",
        "citation_guidance": (
            "Cite the DNAKit software version and selected upstream "
            "workflow. No paper DOI has been assigned by this catalog."
        ),
        "aliases": ["dna-skills"],
        "keywords": ["DNA", "genomics", "search", "workflow", "检索", "入口"],
        "runtime": {
            "python_modules": ["dnakit"],
            "commands": [],
            "environment": [],
            "manual_checks": ["Selected workflows have additional dependencies."],
        },
        "required_skills": [],
        "validation": {"metadata": "passed", "resources": "checksummed", "execution": "not_run"},
        "files": file_manifest(directory),
    }


def plan_entry_record():
    """Record DNAKit's Codex planning skill alongside its workflow catalog."""
    directory = ROOT / "skills" / "dnakit-plan"
    (directory / "LICENSE").write_bytes((ROOT / "LICENSE").read_bytes())
    metadata, _, _ = parse_skill((directory / "SKILL.md").read_bytes())
    return {
        "id": "dnakit-plan",
        "description": metadata["description"],
        "summary_zh": "DNAKit MCP 任务规划与确认后保存计划",
        "category": "workflows",
        "capability": "dnakit-planning",
        "source_id": "dnakit",
        "license": "MIT",
        "upstream_name": "DNAKit planning skill",
        "upstream_url": "https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-plan",
        "upstream_revision": "dnakit-plan-" + PLAN_DATE,
        "attribution": "DNAKit contributors",
        "citation_guidance": (
            "Cite the DNAKit software version and the data, tools, and papers used "
            "when the approved plan is executed. No paper DOI has been assigned to this skill."
        ),
        "aliases": [],
        "keywords": ["DNAKit", "MCP", "Codex", "plan", "planning", "规划", "计划"],
        "runtime": {
            "python_modules": ["dnakit"],
            "commands": ["codex"],
            "environment": [],
            "manual_checks": ["Connect DNAKit MCP for tool discovery when available."],
        },
        "required_skills": [],
        "validation": {"metadata": "passed", "resources": "checksummed", "execution": "not_run"},
        "files": file_manifest(directory),
    }


def disambiguate_aliases(records):
    """A shared name is searchable metadata, never an arbitrary install redirect."""
    identifiers = {entry["id"] for entry in records}
    owners = defaultdict(set)
    for entry in records:
        # Retained alternatives recover shared names after a prior pass removed
        # them from install aliases. This makes incremental adaptation replay
        # preserve the same ambiguity metadata as a fresh upstream import.
        names = set(entry["aliases"]) | {
            alternative["name"]
            for alternative in entry.get("alternatives", [])
            if re.fullmatch(r"[a-z0-9][a-z0-9-]{0,63}", alternative["name"])
            and alternative["name"] != entry["id"]
        }
        for alias in names:
            owners[alias].add(entry["id"])
    ambiguous = {
        name: sorted(values | ({name} if name in identifiers else set()))
        for name, values in owners.items()
        if len(values) > 1 or name in identifiers
    }
    for entry in records:
        shared = set(entry["aliases"]) & ambiguous.keys()
        entry["keywords"] = sorted(set(entry["keywords"]) | shared)
        entry["aliases"] = sorted(set(entry["aliases"]) - shared)
    return ambiguous


def check():
    catalog = json.loads((REGISTRY / "catalog.json").read_text())
    identifiers, capabilities, count = set(), set(), 0
    for entry in catalog["skills"]:
        assert entry["id"] not in identifiers, entry["id"]
        assert entry["capability"] not in capabilities, entry["capability"]
        identifiers.add(entry["id"])
        capabilities.add(entry["capability"])
        directory = ROOT / "skills" / entry["id"]
        metadata, _, _ = parse_skill((directory / "SKILL.md").read_bytes())
        assert metadata["name"] == entry["id"]
        assert {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()} == {
            x["path"] for x in entry["files"]
        }
        for file in entry["files"]:
            path = directory / file["path"]
            assert not path.is_symlink() and path.resolve().is_relative_to(directory.resolve())
            assert hashlib.sha256(path.read_bytes()).hexdigest() == file["sha256"], path
            count += 1
        for file in entry.get("upstream_files", []):
            if file["path"] not in entry.get("modifications", []):
                assert (
                    hashlib.sha256((directory / file["path"]).read_bytes()).hexdigest()
                    == file["sha256"]
                )
        assert any("LICENSE" in x["path"].upper() for x in entry["files"])
    print(
        f"Validated {len(identifiers)} definitions and {count} files; "
        "execution is a separate check."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plan", action="store_true")
    args = parser.parse_args()
    if args.check:
        check()
        return
    sources = json.loads((REGISTRY / "sources.json").read_text())["sources"]
    baseline_sources = [source for source in sources if not is_managed_source(source)]
    selected, decisions = scan(baseline_sources, args.fetch)
    dump(
        WORK / "selection-plan.json",
        [{k: v for k, v in x.items() if k not in {"normalized", "metadata"}} for x in selected],
    )
    dump(WORK / "decisions-plan.json", decisions)
    print(
        json.dumps(
            {
                "selected": len(selected),
                "sources": dict(Counter(x["source"] for x in selected)),
                "decisions": dict(Counter(x["status"] for x in decisions)),
            },
            indent=2,
        ),
        flush=True,
    )
    if args.plan:
        return
    old_catalog = REGISTRY / "catalog.json"
    preflight(json.loads(old_catalog.read_text())["skills"] if old_catalog.exists() else [])
    records = [*import_selected(baseline_sources, selected), entry_record(), plan_entry_record()]
    apply_adaptations(records, decisions)
    records = merge_records(records)
    decisions = merge_decisions(decisions)
    ambiguous_names = disambiguate_aliases(records)
    evidence_path = REGISTRY / "validation.json"
    if evidence_path.exists():
        evidence = json.loads(evidence_path.read_text())
        for entry in records:
            proof = next((x for x in evidence["records"] if x["id"] == entry["id"]), None)
            digest = hashlib.sha256(json.dumps(entry["files"], sort_keys=True).encode()).hexdigest()
            if proof and proof["resource_manifest_sha256"] == digest:
                entry["validation"].update(
                    execution=proof["execution"],
                    evidence="registry/validation.json",
                    scope=evidence["scope"],
                )
    catalog = {
        "schema_version": 1,
        "snapshot_date": DATE,
        "categories": CATEGORIES,
        "query_aliases": QUERY_ALIASES,
        "ambiguous_names": ambiguous_names,
        "sources": [{"id": "dnakit", "repository": "mapengsen/DNAKit"}, *sources],
        "selection": dict(Counter(x["status"] for x in decisions)),
        "skills": sorted(records, key=lambda x: x["id"]),
    }
    extend_metadata(catalog)
    catalog["snapshot_date"] = max(catalog["snapshot_date"], PLAN_DATE)
    dump(REGISTRY / "catalog.json", catalog)
    dump(
        REGISTRY / "decisions.json",
        {"schema_version": 1, "snapshot_date": catalog["snapshot_date"], "decisions": decisions},
    )
    with (REGISTRY / "skills.csv").open("w", encoding="utf-8", newline="") as handle:
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
    dump(
        ROOT / "skills" / "upstream-manifest.json",
        {
            "schema_version": 2,
            "snapshot_date": catalog["snapshot_date"],
            "catalog": "../registry/catalog.json",
            "decisions": "../registry/decisions.json",
            "skills": [
                {"id": x["id"], "source": x["source_id"], "revision": x["upstream_revision"]}
                for x in catalog["skills"]
            ],
        },
    )
    from write_migration import main as write_migration

    write_migration()
    write_expansion_notices()
    check()


if __name__ == "__main__":
    main()
