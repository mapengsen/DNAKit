#!/usr/bin/env python3
"""Audit the user's 335-project workbook without modifying the supplied files.

This version cross-checks Markdown and XLSX, preserves stable project identities,
and records primary-page retrieval separately from software execution and from
verified citation relationships. Fetching reads public documentation only.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urlsplit
from xml.etree import ElementTree as ET

from rich.progress import track

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "registry/user-catalog-335.json"
EVIDENCE = ROOT / "registry/user-catalog-evidence.json"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
FILES = ("CATALOG.md", "FAQ.md", "DNA_catalog.html", "DNA_catalog.xlsx")
SHEETS = ("目录", "文献与存档", "引用与依赖", "功能覆盖", "Skills来源")
COUNTS = dict(zip(SHEETS, (335, 62, 86, 38, 52), strict=True))
URL_ALIASES = {
    "github.com/k-dense-ai/claude-scientific-skills": (
        "github.com/k-dense-ai/scientific-agent-skills"
    ),
    "github.com/weberlab-hhu/helixer": "github.com/usadellab/helixer",
    "github.com/fenderglass/flye": "github.com/mikolmogorov/flye",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def url_key(value: str) -> str:
    """Normalize exact URLs, retaining subdirectories and version distinctions."""
    parts = urlsplit(value)
    key = (parts.netloc + parts.path).casefold().rstrip("/").removesuffix(".git")
    return URL_ALIASES.get(key, key)


def name_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.casefold())


def workbook(path: Path) -> dict[str, list[dict[str, str]]]:
    """Read cached cell text only; never execute spreadsheet formulas/macros."""
    result = {}
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = [
                "".join(item.itertext())
                for item in ET.fromstring(archive.read("xl/sharedStrings.xml"))
            ]
        rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
        rels = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        for sheet in book.findall("m:sheets/m:sheet", NS):
            title = sheet.attrib["name"]
            if title not in SHEETS:
                continue
            target = rels[sheet.attrib[f"{{{rel_ns}}}id"]]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            root = ET.fromstring(archive.read(target))
            rows = []
            for row in root.findall("m:sheetData/m:row", NS):
                values: dict[int, str] = {}
                for cell in row.findall("m:c", NS):
                    number = 0
                    for letter in re.sub(r"\d", "", cell.attrib["r"]):
                        number = number * 26 + ord(letter) - ord("A") + 1
                    value = cell.findtext("m:v", "", NS)
                    if cell.attrib.get("t") == "s":
                        value = strings[int(value)]
                    elif cell.attrib.get("t") == "inlineStr":
                        node = cell.find("m:is", NS)
                        value = "".join(node.itertext()) if node is not None else ""
                    values[number - 1] = value
                rows.append([values.get(i, "") for i in range(max(values, default=-1) + 1)])
            header = rows[2]
            if len(set(header)) != len(header) or "" in header:
                raise ValueError(f"Invalid or duplicate headers in {title}")
            records = [
                {key: row[i] if i < len(row) else "" for i, key in enumerate(header)}
                for row in rows[3:]
                if any(row)
            ]
            if len(records) != COUNTS[title]:
                raise ValueError(f"Unexpected record count in {title}: {len(records)}")
            result[title] = records
    if set(result) != set(SHEETS):
        raise ValueError("The workbook is missing an expected sheet")
    return result


def markdown_rows(path: Path) -> list[tuple[str, str]]:
    result = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|") or "[项目]" not in line:
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if len(cells) != 6:
            raise ValueError("Malformed Markdown project row")
        match = re.search(r"\[项目\]\(([^)]+)\)", cells[3])
        if match is None:
            raise ValueError("Missing Markdown project URL")
        result.append((cells[0], match.group(1)))
    return result


def import_catalog(directory: Path) -> dict[str, Any]:
    sheets = workbook(directory / "DNA_catalog.xlsx")
    md = markdown_rows(directory / "CATALOG.md")
    pairs = [(row["项目 / 接口"], row["项目入口"]) for row in sheets["目录"]]
    if md != pairs or len(pairs) != 335 or len(set(pairs)) != 335:
        raise ValueError("Markdown and workbook project identities/order disagree")
    previous = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else None
    old = json.loads((ROOT / "registry/genomics-survey.json").read_text())["projects"]
    by_url: dict[str, set[str]] = {}
    by_name: dict[str, set[str]] = {}
    for project in old:
        by_name.setdefault(name_key(project["name"]), set()).add(project["id"])
        for link in [project.get("homepage"), project.get("repository")]:
            if link:
                by_url.setdefault(url_key(link), set()).add(project["id"])
    pinned = {r["input_stable_id"]: r for r in previous["entries"]} if previous else {}
    entries = []
    seen_ids: set[str] = set()
    for index, row in enumerate(track(sheets["目录"], description="Cross-check 335 projects"), 1):
        stable_id = row["稳定ID"]
        if stable_id in pinned:
            original = pinned[stable_id]
            if (original["name"], original["project_url"]) != pairs[index - 1]:
                raise ValueError(f"Identity changed for previously reviewed row: {stable_id}")
            project_id = original["project_id"]
            existing = original["existing_project_id"]
            match_method = original["match_method"]
        else:
            candidates = by_url.get(url_key(row["项目入口"]), set())
            match_method = "exact_repository_or_homepage"
            if not candidates:
                candidates = by_name.get(name_key(row["项目 / 接口"]), set())
                match_method = "exact_normalized_name"
            if len(candidates) > 1:
                raise ValueError(f"Ambiguous identity: {row['项目 / 接口']} -> {candidates}")
            existing = next(iter(candidates)) if candidates else None
            project_id = existing or re.sub(
                r"[^a-z0-9]+", "-", row["项目 / 接口"].casefold()
            ).strip("-")
            if not existing:
                match_method = "new_project_record"
        if project_id in seen_ids:
            raise ValueError(f"Multiple workbook rows map to {project_id}")
        seen_ids.add(project_id)
        entries.append(
            {
                "index": index,
                "input_stable_id": stable_id,
                "name": row["项目 / 接口"],
                "project_id": project_id,
                "existing_project_id": existing,
                "match_method": match_method,
                "category": row["方向"],
                "purpose": row["功能"],
                "interface": row["接入类型"],
                "project_url": row["项目入口"],
                "evidence_url": row["项目核验来源"],
                "paper_url": row["论文 / 存档入口"] or None,
                "notes": row["重要限制与备注"],
                "input_fields": row,
                "input_evidence_status": (
                    "user_supplied; independently checked evidence is stored separately"
                ),
            }
        )
    return {
        "schema_version": 1,
        "snapshot_date": "2026-09-20",
        "input_files": [
            {"path": "temp/" + name, "sha256": sha256((directory / name).read_bytes())}
            for name in FILES
        ],
        "baseline_project_count": previous["baseline_project_count"] if previous else len(old),
        "counts": {
            "entries": len(entries),
            "existing": sum(bool(e["existing_project_id"]) for e in entries),
            "new": sum(not e["existing_project_id"] for e in entries),
            **COUNTS,
        },
        "entries": entries,
        "publications": sheets["文献与存档"],
        "relationships": sheets["引用与依赖"],
        "coverage_claims": sheets["功能覆盖"],
        "skill_source_claims": sheets["Skills来源"],
    }


def retrieval_urls(url: str) -> list[str]:
    """Prefer an upstream README; retain subdirectory README paths."""
    parts = urlsplit(url)
    bits = parts.path.strip("/").split("/")
    if parts.netloc.casefold() != "github.com" or len(bits) < 2:
        return [url]
    repo = "/".join(bits[:2])
    if len(bits) > 4 and bits[2] == "blob":
        return [f"https://raw.githubusercontent.com/{repo}/" + "/".join(bits[3:])]
    ref, subdir = "HEAD", ""
    if len(bits) > 4 and bits[2] == "tree":
        ref, subdir = bits[3], "/".join(bits[4:]) + "/"
    return [
        f"https://raw.githubusercontent.com/{repo}/{ref}/{subdir}{name}"
        for name in (
            "README.md",
            "README.rst",
            "README",
            "readme.md",
            "README.txt",
            "README.markdown",
            "readme.MD",
        )
    ]


def fetch_page(url: str, cache: Path) -> dict[str, Any]:
    attempts = []
    for address in retrieval_urls(url):
        try:
            request = urllib.request.Request(
                address, headers={"User-Agent": "DNAKit-source-audit/1.0"}
            )
            with urllib.request.urlopen(request, timeout=20) as response:
                body = response.read(3_000_001)
                if len(body) > 3_000_000:
                    raise ValueError("Source exceeds the 3 MB documentation limit")
                if not body.strip():
                    raise ValueError("Empty documentation response")
                digest = sha256(body)
                (cache / (digest + ".txt")).write_bytes(body)
                return {
                    "url": url,
                    "retrieved_url": response.url,
                    "status": "retrieved",
                    "http_status": response.status,
                    "sha256": digest,
                    "bytes": len(body),
                    "checked_at": datetime.now(timezone.utc).isoformat(),
                    "validation_boundary": (
                        "documentation retrieval; not execution, license clearance "
                        "or scientific validation"
                    ),
                }
        except (OSError, ValueError, urllib.error.URLError) as error:
            attempts.append(
                {
                    "url": address,
                    "error_type": type(error).__name__,
                    "http_status": getattr(error, "code", None),
                }
            )
    return {
        "url": url,
        "status": "unresolved",
        "attempts": attempts,
        "checked_at": datetime.now(timezone.utc).isoformat(),
    }


def fetch_evidence(data: dict[str, Any], cache: Path) -> dict[str, Any]:
    cache.mkdir(parents=True, exist_ok=True)
    urls = {e[key] for e in data["entries"] for key in ("project_url", "evidence_url") if e[key]}
    urls.update(r["证据入口"] for r in data["relationships"])
    urls.add("https://github.com/ChongxiaoLiang/BOLE")
    pages = {}
    with ThreadPoolExecutor(max_workers=10) as pool:
        tasks = {pool.submit(fetch_page, url, cache): url for url in sorted(urls)}
        for task in track(
            as_completed(tasks), total=len(tasks), description="Read primary documentation"
        ):
            pages[tasks[task]] = task.result()
    relationships = []
    for index, relation in enumerate(data["relationships"], 1):
        page = pages[relation["证据入口"]]
        text = ""
        if page["status"] == "retrieved":
            text = html.unescape((cache / (page["sha256"] + ".txt")).read_text(errors="replace"))
        target = relation["对象"]
        needle = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", target).rstrip("/")
        observed = needle.casefold() in unquote(text).casefold()
        # A target appearing on a page does not establish an actual dependency or citation.
        relationships.append(
            {
                "index": index,
                "source": relation["来源"],
                "relation_type": relation["关系类型"],
                "target": target,
                "evidence_url": relation["证据入口"],
                "page_status": page["status"],
                "target_observed": observed,
                "relationship_verified": False,
                "verification": "target_text_observed_relation_requires_context_review"
                if observed
                else "user_supplied_relation_not_confirmed_in_retrieved_page",
            }
        )
    return {
        "schema_version": 1,
        "snapshot_date": data["snapshot_date"],
        "input_sha256": sha256(json.dumps(data, ensure_ascii=False, sort_keys=True).encode()),
        "pages": [pages[url] for url in sorted(pages)],
        "relationships": relationships,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path, default=ROOT / "temp")
    parser.add_argument("--fetch-evidence", action="store_true")
    parser.add_argument("--evidence-cache", type=Path, default=Path("/tmp/dnakit-335-evidence"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check and args.fetch_evidence:
        parser.error("--check is offline and cannot be combined with --fetch-evidence")
    data = import_catalog(args.source_directory)
    content = json.dumps(data, ensure_ascii=False, indent=2) + "\n"
    if args.check:
        if not MANIFEST.exists() or MANIFEST.read_text(encoding="utf-8") != content:
            raise SystemExit("The imported catalog is stale or supplied files have changed")
    else:
        MANIFEST.write_text(content, encoding="utf-8")
    if args.fetch_evidence:
        evidence = fetch_evidence(data, args.evidence_cache)
        EVIDENCE.write_text(
            json.dumps(evidence, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        retrieved = sum(p["status"] == "retrieved" for p in evidence["pages"])
        print(f"Retrieved {retrieved}/{len(evidence['pages'])} primary pages")
    print(
        f"{'Checked' if args.check else 'Imported'} 335 rows: "
        f"{data['counts']['existing']} existing; {data['counts']['new']} new source records"
    )


if __name__ == "__main__":
    main()
