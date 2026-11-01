"""Map each Agent-callable DNAKit function to the place in the documentation that explains it.

Run from the source checkout: python others/docs/build_gui_doc_links.py
--check   fail when src/dnakit/gui/assets/doc_links.json is out of date with docs/
--site D  also verify every link against a site built by ``mkdocs build --site-dir D``

The graphical interface reads the JSON to link Function explorer rows to the English
documentation. A function links to the heading that names it, otherwise to the section
whose text first mentions or imports it; functions the documentation never mentions get no link.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import unicodedata
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
OUTPUT = ROOT / "src/dnakit/gui/assets/doc_links.json"
BASE_URL = "https://mapengsen.github.io/DNAKit/en/"

# Listings that name many functions without explaining them, and pages that are not documentation.
SKIPPED = re.compile(
    r"^(planning/|changelog|api/features/function_tree|api/features/40_public_tool_coverage"
    r"|faq/user_catalog|agent_tools)"
)
QUALIFIED = re.compile(r"\bdnakit(?:\.[A-Za-z_][A-Za-z0-9_]*)+")
IMPORT = re.compile(r"^\s*from\s+(dnakit(?:\.[A-Za-z_][A-Za-z0-9_]*)*)\s+import\s+(.*)$")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
FENCE = re.compile(r"^\s*(`{3,}|~{3,})")
EXPLICIT_ID = re.compile(r"\{[^}]*#([^\s}]+)[^}]*\}\s*$")
HTML_ID = re.compile(r"""\bid=["']([^"']+)["']""")
# A bare function name in inline code, e.g. `toytree_k80_matrix` in an API table.
BARE = re.compile(r"`([a-z][a-z0-9]*_[a-z0-9_]+)(?:\([^`]*\))?`")
CALL = re.compile(r"\b([a-z][a-z0-9]*_[a-z0-9_]+)\(")


def slugify(value: str) -> str:
    """Python-Markdown's default ``toc`` slug, as used by this MkDocs configuration."""

    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode("ascii")
    value = re.sub(r"[^\w\s-]", "", value).strip().lower()
    return re.sub(r"[-\s]+", "-", value)


def heading_text(raw: str) -> str:
    """The plain text a heading renders to (what the slug is computed from)."""

    raw = re.sub(r"\s*\{[^}]*\}\s*$", "", raw)
    raw = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", raw)
    raw = re.sub(r"<[^>]+>", "", raw)
    raw = raw.replace("`", "").replace("**", "").replace("__", "")
    raw = re.sub(r"(?<!\w)[*_](\S[^*_]*?)[*_](?!\w)", r"\1", raw)
    return html.unescape(raw).strip()


def pages() -> dict[str, Path]:
    """English source of every built page, keyed by its URL path below ``/en/``."""

    found: dict[str, Path] = {}
    for path in sorted(DOCS.rglob("*.md")):
        rel = path.relative_to(DOCS).as_posix()
        if rel.endswith(".en.md"):
            key = rel[: -len(".en.md")]
        else:
            key = rel[: -len(".md")]
            if (DOCS / f"{key}.en.md").exists():
                continue  # the English page wins; the Chinese one is the fallback only
        if SKIPPED.match(key):
            continue
        url = key[: -len("index")] if key == "index" or key.endswith("/index") else key + "/"
        found[url] = path
    return found


def page_rank(url: str) -> int:
    """Prefer tutorials, then parameter references, then everything else."""

    if url.startswith("api/features/"):
        return 0
    if url.startswith("faq/api/"):
        return 1
    return 2


def walk(path: Path) -> Iterator[tuple[str, str | None, str, bool]]:
    """Yield ``(section anchor, heading text or None, line, inside a code fence)`` per line."""

    lines = path.read_text(encoding="utf-8").splitlines()
    used = {m.group(1) for line in lines for m in HTML_ID.finditer(line)}
    current = ""
    fence = ""
    for line in lines:
        marker = FENCE.match(line)
        if marker:
            token = marker.group(1)[0]
            fence = "" if fence == token else (fence or token)
        heading = None if fence else HEADING.match(line)
        if heading is None:
            yield current, None, line, bool(fence)
            continue
        raw = heading.group(2)
        explicit = EXPLICIT_ID.search(raw)
        if explicit:
            current = explicit.group(1)
        else:
            base = slugify(heading_text(raw)) or "_"
            current, n = base, 1
            while current in used:
                current, n = f"{base}_{n}", n + 1
        used.add(current)
        yield current, raw, line, False


def page_anchors(path: Path) -> set[str]:
    """Every heading anchor and explicit HTML id of one page."""

    text = path.read_text(encoding="utf-8")
    return set(HTML_ID.findall(text)) | {anchor for anchor, *_ in walk(path)}


def scan(path: Path) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """Return the anchors of (heading hits, qualified or imported mentions, bare mentions).

    Bare names (``toytree_k80_matrix``) are qualified with every ``dnakit.<module>`` the page
    mentions, so they only count on a page about that module.
    """

    headed: dict[str, str] = {}
    mentioned: dict[str, str] = {}
    bare: dict[str, str] = {}
    pending_import: str | None = None
    for anchor, raw, line, fenced in walk(path):
        if raw is not None:
            names = set(QUALIFIED.findall(raw))
            if anchor.startswith("dnakit."):
                names.add(anchor)
            for name in names:
                headed.setdefault(name, anchor)
            continue
        for name in (CALL if fenced else BARE).findall(line):
            bare.setdefault(name, anchor)
        names = set(QUALIFIED.findall(line))
        if pending_import is not None:
            names.update(f"{pending_import}.{n}" for n in re.findall(r"\b(\w+)\b", line))
            if ")" in line:
                pending_import = None
        imported = IMPORT.match(line)
        if imported:
            module, rest = imported.groups()
            rest = rest.split("#")[0]
            names.update(f"{module}.{n}" for n in re.findall(r"\b(\w+)\b(?!\s+as\b)", rest))
            if "(" in rest and ")" not in rest:
                pending_import = module
        for name in names:
            mentioned.setdefault(name, anchor)
    modules = {".".join(name.split(".")[:2]) for name in [*headed, *mentioned]}
    loose = {
        f"{module}.{name}": anchor
        for name, anchor in bare.items()
        for module in sorted(modules)
        if module != "dnakit"
    }
    return headed, mentioned, loose


def names_section(name: str, title: str) -> bool:
    """Whether a heading describes the function, e.g. "Mark duplicate reads" for
    ``mark_bam_duplicates``: most words of the name (by their first five letters) occur in it."""

    words = [w for w in name.rsplit(".", 1)[-1].lower().split("_") if len(w) >= 3]
    title = title.lower()
    return bool(words) and sum(w[:5] in title for w in words) >= 0.6 * len(words)


def build(tools: list[tuple[str, str]]) -> dict[str, object]:
    """Map registry names (and their implementation names) to ``page/#anchor`` links.

    ``mentions`` lists the functions whose link goes to a section about something else that
    only mentions them (the section heading does not name the function).
    """

    # Strength of the evidence: 0 a heading names it, 1 the text names or imports it,
    # 2 the text uses its bare name on a page about its module.
    found: list[dict[str, tuple[int, str]]] = [{}, {}, {}]
    titles: dict[str, str] = {}
    for url, path in sorted(pages().items(), key=lambda item: (page_rank(item[0]), item[0])):
        rank = page_rank(url)
        for anchor, raw, _, _ in walk(path):
            if raw is not None:
                titles[f"{url}#{anchor}"] = heading_text(raw)
        for strength, hits in enumerate(scan(path)):
            for name, anchor in hits.items():
                found[strength].setdefault(name, (rank, f"{url}#{anchor}" if anchor else url))
    links: dict[str, str] = {}
    mentions: list[str] = []
    kinds: Counter[str] = Counter()
    for name, implementation in tools:
        choices = [
            (*hits[key], strength)
            for strength, hits in enumerate(found)
            for key in (name, implementation)
            if key in hits
        ]
        if not choices:
            kinds["undocumented"] += 1
            continue
        _, link, strength = min(choices, key=lambda c: (c[2] > 0, c[0], c[2]))
        links[name] = link
        title = titles.get(link, "")
        if strength == 0:
            kinds["heading"] += 1
        elif names_section(name, title) or names_section(implementation, title):
            kinds["own_section"] += 1
        else:
            kinds["mention"] += 1
            mentions.append(name)
    return {
        "base_url": BASE_URL,
        "source": "Generated by others/docs/build_gui_doc_links.py from docs/; do not edit.",
        "counts": {"functions": len(tools), **dict(sorted(kinds.items()))},
        "links": dict(sorted(links.items())),
        "mentions": sorted(mentions),
    }


def registry_tools() -> list[tuple[str, str]]:
    from dnakit.tools.registry import default_tool_registry

    registry = default_tool_registry(verbose=False)
    return sorted((tool.name, tool.function_name) for tool in registry.tools)


def verify_site(links: dict[str, str], site: Path) -> list[str]:
    """Return links whose page or anchor is missing from a built site."""

    ids: dict[str, set[str]] = {}
    broken = []
    for name, link in links.items():
        page, _, anchor = link.partition("#")
        if page not in ids:
            built = site / "en" / page / "index.html"
            text = built.read_text(encoding="utf-8") if built.is_file() else None
            ids[page] = set(HTML_ID.findall(text)) if text is not None else set()
            if text is None:
                ids[page].add("\0missing")
        if "\0missing" in ids[page] or (anchor and anchor not in ids[page]):
            broken.append(f"{name}: {link}")
    return broken


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--site", type=Path)
    args = parser.parse_args()
    result = build(registry_tools())
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.site:
        broken = verify_site(result["links"], args.site)  # type: ignore[arg-type]
        if broken:
            sys.exit(f"{len(broken)} links missing from the built site:\n" + "\n".join(broken[:40]))
    if args.check:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != text:
            script = Path(__file__).relative_to(ROOT)
            sys.exit(f"Regenerate {OUTPUT.relative_to(ROOT)}: python {script}")
    else:
        OUTPUT.write_text(text, encoding="utf-8")
    print(json.dumps(result["counts"]))


if __name__ == "__main__":
    main()
