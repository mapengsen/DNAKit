"""Write the catalogue pages of third-party MCP tools and Agent Skills into docs/llm_tools.

Run from the source checkout: python others/docs/build_agent_catalog_docs.py
--check   fail when the pages or the mkdocs.yml navigation are out of date

One page per MCP server (docs/llm_tools/mcp_tools/<server>.md) lists every indexed DNA/gene
tool of that server with its arguments and the call through DNAKit; one page per Skill
category (docs/llm_tools/skills_catalog/<category>.md) describes every Skill. Both are built
from the data shipped in the package (dnakit.integrations MCP tool index, dnakit.skills
catalog), in Chinese (``.md``) and English (``.en.md``). Tool and Skill descriptions are
upstream English text in both languages. ``targets()`` gives the page and anchor of every
tool and Skill for the graphical interface (see build_gui_doc_links.py).
"""

# Bilingual Markdown prose intentionally keeps full-width punctuation.
# ruff: noqa: RUF001

from __future__ import annotations

import argparse
import html
import json
import pprint
import re
import sys
from collections import Counter
from importlib import resources
from pathlib import Path
from typing import Any

from dnakit.integrations import describe_mcp_tool
from dnakit.skills import get_skill, list_skill_categories, list_skills

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
MCP_DIR = "llm_tools/mcp_tools"
SKILL_DIR = "llm_tools/skills_catalog"
MKDOCS = ROOT / "mkdocs.yml"
MCP_CALL = "../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found"
LANGUAGES = ("zh", "en")
CHECK_ZH = {"passed": "通过", "checksummed": "已校验", "not_run": "未运行"}
# A server with more indexed tools than this gets one page per tool category.
SPLIT_ABOVE = 200


def tr(zh: str, en: str, language: str) -> str:
    return zh if language == "zh" else en


def esc(text: object) -> str:
    """Upstream text as inline Markdown: no HTML, links, attribute lists or headings."""

    value = html.escape(str(text), quote=False)
    return re.sub(r"([\\\[\]|{}#])", r"\\\1", value)


def one_line(text: object, limit: int = 0) -> str:
    value = " ".join(str(text).split())
    if limit and len(value) > limit:
        value = value[: limit - 1].rstrip() + "…"
    return esc(value)


def paragraphs(text: str, indent: str = "") -> list[str]:
    """Upstream multi-line text as escaped paragraphs (docstring indentation removed)."""

    lines = []
    for block in re.split(r"\n\s*\n", text.strip()):
        rows = [row.strip() for row in block.splitlines() if row.strip()]
        if rows:
            lines.extend(indent + esc(row) for row in rows)
            lines.append("")
    return lines


def schema_type(schema: dict[str, Any]) -> str:
    kind = schema.get("type")
    if isinstance(kind, list):
        kind = " | ".join(str(k) for k in kind)
    elif kind is None and isinstance(schema.get("anyOf"), list):
        kind = " | ".join(str(option.get("type", "?")) for option in schema["anyOf"])
    if kind == "array" and isinstance(schema.get("items"), dict) and "type" in schema["items"]:
        kind = f"array of {schema['items']['type']}"
    return str(kind or "any")


def placeholder(schema: dict[str, Any]) -> Any:
    if "default" in schema:
        return schema["default"]
    if schema.get("enum"):
        return schema["enum"][0]
    kind = schema.get("type")
    if isinstance(kind, list):
        kind = next((k for k in kind if k != "null"), "string")
    return {"integer": 0, "number": 0.0, "boolean": False, "array": [], "object": {}}.get(
        str(kind), "..."
    )


# ---------------------------------------------------------------------- MCP tools


def mcp_index() -> dict[str, Any]:
    path = resources.files("dnakit.integrations").joinpath("mcp_tool_index.json")
    data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return data


def page_slug(server_id: str, category: str, sizes: Counter[str]) -> str:
    """File name (without suffix) of the page that lists one tool."""

    return f"{server_id}-{category}" if sizes[server_id] > SPLIT_ABOVE else server_id


def mcp_pages() -> list[dict[str, Any]]:
    """One entry per server page: server facts, its tools, ``slug`` and ``title``."""

    index = mcp_index()
    sizes = Counter(tool["server_id"] for tool in index["tools"])
    pages = []
    for server in mcp_servers():
        groups: dict[str, list[dict[str, Any]]] = {}
        for tool in server["tools"]:
            slug = page_slug(server["server_id"], tool["category"], sizes)
            groups.setdefault(slug, []).append(tool)
        split = len(groups) > 1 or sizes[server["server_id"]] > SPLIT_ABOVE
        for category in index["categories"] if split else [None]:
            slug = page_slug(server["server_id"], category or "", sizes)
            if slug not in groups:
                continue
            title = f"{server['name']}: {category}" if category else str(server["name"])
            pages.append(
                {**server, "tools": groups[slug], "slug": slug, "title": title, "part": category}
            )
    return pages


def mcp_servers() -> list[dict[str, Any]]:
    """Indexed servers with their facts and tools, ordered by display name."""

    index = mcp_index()
    tools: dict[str, list[dict[str, Any]]] = {}
    for tool in index["tools"]:
        tools.setdefault(tool["server_id"], []).append(
            describe_mcp_tool(tool["server_id"], tool["name"], verbose=False)
        )
    servers = []
    for entry in index["servers"]:
        listed = tools.get(entry["server_id"], [])
        if not listed:
            continue
        servers.append(
            {
                **entry,
                **listed[0]["server"],
                "tools": sorted(listed, key=lambda t: t["tool_name"].casefold()),
            }
        )
    return sorted(servers, key=lambda s: (str(s["name"]).casefold(), s["server_id"]))


def mcp_call(tool: dict[str, Any]) -> str:
    properties = (tool.get("input_schema") or {}).get("properties", {})
    own = {name: placeholder(properties.get(name, {})) for name in tool["required_arguments"]}
    call = tool["call"]
    arguments = call["arguments"]
    if isinstance(arguments, dict):  # dispatched through a generic tool of the server
        arguments = {
            key: own if isinstance(value, str) and "arguments" in value else value
            for key, value in arguments.items()
        }
    else:
        arguments = own
    rendered = pprint.pformat(arguments, width=72, sort_dicts=False)
    rendered = rendered.replace("\n", "\n    ")
    return (
        "from dnakit.integrations import call_mcp_tool\n\n"
        f"result = call_mcp_tool(\n    {tool['server_id']!r},\n    {call['tool_name']!r},\n"
        f"    {rendered},\n)"
    )


def source_link(repository: str, commit: str) -> str:
    if not repository:
        return "—"
    label = repository.removeprefix("https://github.com/")
    if commit:
        return f"[{esc(label)} @ `{commit[:7]}`]({repository}/tree/{commit})"
    return f"[{esc(label)}]({repository})"


def mcp_tool_section(tool: dict[str, Any], categories: dict[str, str], language: str) -> list[str]:
    name = tool["tool_name"]
    schema = tool.get("input_schema") or {}
    properties = schema.get("properties", {})
    required = set(tool["required_arguments"])
    scope = esc(categories.get(tool["category"], ""))
    lines = [f"## `{name}` {{#{name}}}", "", one_line(tool["summary"]), ""]
    facts = [
        tr("**类别：**", "**Category:** ", language)
        + f"`{tool['category']}`"
        + tr(f"（{scope}）", f" ({scope})", language),
        tr("**必填参数：**", "**Required arguments:** ", language)
        + (", ".join(f"`{n}`" for n in tool["required_arguments"]) or tr("无", "none", language)),
    ]
    if tool.get("keywords"):
        facts.append(
            tr("**关键词：**", "**Keywords:** ", language)
            + ", ".join(one_line(k) for k in tool["keywords"])
        )
    requires = tool.get("tool_requires") or {}
    if requires.get("api_keys"):
        facts.append(
            tr("**需要 API 密钥：**", "**Needs API keys:** ", language)
            + ", ".join(f"`{k}`" for k in requires["api_keys"])
        )
    hints = tool.get("annotations") or {}
    if hints.get("read_only_hint") is True:
        facts.append(
            tr(
                "**只读：**是（服务器声明）",
                "**Read-only:** yes (declared by the server)",
                language,
            )
        )
    lines += [f"- {fact}" for fact in facts] + [""]
    if properties:
        lines += [
            tr(
                "| 参数 | 类型 | 必填 | 说明 |",
                "| Argument | Type | Required | Description |",
                language,
            ),
            "| --- | --- | --- | --- |",
        ]
        for arg, spec in properties.items():
            spec = spec if isinstance(spec, dict) else {}
            notes = one_line(spec.get("description", ""), 400)
            if spec.get("enum"):
                values = ", ".join(f"`{v}`" for v in spec["enum"][:12])
                notes += tr(f" 可选值：{values}", f" One of {values}", language)
            if "default" in spec:
                notes += tr(" 默认 ", " Default ", language) + f"`{one_line(spec['default'], 60)}`"
            lines.append(
                f"| `{arg}` | {esc(schema_type(spec))} | "
                f"{tr('是', 'yes', language) if arg in required else tr('否', 'no', language)} | "
                f"{notes.strip()} |"
            )
        lines.append("")
    description = str(tool.get("description") or "").strip()
    if description and " ".join(description.split()) != " ".join(str(tool["summary"]).split()):
        lines.append(
            tr('??? note "上游说明（英文原文）"', '??? note "Upstream description"', language)
        )
        lines.append("")
        lines += paragraphs(description, "    ")
    lines += [
        tr(
            "安装服务器后通过 DNAKit 调用：",
            "Call it through DNAKit once the server is installed:",
            language,
        ),
        "",
        "```python",
        mcp_call(tool),
        "```",
        "",
    ]
    return lines


def mcp_server_page(
    server: dict[str, Any],
    siblings: list[dict[str, Any]],
    categories: dict[str, str],
    language: str,
) -> str:
    tools = server["tools"]
    total = sum(len(page["tools"]) for page in siblings)
    env = server.get("required_environment") or []
    secrets = set(server.get("secret_environment") or [])
    command = " ".join([server.get("command") or "", *(server.get("arguments") or [])]).strip()
    lines = [
        f"# {esc(server['title'])} (`{server['server_id']}`)",
        "",
        tr(
            f"DNAKit 的 MCP 工具索引收录了这个服务器的 **{total}** 个 DNA/基因相关工具"
            f"（服务器共列出 {server['listed_tools']} 个，{server['listed_date']} 读取）。"
            "工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。"
            f"调用方式见 [调用索引中的工具]({MCP_CALL})。下面的工具说明是上游英文原文。",
            f"DNAKit's MCP tool index keeps **{total}** DNA/gene tools of this server "
            f"({server['listed_tools']} listed, read on {server['listed_date']}). They run once "
            "you install the server; DNAKit does not bundle them. See "
            f"[calling a tool that was found]({MCP_CALL}).",
            language,
        ),
        "",
    ]
    if server["part"]:
        others = ", ".join(
            f"[`{page['part']}`]({page['slug']}.md) ({len(page['tools'])})"
            for page in siblings
            if page["slug"] != server["slug"]
        )
        lines += [
            tr(
                f"本页列出类别 `{server['part']}` 的 {len(tools)} 个工具，其他类别：{others}。",
                f"This page lists the {len(tools)} tools in category `{server['part']}`; "
                f"other categories: {others}.",
                language,
            ),
            "",
        ]
    lines += [
        tr("| 项目 | 内容 |", "| Item | Value |", language),
        "| --- | --- |",
        tr("| 源码 | ", "| Source | ", language)
        + source_link(server.get("repository") or "", server.get("source_commit") or "")
        + " |",
        tr("| 许可证 | ", "| License | ", language) + esc(server.get("license") or "?") + " |",
        tr("| 启动命令 | ", "| Command | ", language) + (f"`{command}`" if command else "—") + " |",
        tr("| 传输方式 | ", "| Transports | ", language)
        + (", ".join(server.get("transports") or []) or "—")
        + " |",
        tr("| 安装 | ", "| Installation | ", language)
        + one_line(server.get("installation") or "—")
        + " |",
    ]
    if env:
        lines.append(
            tr("| 必需环境变量 | ", "| Required environment | ", language)
            + ", ".join(
                f"`{e}`" + (tr("（保密）", " (secret)", language) if e in secrets else "")
                for e in env
            )
            + " |"
        )
    if server.get("effects"):
        lines.append(
            tr("| 副作用 | ", "| Effects | ", language) + one_line(server["effects"]) + " |"
        )
    lines += [
        "",
        tr("| 工具 | 类别 | 简介 |", "| Tool | Category | Summary |", language),
        "| --- | --- | --- |",
    ]
    for tool in tools:
        lines.append(
            f"| [`{tool['tool_name']}`](#{tool['tool_name']}) | `{tool['category']}` | "
            f"{one_line(tool['summary'], 160)} |"
        )
    lines.append("")
    for tool in tools:
        lines += mcp_tool_section(tool, categories, language)
    return "\n".join(lines).rstrip() + "\n"


def mcp_index_page(servers: list[dict[str, Any]], language: str) -> str:
    """``servers`` are the server pages (a split server appears once per category)."""

    index = mcp_index()
    categories = index["categories"]
    counts = Counter(tool["category"] for tool in index["tools"])
    lines = [
        tr("# 第三方 MCP 工具目录", "# Third-party MCP tool catalogue", language),
        "",
        tr(
            f"DNAKit 收录了 {len(index['servers'])} 个第三方 MCP 服务器的工具清单，从中保留了 "
            f"**{len(index['tools']):,}** 个与 DNA 和基因相关的工具（{index['built_date']} 建立）。"
            "每个服务器一页，逐个列出工具的简介、参数和通过 DNAKit 调用的写法。"
            "检索与调用的函数见 [MCP 工具索引](../../api/features/111_mcp_tool_index.md)，"
            "也可以在 [网页版](https://dnakit.mapengsen.asia/) 的 MCP tools 标签里搜索。",
            f"DNAKit read the tool lists of {len(index['servers'])} third-party MCP servers and "
            f"kept **{len(index['tools']):,}** DNA and gene tools (built {index['built_date']}). "
            "Each server has a page listing every kept tool with its summary, arguments and "
            "the call through DNAKit. The search and call functions are described in the "
            "[MCP tool index](../../api/features/111_mcp_tool_index.md); the "
            "[web demo](https://dnakit.mapengsen.asia/) searches the same tools in its MCP "
            "tools tab.",
            language,
        ),
        "",
        tr(
            "| 服务器 | 工具数 | 许可证 | 源码 |", "| Server | Tools | License | Source |", language
        ),
        "| --- | ---: | --- | --- |",
    ]
    for server in servers:
        lines.append(
            f"| [{esc(server['title'])}]({server['slug']}.md) | {len(server['tools'])} | "
            f"{esc(server.get('license') or '?')} | "
            f"{source_link(server.get('repository') or '', server.get('source_commit') or '')} |"
        )
    empty = sorted(s["server_id"] for s in index["servers"] if not s["indexed_tools"])
    lines += [
        "",
        tr(
            "以下服务器已列出工具，但没有 DNA/基因相关工具被保留：",
            "These servers were listed but have no DNA/gene tool in the index: ",
            language,
        )
        + ", ".join(f"`{s}`" for s in empty)
        + tr("。", ".", language),
        "",
        tr("## 工具类别 {#categories}", "## Tool categories {#categories}", language),
        "",
        tr("| 类别 | 范围 | 工具数 |", "| Category | Scope | Tools |", language),
        "| --- | --- | ---: |",
    ]
    for key, description in categories.items():
        lines.append(f"| `{key}` | {esc(description)} | {counts[key]} |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------- Skills


def skills_by_category() -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    offset = 0
    while True:
        page = list_skills(offset=offset, limit=100)
        for summary in page["skills"]:
            entry = get_skill(summary["id"])
            grouped.setdefault(entry["category"], []).append(entry)
        offset += len(page["skills"])
        if offset >= page["total"] or not page["skills"]:
            break
    return {key: sorted(value, key=lambda e: e["id"]) for key, value in grouped.items()}


def skill_source(entry: dict[str, Any]) -> str:
    url = entry.get("upstream_url") or ""
    label = esc(entry.get("attribution") or entry.get("source_id") or "")
    revision = entry.get("upstream_revision") or ""
    shown = f"{label} @ `{revision[:7]}`" if revision else label
    return f"[{shown}]({url})" if url else shown


def skill_section(entry: dict[str, Any], language: str) -> list[str]:
    skill = entry["id"]
    runtime = entry.get("runtime") or {}
    lines = [f"## `{skill}` {{#{skill}}}", ""]
    lines += paragraphs(str(entry.get("description") or ""))
    facts = [
        tr("**来源：**", "**Source:** ", language) + skill_source(entry),
        tr("**许可证：**", "**License:** ", language) + f"`{esc(entry.get('license') or '?')}`",
    ]
    needs = []
    if [m for m in runtime.get("python_modules") or [] if m != "dnakit"]:
        needs.append(
            tr("Python 模块 ", "Python modules ", language)
            + ", ".join(f"`{m}`" for m in runtime["python_modules"])
        )
    if runtime.get("commands"):
        needs.append(
            tr("命令 ", "commands ", language) + ", ".join(f"`{c}`" for c in runtime["commands"])
        )
    if runtime.get("environment"):
        needs.append(
            tr("环境变量 ", "environment ", language)
            + ", ".join(f"`{e}`" for e in runtime["environment"])
        )
    if needs:
        facts.append(tr("**运行依赖：**", "**Runs with:** ", language) + "; ".join(needs))
    if entry.get("required_skills"):
        facts.append(
            tr("**依赖的 Skills：**", "**Needs Skills:** ", language)
            + ", ".join(f"`{s}`" for s in entry["required_skills"])
        )
    if entry.get("aliases"):
        facts.append(
            tr("**别名：**", "**Also known as:** ", language)
            + ", ".join(f"`{one_line(a)}`" for a in entry["aliases"])
        )
    validation = {
        key: tr(CHECK_ZH.get(str(value), str(value)), str(value), language)
        for key, value in (entry.get("validation") or {}).items()
    }
    facts.append(
        tr("**核对：**", "**Checks:** ", language)
        + tr(
            f"元数据 {validation.get('metadata', '?')}，资源 {validation.get('resources', '?')}，"
            f"执行 {validation.get('execution', '?')}",
            f"metadata {validation.get('metadata', '?')}, resources "
            f"{validation.get('resources', '?')}, execution {validation.get('execution', '?')}",
            language,
        )
    )
    lines += [f"- {fact}" for fact in facts]
    lines += [
        "",
        "```bash",
        f"dnakit skills info {skill}",
        f"dnakit skills install {skill} --agent claude-code   "
        + tr("# 或 --agent codex / cursor", "# or --agent codex / cursor", language),
        "```",
        "",
    ]
    return lines


def skill_category_page(
    category: str, label: dict[str, str], entries: list[dict[str, Any]], language: str
) -> str:
    lines = [
        f"# {esc(label[language])}",
        "",
        tr(
            f"本类共 **{len(entries)}** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，"
            "分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 "
            "[Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。",
            f"**{len(entries)}** Skills. A Skill is a set of instructions that guides an Agent "
            "through an analysis; the software, reference data and APIs it uses are separate "
            "dependencies. See [Skills](../skills.md) for searching, installing and using them; "
            "the Skills tab of the web demo shows each Skill in full.",
            language,
        ),
        "",
        tr("| Skill | 来源 | 简介 |", "| Skill | Source | Summary |", language),
        "| --- | --- | --- |",
    ]
    for entry in entries:
        lines.append(
            f"| [`{entry['id']}`](#{entry['id']}) | {esc(entry.get('source_id') or '')} | "
            f"{one_line(entry.get('description') or '', 160)} |"
        )
    lines.append("")
    for entry in entries:
        lines += skill_section(entry, language)
    return "\n".join(lines).rstrip() + "\n"


def skill_index_page(grouped: dict[str, list[dict[str, Any]]], language: str) -> str:
    labels = list_skill_categories()
    total = sum(len(v) for v in grouped.values())
    sources = Counter(e.get("source_id") for v in grouped.values() for e in v)
    lines = [
        tr("# Skills 目录", "# Skills catalogue", language),
        "",
        tr(
            f"DNAKit 内置 **{total}** 个 DNA 工作流 Skills，来自 {len(sources)} 个来源，按 "
            f"{len(grouped)} 个类别分页列出。"
            "每个 Skill 给出简介、来源与许可证、运行依赖和安装命令。"
            "使用方法见 [Skills](../skills.md)。",
            f"DNAKit bundles **{total}** DNA workflow Skills from {len(sources)} sources, listed "
            f"by category on {len(grouped)} pages. Each entry gives the summary, source and "
            "license, runtime requirements and the install command. See "
            "[Skills](../skills.md) for how to use them.",
            language,
        ),
        "",
        tr("| 类别 | Skills 数 |", "| Category | Skills |", language),
        "| --- | ---: |",
    ]
    for key, value in labels.items():
        if key in grouped:
            lines.append(f"| [{esc(value[language])}]({key}.md) | {len(grouped[key])} |")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------- output


def build() -> tuple[dict[Path, str], dict[str, list[tuple[str, str]]]]:
    """Return the pages to write and the navigation entries ``(title, path)`` per language."""

    output: dict[Path, str] = {}
    nav: dict[str, list[tuple[str, str]]] = {}
    servers = mcp_pages()
    categories = mcp_index()["categories"]
    grouped = skills_by_category()
    labels = list_skill_categories()
    for language in LANGUAGES:
        suffix = ".md" if language == "zh" else ".en.md"
        output[DOCS / MCP_DIR / f"index{suffix}"] = mcp_index_page(servers, language)
        for server in servers:
            siblings = [s for s in servers if s["server_id"] == server["server_id"]]
            output[DOCS / MCP_DIR / f"{server['slug']}{suffix}"] = mcp_server_page(
                server, siblings, categories, language
            )
        output[DOCS / SKILL_DIR / f"index{suffix}"] = skill_index_page(grouped, language)
        for key in labels:
            if key in grouped:
                output[DOCS / SKILL_DIR / f"{key}{suffix}"] = skill_category_page(
                    key, labels[key], grouped[key], language
                )
        overview = tr("总览", "Overview", language)
        nav[f"mcp-tools-{language}"] = [(overview, f"{MCP_DIR}/index.md")] + [
            (f"{s['title']} ({len(s['tools'])})", f"{MCP_DIR}/{s['slug']}.md") for s in servers
        ]
        nav[f"skills-{language}"] = [(overview, f"{SKILL_DIR}/index.md")] + [
            (f"{labels[k][language]} ({len(grouped[k])})", f"{SKILL_DIR}/{k}.md")
            for k in labels
            if k in grouped
        ]
    return output, nav


def update_nav(text: str, nav: dict[str, list[tuple[str, str]]]) -> str:
    """Replace the lines between ``# BEGIN/END GENERATED NAV: <key>`` markers in mkdocs.yml."""

    for key, entries in nav.items():
        pattern = re.compile(
            rf"^(?P<indent>[ ]*)# BEGIN GENERATED NAV: {key}\n.*?^[ ]*# END GENERATED NAV: {key}\n",
            flags=re.M | re.S,
        )
        match = pattern.search(text)
        if match is None:
            raise SystemExit(f"mkdocs.yml has no '# BEGIN GENERATED NAV: {key}' marker")
        indent = match.group("indent")
        body = "".join(
            f"{indent}- {json.dumps(title, ensure_ascii=False)}: {path}\n"
            for title, path in entries
        )
        block = f"{indent}# BEGIN GENERATED NAV: {key}\n{body}{indent}# END GENERATED NAV: {key}\n"
        text = text[: match.start()] + block + text[match.end() :]
    return text


def targets() -> dict[str, dict[str, str]]:
    """Page and anchor (``page/#anchor``) of every MCP tool (``server/tool``) and Skill."""

    tools = mcp_index()["tools"]
    sizes = Counter(tool["server_id"] for tool in tools)
    mcp = {
        f"{tool['server_id']}/{tool['name']}": (
            f"{MCP_DIR}/{page_slug(tool['server_id'], tool['category'], sizes)}/#{tool['name']}"
        )
        for tool in tools
    }
    skills = {}
    offset = 0
    while True:
        page = list_skills(offset=offset, limit=100)
        for summary in page["skills"]:
            skills[summary["id"]] = f"{SKILL_DIR}/{summary['category']}/#{summary['id']}"
        offset += len(page["skills"])
        if offset >= page["total"] or not page["skills"]:
            break
    return {"mcp_tools": dict(sorted(mcp.items())), "skills": dict(sorted(skills.items()))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output, nav = build()
    output[MKDOCS] = update_nav(MKDOCS.read_text(encoding="utf-8"), nav)
    generated = {DOCS / MCP_DIR, DOCS / SKILL_DIR}
    stale_files = [
        path
        for folder in generated
        if folder.is_dir()
        for path in folder.glob("*.md")
        if path not in output
    ]
    if args.check:
        stale = [
            p for p, c in output.items() if not p.is_file() or p.read_text(encoding="utf-8") != c
        ]
        stale += stale_files
        if stale:
            names = ", ".join(str(p.relative_to(ROOT)) for p in stale[:10])
            sys.exit(f"Regenerate the catalogue pages ({len(stale)} stale): {names}")
    else:
        for path in stale_files:
            path.unlink()
        for path, content in output.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    pages = [p for p in output if p != MKDOCS]
    print(json.dumps({"pages": len(pages), "mode": "check" if args.check else "write"}))


if __name__ == "__main__":
    main()
