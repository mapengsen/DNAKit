"""Prepare reviewable DNAKit adaptations for every pinned bundled skill.

This writes maintained templates and adaptation rules. Applying them is a separate
checksum-guarded step; source URLs, original hashes and license files stay intact.
"""

from __future__ import annotations

import ast
import copy
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAMES = json.loads((ROOT / "registry/skill-names.json").read_text())["canonical"]
PROVIDERS = json.loads((ROOT / "src/dnakit/databases/providers.json").read_text())
OPERATIONS = {item["upstream_name"]: name for name, item in PROVIDERS.items()}
NAME_PATTERN = re.compile(
    r"(?<![\w-])("
    + "|".join(re.escape(old) for old in sorted(NAMES, key=len, reverse=True) if old != NAMES[old])
    + r")(?![\w-])"
)
API_PATTERN = re.compile(r"\b(?:tu|self\.tu)\.tools\.([A-Za-z_][A-Za-z0-9_]*)")


def references(text: str) -> str:
    # Upstream URLs identify the original skill and must never be renamed.
    return "".join(
        part
        if part.startswith(("https://", "http://"))
        else NAME_PATTERN.sub(lambda m: NAMES[m[0]], part)
        for part in re.split(r"(https?://[^\s<>\]\)]+)", text)
    )


def python(text: str) -> str:
    text = text.replace(
        "from tooluniverse.tools.execute_tool import execute_tool",
        "from dnakit.databases import call as execute_tool",
    )
    text = text.replace("from tooluniverse import ToolUniverse", "import dnakit")
    text = re.sub(r"ToolUniverse\([^\n)]*\)", "dnakit.databases", text)
    text = re.sub(
        r"^([ \t]*)(?:tu|self\.tu)\.load_tools\([^\n]*\)[ \t]*(?:#[^\n]*)?$",
        r"\1# DNAKit loads the selected implementation on demand.",
        text,
        flags=re.M,
    )
    text = API_PATTERN.sub(lambda m: "dnakit.databases." + OPERATIONS.get(m[1], m[1]), text)
    text = re.sub(
        r"\btu\.all_tool_dict\b",
        "{item['upstream_name']: item for item in dnakit.databases.list_operations()}",
        text,
    )
    text = re.sub(r"\btu\.(?:run_one_function|run_tool)\b", "dnakit.databases.execute", text)
    # Named calls with keyword arguments become explicit DNAKit callables.
    text = re.sub(
        r"\btu\.(?:execute_tool|use_tools)\(\s*(['\"])(\w+)\1\s*,?\s*",
        lambda m: "dnakit.databases." + OPERATIONS.get(m[2], m[2]) + "(",
        text,
    )
    text = text.replace("from tooluniverse.tools import (", "from dnakit.databases import (")
    # Rewrite imported function identifiers consistently throughout this snippet.
    if "from dnakit.databases import (" in text:
        for block in re.findall(r"from dnakit\.databases import \((.*?)\)", text, re.S):
            for name in re.findall(r"\b\w+\b", block):
                if name in OPERATIONS and name != OPERATIONS[name]:
                    text = re.sub(r"\b" + re.escape(name) + r"\b", OPERATIONS[name], text)
    # Collection-specific bypass flags are not scientific inputs; native schemas are validated.
    if "dnakit.databases" in text:
        text = re.sub(r"(?m)^[ \t]*validate=False,?[ \t]*(?:#[^\n]*)?\n", "", text)
        text = re.sub(r",[ \t]*validate=False(?=[ \t]*\))", "", text)
        text = re.sub(
            r"(?m)^([ \t]*)from tooluniverse_skills import discover_gwas_genes",
            r"\1from python_implementation import discover_gwas_genes",
            text,
        )
        if not re.search(r"(?:^|\n)(?:import dnakit|from dnakit import databases)", text):
            # Keep future imports and a module docstring at the beginning.
            try:
                tree = ast.parse(text)
                last = 0
                for node in tree.body:
                    if (
                        isinstance(node, ast.Expr)
                        and isinstance(node.value, ast.Constant)
                        and isinstance(node.value.value, str)
                    ) or (isinstance(node, ast.ImportFrom) and node.module == "__future__"):
                        last = node.end_lineno
                    else:
                        break
                lines = text.splitlines(keepends=True)
                lines.insert(last, "import dnakit\n")
                text = "".join(lines)
            except SyntaxError:
                text = "import dnakit\n" + text
    return text


def transform(text: str, suffix: str) -> str:
    text = references(text)
    text = text.replace("clawbio.common", "dnakit.workflows.common")
    text = text.replace("encode_connector.client", "dnakit.encode._backend.client")
    text = re.sub(r"\btu\.run_tool\b", "dnakit.databases.execute", text)
    text = re.sub(r"\btu run\s+", "dnakit database call ", text)
    text = re.sub(r"\bpip(?:3)? install tooluniverse\b", "python -m pip install dnakit", text)
    text = text.replace("uv pip install tooluniverse", "uv pip install dnakit")
    text = text.replace(
        "tu.run_batch()", "concurrent.futures.ThreadPoolExecutor.map() with DNAKit functions"
    )
    if suffix == ".py":
        return python(text)
    if suffix == ".md":
        text = re.sub(
            r"(```(?:python|py)\s*\n)(.*?)(\n```)",
            lambda m: m[1] + python(m[2]) + m[3],
            text,
            flags=re.S,
        )
    return text


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    if all((ROOT / "skills" / row["id"] / "DNAKIT.md").exists() for row in catalog["skills"]):
        raise SystemExit(
            "Templates already exist for the full catalog. Edit the maintained templates and "
            "run adaptations.py; do not regenerate them from already adapted skill files."
        )
    existing = {
        rule["id"]: rule
        for rule in json.loads((ROOT / "registry/adaptations.json").read_text())["adaptations"]
    }
    rules, matrix = [], []
    syntax_issues = []
    for entry in catalog["skills"]:
        old, new = entry["id"], NAMES.get(entry["id"], entry["id"])
        template = ROOT / "others/skills/adaptations" / new
        template.mkdir(parents=True, exist_ok=True)
        changed, scripts, referenced = [], [], set()
        for file in entry["files"]:
            relative = file["path"]
            if relative == "UPSTREAM.md" or relative.startswith("licenses/"):
                continue
            path = ROOT / "skills" / old / relative
            if path.suffix not in {
                ".md",
                ".py",
                ".json",
                ".yaml",
                ".yml",
                ".toml",
                ".sh",
                ".txt",
                ".R",
            }:
                continue
            original = path.read_text()
            referenced.update(set(re.findall(r"\b[A-Za-z_]\w*\b", original)) & OPERATIONS.keys())
            text = transform(original, path.suffix)
            if path.suffix.lower() in {".py", ".r", ".sh"} and not any(
                part.startswith("test") for part in Path(relative).parts
            ):
                scripts.append(relative)
            if path.suffix == ".py" and text != original:
                try:
                    ast.parse(text)
                except SyntaxError as exc:
                    syntax_issues.append(
                        {"skill": new, "file": relative, "line": exc.lineno, "message": exc.msg}
                    )
            if relative == "SKILL.md":
                text = re.sub(r"(?m)^name:.*$", "name: " + new, text, count=1)
                closing = text.find("\n---", 4) + len("\n---")
                text = (
                    text[:closing]
                    + (
                        "\n\n## DNAKit execution\n\n"
                        f"Use `dnakit skills doctor {new}` to inspect scientific dependencies. "
                        "Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. "
                        "Keep the workflow, input checks and reporting requirements below.\n"
                    )
                    + text[closing:]
                )
            if text != original:
                output = template / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text(text)
                changed.append(relative)
        runtime = copy.deepcopy(entry.get("runtime", {}))
        modules = runtime.get("python_modules", [])
        runtime["python_modules"] = sorted(
            {
                "dnakit",
                *[
                    m
                    for m in modules
                    if m.split(".")[0] not in {"tooluniverse", "clawbio", "encode_connector"}
                ],
            }
        )
        runtime["commands"] = [
            c
            for c in runtime.get("commands", [])
            if c not in {"tooluniverse", "tu", "clawbio", "encode-toolkit"}
        ]
        kind = (
            "dnakit-database"
            if referenced
            else ("dnakit-bundled-workflow" if scripts else "scientific-backend-retained")
        )
        if entry["source_id"] == "encode-toolkit":
            kind = "dnakit-encode"
        guide = [
            "# DNAKit integration",
            "",
            f"Canonical skill: `{new}`.",
            "",
            "```bash",
            "python -m pip install dnakit",
            f"dnakit skills doctor {new}",
            f"dnakit skills install {new} --agent codex",
            "```",
            "",
            "Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill "
            "and its listed local resources. Real scientific libraries, command-line programs, models, "
            "reference databases and service credentials retain their own requirements.",
            "",
        ]
        if referenced:
            guide.extend(
                [
                    "## Database functions",
                    "",
                    "```python",
                    "import dnakit",
                    "# Inspect the exact parameters before calling a service:",
                    f"schema = dnakit.databases.describe('{OPERATIONS[sorted(referenced)[0]]}')",
                    "```",
                    "",
                    "DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. "
                    "Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.",
                    "",
                    "| Upstream operation | DNAKit function |",
                    "| --- | --- |",
                ]
            )
            guide.extend(
                f"| `{name}` | `dnakit.databases.{OPERATIONS[name]}` |"
                for name in sorted(referenced)
            )
        if scripts:
            guide.extend(
                [
                    "",
                    "## Bundled scripts",
                    "",
                    "Find scripts and run the relevant one with its documented arguments:",
                    "",
                    "```bash",
                    f"dnakit skills scripts {new}",
                    f"dnakit skills run {new} --script {scripts[0]} -- <documented arguments>",
                    "```",
                    "",
                    "The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. "
                    "It forwards standard input/output and exit status. Python callers can use "
                    "`dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.",
                    "",
                ]
            )
            guide.extend(f"- `{script}`" for script in scripts)
        if entry["source_id"] == "encode-toolkit":
            guide.extend(
                [
                    "",
                    "## ENCODE",
                    "",
                    "Use `dnakit.encode` for ENCODE queries, downloads and provenance tracking. "
                    "The API retains the original options; `encode_search_experiments(...)` is now "
                    "`dnakit.encode.search_experiments(...)`. It returns the original JSON text. "
                    "Credential storage uses `dnakit[encode-credentials]`; public searches need no credentials.",
                    "",
                ]
            )
        guide.extend(
            [
                "",
                "## Provenance and validation",
                "",
                "Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). "
                "DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks "
                "do not establish scientific validity or completion of the full workflow.",
                "",
            ]
        )
        (template / "DNAKIT.md").write_text("\n".join(guide))
        upstream_path = ROOT / "skills" / old / "UPSTREAM.md"
        upstream = upstream_path.read_text() if upstream_path.exists() else "# DNAKit-owned skill\n"
        (template / "UPSTREAM.md").write_text(
            upstream
            + "\n\n## DNAKit adaptation — 2026-09-13\n\n"
            + f"Distributed as `{new}`; previous identifier `{old}` remains an alias when different. "
            + "DNAKit changed discovery links, package calls and execution guidance. Original authorship, "
            + "licenses, citations and original file hashes are retained. "
            + "Scientific software/models remain their original implementations. See DNAKIT.md and registry/adaptations.json.\n"
        )
        rule = {
            "id": new,
            "original_id": existing.get(old, {}).get("original_id", old),
            "source_id": entry["source_id"],
            "upstream_revision": entry["upstream_revision"],
            "date": "2026-09-13",
            "implementation": kind,
            "reason": "Use DNAKit discovery, bundled workflows and native collection adapters; retain real scientific backends and upstream attribution.",
            "runtime": runtime,
        }
        if old == "sequence-retrieval":
            rule["implementation"] = "dnakit-native"
        rules.append(rule)
        matrix.append(
            {
                "id": new,
                "original_id": old,
                "source": entry["source_id"],
                "integration": rule["implementation"],
                "native_functions": [OPERATIONS[name] for name in sorted(referenced)],
                "scripts": scripts,
                "retained_runtime": runtime,
                "changed_resources": changed,
                "execution_validation": entry["validation"]["execution"],
            }
        )
    (ROOT / "registry/adaptations.json").write_text(
        json.dumps({"schema_version": 1, "adaptations": rules}, indent=2, ensure_ascii=False) + "\n"
    )
    (ROOT / "registry/migration.json").write_text(
        json.dumps(
            {"schema_version": 1, "skills": matrix, "syntax_issues": syntax_issues},
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    )
    print("Prepared", len(rules), "adaptations; changed-Python syntax issues:", syntax_issues)


if __name__ == "__main__":
    main()
