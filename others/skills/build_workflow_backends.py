"""Import attributed workflow support from reviewed, pinned source archives."""

from __future__ import annotations

import ast
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    sources = json.loads((ROOT / "registry/sources.json").read_text())["sources"]
    for source_id, source_prefix, destination in (
        ("clawbio", "clawbio/common/", "src/dnakit/workflows/common"),
        ("clawbio", "skills/gwas-prs/", "src/dnakit/workflows/prs"),
        ("encode-toolkit", "src/encode_connector/", "src/dnakit/encode/_backend"),
    ):
        source = next(item for item in sources if item["id"] == source_id)
        revision = source["revision"]
        command = ["git", "-C", str(ROOT / source["cache_path"])]

        def read(path: str, command=command, revision=revision) -> bytes:
            return subprocess.check_output([*command, "show", revision + ":" + path])

        paths = subprocess.check_output(
            [*command, "ls-tree", "-r", "--name-only", revision], text=True
        ).splitlines()
        target = ROOT / destination
        target.mkdir(parents=True, exist_ok=True)
        manifest = []
        for path in paths:
            if (
                not path.startswith(source_prefix)
                or not path.endswith((".py", ".json", ".txt"))
                or "/tests/" in path
            ):
                continue
            relative = path.removeprefix(source_prefix)
            if source_id == "encode-toolkit" and relative.endswith("__main__.py"):
                continue
            original = read(path)
            content = original.decode()
            if source_prefix == "skills/gwas-prs/":
                content = content.replace(
                    "from repro_bundle import", "from dnakit.workflows.prs.repro_bundle import"
                )
            content = content.replace("clawbio.common", "dnakit.workflows.common")
            content = content.replace("encode_connector.", "dnakit.encode._backend.")
            content = (
                content.replace('".encode_connector"', '".dnakit-encode"')
                .replace('"encode-connector"', '"dnakit-encode"')
                .replace('".clawbio"', '".dnakit"')
            )
            if source_id == "clawbio" and relative == "__init__.py":
                exports = {}
                for node in ast.parse(content).body:
                    if isinstance(node, ast.ImportFrom):
                        for alias in node.names:
                            if node.module:
                                exports[alias.asname or alias.name] = (node.module, alias.name)
                content = (
                    '"""Lazy shared scientific workflow helpers; see UPSTREAM.json."""\n'
                    "from importlib import import_module\n"
                    f"_EXPORTS = {exports!r}\n"
                    "__all__ = list(_EXPORTS)\n"
                    "def __getattr__(name):\n"
                    "    if name not in _EXPORTS:\n"
                    "        raise AttributeError(name)\n"
                    "    module, attribute = _EXPORTS[name]\n"
                    '    if not module.startswith("dnakit."):\n'
                    '        module = __name__ + "." + module\n'
                    "    return getattr(import_module(module), attribute)\n"
                )
            if source_id == "encode-toolkit" and relative == "server/main.py":
                lines = content.splitlines(keepends=True)
                edits = []
                for node in ast.parse(content).body:
                    remove = isinstance(node, ast.ImportFrom) and (node.module or "").startswith(
                        "mcp"
                    )
                    if isinstance(node, ast.Assign):
                        remove |= any(
                            isinstance(t, ast.Name)
                            and (
                                t.id == "mcp"
                                or t.id
                                in {
                                    "_READONLY_API",
                                    "_READONLY_LOCAL",
                                    "_WRITE_LOCAL",
                                    "_WRITE_LOCAL_API",
                                    "_DOWNLOAD",
                                    "_CREDENTIAL_MGMT",
                                }
                            )
                            for t in node.targets
                        )
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        remove |= node.name == "main"
                        for decorator in node.decorator_list:
                            if ast.unparse(decorator).startswith("mcp.tool("):
                                edits.append((decorator.lineno - 1, decorator.end_lineno, ""))
                    if isinstance(node, ast.If) and "__name__" in ast.unparse(node.test):
                        remove = True
                    if remove:
                        edits.append((node.lineno - 1, node.end_lineno, ""))
                for start, stop, replacement in sorted(edits, reverse=True):
                    lines[start:stop] = [replacement]
                content = "".join(lines).replace(
                    "lifespan(server: FastMCP)", "lifespan(server=None)"
                )
                # Local tracker creation belongs to tracker operations, not read-only HTTP queries.
                content = content.replace(
                    "    _tracker = ExperimentTracker()\n    try:", "    _tracker = None\n    try:"
                )
                content = content.replace(
                    '    if _tracker is None:\n        raise RuntimeError("Experiment tracker not initialized")',
                    "    global _tracker\n    if _tracker is None:\n        _tracker = ExperimentTracker()",
                )
            notice = (
                f"# Adapted from {source['repository']} @ {revision}.\n"
                f"# Original file: {path}. See LICENSE and UPSTREAM.json at the package root.\n"
                "# DNAKit changes: package namespace, isolated cache, lazy imports and native execution.\n"
            )
            output = target / relative
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text((notice if path.endswith(".py") else "") + content)
            manifest.append(
                {
                    "path": relative,
                    "upstream_path": path,
                    "sha256": hashlib.sha256(original).hexdigest(),
                    "adapted_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                }
            )
        if source_prefix == "skills/gwas-prs/":
            (target / "__init__.py").write_text(
                '"""DNAKit polygenic scoring helpers; adapted from ClawBio under MIT."""\n'
            )
        (target / "LICENSE").write_bytes(read("LICENSE"))
        (target / "UPSTREAM.json").write_text(
            json.dumps(
                {
                    "source": source["repository"],
                    "revision": revision,
                    "attribution": source["attribution"],
                    "license": source["content_license"],
                    "files": manifest,
                },
                indent=2,
            )
            + "\n"
        )
        print(source_id, len(manifest), "implementation modules")
        if source_id == "encode-toolkit":
            functions = [
                node
                for node in ast.parse((target / "server/main.py").read_text()).body
                if isinstance(node, ast.AsyncFunctionDef) and node.name.startswith("encode_")
            ]
            output = [
                '"""Native ENCODE query, download, provenance and tracking APIs."""',
                "from __future__ import annotations",
                "from typing import Any, Literal",
                "from ._runtime import invoke as _invoke",
                "",
            ]
            for node in functions:
                name = node.name.removeprefix("encode_")
                parameters = ast.unparse(node.args)
                parameters = parameters + ", " if parameters else ""
                output.extend(
                    [
                        "",
                        f"def {name}({parameters}*, verbose: bool = True) -> str:",
                        "    " + repr(ast.get_docstring(node) or name),
                        "    arguments = dict(locals())",
                        '    arguments.pop("verbose")',
                        f"    return _invoke({name!r}, arguments)",
                    ]
                )
            output.append(
                "\n__all__ = " + repr([node.name.removeprefix("encode_") for node in functions])
            )
            output.extend(
                [
                    "from dnakit._verbosity import register_public as _register_public_verbose",
                    "_register_public_verbose(globals())",
                ]
            )
            (target.parent / "__init__.py").write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    main()
