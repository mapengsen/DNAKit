"""Normalize executable examples and runtime metadata in maintained skill templates."""

from __future__ import annotations

import ast
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "others/skills/adaptations"


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    entries = {entry["id"]: entry for entry in catalog["skills"]}
    names = json.loads((ROOT / "registry/skill-names.json").read_text())["canonical"]
    source = subprocess.check_output(
        [
            "git",
            "-C",
            str(ROOT / "temp/skill-sources/ClawBio"),
            "show",
            "1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d:clawbio/cli.py",
        ],
        text=True,
    )
    aliases = {}
    for node in ast.parse(source).body:
        if not isinstance(node, ast.Assign) or not any(
            isinstance(t, ast.Name) and t.id == "SKILLS" for t in node.targets
        ):
            continue
        for key, value in zip(node.value.keys, node.value.values, strict=True):
            for field, item in zip(value.keys, value.values, strict=True):
                if isinstance(field, ast.Constant) and field.value == "script":
                    parts = re.findall(r"['\"]([^'\"]+)['\"]", ast.unparse(item))
                    if len(parts) >= 2:
                        skill = names.get(parts[0], parts[0])
                        if skill in entries:
                            aliases[key.value] = (skill, "/".join(parts[1:]))

    def script_command(match):
        skill, resource = match[1], match[2]
        canonical = names.get(skill, skill)
        entry = entries.get(canonical)
        if entry and resource in {item["path"] for item in entry["files"]}:
            return f"dnakit skills run {canonical} --script {resource} --"
        return match[0]

    def claw_command(match):
        if match[1] in {"gwas-lookup", "fine-mapping"}:
            return f"dnakit skills support-run {match[1]} --"
        pair = aliases.get(match[1])
        if not pair and match[1] in entries:
            candidate = match[1].replace("-", "_") + ".py"
            if candidate in {row["path"] for row in entries[match[1]]["files"]}:
                pair = (match[1], candidate)
        if not pair:
            raise ValueError(f"Unresolved workflow alias: {match[1]}")
        return f"dnakit skills run {pair[0]} --script {pair[1]} --"

    changed = []
    for skill, entry in entries.items():
        for item in entry["files"]:
            relative = item["path"]
            if (
                Path(relative).suffix not in {".md", ".py", ".sh"}
                or relative == "UPSTREAM.md"
                or relative.startswith("licenses/")
                or any(part.startswith("test") for part in Path(relative).parts)
            ):
                continue
            output = TEMPLATES / skill / relative
            path = output if output.exists() else ROOT / "skills" / skill / relative
            text = original = path.read_text()
            text = re.sub(
                r"\bpython(?:3)?\s+skills/([\w-]+)/([\w./-]+\.(?:py|sh|R))\b",
                script_command,
                text,
            )
            text = re.sub(r"\bpython(?:3)?\s+clawbio\.py\s+run\s+([\w-]+)", claw_command, text)
            if text != original:
                if path.suffix == ".py":
                    ast.parse(text)
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text(text)
                changed.append(f"{skill}/{relative}")
    path = ROOT / "registry/adaptations.json"
    manifest = json.loads(path.read_text())
    providers = json.loads((ROOT / "src/dnakit/databases/providers.json").read_text())
    optional = {
        "chem_tool": ["indigo"],
        "dataset_tool": ["pandas"],
        "variant_fraction_tool": ["pandas"],
        "enrichr_tool": ["networkx"],
        "humanbase_tool": ["networkx"],
        "xml_tool": ["lxml"],
        "cellxgene_census_tool": ["cellxgene_census"],
        "cellxgene_discovery_tool": ["cellxgene_census"],
    }
    for rule in manifest["adaptations"]:
        runtime = rule.get("runtime", {})
        modules = set(runtime.get("python_modules", []))
        for resource in (TEMPLATES / rule["id"]).rglob("*"):
            if resource.suffix in {".md", ".py"} and resource.name != "UPSTREAM.md":
                for name in re.findall(r"dnakit\.databases\.(\w+)", resource.read_text()):
                    modules.update(optional.get(providers.get(name, {}).get("module"), []))
        runtime["python_modules"] = sorted(modules)
        if rule["source_id"] == "encode-toolkit":
            runtime["manual_checks"] = [
                note
                for note in runtime.get("manual_checks", [])
                if "Configure the separate encode-toolkit" not in note
            ]
            notice = "ENCODE adapters run in DNAKit; AGPL-3.0-only attribution is retained."
            if notice not in runtime["manual_checks"]:
                runtime["manual_checks"].append(notice)
        rule["runtime"] = runtime
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(f"Normalized {len(changed)} execution documents and scripts.")


if __name__ == "__main__":
    main()
