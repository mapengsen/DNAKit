"""Reviewed workflow-specific corrections after mechanical namespace migration."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "others/skills/adaptations"


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    mapping = json.loads((ROOT / "registry/skill-names.json").read_text())["canonical"]
    encode = ast.parse((ROOT / "src/dnakit/encode/__init__.py").read_text())
    encode_names = {node.name for node in encode.body if isinstance(node, ast.FunctionDef)}
    for entry in catalog["skills"]:
        target = TEMPLATES / mapping.get(entry["id"], entry["id"])
        if entry["source_id"] == "encode-toolkit":
            for item in entry["files"]:
                relative = item["path"]
                if (
                    not relative.endswith(".md")
                    or relative == "UPSTREAM.md"
                    or relative.startswith("licenses/")
                ):
                    continue
                path = target / relative
                source = path if path.exists() else ROOT / "skills" / entry["id"] / relative
                text = source.read_text()
                for name in encode_names:
                    text = re.sub(r"\bencode_" + name + r"\b", "dnakit.encode." + name, text)
                text = text.replace(
                    "encode_get_tracking_summary", "dnakit.encode.summarize_collection"
                )
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
    # This skill is a setup workflow, so preserve connection/authentication/test steps
    # while switching its server installation to DNAKit's existing MCP interface.
    path = TEMPLATES / "encode-toolkit-setup/SKILL.md"
    text = path.read_text()
    a, b = text.index("## Installation"), text.index("## Authentication")
    text = (
        text[:a]
        + """## Installation and connection

Install DNAKit in the same Python environment that the Agent will launch:

```bash
python -m pip install "dnakit[agent]"
dnakit skills install encode-toolkit-setup --agent codex
codex mcp add dnakit -- /absolute/path/to/dnakit-mcp
codex mcp list
```

Use the actual executable path from that environment. Restart the client after
changing its MCP configuration. Clients using JSON `mcpServers` configuration
can register the same executable as `{"mcpServers":{"dnakit":{"command":"/absolute/path/to/dnakit-mcp"}}}`.

DNAKit exposes compact MCP discovery tools. Search with `search_dnakit_tools`,
inspect the chosen function with `describe_dnakit_tool`, then execute it with
`call_dnakit_tool`. ENCODE functions are named `dnakit.encode.<operation>`.

Verify a public query before downloads or tracking:

```python
import json
import dnakit

result = json.loads(dnakit.encode.search_experiments(
    assay_title="ATAC-seq", organism="Homo sapiens", limit=3,
))
print(result)
```

Check returned errors and actual experiment accessions; do not treat a registered
MCP configuration as proof that a network query succeeded. Verify filter values
with `dnakit.encode.get_metadata` and data availability with `dnakit.encode.get_facets`.

"""
        + text[b:]
    )
    text = text.replace(
        "description: Set up the ENCODE Toolkit server connection. Use when the user needs help installing, configuring, or troubleshooting the ENCODE connector.",
        "description: Set up ENCODE access through DNAKit, including MCP discovery, public queries, credentials and connection troubleshooting.",
    )
    text = text.replace("Windows Credential Manager, or ", "")
    if "Install `dnakit[encode-credentials]`" not in text:
        text = text.replace(
            "## Authentication\n",
            "## Authentication\n\nInstall `dnakit[encode-credentials]` when using persistent credential storage. DNAKit stores its ENCODE tracker and fallback credentials under `~/.dnakit-encode/`.\n",
        )
    path.write_text(text)
    for path in TEMPLATES.rglob("*"):
        if path.suffix not in {".py", ".md"} or path.name == "UPSTREAM.md":
            continue
        original = path.read_text()
        text = original
        for old, new in {
            "PDB_get_structure": "pdb_get_structure",
            "CADD_get_scores": "cadd_get_scores",
            "Orphanet_get_disease_genes": "orphanet_get_disease_genes",
            "STRING_get_functional_enrichment": "string_get_functional_enrichment",
        }.items():
            text = text.replace("dnakit.databases." + old, "dnakit.databases." + new)
        if text != original:
            path.write_text(text)
    path = TEMPLATES / "gwas-trait-to-gene/python_implementation.py"
    path.write_text(
        path.read_text().replace(
            "from tooluniverse_skills import", "from python_implementation import"
        )
    )
    # The PRS stage was deduplicated as a skill, but its executable dependency is
    # still required by WGS-PRS. It now ships as a DNAKit implementation module.
    path = TEMPLATES / "wgs-prs/wgs_prs.py"
    text = path.read_text()
    a = text.index("    def _find_gwas_prs_script(")
    b = text.index("    def _setup_logging(", a)
    text = (
        text[:a]
        + """    def _find_gwas_prs_script(self) -> Optional[Path]:
        if self._clawbio_root:
            candidate = self._clawbio_root / "skills" / "gwas-prs" / "gwas_prs.py"
            if candidate.is_file():
                return candidate
        import dnakit.workflows.prs
        return Path(dnakit.workflows.prs.__file__).with_name("gwas_prs.py")

"""
        + text[b:]
    )
    path.write_text(text)


if __name__ == "__main__":
    main()
