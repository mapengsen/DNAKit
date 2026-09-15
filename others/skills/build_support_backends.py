"""Bundle licensed dependencies of composite skills, without adding duplicate skills."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    source = next(
        item
        for item in json.loads((ROOT / "registry/sources.json").read_text())["sources"]
        if item["id"] == "clawbio"
    )
    command = ["git", "-C", str(ROOT / source["cache_path"])]
    revision = source["revision"]

    def read(path):
        return subprocess.check_output([*command, "show", revision + ":" + path])

    root = ROOT / "src/dnakit/workflows/support"
    root.mkdir(parents=True, exist_ok=True)
    (root / "__init__.py").write_text(
        '"""Licensed support implementations for composite DNAKit workflows."""\n'
    )
    files = subprocess.check_output(
        [*command, "ls-tree", "-r", "--name-only", revision], text=True
    ).splitlines()
    manifest = []
    for name in (
        "clinical-variant-reporter",
        "clinical-trial-finder",
        "omics-target-evidence-mapper",
        "gwas-lookup",
        "fine-mapping",
    ):
        prefix = "skills/" + name + "/"
        document = read(prefix + "SKILL.md").decode()
        metadata = yaml.safe_load(document.split("---", 2)[1])
        if metadata.get("license") != "MIT":
            raise ValueError(f"Review the per-skill license before importing {name}")
        destination = root / name.replace("-", "_")
        destination.mkdir(exist_ok=True)
        for path in files:
            if (
                not path.startswith(prefix)
                or "/tests/" in path
                or Path(path).suffix not in {".py", ".json", ".txt", ".csv", ".tsv", ".vcf", ".md"}
            ):
                continue
            original = read(path)
            text = original.decode()
            text = text.replace("clawbio.common", "dnakit.workflows.common")
            text = text.replace('".clawbio"', '".dnakit"')
            text = text.replace(
                "from datetime import UTC, datetime",
                "from datetime import timezone, datetime\nUTC = timezone.utc",
            )
            target = destination / path.removeprefix(prefix)
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.suffix == ".py":
                text = (
                    f"# Adapted from {source['repository']} @ {revision}, {path}, MIT.\n"
                    "# DNAKit changes: package namespace. See ../LICENSE and ../UPSTREAM.json.\n"
                    + text
                )
            target.write_text(text)
            manifest.append(
                {
                    "path": str(target.relative_to(root)),
                    "original_path": path,
                    "sha256": hashlib.sha256(original).hexdigest(),
                    "adapted_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                }
            )
        (destination / "__init__.py").write_text(
            '"""Support workflow; see package provenance."""\n'
        )
    path = "GENOMEBOOK/PYTHON/04-recombinator.py"
    original = read(path)
    text = original.decode().replace(
        'DATA = BASE / "DATA"', 'DATA = Path.cwd() / "genomebook-data"'
    )
    text = text.replace(
        "def load_registries():\n    with open(TRAIT_REGISTRY) as f:",
        "def load_registries(trait_path=None, disease_path=None):\n    with open(trait_path or TRAIT_REGISTRY) as f:",
    )
    text = text.replace(
        "    with open(DISEASE_REGISTRY) as f:",
        "    with open(disease_path or DISEASE_REGISTRY) as f:",
    )
    text = text.replace(
        'shared_loci = set(parent_a["loci"].keys()) & set(parent_b["loci"].keys())',
        'shared_loci = sorted(set(parent_a["loci"].keys()) & set(parent_b["loci"].keys()))',
    )
    target = root / "recombinator.py"
    target.write_text(
        f"# Adapted from {source['repository']} @ {revision}, {path}, MIT.\n"
        "# DNAKit changes: explicit optional registry paths; runtime data outside installed source.\n"
        + text
    )
    manifest.append(
        {
            "path": "recombinator.py",
            "original_path": path,
            "sha256": hashlib.sha256(original).hexdigest(),
            "adapted_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        }
    )
    (root / "LICENSE").write_bytes(read("LICENSE"))
    (root / "UPSTREAM.json").write_text(
        json.dumps(
            {
                "source": source["repository"],
                "revision": revision,
                "license": "MIT",
                "attribution": source["attribution"],
                "files": manifest,
            },
            indent=2,
        )
        + "\n"
    )
    print("Bundled", len(manifest), "licensed composite-workflow resources")


if __name__ == "__main__":
    main()
