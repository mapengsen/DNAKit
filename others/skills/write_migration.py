"""Derive the migration matrix from the current catalog and maintained templates."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    rows = []
    for entry in catalog["skills"]:
        rule = entry.get("adaptation", {})
        functions = set()
        for path in (ROOT / "others/skills/adaptations" / entry["id"]).rglob("*"):
            if path.suffix in {".py", ".md"} and path.name != "UPSTREAM.md":
                functions.update(
                    re.findall(r"dnakit\.(?:databases|encode|workflows)\.\w+", path.read_text())
                )
        rows.append(
            {
                "id": entry["id"],
                "original_id": rule.get("original_id", entry["id"]),
                "source": entry["source_id"],
                "upstream_revision": entry["upstream_revision"],
                "integration": rule.get("implementation", "scientific-backend-retained"),
                "native_functions": sorted(functions),
                "scripts": [
                    item["path"]
                    for item in entry["files"]
                    if Path(item["path"]).suffix.lower() in {".py", ".r", ".sh"}
                    and not any(p.startswith("test") for p in Path(item["path"]).parts)
                ],
                "retained_runtime": entry.get("runtime", {}),
                "changed_resources": entry.get("modifications", []),
                "execution_validation": entry["validation"]["execution"],
            }
        )
    (ROOT / "registry/migration.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "skills": rows,
                "note": (
                    "Integration, file and interface validation do not establish "
                    "full scientific execution."
                ),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
