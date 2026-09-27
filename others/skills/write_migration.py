"""Derive migration from the catalog, native API manifests and maintained templates."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    catalog = json.loads((ROOT / "registry/catalog.json").read_text())
    rows = []
    total = len(catalog["skills"])
    for index, entry in enumerate(catalog["skills"], 1):
        rule = entry.get("adaptation", {})
        functions = set(entry.get("native_functions", []))
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
        if index % 50 == 0 or index == total:
            print(f"[{index}/{total}] Deriving migration records", file=sys.stderr, flush=True)
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
    print(f"[1/1] Wrote migration records for {len(rows)} skills.", flush=True)


if __name__ == "__main__":
    main()
