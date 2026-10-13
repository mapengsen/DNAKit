"""Release archives must retain exact source and cataloged resource bytes."""

from __future__ import annotations

import hashlib
import io
import json
import tarfile
import zipfile
from pathlib import Path

import pytest
from others.validation.release_archives import audit


def _archives(root: Path, *, tamper: str | None = None) -> Path:
    skill = b"# A cataloged workflow\n"
    catalog = json.dumps(
        {
            "skills": [
                {
                    "id": "example",
                    "files": [
                        {
                            "path": "SKILL.md",
                            "sha256": hashlib.sha256(skill).hexdigest(),
                        }
                    ],
                }
            ]
        }
    ).encode()
    source = {
        "src/dnakit/__init__.py": b"",
        "src/dnakit/py.typed": b"",
        "LICENSE": b"MIT\n",
        "docs/planning/requirements_traceability.csv": b"id,requirement\n",
        "registry/catalog.json": catalog,
        "skills/example/SKILL.md": skill,
    }
    for name, data in source.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    dist = root / "dist"
    dist.mkdir()
    wheel = {
        "dnakit/__init__.py": b"",
        "dnakit/py.typed": b"",
        "dnakit/skills/_data/catalog.json": catalog,
        "dnakit/skills/_data/example/SKILL.md": skill,
        "dnakit-0.2.0.dist-info/licenses/LICENSE": b"MIT\n",
    }
    if tamper == "skill":
        wheel["dnakit/skills/_data/example/SKILL.md"] = b"changed"
    elif tamper == "missing":
        del wheel["dnakit/skills/_data/example/SKILL.md"]
    elif tamper == "code":
        wheel["dnakit/__init__.py"] = b"changed"
    elif tamper == "cache":
        wheel["dnakit/__pycache__/local.pyc"] = b"cache"
    with zipfile.ZipFile(dist / "dnakit-0.2.0-py3-none-any.whl", "w") as archive:
        for name, data in wheel.items():
            archive.writestr(name, data)
    if tamper == "sdist":
        source["skills/example/SKILL.md"] = b"changed"
    with tarfile.open(dist / "dnakit-0.2.0.tar.gz", "w:gz") as archive:
        for name, data in source.items():
            member = tarfile.TarInfo(f"dnakit-0.2.0/{name}")
            member.size = len(data)
            if tamper == "mode" and name == "skills/example/SKILL.md":
                member.mode = 0o755
            archive.addfile(member, io.BytesIO(data))
    return dist


def test_release_archives_preserve_catalog_and_source(tmp_path: Path) -> None:
    audit(_archives(tmp_path), root=tmp_path)


@pytest.mark.parametrize("tamper", ["skill", "missing", "code", "cache", "sdist", "mode"])
def test_release_archives_reject_corruption(tmp_path: Path, tamper: str) -> None:
    with pytest.raises(AssertionError):
        audit(_archives(tmp_path, tamper=tamper), root=tmp_path)
