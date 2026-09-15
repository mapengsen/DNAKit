"""Bundle the canonical root skills directory in wheels without duplicating source files."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = Path(__file__).resolve().parent


class BuildWithSkills(build_py):
    """Copy only cataloged, checksum-verified resources into the installed package."""

    def _skill_outputs(self) -> list[tuple[Path, Path, str | None]]:
        catalog_path = ROOT / "registry" / "catalog.json"
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
        target = Path(self.build_lib) / "dnakit" / "skills" / "_data"
        files = [(catalog_path, target / "catalog.json", None)]
        for skill in catalog["skills"]:
            for resource in skill["files"]:
                relative = Path(skill["id"]) / resource["path"]
                files.append((ROOT / "skills" / relative, target / relative, resource["sha256"]))
        return files

    def run(self) -> None:
        super().run()
        if self.editable_mode:
            return
        bundle = Path(self.build_lib) / "dnakit" / "skills" / "_data"
        if bundle.exists():
            shutil.rmtree(bundle)
        for source, target, expected in self._skill_outputs():
            content = source.read_bytes()
            if expected is not None and hashlib.sha256(content).hexdigest() != expected:
                raise ValueError(f"Skill resource changed without rebuilding its catalog: {source}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)

    def get_outputs(self, include_bytecode: bool = True) -> list[str]:
        outputs = super().get_outputs(include_bytecode)
        if not self.editable_mode:
            outputs.extend(str(target) for _, target, _ in self._skill_outputs())
        return outputs


setup(cmdclass={"build_py": BuildWithSkills})
