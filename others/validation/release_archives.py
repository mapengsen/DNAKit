"""Check release archives against source files and the checksummed Skills catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN = {
    ".cache",
    ".hypothesis",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "dashing_similarity",
    "site",
    ".git",
    ".venv",
    ".dnakit",
}
NATIVE_SUFFIXES = {".a", ".dll", ".dylib", ".exe", ".lib", ".o", ".obj", ".pyd", ".so"}


def audit(dist: Path, *, root: Path = ROOT) -> None:
    """Reject missing code/resources, changed Skills, and unintended archive files."""
    wheels = sorted(dist.glob("dnakit-*.whl"))
    sdists = sorted(dist.glob("dnakit-*.tar.gz"))
    assert len(wheels) == len(sdists) == 1, "Expected exactly one wheel and one sdist."
    catalog_path = root / "registry/catalog.json"
    catalog_bytes = catalog_path.read_bytes()
    catalog = json.loads(catalog_bytes)
    bundle = {"catalog.json": hashlib.sha256(catalog_bytes).hexdigest()}
    for skill in catalog["skills"]:
        for resource in skill["files"]:
            bundle[f"{skill['id']}/{resource['path']}"] = resource["sha256"]
    expected_python = {
        p.relative_to(root / "src").as_posix()
        for p in (root / "src/dnakit").rglob("*")
        if p.is_file() and (p.suffix == ".py" or p.name == "py.typed")
    }
    for archive in (*wheels, *sdists):
        assert archive.stat().st_size <= 100_000_000, f"Archive too large: {archive}"
        wheel = archive.suffix == ".whl"
        contents: dict[str, bytes] = {}
        if wheel:
            with zipfile.ZipFile(archive) as handle:
                entries = handle.infolist()
                names = [entry.filename for entry in entries]
                assert all(
                    not stat.S_ISLNK((entry.external_attr >> 16) & 0xFFFF) for entry in entries
                ), f"Symbolic link in {archive}"
                contents = {entry.filename: handle.read(entry) for entry in entries}
        else:
            with tarfile.open(archive, "r:gz") as handle:
                members = handle.getmembers()
                names = [entry.name for entry in members]
                assert all(entry.isfile() or entry.isdir() for entry in members), (
                    f"Link or special member in {archive}"
                )
                for member in members:
                    if member.isfile():
                        stream = handle.extractfile(member)
                        assert stream is not None
                        contents[member.name] = stream.read()
        assert 0 < len(names) <= 30_000, f"Unexpected archive member count: {archive}"
        assert len(names) == len(set(names)), f"Duplicate archive members: {archive}"
        for name in names:
            path = PurePosixPath(name)
            assert not path.is_absolute() and ".." not in path.parts, f"Unsafe path: {name}"
            assert not FORBIDDEN.intersection(p.lower() for p in path.parts), (
                f"Unintended directory: {name}"
            )
            assert not any(
                path.name.lower().endswith(suffix) or f"{suffix}." in path.name.lower()
                for suffix in NATIVE_SUFFIXES
            ), f"Native binary: {name}"
        if wheel:
            package = {name: data for name, data in contents.items() if name.startswith("dnakit/")}
            bundled = {
                name.removeprefix("dnakit/skills/_data/"): data
                for name, data in package.items()
                if name.startswith("dnakit/skills/_data/")
            }
            assert bundled.keys() == bundle.keys(), "Bundled Skills differ from the catalog."
            for name, data in bundled.items():
                assert hashlib.sha256(data).hexdigest() == bundle[name], (
                    f"Skill checksum mismatch: {name}"
                )
            native = {
                name: data
                for name, data in package.items()
                if not name.startswith("dnakit/skills/_data/")
            }
            for name, data in native.items():
                source = root / "src" / name
                assert source.is_file() and source.read_bytes() == data, (
                    f"Package resource differs from source: {name}"
                )
            assert any(name.endswith("/licenses/LICENSE") for name in contents), (
                "Wheel LICENSE missing."
            )
        else:
            prefix = archive.name.removesuffix(".tar.gz") + "/"
            assert all(name == prefix[:-1] or name.startswith(prefix) for name in names), (
                "Unexpected sdist root."
            )
            relative = {name.removeprefix(prefix): data for name, data in contents.items()}
            for name, data in relative.items():
                if name == "PKG-INFO" or name.startswith("src/dnakit.egg-info/"):
                    continue
                if name == "setup.cfg" and not (root / name).exists():
                    assert data == b"[egg_info]\ntag_build = \ntag_date = 0\n\n", (
                        "Unexpected generated setuptools configuration."
                    )
                    continue
                source = root / name
                assert source.is_file() and source.read_bytes() == data, (
                    f"Sdist file differs from source: {name}"
                )
            for name, digest in bundle.items():
                source_name = (
                    "registry/catalog.json" if name == "catalog.json" else f"skills/{name}"
                )
                assert source_name in relative, f"Sdist skill missing: {source_name}"
                assert hashlib.sha256(relative[source_name]).hexdigest() == digest, (
                    f"Sdist skill checksum mismatch: {source_name}"
                )
            assert {"LICENSE", "docs/planning/requirements_traceability.csv"} <= relative.keys()
            native = {
                name.removeprefix("src/"): data
                for name, data in relative.items()
                if name.startswith("src/dnakit/")
            }
        archived_python = {
            name
            for name in native
            if PurePosixPath(name).suffix == ".py" or PurePosixPath(name).name == "py.typed"
        }
        assert archived_python == expected_python, (
            f"Packaged Python files differ from src/dnakit: {archive}"
        )
        print(f"PASS {archive.name}: {len(contents)} files; {len(bundle) - 1} Skill resources")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dist", type=Path)
    audit(parser.parse_args().dist)
