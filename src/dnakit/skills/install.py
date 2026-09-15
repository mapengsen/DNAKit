"""Explicit, local installation of one bundled skill into an Agent skill directory."""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError

from .catalog import _resource_bytes, get_skill

_CLIENT_DIRS = {
    "codex": ".agents/skills",
    "claude": ".claude/skills",
    "claude-code": ".claude/skills",
    "cursor": ".cursor/skills",
}
_MARKER = ".dnakit-install.json"


def install_skill(
    skill_id: str,
    *,
    agent: str | None = None,
    destination: str | Path | None = None,
    project_dir: str | Path | None = None,
    global_scope: bool = False,
    force: bool = False,
) -> dict[str, Any]:
    """Copy one complete skill, including resources and licenses, without running it.

    ``destination`` is a skill *container*, to which the skill identifier is appended.
    By default, Agent installations are project-scoped. ``force`` replaces only an
    existing directory previously installed by DNAKit under the same identifier.
    """
    entry = get_skill(skill_id)
    if destination is not None:
        if agent is not None or global_scope or project_dir is not None:
            raise ConfigurationError(
                "Use either destination or an Agent/scope selection.",
                code="SKILL_INSTALL_TARGET_INVALID",
            )
        container = Path(destination).expanduser()
    else:
        if agent not in _CLIENT_DIRS:
            raise ConfigurationError(
                "Select --agent codex, claude-code, or cursor; or provide --destination.",
                code="SKILL_INSTALL_TARGET_INVALID",
            )
        if global_scope and project_dir is not None:
            raise ConfigurationError(
                "A global install cannot specify a project.", code="SKILL_INSTALL_TARGET_INVALID"
            )
        base = Path.home() if global_scope else Path(project_dir or Path.cwd()).expanduser()
        container = base / _CLIENT_DIRS[agent]
    container = container.absolute()
    target = container / entry["id"]
    if target.is_symlink():
        raise ConfigurationError(
            "Refusing to replace a symlink skill directory.", code="SKILL_INSTALL_CONFLICT"
        )
    # Validate all bytes before creating or replacing any target files.
    resources = {item["path"]: _resource_bytes(entry, item["path"]) for item in entry["files"]}
    marker = {
        "id": entry["id"],
        "upstream_revision": entry["upstream_revision"],
        "files": entry["files"],
    }
    if target.exists():
        marker_path = target / _MARKER
        try:
            managed = json.loads(marker_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            managed = {}
        if managed.get("id") != entry["id"]:
            raise ConfigurationError(
                "The target is not a DNAKit-managed skill directory.", code="SKILL_INSTALL_CONFLICT"
            )
        identical = all(
            (target / name).is_file()
            and not (target / name).is_symlink()
            and (target / name).read_bytes() == data
            for name, data in resources.items()
        )
        if identical and managed == marker:
            return {"id": entry["id"], "path": str(target), "status": "unchanged"}
        if not force:
            raise ConfigurationError(
                "The installed skill differs from this release.",
                code="SKILL_INSTALL_CONFLICT",
                hint="Review local edits, then use --force to replace this DNAKit-managed skill.",
            )
    container.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{entry['id']}-", dir=container))
    backup: Path | None = None
    try:
        for name, data in resources.items():
            file = staging / name
            file.parent.mkdir(parents=True, exist_ok=True)
            file.write_bytes(data)
            mode = next(item.get("mode", 0o644) for item in entry["files"] if item["path"] == name)
            file.chmod(0o755 if mode & 0o111 else 0o644)
        (staging / _MARKER).write_text(
            json.dumps(marker, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        if target.exists():
            backup = Path(tempfile.mkdtemp(prefix=f".{entry['id']}-backup-", dir=container))
            backup.rmdir()
            target.rename(backup)
        try:
            staging.rename(target)
        except OSError:
            if backup is not None:
                backup.rename(target)
                backup = None
            raise
    finally:
        if staging.exists():
            shutil.rmtree(staging)
        if backup is not None and backup.exists():
            shutil.rmtree(backup)
    return {
        "id": entry["id"],
        "path": str(target),
        "status": "installed",
        "files": len(resources),
        "license": entry["license"],
        "note": "Skill files installed. Configure the declared runtime dependencies separately.",
    }
