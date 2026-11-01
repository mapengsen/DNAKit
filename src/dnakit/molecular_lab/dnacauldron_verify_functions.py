"""Simulate a DNA Cauldron assembly plan from a spreadsheet and write its report.

``dnacauldron.cli.main`` reads an assembly plan (CSV or spreadsheet), imports the part sequences from
a folder and the connector sequences from one sub-folder per collection, simulates every assembly in
the plan and writes the report. ``simulate_dnacauldron_plan`` does the same steps through DNA Cauldron's
own API, returns the simulation statistics and writes the report to a new folder or zip file. The
existing ``simulate_dnacauldron_assembly`` takes a parts mapping instead of a plan and is not used here.
"""

from __future__ import annotations

import importlib.metadata
import os
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnacauldron"
_MAX_PLAN_FILE_BYTES = 50_000_000


def _existing_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    if path.stat().st_size > _MAX_PLAN_FILE_BYTES:
        raise ConfigurationError(f"{name} is larger than {_MAX_PLAN_FILE_BYTES} bytes.")
    return path


def _existing_dir(value: object, name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser()
    if not path.is_dir():
        raise ConfigurationError(f"{name} must be an existing folder of sequence files: {path}")
    return path


def _new_report(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("target_path must be a path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Report already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"


def simulate_dnacauldron_plan(
    assembly_plan_path: str | Path,
    sequences_path: str | Path,
    target_path: str | Path,
    *,
    connectors_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate every assembly in a DNA Cauldron plan and write the simulation report.

    ``assembly_plan_path`` is a CSV or spreadsheet plan (one row per assembly, a ``class: ...`` cell
    naming the assembly class). ``sequences_path`` is a folder of part sequences. ``connectors_path``
    (optional) is a folder whose sub-folders each hold one connector collection. ``target_path`` is a
    new folder or ``.zip`` file for the report. The result lists the outcome counts (valid, cancelled,
    errored) and the report path.
    """
    dnacauldron = require("dnacauldron", _EXTRA)
    plan_file = _existing_file(assembly_plan_path, "assembly_plan_path")
    sequences = _existing_dir(sequences_path, "sequences_path")
    connectors = None if connectors_path is None else _existing_dir(connectors_path, "connectors_path")
    target = _new_report(target_path)
    repository = dnacauldron.SequenceRepository()
    repository.import_records(folder=str(sequences))
    if connectors is not None:
        for subfolder_name in sorted(os.listdir(connectors)):
            sub_path = connectors / subfolder_name
            if sub_path.is_dir():
                repository.import_records(folder=str(sub_path), collection=subfolder_name)
    assembly_plan = dnacauldron.AssemblyPlan.from_spreadsheet(str(plan_file), logger=None)
    simulation = assembly_plan.simulate(sequence_repository=repository)
    stats = {str(key): int(value) for key, value in simulation.compute_stats().items()}
    simulation.write_report(str(target), logger=None)
    return {
        "status": "completed",
        "operation": "simulate_dnacauldron_plan",
        "assembly_plan": plan_file.name,
        "assemblies": len(assembly_plan.assemblies),
        "outcomes": dict(sorted(stats.items())),
        "report_path": str(target),
        "provenance": {
            "package": f"dnacauldron {_version('dnacauldron')}",
            "convention": "assembly plan simulated with AssemblyPlan.simulate; report from write_report",
        },
    }


__all__ = ["simulate_dnacauldron_plan"]

_register_public_verbose(globals())
