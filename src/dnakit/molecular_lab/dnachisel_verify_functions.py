"""Optimise a record whose constraints and objectives are annotated in the file itself (DnaChisel CLI).

``dnachisel.cli.main`` reads a GenBank or SnapGene record, takes its constraints and objectives from the
record's annotations (``DnaOptimizationProblem.from_record``), resolves the constraints, optimises the
objectives and writes either an annotated GenBank file (target ending in ``.gb``) or a DnaChisel report
folder or zip file. ``optimize_dnachisel_source_record`` performs the same steps with DnaChisel's API.
The existing ``optimize_dnachisel_record`` needs its constraints passed explicitly and is not used here.
"""

from __future__ import annotations

import importlib.metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnachisel"


def _source(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("source_path must be a path.")
    path = Path(value).expanduser()
    if not path.is_file():
        raise ConfigurationError(f"source_path must be an existing record file: {path}")
    return path


def _target(value: object) -> Path:
    if not isinstance(value, (str, Path)):
        raise ConfigurationError("target_path must be a path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:  # pragma: no cover
        return "unknown"


def optimize_dnachisel_source_record(
    source_path: str | Path,
    target_path: str | Path,
    *,
    circular: bool = False,
    with_sequence_edits: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Optimise a record using the constraints and objectives annotated in it, and write the result.

    A target ending in ``.gb`` receives an annotated GenBank record (the constraints are resolved first,
    then the objectives optimised). Any other target is a new folder or zip file with DnaChisel's
    optimisation report. ``circular=True`` treats the sequence as circular. ``with_sequence_edits=True``
    annotates the edited positions in the GenBank output. The result gives the constraint and objective
    summaries before and after optimisation.
    """
    dnachisel = require("dnachisel", _EXTRA)
    source = _source(source_path)
    target = _target(target_path)
    problem_class = (
        dnachisel.CircularDnaOptimizationProblem if circular else dnachisel.DnaOptimizationProblem
    )
    problem = problem_class.from_record(str(source), logger=None)
    problem.max_random_iters = 10000
    before = {
        "constraints": problem.constraints_text_summary(),
        "objectives": problem.objectives_text_summary(),
    }
    if str(target).lower().endswith(".gb"):
        problem.resolve_constraints()
        problem.optimize()
        problem.to_record(filepath=str(target), with_sequence_edits=bool(with_sequence_edits))
        output_kind = "genbank"
    else:
        problem.optimize_with_report(str(target), project_name=source.name, file_path=str(source))
        output_kind = "report"
    after = {
        "constraints": problem.constraints_text_summary(),
        "objectives": problem.objectives_text_summary(),
    }
    return {
        "status": "completed",
        "operation": "optimize_dnachisel_source_record",
        "source": source.name,
        "output_path": str(target),
        "output_kind": output_kind,
        "circular": bool(circular),
        "before": before,
        "after": after,
        "provenance": {
            "package": f"dnachisel {_version('dnachisel')}",
            "convention": "DnaOptimizationProblem.from_record, resolve_constraints, optimize, to_record "
            "or optimize_with_report (as the DnaChisel command line)",
        },
    }


__all__ = ["optimize_dnachisel_source_record"]

_register_public_verbose(globals())
