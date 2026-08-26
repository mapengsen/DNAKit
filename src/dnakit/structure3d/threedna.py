"""Parsers for standard 3DNA base-pair steps and DSSR JSON summaries."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from pathlib import Path
from typing import cast

from dnakit.exceptions import ConfigurationError, InputFormatError

from ._shared import threedna_provenance
from .results import BasePairStepParameters, DSSRSummaryResult, ThreeDNAParameterResult

_MAX_PARAMETER_BYTES = 100_000_000
_MAX_BASE_PAIRS = 1_000_000


def _source(path: str | Path) -> tuple[Path, str]:
    if not isinstance(path, (str, Path)):
        raise ConfigurationError("path must be text or Path.", code="INVALID_3DNA_PATH")
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise InputFormatError(
            "3DNA/DSSR input file does not exist.",
            code="THREEDNA_FILE_NOT_FOUND",
            context={"path": str(source)},
        )
    size = source.stat().st_size
    if not 1 <= size <= _MAX_PARAMETER_BYTES:
        raise InputFormatError(
            f"3DNA/DSSR input size must be in [1, {_MAX_PARAMETER_BYTES}] bytes.",
            code="THREEDNA_FILE_SIZE_LIMIT",
        )
    try:
        return source, source.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise InputFormatError(
            "3DNA/DSSR text must be UTF-8 compatible.", code="INVALID_3DNA_ENCODING"
        ) from exc


def read_3dna_bp_step(path: str | Path) -> ThreeDNAParameterResult:
    """Parse the 12 standard rigid-body columns from ``bp_step.par``."""

    source, text = _source(path)
    lines = tuple(line.strip() for line in text.splitlines() if line.strip())
    if not lines:
        raise InputFormatError("3DNA parameter file is empty.", code="EMPTY_3DNA_PARAMETERS")
    try:
        declared_count = int(lines[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise InputFormatError(
            "The first 3DNA bp_step.par line must declare the base-pair count.",
            code="INVALID_3DNA_PARAMETERS",
        ) from exc
    if not 1 <= declared_count <= _MAX_BASE_PAIRS:
        raise InputFormatError(
            "Declared 3DNA base-pair count is outside safety limits.",
            code="INVALID_3DNA_PARAMETERS",
        )
    parameters: list[BasePairStepParameters] = []
    for line in lines[1:]:
        tokens = line.split()
        if len(tokens) != 13 or not any(separator in tokens[0] for separator in ("-", "/")):
            continue
        try:
            values = tuple(float(token) for token in tokens[1:])
        except ValueError:
            continue
        if any(not math.isfinite(value) for value in values):
            raise InputFormatError(
                "3DNA parameters must be finite.", code="INVALID_3DNA_PARAMETERS"
            )
        parameters.append(
            BasePairStepParameters(
                len(parameters) + 1,
                tokens[0],
                *values,
            )
        )
    if len(parameters) != declared_count:
        raise InputFormatError(
            "Parsed 3DNA parameter rows do not match the declared base-pair count.",
            code="THREEDNA_PARAMETER_COUNT_MISMATCH",
            context={"declared": declared_count, "parsed": len(parameters)},
        )
    step_rows = parameters[1:]
    mean_rise = (
        None
        if not step_rows
        else math.fsum(item.rise_angstrom for item in step_rows) / len(step_rows)
    )
    mean_twist = (
        None
        if not step_rows
        else math.fsum(item.twist_degree for item in step_rows) / len(step_rows)
    )
    base_pairs_per_turn = None if mean_twist is None or mean_twist == 0.0 else 360.0 / mean_twist
    pitch = (
        None
        if mean_rise is None or base_pairs_per_turn is None
        else mean_rise * base_pairs_per_turn
    )
    return ThreeDNAParameterResult(
        source_path=str(source),
        base_pair_count=declared_count,
        parameters=tuple(parameters),
        mean_rise_angstrom=mean_rise,
        mean_twist_degree=mean_twist,
        helical_pitch_angstrom=pitch,
        base_pairs_per_turn=base_pairs_per_turn,
        method="3dna-bp-step-12-parameter-parser-v1",
        provenance=threedna_provenance(),
    )


def _count(data: Mapping[str, object], count_key: str, array_key: str) -> int:
    raw_count = data.get(count_key)
    raw_array = data.get(array_key, ())
    if raw_count is None:
        if isinstance(raw_array, list):
            return len(raw_array)
        return 0
    if isinstance(raw_count, bool) or not isinstance(raw_count, int) or raw_count < 0:
        raise InputFormatError(
            f"DSSR {count_key} must be a non-negative integer.",
            code="INVALID_DSSR_JSON",
        )
    if isinstance(raw_array, list) and len(raw_array) != raw_count:
        raise InputFormatError(
            f"DSSR {count_key} disagrees with {array_key} length.",
            code="DSSR_COUNT_MISMATCH",
        )
    return raw_count


def read_dssr_json(path: str | Path) -> DSSRSummaryResult:
    """Parse audited capability counts from an x3dna-dssr ``--json`` result."""

    source, text = _source(path)
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise InputFormatError(
            "DSSR input is not valid JSON.",
            code="INVALID_DSSR_JSON",
            context={"line": exc.lineno, "column": exc.colno},
        ) from exc
    if not isinstance(raw, Mapping):
        raise InputFormatError("DSSR JSON root must be an object.", code="INVALID_DSSR_JSON")
    data = cast(Mapping[str, object], raw)
    version = data.get("dssr_ver", data.get("program", "unknown DSSR version"))
    if not isinstance(version, str) or not version.strip() or len(version) > 1_000:
        raise InputFormatError("DSSR program version is invalid.", code="INVALID_DSSR_JSON")
    torsions_key = "ntParams" if "ntParams" in data else "ntPars"
    raw_torsions = data.get(torsions_key)
    return DSSRSummaryResult(
        source_path=str(source),
        program_version=version,
        nucleotide_count=_count(data, "num_nts", "nts"),
        pair_count=_count(data, "num_pairs", "pairs"),
        helix_count=_count(data, "num_helices", "helices"),
        stem_count=_count(data, "num_stems", "stems"),
        hairpin_count=_count(data, "num_hairpins", "hairpins"),
        hydrogen_bond_count=_count(data, "num_hbonds", "hbonds"),
        backbone_torsion_record_count=len(raw_torsions) if isinstance(raw_torsions, list) else 0,
        method="x3dna-dssr-json-summary-parser-v1",
        provenance=threedna_provenance(),
    )


__all__ = ["read_3dna_bp_step", "read_dssr_json"]
