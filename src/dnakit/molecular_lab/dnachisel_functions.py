"""DnaChisel sequence functions exposed through ``dnakit.molecular_lab``.

Each function imports DnaChisel only when it runs. A missing package raises
``BackendUnavailableError`` with the ``dnakit[bio-dnachisel]`` install hint.
Positions are 0-based; segments and intervals are half-open ``[start, end)``.
Every function returns a JSON-safe dictionary with ``status``, ``operation``,
the result fields and ``provenance``.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-dnachisel"
_MAX_DNA = 200_000
_MAX_LENGTH = 20_000
_MAX_CONSTRAINTS = 20
_MAX_POSITIONS = 100_000
_CONSTRAINTS = frozenset({"AvoidHairpins", "AvoidPattern", "AvoidStopCodons", "EnforceGCContent"})


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _dna(value: object, name: str, *, maximum: int = _MAX_DNA) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = "".join(value.split()).upper()
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - set("ACGT"):
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    return text


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnachisel")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnachisel", "version": installed, **notes}


def _build_constraints(value: object, dnachisel: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("constraints must be a list of specification mappings.")
    if len(value) > _MAX_CONSTRAINTS:
        raise ConfigurationError(f"constraints accepts at most {_MAX_CONSTRAINTS} entries.")
    built: list[Any] = []
    for item in value:
        if not isinstance(item, Mapping) or not isinstance(item.get("name"), str):
            raise ConfigurationError("Each constraint needs a string 'name' key.")
        name = item["name"]
        if name not in _CONSTRAINTS:
            allowed = ", ".join(sorted(_CONSTRAINTS))
            raise ConfigurationError(
                f"Unsupported DnaChisel constraint {name!r}; allowed: {allowed}."
            )
        kwargs: dict[str, Any] = {}
        for key, argument in item.items():
            if key == "name":
                continue
            if key == "location" and isinstance(argument, list):
                argument = tuple(argument)
            kwargs[str(key)] = argument
        try:
            built.append(getattr(dnachisel, name)(**kwargs))
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(f"Invalid DnaChisel constraint {name!r}: {exc}") from exc
    return built


def generate_constrained_dna_sequence(
    sequence_length: int,
    *,
    constraints: Sequence[Mapping[str, Any]] = (),
    seed: int | None = None,
    max_random_iters: int = 5000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a random DNA sequence that satisfies DnaChisel constraints.

    Each constraint is a mapping with a ``name`` from the allow-list and keyword
    arguments, for example ``{"name": "AvoidPattern", "pattern": "GGTCTC"}``.
    The generator is seeded when ``seed`` is given and is reproducible.
    """

    dnachisel = require("dnachisel", _EXTRA)
    length = _integer(sequence_length, "sequence_length", minimum=1, maximum=_MAX_LENGTH)
    specs = _build_constraints(constraints, dnachisel)
    if seed is not None:
        seed = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    iterations = _integer(max_random_iters, "max_random_iters", minimum=1, maximum=100_000)
    try:
        sequence = str(
            dnachisel.random_compatible_dna_sequence(
                sequence_length=length,
                constraints=specs,
                seed=seed,
                max_random_iters=iterations,
                logger=None,
            )
        )
        problem = dnachisel.DnaOptimizationProblem(
            sequence=sequence, constraints=specs, logger=None
        )
        passes = bool(problem.all_constraints_pass())
    except (dnachisel.NoSolutionError, OSError, TypeError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"DnaChisel could not generate a sequence: {exc}") from exc
    return {
        "status": "completed",
        "operation": "generate_constrained_dna_sequence",
        "sequence": sequence,
        "sequence_length": len(sequence),
        "seed": seed,
        "constraint_names": [str(item["name"]) for item in constraints],
        "all_constraints_pass": passes,
        "provenance": _provenance(
            generator="dnachisel.random_compatible_dna_sequence",
            alphabet="ACGT, uniform before constraints",
        ),
    }


def find_sequence_difference_segments(
    sequence_a: str,
    sequence_b: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the contiguous segments where two equal-length DNA sequences differ.

    Adjacent differing positions are merged into one ``[start, end)`` segment.
    The two sequences must have the same length; no alignment is performed.
    """

    module = require("dnachisel.biotools.sequences_differences", _EXTRA)
    first = _dna(sequence_a, "sequence_a")
    second = _dna(sequence_b, "sequence_b")
    if len(first) != len(second):
        raise ConfigurationError("sequence_a and sequence_b must have the same length.")
    raw = module.sequences_differences_segments(first, second)
    segments = sorted([int(start), int(end)] for start, end in raw)
    return {
        "status": "completed",
        "operation": "find_sequence_difference_segments",
        "sequence_length": len(first),
        "segment_count": len(segments),
        "differing_positions": sum(end - start for start, end in segments),
        "segments": segments,
        "provenance": _provenance(method="DnaChisel positional comparison, no alignment"),
    }


def list_common_restriction_enzymes(
    *,
    site_length: int = 6,
    min_suppliers: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """List restriction enzymes with a given site length from DnaChisel's enzyme table.

    ``min_suppliers`` keeps enzymes that at least that many suppliers sell.
    Names are returned in sorted order.
    """

    dnachisel = require("dnachisel", _EXTRA)
    length = _integer(site_length, "site_length", minimum=4, maximum=12)
    suppliers = _integer(min_suppliers, "min_suppliers", minimum=0, maximum=20)
    try:
        names = sorted(
            str(name)
            for name in dnachisel.list_common_enzymes(
                site_length=(length,), min_suppliers=suppliers
            )
        )
    except (TypeError, ValueError, RuntimeError, KeyError) as exc:
        raise ConfigurationError(f"DnaChisel could not list enzymes: {exc}") from exc
    return {
        "status": "completed",
        "operation": "list_common_restriction_enzymes",
        "site_length": length,
        "min_suppliers": suppliers,
        "enzyme_count": len(names),
        "enzymes": names,
        "provenance": _provenance(table="DnaChisel bundled enzyme table"),
    }


def group_nearby_positions(
    positions: Sequence[int],
    *,
    max_gap: int | None = None,
    max_group_spread: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group sorted nucleotide positions that lie close together.

    Neighbouring positions join one group while their difference is below
    ``max_gap``; a group also stays within ``max_group_spread`` of its first
    position. ``None`` disables a rule. Positions must be unique non-negative
    integers; groups and their members are returned sorted.
    """

    module = require("dnachisel.biotools.indices_operations", _EXTRA)
    if isinstance(positions, (str, bytes)) or not isinstance(positions, Sequence):
        raise ConfigurationError("positions must be a list of integers.")
    if len(positions) > _MAX_POSITIONS:
        raise ConfigurationError(f"positions accepts at most {_MAX_POSITIONS} values.")
    values = [_integer(item, "positions", minimum=0, maximum=2**53) for item in positions]
    if len(set(values)) != len(values):
        raise ConfigurationError("positions must be unique.")
    gap = None if max_gap is None else _integer(max_gap, "max_gap", minimum=0, maximum=10**12)
    spread = (
        None
        if max_group_spread is None
        else _integer(max_group_spread, "max_group_spread", minimum=0, maximum=10**12)
    )
    groups = module.group_nearby_indices(sorted(values), max_gap=gap, max_group_spread=spread)
    rows = sorted((sorted(int(item) for item in group) for group in groups), key=lambda g: g[0])
    return {
        "status": "completed",
        "operation": "group_nearby_positions",
        "position_count": len(values),
        "max_gap": gap,
        "max_group_spread": spread,
        "group_count": len(rows),
        "groups": rows,
        "provenance": _provenance(method="DnaChisel group_nearby_indices"),
    }


__all__ = [
    "find_sequence_difference_segments",
    "generate_constrained_dna_sequence",
    "group_nearby_positions",
    "list_common_restriction_enzymes",
]

_register_public_verbose(globals(), tuple(__all__))
