"""Versioned native DNA/DNA nearest-neighbor parameters."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, cast

from dnakit.core._json import to_json_compatible


@dataclass(frozen=True)
class NearestNeighborParameter:
    """One enthalpy/entropy contribution in explicit units."""

    delta_h_kcal_per_mol: float
    delta_s_cal_per_k_mol: float


@dataclass(frozen=True)
class NearestNeighborParameterSet:
    """Immutable identity and values for a published parameter table."""

    name: str
    version: str
    citation: str
    doi: str
    reference_sodium_molar: float
    stacking: MappingProxyType[str, NearestNeighborParameter]
    terminal_at: NearestNeighborParameter
    terminal_gc: NearestNeighborParameter
    symmetry: NearestNeighborParameter

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


_STACKING = {
    "AA": NearestNeighborParameter(-7.9, -22.2),
    "TT": NearestNeighborParameter(-7.9, -22.2),
    "AT": NearestNeighborParameter(-7.2, -20.4),
    "TA": NearestNeighborParameter(-7.2, -21.3),
    "CA": NearestNeighborParameter(-8.5, -22.7),
    "TG": NearestNeighborParameter(-8.5, -22.7),
    "GT": NearestNeighborParameter(-8.4, -22.4),
    "AC": NearestNeighborParameter(-8.4, -22.4),
    "CT": NearestNeighborParameter(-7.8, -21.0),
    "AG": NearestNeighborParameter(-7.8, -21.0),
    "GA": NearestNeighborParameter(-8.2, -22.2),
    "TC": NearestNeighborParameter(-8.2, -22.2),
    "CG": NearestNeighborParameter(-10.6, -27.2),
    "GC": NearestNeighborParameter(-9.8, -24.4),
    "GG": NearestNeighborParameter(-8.0, -19.9),
    "CC": NearestNeighborParameter(-8.0, -19.9),
}

SANTALUCIA_1998 = NearestNeighborParameterSet(
    name="SantaLucia unified DNA/DNA nearest-neighbor",
    version="santalucia1998-v1",
    citation=(
        "SantaLucia J Jr. A unified view of polymer, dumbbell, and "
        "oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 1998."
    ),
    doi="10.1073/pnas.95.4.1460",
    reference_sodium_molar=1.0,
    stacking=MappingProxyType(_STACKING),
    terminal_at=NearestNeighborParameter(2.3, 4.1),
    terminal_gc=NearestNeighborParameter(0.1, -2.8),
    symmetry=NearestNeighborParameter(0.0, -1.4),
)

PARAMETER_SETS = MappingProxyType({SANTALUCIA_1998.version: SANTALUCIA_1998})

__all__ = [
    "PARAMETER_SETS",
    "SANTALUCIA_1998",
    "NearestNeighborParameter",
    "NearestNeighborParameterSet",
]
