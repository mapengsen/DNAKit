"""Validated conditions and safety limits for thermodynamic calculations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, TypeAlias

from dnakit.exceptions import ConfigurationError

SaltModel: TypeAlias = Literal["santalucia1998-monovalent-entropy"]
TmMethod: TypeAlias = Literal["wallace", "nearest_neighbor"]


def _finite_number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"{name} must be a finite number.",
            code="INVALID_THERMODYNAMIC_CONDITION",
            context={"field": name, "value": value},
        )
    return float(value)


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code="INVALID_THERMODYNAMIC_LIMIT",
            context={"field": name, "value": value},
        )
    return value


@dataclass(frozen=True, init=False)
class ThermodynamicConditions:
    """Physical conditions recorded with every temperature-dependent result.

    Concentrations are expressed in mol/L. ``strand_concentration_molar`` is
    the total oligonucleotide strand concentration used by the native Tm
    equation. Native calculations treat Na+ and K+ as an explicit total
    monovalent concentration. Magnesium, dNTP and cosolvent fields are retained
    for backend calculations or explicit empirical corrections.
    """

    temperature_celsius: float
    sodium_molar: float
    potassium_molar: float
    magnesium_molar: float
    dntp_molar: float
    strand_concentration_molar: float
    dmso_percent: float
    dmso_factor_celsius_per_percent: float
    formamide_molar: float
    salt_model: SaltModel

    def __init__(
        self,
        *,
        temperature_celsius: float = 37.0,
        sodium_molar: float = 0.05,
        potassium_molar: float = 0.0,
        magnesium_molar: float = 0.0,
        dntp_molar: float = 0.0,
        strand_concentration_molar: float = 250e-9,
        dmso_percent: float = 0.0,
        dmso_factor_celsius_per_percent: float = 0.6,
        formamide_molar: float = 0.0,
        salt_model: SaltModel = "santalucia1998-monovalent-entropy",
    ) -> None:
        temperature = _finite_number(temperature_celsius, "temperature_celsius")
        sodium = _finite_number(sodium_molar, "sodium_molar")
        potassium = _finite_number(potassium_molar, "potassium_molar")
        magnesium = _finite_number(magnesium_molar, "magnesium_molar")
        dntp = _finite_number(dntp_molar, "dntp_molar")
        concentration = _finite_number(strand_concentration_molar, "strand_concentration_molar")
        dmso = _finite_number(dmso_percent, "dmso_percent")
        dmso_factor = _finite_number(
            dmso_factor_celsius_per_percent, "dmso_factor_celsius_per_percent"
        )
        formamide = _finite_number(formamide_molar, "formamide_molar")
        if not 0.0 <= temperature <= 100.0:
            raise ConfigurationError(
                "temperature_celsius must be in the modeled interval [0, 100] C.",
                code="THERMODYNAMIC_TEMPERATURE_OUT_OF_DOMAIN",
            )
        if not 0.0 <= sodium <= 2.0 or not 0.0 <= potassium <= 2.0:
            raise ConfigurationError(
                "sodium_molar and potassium_molar must each be in [0, 2].",
                code="THERMODYNAMIC_MONOVALENT_OUT_OF_DOMAIN",
            )
        if not 0.0 < sodium + potassium <= 2.0:
            raise ConfigurationError(
                "The total monovalent concentration (Na+ + K+) must be in (0, 2].",
                code="THERMODYNAMIC_MONOVALENT_OUT_OF_DOMAIN",
            )
        if not 0.0 <= magnesium <= 1.0 or not 0.0 <= dntp <= 1.0:
            raise ConfigurationError(
                "magnesium_molar and dntp_molar must be in the interval [0, 1].",
                code="THERMODYNAMIC_DIVALENT_OUT_OF_DOMAIN",
            )
        if not 0.0 < concentration <= 1.0:
            raise ConfigurationError(
                "strand_concentration_molar must be in the interval (0, 1].",
                code="THERMODYNAMIC_CONCENTRATION_OUT_OF_DOMAIN",
            )
        if not 0.0 <= dmso <= 100.0 or not 0.0 <= dmso_factor <= 2.0:
            raise ConfigurationError(
                "dmso_percent must be in [0, 100] and its factor in [0, 2].",
                code="THERMODYNAMIC_COSOLVENT_OUT_OF_DOMAIN",
            )
        if not 0.0 <= formamide <= 30.0:
            raise ConfigurationError(
                "formamide_molar must be in [0, 30].",
                code="THERMODYNAMIC_COSOLVENT_OUT_OF_DOMAIN",
            )
        if salt_model != "santalucia1998-monovalent-entropy":
            raise ConfigurationError(
                "The native implementation supports only the SantaLucia 1998 monovalent model.",
                code="UNSUPPORTED_SALT_MODEL",
                context={"salt_model": salt_model},
            )
        object.__setattr__(self, "temperature_celsius", temperature)
        object.__setattr__(self, "sodium_molar", sodium)
        object.__setattr__(self, "potassium_molar", potassium)
        object.__setattr__(self, "magnesium_molar", magnesium)
        object.__setattr__(self, "dntp_molar", dntp)
        object.__setattr__(self, "strand_concentration_molar", concentration)
        object.__setattr__(self, "dmso_percent", dmso)
        object.__setattr__(self, "dmso_factor_celsius_per_percent", dmso_factor)
        object.__setattr__(self, "formamide_molar", formamide)
        object.__setattr__(self, "salt_model", salt_model)

    @property
    def monovalent_molar(self) -> float:
        """Total Na+ plus K+ concentration used by supported monovalent models."""

        return self.sodium_molar + self.potassium_molar


@dataclass(frozen=True, init=False)
class NearestNeighborConfig:
    """Native model selection and hard sequence-size safety limit."""

    parameter_set: str
    max_sequence_length: int

    def __init__(
        self,
        *,
        parameter_set: str = "santalucia1998-v1",
        max_sequence_length: int = 60,
    ) -> None:
        if parameter_set != "santalucia1998-v1":
            raise ConfigurationError(
                "Unknown native nearest-neighbor parameter set.",
                code="UNKNOWN_THERMODYNAMIC_PARAMETER_SET",
                context={"parameter_set": parameter_set},
            )
        limit = _positive_integer(max_sequence_length, "max_sequence_length")
        if not 2 <= limit <= 60:
            raise ConfigurationError(
                "The native nearest-neighbor applicability limit must be in [2, 60] nt.",
                code="THERMODYNAMIC_APPLICABILITY_LIMIT",
            )
        object.__setattr__(self, "parameter_set", parameter_set)
        object.__setattr__(self, "max_sequence_length", limit)


__all__ = [
    "NearestNeighborConfig",
    "SaltModel",
    "ThermodynamicConditions",
    "TmMethod",
]
