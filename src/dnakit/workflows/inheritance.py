"""Per-locus inheritance simulations over explicit, user-provided registries."""

from __future__ import annotations

import copy
import importlib.util
import random
from pathlib import Path
from typing import Any, cast

from dnakit.exceptions import ConfigurationError


def simulate_inheritance(
    father: dict[str, Any],
    mother: dict[str, Any],
    trait_registry: dict[str, Any],
    disease_registry: dict[str, Any],
    *,
    offspring_count: int = 3,
    generation: int = 1,
    mutation_rate: float = 0.001,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate independent allele inheritance using explicit synthetic or supplied data.

    The attributed Genomebook engine implements per-locus segregation and mutation;
    it does not simulate chromosome-scale crossover or predict real clinical outcomes.
    A local random generator avoids changing randomness elsewhere in the process.
    """
    if (
        isinstance(offspring_count, bool)
        or not isinstance(offspring_count, int)
        or not 1 <= offspring_count <= 10000
    ):
        raise ConfigurationError("offspring_count must be an integer between 1 and 10000.")
    if isinstance(generation, bool) or not isinstance(generation, int) or generation < 1:
        raise ConfigurationError("generation must be a positive integer.")
    if (
        isinstance(mutation_rate, bool)
        or not isinstance(mutation_rate, (int, float))
        or not 0 <= mutation_rate <= 1
    ):
        raise ConfigurationError("mutation_rate must lie between 0 and 1.")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ConfigurationError("seed must be an integer or None.")
    if father.get("sex") != "Male" or mother.get("sex") != "Female":
        raise ConfigurationError(
            "The supplied simulation parent records must use Male/Female sex labels."
        )
    path = Path(__file__).parent / "support/recombinator.py"
    spec = importlib.util.spec_from_file_location("dnakit.workflows._inheritance_run", path)
    if spec is None or spec.loader is None:
        raise ImportError("DNAKit inheritance implementation is missing")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cast(Any, module).random = random.Random(seed)
    module.MUTATION_CONFIG["rate_per_locus"] = mutation_rate
    children = module.breed_pair(
        copy.deepcopy(father),
        copy.deepcopy(mother),
        generation,
        offspring_count,
        copy.deepcopy(trait_registry),
        copy.deepcopy(disease_registry),
    )
    return {
        "offspring": children,
        "seed": seed,
        "mutation_rate": mutation_rate,
        "scope": "per-locus simulation using supplied registries",
        "upstream": "ClawBio Genomebook, MIT",
    }
