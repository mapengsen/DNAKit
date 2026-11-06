"""Codon/mutation-selection simulation and state-frequency estimation (P41)."""

from __future__ import annotations

from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as h
from dnakit.bioinformatics._common import require


def simulate_pyvolve_sequences(
    newick: str,
    *,
    model: str = "gy",
    parameters: dict[str, Any] | None = None,
    size: int = 100,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate nucleotide output with gy/mg/codon/ecm/ecmrest/ecmunrest/mutsel models.

    size is the native number of states (codons for codon models). Parameters are
    the selected model's scientific parameters, never executable code or paths.
    Native model defaults are retained. No automatic simulation files are created.
    """
    kind = h.choice(model, ("gy", "mg", "codon", "ecm", "ecmrest", "ecmunrest", "mutsel"), "model")
    count = h.integer(size, "size", minimum=1, maximum=100000)
    if seed is not None:
        h.integer(seed, "seed", minimum=0)
    module = require("pyvolve", "bio-pyvolve")
    tree = module.read_tree(tree=h.text(newick, "newick"))
    native = module.Model(kind, parameters=dict(parameters or {}))
    partition = module.Partition(models=native, size=count)
    evolver = module.Evolver(tree=tree, partitions=partition)
    evolver(seqfile=None, ratefile=None, infofile=None, seed=seed)
    return {
        "status": "completed",
        "backend": "pyvolve",
        "model": kind,
        "size": count,
        "sequences": evolver.get_sequences(),
    }


def pyvolve_state_frequencies(
    path: str, *, alphabet: str = "nucleotide", file_format: str = "fasta", verbose: bool = True
) -> dict[str, Any]:
    """Estimate A/C/G/T or 61 sense-codon frequencies using native ReadFrequencies."""
    kind = h.choice(alphabet, ("nucleotide", "codon"), "alphabet")
    module = require("pyvolve", "bio-pyvolve")
    freq = module.ReadFrequencies(
        kind,
        file=str(h.existing(path, "path")),
        format=h.choice(
            file_format, ("fasta", "phylip", "phylip-relaxed", "nexus", "stockholm"), "file_format"
        ),
    )
    values = freq.compute_frequencies()
    states = module.Genetics().nucleotides if kind == "nucleotide" else module.Genetics().codons
    return {
        "status": "completed",
        "backend": "pyvolve",
        "alphabet": kind,
        "frequencies": dict(zip(states, values.tolist(), strict=False)),
    }


__all__ = ["pyvolve_state_frequencies", "simulate_pyvolve_sequences"]
_register_public_verbose(globals(), tuple(__all__))
