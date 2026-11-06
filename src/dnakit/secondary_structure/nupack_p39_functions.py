"""NUPACK 4 DNA ensemble utilities and tube design (P39).

The separately licensed NUPACK Python package is required. API contracts follow
https://docs.nupack.org/utilities/, /analysis/ and /design/. No alternative energy
model is substituted when NUPACK is absent.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p40_helpers as _h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError


def _nupack() -> Any:
    try:
        return require("nupack", "nupack")
    except BackendUnavailableError as exc:
        raise BackendUnavailableError(
            "NUPACK 4 Python package is required.",
            hint="Obtain NUPACK 4 from https://www.nupack.org/downloads and follow its Python installation guide.",
        ) from exc


def _model(module: Any, options: Mapping[str, Any] | None) -> Any:
    options = dict(options or {})
    if set(options) - {"material", "ensemble", "celsius", "kelvin", "sodium", "magnesium"}:
        raise ConfigurationError("Unsupported NUPACK physical-model option.")
    options.setdefault("material", "dna")
    return module.Model(**options)


def _result(module: Any, **values: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        "backend": "nupack",
        "version": str(getattr(module, "__version__", "unknown")),
        **values,
    }


def _strands(strands: Sequence[str]) -> list[str]:
    if (
        isinstance(strands, str)
        or not strands
        or any(not isinstance(s, str) or not s for s in strands)
    ):
        raise ConfigurationError("strands must be a nonempty list of nonempty sequences.")
    return list(strands)


def nupack_partition_function(
    strands: Sequence[str], *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Complex ensemble partition function and free energy (kcal/mol); Decimal partition values are strings."""
    module = _nupack()
    partition, energy = module.pfunc(strands=_strands(strands), model=_model(module, model))
    return _result(module, partition_function=str(partition), free_energy=float(energy))


def nupack_structure_energy(
    strands: Sequence[str],
    structure: str,
    *,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Free energy in kcal/mol of a specified dot-parens-plus secondary structure."""
    module = _nupack()
    return _result(
        module,
        free_energy=float(
            module.structure_energy(
                strands=_strands(strands), structure=structure, model=_model(module, model)
            )
        ),
    )


def nupack_loop_energy(
    loop: str, structure: str, *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Calculate one loop's free energy in kcal/mol with the public Model.loop_energy API.

    ``loop`` is the loop sequence (plus signs delimit strand nicks); ``structure``
    is its NUPACK loop structure. Native sequence and structure notation are preserved.
    """
    module = _nupack()
    energy = _model(module, model).loop_energy(
        loop=_h.text(loop, "loop"), structure=_h.text(structure, "structure")
    )
    return _result(module, free_energy=float(energy))


def nupack_loop_stacking_energies(
    loop: str, structure: str, *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Calculate each loop stacking-state free energy in kcal/mol with Model.stack_energies.

    Returns the native state strings (s, b, 5, 3, n) mapped to their energies.
    For loops without stacking states, NUPACK returns the loop energy under state n.
    """
    module = _nupack()
    energies = _model(module, model).stack_energies(
        loop=_h.text(loop, "loop"), structure=_h.text(structure, "structure")
    )
    return _result(
        module, stacking_energies={str(state): float(value) for state, value in energies.items()}
    )


def nupack_sample_structures(
    strands: Sequence[str],
    *,
    count: int = 1,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw Boltzmann secondary structures using NUPACK's native sampling (num_sample)."""
    module = _nupack()
    structures = module.sample(
        strands=_strands(strands),
        num_sample=_h.integer(count, "count", minimum=1),
        model=_model(module, model),
    )
    return _result(module, structures=[str(s) for s in structures])


def nupack_pair_probabilities(
    strands: Sequence[str], *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Equilibrium base-pair probability matrix; diagonal entries denote unpaired probabilities."""
    module = _nupack()
    matrix = module.pairs(strands=_strands(strands), model=_model(module, model)).to_array()
    return _result(module, probabilities=matrix.tolist())


def _energies(values: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "structure": str(x.structure),
            "energy": float(x.energy),
            "stack_energy": float(x.stack_energy),
        }
        for x in values
    ]


def nupack_mfe_structures(
    strands: Sequence[str], *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """MFE proxy secondary structures with proxy and stacking-state free energies."""
    module = _nupack()
    return _result(
        module,
        structures=_energies(module.mfe(strands=_strands(strands), model=_model(module, model))),
    )


def nupack_suboptimal_structures(
    strands: Sequence[str],
    *,
    energy_gap: float = 1.0,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Proxy structures within energy_gap kcal/mol of the minimum stacking-state energy."""
    if energy_gap < 0:
        raise ConfigurationError("energy_gap must be nonnegative.")
    module = _nupack()
    return _result(
        module,
        structures=_energies(
            module.subopt(
                strands=_strands(strands), energy_gap=energy_gap, model=_model(module, model)
            )
        ),
    )


def nupack_ensemble_size(
    strands: Sequence[str], *, model: Mapping[str, Any] | None = None, verbose: bool = True
) -> dict[str, Any]:
    """Number of stacking states, or secondary structures for ensemble='nostacking'."""
    module = _nupack()
    return _result(
        module,
        ensemble_size=str(
            module.ensemble_size(strands=_strands(strands), model=_model(module, model))
        ),
    )


def nupack_design_complex(
    structure: str,
    *,
    constraints: Sequence[str] | None = None,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design one complex to a target structure, optionally with IUPAC sequence constraints."""
    module = _nupack()
    kwargs = {} if constraints is None else {"strands": _strands(constraints)}
    return _result(
        module,
        strands=[
            str(s) for s in module.des(structure=structure, model=_model(module, model), **kwargs)
        ],
    )


def nupack_sequence_distance(first: str, second: str, *, verbose: bool = True) -> dict[str, Any]:
    """NUPACK sequence distance, including its native run-length notation and nick handling."""
    module = _nupack()
    return _result(module, distance=int(module.seq_distance(first, second)))


def nupack_structure_distance(first: str, second: str, *, verbose: bool = True) -> dict[str, Any]:
    """NUPACK secondary-structure distance, accepting its native structure notations."""
    module = _nupack()
    return _result(module, distance=int(module.struc_distance(first, second)))


def _tube(
    module: Any, sequences: Mapping[str, str], concentrations: Mapping[str, float], max_size: int
) -> Any:
    if (
        set(sequences) != set(concentrations)
        or not sequences
        or any(v < 0 for v in concentrations.values())
    ):
        raise ConfigurationError(
            "Sequences and nonnegative molar concentrations must have identical nonempty names."
        )
    strands = {name: module.Strand(seq, name=name) for name, seq in sequences.items()}
    return module.Tube(
        strands={strands[n]: v for n, v in concentrations.items()},
        complexes=module.SetSpec(max_size=_h.integer(max_size, "max_size", minimum=1)),
        name="tube",
    )


def nupack_tube_analysis(
    sequences: Mapping[str, str],
    concentrations: Mapping[str, float],
    *,
    max_size: int = 1,
    save_analysis: str | None = None,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Analyze equilibrium complex concentrations of a DNA test tube; optionally save native analysis for reuse."""
    module = _nupack()
    tube = _tube(module, sequences, concentrations, max_size)
    output = None if save_analysis is None else _h.new_file(save_analysis, "save_analysis")
    result = module.tube_analysis(tubes=[tube], model=_model(module, model))
    if output is not None:
        result.save(str(output))
    values = [
        {"complex": str(c.name), "concentration": float(v)}
        for c, v in result[tube].complex_concentrations.items()
    ]
    return _result(module, concentrations=values, outputs=[] if output is None else [str(output)])


def nupack_complex_concentrations(
    analysis_file: str,
    sequences: Mapping[str, str],
    concentrations: Mapping[str, float],
    *,
    max_size: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recompute concentrations from a trusted native NUPACK AnalysisResult file without recomputing partition functions."""
    module = _nupack()
    data = module.AnalysisResult.load(str(_h.existing(analysis_file, "analysis_file")))
    tube = _tube(module, sequences, concentrations, max_size)
    result = module.complex_concentrations(tube=tube, data=data)
    return _result(
        module,
        concentrations=[
            {"complex": str(c.name), "concentration": float(v)}
            for c, v in result[tube].complex_concentrations.items()
        ],
    )


def _design(
    sequences: Mapping[str, str],
    complexes: list[dict[str, Any]],
    tubes: list[dict[str, Any]],
    model: Mapping[str, Any] | None,
    trials: int | None,
) -> dict[str, Any]:
    module = _nupack()
    strands = {
        n: module.TargetStrand([module.Domain(s, name="domain_" + n)], name=n)
        for n, s in sequences.items()
    }
    targets = {
        c["name"]: module.TargetComplex(
            [strands[n] for n in c["strands"]], c["structure"], name=c["name"]
        )
        for c in complexes
    }
    tubes_native = [
        module.TargetTube(
            {targets[n]: v for n, v in t["targets"].items()},
            off_targets=module.SetSpec(max_size=t.get("max_size", 1)),
            name=t["name"],
        )
        for t in tubes
    ]
    design = module.tube_design(tubes=tubes_native, model=_model(module, model))
    results = (
        [design.evaluate()]
        if trials is None
        else design.run(trials=_h.integer(trials, "trials", minimum=1))
    )
    return _result(
        module,
        results=[
            {
                "strands": {n: str(r.to_analysis(s)) for n, s in strands.items()},
                "ensemble_defect": float(r.defects.ensemble_defect),
                "report": str(r),
            }
            for r in results
        ],
    )


def nupack_design_tubes(
    sequences: Mapping[str, str],
    complexes: list[dict[str, Any]],
    tubes: list[dict[str, Any]],
    *,
    model: Mapping[str, Any] | None = None,
    trials: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Design multiple target tubes from IUPAC strand constraints, target structures and molar concentrations.

    complexes entries: name, strands (names), structure. tubes entries: name,
    targets (complex-name to molarity), max_size (off-target strand count, default 1).
    """
    return _design(sequences, complexes, tubes, model, trials)


def nupack_evaluate_design(
    sequences: Mapping[str, str],
    complexes: list[dict[str, Any]],
    tubes: list[dict[str, Any]],
    *,
    model: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Evaluate fully specified designed DNA strands under the requested physical model and tube targets."""
    return _design(sequences, complexes, tubes, model, None)


__all__ = [
    "nupack_complex_concentrations",
    "nupack_design_complex",
    "nupack_design_tubes",
    "nupack_ensemble_size",
    "nupack_evaluate_design",
    "nupack_loop_energy",
    "nupack_loop_stacking_energies",
    "nupack_mfe_structures",
    "nupack_pair_probabilities",
    "nupack_partition_function",
    "nupack_sample_structures",
    "nupack_sequence_distance",
    "nupack_structure_distance",
    "nupack_structure_energy",
    "nupack_suboptimal_structures",
    "nupack_tube_analysis",
]

_register_public_verbose(globals(), tuple(__all__))
