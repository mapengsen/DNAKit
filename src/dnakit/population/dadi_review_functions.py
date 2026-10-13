"""dadi review wrappers (group G18e): mutation-model rates, ms command lines, SLiM and VCF inputs,
neutral two-population spectra, ancestral misidentification tables, Godambe uncertainty,
nested-model statistics and residual fits.

Every input is local or passed in memory. Model functions are chosen by name from the fixed
family lists in ``dadi_functions``, so no arbitrary code is run. ``ms`` command lines are built
as text only; no simulator is started. Output files are never overwritten. Calls that draw random
numbers save and restore the global Python and NumPy random state. dadi, NumPy and SciPy are
imported only when a function runs.
"""

from __future__ import annotations

import random
from collections.abc import Iterator, Mapping, Sequence
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _finite, _integer, _number, _provenance, _text
from dnakit.population.dadi_functions import (
    _DEMOGRAPHICS_1D,
    _DEMOGRAPHICS_2D,
    _DEMOGRAPHICS_3D,
    _bounds,
    _choice,
    _flag,
    _model_function,
    _nested,
    _new_file,
    _numbers,
    _sizes,
)

_MS_CORE_MODELS = {"split_mig": 4, "IM": 6, "IM_pre": 8}
_MAX_MS_ITERATIONS = 100_000
_MAX_MS_SAMPLE = 1000
_MAX_BOOTSTRAPS = 200
_MAX_SLIM_FILES = 10
_TRI_STATES = 64
_GRID_MIN = 20
_GRID_MAX = 400
_NESTED_TESTS = ("lrt", "wald", "score")
_RATE_CONVENTION = (
    "rate_matrix[i][j] is the rate from state i to state j; the diagonal is ignored; states are "
    "the bases (4), dinucleotides (16) or trinucleotides (64) in dadi's CGTA order"
)


def _dadi_and_numpy() -> tuple[Any, Any]:
    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    return dadi, np


def _none_to_nan(value: Any) -> Any:
    if isinstance(value, list):
        return [_none_to_nan(item) for item in value]
    return float("nan") if value is None else value


def _masked(np: Any, value: object, name: str) -> Any:
    """A nested list of counts as a masked array; None marks an entry dadi masks."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(
            f"{name} must be a nested list of counts; masked entries are None."
        )
    try:
        raw = np.array(_none_to_nan(list(value)), dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular nested list of numbers.") from exc
    present = raw[~np.isnan(raw)]
    if np.isinf(present).any() or (present < 0).any():
        raise ConfigurationError(f"{name} must contain only finite, non-negative values or None.")
    return np.ma.masked_invalid(raw)


def _spectrum(dadi: Any, np: Any, value: object, name: str) -> Any:
    array = _masked(np, value, name)
    if array.ndim not in (1, 2, 3) or min(array.shape) < 3:
        raise ConfigurationError(f"{name} must have 1 to 3 axes, each of length at least 3.")
    return dadi.Spectrum(array)


def _bootstrap_spectra(dadi: Any, np: Any, value: object, shape: tuple[int, ...]) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("bootstraps must be a non-empty list of spectra.")
    if len(value) > _MAX_BOOTSTRAPS:
        raise ConfigurationError(f"bootstraps must list at most {_MAX_BOOTSTRAPS} spectra.")
    spectra = []
    for index, item in enumerate(value):
        spectrum = _spectrum(dadi, np, item, f"bootstraps[{index}]")
        if tuple(spectrum.shape) != shape:
            raise ConfigurationError(f"bootstraps[{index}] must have the shape of data {shape}.")
        spectra.append(spectrum)
    return spectra


def _families(dadi: Any, ndim: int) -> tuple[Any, tuple[str, ...]]:
    if ndim == 1:
        return dadi.Demographics1D, _DEMOGRAPHICS_1D
    if ndim == 2:
        return dadi.Demographics2D, _DEMOGRAPHICS_2D
    return dadi.Demographics3D, _DEMOGRAPHICS_3D


def _extrapolated_model(dadi: Any, observed: Any, model: str) -> Any:
    family, options = _families(dadi, observed.ndim)
    function = _model_function(family, model, options)
    return dadi.Numerics.make_extrap_func(function)


def _grids(value: object) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("grid_points must list one to five grid sizes.")
    grids = [
        _integer(item, "grid_points entry", minimum=_GRID_MIN, maximum=_GRID_MAX) for item in value
    ]
    if not grids or len(grids) > 5:
        raise ConfigurationError("grid_points must list one to five grid sizes.")
    return grids


def _start(value: object, *, positive: bool = True) -> list[float]:
    """Starting parameters; ``positive=False`` allows zeros, as nested models need."""

    start = _numbers(value, "p0")
    if any(item < 0 or (positive and item == 0) for item in start):
        raise ConfigurationError(
            "p0 must contain only positive values." if positive else "p0 must not be negative."
        )
    return start


def _rate_matrix(np: Any, value: object, name: str) -> Any:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a square nested list of rates.")
    rows = [
        list(row) if isinstance(row, Sequence) and not isinstance(row, str) else None
        for row in value
    ]
    size = len(rows)
    if any(row is None or len(row) != size for row in rows):
        raise ConfigurationError(f"{name} must be a square nested list of rates.")
    rates = np.array(
        [[_number(item, f"{name} entry") for item in row] for row in rows], dtype=float
    )
    off = rates[~np.eye(size, dtype=bool)]
    if (off < 0).any():
        raise ConfigurationError(f"{name} must have non-negative off-diagonal rates.")
    return rates


def _distribution(np: Any, value: object, size: int, name: str) -> Any:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != size:
        raise ConfigurationError(f"{name} must give one frequency per state ({size} values).")
    freq = np.array([_number(item, f"{name} entry") for item in value], dtype=float)
    if (freq < 0).any() or abs(float(freq.sum()) - 1.0) > 1e-6:
        raise ConfigurationError(f"{name} must be non-negative and sum to 1.")
    return freq


def _matrix(np: Any, value: Any) -> list[list[float | None]]:
    return [[_finite(float(item)) for item in row] for row in np.asarray(value, dtype=float)]


@contextmanager
def _seeded(np: Any, seed: int) -> Iterator[None]:
    """Seed Python's and NumPy's global generators for one call, then restore both."""

    python_state = random.getstate()
    numpy_state = np.random.get_state()
    try:
        random.seed(seed)
        np.random.seed(seed)
        yield
    finally:
        random.setstate(python_state)
        np.random.set_state(numpy_state)


def _trinucleotide_keys() -> list[str]:
    return [left + center + right for left in "ACGT" for center in "ACGT" for right in "ACGT"]


def total_substitution_rate(
    rate_matrix: Sequence[Sequence[float]],
    stationary: Sequence[float],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Total instantaneous substitution rate of a DNA Markov model, using dadi's helper.

    ``rate_matrix`` is square over 4, 16 or 64 states (bases, dinucleotides or trinucleotides).
    ``stationary`` gives the frequency of each state and must sum to 1. The total is the sum over
    states of frequency times exit rate, as dadi's ``Misc.total_instantaneous_rate()`` computes it.
    """

    dadi, np = _dadi_and_numpy()
    rates = _rate_matrix(np, rate_matrix, "rate_matrix")
    size = rates.shape[0]
    if size not in (4, 16, 64):
        raise ConfigurationError("rate_matrix must have 4, 16 or 64 states.")
    freq = _distribution(np, stationary, size, "stationary")
    try:
        total = float(dadi.Misc.total_instantaneous_rate(rates.copy(), freq))
    except (ValueError, TypeError, IndexError) as exc:
        raise ConfigurationError(f"dadi could not compute the total rate: {exc}") from exc
    return {
        "status": "completed",
        "operation": "total_substitution_rate",
        "num_states": size,
        "total_rate": _finite(total),
        "provenance": _provenance("dadi", _RATE_CONVENTION),
    }


def ms_command_line(
    theta: float,
    sample_sizes: Sequence[int],
    core: str,
    iterations: int,
    *,
    recombination: float = 0.0,
    recombination_sites: int | None = None,
    seeds: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """An ``ms`` command line for coalescent sequence simulation, from dadi's ``Misc.ms_command()``.

    ``core`` is the demography part of the command (for example from ``ms_core_command_2d()``).
    ``sample_sizes`` gives one haploid sample size per population. Recombination sites default to
    ``10 * theta`` when recombination is set. The command is only built; ``ms`` is not run.
    """

    dadi, _ = _dadi_and_numpy()
    scale = _number(theta, "theta", positive=True)
    sizes = [
        _integer(item, "sample size", minimum=1, maximum=_MAX_MS_SAMPLE)
        for item in _sequence(sample_sizes, "sample_sizes")
    ]
    if not sizes:
        raise ConfigurationError("sample_sizes must list at least one population.")
    text = _text(core, "core")
    if "\n" in text or "\r" in text:
        raise ConfigurationError("core must be a single line.")
    count = _integer(iterations, "iterations", minimum=1, maximum=_MAX_MS_ITERATIONS)
    rec = _number(recombination, "recombination", minimum=0.0)
    sites = None
    if recombination_sites is not None:
        sites = _integer(recombination_sites, "recombination_sites", minimum=1)
    seed_values = None
    if seeds is not None:
        seed_list = _sequence(seeds, "seeds")
        if len(seed_list) != 3:
            raise ConfigurationError("seeds must give three integers.")
        seed_values = tuple(_integer(item, "seed", minimum=0) for item in seed_list)
    try:
        command = dadi.Misc.ms_command(
            scale, sizes, text, count, recomb=rec, rsites=sites, seeds=seed_values
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"dadi could not build the ms command: {exc}") from exc
    return {
        "status": "completed",
        "operation": "ms_command_line",
        "command": command,
        "theta": scale,
        "sample_sizes": sizes,
        "iterations": count,
        "provenance": _provenance("dadi", "ms command text; no simulation is run"),
    }


def _sequence(value: object, name: str) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list.")
    return list(value)


def ms_core_command_2d(
    model: str,
    params: Sequence[float],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """The ``ms`` core command of a two-population model, using dadi's ``*_mscore`` helpers.

    ``model`` is ``split_mig`` (params ``nu1, nu2, T, m``), ``IM`` (``s, nu1, nu2, T, m12, m21``)
    or ``IM_pre`` (``nuPre, TPre, s, nu1, nu2, T, m12, m21``). The returned text goes to
    ``ms_command_line()`` as ``core``. Times and rates use dadi's scaling, so they are not
    coalescent units.
    """

    dadi, _ = _dadi_and_numpy()
    name = _choice(model, tuple(_MS_CORE_MODELS), "model")
    values = _numbers(params, "params")
    if len(values) != _MS_CORE_MODELS[name]:
        raise ConfigurationError(f"{name} takes {_MS_CORE_MODELS[name]} parameters.")
    if name in ("IM", "IM_pre"):
        share = values[2] if name == "IM_pre" else values[0]
        if not 0 < share < 1:
            raise ConfigurationError("the ancestral share s must lie strictly between 0 and 1.")
    function = getattr(dadi.Demographics2D, f"{name}_mscore")
    try:
        command = function(tuple(values))
    except (ValueError, ArithmeticError, TypeError) as exc:
        raise ConfigurationError(f"dadi could not build the {name} core command: {exc}") from exc
    return {
        "status": "completed",
        "operation": "ms_core_command_2d",
        "model": name,
        "core": command,
        "provenance": _provenance("dadi", "ms core command text; no simulation is run"),
    }


def expected_neutral_sfs_2d(
    sample_sizes: Sequence[int],
    pts: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected two-population neutral spectrum with no divergence, using dadi's ``snm``.

    The populations never diverge, so each is a constant-size neutral sample. ``pts`` is the
    number of grid points used in the integration. Corners that dadi masks are returned as None.
    """

    dadi, np = _dadi_and_numpy()
    ns = _sizes(sample_sizes, 2, 200)
    points = _integer(pts, "pts", minimum=_GRID_MIN, maximum=_GRID_MAX)
    try:
        spectrum = dadi.Demographics2D.snm(None, ns, points)
    except (ValueError, TypeError, IndexError) as exc:
        raise ConfigurationError(f"dadi could not compute the neutral spectrum: {exc}") from exc
    return {
        "status": "completed",
        "operation": "expected_neutral_sfs_2d",
        "sample_sizes": list(ns),
        "pts": points,
        "spectrum": _nested(np, spectrum),
        "provenance": _provenance("dadi", "neutral two-population spectrum; masks as dadi"),
    }


def trinucleotide_misid_table(
    path: str | Path,
    *,
    expected_substitutions: float,
    rate_matrix: Sequence[Sequence[float]],
    trinucleotide_frequencies: Mapping[str, float],
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the 1 - fux table for ancestral misidentification, using dadi's ``make_fux_table()``.

    ``expected_substitutions`` is the expected number of substitutions per site between the ingroup
    and the outgroup. ``rate_matrix`` is the 64 x 64 trinucleotide transition-rate matrix in dadi's
    CGTA order (rate from state i to state j at ``[i][j]``); the diagonal is ignored. The
    frequencies map each of the 64 trinucleotides over A, C, G and T to its ancestral frequency.
    The table is written to a new file, one line per (first, derived, last, ancestral) context.
    """

    dadi, np = _dadi_and_numpy()
    target = _new_file(path, "path")
    span = _number(expected_substitutions, "expected_substitutions", positive=True)
    rates = _rate_matrix(np, rate_matrix, "rate_matrix")
    if rates.shape != (_TRI_STATES, _TRI_STATES):
        raise ConfigurationError("rate_matrix must be 64 x 64 (trinucleotide states).")
    freqs = _trinucleotide_frequencies(trinucleotide_frequencies)
    try:
        dadi.Misc.make_fux_table(str(target), span, rates.copy(), freqs)
    except (ValueError, TypeError, IndexError, KeyError) as exc:
        raise ConfigurationError(
            f"dadi could not write the misidentification table: {exc}"
        ) from exc
    rows = len(target.read_text(encoding="ascii").splitlines())
    return {
        "status": "completed",
        "operation": "trinucleotide_misid_table",
        "path": str(target),
        "rows": rows,
        "expected_substitutions": span,
        "provenance": _provenance("dadi", "1 - fux table from a trinucleotide rate matrix (HWB)"),
    }


def _trinucleotide_frequencies(value: Mapping[str, float]) -> dict[str, float]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(
            "trinucleotide_frequencies must map trinucleotides to frequencies."
        )
    keys = set(_trinucleotide_keys())
    freqs = {}
    for key in value:
        if not isinstance(key, str) or key not in keys:
            raise ConfigurationError(f"{key!r} is not a trinucleotide over A, C, G and T.")
        freqs[key] = _number(value[key], f"frequency of {key}", minimum=0.0)
    missing = sorted(keys - set(freqs))
    if missing:
        raise ConfigurationError(f"trinucleotide_frequencies is missing {missing[:5]}.")
    return {key: freqs[key] for key in sorted(freqs)}


def slim_output_data_dict(
    paths: Sequence[str | Path],
    *,
    mutation_types: Sequence[str] | None = None,
    prefix: str = "SLIM_",
    verbose: bool = True,
) -> dict[str, Any]:
    """Read SLiM sample outputs, one file per population, into a dadi SNP dictionary.

    Uses dadi's ``Misc.dd_from_SLiM_files()``. Population IDs are ``0, 1, ...`` in file order, and
    each site is keyed ``<prefix><position>.<global mutation id>``. ``mutation_types`` keeps only
    the listed SLiM mutation types. Sites are returned sorted by key, with calls as
    ``[ancestral count, derived count]`` per population.
    """

    dadi, _ = _dadi_and_numpy()
    if isinstance(paths, (str, bytes)) or not isinstance(paths, Sequence) or not paths:
        raise ConfigurationError("paths must list one SLiM output file per population.")
    if len(paths) > _MAX_SLIM_FILES:
        raise ConfigurationError(f"paths must list at most {_MAX_SLIM_FILES} files.")
    sources = [input_path(item) for item in paths]
    types = None
    if mutation_types is not None:
        types = [
            _text(item, "mutation type") for item in _sequence(mutation_types, "mutation_types")
        ]
    label = _text(prefix, "prefix")
    with ExitStack() as stack:
        handles = [stack.enter_context(source.open("r", encoding="utf-8")) for source in sources]
        try:
            data, sizes = dadi.Misc.dd_from_SLiM_files(handles, mut_types=types, chr=label)
        except (ValueError, KeyError, IndexError, UnicodeDecodeError) as exc:
            raise ConfigurationError(f"dadi could not parse the SLiM output: {exc}") from exc
    sites = [
        {
            "key": key,
            "outgroup_allele": int(data[key]["outgroup_allele"]),
            "calls": {
                str(pop): [int(calls[0]), int(calls[1])]
                for pop, calls in sorted(data[key]["calls"].items())
            },
        }
        for key in sorted(data)
    ]
    return {
        "status": "completed",
        "operation": "slim_output_data_dict",
        "source": [str(source) for source in sources],
        "populations": [str(index) for index in range(len(sources))],
        "sample_sizes": [int(size) for size in sizes],
        "mutation_types": types,
        "num_sites": len(sites),
        "sites": sites,
        "provenance": _provenance(
            "dadi", "SLiM sample output parsed by dadi Misc.dd_from_SLiM_files"
        ),
    }


def vcf_subsample_bootstrap_spectra(
    vcf_path: str | Path,
    popinfo_path: str | Path,
    *,
    subsample: Mapping[str, int],
    pop_ids: Sequence[str],
    num_boot: int,
    chunk_size: int,
    seed: int = 0,
    filter_sites: bool = True,
    polarized: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Bootstrap spectra from a VCF, drawing a random subset of individuals each time.

    Uses dadi's ``Misc.bootstraps_subsample_vcf()``. ``subsample`` gives the number of diploid
    individuals to keep in each population, and ``pop_ids`` lists one to three populations from
    the population table ``popinfo_path``. The genome is cut into ``chunk_size`` base-pair chunks
    for resampling. ``seed`` fixes both the individual subsets and the chunk draws; the global
    Python and NumPy random states are restored afterwards.
    """

    dadi, np = _dadi_and_numpy()
    vcf = input_path(vcf_path)
    popinfo = input_path(popinfo_path)
    if isinstance(pop_ids, (str, bytes)) or not isinstance(pop_ids, Sequence) or not pop_ids:
        raise ConfigurationError("pop_ids must list one to three populations.")
    selected = [_text(item, "population ID") for item in pop_ids]
    if len(selected) > 3 or len(set(selected)) != len(selected):
        raise ConfigurationError("pop_ids must list one to three distinct populations.")
    if not isinstance(subsample, Mapping):
        raise ConfigurationError("subsample must map each population to a sample size.")
    sizes = {
        _text(pop, "population ID"): _integer(count, "subsample size", minimum=1, maximum=500)
        for pop, count in subsample.items()
    }
    missing = [pop for pop in selected if pop not in sizes]
    if missing:
        raise ConfigurationError(f"subsample is missing populations {missing}.")
    boots = _integer(num_boot, "num_boot", minimum=1, maximum=_MAX_BOOTSTRAPS)
    size = _integer(chunk_size, "chunk_size", minimum=1)
    seed_value = _integer(seed, "seed", minimum=0)
    flag_filter = _flag(filter_sites, "filter_sites")
    flag_polar = _flag(polarized, "polarized")
    with _seeded(np, seed_value):
        try:
            spectra = dadi.Misc.bootstraps_subsample_vcf(
                str(vcf),
                str(popinfo),
                {pop: sizes[pop] for pop in selected},
                boots,
                size,
                selected,
                filter=flag_filter,
                polarized=flag_polar,
            )
        except (ValueError, KeyError, IndexError, OSError) as exc:
            raise ConfigurationError(f"dadi could not bootstrap {vcf.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "vcf_subsample_bootstrap_spectra",
        "source": str(vcf),
        "pop_ids": selected,
        "subsample": {pop: sizes[pop] for pop in selected},
        "num_boot": boots,
        "chunk_size": size,
        "seed": seed_value,
        "spectra": [_nested(np, spectrum) for spectrum in spectra],
        "provenance": _provenance("dadi", "bootstrap of VCF individuals and chunks; seeded"),
    }


def parameter_uncertainty_gim(
    data: Sequence[Any],
    bootstraps: Sequence[Sequence[Any]],
    *,
    model: str,
    p0: Sequence[float],
    grid_points: Sequence[int] = (20, 30, 40),
    log: bool = False,
    multinom: bool = True,
    eps: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parameter standard deviations from the Godambe information matrix (dadi ``GIM_uncert``).

    ``bootstraps`` are spectra of the same shape as ``data``, for example from
    ``vcf_subsample_bootstrap_spectra()`` or ``bootstrap_snp_spectra()``. The model is chosen by
    name and extrapolated over ``grid_points``. With ``multinom=True`` theta is estimated and its
    uncertainty is the last entry. ``log=True`` gives uncertainties on the log scale.
    """

    dadi, np = _dadi_and_numpy()
    observed = _spectrum(dadi, np, data, "data")
    model_fn = _extrapolated_model(dadi, observed, model)
    boots = _bootstrap_spectra(dadi, np, bootstraps, tuple(observed.shape))
    flag_log = _flag(log, "log")
    start = _start(p0, positive=flag_log)
    grids = _grids(grid_points)
    flag_multi = _flag(multinom, "multinom")
    step = _number(eps, "eps", positive=True)
    if step >= 1.0:
        raise ConfigurationError("eps must be below 1.")
    try:
        uncertainties = dadi.Godambe.GIM_uncert(
            model_fn,
            grids,
            boots,
            start,
            observed,
            log=flag_log,
            multinom=flag_multi,
            eps=step,
        )
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi could not compute the GIM uncertainties: {exc}") from exc
    return {
        "status": "completed",
        "operation": "parameter_uncertainty_gim",
        "model": model,
        "params": start,
        "grid_points": grids,
        "num_boot": len(boots),
        "log": flag_log,
        "multinom": flag_multi,
        "uncertainties": [_finite(float(value)) for value in np.atleast_1d(uncertainties)],
        "provenance": _provenance("dadi", "Godambe information from bootstrap spectra"),
    }


def godambe_matrices(
    data: Sequence[Any],
    bootstraps: Sequence[Sequence[Any]],
    *,
    model: str,
    p0: Sequence[float],
    grid_points: Sequence[int] = (20, 30, 40),
    log: bool = False,
    eps: float = 0.01,
    verbose: bool = True,
) -> dict[str, Any]:
    """Godambe information, Hessian and variability matrices, using dadi's ``get_godambe``.

    The model is evaluated as given, with the Poisson log-likelihood used by dadi's Godambe module,
    so include theta among the parameters if the model should estimate it. The Godambe information
    is ``H J^-1 H`` for the Hessian ``H`` and the bootstrap variability ``J``. ``log=True``
    differentiates with respect to log parameters.
    """

    dadi, np = _dadi_and_numpy()
    observed = _spectrum(dadi, np, data, "data")
    model_fn = _extrapolated_model(dadi, observed, model)
    boots = _bootstrap_spectra(dadi, np, bootstraps, tuple(observed.shape))
    flag_log = _flag(log, "log")
    start = _start(p0, positive=flag_log)
    grids = _grids(grid_points)
    step = _number(eps, "eps", positive=True)
    if step >= 1.0:
        raise ConfigurationError("eps must be below 1.")
    try:
        godambe, hessian, variability, _ = dadi.Godambe.get_godambe(
            model_fn, grids, boots, start, observed, step, log=flag_log
        )
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi could not compute the Godambe matrices: {exc}") from exc
    return {
        "status": "completed",
        "operation": "godambe_matrices",
        "model": model,
        "params": start,
        "grid_points": grids,
        "num_boot": len(boots),
        "log": flag_log,
        "godambe_information": _matrix(np, godambe),
        "hessian": _matrix(np, hessian),
        "variability": _matrix(np, variability),
        "provenance": _provenance("dadi", "Godambe information H J^-1 H from bootstrap spectra"),
    }


def nested_model_statistic(
    data: Sequence[Any],
    bootstraps: Sequence[Sequence[Any]],
    *,
    test: str,
    model: str,
    p0: Sequence[float],
    nested_indices: Sequence[int],
    grid_points: Sequence[int] = (20, 30, 40),
    full_params: Sequence[float] | None = None,
    multinom: bool = True,
    eps: float = 0.01,
    adjusted_and_original: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nested-model test with bootstrap adjustment, using dadi's Godambe module.

    ``test`` selects ``lrt`` (the first-order adjustment factor for the likelihood-ratio test),
    ``wald`` (adjusted Wald statistic) or ``score`` (adjusted score statistic). ``model`` is the
    complex model. ``p0`` gives its best-fit parameters with the nested parameters set to their
    nested values, and ``nested_indices`` gives their positions in the parameter list. The Wald
    test also needs ``full_params``, the complex model's parameter values. ``adjusted_and_original``
    also returns the unadjusted statistic for ``wald`` and ``score``.
    """

    dadi, np = _dadi_and_numpy()
    kind = _choice(test, _NESTED_TESTS, "test")
    observed = _spectrum(dadi, np, data, "data")
    model_fn = _extrapolated_model(dadi, observed, model)
    boots = _bootstrap_spectra(dadi, np, bootstraps, tuple(observed.shape))
    start = _start(p0, positive=False)
    indices = _nested_indices(nested_indices, len(start))
    grids = _grids(grid_points)
    flag_multi = _flag(multinom, "multinom")
    flag_pair = _flag(adjusted_and_original, "adjusted_and_original")
    step = _number(eps, "eps", positive=True)
    if step >= 1.0:
        raise ConfigurationError("eps must be below 1.")
    if flag_pair and kind == "lrt":
        raise ConfigurationError("adjusted_and_original applies to the wald and score tests only.")
    try:
        if kind == "lrt":
            adjustment = dadi.Godambe.LRT_adjust(
                model_fn, grids, boots, start, observed, indices, multinom=flag_multi, eps=step
            )
            result: dict[str, Any] = {"adjustment": _finite(float(adjustment))}
        elif kind == "wald":
            if full_params is None:
                raise ConfigurationError("the wald test needs full_params.")
            full = _numbers(full_params, "full_params")
            value = dadi.Godambe.Wald_stat(
                model_fn,
                grids,
                boots,
                start,
                observed,
                indices,
                full,
                multinom=flag_multi,
                eps=step,
                adj_and_org=flag_pair,
            )
            result = _statistic_pair(value, flag_pair)
        else:
            value = dadi.Godambe.score_stat(
                model_fn,
                grids,
                boots,
                start,
                observed,
                indices,
                multinom=flag_multi,
                eps=step,
                adj_and_org=flag_pair,
            )
            result = _statistic_pair(value, flag_pair)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi could not compute the {kind} statistic: {exc}") from exc
    return {
        "status": "completed",
        "operation": "nested_model_statistic",
        "test": kind,
        "model": model,
        "params": start,
        "nested_indices": indices,
        "grid_points": grids,
        "num_boot": len(boots),
        **result,
        "provenance": _provenance("dadi", "nested-model test with bootstrap moment adjustment"),
    }


def _nested_indices(value: object, count: int) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("nested_indices must list one or more parameter positions.")
    indices = [_integer(item, "nested index", minimum=0, maximum=count - 1) for item in value]
    if len(set(indices)) != len(indices) or len(indices) >= count:
        raise ConfigurationError("nested_indices must be distinct and leave one free parameter.")
    return indices


def _statistic_pair(value: Any, both: bool) -> dict[str, Any]:
    if both:
        adjusted, original = value
        return {"statistic": _finite(float(adjusted)), "unadjusted": _finite(float(original))}
    return {"statistic": _finite(float(value))}


def fit_demographic_model_to_residual(
    data: Sequence[Any],
    target_residual: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    grid_points: Sequence[int] = (20, 30, 40),
    lower_bound: Sequence[float] | None = None,
    upper_bound: Sequence[float] | None = None,
    multinom: bool = True,
    max_iterations: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model so that its residual spectrum matches a target, using dadi's
    ``Inference.optimize_log_resid()``.

    ``target_residual`` is the residual spectrum to match, of the same shape as ``data``, for
    example from a synonymous-site fit. The optimizer is BFGS on log parameters, and it stops after
    ``max_iterations`` steps at most. The result gives the fitted parameters and the objective.
    """

    dadi, np = _dadi_and_numpy()
    observed = _spectrum(dadi, np, data, "data")
    target = _spectrum(dadi, np, target_residual, "target_residual")
    if tuple(target.shape) != tuple(observed.shape):
        raise ConfigurationError("target_residual must have the shape of data.")
    model_fn = _extrapolated_model(dadi, observed, model)
    start = _start(p0)
    grids = _grids(grid_points)
    lower = _bounds(lower_bound, len(start), "lower_bound")
    upper = _bounds(upper_bound, len(start), "upper_bound")
    flag_multi = _flag(multinom, "multinom")
    iterations = _integer(max_iterations, "max_iterations", minimum=1, maximum=10_000)
    try:
        xopt, fopt, _, _, func_calls, grad_calls, warnflag = dadi.Inference.optimize_log_resid(
            start,
            observed,
            model_fn,
            target,
            grids,
            lower_bound=lower,
            upper_bound=upper,
            verbose=0,
            multinom=flag_multi,
            maxiter=iterations,
            full_output=True,
        )
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi could not fit the residual model: {exc}") from exc
    return {
        "status": "completed",
        "operation": "fit_demographic_model_to_residual",
        "model": model,
        "params": [_finite(float(value)) for value in np.atleast_1d(xopt)],
        "objective": _finite(float(fopt)),
        "function_calls": int(func_calls),
        "gradient_calls": int(grad_calls),
        "converged": int(warnflag) == 0,
        "grid_points": grids,
        "multinom": flag_multi,
        "provenance": _provenance("dadi", "BFGS on log parameters against a residual target"),
    }


__all__ = [
    "expected_neutral_sfs_2d",
    "fit_demographic_model_to_residual",
    "godambe_matrices",
    "ms_command_line",
    "ms_core_command_2d",
    "nested_model_statistic",
    "parameter_uncertainty_gim",
    "slim_output_data_dict",
    "total_substitution_rate",
    "trinucleotide_misid_table",
    "vcf_subsample_bootstrap_spectra",
]

_register_public_verbose(globals(), tuple(__all__))
