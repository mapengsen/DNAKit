"""dadi wrappers for family-variant rows that no existing dnakit wrapper reaches (group H1).

The wrappers here cover dadi functions that the earlier ledger called variants of a wrapped
function, although the wrapped function cannot produce their result:

- the log-likelihood family of ``Inference`` (multinomial, per bin, negated) for any model
  spectrum, while ``calculate_sfs_log_likelihood()`` only scores the neutral model;
- the optimally scaled model spectrum;
- the optimizers that work on natural-scale parameters (BFGS), under linear constraints, on a
  parameter grid, and with NLopt, while ``fit_demographic_model()`` only uses log-scale fits.

Model functions are chosen by name from the fixed lists in ``dadi_functions``, and constraints are
linear, so no arbitrary code is run. Every input is in memory. dadi, NumPy, SciPy and NLopt are
imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _finite, _integer, _number, _provenance
from dnakit.population.dadi_functions import (
    _DEMOGRAPHICS_1D,
    _DEMOGRAPHICS_2D,
    _MAX_MAXITER,
    _bounds,
    _choice,
    _counts_array,
    _dadi_and_numpy,
    _flag,
    _model_function,
    _nested,
    _numbers,
    _spectrum_sizes,
)

_NLOPT_ALGORITHMS = ("LN_BOBYQA", "LN_COBYLA", "LN_NELDERMEAD", "LN_SBPLX")
_MAX_GRID_EVALUATIONS = 2000
_MAX_NLOPT_EVALUATIONS = 100_000
_MAX_CONSTRAINTS = 50
_INF = float("inf")


def _same_shape(np: Any, model: object, data: object) -> tuple[Any, Any]:
    model_array = _counts_array(np, model, "model", nonnegative=False)
    data_array = _counts_array(np, data, "data")
    if model_array.shape != data_array.shape:
        raise ConfigurationError("model and data must have the same shape.")
    return model_array, data_array


def calculate_model_log_likelihood(
    model: Sequence[Any],
    data: Sequence[Any],
    *,
    multinom: bool = False,
    per_bin: bool = False,
    negate: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Log-likelihood of observed counts given any model spectrum, using dadi's ``Inference``.

    ``model`` and ``data`` are spectra of the same shape (one to three populations). The default
    is the Poisson log-likelihood ``Inference.ll``. ``multinom=True`` first rescales the model to
    the optimal theta (``ll_multinom``). ``per_bin=True`` returns the value of every entry
    (``ll_per_bin``, ``ll_multinom_per_bin``) instead of their sum, with masked entries as None.
    ``negate=True`` returns the negative total (``minus_ll``, ``minus_ll_multinom``) and cannot
    be combined with ``per_bin``.
    """

    dadi, np = _dadi_and_numpy()
    multi = _flag(multinom, "multinom")
    each = _flag(per_bin, "per_bin")
    flip = _flag(negate, "negate")
    if each and flip:
        raise ConfigurationError("negate applies to the total log-likelihood, not to per_bin.")
    model_array, data_array = _same_shape(np, model, data)
    fs_model = dadi.Spectrum(model_array)
    fs_data = dadi.Spectrum(data_array)
    if each:
        name = "ll_multinom_per_bin" if multi else "ll_per_bin"
    elif flip:
        name = "minus_ll_multinom" if multi else "minus_ll"
    else:
        name = "ll_multinom" if multi else "ll"
    with np.errstate(all="ignore"):
        value = getattr(dadi.Inference, name)(fs_model, fs_data)
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "calculate_model_log_likelihood",
        "dadi_function": name,
        "multinom": multi,
        "per_bin": each,
        "negate": flip,
        "shape": [int(n) for n in data_array.shape],
    }
    if each:
        masked = np.ma.masked_invalid(np.ma.asarray(value, dtype=float))
        result["log_likelihood_per_bin"] = _nested(np, masked)
    else:
        result["log_likelihood"] = _finite(float(value))
    result["provenance"] = _provenance("dadi", f"dadi Inference.{name}; masked entries skipped")
    return result


def scale_model_sfs_to_data(
    model: Sequence[Any],
    data: Sequence[Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a model spectrum by the maximum-likelihood theta for the data, using dadi.

    Returns the scaled spectrum (dadi's ``Inference.optimally_scaled_sfs``) and the factor
    ``theta`` (``Inference.optimal_sfs_scaling``, the data total over the model total on the
    unmasked entries). ``model`` and ``data`` must have the same shape.
    """

    dadi, np = _dadi_and_numpy()
    model_array, data_array = _same_shape(np, model, data)
    fs_model = dadi.Spectrum(model_array)
    fs_data = dadi.Spectrum(data_array)
    with np.errstate(all="ignore"):
        scaled = dadi.Inference.optimally_scaled_sfs(fs_model, fs_data)
        theta = dadi.Inference.optimal_sfs_scaling(fs_model, fs_data)
    return {
        "status": "completed",
        "operation": "scale_model_sfs_to_data",
        "shape": [int(n) for n in data_array.shape],
        "theta": _finite(float(theta)),
        "scaled_spectrum": _nested(np, scaled),
        "provenance": _provenance("dadi", "model spectrum times the optimal theta (Poisson)"),
    }


class _Fit:
    """Observed spectrum, model function and validated options shared by the fit wrappers."""

    def __init__(
        self,
        dadi: Any,
        np: Any,
        data: object,
        model: str,
        p0: object,
        lower_bound: object,
        upper_bound: object,
        multinom: object,
        grid_points: object,
    ) -> None:
        counts = _counts_array(np, data, "data", dims=(1, 2))
        self.observed = dadi.Spectrum(counts)
        self.sizes = tuple(_spectrum_sizes(self.observed))
        family, options = (
            (dadi.Demographics1D, _DEMOGRAPHICS_1D)
            if counts.ndim == 1
            else (dadi.Demographics2D, _DEMOGRAPHICS_2D)
        )
        self.function = _model_function(family, model, options)
        self.model = model
        self.start = _numbers(p0, "p0")
        if any(value <= 0 for value in self.start):
            raise ConfigurationError("p0 must contain only positive values.")
        self.lower = _bounds(lower_bound, len(self.start), "lower_bound")
        self.upper = _bounds(upper_bound, len(self.start), "upper_bound")
        if (
            self.lower is not None
            and self.upper is not None
            and any(lo >= hi for lo, hi in zip(self.lower, self.upper, strict=True))
        ):
            raise ConfigurationError("each lower_bound must be below its upper_bound.")
        self.multinom = _flag(multinom, "multinom")
        self.points = _integer(grid_points, "grid_points", minimum=20, maximum=400)

    def log_likelihood(self, dadi: Any, params: Any) -> float | None:
        spectrum = self.function(tuple(float(value) for value in params), self.sizes, self.points)
        score = dadi.Inference.ll_multinom if self.multinom else dadi.Inference.ll
        return _finite(float(score(spectrum, self.observed)))


def _fit_failure(model: str, exc: Exception) -> ConfigurationError:
    return ConfigurationError(f"dadi could not fit {model!r} to the data: {exc}")


def fit_demographic_model_natural_scale(
    data: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    lower_bound: Sequence[float] | None = None,
    upper_bound: Sequence[float] | None = None,
    multinom: bool = True,
    grid_points: int = 40,
    maxiter: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model by BFGS on the parameters themselves, using ``Inference.optimize``.

    ``fit_demographic_model()`` optimizes the logarithm of the parameters; this function keeps
    the natural scale, so parameters may approach their bounds linearly and ``p0`` need not be
    far from zero. ``data`` is a one- or two-population spectrum, ``model`` a model of that
    family, and ``p0`` positive starting values. With ``multinom=True`` theta is rescaled at every
    step and is not part of ``p0``.
    """

    dadi, np = _dadi_and_numpy()
    fit = _Fit(dadi, np, data, model, p0, lower_bound, upper_bound, multinom, grid_points)
    iterations = _integer(maxiter, "maxiter", minimum=1, maximum=_MAX_MAXITER)
    try:
        output = dadi.Inference.optimize(
            fit.start,
            fit.observed,
            fit.function,
            fit.points,
            lower_bound=fit.lower,
            upper_bound=fit.upper,
            verbose=0,
            multinom=fit.multinom,
            maxiter=iterations,
            full_output=True,
        )
        fitted = output[0]
        log_likelihood = fit.log_likelihood(dadi, fitted)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise _fit_failure(model, exc) from exc
    return {
        "status": "completed",
        "operation": "fit_demographic_model_natural_scale",
        "model": model,
        "multinom": fit.multinom,
        "sample_sizes": list(fit.sizes),
        "grid_points": fit.points,
        "maxiter": iterations,
        "params": [float(value) for value in np.atleast_1d(fitted)],
        "log_likelihood": log_likelihood,
        "function_calls": int(output[4]),
        "converged": int(output[6]) == 0,
        "provenance": _provenance("dadi", "maximum likelihood; dadi Inference.optimize (BFGS)"),
    }


def _linear_constraints(value: object, count: int, name: str) -> list[tuple[list[float], float]]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of constraint mappings.")
    if len(value) > _MAX_CONSTRAINTS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_CONSTRAINTS} constraints.")
    rows = []
    for item in value:
        if not isinstance(item, Mapping) or set(item) - {"coefficients", "offset"}:
            raise ConfigurationError(
                f"each {name} entry needs 'coefficients' and optional 'offset'."
            )
        raw = item.get("coefficients")
        if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence) or len(raw) != count:
            raise ConfigurationError(
                f"{name} coefficients need {count} entries, one per parameter."
            )
        coefficients = [_number(c, f"{name} coefficient", minimum=-_INF) for c in raw]
        offset = _number(item.get("offset", 0.0), f"{name} offset", minimum=-_INF)
        rows.append((coefficients, offset))
    return rows


def _constraint_function(rows: list[tuple[list[float], float]], np: Any) -> Any:
    """``p -> [c . p + offset, ...]``; dadi passes extra arguments, which are ignored."""

    def evaluate(params: Any, *_ignored: Any) -> Any:
        values = np.asarray(params, dtype=float)
        return np.array([float(np.dot(c, values)) + offset for c, offset in rows])

    return evaluate


def fit_demographic_model_constrained(
    data: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    lower_bound: Sequence[float] | None = None,
    upper_bound: Sequence[float] | None = None,
    equality_constraints: Sequence[Mapping[str, Any]] | None = None,
    inequality_constraints: Sequence[Mapping[str, Any]] | None = None,
    multinom: bool = True,
    grid_points: int = 40,
    maxiter: int = 50,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model under linear parameter constraints, using ``optimize_cons``.

    Each constraint is a mapping ``{"coefficients": [c0, c1, ...], "offset": b}`` for the value
    ``c . params + b``. An equality constraint must equal 0 at the solution (for example
    ``nu1 = 2 * nu2`` is ``{"coefficients": [1, -2], "offset": 0}``); an inequality constraint
    must be at least 0 (``nu1 >= nu2`` is ``{"coefficients": [1, -1]}``). The optimizer is
    SLSQP on the natural-scale parameters. The constraint values at the solution are reported.
    """

    dadi, np = _dadi_and_numpy()
    fit = _Fit(dadi, np, data, model, p0, lower_bound, upper_bound, multinom, grid_points)
    iterations = _integer(maxiter, "maxiter", minimum=1, maximum=_MAX_MAXITER)
    count = len(fit.start)
    equalities = _linear_constraints(equality_constraints, count, "equality_constraints")
    inequalities = _linear_constraints(inequality_constraints, count, "inequality_constraints")
    if not equalities and not inequalities:
        raise ConfigurationError("give at least one equality or inequality constraint.")
    try:
        output = dadi.Inference.optimize_cons(
            fit.start,
            fit.observed,
            fit.function,
            fit.points,
            eq_constraint=_constraint_function(equalities, np) if equalities else None,
            ieq_constraint=_constraint_function(inequalities, np) if inequalities else None,
            lower_bound=fit.lower,
            upper_bound=fit.upper,
            verbose=0,
            multinom=fit.multinom,
            maxiter=iterations,
            full_output=True,
        )
        fitted = np.atleast_1d(output[0])
        log_likelihood = fit.log_likelihood(dadi, fitted)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise _fit_failure(model, exc) from exc
    return {
        "status": "completed",
        "operation": "fit_demographic_model_constrained",
        "model": model,
        "multinom": fit.multinom,
        "sample_sizes": list(fit.sizes),
        "grid_points": fit.points,
        "maxiter": iterations,
        "params": [float(value) for value in fitted],
        "log_likelihood": log_likelihood,
        "equality_values": [float(v) for v in _constraint_function(equalities, np)(fitted)]
        if equalities
        else [],
        "inequality_values": [float(v) for v in _constraint_function(inequalities, np)(fitted)]
        if inequalities
        else [],
        "iterations": int(output[2]),
        "converged": int(output[3]) == 0,
        "provenance": _provenance(
            "dadi", "maximum likelihood; dadi Inference.optimize_cons (SLSQP)"
        ),
    }


def _grid_ranges(value: object, count: int) -> list[tuple[float, float, int]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != count:
        raise ConfigurationError(
            f"parameter_ranges must list {count} ranges [start, stop, num_points]."
        )
    ranges = []
    total = 1
    for item in value:
        if isinstance(item, (str, bytes)) or not isinstance(item, Sequence) or len(item) != 3:
            raise ConfigurationError("each range must be [start, stop, num_points].")
        start = _number(item[0], "range start", positive=True)
        stop = _number(item[1], "range stop", positive=True)
        points = _integer(item[2], "range num_points", minimum=2, maximum=_MAX_GRID_EVALUATIONS)
        if stop <= start:
            raise ConfigurationError("each range stop must be above its start.")
        total *= points
        ranges.append((start, stop, points))
    if total > _MAX_GRID_EVALUATIONS:
        raise ConfigurationError(
            f"the grid has {total} points, above the limit of {_MAX_GRID_EVALUATIONS}."
        )
    return ranges


def fit_demographic_model_grid(
    data: Sequence[Any],
    *,
    model: str,
    parameter_ranges: Sequence[Sequence[float]],
    multinom: bool = True,
    grid_points: int = 40,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model by exhaustive search over a parameter grid, using ``optimize_grid``.

    ``parameter_ranges`` has one ``[start, stop, num_points]`` entry per model parameter; the
    grid holds ``num_points`` evenly spaced values from ``start`` to ``stop`` inclusive, and every
    combination is evaluated (at most 2000). The best grid point and its log-likelihood are
    returned. Useful to find starting values before a local fit.
    """

    dadi, np = _dadi_and_numpy()
    counts = _counts_array(np, data, "data", dims=(1, 2))
    observed = dadi.Spectrum(counts)
    family, options = (
        (dadi.Demographics1D, _DEMOGRAPHICS_1D)
        if counts.ndim == 1
        else (dadi.Demographics2D, _DEMOGRAPHICS_2D)
    )
    function = _model_function(family, model, options)
    flag = _flag(multinom, "multinom")
    points = _integer(grid_points, "grid_points", minimum=20, maximum=400)
    if isinstance(parameter_ranges, (str, bytes)) or not isinstance(parameter_ranges, Sequence):
        raise ConfigurationError("parameter_ranges must be a list of [start, stop, num_points].")
    ranges = _grid_ranges(parameter_ranges, len(parameter_ranges))
    grid = tuple(slice(start, stop, complex(0, num)) for start, stop, num in ranges)
    ns = tuple(_spectrum_sizes(observed))
    try:
        best = dadi.Inference.optimize_grid(
            observed, function, points, grid, verbose=0, multinom=flag
        )
        best = [float(value) for value in np.atleast_1d(best)]
        spectrum = function(tuple(best), ns, points)
        score = dadi.Inference.ll_multinom if flag else dadi.Inference.ll
        log_likelihood = _finite(float(score(spectrum, observed)))
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise _fit_failure(model, exc) from exc
    total = 1
    for _start, _stop, num in ranges:
        total *= num
    return {
        "status": "completed",
        "operation": "fit_demographic_model_grid",
        "model": model,
        "multinom": flag,
        "sample_sizes": list(ns),
        "grid_points": points,
        "parameter_ranges": [list(item) for item in ranges],
        "num_evaluations": total,
        "params": best,
        "log_likelihood": log_likelihood,
        "provenance": _provenance("dadi", "maximum likelihood on a grid; dadi optimize_grid"),
    }


def fit_demographic_model_nlopt(
    data: Sequence[Any],
    *,
    model: str,
    p0: Sequence[float],
    lower_bound: Sequence[float] | None = None,
    upper_bound: Sequence[float] | None = None,
    algorithm: str = "LN_BOBYQA",
    multinom: bool = True,
    grid_points: int = 40,
    maxeval: int = 1000,
    ftol_abs: float = 1e-6,
    xtol_abs: float = 1e-6,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a demographic model with a derivative-free NLopt algorithm, using ``Inference.opt``.

    ``algorithm`` is ``LN_BOBYQA`` (default), ``LN_COBYLA``, ``LN_NELDERMEAD`` or ``LN_SBPLX``.
    The search stops after ``maxeval`` model evaluations or when the log-likelihood changes by
    less than ``ftol_abs`` or the parameters by less than ``xtol_abs``. The fit works on the
    natural-scale parameters and needs the ``nlopt`` package.
    """

    dadi, np = _dadi_and_numpy()
    nlopt = require("nlopt", _DADI_EXTRA)
    fit = _Fit(dadi, np, data, model, p0, lower_bound, upper_bound, multinom, grid_points)
    name = _choice(algorithm, _NLOPT_ALGORITHMS, "algorithm")
    evaluations = _integer(maxeval, "maxeval", minimum=1, maximum=_MAX_NLOPT_EVALUATIONS)
    f_tolerance = _number(ftol_abs, "ftol_abs", positive=True)
    x_tolerance = _number(xtol_abs, "xtol_abs", positive=True)
    try:
        fitted, value = dadi.Inference.opt(
            fit.start,
            fit.observed,
            fit.function,
            fit.points,
            multinom=fit.multinom,
            lower_bound=fit.lower,
            upper_bound=fit.upper,
            algorithm=getattr(nlopt, name),
            local_optimizer=getattr(nlopt, name),
            ftol_abs=f_tolerance,
            xtol_abs=x_tolerance,
            maxeval=evaluations,
            verbose=0,
        )
        log_likelihood = fit.log_likelihood(dadi, fitted)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError, RuntimeError) as exc:
        raise _fit_failure(model, exc) from exc
    return {
        "status": "completed",
        "operation": "fit_demographic_model_nlopt",
        "model": model,
        "algorithm": name,
        "multinom": fit.multinom,
        "sample_sizes": list(fit.sizes),
        "grid_points": fit.points,
        "maxeval": evaluations,
        "params": [float(item) for item in np.atleast_1d(fitted)],
        "log_likelihood": log_likelihood,
        "optimizer_value": _finite(float(value)),
        "provenance": _provenance("dadi", f"maximum likelihood; dadi Inference.opt ({name})"),
    }


__all__ = [
    "calculate_model_log_likelihood",
    "fit_demographic_model_constrained",
    "fit_demographic_model_grid",
    "fit_demographic_model_natural_scale",
    "fit_demographic_model_nlopt",
    "scale_model_sfs_to_data",
]

_register_public_verbose(globals(), tuple(__all__))
