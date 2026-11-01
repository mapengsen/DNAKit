"""dadi population-density (phi) wrappers: phi construction, demographic events and spectra.

dadi's diffusion code keeps a population density ``phi`` on a frequency grid, applies events
(population split, admixture, removal) to it, and finally samples an allele-frequency spectrum from
it. These functions expose that pipeline one step at a time. ``phi`` is passed in and returned as a
nested list of floats (one axis per population) together with the number of grid points; dadi's
``Numerics.default_grid(grid_points)`` is the frequency grid on every axis. A large density can be
written to a new ``.npy`` file with ``output_path`` and read back by passing the path as ``phi``.

``spectrum_from_snp_data`` builds a spectrum from a dadi SNP file (dadi reads it into its data
dictionary) and can apply the Hernandez-Williamson-Bustamante ancestral-misidentification
correction. dadi, NumPy are imported only when a function runs. No network is used. Output files
are never overwritten.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA, _integer, _number, _provenance
from dnakit.population.dadi_functions import _load_snp, _nested, _new_file, _populations
from dnakit.population.dadi_review_functions import _seeded, _spectrum

_MIN_GRID = 8
_MAX_GRID = 400
_MAX_INLINE_VALUES = 200_000
_MAX_VALUES = 50_000_000
_MAX_SAMPLE = 200
_ONE_POP_MODELS = ("dominance", "genic", "neutral", "x_linked")
_SPLIT_EVENTS = ("1_to_2", "2_to_3_split_2", "2_to_3_split_1")
_SELECT_OPERATIONS = ("remove", "keep", "reorder")
_PHI_CONVENTION = (
    "phi is a nested list with one axis per population on dadi's default frequency grid of "
    "grid_points points; events follow dadi.PhiManip"
)


def _dadi_and_numpy() -> tuple[Any, Any]:
    return require("dadi", _DADI_EXTRA), require("numpy", _DADI_EXTRA)


def _grid_points(value: object) -> int:
    return _integer(value, "grid_points", minimum=_MIN_GRID, maximum=_MAX_GRID)


def _grid(dadi: Any, points: int) -> Any:
    return dadi.Numerics.default_grid(points)


def _load_phi(np: Any, value: object, points: int, dims: tuple[int, ...] | None = None) -> Any:
    """phi as an array of shape (points,) * ndim, from a nested list or a ``.npy`` path."""

    if isinstance(value, (str, Path)):
        source = input_path(value)
        try:
            array = np.load(source, allow_pickle=False)
        except (OSError, ValueError) as exc:
            raise ConfigurationError(f"phi file {source.name} is not a readable .npy: {exc}") from exc
    elif isinstance(value, Sequence) and not isinstance(value, bytes):
        try:
            array = np.array(value, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("phi must be a rectangular nested list of numbers.") from exc
    else:
        raise ConfigurationError("phi must be a nested list of numbers or the path of a .npy file.")
    array = np.asarray(array, dtype=float)
    if array.ndim < 1 or array.ndim > 5 or any(n != points for n in array.shape):
        raise ConfigurationError(
            f"phi must have 1 to 5 axes of length grid_points={points}; got shape {array.shape}."
        )
    if dims is not None and array.ndim not in dims:
        raise ConfigurationError(f"phi must have {sorted(dims)} populations; got {array.ndim}.")
    if not np.isfinite(array).all():
        raise ConfigurationError("phi must contain only finite values.")
    return array


def _phi_result(
    np: Any, operation: str, phi: Any, points: int, output_path: object, extra: Mapping[str, Any]
) -> dict[str, Any]:
    phi = np.asarray(phi, dtype=float)
    if not np.isfinite(phi).all():
        raise ConfigurationError(f"dadi produced a non-finite phi in {operation}.")
    if phi.size > _MAX_VALUES:
        raise ConfigurationError(f"phi has {phi.size} values; the limit is {_MAX_VALUES}.")
    saved = None
    if output_path is not None:
        target = _new_file(output_path, "output_path")
        if target.suffix.lower() != ".npy":
            raise ConfigurationError("output_path must end with .npy.")
        np.save(target, phi)
        saved = str(target)
    elif phi.size > _MAX_INLINE_VALUES:
        raise ConfigurationError(
            f"phi has {phi.size} values; pass output_path (a new .npy file) to keep it."
        )
    result: dict[str, Any] = {
        "status": "completed",
        "operation": operation,
        "grid_points": points,
        "num_populations": int(phi.ndim),
        "shape": [int(n) for n in phi.shape],
        "phi": phi.tolist() if phi.size <= _MAX_INLINE_VALUES else None,
        "phi_path": saved,
        "phi_sum": float(phi.sum()),
        "provenance": _provenance("dadi", _PHI_CONVENTION),
    }
    result.update(extra)
    return result


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except (ValueError, TypeError, IndexError, KeyError, ZeroDivisionError) as exc:
        raise ConfigurationError(f"dadi {function.__name__} rejected the input: {exc}") from exc


def _fractions(value: object, count: int, name: str) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != count:
        raise ConfigurationError(f"{name} must list {count} admixture proportions.")
    numbers = [_number(item, f"{name} entry") for item in value]
    if any(item > 1 for item in numbers):
        raise ConfigurationError(f"{name} entries must be between 0 and 1.")
    if sum(numbers) > 1 + 1e-12:
        raise ConfigurationError(f"{name} entries must add up to at most 1.")
    return numbers


def phi_one_population(
    *,
    grid_points: int,
    model: str = "dominance",
    nu: float = 1.0,
    theta0: float = 1.0,
    gamma: float = 0.0,
    h: float = 0.5,
    theta: float | None = None,
    beta: float = 1.0,
    alpha: float = 1.0,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the equilibrium population density of one population, using ``dadi.PhiManip``.

    ``model`` selects the dadi routine: ``dominance`` is ``phi_1D`` (numerical quadrature, any
    dominance ``h``, optional ``theta`` for asymmetric mutation), ``genic`` is ``phi_1D_genic``
    (analytic series, ``h`` = 0.5), ``neutral`` is ``phi_1D_snm`` (standard neutral model) and
    ``x_linked`` is ``phi_1D_X`` (X-linked locus, with male/female ``alpha``). ``nu`` is the size,
    ``theta0`` the scaled mutation rate, ``gamma`` the scaled selection coefficient and ``beta``
    the selection-coefficient ratio for the genic variants. Parameters a model does not take are
    ignored, as in dadi.
    """

    dadi, np = _dadi_and_numpy()
    chosen = model if model in _ONE_POP_MODELS else None
    if chosen is None:
        raise ConfigurationError(f"model must be one of {sorted(_ONE_POP_MODELS)}.")
    points = _grid_points(grid_points)
    values = {
        "nu": _number(nu, "nu", positive=True),
        "theta0": _number(theta0, "theta0", positive=True),
        "gamma": _number(gamma, "gamma", minimum=-1e6),
    }
    beta_value = _number(beta, "beta", positive=True)
    xx = _grid(dadi, points)
    manip = dadi.PhiManip
    if chosen == "dominance":
        theta_value = None if theta is None else _number(theta, "theta", positive=True)
        phi = _call(
            manip.phi_1D, xx, nu=values["nu"], theta0=values["theta0"], gamma=values["gamma"],
            h=_number(h, "h"), theta=theta_value, beta=beta_value,
        )
    elif chosen == "genic":
        theta_value = None if theta is None else _number(theta, "theta", positive=True)
        phi = _call(
            manip.phi_1D_genic, xx, nu=values["nu"], theta0=values["theta0"],
            gamma=values["gamma"], theta=theta_value, beta=beta_value,
        )
    elif chosen == "neutral":
        theta_value = None if theta is None else _number(theta, "theta", positive=True)
        phi = _call(
            manip.phi_1D_snm, xx, nu=values["nu"], theta0=values["theta0"], theta=theta_value,
            beta=beta_value,
        )
    else:
        phi = _call(
            manip.phi_1D_X, xx, nu=values["nu"], theta0=values["theta0"], gamma=values["gamma"],
            h=_number(h, "h"), beta=beta_value, alpha=_number(alpha, "alpha", positive=True),
        )
    return _phi_result(np, "phi_one_population", phi, points, output_path, {"model": chosen})


def phi_split_population(
    phi: Sequence[Any] | str | Path,
    *,
    grid_points: int,
    event: str,
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split one population of a density into two, using ``dadi.PhiManip``.

    ``event`` selects the split: ``1_to_2`` (``phi_1D_to_2D``: a one-population density becomes
    a two-population density), ``2_to_3_split_2`` (``phi_2D_to_3D_split_2``: population 2 splits
    into populations 2 and 3) or ``2_to_3_split_1`` (``phi_2D_to_3D_split_1``: population 1
    splits into 1 and 3). The new population starts as a copy of the parent's frequencies.
    """

    dadi, np = _dadi_and_numpy()
    if event not in _SPLIT_EVENTS:
        raise ConfigurationError(f"event must be one of {sorted(_SPLIT_EVENTS)}.")
    points = _grid_points(grid_points)
    ndim_in = 1 if event == "1_to_2" else 2
    density = _load_phi(np, phi, points, (ndim_in,))
    xx = _grid(dadi, points)
    function = {
        "1_to_2": dadi.PhiManip.phi_1D_to_2D,
        "2_to_3_split_2": dadi.PhiManip.phi_2D_to_3D_split_2,
        "2_to_3_split_1": dadi.PhiManip.phi_2D_to_3D_split_1,
    }[event]
    result = _call(function, xx, density)
    return _phi_result(np, "phi_split_population", result, points, output_path, {"event": event})


def phi_create_admixed_population(
    phi: Sequence[Any] | str | Path,
    *,
    grid_points: int,
    proportions: Sequence[float],
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a new population formed by admixture of all existing ones, using ``dadi.PhiManip``.

    A 2-population density gains population 3 (``phi_2D_to_3D_admix``), a 3-population density
    gains population 4 (``phi_3D_to_4D``) and a 4-population density gains population 5
    (``phi_4D_to_5D``). ``proportions`` lists the contribution of populations 1 to n-1 (n = the
    number of populations of ``phi``); the last population contributes the remainder.
    """

    dadi, np = _dadi_and_numpy()
    points = _grid_points(grid_points)
    density = _load_phi(np, phi, points, (2, 3, 4))
    count = density.ndim
    f = _fractions(proportions, count - 1, "proportions")
    grids = [_grid(dadi, points)] * (count + 1)
    function = {
        2: dadi.PhiManip.phi_2D_to_3D_admix,
        3: dadi.PhiManip.phi_3D_to_4D,
        4: dadi.PhiManip.phi_4D_to_5D,
    }[count]
    result = _call(function, density, *f, *grids)
    return _phi_result(
        np, "phi_create_admixed_population", result, points, output_path,
        {"proportions": f},
    )


def phi_admix_into_population(
    phi: Sequence[Any] | str | Path,
    *,
    grid_points: int,
    target_population: int,
    proportions: Sequence[float],
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace one population by a mixture of the other populations, using ``dadi.PhiManip``.

    Covers dadi's ``phi_2D_admix_1_into_2`` / ``phi_2D_admix_2_into_1``, ``phi_3D_admix_*_into_*``,
    ``phi_4D_admix_into_1`` ... ``phi_4D_admix_into_4`` and ``phi_5D_admix_into_1`` ...
    ``phi_5D_admix_into_5``. ``target_population`` is the 1-based population that receives the
    migrants; ``proportions`` lists the share each of the other populations contributes, in
    ascending population order (so for 3 populations and target 3, the entries are the shares of
    populations 1 and 2). The density needs 2 to 5 populations.
    """

    dadi, np = _dadi_and_numpy()
    points = _grid_points(grid_points)
    density = _load_phi(np, phi, points, (2, 3, 4, 5))
    count = density.ndim
    target = _integer(target_population, "target_population", minimum=1, maximum=count)
    f = _fractions(proportions, count - 1, "proportions")
    manip = dadi.PhiManip
    grids = [_grid(dadi, points)] * count
    if count == 2:
        function = manip.phi_2D_admix_1_into_2 if target == 2 else manip.phi_2D_admix_2_into_1
    elif count == 3:
        names = {3: "phi_3D_admix_1_and_2_into_3", 2: "phi_3D_admix_1_and_3_into_2",
                 1: "phi_3D_admix_2_and_3_into_1"}
        function = getattr(manip, names[target])
    else:
        function = getattr(manip, f"phi_{count}D_admix_into_{target}")
    result = _call(function, density, *f, *grids)
    return _phi_result(
        np, "phi_admix_into_population", result, points, output_path,
        {"target_population": target, "proportions": f},
    )


def phi_select_populations(
    phi: Sequence[Any] | str | Path,
    *,
    grid_points: int,
    operation: str,
    populations: Sequence[int],
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove, keep or reorder the populations of a density, using ``dadi.PhiManip``.

    ``operation`` is ``remove`` (``remove_pop``: ``populations`` has the one 1-based population
    to drop), ``keep`` (``filter_pops``: the 1-based populations to keep, marginalising the rest)
    or ``reorder`` (``reorder_pops``: the new order as a permutation of 1..n, where entry i is
    the old population that becomes population i).
    """

    dadi, np = _dadi_and_numpy()
    if operation not in _SELECT_OPERATIONS:
        raise ConfigurationError(f"operation must be one of {sorted(_SELECT_OPERATIONS)}.")
    points = _grid_points(grid_points)
    density = _load_phi(np, phi, points, (1, 2, 3, 4, 5))
    count = density.ndim
    if isinstance(populations, (str, bytes)) or not isinstance(populations, Sequence):
        raise ConfigurationError("populations must be a list of 1-based population numbers.")
    chosen = [_integer(item, "population", minimum=1, maximum=count) for item in populations]
    if len(set(chosen)) != len(chosen) or not chosen:
        raise ConfigurationError("populations must be distinct and not empty.")
    xx = _grid(dadi, points)
    manip = dadi.PhiManip
    if operation == "remove":
        if len(chosen) != 1 or count < 2:
            raise ConfigurationError("remove needs one population and a density of 2+ populations.")
        result = _call(manip.remove_pop, density, xx, chosen[0])
    elif operation == "keep":
        result = _call(manip.filter_pops, density, xx, chosen)
    else:
        if sorted(chosen) != list(range(1, count + 1)):
            raise ConfigurationError(f"reorder needs a permutation of 1..{count}.")
        result = _call(manip.reorder_pops, density, chosen)
    return _phi_result(
        np, "phi_select_populations", result, points, output_path,
        {"select": operation, "populations": chosen},
    )


def _spectrum_result(np: Any, operation: str, fs: Any, extra: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = {
        "status": "completed",
        "operation": operation,
        "spectrum": _nested(np, fs),
        "shape": [int(n) for n in fs.shape],
        "sample_sizes": [int(n) for n in fs.sample_sizes],
        "pop_ids": None if fs.pop_ids is None else list(fs.pop_ids),
        "folded": bool(fs.folded),
        "segregating_sites": float(fs.S()),
        "provenance": _provenance("dadi", "dadi.Spectrum; masked entries are None"),
    }
    result.update(extra)
    return result


def _pop_names(value: object, count: int) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != count:
        raise ConfigurationError(f"pop_ids must list {count} population names.")
    names = [str(item) for item in value]
    if len(set(names)) != count or not all(names):
        raise ConfigurationError("pop_ids must be distinct, non-empty names.")
    return names


def spectrum_from_phi(
    phi: Sequence[Any] | str | Path,
    *,
    grid_points: int,
    sample_sizes: Sequence[int],
    inbreeding_coefficients: Sequence[float] | None = None,
    ploidies: Sequence[int] | None = None,
    mask_corners: bool = True,
    pop_ids: Sequence[str] | None = None,
    admix_props: Sequence[Sequence[float]] | None = None,
    het_ascertained: str | None = None,
    force_direct: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample an allele-frequency spectrum from a population density, using ``dadi.Spectrum``.

    Without ``inbreeding_coefficients`` this is ``Spectrum.from_phi`` (``sample_sizes`` gives the
    number of sampled chromosomes per population, ``admix_props`` and ``het_ascertained`` as in
    dadi). With ``inbreeding_coefficients`` (one per population, 0 to 1) and ``ploidies`` (one
    integer per population) it is ``Spectrum.from_phi_inbreeding``, which samples individuals
    with the given inbreeding; ``force_direct`` is then always true as in dadi. ``mask_corners``
    masks the fixed classes. The result has the spectrum as a nested list (masked entries None).
    """

    dadi, np = _dadi_and_numpy()
    points = _grid_points(grid_points)
    density = _load_phi(np, phi, points)
    count = density.ndim
    if isinstance(sample_sizes, (str, bytes)) or not isinstance(sample_sizes, Sequence) or len(sample_sizes) != count:
        raise ConfigurationError(f"sample_sizes must list {count} integers, one per population.")
    ns = [_integer(item, "sample size", minimum=2, maximum=_MAX_SAMPLE) for item in sample_sizes]
    xxs = [_grid(dadi, points)] * count
    names = _pop_names(pop_ids, count)
    if not isinstance(mask_corners, bool) or not isinstance(force_direct, bool):
        raise ConfigurationError("mask_corners and force_direct must be true or false.")
    kwargs: dict[str, Any] = {"mask_corners": mask_corners, "pop_ids": names}
    if admix_props is not None:
        kwargs["admix_props"] = [[float(v) for v in row] for row in admix_props]
    if het_ascertained is not None:
        kwargs["het_ascertained"] = het_ascertained
    if inbreeding_coefficients is None:
        if ploidies is not None:
            raise ConfigurationError("ploidies is only used together with inbreeding_coefficients.")
        fs = _call(dadi.Spectrum.from_phi, density, ns, xxs, force_direct=force_direct, **kwargs)
        mode = "from_phi"
    else:
        if ploidies is None or len(ploidies) != count or len(inbreeding_coefficients) != count:
            raise ConfigurationError(
                f"inbreeding_coefficients and ploidies must each list {count} values."
            )
        fs_in = [_number(item, "inbreeding coefficient") for item in inbreeding_coefficients]
        if any(item > 1 for item in fs_in):
            raise ConfigurationError("inbreeding coefficients must be between 0 and 1.")
        ploidy = [_integer(item, "ploidy", minimum=1, maximum=8) for item in ploidies]
        fs = _call(dadi.Spectrum.from_phi_inbreeding, density, ns, xxs, fs_in, ploidy, **kwargs)
        mode = "from_phi_inbreeding"
    return _spectrum_result(np, "spectrum_from_phi", fs, {"mode": mode})


def scramble_spectrum_pop_ids(
    spectrum: Sequence[Any],
    *,
    mask_corners: bool = True,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomly reassign the sampled individuals among the populations, using ``dadi.Spectrum``.

    ``Spectrum.scramble_pop_ids`` builds the null spectrum in which population labels carry no
    information (a permutation of individuals between populations, keeping the total sample).
    The spectrum must have at least two populations. The draw uses ``seed`` for NumPy's and
    Python's global generators for this call only, then restores them.
    """

    dadi, np = _dadi_and_numpy()
    seed_value = _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    if not isinstance(mask_corners, bool):
        raise ConfigurationError("mask_corners must be true or false.")
    fs = _spectrum(dadi, np, spectrum, "spectrum")
    if fs.ndim < 2:
        raise ConfigurationError("spectrum needs at least two populations.")
    with _seeded(np, seed_value):
        scrambled = _call(fs.scramble_pop_ids, mask_corners=mask_corners)
    return _spectrum_result(
        np, "scramble_spectrum_pop_ids", scrambled, {"seed": seed_value}
    )


def spectrum_from_snp_data(
    path: str | Path,
    *,
    pop_ids: Sequence[str],
    projections: Sequence[int],
    polarized: bool = True,
    mask_corners: bool = True,
    misidentification_path: str | Path | None = None,
    force_positive: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build an allele-frequency spectrum from a dadi SNP file, using ``dadi.Spectrum``.

    The file is read into dadi's data dictionary (``Misc.make_data_dict``) and projected down to
    ``projections`` chromosomes per population in ``pop_ids``. Without
    ``misidentification_path`` this is ``Spectrum.from_data_dict`` (``polarized`` false folds the
    spectrum). With ``misidentification_path``, a text file of ``<context triplet> <outgroup base>
    <probability>`` lines, it is ``Spectrum.from_data_dict_corrected``, which removes ancestral
    misidentification (Hernandez et al. 2007); ``force_positive`` sets small negative corrected
    entries to zero conserving the total, and ``polarized`` is ignored (the correction is always
    polarized).
    """

    dadi, np = _dadi_and_numpy()
    source = input_path(path)
    data = _load_snp(dadi, source)
    available = _populations(data)
    if isinstance(pop_ids, (str, bytes)) or not isinstance(pop_ids, Sequence) or not pop_ids:
        raise ConfigurationError("pop_ids must be a non-empty list of population IDs.")
    names = [str(item) for item in pop_ids]
    missing = [name for name in names if name not in available]
    if missing or len(set(names)) != len(names):
        raise ConfigurationError(f"pop_ids {names} must be distinct and in the file {available}.")
    if isinstance(projections, (str, bytes)) or not isinstance(projections, Sequence) or len(projections) != len(names):
        raise ConfigurationError("projections must list one integer per population.")
    proj = [_integer(item, "projection", minimum=2, maximum=_MAX_SAMPLE) for item in projections]
    if not all(isinstance(flag, bool) for flag in (polarized, mask_corners, force_positive)):
        raise ConfigurationError("polarized, mask_corners and force_positive must be true or false.")
    if misidentification_path is None:
        fs = _call(
            dadi.Spectrum.from_data_dict, data, names, proj,
            mask_corners=mask_corners, polarized=polarized,
        )
        mode = "from_data_dict"
    else:
        fux = input_path(misidentification_path)
        fs = _call(
            dadi.Spectrum.from_data_dict_corrected, data, names, proj, str(fux),
            force_pos=force_positive, mask_corners=mask_corners,
        )
        mode = "from_data_dict_corrected"
    return _spectrum_result(
        np, "spectrum_from_snp_data", fs, {"mode": mode, "source": str(source)}
    )


__all__ = [
    "phi_admix_into_population",
    "phi_create_admixed_population",
    "phi_one_population",
    "phi_select_populations",
    "phi_split_population",
    "scramble_spectrum_pop_ids",
    "spectrum_from_phi",
    "spectrum_from_snp_data",
]

_register_public_verbose(globals(), tuple(__all__))
