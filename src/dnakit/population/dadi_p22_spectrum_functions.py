"""dadi wrappers for spectrum operations, statistics, low-coverage models and diffusion (P22).

The functions cover dadi's ``Spectrum`` methods (masking, projection, population selection and
reordering, log, unfolding, resampling, summary statistics, ms and sfs_code parsing), the grid
extrapolation wrappers of ``Numerics``, ``Misc.fragment_data_dict``, the ``Godambe`` derivatives,
the diffusion integrations of ``Integration`` for one to five populations (and X-linked loci), and
the GATK multi-sample low-coverage model of ``dadi.LowPass``.

Spectra are nested lists; an entry that dadi masks is ``None``. Random operations take a ``seed``
and leave NumPy's global random state as they found it. dadi and NumPy are imported only when a
function runs. Output files are never overwritten.
"""

from __future__ import annotations

import contextlib
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path
from dnakit.exceptions import ConfigurationError
from dnakit.population.dadi_functions import (
    _DEMOGRAPHICS_1D,
    _DEMOGRAPHICS_2D,
    _DEMOGRAPHICS_3D,
    _load_snp,
)
from dnakit.population.dadi_p22_common import (
    choice,
    dadi_and_numpy,
    existing_file,
    flag,
    guarded,
    integer,
    masked_array,
    nested,
    number,
    reals,
    result,
    scalar,
)

_RESAMPLE_MODES = ("poisson", "fixed_size")
_SIM_FORMATS = ("ms", "sfs_code")
_SFS_CODE_SITES = ("all", "syn", "nonsyn")
_STATISTICS = (
    "S",
    "Watterson_theta",
    "theta_L",
    "Zengs_E",
    "pi",
    "Tajima_D",
    "Fst",
)
_DERIVATIVES = ("hessian", "gradient")
_MAX_SPECTRUM_SAMPLE = 400
_MAX_PHI_GRID = {1: 400, 2: 150, 3: 60, 4: 30, 5: 20}
_MAX_PHI_RETURN = 50_000
_MAX_SIMULATIONS = 200_000
_MAX_LOWPASS_HAPLOTYPES = 24
_MAX_FILE_BYTES = 200_000_000


@contextlib.contextmanager
def _seeded(np: Any, seed: int | None) -> Any:
    """Seed NumPy's global generator for one call and restore its state afterwards."""

    if seed is None:
        yield
        return
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(state)


def _seed(value: object) -> int | None:
    return None if value is None else integer(value, "seed", minimum=0, maximum=2**32 - 1)


def _spectrum(
    dadi: Any, np: Any, value: object, folded: bool, mask_corners: bool, name: str = "spectrum"
) -> Any:
    array = masked_array(np, value, name, (1, 2, 3), nonnegative=False)
    if max(array.shape) - 1 > _MAX_SPECTRUM_SAMPLE:
        raise ConfigurationError(
            f"{name} is larger than {_MAX_SPECTRUM_SAMPLE} chromosomes per population."
        )
    return dadi.Spectrum(
        array, mask_corners=flag(mask_corners, "mask_corners"), data_folded=flag(folded, "folded")
    )


def _describe(np: Any, fs: Any) -> dict[str, Any]:
    return {
        "sample_sizes": [int(n) for n in fs.sample_sizes],
        "folded": bool(fs.folded),
        "spectrum": nested(np, fs),
    }


_CONVENTION = (
    "spectrum[i, j, ...] is the count of sites with i, j, ... derived alleles in the sampled "
    "chromosomes of each population; masked entries are None"
)


def mask_spectrum_corners(
    spectrum: Sequence[Any],
    *,
    folded: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask the entries seen in none or in all samples of a spectrum (``Spectrum.mask_corners``).

    The two corner entries (no derived allele in any population, derived allele in every sample)
    are masked, as dadi does for every spectrum it builds. Other masked entries stay masked.
    ``Spectrum.unmask_all`` is not offered: it raises ``IndexError`` for every input with NumPy 2.
    """

    dadi, np = dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, folded, False)
    guarded("mask_corners", fs.mask_corners)
    return result("mask_spectrum_corners", _CONVENTION, _describe(np, fs))


def project_spectrum(
    spectrum: Sequence[Any],
    sample_sizes: Sequence[int],
    *,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Project a spectrum to smaller sample sizes by hypergeometric sampling (``Spectrum.project``).

    ``sample_sizes`` has one entry per population and none may exceed the current size.
    """

    dadi, np = dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    if isinstance(sample_sizes, (str, bytes)) or not isinstance(sample_sizes, Sequence):
        raise ConfigurationError("sample_sizes must list one integer per population.")
    target = [integer(item, "sample size", minimum=2) for item in sample_sizes]
    projected = guarded("project", lambda: fs.project(target))
    return result(
        "project_spectrum", _CONVENTION, {"projected_to": target, **_describe(np, projected)}
    )


def subset_spectrum_populations(
    spectrum: Sequence[Any],
    *,
    keep: Sequence[int] | None = None,
    marginalize_over: Sequence[int] | None = None,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep some populations of a joint spectrum and sum over the others.

    dadi methods: ``Spectrum.filter_pops``, ``marginalize``.

    Give exactly one of ``keep`` (populations to keep, numbered from 1, as in ``filter_pops``) or
    ``marginalize_over`` (populations to sum over, also numbered from 1; ``Spectrum.marginalize``
    itself counts axes from 0).
    """

    dadi, np = dadi_and_numpy()
    if (keep is None) == (marginalize_over is None):
        raise ConfigurationError("give exactly one of keep or marginalize_over.")
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    corners = flag(mask_corners, "mask_corners")
    chosen = keep if keep is not None else marginalize_over
    if isinstance(chosen, (str, bytes)) or not isinstance(chosen, Sequence) or not chosen:
        raise ConfigurationError(
            "keep and marginalize_over must be non-empty lists of population numbers."
        )
    numbers = [integer(item, "population number", minimum=1, maximum=fs.ndim) for item in chosen]
    if len(set(numbers)) != len(numbers):
        raise ConfigurationError("population numbers must not repeat.")
    if keep is not None:
        out = guarded("filter_pops", lambda: fs.filter_pops(numbers, mask_corners=corners))
    else:
        out = guarded(
            "marginalize",
            lambda: fs.marginalize([item - 1 for item in numbers], mask_corners=corners),
        )
    return result(
        "subset_spectrum_populations",
        _CONVENTION,
        {
            "keep": keep if keep is None else numbers,
            "marginalize_over": marginalize_over if marginalize_over is None else numbers,
            **_describe(np, out),
        },
    )


def log_spectrum(
    spectrum: Sequence[Any],
    *,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Natural logarithm of the entries of a spectrum (``Spectrum.log``).

    Entries that are not positive become masked.
    """

    dadi, np = dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    out = guarded("log", fs.log)
    return result("log_spectrum", _CONVENTION.replace("count", "log count"), _describe(np, out))


def reorder_spectrum_populations(
    spectrum: Sequence[Any],
    new_order: Sequence[int],
    *,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Permute the populations (axes) of a joint spectrum (``Spectrum.reorder_pops``).

    ``new_order`` lists the original population numbers (from 1) in their new order.
    """

    dadi, np = dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    if isinstance(new_order, (str, bytes)) or not isinstance(new_order, Sequence):
        raise ConfigurationError("new_order must list every population number once.")
    order = [integer(item, "population number", minimum=1, maximum=fs.ndim) for item in new_order]
    out = guarded("reorder_pops", lambda: fs.reorder_pops(order))
    return result(
        "reorder_spectrum_populations", _CONVENTION, {"new_order": order, **_describe(np, out)}
    )


def unfold_spectrum(
    spectrum: Sequence[Any],
    *,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unfold a folded spectrum, assuming each allele is equally likely to be ancestral.

    dadi methods: ``Spectrum.unfold``.

    ``spectrum`` must be a folded spectrum (minor-allele counts); the result is
    ``(spectrum + reversed spectrum) / 2``.
    """

    dadi, np = dadi_and_numpy()
    fs = _spectrum(dadi, np, spectrum, True, mask_corners)
    out = guarded("unfold", fs.unfold)
    return result("unfold_spectrum", _CONVENTION, _describe(np, out))


def resample_spectrum(
    spectrum: Sequence[Any],
    mode: str = "poisson",
    *,
    num_sites: int | None = None,
    only_nonmasked: bool = False,
    seed: int | None = None,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Resample a spectrum, using ``Spectrum.sample`` or ``Spectrum.fixed_size_sample``.

    ``mode`` is ``poisson`` (each entry is replaced by a Poisson draw with the entry as mean) or
    ``fixed_size`` (``num_sites`` sites are drawn from the multinomial with the normalized
    spectrum as probabilities; ``only_nonmasked=True`` draws only into unmasked entries).
    ``seed`` makes the draw reproducible.
    """

    dadi, np = dadi_and_numpy()
    way = choice(mode, _RESAMPLE_MODES, "mode")
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    only = flag(only_nonmasked, "only_nonmasked")
    state = _seed(seed)
    if way == "fixed_size":
        if num_sites is None:
            raise ConfigurationError("num_sites is required for mode fixed_size.")
        sites = integer(num_sites, "num_sites", minimum=1, maximum=1_000_000_000)

        def draw() -> Any:
            return fs.fixed_size_sample(sites, only_nonmasked=only)
    else:

        def draw() -> Any:
            return fs.sample()

    with _seeded(np, state):
        out = guarded("resample", draw)
    return result(
        "resample_spectrum",
        _CONVENTION,
        {"mode": way, "num_sites": num_sites, "seed": state, **_describe(np, out)},
    )


def spectrum_from_simulation_file(
    path: str,
    *,
    file_format: str = "ms",
    average: bool = True,
    mask_corners: bool = True,
    pop_assignments: Sequence[int] | None = None,
    sites: str = "all",
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a spectrum from coalescent or forward simulation output.

    dadi methods: ``Spectrum.from_ms_file``, ``from_sfscode_file``.

    ``file_format`` is ``ms`` (output of ms or msms; ``pop_assignments`` optionally gives the number
    of samples per population, otherwise the ms command line is used) or ``sfs_code`` (``sites``
    selects ``all``, ``syn`` or ``nonsyn`` mutations). ``average=True`` averages over the replicate
    runs in the file, ``False`` sums them. Returns the spectrum and the ms command and seeds when
    present.
    """

    dadi, np = dadi_and_numpy()
    kind = choice(file_format, _SIM_FORMATS, "file_format")
    source = existing_file(path, "path")
    if source.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError("the file is larger than the supported size.")
    mean = flag(average, "average")
    corners = flag(mask_corners, "mask_corners")
    header: Any = None
    if kind == "ms":
        assign = None
        if pop_assignments is not None:
            if isinstance(pop_assignments, (str, bytes)) or not isinstance(
                pop_assignments, Sequence
            ):
                raise ConfigurationError(
                    "pop_assignments must list the number of samples per population."
                )
            assign = [integer(item, "pop_assignments entry", minimum=1) for item in pop_assignments]
        fs, header = guarded(
            "from_ms_file",
            lambda: dadi.Spectrum.from_ms_file(
                str(source),
                average=mean,
                mask_corners=corners,
                return_header=True,
                pop_assignments=assign,
            ),
        )
    else:
        selected = choice(sites, _SFS_CODE_SITES, "sites")
        fs, header = guarded(
            "from_sfscode_file",
            lambda: dadi.Spectrum.from_sfscode_file(
                str(source), sites=selected, average=mean, mask_corners=corners, return_header=True
            ),
        )
    command, seeds = (str(header[0]).strip(), str(header[1]).strip()) if header else ("", "")
    return result(
        "spectrum_from_simulation_file",
        _CONVENTION,
        {
            "file_format": kind,
            "path": str(source),
            "average": mean,
            "command": command,
            "seeds": seeds,
            **_describe(np, fs),
        },
    )


def calculate_sfs_statistic(
    spectrum: Sequence[Any],
    statistic: str,
    *,
    folded: bool = False,
    mask_corners: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Population-genetic statistic of a site frequency spectrum, using dadi's ``Spectrum`` methods.

    ``statistic`` is ``S`` (segregating sites), ``Watterson_theta``, ``theta_L`` (Zeng et al. 2006),
    ``Zengs_E`` (Zeng's E), ``pi`` (pairwise differences with the n / (n - 1) correction),
    ``Tajima_D`` (Gillespie's formula) or ``Fst`` (Weir and Cockerham, joint spectrum). All but
    ``Fst`` and ``S`` need a one-population spectrum; ``S`` works for any dimension.
    """

    dadi, np = dadi_and_numpy()
    name = choice(statistic, _STATISTICS, "statistic")
    fs = _spectrum(dadi, np, spectrum, folded, mask_corners)
    value = guarded(name, getattr(fs, name))
    return result(
        "calculate_sfs_statistic",
        "statistic computed from the unmasked entries of the spectrum, as in dadi.Spectrum",
        {
            "statistic": name,
            "sample_sizes": [int(n) for n in fs.sample_sizes],
            "value": scalar(value),
        },
    )


def fragment_snp_data(
    path: str,
    chunk_size: int,
    *,
    include_snp_ids: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split a dadi SNP file into genomic chunks for block bootstrapping.

    dadi methods: ``Misc.fragment_data_dict``.

    ``chunk_size`` is the chunk length in base pairs. SNP IDs must look like
    ``chromosome_position``. The result gives the number of SNPs in every chunk; with
    ``include_snp_ids=True`` also their IDs.
    """

    dadi, _np = dadi_and_numpy()
    source = input_path(path)
    size = integer(chunk_size, "chunk_size", minimum=1, maximum=10**10)
    include = flag(include_snp_ids, "include_snp_ids")
    data = _load_snp(dadi, source)
    chunks = guarded("fragment_data_dict", lambda: dadi.Misc.fragment_data_dict(data, size))
    entries = []
    for index, chunk in enumerate(chunks):
        row: dict[str, Any] = {"index": index, "num_snps": len(chunk)}
        if include:
            row["snp_ids"] = sorted(chunk)
        entries.append(row)
    return result(
        "fragment_snp_data",
        "a chunk holds the SNPs of one chromosome inside one window of chunk_size base pairs",
        {
            "path": str(source),
            "chunk_size": size,
            "num_snps": len(data),
            "num_chunks": len(entries),
            "chunks": entries,
        },
    )


def expected_sfs_extrapolated(
    model: str,
    params: Sequence[float],
    sample_sizes: Sequence[int],
    *,
    grid_points: Sequence[int] = (40, 50, 60),
    log_extrapolation: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected spectrum extrapolated to an infinitely fine grid.

    dadi methods: ``Numerics.make_extrap_func``, ``make_extrap_log_func``.

    ``model`` is a demographic model of the one-, two- or three-population family chosen by the
    number of ``sample_sizes`` (the names accepted by ``expected_sfs_1d``, ``expected_sfs_2d`` and
    ``expected_sfs_3d``). The model is evaluated on each of the ``grid_points`` grid sizes
    (three or more) and Richardson-extrapolated. ``log_extrapolation=True`` extrapolates the
    logarithm of the spectrum (``make_extrap_log_func``), which suits spectra spanning many orders
    of magnitude but fails if any value is not positive.
    """

    dadi, np = dadi_and_numpy()
    if (
        isinstance(sample_sizes, (str, bytes))
        or not isinstance(sample_sizes, Sequence)
        or not 1 <= len(sample_sizes) <= 3
    ):
        raise ConfigurationError("sample_sizes must list one to three integer sample sizes.")
    ns = tuple(
        integer(item, "sample size", minimum=2, maximum=_MAX_SPECTRUM_SAMPLE)
        for item in sample_sizes
    )
    family, options = {
        1: (dadi.Demographics1D, _DEMOGRAPHICS_1D),
        2: (dadi.Demographics2D, _DEMOGRAPHICS_2D),
        3: (dadi.Demographics3D, _DEMOGRAPHICS_3D),
    }[len(ns)]
    name = choice(model, options, "model")
    function = getattr(family, name, None)
    if function is None:
        raise ConfigurationError(f"Installed dadi does not provide model {name!r}.")
    values = reals(params, "params")
    if (
        isinstance(grid_points, (str, bytes))
        or not isinstance(grid_points, Sequence)
        or len(grid_points) < 3
    ):
        raise ConfigurationError("grid_points must list at least three grid sizes.")
    grids = [integer(item, "grid point", minimum=8, maximum=400) for item in grid_points]
    log = flag(log_extrapolation, "log_extrapolation")
    wrapped = (
        dadi.Numerics.make_extrap_log_func(function)
        if log
        else dadi.Numerics.make_extrap_func(function)
    )
    spectrum = guarded(f"extrapolated {name}", lambda: wrapped(tuple(values), ns, grids))
    return result(
        "expected_sfs_extrapolated",
        _CONVENTION + "; extrapolated from the listed grid sizes",
        {
            "model": name,
            "params": values,
            "grid_points": grids,
            "log_extrapolation": log,
            "sample_sizes": list(ns),
            "spectrum": nested(np, spectrum),
        },
    )


def log_likelihood_derivatives(
    data: Sequence[Any],
    model: str,
    params: Sequence[float],
    *,
    derivative: str = "hessian",
    grid_points: Sequence[int] = (20, 30, 40),
    eps: float = 0.01,
    multinom: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hessian or gradient of a model's log-likelihood with respect to its parameters.

    dadi methods: ``Godambe.get_hess``, ``get_grad``.

    The function differentiated is the log-likelihood of ``data`` (one to three populations) under
    the demographic ``model`` extrapolated over ``grid_points``; ``multinom=True`` uses the
    multinomial log-likelihood (theta optimized, so it is not among ``params``), ``False`` the
    Poisson one. Derivatives are finite differences with fractional step ``eps``. ``derivative`` is
    ``hessian`` (matrix of second derivatives) or ``gradient`` (vector of first derivatives).
    """

    dadi, np = dadi_and_numpy()
    way = choice(derivative, _DERIVATIVES, "derivative")
    observed = _spectrum(dadi, np, data, False, True, "data")
    if observed.ndim > 3:
        raise ConfigurationError("data must have at most three populations.")
    family, options = {
        1: (dadi.Demographics1D, _DEMOGRAPHICS_1D),
        2: (dadi.Demographics2D, _DEMOGRAPHICS_2D),
        3: (dadi.Demographics3D, _DEMOGRAPHICS_3D),
    }[observed.ndim]
    name = choice(model, options, "model")
    function = dadi.Numerics.make_extrap_func(getattr(family, name))
    start = reals(params, "params")
    if (
        isinstance(grid_points, (str, bytes))
        or not isinstance(grid_points, Sequence)
        or len(grid_points) < 3
    ):
        raise ConfigurationError("grid_points must list at least three grid sizes.")
    grids = [integer(item, "grid point", minimum=8, maximum=200) for item in grid_points]
    step = number(eps, "eps", positive=True, maximum=0.5)
    multi = flag(multinom, "multinom")
    sizes = tuple(int(n) for n in observed.sample_sizes)
    score = dadi.Inference.ll_multinom if multi else dadi.Inference.ll

    def objective(p: Any, observed_spectrum: Any) -> float:
        return float(score(function(tuple(p), sizes, grids), observed_spectrum))

    derive = dadi.Godambe.get_hess if way == "hessian" else dadi.Godambe.get_grad
    value = guarded(way, lambda: derive(objective, start, step, args=(observed,)))
    matrix = np.asarray(value, dtype=float)
    return result(
        "log_likelihood_derivatives",
        "finite-difference derivatives of the log-likelihood with respect to the parameters",
        {
            "derivative": way,
            "model": name,
            "params": start,
            "grid_points": grids,
            "eps": step,
            "multinom": multi,
            "value": [[scalar(item) for item in row] for row in matrix.reshape(len(start), -1)]
            if way == "hessian"
            else [scalar(item) for item in matrix.reshape(-1)],
        },
    )


def _grid_for(populations: int, grid_points: object) -> int:
    limit = _MAX_PHI_GRID[populations]
    return integer(grid_points, "grid_points", minimum=8, maximum=limit)


def _migration(populations: int, migration: object) -> dict[str, float]:
    kwargs: dict[str, float] = {}
    if migration is None:
        return kwargs
    if (
        isinstance(migration, (str, bytes))
        or not isinstance(migration, Sequence)
        or len(migration) != populations
        or any(
            isinstance(row, (str, bytes))
            or not isinstance(row, Sequence)
            or len(row) != populations
            for row in migration
        )
    ):
        raise ConfigurationError("migration must be a square matrix with one row per population.")
    for i in range(populations):
        for j in range(populations):
            if i != j:
                kwargs[f"m{i + 1}{j + 1}"] = number(migration[i][j], "migration entry", minimum=0.0)
    return kwargs


def _per_pop(value: object, populations: int, name: str, default: float) -> list[float]:
    if value is None:
        return [default] * populations
    return reals(value, name, length=populations)


def integrate_diffusion_phi(
    populations: int,
    duration: float,
    *,
    nu: Sequence[float] | None = None,
    migration: Sequence[Sequence[float]] | None = None,
    gamma: Sequence[float] | None = None,
    h: Sequence[float] | None = None,
    theta0: float = 1.0,
    grid_points: int = 30,
    x_linked: bool = False,
    breeding_ratio: float = 1.0,
    mutation_ratio: float = 1.0,
    sample_sizes: Sequence[int] | None = None,
    return_phi: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Integrate the diffusion density phi over time for one to five populations.

    Uses ``dadi.Integration``.

    The start is the equilibrium density of an ancestral population (``PhiManip.phi_1D`` or
    ``phi_1D_X``), split at time zero into ``populations`` populations (population k + 1 splits
    from population k). The density is then integrated for ``duration`` (units of 2N generations)
    with constant relative sizes ``nu``, selection ``gamma``, dominance ``h`` (one value per
    population, defaults 1, 0 and 0.5) and the migration matrix ``migration[i][j]``, the rate
    into population i+1 from j+1. ``populations`` 1 to 5 select ``one_pop``, ``two_pops``,
    ``three_pops``, ``four_pops`` and ``five_pops``; ``x_linked=True`` (one population only) uses
    ``one_pop_X`` with ``breeding_ratio`` (beta) and ``mutation_ratio`` (alpha). With
    ``sample_sizes`` the expected spectrum is computed from phi (``Spectrum.from_phi``).
    ``return_phi=True`` returns phi itself when it has at most 50,000 entries.
    """

    dadi, np = dadi_and_numpy()
    count = integer(populations, "populations", minimum=1, maximum=5)
    time = number(duration, "duration", minimum=0.0)
    grid = _grid_for(count, grid_points)
    sizes_nu = _per_pop(nu, count, "nu", 1.0)
    if any(item <= 0 for item in sizes_nu):
        raise ConfigurationError("nu entries must be greater than zero.")
    gammas = _per_pop(gamma, count, "gamma", 0.0)
    dominance = _per_pop(h, count, "h", 0.5)
    mutation = number(theta0, "theta0", positive=True)
    beta = number(breeding_ratio, "breeding_ratio", positive=True)
    alpha = number(mutation_ratio, "mutation_ratio", positive=True)
    linked = flag(x_linked, "x_linked")
    if linked and count != 1:
        raise ConfigurationError("x_linked integration is available for one population only.")
    flows = _migration(count, migration)
    include = flag(return_phi, "return_phi")
    xx = dadi.Numerics.default_grid(grid)

    def build() -> Any:
        if linked:
            phi = dadi.PhiManip.phi_1D_X(
                xx,
                nu=sizes_nu[0],
                theta0=mutation,
                gamma=gammas[0],
                h=dominance[0],
                beta=beta,
                alpha=alpha,
            )
            return dadi.Integration.one_pop_X(
                phi,
                xx,
                time,
                nu=sizes_nu[0],
                gamma=gammas[0],
                h=dominance[0],
                beta=beta,
                alpha=alpha,
                theta0=mutation,
            )
        phi = dadi.PhiManip.phi_1D(xx, theta0=mutation, gamma=gammas[0], h=dominance[0])
        if count >= 2:
            phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
        if count >= 3:
            phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
        if count >= 4:
            phi = dadi.PhiManip.phi_3D_to_4D(phi, 0, 0, xx, xx, xx, xx)
        if count >= 5:
            phi = dadi.PhiManip.phi_4D_to_5D(phi, 0, 0, 0, xx, xx, xx, xx, xx)
        kwargs: dict[str, Any] = {"theta0": mutation}
        if count == 1:
            return dadi.Integration.one_pop(
                phi, xx, time, nu=sizes_nu[0], gamma=gammas[0], h=dominance[0], theta0=mutation
            )
        for index in range(count):
            kwargs[f"nu{index + 1}"] = sizes_nu[index]
            kwargs[f"gamma{index + 1}"] = gammas[index]
            kwargs[f"h{index + 1}"] = dominance[index]
        kwargs.update(flows)
        integrate = {
            2: dadi.Integration.two_pops,
            3: dadi.Integration.three_pops,
            4: dadi.Integration.four_pops,
            5: dadi.Integration.five_pops,
        }[count]
        return integrate(phi, xx, time, **kwargs)

    phi = guarded("diffusion integration", build)
    data: dict[str, Any] = {
        "populations": count,
        "x_linked": linked,
        "duration": time,
        "grid_points": grid,
        "phi_shape": [int(n) for n in phi.shape],
        "phi_total": scalar(np.sum(phi)),
        "grid": [float(value) for value in xx],
    }
    if include:
        if phi.size > _MAX_PHI_RETURN:
            raise ConfigurationError(
                f"phi has {phi.size} entries; return_phi is limited to {_MAX_PHI_RETURN}. "
                "Use sample_sizes instead."
            )
        data["phi"] = np.asarray(phi, dtype=float).tolist()
    if sample_sizes is not None:
        if (
            isinstance(sample_sizes, (str, bytes))
            or not isinstance(sample_sizes, Sequence)
            or len(sample_sizes) != count
        ):
            raise ConfigurationError("sample_sizes must list one sample size per population.")
        ns = [
            integer(item, "sample size", minimum=2, maximum=_MAX_SPECTRUM_SAMPLE)
            for item in sample_sizes
        ]
        spectrum = guarded(
            "Spectrum.from_phi", lambda: dadi.Spectrum.from_phi(phi, ns, (xx,) * count)
        )
        data["sample_sizes"] = ns
        data["spectrum"] = nested(np, spectrum)
    return result(
        "integrate_diffusion_phi",
        "phi is the density of sites over derived-allele frequencies on dadi's default grid; "
        "time is in units of 2 Nref generations",
        data,
    )


def _coverage(np: Any, value: object, name: str = "coverage_distribution") -> dict[str, Any]:
    """Map each population to dadi's ``[depths, probabilities]`` array."""

    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError(f"{name} must map each population to its depth probabilities.")
    out: dict[str, Any] = {}
    for pop, probabilities in value.items():
        values = reals(probabilities, f"{name}[{pop}]", minimum_length=2)
        if any(item < 0 for item in values) or abs(sum(values) - 1.0) > 1e-6:
            raise ConfigurationError(f"{name}[{pop}] must be probabilities that sum to 1.")
        out[str(pop)] = np.array([np.arange(len(values)), values], dtype=float)
    return out


def _haplotypes(value: object, name: str) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must list an even number of haplotypes per population.")
    out = [
        integer(item, f"{name} entry", minimum=2, maximum=_MAX_LOWPASS_HAPLOTYPES) for item in value
    ]
    if any(item % 2 for item in out):
        raise ConfigurationError(f"{name} entries must be even.")
    return out


@contextlib.contextmanager
def _lowpass_rng(np: Any, module: Any, seed: int | None) -> Any:
    """Seed both NumPy's global state and the module-level generator of dadi.LowPass."""

    if seed is None:
        yield
        return
    previous = module.rng
    module.rng = np.random.default_rng(seed)
    try:
        with _seeded(np, seed):
            yield
    finally:
        module.rng = previous


def simulate_low_coverage_reads(
    coverage_distribution: Mapping[str, Sequence[float]],
    genotype_partition: Sequence[int],
    sequenced_haplotypes: Sequence[int],
    number_simulations: int,
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simulate read counts across a genotype partition (``dadi.LowPass.simulate_reads``).

    ``genotype_partition`` lists the genotype of every sequenced individual (0 homozygous
    reference, 1 heterozygote, 2 homozygous alternate), population after population;
    ``sequenced_haplotypes`` gives the haplotypes per population (two per individual).
    ``coverage_distribution`` maps each population to the probabilities of read depths 0, 1, 2, ...
    Each simulated locus gets reference and alternate read counts per individual.
    """

    _dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    cov = _coverage(np, coverage_distribution)
    haplotypes = _haplotypes(sequenced_haplotypes, "sequenced_haplotypes")
    if len(cov) != len(haplotypes):
        raise ConfigurationError("coverage_distribution needs one entry per population.")
    individuals = [item // 2 for item in haplotypes]
    if isinstance(genotype_partition, (str, bytes)) or not isinstance(genotype_partition, Sequence):
        raise ConfigurationError("genotype_partition must list a genotype per individual.")
    partition = [integer(item, "genotype", minimum=0, maximum=2) for item in genotype_partition]
    if len(partition) != sum(individuals):
        raise ConfigurationError(
            "genotype_partition must list sequenced_haplotypes / 2 genotypes per population."
        )
    count = integer(number_simulations, "number_simulations", minimum=1, maximum=_MAX_SIMULATIONS)
    state = _seed(seed)
    with _lowpass_rng(np, LowPass, state):
        n_ref, n_alt = guarded(
            "simulate_reads", lambda: LowPass.simulate_reads(cov, partition, individuals, count)
        )
    return result(
        "simulate_low_coverage_reads",
        "rows are simulated loci and columns individuals; counts are reads supporting each allele",
        {
            "number_simulations": count,
            "seed": state,
            "reference_reads": n_ref.tolist(),
            "alternate_reads": n_alt.tolist(),
        },
    )


def subsample_genotypes(
    genotype_calls: Sequence[Sequence[int]],
    subsample_haplotypes: int,
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Subsample genotype calls to fewer haplotypes, skipping loci with too few calls.

    dadi methods: ``LowPass.subsample_genotypes_1D``.

    ``genotype_calls`` has one row per locus and one column per individual, with genotypes 0, 1, 2
    and 99 for a missing call. For ``subsample_haplotypes`` haplotypes, ``subsample_haplotypes / 2``
    called individuals are drawn at random from every locus that has at least that many calls.
    """

    _dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    if not isinstance(genotype_calls, Sequence) or not genotype_calls:
        raise ConfigurationError("genotype_calls must be a non-empty list of rows.")
    try:
        array = np.array(genotype_calls, dtype=int)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("genotype_calls must be a rectangular list of integers.") from exc
    if array.ndim != 2 or not np.isin(array, (0, 1, 2, 99)).all():
        raise ConfigurationError("genotype_calls must be a 2-D list with genotypes 0, 1, 2 or 99.")
    target = integer(
        subsample_haplotypes, "subsample_haplotypes", minimum=2, maximum=2 * array.shape[1]
    )
    if target % 2:
        raise ConfigurationError("subsample_haplotypes must be even.")
    state = _seed(seed)
    with _lowpass_rng(np, LowPass, state):
        out = guarded(
            "subsample_genotypes_1D", lambda: LowPass.subsample_genotypes_1D(array, target)
        )
    return result(
        "subsample_genotypes",
        "rows are the loci that kept enough calls; columns are subsampled individuals",
        {
            "subsample_haplotypes": target,
            "seed": state,
            "genotypes": out.tolist(),
            "num_loci": int(out.shape[0]),
        },
    )


def simulate_gatk_calling_spectrum(
    coverage_distribution: Mapping[str, Sequence[float]],
    allele_counts: Sequence[int],
    sequenced_haplotypes: Sequence[int],
    subsample_haplotypes: Sequence[int],
    number_simulations: int,
    *,
    inbreeding: Sequence[float] | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Spectrum produced by GATK multi-sample calling of an allele at given frequencies.

    dadi methods: ``LowPass.simulate_GATK_multisample_calling``.

    ``allele_counts`` is the true number of derived haplotypes in each population's sequenced
    sample. Loci are simulated with the read depths of ``coverage_distribution``, called as GATK
    would (variants need at least two alternate reads), and subsampled to
    ``subsample_haplotypes``. The result is the probability of each observed derived-allele count
    (0 includes the loci that were not called). ``inbreeding`` gives the coefficient per
    population (default 0).
    """

    _dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    cov = _coverage(np, coverage_distribution)
    sequenced = _haplotypes(sequenced_haplotypes, "sequenced_haplotypes")
    kept = _haplotypes(subsample_haplotypes, "subsample_haplotypes")
    if len(cov) != len(sequenced) or len(kept) != len(sequenced):
        raise ConfigurationError(
            "coverage_distribution and the haplotype lists need one entry per population."
        )
    if any(k > s for k, s in zip(kept, sequenced, strict=True)):
        raise ConfigurationError("subsample_haplotypes cannot exceed sequenced_haplotypes.")
    if (
        isinstance(allele_counts, (str, bytes))
        or not isinstance(allele_counts, Sequence)
        or len(allele_counts) != len(sequenced)
    ):
        raise ConfigurationError("allele_counts needs one count per population.")
    counts = [
        integer(item, "allele count", minimum=0, maximum=s)
        for item, s in zip(allele_counts, sequenced, strict=True)
    ]
    coefficients = _per_pop(inbreeding, len(sequenced), "inbreeding", 0.0)
    if any(not 0 <= item < 1 for item in coefficients):
        raise ConfigurationError("inbreeding coefficients must be in [0, 1).")
    sims = integer(number_simulations, "number_simulations", minimum=1, maximum=_MAX_SIMULATIONS)
    state = _seed(seed)
    with _lowpass_rng(np, LowPass, state):
        out = guarded(
            "simulate_GATK_multisample_calling",
            lambda: LowPass.simulate_GATK_multisample_calling(
                cov, counts, sequenced, kept, sims, coefficients
            ),
        )
    return result(
        "simulate_gatk_calling_spectrum",
        "entry i, j, ... is the probability of observing i, j, ... derived alleles after "
        "calling and subsampling",
        {
            "allele_counts": counts,
            "sequenced_haplotypes": sequenced,
            "subsample_haplotypes": kept,
            "number_simulations": sims,
            "seed": state,
            "probabilities": nested(np, out),
        },
    )


def gatk_no_call_probability(
    coverage_distribution: Sequence[float],
    sequenced_haplotypes: int,
    *,
    inbreeding: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability that GATK multi-sample calling makes no genotype call, for every allele count.

    dadi methods: ``LowPass.probability_of_no_call_1D_GATK_multisample``.

    ``coverage_distribution`` lists the probabilities of read depths 0, 1, 2, ... for one
    population of ``sequenced_haplotypes`` haplotypes. Entry k of the result is the probability
    for a site with k derived haplotypes (k = 0 .. n).
    """

    _dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    cov = _coverage(np, {"pop": coverage_distribution})["pop"]
    n = _haplotypes([sequenced_haplotypes], "sequenced_haplotypes")[0]
    coefficient = number(inbreeding, "inbreeding", minimum=0.0, maximum=0.999999)
    out = guarded(
        "no-call probability",
        lambda: LowPass.probability_of_no_call_1D_GATK_multisample(cov, n, coefficient),
    )
    return result(
        "gatk_no_call_probability",
        "entry k is the probability of no genotype call at a site with k derived haplotypes",
        {
            "sequenced_haplotypes": n,
            "inbreeding": coefficient,
            "probabilities": [scalar(item) for item in out],
        },
    )


def probability_enough_individuals_covered(
    coverage_distribution: Sequence[float],
    sequenced_haplotypes: int,
    subsample_haplotypes: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Probability that enough individuals are covered to obtain the subsample.

    dadi methods: ``LowPass.probability_enough_individuals_covered``.

    For a site already called as variant, the chance that at least ``subsample_haplotypes / 2``
    of the ``sequenced_haplotypes / 2`` individuals have one or more reads.
    """

    _dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    cov = _coverage(np, {"pop": coverage_distribution})["pop"]
    n = _haplotypes([sequenced_haplotypes], "sequenced_haplotypes")[0]
    kept = _haplotypes([subsample_haplotypes], "subsample_haplotypes")[0]
    if kept > n:
        raise ConfigurationError("subsample_haplotypes cannot exceed sequenced_haplotypes.")
    value = guarded(
        "probability_enough_individuals_covered",
        lambda: LowPass.probability_enough_individuals_covered(cov, n, kept),
    )
    return result(
        "probability_enough_individuals_covered",
        "probability that the number of covered individuals reaches subsample_haplotypes / 2",
        {"sequenced_haplotypes": n, "subsample_haplotypes": kept, "probability": scalar(value)},
    )


def expected_sfs_low_pass(
    model: str,
    params: Sequence[float],
    coverage_distribution: Mapping[str, Sequence[float]],
    sequenced_haplotypes: Sequence[int],
    subsample_haplotypes: Sequence[int],
    *,
    grid_points: Sequence[int] = (30, 40, 50),
    inbreeding: Sequence[float] | None = None,
    simulation_threshold: float = 1e-2,
    simulations_per_frequency: int = 1000,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected spectrum of a demographic model distorted by low-coverage GATK multi-sample calling.

    dadi methods: ``LowPass.make_low_pass_func_GATK_multisample``.

    The model (a one-, two- or three-population name of the family given by the number of
    populations) is computed for ``sequenced_haplotypes`` chromosomes, extrapolated over
    ``grid_points``, and transformed to ``subsample_haplotypes`` chromosomes by the probability of
    no call, genotype-calling errors and projection. Entries below ``simulation_threshold`` of the
    analytic transformation are replaced by simulations (``simulations_per_frequency`` loci per
    frequency; 0 threshold always simulates, 1 never).
    """

    dadi, np = dadi_and_numpy()
    from dadi.LowPass import LowPass

    cov = _coverage(np, coverage_distribution)
    pop_ids = list(cov)
    sequenced = _haplotypes(sequenced_haplotypes, "sequenced_haplotypes")
    kept = _haplotypes(subsample_haplotypes, "subsample_haplotypes")
    count = len(sequenced)
    if count > 3 or len(kept) != count or len(cov) != count:
        raise ConfigurationError(
            "give one to three populations, with one entry per population in every list."
        )
    if any(k > s for k, s in zip(kept, sequenced, strict=True)):
        raise ConfigurationError("subsample_haplotypes cannot exceed sequenced_haplotypes.")
    family, options = {
        1: (dadi.Demographics1D, _DEMOGRAPHICS_1D),
        2: (dadi.Demographics2D, _DEMOGRAPHICS_2D),
        3: (dadi.Demographics3D, _DEMOGRAPHICS_3D),
    }[count]
    name = choice(model, options, "model")
    values = reals(params, "params")
    if (
        isinstance(grid_points, (str, bytes))
        or not isinstance(grid_points, Sequence)
        or len(grid_points) < 3
    ):
        raise ConfigurationError("grid_points must list at least three grid sizes.")
    grids = [integer(item, "grid point", minimum=8, maximum=200) for item in grid_points]
    coefficients = _per_pop(inbreeding, count, "inbreeding", 0.0)
    if any(not 0 <= item < 1 for item in coefficients):
        raise ConfigurationError("inbreeding coefficients must be in [0, 1).")
    threshold = number(simulation_threshold, "simulation_threshold", minimum=0.0, maximum=1.0)
    sims = integer(
        simulations_per_frequency, "simulations_per_frequency", minimum=1, maximum=100_000
    )
    state = _seed(seed)
    function = dadi.Numerics.make_extrap_func(getattr(family, name))

    def compute() -> Any:
        wrapped = LowPass.make_low_pass_func_GATK_multisample(
            function,
            cov,
            pop_ids,
            sequenced,
            kept,
            sim_threshold=threshold,
            Fx=coefficients,
            nsim=sims,
        )
        return wrapped(tuple(values), None, grids)

    with _lowpass_rng(np, LowPass, state):
        spectrum = guarded("low-pass model", compute)
    return result(
        "expected_sfs_low_pass",
        _CONVENTION + "; after low-coverage calling and subsampling",
        {
            "model": name,
            "params": values,
            "sequenced_haplotypes": sequenced,
            "subsample_haplotypes": kept,
            "grid_points": grids,
            "seed": state,
            "sample_sizes": [int(n) for n in spectrum.sample_sizes],
            "spectrum": nested(np, spectrum),
        },
    )


__all__ = [
    "calculate_sfs_statistic",
    "expected_sfs_extrapolated",
    "expected_sfs_low_pass",
    "fragment_snp_data",
    "gatk_no_call_probability",
    "integrate_diffusion_phi",
    "log_likelihood_derivatives",
    "log_spectrum",
    "mask_spectrum_corners",
    "probability_enough_individuals_covered",
    "project_spectrum",
    "reorder_spectrum_populations",
    "resample_spectrum",
    "simulate_gatk_calling_spectrum",
    "simulate_low_coverage_reads",
    "spectrum_from_simulation_file",
    "subsample_genotypes",
    "subset_spectrum_populations",
    "unfold_spectrum",
]

_register_public_verbose(globals(), tuple(__all__))
