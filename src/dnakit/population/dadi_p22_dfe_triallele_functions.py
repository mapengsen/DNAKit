"""dadi wrappers for DFE-integrated spectra and the triallelic density model (pack P22).

``dfe_integrated_spectrum`` integrates cached selection spectra over a distribution of fitness
effects (the ``Cache1D.integrate*`` and ``Cache2D.integrate*`` methods of ``dadi.DFE``). The
triallelic functions cover ``dadi.Triallele``: the demographic models, the density integration and
sampling, the alternate-mutation spectrum, spectrum folding/unfolding/log/masking and the number of
segregating sites.

dadi's triallelic demographic models (``demographics.equilibrium``, ``two_epoch``, ``three_epoch``,
``bottlegrowth``) call ``np.float``, which NumPy 1.24 removed, and the forward integration they rely
on (``integration.advance``) returns NaN or zeros for every input with the installed SciPy and
Cython modules (the C tridiagonal solver receives non-contiguous array views). Only the exact
neutral equilibrium therefore works; ``triallelic_neutral_equilibrium`` runs it with dadi's own
``equilibrium_neutral_exact`` and ``sample``. Ancestral misidentification, projection and optimal
scaling of triallelic spectra fail upstream as well and are not offered.

Spectra are nested lists; an entry that dadi masks is ``None``. dadi and NumPy are imported only
when a function runs.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import _DADI_EXTRA
from dnakit.population.dadi_p22_common import (
    choice,
    dadi_and_numpy,
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

_DFE_KINDS = (
    "negative_1d",
    "point_pos_1d",
    "bivariate_2d",
    "bivariate_point_pos_2d",
    "bivariate_symmetric_point_pos_2d",
)
_DFE_1D_DISTRIBUTIONS = ("gamma", "lognormal", "exponential", "beta")
_DFE_2D_DISTRIBUTIONS = ("biv_lognormal", "biv_ind_gamma")
_DFE_1D_MODELS = (
    "equil",
    "two_epoch_sel",
    "three_epoch_sel",
    "growth_sel",
    "bottlegrowth_1d_sel",
    "IM_pre_sel_single_gamma",
    "IM_sel_single_gamma",
    "split_mig_sel_single_gamma",
    "split_asym_mig_sel_single_gamma",
    "split_delay_mig_sel_single_gamma",
    "bottlegrowth_2d_sel_single_gamma",
    "bottlegrowth_split_sel_single_gamma",
    "bottlegrowth_split_mig_sel_single_gamma",
)
_DFE_2D_MODELS = (
    "split_mig_sel",
    "split_asym_mig_sel",
    "split_delay_mig_sel",
    "IM_sel",
    "IM_pre_sel",
    "bottlegrowth_2d_sel",
    "bottlegrowth_split_sel",
    "bottlegrowth_split_mig_sel",
)
_MAX_DFE_SAMPLE = 40
_MAX_GAMMA_POINTS = 60
_MAX_GRID = 200
_TRIALLELIC_TRANSFORMS = ("mask_infeasible", "unfold", "fold_major", "fold_ancestral", "log")
_MAX_TRI_GRID = 120
_MAX_TRI_SAMPLE = 100


def _grid_triplet(value: object) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != 3:
        raise ConfigurationError("grid_points must list three integer grid sizes.")
    return [integer(item, "grid point", minimum=8, maximum=_MAX_GRID) for item in value]


def dfe_integrated_spectrum(
    kind: str,
    model: str,
    demographic_params: Sequence[float],
    sample_sizes: Sequence[int],
    dfe_params: Sequence[float],
    *,
    distribution: str | None = None,
    theta: float = 1.0,
    n_positive: int = 1,
    grid_points: Sequence[int] = (20, 25, 30),
    gamma_points: int = 20,
    gamma_max: float = 2000.0,
    exterior_int: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected spectrum integrated over a distribution of fitness effects (``dadi.DFE`` caches).

    A selection cache (spectra for a grid of negative gammas plus the point-mass gammas the kind
    needs) is built for the demographic ``model`` of ``dadi.DFE.DemogSelModels`` and integrated
    over the DFE. ``kind`` selects the dadi method: ``negative_1d`` (``Cache1D.integrate``,
    negative-gamma distribution, ``dfe_params`` are the distribution parameters), ``point_pos_1d``
    (``Cache1D.integrate_point_pos``; the last ``2 * n_positive`` values are proportion and gamma
    of each positive point mass), ``bivariate_2d`` (``Cache2D.integrate``, bivariate negative
    distribution), ``bivariate_point_pos_2d`` (``Cache2D.integrate_point_pos``; last four values
    ppos1, gamma1, ppos2, gamma2) and ``bivariate_symmetric_point_pos_2d``
    (``Cache2D.integrate_symmetric_point_pos``; last two values ppos, gamma).

    One-dimensional kinds use single-gamma models (``equil``, ``two_epoch_sel``,
    ``three_epoch_sel``, ``growth_sel``, ``bottlegrowth_1d_sel``, ``split_mig_sel_single_gamma``
    and the other
    ``*_single_gamma`` models) and ``distribution`` ``gamma`` (default), ``lognormal``,
    ``exponential`` or ``beta``. Two-dimensional kinds use two-gamma models (``split_mig_sel``,
    ``IM_sel``, ``split_asym_mig_sel``, ...) and ``biv_lognormal`` (default; mu, sigma, rho) or
    ``biv_ind_gamma``. ``gamma_points`` log-spaced gammas up to ``gamma_max`` are cached;
    ``exterior_int=False`` skips the density outside that range. The result is multiplied by
    ``theta``.
    """

    _dadi, np = dadi_and_numpy()
    name = choice(kind, _DFE_KINDS, "kind")
    one_dim = name.endswith("_1d")
    options = _DFE_1D_MODELS if one_dim else _DFE_2D_MODELS
    model_name = choice(model, options, "model")
    demo = (
        []
        if model_name == "equil" and not demographic_params
        else reals(demographic_params, "demographic_params")
    )
    if isinstance(sample_sizes, (str, bytes)) or not isinstance(sample_sizes, Sequence):
        raise ConfigurationError("sample_sizes must list one or two integer sample sizes.")
    if not 1 <= len(sample_sizes) <= 2:
        raise ConfigurationError("sample_sizes must list one or two integer sample sizes.")
    ns = [integer(item, "sample size", minimum=2, maximum=_MAX_DFE_SAMPLE) for item in sample_sizes]
    values = reals(dfe_params, "dfe_params")
    pts = _grid_triplet(grid_points)
    gammas = integer(gamma_points, "gamma_points", minimum=5, maximum=_MAX_GAMMA_POINTS)
    ceiling = number(gamma_max, "gamma_max", positive=True)
    scale = number(theta, "theta", positive=True)
    exterior = flag(exterior_int, "exterior_int")
    positives = integer(n_positive, "n_positive", minimum=1, maximum=10)
    distributions = _DFE_1D_DISTRIBUTIONS if one_dim else _DFE_2D_DISTRIBUTIONS
    dist_name = choice(
        distribution if distribution is not None else distributions[0],
        distributions,
        "distribution",
    )
    cache1d_mod = require("dadi.DFE.Cache1D_mod", _DADI_EXTRA)
    cache2d_mod = require("dadi.DFE.Cache2D_mod", _DADI_EXTRA)
    selection = require("dadi.DFE.DemogSelModels", _DADI_EXTRA)
    pdfs = require("dadi.DFE.PDFs", _DADI_EXTRA)
    function = getattr(selection, model_name)
    density = getattr(pdfs, dist_name)
    if name == "point_pos_1d":
        if len(values) < 2 * positives:
            raise ConfigurationError(
                "dfe_params need the distribution parameters and 2 * n_positive values."
            )
        point_gammas = values[-2 * positives + 1 :: 2]
    elif name == "bivariate_point_pos_2d":
        if len(values) < 5:
            raise ConfigurationError("dfe_params end with ppos1, gamma1, ppos2, gamma2.")
        point_gammas = [values[-3], values[-1]]
    elif name == "bivariate_symmetric_point_pos_2d":
        if len(values) < 3:
            raise ConfigurationError("dfe_params end with ppos and gamma.")
        point_gammas = [values[-1]]
    else:
        point_gammas = []
    point_gammas = list(dict.fromkeys(point_gammas))  # a repeated gamma must be cached once
    common = {
        "gamma_bounds": (1e-4, ceiling),
        "gamma_pts": gammas,
        "additional_gammas": point_gammas,
        "cpus": 1,
        "verbose": False,
    }

    def compute() -> Any:
        if one_dim:
            cache = cache1d_mod.Cache1D(tuple(demo), list(ns), function, pts, **common)
            if name == "negative_1d":
                return cache.integrate(values, None, density, scale, None, exterior_int=exterior)
            return cache.integrate_point_pos(
                values, None, density, scale, function, positives, None, exterior_int=exterior
            )
        cache = cache2d_mod.Cache2D(tuple(demo), list(ns), function, pts, **common)
        if name == "bivariate_2d":
            return cache.integrate(values, None, density, scale, None, exterior_int=exterior)
        if name == "bivariate_point_pos_2d":
            rho = values[-5] if dist_name == "biv_lognormal" else 0.0
            return cache.integrate_point_pos(values, None, density, scale, rho, None)
        return cache.integrate_symmetric_point_pos(values, None, density, scale, None)

    spectrum = guarded(f"{name} DFE integration", compute)
    return result(
        "dfe_integrated_spectrum",
        "spectrum[i] is the expected count in derived-allele class i; cached selection spectra "
        "weighted by the DFE density",
        {
            "kind": name,
            "model": model_name,
            "distribution": dist_name,
            "demographic_params": demo,
            "sample_sizes": ns,
            "dfe_params": values,
            "theta": scale,
            "grid_points": pts,
            "gamma_points": gammas,
            "gamma_max": ceiling,
            "exterior_int": exterior,
            "spectrum": nested(np, spectrum),
        },
    )


def _tri_spectrum(
    dadi: Any,
    np: Any,
    value: object,
    name: str,
    folded_major: bool,
    folded_ancestral: bool,
    infeasible: bool,
) -> Any:
    from dadi.Triallele import TriSpectrum_mod

    array = masked_array(np, value, name, (2,), nonnegative=False, square=True)
    return TriSpectrum_mod.TriSpectrum(
        np.ma.getdata(array),
        mask=np.ma.getmaskarray(array),
        mask_infeasible=infeasible,
        data_folded_major=folded_major,
        data_folded_ancestral=folded_ancestral,
        check_folding_major=False,
        check_folding_ancestral=False,
    )


def _tri_description(np: Any, fs: Any) -> dict[str, Any]:
    return {
        "sample_size": int(fs.shape[0]) - 1,
        "folded_major": bool(fs.folded_major),
        "folded_ancestral": bool(fs.folded_ancestral),
        "spectrum": nested(np, fs),
    }


def triallelic_neutral_equilibrium(
    grid_points: int,
    sample_size: int | None = None,
    *,
    folded: bool = False,
    phi: Sequence[Sequence[float]] | None = None,
    return_density: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Neutral equilibrium triallelic density and spectrum, using ``dadi.Triallele``.

    The density is the exact neutral solution of Jenkins et al. (2014) on the grid
    ``x = linspace(0, 1, grid_points + 1)`` with both mutation rates equal to 1
    (``integration.equilibrium_neutral_exact``). With ``sample_size`` it is sampled to the expected
    triallelic spectrum of that many chromosomes (``numerics.sample``), which is the
    ``dadi.Triallele.demographics.equilibrium`` spectrum for neutral alleles; ``folded=True``
    folds it by major and minor derived allele. Pass ``phi`` (a square list with one row per
    grid point) to sample your own density instead; ``grid_points`` is then ignored.
    ``return_density=True`` adds the density to the result.

    Only the neutral equilibrium is offered: dadi's forward integration of the triallelic density
    (``integration.advance``, used for selection, other mutation rates and the epoch models) returns
    NaN or zeros for every input with the installed SciPy and Cython modules.
    """

    _dadi, np = dadi_and_numpy()
    fold = flag(folded, "folded")
    include = flag(return_density, "return_density")
    n = (
        None
        if sample_size is None
        else integer(sample_size, "sample_size", minimum=3, maximum=_MAX_TRI_SAMPLE)
    )
    from dadi.Triallele import integration, numerics

    if phi is not None:
        if n is None:
            raise ConfigurationError("sample_size is required when phi is given.")
        density = masked_array(np, phi, "phi", (2,), nonnegative=False, square=True)
        grid = int(density.shape[0]) - 1
        if grid > _MAX_TRI_GRID:
            raise ConfigurationError("phi is larger than the supported grid.")
        phi_array = np.ma.getdata(density)
    else:
        grid = integer(grid_points, "grid_points", minimum=10, maximum=_MAX_TRI_GRID)
        phi_array = integration.equilibrium_neutral_exact(np.linspace(0, 1, grid + 1))
    x = np.linspace(0, 1, grid + 1)
    data: dict[str, Any] = {
        "grid_points": grid,
        "density_source": "given" if phi is not None else "neutral_exact",
        "density_total": scalar(np.sum(phi_array)),
    }
    if include:
        data["phi"] = np.asarray(phi_array, dtype=float).tolist()
    if n is not None:

        def sample() -> Any:
            spectrum = numerics.sample(phi_array, n, x)
            return spectrum.fold_major() if fold else spectrum

        data.update(_tri_description(np, guarded("triallelic density sampling", sample)))
    elif fold:
        raise ConfigurationError("folded needs sample_size.")
    return result(
        "triallelic_neutral_equilibrium",
        "phi[i, j] is the density of sites with derived-allele frequencies x[i], x[j]; "
        "spectrum[i, j] is the expected number of sites with i and j copies of the two derived "
        "alleles among n chromosomes (infeasible entries are None)",
        data,
    )


def alternate_mutation_triallelic_spectrum(
    sample_size: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Exact neutral triallelic spectrum when two mutations enter at one site (Jenkins et al. 2014).

    Both derived alleles arise at the same site (the alternate mutation mechanism of
    ``dadi.Triallele.integration.alt_mut_mech_sample_spectrum``); population size changes do not
    alter this neutral spectrum. ``sample_size`` is the number of sampled chromosomes (at least 3).
    """

    _dadi, np = dadi_and_numpy()
    n = integer(sample_size, "sample_size", minimum=3, maximum=_MAX_TRI_SAMPLE)
    from dadi.Triallele import integration

    spectrum = guarded(
        "alternate mutation spectrum", lambda: integration.alt_mut_mech_sample_spectrum(n)
    )
    return result(
        "alternate_mutation_triallelic_spectrum",
        "spectrum[i, j] is the expected number of triallelic sites with i and j copies of the two "
        "derived alleles; infeasible entries are None",
        {"sample_size": n, "spectrum": nested(np, spectrum)},
    )


def triallelic_spectrum_transform(
    spectrum: Sequence[Sequence[float | None]],
    operation: str,
    *,
    folded_major: bool = False,
    folded_ancestral: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask, unfold, fold or log-transform a triallelic spectrum.

    Uses ``dadi.Triallele.TriSpectrum``.

    ``operation`` is ``mask_infeasible`` (mask entries with a zero count or i + j >= n),
    ``unfold`` (undo the major/minor fold), ``fold_major`` (fold by major and minor derived allele),
    ``fold_ancestral`` (fold when the ancestral state is unknown, tracking the two minor allele
    counts) or ``log`` (natural log of the entries). ``spectrum`` is a square nested list with
    ``None`` for masked entries; ``folded_major`` and ``folded_ancestral`` state its current
    folding.
    """

    dadi, np = dadi_and_numpy()
    name = choice(operation, _TRIALLELIC_TRANSFORMS, "operation")
    fs = _tri_spectrum(
        dadi,
        np,
        spectrum,
        "spectrum",
        flag(folded_major, "folded_major"),
        flag(folded_ancestral, "folded_ancestral"),
        False,
    )
    method = getattr(fs, name)
    changed = guarded(f"triallelic {name}", method)
    if name == "mask_infeasible":
        changed = fs
    return result(
        "triallelic_spectrum_transform",
        "spectrum[i, j] is the count with i and j copies of the two derived alleles; masked "
        "entries are None",
        {"transform": name, **_tri_description(np, changed)},
    )


def triallelic_segregating_sites(
    spectrum: Sequence[Sequence[float | None]],
    *,
    folded_major: bool = False,
    folded_ancestral: bool = False,
    apply_infeasible_mask: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of triallelic segregating sites: the sum of the unmasked entries.

    dadi methods: ``TriSpectrum.S``.

    ``None`` entries of ``spectrum`` are masked. ``apply_infeasible_mask=True`` additionally masks
    the infeasible entries (zero counts, i + j >= n) as dadi does when it builds a spectrum.
    """

    dadi, np = dadi_and_numpy()
    fs = _tri_spectrum(
        dadi,
        np,
        spectrum,
        "spectrum",
        flag(folded_major, "folded_major"),
        flag(folded_ancestral, "folded_ancestral"),
        flag(apply_infeasible_mask, "apply_infeasible_mask"),
    )
    total = guarded("triallelic segregating sites", fs.S)
    return result(
        "triallelic_segregating_sites",
        "sum of the unmasked entries of the triallelic spectrum",
        {"sample_size": int(fs.shape[0]) - 1, "segregating_sites": scalar(total)},
    )


__all__ = [
    "alternate_mutation_triallelic_spectrum",
    "dfe_integrated_spectrum",
    "triallelic_neutral_equilibrium",
    "triallelic_segregating_sites",
    "triallelic_spectrum_transform",
]

_register_public_verbose(globals(), tuple(__all__))
